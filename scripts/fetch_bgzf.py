"""Fetch BGZF files over HTTP(S) in byte ranges, checking every block, until each matches its md5.

Written for the 1000 Genomes panel of W1 amendment 1 (`docs/design/w1-ancestry-plan.md`). On 15
October whole-file downloads from EBI came back at the right size with a few corrupt BGZF
blocks (chr8: 8 of 828,019, each 2 to 4 kB, scattered), and many transfers stalled at 0 B/s.
Here each file is assembled in place from 8 MB ranges, fetched by up to `--streams` `curl`
processes at once, each aborted and retried when it stalls. Every block's CRC32 and length are
then checked; the ranges holding a bad block are fetched again, and the file is checked again.
A file is renamed to its final name only when its md5 matches the manifest. A local copy of
the file, complete or partial, seeds the work file, so its good blocks are kept.

Usage: `python scripts/fetch_bgzf.py --base URL --manifest MANIFEST --out DIR NAME...`
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import struct
import subprocess
import sys
import threading
import time
import zlib
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from pathlib import Path

CHUNK = 8 << 20
ROUNDS = 10


def scan(path: Path, size: int) -> list:
    """Byte ranges `(start, stop)` of BGZF blocks of `path` that fail their CRC32 or length,
    over the first `size` bytes. Where a header does not parse, the range runs from the last
    good block's start to one chunk past the failure, since the bad bytes may be in that
    block's size field."""
    bad, offset, last = [], 0, 0
    with open(path, "rb") as h:
        while offset < size:
            h.seek(offset)
            head = h.read(18)
            if len(head) < 18 or head[:4] != b"\x1f\x8b\x08\x04" or head[12:14] != b"BC":
                return bad + [(last, min(size, offset + CHUNK))]
            length = struct.unpack("<H", head[16:18])[0] + 1
            body = h.read(length - 18)
            ok = len(body) == length - 18
            if ok:
                try:
                    data = zlib.decompress(body[:-8], -15)
                    crc, isize = struct.unpack("<II", body[-8:])
                    ok = zlib.crc32(data) == crc and len(data) == isize
                except zlib.error:
                    ok = False
            if not ok:
                bad.append((offset, offset + length))
            last, offset = offset, offset + length
    return bad if offset == size else bad + [(last, size)]


def chunks(start: int, stop: int, chunk: int = CHUNK) -> range:
    """Indices of the chunks that overlap the byte range `[start, stop)`."""
    return range(start // chunk, (stop - 1) // chunk + 1)


def md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 24), b""):
            h.update(block)
    return h.hexdigest()


class Fetch:
    """One file: its work file, the chunks written so far, and their fetching."""

    def __init__(self, name: str, url: str, want: str, out: Path, streams: threading.Semaphore,
                 scans: ProcessPoolExecutor):
        self.name, self.url, self.want, self.out, self.streams = name, url, want, out, streams
        self.scans = scans  # block checks run in processes: one Python thread cannot do 13 files
        self.final = out / name
        self.work = out / f"{name}.work"
        self.state = out / f"{name}.state.json"
        self.lock = threading.Lock()
        self.fetches: dict = {}  # chunk -> md5s of the bytes fetched for it, in order

    def log(self, text: str) -> None:
        print(f"{time.strftime('%H:%M:%S')} {self.name}: {text}", flush=True)

    def size(self) -> int:
        out = subprocess.run(["curl", "-sSfI", "--max-time", "60", self.url],
                             capture_output=True, text=True, check=True).stdout
        return int(next(line.split()[1] for line in out.splitlines()
                        if line.lower().startswith("content-length")))

    def prepare(self, size: int) -> set:
        """The set of chunks already in the work file, seeding it from a local copy the first
        time."""
        if self.state.exists() and self.work.exists():
            return set(json.loads(self.state.read_text())["done"])
        seeds = [p for p in (self.final, self.out / f"{self.name}.part") if p.exists()]
        seed = max(seeds, key=lambda p: p.stat().st_size) if seeds else None
        with open(self.work, "wb") as w:
            if seed is not None:
                with open(seed, "rb") as s:
                    shutil.copyfileobj(s, w, 1 << 24)
            w.truncate(size)
        have = min(size, seed.stat().st_size) if seed is not None else 0
        done = set(range(have // CHUNK)) | ({size // CHUNK} if have == size else set())
        self.save(done)
        self.log(f"size {size}, seeded {have} bytes from {seed.name if seed else 'nothing'}")
        return done

    def save(self, done: set) -> None:
        tmp = self.state.with_suffix(".tmp")
        tmp.write_text(json.dumps({"done": sorted(done)}))
        tmp.replace(self.state)

    def fetch_chunk(self, k: int, size: int, done: set) -> None:
        start, stop = k * CHUNK, min(size, (k + 1) * CHUNK)
        tmp = self.out / f"{self.name}.{k}.tmp"
        for attempt in range(12):
            with self.streams:
                r = subprocess.run(["curl", "-sS", "-f", "-r", f"{start}-{stop - 1}",
                                    "--speed-limit", "10000", "--speed-time", "60",
                                    "--max-time", "1800", "-o", str(tmp), self.url],
                                   capture_output=True, text=True)
            got = tmp.stat().st_size if tmp.exists() else -1
            if r.returncode == 0 and got == stop - start:
                data = tmp.read_bytes()
                tmp.unlink()
                fd = os.open(self.work, os.O_WRONLY)
                try:
                    os.pwrite(fd, data, start)
                finally:
                    os.close(fd)
                with self.lock:
                    self.fetches.setdefault(k, []).append(hashlib.md5(data).hexdigest())
                    done.add(k)
                    self.save(done)
                return
            if got > stop - start:
                raise RuntimeError(f"{self.name}: the server ignored the range request")
            tmp.unlink(missing_ok=True)
            time.sleep(min(300, 15 * 2 ** attempt))
        raise RuntimeError(f"{self.name}: chunk {k} failed 12 times")

    def run(self) -> bool:
        try:
            return self._run()
        except Exception as e:  # noqa: BLE001  one file's failure must not stop the others
            self.log(f"failed: {type(e).__name__}: {e}")
            return False

    def _run(self) -> bool:
        if self.final.exists() and md5(self.final) == self.want:
            self.log("already verified")
            return True
        size = self.size()
        done = self.prepare(size)
        total = (size + CHUNK - 1) // CHUNK
        for round_ in range(ROUNDS):
            todo = [k for k in range(total) if k not in done]
            if todo:
                self.log(f"round {round_ + 1}: fetching {len(todo)} of {total} chunks")
                with ThreadPoolExecutor(min(32, len(todo))) as pool:
                    list(pool.map(lambda k: self.fetch_chunk(k, size, done), todo))
            bad = self.scans.submit(scan, self.work, size).result()
            if not bad:
                break
            again = sorted({k for a, b in bad for k in chunks(a, b)})
            same = [k for k in again if len(self.fetches.get(k, [])) >= 2
                    and len(set(self.fetches[k][-2:])) == 1]
            self.log(f"round {round_ + 1}: {len(bad)} bad blocks in {len(again)} chunks"
                     + (f"; chunks fetched twice with the same bytes: {same}" if same else ""))
            with self.lock:
                done.difference_update(again)
                self.save(done)
        else:
            self.log(f"still bad after {ROUNDS} rounds")
            return False
        got = self.scans.submit(md5, self.work).result()
        if got != self.want:
            self.log(f"every block passes but md5 {got} is not the manifest's {self.want}")
            return False
        self.work.replace(self.final)
        self.state.unlink()
        (self.out / f"{self.name}.part").unlink(missing_ok=True)
        refetched = sum(len(v) - 1 for v in self.fetches.values())
        self.log(f"md5 matches; {refetched} chunks were fetched more than once")
        return True


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", required=True)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--streams", type=int, default=16)
    parser.add_argument("names", nargs="+")
    args = parser.parse_args()
    want = dict(line.split()[:2] for line in args.manifest.read_text().splitlines() if line.strip())
    streams = threading.Semaphore(args.streams)
    with ProcessPoolExecutor(8) as scans:
        jobs = [Fetch(n, f"{args.base}/{n}", want[n], args.out, streams, scans)
                for n in args.names]
        with ThreadPoolExecutor(len(jobs)) as pool:
            ok = list(pool.map(Fetch.run, jobs))
    failed = [j.name for j, good in zip(jobs, ok) if not good]
    print(f"{len(jobs) - len(failed)} of {len(jobs)} files verified"
          + (f"; failed: {' '.join(failed)}" if failed else ""), flush=True)
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
