"""Block checks of `scripts/fetch_bgzf.py` on a small BGZF file with one corrupt block."""

import importlib.util
import struct
import zlib
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent.parent / "scripts" / "fetch_bgzf.py"


def _module():
    spec = importlib.util.spec_from_file_location("fetch_bgzf", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def block(payload: bytes) -> bytes:
    """One BGZF block, as htslib writes it."""
    c = zlib.compressobj(6, zlib.DEFLATED, -15)
    body = c.compress(payload) + c.flush()
    size = 18 + len(body) + 8
    head = (b"\x1f\x8b\x08\x04" + b"\0" * 4 + b"\0\xff" + b"\x06\0" + b"BC" + b"\x02\0"
            + struct.pack("<H", size - 1))
    return head + body + struct.pack("<II", zlib.crc32(payload), len(payload))


def test_scan_finds_the_corrupt_block_and_only_it(tmp_path):
    m = _module()
    a, b, eof = block(b"ACGT" * 500), block(b"chr8\t1\t.\tA\tG\n" * 300), block(b"")
    path = tmp_path / "x.vcf.gz"
    path.write_bytes(a + b + eof)
    assert m.scan(path, path.stat().st_size) == []
    data = bytearray(a + b + eof)
    data[len(a) + 30] ^= 0x40  # one bit inside the second block's deflate data
    path.write_bytes(bytes(data))
    assert m.scan(path, len(data)) == [(len(a), len(a) + len(b))]


def test_scan_covers_a_broken_header_from_the_last_good_block(tmp_path):
    m = _module()
    a, b, eof = block(b"x" * 3000), block(b"y" * 3000), block(b"")
    data = bytearray(a + b + eof)
    data[len(a)] = 0  # the second block's magic
    path = tmp_path / "y.vcf.gz"
    path.write_bytes(bytes(data))
    bad = m.scan(path, len(data))
    assert len(bad) == 1 and bad[0][0] == 0 and bad[0][1] >= len(a) + len(b)


def test_scan_flags_a_file_cut_short(tmp_path):
    m = _module()
    a, b = block(b"x" * 3000), block(b"y" * 3000)
    path = tmp_path / "z.vcf.gz"
    path.write_bytes(a + b[:20])
    assert m.scan(path, len(a) + 20)[-1][1] == len(a) + 20


def test_chunks_overlapping_a_range():
    m = _module()
    assert list(m.chunks(0, 10, chunk=4)) == [0, 1, 2]
    assert list(m.chunks(4, 8, chunk=4)) == [1]
    assert list(m.chunks(7, 9, chunk=4)) == [1, 2]
