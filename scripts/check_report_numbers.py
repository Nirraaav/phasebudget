"""Check that every headline number in the report still appears in its run file.

Reads docs/presentation/numbers.tsv. For each row it opens run_file and looks
for evidence_substring on the stated line. If the text has moved to another
line it prints "moved" with the new line; if it is gone it prints "MISSING".
Rows whose run_file is NOT FOUND (numbers computed from others, or recorded
only in WORKLOG.md) are listed as "skip".

Run from the repository root:

    python3 scripts/check_report_numbers.py

Exits 1 if any row is MISSING, else 0. Standard library only.
"""
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TABLE = ROOT / "docs" / "presentation" / "numbers.tsv"


def main() -> int:
    missing = 0
    counts = {"ok": 0, "moved": 0, "MISSING": 0, "skip": 0}
    with open(TABLE, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE))
    for row in rows:
        label = f"[{row['section']}] {row['value_in_report']}: {row['claim']}"
        if row["run_file"] == "NOT FOUND":
            status, where = "skip", row["evidence_substring"]
        else:
            path = ROOT / row["run_file"]
            lines = path.read_text(errors="replace").splitlines() if path.exists() else []
            sub, n = row["evidence_substring"], int(row["line"])
            hits = [i + 1 for i, text in enumerate(lines) if sub in text]
            if n <= len(lines) and sub in lines[n - 1]:
                status, where = "ok", f"{row['run_file']}:{n}"
            elif hits:
                status, where = "moved", f"{row['run_file']}:{hits[0]} (table says {n})"
            else:
                status, where = "MISSING", f"{row['run_file']}:{n}"
                missing += 1
        counts[status] += 1
        print(f"{status:7} {label} -> {where}")
    print(", ".join(f"{k} {v}" for k, v in counts.items()))
    return 1 if missing else 0


if __name__ == "__main__":
    sys.exit(main())
