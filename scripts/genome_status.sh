#!/usr/bin/env bash
# Progress of the genome-wide run and the trio benchmark.
cd "$(dirname "$0")/.."
status=data/genome/status.json
if [ "${1:-}" = ont ]; then status=data/genome_ont/status.json; fi
for unit in phasebudget-genome phasebudget-genome-ont phasebudget-tools phasebudget-trio; do
  printf '%-20s %s\n' "$unit" "$(systemctl --user is-active $unit 2>/dev/null)"
done
if [ -f "$status" ]; then
  STATUS="$status" python3 - <<'PY'
import json, os
s = json.load(open(os.environ["STATUS"]))
print(f"updated {s['time']}: {s['done']} jobs done, {s.get('pending', 0)} pending, "
      f"{len(s.get('running', []))} running, {len(s['failed'])} failed")
for name in s.get("running", []):
    print("  running", name)
for name in s["failed"][:20]:
    print("  FAILED ", name)
if "free_disk_gb" in s:
    print(f"free disk {s['free_disk_gb']} GB, available memory {s['mem_available_gb']} GB")
PY
fi
tail -3 tmp/t3_trio_budget.txt 2>/dev/null
