"""Do the port's selections name the same reads as `whatshap phase --output-read-list`?

For each depth and cap given, runs `select_whatshap` on the extracted
fragments (the reads calling two or more variants, in WhatsHap's order) and
compares the selected read names with the list WhatsHap wrote. Prints the
counts and the size of the symmetric difference.
"""

from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path

from phasebudget.selection import select_whatshap


def _extractor():
    path = Path(__file__).with_name("extract_fragments.py")
    spec = importlib.util.spec_from_file_location("extract_fragments", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fragments")
    parser.add_argument("read_list")
    parser.add_argument("--cap", type=int, required=True)
    args = parser.parse_args()
    matrix, positions, names, _ = _extractor().load_fragments(args.fragments)
    informative = matrix.subset([f.read_id for f in matrix if f.n_observed >= 2])
    ours = {names[r] for r in select_whatshap(informative, args.cap, positions=positions).selected}
    with open(args.read_list) as handle:
        theirs = {line.split("\t")[0] for line in handle if not line.startswith("#")}
    print(f"{Path(args.fragments).name} cap {args.cap}: port {len(ours)}, whatshap phase "
          f"{len(theirs)}, symmetric difference {len(ours ^ theirs)}")


if __name__ == "__main__":
    main()
