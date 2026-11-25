"""Round 4, W3b: train the neural per-pair score on a GPU and write its estimate for every pair.

`docs/design/w3b-neural-score.md` states the model. For one depth and phaser,
it reads `w3_features.py`'s and `w3_reads.py`'s arrays for the seven genomes of
`w3_pooled.py`, trains on chromosomes 1 to 6 of all of them, and writes the
estimated probability that each pair is wrong to
`data/w3_neural/<platform>_<depth>_<phaser>.npy`, in the features' order.
`w3_pooled.py` then calibrates and scores it.

Usage: `PYTHONPATH=src:benchmarks .venv-gpu/bin/python benchmarks/w3_neural.py --depth 30
--phaser whatshap --device cuda:0`
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

OUT = HERE.parent / "data" / "w3_neural"
EPOCHS, BATCH, LR, MAX_WEIGHT, SEED = 3, 4096, 1e-3, 100.0, 0


def load(platform: str, depth: str, phaser: str, stratum):
    """`(pair features, read tensor, read mask, chromosomes, wrong)` for one genome."""
    import numpy as np
    import torch

    import w3_gate as w3
    from w3_reads import OUT as READS

    x, chrom, pos1, wrong, _ = w3.load(depth, phaser, stratum, platform)
    r = np.load(READS / f"{platform}_{depth}_{phaser}.npz")
    assert (r["pos1"] == pos1).all() and (r["chrom"] == chrom).all(), platform
    x = np.sign(x) * np.log1p(np.abs(x))  # counts and distances span orders of magnitude
    calls = torch.from_numpy(r["calls"]).float()
    weights = torch.from_numpy(r["weights"].astype(np.float32)) / 60.0
    called = torch.log1p(torch.from_numpy(r["read_calls"].astype(np.float32))).unsqueeze(-1)
    reads = torch.cat([calls, calls.abs() * weights, called], dim=-1)  # pairs x 32 x 17
    mask = torch.from_numpy(r["read_calls"] > 0).float()
    return torch.from_numpy(x).float(), reads, mask, chrom, torch.from_numpy(wrong).float()


def model(n_pair: int):
    import torch
    from torch import nn

    class Net(nn.Module):
        def __init__(self):
            super().__init__()
            self.read = nn.Sequential(nn.Linear(17, 64), nn.ReLU(), nn.Linear(64, 64), nn.ReLU())
            self.head = nn.Sequential(nn.Linear(128 + n_pair, 64), nn.ReLU(), nn.Linear(64, 1))

        def forward(self, pair, reads, mask):
            h = self.read(reads)
            m = mask.unsqueeze(-1)
            mean = (h * m).sum(1) / m.sum(1).clamp(min=1)
            top = torch.where(m > 0, h, torch.full_like(h, -1e4)).max(1).values
            top = torch.where(m.sum(1) > 0, top, torch.zeros_like(top))
            return self.head(torch.cat([mean, top, pair], dim=-1)).squeeze(-1)

    return Net()


def main() -> None:
    import numpy as np
    import torch

    import genome_common as common
    import w3_gate as w3
    from stratified_gate import strata
    from w3_pooled import GENOMES

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--depth", required=True)
    parser.add_argument("--phaser", required=True)
    parser.add_argument("--device", default="cuda:0")
    args = parser.parse_args()
    torch.manual_seed(SEED)
    stratum = strata(list(common.CONFIRMATORY))
    data = {g: load(g, args.depth, args.phaser, stratum) for g in GENOMES}
    train = [np.isin(d[3], w3.TRAIN) for d in data.values()]
    pair = torch.cat([d[0][t] for d, t in zip(data.values(), train)])
    mean, std = pair.mean(0), pair.std(0).clamp(min=1e-6)
    pair = (pair - mean) / std
    reads = torch.cat([d[1][t] for d, t in zip(data.values(), train)])
    mask = torch.cat([d[2][t] for d, t in zip(data.values(), train)])
    wrong = torch.cat([d[4][t] for d, t in zip(data.values(), train)])
    positives = float(wrong.sum())
    weight = min(MAX_WEIGHT, (len(wrong) - positives) / max(1.0, positives))
    print(f"{args.depth}x {args.phaser}: train {len(wrong)} pairs, {int(positives)} wrong, "
          f"positive weight {weight:.1f}, device {args.device}", flush=True)
    net = model(pair.shape[1]).to(args.device)
    opt = torch.optim.Adam(net.parameters(), lr=LR)
    loss_fn = torch.nn.BCEWithLogitsLoss(pos_weight=torch.tensor(weight, device=args.device))
    generator = torch.Generator().manual_seed(SEED)
    for epoch in range(EPOCHS):
        order, total = torch.randperm(len(wrong), generator=generator), 0.0
        for k in range(0, len(order), BATCH):
            idx = order[k:k + BATCH]
            out = net(pair[idx].to(args.device), reads[idx].to(args.device),
                      mask[idx].to(args.device))
            loss = loss_fn(out, wrong[idx].to(args.device))
            opt.zero_grad()
            loss.backward()
            opt.step()
            total += float(loss) * len(idx)
        print(f"  epoch {epoch + 1}: mean loss {total / len(order):.4f}", flush=True)
    net.eval()
    OUT.mkdir(parents=True, exist_ok=True)
    with torch.no_grad():
        for g, (x, r, m, _, _) in data.items():
            x = (x - mean) / std
            p = torch.cat([torch.sigmoid(net(x[k:k + BATCH].to(args.device),
                                             r[k:k + BATCH].to(args.device),
                                             m[k:k + BATCH].to(args.device))).cpu()
                           for k in range(0, len(x), BATCH)])
            np.save(OUT / f"{g}_{args.depth}_{args.phaser}.npy", p.numpy().astype(np.float32))
    print(f"  wrote {len(data)} genomes", flush=True)


if __name__ == "__main__":
    main()
