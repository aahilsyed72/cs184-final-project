#!/usr/bin/env python3
"""Summarize results/timing.csv (median over repeats) and results/accuracy.csv
(pixel mismatch vs. the fine-step reference). Prints markdown tables and saves
them to results/bench_summary.md."""
import csv
import statistics
from collections import defaultdict
from pathlib import Path

import matplotlib.image as mpimg
import numpy as np

root = Path(__file__).resolve().parent.parent
res_dir = root / "results"
out = []

# ---- timing
groups = defaultdict(list)
for row in csv.DictReader(open(res_dir / "timing.csv")):
    key = (int(row["width"]), int(row["height"]), row["step_mode"], row["step_param"], int(row["threads"]))
    groups[key].append(row)

out.append("### Render timing (median of repeats, 1 spp, checker sky, camera at r = 20)\n")
out.append("| resolution | step mode | threads | render time (s) | speedup vs 1 thread | total RK4 steps | mean steps/ray | steps/s |")
out.append("|---|---|---|---|---|---|---|---|")
for key in sorted(groups, key=lambda k: (k[0], k[2], k[4])):
    w, h, mode, param, thr = key
    rows = groups[key]
    t = statistics.median(float(r["render_s"]) for r in rows)
    t1 = statistics.median(float(r["render_s"]) for r in groups[(w, h, mode, param, 1)])
    steps = int(rows[0]["total_steps"])
    out.append(f"| {w}x{h} | {mode} {param} | {thr} | {t:.3f} | {t1 / t:.2f}x | {steps:,} | "
               f"{float(rows[0]['mean_steps_per_ray']):.0f} | {steps / t:.3g} |")

# ---- accuracy
acc_csv = res_dir / "accuracy.csv"
if acc_csv.exists():
    rows = list(csv.DictReader(open(acc_csv)))
    ref = mpimg.imread(res_dir / "accuracy" / "ref.png")[..., :3]
    ref_steps = int(rows[0]["total_steps"])
    out.append("\n### Accuracy vs. cost (640x360, vs. adaptive step 0.002 reference)\n")
    out.append("| step mode | mean steps/ray | render time (s) | captured % | pixels differing from reference |")
    out.append("|---|---|---|---|---|")
    for r in rows:
        mode, p = r["step_mode"], r["step_param"]
        name = "ref" if (mode == "adaptive" and float(p) == 0.002) else f"{mode}_{p}"
        img = mpimg.imread(res_dir / "accuracy" / f"{name}.png")[..., :3]
        # A pixel "differs" if any channel is off by more than 2/255.
        diff = (np.abs(img - ref) > 2.5 / 255).any(axis=2).mean() * 100
        label = f"{mode} {p}" + (" (reference)" if name == "ref" else "")
        out.append(f"| {label} | {float(r['mean_steps_per_ray']):.0f} | {float(r['render_s']):.3f} | "
                   f"{float(r['captured_pct']):.3f} | {diff:.3f}% |")

text = "\n".join(out) + "\n"
(res_dir / "bench_summary.md").write_text(text)
print(text)
