"""Regenerate Fig. 6 and Table 2 from the DataScope exports."""
import csv
from pathlib import Path
from statistics import mean, variance

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
FILES = {"PD": "pid.csv", "Adaptive": "adaptive.csv",
         "Fixed-learning": "selflearning.csv", "FSLC": "slc.csv"}
PUBLISHED = {"PD": (-0.8749, 1.8725, 1.2183, 1.8533),
             "Adaptive": (-0.8484, 0.4395, 1.2526, 0.1873),
             "Fixed-learning": (0.8450, 0.2190, 1.1239, 0.1954),
             "FSLC": (0.8004, 0.1209, 0.9094, 0.0979)}


def main():
    output = ROOT / "results"
    output.mkdir(exist_ok=True)
    fig, axes = plt.subplots(2, 1, sharex=True, figsize=(8, 5))
    stats = []
    for method, filename in FILES.items():
        with (ROOT / "data" / "hardware" / filename).open(encoding="utf-8-sig", newline="") as f:
            rows = {row[0].strip(): row[1:] for row in csv.reader(f) if row}
        ids = [int(i) for i in rows["DataScope"][299:800]]
        assert ids == list(range(300, 801)), filename
        values = []
        for ax, channel in zip(axes, ("CH1", "CH2")):
            angles = [float(v) for v in rows[channel + ":"][299:800]]
            assert len(angles) == 501, filename
            ax.plot([i / 100 for i in ids], angles, label=method, linewidth=1)
            ax.set_ylabel(f"{channel} (deg)")
            values.extend((mean(angles), variance(angles)))
        if tuple(round(v, 4) for v in values) != PUBLISHED[method]:
            raise ValueError(f"{filename} differs from Table 2")
        stats.append([method, *values])
    axes[0].legend(ncol=4, fontsize=8)
    axes[-1].set_xlabel("Time (s)")
    axes[-1].set_xlim(3, 8)
    fig.tight_layout()
    fig.savefig(output / "measured_tilt.png", dpi=220)
    plt.close(fig)
    with (output / "table2.csv").open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Method", "CH1Mean_deg", "CH1Variance_deg2",
                         "CH2Mean_deg", "CH2Variance_deg2"])
        writer.writerows(reversed(stats))
    print("All 16 entries match Table 2 at four decimal places.")


if __name__ == "__main__":
    main()
