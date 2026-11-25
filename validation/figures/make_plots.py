"""Turn validation_results.csv (written by validation.exe) into thesis figures.
Usage:  python make_plots.py [validation_results.csv]      -> fig_uniaxial.png, fig_cantilever.png, fig_patch.png, fig_validation_all.png
Needs: python, matplotlib."""
import csv, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

path = sys.argv[1] if len(sys.argv) > 1 else "validation_results.csv"
tests = {"uniaxial": {}, "cantilever": {}, "patch_indentation": {}}
with open(path, newline="") as f:
    for row in csv.DictReader(f):
        t = row["test"]
        if t not in tests:
            continue
        nu = row["case"].split("=")[1].split(";")[0]
        tests[t].setdefault(nu, []).append((float(row["h_mm_or_info"]), float(row["error_percent"])))

titles = {
    "uniaxial": ("Test 1: clamped bar in compression", "error of end displacement [%]"),
    "cantilever": ("Test 2: cantilever beam, tip load", "error of tip deflection [%]"),
    "patch_indentation": ("Test 3: half-space under the proxy load patch", "error of centre displacement [%]"),
}
colors = {"0.30": "#1f77b4", "0.45": "#d62728"}

def draw(ax, key):
    for nu, pts in sorted(tests[key].items()):
        pts.sort()
        ax.plot([p[0] for p in pts], [p[1] for p in pts], "o-", color=colors.get(nu, "k"), label=f"$\\nu$ = {nu}")
    ax.axhline(0, color="0.4", lw=0.8)
    ax.set_xlabel("tetrahedral cell size h [mm]")
    ax.set_ylabel(titles[key][1])
    ax.set_title(titles[key][0], fontsize=10)
    ax.grid(alpha=0.3)
    ax.legend()

names = {"uniaxial": "fig_uniaxial.png", "cantilever": "fig_cantilever.png", "patch_indentation": "fig_patch.png"}
for key, fn in names.items():
    fig, ax = plt.subplots(figsize=(5, 3.6))
    draw(ax, key)
    fig.tight_layout(); fig.savefig(fn, dpi=200); plt.close(fig)

fig, axes = plt.subplots(1, 3, figsize=(14, 3.8))
for ax, key in zip(axes, names):
    draw(ax, key)
fig.tight_layout(); fig.savefig("fig_validation_all.png", dpi=200)
print("written:", ", ".join(list(names.values()) + ["fig_validation_all.png"]))
