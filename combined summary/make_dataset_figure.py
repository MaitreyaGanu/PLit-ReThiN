"""Combined ARI/NMI figure for one dataset (Poisson + NB CSVs).
usage: python3 make_dataset_figure.py <poisson.csv> <nb.csv> <out.pdf>"""
import sys, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

P, N = pd.read_csv(sys.argv[1]), pd.read_csv(sys.argv[2])
df = pd.concat([P, N[N.Method == "PLit_NB"]])

INK, INK2, GRID, GREY, RAND = "#0b0b0b", "#52514e", "#e6e5e1", "#8f8e89", "#b5b4ae"
# (csv name, label, colour, marker, linewidth, linestyle, zorder)
# fixed entity -> hue mapping (validated 8-slot categorical palette, same in every dataset figure)
SPEC = [("PLit",              "PLit (Poisson)$^\\dagger$", "#eb6834", "s", 2.0, "-",  6),
        ("PLit_NB",           "PLit (NB)$^\\dagger$",      "#1baf7a", "^", 2.0, "-",  6),
        ("ReThiN",            "ReThiN$^\\dagger$",         "#2a78d6", "o", 2.0, "-",  6),
        ("M3Drop",            "M3Drop",                    "#eda100", "D", 1.2, "-",  4),
        ("Pearson_Residuals", "Pearson residuals",         "#e87ba4", "v", 1.2, "-",  4),
        ("scran_HVG",         "scran HVG",                 "#008300", "P", 1.2, "-",  4),
        ("scry_Deviance",     "scry Deviance",             "#4a3aa7", "X", 1.2, "-",  4),
        ("Seurat_VST",        "Seurat VST",                "#e34948", "*", 1.2, "-",  4),
        ("Random_Baseline",   "Random",                    GREY,      None, 1.0, "--", 2)]

plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "pdf.fonttype": 42})
fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.1), constrained_layout=True)
for ax, met in zip(axes, ["ARI", "NMI"]):
    for m, lab, c, mk, lw, ls, z in SPEC:
        d = df[df.Method == m].sort_values("K")
        ax.errorbar(d.K, d[met], yerr=d[met + "_sd"], color=c, marker=mk, ms=6 if mk != "*" else 8,
                    mec="white", mew=0.6, lw=lw, ls=ls, elinewidth=0.8, capsize=2, zorder=z, label=lab)
    ax.set_xscale("log"); ax.set_xticks([100, 200, 500, 1000]); ax.set_xticklabels(["100", "200", "500", "1000"])
    ax.minorticks_off()
    ax.set_xlabel("Genes selected ($K$)"); ax.set_ylabel(f"{met} (mean $\\pm$ SD)")
    ax.grid(axis="y", color=GRID, lw=0.6); ax.set_axisbelow(True)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
from matplotlib.lines import Line2D
h = [Line2D([], [], color=c, marker=mk, ms=6 if mk != "*" else 8, mec="white", mew=0.6, lw=lw, ls=ls)
     for _, _, c, mk, lw, ls, _ in SPEC]
l = [lab for _, lab, *_ in SPEC]
order = [0, 5, 1, 6, 2, 7, 3, 8, 4]          # column-major fill -> row 1: proposed + M3Drop, Pearson
fig.legend([h[i] for i in order], [l[i] for i in order], loc="outside upper center", ncol=5, frameon=False, fontsize=8,
           handlelength=2.2, columnspacing=1.2, labelcolor=INK)
fig.savefig(sys.argv[3]); fig.savefig(sys.argv[3].replace(".pdf", ".png"), dpi=200)
