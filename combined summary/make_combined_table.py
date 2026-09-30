import sys, pandas as pd, numpy as np

ROWS = [("PLit", "PLit (Poisson)$^\\dagger$", "P"), ("PLit_NB", "PLit (NB)$^\\dagger$", "NB"),
        ("ReThiN", "ReThiN$^\\dagger$", "P"), None,
        ("M3Drop", "M3Drop", "P"), ("Pearson_Residuals", "Pearson Residuals", "P"),
        ("scran_HVG", "scran HVG", "P"), ("scry_Deviance", "scry Deviance", "P"),
        ("Seurat_VST", "Seurat VST", "P"), None,
        ("Random_Baseline", "Random Baseline", "P")]
KS = [100, 200, 500, 1000]

def main(pois_csv, nb_csv, dataset, label):
    P, N = pd.read_csv(pois_csv), pd.read_csv(nb_csv)
    # baseline consistency check
    for m in set(P.Method) & set(N.Method):
        a = P[P.Method == m].sort_values("K")[["ARI", "NMI"]].round(3).values
        b = N[N.Method == m].sort_values("K")[["ARI", "NMI"]].round(3).values
        if a.shape != b.shape or not np.allclose(a, b):
            print(f"% WARNING: {m} differs between Poisson and NB CSVs", file=sys.stderr)
    src = {"P": P, "NB": N}
    val = {}
    for r in ROWS:
        if r is None: continue
        m, _, s = r
        d = src[s][src[s].Method == m].set_index("K")
        for met in ["ARI", "NMI"]:
            for k in KS:
                val[(m, met, k)] = (round(d.loc[k, met], 3), round(d.loc[k, met + "_sd"], 3))
    top = {(met, k): sorted({val[(r[0], met, k)][0] for r in ROWS if r and r[0] != "Random_Baseline"}, reverse=True)[:3]
           for met in ["ARI", "NMI"] for k in KS}
    out = ["\\begin{table}[H]", "\\centering", "\\renewcommand{\\arraystretch}{1.2}",
           f"\\caption{{Results on the {dataset} dataset (mean $\\pm$ SD over 5 subsampling seeds). Per column, excluding Random Baseline: best in \\textbf{{bold}}, second \\uline{{underlined}}, third \\uuline{{double-underlined}}. $^\\dagger$Proposed method.}}", f"\\label{{tab:{label}_results}}",
           "\\resizebox{\\textwidth}{!}{%", "\\begin{tabular}{lcccccccc}", "\\toprule",
           "& \\multicolumn{4}{c}{\\textbf{ARI} (mean $\\pm$ SD)} & \\multicolumn{4}{c}{\\textbf{NMI} (mean $\\pm$ SD)} \\\\",
           "\\cmidrule(lr){2-5}\\cmidrule(lr){6-9}",
           "\\textbf{Method} & $K=100$ & $K=200$ & $K=500$ & $K=1000$ & $K=100$ & $K=200$ & $K=500$ & $K=1000$ \\\\",
           "\\midrule"]
    for r in ROWS:
        if r is None: out.append("\\midrule"); continue
        m, name, _ = r
        cells = []
        for met in ["ARI", "NMI"]:
            for k in KS:
                mu, sd = val[(m, met, k)]
                rk = top[(met, k)].index(mu) + 1 if m != "Random_Baseline" and mu in top[(met, k)] else 0
                if rk == 1:   cells.append(f"$\\mathbf{{{mu:.3f}}} \\pm {sd:.3f}$")
                elif rk == 2: cells.append(f"\\uline{{${mu:.3f}$}}${{}}\\pm {sd:.3f}$")
                elif rk == 3: cells.append(f"\\uuline{{${mu:.3f}$}}${{}}\\pm {sd:.3f}$")
                else:         cells.append(f"${mu:.3f} \\pm {sd:.3f}$")
        out.append(f"{name} & " + " & ".join(cells) + " \\\\")
    out += ["\\bottomrule",
            "\\end{tabular}%", "}", "\\end{table}"]
    print("\n".join(out))

if __name__ == "__main__":
    main(*sys.argv[1:5])
