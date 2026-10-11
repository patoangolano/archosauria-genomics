"""Step 2 - karyotype structure, microchromosome counts and GC compartmentation.

Reads data/ncbi_sequence_reports.json.gz; writes results/tabela2_cariotipo_gc.csv.
"""
import gzip, json
import numpy as np, pandas as pd
from scipy.stats import spearmanr

MICRO_MAX_BP, MACRO_MIN_BP = 20e6, 40e6
seqreps = json.load(gzip.open("../data/ncbi_sequence_reports.json.gz", "rt"))


def chrom_table(records):
    recs = [r for r in records if r.get("role") == "assembled-molecule"
            and r.get("assigned_molecule_location_type") == "Chromosome"]
    t = pd.DataFrame([{"chr": r["chr_name"], "length": r["length"],
                       "gc": r.get("gc_percent")} for r in recs])
    return t.sort_values("length", ascending=False).reset_index(drop=True)


rows = []
for sp, records in seqreps.items():
    t = chrom_table(records)
    auto = t[~t.chr.isin(["Z", "W", "X", "Y", "MT"])]
    micro, macro = auto[auto.length < MICRO_MAX_BP], auto[auto.length >= MACRO_MIN_BP]
    a = auto.dropna(subset=["gc"])
    rho, p = spearmanr(np.log10(a.length), a.gc) if len(a) > 4 else (np.nan, np.nan)
    rows.append({
        "species": sp, "n_assembled_chr": len(t), "n_autosomes": len(auto),
        "largest_Mb": round(t.length.max() / 1e6, 1),
        "smallest_auto_Mb": round(auto.length.min() / 1e6, 2),
        "n_micro_lt20Mb": len(micro),
        "pct_genome_in_micro": round(100 * micro.length.sum() / t.length.sum(), 1),
        "GC_macro_ge40Mb": round(float(np.average(macro.gc, weights=macro.length)), 2),
        "GC_micro_lt20Mb": round(float(np.average(micro.gc, weights=micro.length)), 2)
                            if len(micro) else None,
        "spearman_rho_len_gc": round(float(rho), 3), "spearman_P": f"{p:.2e}"})

df = pd.DataFrame(rows)
df["GC_delta"] = (df.GC_micro_lt20Mb - df.GC_macro_ge40Mb).round(2)
df.to_csv("../results/tabela2_cariotipo_gc.csv", index=False)
print(df.to_string(index=False))
