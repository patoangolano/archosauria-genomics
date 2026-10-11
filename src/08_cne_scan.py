"""Cross-species BLAST of Hox clusters -> conserved blocks, split coding vs non-coding (CNE)."""
import json, os, subprocess, sys
from collections import defaultdict

CLUSTERS = ["A", "B", "C", "D"]
QUERY = "Gallus_gallus"
TARGETS = ["Alligator_mississippiensis", "Crocodylus_porosus"]
MIN_LEN, MIN_PID = 50, 65.0
os.makedirs("../data/blast", exist_ok=True)

def run_blast(q, s, out, task="dc-megablast"):
    subprocess.run(["makeblastdb", "-in", s, "-dbtype", "nucl", "-out", out + ".db"],
                   check=True, capture_output=True)
    cols = "qseqid sseqid pident length qstart qend sstart send evalue bitscore"
    r = subprocess.run(["blastn", "-task", task, "-query", q, "-db", out + ".db",
                        "-outfmt", "6 " + cols, "-evalue", "1e-5",
                        "-num_threads", "8", "-max_target_seqs", "5000"],
                       check=True, capture_output=True, text=True)
    hits = []
    for line in r.stdout.strip().splitlines():
        if not line: continue
        f = line.split("\t")
        hits.append({"pident": float(f[2]), "length": int(f[3]),
                     "qstart": int(f[4]), "qend": int(f[5]),
                     "sstart": int(f[6]), "send": int(f[7]),
                     "evalue": float(f[8]), "bitscore": float(f[9])})
    return hits

def merge(ivs):
    if not ivs: return []
    ivs = sorted((min(a, b), max(a, b)) for a, b in ivs)
    out = [list(ivs[0])]
    for a, b in ivs[1:]:
        if a <= out[-1][1] + 1: out[-1][1] = max(out[-1][1], b)
        else: out.append([a, b])
    return [tuple(x) for x in out]

def overlap(iv, ivs):
    a, b = iv
    return sum(max(0, min(b, y) - max(a, x) + 1) for x, y in ivs)

results, cne_table = {}, []
for cl in CLUSTERS:
    qf = f"../data/hox/{QUERY}_HOX{cl}.fa"
    qlen = sum(len(l.strip()) for l in open(qf) if not l.startswith(">"))
    # chicken coding exons inside the chicken region
    ex = run_blast(f"../data/hox/{QUERY}_cds.fa", qf, f"../data/blast/{cl}_exon", task="megablast")
    exons = merge([(h["sstart"], h["send"]) for h in ex if h["pident"] >= 95 and h["length"] >= 30])
    exon_bp = sum(b - a + 1 for a, b in exons)
    for tg in TARGETS:
        hits = run_blast(qf, f"../data/hox/{tg}_HOX{cl}.fa", f"../data/blast/{cl}_{tg[:3]}")
        hits = [h for h in hits if h["length"] >= MIN_LEN and h["pident"] >= MIN_PID]
        blocks = merge([(h["qstart"], h["qend"]) for h in hits])
        cons_bp = sum(b - a + 1 for a, b in blocks)
        cnes = [bl for bl in blocks if overlap(bl, exons) < 0.2 * (bl[1] - bl[0] + 1)]
        cne_bp = sum(b - a + 1 for a, b in cnes)
        # identity of hits whose query block is non-coding
        ncs = [h for h in hits if overlap((h["qstart"], h["qend"]), exons)
               < 0.2 * (h["qend"] - h["qstart"] + 1)]
        pid_nc = sum(h["pident"] * h["length"] for h in ncs) / max(1, sum(h["length"] for h in ncs))
        pid_all = sum(h["pident"] * h["length"] for h in hits) / max(1, sum(h["length"] for h in hits))
        lens = sorted(b - a + 1 for a, b in cnes)
        cne_table.append({
            "cluster": f"HOX{cl}", "target": tg, "region_kb": round(qlen / 1000, 1),
            "exon_bp": exon_bp, "n_blocks": len(blocks), "cons_bp": cons_bp,
            "pct_region_conserved": round(100 * cons_bp / qlen, 2),
            "n_CNE": len(cnes), "CNE_bp": cne_bp,
            "pct_region_CNE": round(100 * cne_bp / qlen, 2),
            "CNE_median_bp": lens[len(lens) // 2] if lens else 0,
            "CNE_max_bp": lens[-1] if lens else 0,
            "mean_pid_all": round(pid_all, 1), "mean_pid_noncoding": round(pid_nc, 1)})
        results[f"HOX{cl}|{tg}"] = {"blocks": blocks, "cnes": cnes, "exons": exons,
                                    "qlen": qlen, "n_hits": len(hits)}

json.dump({"table": cne_table, "detail": {k: {"cnes": v["cnes"], "qlen": v["qlen"],
           "exons": v["exons"], "blocks": v["blocks"]} for k, v in results.items()}},
          open("../data/cne_results.json", "w"))
import csv
with open("../results/cne_hox_summary.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(cne_table[0].keys())); w.writeheader(); w.writerows(cne_table)
for r in cne_table:
    print(f"{r['cluster']} vs {r['target'][:12]:12s} region={r['region_kb']:6.1f}kb "
          f"cons={r['pct_region_conserved']:5.2f}% CNEs={r['n_CNE']:4d} "
          f"CNEbp={r['CNE_bp']:6d} ({r['pct_region_CNE']:5.2f}%) "
          f"medCNE={r['CNE_median_bp']:4d} maxCNE={r['CNE_max_bp']:5d} "
          f"pid_nc={r['mean_pid_noncoding']:.1f}%")
