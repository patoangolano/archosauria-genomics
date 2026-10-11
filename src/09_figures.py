"""Step 9 - the three publication figures.

Writes figures/fig1_arquitetura_genomica.png, fig2_taxas_moleculares.png
and fig3_hox_cnes.png.
"""
import gzip, json
import numpy as np, pandas as pd
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.stats import spearmanr

mpl.rcParams.update({"font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
                     "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 7,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 150, "savefig.dpi": 300})
AVE, CROC, TURT, GREY = "#1f6f8b", "#c4761a", "#6b8f4e", "#9a9a9a"
STYLE = {"Gallus gallus": (AVE, "o", 1.0),
         "Alligator mississippiensis": (CROC, "s", 1.0),
         "Chrysemys picta": (TURT, "^", 0.75)}


def panel_letter(ax, letter):
    ax.text(-0.18, 1.02, letter, transform=ax.transAxes, fontsize=10,
            fontweight="bold", va="bottom", ha="left")


def chrom_table(records):
    recs = [r for r in records if r.get("role") == "assembled-molecule"
            and r.get("assigned_molecule_location_type") == "Chromosome"]
    t = pd.DataFrame([{"chr": r["chr_name"], "length": r["length"],
                       "gc": r.get("gc_percent")} for r in recs])
    return t.sort_values("length", ascending=False).reset_index(drop=True)


seqreps = json.load(gzip.open("../data/ncbi_sequence_reports.json.gz", "rt"))
karyo = {sp: chrom_table(v) for sp, v in seqreps.items()}

# ---------------- Figure 1 ----------------
fig1, axs = plt.subplots(1, 2, figsize=(7.2, 3.1))
ax = axs[0]
for sp, (c, m, al) in STYLE.items():
    t = karyo[sp]
    t = t[~t.chr.isin(["W", "MT"])].sort_values("length", ascending=False)
    ax.plot(range(1, len(t) + 1), t.length / 1e6, marker=m, ms=3.2, lw=1.1, color=c, alpha=al)
    ax.text(len(t) + 1.2, t.length.min() / 1e6, sp.split()[0], color=c, fontsize=6,
            va="center", style="italic")
ax.axhline(20, color=GREY, lw=0.8, ls=":", zorder=0)
ax.text(51, 0.42, "limiar de microcromossomo: 20 Mb", color=GREY, fontsize=6,
        ha="right", va="center")
ax.set_yscale("log"); ax.set_xlim(0, 52)
ax.set_xlabel("cromossomo (ordenado por tamanho)"); ax.set_ylabel("comprimento (Mb)")
ax.set_yticks([1, 10, 100]); ax.set_yticklabels(["1", "10", "100"])
ax.set_title("A galinha tem 28 microcromossomos;\no jacaré, nenhum")
panel_letter(ax, "a")

ax = axs[1]
RHO_AT = {"Gallus gallus": (0.17, 40.0), "Alligator mississippiensis": (26, 57.5),
          "Chrysemys picta": (26, 55.2)}
for sp, (c, m, al) in STYLE.items():
    t = karyo[sp].dropna(subset=["gc"]); t = t[~t.chr.isin(["W", "MT"])]
    ax.scatter(t.length / 1e6, t.gc, s=13, color=c, marker=m, alpha=al,
               edgecolor="white", linewidth=0.3, zorder=3)
    x = np.log10(t.length / 1e6)
    slope, intercept = np.polyfit(x, t.gc, 1)
    xs = np.linspace(x.min(), x.max(), 50)
    ax.plot(10 ** xs, intercept + slope * xs, color=c, lw=1.1, alpha=0.85, zorder=2)
    rho, _ = spearmanr(x, t.gc)
    ax.text(*RHO_AT[sp], f"{sp.split()[0]}  \u03c1 = {rho:.2f}", color=c, fontsize=6,
            style="italic", ha="left", va="center", zorder=5)
ax.set_xscale("log"); ax.set_ylim(37.5, 66)
ax.set_xlabel("comprimento do cromossomo (Mb)"); ax.set_ylabel("conteúdo GC (%)")
ax.set_xticks([1, 10, 100]); ax.set_xticklabels(["1", "10", "100"])
ax.set_title("O gradiente GC–tamanho é comum\na aves, crocodilianos e quelônios")
panel_letter(ax, "b")
fig1.tight_layout(); fig1.savefig("../figures/fig1_arquitetura_genomica.png", bbox_inches="tight")

# ---------------- Figure 2 ----------------
def parse_newick(s):
    s = s.strip().rstrip(";"); pos = [0]

    def node():
        if s[pos[0]] == "(":
            pos[0] += 1
            ch = [node()]
            while s[pos[0]] == ",":
                pos[0] += 1
                ch.append(node())
            pos[0] += 1
            nm = ""
            while pos[0] < len(s) and s[pos[0]] not in ",():":
                nm += s[pos[0]]; pos[0] += 1
        else:
            nm, ch = "", []
            while s[pos[0]] not in ",():":
                nm += s[pos[0]]; pos[0] += 1
        bl = 0.0
        if pos[0] < len(s) and s[pos[0]] == ":":
            pos[0] += 1; b = ""
            while pos[0] < len(s) and s[pos[0]] not in ",()":
                b += s[pos[0]]; pos[0] += 1
            bl = float(b)
        return {"name": nm, "bl": bl, "ch": ch}
    return node()


PT = {"Gallus_gallus": ("Gallus gallus", AVE),
      "Alligator_mississippiensis": ("Alligator mississippiensis", CROC),
      "Crocodylus_porosus": ("Crocodylus porosus", CROC),
      "Gavialis_gangeticus": ("Gavialis gangeticus", CROC),
      "Chrysemys_picta": ("Chrysemys picta", TURT),
      "Anolis_carolinensis": ("Anolis carolinensis", GREY),
      "Homo_sapiens": ("Homo sapiens", GREY)}
ORDER = list(PT)
tree = parse_newick(open("../alignments/trees/p3_fix.treefile").read())
ypos, segs = {}, []


def layout(n, x0):
    x = x0 + n["bl"]
    if not n["ch"]:
        y = ORDER.index(n["name"]); ypos[n["name"]] = y
        segs.append((x0, x, y, y, n["name"])); return y
    ys = [layout(c, x) for c in n["ch"]]
    y = sum(ys) / len(ys)
    segs.append((x0, x, y, y, None)); segs.append((x, x, min(ys), max(ys), "V"))
    return y


layout(tree, 0.0)
rates = json.load(open("../data/rate_results.json"))["rates"]
fig2, axs = plt.subplots(1, 2, figsize=(7.2, 3.0), gridspec_kw={"width_ratios": [1.25, 1]})
ax = axs[0]
for x0, x1, y0, y1, nm in segs:
    ax.plot([x0, x1], [y0, y1], color=PT[nm][1] if nm in PT else "#444444",
            lw=1.4 if nm in PT else 1.0, solid_capstyle="butt")
for sp, y in ypos.items():
    xt = max(x1 for x0, x1, ya, yb, nm in segs if nm == sp)
    ax.text(xt + 0.012, y, PT[sp][0], fontsize=6.4, va="center", style="italic", color=PT[sp][1])
ax.set_ylim(-0.7, 6.7); ax.invert_yaxis(); ax.set_yticks([]); ax.set_xlim(-0.02, 1.08)
ax.spines[["left", "top", "right"]].set_visible(False)
ax.set_xlabel("substituições por sítio (3ª posição de códon)")
ax.set_title("Na 3ª posição de códon, o ramo da galinha\né ~2× mais longo que o dos crocodilianos")
panel_letter(ax, "a")

ax = axs[1]
gal = [r["Gallus"] for r in rates]
cro = [r["crocodylia_mean"] for r in rates]
xs, w = np.arange(3), 0.34
ax.bar(xs - w / 2, gal, w, color=AVE, label="Gallus gallus (ave)")
ax.bar(xs + w / 2, cro, w, color=CROC, label="Crocodylia (média de 3 spp.)")
for i, (g, c) in enumerate(zip(gal, cro)):
    ax.text(i, max(g, c) + 0.012, f"{g / c:.2f}×", ha="center", fontsize=6.5, color="#333333")
ax.set_xticks(xs); ax.set_xticklabels(["Proteína", "Códon 1+2", "Códon 3"])
ax.set_ylabel("distância desde o ancestral\ndos arcossauros (subst./sítio)")
ax.set_title("A aceleração aviária concentra-se\nnos sítios quase neutros")
ax.legend(frameon=False, fontsize=6, loc="upper left"); ax.margins(y=0.18)
panel_letter(ax, "b")
fig2.tight_layout(); fig2.savefig("../figures/fig2_taxas_moleculares.png", bbox_inches="tight")

# ---------------- Figure 3 ----------------
cne = json.load(open("../data/cne_results.json"))
hox = pd.DataFrame(json.load(open("../data/hox_locations.json"))).drop_duplicates(
    subset=["species", "symbol"])
hox["cluster"] = hox.symbol.str[3]
spans = pd.read_csv("../results/tabela5_hox_sintenia.csv")
FLANK = 20_000

fig3, axs3 = plt.subplots(4, 1, figsize=(7.2, 5.0))
for k, cl in enumerate("ABCD"):
    ax = axs3[k]
    det = cne["detail"][f"HOX{cl}|Alligator_mississippiensis"]
    L = det["qlen"]
    start = int(spans[(spans.species == "Gallus_gallus") &
                      (spans.cluster == cl)].start.iloc[0]) - FLANK
    ax.add_patch(plt.Rectangle((0, 0.30), L / 1000, 0.14, color="#e9e9e9", lw=0))
    for a, b in det["exons"]:
        ax.add_patch(plt.Rectangle((a / 1000, 0.28), (b - a + 1) / 1000, 0.18,
                                   color="#555555", lw=0))
    for a, b in det["cnes"]:
        ax.add_patch(plt.Rectangle((a / 1000, 0.52), max((b - a + 1) / 1000, 0.35), 0.24,
                                   color=AVE, lw=0))
    sub = hox[(hox.species == "Gallus_gallus") & (hox.cluster == cl)].sort_values("begin")
    prev_x, lo = -99, True
    for _, g in sub.iterrows():
        xm = ((g.begin + g.end) / 2 - start + 1) / 1000
        lo = True if (xm - prev_x) > L / 1000 * 0.045 else not lo
        ax.text(xm, 0.045 if lo else -0.085, g.symbol.replace("HOX" + cl, ""),
                fontsize=5.4, ha="center", color="#333333")
        prev_x = xm
    bp = sum(b - a + 1 for a, b in det["cnes"])
    ax.text(L / 1000, 0.92, f"{len(det['cnes'])} CNEs · {bp / 1000:.1f} kb",
            ha="right", va="center", fontsize=6, color=AVE)
    ax.set_xlim(-2, L / 1000 + 2); ax.set_ylim(-0.17, 1.05); ax.set_yticks([])
    ax.spines[["left", "top", "right"]].set_visible(False)
    ax.set_ylabel(f"HOX{cl}", rotation=0, ha="right", va="center", fontsize=7.5)
axs3[3].set_xlabel("posição no cluster da galinha (kb)")
axs3[0].set_title("Elementos não-codificantes conservados entre galinha e jacaré\n"
                  "persistem nos quatro clusters Hox")
fig3.legend(handles=[Patch(color=AVE, label="CNE (conservado, não-codificante)"),
                     Patch(color="#555555", label="éxon codificante Hox")],
            frameon=False, fontsize=6.5, ncol=2, loc="lower center", bbox_to_anchor=(0.55, -0.025))
fig3.tight_layout(rect=[0, 0.035, 1, 1])
fig3.savefig("../figures/fig3_hox_cnes.png", bbox_inches="tight")
print("figuras escritas em ../figures/")
