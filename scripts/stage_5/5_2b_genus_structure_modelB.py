import re
import sqlite3
import sys
import platform
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from scipy.stats import binomtest

rng = np.random.default_rng(11)

plt.rcParams.update({
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "axes.grid": False,
    "font.size": 12,
})


MODEL = "B"

INPUT_FILE = "out/final_feature_table.csv"

FOCUS_FEATURE = "nirK_aa_frac_W" if MODEL == "A" else "pcuac_aa_frac_W"

# Protein length column, used to convert the fraction back into a residue
# count. Must be the same length used as the denominator in 4.2.
LENGTH_COL = "nirK_length" if MODEL == "A" else "pcuac_length"

# The class shared by both models is "both". The contrasting class differs.
NEG_CLASS = "nirK_only" if MODEL == "A" else "pcuac_only"
POS_CLASS = "both"

# Readable labels for figures. Avoid raw column names in anything a
# non-specialist reader will see.
NEG_LABEL = "NirK only" if MODEL == "A" else "PCuAC only"
POS_LABEL = "Both proteins"
PROTEIN = "NirK" if MODEL == "A" else "PCuAC"
PARTNER = "PCuAC" if MODEL == "A" else "NirK"

# Single-letter residue code, taken from the feature name
RESIDUE = FOCUS_FEATURE.rsplit("_", 1)[-1]
RESIDUE_NAME_MAP = {
    "W": "tryptophan", "C": "cysteine", "A": "alanine", "M": "methionine",
    "H": "histidine", "Y": "tyrosine", "F": "phenylalanine", "G": "glycine",
    "P": "proline",
}
RESIDUE_NAME = RESIDUE_NAME_MAP.get(RESIDUE, f"{RESIDUE} residue")

OUTPUT_DIR = Path("out") / "genus_structure" / f"model_{MODEL}"
PLOTS_DIR = OUTPUT_DIR / "plots"
TABLES_DIR = OUTPUT_DIR / "tables"
PLOTS_DIR.mkdir(parents=True, exist_ok=True)
TABLES_DIR.mkdir(parents=True, exist_ok=True)

# Analysis thresholds
MIN_PER_CLASS = 5     # per-class minimum for a within-genus comparison
MIN_SPECIES = 8        # minimum genus size for the strip figure
MAX_GENERA = 25        # genera shown in the strip figure
MIN_COUNT_GROUP = 20   # minimum species before a count group is plotted
N_PERM = 500            # permutation replicates for the variance null
TOP_N_COMPOSITION = 15  # largest genera shown in the composition figure

USE_TAXONOMIZR = False
TAXONOMIZR_DB = "accessionTaxa.sql"

QUALIFIERS = {
    "uncultured", "unclassified", "candidatus", "unidentified",
    "endosymbiont", "bacterium", "bacteria", "archaeon", "organism",
}

PALETTE = {NEG_LABEL: "#B2182B", "Mixed": "#999999", POS_LABEL: "#2166AC"}
CLASS_COLOURS = {NEG_LABEL: "#B2182B", POS_LABEL: "#2166AC"}


def add_title(ax, title, subtitle=None):
    """Bold left-aligned title with a smaller grey subtitle underneath."""
    if subtitle:
        ax.text(0.0, 1.14, title, transform=ax.transAxes, fontsize=13,
                fontweight="bold", ha="left", va="bottom")
        ax.text(0.0, 1.03, subtitle, transform=ax.transAxes, fontsize=9,
                color="0.30", ha="left", va="bottom")
    else:
        ax.text(0.0, 1.03, title, transform=ax.transAxes, fontsize=13,
                fontweight="bold", ha="left", va="bottom")


def strip_top_right(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def derive_genus(species_key: pd.Series) -> pd.Series:
    """Genus is the first token of the species key, ignoring qualifiers."""
    def _first_token(s):
        if pd.isna(s):
            return np.nan
        cleaned = re.sub(r"[\[\]]", "", s)
        tokens = [t for t in re.split(r"\s+", cleaned)
                  if t != "" and t not in QUALIFIERS]
        return tokens[0] if tokens else np.nan
    return species_key.apply(_first_token)


df = pd.read_csv(INPUT_FILE)
assert FOCUS_FEATURE in df.columns, f"{FOCUS_FEATURE} not found in {INPUT_FILE}"
assert LENGTH_COL in df.columns, f"{LENGTH_COL} not found in {INPUT_FILE}"

model_df = df[df["class"].isin([NEG_CLASS, POS_CLASS])].copy()
model_df["genus"] = derive_genus(model_df["species"])
model_df["value"] = model_df[FOCUS_FEATURE]
model_df["value_pct"] = 100 * model_df["value"]
model_df["prot_length"] = model_df[LENGTH_COL]
model_df["class_label"] = pd.Categorical(
    np.where(model_df["class"] == POS_CLASS, POS_LABEL, NEG_LABEL),
    categories=[NEG_LABEL, POS_LABEL],
)
model_df = model_df[model_df["value"].notna()]

n_before = len(model_df)
model_df = model_df[model_df["genus"].notna()].copy()
n_unresolved = n_before - len(model_df)

print(f"Genus structure for Model {MODEL}")
print(f"Focus feature: {FOCUS_FEATURE} ( {RESIDUE_NAME} )")
print(f"Species analysed: {len(model_df)}")
print(f"Dropped (no resolvable genus): {n_unresolved}")
print(f"Distinct genera: {model_df['genus'].nunique()}")
print(model_df["class_label"].value_counts())
print()

def _lineage_rank(cur, taxid, rank, cache):
    """ NCBI taxonomy tree from taxid to the ancestor at `rank`."""
    key = (taxid, rank)
    if key in cache:
        return cache[key]
    seen = set()
    current = taxid
    result = np.nan
    while current is not None and current not in seen:
        seen.add(current)
        cur.execute("SELECT parent, rank FROM nodes WHERE id = ?", (current,))
        row = cur.fetchone()
        if row is None:
            break
        parent, node_rank = row
        if node_rank == rank:
            cur.execute(
                "SELECT name FROM names WHERE id = ? AND class = 'scientific name'",
                (current,),
            )
            name_row = cur.fetchone()
            result = name_row[0] if name_row else np.nan
            break
        if parent == current:
            break
        current = parent
    cache[key] = result
    return result


def assign_phylum(model_df, db_path, tables_dir):
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()
    cache = {}

    first_taxid = model_df["taxid"].astype(str).str.split(";").str[0].astype(int)

    model_df["phylum"] = first_taxid.apply(lambda t: _lineage_rank(cur, t, "phylum", cache))
    # `class` is both a taxonomic rank and the outcome variable, so only
    # phylum and genus are taken from the lineage lookup.
    ncbi_genus = first_taxid.apply(lambda t: _lineage_rank(cur, t, "genus", cache))
    model_df["ncbi_genus"] = ncbi_genus.astype(str).str.lower()
    conn.close()

    n_no_phylum = model_df["phylum"].isna().sum()
    print(f"Phylum assigned for {len(model_df) - n_no_phylum} species")
    print(f"No phylum assignment: {n_no_phylum}")
    print(model_df["phylum"].value_counts(dropna=False).sort_values(ascending=False))

    disagreement = model_df[
        model_df["ncbi_genus"].notna() & (model_df["ncbi_genus"] != model_df["genus"])
    ][["species", "genus", "ncbi_genus", "phylum", "class"]]
    disagreement.to_csv(tables_dir / "genus_parsing_disagreements.csv", index=False)
    print(f"\nGenus parsing disagreements (MANUAL CHECK): {len(disagreement)}\n")

    phylum_by_class = (
        model_df[model_df["phylum"].notna()]
        .groupby("phylum")
        .apply(lambda g: pd.Series({
            "n_species": len(g),
            "n_positive": int((g["class"] == POS_CLASS).sum()),
        }))
        .reset_index()
    )
    phylum_by_class["prop_positive"] = (
        phylum_by_class["n_positive"] / phylum_by_class["n_species"]
    ).round(3)
    phylum_by_class = phylum_by_class.sort_values("n_species", ascending=False)
    phylum_by_class.to_csv(tables_dir / "phylum_by_class.csv", index=False)
    print(phylum_by_class)
    print()
    return True


has_phylum = False

if USE_TAXONOMIZR:
    if not Path(TAXONOMIZR_DB).exists():
        print(f"Taxonomy database {TAXONOMIZR_DB} not found; skipping phylum.")
        print("Run taxonomizr::prepareDatabase() once (in R) to create it.\n")
    else:
        has_phylum = assign_phylum(model_df, TAXONOMIZR_DB, TABLES_DIR)


genus_stats = (
    model_df.groupby("genus")
    .apply(lambda g: pd.Series({
        "n_species": len(g),
        "n_negative": int((g["class"] == NEG_CLASS).sum()),
        "n_positive": int((g["class"] == POS_CLASS).sum()),
        "q25": np.quantile(g["value_pct"], 0.25),
        "median_pct": g["value_pct"].median(),
        "q75": np.quantile(g["value_pct"], 0.75),
    }))
    .reset_index()
)
genus_stats["prop_positive"] = genus_stats["n_positive"] / genus_stats["n_species"]
genus_stats["iqr_pct"] = genus_stats["q75"] - genus_stats["q25"]


def _composition(p):
    if p == 0:
        return NEG_LABEL
    if p == 1:
        return POS_LABEL
    return "Mixed"


genus_stats["composition"] = pd.Categorical(
    genus_stats["prop_positive"].apply(_composition),
    categories=[NEG_LABEL, "Mixed", POS_LABEL],
)
genus_stats = genus_stats.sort_values("n_species", ascending=False)
genus_stats.to_csv(TABLES_DIR / "genus_stats.csv", index=False)

# For a species in a class-pure genus, the genus determines the label
# outright. No sequence feature is required to predict it.
purity = (
    genus_stats.groupby("composition", observed=False)
    .agg(n_genera=("genus", "count"), n_species=("n_species", "sum"))
    .reset_index()
)
purity.to_csv(TABLES_DIR / "genus_purity.csv", index=False)
print(purity)

n_pure_genera = int(purity.loc[purity["composition"] != "Mixed", "n_genera"].sum())
n_pure_species = int(purity.loc[purity["composition"] != "Mixed", "n_species"].sum())

print(f"\nClass-pure genera: {n_pure_genera} of {len(genus_stats)} "
      f"({100 * n_pure_genera / len(genus_stats):.0f}%)")
print(f"Species in class-pure genera: {n_pure_species} of {len(model_df)} "
      f"({100 * n_pure_species / len(model_df):.0f}%)\n")


comp_order = [NEG_LABEL, "Mixed", POS_LABEL]
purity_ord = purity.set_index("composition").reindex(comp_order).reset_index()

fig1, ax1 = plt.subplots(figsize=(7, 5))
colours = [PALETTE[c] for c in purity_ord["composition"]]
bars = ax1.bar(purity_ord["composition"].astype(str), purity_ord["n_genera"],
                width=0.65, color=colours)
for bar, n_genera, n_species in zip(bars, purity_ord["n_genera"], purity_ord["n_species"]):
    ax1.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
              f"{int(n_genera)} genera\n{int(n_species)} species",
              ha="center", va="bottom", fontsize=9)
ax1.set_ylim(0, purity_ord["n_genera"].max() * 1.18)
ax1.set_ylabel("Number of genera")
strip_top_right(ax1)
add_title(ax1, "Most genera contain only one class",
          "Co-occurrence is usually fixed across a genus, so for most\n"
          "species the genus alone determines the label")
fig1.tight_layout(rect=[0, 0, 1, 0.88])
fig1.savefig(PLOTS_DIR / "fig1_genus_purity.png", dpi=300, bbox_inches="tight")
fig1.savefig(PLOTS_DIR / "fig1_genus_purity.pdf", bbox_inches="tight")
plt.show()


top_composition = (
    genus_stats.sort_values("n_species", ascending=False)
    .head(TOP_N_COMPOSITION)
    .sort_values("n_species", ascending=True)  # ascending so barh puts the largest at the top
    .reset_index(drop=True)
)

y_pos2 = np.arange(len(top_composition))
fig2, ax2 = plt.subplots(figsize=(9, 6.5))
ax2.barh(y_pos2, top_composition["n_positive"], color=CLASS_COLOURS[POS_LABEL], label=POS_LABEL)
ax2.barh(y_pos2, top_composition["n_negative"], left=top_composition["n_positive"],
         color=CLASS_COLOURS[NEG_LABEL], label=NEG_LABEL)
ax2.set_yticks(y_pos2)
ax2.set_yticklabels(top_composition["genus"], fontstyle="italic", fontsize=10)
ax2.set_xlabel("Species")
strip_top_right(ax2)
ax2.legend(handles=[mpatches.Patch(color=CLASS_COLOURS[NEG_LABEL], label=NEG_LABEL),
                     mpatches.Patch(color=CLASS_COLOURS[POS_LABEL], label=POS_LABEL)],
           loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2, frameon=False)
add_title(ax2, f"Class composition of the {TOP_N_COMPOSITION} largest genera",
          f"Model {MODEL} population. Qualifier tokens (e.g. uncultured, unclassified) "
          "are excluded as they are not genera.")
fig2.subplots_adjust(left=0.15, right=0.95, top=0.82, bottom=0.16)
fig2.savefig(PLOTS_DIR / "fig2_genus_composition.png", dpi=300, bbox_inches="tight")
fig2.savefig(PLOTS_DIR / "fig2_genus_composition.pdf", bbox_inches="tight")
plt.show()


def variance_explained(value, group):
    value = np.asarray(value, dtype=float)
    group = np.asarray(group, dtype=object)
    keep = ~pd.isna(group)
    value = value[keep]
    group = group[keep]
    grand = value.mean()
    ss_total = np.sum((value - grand) ** 2)
    if ss_total == 0:
        return np.nan
    tmp = pd.DataFrame({"value": value, "group": group})
    means = tmp.groupby("group")["value"].mean()
    ns = tmp.groupby("group")["value"].size()
    ss_between = np.sum(ns * (means - grand) ** 2)
    return ss_between / ss_total


def perm_null(value, group, n_perm):
    group = np.asarray(group, dtype=object)
    out = np.empty(n_perm)
    for i in range(n_perm):
        out[i] = variance_explained(value, rng.permutation(group))
    return out


groupings = {"genus": model_df["genus"].values, "class": model_df["class"].values}
if has_phylum:
    groupings["phylum"] = model_df["phylum"].values

rows = []
for g_name, grp in groupings.items():
    obs = variance_explained(model_df["value"].values, grp)
    null = perm_null(model_df["value"].values, grp, N_PERM)
    rows.append({
        "grouping": g_name,
        "n_groups": pd.Series(grp).dropna().nunique(),
        "r2_observed": obs,
        "r2_null_mean": null.mean(),
        "r2_excess": obs - null.mean(),
        "p_permutation": (np.sum(null >= obs) + 1) / (N_PERM + 1),
    })
decomposition = pd.DataFrame(rows)
num_cols = decomposition.select_dtypes(include=[np.number]).columns
decomposition[num_cols] = decomposition[num_cols].round(4)
decomposition.to_csv(TABLES_DIR / "variance_decomposition.csv", index=False)
print(f"Variance in {FOCUS_FEATURE} explained by each grouping")
print(decomposition)
print("\nr2_excess is the observed value minus the permutation mean: the")
print("structure remaining once inflation from group count is removed.\n")

overall_iqr_pct = 100 * (np.quantile(model_df["value"], 0.75) - np.quantile(model_df["value"], 0.25))
median_within_iqr_pct = genus_stats.loc[
    genus_stats["n_species"] >= MIN_PER_CLASS, "iqr_pct"
].median()

print(f"Overall IQR: {overall_iqr_pct:.4f}%   "
      f"median within-genus IQR: {median_within_iqr_pct:.4f}% "
      f"({100 * median_within_iqr_pct / overall_iqr_pct:.1f}% of overall)\n")


within_genus_rows = []
for genus, g in model_df.groupby("genus"):
    n_neg = int((g["class"] == NEG_CLASS).sum())
    n_pos = int((g["class"] == POS_CLASS).sum())
    if n_neg < MIN_PER_CLASS or n_pos < MIN_PER_CLASS:
        continue
    med_neg = g.loc[g["class"] == NEG_CLASS, "value"].median()
    med_pos = g.loc[g["class"] == POS_CLASS, "value"].median()
    within_genus_rows.append({
        "genus": genus,
        "n_negative": n_neg,
        "n_positive": n_pos,
        "median_negative": med_neg,
        "median_positive": med_pos,
        "difference": med_pos - med_neg,
    })

within_genus = pd.DataFrame(within_genus_rows)
if len(within_genus):
    within_genus["_n_total"] = within_genus["n_negative"] + within_genus["n_positive"]
    within_genus = (within_genus.sort_values("_n_total", ascending=False)
                     .drop(columns="_n_total"))
within_genus.to_csv(TABLES_DIR / "within_genus_comparison.csv", index=False)

overall_difference = (
    model_df.loc[model_df["class"] == POS_CLASS, "value"].median()
    - model_df.loc[model_df["class"] == NEG_CLASS, "value"].median()
)

n_tested = len(within_genus)
if n_tested:
    n_same = int((np.sign(within_genus["difference"]) == np.sign(overall_difference)).sum())
    n_opposite = int((np.sign(within_genus["difference"]) == -np.sign(overall_difference)).sum())
    n_tied = int((within_genus["difference"] == 0).sum())
else:
    n_same = n_opposite = n_tied = 0

print(f"Within-genus comparison on {FOCUS_FEATURE}")
print(f"Minimum species per class: {MIN_PER_CLASS}")
print(f"Genera testable: {n_tested}")
print(f"Pooled median difference: {overall_difference:.3g}")

sign_test = None
if n_tested > 0:
    print(f"  same direction:     {n_same}")
    print(f"  opposite direction: {n_opposite}")
    print(f"  tied (exactly 0):   {n_tied}")

    n_nontied = n_same + n_opposite
    if n_nontied > 0:
        sign_test = binomtest(n_same, n_nontied, p=0.5)
        ci = sign_test.proportion_ci(confidence_level=0.95, method="exact")
        print(f"Sign test on {n_nontied} non-tied genera: p = {sign_test.pvalue:.3g}")
        print(f"  95% CI: {ci.low:.3g} - {ci.high:.3g}")
        print("\nNOTE: with few testable genera this cannot distinguish 'no effect'")
        print("from 'insufficient power'. Report the interval, not just the p-value.")
print()


counted = model_df.copy()
counted["raw_count"] = counted["value"] * counted["prot_length"]
counted["res_count"] = counted["raw_count"].round().astype(int)
counted["residual"] = (counted["raw_count"] - counted["res_count"]).abs()

max_residual = counted["residual"].max()
print(f"Reconstructed {RESIDUE_NAME} counts")
print(f"  max distance from an integer: {max_residual:.3g}")
if max_residual > 0.01:
    print(f"  WARNING: values are not close to integers. {LENGTH_COL} may not")
    print("  be the denominator used in 4.2; treat the count figures as approximate.")
print(pd.crosstab(counted["res_count"], counted["class_label"]))
print()

count_prop = (
    counted.groupby("res_count")
    .agg(n_species=("class", "size"),
         n_positive=("class", lambda s: int((s == POS_CLASS).sum())))
    .reset_index()
)
count_prop["prop_positive"] = count_prop["n_positive"] / count_prop["n_species"]
count_prop.to_csv(TABLES_DIR / "count_vs_class.csv", index=False)
print(f"{PARTNER} co-occurrence by {RESIDUE_NAME} count")
print(count_prop.assign(prop_positive=count_prop["prop_positive"].round(3)))
print()


mosaic_df = count_prop[count_prop["n_species"] >= MIN_COUNT_GROUP].copy()
mosaic_df["xmax"] = mosaic_df["n_species"].cumsum()
mosaic_df["xmin"] = mosaic_df["xmax"] - mosaic_df["n_species"]
mosaic_df["xmid"] = (mosaic_df["xmin"] + mosaic_df["xmax"]) / 2

overall_prop = (counted["class"] == POS_CLASS).mean()

fig3, ax3 = plt.subplots(figsize=(9, 6))
for _, row in mosaic_df.iterrows():
    width = row["xmax"] - row["xmin"]
    ax3.add_patch(mpatches.Rectangle(
        (row["xmin"], 0), width, row["prop_positive"],
        facecolor=CLASS_COLOURS[POS_LABEL], edgecolor="white", linewidth=1))
    ax3.add_patch(mpatches.Rectangle(
        (row["xmin"], row["prop_positive"]), width, 1 - row["prop_positive"],
        facecolor=CLASS_COLOURS[NEG_LABEL], edgecolor="white", linewidth=1))
    ax3.text(row["xmid"], row["prop_positive"] + 0.045,
              f"{100 * row['prop_positive']:.0f}%",
              ha="center", fontsize=11, fontweight="bold", color="0.15")
    ax3.text(row["xmid"], -0.05, f"{int(row['res_count'])}",
              ha="center", fontsize=11, fontweight="bold", color="0.20")
    ax3.text(row["xmid"], -0.105, f"n={int(row['n_species'])}",
              ha="center", fontsize=8.5, color="0.45")

ax3.axhline(overall_prop, linestyle="--", color="0.20", linewidth=0.5)
ax3.text(0, overall_prop + 0.035, f"dataset average {100 * overall_prop:.0f}%",
          ha="left", fontsize=8, color="0.20")

ax3.set_xlim(0, mosaic_df["xmax"].max())
ax3.set_ylim(-0.14, 1.05)
ax3.set_yticks(np.arange(0, 1.01, 0.25))
ax3.set_yticklabels([f"{int(t * 100)}%" for t in np.arange(0, 1.01, 0.25)])
ax3.set_xticks([])
ax3.set_ylabel(f"Species carrying {PARTNER}")
legend_handles = [mpatches.Patch(color=CLASS_COLOURS[POS_LABEL], label=POS_LABEL),
                   mpatches.Patch(color=CLASS_COLOURS[NEG_LABEL], label=NEG_LABEL)]
ax3.legend(handles=legend_handles, loc="upper center", bbox_to_anchor=(0.5, -0.02),
            ncol=2, frameon=False)
strip_top_right(ax3)
add_title(ax3, f"Fewer {RESIDUE_NAME}s, more {PARTNER}",
          f"Bar width = number of species; height = share carrying {PARTNER}.\n"
          f"Numbers below each bar are the {RESIDUE_NAME} count.")
fig3.tight_layout(rect=[0, 0, 1, 0.86])
fig3.savefig(PLOTS_DIR / "fig3_count_vs_cooccurrence.png", dpi=300, bbox_inches="tight")
fig3.savefig(PLOTS_DIR / "fig3_count_vs_cooccurrence.pdf", bbox_inches="tight")
plt.show()


fig4, ax4 = plt.subplots(figsize=(9, 6))
size_range = mosaic_df["n_species"].max() - mosaic_df["n_species"].min()
sizes = 250 + 900 * (mosaic_df["n_species"] - mosaic_df["n_species"].min()) / (size_range if size_range > 0 else 1)

ax4.plot(mosaic_df["res_count"], 100 * mosaic_df["prop_positive"],
          color=CLASS_COLOURS[POS_LABEL], linewidth=2, zorder=2)
sc = ax4.scatter(mosaic_df["res_count"], 100 * mosaic_df["prop_positive"],
                  s=sizes, c=mosaic_df["n_species"], cmap="Blues",
                  edgecolor="black", linewidth=0.7, zorder=3)
for _, row in mosaic_df.iterrows():
    ax4.annotate(f"{100 * row['prop_positive']:.0f}%",
                 (row["res_count"], 100 * row["prop_positive"]),
                 textcoords="offset points", xytext=(0, 14),
                 ha="center", fontsize=11, fontweight="bold")

ax4.axhline(100 * overall_prop, linestyle="--", color=CLASS_COLOURS[POS_LABEL],
            linewidth=1, alpha=0.7)
ax4.text(mosaic_df["res_count"].min(), 100 * overall_prop - 4,
          f"Dataset average: {100 * overall_prop:.0f}%", fontsize=9, color="0.3")

cbar = fig4.colorbar(sc, ax=ax4)
cbar.set_label("Species count")
ax4.set_xlabel(f"{RESIDUE_NAME.capitalize()} count")
ax4.set_ylabel(f"{PARTNER} (%)")
ax4.set_ylim(0, max(100 * mosaic_df["prop_positive"].max() * 1.35, 20))
strip_top_right(ax4)
add_title(ax4, f"{PARTNER} frequency drops as {RESIDUE_NAME} count increases",
          "Same relationship as Figure 3. Marker size/colour = species count "
          f"at that {RESIDUE_NAME} count.")
fig4.tight_layout(rect=[0, 0, 1, 0.86])
fig4.savefig(PLOTS_DIR / "fig4_count_vs_pcuac_bubble.png", dpi=300, bbox_inches="tight")
fig4.savefig(PLOTS_DIR / "fig4_count_vs_pcuac_bubble.pdf", bbox_inches="tight")
plt.show()


strip_genera = (
    genus_stats[genus_stats["n_species"] >= MIN_SPECIES]
    .sort_values("n_species", ascending=False)
    .head(MAX_GENERA)["genus"]
    .tolist()
)
strip_df = counted[counted["genus"].isin(strip_genera)].copy()

order_stats = (
    strip_df.groupby("genus")
    .apply(lambda g: pd.Series({
        "m": g["res_count"].median(),
        "p": (g["class"] == POS_CLASS).mean(),
    }))
    .reset_index()
    .sort_values(["m", "p"], ascending=[True, False])
)
genus_order = order_stats["genus"].tolist()
strip_df["genus"] = pd.Categorical(strip_df["genus"], categories=genus_order, ordered=True)

if has_phylum:
    strip_df["phylum"] = strip_df["phylum"].fillna("Unassigned")


def plot_strip(ax, sub_df, genus_order_sub):
    y_positions = {g: i for i, g in enumerate(genus_order_sub)}
    for cls_label, colour in CLASS_COLOURS.items():
        sel = sub_df[sub_df["class_label"] == cls_label]
        if sel.empty:
            continue
        y = sel["genus"].map(y_positions).astype(float).values
        y_jit = y + rng.uniform(-0.22, 0.22, size=len(y))
        x_jit = sel["res_count"].values + rng.uniform(-0.18, 0.18, size=len(y))
        ax.scatter(x_jit, y_jit, s=18, alpha=0.65, color=colour,
                   edgecolors="none", label=cls_label)
    ax.set_yticks(range(len(genus_order_sub)))
    ax.set_yticklabels(genus_order_sub, fontstyle="italic", fontsize=9)
    ax.set_ylim(-0.5, len(genus_order_sub) - 0.5)


max_count = int(strip_df["res_count"].max())

if has_phylum:
    # Phyla in the order their genera first appear along genus_order
    phylum_lookup = strip_df.drop_duplicates("genus").set_index("genus")["phylum"]
    phylum_order, seen = [], set()
    for g in genus_order:
        p = phylum_lookup[g]
        if p not in seen:
            seen.add(p)
            phylum_order.append(p)

    genus_by_phylum = {
        p: [g for g in genus_order if phylum_lookup[g] == p]
        for p in phylum_order
    }
    heights = [len(genus_by_phylum[p]) for p in phylum_order]

    fig5, axes5 = plt.subplots(
        len(phylum_order), 1, figsize=(9, 8), sharex=True,
        gridspec_kw={"height_ratios": heights, "hspace": 0.15},
    )
    axes5 = np.atleast_1d(axes5)
    for ax, p in zip(axes5, phylum_order):
        genera_p = genus_by_phylum[p]
        sub_df = strip_df[strip_df["genus"].isin(genera_p)]
        plot_strip(ax, sub_df, genera_p)
        ax.set_ylabel(textwrap.fill(p, 14), rotation=0, ha="right", va="center",
                      fontsize=9, color="0.25")
        strip_top_right(ax)
    axes5[-1].set_xlabel(f"Number of {RESIDUE_NAME} residues in {PROTEIN}")
else:
    fig5, ax5 = plt.subplots(figsize=(9, 8))
    plot_strip(ax5, strip_df, genus_order)
    ax5.set_xlabel(f"Number of {RESIDUE_NAME} residues in {PROTEIN}")
    strip_top_right(ax5)
    axes5 = np.array([ax5])

axes5[0].set_xticks(range(0, max_count + 1))
handles = [Line2D([0], [0], marker="o", linestyle="", color=CLASS_COLOURS[label], label=label)
           for label in CLASS_COLOURS]
fig5.legend(handles=handles, loc="lower center", ncol=2, frameon=False,
            bbox_to_anchor=(0.5, -0.02))

add_title(axes5[0],
          f"{RESIDUE_NAME.capitalize()} count is fixed within lineages",
          f"One point per species ({MAX_GENERA} largest genera). Each genus occupies a single count.\n"
          f"Genera at low counts carry {PARTNER} far more often.")

fig5.tight_layout(rect=[0, 0.03, 1, 0.94])
fig5.savefig(PLOTS_DIR / "fig5_count_by_genus.png", dpi=300, bbox_inches="tight")
fig5.savefig(PLOTS_DIR / "fig5_count_by_genus.pdf", bbox_inches="tight")
plt.show()

def r2_of(g):
    return decomposition.loc[decomposition["grouping"] == g, "r2_observed"].iloc[0]


def excess_of(g):
    return decomposition.loc[decomposition["grouping"] == g, "r2_excess"].iloc[0]


low_counts = count_prop[(count_prop["n_species"] >= MIN_COUNT_GROUP) & (count_prop["res_count"] <= 2)]
high_counts = count_prop[(count_prop["n_species"] >= MIN_COUNT_GROUP) & (count_prop["res_count"] >= 3)]

summary_rows = [
    ("Model", MODEL),
    ("Focus feature", FOCUS_FEATURE),
    ("Species analysed", str(len(model_df))),
    ("Dropped (no resolvable genus)", str(n_unresolved)),
    ("Distinct genera", str(model_df["genus"].nunique())),
    ("Class-pure genera", str(n_pure_genera)),
    ("Species in class-pure genera", str(n_pure_species)),
    ("Variance explained by genus", str(r2_of("genus"))),
    ("excess over permutation null", str(excess_of("genus"))),
    ("Variance explained by class", str(r2_of("class"))),
    ("excess over permutation null", str(excess_of("class"))),
    ("Variance explained by phylum", str(r2_of("phylum")) if has_phylum else "not computed"),
    ("Overall IQR (%)", f"{overall_iqr_pct:.4g}"),
    ("Median within-genus IQR (%)", f"{median_within_iqr_pct:.4g}"),
    ("Within-genus spread, % of overall", f"{100 * median_within_iqr_pct / overall_iqr_pct:.1f}"),
    ("Genera testable within-genus", str(n_tested)),
    ("same direction", str(n_same)),
    ("opposite direction", str(n_opposite)),
    ("tied", str(n_tied)),
    ("Sign test p", f"{sign_test.pvalue:.3g}" if sign_test is not None else "not performed"),
    (f"1-2 {RESIDUE}: % carrying {PARTNER}",
     f"{100 * low_counts['n_positive'].sum() / low_counts['n_species'].sum():.1f}"),
    (f"3+ {RESIDUE}: % carrying {PARTNER}",
     f"{100 * high_counts['n_positive'].sum() / high_counts['n_species'].sum():.1f}"),
]
summary_table = pd.DataFrame(summary_rows, columns=["metric", "value"])
summary_table.to_csv(TABLES_DIR / "genus_structure_summary.csv", index=False)
print(summary_table.to_string(index=False))

session_info = [
    f"Python {sys.version}",
    f"Platform: {platform.platform()}",
    f"numpy {np.__version__}",
    f"pandas {pd.__version__}",
    f"scipy {scipy.__version__}",
    f"matplotlib {matplotlib.__version__}",
]
(TABLES_DIR / "5.2_sessionInfo.txt").write_text("\n".join(session_info))

print(f"\nOutputs written to {OUTPUT_DIR}\n")
