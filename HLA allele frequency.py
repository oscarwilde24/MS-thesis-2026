# ============================================================
# HLA Allele Frequency Analysis & Visualization
# Global → India → Regional (Color-Stable, High-Contrast)
# ============================================================

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from pathlib import Path

sns.set(style="whitegrid")

# ============================================================
# 0. INPUT PATHS & CONSTANTS
# ============================================================

BASE_DIR = Path(
    r"D:\OneDrive - Jawaharlal Nehru Centre for Advanced Scientific Research"
    r"\Academic\Year 3\Codes\Allele frequency"
)

HLA_FILES = {
    "A": BASE_DIR / "A" / "HLA_A_AlleleFrequencies_AllPages_FIXED.csv",
    "B": BASE_DIR / "B" / "HLA_B_AlleleFrequencies_AllPages_FIXED.csv",
    "C": BASE_DIR / "C" / "HLA_C_AlleleFrequencies_AllPages_FIXED.csv",
}

TARGET_COVERAGE = 0.80

TOP_COLORS = ["#C40C0C", "#FF6500", "#CC561E"]  # fixed global top-3
GLOBAL_PALETTE = "Blues"
INDIA_PALETTE = "YlGn"
REGION_PALETTE = "BuPu"

FALLBACK_COLOR = "#7A7A7A"

# ============================================================
# 1. STRONG (NON-PASTEL) PALETTE GENERATOR
# ============================================================

def strong_palette(name, n, low=0.45, high=0.95):
    """
    Extract only dark, high-contrast colors from a colormap.
    Avoids washed-out tones on white background.
    """
    cmap = cm.get_cmap(name)
    values = np.linspace(low, high, n)
    return [cmap(v) for v in values]

# ============================================================
# 2. INDIA REGION CLASSIFIER
# ============================================================

def classify_region(population: str) -> str:
    pop = str(population).lower()

    if any(k in pop for k in ["karnataka", "kannada"]):
        return "Karnataka"

    if any(k in pop for k in [
        "andhra", "telugu", "telangana",
        "tamil", "tamil nadu",
        "kerala", "malayalam",
        "ezhava", "nair", "nadar",
        "malabar", "syrian christian"
    ]):
        return "South_India"

    if any(k in pop for k in [
        "delhi", "north", "central",
        "gujarat", "mumbai", "maratha",
        "bengal", "jalpaiguri",
        "northeast", "north east",
        "ucbb"
    ]):
        return "North_India"

    return "Other"

# ============================================================
# 3. LOAD & CLEAN DATA
# ============================================================

def load_and_clean(path: Path, locus: str) -> pd.DataFrame:
    df = pd.read_csv(path, encoding="latin1")

    df = df[["Allele", "Region", "Freq", "Sample Size"]]
    df.columns = ["Allele", "Population", "Freq", "Sample_Size"]

    df["Freq"] = (
        df["Freq"].astype(str)
        .str.replace(r"[^\d.]", "", regex=True)
        .astype(float)
    )

    df["Sample_Size"] = (
        df["Sample_Size"].astype(str)
        .str.replace(",", "", regex=False)
        .astype(float)
    )

    df.dropna(subset=["Freq", "Sample_Size"], inplace=True)

    pattern = rf"^{locus}\*\d{{2}}:\d{{2}}"
    df = df[df["Allele"].str.match(pattern)].copy()

    df["Allele"] = df["Allele"].str.extract(
        rf"^({locus}\*\d{{2}}:\d{{2}})"
    )

    df["Allele_Copies"] = df["Freq"] * 2 * df["Sample_Size"]
    df["Total_Alleles"] = 2 * df["Sample_Size"]

    return df

# ============================================================
# 4. COLLAPSE TO POPULATION LEVEL
# ============================================================

def collapse_population(df: pd.DataFrame, keep_region=False) -> pd.DataFrame:
    cols = ["Allele", "Population"]
    if keep_region:
        cols.append("Region_Group")

    out = (
        df.groupby(cols, as_index=False)
        .agg(
            Allele_Copies=("Allele_Copies", "sum"),
            Total_Alleles=("Total_Alleles", "sum"),
            Sample_Size=("Sample_Size", "sum")
        )
    )

    out["Allele_Frequency"] = out["Allele_Copies"] / out["Total_Alleles"]
    return out

# ============================================================
# 5. GLOBAL STATISTICS
# ============================================================

def global_stats(df: pd.DataFrame):
    af = (
        df.groupby("Allele")["Allele_Copies"].sum()
        / df.groupby("Allele")["Total_Alleles"].sum()
    )

    carrier = 1 - (1 - af) ** 2
    homo = af ** 2

    return af.sort_values(ascending=False), carrier, homo

def minimal_coverage(af, carrier, homo, target):
    cumulative_not = 1.0
    rows = []

    for allele in carrier.sort_values(ascending=False).index:
        cumulative_not *= (1 - carrier[allele])
        coverage = 1 - cumulative_not

        rows.append({
            "Allele": allele,
            "Global_AF": af[allele],
            "Carrier_%": carrier[allele] * 100,
            "Homozygous_%": homo[allele] * 100,
            "Cumulative_Coverage_%": coverage * 100
        })

        if coverage >= target:
            break

    return pd.DataFrame(rows)

# ============================================================
# 6. COLOR MAP GOVERNANCE (STRICT)
# ============================================================

def build_global_color_map(global_af):
    cmap = {}

    top3 = global_af.head(3).index.tolist()
    for allele, color in zip(top3, TOP_COLORS):
        cmap[allele] = color

    remaining = global_af.index.difference(top3)
    colors = strong_palette(GLOBAL_PALETTE, len(remaining))

    for allele, color in zip(remaining, colors):
        cmap[allele] = color

    return cmap, set(global_af.index)

def extend_color_map(base_cmap, used, af, palette):
    cmap = base_cmap.copy()
    new = af.index.difference(used)

    if len(new) == 0:
        return cmap, set()

    colors = strong_palette(palette, len(new))

    for allele, color in zip(new, colors):
        cmap[allele] = color

    return cmap, set(new)

# ============================================================
# 7. TOP-10 BAR PLOT (COLOR-SAFE)
# ============================================================

def plot_top10(series, cmap, title, outfile):
    df = (
        series.sort_values(ascending=False)
        .head(10)
        .rename("Frequency")
        .reset_index()
    )

    df["Color"] = df["Allele"].map(cmap).fillna(FALLBACK_COLOR)

    plt.figure(figsize=(10, 6))
    ax = sns.barplot(
        data=df,
        x="Allele",
        y="Frequency",
        hue="Allele",
        palette=dict(zip(df["Allele"], df["Color"])),
        legend=False
    )

    ax.set_title(title, fontsize=18)
    ax.set_xlabel("Allele", fontsize=14)
    ax.set_ylabel("Allele Frequency", fontsize=14)
    ax.tick_params(axis="x", rotation=45)

    plt.tight_layout()
    plt.savefig(outfile, dpi=300)
    plt.close()

# ============================================================
# 8. MAIN ANALYSIS LOOP
# ============================================================

for locus, path in HLA_FILES.items():
    print(f"\n🧬 Processing HLA-{locus}")

    raw = load_and_clean(path, locus)
    collapsed = collapse_population(raw)

    af, carrier, homo = global_stats(collapsed)
    minimal_coverage(af, carrier, homo, TARGET_COVERAGE)\
        .to_csv(f"HLA_{locus}_Global_80pct_Coverage.csv", index=False)

    cmap, global_set = build_global_color_map(af)

    plot_top10(
        af, cmap,
        f"HLA-{locus} Top 10 Alleles (Global)",
        f"HLA_{locus}_Top10_Global.png"
    )

    india = raw[raw["Population"].str.startswith("India", na=False)].copy()
    india["Region_Group"] = india["Population"].apply(classify_region)
    india_collapsed = collapse_population(india, keep_region=True)

    india_af = (
        india_collapsed.groupby("Allele")["Allele_Copies"].sum()
        / india_collapsed.groupby("Allele")["Total_Alleles"].sum()
    )

    cmap, india_new = extend_color_map(cmap, global_set, india_af, INDIA_PALETTE)
    used = global_set | india_new

    plot_top10(
        india_af, cmap,
        f"HLA-{locus} Top 10 – All India",
        f"HLA_{locus}_Top10_India_All.png"
    )

    for region in ["South_India", "North_India", "Karnataka"]:
        sub = india_collapsed[india_collapsed["Region_Group"] == region]
        if sub.empty:
            continue

        af_r = (
            sub.groupby("Allele")["Allele_Copies"].sum()
            / sub.groupby("Allele")["Total_Alleles"].sum()
        )

        cmap, _ = extend_color_map(cmap, used, af_r, REGION_PALETTE)

        plot_top10(
            af_r, cmap,
            f"HLA-{locus} Top 10 – {region.replace('_', ' ')}",
            f"HLA_{locus}_Top10_{region}.png"
        )

    print(f"✅ Global color assignments locked for HLA-{locus}")
