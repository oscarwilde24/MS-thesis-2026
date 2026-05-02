import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from pathlib import Path
from tqdm import tqdm

# Set visual style for plots - "ticks" provides a clean look with axis markers
sns.set(style="ticks")

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

TOP_COLORS = ["#C40C0C", "#FF6500", "#CC561E"]  # Colors for Top 3 alleles
GLOBAL_PALETTE = "Blues"                         # Palette for remaining global alleles
INDIA_PALETTE = "YlGn"                          # Palette for unique Indian alleles
FALLBACK_COLOR = "#7A7A7A"                      # Default color for missing values

# ============================================================
# 1. HELPERS
# ============================================================
def strong_palette(name, n, low=0.45, high=0.95):
    """Generates a high-contrast color list from a Matplotlib colormap."""
    cmap = cm.get_cmap(name)
    return [cmap(v) for v in np.linspace(low, high, n)]

def load_and_clean(path: Path, locus: str) -> pd.DataFrame:
    """Reads CSV, cleans string-based numbers, and filters for specific HLA locus."""
    df = pd.read_csv(path, encoding="latin1")[["Allele", "Region", "Freq", "Sample Size"]]
    df.columns = ["Allele", "Population", "Freq", "Sample_Size"]
    
    df["Freq"] = df["Freq"].astype(str).str.replace(r"[^\d.]", "", regex=True).astype(float)
    df["Sample_Size"] = df["Sample_Size"].astype(str).str.replace(",", "", regex=False).astype(float)
    df.dropna(subset=["Freq", "Sample_Size"], inplace=True)

    pattern = rf"^{locus}\*\d{{2}}:\d{{2}}"
    df = df[df["Allele"].str.match(pattern)].copy()
    df["Allele"] = df["Allele"].str.extract(rf"^({locus}\*\d{{2}}:\d{{2}})")

    df["Allele_Copies"] = df["Freq"] * 2 * df["Sample_Size"]
    df["Total_Alleles"] = 2 * df["Sample_Size"]
    return df

def collapse_population(df: pd.DataFrame) -> pd.DataFrame:
    """Aggregates data across populations to calculate mean allele frequency."""
    out = df.groupby(["Allele"], as_index=False).agg(
        Allele_Copies=("Allele_Copies", "sum"),
        Total_Alleles=("Total_Alleles", "sum")
    )
    out["Allele_Frequency"] = out["Allele_Copies"] / out["Total_Alleles"]
    return out

# ============================================================
# 2. MERGED PLOTTING FUNCTION
# ============================================================
def plot_merged(ax, af_global, af_india, cmap, locus):
    """Creates a side-by-side bar comparison of Global vs India within one axis."""
    top_global = af_global.sort_values(ascending=False).head(10)
    top_india = af_india.sort_values(ascending=False).head(10)

    gap = 1
    x_global = np.arange(len(top_global))
    x_india = np.arange(len(top_india)) + len(top_global) + gap

    ax.bar(x_global, top_global.values,
           color=[cmap.get(a, FALLBACK_COLOR) for a in top_global.index])

    ax.bar(x_india, top_india.values,
           color=[cmap.get(a, FALLBACK_COLOR) for a in top_india.index])

    xticks = list(x_global) + list(x_india)
    xlabels = list(top_global.index) + list(top_india.index)
    ax.set_xticks(xticks)
    
    # Configure tick label font sizes
    ax.set_xticklabels(xlabels, rotation=45, ha="right", fontsize=20)
    ax.tick_params(axis='y', labelsize=20)

    # Locus-specific Y-axis ranges
    if locus == "B":
        ax.set_ylim(0, 0.15)
    else:
        ax.set_ylim(0, 0.25)

    # Remove the bounding box (top and right spines)
    sns.despine(ax=ax)

    ax.set_xlabel("")
    ax.set_ylabel("")

# ============================================================
# 3. MAIN
# ============================================================
fig, axes = plt.subplots(3, 1, figsize=(16, 18))
loci = ["A", "B", "C"]

for i, locus in enumerate(tqdm(loci, desc="Generating Merged Figure")):
    raw = load_and_clean(HLA_FILES[locus], locus)

    collapsed = collapse_population(raw)
    af_global = (collapsed.set_index("Allele")["Allele_Copies"] /
                 collapsed.set_index("Allele")["Total_Alleles"])

    india = raw[raw["Population"].str.startswith("India", na=False)].copy()
    india_collapsed = collapse_population(india)
    af_india = (india_collapsed.set_index("Allele")["Allele_Copies"] /
                india_collapsed.set_index("Allele")["Total_Alleles"])

    cmap = {}
    top3 = af_global.sort_values(ascending=False).head(3).index
    for allele, color in zip(top3, TOP_COLORS):
        cmap[allele] = color

    remaining = af_global.index.difference(top3)
    for allele, color in zip(remaining, strong_palette(GLOBAL_PALETTE, len(remaining))):
        cmap[allele] = color

    new_alleles = af_india.index.difference(af_global.index)
    for allele, color in zip(new_alleles, strong_palette(INDIA_PALETTE, len(new_alleles))):
        cmap[allele] = color

    # Pass the locus to the plotting function for specific Y-limits
    plot_merged(axes[i], af_global, af_india, cmap, locus)

plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig("HLA_Merged_Global_India.png", dpi=300, bbox_inches='tight')

print("\n✅ Saved: HLA_Merged_Global_India.png")
