import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

def smart_read(file):
    try:
        return pd.read_excel(file, engine="openpyxl")
    except Exception:
        try:
            return pd.read_csv(file, sep="\t")
        except:
            return pd.read_csv(file, sep="\s+")

def extract_total(df):
    # Look for a column likely containing totals
    col_candidates = [c for c in df.columns if any(k in c.lower() for k in ["total", "count", "num", "number"])]
    if col_candidates:
        col = col_candidates[0]
        # Try to find "Total" row if it exists
        total_row = df[df.astype(str).apply(lambda x: x.str.contains("total", case=False)).any(axis=1)]
        if not total_row.empty:
            val = total_row[col].astype(float).values[0]
        else:
            val = df[col].astype(float).sum()
        return int(val)
    else:
        # fallback: last numeric value in the file
        nums = df.select_dtypes(include="number").values.flatten()
        return int(nums[-1]) if len(nums) else 0

# --- Load files ---
snp = smart_read("03.snp.stat.xls")
indel = smart_read("04indel.stat.xls")

# --- Extract totals ---
snp_total = extract_total(snp)
indel_total = extract_total(indel)

# --- Sanity check ---
if snp_total > 5_000_000:  # unrealistic for GBS
    print(f"⚠️ Warning: SNP count unusually high ({snp_total}). Double-check source file format.")
if indel_total > 5_000_000:
    print(f"⚠️ Warning: INDEL count unusually high ({indel_total}). Double-check source file format.")

# --- Create summary table ---
summary = pd.DataFrame({
    "Variant Type": ["SNPs", "INDELs"],
    "Count": [snp_total, indel_total]
})

# --- Save ---
summary.to_csv("Faba_variant_summary.tsv", sep="\t", index=False)
print("✅ Corrected Faba Bean Variant Summary")
print(summary.to_string(index=False))

# --- Plot ---
plt.figure(figsize=(6,4))
sns.barplot(x="Variant Type", y="Count", hue="Variant Type", data=summary,
            palette=["#1f77b4", "#2ca02c"], legend=False)
plt.title("Faba Bean Variant Composition", fontsize=12, fontweight="bold")
plt.ylabel("Variant Count", fontsize=11, fontweight="bold")
plt.xlabel("")
plt.xticks(fontsize=10, fontweight="bold")
plt.yticks(fontsize=9, fontweight="bold")

for i, row in summary.iterrows():
    plt.text(i, row["Count"], f"{row['Count']:,}", ha='center', va='bottom', fontweight='bold', fontsize=10)

plt.tight_layout()
plt.savefig("Faba_variant_summary.png", dpi=300, bbox_inches="tight")
plt.close()
print("\n📄 Updated outputs:")
print(" - Faba_variant_summary.tsv")
print(" - Faba_variant_summary.png")
