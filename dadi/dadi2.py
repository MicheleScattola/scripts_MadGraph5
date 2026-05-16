import json
import os

import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as stats
import seaborn as sns

FOLDER = "/home/mike/Physics/labParticles/modulo3"
FILE = os.path.join(FOLDER, "datasetDadi.json")
P_THRESHOLD = 0.05
OUTPUT_FOLDER = '/home/mike/Physics/labParticles/scripts/dadi'

def build_summary(df: pd.DataFrame, ntrials: int) -> pd.DataFrame:
  meta = df[["name", "color", "weight", "texture", "material"]].copy()

  exploded = df[["name", "results"]].explode("results")
  counts = (
    exploded.groupby(["name", "results"])
      .size() # counts how many elements
      .unstack("results",fill_value=0) # pivots table from multi-index list to matrix
      .reindex(columns=range(1, 7), fill_value=0) # re-arrange columns and fill possible gaps
  )

  n = counts.sum(axis=1)
  expected = n / 6.0
  chi2 = ((counts.sub(expected, axis=0) ** 2).div(expected, axis=0)).sum(axis=1) # axis=0 sum/sub over columns, axis=1 over rows
  chi2_values = chi2.to_numpy()
  pvalue_values = stats.chi2.sf(chi2_values, df=5) # sf is the survival function P(X>=x)
  pval_pre = pd.Series(pvalue_values, index=chi2.index)
  pval_post = 1. - pow(1. - pval_pre, ntrials)

  stats_df = pd.DataFrame(
    {
      "chi2": chi2.values,
      "pvalue": pval_pre.values,
      "pvalue_adj": pval_post.values,
    },
    index=chi2.index,
  ).reset_index(names="name")
  return meta.merge(stats_df, on="name", how="inner")

def plot_correlations(df: pd.DataFrame) -> None:
  fig, axes = plt.subplots(2, 2, figsize=(13, 10))

  texture_color = pd.crosstab(df["texture"], df["color"])
  sns.heatmap(texture_color, cmap="flare", ax=axes[0, 0])
  axes[0, 0].set_title("Texture vs Color")

  material_color = pd.crosstab(df["material"], df["color"])
  sns.heatmap(material_color, cmap="flare", ax=axes[0, 1])
  axes[0, 1].set_title("Material vs Color")

  texture_material = pd.crosstab(df["texture"], df["material"])
  sns.heatmap(texture_material, cmap="flare", ax=axes[1, 0])
  axes[1, 0].set_title("Texture vs Material")

  weight_bins = pd.cut(df["weight"], bins=12)
  texture_weight = pd.crosstab(df["texture"], weight_bins)
  sns.heatmap(texture_weight, cmap="flare", ax=axes[1, 1])
  axes[1, 1].set_title("Texture vs Weight")
  axes[1, 1].set_xlabel("weight bin")

  fig.tight_layout()
  plt.savefig(os.path.join(OUTPUT_FOLDER,'plots/corr.png'))
  #plt.show() 


if __name__ == "__main__":
  with open(FILE, "r", encoding="utf-8") as file:
    data = json.load(file)

  df = pd.DataFrame(data)
  print(f"columns: {list(df.columns)}")
  print(f"unique results sizes: {sorted(df['results'].map(len).unique().tolist())}")

  summary = build_summary(df,len(df)).sort_values("pvalue", ascending=True)
  filtered = summary[summary["pvalue"] <= P_THRESHOLD].copy()

  print(f"{len(summary)} dices before cleanup ...")
  print(f"{len(filtered)} dices after cleanup ...")
  top_n = min(6, len(filtered))
  print([f"{v:.2e}" for v in filtered["pvalue"].head(top_n)])
  print(filtered["color"].head(top_n).tolist())

  if filtered.empty:
    print("No dice below threshold; skip plotting.")
  else:
    sns.set_theme(style="whitegrid")
    plot_correlations(filtered)

