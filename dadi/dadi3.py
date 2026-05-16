import json
import os
import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
import seaborn as sns

FOLDER = "/home/mike/Physics/labParticles/modulo3"
FILE = os.path.join(FOLDER, "datasetDadi.json")
P_THRESHOLD = 0.05
OUTPUT_FOLDER = '/home/mike/Physics/labParticles/scripts/dadi'


def analyze_feature(data: pd.DataFrame, feature: str, ntrials: int):
  print(f"\n--- Analyzing feature: {feature} ---\n")
  
  for key, df in data.groupby(feature):
    nentries = len(df)
    print(f"key: {key}, entries = {nentries}")
    if nentries == 0:
      continue

    # check values stats
    counts = df["results"].value_counts().sort_index().reindex(range(1, 7), fill_value=0)
    #print(f"{counts}\n")
    counts_err = [math.sqrt(int(c)) for c in counts]
    # check chisquare
    expected = nentries / 6.0
    chi = stats.chisquare(f_obs=counts.values, f_exp=[expected] * 6)
    pvalue_pre = chi.pvalue
    print(f"Chi-square: statistic={chi.statistic:.3}, pvalue={pvalue_pre:.3}")

    pvalue_penalized = 1.0-pow(1.0-pvalue_pre,ntrials)
    if(pvalue_penalized < P_THRESHOLD):
      print(f"!!! Outlier found !!!")
      print(f"p-value (penalized) = {pvalue_penalized:.3e}")

      out_dir = os.path.join(OUTPUT_FOLDER, f'plots/{feature}')
      os.makedirs(out_dir, exist_ok=True)

      plt.figure()
      x = np.arange(1, 7)
      y = [c/float(nentries)*6-0 for c in counts.values]
      yerr = [math.sqrt(c)/float(nentries)*6.0 for c in counts.values]
      # plot dots with error bars
      plt.errorbar(x, y, yerr=yerr, fmt='o', color='C0', ecolor='C0', capsize=3,label=r'data$\pm\sqrt{counts}$')
      # expected proportion line at 1/6
      expected_prop = 1.0 / 6.0 *6.0
      plt.hlines(expected_prop, xmin=0.5, xmax=6.5, colors='red', linestyles='dotted',label='expected 1/6')
      plt.xticks(x)
      plt.title(f"{feature}:{key} (n={nentries}, p-value_pen. = {pvalue_penalized:.2e})")
      plt.xlabel('Dice outcome')
      plt.ylabel('Frequency (normalized)')
      plt.legend()
      plt.savefig(os.path.join(out_dir, f'{feature}_{key}.png'))
      plt.close()

    print()
    



if __name__ == "__main__":
  with open(FILE, "r", encoding="utf-8") as file:
    data = json.load(file)

  df = pd.DataFrame(data)
  print(f"columns: {list(df.columns)}")
  df = df.explode("results")

  analysis = ["color","texture"]
  print(f"Analyzing features: {analysis}")
  ntrials = 0
  for feature in analysis:
    nunique = df[feature].nunique()
    ntrials += nunique
    print(f"Feature {feature} with {nunique} classes")

  print(f"Using p-value penalization with n_trials = {ntrials}")

  for feature in analysis:
    analyze_feature(df,feature,ntrials)


  
