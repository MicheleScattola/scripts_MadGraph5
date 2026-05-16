import os
import json
import pandas as pd 
import scipy.stats as stats
from collections import Counter
import matplotlib.pyplot as plt 
FOLDER = '/home/mike/Physics/labParticles/modulo3'
OUTPUT_FOLDER = os.path.join(FOLDER,'plots')
FILE = os.path.join(FOLDER,'datasetDadi.json')
EXP_MEAN = (sum(range(1,7))/6)

class dado:

  def __init__(self,name,color,weight,texture,material,values):
    self.name = name
    self.color = color
    self.weight = weight
    self.texture = texture
    self.material = material
    self.values = values
    self.chi2 = None
    self.pvalue = None

    
def cleanup(dadi,threshold):
  return [d for d in dadi if d.pvalue <= threshold]

# MAIN
if __name__ == "__main__":

  with open(FILE,'r') as file:
    
    data = json.load(file)

    print(f"columns: {list(data[0].keys()) if data else []}")
    df = pd.DataFrame(data=data)
    lengths = df['results'].map(len)
    print(f"unique results sizes: {sorted(lengths.unique().tolist())}")
    
    dadi = []

    for name, db in df.groupby('name'):
      values = db['results'].explode().tolist()
      first = db.iloc[0]
      d = dado(
        name=name,
        color=first['color'],
        weight=first['weight'],
        texture=first['texture'],
        material=first['material'],
        values=values
      )
      dadi.append(d)

    print(f"initialized {len(dadi)} dado objects")

    for d in dadi:
      frequencies = Counter(d.values)
      observed = [frequencies.get(face, 0) for face in range(1, 7)] # get the value of 'face' , in case of error return 0
      expected = [len(d.values) / 6.0] * 6 # [] of size *6
      chi2, pvalue = stats.chisquare(f_obs=observed, f_exp=expected)
      d.chi2 = float(chi2)
      d.pvalue = float(pvalue)

    print(f"computed chi-square and p-value for {len(dadi)} objects")

    dadi.sort(key=lambda x: x.pvalue)
    print(f"{len(dadi)} dices before cleanup ...")
    dadi = cleanup(dadi,0.05)

    print(f"{len(dadi)} dices after cleanup ...")

    print([f"{dadi[i].pvalue:.2e}" for i in range(6)])
    print([f"{dadi[i].color}" for i in range(6)])

    plt.figure()
    fig, ax = plt.subplots(2,3)
    ax[0,0].hist2d([d.texture for d in dadi],[d.color for d in dadi])
    ax[0,1].hist2d([d.material for d in dadi],[d.color for d in dadi])
    ax[1,0].hist2d([d.texture for d in dadi],[d.weight for d in dadi])
    ax[1,1].hist2d([d.material for d in dadi],[d.weight for d in dadi])
    ax[2,0].hist2d([d.texture for d in dadi],[d.material for d in dadi])
    ax[2,1].hist2d([d.color for d in dadi],[d.weight for d in dadi])

    plt.show()

