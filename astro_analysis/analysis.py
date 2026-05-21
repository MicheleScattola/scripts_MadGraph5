import os
import matplotlib.pyplot as plt
import numpy as np
from math import sin,cos,asin,acos
import healpy

def plot(data,xcol,ycol,marker='o',color='red',label=None):
  plt.plot(data[xcol], data[ycol], 'or',color=color,marker=marker,label=label)
  plt.xlabel(xcol)
  plt.ylabel(ycol)
  plt.show()

def angular_distance(ra1,dec1,ra2,dec2):
  return np.degrees(acos(sin(dec1)*sin(dec2)+cos(dec1)*cos(dec2)*cos(ra1-ra2)))


if __name__ == "__main__":
  FILE = '/home/mike/Physics/labParticles/data_files/auger.txt'
  SAG_RA = 266.417
  SAG_DEC = -29.008

  if not os.path.exists(FILE):
    raise FileNotFoundError('Input file not found!')

  data = np.genfromtxt(FILE, names=True)
  labels = data.dtype.names
  #print(labels)

  source_ra = []
  source_dec = []
  theta_max = 30.
  for obj in data:
    theta = angular_distance(np.radians(obj['RA']),np.radians(obj['Dec']),np.radians(SAG_RA),np.radians(SAG_DEC))
    if(abs(theta)<theta_max and abs(theta)>1.):
      source_ra.append(obj['RA'])
      source_dec.append(obj['Dec'])

    
  #print(sources)
  plt.plot(data['RA'],data['Dec'],marker='o',color='blue',alpha=0.3,linestyle='None',label='all sources')
  plt.plot(source_ra,source_dec,marker='o',color='red',alpha=0.7,linestyle='None',label=f'sources d<{theta_max}°')
  plt.plot(SAG_RA,SAG_DEC,marker='*',color='green',label='SAG A*',linestyle='None')
  plt.xlabel('Right ascension')
  plt.ylabel('Declination')
  plt.title('RA vs Dec')
  plt.legend()
  plt.show()  