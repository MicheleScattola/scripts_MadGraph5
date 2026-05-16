#!/bin/bash

# Stop on error
set -e 


echo ">>> Running on MC_1 : ggH"
python3.12 rdf.py -i /mnt/share/physics_data/4lep/MC/mc_345060.ggH125_ZZ4lep.4lep.root -o ggH

echo ">>> Running on MC_2 : VBF"
python3.12 rdf.py -i /mnt/share/physics_data/4lep/MC/mc_344235.VBFH125_ZZ4lep.4lep.root -o VBFH