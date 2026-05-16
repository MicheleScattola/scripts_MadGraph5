## $H \to Z Z^* \to 4l$

#### Requirements:
-  one vertex with two associated ID tracks with transverse momentum $p_T > 500 MeV$
-  muons: $p_T > 5 GeV$ , $|\eta| < 2.7$, **except** calorimeter-tagged: $p_T>15GeV$, $|\eta|<0.1$
-  muons: $|d_O|/\sigma_{d_0} < 5$
-  electrons: $|d_0|/\sigma_{d_0} < 3$
-  electrons: $E_T>7GeV$, $|\eta|<2.47$
-  jets (anti-kt): $R=0.4$, $p_T>30GeV$, $|\eta|<4.5$
-  b-tagging: works for $|\eta|<2.5$ with 70% efficiency

#### Selection:
Same flavor - opposite charge (SFOC) lepton pairs $\Rightarrow$ two pairs ordered by invariant mass closest to Higgs: $m_{1,2}$ *leading* , $m_{3,4}$ *sub-leading*.

Three leading leptons are required to satisfy $$p_T>20,15,10 \text{ GeV}$$

- n_leptons = 4
- lepton sum charge*flavor = 0
- DR(l1,l2) > 0.10

#### Mar 5 Mag
weight = evt_weight * 1000 * lumi * xsec / sum_weight