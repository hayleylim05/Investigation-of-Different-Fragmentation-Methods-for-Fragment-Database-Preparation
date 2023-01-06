cdk2mols = Chem.SDMolSupplier('fda.sdf')

from rdkit import Chem
from rdkit.Chem import BRICS
import pandas as pd
from rdkit import DataStructs
from rdkit.Chem import Descriptors
from rdkit.Chem.Descriptors import ExactMolWt
from rdkit.Chem.Descriptors import HeavyAtomMolWt
import numpy 
from rdkit.Chem import Lipinski
from rdkit.Chem import MACCSkeys

# generation of BRICS fragments 

frags_brics_rdkit = []

for m in cdk2mols:
    pieces = BRICS.BRICSDecompose(m)
    frags_brics_rdkit.append(pieces)

# analysis of BRICS fragments 

# Calculation of no. of fragments, avg molecular weight, standard deviation and plot histogram 

weight_each_brics = []

for r in frags_brics_rdkit:
    weight_each_brics.append(HeavyAtomMolWt(r)) 

print("Number of fragments:", len(frags_brics_rdkit))
print("Average molecular weight of fragments:", sum(weight_each_brics)/len(weight_each_brics))
print("Standard deviation:", numpy.std(weight_each_brics))

import matplotlib.pyplot as plt

plt.style.use('ggplot')
plt.hist(weight_each_brics, bins=10)
plt.show() # plotting of the histogram 


# Calculation of no. of fragments with < 3 hydrogen donors 
h_donors_brics = []

for x in frags_brics_rdkit:
    if Chem.Lipinski.NumHDonors(x) < 3:
        h_donors_brics.append(Chem.Lipinski.NumHDonors(x))

print('No. of fragments with H donors < 3:',len(h_donors_brics))
print('Fraction:', len(h_donors_brics)/len(frags_brics_rdkit))

# Calculation of no. of fragments with < 3 hydrogen acceptors 
h_acceptors_brics = []

for x in frags_brics_rdkit:
    if Chem.Lipinski.NumHAcceptors(x) < 3:
        h_acceptors_brics.append(Chem.Lipinski.NumHAcceptors(x))

print('No. of fragments with H acceptors < 3:',len(h_acceptors_brics))
print('Fraction:', len(h_acceptors_brics)/len(frags_brics_rdkit))

# Calculation of no. of fragments with logP < 3 
logp_brics = []

for x in frags_brics_rdkit:
    if Chem.Crippen.MolLogP(x) < 3:
        logp_brics.append(Chem.Crippen.MolLogP(x))
        
print('No. of fragments with logP < 3:',len(logP_brics))
print('Fraction:', len(logp_brics)/len(frags_brics_rdkit))

# Calculation of no. of fragments with mw < 300
mw_brics = []

for x in frags_brics_rdkit:
    if HeavyAtomMolWt(x) < 300:
        mw_brics.append(HeavyAtomMolWt(x))
        
print('No. of fragments with mw < 300:',len(mw_brics))
print('Fraction:', len(mw_brics)/len(frags_brics_rdkit))

# similarity fragments 

fps_brics = [MACCSkeys.GenMACCSKeys(x) for x in frags_brics_rdkit] # generation of MACCS keys 

# creating a matrix to calculate the similarity of all RECAP fragments against all other RECAP fragments 

sim_brics = [[] for c in fps_brics]
for c in range(len(fps_brics)):
    for cc in range(len(fps_brics)):
        if c == cc:
            sim_brics[c].append(1)
        else:
            sim_brics[c].append(DataStructs.FingerprintSimilarity(fps_brics[c], fps_brics[cc]))

# creating a dictionary to print out what pairs of molecules had what tanimoto similarity 

sim_dict_brics = {}

for d in range(len(sim_brics)):
    data_buffer_brics = []
    for dd in range(len(sim_brics)):
        
        if d == dd:
            continue
        
        if sim_brics[d][dd] < 0.001:  
            data_buffer_brics.append(dd)
            
    sim_dict_brics[d] = data_buffer_brics
    
    if len(data_buffer_brics) > 0:
        print("{} | {}".format(d, [str(ddf) for ddf in data_buffer_brics]))