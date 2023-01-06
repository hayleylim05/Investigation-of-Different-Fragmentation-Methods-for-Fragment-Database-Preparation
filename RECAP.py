from rdkit import Chem
from rdkit.Chem import Recap
import pandas as pd
from rdkit import DataStructs
from rdkit.Chem import Descriptors
from rdkit.Chem.Descriptors import ExactMolWt
from rdkit.Chem.Descriptors import HeavyAtomMolWt
import numpy 
from rdkit.Chem import Lipinski
from rdkit.Chem import MACCSkeys

cdk2mols = Chem.SDMolSupplier('fda.sdf')

# recap fragments + fingerprints 

frags_recap_smiles = []
frags_recap_rdkit = []

for mol in cdk2mols:
    frags_recap = Recap.RecapDecompose(mol)
    frags_recap_smiles.append(sorted(frags_recap.GetLeaves()))
    for x in sorted(frags_recap.GetLeaves()):
        frags_recap_rdkit.append(Chem.MolFromSmiles(x))

# analysis of fragments 

# Calculation of no. of fragments, avg molecular weight, standard deviation and plot histogram 

weight_each_recap = []

for r in frags_recap_rdkit:
    weight_each_recap.append(HeavyAtomMolWt(r)) # creates a list of all individual weights of fragments

print("Number of fragments:", len(frags_recap_rdkit))
print("Average molecular weight of fragments:", sum(weight_each_recap)/len(weight_each_recap))
print("Standard deviation:", numpy.std(weight_each_recap))

import matplotlib.pyplot as plt

plt.style.use('ggplot')
plt.hist(weight_each_recap, bins=10)
plt.show() # displays the histogram 

# Calculation of no. of fragments with < 3 hydrogen donors 
h_donors_recap = []

for x in frags_recap_rdkit:
    if Chem.Lipinski.NumHDonors(x) < 3:
        h_donors_recap.append(Chem.Lipinski.NumHDonors(x))

print('No. of fragments with H donors < 3:',len(h_donors_recap))
print('Fraction:', len(h_donors_recap)/len(frags_recap_rdkit))

# Calculation of no. of fragments with < 3 hydrogen acceptors 
h_acceptors_recap = []

for x in frags_recap_rdkit:
    if Chem.Lipinski.NumHAcceptors(x) < 3:
        h_acceptors_recap.append(Chem.Lipinski.NumHAcceptors(x))

print('No. of fragments with H acceptors < 3:',len(h_acceptors_recap))
print('Fraction:', len(h_acceptors_recap)/len(frags_recap_rdkit))

# Calculation of no. of fragments with logP < 3 
logp_recap = []

for x in frags_recap_rdkit:
    if Chem.Crippen.MolLogP(x) < 3:
        logp_recap.append(Chem.Crippen.MolLogP(x))
        
print('No. of fragments with logP < 3:',len(logP_recap))
print('Fraction:', len(logp_recap)/len(frags_recap_rdkit))

# Calculation of no. of fragments with mw < 300
mw_recap = []

for x in frags_recap_rdkit:
    if HeavyAtomMolWt(x) < 300:
        mw_recap.append(HeavyAtomMolWt(x))
        
print('No. of fragments with mw < 300:',len(mw_recap))
print('Fraction:', len(mw_recap)/len(frags_recap_rdkit))

# similarity analysis 

fps_recap = [MACCSkeys.GenMACCSKeys(x) for x in frags_recap_rdkit] # generation of MACCS keys 

# creating a matrix to calculate the similarity of all RECAP fragments against all other RECAP fragments 

sim_recap = [[] for c in fps_recap]
for c in range(len(fps_recap)):
    for cc in range(len(fps_recap)):
        if c == cc:
            sim_recap[c].append(1)
        else:
            sim_recap[c].append(DataStructs.FingerprintSimilarity(fps_recap[c], fps_recap[cc]))

# creating a dictionary to print out what pairs of molecules had what tanimoto similarity 

sim_dict_recap = {}

for d in range(len(sim_recap)):
    data_buffer_recap = []
    for dd in range(len(sim_recap)):
        
        if d == dd:
            continue
        
        if sim_recap[d][dd] < 0.001: 
            data_buffer_recap.append(dd)
            
    sim_dict_recap[d] = data_buffer_recap
    
    if len(data_buffer_recap) > 0:
        print("{} | {}".format(d, [str(ddf) for ddf in data_buffer_recap]))



