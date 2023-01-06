from rdkit.Chem import Descriptors
from rdkit.Chem.Descriptors import HeavyAtomMolWt
import numpy

# generation of SPRESSO fragments 

for i in $(ls fda_*.sdf)
do
/home/zhenyu/decomposition/spresso_decompose -l "$i" -f ${i:0:-4}_f.sdf -o ${i:0:-4}_o.sdf # done by supervisor as my laptop did not have the right specifications 
done 

# extraction of SPRESSO fragments from file sent by supervisor 
import os

path = "./"

frags_spresso_rdkit = []

for dire in os.listdir():
    if ".sdf" in dire:
        frags_spresso_rdkit.append(Chem.SDMolSupplier(dire))

      
# analysis of SPRESSO fragments 
weight_each_spresso = []


for r in frags_spresso_rdkit:
    weight_each_spresso.append(HeavyAtomMolWt(r))

print("Number of fragments:", len(frags_spresso_rdkit))
print("Average molecular weight of fragments:", sum(weight_each_spresso)/len(frags_spresso_rdkit))
print("Standard deviation:", numpy.std(weight_each_spresso))

import matplotlib.pyplot as plt

plt.style.use('ggplot')
plt.hist(weight_each_spresso, bins=10)
plt.show()

# Calculation of no. of fragments with < 3 hydrogen donors 
h_donors_bm = []

for x in frags_bm_rdkit:
    if Chem.Lipinski.NumHDonors(x) < 3:
        h_donors_bm.append(Chem.Lipinski.NumHDonors(x))

print('No. of fragments with H donors < 3:',len(h_donors_bm))
print('Fraction:', len(h_donors_bm)/len(frags_bm_rdkit))

# Calculation of no. of fragments with < 3 hydrogen acceptors 
h_acceptors_spresso = []

for x in frags_spresso_rdkit:
    if Chem.Lipinski.NumHAcceptors(x) < 3:
        h_acceptors_spresso.append(Chem.Lipinski.NumHAcceptors(x))

print('No. of fragments with H acceptors < 3:',len(h_acceptors_spresso))
print('Fraction:', len(h_acceptors_spresso)/len(frags_spresso_rdkit))

# Calculation of no. of fragments with logP < 3 
logp_spresso = []

for x in frags_spresso_rdkit:
    if Chem.Crippen.MolLogP(x) < 3:
        logp_spresso.append(Chem.Crippen.MolLogP(x))
        
print('No. of fragments with logP < 3:',len(logp_spresso))
print('Fraction:', len(logp_spresso)/len(frags_spresso_rdkit))

# Calculation of no. of fragments with mw < 300
mw_spresso= []

for x in frags_spresso_rdkit:
    if HeavyAtomMolWt(x) < 300:
        mw_spresso.append(HeavyAtomMolWt(x))
        
print('No. of fragments with mw < 300:',len(mw_spresso))
print('Fraction:', len(mw_spresso)/len(frags_spresso_rdkit))

# similarity analysis 

fps_spresso = [MACCSkeys.GenMACCSKeys(x) for x in frags_spresso_rdki]

sim_spresso_maccs = [[] for c in fps_spresso]
for c in range(len(maccs_spresso)):
    for cc in range(len(maccs_spresso)):
        if c == cc:
            sim_spresso_maccs[c].append(1)
        else:
            sim_spresso_maccs[c].append(DataStructs.FingerprintSimilarity(fps_spresso[c], fps_spresso[cc]))

sim_dict_spresso_maccs = {}

for d in range(len(sim_spresso_maccs)):
    data_buffer_spresso_maccs = []
    for dd in range(len(sim_spresso_maccs)):
        # to keep the data that would be above threshold
        
        if d == dd:
            continue
        
        if sim_spresso_maccs[d][dd] < 0.000001:
            data_buffer_spresso_maccs.append(dd)
            
    sim_dict_spresso_maccs[d] = data_buffer_spresso_maccs
    
    if len(data_buffer_spresso_maccs) > 0:
        print("{} | {}".format(d, [str(ddf) for ddf in data_buffer_spresso_maccs]))

