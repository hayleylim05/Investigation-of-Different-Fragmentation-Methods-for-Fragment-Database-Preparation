from rdkit import Chem
import scaffoldgraph as sg
from rdkit import DataStructs
from rdkit.Chem import Descriptors
from rdkit.Chem.Descriptors import HeavyAtomMolWt
import numpy 
from rdkit.Chem import Lipinski 
from rdkit.Chem import MACCSkeys
import matplotlib.pyplot as plt

cdk2mols = Chem.SDMolSupplier('fda.sdf') 

# generation of Bemis Murcko fragments 
frags_bm_rdkit = []

for mol in cdk2mols:
    frags_bm_rdkit.append(sg.get_all_murcko_fragments(mol))

fps_bm = [Chem.RDKFingerprint(b) for b in frags_bm_rdkit]

frags_bm_rdkit = set()

for mol in cdk2mols:
    frags_bm_rdkit.update(sg.get_all_murcko_fragments(mol))
  
# analysis of Bemis Murcko fragments 
  
weight_each_bm = []


for r in frags_bm_rdkit:
    weight_each_bm.append(HeavyAtomMolWt(r))

# Calculation of no. of fragments, avg mw of fragments and the standard deviation
print("Number of fragments:", len(frags_bm_rdkit))
print("Average molecular weight of fragments:", sum(weight_each_bm)/len(weight_each_bm))
print("Standard deviation:", numpy.std(weight_each_bm))

plt.style.use('ggplot')
plt.hist(weight_each_bm, bins=10)
plt.show() # histogram of mw of fragments 

# Calculation of no. of fragments with < 3 hydrogen donors 
h_donors_bm = []

for x in frags_bm_rdkit:
    if Chem.Lipinski.NumHDonors(x) < 3:
        h_donors_bm.append(Chem.Lipinski.NumHDonors(x))

print('No. of fragments with H donors < 3:',len(h_donors_bm))
print('Fraction:', len(h_donors_bm)/len(frags_bm_rdkit))

# Calculation of no. of fragments with < 3 hydrogen acceptors 
h_acceptors_recap = []

for x in frags_bm_rdkit:
    if Chem.Lipinski.NumHAcceptors(x) < 3:
        h_acceptors_recap.append(Chem.Lipinski.NumHAcceptors(x))

print('No. of fragments with H acceptors < 3:',len(h_acceptors_recap))
print('Fraction:', len(h_acceptors_recap)/len(frags_bm_rdkit))

# Calculation of no. of fragments with logP < 3 
logp_bm = []

for x in frags_bm_rdkit:
    if Chem.Crippen.MolLogP(x) < 3:
        logp_bm.append(Chem.Crippen.MolLogP(x))
        
print('No. of fragments with logP < 3:',len(logp_bm_bm))
print('Fraction:', len(logp_bm)/len(frags_bm_rdkit))

# Calculation of no. of fragments with mw < 300
mw_bm= []

for x in frags_bm_rdkit:
    if HeavyAtomMolWt(x) < 300:
        mw_bm.append(HeavyAtomMolWt(x))
        
print('No. of fragments with mw < 300:',len(mw_bm))
print('Fraction:', len(mw_bm)/len(frags_bm_rdkit))

      
# similarity analysis 

maccs_bm = [MACCSkeys.GenMACCSKeys(x) for x in frags_bm_rdkit] # generation of MACCS keys 

# generation of matrix to calculate Tanimoto similarity of one Bemis Murcko fragment against all other Bemis Murcko fragments 
sim_bm_maccs = [[] for a in maccs_bm]
for a in range(len(maccs_bm)):
    for aa in range(len(maccs_bm)):
        if a == aa:
            sim_bm_maccs[a].append(1)
        else:
            sim_bm_maccs[a].append(DataStructs.FingerprintSimilarity(maccs_bm[a], maccs_bm[aa]))

sim_dict_bm_maccs = {} # creation of dictionary to find out what pairs of fragments fall above the threshold of 0.85

for d in range(len(sim_bm_maccs)):
    data_buffer_bm_maccs = []
    for dd in range(len(sim_bm_maccs)):
        # to keep the data that would be above threshold
        
        if d == dd:
            continue
        
        if sim_bm_maccs[d][dd] > 0.85: 
            data_buffer_bm_maccs.append(dd)
            
    sim_dict_bm_maccs[d] = data_buffer_bm_maccs
    
    if len(data_buffer_bm_maccs) > 0:
        print("{} | {}".format(d, [str(ddf) for ddf in data_buffer_bm_maccs]))