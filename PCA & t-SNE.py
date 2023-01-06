recap_smiles = []
brics_smiles = []
bm_smiles = []
spresso_smiles = []

for a in frags_recap_rdkit:
    recap_smiles.append(Chem.MolToSmiles(a))

for b in frags_brics_rdkit:
    brics_smiles.append(Chem.MolToSmiles(b))

for c in frags_bm_rdkit:
    bm_smiles.append(Chem.MolToSmiles(c))

for d in frags_spresso_rdkit:
    spresso_smiles.append(Chem.MolToSmiles(d))

# creation of df 
  
master_list = [ ]


all_frags = [['recap', recap_smiles],['brics', brics_smiles],['bm', bm_smiles],['spresso', spresso_smiles]]


for frag_name, frags in all_frags:
    for frag in frags:
        master_list.append( [ frag, frag_name ] )


import pandas as pd

df = pd.DataFrame(master_list, columns = [ 'Fragments', 'Method' ])

# taken from https://github.com/PatWalters/workshop.git 
# (Walters P, 2_visualizing_chemical_space.ipynb, 2019, https://github.com/PatWalters/workshop/blob/95513578c5b46cb30c3e450e38dba8f5bd84507b/predictive_models/2_visualizing_chemical_space.ipynb)
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
import numpy as np
from tqdm import tqdm
import seaborn as sns
import matplotlib.pyplot as plt


def fp_list_from_smiles_list(smiles_list,n_bits=2048):
    fp_list = []
    for smiles in tqdm(smiles_list):
        mol = Chem.MolFromSmiles(smiles)
        fp_list.append(fp_as_array(mol,n_bits))
    return fp_list

def fp_as_array(mol,n_bits=2048):
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=n_bits)
    arr = np.zeros((1,), int)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr

fp_list = fp_list_from_smiles_list(df.Fragments)

# plotting of PCA 
from sklearn.decomposition import PCA
pca = PCA(n_components=2)
crds = pca.fit_transform(fp_list) 

ax = sns.scatterplot(data=crds_df.query("Method == 'recap'"),x="PC_1",y="PC_2",color='red')
ax = sns.scatterplot(data=crds_df.query("Method == 'brics'"),x="PC_1",y="PC_2",color='blue')
ax = sns.scatterplot(data=crds_df.query("Method == 'bm'"),x="PC_1",y="PC_2",color='brown')
ax = sns.scatterplot(data=crds_df.query("Method == 'spresso'"),x="PC_1",y="PC_2",color='orange')


_ = plt.legend(labels=['RECAP', 'BRICS','BM','SPRESSO'])

# plotting of t-SNE
pca = PCA(n_components=50)
crds = pca.fit_transform(fp_list) 

from sklearn.manifold import TSNE
crds_embedded = TSNE(n_components=2).fit_transform(crds)

tsne_df = pd.DataFrame(crds_embedded,columns=["X","Y"])
tsne_df['Method'] = df['Method']


ax = sns.scatterplot(data=tsne_df.query("Method == 'spresso'"),x="X",y="Y",color='orange')
ax = sns.scatterplot(data=tsne_df.query("Method == 'bm'"),x="X",y="Y",color='brown')
ax = sns.scatterplot(data=tsne_df.query("Method == 'brics'"),x="X",y="Y",color='blue')
ax = sns.scatterplot(data=tsne_df.query("Method == 'recap'"),x="X",y="Y",color='red')

_ = plt.legend(labels=['SPRESSO','BM','BRICS','RECAP'])