"""Original NumPy SVD PCA algorithm, retained for manuscript rendering precision."""
from pathlib import Path
import numpy as np
import pandas as pd
root=Path(__file__).resolve().parents[1]
meta=pd.read_csv(root/'data/metadata/samples.csv'); variance=[]
for model in ['OVCAR5','OVCAR8','ID8']:
    mm=meta[meta.cell_line.eq(model)].reset_index(drop=True)
    v=pd.read_csv(root/f'results/differential_expression/{model}/VST.csv').set_index('gene_id')[mm.sample_id].to_numpy()
    sel=np.argsort(-v.var(1,ddof=1))[:1000]
    z=v[sel].T;z-=z.mean(0)
    u,s,vt=np.linalg.svd(z,full_matrices=False);sc=u*s;ev=s*s/(s*s).sum()
    for k in range(len(s)):
        if vt[k,np.argmax(abs(vt[k]))]<0:sc[:,k]*=-1
    pp=mm[['sample_id','label','batch','treatment']].copy();pp['PC1']=sc[:,0];pp['PC2']=sc[:,1]
    pp.to_csv(root/f'results/plot_data/{model}_paired_VST_PCA_scores.csv',index=False)
    variance.append(dict(cell_line=model,PC1_percent=100*ev[0],PC2_percent=100*ev[1]))
pd.DataFrame(variance).to_csv(root/'results/plot_data/VST_PCA_variance.csv',index=False)
(root/'environment/python_versions.txt').write_text(f'NumPy {np.__version__}\npandas {pd.__version__}\n')
