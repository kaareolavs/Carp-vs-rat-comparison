import sys, warnings; warnings.filterwarnings('ignore')
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
from scipy import stats
from data import Tissue; from stats import sex_comparison, bh
from sampling import load
S=load()
res={}
for tn in ['Heart','Brain','Liver']:
    x=Tissue(tn); R=sex_comparison(x)
    comp=[c for c in x.compounds if R['keep'][c]]
    Y=x.log2_imputed(comp)
    info=S[(S.organ==tn)&S.hmt.notna()].set_index('hmt')
    males=[s for s in x.cols('Rat') if x.samples.loc[s,'sex']=='M' and s in info.index and pd.notna(info.loc[s,'weight'])]
    g=x.samples.loc[males,'group']
    out=[]
    for proxy in ['num','weight']:
        v=info.loc[males,proxy].astype(float); vr=v-v.groupby(g).transform('mean')
        sl=[];pp=[]
        for k in comp:
            y=Y.loc[k,males]; yr=y-y.groupby(g).transform('mean')
            r,p=stats.pearsonr(vr,yr); sl.append(r); pp.append(p)
        res[(tn,proxy)]=pd.DataFrame({'r':sl,'p':pp},index=comp)
    fm0=R['sex_at']['Ctr'].loc[comp,'effect']; fm30=R['sex_at']['30min'].loc[comp,'effect']
    print(f'\n## {tn}: males n={len(males)}, groups {g.value_counts().to_dict()}')
    for proxy in ['num','weight']:
        a=res[(tn,proxy)]
        print(f'  {proxy:6}: metabolites with p<0.05 vs {proxy}: {(a.p<0.05).sum()}/{len(a)} (expected ~{0.05*len(a):.0f}); q<0.1: {(bh(a.p)<0.1).sum()}')
        for lab,fm in [('F/M ctrl',fm0),('F/M 30min',fm30)]:
            ok=fm.notna()&a.r.notna()
            rho,p=stats.spearmanr(a.r[ok],fm[ok]); print(f'     corr(male {proxy}-trend, {lab}) rho={rho:+.2f} p={p:.3f}')
    a=res[(tn,'num')]
    top=a[a.p<0.05].sort_values('r')
    print('  sampling-order trend p<0.05:', ', '.join(f"{x.meta.loc[k,'label']}({r:+.2f}; F/M ctrl {fm0[k]:+.1f})" for k,r in top.r.items()))
    for k in ['Hydroxyproline','Ornithine','Carnitine','Val','Phe','Tyr']:
        if k in comp: print(f'    {k:15} r(num)={res[(tn,"num")].loc[k,"r"]:+.2f} p={res[(tn,"num")].loc[k,"p"]:.2f}  r(wt)={res[(tn,"weight")].loc[k,"r"]:+.2f}  F/M ctrl {2**fm0[k]:.2f} 30 {2**fm30[k]:.2f}')
