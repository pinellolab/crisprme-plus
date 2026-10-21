#!/usr/bin/env python3
"""Fit the DEPLOYMENT off-target probability calibration artifacts (v2).

Adds, beyond the point map + bootstrap ensemble + per-bin phi:
  * per-REPORT-BAND overdispersion phi (incl. the >=0.2 AGGREGATE) so aggregate
    intervals are computed on the combined band directly, not by summing
    independent sub-bands (cross-band guide-promiscuity correlation is ~+0.66;
    independent combination wrongly inflates the aggregate's lower bound).
  * per-band guide-positive counts + #contributing guides (IND resolution).
  * guide-cluster-bootstrap AUPR CIs (the point AUPR is high-variance).
  * per-held-out-guide (=one novel guide) TOTAL-count coverage per report band
    (the valid single-guide coverage number).

Writes <scorer>_calibration.npz with grid, p_point, p_ens, phi (per-bin, legacy),
and meta carrying phi_bands / band_counts / aupr / loo_coverage. numpy/scipy only."""
import csv, json, numpy as np

EPS=1e-6; B=2000; NBOOT=1000; NBOOT_AUPR=300; COVER_INNER=400
CENTER=(np.arange(B)+0.5)/B
BIN_BAND=np.digitize(CENTER,[0.1,0.2,0.5])             # 0..3 fine bins
# report bands the card uses: label -> (lo, hi)
REPORT_BANDS={"0-0.1":(0.0,0.1),"0.1-0.2":(0.1,0.2),"0.2-0.5":(0.2,0.5),
              "0.5-1":(0.5,1.0),"ge0.2":(0.2,1.0)}

def load():
    g=[];cfd=[];cb=[];lab=[];cbh=[]
    with open("scored.tsv") as s,open("benchmark.tsv") as b:
        rs,rb=csv.reader(s,delimiter="\t"),csv.reader(b,delimiter="\t");next(rs);next(rb)
        for a,c in zip(rs,rb):
            lab.append(int(a[5]));cfd.append(float(a[6]));cb.append(float(a[7]))
            g.append(c[0]);cbh.append(c[11]=="1")
    return np.array(g),np.array(cfd),np.array(cb),np.array(lab,dtype=np.int64),np.array(cbh)

def binof(s): return np.clip((s*B).astype(int),0,B-1)
def tables(scores,labels,guides):
    gids=np.unique(guides); gi={x:i for i,x in enumerate(gids)}; G=len(gids)
    bi=binof(scores); gidx=np.array([gi[x] for x in guides])
    N=np.zeros((G,B));K=np.zeros((G,B)); np.add.at(N,(gidx,bi),1.0); np.add.at(K,(gidx,bi),labels.astype(float))
    return gids,N,K,gidx

def fit_iso(n,k):
    nz=n>0
    if nz.sum()==0: return np.full(B,EPS)
    c=CENTER[nz];v=k[nz]/n[nz];w=n[nz];bv=[];bw=[];bc=[]
    for i in range(len(v)):
        cv=v[i];cw=w[i];cnt=1
        while bv and bv[-1]>=cv:
            pv=bv.pop();pw=bw.pop();pc=bc.pop();cv=(pv*pw+cv*cw)/(pw+cw);cw=pw+cw;cnt+=pc
        bv.append(cv);bw.append(cw);bc.append(cnt)
    f=np.empty(len(v));p=0
    for cv,cw,cnt in zip(bv,bw,bc): f[p:p+cnt]=cv;p+=cnt
    return np.clip(np.interp(CENTER,c,np.clip(f,EPS,1-EPS)),EPS,1-EPS)

def bin_in_range(lo,hi): return (CENTER>=lo)&(CENTER<hi) if hi<1.0 else (CENTER>=lo)
def phi_mom(e,o):
    use=e>1e-6; e=e[use];o=o[use]; d=np.sum(e**2)
    return max(0.0,float(np.sum((o-e)**2-e)/d)) if d>0 and len(e)>1 else 0.0

def aupr_exact(scores,labels):
    o=np.argsort(-scores,kind="mergesort"); y=labels[o]; P=y.sum()
    if P==0: return 0.0
    tp=np.cumsum(y); fp=np.cumsum(1-y); rec=tp/P; prec=tp/(tp+fp)
    return float(np.trapz(prec,rec))

def fit_scorer(name, scores, labels, guides, provenance):
    gids,N,K,gidx=tables(scores,labels,guides); G=len(gids)
    Ntot=N.sum(0); Ktot=K.sum(0)
    p_point=fit_iso(Ntot,Ktot)
    rng=np.random.default_rng(20261020)
    # bootstrap ensemble of maps (guide resample + refit)
    ens=np.empty((NBOOT,B))
    for b in range(NBOOT):
        s=rng.integers(0,G,G); ens[b]=fit_iso(N[s].sum(0),K[s].sum(0))
    # per-guide LOGO per-bin expected contribution (for per-range phi, counts, coverage)
    Econ=np.zeros((G,B))
    for q in range(G):
        f=fit_iso(Ntot-N[q],Ktot-K[q]); Econ[q]=N[q]*f
    # per-report-band phi + guide-positive counts + #contributing guides
    phi_bands={}; band_counts={}
    for lbl,(lo,hi) in REPORT_BANDS.items():
        m=bin_in_range(lo,hi)
        e=Econ[:,m].sum(1); o=K[:,m].sum(1)
        phi_bands[lbl]=phi_mom(e,o)
        band_counts[lbl]=dict(positives=int(o.sum()),
                              guides_with_positive=int((o>0).sum()),
                              guides_contributing=int((e>1e-6).sum()))
    # legacy per-bin phi (fine bins) kept for back-compat
    phi_bin=np.zeros(4)
    for bd in range(4):
        m=BIN_BAND==bd; phi_bin[bd]=phi_mom(Econ[:,m].sum(1),K[:,m].sum(1))
    # guide-cluster-bootstrap AUPR CI
    gsite={i:np.where(gidx==i)[0] for i in range(G)}
    auprs=[]
    for b in range(NBOOT_AUPR):
        s=rng.integers(0,G,G); idx=np.concatenate([gsite[i] for i in s])
        auprs.append(aupr_exact(scores[idx],labels[idx]))
    aupr=dict(point=aupr_exact(scores,labels),
              ci=[float(np.percentile(auprs,5)),float(np.percentile(auprs,95))])
    # per-held-out-guide (one novel guide) TOTAL-count coverage per report band
    loo_cov={}
    for lbl,(lo,hi) in REPORT_BANDS.items():
        m=bin_in_range(lo,hi); ph=phi_bands[lbl]; cov=0
        for q in range(G):
            others=np.array([i for i in range(G) if i!=q])
            draws=np.empty(COVER_INNER)
            for b in range(COVER_INNER):
                s=others[rng.integers(0,len(others),len(others))]
                f=fit_iso(N[s].sum(0),K[s].sum(0)); lam=(N[q][m]*f[m]).sum()
                if ph>1e-9 and lam>0: lam=rng.gamma(1/ph,lam*ph)
                draws[b]=rng.poisson(max(lam,0))
            olo,ohi=np.percentile(draws,[5,95]); o=K[q][m].sum()
            if olo<=o<=ohi: cov+=1
        loo_cov[lbl]=round(cov/G,3)
    prevalence=float(labels.sum())/len(labels)
    meta=dict(scorer=name,n_guides=int(G),n_positives=int(labels.sum()),n_sites=int(len(labels)),
              prevalence=prevalence,band_edges=[0.1,0.2,0.5],report_bands=REPORT_BANDS,
              grid_B=B,n_boot=NBOOT,phi_bands=phi_bands,band_counts=band_counts,
              aupr=aupr,loo_total_count_coverage_nominal90=loo_cov,provenance=provenance,
              note="interval per report band = two-level guide bootstrap over p_ens + Gamma(mean1,var=phi_bands[band]) + Poisson; aggregates use the band's OWN phi (not summed sub-bands). Counts at benchmark prevalence; scale to deployed prior explicitly.")
    np.savez_compressed(f"{name}_calibration.npz",grid=CENTER,p_point=p_point,p_ens=ens,
                        phi=phi_bin,meta=json.dumps(meta))
    print(f"[{name}] G={G} pos={int(labels.sum())} AUPR={aupr['point']:.3f} CI{aupr['ci']}")
    print(f"   phi_bands: "+", ".join(f"{k}={v:.3f}" for k,v in phi_bands.items()))
    print(f"   LOO total-count coverage(nom90): "+", ".join(f"{k}:{int(v*100)}%" for k,v in loo_cov.items()))
    print(f"   band_counts ge0.2: {band_counts['ge0.2']}")

def main():
    g,cfd,cb,lab,cbh=load()
    fit_scorer("cfd",cfd,lab,g,"GUIDE-seq/validated external-eval CSVs (Tsai2015+RefinedTrueOT+Chen2017+Listgarten2018); 359 pos/46 guides; <=1-gap; CFD scores all")
    m=cbh&(cb>=0)
    fit_scorer("crispr_bulge",cb[m],lab[m],g[m],"same benchmark, CRISPR-Bulge held-out in-domain guides; 287 pos/40 guides")

if __name__=="__main__": main()
