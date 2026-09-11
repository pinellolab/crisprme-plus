#!/usr/bin/env python
"""One-off extractor: build a tiny self-contained indel-SNP CI fixture from the
e2e chr22 data. Produces a few-KB JSON the brute-force oracle test consumes with
NO large files and NO crispritz. Run from /Users/lp698/Projects/test_crisprme.
"""
import sys, bisect, json, gzip
sys.path.insert(0, '/Users/lp698/Projects/crisprme-plus/PostProcess')
import overlay_indel_snps as ov
import indel_snp_cis as isc

BASE = '/Users/lp698/Projects/test_crisprme/e2e'
ENR  = f'{BASE}/Genomes/hg38_chr22+hg38_1000G/chr22.enriched.fa'
FAKE = f'{BASE}/Genomes/hg38_chr22+hg38_1000G_INDELS/fake_hg38_1000G_chr22/fakechr22.fa'
LOG  = f'{BASE}/Dictionaries/log_indels_hg38_1000G/logchr22.txt'
DICT = f'{BASE}/Dictionaries/dictionaries_hg38_1000G/my_dict_chr22.json'
VCF  = f'{BASE}/VCFs/hg38_1000G/ALL.chr22.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz'
TGT  = f'{BASE}/Results/crisprme-test-out/crispritz_targets/indels_hg38_chr22+hg38_1000G_20bp-NGG-SpCas9.txt_guides.txt_4_1_1.targets.txt'
IUP = {'M':'AC','R':'AG','W':'AT','S':'CG','Y':'CT','K':'GT'}

E = ov.read_enriched_chromosome(ENR)
AMB = ov._IUPAC_AMBIG
with ov._open_text(FAKE) as fh:
    fh.readline(); FSEQ = "".join(l.strip() for l in fh)

# indel records keyed by fake span
recs=[]
with open(LOG) as fh:
    fh.readline()
    for line in fh:
        c=line.rstrip("\n").split("\t")
        if len(c)<6: continue
        fs,fe=map(int,c[5].split(","))
        sp=int(c[0].split("_")[1].split("-")[0])
        p=c[4].split("_"); ref,alt=p[-2],p[-1]; pos=int(p[-2-1]) if False else int(c[4].split("_")[1])
        recs.append((fs,fe,sp,ref,alt,pos,c[4]))
recs.sort(); starts=[r[0] for r in recs]
def find(P):
    i=bisect.bisect_right(starts,P)-1
    return recs[i] if i>=0 and recs[i][0]<=P<recs[i][1] else None

# group targets by indel; pick indels that (a) include the known co-occurrence,
# (b) span both strands, (c) include an indel whose targets have NO usable SNP.
targets_by_indel={}
with open(TGT) as fh:
    fh.readline()
    for row in fh:
        line=row.rstrip("\n").split("\t")
        r=find(int(line[4]))
        if r is None: continue
        targets_by_indel.setdefault(r, []).append(line)

def used_cols(line, r):
    fs,fe,sp,ref,alt,pos,desc=r
    dna,guide,strand=line[2],line[1],line[6]
    real_of=isc.build_offset_to_real(dna,int(line[4]),strand,fs,sp,ref,alt)
    cols=[]
    for j,ch in enumerate(dna):
        gj=guide[j]
        if ch=='-' or gj in '-N' or real_of[j] is None: continue
        if not (0<=real_of[j]<len(E)): continue
        if E[real_of[j]] in AMB: cols.append((j, real_of[j]))
    return cols

# select up to 6 indels: prioritise ones with SNP-overlapping targets + strand variety
pri=[]
for r,tl in targets_by_indel.items():
    has_snp=any(used_cols(t,r) for t in tl)
    strands={t[6] for t in tl}
    pri.append((has_snp, len(strands), r, tl))
pri.sort(key=lambda x:(-x[0], -x[1]))
chosen=[]
seen_plus=seen_minus=False
for has_snp, nstr, r, tl in pri:
    if len(chosen)>=6: break
    chosen.append((r,tl))
# ensure the known co-occurrence indel (pos 22880651) is included
if not any(r[5]==22880651 for r,_ in chosen):
    for r,tl in targets_by_indel.items():
        if r[5]==22880651: chosen.insert(0,(r,tl)); break
chosen=chosen[:6]

# needed VCF indel positions (1-based POS = the log indel POS)
need_pos={r[5] for r,_ in chosen}
maxpos=max(need_pos)
# targeted VCF pass: phased indel GT for the chosen indels
indel_gt={}
sample_cols=None
with gzip.open(VCF,'rt') as fh:
    for line in fh:
        if line.startswith('##'): continue
        if line.startswith('#CHROM'):
            sample_cols=line.rstrip("\n").split("\t")[9:]; continue
        c=line.split("\t")
        pos=int(c[1])
        if pos>maxpos: break
        if pos not in need_pos: continue
        ref=c[3]
        for ai,alt in enumerate(c[4].split(","), start=1):
            if len(ref)==len(alt): continue  # SNP alt at same pos, skip
            key=f"{pos}_{ref}_{alt}"
            gts={}
            for col,g in enumerate(c[9:]):
                import build_indel_genotypes as big
                ng=big.normalize_gt_for_alt(g, ai)
                if ng is not None: gts[sample_cols[col]]=ng
            indel_gt[key]=gts

# needed dict SNP keys (real+1) across chosen indels' target windows
need_dict=set()
for r,tl in chosen:
    for t in tl:
        for j,real in used_cols(t,r):
            need_dict.add(real+1)
# extract those dict entries (stream the 1.5GB json once)
snp_gt={}
with open(DICT) as fh:
    dd=json.load(fh)
for k in list(dd):
    p=int(k.split(",")[1])
    if p in need_dict: snp_gt[p]=dd[k]
del dd

# assemble the fixture
fixture={"chrom":"chr22","indels":[]}
for r,tl in chosen:
    fs,fe,sp,ref,alt,pos,desc=r
    # enriched window for this fake segment's real positions (+ a small margin)
    reals={}
    for k in range(fe-fs):
        rp=ov.map_fake_offset_to_real(k, sp, len(ref), len(alt))
        if rp is not None and 0<=rp<len(E): reals[rp]=E[rp]
    # SNP ref/alt (forward) at each overlapping dict position + genotypes
    snps={}
    for t in tl:
        for j,real in used_cols(t,r):
            dictpos=real+1
            if dictpos not in snp_gt: continue
            fref=FSEQ[fs + (real - sp)] if 0<=fs+(real-sp)<len(FSEQ) else None  # plain fake ref
            # derive ref/alt from the IUPAC + the plain fake base
            code=E[real]; pair=IUP.get(code.upper())
            if not pair: continue
            # forward ref = plain fake base at this fake offset
            fake_off = None
            # find the fake offset within [fs,fe) that maps to this real (search)
            for kk in range(fe-fs):
                if ov.map_fake_offset_to_real(kk, sp, len(ref), len(alt))==real:
                    fake_off=fs+kk; break
            if fake_off is None: continue
            fbase=FSEQ[fake_off].upper()
            aalt = pair[0] if pair[1]==fbase else (pair[1] if pair[0]==fbase else None)
            if aalt is None: continue
            snps[dictpos]={"ref":fbase,"alt":aalt,"gt":snp_gt[dictpos]}
    fixture["indels"].append({
        "start_position":sp, "ref":ref, "alt":alt, "pos":pos, "desc":desc,
        "fake_start":fs, "fake_end":fe, "fake_seq":FSEQ[fs:fe],
        "enriched":{str(k):v for k,v in reals.items()},
        "snps":{str(k):v for k,v in snps.items()},
        "indel_gt":indel_gt.get(f"{pos}_{ref}_{alt}", {}),
        "targets":[{"bulge_type":t[0],"guide":t[1],"dna":t[2],"chrom":t[3],
                    "fake_pos":int(t[4]),"strand":t[6]} for t in tl],
    })

out='/Users/lp698/Projects/crisprme-plus/PostProcess/fixtures/indel_snp_chr22_fixture.json'
import os; os.makedirs(os.path.dirname(out), exist_ok=True)
with open(out,'w') as f: json.dump(fixture,f,indent=1)
sz=os.path.getsize(out)
print(f"fixture: {len(fixture['indels'])} indels, {sum(len(i['targets']) for i in fixture['indels'])} targets, "
      f"{sum(len(i['snps']) for i in fixture['indels'])} SNPs, {sum(len(i['indel_gt']) for i in fixture['indels'])} indel-GT rows")
print(f"size: {sz/1024:.1f} KB -> {out}")
