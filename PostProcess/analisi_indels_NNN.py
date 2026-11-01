#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Semi-Brute Force Calculation of Samples
For each target:
- If it belongs to a new cluster, save the previous cluster, then compute the 
    possible existing decompositions and move on to the next target.
- If it belongs to the same cluster:
    - If the target falls into a category that has already been analyzed 
        (e.g., X0, DNA1, DNA2, RNA1, RNA2, where the number indicates the number 
        of bulges present), use the decomposition data already computed to avoid 
        recalculating it.
    - If the category has not been analyzed yet, perform a full decomposition 
        (similar to top1) and store the data for subsequent targets.
TODO: Parallelize the analysis.
NOTE: Verify if it is compatible with PAM at the beginning.

"""

# argv1: BED file containing annotations.
# argv2: .cluster.txt file, sorted by chromosome. As of 03/03, it was the top1 file sorted by chromosome.
# argv3: Output file name.
# argv4: Directory containing dictionaries.
# argv5: PAM file.
# argv6: Maximum allowed mismatches.
# argv7: Genome reference directory (e.g., ../../Genomes/hg38_ref).
# argv8: Guide file.
# argv9: Maximum allowed DNA bulges.
# argv10: Maximum allowed RNA bulges.
# argv11 (optional): the search's cap on total edits, mismatches + bulges (complete-search
#     --max-total-edits; -1 or absent: none). A target resolved per carrier haplotype
#     can carry more mismatches than the engine counted on the IUPAC-coded genome, and
#     is dropped past this cap as past --mm.

import pickle
from typing import final  # to read CFD matrices
from supportFunctions.loadSample import associateSample
import re
import multiprocessing
import string
import numpy as np
import pandas as pd
import subprocess
import concurrent.futures
from intervaltree import Interval, IntervalTree
import bisect
import os
import itertools
import time
import json
import gzip
import sys
import warnings

# Position-aware CFD bulge penalty (2.8.0); pass-through fallback on a very old deploy.
try:
    from bulge_cfd_penalty import cfd_score_bulge_aware as _cfd_bulge_aware
except Exception:  # pragma: no cover - module always shipped alongside
    _cfd_bulge_aware = None
# Clean-room CFD base scorer + its JSON tables (own implementation of Doench et al.,
# Nat. Biotechnol. 2016), replacing the inherited calc_cfd/revcom/get_mm_pam_scores.
from cfd_score import cfd as _cfd_clean, cfd_tables as _cfd_tables

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)
SIZE_DOENCH = 10000
N_THR = 3


def doenchForIupac(sequence_doench, guide_seq, genome_type):
    pos_iupac = []
    var = []
    for pos, c in enumerate(sequence_doench):
        if c in iupac_code:
            pos_iupac.append(pos)
            var.append(iupac_code[c])

    if var:
        for i in itertools.product(*var):
            t = list(sequence_doench)
            for p, el in enumerate(pos_iupac):
                t[el] = i[p]
            targets_for_doench[guide_seq][genome_type].append("".join(t))
    else:
        targets_for_doench[guide_seq][genome_type].append(sequence_doench)




class reversor:
    """
    Nel caso debba ordinare più campi però con reverse diversi, eg uno True e l'altro False, posso usare questa classe nella chiave per
    simulare il contrario del reverse applicato
    """

    def __init__(self, obj):
        self.obj = obj

    def __eq__(self, other):
        return other.obj == self.obj

    def __lt__(self, other):
        return other.obj < self.obj


def convert_fasta_to_dict(f):
    fasta = {}
    with open(f) as file_one:
        for line in file_one:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                active_sequence_name = line[1:]
                if active_sequence_name not in fasta:
                    fasta[active_sequence_name] = ""
                continue
            sequence = line
            fasta[active_sequence_name] = sequence
    return fasta


def alignRefFromVar(line):  # chr_fake, start_pos, len_guide, bulge):
    t = line.copy()
    # chr_fake = t[10].split('_')
    len_guide = len(t[2])
    start_pos = int(t[4])
    true_chr = t[3]
    good_chr_fake = true_chr + ":" + str(start_pos) + "-" + str(start_pos + len_guide)
    # dict_ref_seq[good_chr_fake].upper() #file_fasta.readline().strip().upper()
    sequence = genomeStr[start_pos : start_pos + len_guide].upper()
    # print('ref sequence', sequence, len(sequence), len_guide)

    if t[6] == "-":  # if negative strand reverse ref sequence
        # target = t[2][::-1]
        sequence = reverse_complement_table(sequence)
        target = t[2]
    else:
        target = t[2]

    if t[0] == "RNA":  # if RNA bulges in sequences find gaps and update ref sequence
        tmp_gap_position = [g.start() for g in re.finditer("-", target)]
        sequence = list(sequence)
        # cut sequence to extract only the part interested by counting the RNA bulges
        sequence = sequence[len(tmp_gap_position) :]
        for tmp_g_p in tmp_gap_position:
            # sequence[tmp_g_p] = '-'
            sequence.insert(tmp_g_p, "-")
        # if t[6] == '-':
        #     sequence = reverse_complement_table(''.join(sequence))
        # else:

        # sequence = sequence[len(tmp_gap_position):]
        sequence = "".join(sequence)
        # if t[6] == '+':
        #     sequence = sequence[0:len_guide]
        # else:
        #     sequence = reverse_complement_table(''.join(sequence))
        #     sequence = sequence[len(tmp_gap_position):]
    # elif t[6] == '-':
    #     sequence = reverse_complement_table(sequence)
    guide_no_pam = t[1][pos_beg:pos_end]
    list_t = list(sequence)

    # print('ref seq with bulges', sequence)
    # align ref sequence with var
    for position_t, char_t in enumerate(sequence[pos_beg:pos_end]):
        if char_t.upper() != guide_no_pam[position_t]:
            if guide_no_pam[position_t] != "-":
                list_t[sum_for_mms + position_t] = char_t.lower()
    sequence = "".join(list_t)

    # print('ref seq realigned', sequence)
    return sequence


# print('ESECUZIONE DI ANNOTATION E CALC SAMPLE INSIEME')
# print('TEST PER ANNOTAZIONE COMPLETA: I TARGET SENZA ANNOTAZIONE SONO SALVATI COME \"n\"')
# print('SE UN  TARGET HA 1+ ANNOTAZIONI, LE SALVA IN SINGOLA UNICA RIGA')
# print('RIMOZIONE DEI TARGET CHE NON HANNO SAMPLES')
# print('CALCOLO SCORES')
# print('SOSTITUZIONE IUPAC DI TUTTI I TARGET CON CHAR DEL TOP1SCOMPOSTO')
# print("READING INPUT FILES")
# Dictionaries for annotating samples

# Dict for populations
# pop_file = pd.read_excel(os.path.dirname(os.path.realpath(__file__)) + '/20130606_sample_info.xlsx')
# all_samples = pop_file.Sample.to_list()
# all_pop = pop_file.Population.to_list()
# dict_sample_to_pop = dict()
# for  pos, i in enumerate(all_samples):
#     try:
#         dict_sample_to_pop[i] = all_pop[pos]        #{'S1':'POP1', 'S2':'POP1', ...}
#     except:
#         dict_sample_to_pop[i] = all_pop[pos]

# #Dict for superpopulation
# dict_pop_to_sup = {'CHB':'EAS', 'JPT':'EAS', 'CHS':'EAS', 'CDX':'EAS', 'KHV':'EAS',
#                     'CEU':'EUR', 'TSI':'EUR', 'FIN':'EUR', 'GBR':'EUR', 'IBS':'EUR',
#                     'YRI':'AFR', 'LWK':'AFR', 'GWD':'AFR', 'MSL':'AFR', 'ESN':'AFR', 'ASW':'AFR', 'ACB':'AFR',
#                     'MXL':'AMR', 'PUR':'AMR', 'CLM':'AMR', 'PEL':'AMR',
#                     'GIH':'SAS', 'PJL':'SAS', 'BEB':'SAS', 'STU':'SAS', 'ITU':'SAS'
# }
# superpopulation = ['EAS', 'EUR', 'AFR', 'AMR','SAS']

# dict_sample_to_pop, dict_pop_to_sup, dict_superpop_to_pop, dict_pop_to_sample, all_samples, all_pop, superpopulation, gender_sample = associateSample.loadSampleAssociation(sys.argv[11])


# READ INPUT FILES
annotationFile = sys.argv[1]  # file with annotation
resultsFile = sys.argv[2]  # file with results from search
outputFile = sys.argv[3]  # file with annotated results

# Get pam and guide length for new count mismatch samples
pam_at_beginning = False
with open(sys.argv[5]) as pam:
    line = pam.read().strip()
    pam = line.split(" ")[0]
    len_pam = int(line.split(" ")[1])
    guide_len = len(pam) - len_pam
    pos_beg = 0
    pos_end = None
    pam_begin = 0
    pam_end = len_pam * (-1)
    if len_pam < 0:
        guide_len = len(pam) + len_pam
        pam = pam[: (len_pam * (-1))]
        len_pam = len_pam * (-1)
        pos_beg = len_pam
        pos_end = None
        pam_begin = 0
        pam_end = len_pam
        pam_at_beginning = True
    else:
        pam = pam[(len_pam * (-1)) :]
        pos_beg = 0
        pos_end = len_pam * (-1)
        pam_begin = len_pam * (-1)
        pam_end = None

# Second (ML) score column = CRISPR-Bulge (the legacy RandomForest scorer was retired); env var kept only so the
# report/web can label the column.
_SCORER_SELECT = os.environ.get("CRISPRME_SCORER_SELECT", "crispr-bulge").lower()

do_scores = True
if guide_len != 20 or len_pam != 3 or pam_at_beginning:
    # sys.stderr.write('CFD SCORE IS NOT CALCULATED WITH GUIDES LENGTH != 20 OR PAM LENGTH !=3 OR UPSTREAM PAM')
    do_scores = False
    with open(outputFile + ".acfd.txt", "w+") as result:
        result.write("NO SCORES")

iupac_code_set = {
    "R": {"A", "G"},
    "Y": {"C", "T"},
    "S": {"G", "C"},
    "W": {"A", "T"},
    "K": {"G", "T"},
    "M": {"A", "C"},
    "B": {"C", "G", "T"},
    "D": {"A", "G", "T"},
    "H": {"A", "C", "T"},
    "V": {"A", "C", "G"},
    "r": {"A", "G"},
    "y": {"C", "T"},
    "s": {"G", "C"},
    "w": {"A", "T"},
    "k": {"G", "T"},
    "m": {"A", "C"},
    "b": {"C", "G", "T"},
    "d": {"A", "G", "T"},
    "h": {"A", "C", "T"},
    "v": {"A", "C", "G"},
    "A": {"A"},
    "T": {"T"},
    "C": {"C"},
    "G": {"G"},
    "a": {"a"},
    "t": {"t"},
    "c": {"c"},
    "g": {"g"},
    "N": {"A", "T", "G", "C"},
}


# OPEN INPUT FILES AND PREPARE OUTPUT FILE
inResult = open(resultsFile, "r")  # resultfile open
# inAnnotationFile = open(annotationFile, "r")  # file with annotations open
# outFileSampleAll = open(outputFile + '.samples.all.annotation.txt', 'w')  # outfile open (file with IUPAC targets and associated samples and annotation)

count_removed_target = 0
process = subprocess.Popen(
    ["wc", "-l", resultsFile], stdout=subprocess.PIPE, stderr=subprocess.PIPE
)
out, err = process.communicate()
total_line = int(out.decode("UTF-8").split(" ")[0])
if total_line < 1:
    print("WARNING! Input file has no targets")
    sys.exit()
if total_line < 10:
    mod_tot_line = 1
else:
    mod_tot_line = int(total_line / 10)
# VARIABLE INIT
guideDict = {}
totalDict = {}

start_time = time.time()

# print("EXECUTING PRELIMINARY OPERATIONS")

# annotationsTree = IntervalTree()
# annotationsSet = set()
# #guidesSet = set()       #NOTE/BUG if guide finds 0 targets, it will not be annotated

# for line in inAnnotationFile:
#     x = line.split('\t')
#     x[3] = str(x[3]).rstrip("\n")
#     annotationsTree[int(x[1]):int(x[2])] = str(x[0])+'\t'+str(x[3])
#     annotationsSet.add(str(x[3]))

# totalDict['targets'] = [0]*10
# for item in annotationsSet:
#     totalDict[item] = [0]*10

# print("PRELIMINARY OPERATIONS COMPLETED IN: %s seconds" % (time.time() - start_time))

start_time = time.time()

# print("EXECUTING ANNOTATION")

with open(resultsFile, "r") as resFile:
    header_len = len(resFile.readline().strip().split("\t"))

# if header_len == 15:    #'Both' case : comparison variant/ref is active
# header = '#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation Type\tReal Guide'
header = "#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\t#Seq_in_cluster\tReference"


mm_pos = 7  # position of mismatch column
bulge_pos = 8
max_dna_bulges = int(sys.argv[9])
max_rna_bulges = int(sys.argv[10])
max_total_edits = int(sys.argv[11]) if len(sys.argv) > 11 else -1
max_bulges = max_dna_bulges
if max_rna_bulges > max_bulges:
    max_bulges = max_rna_bulges

blank_add_begin = " "  # Needed for replacing IUPAC in cluster targets
blank_add_end = ""
pam_multiplier = 1
pam_multiplier_negative = 0
start_sample_for_cluster = 0
cluster_step = 1  # If PAM end, go left to right
# when updatig lowercase for nem_mm, this value represents the offset for the pam position (mainly needed only if pam at beginning)
sum_for_mms = 0
# Values to check new iupac when working on cluster targets
end_sample_for_cluster = max_dna_bulges + max_rna_bulges
if pam_at_beginning:
    blank_add_begin = ""
    blank_add_end = " "
    pam_multiplier = 0  # Since ' ' are at end, and '-' to reinsert are before the ' ', need to put max_dna_bulges and rna_bulges of target to 0
    pam_multiplier_negative = 1
    # For PAM at beginning, start from last nucleotide and go to left
    end_sample_for_cluster = len_pam + guide_len - max_rna_bulges
    start_sample_for_cluster = len_pam + guide_len + max_dna_bulges
    cluster_step = -1  # If PAM beginning, go right to left
    sum_for_mms = len_pam
# outFileSampleAll.write(header + '\n')
summary_samples = True

header_list = header.strip().split("\t")
# Variables for summary samples code
"""
{
    GUIDE1 -> {
        SAMPLE/POP/SUPERPOP1 ->{
            targets -> [0 0 0 0 0 0 0 0 0],
            ann1 -> [0 0 0 0 0 0 0 0 0],
            ann2 -> [0 0 0 0 0 0 0 0 0],
        },
        SAMPLE/POP/SUPERPOP2 ->{
            targets -> [0 0 0 0 0 0 0 0 0],
            ann1 -> [0 0 0 0 0 0 0 0 0],
            ann2 -> [0 0 0 0 0 0 0 0 0],
        }
    }
    GUIDE2 -> {
        SAMPLE/POP/SUPERPOP1 ->{
            targets -> [0 0 0 0 0 0 0 0 0],
            ann1 -> [0 0 0 0 0 0 0 0 0],
            ann2 -> [0 0 0 0 0 0 0 0 0],
        },
        SAMPLE/POP/SUPERPOP2 ->{
            targets -> [0 0 0 0 0 0 0 0 0],
            ann1 -> [0 0 0 0 0 0 0 0 0],
            ann2 -> [0 0 0 0 0 0 0 0 0],
        }
    }
}

Per pop e superpop, se ho due sample stessa famiglia stesso target, conto solo una volta (visited_pop and visited_superpop array)
"""
count_sample = (
    dict()
)  # NOTE cout_sample -> GUIDE -> SAMPLE -> has targets + ann1 + ann2 ... + refposition
# refposition is a key unrelated to the other keys (targets ann1 ann2 ...) and it's used to classify the sample (0 0+ 1 1+).
# it's put in here just to avoid to duplicate the entire guide -> sample ->     structure
# refposition -> [class , number of specific VAR on target to add/remove]  #Save class (0 at start) and number of ontarget var specific
# for that sample.
# Count number of REF target and REF part of semicommon
ontarget_reference_count = dict()
count_pop = dict()
# NOTE added key 'distributions' for population distribution images
count_superpop = dict()
# count_superpop-> GUIDE -> SUPERPOP -> targets ann1 ann2 ... distributions
# distributions is an array of len mms+bulge, each position contains an array [0,0,0] of len bulge+1 (indicating no bulge, 1 bulge, 2bulge ...)

# Create -Summary_total for a file ref.Annotation.summary.txt from the y and n values of Var_uniq column
summary_barplot_from_total = False
if "Var_uniq" in header:
    summary_barplot_from_total = True
    vu_pos = header_list.index("Var_uniq")
# count_unique = dict()
# count_unique['targets'] = [0]*10
# count_unique_for_guide = dict()
# for item in annotationsSet:
#     count_unique[item] = [0]*10

# Variables for samples calculation
total_error = 0


current_chr = "none"
chr_name = "none"


def rev_comp(a):
    if a == "A" or a == "a":
        return "T"
    if a == "T" or a == "t":
        return "A"
    if a == "C" or a == "c":
        return "G"
    return "C"


def _cfd_from_alignment(aligned_crrna, aligned_proto):
    """CFD for one in-place display alignment (PAM included, gaps '-' in place) with the
    position-aware bulge penalty (2.8.0); bulge-agnostic clean-room fallback if the helper
    is unavailable. Mirrors new_simple_analysis._cfd_from_alignment."""
    if _cfd_bulge_aware is not None:
        return _cfd_bulge_aware(
            aligned_crrna, aligned_proto, mm_scores, pam_scores, do_scores=do_scores
        )
    up = aligned_proto.upper()
    return _cfd_clean(aligned_crrna, up[:-3], up[-2:], (mm_scores, pam_scores))


def preprocess_CFD_score(target):
    # preprocess target then calculate CFD score (unified bulge-aware path -- see
    # new_simple_analysis.preprocess_CFD_score).
    if do_scores:
        cfd_score = _cfd_from_alignment(target[1], target[2])
        target.append("{:.3f}".format(cfd_score))
        if target[-3] == 55:
            target[-3] = "{:.3f}".format(cfd_score)
        if target[-3] == 33:
            cfd_ref_score = _cfd_from_alignment(target[1], target[-4])
            target[-3] = "{:.3f}".format(cfd_ref_score)
    else:
        cfd_score = -1
        target.append("{:.3f}".format(cfd_score))
        target[-3] = "{:.3f}".format(cfd_score)
    # print(target)
    # print('cfd', cfd_score)
    return target


def preprocess_CRISPR_BULGE_score(cluster_targets):
    """CRISPR-Bulge scoring (analogue of the retired legacy-scorer preprocessing path, indel path).

    Identical output contract + alt/ref two-pass + null semantics as the CRISPR-Bulge path,
    but scores the aligned (sgRNA, off-target) pair only (NO 29-nt genomic window) via
    the scorer-runner (dedicated cbulge env; graceful -1.0 if absent). Unscoreable rows
    (N in aligned DNA, or length mismatch) get an equal-length dummy then are nulled by
    index, so CRISPR-Bulge's per-row len(sg)==len(off) requirement always holds.
    """
    import scorer_runner  # lazy: only when CRISPR-Bulge is selected

    cluster_scored = list()
    index_to_null = list()

    if not do_scores:
        for target in cluster_targets:
            t = target.copy()
            t[-2] = "{:.6g}".format(-1)
            t.append("{:.6g}".format(-1))
            cluster_scored.append(t)
        return cluster_scored

    # seq_len=24 model: faithfully scores a protospacer+PAM alignment with AT MOST ONE 1-bp
    # bulge. Score ONLY what the model represents and NULL (-1) everything else so a value is
    # always trustworthy: N in the aligned off, sg/off length mismatch, >1 gap column (>=2
    # bulges OR a single >=2-bp bulge), or over-length (>24). Nulled rows fall back to CFD
    # (PRIMARY) + mismatches/bulges. We deliberately do NOT collapse/reduce multi-bulge
    # alignments to one bulge -- an adversarial review showed the naive reduction corrupts the
    # alignment and its "conservative" direction is unprovable for a non-monotonic GRU (silent
    # under-flag risk). Multi-bulge OTs are the rare, heavily-disrupted tail; CFD + edit distance
    # carry them. A scoreable equal-length dummy is substituted for nulled rows.
    _MAXLEN = 24
    _SG_DUMMY = "A" * 20 + "NGG"
    _OFF_DUMMY = "A" * 20 + "AGG"

    def _pair(target, off_field, index):
        sg = str(target[1])[: len(str(target[1])) - 3] + "NGG"
        off = str(off_field)
        total_gaps = sg.count("-") + off.count("-")
        if (("N" in off) or ("n" in off) or (len(off) != len(sg))
                or total_gaps > 1 or len(sg) > _MAXLEN):
            index_to_null.append(index)
            return _SG_DUMMY, _OFF_DUMMY
        return sg, off

    pairs_alt = [_pair(t, t[2], i) for i, t in enumerate(cluster_targets)]
    scores_alt = scorer_runner.CRISPR_BULGE_predict_list(
        [p[0] for p in pairs_alt], [p[1] for p in pairs_alt])

    pairs_ref = []
    for i, target in enumerate(cluster_targets):
        ref = str(target[-3]) if "n" not in str(target[-3]) else str(target[2])
        pairs_ref.append(_pair(target, ref, i))
    scores_ref = scorer_runner.CRISPR_BULGE_predict_list(
        [p[0] for p in pairs_ref], [p[1] for p in pairs_ref])

    # CRISPR-Bulge scores at FULL precision ("{:.6g}") -- the 95%-recall sensitive
    # default (1.9e-4) is below the old 3-decimal floor (see cb_operating_point).
    for index, target in enumerate(cluster_targets):
        t = target.copy()
        if index in index_to_null:
            s = -1
            t[-2] = "{:.6g}".format(s)
            t.append("{:.6g}".format(s))
        else:
            if t[-2] == 55:
                t[-2] = "{:.6g}".format(scores_alt[index])
                t.append("{:.6g}".format(scores_alt[index]))
            if t[-2] == 33:
                t[-2] = "{:.6g}".format(scores_ref[index])
                t.append("{:.6g}".format(scores_alt[index]))
        cluster_scored.append(t)

    return cluster_scored


def calculate_scores(cluster_to_save):
    # function to calculate score for each input target
    # input is target line splitted in list format
    # list of functions to calculate specific score (to add a score, simply add your function to this call and update the clusters list in return)
    cluster_with_CFD_score = list()
    cluster_with_CRISPR_BULGE_score = list()

    # time_start = time.time()
    # print('CALCULATING CFD SCORE')

    for target in cluster_to_save:  # calculate CFD score for each target
        target_CFD = target.copy()
        cluster_with_CFD_score.append(preprocess_CFD_score(target_CFD))

    # second (ML) score column = CRISPR-Bulge (the legacy RandomForest scorer was retired); same .bestCRISPR_BULGE.txt /
    # CRISPR_BULGE_score column names kept as stable identifiers (report/web label it CRISPR-Bulge).
    cluster_with_CRISPR_BULGE_score = preprocess_CRISPR_BULGE_score(cluster_to_save)

    # REMOVED TO CHECK IF FILE IS RETURN WITH IDENTICAL ROWS COUNT

    # analyze CFD scored targets, returning for each guide,chr,cluster_pos the highest scoring target (or multiple in case of equal)
    # df_CFD = pd.DataFrame(cluster_with_CFD_score, columns=['Bulge_type', 'crRNA', 'DNA', 'Chromosome',
    #    'Position', 'Cluster_Position', 'Direction', 'Mismatches',
    #    'Bulge_Size', 'Total', 'PAM_gen', 'Var_uniq', 'Samples', 'Annotation_Type',
    #    'Real_Guide', 'rsID', 'AF', 'SNP', 'Reference_target', 'CFD',
    #    'Seq_in_cluster', 'CFD_ref'])
    # group by over real_guide,chr,cluster_pos to avoid mixing targets
    # select lowest count of mm+bul
    # idx_fewest_mm_bul = df_CFD.groupby(['Real_Guide', 'Chromosome', 'Cluster_Position', 'SNP', 'Samples'])[
    #     'Total'].transform(min) == df_CFD['Total']
    # df_CFD_fewest = df_CFD[idx_fewest_mm_bul]
    # df_CFD_fewest.drop_duplicates(
    #     ['Real_Guide', 'Chromosome', 'Cluster_Position', 'SNP', 'Samples', 'Mismatches', 'Bulge_Size'], inplace=True)
    # # select highest score
    # idx_max_score = df_CFD.groupby(['Real_Guide', 'Chromosome', 'Cluster_Position', 'SNP', 'Samples'])[
    #     'CFD'].transform(max) == df_CFD['CFD']
    # df_CFD_best_score = df_CFD[idx_max_score]
    # # remove duplicate rows (possible due to haplotypes and variants)
    # df_CFD_best_score.drop_duplicates(['Real_Guide', 'Chromosome',
    #                                    'Cluster_Position', 'SNP', 'Samples', 'CFD'], inplace=True)
    # frames = [df_CFD_fewest, df_CFD_best_score]
    # df_CFD = pd.concat(frames)
    # cluster_with_CFD_score = df_CFD.values.tolist()

    # df_CRISPR_BULGE = pd.DataFrame(cluster_with_CRISPR_BULGE_score, columns=['Bulge_type', 'crRNA', 'DNA', 'Chromosome',
    #                                                              'Position', 'Cluster_Position', 'Direction', 'Mismatches',
    #                                                              'Bulge_Size', 'Total', 'PAM_gen', 'Var_uniq', 'Samples', 'Annotation_Type',
    #                                                              'Real_Guide', 'rsID', 'AF', 'SNP', 'Reference_target', 'CFD',
    #                                                              'Seq_in_cluster', 'CFD_ref'])
    # # select lowest count of mm+bul
    # idx_fewest_mm_bul = df_CRISPR_BULGE.groupby(['Real_Guide', 'Chromosome', 'Cluster_Position', 'SNP', 'Samples'])[
    #     'Total'].transform(min) == df_CRISPR_BULGE['Total']
    # df_CRISPR_BULGE_fewest = df_CRISPR_BULGE[idx_fewest_mm_bul]
    # df_CRISPR_BULGE_fewest.drop_duplicates(
    #     ['Real_Guide', 'Chromosome', 'Cluster_Position', 'SNP', 'Samples', 'Mismatches', 'Bulge_Size'], inplace=True)
    # # select highest score
    # idx_max_score = df_CRISPR_BULGE.groupby(['Real_Guide', 'Chromosome', 'Cluster_Position', 'SNP', 'Samples'])[
    #     'CFD'].transform(max) == df_CRISPR_BULGE['CFD']
    # df_CRISPR_BULGE_best_score = df_CRISPR_BULGE[idx_max_score]
    # # remove duplicate rows (possible due to haplotypes and variants)
    # df_CRISPR_BULGE_best_score.drop_duplicates(['Real_Guide', 'Chromosome',
    #                                       'Cluster_Position', 'SNP', 'Samples', 'CFD'], inplace=True)
    # frames = [df_CRISPR_BULGE_fewest, df_CRISPR_BULGE_best_score]
    # df_CRISPR_BULGE = pd.concat(frames)
    # cluster_with_CRISPR_BULGE_score = df_CRISPR_BULGE.values.tolist()

    return [cluster_with_CFD_score, cluster_with_CRISPR_BULGE_score]


iupac_code = {
    "R": ("A", "G"),
    "Y": ("C", "T"),
    "S": ("G", "C"),
    "W": ("A", "T"),
    "K": ("G", "T"),
    "M": ("A", "C"),
    "B": ("C", "G", "T"),
    "D": ("A", "G", "T"),
    "H": ("A", "C", "T"),
    "V": ("A", "C", "G"),
    "r": ("A", "G"),
    "y": ("C", "T"),
    "s": ("G", "C"),
    "w": ("A", "T"),
    "k": ("G", "T"),
    "m": ("A", "C"),
    "b": ("C", "G", "T"),
    "d": ("A", "G", "T"),
    "h": ("A", "C", "T"),
    "v": ("A", "C", "G"),
    "N": ("A", "T", "C", "G"),
}


# The IUPAC codes of an indel target's DNA column (the SNP overlay of the fake-indel
# genome, CRISPRME_INDEL_SNP) are resolved per carrier haplotype by _overlaid_rows,
# below: never column by column, which would pair one sample's indel with SNP alleles
# only other samples carry.

# Whole-genome IUPAC -> first-listed-allele map (fast str.translate), applied to
# genomeStr once at load. Mirrors iupac_code's first entry: R->A Y->C S->G W->A K->G
# M->A B->C D->A H->A V->A. 'N' is absent, so it survives for the existing N-guard.
# genomeStr is the PLAIN reference chromosome (post_analisi_indel.sh passes
# <ref>/<chrom>.fa), which holds no IUPAC code, so this is a no-op kept as a guard.
_IUPAC_FIRST_ALLELE_TAB = str.maketrans(
    "RYSWKMBDHVryswkmbdhv", "ACGAGACAAAacgagacaaa"
)


# For scoring of CFD And Doench
tab = str.maketrans("ACTGRYSWMKHDBVactgryswmkhdbv", "TGACYRSWKMDHVBtgacyrswkmdhvb")


def reverse_complement_table(seq):
    return seq.translate(tab)[::-1]


mm_scores, pam_scores = _cfd_tables()  # clean-room JSON tables
guides_dict = dict()  # For CFD score
guides_dict_doench = dict()
targets_for_doench = dict()

N_THR = multiprocessing.cpu_count() // 2
refgenomedir = sys.argv[7]

max_doench = 0
sum_cfd = 0
cfd_scores = []

# open fasta to read ref chromosome
inFasta = open(sys.argv[7], "r")
current_chr = inFasta.readline().strip().replace(">", "")  # lettura fasta del chr
genomeStr = inFasta.readlines()  # lettura fasta del chr
genomeStr = "".join(genomeStr).upper()
# string of the whole chromosome on single line
genomeStr = genomeStr.replace("\n", "")
# Keep every genome-derived output column (ref alignment, CRISPR-Bulge context,
# motif/radar counts) pure A/C/G/T should an IUPAC code ever reach genomeStr; a no-op
# on the plain reference it is read from. 'N' is preserved for the existing N-guards.
genomeStr = genomeStr.translate(_IUPAC_FIRST_ALLELE_TAB)


start_time_total = time.time()
lines_processed = 0
allowed_mms = int(sys.argv[6])
current_guide_chr_pos_direction = "no"

# save best files
cfd_best = open(outputFile + ".bestCFD_INDEL.txt", "w")
cfd_best.write(header + "\tCFD\n")  # Write header

crispr_bulge_best = open(outputFile + ".bestCRISPR_BULGE_INDEL.txt", "w")
crispr_bulge_best.write(header + "\tCFD\n")

mmblg_best = open(outputFile + ".bestmmblg_INDEL.txt", "w")
mmblg_best.write(header + "\tCFD\n")  # Write header


INDELS_tree = IntervalTree()
# prefer a gzip-compressed indel log (log<chrom>.txt.gz) when present so batteries
# installs keep the logs compressed on disk and read them on the fly; fall back to plain
_logbase = os.path.realpath(sys.argv[4]) + "/log" + current_chr + ".txt"
_logpath = _logbase + ".gz" if os.path.exists(_logbase + ".gz") else _logbase
with (gzip.open(_logpath, "rt") if _logpath.endswith(".gz") else open(_logpath, "r")) as log:
    # print('indel processing:', current_chr)
    log.readline()
    for entry in log:
        splitted = entry.strip().split("\t")
        fake_pos = splitted[5].strip().split(",")  # inizio,fine del fakechr
        # start       end fakechr          #pos_vera,samples,rsID,MAF,INDELinfo,start_fake,ref_seq
        INDELS_tree[int(fake_pos[0]) : int(fake_pos[1])] = [
            splitted[0],
            splitted[1],
            splitted[2],
            splitted[3],
            splitted[4],
            int(fake_pos[0]),
            splitted[6],
        ]
# datastore = datastore.to_dict(orient='index')

# [indel-snp] SNP+indel co-occurrence (default-on; CRISPRME_INDEL_SNP=0 opts out).
# The SNP registry/genotype tiers are SIBLINGS of the log_indels dir (sys.argv[4])
# under Dictionaries/, so no new CLI/shell args are needed. Each one is optional:
#   * the phased indel genotype store and the Tier-1 SNP genotypes (or, on a classic
#     dictionary install, the per-sample SNP dictionary) resolve an overlaid target
#     per carrier haplotype (_overlaid_rows); without them each carrier gets the
#     worst case its own alleles allow, or the code's alleles allow, PUTATIVE;
#   * the Tier-0 registry names the SNPs (rsID, AF) and gives the panel size for the
#     joint allele frequency of the SNP+indel co-occurrence companion TSV, which is
#     written only with it.
# A target with no IUPAC code (every target of a classic, non-overlaid fake-indel
# genome) never reaches any of this: its row is the one the engine reported.
import indel_snp_cis as _isc

_indel_snp = os.environ.get("CRISPRME_INDEL_SNP", "1") in ("1", "true", "True", "yes")
_reg = _gt = _indelgt = _cooc_out = _dict_snps = None
_droot = os.path.dirname(os.path.realpath(sys.argv[4]))
_vcf = os.path.basename(os.path.realpath(sys.argv[4]))[len("log_indels_"):]
if _indel_snp:
    try:
        import tier0_registry as _t0
        _rb = os.path.join(_droot, "registry_" + _vcf, "reg_" + current_chr + ".bin")
        _ri = os.path.join(_droot, "registry_" + _vcf, "reg_" + current_chr + ".idx")
        if os.path.isfile(_rb) and os.path.isfile(_ri):
            _reg = _t0.RegistryReader(_rb, _ri)
    except Exception:  # noqa: BLE001 - registry optional
        _reg = None
    try:
        import tier1_genotypes as _t1
        _gb = os.path.join(_droot, "genotypes_" + _vcf, "gt_" + current_chr + ".bin")
        _gi = os.path.join(_droot, "genotypes_" + _vcf, "gt_" + current_chr + ".idx")
        if os.path.isfile(_gb) and os.path.isfile(_gi):
            _gt = _t1.GenotypeReader(_gb, _gi)
    except Exception:  # noqa: BLE001 - genotype tier optional
        _gt = None
    try:
        import build_indel_genotypes as _big
        _ig = os.path.join(_droot, "indel_genotypes_" + _vcf,
                           "gt_indel_" + current_chr + ".tsv.gz")
        if os.path.isfile(_ig):
            _indelgt = _big.IndelGenotypeReader(_ig)
    except Exception:  # noqa: BLE001 - phased indel GT optional (-> PUTATIVE)
        _indelgt = None
    if _reg is not None:
        try:
            _cooc_out = open(outputFile + ".indel_snp_cooc.tsv", "w")
            _cooc_out.write(
                "chrom\tindel_pos\tindel_ref\tindel_alt\tofftarget_start\tstrand\t"
                "snp_dictpos\tsnp_rsid\tphase\tjoint_af\tn_cis\tcis_samples\n"
            )
        except Exception as _isp_err:  # noqa: BLE001 - the companion is optional
            print(f"WARNING [indel-snp]: SNP+indel co-occurrence companion disabled "
                  f"({_isp_err})", flush=True)
            _cooc_out = None

# [mega] sites-only INDEL per-dataset AF sidecar (INDEPENDENT of indel+SNP cooc):
# a merged all-source "mega" index has no genotypes, so indel off-targets cannot get
# their per-dataset AF from the fake-indel MAF (which is empty). When present, the
# build_mega_indel_af.py sidecar (indel_af_<vcf>/indel_af_<chrom>.tsv.gz, a sibling
# of registry_<vcf>/) supplies per-dataset AF + AF_max, emitted to a companion TSV.
# Fully guarded: absent sidecar or any error -> no-op; classic indel row untouched.
_indel_af = _indel_af_out = None
_indel_af_errs = [0]  # emission-error counter (warn once, never break classic path)
try:
    import build_mega_indel_af as _bia
    _root2 = os.path.dirname(os.path.realpath(sys.argv[4]))
    _vcf2 = os.path.basename(os.path.realpath(sys.argv[4]))[len("log_indels_"):]
    _iaf = os.path.join(_root2, "indel_af_" + _vcf2,
                        "indel_af_" + current_chr + ".tsv.gz")
    if os.path.isfile(_iaf):
        _indel_af = _bia.IndelAfReader(_iaf)
        _indel_af_out = open(outputFile + ".indel_af.tsv", "w")
        _indel_af_out.write(
            "chrom\tofftarget_start\tindel_pos\tindel_ref\tindel_alt\t"
            + "\t".join("AF_" + _l for _l in _indel_af._labels) + "\tAF_max\n")
except Exception as _iaf_err:  # noqa: BLE001 - never break the classic path
    print(f"WARNING [mega-indel-af]: sidecar disabled ({_iaf_err})", flush=True)
    _indel_af = _indel_af_out = None


def _indel_window(indel_data):
    """(pos, REF, ALT, real 0-based start of the fake window) of a log entry."""
    _ip = indel_data[4].split("_")  # chrN_pos_REF_ALT
    return (int(_ip[-3]), _ip[-2], _ip[-1],
            int(indel_data[0].split("_")[1].split("-")[0]))


def _overlaid_positions(path):
    """The 1-based positions of every IUPAC code of this chromosome's overlaid indel
    targets: the per-sample dictionary entries the per-haplotype resolution reads."""
    out = set()
    with open(path) as fh:
        fh.readline()  # header
        for raw in fh:
            ln = raw.rstrip("\n").split("\t")
            if len(ln) < 7 or ln[3] != "fake" + current_chr:
                continue
            amb = _isc.ambiguous_columns(ln[2])
            if not amb:
                continue
            try:
                idata = sorted(INDELS_tree[int(ln[4])])[0].data
            except Exception:  # noqa: BLE001 - no indel here: skipped by the loop too
                continue
            _ipos, _iref, _ialt, _istart = _indel_window(idata)
            real_of = _isc.build_offset_to_real(
                ln[2], int(ln[4]), ln[6], int(idata[5]), _istart, _iref, _ialt)
            for j in amb:
                if real_of[j] is not None:
                    out.add(real_of[j] + 1)
    return out


# A classic dictionary install (an index built without --samplesID) has no Tier-1
# store: the per-sample genotypes live in the SNP dictionary, read here for the
# positions the overlaid targets need.
if _indel_snp and _gt is None:
    _dict_path = os.path.join(_droot, "dictionaries_" + _vcf,
                              "my_dict_" + current_chr + ".json")
    if os.path.exists(_dict_path) or os.path.exists(_dict_path + ".gz"):
        try:
            _dict_snps = _isc.load_dict_snps(_dict_path, current_chr,
                                             _overlaid_positions(resultsFile))
        except Exception as _dict_err:  # noqa: BLE001 - resolve without genotypes
            print(f"WARNING [indel-snp]: SNP dictionary {_dict_path} unreadable "
                  f"({_dict_err}); overlaid indel targets are resolved without "
                  "per-sample genotypes (PUTATIVE)", flush=True)
            _dict_snps = None

_snp_cache = {}    # pos1 -> {alt: {sample: genotype}} from the genotype source
_igt_cache = {}    # indel token -> its carriers (indel_snp_cis.IndelCarriers)
# A common variant has thousands of carriers, so each cache is bounded by the carriers
# it holds, not by its entries: once over the bound it is emptied and refilled.
_CACHE_CARRIERS = 1000000
_cache_held = {}


def _cache_put(cache, key, value, n):
    """Store ``value`` (``n`` carriers) in ``cache``, emptying it first when it would
    hold more than _CACHE_CARRIERS carriers."""
    held = _cache_held.get(id(cache), 0)
    if held + n > _CACHE_CARRIERS:
        cache.clear()
        held = 0
    cache[key] = value
    _cache_held[id(cache)] = held + n
    return value
# how the overlaid targets were resolved (printed at the end, on stdout)
_ovl_stats = {"targets": 0, "rows": 0, "dropped": 0, "putative": 0, "unmapped": 0,
              "second_indel": 0}


def _snp_record(pos1):
    """{forward alt: {sample: genotype}} the genotype source records at ``pos1``
    (Tier-1, else the per-sample dictionary; {} where it records none), or None when
    there is no genotype source."""
    rec = _snp_cache.get(pos1)
    if rec is not None:
        return rec
    if _gt is not None:
        rec = {}
        for alt in _gt.alts_at(pos1):
            if len(alt) == 1 and alt in "ACGT":
                rec[alt] = dict(tok.partition(":")[::2]
                                for tok in (_gt.carrier_tokens(pos1, alt) or []))
    elif _dict_snps is not None:
        rec = {alt: carriers for (_r, alt, _rs, _af, carriers) in _dict_snps.get(pos1, ())}
    else:
        return None
    return _cache_put(_snp_cache, pos1, rec, 1 + sum(len(c) for c in rec.values()))


def _snp_name(pos1, alt):
    """(rsID, MAF) of a SNP alt, formatted as the SNP path writes them: the
    registry's, else the per-sample dictionary's, else ('.', '0')."""
    if _reg is not None:
        rs = _reg.rsid(pos1, alt) or "."
        try:
            grp = _reg.lookup(pos1, alt)
            af = grp["global"].allele_freq() if grp and "global" in grp else 0.0
        except Exception:  # noqa: BLE001
            af = 0.0
        return rs, "%.6g" % af
    if _dict_snps is not None:
        for (_r, a, rs, af, _c) in _dict_snps.get(pos1, ()):
            if a == alt:
                try:
                    return rs, "%.6g" % float(af)
                except ValueError:
                    return rs, "0"
    return ".", "0"


def _snp_global(pos1, alt):
    """The registry's global group (Counts) of a SNP alt, or None."""
    try:
        grp = _reg.lookup(pos1, alt)
        return grp["global"] if grp and "global" in grp else None
    except Exception:  # noqa: BLE001
        return None


_indel_index = None  # [(pos, REF, ALT)] of the indel store, by position (built lazily)
_indel_starts = []
_indel_maxref = 0
_jgt_cache = {}      # (pos, REF, ALT) -> {sample: genotype} of an indel of the store


def _other_indels(lo, hi, ipos, iref, ialt):
    """Carriers ({sample: genotype}) of each OTHER indel of the phased store that
    changes the reference bases lo..hi (0-based) of a haplotype: a slot carrying one
    of them does not have the fake-indel window's sequence there. An indel whose
    REF overlaps this one's cannot share its haplotype and is skipped."""
    global _indel_index, _indel_starts, _indel_maxref
    if _indelgt is None:
        return ()
    if _indel_index is None:
        _indel_index = _indelgt.indels()
        _indel_starts = [p for p, _r, _a in _indel_index]
        _indel_maxref = max((len(r) for _p, r, _a in _indel_index), default=1)
    i0, i1 = ipos - 1, ipos - 1 + len(iref)  # this indel's REF, 0-based half-open
    out = []
    k = bisect.bisect_left(_indel_starts, lo - _indel_maxref)
    while k < len(_indel_index) and _indel_index[k][0] - 1 <= hi:
        p, r, a = _indel_index[k]
        k += 1
        if (p, r, a) == (ipos, iref, ialt) or (p - 1 < i1 and i0 < p - 1 + len(r)):
            continue
        if not _isc.changes_span(p, r, a, lo, hi):
            continue
        gts = _jgt_cache.get((p, r, a))
        if gts is None:
            gts = _indelgt.carriers_dict(p, r, a)
            _cache_put(_jgt_cache, (p, r, a), gts, 1 + len(gts))
        out.append(gts)
    return out


def _indel_carriers(indel_data, ipos, iref, ialt):
    """The indel's carriers (indel_snp_cis.IndelCarriers): its phased genotypes from
    the indel store, else the indel log's carriers with their phase unknown."""
    key = indel_data[4]
    carriers = _igt_cache.get(key)
    if carriers is None:
        igt = _indelgt.carriers_dict(ipos, iref, ialt) if _indelgt is not None else {}
        if igt:
            carriers = _isc.IndelCarriers(igt)
        else:
            carriers = _isc.IndelCarriers(
                {s: "0/1" for s in indel_data[1].split(",") if s}, phase_known=False)
        _cache_put(_igt_cache, key, carriers, 1 + len(carriers.parsed))
    return carriers


def _overlaid_rows(line, final_result, indel_data, amb):
    """Resolve an overlaid indel target per carrier haplotype (indel_snp_cis).

    ``amb`` holds the columns of the engine's DNA (line[2]) that carry an IUPAC code.
    Each indel carrier slot gets its own allele there (the SNP alt it carries, else
    the reference base), and one row is returned per resulting sequence that still
    has the PAM and fits the search budget, listing only the carriers that have it
    and naming the SNP alts it holds next to the indel (Variant info, rsID, MAF).
    A carrier the genotypes cannot place on a slot gets the strongest target its own
    alleles allow (PUTATIVE); with no genotypes at all, the indel log's carriers get
    the strongest target the codes allow (PUTATIVE). With CRISPRME_INDEL_SNP=0 every
    code reads as the reference base, as on a fake-indel genome without the overlay.
    The SNP+indel co-occurrence companion gets one line per SNP set a row requires.
    Returns [] when no carrier has a target here."""
    _ipos, _iref, _ialt, _istart = _indel_window(indel_data)
    strand = line[6]
    real_of = _isc.build_offset_to_real(
        line[2], int(line[4]), strand, int(indel_data[5]), _istart, _iref, _ialt)
    have_gt = _indel_snp and (_gt is not None or _dict_snps is not None)
    cols = []
    for j in amb:
        real = real_of[j]
        code = line[2][j]
        if real is None or not 0 <= real < len(genomeStr):
            # never expected (the overlay codes only reference positions): read the
            # code as its first allele, unresolved, and say so at the end
            _ovl_stats["unmapped"] += 1
            first = _isc._CODE[code.upper()][0]
            ref_fwd = _isc.complement(first) if strand == "-" else first
            cols.append(_isc.SnpColumn(j, -1, ref_fwd, None, code, strand, known=False))
            continue
        ref_fwd = genomeStr[real]
        if not _indel_snp:
            cols.append(_isc.SnpColumn(j, real, ref_fwd, None, ref_fwd if strand == "+"
                                       else _isc.complement(ref_fwd), strand,
                                       known=False))
            continue
        rec = _snp_record(real + 1) if have_gt else None
        col = _isc.SnpColumn(j, real, ref_fwd, rec, code, strand, known=rec is not None)
        if col.ref not in col.code:
            _ovl_stats["unmapped"] += 1  # the code does not hold the reference base
            col.known = False
        cols.append(col)
    carriers = None
    others = ()
    if have_gt:
        # the indel's phased carriers (the store), else the log's, phase unknown
        carriers = _indel_carriers(indel_data, _ipos, _iref, _ialt)
        if carriers.phase_known:
            _span = [r for r in real_of if r is not None]
            if _span:
                others = _other_indels(min(_span), max(_span), _ipos, _iref, _ialt)
    if do_scores:
        def score_fn(bases, _cols=cols):
            return _cfd_from_alignment(line[1], _isc.realize(line[2], line[1], _cols, bases))
    else:
        score_fn = None
    rows = _isc.resolve_target(
        line[1], line[2], cols, carriers,
        _isc.pam_columns(len(line[2]), pam, pam_at_beginning),
        range(len(line[2]))[pos_beg:pos_end], allowed_mms, max_total_edits,
        int(line[8]), score_fn=score_fn,
        max_dna_bulges=max_dna_bulges, max_rna_bulges=max_rna_bulges, others=others)
    _ovl_stats["targets"] += 1
    _ovl_stats["rows"] += len(rows)
    if not rows:
        _ovl_stats["dropped"] += 1
    cooc = {}  # (phase, SNP alts the row requires) -> [samples, slots]
    out = []
    for r in rows:
        if r["phase"] == _isc.PUTATIVE:
            _ovl_stats["putative"] += 1
        fr = final_result.copy()
        fr[2] = r["seq"]
        fr[7] = str(r["mm"])
        fr[9] = str(r["mm"] + int(line[8]))
        if r["samples"] is not None:
            fr[12] = ",".join(r["samples"])
        # the SNP alts this haplotype holds, forward strand, by position; the indel
        # token first (resultIntegrator strand-corrects each comma-separated token)
        alts = sorted(
            (c.pos1, c.ref_fwd, _isc.complement(b) if strand == "-" else b)
            for c, b in zip(cols, r["bases"]) if b != c.ref and c.pos1 > 0)
        if alts:
            names = [_snp_name(p, a) for p, _rf, a in alts]
            fr[17] = ",".join([indel_data[4]]
                              + ["%s_%d_%s_%s" % (current_chr, p, rf, a) for p, rf, a in alts])
            fr[15] = ",".join([indel_data[2]] + [n[0] for n in names])
            fr[16] = ",".join([indel_data[3]] + [n[1] for n in names])
        used = tuple(sorted(
            (cols[k].pos1, _isc.complement(r["bases"][k]) if strand == "-" else r["bases"][k])
            for k in r["enabling"] if cols[k].pos1 > 0))
        if used and _cooc_out is not None:
            ent = cooc.setdefault((r["phase"], used, r["samples"] is None), [set(), 0])
            ent[0] |= set(r["samples"] or ())
            ent[1] += r["slots"]
        out.append(fr)
    for (phase, used, sites_only), (samples, slots) in sorted(cooc.items()):
        head = (current_chr + "\t" + str(_ipos) + "\t" + _iref + "\t" + _ialt + "\t"
                + final_result[4] + "\t" + strand + "\t"
                + ";".join(str(p) for p, _a in used) + "\t"
                + ";".join(_snp_name(p, a)[0] for p, a in used) + "\t" + phase + "\t")
        if sites_only:
            # the 2.5.2 fallback: no genotypes to count cis carriers with, so the joint
            # frequency is bounded by the rarest participating allele
            _iaf = None
            if _indel_af is not None:
                try:
                    _iar = _indel_af.lookup(_ipos, _iref, _ialt)
                    _iaf = _iar.get("AF_max") if _iar else None
                except Exception:  # noqa: BLE001
                    _iaf = None
            _snp_afs = []
            for p, a in used:
                g = _snp_global(p, a)
                _snp_afs.append(g.allele_freq() if g is not None else None)
            _cooc_out.write(head + ("%.6g" % _isc.putative_joint_af([_iaf] + _snp_afs))
                            + "\tNA\tNA\n")
            continue
        g = _snp_global(*used[0])
        _cooc_out.write(head + ("%.6g" % _isc.joint_af(slots, g.AN if g is not None else 0))
                        + "\t" + str(len(samples)) + "\t" + ",".join(sorted(samples)[:50])
                        + "\n")
    return out


def _second_indel_rows(line, final_result, indel_data):
    """A target with no IUPAC code is the same sequence on every slot that carries
    its indel -- except a slot that also carries a second indel inside the target,
    which the one-indel fake-indel window cannot hold. With the phased indel store,
    such carriers leave the row (the log's list, in its order); a row left with no
    carrier is dropped. Without a second indel there (nearly always) the row is the
    one the engine reported, unchanged."""
    if _indelgt is None or not _indel_snp:
        return [final_result]
    _ipos, _iref, _ialt, _istart = _indel_window(indel_data)
    real_of = _isc.build_offset_to_real(
        line[2], int(line[4]), line[6], int(indel_data[5]), _istart, _iref, _ialt)
    span = [r for r in real_of if r is not None]
    others = _other_indels(min(span), max(span), _ipos, _iref, _ialt) if span else ()
    if not others:
        return [final_result]
    carriers = _indel_carriers(indel_data, _ipos, _iref, _ialt)
    exact, putative = _isc.resolve_carriers([], carriers, others=others)
    kept = set().union(*[e[0] for e in exact.values()], *putative.values())
    gone = set(carriers.parsed) - kept
    if not gone:
        return [final_result]
    _ovl_stats["second_indel"] += 1
    left = [s for s in final_result[12].split(",") if s and s not in gone]
    if not left:
        return []
    fr = final_result.copy()
    fr[12] = ",".join(left)
    return [fr]


def _report_overlaid():
    """One stdout line on how the overlaid targets were resolved (stdout only:
    the pipeline treats any stderr write as a fatal stage failure)."""
    if _ovl_stats["second_indel"]:
        print("indel-snp: %d target(s) without a SNP code lost the carriers whose "
              "haplotype holds a second indel inside them"
              % _ovl_stats["second_indel"], flush=True)
    if not _ovl_stats["targets"]:
        return
    print("indel-snp: %s overlaid target(s) resolved per carrier haplotype into %d "
          "row(s) (%d PUTATIVE); %d had no carrier with a target there"
          % (_ovl_stats["targets"], _ovl_stats["rows"], _ovl_stats["putative"],
             _ovl_stats["dropped"]), flush=True)
    if _ovl_stats["unmapped"]:
        print("WARNING [indel-snp]: %d IUPAC column(s) of %s did not map onto the "
              "reference base they code; they were resolved without genotypes"
              % (_ovl_stats["unmapped"], current_chr), flush=True)


def _queue_row(final_result):
    """Reference alignment, PAM creation and the scoring sentinels for one indel
    row, then queue it for scoring."""
    ref_aligned = alignRefFromVar(final_result)
    # PAM creation on the INDEL path. Mirrors the SNP path's predicate in
    # new_simple_analysis.py: if ANY reference base in the PAM window is
    # incompatible with the PAM motif, the reference could not have carried a
    # usable PAM at this locus, so the variant created it. The stored value is the
    # ALT PAM string, exactly as the SNP path stores it.
    #
    # PAM_gen is index 10, the same slot as on the SNP path. The header that
    # analisi_indels_NNN.sh seds into the cluster file lists an extra
    # Chromosome_fake column at 10, which would put PAM_gen at 11 -- but that
    # header is misaligned with the rows the analyzer actually reads. The rows
    # here are 18 fields with no Chromosome_fake, which is also what this
    # function's own assignments assume: [12] samples, [15] rsID, [16] MAF,
    # [17] indel info all line up with the SNP layout and would each be off by
    # one against the sed header. Do not re-derive this index from that header.
    #
    # pam_begin/pam_end already encode PAM orientation (trailing window for a 3'
    # PAM, leading for a 5' one) and `pam` is the bare motif aligned with that
    # window, so this is correct for both. Left as the existing placeholder when
    # the window holds a character the IUPAC table does not define -- e.g. a '-'
    # from a bulge falling inside the PAM -- so this can never raise where the
    # column was previously simply unset.
    # .upper(): iupac_code_set maps lowercase keys to LOWERCASE sets ("a": {"a"}),
    # so a lowercase reference base would never intersect an uppercase motif set and
    # would be scored as a creation event. In practice the PAM window is always
    # uppercase (the guide carries N there, so those positions always "match" and
    # alignRefFromVar leaves them uppercase), making this defensive only -- but it
    # removes the false-positive class outright rather than relying on that.
    ref_pam = ref_aligned[pam_begin:pam_end].upper()
    if len(ref_pam) == len(pam) and all(c in iupac_code_set for c in ref_pam):
        for pos_pam, pam_char in enumerate(ref_pam):
            if not iupac_code_set[pam[pos_pam]] & iupac_code_set[pam_char]:
                final_result[10] = final_result[2][pam_begin:pam_end].upper()
                break
    final_result.append(ref_aligned)
    # number to activate ref score calculation (active if target is alternative)
    final_result.append(33)
    # position of tmp_mms (removed later after processing)
    final_result.append(0)
    cluster_to_save.append(final_result)


print("Analysis of " + current_chr)

save_cluster_targets = True
remove_iupac = False
save_total_general_table = False
# chiave ['add'] = target semicommon da aggiungere alla tab generale delle guide, metto i valori di total presi dal primo target REF nel cluster di un semicommon che esiste
add_to_general_table = dict()
# chiave ['distribution'] = array len total, ogni cella divisa per mm,1B,2B..., per populationDistribution
last_annotation = ""  # in un cluster, l'annotazione è la stessa, quindi la calcolo solo una volta e poi la riscrivo per gli altri target del cluster
last_samples = (
    []
)  # se due target hanno la stessa scomposizione, hanno gli stessi samples, quindi non li ricalcolo
next(inResult)  # Skip header
# contains all targets for each cluster, it's then sorted by CFD. Only the target with highest CFD is saved
# cluster_to_save = []
# sum number times a cfd value appears, TODO aggiungere anche il conteggio per ogni sample
cfd_for_graph = {"ref": [0] * 101, "var": [0] * 101}

# Categorie per la scomposizione, chiavi: 'BulgeType' + 'BulgeSize'
dictionary_entries = ["X0"]
for i in range(int(max_dna_bulges)):
    dictionary_entries.append("DNA" + str(i + 1))
for i in range(int(max_rna_bulges)):
    dictionary_entries.append("RNA" + str(i + 1))


# dict_ref_seq = convert_fasta_to_dict(sys.argv[12])

cluster_class = None
datastore = None


global_start = time.time()

# list saving all the targets reported in chr##
cluster_to_save = list()
# rows carrying another contig's label (see the check in the loop); counted so the
# run reports them instead of dropping them silently
foreign_rows = 0

for line in inResult:
    # print("line in INDEL analysis", line)

    line = line.strip().split("\t")
    # print(line)
    guide_no_bulge = line[1].replace("-", "")
    # copy line to avoid problem with memory managament and intrisic pointers
    final_result = line.copy()
    # The indel search runs per fake contig and post_analisi_indel.sh subsets the
    # merged targets by "fake${chrom}", so every row reaching THIS invocation must
    # carry THIS chromosome's contig label. Verify it before trusting the lookup
    # below: the fake coordinate space is shared across chromosomes (each starts at
    # 0 with near-uniform windows), so a row from another contig resolves to a REAL
    # indel here and would be reported with that indel's rsID, AF and carrier
    # samples rather than rejected.
    if line[3] != "fake" + current_chr:
        foreign_rows += 1
        continue
    # extract indel data from the INDEL tree
    try:
        indel_data = sorted(INDELS_tree[int(line[4])])[0].data
    except:
        # target found has no INDEL, so it's REF and must be skipped
        continue
    # assign extracted data to final target
    final_result[3] = current_chr  # current analyzed chr
    final_result[12] = indel_data[1]  # samples
    final_result[15] = indel_data[2]  # rsID
    final_result[16] = indel_data[3]  # MAF
    final_result[17] = indel_data[4]  # INDELinfo

    # correct position if PAM at beginning and strand
    if pam_at_beginning:
        if line[0] == "RNA" and line[6] == "-":
            # cluster_position
            fake_start_target = int(line[5]) - int(indel_data[5])
        else:
            # real_position
            fake_start_target = int(line[4]) - int(indel_data[5])
    else:
        if line[0] == "RNA" and line[6] == "+":
            # cluster_position
            fake_start_target = int(line[5]) - int(indel_data[5])
        else:
            # real_position
            fake_start_target = int(line[4]) - int(indel_data[5])

    # correct start with real chr spatial info and length of the target
    true_start_target = (
        int(indel_data[0].split("_")[1].split("-")[0]) + fake_start_target
    )
    # correct cluster position
    diff_pos_clus = int(line[4]) - int(line[5])
    # real_start
    final_result[4] = str(true_start_target)
    # real_start_cluster
    final_result[5] = str(true_start_target - diff_pos_clus)

    # [mega] per-dataset INDEL AF sidecar lookup -> companion row (guarded no-op if
    # the sidecar is absent). indel_data[4] is "chrN_pos_REF_ALT"; the SNP-only
    # Tier-0 registry cannot hold the indel, so its per-dataset AF comes from here.
    if _indel_af is not None:
        try:
            _iap = indel_data[4].split("_")  # chrN_pos_REF_ALT
            _iarec = _indel_af.lookup(int(_iap[-3]), _iap[-2], _iap[-1])
            if _iarec is not None:
                _iavals = [_iarec.get(_l) for _l in _indel_af._labels]
                _indel_af_out.write(
                    current_chr + "\t" + final_result[4] + "\t" + _iap[-3] + "\t"
                    + _iap[-2] + "\t" + _iap[-1] + "\t"
                    + "\t".join("." if _v is None else ("%.6g" % _v) for _v in _iavals)
                    + "\t" + ("." if _iarec.get("AF_max") is None
                              else ("%.6g" % _iarec["AF_max"])) + "\n")
        except Exception as _iae:  # noqa: BLE001 - never break the classic indel row
            _indel_af_errs[0] += 1
            if _indel_af_errs[0] == 1:
                print("WARNING [mega-indel-af]: emission error (further "
                      "occurrences suppressed): %s" % _iae, flush=True)

    # [indel-snp] An overlaid target (IUPAC codes in line[2], the SNP overlay of the
    # fake-indel genome) stands for every combination of its SNP alleles at once:
    # resolve it per carrier haplotype, one row per sequence a carrier really has.
    # A target with no code is the same sequence for every carrier of the indel and is
    # reported as the engine found it, less any carrier whose haplotype holds a second
    # indel inside it (_second_indel_rows; a classic, non-overlaid build has no phased
    # indel store, so its rows are untouched).
    _amb = _isc.ambiguous_columns(line[2])
    for _row in (_overlaid_rows(line, final_result, indel_data, _amb) if _amb
                 else _second_indel_rows(line, final_result, indel_data)):
        _queue_row(_row)

    if len(cluster_to_save) >= 100000:
        # after reading 100k lines from file and creating the cluster, start processing it
        clusters_with_scores = calculate_scores(cluster_to_save)

        for count, cluster in enumerate(clusters_with_scores):
            for target in cluster:
                if count == 0:  # CFD target
                    # remove count of tmp_mms
                    target.pop(-2)
                    # save CFD targets
                    cfd_best.write("\t".join(target) + "\t" + str(0) + "\n")
                    # save mm-bul targets
                    mmblg_best.write("\t".join(target) + "\t" + str(0) + "\n")
                if count == 1:  # CRISPR-Bulge target
                    # remove count of tmp_mms
                    target.pop(-2)
                    # save CRISPR-Bulge targets
                    crispr_bulge_best.write("\t".join(target) + "\t" + str(0) + "\n")
        cluster_to_save = list()

if len(cluster_to_save):
    pass
else:
    # skip processing and return empty file
    # close open files
    cfd_best.close()
    mmblg_best.close()
    crispr_bulge_best.close()
    if _cooc_out is not None:
        _cooc_out.close()
    if _indel_af_out is not None:
        _indel_af_out.close()
    # update header
    os.system(
        "sed -i '1s/.*/#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\tReference\tCFD_ref\tCFD\t#Seq_in_cluster/' "
        + outputFile
        + ".bestCFD_INDEL.txt"
    )
    os.system(
        "sed -i '1s/.*/#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\tReference\tCFD_ref\tCFD\t#Seq_in_cluster/' "
        + outputFile
        + ".bestCRISPR_BULGE_INDEL.txt"
    )
    os.system(
        "sed -i '1s/.*/#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\tReference\tCFD_ref\tCFD\t#Seq_in_cluster/' "
        + outputFile
        + ".bestmmblg_INDEL.txt"
    )
    _report_overlaid()
    # print exit messagge and exit 0
    print("Done", current_chr, time.time() - start_time)
    print("ANALYSIS COMPLETE IN", time.time() - global_start)
    exit(0)

# start processing scores after the whole target file is read
clusters_with_scores = calculate_scores(cluster_to_save)

for count, cluster in enumerate(clusters_with_scores):
    for target in cluster:
        # print(target)
        if count == 0:  # CFD target
            # remove count of tmp_mms
            target.pop(-2)
            # save CFD targets
            # print('CFD', target)
            cfd_best.write("\t".join(target) + "\t" + str(0) + "\n")
            # save mm-bul targets
            mmblg_best.write("\t".join(target) + "\t" + str(0) + "\n")
        if count == 1:  # CRISPR-Bulge target
            # remove count of tmp_mms
            target.pop(-2)
            # print('CRISPR-Bulge', target)
            # save CRISPR-Bulge targets
            crispr_bulge_best.write("\t".join(target) + "\t" + str(0) + "\n")

cfd_best.close()
mmblg_best.close()
crispr_bulge_best.close()
if _cooc_out is not None:
    _cooc_out.close()
if _indel_af_out is not None:
    _indel_af_out.close()

os.system(
    "sed -i '1s/.*/#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\tReference\tCFD_ref\tCFD\t#Seq_in_cluster/' "
    + outputFile
    + ".bestCFD_INDEL.txt"
)
os.system(
    "sed -i '1s/.*/#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\tReference\tCFD_ref\tCFD\t#Seq_in_cluster/' "
    + outputFile
    + ".bestCRISPR_BULGE_INDEL.txt"
)
os.system(
    "sed -i '1s/.*/#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\tReference\tCFD_ref\tCFD\t#Seq_in_cluster/' "
    + outputFile
    + ".bestmmblg_INDEL.txt"
)


if foreign_rows:
    # STDOUT only -- this pipeline treats any stderr write as a fatal stage
    # failure (`[ -s $logerror ]`), so warning on stderr would abort the search.
    print(
        "WARNING: %d indel target row(s) in the %s subset carried another "
        "contig's label and were dropped; they would otherwise have been "
        "reported against %s indels they do not belong to."
        % (foreign_rows, current_chr, current_chr)
    )
_report_overlaid()

print("Done", current_chr, time.time() - start_time)

print("ANALYSIS COMPLETE IN", time.time() - global_start)
