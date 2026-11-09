#!/usr/bin/env python3
"""An independent check of CRISPRme's indel (and MNV) off-target rows on real data.

    check_indel_rows.py --results <results dir> --reference <genome dir or FASTA>
                        --vcf <bgzipped VCF> [--vcf ...] | --vcf-dir <dir>
                        --pam <PAM file> --guide <guide or guide file>
                        --mm 4 --bDNA 1 --bRNA 1 [--max-total-edits 6]
                        [--samples <samplesID file> ...]
                        [--chrom chr22 | --region chr22:1-500000]
                        [--mode soundness|completeness|both] [--snp-rows report|off]
                        [--jobs N] [--known known.tsv] [--out failures.tsv]

WHAT IT IS. A second opinion on the indel path, built only from what a user has: the
reference FASTA, the VCF(s) the search read, the samples list, the guide, the PAM file,
the search budget and the result folder. It imports nothing but Python's standard
library and no code of CRISPRme's variant or indel machinery
(``PostProcess/indel_snp_cis.py``, ``site_alleles.py``, ``analisi_indels_NNN.py``,
``overlay_indel_snps*``, ``resultIntegrator``, ``remove_contiguous_samples``, the CRISPRitz
enricher, the dictionaries, the indel log or the fake contigs), and none of it is copied
here. It reads the VCF with its own BGZF + tabix reader, writes each sample's haplotypes
itself from the VCF 4.x definition of a record, keeps its own map from haplotype bases to
reference coordinates and aligns the guide with its own exhaustive aligner. The CFD
tables CRISPRme ships (``PostProcess/cfd_*_scores.json``) are read as data.

``crisprme.py validate-test`` compares the engine's SITES (chromosome, start, strand) on
the reference and on the SNV-only 1000 Genomes chr22 genome, so an indel row's
coordinate, its aligned strings, its REF/ALT text, its variants, its carriers and the
alternative-alignment table are never checked against an independent truth on real data.
This is that check.

THE TRUTH IS THE VCF, AS THE SEARCH READ IT. A sample's haplotype is the reference with
every ALT allele its genotype puts on that slot, each in place of its REF bases (VCF 4.x),
of the records CRISPRme reads; a record is read with its REF and ALT in upper case, and so
is a row's Variant_info (var_text: a row writes a record as its VCF does, some HPRC records
in lower case). CRISPRitz's enricher skips every line whose FILTER is
neither PASS nor '.' (those exact strings): no such line is in the enriched genome, its
dictionaries, registry and genotype store, the indel log or the fake-indel windows. The
phased indel store (build_indel_genotypes) holds the same records since CRISPRme 2.10.0;
one written before streamed the VCF without reading FILTER and kept every indel of it, and
CRISPRme's carrier resolution reads that store: a second indel or an overlapping one on a
carrier's slot takes it off a target, and Track B drops a SNP row's carrier for an indel
of its slot, whatever that indel's FILTER (the oracle's sites O, Q and Y1 with their
indels marked LowQual are resolved so by a store of 2.8.x). So of a genotyped dataset such
a line is read as the reference, but for its indels where the run read a store written
before 2.10.0 (RELEASES, below), which stay on their carriers' haplotypes as the store has
them; either way they are no record of the search (no fake-indel window holds one: a row
naming one fails, variant-filtered, and completeness demands none of their targets,
limit-filtered-record), and its SNVs and MNVs are read as the reference. Of a sites-only
dataset such a line keeps its SNVs for its SNP rows: the hg38_mega index's SNP rows read
them whatever their FILTER (its registry holds every SNV of its VCF, and its SNP rows name
some: chr22:11574550 C>G, LOW_VQSLOD), while an index CRISPRitz's enricher builds from a
sites-only VCF does not hold them, and a run does not say which it read -- so a SNP row
through one is a row CRISPRme may write or not (below). Its fake-indel windows code what
CRISPRitz's overlay codes from the enriched genome, its SNVs whose FILTER passes: every
allele of a split site on a genome of CRISPRitz 2.10.0, the last record of one on a genome
of 2.9.1 or earlier (the published mega's indel rows name no other: at chr22:25093750 the
PASS C>G, before the LOW_VQSLOD C>T), so a row of it (no carrier) names SNVs of its own
VCF whose FILTER passes (else variant-not-in-dataset, variant-filtered; on a genome of
2.9.1, an earlier record of a split site is the note variant-shadowed). Its indels and
MNVs outside PASS and '.' are dropped (hg38_mega's indel log holds exactly the PASS indels
of its VCF). A slot holds a row's target when it reads the target's bases from the row's
start THROUGH the record: each base of the target the record's ALT puts there at its place
(and the base after them), as the row's target has it -- a second insertion between the
start and the record that lengthens the stretch and pushes some of the record's bases out
of it makes another target, which needs two indel-path records. Two alleles of one slot
whose REF spans meet (a SNP at an insertion's own position, a SNP on a base a deletion
removes, a complex allele split into overlapping records, two ALTs at one position phased
onto one haplotype) describe one stretch twice: the VCF gives no sequence there. Where
every allele that meets the record's REF is one METHODS section 4 reads the record's slot
without -- another allele at the record's own position (a phasing artefact), an allele
only on bases the record deletes, an MNV (in no store CRISPRme reads beside the record) --
and no allele kept meets one left out, the slot is read without them: the one reading of
it that holds the record, so the one a row naming the record (and not them) can mean. A
listed carrier must hold the target there, and is a note when it does
(overlap-reading-agrees); a carrier left off whose slot holds it there is counted
(limit-overlapping-records: METHODS says such carriers stay listed; not turned into a
failure). Any other overlap of the record (an indel or a SNV at another position: a
complex allele, which METHODS takes off, or an undocumented one) leaves the slot no
reading. Two other alleles of the slot whose REFs meet (an overlap zone) leave it no
sequence only where they can change what it says of the target -- one of them is named, or
changes a base read from the target's start (no margin: a zone wholly before or after that
stretch changes none of its bases) --: every reading of such zones (any of their alleles
no two of which meet, none included) is tried, and where all of them agree the slot is
judged so (a carrier whose slot holds the target on every reading is required, one whose
slot holds it on none fails when listed); where they part, the slot has no reading. A slot
with no reading gives no truth there, so CRISPRme's documented rules are the only model: a
listed carrier they take off this target for a documented cause (a complex allele, a
second indel the slot cannot be read through, the last record of a split site, the
PUTATIVE rule) fails (carrier-off-by-rule), and any other is counted
(limit-overlapping-records). Whatever its overlaps, a slot that holds a named variant on
no reading (neither on it nor of unknown phase there) lacks the target. A sample's alleles
are on known slots when its calls are phased ('|', or haploid). Where a call that touches
the target is unphased, has a missing allele, or has fewer copies than the record's (a
haploid call against a diploid one), which of its alleles share the record's haplotype is
unknown -- and so is it where the record's slot holds two SNV ALTs of one position that
CRISPRme's stores both read (every allele of a split site coded; split records phased one by
one), which have no one sequence in the VCF and CRISPRme resolves PUTATIVE (decision B4): when
the record's slot is known and no placement of those calls lets it hold the
target (a phased second indel inside it, say), or every placement does, the slot is judged
so; otherwise the truth is CRISPRme's documented PUTATIVE rule (METHODS section 4): the
sample is listed for the strongest target its own alleles allow on the row's alignment --
the guide's base wherever it carries it and the alignment pairs the position with a guide
base; elsewhere (a DNA bulge included), and at the PAM jointly, the alleles that maximize
the alignment's CFD among those the motif allows (a product: where it is 0 for every
choice, a column no choice changes scoring 0, every choice is a strongest one; METHODS
orders no CFD tie, so a tie is left open); a search that scores no CFD (a 5' PAM, or a
guide or PAM not 20 and 3 nt long) takes there the reference base where the alleles (and,
at the PAM, the motif) allow it, else the first of A, C, G, T they allow, read in the
target's orientation -- and never for an allele it does not carry; its copies are counted
against the record's call (a haploid ALT call beside a diploid record leaves the reference
possible). The rule chooses only at the fake-indel window's code columns: a reference base
of the target off the record's own REF span. The window writes the ALT verbatim, the
anchor included, and its overlay codes flank bases only, so a base on the record's REF
span is never a choice, whatever SNV the VCF holds there: a call of unknown phase whose
REF meets the record's would leave the record's slot with no sequence, and the one reading
that gives it one (the documented multi-allelic heterozygote, one allele per haplotype)
puts it on the other slot.

THE DOCUMENTED VIEW ONLY CLASSIFIES. Where CRISPRme's output and the truth part, the
rules CRISPRme documents for reading a carrier (METHODS section 4, CHANGELOG 2.8.1 Known
limits), re-implemented here from that text, say whether the difference is a documented
limit (counted and listed under its own class, non-gating) or not (gating):
  * the fake-indel window holds one indel-path record (an indel or an MNV) and the SNPs of
    the enriched genome -- every allele of a split site on a genome of CRISPRitz 2.10.0,
    where a slot holding two ALTs of one position has no sequence and is resolved PUTATIVE
    (the truth there, as for a call of unknown phase: THE TRUTH), and two listings of one
    REF>ALT at a position are one allele, on a line of several ALTs too; on one of 2.9.1 or
    earlier, only the last of several SNP records at one position, of the same REF>ALT too,
    or a last line of a '*' alone (limit-last-alt-wins) (RELEASES) --, and 25 reference
    bases before the record's first
    REF base and 27 after its last. An index keeps a site, a PAM and the guide's length of
    bases from it (its guide region), not an alignment: a target that runs past the window,
    or whose site's guide region does (one with an RNA bulge reaches past the target's own
    end), is held by no index (limit-window-reach; a 27-nt target with a DNA bulge,
    SaCas9's say), and one whose site's guide region ends within the INDELS index's budget
    of the window's edge (the pad the budget's DNA bulges use) by none of format 1
    (limit-window-edge, an index of CRISPRitz 2.9.1 or earlier: 2.10.0's format 2 holds
    it);
  * the window's target is resolved with each slot's own SNP alleles (a row lists the
    haplotypes left with its sequence and names the SNP alts they hold), so a slot is a
    candidate for a row only where its own alleles give the row's target; a second phased
    indel inside the target is then read on the slot's bases, from the target's first
    base on or up to its last base, through the record (each base the record puts in the
    target at its place, and the base after them, in the stretch read: an indel that
    lengthens the stretch before the record pushes them out), before any call of unknown
    phase sends the carrier to the PUTATIVE rule; an unreadable or different reading
    takes the slot off, and a target that needs two indel-path records is never reported
    (limit-two-indels) -- unless CRISPRme's reading of the slot against the window's
    target of the same bases, start and variants holds (a repeat unit gained or lost at
    the target's edge puts them back), where it keeps the carrier and completeness
    demands that row;
  * another allele at the record's own position is a phasing artefact, a variant only on
    bases the record deletes changes nothing, and an overlapping indel at another position
    takes the slot off (limit-overlapping-records; a listed carrier the rule takes off is
    carrier-off-by-rule, gating) -- the truth reads the first two as the rule does, so an
    overlap explains a disagreement only on its own slot, as the clash it is; a reading
    of the slot through a second indel that meets two of its variants that overlap (or a
    call of unknown phase) cannot be read, either way;
  * the carriers of an MNV, of an alt of a multi-ALT SNP line whose number the genotypes
    cannot tell (the number every carrier of the alt holds, less those the position's other
    alts of one character of a line of several ALTs settle -- a '*' among them, an indel's
    ALT in no dictionary entry, an alt whose carriers hold only 0 and 1 settling none: its
    1 is its own line's --, until nothing changes: none or several left), and of a second
    indel with an
    unphased call are resolved PUTATIVE (limit-putative-carrier where that differs from
    their known haplotype -- the target a known slot holds listed only under the
    alignments where the rule picks it, worse than its best, by the best of those -- or
    where the rule ties the target with another and the row leaves the carrier off);
  * an MNV beside the record on the slot is in neither store CRISPRme reads (the phased
    indel store keeps records whose REF and ALT lengths differ, the SNP tiers SNVs): a
    target that needs it needs two indel-path records (limit-two-indels), but a carrier
    listed for a sequence such an MNV changes is listed for a sequence it does not have,
    which nothing documents (variant-missing, gating).
  A disagreement no documented rule explains gates (unexplained-disagreement: a defect of
  CRISPRme or of this checker).

TWO MODES.

  SOUNDNESS (every reported indel/MNV row is real). Each row block of the integrated
  results (highest-CFD, fewest mm+b, highest CRISPR-Bulge) and of the alternative
  alignments whose variant list names an indel or an MNV:
    * is labelled REF/ALT_origin alt, names exactly one indel-path record (an indel or an
      MNV) and SNPs, every one a record of the VCF the search read (FILTER: a record it
      skips is variant-filtered), forward strand, on the row's chromosome;
    * its aligned spacer, gaps removed, is the guide (--guide, never the table's own
      column) followed (or, for a 5' PAM, preceded) by the PAM; its mismatches, bulge
      types and sizes and their total recount exactly, within the budget; no gap falls in
      the PAM; the PAM matches the searched motif; no IUPAC code is left;
    * its target (gaps removed, PAM included) is the reference with exactly the named
      variants applied, read from the reported start under the convention below, and
      holds a base of its indel; every named SNP changes a base of the target;
    * its REF text is the reference read from the start, gapped as the target;
    * every listed carrier is a sample of the panel (--samples) and holds it, and no
      carrier is missing (CARRIERS).
  An ALT-origin block that names no variant is reported (alt-without-variant). A row that
  names an earlier SNP record of a split site, which the documented limit says CRISPRme's
  stores do not hold, is a note (variant-shadowed, non-gating): the record is real and the
  target is checked with it.

  COMPLETENESS (no real indel off-target is missed -- the real-data brute force). For
  every indel and MNV record of the VCF in the checked region: each carrier's haplotypes
  around it are written (a known slot with its known alleles: a target that touches a
  call of unknown phase is the PUTATIVE rule's, below, and no other depends on where such
  a call is), identical sequences are aligned once, every in-budget alignment on both
  strands with the PAM that holds a base the record changes (the junction rule) is
  enumerated and mapped to the reference, and each (sample, target) must be reported for
  that sample at that site
  by its best alignment -- as an integrated row or an alternative alignment whose rows of
  the same site, variant list and carrier set are no worse than it -- or merged into a row
  of the same site by the documented merge into a row no worse than it. A target that
  touches a call of unknown phase is demanded as the strongest target the sample's alleles
  allow there (the PUTATIVE rule; a base on the record's own REF span is no choice: the
  window writes the ALT verbatim, the anchor included, and an allele there would describe
  the record's bases twice, or lie on bases it deletes).
  A run that writes no alternative alignments (.alt_mode = no-alt-alignments: by default
  a population-level run, as the documented chr22 self-test's, which --ci-cd-test does
  not make per-sample) keeps one row per merge cluster, so there each target's cluster
  must be reported with a fewest-mm+b block no worse than it and a highest-CFD block
  scoring at least its CFD. A carrier's slot whose overlap zones change a target on some
  reading and not on another demands nothing there: such a target is counted
  (limit-overlapping-records), never demanded.

  NOT CHECKED: the CFD_score_REF/ALT, CFD_risk_score and CRISPR-Bulge columns (and which
  row the CRISPR-Bulge block keeps), Variant_info_spacer+PAM, PAM_creation, the seed
  counts, Variant_MAF and Variant_rsID, Other_motifs, Not_found_in_REF, the annotations;
  a row's own CFD_score is not compared (the merge checks recompute each fold row's CFD
  from its aligned strings).

DEFINITIONS (fix/282-indel-alignment; METHODS section 4, "SNP+indel off-targets", and
the merge of ``remove_contiguous_samples.py``).

  Coordinate. A row starts at the 0-based reference coordinate of its target's FIRST
  ALIGNED BASE, the leftmost on either strand. A base an insertion adds has none: such a
  row starts at the reference base the added bases follow, the insertion's ANCHOR (its
  VCF POS - 1, 0-based). An ALT allele's bases are laid on the reference after the bases
  REF and ALT share at either end: the j-th remaining ALT base stands in the place of the
  j-th remaining REF base while there is one, and the rest are added after the last of
  them (after the shared base before them when none is replaced).

  REF text. ``Aligned_protospacer+PAM_REF`` is the reference read from the start for as
  many bases as the target holds, in the target's orientation, with the target's gap
  columns put back at the same indices (letter case is not compared).

  Bulges. A gap in the aligned spacer is a DNA bulge, a gap in the aligned target an RNA
  bulge. None falls in the PAM columns; the spacer's guide end away from the PAM is never
  a gap (a leading RNA bulge, and a DNA bulge right next to the PAM, are allowed).

  Junction rule. A target holds its record when it holds a base the ALT adds or
  replaces; for a pure deletion, when it holds the bases on both sides of the deleted
  ones, side by side. A target that holds none lies in a flank: reference or SNP sequence,
  which the reference and SNP searches report; it is no indel row, and completeness does
  not ask for it.

  CARRIERS (soundness). The panel is the samples of --samples (the samplesID files the
  search read), else of the run's .sampleID.txt (CRISPRme writes there the samples of
  those files), else every sample of the VCF(s): CRISPRme lists no carrier outside it. The
  pool is the listed samples and every panel carrier of the row's record. A listed sample
  outside the panel fails (carrier-not-in-panel). A listed sample whose haplotypes do not
  hold the target at its start with every named variant, or that is listed for a weaker
  target than its strongest (PUTATIVE), fails unless the documented view lists it (then a
  limit, with its cause); a sample whose haplotype holds the target at the start with the
  named variants, or whose strongest PUTATIVE target it is, fails when the row omits it,
  unless the documented view leaves it off (a limit). A slot whose alleles overlap the
  record is judged as THE TRUTH says: read without the alleles METHODS reads it without,
  where it gives that reading, and held to its named variants always. The merge never
  moves a carrier off its sequence's row (same position, variants and aligned target
  merge carriers; equal carrier sets fold to the best alignment), so the row itself must
  list it. A
  SITES-ONLY VCF (eight columns, no sample) genotypes nobody, and a sites-only dataset's
  rows list no carrier (CHANGELOG 2.8.1; its post-analysis is the population-level one, its
  rows the worst case of each window): so a row that lists no carrier for a record a
  sites-only VCF holds is that dataset's, and its strings, target and variants are checked,
  its carriers not (limit-no-genotypes, one line per row block). The decision is the
  record's, not the run's: beside a genotyped dataset, as the default index with
  hg38_mega, the genotyped dataset's rows list their carriers and are checked as above, and
  a record both hold is demanded for its genotyped carriers by completeness. A row that
  lists no carrier for a record no sites-only VCF holds fails (no-carrier-listed), and a
  sample listed for a record only sites-only VCFs hold carries nothing
  (carrier-without-indel).

  Merge (completeness). The pipeline clusters a guide's rows on a chromosome with Cluster
  Positions within --merge bp (3) of the previous, chained, whatever the strand.
  CRISPRitz's Cluster Position is the start moved to the PAM's fixed end (less the RNA
  bulges and plus the DNA bulges on the + strand with a 3' PAM, the - strand with a 5'
  one), which every alignment of one PAM site shares. Inside a cluster, variant rows of
  one (Cluster Position, variant list) -- and for an indel/MNV row one aligned target --
  merge their carriers, and the rows of one site, variant list and carrier set fold into
  the first of them in each file's order: the fewest file sorts on the sorting criteria in
  turn (--sorting-criteria, crisprme.py's default mm+bulges,mm), the CFD file on the CFD
  (as the analysis writes it, to three decimals) and then on the same criteria
  (merge_close_targets_cfd.sh hands both files --sorting-criteria; a later CRISPRme branch
  hands the CFD file --sorting-criteria-scoring, so where a run's command line names one, a
  CFD tie no worse under either is taken as no worse). Two alignments of one
  target share its start, variants, carriers and site, so they always fold. The site is the
  cluster position on the haplotype carrying the record: the Cluster Position plus the
  rank of the row's first base among the bases an insertion adds (0 on any other base), as
  the fix branch's merge compares rows that start inside one insertion (c3a75fe; at
  d7a6fff it compared the Cluster Position alone, which folds two targets of one
  insertion). CRISPRme's fold key holds no strand, while the alignments of one target
  share theirs: rows of the other strand are other targets, which a fold across strands
  would lose, so a fold group is taken on one strand. A fold group is NO WORSE than a
  target when its fewest-mm+b block's key under the sorting criteria is at most the least
  key of the target's alignments, and its highest-CFD block, its CFD recomputed from its
  own aligned strings, scores at least the target's best CFD at the three decimals the
  merge sorts on -- and, where the two surely print the same three decimals, has a key no
  greater than the least key of the target's alignments of that CFD. (Not compared: the
  merge's keys after the criteria -- the Cluster Position, the carrier count, the aligned
  text -- which only order alignments equal in all of these.) So a (sample, target) is
  satisfied by a row reporting
  exactly it with its variants whose fold group (the rows of its site naming those
  variants with its carrier set) is no worse than it, or by the rows of its site naming
  the same variants whose carrier set is the target's own as CRISPRme writes it and that
  are no worse than it: every carrier both the documented view and its own haplotypes
  list for it, and nobody else but carriers either lists or the view leaves open
  (PUTATIVE) -- the view only classifies, so one carrier it lists against the truth makes
  no veto; a carrier whose slot overlaps the record is never required, and allowed where
  its slot, read as METHODS reads it (without the overlapping alleles), holds the target,
  or, a slot with no reading (no one sequence in the VCF), where the view lists it or
  leaves it open. A PUTATIVE target is the carrier's own: it is
  held to the alignments under which the rule picks it for that carrier (the choice is
  the alignment's highest CFD; CRISPRme's merge keys an indel or MNV row by its aligned
  spacer as well as its target, so a PUTATIVE carrier stays on the alignment it was
  resolved on), and one the rule ties with others to its own best. The
  verdicts depend on nothing but the input (no set's order). In a run with one row per
  cluster, the cluster of
  a target is the integrated row with a block within --merge bp of its Cluster Position,
  or of one the rows CRISPRme writes chain to it, each within --merge bp of the next: the
  reference rows, the SNP rows (the mode is the one the run ran, its .search_mode: a
  population-level run's the worst case of each alignment over each dataset's SNP
  alleles; a per-sample run's only the haplotypes individuals hold, a call with fewer
  copies than the sample's others of unknown phase; the alleles its stores hold, the last
  record of a split site) and the targets some carrier of an indel-path record holds (the
  per-carrier resolution gives no row for a raw target no carrier holds, nor for one that
  needs two indel-path records; a slot whose overlapping alleles METHODS reads it without
  is read so, as a clean slot: its exact targets and the window's targets a second
  indel-path record leaves readable; one it takes off -- a complex allele -- gives none;
  a sites-only dataset's rows are the worst case of each alignment of its window, the
  record's ALT written verbatim and an SNV coded only off its REF span, and a record a
  genotyped and a sites-only dataset both hold gives the rows of both), every
  record the chain reaches read, up to 3 kb either side (a chain still going there with no
  row as good as the target is this checker's coverage gap, coverage-chain-limit, not a
  miss); there is no other fallback. The rows CRISPRme may write or not, which no
  documented rule settles -- the one-record reading of a slot whose overlap METHODS does
  not cover, a target an overlap zone of two other alleles of a slot changes (on the
  reading with none of them, or with some written: the window reads the stores' alleles
  there), a target of a slot read without its overlap that touches a call of unknown phase
  (resolved PUTATIVE), a target through an earlier record of a split site or an untied
  multi-ALT alt (its stores read another allele there), a per-sample run's capped
  windows' greedy worst-case rows (METHODS section 5's known limit: only the classic
  dictionary path writes them, and a run does not say which path it took), and a
  per-sample SNP target of a slot that also holds an indel of the phased store touching
  it (since 2.10.0 CRISPRme drops such a carrier from a SNP row unless its reading of the
  slot through a repeat keeps the target's bases, which the chain does not model) -- are
  followed only when the rows it writes reach no row as good as the target, and never one
  that needs a second indel-path record: a cluster that reaches one only through them is
  counted under the class of their cause (limit-overlapping-records, limit-last-alt-wins,
  limit-putative-carrier, limit-capped-window, coverage-snp-row-indel), neither passed nor
  failed. That row must be
  as good as the target in both blocks: no more edits by the first sorting criterion (the
  merge picks between a cluster's reference and variant rows on mm+bulges, or the CFD,
  alone) and at least its CFD. A row's CFD is always recomputed from its aligned strings:
  the merge sorts on the CFD of a row's own alignment (CFD_score_ALT), while CFD_score
  prints the higher of its REF and ALT scores.

LIMITS (non-gating, counted per class, one TSV line per row block or per (record,
target)): limit-last-alt-wins, limit-two-indels, limit-overlapping-records,
limit-putative-carrier (all above), limit-window-reach (in completeness, per record: the
carriers of a target in budget that runs past its fake-indel window, or whose site's guide
region does, held by no index), limit-window-edge (the same, of a target whose site an
INDELS index of format 1 does not hold, at the window's edge), limit-max-degeneracy (a
cluster reached only through a SNP row over the run's --max-degeneracy),
limit-capped-window (a cluster of a per-sample run of 2.8.x without alternative alignments
reached only through a capped window's greedy row, above),
limit-contig-end (a record within one target of a
contig's end, or a row whose REF text runs past it), limit-n-bases (a stretch of N in a
record's window: no target is aligned across it), limit-no-genotypes (a row of a
sites-only dataset, which lists no carrier: the row itself is checked, one line per
block; and in completeness one line per chromosome counting the records only sites-only
VCFs hold, which no sample carries), limit-filtered-record (in completeness one line per
chromosome counting the indel and MNV records whose FILTER the search skips, none of whose
targets is demanded; and a cluster reached only through a sites-only dataset's SNP row
through one of its SNVs outside PASS and '.', above), and in the SNP-row report
limit-snp-row-mnv (a
listed carrier whose slot of the row's SNPs holds an MNV that changes the target: the
phased indel store CRISPRme reads a SNP row's carriers against holds no MNV) and
limit-snp-row-unphased; variant-shadowed and overlap-reading-agrees (a listed carrier
whose slot, read as METHODS section 4 reads a slot whose alleles overlap the record,
holds the target: the row agrees with it) are notes. Five non-gating classes are this
checker's own COVERAGE GAPS, not CRISPRme limits: coverage-overlap-cap (a carrier whose
slot's overlap zones within a target's reach of a record have more than 1024 readings:
only each of their alleles alone is read for the targets only such readings hold, which
are counted, never demanded, and are rows CRISPRme may write or not that a merge cluster
can chain through), coverage-putative-cap (a carrier the
PUTATIVE rule resolves with more than 10 positions within a target's reach of a record,
off its REF span, where its alleles allow two bases: its targets that touch a call of
unknown phase are not enumerated, 2^k choices each), coverage-chain-limit (a
population-level cluster that runs past the 3 kb either side the chain follows, and
holds no row as good as the target within them: its row may lie further),
coverage-snp-row-indel (a cluster of a per-sample run without alternative alignments that
reaches a row as good as the target only through a SNP target of a slot that also holds
an indel touching it, a row CRISPRme writes only where its reading of the slot through a
repeat keeps the carrier, above) and coverage-unindexed-long-ref
(a VCF read without a .tbi holds a record whose REF is longer than the 20 kb an unindexed
read looks back: a deletion that long, starting further before a region, is not read
there; give the VCF a .tbi).
A target all of whose carriers are such limits is demanded of no one, so a row that lists
only them is not missed when it is lost: each completeness limit line says how many of its
carriers a row naming the record lists, and that count is what changes then (CRISPRme's
reading of those slots -- the alleles at a record's own position, the record a split
site's store keeps -- is documented only as 'stays listed' and 'the last record', and is
not modeled into a gating demand). Rows on a chromosome no VCF covers, or with no reference
sequence, and table lines shorter than their header are reported (gating): nothing is
skipped silently.

OUTPUT. A TSV of every finding with its reason class, gating flag and evidence, a
per-class summary on stdout, exit status 1 on any gating finding (2 on a usage error or
a configuration it does not model, below: nothing is checked then). The check's inputs
must be the search's: where the run's own record of them says otherwise -- its
.Params.txt (Mismatches, DNA, RNA, Max_total_edits, Pam, and Genome_idx, the variant
datasets it searched: one VCF per dataset must hold each chromosome checked, no fewer and
no more; a search without bulges records 'Genome_idx None' whatever it searched, and its
datasets are then the lines of the --vcf list its command line names, read in its
working directory, the folder that holds Results/<output>), its .guides.txt (a guide it
did not search; some of its guides are a check of
those, and the others are counted, one guide-not-given line each: their rows are not
checked and completeness demands none of their targets), its .sampleID.txt (against
--samples) and its .search_mode, else its command line's mode flag (against --analysis)
-- the run is refused (exit 2): it would be modeled partially. So is a run whose own
records contradict each other (a .search_mode per-sample beside --population-level, an
.alt_mode against --alt-alignments or --no-alt-alignments: no search writes them
together), a search of a genotyped dataset without its genotype store that says so
(its log_verbose.txt, or a --per-sample search that ran population-level: PRECONDITION),
and one that does not record what the check needs to model it: completeness
needs whether the run wrote alternative alignments (its .alt_mode, else its command line)
and, where it wrote none, the mode it ran (its .search_mode, else --per-sample or
--population-level on its command line; --ci-cd-test chooses no mode: complete-test's
searches take complete-search's context-aware default, population-level on the
self-test's index), which --analysis can give; the SNP-row report needs the mode too.
``--known`` lists findings already reported upstream: they turn non-gating
(known-defect), and a listed one that no longer shows is reported (known-not-seen).
SNP-only rows of a per-sample run can be checked in report-only mode
(``--snp-rows report``: their strings, their target and that every named variant is in
it, their listed carriers, and the phased carriers a SNP row of their site leaves out). A
carrier is read as CRISPRme's SNP rows read it (METHODS section 4): one with an SNV call
of unknown phase in the target ('/', homozygous too, or haploid) is listed for every
non-empty subset of its variants there (limit-snp-row-unphased); the merge joins the
carriers of a site's SNP rows of one variant list (its key for a SNP row is the Cluster
Position and the variants, no alignment), so a carrier may be listed for another
alignment of the site that its slot holds; and since 2.10.0 a carrier whose slot of the
row's SNPs holds an indel of the phased store that changes the target is dropped, unless
its own bases, read from the target's first base on or up to its last through the row's
SNPs, still are the target's, or the indel's call cannot be tied to the slot (kept,
PUTATIVE), or it sits at a named SNP's position (a phasing artefact) -- a slot whose
bases cannot be read (an indel near whose slot is unknown, two records of one base)
leaves the row. Listing a carrier the rule drops is reported (snp-row-indel-carrier); an
MNV the store does not hold, limit-snp-row-mnv. Population-level SNP rows are a documented
worst-case union and are out of scope, as are a sites-only dataset's SNP rows, which are
population-level whatever the search's mode (one limit-no-genotypes line each). Without
--max-total-edits the cap is the run's own
(its .Params.txt, else its command line), else complete-search's 6; without --samples the
panel is the run's
.sampleID.txt, else every sample of the VCF(s). Rows are keyed by guide: one guide's row
stands for no demand of another.

RELEASES. Where CRISPRme's and CRISPRitz's releases read the data otherwise, the run's own
records say which read it (Semantics, run_semantics), and a run whose records do not tell,
where its data makes it matter, is REFUSED (exit 2). A run of CRISPRme 2.10.0 or later
writes .indexes_searched.txt (the indexes that served it; and .carrier_gate.json); a run
without them is an earlier CRISPRme's, whose CRISPRitz (2.9.1 or earlier) coded the last
record of a split SNP site and built indexes of format 1, and whose phased indel store
kept the indels of every FILTER. Of a run of 2.10.0: the variant genome codes every allele
of a split site where the variant index the record names flags Q in its .bin headers, or
where its enriched genome (Genomes/<genome>+<dataset>/<chrom>.enriched.fa of the run's
working directory) codes an earlier record's ALT at a split site of the VCF; the last
record alone where that index is of format 1, or the enriched genome holds the last ALT
alone there, or an index of format 2 flags no Q while the VCF has a position of four bases
(refused where nothing tells and the VCF has a split site, two or more lines of an SNV ALT
at one position, a line of several ALTs among them, or such a line before a last line of a
'*' alone; a line whose ALT is a '*' is no record of one that tells: no code of the genome
holds it); every allele coded, two listings of one SNV
REF>ALT at a position are one allele, on a line of several ALTs too, their calls merged as
CRISPRitz 2.10.0's enricher merges them into the allele's one dictionary entry, which
CRISPRme reads (its CHANGELOG and its docs/MULTIALLELIC_ENRICHMENT.md section 2.3: a line of
several ALTs' calls of it first written as a line of one ALT writes them; phased alike --
as many alleles, both phased -- a slot holds the ALT where either line puts it; otherwise
the call with more copies of the ALT, the earlier line's on a tie, as it is), and a slot
holding two ALTs of one position is
resolved PUTATIVE (METHODS section 4), the truth there; with the last record alone coded, an
earlier line of the same REF>ALT is read as shadowed too (the dictionary entry was the last
line's, the earlier line's carriers lost: limit-last-alt-wins), and so is every line before
a last line of a '*' alone (2.9.1's dictionary entry of the position is the '*''s, while its
genome keeps the code the earlier lines wrote: a '*' adds no base, and one base writes no
code). A run of every allele whose VCF holds, near a row or a record it checks, two or more
lines that each list several ALTs at one position is refused (exit 2): each entry of the
position's dictionary value keeps its line's raw allele numbers, which CRISPRme's readers
cannot tell apart (METHODS section 4), and what CRISPRme then lists is not modeled. Two
lines of one indel or MNV REF>ALT at a position are one allele on every release,
their calls merged by the same rule (decision B2: CRISPRme 2.10.0's readers of the phased
indel store merge them so; one before read the later line alone, a defect the truth
catches). The INDELS
index's format is that of its .bin files, else its log's format-1 NOTE of it (refused
where neither tells and a target of the guide, PAM and budget can lie at a window's edge).
The phased indel store holds only the records whose FILTER passes but where the run's
log_verbose.txt notes a store written before 2.10.0, for its dataset and chromosome
(refused without the log, where a genotyped VCF holds an indel record of another FILTER).
A per-sample run of 2.10.0 reports a window over CRISPRME_IUPAC_CAP by its carriers'
haplotypes on the classic dictionaries too, so no greedy row of a capped window chains a
cluster (limit-capped-window is 2.8.x's). --max-degeneracy K (2.10.0; its .Params.txt
Max_degeneracy) cuts back or leaves out a SNP row resting on more than K ambiguous
columns, so in a cluster chain an alignment over more than K positions of several options
is a row CRISPRme may write or not (limit-max-degeneracy); an indel row is not cut. The
carrier gate (2.10.0) ranks a population-level row and labels it in its companion and the
report, leaving the integrated results and their alternative alignments as they are: it
changes nothing this checker reads. An alt of a VCF line that lists several ALTs whose
number no genotype ties is resolved PUTATIVE (untied), as by 2.8.x's indel path, now by
every reader of CRISPRme 2.10.0, which ties numbers over every alt of one character the
dictionary lists at the position (a '*' too). A run whose log says that the SNP+indel
co-occurrence was off (CRISPRME_INDEL_SNP=0, the search's NOTE since 2.10.0) is REFUSED:
every SNP code of its fake-indel genomes read as the reference, its indel rows name no SNP
and list every carrier of their indel, and its SNP rows keep a carrier whose haplotype
holds an indel in the target -- none of the per-carrier resolution modeled here (a run of
an earlier CRISPRme made so records nothing that tells it: check only runs made with the
default).

PRECONDITION. The search read each genotyped dataset's genotypes (its Tier-1 genotype
store, or classic dictionaries), as every index CRISPRme installs or builds with
--samplesID does: it resolves each carrier of an indel row on its own haplotypes, which
is what this checker models. An install without them (download --what index
--no-genotypes) searches population-level and lists every carrier of the indel on a
SNP+indel row, with or without its SNPs (CHANGELOG 2.8.1). A search of it made with
--per-sample records that (its .search_mode population-level beside --per-sample on its
command line: complete-search writes that pair only where every dataset it searches is
sites-only or where a genotyped one has no genotype store), and with a genotyped VCF
given it is REFUSED (exit 2). In any mode the run's log_verbose.txt, the pipeline's
stdout kept in the result folder, records it too: its SNP post-analysis says
'Registry-only install for <chrom>: ...' for such a dataset (hg38_mega's aggregate
registry says so as well, marked '(aggregate info_af registry'; not counted), and with a
genotyped VCF given the run is REFUSED (exit 2). Without that line (a result folder
without its log_verbose.txt, or a dataset its SNP post-analysis never read), a run made
with --population-level or no mode flag records nothing else that tells it apart (its
.search_mode, .command_line.txt and .Params.txt read as a search that read the store;
complete-search's warning goes to the terminal only): check only runs made with the
genotypes the VCF holds.

VCF ACCESS. A .tbi beside a VCF makes every read a seek. Without one, each VCF is read
once into an in-memory block index (nothing is written), which costs one decompression of
the file; a canonical install (one VCF per chromosome, 28 GB) is then feasible with
--jobs. An unindexed read starts 20 kb before its region, so give a .tbi for a panel
with deletions longer than that (HPRC's chr22 has 582; a check lists each VCF whose reads
met one, coverage-unindexed-long-ref). Several VCFs of one chromosome are several
datasets, and each sample's haplotypes are read from the one that genotypes it: a
genotyped one beside a sites-only one is the default index with hg38_mega, and genotyped
datasets with no sample in common are checked so. Two genotyped VCFs of one chromosome
that share a sample of the panel are REFUSED (exit 2): CRISPRme resolves each dataset on
its own genotypes and its merge joins their rows, so a shared sample can be listed by
either dataset's haplotypes, which this checker does not model (run a search of each
dataset alone and check each: such a run given one of its VCFs is refused too, OUTPUT). A
sites-only VCF's split sites are not read as shadowed by its SNP rows: CRISPRme's mega
index holds earlier records of some (its SNP rows name them, as chr22:11337076 G>C before
G>T, both PASS), so which it holds cannot be read off the VCF's order; its fake-indel
windows hold every record whose FILTER passes on a genome of CRISPRitz 2.10.0, the last such
record on one of 2.9.1 or earlier (above).
``--jobs N`` splits the work by region over N processes.
"""

import argparse
import bisect
import csv
import glob
import gzip
import io
import itertools
import json
import os
import re
import struct
import sys
import zipfile
import zlib
from collections import defaultdict

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
# the FILTER values CRISPRitz's enricher reads (Enrichment/enricher.py and enricher.cpp:
# `if line[6] not in ('PASS', '.')`, the line skipped; PostProcess/validate_inputs.py's
# ENRICHER_PASS_VALUES): no line outside them is in the enriched genome, its dictionaries,
# its registry and genotype store, its indel log or its fake-indel windows
PASS_FILTERS = ("PASS", ".")
IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T", "U": "T",
    "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}
BLOCKS = ("highest_CFD", "fewest_mm+b", "highest_CRISPR_BULGE")
# the fake-indel window CRISPRitz's enricher writes for a record (enricher.py: the bases
# genome[POS-26 : POS+26+len(REF)], the REF in it replaced by the ALT): 25 reference bases
# before the record's first REF base and 27 after its last, a separator N on either side
WINDOW_BEFORE, WINDOW_AFTER = 25, 27
PUTATIVE_CAP = 10  # positions of unknown phase enumerated per carrier and record
IUPAC_CAP = 10     # CRISPRME_IUPAC_CAP's default: a window over it is capped (METHODS 5)
NO_VARIANT = ("", "NA", "n")


def rc(s):
    return s.translate(COMP)[::-1]


def comp(b):
    return b.translate(COMP)


def iupac_match(motif_base, base):
    return base.upper() in IUPAC.get(motif_base.upper(), "")


# the target bases each spacer character (a guide base or an IUPAC PAM code, upper case)
# matches, for the aligner's inner loop
_MATCHES = {m: frozenset(v) for m, v in IUPAC.items()}


# --------------------------------------------------------------------------- #
# BGZF + tabix: read a VCF region without pysam (an index is optional)
# --------------------------------------------------------------------------- #
class Bgzf(object):
    """Just enough BGZF to start reading at a virtual offset. A BGZF file is a series of
    gzip members, each carrying its size in a BC extra field."""

    def __init__(self, path):
        self.path = path
        self.fh = open(path, "rb")
        head = self.fh.read(18)
        self.fh.seek(0)
        self.bgzf = (len(head) >= 18 and head[:2] == b"\x1f\x8b" and (head[3] & 4)
                     and head[12:14] == b"BC")

    def close(self):
        self.fh.close()

    def lines(self, voffset=0):
        """Decoded lines from a virtual offset (coffset << 16 | uoffset). A plain gzip
        file has no blocks to seek to, so it is streamed from its start."""
        if not self.bgzf:
            with gzip.open(self.path, "rt") as fh:
                for line in fh:
                    yield line.rstrip("\n")
            return
        self.fh.seek(voffset >> 16)
        skip = voffset & 0xFFFF
        tail = b""
        while True:
            block = self._member()
            if block is None:
                break
            if skip:
                block = block[skip:]
                skip = 0
            buf = tail + block
            parts = buf.split(b"\n")
            tail = parts.pop()
            for p in parts:
                yield p.decode("utf-8", "replace")
        if tail:
            yield tail.decode("utf-8", "replace")

    def lines_at(self, voffset):
        """(virtual offset, line) of each line from a virtual offset on: a line's virtual
        offset is its first byte's (the block's compressed offset << 16 | the byte's
        offset in the block), as a tabix chunk names its ends."""
        self.fh.seek(voffset >> 16)
        pos = voffset & 0xFFFF
        tail, tail_v = b"", None
        while True:
            at = self.fh.tell()
            block = self._member()
            if block is None:
                break
            while True:
                nl = block.find(b"\n", pos)
                if nl < 0:
                    if pos < len(block):
                        if not tail:
                            tail_v = (at << 16) | pos
                        tail += block[pos:]
                    break
                if tail:
                    line, v = tail + block[pos:nl], tail_v
                    tail, tail_v = b"", None
                else:
                    line, v = block[pos:nl], (at << 16) | pos
                yield v, line.decode("utf-8", "replace")
                pos = nl + 1
            pos = 0
        if tail:
            yield tail_v, tail.decode("utf-8", "replace")

    def blocks(self):
        """(compressed offset, decompressed bytes) of every block, in order."""
        self.fh.seek(0)
        while True:
            at = self.fh.tell()
            block = self._member()
            if block is None:
                return
            yield at, block

    def _member(self):
        head = self.fh.read(12)
        if len(head) < 12:
            return None
        if head[:2] != b"\x1f\x8b" or not (head[3] & 4):
            raise ValueError("not a BGZF member in %s" % self.path)
        xlen = struct.unpack("<H", head[10:12])[0]
        extra = self.fh.read(xlen)
        bsize = None
        i = 0
        while i + 4 <= xlen:
            slen = struct.unpack("<H", extra[i + 2:i + 4])[0]
            if extra[i] == 66 and extra[i + 1] == 67 and slen == 2:
                bsize = struct.unpack("<H", extra[i + 4:i + 6])[0] + 1
            i += 4 + slen
        if bsize is None:
            raise ValueError("a gzip member without BGZF's BC field in %s" % self.path)
        payload = self.fh.read(bsize - 12 - xlen)
        return zlib.decompressobj(-15).decompress(payload[:-8])


def _reg2bins(beg, end):
    """The tabix bins that can hold a feature overlapping [beg, end) (0-based)."""
    out = [0]
    end -= 1
    for sh, start in ((26, 1), (23, 9), (20, 73), (17, 585), (14, 4681)):
        out.extend(range(start + (beg >> sh), start + (end >> sh) + 1))
    return out


class Tabix(object):
    """A tabix (.tbi) index, read here rather than through pysam."""

    def __init__(self, path):
        with gzip.open(path, "rb") as fh:
            buf = fh.read()
        if buf[:4] != b"TBI\x01":
            raise ValueError("not a tabix index: %s" % path)
        n_ref, _fmt, _cs, _cb, _ce, _meta, _skip, l_nm = struct.unpack("<8i", buf[4:36])
        names = buf[36:36 + l_nm].split(b"\x00")
        self.refs = {}
        off = 36 + l_nm
        for k in range(n_ref):
            n_bin = struct.unpack("<i", buf[off:off + 4])[0]
            off += 4
            bins = {}
            for _ in range(n_bin):
                b, n_chunk = struct.unpack("<Ii", buf[off:off + 8])
                off += 8
                chunks = []
                for _c in range(n_chunk):
                    cb, ce = struct.unpack("<QQ", buf[off:off + 16])
                    off += 16
                    chunks.append((cb, ce))
                bins[b] = chunks
            n_intv = struct.unpack("<i", buf[off:off + 4])[0]
            off += 4
            intv = struct.unpack("<%dQ" % n_intv, buf[off:off + 8 * n_intv]) if n_intv else ()
            off += 8 * n_intv
            self.refs[names[k].decode()] = (bins, intv)

    def chunks(self, chrom, beg, end):
        """Merged (begin, end) virtual-offset chunks for the region, 0-based [beg, end)."""
        if chrom not in self.refs:
            return []
        bins, intv = self.refs[chrom]
        cand = []
        for b in _reg2bins(beg, end):
            cand.extend(bins.get(b, ()))
        if not cand:
            return []
        lo = intv[min(beg >> 14, len(intv) - 1)] if intv else 0
        cand = sorted((cb, ce) for cb, ce in cand if ce > lo)
        out = []
        for cb, ce in cand:
            cb = max(cb, lo)
            if out and cb <= out[-1][1]:
                out[-1] = (out[-1][0], max(out[-1][1], ce))
            else:
                out.append((cb, ce))
        return out


class Vcf(object):
    """A bgzipped VCF read by region, with its genotypes. Our own reader."""

    LONGEST = 20000  # how far before a region an unindexed read starts (long deletions)

    def __init__(self, path):
        self.path = path
        self.long_ref = 0      # the longest REF an unindexed read met (0: none longer)
        self.samples = []
        self.contigs = set()
        self.tbi = None
        self.first_chrom = None
        self.blockidx = None  # {chrom: ([first POS of a block], [block offset])}
        if os.path.exists(path + ".tbi"):
            try:
                self.tbi = Tabix(path + ".tbi")
            except Exception:  # noqa: BLE001 - an unreadable index: read without one
                self.tbi = None
        bg = Bgzf(path)
        self.bgzf = bg.bgzf
        try:
            for line in bg.lines(0):
                if line.startswith("##"):
                    if line.startswith("##contig=<ID="):
                        self.contigs.add(line[13:].split(",")[0].split(">")[0])
                    continue
                if line.startswith("#CHROM"):
                    self.samples = line.rstrip("\n").split("\t")[9:]
                    continue
                if line:
                    self.first_chrom = line.split("\t", 1)[0]
                break
        finally:
            bg.close()

    def build_blockidx(self):
        """One pass over the file, keeping the CHROM and POS of the first whole line of
        each block: an index in memory (nothing is written next to a shared file)."""
        if self.blockidx is not None or self.tbi is not None or not self.bgzf:
            return
        idx = {}
        bg = Bgzf(self.path)
        try:
            for at, block in bg.blocks():
                nl = block.find(b"\n")
                while nl >= 0 and nl + 1 < len(block) and block[nl + 1:nl + 2] == b"#":
                    nl = block.find(b"\n", nl + 1)
                if nl < 0 or nl + 1 >= len(block):
                    continue
                end = block.find(b"\t", nl + 1)
                end2 = block.find(b"\t", end + 1) if end >= 0 else -1
                if end < 0 or end2 < 0:
                    continue
                chrom = block[nl + 1:end].decode()
                try:
                    pos = int(block[end + 1:end2])
                except ValueError:
                    continue
                ent = idx.setdefault(chrom, ([], []))
                ent[0].append(pos)
                ent[1].append(at)
                # a block that starts on one chromosome and holds the next: the next one's
                # first line is only found by scanning the block, which the reader does
                last = block.rfind(b"\n", 0, len(block) - 1)
                if last > nl:
                    e = block.find(b"\t", last + 1)
                    if e > 0:
                        lc = block[last + 1:e].decode()
                        if lc != chrom:
                            idx.setdefault(lc, ([], []))
                            idx[lc][0].append(0)
                            idx[lc][1].append(at)
        finally:
            bg.close()
        self.blockidx = idx

    def contig(self, chrom):
        """What this VCF calls ``chrom`` (chr22 or 22)."""
        names = [chrom, chrom[3:] if chrom.startswith("chr") else "chr" + chrom]
        known = (set(self.tbi.refs) if self.tbi is not None
                 else set(self.blockidx) if self.blockidx is not None
                 else self.contigs | ({self.first_chrom} if self.first_chrom else set()))
        for n in names:
            if n in known:
                return n
        return chrom

    def chroms(self):
        """The chromosomes this VCF holds (as it names them)."""
        if self.tbi is not None:
            return list(self.tbi.refs)
        if self.bgzf:
            self.build_blockidx()
            return list(self.blockidx)
        return [self.first_chrom] if self.first_chrom else []

    def covers(self, chrom):
        """True / False when this VCF can tell; None for a plain-gzip file that has to be
        read to know."""
        name = self.contig(chrom)
        if self.tbi is not None:
            return name in self.tbi.refs
        if self.bgzf:
            self.build_blockidx()
            return name in self.blockidx
        return None if self.first_chrom != name else True

    def holds(self, chrom):
        """Whether this VCF holds a record of ``chrom`` (a plain-gzip file that cannot
        tell from its first line is read until it finds one)."""
        known = self.covers(chrom)
        if known is not None:
            return known
        name = self.contig(chrom)
        with gzip.open(self.path, "rt") as fh:
            return any(line.startswith(name + "\t") for line in fh)

    def records(self, chrom, beg, end):
        """Every record whose REF span overlaps [beg, end) (0-based), as a dict."""
        chrom = self.contig(chrom)
        bg = Bgzf(self.path)
        try:
            if self.tbi is not None:
                # each chunk is read from its first line to its end (a virtual offset, the
                # end exclusive), and no line twice: the chunks come sorted, and a line an
                # earlier chunk's read reached is not yielded again
                done = -1
                for cb, ce in self.tbi.chunks(chrom, beg, end):
                    for v, line in bg.lines_at(cb):
                        if v >= ce:
                            break
                        if v <= done:
                            continue
                        done = v
                        r = self._parse(line, chrom, beg, end)
                        if r == "past":
                            break
                        if r:
                            yield r
                return
            start = 0
            if self.bgzf:
                self.build_blockidx()
                ent = self.blockidx.get(chrom)
                if not ent:
                    return
                i = bisect.bisect_right(ent[0], beg - self.LONGEST) - 1
                start = ent[1][max(i, 0)] << 16
            seen = False
            for line in bg.lines(start):
                if not line or line[0] == "#":
                    continue
                if not line.startswith(chrom + "\t"):
                    if seen:
                        break  # past the chromosome (a VCF is sorted)
                    continue
                seen = True
                r = self._parse(line, chrom, beg, end)
                if r == "past":
                    break
                if r:
                    if len(r["ref"]) > self.LONGEST:
                        # such a deletion is not read where a region starts more than
                        # LONGEST after it: noted (coverage-unindexed-long-ref)
                        self.long_ref = max(self.long_ref, len(r["ref"]))
                    yield r
        finally:
            bg.close()

    def _parse(self, line, chrom, beg, end):
        """A record as a dict; the genotype columns stay unsplit (only the lines a window
        keeps are split). A sites-only line (eight columns, no FORMAT) is a record with no
        genotypes. Its FILTER is kept as written (Window reads it: PASS_FILTERS)."""
        if not line or line.startswith("#"):
            return None
        head, sep, tail = line.rstrip("\n").partition("\t")
        if head != chrom or not sep:
            return None
        f = tail.split("\t", 8)
        if len(f) < 7:
            return None
        pos1 = int(f[0])
        if pos1 - 1 >= end:
            return "past"
        ref = f[2].upper()
        if pos1 - 1 + len(ref) <= beg:
            return None
        return {"chrom": head, "pos1": pos1, "id": f[1], "ref": ref,
                "alts": [a.upper() for a in f[3].split(",")], "filter": f[5],
                "info": f[6], "format": f[7] if len(f) > 7 else "",
                "tail": f[8] if len(f) > 8 else ""}


# --------------------------------------------------------------------------- #
# reference FASTA, read without writing an .fai
# --------------------------------------------------------------------------- #
class Reference(object):
    """A per-chromosome FASTA folder or one FASTA of any number of contigs. Only the
    named contig is read (chr22 or 22), and nothing is written beside it."""

    def __init__(self, path):
        self.path = path
        self.cache = {}

    def _names(self, chrom):
        return [chrom, chrom[3:] if chrom.startswith("chr") else "chr" + chrom]

    def _file(self, chrom):
        if os.path.isdir(self.path):
            for n in self._names(chrom):
                for name in ("%s.fa" % n, "%s.fasta" % n, "%s.fa.gz" % n, "%s.fasta.gz" % n):
                    p = os.path.join(self.path, name)
                    if os.path.exists(p):
                        return p
            raise FileNotFoundError("no FASTA for %s under %s" % (chrom, self.path))
        return self.path

    def has(self, chrom):
        """Whether the reference holds ``chrom`` (a file of the folder, or a header of
        the one FASTA, read once)."""
        if os.path.isdir(self.path):
            try:
                self._file(chrom)
                return True
            except FileNotFoundError:
                return False
        if not hasattr(self, "_names_in"):
            op = gzip.open if self.path.endswith(".gz") else open
            with op(self.path, "rt") as fh:
                self._names_in = {line[1:].split()[0] for line in fh
                                  if line.startswith(">") and line[1:].split()}
        return bool(set(self._names(chrom)) & self._names_in)

    def seq(self, chrom):
        if chrom in self.cache:
            return self.cache[chrom]
        p = self._file(chrom)
        names = set(self._names(chrom))
        op = gzip.open if p.endswith(".gz") else open
        parts = []
        found = False
        with op(p, "rt") as fh:
            want = False
            for line in fh:
                if line.startswith(">"):
                    name = line[1:].split()[0] if line[1:].split() else ""
                    if want:
                        break  # the contig is read whole
                    want = name in names
                    found = found or want
                    continue
                if want:
                    parts.append(line.strip())
        if not found:
            raise KeyError("no contig %s in %s" % (chrom, p))
        s = "".join(parts).upper()
        if not getattr(self, "keep_all", False):
            self.cache = {}
        self.cache[chrom] = s
        return s


# --------------------------------------------------------------------------- #
# ALT alleles and genotypes, as VCF 4.x defines them
# --------------------------------------------------------------------------- #
class Allele(object):
    """One ALT allele of a VCF line. On a haplotype that carries it, ALT stands in place
    of the REF bases [lo, hi) (0-based). The bases REF and ALT share at either end stay
    where they are, so the allele's edit is what lies between them: it replaces the
    reference bases [edit_lo, edit_hi) with ``body``."""

    __slots__ = ("key", "lo", "hi", "line", "number", "edit_lo", "edit_hi", "body",
                 "kind", "_bases", "searched", "filtered")

    def __init__(self, pos1, ref, alt, line, number):
        self.searched = True    # False: an indel CRISPRme's search skips (FILTER) and only
                                # its phased indel store holds (Window)
        self.filtered = None    # the FILTER of a line outside PASS_FILTERS that is kept
        self.key = (pos1, ref, alt)
        self.lo = pos1 - 1
        self.hi = self.lo + len(ref)
        self.line = line        # the line's id in its Window
        self.number = number    # the allele's number in the GT field (1 = the first ALT)
        head = len(os.path.commonprefix([ref, alt]))
        back = len(os.path.commonprefix([ref[head:][::-1], alt[head:][::-1]]))
        self.edit_lo = self.lo + head
        self.edit_hi = self.hi - back
        self.body = alt[head:len(alt) - back]
        if len(ref) == 1 and len(alt) == 1:
            self.kind = "snv"
        elif len(ref) == len(alt):
            self.kind = "mnv"
        else:
            self.kind = "indel"
        self._bases = None

    def token(self, chrom):
        return "%s_%d_%s_%s" % (chrom, self.key[0], self.key[1], self.key[2])

    def bases(self):
        """[(base, (coordinate, k))] of the edit's ALT bases: the j-th stands at reference
        coordinate edit_lo + j while a replaced base is left to stand in (k = 0); the rest
        are the 1st, 2nd... base added after the last replaced base, or after the shared
        base before the edit when it replaces none (an insertion's anchor)."""
        return list(zip(self.body, self.places()))

    def places(self):
        """The places of the edit's ALT bases (see bases())."""
        if self._bases is None:
            width = self.edit_hi - self.edit_lo
            stand = min(width, len(self.body))
            out = [(self.edit_lo + j, 0) for j in range(stand)]
            anchor = self.edit_lo + width - 1
            out.extend((anchor, k) for k in range(1, len(self.body) - stand + 1))
            self._bases = out
        return self._bases

    def window(self, lo, hi):
        """(first, last) indexes of the edit's ALT bases whose place lies in [lo, hi):
        the bases standing there, and the added ones when their anchor does."""
        width = self.edit_hi - self.edit_lo
        stand = min(width, len(self.body))
        a = max(0, lo - self.edit_lo)
        b = max(a, min(stand, hi - self.edit_lo))
        if len(self.body) > stand and lo <= self.edit_lo + width - 1 < hi:
            b = len(self.body)  # the anchor is in the window, and so are its added bases
        return a, b

    def adds(self):
        return len(self.body) > self.edit_hi - self.edit_lo


def touches(al, a, b):
    """Whether an allele changes what a haplotype holds between reference coordinates a
    and b (inclusive): it replaces or deletes one of them, or adds bases between two."""
    if al.edit_hi > al.edit_lo and al.edit_lo <= b and al.edit_hi - 1 >= a:
        return True
    if al.adds():
        anchor = al.edit_hi - 1 if al.edit_hi > al.edit_lo else al.edit_lo - 1
        return a <= anchor < b
    return False


def changes_target(al, t):
    """Whether an allele changes target ``t``: the target holds a base it wrote (the
    target was read on a haplotype that carries it -- one ending inside an insertion
    anchored at its last coordinate included), or writing it would change the reference
    stretch the target covers (touches())."""
    return any(x is al for x in t.src) or touches(al, t.span[0], t.span[1])


def parse_gt(gt):
    """(allele numbers, phased) of a GT value: None for a missing allele; '|' phased,
    '/' unphased, a single allele (haploid) phased."""
    if "/" in gt:
        parts, phased = gt.replace("|", "/").split("/"), False
    elif "|" in gt:
        parts, phased = gt.split("|"), True
    else:
        parts, phased = [gt], True
    out = []
    for p in parts:
        try:
            out.append(int(p))
        except ValueError:
            out.append(None)
    return tuple(out), phased


_HOMREF = frozenset(("0|0", "0/0", "0", "0|0|0", "0/0/0"))


def merge_call(a, b):
    """One sample's call of one SNV ALT (1, every other called allele 0) that two lines
    give, ``a`` the earlier line's: the one call CRISPRitz 2.10.0's enricher writes into the
    allele's one dictionary entry, which CRISPRme reads (CRISPRitz's CHANGELOG and its
    docs/MULTIALLELIC_ENRICHMENT.md section 2.3). Phased alike -- as many alleles, both
    phased ('|', or haploid) --, a slot holds the ALT where either call puts it there, else
    a called allele (0) where either has one; otherwise, as unphased slots cannot be
    matched up, the call with more copies of the ALT, the earlier one on a tie, as it is."""
    (xa, pa), (xb, pb) = a, b
    if pa and pb and len(xa) == len(xb):
        return tuple(1 if 1 in (x, y) else (0 if (x, y) != (None, None) else None)
                     for x, y in zip(xa, xb)), True
    return b if xb.count(1) > xa.count(1) else a


def merge_calls(into, more):
    """Merge ``more`` ({sample: call} of a later line of one SNV REF>ALT) into ``into``."""
    for s, c in more.items():
        into[s] = merge_call(into[s], c) if s in into else c


def split_call(call, number):
    """A call of a line of several ALTs as a line of one ALT writes it for its allele
    ``number`` (CRISPRitz's enricher's form of a dictionary entry it merges): that allele
    1, every other called allele 0, a missing one missing."""
    alle, phased = call
    return tuple(None if x is None else (1 if x == number else 0) for x in alle), phased


def parse_calls(fmt, tail, names):
    """{sample: (allele numbers, phased)} of every sample whose call is not all reference
    (GT must be the first FORMAT key, as VCF 4.x requires when it is present). A sample
    whose name is None (outside the panel) is left out."""
    out = {}
    if not tail or not fmt.startswith("GT"):
        return out
    fields = tail.split("\t")
    if fmt != "GT":
        fields = [g.partition(":")[0] for g in fields]
    homref = _HOMREF
    for i in [i for i, g in enumerate(fields[:len(names)]) if g not in homref]:
        if names[i] is None:
            continue
        alle, phased = parse_gt(fields[i])
        if any(x != 0 for x in alle):
            out[names[i]] = (alle, phased)
    return out


class Window(object):
    """The VCF lines of one region: their ALT alleles and the samples' calls.

    ``shadowed`` and ``untied`` belong to the documented view only: the SNV alleles whose
    carriers CRISPRme's SNP stores never see (an earlier line of a split multi-allelic
    site: only the last is kept, CHANGELOG 2.8.1 Known limits -- its entry replaced the
    earlier's in CRISPRitz 2.9.1's dictionary, of the same REF>ALT too, whose ALT the
    genome still codes), and the alts of a line of several ALTs whose allele number the
    genotypes cannot tell CRISPRme (METHODS section 4: their carriers are resolved as
    unphased ones; _untie).

    ``panel`` (a set, or None for every sample): the samples of the search's panel; the
    calls of any other sample are not read (CRISPRme lists no carrier outside it).

    ``sem`` (Semantics, or None for the readings of CRISPRme and CRISPRitz before 2.10.0)
    says, by VCF, how the run read the lines: with every allele of a split SNP site coded
    (CRISPRitz 2.10.0) no earlier record is shadowed, and two listings of one SNV REF>ALT at
    a position are one allele, on a line of several ALTs too, their calls merged as
    CRISPRitz 2.10.0's enricher merges them into the allele's one dictionary entry
    (merge_call; the allele then has a line of its own, its calls written as a line of one
    ALT writes them); with the phased indel store of CRISPRme 2.10.0 a genotyped line
    outside PASS_FILTERS keeps no indel either. Two listings of one indel or MNV REF>ALT at a
    position (lines whose FILTER passes) are one allele whatever the release (decision B2),
    merged so too."""

    def __init__(self, refseq, chrom, beg, end, vcfs, needed=None, panel=None, sem=None):
        self.chrom = chrom
        self.sem = sem
        self.refseq = refseq
        self.beg = max(0, beg)
        self.end = end
        self.needed = needed
        self.panel = panel
        self.vcfs = vcfs
        self.samples = []
        self.alleles = []
        self.by_key = defaultdict(list)
        self.line_lo = []       # line id -> 0-based start of its REF
        self.line_hi = []
        self.line_alleles = []  # line id -> [Allele or None per ALT]
        self.calls = []         # line id -> {sample: (allele numbers, phased)}
        self.line_sites_only = []  # line id -> whether its VCF has no sample columns
        self.line_vcf = []      # line id -> the index of its VCF in ``vcfs`` (a dataset)
        self.enum = {}          # (guide, record key) -> its enumeration (record_demands)
        self.cps = {}           # (guide, record key) -> the Cluster Positions of its rows
        self.snpcps = {}        # (guide, piece) -> those of the reference and SNP rows, and
                                # of a per-sample run's capped windows (snp_row_cps)
        self.shadowed = set()
        self.so_shadowed = set()  # a sites-only VCF's earlier records of a split site,
                                  # which its fake-indel windows do not hold (below)
        self.untied = set()
        self.mixed = set()   # 1-based positions whose dictionary lines (every allele coded)
                             # are two or more lines of several ALTs (_dictionary)
        self.snv_at = defaultdict(list)  # reference coordinate -> SNV alleles there
        self.filtered_keys = {}  # record key -> the FILTER of a line of it outside
                                 # PASS_FILTERS (kept or not)
        self.dropped = []        # (0-based start, key, kind, sites-only) of the alleles of
                                 # such lines that no store CRISPRme reads holds
        self._order = None
        seen = set()
        for vi, v in enumerate(vcfs):
            all_alleles, store_filter = (False, False) if sem is None else (
                sem.all_alleles(v, chrom), sem.store_filter(v, chrom))
            names = [s if panel is None or s in panel else None for s in v.samples]
            for s in names:
                if s is not None and s not in seen:
                    seen.add(s)
                    self.samples.append(s)
            snv_lines = defaultdict(list)  # POS -> line ids of a sites-only VCF's SNV-only
                                           # lines, file order
            at = None       # the POS of ``group`` and ``listed``
            group = []      # its dictionary lines (in_dict), file order: (line id or None, r)
            listed = {}     # (POS, REF, ALT) -> the Allele of its first listing there
            for r in v.records(chrom, self.beg, self.end):
                lo = r["pos1"] - 1
                hi = lo + len(r["ref"])
                if not self._wanted(lo, hi):
                    continue
                if r["pos1"] != at:
                    if v.samples:
                        self._dictionary(group, all_alleles, v.samples)
                    at, group, listed = r["pos1"], [], {}
                passing = r["filter"] in PASS_FILTERS
                # a line CRISPRitz's enricher lists in the per-sample dictionary: its REF one
                # base, an ALT of one character at least (a '*' too), its FILTER read
                in_dict = passing and len(r["ref"]) == 1 and any(len(x) == 1
                                                                 for x in r["alts"])
                calls = None
                merged = set()   # the ALT numbers of the line its first listing holds
                if passing and v.samples:
                    # a REF>ALT listed again at its position, on a line of one ALT or of
                    # several, is one allele (decision B2): this line's calls of it go into
                    # its first listing's. An SNV so where every allele is coded (CRISPRitz
                    # 2.10.0's dictionary merges the two entries; the last record's genome
                    # kept the later one: shadowed), an indel or an MNV always (CRISPRme
                    # 2.10.0's readers of the phased indel store merge the two lines, as
                    # CRISPRitz merges an SNV's; one before read the later line alone, a
                    # defect: the carriers only the earlier line listed went missing)
                    for k, alt in enumerate(r["alts"], 1):
                        first = listed.get((r["pos1"], r["ref"], alt))
                        if first is None or (first.kind == "snv" and not (
                                all_alleles and in_dict)):
                            continue
                        if calls is None:
                            calls = parse_calls(r["format"], r["tail"], names)
                        merged.add(k)
                        self._merge_into(first, {s: split_call(c, k)
                                                 for s, c in calls.items() if k in c[0]})
                lid = len(self.line_lo)
                als = []
                for k, alt in enumerate(r["alts"], 1):
                    if (k in merged or not alt or alt in (".", "*") or alt.startswith("<")
                            or set(alt) - set("ACGTN")):
                        als.append(None)
                        continue
                    a = Allele(r["pos1"], r["ref"], alt, lid, k)
                    if not passing:
                        # CRISPRitz's enricher skips the line (PASS_FILTERS). Of a
                        # genotyped dataset the indels stay, on their carriers' haplotypes,
                        # where its phased indel store keeps every indel of the VCF (as up
                        # to 2.8.1; since 2.10.0 it skips them too: store_filter); of a
                        # sites-only dataset the SNVs, which hg38_mega's SNP rows read
                        # whatever their FILTER (its fake-indel windows do not:
                        # sites_only_cps). The rest is in no store CRISPRme reads: read as
                        # the reference
                        keep = (a.kind == "indel" and not store_filter) if v.samples \
                            else a.kind == "snv"
                        self.filtered_keys.setdefault(a.key, r["filter"])
                        if not keep:
                            self.dropped.append((a.lo, a.key, a.kind, not v.samples))
                            als.append(None)
                            continue
                        a.filtered = r["filter"]
                        a.searched = not v.samples
                    als.append(a)
                    if (a.kind != "snv" and passing and v.samples) or (
                            a.kind == "snv" and in_dict and all_alleles):
                        listed.setdefault(a.key, a)
                kept = any(a is not None for a in als)
                if in_dict:
                    group.append((lid if kept else None, r))
                if not kept:
                    continue
                if calls is None:
                    calls = parse_calls(r["format"], r["tail"], names)
                self.line_lo.append(lo)
                self.line_hi.append(hi)
                self.line_alleles.append(als)
                self.calls.append(calls)
                self.line_sites_only.append(not v.samples)
                self.line_vcf.append(vi)
                for a in als:
                    if a is None:
                        continue
                    self.alleles.append(a)
                    self.by_key[a.key].append(a)
                    if a.kind == "snv":
                        self.snv_at[a.lo].append(a)
                if not v.samples and len(r["ref"]) == 1 and all(
                        a.kind == "snv" for a in als if a is not None):
                    snv_lines[r["pos1"]].append(lid)
            if v.samples:
                self._dictionary(group, all_alleles, v.samples)
                continue
            for lids in snv_lines.values():
                if len(lids) < 2 or all_alleles:
                    continue   # (every allele of the site coded: none shadowed)
                # a sites-only VCF's split sites are not read as shadowed by its SNP rows:
                # the mega index CRISPRme ships holds earlier records of some (its SNP rows
                # name them, as chr22:11337076 G>C before G>T), so which it holds cannot be
                # read off the order of the VCF given. Its fake-indel windows hold what
                # CRISPRitz's overlay codes from the enriched genome, the last record whose
                # FILTER passes (the mega's indel rows name it at all 280 split sites they
                # name on chr22, the last PASS one at chr22:25093750, whose last record is
                # LOW_VQSLOD)
                ok = [lid for lid in lids if all(
                    a is None or a.filtered is None for a in self.line_alleles[lid])]
                if len(ok) > 1:
                    last = {a.key for a in self.line_alleles[ok[-1]] if a is not None}
                    for lid in ok[:-1]:
                        for a in self.line_alleles[lid]:
                            if a is not None and a.key not in last:
                                self.so_shadowed.add(a)

    def _wanted(self, lo, hi):
        if self.needed is None:
            return True
        i = bisect.bisect_right(self.needed, (lo, float("inf")))
        for j in (i - 1, i):
            if 0 <= j < len(self.needed):
                a, b = self.needed[j]
                if lo < b and hi > a:
                    return True
        return False

    def _merge_into(self, a, more):
        """Merge ``more`` ({sample: call of allele ``a``, 1 the allele}, a later line's) into
        ``a``'s calls, as CRISPRitz 2.10.0's enricher merges two listings of one REF>ALT into
        the allele's one dictionary entry (merge_call). An allele of a line of several ALTs
        is first given a line of its own, its calls written as a line of one ALT writes them
        (split_call), as the enricher writes the merged entry; the line keeps its other
        ALTs."""
        lid = a.line
        if len(self.line_alleles[lid]) > 1:
            own = len(self.line_lo)
            self.line_lo.append(self.line_lo[lid])
            self.line_hi.append(self.line_hi[lid])
            self.line_alleles.append([a])
            self.calls.append({s: split_call(c, a.number) for s, c in self.calls[lid].items()
                               if a.number in c[0]})
            self.line_sites_only.append(self.line_sites_only[lid])
            self.line_vcf.append(self.line_vcf[lid])
            self.line_alleles[lid][a.number - 1] = None
            a.line, a.number = own, 1
        merge_calls(self.calls[a.line], more)

    def _dictionary(self, group, all_alleles, samples):
        """Read the dictionary lines ``group`` of one position of a genotyped VCF ([(line id
        or None, record)], file order) as CRISPRitz's per-sample dictionary holds them and
        CRISPRme reads it: on a genome of the last record (CRISPRitz 2.9.1 and earlier) the
        last line's entries replaced the earlier ones, whose SNVs are shadowed -- a line of
        '*' alone too, which holds no allele here: 2.9.1's dictionary entry of the position
        is its '*' then, and its genome keeps the code the earlier lines wrote (a '*' adds
        no base, and one base writes no code), so their ALTs are searched and their
        carriers are in no store, as at any other split site; with every allele coded,
        every line's (two listings of one REF>ALT one entry: Window). Then the alts no
        allele number ties are untied (_untie)."""
        if not all_alleles:
            for lid, _r in group[:-1]:
                for a in (self.line_alleles[lid] if lid is not None else ()):
                    if a is not None and a.kind == "snv":
                        self.shadowed.add(a)
            group = group[-1:]
        elif sum(1 for _lid, r in group if len(r["alts"]) > 1) > 1:
            # two lines of several ALTs: each entry keeps its line's raw allele numbers, and
            # CRISPRme's readers cannot tell the two numberings apart (METHODS section 4: an
            # alt can be tied to the other line's number), which is not modeled (main
            # refuses the run)
            self.mixed.add(group[0][1]["pos1"])
        if any(len(r["alts"]) > 1 for _lid, r in group):
            self._untie(group, samples)

    def _untie(self, group, samples):
        """Mark the alts of a line of several ALTs whose allele number CRISPRme's readers
        cannot tell (documented view), by METHODS section 4's rule, read from its text:
        each alt of one character that the dictionary lists at the position (an indel's
        ALT is in no entry; a '*' is listed; two listings of one REF>ALT are one entry,
        written as a line of one ALT writes it) has as number the one every carrier of it
        holds; where that leaves several, the numbers the position's other alts of a line
        of several ALTs settle are taken off, again and again, until none is left to take
        off. An alt whose carriers hold no number but 0 and 1 (a line of one ALT, or the
        entry merged from two lines) settles nothing: its 1 is its own line's. An alt left
        with no number, or several, is untied: CRISPRme resolves its carriers as unphased
        ones. The calls are those of every sample of the VCF (``samples``), as the
        dictionary lists them, whatever the panel."""
        alts = {}   # (REF, ALT) -> [allele or None, [each carrier's numbers], one-ALT form]
        for lid, r in group:
            calls = parse_calls(r["format"], r["tail"], samples)
            for k, alt in enumerate(r["alts"], 1):
                if len(alt) != 1:
                    continue
                held = [{n for n in alle if n} for alle, _p in calls.values() if k in alle]
                slot = alts.get((r["ref"], alt))
                if slot is None:
                    allele = None if lid is None else self.line_alleles[lid][k - 1]
                    alts[(r["ref"], alt)] = [allele, held, False]
                else:
                    slot[1].extend(held)
                    slot[2] = True   # merged: the entry is written in one-ALT form

        def numbers(held, one_alt_form):
            if not held:
                return None   # no carrier: nothing to tie
            return {1} if one_alt_form else set.intersection(*held)

        state = [(numbers(held, one), not one and any(n > 1 for h in held for n in h))
                 for _a, held, one in alts.values()]
        while True:
            off = {min(c) for c, settles in state if settles and c and len(c) == 1}
            after = [((c - off) if c is not None and len(c) > 1 else c, settles)
                     for c, settles in state]
            if after == state:
                break
            state = after
        for (allele, _h, _o), (c, _s) in zip(alts.values(), state):
            if allele is not None and c is not None and len(c) != 1:
                self.untied.add(allele)

    def lines_in(self, lo, hi):
        """The line ids whose REF span meets [lo, hi), in position order."""
        if self._order is None:
            self._order = sorted(range(len(self.line_lo)), key=lambda i: self.line_lo[i])
            self._starts = [self.line_lo[i] for i in self._order]
            self._widest = max([h - l for l, h in zip(self.line_lo, self.line_hi)] or [1])
        i = bisect.bisect_left(self._starts, lo - self._widest)
        j = bisect.bisect_left(self._starts, hi)
        return [k for k in self._order[i:j] if self.line_hi[k] > lo and self.line_lo[k] < hi]

    def call(self, allele, sample):
        return self.calls[allele.line].get(sample)

    def sites_only(self, key):
        """Whether a sites-only VCF (no sample columns) holds the record ``key``: the
        rows a sites-only dataset gives of it list no carrier (CHANGELOG 2.8.1)."""
        return any(self.line_sites_only[a.line] for a in self.by_key.get(key, ()))

    def genotyped(self, key):
        """Whether a VCF with sample columns holds the record ``key`` on a line
        CRISPRme's search reads (PASS_FILTERS)."""
        return any(not self.line_sites_only[a.line] and a.searched
                   for a in self.by_key.get(key, ()))

    def searched(self, key):
        """Whether CRISPRme's search reads the record ``key`` (a line of it passes its
        FILTER, PASS_FILTERS): an indel only the phased indel store holds is on its
        carriers' haplotypes, and no fake-indel window holds it."""
        return any(a.searched for a in self.by_key.get(key, ()))

    def carriers(self, allele):
        """The samples whose call holds the allele on some slot (none for an allele
        CRISPRme's search does not read: no row of it lists a carrier)."""
        if not allele.searched:
            return []
        return [s for s, (alle, _p) in self.calls[allele.line].items()
                if allele.number in alle]

    def signature(self, sample, lids):
        return tuple((lid, self.calls[lid][sample]) for lid in lids
                     if sample in self.calls[lid])


# --------------------------------------------------------------------------- #
# haplotypes: the reference with a slot's alleles in place
# --------------------------------------------------------------------------- #
class Hap(object):
    """A stretch of one haplotype: its bases, each base's place -- (reference coordinate,
    0) for a base standing there, (coordinate, k) for the k-th base an allele adds after
    it -- the allele that wrote it (None: the reference), and the spans where the alleles
    given overlap (no sequence: ``clash``)."""

    __slots__ = ("seq", "place", "src", "clash")

    def __init__(self, seq, place, src, clash):
        self.seq = seq
        self.place = place
        self.src = src
        self.clash = clash

    def span_at(self, coord):
        """(first, end) indexes of the bases whose place has reference coordinate
        ``coord``: the places run in order along a haplotype, so a bisection finds them."""
        i = bisect.bisect_left(self.place, (coord, 0))
        j = bisect.bisect_left(self.place, (coord + 1, 0), i)
        return i, j

    def at(self, coord):
        """Indexes of the bases whose place has reference coordinate ``coord``."""
        i, j = self.span_at(coord)
        return range(i, j)

    def index(self, place):
        i = bisect.bisect_left(self.place, place)
        return i if i < len(self.place) and self.place[i] == place else None

    def clashes(self, a, b):
        """The overlap zones (lo, hi, alleles) that meet reference coordinates a..b."""
        return [z for z in self.clash if z[0] <= b and z[1] - 1 >= a]


def make_hap(refseq, alleles, lo, hi):
    """[lo, hi) of the haplotype that holds ``alleles``. Alleles whose REF spans meet
    describe one stretch twice: none of them is written and the stretch is in ``clash``.
    Only the window is written: an allele reaching into it from either side puts in the
    bases of its edit that stand (or are added after a base) inside it."""
    hi = min(hi, len(refseq))
    use = sorted((a for a in alleles if a.hi > lo and a.lo < hi), key=lambda a: (a.lo, a.hi))
    runs = []
    for a in use:
        if runs and a.lo < runs[-1][1]:
            runs[-1][0].append(a)
            runs[-1][1] = max(runs[-1][1], a.hi)
        else:
            runs.append([[a], a.hi])
    seq, place, src, clash = [], [], [], []
    p = lo

    def ref_upto(stop):
        stop = min(stop, hi)
        if stop > p:
            seq.append(refseq[p:stop])
            place.extend([(q, 0) for q in range(p, stop)])
            src.extend([None] * (stop - p))
        return max(p, stop)

    for members, end in runs:
        if len(members) > 1:
            clash.append((members[0].lo, end, tuple(members)))
            continue
        a = members[0]
        p = ref_upto(a.edit_lo)
        i, j = a.window(lo, hi)
        if j > i:
            seq.append(a.body[i:j])
            place.extend(a.places()[i:j])
            src.extend([a] * (j - i))
        p = max(p, a.edit_hi)
    ref_upto(hi)
    return Hap("".join(seq), place, src, clash)


def holds_edit(hap, i, n, al):
    """The junction rule on hap[i:i + n]: it holds a base the allele adds or replaces, or,
    for a pure deletion, the bases either side of the deleted ones, side by side."""
    if al.body:
        return any(hap.src[k] is al for k in range(i, i + n))
    before, after = (al.edit_lo - 1, 0), (al.edit_hi, 0)
    for k in range(i, i + n - 1):
        if hap.place[k] == before and hap.place[k + 1] == after:
            return True
    return False


def reads_at(hap, coord, probe, al=None):
    """The indexes where a base of reference coordinate ``coord`` starts ``probe`` (and,
    with ``al``, the stretch holds that allele by the junction rule)."""
    n = len(probe)
    i, j = hap.span_at(coord)
    out = []
    k = hap.seq.find(probe, i, j - 1 + n)
    while 0 <= k < j:
        if al is None or holds_edit(hap, k, n, al):
            out.append(k)
        k = hap.seq.find(probe, k + 1, j - 1 + n)
    return out


def span_of(places):
    cs = [pl[0] for pl in places]
    return min(cs), max(cs)


def meets_ref(a, al):
    """Whether allele ``a``'s REF span shares a base with record ``al``'s: on one slot the
    two would describe those bases twice (no sequence there)."""
    return a.lo < al.hi and a.hi > al.lo


# the overlaps of a record that METHODS section 4 reads the record's slot without (the
# record stays listed): another allele at its own position (a phasing artefact, or the
# other allele of a multi-allelic heterozygote), an allele only on bases it deletes, an
# MNV (in no store CRISPRme reads beside the record)
LEFT_OUT = frozenset(("artefact", "deleted", "mnv"))
# the documented causes an overlap gives (the 'mnv' overlap gives none: it is not read)
OVERLAP_CAUSES = frozenset(("artefact", "deleted", "complex", "undocumented"))


def overlap_kind(a, al, phase_known=True):
    """How METHODS section 4 reads allele ``a``, on record ``al``'s slot, where its REF
    meets the record's: 'artefact' (another allele at the record's own position -- an
    alternative to it, so a slot phased with both is a phasing artefact -- or, its phase
    unknown, the other allele of a multi-allelic heterozygote, one allele per haplotype),
    'deleted' (only on bases the record, a pure deletion, removes: it changes nothing),
    'mnv' (an MNV, in no store CRISPRme reads beside the record: not read), 'complex' (an
    indel at another position, one complex allele split into overlapping records: the
    slot is taken off the target's rows whatever its bases) or 'undocumented' (a SNV at
    another position, which METHODS does not cover: read as taking the slot off). The
    first three leave the record's slot listed, read without them."""
    if a.key[0] == al.key[0]:
        return "artefact"
    if not al.body and not a.adds() and a.lo >= al.edit_lo and a.hi <= al.edit_hi:
        return "deleted"
    if not phase_known:
        return "artefact"
    if a.kind == "snv":
        return "undocumented"
    if a.kind == "indel":
        return "complex"
    return "mnv"


def slot_reading(al, present):
    """(the alleles a slot holding record ``al`` and ``present`` is read with, the
    overlapping alleles left out, the kinds of the overlaps; overlap_kind). Where every
    allele whose REF meets the record's is one METHODS section 4 reads the slot without
    (LEFT_OUT) and no allele kept meets one left out, the slot is read without them: the
    one reading of the slot that holds the record (a row naming the record and one of
    them names two records of one base). Otherwise the record keeps its overlap: the VCF
    gives the slot no sequence there (a complex or undocumented overlap, which METHODS
    reads as taking the slot off, or an allele kept that meets one left out)."""
    over = [a for a in present if meets_ref(a, al)]
    if not over:
        return present, [], frozenset()
    kinds = frozenset(overlap_kind(a, al) for a in over)
    rest = [a for a in present if not meets_ref(a, al)]
    if not kinds <= LEFT_OUT or any(meets_ref(b, a) for a in over for b in rest):
        return present, [], kinds
    return rest, over, kinds


def record_clash(al, use):
    """Whether the record still overlaps an allele of a slot read with ``use``: no
    reading of the slot holds it (slot_reading)."""
    return any(meets_ref(a, al) for a in use if a is not al)


READINGS_CAP = 1024  # readings of a slot's overlap zones tried for one target (or record)


def zone_readings(members):
    """The readings of one overlap zone of a slot (alleles whose REF spans meet: the VCF
    gives the stretch no sequence, so no reading is the truth): every set of its alleles
    no two of which meet, none of them included."""
    out = [()]
    for a in members:
        out.extend([x + (a,) for x in out if not any(meets_ref(a, b) for b in x)])
    return out


def slot_reach(tgt, alleles):
    """The furthest reference coordinate a target's bases, read from its start on a slot
    holding ``alleles``, can reach: its own last one, or as many bases on from its start
    as it holds plus every base those alleles delete there."""
    gone = sum(max(0, (a.edit_hi - a.edit_lo) - len(a.body)) for a in alleles
               if a.edit_hi > tgt.start and a.edit_lo < tgt.hi)
    return max(tgt.span[1], tgt.start + tgt.n - 1 + gone)


def zones_for(hap, tgt, named_keys, reach):
    """The overlap zones of ``hap`` that can change what the slot says of target ``tgt``:
    one of their alleles is a variant the row names, or changes a base the slot reads
    from the target's start up to ``reach`` (touches()). A zone none of whose alleles
    does leaves the target's bases, and so every verdict on it, the same on every reading
    of the zone: it is no reason to leave the slot unjudged."""
    return [z for z in hap.clash
            if any(a.key in named_keys or touches(a, tgt.start, reach) for a in z[2])]


def zone_alleles(hap):
    return [a for z in hap.clash for a in z[2]]


def record_core(place, src, al):
    """The places a reading of a target must hold to be read THROUGH record ``al`` (METHODS
    section 4: the target's bases 'through the indel ... not in a stretch that another of
    its variants lengthened, pushing the indel out of it'): the places of the target's
    bases the record's ALT puts there -- each base it adds or replaces and each base its
    REF and ALT share -- and, where the target goes on past the last of them, the base
    right after it (for a deletion, its junction)."""
    idx = [i for i, (pl, s) in enumerate(zip(place, src))
           if (s is not None and s.key == al.key)
           or (pl[1] == 0 and al.lo <= pl[0] < al.hi)]
    if not idx:
        return frozenset()
    out = {place[i] for i in idx}
    if idx[-1] + 1 < len(place):
        out.add(place[idx[-1] + 1])
    return frozenset(out)


def reading_holds(hap, probe, place, core, unread=()):
    """CRISPRme's documented reading of a slot (METHODS section 4): whether ``hap`` holds
    the bases ``probe`` at the target's own place, read from the place of its first base
    on or up to the place of its last base, through its record (every place of ``core``
    in the stretch read), and the stretch meets no allele of ``unread`` (a call of unknown
    phase: those bases cannot be read)."""
    n = len(probe)
    for k in (hap.index(place[0]), None if hap.index(place[-1]) is None
              else hap.index(place[-1]) - n + 1):
        if k is None or k < 0 or hap.seq[k:k + n] != probe:
            continue
        stretch = hap.place[k:k + n]
        if not core <= set(stretch):
            continue
        sp = span_of(stretch)
        if not any(touches(a, sp[0], sp[1]) for a in unread):
            return True
    return False


# --------------------------------------------------------------------------- #
# CFD, from the published tables CRISPRme ships as data
# --------------------------------------------------------------------------- #
class Cfd(object):
    """The CFD score (Doench et al. 2016) as METHODS section 8 documents it: one factor
    per spacer column at its guide position (a DNA bulge's extra base takes no guide
    position), a bulge's gap factor at its position, the two undefined terminal columns
    extrapolated -- lenient at the PAM-distal 5' tip (0.9, the value CRISPRme's scorer
    uses), abolished at the PAM junction -- times the factor of the PAM's last two bases.
    Defined only for a 20-nt guide with a 3-nt 3' PAM, as CRISPRme computes it."""

    TIP, JUNCTION = 0.9, 0.0

    def __init__(self, folder):
        with open(os.path.join(folder, "cfd_mismatch_scores.json")) as fh:
            self.mm = json.load(fh)
        with open(os.path.join(folder, "cfd_pam_scores.json")) as fh:
            self.pam = json.load(fh)

    @staticmethod
    def applies(guide_len, pam_len, pam_5prime):
        return guide_len == 20 and pam_len == 3 and not pam_5prime

    def factor(self, g, d, pos):
        """One spacer column: guide base (or '-') g against target base (or '-') d."""
        g = g.upper().replace("T", "U")
        d = d.upper()
        if g == "-":
            if pos <= 1:
                return self.TIP
            if pos > 20:
                return self.JUNCTION
            return self.mm.get("r-:d%s,%d" % (comp(d), pos), self.JUNCTION)
        if d == "-":
            if pos <= 1:
                return self.TIP
            if pos > 20:
                return self.JUNCTION
            return self.mm.get("r%s:d-,%d" % (g, pos), 1.0)
        if d not in "ACGT":
            return 0.0
        if g == d.replace("T", "U"):
            return 1.0
        return self.mm.get("r%s:d%s,%d" % (g, comp(d), pos), 1.0)

    @staticmethod
    def positions(crrna, n_spacer):
        """The guide position of each spacer column."""
        out = []
        dna = 0
        for c in range(n_spacer):
            out.append(c + 1 - dna)
            if crrna[c] == "-":
                dna += 1
        return out

    def score(self, crrna, target, pam_len=3):
        n = len(crrna) - pam_len
        pos = self.positions(crrna, n)
        s = 1.0
        for c in range(n):
            s *= self.factor(crrna[c], target[c], pos[c])
        return s * self.pam.get(target[-2:].upper(), 0.0)


# --------------------------------------------------------------------------- #
# our own exhaustive aligner
# --------------------------------------------------------------------------- #
_ALIGN_CACHE = {}


def align_pairs(spacer, target, max_mm, max_dna, max_rna, max_total, pam_len,
                pam_5prime):
    """Every alignment of ``spacer`` (guide + PAM motif, IUPAC) to ``target`` (plain
    bases) within the budget. A gap in the spacer is a DNA bulge, a gap in the target an
    RNA bulge; the PAM holds no gap, the spacer's first and last columns never one.
    Yields (aligned spacer, aligned target, mismatches, dna bulges, rna bulges)."""
    ck = (spacer, target, max_mm, max_dna, max_rna, max_total, pam_len, pam_5prime)
    hit = _ALIGN_CACHE.get(ck)
    if hit is not None:
        return hit
    n, m = len(spacer), len(target)
    target = target.upper()
    matches = [_MATCHES.get(c.upper(), frozenset()) for c in spacer]
    if pam_5prime:
        pam_first, pam_last = 0, pam_len - 1
    else:
        pam_first, pam_last = n - pam_len, n - 1
    out = []
    if not _feasible(matches, target, n, m, max_mm, max_dna, max_rna, max_total,
                     pam_first, pam_last, pam_5prime):
        _ALIGN_CACHE[ck] = out
        return out

    def walk(i, j, mm, dn, rn, sa, ta):
        if mm > max_mm or dn > max_dna or rn > max_rna or mm + dn + rn > max_total:
            return
        if n - i > m - j + max_rna or m - j > n - i + max_dna:
            return
        if i == n:
            if j == m and sa and sa[0] != "-" and sa[-1] != "-":
                out.append(("".join(sa), "".join(ta), mm, dn, rn))
            return
        in_pam = pam_first <= i <= pam_last
        if j < m:
            ok = target[j] in matches[i]
            if ok or not in_pam:
                sa.append(spacer[i])
                ta.append(target[j])
                walk(i + 1, j + 1, mm + (0 if ok else 1), dn, rn, sa, ta)
                sa.pop()
                ta.pop()
        dna_ok = (i <= pam_first) if not pam_5prime else (i > pam_last)
        if dna_ok and dn < max_dna and j < m and len(sa) > 0:
            sa.append("-")
            ta.append(target[j])
            walk(i, j + 1, mm, dn + 1, rn, sa, ta)
            sa.pop()
            ta.pop()
        if not in_pam and rn < max_rna:
            sa.append(spacer[i])
            ta.append("-")
            walk(i + 1, j, mm, dn, rn + 1, sa, ta)
            sa.pop()
            ta.pop()

    walk(0, 0, 0, 0, 0, [], [])
    if len(_ALIGN_CACHE) > 200000:
        _ALIGN_CACHE.clear()
    _ALIGN_CACHE[ck] = out
    return out


def _feasible(matches, target, n, m, max_mm, max_dna, max_rna, max_total, pam_first,
              pam_last, pam_5prime):
    """Whether any alignment of the spacer to ``target`` fits the budget: the fewest
    mismatches for each count of DNA and RNA bulges, by dynamic programming over the
    spacer's positions (a DNA bulge takes a target base, an RNA bulge a spacer base, so the
    target position follows). It only rules out pieces with no in-budget alignment, so the
    exhaustive walk that follows finds exactly what it found without it."""
    inf = max_mm + 1
    # cost[(dn, rn)] at spacer position i: fewest mismatches so far
    cost = {(0, 0): 0}
    for i in range(n + 1):
        # DNA bulges at this spacer position (target bases with no spacer base)
        dna_ok = (i <= pam_first) if not pam_5prime else (i > pam_last)
        if dna_ok and i > 0 and max_dna:
            for dn in range(max_dna):
                for rn in range(max_rna + 1):
                    c = cost.get((dn, rn))
                    j = i + dn - rn
                    if c is None or j >= m:
                        continue
                    if c < cost.get((dn + 1, rn), inf):
                        cost[(dn + 1, rn)] = c
        if i == n:
            break
        nxt = {}
        in_pam = pam_first <= i <= pam_last
        for (dn, rn), c in cost.items():
            j = i + dn - rn
            if j < m:
                ok = target[j] in matches[i]
                if ok or not in_pam:
                    c2 = c + (0 if ok else 1)
                    if c2 <= max_mm and c2 < nxt.get((dn, rn), inf):
                        nxt[(dn, rn)] = c2
            if not in_pam and rn < max_rna and c < nxt.get((dn, rn + 1), inf):
                nxt[(dn, rn + 1)] = c
        cost = nxt
        if not cost:
            return False
    for (dn, rn), c in cost.items():
        if i + dn - rn == m and c <= max_mm and c + dn + rn <= max_total:
            return True
    return False


def pam_sites(motif, seq, pam_len):
    out = []
    for q in range(0, len(seq) - pam_len + 1):
        if all(iupac_match(m, b) for m, b in zip(motif, seq[q:q + pam_len])):
            out.append(q)
    return out


def scan_sequence(spacer, seq, max_mm, max_dna, max_rna, max_total, pam_len, pam_5prime,
                  skipped=None):
    """Every in-budget alignment on ``seq``, both strands, anchored on the PAM. Yields
    (strand, offset of the leftmost aligned base on the + reading, aligned spacer,
    aligned target, mm, dna bulges, rna bulges). A piece holding an N is not aligned (no
    base pairs with N); ``skipped`` (a list) gets one entry per such piece."""
    n = len(spacer)
    motif = spacer[:pam_len] if pam_5prime else spacer[n - pam_len:]
    proto = n - pam_len
    for strand, s in (("+", seq), ("-", rc(seq))):
        for q in pam_sites(motif, s, pam_len):
            for width in range(proto - max_rna, proto + max_dna + 1):
                beg = (q - width) if not pam_5prime else q
                end = (q + pam_len) if not pam_5prime else (q + pam_len + width)
                if beg < 0 or end > len(s):
                    continue
                piece = s[beg:end]
                if "N" in piece:
                    if skipped is not None:
                        skipped.append((strand, beg))
                    continue
                for sa, ta, mm, dn, rn in align_pairs(
                        spacer, piece, max_mm, max_dna, max_rna, max_total, pam_len,
                        pam_5prime):
                    off = beg if strand == "+" else len(s) - end
                    yield strand, off, sa, ta, mm, dn, rn


# --------------------------------------------------------------------------- #
# results tables
# --------------------------------------------------------------------------- #
def open_table(path):
    if path.endswith(".zip"):
        with zipfile.ZipFile(path) as z:
            names = [n for n in z.namelist() if n.endswith(".tsv")]
            if not names:
                return io.StringIO("")
            return io.StringIO(z.read(names[0]).decode())
    return open(path, newline="")


def parse_token(tok):
    p = tok.rsplit("_", 3)
    return p[0], int(p[1]), p[2].upper(), p[3].upper()


def var_text(text):
    """A row's Variant_info with each token's REF and ALT in upper case, as Vcf._parse
    reads every record (so an allele's token): a row writes a record as its VCF does, and
    HPRC writes 1,546 indel records of chr22 in lower case (chr22_36126395_ccctctccctcctttc
    ctTCCTT_c), whose rows completeness looked up by their raw tokens and never found."""
    if text in NO_VARIANT:
        return text
    out = []
    for t in text.split(","):
        p = t.rsplit("_", 3)
        out.append("_".join(p[:2] + [p[2].upper(), p[3].upper()]) if len(p) == 4 else t)
    return ",".join(out)


def is_indel_token(tok):
    p = tok.rsplit("_", 3)
    return len(p) == 4 and (len(p[2]) > 1 or len(p[3]) > 1)


def cluster_of(start, strand, probe_len, full_len, pam_5prime):
    """CRISPRitz's Cluster Position of an alignment: its start, moved to the PAM's fixed
    end -- less the RNA bulges and plus the DNA bulges on the + strand with a 3' PAM (and
    the - strand with a 5' PAM), the start itself otherwise. Every alignment of one PAM
    site shares it, which is what the merge folds on."""
    if (strand == "+") != pam_5prime:
        return start + probe_len - full_len
    return start


def guide_of(spacer_col, pam_len, pam_5prime):
    """The guide of a row, from its Spacer+PAM column (the guide and N's over the PAM)."""
    spc = spacer_col.upper()
    return spc[pam_len:] if pam_5prime else spc[:len(spc) - pam_len]


def row_cfd(cfd_model, crrna, aligned, pam_len):
    """The CFD of a row's alignment, from its own aligned strings (None without tables):
    what the analysis writes in the column the merge sorts on (CFD_score_ALT of the
    results; their CFD_score is the higher of the REF and ALT scores)."""
    if cfd_model is None or not crrna or aligned in NO_VARIANT or len(crrna) != len(
            aligned):
        return None
    return cfd_model.score(crrna, aligned.upper(), pam_len)


def read_rows(path, table, chrom=None, lo=None, hi=None, margin=0, pam_5prime=False,
              loci=None, bad=None, pam_len=3, cfd_model=None, novar=None):
    """Every block (dict) of a results table that names a variant. Rows whose start is
    within [lo - margin, hi + margin) are read; those within [lo, hi) are marked for the
    soundness check ('own'). ``loci`` (a dict) gets, per (guide, chromosome), every
    integrated line's (guide, Cluster Position of each block, the fewest-mm+b block's
    edits under each sorting criterion, the highest-CFD block's CFD -- recomputed from
    its aligned strings when ``cfd_model`` is given, the value the merge sorts on, else
    as printed --, line) for the per-cluster completeness check; ``bad`` (a list) gets
    the lines shorter than the header, ``novar`` (a list) the ALT-origin blocks of
    [lo, hi) that name no variant."""
    csv.field_size_limit(1 << 30)
    out = []
    with open_table(path) as fh:
        r = csv.reader(fh, delimiter="\t")
        try:
            head = next(r)
        except StopIteration:
            return out
        col = {n: i for i, n in enumerate(head)}
        if "Chromosome" not in col:
            return out
        for line_no, row in enumerate(r, start=2):
            if not row:
                continue
            if len(row) < len(head):
                if bad is not None:
                    bad.append((table, line_no, len(row), len(head)))
                continue
            ch = row[col["Chromosome"]]
            if chrom and ch != chrom:
                continue
            full_len = len(row[col["Spacer+PAM"]]) if "Spacer+PAM" in col else 23
            guide = (guide_of(row[col["Spacer+PAM"]], pam_len, pam_5prime)
                     if "Spacer+PAM" in col else "")
            line_blocks = []
            for b in BLOCKS:
                try:
                    start = int(row[col["Start_coordinate_(%s)" % b]])
                except (KeyError, ValueError):
                    continue
                if lo is not None and not (lo - margin <= start < hi + margin):
                    continue
                strand = row[col["Strand_(%s)" % b]]
                alt = row[col["Aligned_protospacer+PAM_ALT_(%s)" % b]]
                ref = row[col["Aligned_protospacer+PAM_REF_(%s)" % b]]
                shown = alt if alt not in ("", "NA", "n") else ref
                width = len(shown.replace("-", ""))
                line_blocks.append((b, cluster_of(start, strand, width, full_len,
                                                  pam_5prime)))
                origin = row[col["REF/ALT_origin_(%s)" % b]]
                var = var_text(row[col["Variant_info_genome_(%s)" % b]])
                cfd_col = col.get("CFD_score_(%s)" % b)
                block = {
                    "table": table, "line": line_no, "block": b, "chrom": ch,
                    "spacer_col": row[col["Spacer+PAM"]], "guide": guide,
                    "cfd": row[cfd_col] if cfd_col is not None else "",
                    "start": start,
                    "own": lo is None or lo <= start < hi,
                    "strand": strand,
                    "crrna": row[col["Aligned_spacer+PAM_(%s)" % b]],
                    "ref": ref,
                    "alt": alt,
                    "pam": row[col["PAM_(%s)" % b]],
                    "mm": row[col["Mismatches_(%s)" % b]],
                    "bulges": row[col["Bulges_(%s)" % b]],
                    "total": row[col["Mismatches+bulges_(%s)" % b]],
                    "btype": row[col["Bulge_type_(%s)" % b]],
                    "origin": origin,
                    "var": var,
                    "samples": row[col["Variant_samples_(%s)" % b]],
                }
                if var in NO_VARIANT:
                    # an ALT row is made by variants and must name them
                    if origin == "alt" and novar is not None and block["own"]:
                        novar.append(block)
                    continue
                out.append(block)
            if loci is not None and line_blocks:
                try:
                    crit = edits(int(row[col["Mismatches_(fewest_mm+b)"]]),
                                 int(row[col["Bulges_(fewest_mm+b)"]]))
                except (KeyError, ValueError):
                    crit = {}
                try:
                    cfd = float(row[col["CFD_score_(highest_CFD)"]])
                except (KeyError, ValueError):
                    cfd = None
                try:
                    # the merge sorts on the CFD of the row's own alignment (its
                    # CFD_score_ALT); CFD_score prints the higher of REF and ALT
                    b = "highest_CFD"
                    alt = row[col["Aligned_protospacer+PAM_ALT_(%s)" % b]]
                    shown = (alt if alt not in NO_VARIANT
                             else row[col["Aligned_protospacer+PAM_REF_(%s)" % b]])
                    calc = row_cfd(cfd_model, row[col["Aligned_spacer+PAM_(%s)" % b]],
                                   shown, pam_len)
                    if cfd is not None and calc is not None:
                        cfd = calc
                except KeyError:
                    pass
                loci.setdefault((guide, ch), []).append(
                    (guide, [c for _b, c in line_blocks], crit, cfd, line_no))
    return out


def result_tables(resdir):
    """(path, label) of the integrated results and the alternative alignments."""
    out = []
    for p in sorted(glob.glob(os.path.join(resdir, "*integrated_results.tsv"))):
        out.append((p, "integrated"))
    if not out:
        for p in sorted(glob.glob(os.path.join(resdir, "*integrated_results.zip"))):
            out.append((p, "integrated"))
    stem = os.path.join(resdir, "*all_results_with_alternative_alignments")
    alt = sorted(glob.glob(stem + ".tsv")) or sorted(glob.glob(stem + ".zip"))
    for p in alt:
        out.append((p, "alternative"))
    return out


def _flag(text, flag):
    """Whether a command line holds ``flag`` as a word of its own."""
    return re.search(r"(^|\s)%s(\s|$)" % re.escape(flag), text) is not None


def run_kind(resdir):
    """(analysis mode, alternative alignments written, merge window, the sorting
    criteria of the fewest block (--sorting-criteria, crisprme.py's default
    'mm+bulges,mm' when the command line gives none), the scoring criteria the command
    line gives (--sorting-criteria-scoring, else None), total-edits cap or None, where
    the mode was read, whether the command line is there and names no alternative-
    alignment flag) of a run, from its own markers.

    The mode is the one the search ran, which complete-search writes to .search_mode
    ('per-sample' or 'population-level', since 2.8.0) whatever chose it: --per-sample,
    --population-level, or its context-aware default (per-sample on a phased genotyped
    index, population-level otherwise) -- and --per-sample itself runs population-level
    where what it searches has no genotypes or no genotype store. Without .search_mode,
    the command line's --population-level (which nothing overrides) or --per-sample;
    else None. --ci-cd-test chooses no mode (it keeps the intermediate files):
    complete-test's searches take the default, population-level on the self-test's
    index. Whether the run wrote alternative alignments: .alt_mode, else the command
    line's --alt-alignments or --no-alt-alignments, else the mode's default (written per
    sample, not population-level) where the command line is there to say no flag
    overrode it; else None. A run without them keeps one row per merge cluster. The cap
    is the one the search used (.Params.txt's Max_total_edits, else --max-total-edits in
    .command_line.txt)."""
    mode = alt = cap = scoring = None
    source = None
    merge, criterion = 3, "mm+bulges,mm"
    p = os.path.join(resdir, ".search_mode")
    if os.path.exists(p):
        with open(p) as fh:
            v = fh.read().strip()
        if v in ("per-sample", "population-level"):
            mode, source = v, ".search_mode"
    p = os.path.join(resdir, ".alt_mode")
    if os.path.exists(p):
        with open(p) as fh:
            v = fh.read().strip()
        alt = {"no-alt-alignments": False, "alt-alignments": True}.get(v)
    p = os.path.join(resdir, ".command_line.txt")
    cmdline = os.path.exists(p)
    alt_flag = None
    text = ""
    if cmdline:
        with open(p) as fh:
            text = fh.read()
        if mode is None:
            if _flag(text, "--population-level"):
                mode, source = "population-level", "command line"
            elif _flag(text, "--per-sample"):
                mode, source = "per-sample", "command line"
        if _flag(text, "--no-alt-alignments"):
            alt_flag = False
        elif _flag(text, "--alt-alignments"):
            alt_flag = True
        m = re.search(r"--merge\s+(\d+)", text)
        if m:
            merge = int(m.group(1))
        m = re.search(r"--sorting-criteria\s+(\S+)", text)
        if m:
            criterion = m.group(1)
        m = re.search(r"--sorting-criteria-scoring\s+(\S+)", text)
        if m:
            scoring = m.group(1)
        m = re.search(r"--max-total-edits\s+(-?\d+)", text)
        if m:
            cap = int(m.group(1))
    for name in (".Params.txt", "Params.txt"):
        p = os.path.join(resdir, name)
        if os.path.exists(p):
            with open(p) as fh:
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    if len(f) == 2 and f[0] == "Max_total_edits":
                        try:
                            cap = int(f[1])
                        except ValueError:
                            pass
            break
    if alt is not None and alt_flag is not None and alt != alt_flag:
        raise ValueError(".alt_mode says %s, while the command line says %s"
                         % ("alt-alignments" if alt else "no-alt-alignments",
                            "--alt-alignments" if alt_flag else "--no-alt-alignments"))
    if source == ".search_mode" and mode == "per-sample" and \
            _flag(text, "--population-level"):
        # --population-level forces the population-level analysis: nothing makes it
        # per-sample (the other way, --per-sample on what has no genotypes, is real)
        raise ValueError(".search_mode says per-sample, while the command line says "
                         "--population-level")
    if alt is None:
        alt = alt_flag
    # the command line is there and gives no --alt-alignments / --no-alt-alignments: the
    # mode's default holds (a mode given with --analysis, where the run records none)
    flagless = cmdline and alt_flag is None
    if alt is None and flagless and mode is not None:
        alt = mode == "per-sample"
    return mode, alt, merge, criterion, scoring, cap, source, flagless


def per_sample_downgraded(resdir):
    """Whether the run asked for --per-sample (its .command_line.txt) and ran
    population-level (its .search_mode). complete-search writes that pair in two cases
    only (crisprme.py): every dataset it searches is sites-only, with no genotypes to
    resolve; or a genotyped one has no genotype store (Dictionaries/genotypes_<vcf>/) and
    no classic per-sample dictionaries -- an install made with download --what index
    --no-genotypes --, whose SNP+indel rows list every carrier of the indel. With a
    genotyped VCF of the search given, it is the second."""
    p = os.path.join(resdir, ".search_mode")
    q = os.path.join(resdir, ".command_line.txt")
    if not (os.path.isfile(p) and os.path.isfile(q)):
        return False
    with open(p) as fh:
        mode = fh.read().strip()
    with open(q) as fh:
        text = fh.read()
    return mode == "population-level" and _flag(text, "--per-sample")


REGISTRY_ONLY = "Registry-only install for "


def registry_only_logged(resdir):
    """The chromosomes for which the run's SNP post-analysis logged, in its
    log_verbose.txt (the pipeline's stdout, kept in the result folder), that a dataset it
    read had a Tier-0 registry and neither a genotype store nor classic per-sample
    dictionaries (new_simple_analysis.py: 'Registry-only install for <chrom>: ... no
    genotype tier to resolve carriers.') -- an install made with download --what index
    --no-genotypes --, leaving out an aggregate registry's line (the sites-only
    hg38_mega's, which says '(aggregate info_af registry; ...' on it): a sorted list, or
    None when the folder holds no log_verbose.txt."""
    p = os.path.join(resdir, "log_verbose.txt")
    if not os.path.isfile(p):
        return None
    out = set()
    with open(p, errors="replace") as fh:
        for line in fh:
            if line.startswith(REGISTRY_ONLY) and "aggregate info_af registry" not in line:
                out.add(line[len(REGISTRY_ONLY):].split(":", 1)[0].strip())
    return sorted(out)


INDEXES_SEARCHED = ".indexes_searched.txt"
# the search's own record that the SNP+indel co-occurrence was off (the search script,
# 2.10.0)
OPT_OUT = "NOTE: the SNP+indel co-occurrence is off (CRISPRME_INDEL_SNP=0)"
# the indel post-analysis' record of a phased indel store written before 2.10.0 (one line
# per dataset and chromosome), which holds the indels of every FILTER
OLD_STORE = re.compile(r"NOTE \[indel-snp\]: the phased indel genotype store of (\S+) on "
                       r"(\S+) was written before CRISPRme 2\.10\.0")
# the search script's record of an index of format 1 it took (CRISPRitz 2.10.0 or later)
FORMAT1_NOTE = re.compile(r"^NOTE: (\S+) is in index format 1 \(built by CRISPRitz 2\.9\.1 "
                          r"or earlier\)", re.M)
CZI_MAGIC = b"CZI\xfe"   # the first bytes of a .bin of index format 2 (CRISPRitz 2.10.0)
FLAG_FOUR_ALLELE = 4     # its header flag: a target of the file holds Q


class Semantics(object):
    """How the CRISPRme and CRISPRitz that made a run read its data, where their releases
    differ, read from the run's own records (run_semantics), by dataset and chromosome:

      all_alleles  the variant genome codes every allele of a split SNP site (CRISPRitz
                   2.10.0: no earlier record is read as the reference, two listings of one
                   REF>ALT at a position are one allele, on a line of several ALTs too,
                   their calls merged as its enricher merges them,
                   and a slot holding two ALTs of one position is resolved PUTATIVE,
                   METHODS section 4), or only the last record's (CRISPRitz 2.9.1 and
                   earlier, every index published until the rebuild for 2.10.0, whose
                   dictionary kept the last record's entry of a position, of the same
                   REF>ALT too: limit-last-alt-wins);
      store_filter the phased indel genotype store holds only the indels of records whose
                   FILTER CRISPRitz's enricher reads (CRISPRme 2.10.0), or every indel of
                   the VCF (a store written before);
      indels       the INDELS index's format and pad (its budget, bMax): an index of format
                   1 does not hold a target whose site's guide region ends within bMax
                   bases of its fake-indel window's edge (the separator N in its pad),
                   which one of format 2 holds (CRISPRitz 2.10.0; window_cut);
      capped       whether a per-sample run may report a window over CRISPRME_IUPAC_CAP by
                   its greedy worst-case row (the classic dictionaries up to 2.8.1; 2.10.0
                   enumerates the carriers' haplotypes there).

    ``values`` maps (dataset or None, chromosome or None) to a dict of those, None for one
    the records do not tell; ``dataset_of`` maps a VCF path to its dataset."""

    def __init__(self, release2100, values, dataset_of):
        self.release2100 = release2100
        self.values = values
        self.dataset_of = dataset_of

    def _get(self, vcf, chrom, what):
        ds = self.dataset_of.get(getattr(vcf, "path", vcf))
        for key in ((ds, chrom), (ds, None), (None, chrom), (None, None)):
            v = self.values.get(key)
            if v is not None and what in v:
                return v[what]
        return None

    def all_alleles(self, vcf, chrom):
        return bool(self._get(vcf, chrom, "all_alleles"))

    def store_filter(self, vcf, chrom):
        return bool(self._get(vcf, chrom, "store_filter"))

    def indels(self, vcf, chrom):
        """(format, pad) of the INDELS index of the VCF's dataset on ``chrom``."""
        return self._get(vcf, chrom, "indels") or (None, None)

    @property
    def capped(self):
        return not self.release2100


def _read_log(resdir):
    p = os.path.join(resdir, "log_verbose.txt")
    if not os.path.isfile(p):
        return None
    with open(p, errors="replace") as fh:
        return fh.read()


def indexes_searched(resdir):
    """{dataset: (variant index directory or None, INDELS index name)} of the run's
    .indexes_searched.txt (CRISPRme 2.10.0: 'Variant_index:<dataset>\t<name> (in <dir>):
    built for ...; its indels: <name>', or 'none: a search with no bulges ...'), None when
    the run has none (a CRISPRme before 2.10.0)."""
    p = os.path.join(resdir, INDEXES_SEARCHED)
    if not os.path.isfile(p):
        return None
    out = {}
    with open(p, errors="replace") as fh:
        for line in fh:
            if not line.startswith("Variant_index:"):
                continue
            head, _tab, rest = line.rstrip("\n").partition("\t")
            ds = head[len("Variant_index:"):]
            indels = rest.rpartition("; its indels: ")[2].strip() or None
            m = re.match(r"(\S+) \(in (.*?)\): built for ", rest)
            out[ds] = (os.path.join(m.group(2), m.group(1)) if m else None, indels)
    return out


def bin_header(path):
    """(index format, four-allele flag) of a CRISPRitz .bin, read from its first bytes
    (CRISPRitz's docs/INDEX_FORMAT.md: format 2 opens with the magic 'CZI'+0xFE, then the
    version, the header length and the flags; format 1 with its leaf count), None when it
    cannot be read."""
    try:
        with open(path, "rb") as fh:
            head = fh.read(16)
    except OSError:
        return None
    if head[:4] == CZI_MAGIC and len(head) >= 16:
        version, _hlen, flags = struct.unpack("=3i", head[4:16])
        return version, bool(flags & FLAG_FOUR_ALLELE)
    return (1, False) if len(head) >= 8 else None


def _bins_of(index_dir, chrom):
    """The .bin files of ``index_dir`` (or of its sub-indexes, an INDELS index's) for
    ``chrom``: <PAM>_<chrom>_<n>.bin, or the sub-index of <chrom>'s fake-indel genome."""
    out = []
    if not index_dir or not os.path.isdir(index_dir):
        return out
    for name in sorted(os.listdir(index_dir)):
        path = os.path.join(index_dir, name)
        if os.path.isdir(path) and name.endswith("fake" + chrom):
            out.extend(os.path.join(path, x) for x in sorted(os.listdir(path))
                       if x.endswith(".bin") and not x.startswith("."))
        elif name.endswith(".bin") and not name.startswith(".") and \
                "_%s_" % chrom in name:
            out.append(path)
    return out


def _budget_of(name):
    """The budget (bMax) in an index's name (<guide>-<PAM>_<N>_<genome>, <PAM>_<N>_...)."""
    m = re.search(r"_(\d+)_", os.path.basename(name or ""))
    return int(m.group(1)) if m else None


def _enriched_reads_all(path, sites):
    """Whether the variant genome at ``path`` (a FASTA of one sequence) codes every allele
    of the split SNP sites ``sites`` (split_sites): True where a site's code holds an
    earlier line's ALT beside the last line's, False where every site's holds the last
    line's alone, None where no site tells them apart."""
    try:
        with open(path, "rb") as fh:
            header = fh.readline()
            first = fh.readline()
            width = len(first.rstrip(b"\r\n"))
            eol = len(first) - width
            if width <= 0:
                return None
            verdict = None
            for q, ref, alts, last in sites:
                fh.seek(len(header) + (q // width) * (width + eol) + q % width)
                code = fh.read(1).decode("ascii", "replace").upper()
                got = set(IUPAC.get(code, code)) if code != "Q" else set("ACGT")
                earlier = alts - last
                if not earlier:
                    continue
                if earlier & got:
                    return True
                if last <= got:
                    verdict = False
            return verdict
    except OSError:
        return None


def split_sites(vcf, chrom, most=12, telling=True):
    """Up to ``most`` split SNP sites of ``vcf`` on ``chrom``: [(0-based position, REF,
    {ALT}, {the last line's ALTs})], two or more lines at one position whose FILTER passes,
    each with a REF of one base and an ALT of one of A, C, G and T (a line of several ALTs
    too): the positions a genome of the last record (CRISPRitz 2.9.1 and earlier, whose
    dictionary kept the last line's entries there, and its code the last line's ALTs) and
    one of every allele read otherwise. ``telling``: only those whose earlier lines hold an
    ALT the last line does not, which the enriched genome's code tells apart. A '*' (a
    spanning deletion's) is no such ALT: no IUPAC code of the enriched genome holds it and
    CRISPRitz codes nothing for it, so it tells nothing of how the genome was coded. Not
    ``telling`` (every position the two genomes read otherwise), a line whose ALTs of one
    character are no base (a '*' alone) counts as a line too, where it is the last: 2.9.1's
    dictionary entry of the position is its '*' then, and the earlier lines' carriers are
    lost, while the code they wrote stays (Window._dictionary)."""
    out = []
    pos = ref = None
    group = []   # the ALT sets of the position's lines, file order

    def flush():
        alts = set().union(*group)
        if telling:
            if len(group) > 1 and alts - group[-1]:
                out.append((pos - 1, ref, alts, group[-1]))
        elif any(group[:-1]):
            out.append((pos - 1, ref, alts, group[-1]))

    for r in vcf.records(chrom, 0, 1 << 31):
        alts = {a.upper() for a in r["alts"]} & set("ACGT")
        if r["filter"] not in PASS_FILTERS or len(r["ref"]) != 1:
            continue
        if not alts and (telling or not any(len(a) == 1 for a in r["alts"])):
            continue
        if r["pos1"] != pos:
            if group:
                flush()
                if len(out) >= most:
                    return out
            pos, ref, group = r["pos1"], r["ref"], []
        group.append(alts)
    if group:
        flush()
    return out


def run_semantics(resdir, vcfs, chroms, datasets, genome):
    """(Semantics, why not) of a run, from its own records -- or (None, why) for one whose
    opt-out of the SNP+indel co-occurrence they record (refused). A run of CRISPRme 2.10.0
    or later writes .indexes_searched.txt (and .carrier_gate.json); a run without either is
    an earlier CRISPRme's, which read every index with CRISPRitz 2.9.1 or earlier (format
    1, the last record of a split site) and a phased indel store of every FILTER. Of a run
    of 2.10.0: the INDELS index's format is its .bin files' (the index its
    .indexes_searched.txt names), else its log's format-1 NOTE; the variant genome codes
    every allele where its variant index's .bin files flag Q, or where its enriched genome
    (Genomes/<genome>+<dataset>/<chrom>.enriched.fa in the run's working directory) holds
    an earlier record's ALT at a split site, and the last only where an index of format 1
    served it, or where the enriched genome holds the last ALT alone there, or where its
    index of format 2 flags no Q while its VCF has four bases at a position; the phased
    indel store applies the FILTER rule but where the run's log notes an older store
    (none where the log is gone). What none of them tells is None: main refuses the run
    where it matters."""
    log = _read_log(resdir)
    if log is not None and OPT_OUT in log:
        return None, "opt-out"
    release2100 = (os.path.isfile(os.path.join(resdir, INDEXES_SEARCHED))
                   or os.path.isfile(os.path.join(resdir, ".carrier_gate.json")))
    workdir = os.path.dirname(os.path.dirname(os.path.abspath(resdir)))
    dataset_of = {}
    for v in vcfs:
        parent = os.path.basename(os.path.dirname(os.path.abspath(v.path)))
        if datasets and parent in datasets:
            dataset_of[v.path] = parent
        elif datasets and len(datasets) == 1:
            dataset_of[v.path] = datasets[0]
    values = {}
    if not release2100:
        values[(None, None)] = {"all_alleles": False, "store_filter": False,
                                "indels": (1, None)}
        return Semantics(False, values, dataset_of), None
    old = set()
    format1 = set()
    if log is not None:
        old = {(m.group(1), m.group(2)) for m in OLD_STORE.finditer(log)}
        format1 = {os.path.basename(m.group(1).rstrip("/"))
                   for m in FORMAT1_NOTE.finditer(log)}
    searched = indexes_searched(resdir) or {}
    by_ds = defaultdict(list)
    for v in vcfs:
        by_ds[dataset_of.get(v.path)].append(v)
    libs = [os.path.join(workdir, "genome_library")]
    p = os.path.join(resdir, ".command_line.txt")
    if os.path.isfile(p):
        with open(p, errors="replace") as fh:
            m = re.search(r"--index-path\s+(\S+)", fh.read())
        if m:
            libs.append(m.group(1) if os.path.isabs(m.group(1))
                        else os.path.join(workdir, m.group(1)))
    for ds, vs in by_ds.items():
        var_dir, indels_name = searched.get(ds, (None, None)) if ds else (None, None)
        indels_dir = None
        if indels_name:
            where = [os.path.dirname(var_dir)] if var_dir else []
            for lib in where + libs:
                if os.path.isdir(os.path.join(lib, indels_name)):
                    indels_dir = os.path.join(lib, indels_name)
                    break
        for ch in chroms:
            val = {}
            val["store_filter"] = None if log is None else (ds, ch) not in old
            fmt = None
            heads = [bin_header(b) for b in _bins_of(indels_dir, ch)]
            heads = [h for h in heads if h is not None]
            if heads:
                fmt = 2 if all(h[0] == 2 for h in heads) else 1
            elif indels_name and indels_name in format1:
                fmt = 1
            val["indels"] = (fmt, _budget_of(indels_name))
            alls = None
            var_heads = [h for h in (bin_header(b) for b in _bins_of(var_dir, ch)) if h]
            if var_heads and all(h[0] == 1 for h in var_heads):
                alls = False
            elif var_heads and any(h[1] for h in var_heads):
                alls = True
            if alls is None and ds:
                fa = os.path.join(workdir, "Genomes", "%s+%s" % (genome, ds),
                                  ch + ".enriched.fa")
                if os.path.isfile(fa):
                    sites = []
                    for v in vs:
                        if v.covers(ch) is not False:
                            sites.extend(split_sites(v, ch))
                    alls = _enriched_reads_all(fa, sites)
            if alls is None and var_heads and all(h[0] == 2 for h in var_heads) and \
                    any(four_base_position(v, ch) for v in vs):
                alls = False   # format 2 with no Q while the VCF has four bases somewhere
            val["all_alleles"] = alls
            values[(ds, ch)] = val
    return Semantics(True, values, dataset_of), None


def four_base_position(vcf, chrom):
    """Whether ``vcf`` has a position on ``chrom`` whose reference base and the SNV ALTs of
    its lines that pass their FILTER hold A, C, G and T (CRISPRitz 2.10.0 codes it Q)."""
    pos, held = None, set()
    for r in vcf.records(chrom, 0, 1 << 31):
        if r["filter"] not in PASS_FILTERS or len(r["ref"]) != 1:
            continue
        if r["pos1"] != pos:
            pos, held = r["pos1"], {r["ref"].upper()}
        held |= {a.upper() for a in r["alts"] if len(a) == 1}
        if held >= set("ACGT"):
            return True
    return False


def run_params(resdir):
    """The run's .Params.txt as {key: value} (empty when it has none)."""
    out = {}
    for name in (".Params.txt", "Params.txt"):
        p = os.path.join(resdir, name)
        if os.path.exists(p):
            with open(p) as fh:
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    if len(f) == 2:
                        out[f[0]] = f[1]
            break
    return out


def run_guides(resdir, pam_motif, pam_5prime):
    """The guides the run searched (its .guides.txt), or None when it does not say."""
    p = os.path.join(resdir, ".guides.txt")
    if not os.path.isfile(p):
        return None
    try:
        with open(p) as fh:
            return set(read_guides([x.strip() for x in fh if x.strip()], len(pam_motif),
                                   pam_5prime))
    except (OSError, ValueError):
        return None


def run_mismatch(resdir, args, pam_motif, pam_5prime, guides):
    """What the check was given that the run's own record of its search contradicts (a
    list of texts; empty when nothing does): the budget (--mm, --bDNA, --bRNA,
    --max-total-edits against .Params.txt's Mismatches, DNA, RNA, Max_total_edits), the
    PAM (its Pam), a guide the run did not search (its .guides.txt; some of its guides
    are a check of those guides: rows are keyed by guide) and the samples list (--samples
    against the run's .sampleID.txt, which CRISPRme writes from the search's samplesID
    files). A check of other inputs than the search's models the run partially."""
    out = []
    params = run_params(resdir)
    for key, flag, given in (("Mismatches", "--mm", args.mm), ("DNA", "--bDNA", args.bdna),
                             ("RNA", "--bRNA", args.brna),
                             ("Max_total_edits", "--max-total-edits", args.max_total_edits)):
        try:
            v = int(params[key])
        except (KeyError, ValueError):
            continue
        if given is not None and v != given:
            out.append("%s %d, while the run searched with %s %d" % (flag, given, key, v))
    if params.get("Pam") and params["Pam"].upper() != pam_motif.upper():
        out.append("the PAM file's motif is %s, while the run searched %s"
                   % (pam_motif, params["Pam"]))
    ran = run_guides(resdir, pam_motif, pam_5prime)
    if ran is not None and set(guides) - ran:
        out.append("--guide gives %s, which the run did not search (it searched %s)"
                   % (",".join(sorted(set(guides) - ran)), ",".join(sorted(ran))))
    p = os.path.join(resdir, ".sampleID.txt")
    if args.samples and os.path.isfile(p):
        try:
            given, ran = read_panel(args.samples), read_panel([p])
        except OSError:
            given = ran = None
        if given is not None and given != ran:
            out.append("--samples gives %d sample(s), %d of them not in the run's "
                       ".sampleID.txt, which lists %d not given"
                       % (len(given), len(given - ran), len(ran - given)))
    return out


def run_datasets(resdir):
    """(datasets, where read): the variant datasets the run searched, or (None, why not)
    when the run does not say. A search with bulges records them in its .Params.txt
    (Genome_idx: one '<PAM>_<bulges>_<genome>+<dataset>' item per dataset,
    comma-separated; an item with no '+' is the reference alone). A search without bulges
    reads no index and records 'Genome_idx None' whatever it searched (crisprme.py): its
    datasets are then the lines of the --vcf list its command line names (a dataset
    folder per line), read where complete-search read it -- its working directory, which
    holds Results/<output> (an absolute path as given); none when the command line gives
    no --vcf (the reference alone)."""
    for name in (".Params.txt", "Params.txt"):
        p = os.path.join(resdir, name)
        if not os.path.exists(p):
            continue
        with open(p) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) == 2 and f[0] == "Genome_idx":
                    if f[1].strip() == "None":
                        return _datasets_of_command_line(resdir)
                    return ([x.split("+", 1)[1] for x in f[1].split(",") if "+" in x],
                            "Genome_idx in its .Params.txt")
        return None, "its .Params.txt has no Genome_idx"
    return None, "it has no .Params.txt"


def _datasets_of_command_line(resdir):
    """run_datasets for a run whose .Params.txt says 'Genome_idx None' (no bulges)."""
    p = os.path.join(resdir, ".command_line.txt")
    if not os.path.isfile(p):
        return None, "Genome_idx None (a search without bulges) and no .command_line.txt"
    with open(p) as fh:
        text = fh.read()
    m = re.search(r"(?:^|\s)--vcf\s+(\S+)", text)
    if not m:
        return [], "its command line gives no --vcf"
    given = m.group(1)
    path = given if os.path.isabs(given) else os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(resdir))), given)
    if not os.path.isfile(path):
        return None, ("Genome_idx None (a search without bulges) and its --vcf list %s is "
                      "not at %s" % (given, path))
    names = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line:
                names.append(os.path.basename(line.rstrip("/")))
    return names, "the --vcf list %s its command line names (Genome_idx None: a search " \
                  "without bulges)" % path


def samples_of(row):
    return [s for s in row["samples"].split(",") if s and s not in ("NA", "n", ".")]


# --------------------------------------------------------------------------- #
# findings
# --------------------------------------------------------------------------- #
COLS = ["gating", "class", "table", "line", "block", "chrom", "start", "strand",
        "spacer", "target", "variants", "samples", "detail"]


class Report(object):
    def __init__(self):
        self.rows = []

    def add(self, gating, cls, row, detail, samples=None):
        self.rows.append({
            "gating": "yes" if gating else "no", "class": cls,
            "table": row.get("table", ""), "line": row.get("line", ""),
            "block": row.get("block", ""), "chrom": row.get("chrom", ""),
            "start": row.get("start", ""), "strand": row.get("strand", ""),
            "spacer": row.get("crrna", ""), "target": row.get("alt", ""),
            "variants": row.get("var", ""),
            "samples": samples if samples is not None else row.get("samples", ""),
            "detail": detail,
        })

    def gating(self):
        return sum(1 for r in self.rows if r["gating"] == "yes")

    def counts(self):
        out = defaultdict(int)
        for r in self.rows:
            out[r["class"]] += 1
        return out

    def write(self, path):
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n",
                               extrasaction="ignore")
            w.writeheader()
            for r in self.rows:
                w.writerow(r)


def recount(crrna, target, pam_len, pam_5prime):
    """(mismatches, dna bulges, rna bulges, problem) of an aligned pair as CRISPRme
    writes it: the spacer's gaps are DNA bulges, the target's RNA bulges, case ignored,
    no gap in the PAM, and the spacer's guide end away from the PAM never a gap."""
    if len(crrna) != len(target):
        return 0, 0, 0, "aligned spacer and target differ in length"
    n = len(crrna)
    pam_cols = set(range(0, pam_len)) if pam_5prime else set(range(n - pam_len, n))
    far = n - 1 if pam_5prime else 0
    if crrna[far] == "-":
        return 0, 0, 0, "a gap at the spacer's guide end"
    mm = dn = rn = 0
    for i, (g, d) in enumerate(zip(crrna, target)):
        if g == "-" and d == "-":
            return 0, 0, 0, "a column with a gap in both strings"
        if i in pam_cols and (g == "-" or d == "-"):
            return 0, 0, 0, "a gap in the PAM"
        if g == "-":
            dn += 1
            continue
        if d == "-":
            rn += 1
            continue
        if i in pam_cols:
            continue
        if not iupac_match(g, d):
            mm += 1
    return mm, dn, rn, None


def base_columns(alt_aligned, strand):
    """For each base of the target in + order, its column in the aligned strings."""
    cols = [j for j, ch in enumerate(alt_aligned) if ch != "-"]
    return cols if strand == "+" else cols[::-1]


class Ctx(object):
    """What every check needs: the budget, the PAM, the guides, the CFD tables."""

    def __init__(self, args, pam_motif, pam_len, pam_5prime, guides, cfd=None):
        self.args = args
        self.mm = args.mm
        self.bdna = args.bdna
        self.brna = args.brna
        self.cap = args.max_total_edits if args.max_total_edits is not None else 99
        self.pam_motif = pam_motif
        self.pam_len = pam_len
        self.pam_5prime = pam_5prime
        self.guides = guides  # [guide without PAM]
        # whether the search scores CFD (CRISPRme's indel path does only for a 20-nt guide
        # and a 3-nt 3' PAM: elsewhere its PUTATIVE rule takes the reference base, else
        # the first of A, C, G, T); the tables missing here leave such a choice open
        self.cfd_applies = bool(guides) and Cfd.applies(len(guides[0]), pam_len,
                                                        pam_5prime)
        self.cfd = cfd if cfd is not None and self.cfd_applies else None
        self.panel = None     # the search's samples (a set), None for every VCF sample
        # the SNP rows a merge cluster can hold (snp_row_sets): a population-level run's are
        # the worst case of each alignment (every allele an option), a per-sample run's the
        # haplotypes individuals hold
        self.snp_union = True
        # how the run's CRISPRme and CRISPRitz read its data (run_semantics), or None for
        # the readings before 2.10.0
        self.sem = None
        # the run's --max-degeneracy (its .Params.txt's Max_degeneracy, 2.10.0), or None
        self.max_degeneracy = None

    def full(self, guide):
        return (self.pam_motif + guide) if self.pam_5prime else (guide + self.pam_motif)

    def pam_cols(self, n):
        return range(0, self.pam_len) if self.pam_5prime else range(n - self.pam_len, n)


# --------------------------------------------------------------------------- #
# a row's own strings
# --------------------------------------------------------------------------- #
def alignment_checks(row, ctx, rep, gating):
    """The checks of a row's own strings: its guide, its aligned spacer, the recount, the
    bulge type, the budget, the PAM. False when the row cannot be read further."""
    crrna, alt = row["crrna"], row["alt"]
    stripped = alt.replace("-", "").upper()
    if set(stripped) - set("ACGT"):
        rep.add(gating, "iupac-left", row, "the target holds %s"
                % "".join(sorted(set(stripped) - set("ACGT"))))
        return False
    pam_len, p5 = ctx.pam_len, ctx.pam_5prime
    spc = row["spacer_col"].upper()
    sguide = spc[pam_len:] if p5 else spc[:len(spc) - pam_len]
    guide = next((g for g in ctx.guides if g == sguide), None)
    if guide is None:
        rep.add(gating, "guide-unknown", row,
                "the row's Spacer+PAM %s is none of the given guides" % row["spacer_col"])
        return False
    nogap = crrna.replace("-", "").upper()
    g_part = nogap[pam_len:] if p5 else nogap[:len(nogap) - pam_len]
    p_part = nogap[:pam_len] if p5 else nogap[len(nogap) - pam_len:]
    if g_part != guide or not (p_part == "N" * pam_len or p_part == ctx.pam_motif):
        rep.add(gating, "spacer-not-guide", row,
                "the aligned spacer, gaps removed, is %s, not the guide %s with its PAM"
                % (nogap, guide))
        return False
    mm, dn, rn, problem = recount(crrna, alt, pam_len, p5)
    if problem:
        rep.add(gating, "alignment-invalid", row, problem)
        return False
    try:
        r_mm, r_bul, r_tot = int(row["mm"]), int(row["bulges"]), int(row["total"])
    except ValueError:
        rep.add(gating, "alignment-invalid", row, "unreadable counts")
        return False
    if (mm, dn + rn) != (r_mm, r_bul) or r_tot != r_mm + r_bul:
        rep.add(gating, "counts-mismatch", row,
                "reported %d mm + %d bulges = %d; the alignment has %d mm, %d DNA and "
                "%d RNA bulges" % (r_mm, r_bul, r_tot, mm, dn, rn))
    want = ",".join([t for t, c in (("RNA", rn), ("DNA", dn)) if c]) or "X"
    if row["btype"] != want and not (want == "X" and row["btype"] in ("", "NA")):
        rep.add(gating, "bulge-type-mismatch", row,
                "reported %s, the alignment has %d DNA and %d RNA bulges"
                % (row["btype"], dn, rn))
    if mm > ctx.mm or dn > ctx.bdna or rn > ctx.brna or mm + dn + rn > ctx.cap:
        rep.add(gating, "budget-exceeded", row,
                "%d mm, %d DNA, %d RNA bulges against %d/%d/%d, cap %s"
                % (mm, dn, rn, ctx.mm, ctx.bdna, ctx.brna, ctx.args.max_total_edits))
    pam_seq = stripped[:pam_len] if p5 else stripped[-pam_len:]
    if not all(iupac_match(m, b) for m, b in zip(ctx.pam_motif, pam_seq)):
        rep.add(gating, "pam-mismatch", row, "PAM %s against the motif %s"
                % (pam_seq, ctx.pam_motif))
    if row["pam"] and row["pam"] != "NA" and row["pam"].upper() != pam_seq:
        rep.add(gating, "pam-column-mismatch", row,
                "the PAM column says %s, the target has %s" % (row["pam"], pam_seq))
    return True


def ref_text_check(row, win, rep, gating):
    stripped = row["alt"].replace("-", "").upper()
    span = win.refseq[row["start"]:row["start"] + len(stripped)]
    if len(span) != len(stripped):
        rep.add(False, "limit-contig-end", row,
                "the REF text cannot be read: the target runs past the contig's end")
        return
    text = list(span if row["strand"] == "+" else rc(span))
    for i, c in enumerate(row["alt"]):
        if c == "-":
            text.insert(i, "-")
    if row["ref"].upper() != "".join(text).upper():
        rep.add(gating, "ref-text-mismatch", row,
                "REF %s, the reference at %d gapped as the target is %s"
                % (row["ref"], row["start"], "".join(text)))


def named_alleles(row, win, rep, gating):
    """The row's variants as Allele objects (the first line holding each), or None after
    a finding."""
    named = []
    for t in [x for x in row["var"].split(",") if x and x != "NA"]:
        try:
            c, p, r, a = parse_token(t)
        except (ValueError, IndexError):
            rep.add(gating, "variant-unparseable", row, t)
            return None
        found = win.by_key.get((p, r, a)) if c == row["chrom"] else None
        if found and not any(x.searched for x in found):
            rep.add(gating, "variant-filtered", row,
                    "%s has FILTER %s, which CRISPRitz's enricher skips (only PASS and '.' "
                    "are read): no fake-indel window, indel log or SNP store holds it (the "
                    "phased indel store keeps it, for its carriers' other targets), so no "
                    "row names it" % (t, found[0].filtered))
            return None
        if not found:
            if c == row["chrom"] and (p, r, a) in win.filtered_keys:
                rep.add(gating, "variant-filtered", row,
                        "%s has FILTER %s, which CRISPRitz's enricher skips (only PASS and "
                        "'.' are read): no store CRISPRme's search reads holds it"
                        % (t, win.filtered_keys[(p, r, a)]))
            else:
                rep.add(gating, "variant-not-in-vcf", row,
                        "%s is not a record of the VCF (forward strand)" % t)
            return None
        if any(x in win.shadowed for x in found) and not any(
                x not in win.shadowed for x in found):
            # the record is real, so the row is not wrong for naming it (its target is
            # checked against it); what it contradicts is a documented limit: a note
            rep.add(False, "variant-shadowed", row,
                    "%s is an earlier SNP record of a split site, which the documented "
                    "limit says CRISPRme's stores do not hold (CHANGELOG 2.8.1 Known "
                    "limits): this index holds it" % t)
        # (of two lines of one REF>ALT, the one CRISPRme's stores read)
        named.append(next((x for x in found if x.searched and x not in win.shadowed),
                          next((x for x in found if x.searched), found[0])))
    return named


# --------------------------------------------------------------------------- #
# a row's target: the reference with its named variants
# --------------------------------------------------------------------------- #
class Target(object):
    """A target holding ``main`` (the row's indel-path record): its + reading ``probe``,
    the places of its bases (on the sequence that makes it), the reference coordinates it
    spans, its SNV code columns (base index, coordinate), the window the checks write
    haplotypes over and the row's aligned strings (None in completeness when unknown).

    The code columns are those of CRISPRme's fake-indel window, the only places its
    PUTATIVE rule chooses a base at: a reference base of the target, off the record's own
    REF span, where an SNV of the VCF lies. The window writes the record's ALT verbatim and
    its overlay codes flank bases only (METHODS section 4), so a base on the record's REF
    span -- an indel's anchor, the bases an unnormalized record keeps -- is never a choice,
    whatever SNV the VCF holds there (a multi-allelic heterozygote's other allele, one per
    haplotype)."""

    __slots__ = ("main", "named", "probe", "strand", "start", "place", "span", "codes",
                 "lo", "hi", "crrna", "aligned", "n", "src")

    def __init__(self, win, main, named, probe, strand, place, src, crrna, aligned):
        self.main = main
        self.named = named
        self.probe = probe
        self.strand = strand
        self.place = place
        self.src = src
        self.n = len(probe)
        self.start = place[0][0]
        self.span = span_of(place)
        self.codes = [(i, pl[0]) for i, (pl, s) in enumerate(zip(place, src))
                      if pl[1] == 0 and s is not main and win.snv_at.get(pl[0])
                      and not (main is not None and main.lo <= pl[0] < main.hi)]
        self.lo = max(0, min(self.span[0], main.lo) - 2 * self.n - 60)
        self.hi = min(len(win.refseq), max(self.span[1], main.hi) + 2 * self.n + 60)
        self.crrna = crrna
        self.aligned = aligned


def row_target(win, row, named, main):
    """The row's target on the reference with exactly its named variants: (Target,
    problem class, detail)."""
    stripped = row["alt"].replace("-", "").upper()
    probe = stripped if row["strand"] == "+" else rc(stripped)
    n = len(probe)
    # a row of the record starts within one target of its REF span (inserted bases are
    # written after their anchor whatever the window's size)
    width = (main.hi - main.lo) if main else 1
    anchor = main.lo if main else row["start"]
    lo = max(0, min(row["start"], anchor) - n - width - 60)
    hi = min(len(win.refseq), max(row["start"], anchor + width) + 2 * n + width + 60)
    hap = make_hap(win.refseq, named, lo, hi)
    if hap.clash:
        return None, "variants-overlap", (
            "the named variants cannot all be on one haplotype: two of them describe the "
            "reference bases %s" % ",".join("%d..%d" % (z[0], z[1] - 1) for z in hap.clash))
    here = reads_at(hap, row["start"], probe)
    for i in here:
        if main is None or holds_edit(hap, i, n, main):
            return Target(win, main, named, probe, row["strand"], hap.place[i:i + n],
                          hap.src[i:i + n], row["crrna"], row["alt"]), None, None
    through = []
    k = hap.seq.find(probe)
    while k >= 0:
        if main is None or holds_edit(hap, k, n, main):
            through.append(hap.place[k][0])
        k = hap.seq.find(probe, k + 1)
    if through:
        return None, "start-not-on-variant-sequence", (
            "the target is the reference with its named variants%s starting at %s, not at "
            "%d" % (", through its record," if main else "",
                    ",".join(map(str, sorted(set(through))[:6])), row["start"]))
    if here:
        return None, "target-holds-no-indel-base", (
            "the target at %d holds no base its record %s changes (a flank target)"
            % (row["start"], main.token(win.chrom)))
    return None, "target-not-on-variant-sequence", (
        "the target is not the reference with its named variants near %d: a base no "
        "named variant explains" % row["start"])


# --------------------------------------------------------------------------- #
# the truth: the sample's haplotypes, as its calls give them
# --------------------------------------------------------------------------- #
def main_call(win, key, sample):
    """(the line's allele of ``key``, the sample's call) where the sample carries it on a
    line CRISPRme's search reads (a row's record: the dataset that skips the line lists
    none of its carriers)."""
    for a in win.by_key.get(key, ()):
        if not a.searched:
            continue
        c = win.call(a, sample)
        if c is not None and a.number in c[0]:
            return a, c
    return None, None


def record_slots(call, number):
    """(the slots holding allele ``number``, whether the call says which slots they
    are). An unphased heterozygous call holds it on one slot of unknown place (None)."""
    alle, phased = call
    hs = [h for h, x in enumerate(alle) if x == number]
    if not hs:
        return [], True
    if phased or len(alle) == 1 or all(x == number for x in alle):
        return hs, True
    return [None], False


def on_slot(call_main, h, call, number):
    """Whether slot h of the record's call also holds allele ``number`` of another line
    ('yes' / 'no' / '?'): known when both calls are phased alike, or the other call is
    homozygous or lacks the allele; unknown otherwise (an unphased call, a missing allele,
    a haploid call against a diploid one, or the record's own slot unknown)."""
    if call is None:
        return "no"
    alle, phased = call
    if number not in alle and None not in alle:
        return "no"
    ploidy = len(call_main[0])
    if len(alle) == ploidy and all(x == number for x in alle):
        return "yes"
    if h is None or not phased or not call_main[1] or len(alle) != ploidy:
        return "?"
    x = alle[h]
    if x is None:
        return "?"
    return "yes" if x == number else "no"


def one_per_key(win, alleles):
    """``alleles`` with each REF>ALT once: of two lines of one (a split site read as its
    last record, on a genome of CRISPRitz 2.9.1), the one CRISPRme's stores read (not
    shadowed). The VCF says the same allele twice there, not two alleles of one base."""
    best = {}
    for a in alleles:
        b = best.get(a.key)
        if b is None or (b in win.shadowed and a not in win.shadowed):
            best[a.key] = a
    return [a for a in alleles if best[a.key] is a]


def slot_alleles(win, sample, main_al, call_main, h, lids):
    """(alleles the slot holds, alleles it may hold) among the lines ``lids``, each REF>ALT
    once (one_per_key). Two or more SNV ALTs of one position that the slot holds (split
    records phased one by one), each read by CRISPRme's stores (none shadowed: every allele
    coded) and none on the record's REF span, are alleles it may hold: the slot has no one
    sequence in the VCF there, and CRISPRme resolves such a carrier PUTATIVE (METHODS section
    4, decision B4), the rule it resolves a call of unknown phase by -- which so is the
    truth there."""
    present, unknown = [], []
    for lid in lids:
        call = win.calls[lid].get(sample)
        if call is None:
            continue
        for a in win.line_alleles[lid]:
            if a is None or a is main_al:
                continue
            if lid == main_al.line:
                continue  # one allele of a line per slot: the slot's is the record's
            r = on_slot(call_main, h, call, a.number)
            if r == "yes":
                present.append(a)
            elif r == "?":
                unknown.append(a)
    present = one_per_key(win, present)
    held = {a.key for a in present}
    unknown = [a for a in one_per_key(win, unknown) if a.key not in held]
    at = defaultdict(list)
    for a in present:
        if a.kind == "snv" and a not in win.shadowed and not meets_ref(a, main_al):
            at[a.lo].append(a)
    two = {id(a) for alts in at.values() if len(alts) > 1 for a in alts}
    if two:
        unknown = unknown + [a for a in present if id(a) in two]
        present = [a for a in present if id(a) not in two]
    return present, unknown


def truth(win, tgt, sample, lids):
    """What the sample's own haplotypes say of the row's target:
      ('none',)                       it carries no allele of the row's record
      ('exact', holds, info)          its alleles where the target reads are on known slots
      ('resolved', info)              only a slot read without alleles that overlap the
                                      record (slot_reading) holds it
      ('clash', info)                 a slot of the record has no one sequence there
      ('unknown', info)               which alleles share the record's slot there is unknown
    ``holds``: some slot of the record reads the target from the row's start, through the
    record -- the record's bases the target holds at their places, as the row's target has
    them (record_core) -- and holds every named variant. A slot holding an allele whose REF
    meets the record's has no sequence there: where METHODS section 4 reads it without that
    allele (a phasing artefact at the record's position, an allele only on bases it
    deletes, an MNV; slot_reading) the slot is read so -- the one reading of it that holds
    the record, and so the one a row naming the record and not that allele can mean -- and
    otherwise it is a clash, whatever its calls of unknown phase. Two other alleles of the
    slot whose REFs meet (a zone) leave it no sequence only where they can change what it
    says of the target (zones_for: one of them is named, or changes a base read from the
    target's start): every reading of such zones -- any of their alleles no two of which
    meet, or none -- is tried, and where all of them say the same, with no call of
    unknown phase in the way, the slot is judged so; where they part it is a clash. A slot
    that holds a named variant on no reading (neither on it nor of unknown phase there)
    lacks the target whatever its overlaps. A call of unknown phase whose REF meets the
    record's would leave the record's slot with no sequence: the reading under which the
    slot has one (a multi-allelic heterozygote, one allele per haplotype, as the
    documented rule reads it) puts it on the other slot, so it changes nothing here."""
    main_al, call = main_call(win, tgt.main.key, sample)
    if main_al is None:
        return ("none",)
    slots, _known = record_slots(call, main_al.number)
    named = [a for a in tgt.named if a.key != tgt.main.key]
    named_keys = {a.key for a in named}
    core = record_core(tgt.place, tgt.src, tgt.main)
    info = []
    verdicts = []
    for h in slots:
        present, unknown = slot_alleles(win, sample, main_al, call, h, lids)
        use, left, _kinds = slot_reading(main_al, present)
        keys = {a.key for a in use}
        can_named = all(k in keys or any(u.key == k for u in unknown) for k in named_keys)
        hap = make_hap(win.refseq, [main_al] + use, tgt.lo, tgt.hi)
        extra = [a for a in use if touches(a, tgt.span[0], tgt.span[1])
                 and a.key not in named_keys]
        if not can_named:
            # a named variant the slot holds on no reading: it never holds the target
            info.append({"slot": h, "holds": False, "named": False, "extra": extra,
                         "unknown": [], "clash": [], "left": left})
            verdicts.append((False, False, False, False))
            continue
        if record_clash(main_al, use):
            # the record overlaps an allele no documented reading leaves out: no sequence,
            # and no call of unknown phase gives it one
            info.append({"slot": h, "holds": False, "named": True, "extra": extra,
                         "unknown": [], "clash": hap.clashes(main_al.lo, main_al.hi - 1),
                         "left": left})
            verdicts.append((False, False, True, False))
            continue
        # two other alleles of the slot that overlap each other (a zone: no sequence there)
        # leave the slot unjudged only where they can change what it says of the target:
        # every reading of such zones (any of their alleles no two of which meet, or none)
        # is tried, and the slot is judged where all of them agree
        zones = zones_for(hap, tgt, named_keys, slot_reach(tgt, [main_al] + use))
        readings = [use]
        if zones:
            choices = [zone_readings(z[2]) for z in zones]
            count = 1
            for c in choices:
                count *= len(c)
            if count > READINGS_CAP:
                info.append({"slot": h, "holds": False, "named": True, "extra": extra,
                             "unknown": [], "clash": zones, "left": left})
                verdicts.append((False, False, True, bool(left)))
                continue
            inzone = {id(a) for z in zones for a in z[2]}
            fixed = [a for a in use if id(a) not in inzone]
            readings = [fixed + [a for x in combo for a in x]
                        for combo in itertools.product(*choices)]
        results = []
        for r in readings:
            hr = hap if r is use else make_hap(win.refseq, [main_al] + r, tgt.lo, tgt.hi)
            rkeys = {a.key for a in r}
            hits = [k for k in reads_at(hr, tgt.start, tgt.probe, main_al)
                    if core <= set(hr.place[k:k + tgt.n])]
            span = span_of(hr.place[hits[0]:hits[0] + tgt.n]) if hits else tgt.span
            lo_, hi_ = min(span[0], tgt.span[0]), max(span[1], tgt.span[1])
            unk = [a for a in unknown if touches(a, lo_, hi_) and not meets_ref(a, main_al)]
            has_named = all(k in rkeys for k in named_keys)
            holds = bool(hits) and has_named
            if unk and h is not None:
                # the record's slot is known: what its known alleles say comes first. When
                # no placement of the calls of unknown phase lets the slot hold the target
                # (a phased second indel inside it, say), or every placement does, the
                # slot is judged; only when they decide is the PUTATIVE rule the truth
                can, must = placements(win, tgt, main_al, r, unk, named)
                if can is not None and (not can or must):
                    holds, unk = bool(must), []
            results.append((holds, has_named, unk))
        if len(results) > 1 and (any(u for _h, _n, u in results)
                                 or len({hd for hd, _n, _u in results}) > 1):
            # the readings of the zones part (or meet a call of unknown phase): the VCF
            # gives the slot no one sequence where the target reads
            info.append({"slot": h, "holds": False, "named": True, "extra": extra,
                         "unknown": [], "clash": zones, "left": left})
            verdicts.append((False, False, True, bool(left)))
            continue
        holds, _n, unk = results[0]
        has_named = any(n for _h, n, _u in results)
        info.append({"slot": h, "holds": holds, "named": has_named, "extra": extra,
                     "unknown": unk, "clash": [], "left": left})
        verdicts.append((holds, bool(unk), False, bool(left)))
    if any(hd and not u and not c and not lf for hd, u, c, lf in verdicts):
        return ("exact", True, info)
    if any(hd and not u and not c for hd, u, c, _lf in verdicts):
        return ("resolved", info)
    if any(u for _hd, u, _c, _lf in verdicts):
        return ("unknown", info)
    if any(c for _hd, _u, c, _lf in verdicts):
        return ("clash", info)
    return ("exact", False, info)


def placements(win, tgt, main_al, present, unk, named, cap=PUTATIVE_CAP):
    """(can, must) for a known slot of the record: whether some placement of the calls of
    unknown phase ``unk`` (each of a line's alleles on the slot, or none) lets the slot
    read the target at its start through the record (record_core) with every named
    variant, and whether every placement does. (None, None) when a placement gives
    overlapping alleles that can change what the slot says of the target (zones_for: no
    sequence) or there are too many lines to try. A call whose REF meets the record's is
    left out (truth())."""
    by_line = defaultdict(list)
    for a in unk:
        if not meets_ref(a, main_al):
            by_line[a.line].append(a)
    if len(by_line) > cap:
        return None, None
    want = {a.key for a in named}
    core = record_core(tgt.place, tgt.src, tgt.main)
    can, must = False, True
    for combo in itertools.product(*[[None] + v for v in by_line.values()]):
        chosen = [a for a in combo if a is not None]
        hap = make_hap(win.refseq, [main_al] + present + chosen, tgt.lo, tgt.hi)
        if zones_for(hap, tgt, want, slot_reach(tgt, [main_al] + present + chosen)):
            return None, None
        keys = {a.key for a in present} | {a.key for a in chosen}
        ok = want <= keys and any(core <= set(hap.place[k:k + tgt.n])
                                  for k in reads_at(hap, tgt.start, tgt.probe, main_al))
        can = can or ok
        must = must and ok
    return can, must


def options(win, sample, q, visible_only=False, ploidy=0):
    """The bases the sample's own alleles allow at reference coordinate q: the ALT of every
    SNV allele there it carries, and the reference base unless its ALT copies there are as
    many as its slots (a missing allele may be the reference): its copies counted against
    the record's call, as METHODS section 4 counts them, so two ALTs of one position phased
    onto one slot, which the VCF gives no sequence, leave the reference out as they do in
    CRISPRme. Its slots are as many as the most copies any of its calls there has, or
    ``ploidy`` (the record's call, or the most its calls around have) when that is more:
    METHODS section 4 resolves 'a haploid call against a diploid one' PUTATIVE, so a
    haploid ALT call beside a diploid record leaves the reference possible on the other
    copy. ``visible_only``: as CRISPRme's stores see them (an earlier line of a split site
    reads as the reference). Two lines of one REF>ALT (a split site read as its last record)
    give one allele, its copies those of their calls merged (merge_call), not their sum."""
    alts = set()
    copies = 0
    missing = False
    calls = {}   # REF>ALT -> (allele, the sample's call of it: 1 the ALT, 0 another allele)
    for a in win.snv_at.get(q, ()):
        call = win.call(a, sample)
        if call is None:
            continue
        alle = call[0]
        ploidy = max(ploidy, len(alle))
        if visible_only and a in win.shadowed:
            continue
        mine = (tuple(None if x is None else int(x == a.number) for x in alle), call[1])
        seen = calls.get(a.key)
        calls[a.key] = (a, mine if seen is None else merge_call(seen[1], mine))
    for a, (alle, _phased) in calls.values():
        for x in alle:
            if x == 1:
                alts.add(a.body)
                copies += 1
            elif x is None:
                missing = True
    if missing or copies < max(ploidy, 1):
        alts.add(win.refseq[q])
    return alts


def strongest(ctx, tgt, opts_by_code, refseq=None):
    """The documented PUTATIVE choice (METHODS section 4) over the target's code columns:
    {base index: set of allowed bases (+ orientation)} for the spacer columns, and the set
    of allowed tuples for the PAM's code columns (jointly), plus whether every choice is
    determined. The guide's base wins wherever it is an option and the alignment pairs
    the column with a guide base; elsewhere (a DNA bulge included), and at the PAM, the
    options that maximize the CFD of the whole alignment. CFD is a product of one factor
    per spacer column and the PAM's factor, so where the best product is above 0 each
    column takes its own best options; where it is 0 (a column no choice changes scores
    0, as a DNA bulge whose base the table gives 0, or every option of a column does, or
    every PAM the motif allows) every choice ties at 0 and each is a strongest one (METHODS
    documents no order among CFD ties: the choice is left open). A search that scores no
    CFD (a 5' PAM, or a guide or PAM not 20 and 3 nt long) takes there the reference base
    where the carrier's alleles allow it (at a PAM column, and the motif), else the first
    of A, C, G, T they allow, read in the target's orientation (METHODS section 4): one
    choice. ``refseq``: the reference, for that rule (without it every option ties). A
    search that scores CFD checked without the CFD tables leaves its choices open."""
    if tgt.crrna is None or tgt.aligned is None:
        return None
    if getattr(ctx, "cfd_applies", True):
        refseq = None   # a CFD search: without the tables, its choices are left open
    n_cols = len(tgt.crrna)
    cols = base_columns(tgt.aligned, tgt.strand)
    pam = set(ctx.pam_cols(n_cols))
    n_spacer = n_cols - ctx.pam_len
    positions = Cfd.positions(tgt.crrna, n_spacer)
    ori = (lambda b: b) if tgt.strand == "+" else comp
    cfd = ctx.cfd
    allowed = {}
    scores = {}       # base index -> {option: its column's CFD factor} (a choice to make)
    pam_codes = []
    code_cols = set()
    fixed = 1.0       # the factors no choice changes
    def first(q, opts, ok=lambda b: True):
        # no CFD: the reference base where it is allowed, else the first of A, C, G, T
        # (target orientation); None when no option is allowed
        good = [b for b in opts if ok(b)]
        if refseq is not None and refseq[q].upper() in good:
            return refseq[q].upper()
        return min(good, key=lambda b: "ACGT".find(ori(b))) if good else None

    at_q = dict(tgt.codes)
    for i, _q in tgt.codes:
        opts = opts_by_code[i]
        col = cols[i]
        code_cols.add(col)
        if col in pam:
            pam_codes.append(i)
            continue
        g = tgt.crrna[col].upper()
        o_opts = {ori(b): b for b in opts}
        if g != "-" and g in o_opts:
            allowed[i] = {o_opts[g]}   # the guide's base wins (a factor of 1)
            continue
        if cfd is None:
            allowed[i] = {first(at_q[i], opts)} if refseq is not None else set(opts)
            continue
        scores[i] = {b: cfd.factor(g, ori(b), positions[col]) for b in opts}
    if cfd is not None:
        for col in range(n_spacer):
            if col not in code_cols:
                fixed *= cfd.factor(tgt.crrna[col], tgt.aligned[col], positions[col])
    pam_ok = None
    pam_top = 1.0
    combos = []
    if pam_codes:
        target_bases = list(tgt.aligned.replace("-", "").upper())
        # target_bases is in the target's orientation: base index i (+ order) sits at
        # position i on the + strand, n-1-i on the - strand
        def at(i):
            return i if tgt.strand == "+" else len(target_bases) - 1 - i
        for combo in itertools.product(*[sorted(opts_by_code[i]) for i in pam_codes]):
            tb = list(target_bases)
            for i, b in zip(pam_codes, combo):
                tb[at(i)] = ori(b)
            pseq = "".join(tb[:ctx.pam_len] if ctx.pam_5prime else tb[-ctx.pam_len:])
            if not all(iupac_match(m, x) for m, x in zip(ctx.pam_motif, pseq)):
                continue
            score = cfd.pam.get(pseq[-2:], 0.0) if cfd is not None else 0.0
            combos.append((score, combo))
        pam_top = max((s for s, _c in combos), default=0.0)
        if cfd is None and refseq is not None:
            # no CFD: each PAM column on its own, the reference base where its alleles and
            # the motif allow it, else the first of A, C, G, T they allow
            n_t = len(target_bases)

            def motif_ok(i):
                k = at(i) if ctx.pam_5prime else at(i) - (n_t - ctx.pam_len)
                return lambda b: iupac_match(ctx.pam_motif[k], ori(b))
            pick = tuple(first(at_q[i], opts_by_code[i], motif_ok(i)) for i in pam_codes)
            combos = [(0.0, pick)] if None not in pick and any(
                c == pick for _s, c in combos) else []
    elif cfd is not None:
        tail = tgt.aligned.replace("-", "").upper()
        pam_top = cfd.pam.get((tail[-ctx.pam_len:])[-2:], 0.0)
    best = fixed * pam_top
    for sc in scores.values():
        best *= max(sc.values())
    tie = (cfd is None and refseq is None) or (cfd is not None and best <= 0.0)
    for i, sc in scores.items():
        top = max(sc.values())
        allowed[i] = set(sc) if tie else {b for b, v in sc.items() if v >= top - 1e-12}
    if pam_codes:
        pam_ok = ({c for _s, c in combos} if tie
                  else {c for s, c in combos if s >= pam_top - 1e-12})
    determined = all(len(v) == 1 for v in allowed.values()) and (
        not pam_codes or len(pam_ok) == 1)
    return allowed, pam_codes, pam_ok, determined


def record_ploidy(win, al, sample):
    """How many copies the sample's call of record ``al`` has (0 when it carries none)."""
    c = main_call(win, al.key, sample)[1]
    return len(c[0]) if c is not None else 0


def putative_check(win, ctx, tgt, sample, visible_only=False):
    """(allowed, strongest, determined) of the row's bases for a sample read by the
    PUTATIVE rule: every base at a code column is an allele of the sample; the bases are
    its strongest choice; that choice is one only (no tie the rule leaves open)."""
    ploidy = record_ploidy(win, tgt.main, sample)
    opts = {i: options(win, sample, q, visible_only, ploidy) for i, q in tgt.codes}
    row = {i: tgt.probe[i] for i, _q in tgt.codes}
    allowed = all(row[i] in opts[i] for i in opts)
    st = strongest(ctx, tgt, opts, win.refseq)
    if st is None:
        return allowed, allowed, False
    per, pam_codes, pam_ok, determined = st
    strong = all(row[i] in per[i] for i in per)
    if pam_codes:
        strong = strong and tuple(row[i] for i in pam_codes) in pam_ok
    return allowed, allowed and strong, determined


# --------------------------------------------------------------------------- #
# CRISPRme's documented view of a carrier (classification only)
# --------------------------------------------------------------------------- #
def documented(win, ctx, tgt, sample, lids):
    """What CRISPRme's documented rules make of a carrier of the row's record for the
    row's target (METHODS section 4; CHANGELOG 2.8.1 Known limits), re-implemented from
    that text: (verdict, causes), verdict 'lists', 'off' or 'maybe' (PUTATIVE, and the
    rule leaves the choice open), causes a set naming what decided it."""
    main_al, call = main_call(win, tgt.main.key, sample)
    if main_al is None:
        return "off", {"no-record"}
    causes = set()
    if main_al.kind == "mnv":
        causes.add("putative-mnv")  # MNVs are not in the phased indel store
        return _documented_putative(win, ctx, tgt, sample, causes)
    slots, known = record_slots(call, main_al.number)
    if not known:
        causes.add("putative")
        return _documented_putative(win, ctx, tgt, sample, causes)
    s0, s1 = tgt.span
    # the reading must be through the record: the record's bases the row's target holds,
    # at their places (and the base past them), in the stretch read
    core = record_core(tgt.place, tgt.src, tgt.main)
    verdicts = []     # per slot: (the rules list it, an MNV on it changes the target)
    put_mnv = []      # the same flag for each slot the rules resolve PUTATIVE
    putative = None   # a call of unknown phase (or an untied alt) on a slot the rules keep
    for h in slots:
        snvs, seconds = [], []
        off = None
        unsure = puts = None
        # an MNV is in no store CRISPRme reads for an indel's carrier (the phased indel
        # store keeps records whose REF and ALT lengths differ, the SNP tiers SNVs): the
        # slot is read without it, whatever it changes
        mnv_here = False
        for lid in lids:
            c = win.calls[lid].get(sample)
            if c is None or lid == main_al.line:
                continue
            for a in win.line_alleles[lid]:
                if a is None:
                    continue
                r = on_slot(call, h, c, a.number)
                if r == "no":
                    continue
                if meets_ref(a, main_al):
                    # overlap_kind: an allele at the record's own position (or, its phase
                    # unknown, a multi-allelic heterozygote) and one only on bases it
                    # deletes stay out of the reading; an MNV is not read; an indel at
                    # another position takes the slot off, and so (undocumented) does a SNV
                    ok = overlap_kind(a, main_al, phase_known=(r != "?"))
                    if ok == "complex":
                        off = "complex"
                    elif ok == "undocumented":
                        off = off or "undocumented"
                    elif ok != "mnv":
                        causes.add(ok)
                    continue
                if a.kind == "snv":
                    if not touches(a, s0, s1):
                        continue
                    if r == "?":
                        puts = puts or "putative"
                        continue
                    if a in win.untied:
                        puts = puts or "putative-untied"
                        continue
                    if a in win.shadowed:
                        causes.add("last-alt-wins")
                        continue
                    if any(x.lo == a.lo for x in snvs):
                        # two ALTs of one position on the slot (split records phased one
                        # by one): no haplotype, and CRISPRme resolves the carrier
                        # PUTATIVE (METHODS section 4; every allele coded since 2.10.0)
                        puts = puts or "putative-two-alts"
                        snvs = [x for x in snvs if x.lo != a.lo]
                        continue
                    snvs.append(a)
                elif a.kind == "indel":
                    if r == "?":
                        if changes_target(a, tgt):
                            unsure = "putative"
                        continue
                    seconds.append(a)
                elif changes_target(a, tgt):  # an MNV, not read
                    if r == "yes":
                        mnv_here = True
                    else:
                        causes.add("second-indel")
        if unsure:
            # a second indel inside the target with an unphased call: listed PUTATIVE,
            # even when a phased second indel takes the slot off (METHODS section 4)
            causes.add(unsure)
            return _documented_putative(win, ctx, tgt, sample, causes)
        # the slot's phased variants come first (CRISPRme takes a slot off before it
        # resolves anything PUTATIVE): an overlapping record at another position, and a
        # second indel inside the target whose reading fails
        if off:
            causes.add(off)
            verdicts.append((False, mnv_here))
            continue
        inside = [a for a in seconds if changes_target(a, tgt)]
        # the window's target as CRISPRme resolves it for this slot: the record with the
        # slot's own SNP alleles (METHODS section 4: a row lists the haplotypes left with
        # its sequence and names the SNP alts they hold), read at the row's target's place
        window = make_hap(win.refseq, [main_al] + snvs, tgt.lo, tgt.hi)
        own = any(core <= set(window.place[k:k + tgt.n])
                  for k in reads_at(window, tgt.start, tgt.probe, main_al))
        if not inside:
            if puts:
                putative = putative or puts
                put_mnv.append(mnv_here)
                continue
            verdicts.append((own, mnv_here))
            continue
        causes.add("second-indel")
        if not puts and not own:
            # the slot's own SNP alleles make another sequence of the window (its row, if
            # any): a second indel that happens to give this row's bases (a repeat) does
            # not list it here, for SNP alts it does not hold
            verdicts.append((False, mnv_here))
            continue
        # the slot's bases, with its phased store indels and its own SNP alleles around
        snv_all = []
        unread = []
        for lid in lids:
            c = win.calls[lid].get(sample)
            if c is None or lid == main_al.line:
                continue
            for a in win.line_alleles[lid]:
                if a is None or a.kind != "snv" or a in win.shadowed:
                    continue
                if a.lo < main_al.hi and a.hi > main_al.lo:
                    continue
                r = on_slot(call, h, c, a.number)
                if r == "yes":
                    snv_all.append(a)
                elif r == "?":
                    unread.append(a)
        hap = make_hap(win.refseq, [main_al] + snv_all + seconds, tgt.lo, tgt.hi)
        # a way of reading that meets two of its variants that overlap (two records of one
        # base: no sequence) cannot be read, as one that meets a call of unknown phase
        ok = reading_holds(hap, tgt.probe, tgt.place, core, unread + zone_alleles(hap))
        if ok and puts:
            putative = putative or puts   # read through, then resolved PUTATIVE
            put_mnv.append(mnv_here)
            continue
        verdicts.append((ok, mnv_here))
    # 'mnv-unseen': the rules list the carrier only through slots read without an MNV
    # that changes the target there -- a sequence the carrier does not have, which no
    # documented rule lists (gating); a slot the rules leave off whose MNV changes the
    # target needs a second indel-path record for it (the two-indel limit)
    if putative is not None:
        # a slot the rules keep has a call of unknown phase in the target: PUTATIVE
        causes.add(putative)
        if all(put_mnv) and not any(v for v, _m in verdicts):
            causes.add("mnv-unseen")
        return _documented_putative(win, ctx, tgt, sample, causes)
    listing = [m for v, m in verdicts if v]
    if listing:
        if all(listing):
            causes.add("mnv-unseen")
        return "lists", causes
    if any(m for _v, m in verdicts):
        causes.add("second-indel")
    return "off", causes


def _documented_putative(win, ctx, tgt, sample, causes):
    allowed, strong, determined = putative_check(win, ctx, tgt, sample, visible_only=True)
    if any(a in win.shadowed and win.call(a, sample) is not None
           for _i, q in tgt.codes for a in win.snv_at.get(q, ())):
        causes.add("last-alt-wins")
    if not (allowed and strong):
        return "off", causes
    return ("lists" if determined else "maybe"), causes


LIMIT_OF = (("last-alt-wins", "limit-last-alt-wins"),
            ("second-indel", "limit-two-indels"),
            ("mnv-unseen", "limit-two-indels"),   # only where the carrier is left off
            ("putative-mnv", "limit-putative-carrier"),
            ("putative-untied", "limit-putative-carrier"),
            ("putative-two-alts", "limit-putative-carrier"),
            ("putative", "limit-putative-carrier"),
            ("complex", "limit-overlapping-records"),
            ("undocumented", "limit-overlapping-records"),
            ("artefact", "limit-overlapping-records"),
            ("deleted", "limit-overlapping-records"))


def limit_of(causes):
    for c, cls in LIMIT_OF:
        if c in causes:
            return cls
    return None


# --------------------------------------------------------------------------- #
# soundness
# --------------------------------------------------------------------------- #
CARRIER_TEXT = {
    "carrier-without-indel": "listed without carrying %s",
    "carrier-without-target": "listed, and no haplotype of theirs holds the target at its "
                              "start with %s",
    "carrier-without-variants": "listed, and no slot of theirs holding %s holds every "
                                "named variant",
    "variant-missing": "listed, and their slot of %s holds a variant inside the target "
                       "the row does not name",
    "carrier-without-allele": "listed (phase unknown there) for an allele they do not "
                              "carry, with %s",
    "carrier-not-strongest": "listed (phase unknown there) for a weaker target than the "
                             "strongest their alleles allow, with %s (METHODS section 4: "
                             "the guide's base wherever they carry it, else the highest "
                             "CFD; a search that scores no CFD the reference base, else "
                             "the first of A, C, G, T)",
    "carrier-off-by-rule": "listed, while their slot of %s has no one sequence in the VCF "
                           "where the target reads (alleles that overlap) and CRISPRme's "
                           "documented rules take it off this target (METHODS section 4: "
                           "a complex allele, a second indel the slot cannot be read "
                           "through, the last record of a split site, the PUTATIVE "
                           "rule)",
    "carrier-missing": "hold the target at its start with %s (or it is their strongest), "
                       "are listed by CRISPRme's documented rules, and are not on the row",
    "carrier-not-in-panel": "listed with %s, and not samples of the search's panel "
                            "(--samples, else the run's .sampleID.txt)",
    "unexplained-disagreement": "their haplotypes and CRISPRme's documented rules part on "
                                "this row of %s for no documented reason",
    "limit-last-alt-wins": "differ from the row by an earlier SNP record of a split site, "
                           "which CRISPRme does not read, with %s",
    "limit-two-indels": "need a second indel-path record on the slot of %s, or have one "
                        "read by the documented rule",
    "limit-overlapping-records": "hold alleles that overlap %s or each other where the "
                                 "target reads (no one sequence in the VCF there)",
    "overlap-reading-agrees": "are listed, and their slot of %s, whose alleles overlap the "
                              "record, holds the target as METHODS section 4 reads it "
                              "(without the alleles it reads such a slot without)",
    "limit-putative-carrier": "are carriers of %s that CRISPRme resolves PUTATIVE (an MNV, "
                              "an untied multi-ALT alt, an unphased call) while their "
                              "haplotype says otherwise",
}


MNV_UNSEEN = ("an MNV on their slot changes the target; CRISPRme reads no MNV beside the "
              "record (its phased indel store keeps records whose REF and ALT lengths "
              "differ, its SNP tiers SNVs), and no documented rule lists a carrier for a "
              "sequence it does not have")


def carrier_verdict(win, ctx, tgt, sample, listed, lids):
    """(gating, class[, detail]) for one sample of a row, or None when it agrees."""
    tv = truth(win, tgt, sample, lids)
    kind = tv[0]
    if kind == "none":
        if listed and not win.genotyped(tgt.main.key):
            return (True, "carrier-without-indel",
                    "only sites-only VCF(s) hold the record: nobody's genotype carries it, "
                    "and a sites-only dataset's rows list no carrier")
        return (True, "carrier-without-indel") if listed else None
    if kind == "exact":
        holds, info = tv[1], tv[2]
        if holds == listed:
            return None
        verdict, causes = documented(win, ctx, tgt, sample, lids)
        if holds:  # on its haplotype, not on the row
            if verdict == "lists":
                return True, "carrier-missing"
            if verdict == "maybe":
                # CRISPRme resolves it PUTATIVE and its rule ties this target with
                # another: it may be listed on the other (completeness looks there)
                return False, "limit-putative-carrier", (
                    "documented: %s; a PUTATIVE tie with another target"
                    % ",".join(sorted(causes)))
            return unexplained_or(causes)
        if verdict in ("lists", "maybe"):
            if "mnv-unseen" in causes:
                return True, "variant-missing", MNV_UNSEEN
            return unexplained_or(causes)
        if "complex" in causes:
            return True, "carrier-off-by-rule"
        return True, diagnose(info)
    if kind == "resolved":
        # the record's slot holds an allele that overlaps it, which METHODS section 4
        # reads the slot without (it stays listed); read so, the slot holds the target
        left = sorted({a.token(win.chrom) for x in tv[1] for a in x["left"]})
        if listed:
            # the row agrees with the slot as METHODS reads it: a note of its own, apart
            # from the slots this checker leaves unjudged
            return False, "overlap-reading-agrees", "read without %s" % ", ".join(left)
        return False, "limit-overlapping-records", (
            "not listed, while their slot, read without the alleles that overlap the "
            "record as METHODS section 4 reads it (%s; such a carrier stays listed), holds "
            "the target: not modeled into a failure" % ", ".join(left))
    if kind == "clash":
        # no one sequence of the record's slot where the target reads (a complex or
        # undocumented overlap of the record, or a zone of two other alleles whose
        # readings part), and a slot that can hold every named variant (truth: one that
        # cannot lacks the target whatever its overlaps). The VCF gives no truth there, so
        # CRISPRme's documented rules are the only model: a listed carrier they take off
        # this target for a documented cause is listed against them. One an undocumented
        # overlap took off (a SNV at another position on the record's REF) is not judged
        verdict, causes = documented(win, ctx, tgt, sample, lids)
        # (a row that names an earlier record of a split site shows the run's stores hold
        # it, against the documented limit: the rules are no model of that run there, and
        # such a row is a note, variant-shadowed)
        if listed and verdict == "off" and "undocumented" not in causes and \
                not any(a in win.shadowed for a in tgt.named):
            return True, "carrier-off-by-rule", "documented: %s" % (
                ",".join(sorted(causes)) or "the window's target of the slot's own alleles "
                                            "is another")
        zones = [z for x in tv[1] for z in x["clash"]]
        return False, "limit-overlapping-records", "; ".join(sorted(clash_text(win, zones)))
    # unknown phase where the target reads: the documented PUTATIVE rule is the truth
    allowed, strong, determined = putative_check(win, ctx, tgt, sample)
    if listed:
        if strong:
            return None
        verdict, causes = documented(win, ctx, tgt, sample, lids)
        if verdict in ("lists", "maybe"):
            if "mnv-unseen" in causes:
                return True, "variant-missing", MNV_UNSEEN
            return unexplained_or(causes)
        return True, ("carrier-without-allele" if not allowed else "carrier-not-strongest")
    if not (strong and determined):
        return None
    verdict, causes = documented(win, ctx, tgt, sample, lids)
    if verdict == "lists":
        return True, "carrier-missing"
    if verdict == "maybe":
        return None
    return unexplained_or(causes)


def unexplained_or(causes):
    """The limit a documented cause makes of a disagreement; with none, the disagreement
    is unexplained (gating: a defect of CRISPRme or of this checker). An overlap of the
    record explains nothing here: the truth reads a slot without the alleles METHODS
    section 4 reads it without, as the documented view does, and a slot with an overlap
    no reading leaves out is a clash, judged apart (carrier_verdict) -- so a cause one
    gives is another slot's, or the same reading's, and leaves the disagreement as it is."""
    cls = limit_of(set(causes) - OVERLAP_CAUSES)
    if cls:
        return False, cls, "documented: %s" % ",".join(sorted(causes))
    return True, "unexplained-disagreement"


def diagnose(info):
    """The class of a listed carrier whose known slots of the record lack the target."""
    if any(x["named"] for x in info):
        if any(x["named"] and x["extra"] for x in info):
            return "variant-missing"
        return "carrier-without-target"
    return "carrier-without-variants"


def local_lines(win, tgt):
    return win.lines_in(tgt.lo, tgt.hi)


def check_row(row, win, ctx, rep, gating=True):
    """Every soundness check of one indel/MNV row block."""
    if row["origin"] != "alt":
        rep.add(gating, "origin-not-alt", row,
                "REF/ALT_origin is %s on a row naming an indel-path record" % row["origin"])
    if not alignment_checks(row, ctx, rep, gating):
        return
    ref_text_check(row, win, rep, gating)
    named = named_alleles(row, win, rep, gating)
    if named is None:
        return
    mains = [a for a in named if a.kind != "snv"]
    if len(mains) != 1:
        rep.add(gating, "variants-not-one-indel", row,
                "the row names %d indel/MNV records; a fake-indel window holds one"
                % len(mains))
        return
    main = mains[0]
    tgt, cls, detail = row_target(win, row, named, main)
    if tgt is None:
        rep.add(gating, cls, row, detail)
        return
    for a in named:
        if a is not main and not any(s is a for s in tgt.src):
            rep.add(gating, "variant-outside-target", row,
                    "%s changes no base of the target (%d..%d)"
                    % (a.token(win.chrom), tgt.span[0], tgt.span[1]))
    samples = samples_of(row)
    tok = main.token(win.chrom)
    if not samples:
        if win.sites_only(main.key):
            # a row of a sites-only dataset (no genotypes): CRISPRme documents no carrier
            # on it, so its strings, target and variants are all there is to check. The
            # decision is the record's, not the run's: beside a genotyped dataset (the
            # default index with hg38_mega) that dataset's rows still list carriers
            rep.add(False, "limit-no-genotypes", row,
                    "%s is a record of a sites-only VCF (no sample columns), whose "
                    "dataset's rows list no carrier (CHANGELOG 2.8.1): the row's strings, "
                    "target and variants are checked, no carrier is%s"
                    % (tok, "; the record's carriers in the genotyped VCF(s) are "
                            "demanded by completeness" if win.genotyped(main.key) else ""))
            # its SNVs are those its fake-indel window codes: the dataset's own, whose
            # FILTER passes, the last record of a split site (sites_only_cps)
            vis = {win.line_vcf[a.line] for a in win.by_key.get(main.key, ())
                   if win.line_sites_only[a.line]}
            for a in named:
                if a.kind != "snv":
                    continue
                own = [x for x in win.by_key.get(a.key, ())
                       if win.line_sites_only[x.line] and win.line_vcf[x.line] in vis]
                if not own:
                    rep.add(gating, "variant-not-in-dataset", row,
                            "%s is no record of the sites-only VCF that holds %s, whose "
                            "fake-indel window codes that dataset's SNVs only"
                            % (a.token(win.chrom), tok))
                elif all(x.filtered is not None for x in own):
                    rep.add(gating, "variant-filtered", row,
                            "%s has FILTER %s in the sites-only VCF that holds %s: its "
                            "fake-indel window codes the enriched genome's SNVs, which "
                            "CRISPRitz's enricher keeps only where FILTER is PASS or '.'"
                            % (a.token(win.chrom), own[0].filtered, tok))
                elif all(x in win.so_shadowed for x in own if x.filtered is None):
                    rep.add(False, "variant-shadowed", row,
                            "%s is an earlier SNP record of a split site of the sites-only "
                            "VCF that holds %s, whose fake-indel window codes the last "
                            "(CHANGELOG 2.8.1 Known limits): this index holds it"
                            % (a.token(win.chrom), tok))
            return
        rep.add(gating, "no-carrier-listed", row, "the row lists no carrier")
        return
    lids = local_lines(win, tgt)
    listed = set(samples)
    out = defaultdict(list)
    why = defaultdict(set)
    if ctx.panel is not None:
        # CRISPRme lists no carrier outside the search's panel
        for s in sorted(listed - ctx.panel):
            out[(True, "carrier-not-in-panel")].append(s)
        listed &= ctx.panel
    pool = set(listed)
    for a in win.by_key.get(main.key, ()):
        pool.update(win.carriers(a))
    cache = {}
    for s in sorted(pool):
        sig = (s in listed, win.signature(s, lids))
        v = cache.get(sig)
        if v is None:
            v = carrier_verdict(win, ctx, tgt, s, s in listed, lids) or ("ok",)
            cache[sig] = v
        if v != ("ok",):
            out[v[:2]].append(s)
            if len(v) > 2 and v[2]:
                why[v[:2]].add(v[2])
    for (g, cls), who in sorted(out.items(), key=lambda kv: kv[0][1]):
        extra = sorted(why.get((g, cls), ()))
        detail = "%d sample(s) %s%s: %s" % (
            len(who), CARRIER_TEXT.get(cls, "%s") % tok,
            " (%s%s)" % ("; ".join(extra[:4]), "; ..." if len(extra) > 4 else "")
            if extra else "", ",".join(who[:12]))
        rep.add(g and gating, cls, row, detail, samples=",".join(who))


# the reasons a listed carrier of a SNP row lacks its target, a documented limit first
SNP_ROW_ORDER = ("limit-snp-row-unphased", "limit-snp-row-mnv", "limit-overlapping-records",
                 "limit-last-alt-wins", "snp-row-indel-carrier")


def snp_unphased(win, calls, span, kinds=("snv",)):
    """Whether a carrier's calls at a SNP row's target spanning ``span`` leave its
    sequence there of unknown phase. ``kinds`` ('snv'): an SNV alt it carries inside the
    target written unphased ('/', homozygous too) or as one allele (a haploid call) --
    CRISPRme's SNP rows then take it for an unphased individual, listed for every
    non-empty subset of the variants it carries there (METHODS section 4, "Behavior by
    phasing status"). ('indel', 'mnv'): an indel or MNV it carries inside the target with
    an unphased heterozygous call, which cannot be tied to a slot (an indel there keeps it
    listed, PUTATIVE: METHODS section 8). A missing allele of a phased call leaves its slot
    without the alt. Phase sets (PS) are not read: a '|' call is phased along the
    chromosome, as CRISPRme reads one with no PS. ('snv') A phased slot holding two ALTs of
    one position (every allele of a split site coded, 2.10.0) is no haplotype either: the
    observed path resolves its sample PUTATIVE, every combination of one ALT a position;
    and an alt of a line of several ALTs whose number CRISPRme's readers cannot tie
    (Window.untied) is read as one copy on an unknown slot, so its carrier is an unphased
    individual there too (METHODS section 4)."""
    held = defaultdict(set)   # (position, slot) -> the SNV ALTs a phased call puts there
    for lid, c in calls.items():
        if win.line_lo[lid] > span[1] or win.line_hi[lid] - 1 < span[0]:
            continue
        for a in win.line_alleles[lid]:
            if a is None or a.kind not in kinds or a.number not in c[0] or \
                    not touches(a, span[0], span[1]):
                continue
            if a.kind == "snv" and (not c[1] or len(c[0]) == 1 or a in win.untied):
                return True
            if a.kind != "snv" and not c[1] and len(set(c[0])) > 1:
                return True
            if a.kind == "snv" and a not in win.shadowed:
                for h, x in enumerate(c[0]):
                    if x == a.number:
                        held[(a.lo, h)].add(a.body)
    return any(len(v) > 1 for v in held.values())


def snp_slot_on(win, a, sample, h, missing=None):
    """Whether slot h of the sample holds allele ``a``, as CRISPRme's SNP path reads a
    genotype for a slot (METHODS section 4): True where every copy of its call is the
    allele (homozygous, or haploid) or a phased call puts it on h; False where it does not
    carry it, or a phased call puts another allele on h; ``missing`` where a phased call
    has a missing allele on h (an SNV's allele there is unknown, None; the phased indel
    store reads it as not the indel, False); None where the call cannot be tied to a slot
    (unphased)."""
    c = win.calls[a.line].get(sample)
    if c is None:
        return False
    alle, phased = c
    if a.number not in alle:
        return False
    if all(x == a.number for x in alle):
        return True
    if phased and h < len(alle):
        if alle[h] is None:
            return missing
        return alle[h] == a.number
    return None


def snp_row_keeps(win, sample, h, start, probe, named, lids):
    """Whether CRISPRme keeps slot h of the sample on a per-sample SNP row whose target is
    ``probe`` (+ reading) at ``start`` and names the SNPs ``named``, by the rule METHODS
    section 4 gives since 2.10.0 (CRISPRME_INDEL_SNP=1, the default). The slot is read
    against the phased indel store (the records whose REF and ALT lengths differ: no
    MNV): with no indel of it on the slot that changes the target's bases -- deletes or
    replaces a base of its span or inserts between two of them; one at a named SNP's
    position is another allele of that base, a phasing artefact, and changes nothing --
    it stays (also where such an indel's call cannot be tied to the slot: PUTATIVE).
    Otherwise it stays only where its own bases still read the target at its place, from
    the target's first base on or up to its last (a repeat unit gained or lost at the edge
    leaves them), the base of every named SNP in the stretch: written with every store
    indel on the slot that reaches the stretch read, over the target's length either side
    and the bases those indels add or lose, and with the slot's own SNP alleles outside
    the target (the row's inside it). A way that meets an indel whose call cannot be tied
    to the slot, two records of one base (two indels whose REFs share a base, a SNP of the
    slot on a base of an indel's REF), or a base whose allele the slot's calls leave
    unknown, reads nothing; with neither way holding, the carrier leaves the row."""
    n = len(probe)
    lo, hi = start, start + n - 1
    row = {a.lo: a.body for a in named}

    def indels_in(a, b):
        return [x for lid in win.lines_in(max(0, a), b + 1) for x in win.line_alleles[lid]
                if x is not None and x.kind == "indel" and x.lo not in row]
    mine = [a for a in indels_in(lo - 1, hi + 1)
            if touches(a, lo, hi) and snp_slot_on(win, a, sample, h, False)]
    if not mine:
        return True
    margin = n + sum(abs((a.hi - a.lo) - len(a.key[2])) for a in mine)
    first, last = lo - margin, hi + margin
    store = indels_in(first, last)

    def value(q):
        # the slot's base at q as CRISPRme's SNP stores give it: the row's allele inside
        # the target; outside, the slot's own SNP allele (the last record of a split site;
        # None where its slot cannot be told, an untied multi-ALT alt included)
        if lo <= q <= hi:
            return row.get(q, win.refseq[q])
        got = None
        for x in win.snv_at.get(q, ()):
            if x in win.shadowed:
                continue
            o = snp_slot_on(win, x, sample, h, None)
            if o is not False and x in win.untied:
                return None
            if o is None:
                return None
            if o:
                if got is not None:
                    return None
                got = x.body
        return got or win.refseq[q]

    for side in (0, 1):
        held, bad = list(mine), False
        for a in store:
            if a in mine or not (a.hi > first and a.lo <= last):
                continue
            if not ((a.hi - 1 >= lo) if side == 0 else (a.lo <= hi)):
                continue   # a way reads the indels that reach its stretch
            o = snp_slot_on(win, a, sample, h, False)
            if o is None:
                bad = True
                break
            if o:
                held.append(a)
        if bad:
            continue
        if any(value(q) != win.refseq[q] for a in held for q in range(a.lo, a.hi)):
            continue   # a SNP allele on a base of an indel's REF: two records of one base
        spans = sorted((a.lo, a.hi) for a in held)
        if any(y[0] < x[1] for x, y in zip(spans, spans[1:])):
            continue   # two indels whose REFs share a base
        m2 = n + sum(abs((a.hi - a.lo) - len(a.key[2])) for a in held)
        hap = make_hap(win.refseq, held, max(0, lo - m2), hi + m2 + 1)
        if side == 0:
            a0 = bisect.bisect_left(hap.place, (lo, 0))
        else:
            a0 = bisect.bisect_right(hap.place, (hi, 0)) - n
        if a0 < 0 or a0 + n > len(hap.place):
            continue
        places = hap.place[a0:a0 + n]
        if not all((q, 0) in places for q in row):
            continue   # a named SNP pushed out of the stretch, or its base deleted
        good = True
        for k in range(n):
            pl, src = places[k], hap.src[a0 + k]
            base = hap.seq[a0 + k] if src is not None else (
                value(pl[0]) if pl[1] == 0 else None)
            if base is None or base.upper() != probe[k].upper():
                good = False
                break
        if good:
            return True
    return False


def snp_merged_on(win, ctx, sample, h, keys, row, cp, lo, hi, lids):
    """Whether the slot of the sample with the alleles ``keys`` holds another alignment of
    the SNP row's site -- the same Cluster Position, an in-budget alignment of the row's
    guide whose target names exactly the row's variants (the SNP alts in it) -- that
    CRISPRme keeps it on (snp_row_keeps): the merge joins the carriers of a site's SNP rows
    of one variant list (its key for a SNP row is the Cluster Position and the variants,
    no alignment), so the row lists it as rightly as that alignment's row would."""
    guide = row.get("guide") or ""
    if not guide:
        return False
    full = ctx.full(guide)
    var = {t for t in row["var"].split(",") if t and t != "NA"}
    snvs = [a for a in keys if a.kind == "snv"]
    one = len(full) + max(ctx.bdna, ctx.brna) + 2
    hv = make_hap(win.refseq, snvs, max(0, lo), min(len(win.refseq), hi))
    seg_lo = bisect.bisect_left(hv.place, (cp - 2 * one, 0))
    seg_hi = bisect.bisect_left(hv.place, (cp + 2 * one, 0))
    seg = hv.seq[seg_lo:seg_hi]
    for strand, off, sa, ta, _mm, _dn, _rn in scan_sequence(
            full, seg, ctx.mm, ctx.bdna, ctx.brna, ctx.cap, ctx.pam_len, ctx.pam_5prime):
        w = len(ta.replace("-", ""))
        i = seg_lo + off
        place = hv.place[i:i + w]
        if cluster_of(place[0][0], strand, w, len(full), ctx.pam_5prime) != cp:
            continue
        got = [x for x in hv.src[i:i + w] if x is not None]
        if {x.token(win.chrom) for x in got} != var:
            continue
        if snp_row_keeps(win, sample, h, place[0][0], hv.seq[i:i + w], got, lids):
            return True
    return False


def snp_slot_reads(hap, start, probe, want):
    """Whether a slot's haplotype ``hap`` has a SNP row's target as CRISPRme reads a SNP
    row's carrier since 2.10.0 (METHODS section 4): the slot's bases read from the
    target's first base on, or up to its last, are the target's -- a repeat unit gained or
    lost at the target's edge leaves them as they are -- with the base of every SNP the row
    names (``want``: their keys) in the stretch read (an indel that pushes a named SNP
    out, or deletes its base, gives the target thanks to the indel, not the SNP)."""
    n = len(probe)

    def through(i):
        got = {x.key for x in hap.src[i:i + n] if x is not None}
        return want <= got

    if any(through(i) for i in reads_at(hap, start, probe)):
        return True
    j = hap.index((start + n - 1, 0))
    return (j is not None and j - n + 1 >= 0 and hap.seq[j - n + 1:j + 1] == probe
            and through(j - n + 1))


def check_snp_row(row, win, ctx, rep, snp_listed=None):
    """The same checks on a SNP-only row of a per-sample run, report-only: its own
    strings, its target, every named variant inside it, every listed carrier, and the
    carriers it leaves out -- a sample whose phased slot holds exactly the named SNPs
    inside the target, and no indel that changes it, is listed on a row of this site
    naming them (``snp_listed``: (guide, chromosome, Cluster Position, variants) -> the
    samples SNP rows there list; such rows merge their carriers). Since 2.10.0 a
    per-sample SNP row does not list a carrier whose haplotype holding its SNPs also holds
    an indel of the phased indel store that changes the target's bases, read through a
    repeat at its edge (METHODS section 4): listing one is a defect
    (snp-row-indel-carrier); an MNV there (in no store), an indel at a named SNP's
    position (a phasing artefact) and an unphased call are the documented exceptions. A
    row that lists no carrier and names only records a sites-only VCF holds is a
    sites-only dataset's, whose SNP post-analysis is the population-level one whatever
    the search's mode (CHANGELOG 2.8.1): a worst-case union, out of scope, so one line
    says so."""
    if not alignment_checks(row, ctx, rep, False):
        return
    ref_text_check(row, win, rep, False)
    named = named_alleles(row, win, rep, False)
    if named is None:
        return
    if not samples_of(row) and named and all(win.sites_only(a.key) for a in named):
        rep.add(False, "limit-no-genotypes", row,
                "a SNP row of a sites-only dataset (no carrier; %s held by a sites-only "
                "VCF): its post-analysis is the population-level one, a documented "
                "worst-case union, which the SNP-row report leaves out"
                % ",".join(a.token(win.chrom) for a in named))
        return
    hap = None
    stripped = row["alt"].replace("-", "").upper()
    probe = stripped if row["strand"] == "+" else rc(stripped)
    n = len(probe)
    lo = max(0, row["start"] - 2 * n - 60)
    hi = min(len(win.refseq), row["start"] + 3 * n + 60)
    hap = make_hap(win.refseq, named, lo, hi)
    hits = reads_at(hap, row["start"], probe)
    if not hits:
        rep.add(False, "target-not-on-variant-sequence", row,
                "the target is not the reference with its named variants at %d"
                % row["start"])
        return
    span = span_of(hap.place[hits[0]:hits[0] + n])
    src = hap.src[hits[0]:hits[0] + n]
    for a in named:
        if not any(s is a for s in src):
            rep.add(False, "variant-outside-target", row,
                    "%s changes no base of the target (%d..%d)"
                    % (a.token(win.chrom), span[0], span[1]))
    # the target's own SNPs: a slot holds the target with them (a named variant outside
    # it is the row's own fault, reported above)
    want = {a.key for a in named if any(s is a for s in src)}
    lids = win.lines_in(lo, hi)
    out = defaultdict(list)
    var = frozenset(x for x in row["var"].split(",") if x and x != "NA")
    cp = cluster_of(row["start"], row["strand"], n, len(row["spacer_col"]),
                    ctx.pam_5prime)
    inside = [a for a in named if a.key in want]
    for s in samples_of(row):
        if ctx.panel is not None and s not in ctx.panel:
            out["carrier-not-in-panel"].append(s)
            continue
        calls = {lid: win.calls[lid][s] for lid in lids if s in win.calls[lid]}
        if not calls:
            out["carrier-without-target"].append(s)
            continue
        if snp_unphased(win, calls, span):
            out["limit-snp-row-unphased"].append(s)
            continue
        ploidy = max(len(c[0]) for c in calls.values())
        ok, whys = False, set()
        for h in range(ploidy):
            keys = [a for lid, c in calls.items() for a in win.line_alleles[lid]
                    if a is not None and len(c[0]) > h and c[0][h] == a.number]
            if not all(any(a is k for k in keys) or any(a.key == k.key for k in keys)
                       for a in named):
                continue
            hh = make_hap(win.refseq, keys, lo, hi)
            hv = make_hap(win.refseq, [a for a in keys if a.kind == "snv"], lo, hi)
            if snp_slot_reads(hv, row["start"], probe, want):
                # the slot's SNP alleles give the target (the haplotype SNP rows list a
                # carrier for); since 2.10.0 CRISPRme reads the slot against the phased
                # indel store as METHODS section 4 says (snp_row_keeps), which is then
                # the truth of its bases -- but for an MNV, in no store CRISPRme reads,
                # and an indel or MNV inside the target whose call cannot be tied to the
                # slot (kept, PUTATIVE: its sequence there is not known)
                if snp_row_keeps(win, s, h, row["start"], probe, inside, lids):
                    if snp_unphased(win, calls, span, ("indel", "mnv")):
                        whys.add("limit-snp-row-unphased")
                        continue
                    if snp_slot_reads(hh, row["start"], probe, want):
                        ok = True
                        break
                    if any(a.kind == "mnv" and touches(a, span[0], span[1]) for a in keys):
                        whys.add("limit-snp-row-mnv")
                    elif any(touches(a, span[0], span[1]) for a in zone_alleles(hh)):
                        # two records of one base (an indel at a named SNP's position: a
                        # phasing artefact, which the rule keeps listed)
                        whys.add("limit-overlapping-records")
                    else:
                        ok = True   # read through a repeat at the edge, as the rule reads
                        break
                    continue
                if snp_merged_on(win, ctx, s, h, keys, row, cp, lo, hi, lids):
                    ok = True   # listed through another alignment of the site (below)
                    break
                whys.add("snp-row-indel-carrier")
                continue
            if snp_merged_on(win, ctx, s, h, keys, row, cp, lo, hi, lids):
                # the merge joins the carriers of a site's SNP rows of one variant list
                # (its key for a SNP row is the Cluster Position and the variants, no
                # alignment): the slot holds another alignment of this site with these
                # SNPs, which it is listed for
                ok = True
                break
            extra = [a for a in keys if touches(a, span[0], span[1])
                     and all(a.key != x.key for x in named)]
            if any(a in win.shadowed for a in extra):
                whys.add("limit-last-alt-wins")
            elif any(touches(a, span[0], span[1]) for a in zone_alleles(hh)):
                whys.add("limit-overlapping-records")
        if not ok:
            # a documented limit on some slot explains the listing; else a defect
            out[next((c for c in SNP_ROW_ORDER if c in whys),
                     "carrier-without-target")].append(s)
    # the carriers it leaves out (phased calls only: an unphased carrier is listed for
    # every subset it may hold, METHODS section 4, which is not checked here)
    site = (snp_listed or {}).get((row.get("guide", ""), row["chrom"], cp, var), set())
    pool = set()
    for a in named:
        for x in win.by_key.get(a.key, ()):
            pool.update(win.carriers(x))
    pool -= set(samples_of(row)) | set(site)
    for s in sorted(pool):
        calls = {lid: win.calls[lid][s] for lid in lids if s in win.calls[lid]}
        if not calls or snp_unphased(win, calls, span):
            continue
        touching = [lid for lid in calls if win.line_lo[lid] <= span[1]
                    and win.line_hi[lid] - 1 >= span[0]]
        ploidy = max(len(c[0]) for c in calls.values())
        for h in range(ploidy):
            snvs = [a for lid in touching for a in win.line_alleles[lid]
                    if a is not None and a.kind == "snv" and len(calls[lid][0]) > h
                    and calls[lid][0][h] == a.number and touches(a, span[0], span[1])]
            if {a.key for a in snvs} != want or any(a in win.shadowed for a in snvs):
                continue
            if not snp_slot_reads(make_hap(win.refseq, snvs, lo, hi), row["start"], probe,
                                  want):
                continue
            # CRISPRme (2.10.0) drops a carrier whose slot holds an indel of the phased
            # store that changes the target, as METHODS section 4 reads the slot
            if not snp_row_keeps(win, s, h, row["start"], probe, inside, lids):
                continue
            out["carrier-missing"].append(s)
            break
    text = {"snp-row-indel-carrier": "hold, on the slot of the row's SNPs, an indel that "
                                     "changes the target's bases (read through a repeat at "
                                     "its edge): their sequence there is the indel's, and "
                                     "since 2.10.0 a per-sample SNP row does not list such a "
                                     "carrier (METHODS section 4; report-only)",
            "limit-snp-row-mnv": "hold, on the slot of the row's SNPs, an MNV that changes "
                                 "the target: the phased indel store CRISPRme reads a SNP "
                                 "row's carriers against holds no MNV (METHODS section 8)",
            "carrier-missing": "are not listed, nor on another row of this site naming "
                               "these SNPs, and their phased slot holds exactly these SNPs "
                               "inside the target, which it reads, with no indel that "
                               "changes it (report-only)",
            "limit-last-alt-wins": "hold an earlier SNP record of a split site inside the "
                                   "target (CRISPRme reads the last record)",
            "limit-overlapping-records": "hold overlapping alleles inside the target (an "
                                         "indel at a named SNP's position on its slot is a "
                                         "phasing artefact, which CRISPRme keeps listed)",
            "limit-snp-row-unphased": "have a call of unknown phase inside the target: a "
                                      "SNP row lists every subset a carrier may hold "
                                      "(METHODS section 4), not checked here",
            "carrier-without-target": "hold no haplotype with this target and every named "
                                      "variant",
            "carrier-not-in-panel": "are not samples of the search's panel (--samples, "
                                    "else the run's .sampleID.txt)"}
    for cls, who in sorted(out.items()):
        rep.add(False, cls, row, "%d %ssample(s) %s: %s"
                % (len(who), "" if cls == "carrier-missing" else "listed ", text[cls],
                   ",".join(who[:12])), samples=",".join(who))


# --------------------------------------------------------------------------- #
# completeness
# --------------------------------------------------------------------------- #
class Best(object):
    """The best of a target's alignments: ``best`` the fewest mismatches+bulges,
    mismatches and bulges (each on its own), ``cfd`` the highest CFD (None without CFD),
    and ``keys`` every alignment's (mismatches, bulges, CFD), which the merge's sort keys
    are read from: what the rows standing for the target after the merge must match.
    ``real`` is False for a summary of several alignments (its minima are no one
    alignment's, so they give no key)."""

    __slots__ = ("best", "cfd", "keys")

    def __init__(self, best, cfd, real=True):
        self.best = dict(best)
        self.cfd = cfd
        self.keys = {(best["mm"], best["bulges"], cfd)} if real else set()

    def add(self, best, cfd, real=True):
        for c, v in best.items():
            if v < self.best[c]:
                self.best[c] = v
        if cfd is not None and (self.cfd is None or cfd > self.cfd):
            self.cfd = cfd
        if real:
            self.keys.add((best["mm"], best["bulges"], cfd))


class Demand(Best):
    """A target a carrier must be reported for: exact (its known haplotype) or the
    strongest of its PUTATIVE choices (any one of ``alts`` when the rule leaves a tie).
    ``best`` and ``cfd`` are the best over its alignments (and its alternatives), ``own``
    each alternative's own Best (by strand, start, + reading and variants), ``aligned``
    each alternative's alignments (a Target each), and ``k`` the rank of its first base
    among the bases an insertion adds (0 on any other base)."""

    __slots__ = ("kind", "alts", "k", "own", "aligned")

    def __init__(self, kind, alts, best, cfd, k):
        Best.__init__(self, best, cfd, real=False)   # add_alt notes each alignment
        self.kind = kind
        self.alts = alts  # [(Target, variant tokens frozenset)]
        self.k = k
        self.own = {}
        self.aligned = {}

    @property
    def total(self):
        return self.best["mm+bulges"]

    def add(self, best, cfd, k=None):
        Best.add(self, best, cfd, real=False)
        if k is not None:
            self.k = min(self.k, k)

    def add_alt(self, t, toks, best, cfd):
        """One alignment of alternative ``t`` (its own best, and the demand's)."""
        key = (t.strand, t.start, t.probe, toks)
        self.keys.add((best["mm"], best["bulges"], cfd))
        b = self.own.get(key)
        if b is None:
            self.own[key] = Best(best, cfd)
        else:
            b.add(best, cfd)
        ts = self.aligned.setdefault(key, [])
        if not any(x.crrna == t.crrna and x.aligned == t.aligned for x in ts):
            ts.append(t)

    def alignments(self, t, toks):
        """The alignments of alternative ``t`` (Targets: a PUTATIVE choice, and so
        CRISPRme's documented view of a carrier, can depend on the alignment)."""
        return self.aligned.get((t.strand, t.start, t.probe, toks)) or [t]

    def best_of(self, t, toks):
        """The Best a row standing for alternative ``t`` must match: its own alignments'
        (a PUTATIVE target the rule ties with others is compared with itself)."""
        return self.own.get((t.strand, t.start, t.probe, toks), self)


def edits(mm, bulges):
    """An alignment's value under each sorting criterion of the fewest block."""
    return {"mm+bulges": mm + bulges, "mm": mm, "bulges": bulges}


def completeness(win, ctx, rep, reported, region, popmode, stats=None):
    """The real-data brute force over the indel and MNV records of ``region``. ``stats``
    counts the records only sites-only VCFs hold (no sample carries them: nothing to
    demand, and they are counted, not skipped in silence) and those CRISPRme's search
    skips for their FILTER (PASS_FILTERS: no fake-indel window holds them, so no row of
    theirs is demanded; counted too)."""
    lo_r, hi_r = region
    if stats is not None:
        skipped = {key for lo, key, kind, _so in win.dropped
                   if kind != "snv" and lo_r <= lo < hi_r}
        skipped |= {a.key for a in win.alleles
                    if a.kind != "snv" and not a.searched and lo_r <= a.lo < hi_r}
        n = sum(1 for key in skipped if not win.searched(key))
        if n:
            stats["filtered"] = stats.get("filtered", 0) + n
    for gi, guide in enumerate(ctx.guides):
        full = ctx.full(guide)
        one = len(full) + max(ctx.bdna, ctx.brna) + 2   # one target and its bulges
        reach = 2 * one + 20
        for al in list(win.alleles):
            if al.kind == "snv" or not (lo_r <= al.lo < hi_r):
                continue
            if win.by_key[al.key][0] is not al:
                continue  # one record per key (a second line of it is the same record)
            if not win.searched(al.key):
                # an indel only the phased indel store holds (counted above): on its
                # carriers' haplotypes, a second indel-path record of their other
                # targets; no row of its own
                continue
            carriers = set()
            for a in win.by_key[al.key]:
                carriers.update(win.carriers(a))
            if not carriers:
                if gi == 0 and stats is not None and not win.genotyped(al.key):
                    stats["sites_only"] = stats.get("sites_only", 0) + 1
                continue
            lo = max(0, al.lo - reach)
            hi = min(len(win.refseq), al.hi + reach)
            if al.lo - one < 0 or al.hi + one > len(win.refseq):
                rep.add(False, "limit-contig-end", {
                    "table": "completeness", "chrom": win.chrom, "start": al.lo,
                    "var": al.token(win.chrom)},
                    "the record is within one target of the contig's end")
                continue
            comp_record(win, ctx, rep, reported, al, carriers, guide, full, lo, hi,
                        popmode, (al.edit_lo - one, al.edit_hi + one))


def scan_cached(cache, ctx, full, seq, skipped):
    hits = cache.get(seq)
    if hits is None:
        hits = list(scan_sequence(full, seq, ctx.mm, ctx.bdna, ctx.brna, ctx.cap,
                                  ctx.pam_len, ctx.pam_5prime, skipped))
        cache[seq] = hits
    return hits


def demand_tokens(win, al, place, src, present):
    """The variants CRISPRme names for a target of the record: the record and the SNV
    alleles whose bases it holds; and the other indel-path alleles that change it (one
    that wrote a base of it, or whose edit lies inside it)."""
    toks = {al.token(win.chrom)}
    others = []
    s0, s1 = span_of(place)
    for s in src:
        if s is not None and s is not al and s.kind == "snv":
            toks.add(s.token(win.chrom))
    for a in present:
        if a is not al and a.kind != "snv" and (any(x is a for x in src)
                                                or touches(a, s0, s1)):
            others.append(a)
    return frozenset(toks), others


def one_record_reading(win, main_al, slot_hap, alleles, t, toks, lo, hi, unknown=()):
    """The fake-indel window's target that CRISPRme lists a slot's carrier for where the
    slot's own target ``t`` needs a second indel-path record (METHODS section 4), or None.
    The window's reading of the slot -- the record and the slot's SNV alleles, no second
    indel-path record -- must hold t's bases at t's start with t's variants, and the slot
    must hold that reading's bases at its place, read from its first base on or up to its
    last base, through the record (record_core: a unit a tandem repeat gains or loses at
    the target's edge leaves them as they are; one gained between the target's start and
    the record pushes the record's bases out of the stretch, and such a carrier is not
    listed), and the stretch read meets no call of unknown phase and no allele of an
    overlap zone of the slot (``unknown``: either way, bases that cannot be read)."""
    hap = make_hap(win.refseq, [main_al] + [a for a in alleles if a.kind == "snv"], lo, hi)
    zone = zone_alleles(hap)
    for i in reads_at(hap, t.start, t.probe, main_al):
        place, src = hap.place[i:i + t.n], hap.src[i:i + t.n]
        sp = span_of(place)
        if any(touches(a, sp[0], sp[1]) for a in zone):
            continue   # two SNVs of one base where the window's target reads
        if demand_tokens(win, main_al, place, src, [])[0] != toks:
            continue
        if reading_holds(slot_hap, t.probe, place, record_core(place, src, main_al),
                         unknown):
            return Target(win, main_al, [main_al], t.probe, t.strand, place, src, t.crrna,
                          t.aligned)
    return None


def window_cut(win, ctx, main_al, place, strand, full_n):
    """Why CRISPRme's search cannot hold a target of record ``main_al`` whose bases lie at
    ``place`` on a haplotype carrying it, or None. A CRISPRitz index keeps a target's site,
    its PAM, not its alignment (CRISPRitz's docs/INDEX_FORMAT.md, sections 2 and 3): the
    site's guide region is the guide's length of bases from the PAM (``full_n``: the guide's
    length plus the PAM's), so it reaches past a target with RNA bulges and stops short of
    one with DNA bulges, by the difference of their lengths. 'reach' for one that runs past
    the record's fake-indel window (WINDOW_BEFORE reference bases before its first REF base,
    WINDOW_AFTER after its last: the separator N lies beyond, whatever the index), its own
    bases or its site's guide region (an index of either format drops a site whose guide
    region holds an N); 'edge' for one an INDELS index of format 1 leaves out: its site's
    pad, the index's budget of bases past the guide region's far end (the end away from
    the PAM), holds the separator N (CRISPRitz 2.9.1 drops such a site at build time;
    2.10.0, format 2, keeps it, its pad allowed to run into an N). Only where the run's
    records say the format (main refuses a run where they do not and it matters)."""
    lo_w = main_al.lo - WINDOW_BEFORE
    hi_w = main_al.hi + WINDOW_AFTER     # exclusive
    first, last = place[0], place[-1]
    if (first[1] == 0 and first[0] < lo_w) or last[0] >= hi_w:
        return "reach"
    left = (strand == "+") != ctx.pam_5prime      # the guide's far end on the left
    if left:
        d = first[0] - lo_w if first[1] == 0 else WINDOW_BEFORE
    else:
        d = hi_w - 1 - last[0]
    # the window's bases past the site's guide region (beyond the target's own end: the
    # target's bases less the site's)
    d += len(place) - full_n
    if d < 0:
        return "reach"
    sem = win.sem
    if sem is None:
        return None
    fmt, pad = sem.indels(win.vcfs[win.line_vcf[main_al.line]], win.chrom)
    if fmt != 1:
        return None
    pad = pad or (max(ctx.bdna, ctx.brna) + 1)
    return "edge" if d < pad else None


def slot_targets(win, ctx, full, cache, main_al, use, unknown, lo, hi, reach, skipped):
    """([(kind, Cluster Position, payload)], capped): the in-budget targets of record
    ``main_al`` on a slot of its carriers read with ``use`` (the slot's other known
    alleles, none of which meets the record's REF) beside the calls of unknown phase
    ``unknown``, both strands, each holding the record by the junction rule:
      'putative'  it touches a call of unknown phase: the PUTATIVE rule's (payload None)
      'zone'      an overlap zone of the slot (two of its alleles whose REFs meet: no
                  sequence there) changes it -- on the reading with no allele of the zone,
                  or on a reading of a zone in reach with some written (zone_readings) --
                  and it needs no other indel-path record of the slot: a row CRISPRme may
                  write or not (the window reads the stores' alleles there: the last record
                  of a split site, say; payload: the zones)
      'zone-noted' the same, of a target that needs another indel-path record: no row
                  (METHODS section 4 reads such a slot through that record on its bases,
                  which meets the two variants that overlap, and lists the carrier on no
                  row of it), but its carrier is counted with the others
      'exact'     the slot holds it and no other indel-path record of the slot changes it
      'window'    another indel-path record of the slot changes it, and CRISPRme's reading
                  of the slot against the window's target of the same bases, start and
                  variants holds (one_record_reading): the target is that one
      'second'    another indel-path record changes it and that reading does not hold:
                  never a row of CRISPRme's (one record per fake-indel window)
      'cut'       CRISPRme's search cannot hold it (window_cut): past the fake-indel
                  window, or at its edge in an INDELS index of format 1 (payload: why)
    (payload of 'exact', 'window' and 'second': Target, variant tokens, strand, + reading,
    mm, dna, rna, aligned spacer, aligned target). ``capped``: the zones in reach had more
    readings than READINGS_CAP, and only each of their alleles alone was read (no demand
    depends on the other readings: their targets are counted, and are rows a merge cluster
    may chain through)."""
    out = []
    hap = make_hap(win.refseq, [main_al] + use, lo, hi)
    zone = zone_alleles(hap)
    readings = [(use, hap, None)]
    # the other readings of the zones in reach: each with some of their alleles written
    near = [z for z in hap.clash
            if any(touches(a, reach[0], reach[1]) or (a.lo <= reach[1] and a.hi > reach[0])
                   for a in z[2])]
    capped = False
    if near:
        choices = [zone_readings(z[2]) for z in near]
        count = 1
        for c in choices:
            count *= len(c)
        if count > READINGS_CAP:
            capped = True
            combos = [[a] for z in near for a in z[2]]
        else:
            combos = [[a for x in combo for a in x] for combo in itertools.product(*choices)]
        inzone = {id(a) for z in near for a in z[2]}
        fixed = [a for a in use if id(a) not in inzone]
        for chosen in combos:
            if chosen:
                readings.append((fixed + chosen, make_hap(win.refseq, [main_al] + fixed + chosen,
                                                          lo, hi), chosen))
    for alleles, hr, chosen in readings:
        for strand, off, sa, ta, mm, dn, rn in scan_cached(cache, ctx, full, hr.seq,
                                                           skipped):
            width = len(ta.replace("-", ""))
            if off + width > len(hr.seq) or not holds_edit(hr, off, width, main_al):
                continue
            place = hr.place[off:off + width]
            src = hr.src[off:off + width]
            s0, s1 = span_of(place)
            hit = [z for z in hr.clash if any(touches(a, s0, s1) for a in z[2])]
            if chosen is not None:
                if not any(touches(a, s0, s1) for a in chosen):
                    continue   # the same target as on the first reading
                hit = hit or near
            cp = cluster_of(place[0][0], strand, width, len(full), ctx.pam_5prime)
            if chosen is None and any(touches(a, s0, s1) for a in unknown):
                out.append(("putative", cp, None))
                continue
            probe = hr.seq[off:off + width]
            toks, others = demand_tokens(win, main_al, place, src, alleles)
            # (a target no other indel-path record of the slot changes is read on the
            # record's own window, as are its places; one that needs another is no row of
            # it whatever its reach: the two-indel limit)
            cut = None if others else window_cut(win, ctx, main_al, place, strand,
                                                  len(full))
            if cut:
                out.append(("cut", None, cut))
                continue
            if hit:
                out.append(("zone-noted" if others else "zone", cp, hit))
                continue
            t = Target(win, main_al, [main_al], probe, strand, place, src, sa, ta)
            kind = "exact"
            if others:
                w = one_record_reading(win, main_al, hr, alleles, t, toks, lo, hi,
                                       list(unknown) + zone)
                if w is None:
                    kind = "second"
                else:
                    kind, t = "window", w
                    cut = window_cut(win, ctx, main_al, w.place, strand, len(full))
                    if cut:
                        out.append(("cut", None, cut))
                        continue
            out.append((kind, cp, (t, toks, strand, probe, mm, dn, rn, sa, ta)))
    return out, capped


def record_demands(win, ctx, al, carriers, full, lo, hi, reach):
    """The demands of one record: every target its carriers' haplotypes hold with it, by
    the carrier group (a call signature). ``reach``: the reference coordinates a target
    holding the record's edit can cover. A dict: 'demands' ((kind, key) -> samples),
    'demand_of' (-> Demand), 'clash_who'/'clash_how' (samples whose slot's alleles overlap
    where a target reads, and how), 'capped' (the PUTATIVE rule's coverage gap),
    'zone_capped' (carriers whose slot's overlap zones near the record had more readings
    than READINGS_CAP), 'skipped' (windows holding N), 'lids', 'cps': the Cluster
    Positions of the targets CRISPRme writes rows for -- the targets some carrier is
    demanded for (exact, the window's target a second indel-path record leaves readable,
    or PUTATIVE; never one that needs a second indel-path record), and those a slot whose
    overlapping alleles METHODS section 4 reads it without (a phasing artefact at the
    record's position, an allele only on bases it deletes, an MNV; it stays listed) holds
    so read, enumerated as a clean slot's: exact, or the window's target -- and
    'cps_open': those of the rows CRISPRme may write or not, which no documented rule
    settles, each with why (OPEN_CLASS): 'overlap' (the one-record reading of a slot whose
    record meets a SNV at another position, or an allele kept that meets one left out; a
    target an overlap zone of two other alleles of the slot changes, or that only a
    reading of such a zone holds; a target of a slot read without its overlap that touches
    a call of unknown phase, which CRISPRme resolves PUTATIVE), 'last-alt-wins' and
    'putative-untied' (a demanded target through a SNV CRISPRme's stores read otherwise:
    an earlier record of a split site, an alt of a multi-ALT line no genotype ties). No
    target that needs a second indel-path record is a row either way. A slot whose record
    meets an indel at another position (a complex allele) gives none: METHODS takes it off
    the target's rows."""
    lids = win.lines_in(lo, hi)
    groups = defaultdict(list)
    for s in sorted(carriers):   # a set: sorted, so that nothing depends on hash order
        _a, call = main_call(win, al.key, s)
        groups[(call, win.signature(s, lids))].append(s)
    cache = {}
    skipped = []
    demands = defaultdict(set)   # (kind, key of the demand) -> samples
    demand_of = {}
    clash_who = set()
    clash_how = set()
    clash_cps = set()     # rows CRISPRme writes for a slot read without its overlaps
    open_cps = defaultdict(set)   # rows it may write or not -> why (OPEN_CLASS)
    capped = set()
    zone_capped = set()
    cut = {}   # why (window_cut) -> the samples one of whose targets CRISPRme cannot hold

    def unread_overlap(members, main_al, present, unknown, kinds):
        """A slot whose record overlaps an allele no documented reading leaves out: no
        sequence. Its carriers are counted where a target of the record reads on it with
        the overlapping alleles left out (the one-record reading CRISPRme's window makes);
        its rows chain a cluster only as METHODS reads the slot -- none where it takes the
        slot off (a complex allele), left open where it says nothing -- and never one that
        needs a second indel-path record."""
        full_hap = make_hap(win.refseq, [main_al] + present, lo, hi)
        own = full_hap.clashes(main_al.lo, main_al.hi - 1)
        over = {id(a) for z in own for a in z[2] if a is not main_al}
        rest = [a for a in present if id(a) not in over]
        got, zc = slot_targets(win, ctx, full, cache, main_al, rest, unknown, lo, hi,
                               reach, skipped)
        for kind, _c, why in got:
            if kind == "cut":
                cut.setdefault(why, set()).update(members)
        got = [x for x in got if x[0] != "cut"]
        if zc:
            zone_capped.update(members)
        if got:
            clash_who.update(members)
            clash_how.update(clash_text(win, own))
            if "complex" not in kinds:
                for kind, c, _x in got:
                    if kind not in ("second", "zone-noted"):
                        open_cps[c].add("overlap")

    for (call, _sig), members in groups.items():
        s = members[0]
        main_al, call = main_call(win, al.key, s)
        slots, known = record_slots(call, main_al.number)
        need_putative = not known
        resolved_unknown = False
        if not known:
            # the record's slot is unknown, but an allele every slot of the carrier holds
            # (a homozygous call) is on it: one that overlaps the record and that no
            # documented reading leaves out leaves it no sequence, so it demands nothing
            # (counted where a target of the record reads without that allele)
            both, both_unknown = slot_alleles(win, s, main_al, call, None, lids)
            use, _left, kinds = slot_reading(main_al, both)
            if record_clash(main_al, use):
                unread_overlap(members, main_al, both,
                               [a for a in both_unknown if not meets_ref(a, main_al)], kinds)
                need_putative = False
        for h in (slots if known else ()):
            present, unknown = slot_alleles(win, s, main_al, call, h, lids)
            # a call of unknown phase whose REF meets the record's would leave the slot
            # with no sequence; the reading that gives it one puts it on the other slot
            # (truth()), so it changes nothing here
            unknown = [a for a in unknown if not meets_ref(a, main_al)]
            # the slot is written with its known alleles only: a target that touches a
            # call of unknown phase is the PUTATIVE rule's (putative_targets enumerates
            # every choice there, and those a known slot can hold), and no other target
            # depends on where those calls are -- none of them changes its bases
            in_reach = any(touches(a, reach[0], reach[1])
                           or (a.lo <= reach[1] and a.hi > reach[0]) for a in unknown)
            use, left, kinds = slot_reading(main_al, present)
            if record_clash(main_al, use):
                unread_overlap(members, main_al, present, unknown, kinds)
                continue
            need_putative = need_putative or (in_reach and not left)
            got, zc = slot_targets(win, ctx, full, cache, main_al, use, unknown, lo, hi,
                                   reach, skipped)
            for kind, _c, why in got:
                if kind == "cut":
                    cut.setdefault(why, set()).update(members)
            got = [x for x in got if x[0] != "cut"]
            if zc:
                zone_capped.update(members)
            if left:
                # METHODS section 4 reads the slot without its overlapping alleles (it
                # stays listed): read so, it is a clean slot, whose targets are counted,
                # not demanded (the reading is not turned into a gating demand), and chain
                # a cluster as CRISPRme writes them -- exact ones and the window's targets
                # a second indel-path record leaves readable; one that touches a call of
                # unknown phase (resolved PUTATIVE) or an overlap zone is a row it may
                # write or not, one that needs a second indel-path record none
                resolved_unknown = resolved_unknown or in_reach
                if got:
                    clash_who.update(members)
                    clash_how.update(clash_text(win, make_hap(
                        win.refseq, [main_al] + present, lo, hi).clashes(
                            main_al.lo, main_al.hi - 1)))
                for kind, c, _x in got:
                    if kind in ("exact", "window"):
                        clash_cps.add(c)
                    elif kind in ("putative", "zone"):
                        open_cps[c].add("overlap")
                continue
            for kind, c, x in got:
                if kind == "putative":
                    continue   # the PUTATIVE rule's (need_putative is set)
                if kind in ("zone", "zone-noted"):
                    # an overlap zone of two other alleles of the slot changes it (on some
                    # reading of the zone): no one sequence, counted, not demanded; a row
                    # CRISPRme may write or not unless it needs a second indel-path record
                    clash_who.update(members)
                    clash_how.update(clash_text(win, x))
                    if kind == "zone":
                        open_cps[c].add("overlap")
                    continue
                t, toks, strand, probe, mm, dn, rn, sa, ta = x
                # a target that needs a second indel-path record on the slot is demanded
                # too: a row reporting it satisfies it, and a miss is the documented limit
                # (one record per fake-indel window) -- unless CRISPRme's reading of the
                # slot against the window's target of the same bases, start and variants
                # holds (one_record_reading), where it documents listing the carrier: that
                # target is demanded exactly
                key = ("second" if kind == "second" else "exact", strand, t.start, probe,
                       toks)
                demands[key].update(members)
                keep_demand(demand_of, key, t, toks, edits(mm, dn + rn),
                            target_cfd(ctx, sa, ta))
        if need_putative:
            putative_targets(win, ctx, None, al, main_al, s, members, lids, lo, hi, full,
                             known, cache, skipped, demands, demand_of, capped, reach,
                             cut=cut)
        if resolved_unknown:
            # the PUTATIVE choices of a slot read without its overlap: rows CRISPRme may
            # write or not (its reading of such a slot is not turned into a demand)
            scratch, scratch_of = defaultdict(set), {}
            putative_targets(win, ctx, None, al, main_al, s, members, lids, lo, hi, full,
                             known, cache, skipped, scratch, scratch_of, capped, reach,
                             resolved=True)
            for d in scratch_of.values():
                for t, _x in d.alts:
                    open_cps[cluster_of(t.start, t.strand, t.n, len(full),
                                        ctx.pam_5prime)].add("overlap")
    cps = set(clash_cps)
    for key, d in demand_of.items():
        if d.kind != "second":
            for t, _x in d.alts:
                c = cluster_of(t.start, t.strand, t.n, len(full), ctx.pam_5prime)
                # a target holding a SNV CRISPRme's stores read otherwise (an earlier
                # record of a split site, an alt of a multi-ALT line no genotype ties) is
                # one it may write or not: its row reads the stores' allele there
                why = {("last-alt-wins" if x in win.shadowed else "putative-untied")
                       for x in t.src if x is not None and x.kind == "snv"
                       and (x in win.shadowed or x in win.untied)}
                if why:
                    open_cps[c] |= why
                else:
                    cps.add(c)
    return {"demands": demands, "demand_of": demand_of, "clash_who": clash_who,
            "clash_how": clash_how, "capped": capped, "zone_capped": zone_capped,
            "skipped": skipped, "lids": lids, "cps": frozenset(cps), "cut": cut,
            "cps_open": {c: frozenset(w) for c, w in open_cps.items() if c not in cps}}


def comp_record(win, ctx, rep, reported, al, carriers, guide, full, lo, hi, popmode,
                reach=None):
    """The demands of one record (record_demands, or the enumeration the cluster chain
    made of it on this window), checked against the tables."""
    reach = reach or (lo, hi)
    chrom = win.chrom
    e = win.enum.pop((guide, al.key), None)
    if e is None:
        e = record_demands(win, ctx, al, carriers, full, lo, hi, reach)
        # the record's rows as a cluster chain reads them: its carriers' targets and, for
        # a record a sites-only dataset holds too, that dataset's worst-case rows (kept
        # whole, so that no chain depends on which record the work met first)
        record_cps(win, ctx, guide, al, e=e)
    demands, demand_of, lids = e["demands"], e["demand_of"], e["lids"]
    clash_who, clash_how, capped, skipped = (e["clash_who"], e["clash_how"], e["capped"],
                                             e["skipped"])
    tok = al.token(chrom)
    resolve(win, ctx, rep, reported, al, carriers, lids, demands, demand_of, popmode,
            guide)
    if clash_who:
        who = sorted(clash_who)
        how = sorted(clash_how)
        rep.add(False, "limit-overlapping-records", {
            "table": "completeness", "chrom": chrom, "start": al.lo, "var": tok},
            "%d sample(s) hold %s on a slot whose alleles overlap where a target of it "
            "reads (one in budget reads there with the overlapping alleles left out), so "
            "the VCF gives no one sequence there and they are not demanded (overlapping "
            "on one slot: %s%s; %d of them listed on a row naming the record): %s"
            % (len(who), tok, "; ".join(how[:4]), "; ..." if len(how) > 4 else "",
               reported.on_rows(guide, who, tok), ",".join(who[:12])),
            samples=",".join(who))
    if capped:
        who = sorted(capped)
        rep.add(False, "coverage-putative-cap", {
            "table": "completeness", "chrom": chrom, "start": al.lo, "var": tok},
            "this checker's own coverage gap, not a CRISPRme limit: %d sample(s) have more "
            "than %d positions of unknown phase within a target's reach of %s, and their "
            "targets there that touch one are not enumerated (2^k haplotypes each): %s"
            % (len(who), PUTATIVE_CAP, tok, ",".join(who[:12])), samples=",".join(who))
    zone_capped = e.get("zone_capped")
    if zone_capped:
        who = sorted(zone_capped)
        rep.add(False, "coverage-overlap-cap", {
            "table": "completeness", "chrom": chrom, "start": al.lo, "var": tok},
            "this checker's own coverage gap, not a CRISPRme limit: %d sample(s) hold %s on "
            "a slot whose overlap zones (alleles whose REFs meet) within a target's reach "
            "have more than %d readings, and only each of their alleles alone was read for "
            "the targets only such readings hold -- counted, never demanded, and rows "
            "CRISPRme may write or not, which a merge cluster can chain through: %s"
            % (len(who), tok, READINGS_CAP, ",".join(who[:12])), samples=",".join(who))
    if skipped:
        rep.add(False, "limit-n-bases", {
            "table": "completeness", "chrom": chrom, "start": al.lo, "var": tok},
            "%d window(s) around %s hold an N and are not aligned (no base pairs with N)"
            % (len(skipped), tok))
    cut = e.get("cut") or {}
    if cut.get("reach"):
        who = sorted(cut["reach"])
        rep.add(False, "limit-window-reach", {
            "table": "completeness", "chrom": chrom, "start": al.lo, "var": tok},
            "%d sample(s) hold a target of %s in budget that runs past its fake-indel "
            "window (%d bases before the record's first REF base, %d after its last), or "
            "whose site's guide region does (a target of a %d-nt guide and PAM with a DNA "
            "bulge, or with an RNA bulge at the window's first base, say), which no index "
            "holds (METHODS section 4): not demanded (%d of them listed on a row naming "
            "the record): %s" % (len(who), tok, WINDOW_BEFORE, WINDOW_AFTER, len(full),
                                 reported.on_rows(guide, who, tok), ",".join(who[:12])),
            samples=",".join(who))
    if cut.get("edge"):
        who = sorted(cut["edge"])
        rep.add(False, "limit-window-edge", {
            "table": "completeness", "chrom": chrom, "start": al.lo, "var": tok},
            "%d sample(s) hold a target of %s whose site's guide region (the guide's "
            "length of bases from its PAM) ends within the INDELS index's budget of its "
            "fake-indel window's edge, which an index of format 1 (CRISPRitz 2.9.1 or "
            "earlier, the run's) does not hold and one of format 2 does (METHODS section "
            "9): not demanded (%d of them listed on a row naming the record): %s"
            % (len(who), tok, reported.on_rows(guide, who, tok), ",".join(who[:12])),
            samples=",".join(who))


def keep_demand(demand_of, key, t, toks, best, cfd):
    """Record one alignment of a demanded target: the target is one sequence at one
    start (or two copies of it one insertion holds), and its many alignments differ in
    edits and CFD, so the demand keeps the best of each (what the row standing for it
    must match)."""
    d = demand_of.get(key)
    if d is None:
        d = demand_of[key] = Demand(key[0], [(t, toks)], best, cfd, t.place[0][1])
    else:
        d.add(best, cfd, t.place[0][1])
    d.add_alt(t, toks, best, cfd)


def clash_text(win, zones):
    return {"+".join(a.token(win.chrom) for a in z[2]) for z in zones}


def target_cfd(ctx, sa, ta):
    if ctx.cfd is None:
        return None
    return ctx.cfd.score(sa, ta.upper(), ctx.pam_len)


def slot_allows(win, t, main_al, named, slot_sets):
    """Whether some known slot of the record can hold target ``t`` for some placement of
    its calls of unknown phase (placements); True when that cannot be told."""
    for present, unknown in slot_sets:
        unk = [a for a in unknown if touches(a, t.span[0], t.span[1])]
        can, _must = placements(win, t, main_al, present, unk, named)
        if can is None or can:
            return True
    return False


def allele_of(win, sample, q, body):
    """The SNV allele at reference coordinate q whose ALT is ``body`` that the sample's call
    holds: of two lines of one REF>ALT (a split site read as its last record), the one
    CRISPRme's stores read where the sample is on it, else the line it is on."""
    alts = [a for a in win.snv_at[q] if a.body == body]
    held = [a for a in alts if a.number in (win.call(a, sample) or ((), False))[0]]
    pool = held or alts
    return next((a for a in pool if a not in win.shadowed), pool[0])


def putative_targets(win, ctx, _unused, al, main_al, s, members, lids, lo, hi, full, known,
                     cache, skipped, demands, demand_of, capped, reach=None, resolved=False,
                     cut=None):
    """The PUTATIVE demands of one carrier group: every in-budget target its own alleles
    allow with the record (any choice at each SNV position: the documented rule ignores
    phase), grouped by alignment; the strongest choice of each is demanded. With the
    record's slot known, only alignments that touch a call of unknown phase are, and only
    choices some placement of those calls lets a slot of the record hold (its known
    alleles, a phased second indel say, can rule a choice out). A position on the
    record's own REF span is no choice: an allele there would describe the record's bases
    twice on its slot (no sequence), and one only on bases it deletes changes nothing
    (METHODS section 4); a call of unknown phase there is off the record's slot
    (truth()). A known slot whose alleles overlap the record demands nothing here, as it
    demands no exact target (record_demands counts its carriers): one with no reading
    holds no target, and one read without the overlapping alleles is CRISPRme's documented
    reading only, which this checker does not turn into a gating demand. ``resolved``:
    the same choices for the slots METHODS reads without their overlap, and those only,
    each read so -- the rows CRISPRme may write for them (record_demands leaves them
    open), not demands."""
    main_call_ = main_call(win, al.key, s)[1]
    unknown_all = []
    slot_sets = []
    if known:
        for h in record_slots(main_call_, main_al.number)[0]:
            present_h, unknown_h = slot_alleles(win, s, main_al, main_call_, h, lids)
            if any(meets_ref(a, main_al) for a in present_h):
                if not resolved:
                    continue
                use_h, left_h, _k = slot_reading(main_al, present_h)
                if not left_h or record_clash(main_al, use_h):
                    continue
                present_h = use_h
            elif resolved:
                continue
            unknown_h = [a for a in unknown_h if not meets_ref(a, main_al)]
            unknown_all.extend(unknown_h)
            slot_sets.append((present_h, unknown_h))
        if not slot_sets:
            return
    r_lo, r_hi = reach or (lo, hi)
    a_lo, a_hi = max(lo, r_lo), min(hi - 1, r_hi)
    positions = sorted({a.lo for lid in win.lines_in(a_lo, a_hi + 1)
                        for a in win.line_alleles[lid]
                        if a is not None and a.kind == "snv" and a_lo <= a.lo <= a_hi
                        and not main_al.lo <= a.lo < main_al.hi})
    ploidy = len(main_call_[0])
    choice = []
    for q in positions:
        opts = options(win, s, q, ploidy=ploidy)
        if len(opts) > 1:
            choice.append((q, sorted(opts)))
    if len(choice) > PUTATIVE_CAP:
        capped.update(members)
        return
    fixed = []
    for q in positions:
        opts = options(win, s, q, ploidy=ploidy)
        if len(opts) == 1:
            b = next(iter(opts))
            if b != win.refseq[q]:
                fixed.append(allele_of(win, s, q, b))
    groups = defaultdict(dict)  # alignment -> {probe: (Target, mm, bulges)}
    for combo in itertools.product(*[opts for _q, opts in choice]):
        chosen = list(fixed)
        for (q, _opts), b in zip(choice, combo):
            if b != win.refseq[q]:
                chosen.append(allele_of(win, s, q, b))
        hap = make_hap(win.refseq, [main_al] + chosen, lo, hi)
        for strand, off, sa, ta, mm, dn, rn in scan_cached(cache, ctx, full, hap.seq,
                                                           skipped):
            width = len(ta.replace("-", ""))
            if off + width > len(hap.seq) or not holds_edit(hap, off, width, main_al):
                continue
            place = hap.place[off:off + width]
            s0, s1 = span_of(place)
            if known and not any(touches(a, s0, s1) for a in unknown_all):
                continue
            why = window_cut(win, ctx, main_al, place, strand, len(full))
            if why:
                if cut is not None:
                    cut.setdefault(why, set()).update(members)
                continue
            probe = hap.seq[off:off + width]
            gap = "".join("-" if c == "-" else "x" for c in ta)
            g = groups[(strand, tuple(place), sa, gap)]
            if probe not in g:
                t = Target(win, main_al, [main_al], probe, strand, place,
                           hap.src[off:off + width], sa, ta)
                g[probe] = (t, mm, dn + rn)
    # the strongest choice depends on the alignment (its CFD), so the alignments under
    # which a target is this group's strongest are the group's own: a demand is the
    # targets the rule picks (one, or a tie) together with every alignment that picks
    # them, and carriers share it only when their alignments are the same (pooling
    # another group's alignments would hold this one to a target it is not listed for)
    picks = defaultdict(list)   # (strand, start, + readings, variants) -> [(aid, alts)]
    for aid, members_t in groups.items():
        alts = []
        for probe, (t, mm, bul) in members_t.items():
            if not putative_check(win, ctx, t, s)[1]:
                continue
            snvs = [x for x in t.src if x is not None and x is not main_al
                    and x.kind == "snv"]
            if known and not slot_allows(win, t, main_al, snvs, slot_sets):
                continue
            toks = frozenset([al.token(win.chrom)] + [x.token(win.chrom) for x in snvs])
            alts.append((t, toks, mm, bul))
        if not alts:
            continue
        picks[(aid[0], alts[0][0].start, tuple(sorted(t.probe for t, *_x in alts)),
               tuple(sorted(tuple(sorted(x)) for _t, x, *_y in alts)))].append((aid, alts))
    for (strand, start, probes, toks_all), found in picks.items():
        key = ("putative", strand, start, probes, toks_all,
               frozenset(aid for aid, _alts in found))
        demands[key].update(members)
        for _aid, alts in found:
            best = {c: min(edits(mm, bul)[c] for _t, _x, mm, bul in alts)
                    for c in ("mm+bulges", "mm", "bulges")}
            best_cfd = (max((target_cfd(ctx, t.crrna, t.aligned) or 0.0)
                            for t, *_x in alts) if ctx.cfd else None)
            k = min(t.place[0][1] for t, *_x in alts)
            if key not in demand_of:
                demand_of[key] = Demand("putative", [(t, x) for t, x, *_y in alts], best,
                                        best_cfd, k)
            else:
                demand_of[key].add(best, best_cfd, k)
            d = demand_of[key]
            for t, x, mm, bul in alts:
                d.add_alt(t, x, edits(mm, bul), target_cfd(ctx, t.crrna, t.aligned))


def resolve(win, ctx, rep, reported, al, carriers, lids, demands, demand_of, popmode,
            guide):
    """Each demand against the tables; a miss is classified by the documented view."""
    chrom = win.chrom
    full_len = len(ctx.full(guide))
    fold_sets = {}
    tok0 = al.token(chrom)

    def first_k(r):
        return reported.first_k(win, r)

    for key in sorted(demands, key=lambda k: (k[0], k[1], k[2], str(k[3]))):
        d = demand_of[key]
        who = sorted(demands[key])
        t0, toks0 = d.alts[0]
        cl = cluster_of(t0.start, t0.strand, t0.n, full_len, ctx.pam_5prime)
        row = {"table": "completeness", "chrom": chrom, "start": t0.start,
               "strand": t0.strand, "crrna": t0.crrna, "alt": t0.probe,
               "var": ",".join(sorted(toks0))}
        if popmode:
            # one row per cluster: the cluster must hold a row as good as the target, for
            # every carrier CRISPRme's documented rules report it for
            ok, why, cut, via = reported.locus(
                guide, chrom, cl, d, lambda wide: cluster_chain(
                    win, ctx, guide, cl, reported.merge, wide))
            if ok and via:
                # its cluster reaches a row as good as it only through rows CRISPRme may
                # write or not, which no documented rule settles: counted, not judged
                cls = next(c for k, c in OPEN_CLASS if k in via)
                rep.add(False, cls, row,
                        "%d sample(s) hold this target, and the rows CRISPRme writes chain "
                        "its cluster to no row as good as it (%s); it reaches one only "
                        "through rows CRISPRme may write or not: %s (%d of them listed on "
                        "a row naming the record): %s"
                        % (len(who), why, "; ".join(OPEN_TEXT[k] for k, _c in OPEN_CLASS
                                                     if k in via),
                           reported.on_rows(guide, who, tok0), ",".join(who[:12])),
                        samples=",".join(who))
                continue
            if not ok and cut:
                # the cluster runs past the 3 kb either side the chain follows: its row
                # may lie further, so this checker cannot tell (its own coverage gap)
                rep.add(False, "coverage-chain-limit", row,
                        "this checker's own coverage gap, not a CRISPRme limit: %d "
                        "sample(s) hold this target, whose merge cluster runs past the 3 "
                        "kb either side its chain follows, and %s: its row may lie "
                        "further (%d of them listed on a row naming the record): %s"
                        % (len(who), why, reported.on_rows(guide, who, tok0),
                           ",".join(who[:12])), samples=",".join(who))
                continue
            if d.kind == "second":
                if not ok:
                    rep.add(False, "limit-two-indels", row,
                            "%d sample(s) hold this target only with a second "
                            "indel-path record on the slot (one record per fake-indel "
                            "window), and its cluster %s (%d of them listed on a row "
                            "naming the record): %s"
                            % (len(who), why, reported.on_rows(guide, who, tok0),
                               ",".join(who[:12])), samples=",".join(who))
                continue
            if ok:
                continue
            # CRISPRme must report the target when its documented rules list any one of
            # the carriers for it (one call signature is enough to ask)
            causes, reach = set(), None
            seen = set()
            for s in who:
                sig = win.signature(s, lids)
                if sig in seen:
                    continue
                seen.add(sig)
                verdict, cs = documented_for(win, ctx, d, s, lids)
                causes |= cs
                if verdict != "off":
                    reach = s
                    break
            if reach is None:
                g, cls = unexplained_or(causes)[:2]
                rep.add(g, cls, row,
                        "%d sample(s): their cluster %s, and CRISPRme's documented rules "
                        "do not report this target for them (documented: %s; %d of them "
                        "listed on a row naming the record): %s"
                        % (len(who), why, ",".join(sorted(causes)),
                           reported.on_rows(guide, who, tok0), ",".join(who[:12])),
                        samples=",".join(who))
                continue
            rep.add(True, "missing-locus", row,
                    "%d sample(s) hold this target (%d edits%s) at Cluster Position %d "
                    "(%s among them, listed by CRISPRme's documented rules), and %s"
                    % (len(who), d.total, "" if d.cfd is None else ", CFD %.4f" % d.cfd,
                       cl, reach, why), samples=",".join(who))
            continue
        missing = []
        worse = {}
        listed_worse = set()   # listed for the target exactly, but not by its best
        for s in who:
            # a target that needs two indel-path records is reported, for the sample,
            # when a row lists it there under either record's window (the documented
            # limit is that it is not: its alignment is not compared)
            if d.kind == "second" and any(reported.has(guide, s, t.strand, t.start,
                                                       t.probe) for t, _x in d.alts):
                continue
            hit = False
            # the rows that report an alternative exactly for the sample, and whether
            # their fold groups (the rows of their site naming those variants with those
            # carriers, which the merge folds into the best by each block's order: two
            # alignments of one target with one carrier set always fold) are worse
            for t, toks in d.alts:
                listed, bad = reported.exact_for(guide, s, t, toks, d.best_of(t, toks),
                                                 first_k)
                if not listed:
                    continue
                if bad is None:
                    hit = True
                    break
                listed_worse.add(s)
                worse.setdefault(s, "the rows of its site naming its variants with the "
                                    "carriers of the rows that list them for it: %s" % bad)
            if hit:
                continue
            for t, toks in d.alts:
                hcl = cluster_of(t.start, t.strand, t.n, full_len, ctx.pam_5prime) + (
                    t.place[0][1] if d.kind == "putative" else d.k)
                # the carriers the target's row must list and may list, under each of its
                # alignments (a PUTATIVE choice can depend on it): documented_carriers
                sets = []
                for ta in d.alignments(t, toks):
                    fk = (ta.strand, ta.start, ta.probe, toks, ta.crrna, ta.aligned)
                    if fk not in fold_sets:
                        fold_sets[fk] = documented_carriers(win, ctx, ta, carriers, lids)
                    if fold_sets[fk] not in sets:
                        sets.append(fold_sets[fk])
                for req, opt in sets:
                    ok, why = reported.folded(guide, s, hcl, t.strand, toks, req,
                                              opt | set(who), d.best_of(t, toks), first_k)
                    if ok:
                        hit = True
                        break
                    if why and s not in worse:
                        worse[s] = why
                if hit:
                    break
            if not hit:
                missing.append(s)
        if not missing:
            continue
        out = defaultdict(list)
        if d.kind == "second":
            rep.add(False, "limit-two-indels", row,
                    "%d sample(s) hold this target only with a second indel-path record "
                    "on the slot, and no row reports it for them (one record per "
                    "fake-indel window; %d of them listed on a row naming the record): %s"
                    % (len(missing), reported.on_rows(guide, missing, tok0),
                       ",".join(missing[:12])), samples=",".join(missing))
            continue
        for s in missing:
            verdict, causes = documented_for(win, ctx, d, s, lids)
            if verdict in ("lists", "maybe") and s in listed_worse:
                if any(c.startswith("putative") for c in causes) and \
                        listed_by_the_rule(win, ctx, reported, guide, d, s, first_k):
                    # resolved PUTATIVE against its known slot (an untied multi-ALT alt,
                    # an MNV), the carrier is listed for this target only under the
                    # alignments where the rule picks it, and a row reports it no worse
                    # than the best of those
                    out[(False, "limit-putative-carrier",
                         "%s; the PUTATIVE rule picks this target only under alignments "
                         "worse than its best, and a row lists it by the best of those"
                         % ",".join(sorted(causes)))].append(s)
                    continue
                # listed for the target, by an alignment worse than its best: the merge
                # kept the worse of its alignments (or lost the best)
                out[(True, "missing-best-alignment", "")].append(s)
            elif verdict == "lists":
                out[(True, "missing-row" if d.kind == "exact" else "missing-putative-row",
                     "")].append(s)
            elif verdict == "maybe":
                if reported.listed_at(guide, s, cl, al.token(chrom)):
                    # CRISPRme's rule ties this target with another and the sample is
                    # listed at the site naming the record: likely on the other one
                    out[(False, "limit-putative-carrier",
                         "%s; a PUTATIVE tie, listed at the site on another row"
                         % ",".join(sorted(causes)))].append(s)
                    continue
                out[(True, "missing-putative-row", "")].append(s)
            else:
                g, cls = unexplained_or(causes)[:2]
                out[(g, cls, ",".join(sorted(causes)))].append(s)
        for (g, cls, why), ss in sorted(out.items(), key=lambda kv: kv[0][1]):
            text = ("hold%s this target at %d (+ reading %s, best %d edits%s)%s and no row "
                    "reports it for them no worse than that: neither one of exactly this "
                    "target and variants whose rows of its site with its carriers are no "
                    "worse, nor a fold into rows of its site naming the same variants with "
                    "the target's carriers%s"
                    % ("s" if g else "", t0.start, t0.probe, d.total,
                       "" if d.cfd is None else ", CFD %.4f" % d.cfd,
                       " as the strongest PUTATIVE target" if d.kind == "putative" else "",
                       "" if g else " (documented: %s)" % why))
            if g:  # a failure is one line per (sample, target)
                for x in ss:
                    rep.add(True, cls, row, "%s %s%s" % (
                        x, text, " (%s)" % worse[x] if x in worse else ""), samples=x)
            else:
                rep.add(False, cls, row, "%d sample(s) %s (%d of them listed on a row "
                        "naming the record): %s"
                        % (len(ss), text, reported.on_rows(guide, ss, tok0),
                           ",".join(ss[:12])), samples=",".join(ss))


def documented_for(win, ctx, d, s, lids):
    """The documented verdict on a demand: 'lists' if CRISPRme's rules report any of its
    alternatives for the sample under any of its alignments, 'maybe' if they leave the
    choice open."""
    seen = set()
    best = "off"
    for t0, toks in d.alts:
        for t in d.alignments(t0, toks):
            t.named = [t.main] + [x for x in t.src if x is not None and x is not t.main
                                  and x.kind == "snv"]
            verdict, causes = documented(win, ctx, t, s, win.lines_in(t.lo, t.hi))
            seen |= causes
            if verdict == "lists":
                return "lists", seen
            if verdict == "maybe":
                best = "maybe"
    return best, seen


def listed_by_the_rule(win, ctx, reported, guide, d, s, first_k):
    """Whether a row reports the sample for an alternative of demand ``d`` no worse than
    the best of the alignments under which CRISPRme's documented rules list it there
    (documented_best): a PUTATIVE carrier's target, listed where the rule picks it."""
    for t0, toks in d.alts:
        db = documented_best(win, ctx, d, t0, toks, s)
        if db is not None and reported.exact_for(guide, s, t0, toks, db,
                                                 first_k) == (True, None):
            return True
    return False


def documented_best(win, ctx, d, t0, toks, s):
    """The Best of alternative ``t0`` of a demand over the alignments under which
    CRISPRme's documented rules list the sample for it (a PUTATIVE choice depends on the
    alignment), or None when there are none."""
    best = None
    for t in d.alignments(t0, toks):
        if t.crrna is None or t.aligned is None:
            continue
        t.named = [t.main] + [x for x in t.src if x is not None and x is not t.main
                              and x.kind == "snv"]
        if documented(win, ctx, t, s, win.lines_in(t.lo, t.hi))[0] == "off":
            continue
        mm, dn, rn, _p = recount(t.crrna, t.aligned, ctx.pam_len, ctx.pam_5prime)
        e, cfd = edits(mm, dn + rn), target_cfd(ctx, t.crrna, t.aligned)
        if best is None:
            best = Best(e, cfd)
        else:
            best.add(e, cfd)
    return best


def truth_lists(win, ctx, t, s, lids):
    """What the sample's own haplotypes say of its listing for target ``t`` (one
    alignment, t.named set): 'yes' (a slot of the record holds it through the record with
    its variants, or it is the sample's strongest PUTATIVE target, the rule's choice
    fixed), 'maybe' (a tie the PUTATIVE rule leaves), 'resolved' (only a slot whose
    alleles overlap the record holds it, read without them as METHODS section 4 reads
    it), 'clash' (a slot with no reading there may) or 'no'."""
    tv = truth(win, t, s, lids)
    if tv[0] == "exact":
        return "yes" if tv[1] else "no"
    if tv[0] in ("clash", "resolved"):
        return tv[0]
    if tv[0] == "unknown":
        _allowed, strong, determined = putative_check(win, ctx, t, s)
        return ("yes" if determined else "maybe") if strong else "no"
    return "no"


def documented_carriers(win, ctx, t, carriers, lids):
    """(required, allowed): the carriers the rows of target ``t``'s fold group (under its
    alignment ``t``) must list -- those CRISPRme's documented rules list for it whose own
    haplotypes hold it too -- and those they may list besides: the ones the rules list or
    leave open (PUTATIVE) and the ones whose haplotypes hold it while the rules leave them
    off. A carrier the rules list for a sequence its haplotypes do not hold is allowed and
    not required: those rules are this checker's re-implementation, which only
    classifies, and one carrier it gets wrong must not veto the fold of every true
    carrier of the target. A carrier whose slot's alleles overlap the record is never
    required (the reading of such slots is not turned into a gating demand); it is
    allowed where the slot, read without them as METHODS section 4 reads it, holds the
    target ('resolved': the one reading that holds the record), and a slot with no reading
    ('clash': the VCF gives no sequence) only where the documented view lists it or leaves
    it open."""
    req, opt = set(), set()
    lids = win.lines_in(t.lo, t.hi)
    t.named = [t.main] + [x for x in t.src if x is not None and x is not t.main
                          and x.kind == "snv"]
    cache = {}
    for s in sorted(carriers):
        sig = win.signature(s, lids)
        v = cache.get(sig)
        if v is None:
            v = cache[sig] = (documented(win, ctx, t, s, lids)[0],
                              truth_lists(win, ctx, t, s, lids))
        verdict, tr = v
        if tr == "resolved":
            opt.add(s)
        elif tr == "clash":
            if verdict != "off":
                opt.add(s)
        elif verdict == "lists" and tr == "yes":
            req.add(s)
        elif verdict != "off" or tr != "no":
            opt.add(s)
    return frozenset(req), frozenset(opt)


# --------------------------------------------------------------------------- #
# a merge cluster's extent (one row per cluster): the alignments CRISPRme chains
# --------------------------------------------------------------------------- #
def _feasible_sets(matches, target, n, m, max_mm, max_dna, max_rna, max_total, pam_first,
                   pam_last, pam_5prime):
    """_feasible on a target whose bases are sets (IUPAC: any allele of a position): an
    alignment exists if some choice of each position's alleles has one in budget. Used
    only to find where a merge cluster can extend, so it is kept apart from the aligner's
    hot path."""
    inf = max_mm + 1
    cost = {(0, 0): 0}
    i = 0
    for i in range(n + 1):
        dna_ok = (i <= pam_first) if not pam_5prime else (i > pam_last)
        if dna_ok and i > 0 and max_dna:
            for dn in range(max_dna):
                for rn in range(max_rna + 1):
                    c = cost.get((dn, rn))
                    j = i + dn - rn
                    if c is None or j >= m:
                        continue
                    if c < cost.get((dn + 1, rn), inf):
                        cost[(dn + 1, rn)] = c
        if i == n:
            break
        nxt = {}
        in_pam = pam_first <= i <= pam_last
        for (dn, rn), c in cost.items():
            j = i + dn - rn
            if j < m:
                ok = not target[j].isdisjoint(matches[i])
                if ok or not in_pam:
                    c2 = c + (0 if ok else 1)
                    if c2 <= max_mm and c2 < nxt.get((dn, rn), inf):
                        nxt[(dn, rn)] = c2
            if not in_pam and rn < max_rna and c < nxt.get((dn, rn + 1), inf):
                nxt[(dn, rn + 1)] = c
        cost = nxt
        if not cost:
            return False
    for (dn, rn), c in cost.items():
        if i + dn - rn == m and c <= max_mm and c + dn + rn <= max_total:
            return True
    return False


CHUNK = 64   # the cluster chain reads the reference's SNP rows a piece this size at a time


def snp_row_sets(win, ctx, lo, hi):
    """[lo, hi) of the reference as the reference and SNP rows CRISPRme can write there see
    it: base-set sequences (snp_row_seqs, whatever indels may take a per-sample row's
    carriers off)."""
    return [sets for sets, _drops in snp_row_seqs(win, ctx, lo, hi)]


def snp_row_seqs(win, ctx, lo, hi):
    """[(sets, drops)]: [lo, hi) of the reference as the reference and SNP rows CRISPRme can
    write there see it: base-set sequences (any base of a set, as CRISPRitz reads an IUPAC
    code). A population-level run's SNP rows are the worst case of each engine alignment
    over the SNP alleles of the enriched genome (new_simple_analysis's worst-case
    representative), so every allele of a dataset is an option at its position -- one
    dataset at a time, as a search enriches each apart; a per-sample run reports only a
    haplotype some individual holds (its PER-SAMPLE CONTRACT): a phased slot's own
    alleles, and at a call of unknown phase either of the sample's alleles (an unphased
    individual is listed for every subset of its alleles), beside the reference rows; a
    sites-only dataset's SNP rows are population-level whatever the mode (CHANGELOG
    2.8.1). Only the alleles CRISPRme's stores hold count (not an earlier record of a split
    site). ``drops``: the indels of a per-sample slot's own dataset that its call puts on
    that slot (phased, or homozygous): since 2.10.0 a per-sample SNP row does not list a
    carrier whose haplotype holding its SNPs holds an indel of the phased indel store that
    changes the target's bases (METHODS section 4; one whose indel call is unphased
    stays), so the slot's targets such an indel touches are rows CRISPRme may write or not
    (a repeat read through keeps the carrier), and the rest rows it writes. A sequence
    some slot gives with no such indel has none (empty)."""
    base = [frozenset((b,)) for b in win.refseq[lo:hi]]
    out = {tuple(base): {frozenset()}}
    by_vcf = defaultdict(list)
    indels = defaultdict(list)   # dataset -> its indels whose REF meets [lo - 1, hi)
    for lid in win.lines_in(max(0, lo - 1), hi):
        for a in win.line_alleles[lid]:
            if a is None:
                continue
            if a.kind == "snv" and lo <= a.lo < hi and a not in win.shadowed:
                by_vcf[win.line_vcf[lid]].append(a)
            elif a.kind == "indel" and a.lo < hi and a.hi > lo - 1:
                indels[win.line_vcf[lid]].append(a)
    for vi, als in sorted(by_vcf.items()):
        if ctx.snp_union or not win.vcfs[vi].samples:
            sets = list(base)
            # a sites-only dataset's SNV outside PASS_FILTERS is a row CRISPRme may write or
            # not (filtered_sets): hg38_mega's index holds such SNVs, one CRISPRitz's
            # enricher builds does not
            for a in [a for a in als if a.filtered is None]:
                sets[a.lo - lo] = sets[a.lo - lo] | {a.body}
            out.setdefault(tuple(sets), set()).add(frozenset())
            continue
        at = defaultdict(list)    # position -> the dataset's alleles there
        for a in als:
            at[a.lo].append(a)
        samples = set()
        for a in als:
            samples.update(win.calls[a.line])
        for s in sorted(samples):
            # the sample's calls in the piece; a call with fewer copies than its others
            # (a haploid call against a diploid one) is of unknown phase (METHODS 4)
            cs_at = {}
            for q, here in at.items():
                cs = [(a, win.calls[a.line].get(s)) for a in here]
                cs = [(a, c) for a, c in cs if c is not None]
                if cs:
                    cs_at[q] = cs
            ploidy = max([len(c[0]) for cs in cs_at.values() for _a, c in cs] or [1])
            known, open_ = {}, {}
            # a phased slot holding two ALTs of one position (every allele coded) is no
            # haplotype: the sample is PUTATIVE in the piece, each of its positions open
            two = any(sum(1 for a, c in cs
                          if c[1] and len(c[0]) > h and c[0][h] == a.number) > 1
                      for cs in cs_at.values() for h in range(ploidy))
            for q, cs in cs_at.items():
                if not two and all(c[1] and None not in c[0] and len(c[0]) == ploidy
                                   for _a, c in cs):
                    known[q] = {h: a.body for a, c in cs for h, x in enumerate(c[0])
                                if x == a.number}
                else:
                    open_[q] = frozenset(options(win, s, q, visible_only=True,
                                                 ploidy=ploidy))
            for h in range(ploidy):
                sets = list(base)
                for q, by_slot in known.items():
                    if h in by_slot:
                        sets[q - lo] = frozenset((by_slot[h],))
                for q, opts in open_.items():
                    sets[q - lo] = opts
                # the slot's own indels (Track B): on it by a phased call of the sample's
                # ploidy, or a homozygous one; an unphased call keeps the carrier listed
                drops = []
                for a in indels.get(vi, ()):
                    c = win.calls[a.line].get(s)
                    if c is None or len(c[0]) != ploidy or None in c[0]:
                        continue
                    if all(x == a.number for x in c[0]) or (c[1] and c[0][h] == a.number):
                        drops.append(a)
                out.setdefault(tuple(sets), set()).add(frozenset(drops))
    res = []
    for x in sorted(out, key=lambda x: [sorted(y) for y in x]):
        ds = out[x]
        if frozenset() in ds:
            res.append((list(x), frozenset()))
        else:
            for d in sorted(ds, key=lambda d: sorted(a.key for a in d)):
                res.append((list(x), d))
    return res


def capped_sets(win, ctx, lo, hi):
    """[lo, hi) of the reference as a capped window's greedy worst-case row sees it, for
    each genotyped dataset of a per-sample run: every allele of the dataset an option at
    its position (the representative takes, at each column, the allele that most lowers
    the mismatch count, METHODS section 5). Under --per-sample a capped window emits it
    only on the classic dictionary path (an index built without --samplesID, as the
    benchmark's complete-test builds; a documented known limit), which a run does not say
    it took: these rows are ones CRISPRme may write or not. A run of CRISPRme 2.10.0 writes
    none (Semantics.capped): its classic dictionaries' per-sample windows over the cap are
    the haplotypes their carriers have, as on the tiers."""
    if ctx.snp_union or (ctx.sem is not None and not ctx.sem.capped):
        return []
    by_vcf = defaultdict(list)
    for lid in win.lines_in(lo, hi):
        for a in win.line_alleles[lid]:
            if a is not None and a.kind == "snv" and lo <= a.lo < hi and \
                    a not in win.shadowed and win.vcfs[win.line_vcf[lid]].samples:
                by_vcf[win.line_vcf[lid]].append(a)
    out = []
    for _vi, als in sorted(by_vcf.items()):
        if len({a.lo for a in als}) <= IUPAC_CAP:
            continue   # no window of the piece can hold more than the cap
        sets = [frozenset((b,)) for b in win.refseq[lo:hi]]
        for a in als:
            sets[a.lo - lo] = sets[a.lo - lo] | {a.body}
        out.append(sets)
    return out


def filtered_sets(win, ctx, lo, hi):
    """[(sets, alleles)]: [lo, hi) of the reference as a sites-only dataset's SNP rows see
    it with its SNVs whose FILTER is outside PASS_FILTERS (``alleles``) among the options:
    the hg38_mega index holds such SNVs (its registry holds every one of them, and its rows
    name them: chr22:11574550 C>G, LOW_VQSLOD), an index CRISPRitz's enricher builds from a
    sites-only VCF does not, and a run does not say which it read -- so an alignment
    through one is a row CRISPRme may write or not."""
    by_vcf = defaultdict(list)
    for lid in win.lines_in(lo, hi):
        if not win.line_sites_only[lid]:
            continue
        for a in win.line_alleles[lid]:
            if a is not None and a.kind == "snv" and lo <= a.lo < hi:
                by_vcf[win.line_vcf[lid]].append(a)
    out = []
    for _vi, als in sorted(by_vcf.items()):
        filt = [a for a in als if a.filtered is not None]
        if not filt:
            continue
        sets = [frozenset((b,)) for b in win.refseq[lo:hi]]
        for a in als:
            sets[a.lo - lo] = sets[a.lo - lo] | {a.body}
        out.append((sets, filt))
    return out


def snp_row_cps(win, ctx, guide, full, a, b):
    """(certain, capped, dropped, filtered): the Cluster Positions of the in-budget alignments of the
    guide, both strands, on [a, b) of the reference as its reference and SNP rows see it
    (snp_row_seqs), those of a per-sample run's capped windows' greedy worst-case rows
    (capped_sets: alignments over more than CRISPRME_IUPAC_CAP variant positions), and
    those of a per-sample slot's SNP targets that an indel of the slot touches (rows
    CRISPRme writes only where its reading of the slot through a repeat keeps the
    carrier), and those of a sites-only dataset's SNP rows through one of its SNVs whose
    FILTER is outside PASS_FILTERS (filtered_sets), read a piece at a time -- each piece
    reaching one target and its bulges into the next, so that every alignment lies whole
    in the piece where it starts -- and kept on the window."""
    span = len(full) + max(ctx.bdna, ctx.brna) + 2
    out, capped, dropped, filtered, degenerate = set(), set(), set(), set(), set()
    for c in range(max(0, a) // CHUNK, max(a, b - 1) // CHUNK + 1):
        got = win.snpcps.get((guide, c))
        if got is None:
            lo = c * CHUNK
            hi = min(len(win.refseq), lo + CHUNK + span)
            places = [(q, 0) for q in range(lo, hi)]
            rows, maybe, degen = set(), set(), set()
            k = ctx.max_degeneracy
            for sets, drops in snp_row_seqs(win, ctx, lo, hi):
                # --max-degeneracy K: a SNP row resting on more than K ambiguous columns is
                # cut back or left out (METHODS section 5), so an alignment over more than
                # K positions of several options is a row CRISPRme may write or not
                cut = k is not None and any(len(x) > 1 for x in sets)
                if drops:
                    sure, touched = feasible_cps(ctx, full, sets, places, drops=drops)
                    if cut:
                        keep, _t = feasible_cps(ctx, full, sets, places, drops=drops,
                                                sparse=k)
                        degen |= sure - keep
                        sure = keep
                    rows |= sure
                    maybe |= touched
                else:
                    got_rows = feasible_cps(ctx, full, sets, places)
                    if cut:
                        keep = feasible_cps(ctx, full, sets, places, sparse=k)
                        degen |= got_rows - keep
                        got_rows = keep
                    rows |= got_rows
            dense = set()
            for sets in capped_sets(win, ctx, lo, hi):
                dense |= feasible_cps(ctx, full, sets, places, dense=IUPAC_CAP)
            filt = set()
            for sets, als in filtered_sets(win, ctx, lo, hi):
                filt |= feasible_cps(ctx, full, sets, places, drops=als)[1]
            got = win.snpcps[(guide, c)] = (rows, dense - rows, maybe - rows, filt - rows,
                                            degen - rows)
        out |= got[0]
        capped |= got[1]
        dropped |= got[2]
        filtered |= got[3]
        degenerate |= got[4]
    return out, capped - out, dropped - out, filtered - out, degenerate - out


def sites_only_cps(win, ctx, full, al, lo, hi):
    """The Cluster Positions of the rows a sites-only dataset that holds record ``al``
    writes for it: the worst case of each alignment of its fake-indel window over that
    dataset's SNP alleles (CHANGELOG 2.8.1), the window being the reference over [lo, hi)
    with the record's ALT written verbatim and an SNV's alleles coded at a flank base --
    a reference base off the record's REF span: the overlay codes no base of the ALT, the
    anchor included (METHODS section 4), so an SNV at the record's own position is in no
    row of the window. The overlay codes the enriched genome's SNVs: those whose FILTER
    passes (PASS_FILTERS): every allele of a split site on a genome of CRISPRitz 2.10.0, the
    last record of one on a genome of 2.9.1 or earlier (so_shadowed)."""
    out = set()
    for vi in sorted({win.line_vcf[a.line] for a in win.by_key.get(al.key, ())
                      if win.line_sites_only[a.line]}):
        hap = make_hap(win.refseq, [al], lo, hi)
        sets = []
        for b, pl, src in zip(hap.seq, hap.place, hap.src):
            opts = {b}
            if src is None and pl[1] == 0 and not al.lo <= pl[0] < al.hi:
                opts.update(x.body for x in win.snv_at.get(pl[0], ())
                            if win.line_vcf[x.line] == vi and x not in win.shadowed
                            and x.filtered is None and x not in win.so_shadowed)
            sets.append(frozenset(opts))
        out |= feasible_cps(ctx, full, sets, hap.place)
    return out


def record_cps(win, ctx, guide, al, e=None):
    """(certain, open): the Cluster Positions of the rows CRISPRme writes for record
    ``al``, whatever the mode -- a genotyped dataset's are the targets its carriers'
    haplotypes hold (record_demands: the per-carrier resolution gives no row for a raw
    target no carrier holds, METHODS section 4; a slot whose alleles overlap the record
    gives its targets as METHODS reads it); a sites-only dataset's, the worst case of each
    alignment of its fake-indel window over that dataset's SNP alleles (sites_only_cps),
    both for a record both datasets hold -- and {Cluster Position: why} of the rows it may
    write or not, which no documented rule settles (record_demands' 'cps_open'). ``e``:
    record_demands' enumeration of the record, when completeness made it already (the
    value kept for the record holds both parts whoever computed it first). A record
    within one target of the contig's end has none (completeness lists it)."""
    key = (guide, al.key)
    got = win.cps.get(key)
    if got is not None:
        return got
    full = ctx.full(guide)
    one = len(full) + max(ctx.bdna, ctx.brna) + 2
    reach = 2 * one + 20
    out, opened = set(), defaultdict(set)
    if al.lo - one >= 0 and al.hi + one <= len(win.refseq):
        lo, hi = max(0, al.lo - reach), min(len(win.refseq), al.hi + reach)
        out |= sites_only_cps(win, ctx, full, al, lo, hi)
        if e is None:
            carriers = set()
            for a in win.by_key.get(al.key, ()):
                carriers.update(win.carriers(a))
            if carriers:
                w, x = win, al
                if lo < win.beg or hi > win.end:
                    w = Window(win.refseq, win.chrom, max(0, lo - 1000), hi + 1000,
                               win.vcfs, panel=win.panel, sem=win.sem)
                    x = w.by_key[al.key][0]
                e = record_demands(w, ctx, x, carriers, full, lo, hi,
                                   (x.edit_lo - one, x.edit_hi + one))
                if w is win:
                    win.enum[key] = e
        if e is not None:
            out |= e["cps"]
            for c, why in e["cps_open"].items():
                opened[c] |= why
    win.cps[key] = (frozenset(out), {c: frozenset(w) for c, w in opened.items()
                                     if c not in out})
    return win.cps[key]


# why CRISPRme may write a row or not, and the class a cluster that reaches a row as good as
# a target only through such rows is counted under: a slot whose overlap no documented rule
# settles; a target through an earlier record of a split site, or an untied multi-ALT alt,
# which its stores read otherwise; a per-sample run's capped window, whose greedy
# worst-case row only the classic dictionary path writes (METHODS section 5, known limit)
OPEN_CLASS = (("overlap", "limit-overlapping-records"),
              ("last-alt-wins", "limit-last-alt-wins"),
              ("putative-untied", "limit-putative-carrier"),
              ("capped", "limit-capped-window"),
              ("filtered", "limit-filtered-record"),
              ("degeneracy", "limit-max-degeneracy"),
              ("snp-row-indel", "coverage-snp-row-indel"))
OPEN_TEXT = {
    "overlap": "the rows of a slot whose alleles overlap, which no documented rule says "
               "CRISPRme writes or not (the one-record reading of a slot whose record meets "
               "a SNV at another position; a target an overlap zone of two other alleles "
               "of a slot changes, on any reading of the zone; a target of a slot read "
               "without its overlap that touches a call of unknown phase)",
    "last-alt-wins": "a target through an earlier record of a split SNV site, which "
                     "CRISPRme's stores read as the last",
    "putative-untied": "a target through an alt of a multi-ALT SNV line no genotype ties, "
                       "which CRISPRme resolves PUTATIVE",
    "capped": "a capped window's greedy worst-case row, which a per-sample run writes only "
              "on the classic dictionary path (an index built without --samplesID; METHODS "
              "section 5, known limit), and the run does not say which it took",
    "filtered": "a row of a sites-only dataset through one of its SNVs whose FILTER is "
                "not PASS or '.': the hg38_mega index holds such SNVs, an index CRISPRitz's "
                "enricher builds does not, and a run does not say which it read",
    "degeneracy": "a SNP row resting on more ambiguous columns than the run's "
                  "--max-degeneracy, which CRISPRme cuts back to one within it or leaves "
                  "out (METHODS section 5)",
    "snp-row-indel": "a per-sample SNP row of a haplotype whose slot also holds an indel "
                     "that touches the target: CRISPRme drops such a carrier unless its "
                     "reading of the slot through a repeat keeps the target's bases "
                     "(METHODS section 4, 2.10.0), which this checker's cluster chain does "
                     "not model (its own coverage gap)"}


def chain_cps(win, ctx, guide, full, a, b):
    """(rows, open) on [a, b): the Cluster Positions of the rows CRISPRme writes there --
    the reference and SNP rows (snp_row_cps) and the rows of every indel-path record that
    touches it (record_cps; a merge cluster holds the rows of every record, and the
    alleles of one multi-allelic site put their targets a few bases apart) -- and
    {Cluster Position: why} of the rows it may write or not (OPEN_CLASS)."""
    rows, capped, dropped, filtered, degenerate = snp_row_cps(win, ctx, guide, full, a, b)
    out = set(rows)
    opened = defaultdict(set)
    for c in degenerate:
        opened[c].add("degeneracy")
    for c in capped:
        opened[c].add("capped")
    for c in dropped:
        opened[c].add("snp-row-indel")
    for c in filtered:
        opened[c].add("filtered")
    seen = set()
    for lid in win.lines_in(a, b):
        for x in win.line_alleles[lid]:
            if x is not None and x.kind != "snv" and x.key not in seen:
                seen.add(x.key)
                certain, opn = record_cps(win, ctx, guide, win.by_key[x.key][0])
                out |= certain
                for c, why in opn.items():
                    opened[c] |= why
    return out, {c: w for c, w in opened.items() if c not in out}


def feasible_cps(ctx, full, sets, places, dense=None, drops=None, sparse=None):
    """The Cluster Positions of the PAM sites of ``sets`` (both strands) where some
    alignment of the guide is in budget (with ``dense``, only alignments over more than
    that many positions of several options: a capped window; with ``sparse``, only those
    over at most that many: a row --max-degeneracy cannot cut). With ``drops`` (alleles),
    (those of alignments whose target no allele of ``drops`` touches, those of the rest)."""
    n = len(full)
    pam_len, p5 = ctx.pam_len, ctx.pam_5prime
    motif = full[:pam_len] if p5 else full[n - pam_len:]
    proto = n - pam_len
    matches = [_MATCHES.get(c.upper(), frozenset()) for c in full]
    pam_first, pam_last = (0, pam_len - 1) if p5 else (n - pam_len, n - 1)
    motif_sets = [_MATCHES.get(c.upper(), frozenset()) for c in motif]
    size = len(sets)
    rc_sets = [frozenset(comp(x) for x in s) for s in reversed(sets)]
    out, touched = set(), set()
    for strand, ss in (("+", sets), ("-", rc_sets)):
        many = [0]
        for x in ss:
            many.append(many[-1] + (len(x) > 1))
        for q in range(0, size - pam_len + 1):
            if any(ss[q + j].isdisjoint(motif_sets[j]) for j in range(pam_len)):
                continue
            for width in range(proto - ctx.brna, proto + ctx.bdna + 1):
                beg = (q - width) if not p5 else q
                end = (q + pam_len) if not p5 else (q + pam_len + width)
                if beg < 0 or end > size:
                    continue
                if dense is not None and many[end] - many[beg] <= dense:
                    continue
                if sparse is not None and many[end] - many[beg] > sparse:
                    continue
                piece = ss[beg:end]
                if any("N" in x for x in piece):
                    continue
                if not _feasible_sets(matches, piece, n, end - beg, ctx.mm, ctx.bdna,
                                      ctx.brna, ctx.cap, pam_first, pam_last, p5):
                    continue
                off = beg if strand == "+" else size - end
                cp = cluster_of(places[off][0], strand, end - beg, n, p5)
                if drops and any(touches(x, places[off][0], places[off + end - beg - 1][0])
                                 for x in drops):
                    touched.add(cp)
                else:
                    out.add(cp)
    if drops is not None:
        return out, touched - out
    return out


def chain_from(c0, cps, merge):
    """The run of Cluster Positions, each within ``merge`` of the next, that holds c0."""
    pts = sorted(set(cps) | {c0})
    i = pts.index(c0)
    lo = hi = i
    while lo > 0 and pts[lo] - pts[lo - 1] <= merge:
        lo -= 1
    while hi + 1 < len(pts) and pts[hi + 1] - pts[hi] <= merge:
        hi += 1
    return pts[lo:hi + 1]


def cluster_chain(win, ctx, guide, c0, merge, wide=False, step=150, most=3000):
    """The Cluster Positions CRISPRme's merge can chain with c0: those of the rows it
    writes (chain_cps: the reference and SNP rows the mode allows, and each indel-path
    record's targets some carrier holds) -- and, ``wide``, those of the rows it may write
    or not --, each within ``merge`` bp of the next, from c0 on (CRISPRme clusters a
    guide's rows on a chromosome so, whatever their strand). The window grows while the
    chain reaches its edge (up to ``most`` bp either side). Returns (members, cut, why):
    ``cut`` when the chain still reached the edge at that limit, so the cluster may hold
    rows past what was followed; ``why``, the reasons (OPEN_CLASS) of the open rows among
    the members."""
    full = ctx.full(guide)
    lo, hi = c0 - step, c0 + step
    grown = False
    while True:
        a = max(0, lo - len(full) - 10)
        b = min(len(win.refseq), hi + len(full) + 10)
        if not grown and (win.needed is not None or a < win.beg or b > win.end):
            # the chain reaches past what the window read (a tile's reads a pad around
            # the tile): every record it can still reach is read, once
            grown = True
            win = Window(win.refseq, win.chrom,
                         min(win.beg, max(0, c0 - most - step - len(full) - 10)),
                         max(win.end, c0 + most + step + len(full) + 10), win.vcfs,
                         panel=win.panel, sem=win.sem)
        rows, opened = chain_cps(win, ctx, guide, full, a, b)
        cps = (rows | set(opened)) if wide else rows
        members = chain_from(c0, cps, merge)
        inside = members[0] - merge > lo and members[-1] + merge < hi
        if inside or hi - lo >= 2 * most:
            why = set()
            if wide:
                for c in members:
                    why |= opened.get(c, set())
            return members, not inside, why
        lo -= step
        hi += step


class Entry(object):
    """One reported row block, for completeness: its carriers, Cluster Position, + reading,
    edits under each sorting criterion and CFD -- recomputed from its own aligned strings
    when the CFD tables apply (``printed`` keeps the table's value)."""

    __slots__ = ("row", "carriers", "cp", "block", "crit", "cfd", "printed", "probe")

    def __init__(self, row, carriers, cp, probe=None, cfd_model=None, pam_len=3):
        self.row = row
        self.carriers = carriers
        self.cp = cp
        self.probe = probe
        self.block = row["block"]
        self.crit = {}
        try:
            mm, bul = int(row["mm"]), int(row["bulges"])
            self.crit = edits(mm, bul)
        except (KeyError, ValueError):
            pass
        try:
            self.printed = float(row.get("cfd", ""))
        except ValueError:
            self.printed = None
        calc = row_cfd(cfd_model, row.get("crrna", ""), row.get("alt", ""), pam_len)
        self.cfd = calc if calc is not None else self.printed


def cfd3(x):
    """A CFD as the analysis writes it ("{:.3f}"), the value the merge sorts on."""
    return float("%.3f" % x)


def cfd_below(a, b):
    """Whether CFD ``a`` (a row's) sorts below CFD ``b`` (a target's) in the merge, which
    sorts on three decimals. A difference of float rounding alone is none: a value within
    1e-9 of a rounding boundary (x.xxx5) may print either way, as products taken in
    another order do, so the row's is rounded up there and the target's down."""
    return cfd3(a + 1e-9) < cfd3(b - 1e-9) and b - a > 1e-9


def cfd_same(a, b):
    """Whether CFDs ``a`` and ``b`` surely print the same three decimals (neither within
    1e-9 of a rounding boundary): only then does the merge's CFD file go on to its next
    sort key for them."""
    return len({cfd3(a - 1e-9), cfd3(a + 1e-9), cfd3(b - 1e-9), cfd3(b + 1e-9)}) == 1


def crit_key(crits, mm, bulges):
    """An alignment's key in the merge's sort on the sorting criteria ``crits``."""
    e = edits(mm, bulges)
    return tuple(e[c] for c in crits)


def crit_text(crits, key):
    return ", ".join("%d %s" % (v, c) for v, c in zip(key, crits))


def row_first_k(win, row):
    """The rank of a reported row's first base among the bases its insertion adds (0 when
    it is any other base), read from the row's own target on the reference with its named
    variants: its cluster position on the haplotype is its Cluster Position plus that.
    Two places hold the same target only in a repeat it lies wholly in; the first is
    taken, as CRISPRme takes it (indel_snp_cis.row_first_base, on the fix branch). None
    when the row's target does not read at its own start through its record (a row
    soundness fails, as 2.8.1's rows placed past an insertion's anchor): such a row stands
    for no other alignment."""
    named = []
    for t in row["var"].split(","):
        try:
            _c, p, r, a = parse_token(t)
        except (ValueError, IndexError):
            return None
        found = win.by_key.get((p, r, a))
        if not found:
            return None
        named.append(found[0])
    mains = [a for a in named if a.kind != "snv"]
    if len(mains) != 1:
        return None
    main = mains[0]
    stripped = row["alt"].replace("-", "").upper()
    probe = stripped if row["strand"] == "+" else rc(stripped)
    n = len(probe)
    width = main.hi - main.lo
    lo = max(0, min(row["start"], main.lo) - n - width - 60)
    hi = min(len(win.refseq), max(row["start"], main.hi) + 2 * n + width + 60)
    hap = make_hap(win.refseq, named, lo, hi)
    if hap.clash:
        return None
    hits = reads_at(hap, row["start"], probe, main)
    return hap.place[hits[0]][1] if hits else None


class Reported(object):
    """What the tables report, for completeness: every row block that names a variant, by
    guide -- its carriers, its Cluster Position, its block's edits and CFD -- and every
    integrated line's clusters, edits and CFD (runs with one row per cluster)."""

    def __init__(self, rows, ctx, loci=None, merge=3, criterion="mm+bulges,mm",
                 scoring=None):
        self.exact = set()
        self.exact_any = set()           # (guide, sample, strand, start, + reading)
        self.by_var = defaultdict(list)  # (guide, variants) -> [Entry]
        self.by_target = defaultdict(list)  # (guide, strand, start, + reading, variants)
        self.at = set()                  # (guide, sample, cluster, record token)
        self.on_record = defaultdict(set)  # (guide, record token) -> samples listed on it
        self.ctx = ctx
        self.merge = merge
        # the sorting criteria, in order: the merge sorts the fewest file on them and the
        # CFD file on the CFD, then on them (merge_close_targets_cfd.sh hands both files
        # --sorting-criteria); the first alone is compared where the reference row and
        # the variant row of one cluster are chosen between (on mm+bulges or CFD only)
        self.criteria = tuple(c for c in criterion.split(",")
                              if c in ("mm+bulges", "mm", "bulges")) or ("mm+bulges", "mm")
        self.criterion = self.criteria[0]
        # --sorting-criteria-scoring: a CRISPRme whose merge gives it to the CFD file (a
        # later branch's does; 2.8.x's never reads it) breaks a CFD tie on it instead, so
        # a tie no worse under either is taken as no worse
        self.scoring = tuple(c for c in (scoring or "").split(",")
                             if c in ("mm+bulges", "mm", "bulges")) or None
        self._k = {}
        self._groups = {}   # (guide, variants) -> {(strand, carriers, site): [Entry]}
        self.loci = {}
        for gk, lines in (loci or {}).items():
            pos = sorted((c, k) for k, ln in enumerate(lines) for c in ln[1])
            self.loci[gk] = (pos, lines)
        for r in rows:
            t = r["alt"].replace("-", "").upper()
            probe = t if r["strand"] == "+" else rc(t)
            carriers = frozenset(samples_of(r))
            var = frozenset(x for x in r["var"].split(",") if x)
            g = r.get("guide", "")
            cl = cluster_of(r["start"], r["strand"], len(probe), len(r["spacer_col"]),
                            ctx.pam_5prime)
            for s in carriers:
                self.exact.add((g, s, r["strand"], r["start"], probe, var))
                self.exact_any.add((g, s, r["strand"], r["start"], probe))
                for v in var:
                    if is_indel_token(v):
                        self.at.add((g, s, cl, v))
                        self.on_record[(g, v)].add(s)
            e = Entry(r, carriers, cl, probe, ctx.cfd, ctx.pam_len)
            self.by_var[(g, var)].append(e)
            self.by_target[(g, r["strand"], r["start"], probe, var)].append(e)

    def has(self, guide, s, strand, start, probe, var=None):
        """Whether the sample is listed for the target with the variants ``var`` (any
        variants: None)."""
        if var is None:
            return (guide, s, strand, start, probe) in self.exact_any
        return (guide, s, strand, start, probe, var) in self.exact

    def first_k(self, win, row):
        key = (row["table"], row["line"], row["block"])
        if key not in self._k:
            self._k[key] = row_first_k(win, row)
        return self._k[key]

    def group(self, guide, var, carriers, site, strand, first_k):
        """The rows that stand for one target's alignments after the merge's fold: those
        naming ``var`` with the carrier set ``carriers`` on ``strand`` whose target starts
        at haplotype cluster position ``site`` (a row whose target is not at its own start
        belongs to none). CRISPRme's fold key holds no strand, but the alignments of one
        target share it: a row of the other strand is another target, which a fold would
        lose, so it stands for none of this one's."""
        groups = self._groups.get((guide, var))
        if groups is None:
            groups = self._groups[(guide, var)] = defaultdict(list)
            for e in self.by_var.get((guide, var), ()):
                k = first_k(e.row)
                if k is not None:
                    groups[(e.row["strand"], e.carriers, e.cp + k)].append(e)
        return groups.get((strand, carriers, site), [])

    def exact_for(self, guide, sample, t, toks, best, first_k):
        """(listed, why) of target ``t`` (its strand, start, + reading and the variants
        ``toks``) for ``sample``: whether a row block reports it exactly for the sample,
        and why the fold groups of those blocks -- every row of its site naming those
        variants with one block's carrier set, which the merge folds into the best by each
        block's order (METHODS section 4: one entry per carrier set) -- are worse than
        ``best``, taken together (None: no worse). A PUTATIVE co-carrier listed on some of
        the target's alignments only splits them between two carrier sets: the alignment
        with the fewest edits can stand in one group and the one of the highest CFD in the
        other, each group's block being its best there, so the sample's best of each block
        is the best of its groups' (judged one group at a time, both alignments reported
        for the sample failed, missing-best-alignment)."""
        rows = []
        listed = False
        for e in self.by_target.get((guide, t.strand, t.start, t.probe, toks), ()):
            if sample not in e.carriers:
                continue
            listed = True
            k = first_k(e.row)
            if k is not None:
                rows.extend(self.group(guide, toks, e.carriers, e.cp + k, t.strand, first_k))
        if not listed:
            return False, None
        if not rows:
            return True, "its row's target is not at its own start"
        return True, self.worse(list({id(e): e for e in rows}.values()), best)

    def folded(self, guide, sample, hcl, strand, var, required, optional, d, first_k):
        """(ok, why): rows of the same site naming the same variants whose carriers are the
        target's own -- every carrier of ``required`` and nobody outside it and
        ``optional`` (documented_carriers) -- and no worse than it in each block: the merge
        folds those of one carrier set into the first of them by each file's order
        (_fold_same_carriers), which worse() compares on the merge's sort keys, and the
        groups of several such carrier sets (a PUTATIVE co-carrier listed on some of the
        target's alignments only) are judged together, as in exact_for. The site is the
        cluster position on the haplotype carrying the record, ``hcl``: the Cluster
        Position plus the rank of the first base among an insertion's added bases
        (first_k of a row), which is the Cluster Position itself on any row that does not
        start inside an insertion. (On d7a6fff the merge compares rows at their Cluster
        Position alone, so it folds two targets that start on different bases of one
        insertion; the fix branch's c3a75fe compares them on the haplotype.) The rows are
        on the target's ``strand``: CRISPRme's fold key holds none, and a fold across
        strands loses another target."""
        allowed = set(required) | set(optional)
        rows = []
        for e in self.by_var.get((guide, var), ()):
            if e.cp > hcl or sample not in e.carriers or e.row["strand"] != strand:
                continue
            if not (required <= e.carriers <= allowed):
                continue
            k = first_k(e.row)
            if k is None or e.cp + k != hcl:
                continue   # another site, or a row whose target is not at its start
            rows.append(e)
        if not rows:
            return False, None
        why = self.worse(rows, d)
        return why is None, why

    def worse(self, es, d):
        """None when one fold group's rows ``es`` are no worse than ``d`` (a Best) by the
        merge's sort keys in each block, else why. The fewest file sorts on the sorting
        criteria in order (mm+bulges, then mm, by default); the CFD file on the CFD at
        the three decimals the analysis writes, then on the same criteria. So the
        fewest-mm+b block's key must be at most the least key of the target's alignments,
        and the highest-CFD block's CFD at least the target's best -- and, where the two
        surely print the same three decimals, its key at most the least key of the
        target's alignments of that CFD. A group the run lists has a row in each block
        (each file keeps one per group). A row's CFD is recomputed from its aligned
        strings. Not compared: the keys after the criteria (the Cluster Position, the
        carrier count, the aligned text), which only order alignments equal in all of
        these."""
        crits = self.criteria
        few_rows = [e for e in es if e.block == "fewest_mm+b"]
        cfd_rows = [e for e in es if e.block == "highest_CFD"]
        if not few_rows or not cfd_rows:
            return "no %s block of its site lists its carriers" % (
                "fewest-mm+b" if not few_rows else "highest-CFD")
        keys = sorted(d.keys, key=lambda k: (k[0], k[1], -1.0 if k[2] is None else k[2]))
        few = [tuple(e.crit[c] for c in crits) for e in few_rows
               if all(c in e.crit for c in crits)]
        if few and keys:
            rk = min(few)
            tk = min(crit_key(crits, mm, b) for mm, b, _c in keys)
            if rk > tk:
                if rk[0] > tk[0]:
                    return ("the row of its site with its carriers has a fewest-mm+b block "
                            "of %d (%s), the target %d" % (rk[0], crits[0], tk[0]))
                return ("the row of its site with its carriers has a fewest-mm+b block of "
                        "%s, the target %s (worse by the merge's next sorting criterion, "
                        "--sorting-criteria %s)"
                        % (crit_text(crits, rk), crit_text(crits, tk), ",".join(crits)))
        cfd = [e for e in cfd_rows if e.cfd is not None]   # none: no CFD column
        if cfd and d.cfd is not None:
            top = max(cfd, key=lambda e: e.cfd)
            if cfd_below(top.cfd, d.cfd):
                return ("the row of its site with its carriers has a highest-CFD block of "
                        "%.3f (%.6f, printed %s), the target %.3f (%.6f)"
                        % (cfd3(top.cfd), top.cfd, top.row.get("cfd", "") or "nothing",
                           cfd3(d.cfd), d.cfd))
            why = None
            for crs in [crits] + ([self.scoring] if self.scoring and
                                  self.scoring != crits else []):
                tied = [crit_key(crs, mm, b) for mm, b, c in keys
                        if c is not None and cfd_same(c, d.cfd)]
                if not tied:
                    return None
                tk = min(tied)
                ok = False
                for e in cfd:
                    if cfd_below(e.cfd, d.cfd):
                        continue
                    if (not cfd_same(e.cfd, d.cfd) or not all(c in e.crit for c in crs)
                            or tuple(e.crit[c] for c in crs) <= tk):
                        ok = True   # a higher (or not surely equal) CFD, or no worse key
                        break
                if ok:
                    return None
                rk = min(tuple(e.crit[c] for c in crs) for e in cfd
                         if cfd_same(e.cfd, d.cfd))
                why = why or ("the row of its site with its carriers has a highest-CFD block "
                              "of %.3f with %s, the target %.3f with %s (a tie at three "
                              "decimals, which the merge breaks on --sorting-criteria %s)"
                              % (cfd3(top.cfd), crit_text(crs, rk), cfd3(d.cfd),
                                 crit_text(crs, tk), ",".join(crs)))
            return why
        return None

    def on_rows(self, guide, samples, tok):
        """How many of ``samples`` some row of the tables lists with the record ``tok``:
        a limit line says it, so that a row that lists only such carriers does not leave
        the findings as they are when it is lost."""
        got = self.on_record.get((guide, tok), ())
        return sum(1 for s in samples if s in got)

    def listed_at(self, guide, sample, cluster, tok):
        return (guide, sample, cluster, tok) in self.at

    def near(self, pos, cps):
        out = set()
        for c in cps:
            i = bisect.bisect_left(pos, (c - self.merge, -1))
            while i < len(pos) and pos[i][0] <= c + self.merge:
                out.add(pos[i][1])
                i += 1
        return out

    def locus(self, guide, chrom, c, d, chain):
        """(ok, why, cut, via) for a target of Cluster Position c in a run with one row per
        merge cluster. The cluster's row is an integrated line with a block within
        --merge bp of c, or of a Cluster Position chained to c (chain(wide): those of the
        rows CRISPRme writes, each within --merge bp of the next); none is a miss. It
        must have a fewest-mm+b block no worse than the target by the first sorting
        criterion and a highest-CFD block of at least its CFD at the three decimals the
        merge sorts on (the row's CFD recomputed from its aligned strings). ``cut``: the
        chain stopped at its limit (cluster_chain). ``via``: the reasons (OPEN_CLASS) of
        the rows CRISPRme may write or not that the cluster reaches such a row only
        through (empty when it needs none) -- the verdict then rests on what no
        documented rule settles."""
        pos, lines = self.loci.get((guide, chrom), ([], []))
        near = self.near(pos, [c])
        if near:
            ok, why = self.good(near, lines, d)
            return ok, why, False, set()
        members, cut, _w = chain(False)
        near = self.near(pos, members)
        ok, why = self.good(near, lines, d) if near else (False, None)
        if ok:
            return True, "", cut, set()
        if why is None:
            why = ("no integrated row lies in its cluster (Cluster Positions %d..%d: the "
                   "rows CRISPRme writes there -- reference, SNP and the targets some "
                   "carrier of an indel-path record holds -- that chain with it within %d "
                   "bp)" % (members[0], members[-1], self.merge))
        wider, _cut, via = chain(True)
        far = self.near(pos, wider)
        if via and far and self.good(far, lines, d)[0]:
            return True, why, cut, via
        return False, why, cut, set()

    def good(self, near, lines, d):
        """(ok, why): whether one of the integrated lines ``near`` is as good as the
        target in both blocks (locus)."""
        why = []
        crit = self.criterion
        for k in sorted(near):
            _g, _clusters, r_crit, r_cfd, ln = lines[k]
            v = r_crit.get(crit) if r_crit else None
            if v is not None and v > d.best[crit]:
                why.append("the row at line %d has a fewest-mm+b block of %d (%s)"
                           % (ln, v, crit))
                continue
            if d.cfd is not None and r_cfd is not None and cfd_below(r_cfd, d.cfd):
                why.append("the row at line %d has a highest-CFD block of %.3f (%.6f), "
                           "the target %.3f" % (ln, cfd3(r_cfd), r_cfd, cfd3(d.cfd)))
                continue
            return True, ""
        return False, "; ".join(why)

# --------------------------------------------------------------------------- #
# driver
# --------------------------------------------------------------------------- #
def read_panel(paths):
    """The sample IDs of samplesID files (the first column of each line; '#' lines and a
    SAMPLE_ID header are skipped)."""
    out = set()
    for path in paths:
        op = gzip.open if path.endswith(".gz") else open
        with op(path, "rt") as fh:
            for line in fh:
                s = line.strip()
                if not s or s.startswith("#"):
                    continue
                sid = s.split()[0]
                if sid != "SAMPLE_ID":
                    out.add(sid)
    return out


def read_pam(path):
    """(motif, length, 5-prime?) from a CRISPRme PAM file, e.g. 'NNN...NRG 3'."""
    with open(path) as fh:
        text = fh.read().split()
    motif, k = text[0].upper(), int(text[1])
    if k < 0:
        return motif[:-k], -k, True
    return motif[-k:], k, False


def read_guides(values, pam_len, pam_5prime):
    """The guides (without PAM) from --guide values: sequences or files of one guide per
    line, as CRISPRme reads them (a guide may carry N's over the PAM)."""
    out = []
    for v in values:
        items = [v]
        if os.path.isfile(v):
            with open(v) as fh:
                items = [line.strip() for line in fh if line.strip()]
        for g in items:
            g = g.upper()
            if pam_5prime and g.startswith("N" * pam_len):
                g = g[pam_len:]
            elif not pam_5prime and g.endswith("N" * pam_len):
                g = g[:len(g) - pam_len]
            if set(g) - set("ACGTU"):
                raise ValueError("not a guide: %s" % g)
            out.append(g.replace("U", "T"))
    return out


def load_known(path):
    out = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[0] == "class":
                continue
            f += [""] * (7 - len(f))
            out.append({"class": f[0], "chrom": f[1], "start": f[2], "strand": f[3],
                        "variants": f[4], "samples": f[5], "reason": f[6]})
    return out


def apply_known(rep, known, chroms, modes):
    used = set()
    for r in rep.rows:
        for i, k in enumerate(known):
            if (r["class"] == k["class"] and r["chrom"] == k["chrom"]
                    and str(r["start"]) == k["start"] and r["strand"] == k["strand"]
                    and r["variants"] == k["variants"]
                    and (not k["samples"] or r["samples"] == k["samples"])):
                used.add(i)
                r["detail"] = "%s: %s | known: %s" % (r["class"], r["detail"], k["reason"])
                r["class"] = "known-defect"
                r["gating"] = "no"
                break
    for i, k in enumerate(known):
        if i in used or k["chrom"] not in chroms:
            continue
        mode = "completeness" if k["class"].startswith("missing") else "soundness"
        if mode not in modes:
            continue
        rep.rows.append({"gating": "no", "class": "known-not-seen", "table": "",
                         "line": "", "block": "", "chrom": k["chrom"], "start": k["start"],
                         "strand": k["strand"], "spacer": "", "target": "",
                         "variants": k["variants"], "samples": k["samples"],
                         "detail": "listed in --known (%s) and not found: fixed, or "
                                   "outside this run?" % k["reason"]})


G = {}  # what the work units share (set before the processes fork)


def _blockidx_of(path):
    v = Vcf(path)
    v.build_blockidx()
    return v.blockidx


def row_span(r, pad):
    """The reference stretch a soundness check of row ``r`` reads records over: its start
    and target and the REF span of every record it names, ``pad`` either side (a long
    deletion's far junction, and the carriers' alleles there, lie past the start by the
    deletion's length)."""
    lo = r["start"]
    hi = lo + len(r["alt"].replace("-", ""))
    for t in r["var"].split(","):
        try:
            _c, p, ref, _a = parse_token(t)
        except (ValueError, IndexError):
            continue
        lo = min(lo, p - 1)
        hi = max(hi, p - 1 + len(ref))
    return max(0, lo - pad), hi + pad


def run_unit(unit):
    """One work unit: ('sound' | 'comp', chrom, lo, hi, rows). Returns (the findings, the
    unit's counts: {'sites_only': records only sites-only VCFs hold, 'filtered': records
    the search skips for their FILTER})."""
    kind, ch, tlo, thi, rows = unit
    ctx = G["ctx"]
    rep = Report()
    stats = {"chrom": ch}
    refseq = G["ref"].seq(ch)
    pad = 2 * (len(ctx.guides[0]) + ctx.pam_len) + 1000
    vcfs = G["vcfs_for"][ch]
    if kind == "sound":
        needed = []
        for a, b in sorted(row_span(r, pad) for r in rows):
            if needed and a <= needed[-1][1]:
                needed[-1] = (needed[-1][0], max(needed[-1][1], b))
            else:
                needed.append((a, b))
        win = Window(refseq, ch, needed[0][0], needed[-1][1], vcfs, needed=needed,
                     panel=ctx.panel, sem=ctx.sem)
        for r in sorted(rows, key=lambda x: x["start"]):
            if any(is_indel_token(x) for x in r["var"].split(",")):
                check_row(r, win, ctx, rep)
            else:
                check_snp_row(r, win, ctx, rep, G.get("snp_listed", {}))
    else:
        win = Window(refseq, ch, max(0, tlo - pad), thi + pad, vcfs, panel=ctx.panel,
                     sem=ctx.sem)
        # a record of the tile whose REF runs far (a long deletion) has its targets'
        # haplotypes written past what was read: read again, to its far end
        one = max(len(ctx.full(g)) for g in ctx.guides) + max(ctx.bdna, ctx.brna) + 2
        far = max([a.hi for a in win.alleles if a.kind != "snv" and tlo <= a.lo < thi]
                  or [0]) + 2 * one + 20   # the far end of what comp_record writes
        if far > win.end:
            win = Window(refseq, ch, max(0, tlo - pad), far + pad, vcfs, panel=ctx.panel,
                         sem=ctx.sem)
        completeness(win, ctx, rep, G["reported"][ch], (tlo, thi), G["popmode"], stats)
    for v in vcfs:
        if v.long_ref:
            stats.setdefault("long_ref", {})[v.path] = v.long_ref
    stats["mixed"] = sorted(win.mixed)
    return rep.rows, stats


def default_cfd_dir():
    here = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(os.path.dirname(here), "PostProcess")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--results", required=True, help="a CRISPRme result folder")
    p.add_argument("--reference", required=True,
                   help="the reference genome: a folder of per-chromosome FASTAs or one "
                        "FASTA (any number of contigs)")
    p.add_argument("--vcf", action="append", default=[], help="a bgzipped VCF (repeatable)")
    p.add_argument("--vcf-dir", action="append", default=[], help="a folder of VCFs")
    p.add_argument("--pam", required=True, help="the PAM file the search used")
    p.add_argument("--guide", action="append", required=True,
                   help="the guide (or the guide file the search read); repeatable")
    p.add_argument("--mm", type=int, required=True)
    p.add_argument("--bDNA", dest="bdna", type=int, required=True)
    p.add_argument("--bRNA", dest="brna", type=int, required=True)
    p.add_argument("--max-total-edits", type=int, default=None,
                   help="the cap on mismatches + bulges the search used (default: the "
                        "run's own, from its .Params.txt or .command_line.txt, else 6, "
                        "complete-search's default; a negative value: no cap)")
    p.add_argument("--samples", action="append", default=[],
                   help="a samplesID file the search read (repeatable: the panel is their "
                        "union; the sample ID is the first column, '#' lines are skipped). "
                        "Default: the run's .sampleID.txt (CRISPRme writes there the "
                        "samples of those files), else every sample of the VCF(s)")
    p.add_argument("--chrom", default=None)
    p.add_argument("--region", default=None, help="chrom:start-end (1-based, inclusive)")
    p.add_argument("--mode", default="both", choices=("soundness", "completeness", "both"))
    p.add_argument("--analysis", default="auto",
                   choices=("auto", "per-sample", "population-level"),
                   help="the run's analysis mode (auto: read from the result folder; "
                        "given, it must be the run's own where its command line says)")
    p.add_argument("--snp-rows", default="off", choices=("off", "report"),
                   help="also check SNP-only rows of a per-sample run, report-only")
    p.add_argument("--known", default=None,
                   help="a TSV of findings already reported upstream (columns: class, "
                        "chrom, start, strand, variants, samples (empty: any), reason): "
                        "a matching finding turns non-gating as known-defect, and an "
                        "entry no longer found is reported as known-not-seen")
    p.add_argument("--cfd-tables", default=None,
                   help="the folder of cfd_mismatch_scores.json and cfd_pam_scores.json "
                        "(default: the repository's PostProcess)")
    p.add_argument("--jobs", type=int, default=1, help="processes to split the work over")
    p.add_argument("--out", default=None, help="the TSV of findings")
    p.add_argument("--tile", type=int, default=50000,
                   help="how many bases of a chromosome one work unit covers")
    p.add_argument("--progress", type=int, default=0, help="a line every N work units")
    args = p.parse_args(argv)

    pam_motif, pam_len, pam_5prime = read_pam(args.pam)
    try:
        guides = read_guides(args.guide, pam_len, pam_5prime)
    except ValueError as e:
        sys.stderr.write("ERROR: %s\n" % e)
        return 2
    vcf_paths = list(args.vcf)
    for d in args.vcf_dir:
        vcf_paths.extend(sorted(glob.glob(os.path.join(d, "*.vcf.gz"))))
    seen_paths = set()
    vcf_paths = [v for v in vcf_paths if not (os.path.realpath(v) in seen_paths
                                              or seen_paths.add(os.path.realpath(v)))]
    if not vcf_paths:
        sys.stderr.write("ERROR: no VCF given (--vcf / --vcf-dir)\n")
        return 2
    chrom = args.chrom
    lo = hi = None
    if args.region:
        m = re.match(r"^([^:]+):(\d+)-(\d+)$", args.region)
        if not m:
            sys.stderr.write("ERROR: --region wants chrom:start-end\n")
            return 2
        chrom, lo, hi = m.group(1), int(m.group(2)) - 1, int(m.group(3))
    tables = result_tables(args.results)
    if not tables:
        sys.stderr.write("ERROR: no results table under %s\n" % args.results)
        return 2
    try:
        (mode_run, alt_written, merge, criterion, scoring, run_cap, mode_from,
         flagless) = run_kind(args.results)
    except ValueError as e:
        # a run's own markers that no CRISPRme run writes together: which one holds is
        # unknown, so the run would be modeled on a guess
        sys.stderr.write("ERROR: the run's own records contradict each other (%s): the run "
                         "is not checked\n" % e)
        return 2
    # the check's inputs must be the search's: where the run's own record of them says
    # otherwise (its .Params.txt, .guides.txt, .sampleID.txt, .search_mode), the run would
    # be modeled partially, so it is refused
    wrong = run_mismatch(args.results, args, pam_motif, pam_5prime, guides)
    if args.analysis != "auto" and mode_run is not None and args.analysis != mode_run:
        wrong.append("--analysis %s, while the run's %s says %s"
                     % (args.analysis, mode_from, mode_run))
    if wrong:
        sys.stderr.write("ERROR: the check was given other inputs than the search used (%s): "
                         "give the search's own (budget, PAM, guides, samples list, mode); "
                         "the run is not checked\n" % "; ".join(wrong))
        return 2
    if args.analysis != "auto" and mode_run is None:
        # the run records no mode: the one given is taken, and with it the mode's default
        # for the alternative alignments where the command line names no flag for them
        mode_run, mode_from = args.analysis, "--analysis"
        if alt_written is None and flagless:
            alt_written = mode_run == "per-sample"
    # what the run records must be enough to model it: completeness needs whether it
    # wrote alternative alignments (one row per merge cluster without them), and then the
    # mode, whose SNP rows the cluster chains through; the SNP-row report needs the mode
    # (a population-level run's SNP rows are a worst-case union). --ci-cd-test says
    # neither: complete-test's searches take complete-search's default mode
    short = []
    if args.mode in ("completeness", "both"):
        if alt_written is None:
            short.append("completeness needs whether the run wrote alternative alignments "
                         "(its .alt_mode; else its command line's --alt-alignments, "
                         "--no-alt-alignments or mode), and the run records none")
        elif alt_written is False and mode_run is None:
            short.append("completeness of a run without alternative alignments (one row per "
                         "merge cluster) needs the run's analysis mode, whose SNP rows the "
                         "cluster chains through (its .search_mode; else its command line's "
                         "--per-sample or --population-level; --ci-cd-test chooses none), "
                         "and the run records none: give --analysis")
    if args.snp_rows == "report" and mode_run is None:
        short.append("--snp-rows report needs the run's analysis mode (its .search_mode; else "
                     "its command line's --per-sample or --population-level), and the run "
                     "records none: give --analysis")
    if short:
        sys.stderr.write("ERROR: %s; the run is not checked\n" % "; ".join(short))
        return 2
    # the total-edits cap: the one given, else the run's own, else complete-search's 6
    if args.max_total_edits is None:
        if run_cap is None:
            sys.stderr.write("check_indel_rows: the run does not say its --max-total-edits "
                             "(no .Params.txt): taking 6, complete-search's default\n")
            run_cap = 6
        args.max_total_edits = run_cap
    if args.max_total_edits < 0:
        args.max_total_edits = None   # no cap
    # one row per cluster (no alternative alignments): completeness is checked per cluster
    popmode = alt_written is False
    cfd = None
    try:
        cfd = Cfd(args.cfd_tables or default_cfd_dir())
    except (OSError, ValueError) as e:
        sys.stderr.write("check_indel_rows: no CFD tables (%s): the PUTATIVE choices CFD "
                         "decides are left open and the per-cluster check compares edits "
                         "only\n" % e)
    ctx = Ctx(args, pam_motif, pam_len, pam_5prime, guides, cfd)
    # a merge cluster's SNP rows: a per-sample run's are the haplotypes individuals hold,
    # a population-level run's the worst case of each alignment (only read per cluster,
    # where the mode is known)
    ctx.snp_union = mode_run != "per-sample"
    # a target of a record inside the region starts at most one target (and its bulges)
    # before the record and may hold a long REF: rows are read with that margin, so that
    # completeness sees the rows of the region's records; soundness checks the region's own
    full_len = len(ctx.full(guides[0]))
    margin = full_len + max(ctx.bdna, ctx.brna) + Vcf.LONGEST
    rows = []
    loci = {}
    bad = []
    novar = []
    for path, label in tables:
        rows.extend(read_rows(path, label, chrom, lo, hi, margin, pam_5prime,
                              loci=loci if label == "integrated" else None, bad=bad,
                              pam_len=pam_len, cfd_model=ctx.cfd, novar=novar))
    indel_rows = [r for r in rows if any(is_indel_token(t) for t in r["var"].split(","))]
    snp_rows = [r for r in rows if r["origin"] == "alt"
                and not any(is_indel_token(t) for t in r["var"].split(","))]
    sys.stderr.write("check_indel_rows: %d indel/MNV row block(s), %d SNP-only, %s run%s%s, "
                     "from %s\n" % (sum(1 for r in indel_rows if r["own"]),
                                    sum(1 for r in snp_rows if r["own"]),
                                    mode_run or "unknown-mode",
                                    " (its %s)" % mode_from if mode_from else "",
                                    " without alternative alignments" if popmode else
                                    "" if alt_written else
                                    " (whether it wrote alternative alignments unknown)"
                                    if alt_written is None else "",
                                    ", ".join(os.path.basename(p) for p, _l in tables)))
    if mode_run == "population-level" and args.snp_rows == "report":
        sys.stderr.write("check_indel_rows: a population-level run's SNP rows are a "
                         "worst-case union; --snp-rows report is ignored\n")
        args.snp_rows = "off"
    ref = Reference(args.reference)
    ref.keep_all = args.jobs > 1  # loaded once in this process, shared with the workers
    vcfs_all = [Vcf(p) for p in vcf_paths]
    no_genotypes = not any(v.samples for v in vcfs_all)
    logged = registry_only_logged(args.results) if not no_genotypes else None
    if logged:
        # a genotyped dataset searched without its genotype store, in any mode: the run's
        # own log says so (PRECONDITION), and it is refused
        sys.stderr.write(
            "ERROR: the run's log_verbose.txt says that a dataset it searched had neither a "
            "genotype store (Dictionaries/genotypes_<vcf>/) nor classic per-sample "
            "dictionaries ('%s%s: ...', not an aggregate registry's), as after download "
            "--what index --no-genotypes, and the VCF(s) given hold genotypes. Such a "
            "search lists every carrier of an indel on its SNP+indel rows, which this "
            "checker does not model (it models a search that read the genotypes the VCF "
            "holds: install the index with its genotype store and search again); the run "
            "is not checked\n" % (REGISTRY_ONLY, ", ".join(logged)))
        return 2
    if not no_genotypes and per_sample_downgraded(args.results):
        # --per-sample on a genotyped dataset without its genotype store: the run records it
        # (the PRECONDITION's one case a result folder tells apart), and it is refused
        sys.stderr.write(
            "ERROR: the run asked for --per-sample and ran population-level (its "
            ".search_mode), which complete-search does where every dataset it searches is "
            "sites-only or where a genotyped one has no genotype store "
            "(Dictionaries/genotypes_<vcf>/) and no classic per-sample dictionaries (an "
            "install made with download --what index --no-genotypes); the VCF(s) given hold "
            "genotypes, so it is the second. Such a run lists every carrier of an indel on "
            "its SNP+indel rows, which this checker does not model (it models a search that "
            "read the genotypes the VCF holds: install the index with its genotype store "
            "and search again); the run is not checked\n")
        return 2
    # the panel: the samplesID files given, else the run's .sampleID.txt, else every
    # sample of the VCF(s); CRISPRme lists no carrier outside it
    panel_from = args.samples
    if not panel_from and os.path.isfile(os.path.join(args.results, ".sampleID.txt")):
        panel_from = [os.path.join(args.results, ".sampleID.txt")]
    if panel_from and not no_genotypes:
        try:
            ctx.panel = read_panel(panel_from)
        except OSError as e:
            sys.stderr.write("ERROR: cannot read the samples list: %s\n" % e)
            return 2
        in_vcf = {s for v in vcfs_all for s in v.samples}
        if not ctx.panel & in_vcf:
            sys.stderr.write("ERROR: no sample of %s is in the VCF(s)\n"
                             % ", ".join(panel_from))
            return 2
        sys.stderr.write("check_indel_rows: the panel is the %d sample(s) of %s (%d of the "
                         "VCF's %d)\n" % (len(ctx.panel), ", ".join(panel_from),
                                          len(ctx.panel & in_vcf), len(in_vcf)))
    elif not no_genotypes:
        sys.stderr.write("check_indel_rows: no samples list (--samples, or the run's "
                         ".sampleID.txt): every sample of the VCF(s) is a panel sample\n")
    if no_genotypes:
        sys.stderr.write("check_indel_rows: the VCF(s) hold no genotypes (a sites-only "
                         "panel): a row's carriers are not checked, only the row itself, "
                         "and completeness has no haplotype to enumerate\n")
    elif any(not v.samples for v in vcfs_all):
        sys.stderr.write("check_indel_rows: %d of the %d VCF(s) hold no genotypes (a "
                         "sites-only dataset beside a genotyped one): a row that lists no "
                         "carrier for a record such a VCF holds is that dataset's, and its "
                         "strings, target and variants are checked, not its carriers\n"
                         % (sum(1 for v in vcfs_all if not v.samples), len(vcfs_all)))
    todo = [v for v in vcfs_all if v.tbi is None and v.bgzf]
    if todo:
        sys.stderr.write("check_indel_rows: %d VCF(s) without a .tbi: reading each once "
                         "into an in-memory block index\n" % len(todo))
        if args.jobs > 1 and len(todo) > 1:
            import multiprocessing as mp
            with mp.get_context("fork").Pool(min(args.jobs, len(todo))) as pool:
                for v, idx in zip(todo, pool.map(_blockidx_of, [v.path for v in todo])):
                    v.blockidx = idx
        else:
            for v in todo:
                v.build_blockidx()

    rep = Report()
    for table, line_no, nf, nh in bad:
        rep.add(True, "row-unreadable", {"table": table, "line": line_no},
                "the line has %d fields, the header %d: not read" % (nf, nh))
    for r in (novar if args.mode in ("soundness", "both") else ()):
        rep.add(True, "alt-without-variant", r,
                "REF/ALT_origin is alt and Variant_info_genome is %r: an ALT row is made "
                "by variants and must name them (the block is not checked further)"
                % r["var"])
    check_rows = [r for r in indel_rows + (snp_rows if args.snp_rows == "report" else [])
                  if r["own"]]
    # the rows of a guide the run searched and --guide does not give are another check's
    # (rows are keyed by guide): counted, one line per guide, not checked
    ran = run_guides(args.results, pam_motif, pam_5prime)
    if ran is not None:
        other = defaultdict(int)
        if args.mode in ("soundness", "both"):
            keep = []
            for r in check_rows:
                if r["guide"] in ran and r["guide"] not in guides:
                    other[r["guide"]] += 1
                else:
                    keep.append(r)
            check_rows = keep
        # one line per guide the run searched and --guide leaves out, in either mode:
        # completeness demands none of its targets either
        for g in sorted(set(ran) - set(guides)):
            parts = []
            if args.mode in ("soundness", "both"):
                parts.append("%d row block(s) of the guide %s, which the run searched and "
                             "--guide does not give, are not checked" % (other[g], g))
            if args.mode in ("completeness", "both"):
                parts.append("and completeness demands none of its targets" if parts else
                             "completeness demands none of the targets of the guide %s, "
                             "which the run searched and --guide does not give" % g)
            rep.add(False, "guide-not-given", {"table": "", "crrna": g},
                    "%s (give it to check them)" % ", ".join(parts))
    by_chrom = defaultdict(list)
    for r in check_rows:
        by_chrom[r["chrom"]].append(r)
    chroms = set(by_chrom) if args.mode in ("soundness", "both") else set()
    if args.mode in ("completeness", "both"):
        if chrom:
            chroms.add(chrom)
        else:
            # every chromosome a VCF holds (and the reference has), not only those with
            # rows: a chromosome whose indel rows were all lost is checked too
            prefix = any(r["chrom"].startswith("chr") for r in rows) or not rows
            for v in vcfs_all:
                for n in v.chroms():
                    name = n if (n.startswith("chr") or not prefix) else "chr" + n
                    if ref.has(name):
                        chroms.add(name)
            chroms |= set(by_chrom)
            if not chroms:
                sys.stderr.write("ERROR: completeness found no chromosome to check; give "
                                 "--chrom or --region\n")
                return 2
    vcfs_for = {}
    lengths = {}
    for ch in sorted(chroms):
        cov = [v for v in vcfs_all if v.covers(ch) is not False]
        vcfs_for[ch] = cov
        why = None
        if not cov:
            why = "no VCF given covers %s" % ch
        else:
            try:
                lengths[ch] = len(ref.seq(ch))
            except (KeyError, FileNotFoundError) as e:
                why = "no reference sequence for %s (%s)" % (ch, e)
        if why is None:
            continue
        n = 0
        for r in by_chrom.get(ch, ()):
            n += 1
            rep.add(True, "row-not-checked", r, why + ": this row is not checked")
        if args.mode in ("completeness", "both") and chrom == ch:
            rep.add(True, "row-not-checked", {"chrom": ch, "table": "completeness"},
                    why + ": completeness not run")
        sys.stderr.write("check_indel_rows: %s (%d row block(s) not checked)\n" % (why, n))
    # the VCFs given must be the search's: one per variant dataset it read (.Params.txt),
    # on every chromosome checked. With fewer, the rows of a dataset left out are read
    # against the others' panel (a run of two datasets that share samples, given one of
    # them, passes or fails by which); with more, a dataset the search never read is
    # demanded of it. Either way the run would be modeled partially: it is refused
    datasets, datasets_from = run_datasets(args.results)
    if datasets is None:
        sys.stderr.write("check_indel_rows: the run does not say which variant datasets it "
                         "searched (%s): the VCF(s) given are taken for the search's own\n"
                         % datasets_from)
    if datasets is not None:
        for ch in sorted(chroms):
            here = [v for v in vcfs_for.get(ch, ()) if v.holds(ch)]
            if not here or len(here) == len(datasets):
                continue   # (a chromosome no VCF holds is reported above, row-not-checked)
            sys.stderr.write(
                "ERROR: the run searched %d variant dataset(s) (%s: %s), and %d of the "
                "VCF(s) given hold %s (%s). Give the VCF "
                "of %s of every dataset the search read, one per dataset, and no other: "
                "%s, so the run would be modeled partially, and it is not checked (if a "
                "dataset of the search has no VCF for %s, check the other chromosomes with "
                "--chrom).\n"
                % (len(datasets), datasets_from, ", ".join(datasets) or "the reference alone",
                   len(here), ch, ", ".join(v.path for v in here), ch,
                   "with fewer, the rows of a dataset left out are read against the "
                   "others' panel" if len(here) < len(datasets) else
                   "with more, a dataset the search never read is demanded of it", ch))
            return 2
    # two genotyped datasets that share samples: CRISPRme resolves each dataset on its own
    # genotypes and its merge joins the rows of all of them, so a shared sample can be
    # listed by either dataset's haplotypes; this checker reads one haplotype set per
    # sample, so it would model such a run partially. It is refused, never half-checked
    for ch in sorted(chroms):
        typed = [v for v in vcfs_for.get(ch, ()) if v.samples]
        for i, v in enumerate(typed):
            for w in typed[i + 1:]:
                both = set(v.samples) & set(w.samples)
                if ctx.panel is not None:
                    both &= ctx.panel
                if not both or not (v.holds(ch) and w.holds(ch)):
                    continue
                sys.stderr.write(
                    "ERROR: %s and %s both genotype %d sample(s) of the panel on %s "
                    "(%s%s): two genotyped datasets that share samples are not modeled -- "
                    "CRISPRme resolves each dataset on its own genotypes and its merge "
                    "joins their rows, while this checker reads one haplotype set per "
                    "sample -- so the run is not checked. Run a search of each dataset "
                    "alone and check each (giving this run one of its VCFs is refused "
                    "too); a sites-only dataset beside a genotyped one, or datasets with "
                    "no sample in common, are checked.\n"
                    % (v.path, w.path, len(both), ch, ",".join(sorted(both)[:6]),
                       ", ..." if len(both) > 6 else ""))
                return 2
    # how the run's CRISPRme and CRISPRitz read its data, where their releases differ
    # (Semantics): read from its own records, and refused where they do not tell and the
    # data makes it matter
    genome = os.path.basename(os.path.normpath(run_params(args.results).get(
        "Genome_ref", run_params(args.results).get("Genome_selected", "")))) or "hg38"
    sem, why = run_semantics(args.results, vcfs_all, sorted(chroms), datasets, genome)
    if sem is None:
        sys.stderr.write(
            "ERROR: the run's log_verbose.txt says the SNP+indel co-occurrence was off "
            "('%s'): every SNP code of its fake-indel genomes read as the reference base, "
            "so its indel rows name no SNP and list every carrier of their indel, its SNP "
            "rows keep a carrier whose haplotype holds an indel in the target, and none of "
            "the per-carrier resolution this checker models applies (METHODS section 4); "
            "the run is not checked\n" % OPT_OUT)
        return 2
    ctx.sem = sem
    try:
        ctx.max_degeneracy = int(run_params(args.results)["Max_degeneracy"])
    except (KeyError, ValueError):
        m = None
        p = os.path.join(args.results, ".command_line.txt")
        if os.path.isfile(p):
            with open(p) as fh:
                m = re.search(r"--max-degeneracy\s+(\d+)", fh.read())
        ctx.max_degeneracy = int(m.group(1)) if m else None
    full_n = max(len(ctx.full(g)) for g in guides)   # the longest site of a guide
    unknown = []
    for ch in sorted(chroms):
        for v in vcfs_for.get(ch, ()):
            ds = sem.dataset_of.get(v.path)
            val = sem.values.get((ds, ch)) or sem.values.get((None, None)) or {}
            if v.samples and val.get("store_filter") is None and any(
                    r["filter"] not in PASS_FILTERS and any(
                        len(a) != len(r["ref"]) for a in r["alts"])
                    for r in v.records(ch, 0, 1 << 31)):
                unknown.append("which indels the phased indel genotype store of %s held "
                               "on %s (since 2.10.0 only those of records whose FILTER is "
                               "PASS or '.', which its log tells; the run has no "
                               "log_verbose.txt), and %s holds indel records of another "
                               "FILTER there" % (ds or v.path, ch, v.path))
            if val.get("all_alleles") is None and split_sites(v, ch, most=1, telling=False):
                unknown.append("whether the variant genome of %s on %s codes every allele "
                               "of a split SNP site (CRISPRitz 2.10.0) or the last "
                               "record's (2.9.1 and earlier) -- neither its variant index, "
                               "nor its enriched genome in the working directory, tells "
                               "-- and %s has such sites there"
                               % (ds or v.path, ch, v.path))
            fmt, pad = val.get("indels") or (None, None)
            # the nearest a site holding a base of its record lies to the window's first
            # base: its PAM ending on the record's first REF base (an MNV's first replaced
            # base), WINDOW_BEFORE bases in, its guide region from full_n - 1 bases before
            if fmt is None and v.samples is not None and \
                    WINDOW_BEFORE + 1 - full_n < (pad or max(args.bdna, args.brna) + 1):
                unknown.append("the index format of the INDELS index of %s on %s (its "
                               "%s names no index found, and its log no format-1 NOTE of "
                               "it), and a target of this guide, PAM and budget can lie at "
                               "a fake-indel window's edge, which an index of format 1 "
                               "drops and one of format 2 holds" % (ds or v.path, ch,
                                                                      INDEXES_SEARCHED))
    if unknown:
        sys.stderr.write("ERROR: the run's own records do not say %s: the run is not "
                         "checked\n" % "; nor ".join(dict.fromkeys(unknown)))
        return 2
    sys.stderr.write("check_indel_rows: read as a run of CRISPRme %s (%s)\n" % (
        "2.10.0 or later" if sem.release2100 else "2.8.1 or earlier",
        "; ".join("%s%s: %s, store %s, INDELS format %s" % (
            k[0] or "every dataset", " on %s" % k[1] if k[1] else "",
            {True: "every allele", False: "the last record of a split site"}.get(
                v.get("all_alleles"), "no split SNP site to tell"),
            "of PASS and '.'" if v.get("store_filter") else "of every FILTER",
            (v.get("indels") or (None,))[0]) for k, v in sorted(
                sem.values.items(), key=lambda kv: (str(kv[0][0]), str(kv[0][1]))))
        or "no chromosome to read it for"))
    if not ref.keep_all:
        ref.cache = {}
    reported = {}
    for ch in chroms:
        here = [r for r in indel_rows if r["chrom"] == ch]
        reported[ch] = Reported(here, ctx, loci={k: v for k, v in loci.items()
                                                  if k[1] == ch},
                                merge=merge, criterion=criterion, scoring=scoring)
    # the samples SNP rows list per site and variant list (they merge their carriers
    # there), for the SNP-row report's missing carriers
    snp_listed = defaultdict(set)
    if args.snp_rows == "report":
        for r in snp_rows:
            t = r["alt"].replace("-", "")
            snp_listed[(r["guide"], r["chrom"], cluster_of(
                r["start"], r["strand"], len(t), len(r["spacer_col"]), pam_5prime),
                frozenset(x for x in r["var"].split(",") if x and x != "NA"))].update(
                    samples_of(r))

    units = []
    tile = args.tile
    for ch in sorted(chroms):
        if ch not in lengths:
            continue
        if args.mode in ("soundness", "both"):
            tiles = defaultdict(list)
            for r in by_chrom.get(ch, ()):
                tiles[r["start"] // tile].append(r)
            for t in sorted(tiles):
                units.append(("sound", ch, t * tile, (t + 1) * tile, tiles[t]))
        if args.mode in ("completeness", "both"):
            clo = lo if lo is not None else 0
            chi = min(hi, lengths[ch]) if hi is not None else lengths[ch]
            t = clo
            while t < chi:
                units.append(("comp", ch, t, min(chi, t + tile), []))
                t += tile
    G.update(ctx=ctx, ref=ref, vcfs_for=vcfs_for, reported=reported, popmode=popmode,
             snp_listed=snp_listed)
    done = 0
    # each unit's findings are appended to <out>.partial as it ends, so a long run that is
    # stopped keeps what it found; the sorted TSV replaces it at the end
    partial = None
    if args.out:
        partial = open(args.out + ".partial", "w", newline="")
        pw = csv.DictWriter(partial, fieldnames=COLS, delimiter="\t", lineterminator="\n",
                            extrasaction="ignore")
        pw.writeheader()

    sites_only = defaultdict(int)   # chromosome -> records only sites-only VCFs hold
    filtered = defaultdict(int)     # chromosome -> records the search skips for FILTER
    long_ref = {}                   # unindexed VCF -> the longest REF its reads met
    mixed = set()                   # (chromosome, POS) of two lines of several ALTs

    def keep(result):
        found, stats = result
        rep.rows.extend(found)
        mixed.update((stats["chrom"], p) for p in stats.get("mixed", ()))
        sites_only[stats["chrom"]] += stats.get("sites_only", 0)
        filtered[stats["chrom"]] += stats.get("filtered", 0)
        for path, n in stats.get("long_ref", {}).items():
            long_ref[path] = max(n, long_ref.get(path, 0))
        if partial is not None:
            for r in found:
                pw.writerow(r)
            partial.flush()
    if args.jobs > 1 and len(units) > 1:
        import multiprocessing as mp
        ctx_mp = mp.get_context("fork")
        with ctx_mp.Pool(args.jobs) as pool:
            for found in pool.imap_unordered(run_unit, units, chunksize=1):
                keep(found)
                done += 1
                if args.progress and done % args.progress == 0:
                    sys.stderr.write("  %d/%d unit(s), %d finding(s)\n"
                                     % (done, len(units), len(rep.rows)))
                    sys.stderr.flush()
    else:
        for u in sorted(units, key=lambda x: (x[1], x[2])):
            keep(run_unit(u))
            done += 1
            if args.progress and done % args.progress == 0:
                sys.stderr.write("  %s %d/%d unit(s), %d finding(s)\n"
                                 % (u[1], done, len(units), len(rep.rows)))
                sys.stderr.flush()
    if mixed:
        if partial is not None:
            partial.close()
            os.remove(args.out + ".partial")
        sys.stderr.write(
            "ERROR: %s hold%s two or more VCF lines that each list several ALTs (a REF of "
            "one base, one ALT of one character at least, FILTER PASS or '.'), on a variant "
            "genome that codes every allele (CRISPRitz 2.10.0): its dictionary lists the "
            "one-base ALTs of every such line in the position's one value, each entry with "
            "its line's own allele numbers (1|2), which CRISPRme's readers cannot tell apart "
            "(METHODS section 4: an ALT can be tied to the other line's number and read on "
            "another haplotype slot), and this checker does not model what CRISPRme then "
            "lists: the run is not checked. Split such lines (bcftools norm -m -any, as the "
            "published panels are built) or join them into one, and search again.\n"
            % (", ".join("%s:%d" % m for m in sorted(mixed)[:5])
               + (" and %d more" % (len(mixed) - 5) if len(mixed) > 5 else ""),
               "s" if len(mixed) == 1 else ""))
        return 2
    for ch in sorted(c for c, n in sites_only.items() if n):
        rep.add(False, "limit-no-genotypes", {"table": "completeness", "chrom": ch},
                "%d indel/MNV record(s) of %s%s are held only by sites-only VCF(s) (no "
                "sample columns): no sample carries them, so there is no haplotype to "
                "enumerate and nothing to demand (a sites-only dataset's rows list no "
                "carrier; soundness checks those rows themselves)"
                % (sites_only[ch], ch, "" if lo is None else " in the region"))
    for ch in sorted(c for c, n in filtered.items() if n):
        # what a genotyped dataset's such indels are, by the phased indel store the run read
        rules = {sem.store_filter(v, ch) for v in vcfs_for.get(ch, ()) if v.samples}
        store = {
            frozenset([True]): " A genotyped dataset's indels among them are read as the "
                               "reference too: the phased indel store of CRISPRme 2.10.0 "
                               "skips them, as the enricher does.",
            frozenset([False]): " A genotyped dataset's indels among them stay on their "
                                "carriers' haplotypes, as the phased indel store the run read "
                                "(written before CRISPRme 2.10.0), which keeps every indel of "
                                "the VCF, reads them (a second indel-path record of the "
                                "targets they change).",
            frozenset([True, False]): " A genotyped dataset's indels among them are read as "
                                      "the reference where the run read a phased indel store "
                                      "of CRISPRme 2.10.0, and stay on their carriers' "
                                      "haplotypes where it read one written before.",
        }.get(frozenset(rules), "")
        rep.add(False, "limit-filtered-record", {"table": "completeness", "chrom": ch},
                "%d indel/MNV record(s) of %s%s have a FILTER other than PASS or '.', which "
                "CRISPRitz's enricher skips: no fake-indel window or indel log holds them, "
                "so no row of theirs is demanded.%s" % (
                    filtered[ch], ch, "" if lo is None else " in the region", store))
    for path in sorted(long_ref):
        rep.add(False, "coverage-unindexed-long-ref", {"table": ""},
                "this checker's own coverage gap, not a CRISPRme limit: %s has no .tbi "
                "and holds a record whose REF spans %d bp; a read without an index starts "
                "%d bp before its region, so a deletion longer than that which starts "
                "further before a region is not read there (its carriers' haplotypes there "
                "lack it). Index the VCF (tabix -p vcf, on a copy or a link in a folder of "
                "your own) and run again" % (path, long_ref[path], Vcf.LONGEST))
    if args.known:
        modes = {"soundness", "completeness"} if args.mode == "both" else {args.mode}
        apply_known(rep, load_known(args.known), chroms, modes)
    # every column in the key: the same input gives the same TSV, whatever order the work
    # units ended in
    rep.rows.sort(key=lambda r: (r["chrom"], str(r["start"]).zfill(12), r["class"],
                                 r["table"], str(r["line"]), r["block"], r["strand"],
                                 r["target"], r["variants"], r["samples"], r["detail"],
                                 r["gating"], r["spacer"]))
    if args.out:
        rep.write(args.out)
        partial.close()
        os.remove(args.out + ".partial")
    counts = rep.counts()
    print("check_indel_rows: %d finding(s), %d gating" % (len(rep.rows), rep.gating()))
    for cls in sorted(counts):
        g = sum(1 for r in rep.rows if r["class"] == cls and r["gating"] == "yes")
        print("  %-34s %7d%s" % (cls, counts[cls], "" if g else "  (non-gating)"))
    print("check_indel_rows: %s" % ("FAIL" if rep.gating() else "PASS"))
    return 1 if rep.gating() else 0


if __name__ == "__main__":
    sys.exit(main())
