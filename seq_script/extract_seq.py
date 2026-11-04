#Use bedtool to extract sequence
import subprocess
import tempfile
import os
from os.path import isfile, isdir,join      #for getting lst of chr to know file extension and if enriched
from os import listdir

#Input chr1:11,130,540-11,130,751
def extractSequence(name, input_range, genome_selected):
    """The sequence at input_range (chrom:start-end) of genome genome_selected, read by
    bedtools getfasta from the FASTA of chrom in Genomes/<genome_selected>/ of the data
    folder (the working directory).

    The search form passes on what was typed, and a crafted request any text: name (the
    sequence's '>' header) is not used, genome_selected must be a genome folder of
    Genomes/ and chrom one of its FASTAs (<chrom>.fa or <chrom>.fasta), else
    FileNotFoundError. bedtools runs with no shell, on a BED and a link to the FASTA in a
    private temporary folder, where it writes the FASTA's index too: nothing is written
    in the data folder, and no index there is removed."""
    chrom = input_range.split(':')[0]
    start_position = input_range.split(':')[1].split('-')[0].replace(',','').replace('.','').replace(' ','')
    end_position = input_range.split(':')[1].split('-')[1].replace(',','').replace('.','').replace(' ','')

    genomes = join(os.getcwd(), 'Genomes')
    genome_dir = join(genomes, genome_selected)
    if genome_selected.startswith('.') or genome_selected not in listdir(genomes) or not isdir(genome_dir):
        raise FileNotFoundError(f'No genome {genome_selected} in {genomes}')
    fasta = next((f for f in (chrom + '.fa', chrom + '.fasta')
                  if f in listdir(genome_dir) and isfile(join(genome_dir, f))), None)
    if fasta is None:
        raise FileNotFoundError(f'No FASTA of {chrom} in {genome_dir}')
    with tempfile.TemporaryDirectory() as tmp:
        bed = join(tmp, 'region.bed')
        with open(bed, 'w') as b:
            b.write(chrom + '\t' + start_position + '\t' + end_position)
        os.symlink(join(genome_dir, fasta), join(tmp, fasta))
        output_extract = subprocess.check_output(
            ['bedtools', 'getfasta', '-fi', join(tmp, fasta), '-bed', bed]).decode("utf-8")
    ret_string = output_extract.split('\n')[1].strip()
    return ret_string
