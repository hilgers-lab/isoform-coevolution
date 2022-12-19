import os
import glob
from os.path import join as pjoin
import pandas as pd

MAFS = glob.glob('input/mafs/*.maf')
LABELS, = glob_wildcards('input/gff/{label}GenesRegions.gff')

BOX_CSV = config['gene_regions_csv']
box_df = pd.read_csv(BOX_CSV, index_col=False, 
                        header=0)

box_genes = list(box_df['gene_id'].unique())

box_plot_dirs = expand(
    "output/coevo-box/label-{label}_gene-{gene}/",
    label=LABELS,
    gene=box_genes,
)

rule flags: 
    input: 
        expand('flags/{label}.done', label=LABELS)
    output: 
        'flags/all.done'

rule all: 
    input: 
        box_plot_dirs


# output a MAF file for each gene
checkpoint gff_to_mafs:
    input: 
        gff='input/gff/{label}GenesRegions.gff', 
        mafs=MAFS
    output:
        tmp_bed=temp('mafs/{label}_tmp.bed'),
        maf_dir=directory('mafs/{label}')
    conda: 
        "envs/bedops_and_ucsc_mafs.yaml"
    shell:
        "gff2bed < {input.gff} | " # turn the gff into a NEWT
        "sed 's/^/chr/' | " # add chr to each
        """awk -F '\\t' 'OFS="\\t" {{print $1, $2, $3, $10}}' | """ # use gene_id as 4th column
        "sed 's/gene_id //g' | " # keep only the actual gene_id
        "sed 's/\\\"//g' > {output.tmp_bed}; " # remove pesky quotes
        "mafsInRegion {output.tmp_bed} -outDir {output.maf_dir} {input.mafs}"

   
def get_file_names(wc):
    "Here be magic"
    ck_output = checkpoints.gff_to_mafs.get(**wc).output['maf_dir']
    GENES, = glob_wildcards(pjoin(ck_output, '{gene}.maf'))
    return expand('output/coevo-box/label-{label}_gene-{gene}', gene=GENES, label=wc.label)

rule almost_done: 
    input: get_file_names
    output: "flags/{label}.done"
    shell: "touch {output}"


rule maf_to_fa: 
    input: 
        'mafs/{label}/{gene}.maf'
    output: 
        'fasta/{label}/{gene}.fasta'
    conda: 
        'envs/bx-python.yaml'
    shell: 
        "export SPECIES=$(maf_species_in_all_files.py {input}); "
        "maf_to_concat_fasta.py $SPECIES < {input} > {output}"



# temporarily store the matrix since we read it multiple times
rule dna_coevo_matrix: 
    input: 
        'fasta/{label}/{gene}.fasta'
    output: 
        raw_mat='matrix/{label}/{gene}_raw.mat',
        refined_mat='matrix/{label}/{gene}_refined.mat'
    conda: 
        'envs/coevo.yaml'
    benchmark: 
        'benchmarks/dna_coevo_matrix/label-{label}_gene-{gene}.txt'
    resources: 
        mem_mb=40000
    shell: 
        "python scripts/DNAcoevolution.py "
        "-i {input} "
        "-s {output.raw_mat} "
        "-v {output.refined_mat} "

rule dna_coevo_plain: 
    input: 
        raw_mat='matrix/{label}/{gene}_raw.mat',
        refined_mat='matrix/{label}/{gene}_refined.mat'
    output: 
        directory('output/coevo-plain/label-{label}_gene-{gene}')
    conda: 
        'envs/coevo.yaml'
    resources: 
        mem_mb=40000
    benchmark: 
        'benchmarks/dna_coevo_plain/label-{label}_gene-{gene}.txt'
    shell: 
        "python scripts/DNAcoevolution.py "
        "-a {input.raw_mat} "
        "-r {input.refined_mat} "
        "-o {output} "


rule dna_coevo_box: 
    input: 
        fasta='fasta/{label}/{gene}.fasta',
        csv='config/all_box_genes.csv',
    output: 
        directory('output/coevo-box/label-{label}_gene-{gene}')
    conda: 
        'envs/coevo.yaml'
    resources: 
        mem_mb=80000,
        disk_mb=40000
    benchmark: 
        'benchmarks/dna_coevo_box/label-{label}_gene-{gene}.txt'
    shell: 
        "python scripts/DNAcoevolution2.py "
        "-i {input.fasta} "
        "-c {input.csv} "
        "-o {output} "
        "-g {wildcards.gene} "

# rule covariance_gene_level:
#     input: 
#     shell: 
