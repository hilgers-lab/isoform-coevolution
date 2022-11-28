import os
import glob
from os.path import join as pjoin

MAFS = glob.glob('input/mafs/*.maf')
LABELS, = glob_wildcards('input/gff/{label}GenesRegions.gff')

rule all: 
    input: expand('flags/{label}.done', label=LABELS)


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
    return expand('figures/{label}/{gene}_coevo.png', gene=GENES, label=wc.label)


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

# if there is a TSV with the following columns: 
# chr, start, end, feature at the given path tsv/label/gene.tsv
# then color in the respective rectangles in the 2D array and output the 
# average for each rectangle
def _covariance_gene_input(wc):
    tsv_path = 'tsv/{}/{}.tsv'.format(wc.label, wc.gene):
    return '-t {}'.format(tsv_path) if os.path.exists(tsv_path) else ''

rule dna_coevo: 
    input: 
        'fasta/{label}/{gene}.fasta'
    output: 
        directory('output/{label}/{gene})'
    params: 
        tsv=_covariance_gene_tsv
    conda: 
        'envs/coevo.yaml'
    shell: 
        "python scripts/DNAcoevolution.py "
        "-i {input.fasta} "
        "{params.tsv} "
        "-o {output} "

checkpoint tsv_genes:
    input: 
        config['input/gene_regions.tsv']
    output: 
        touch('flags/tsv_genes.checkpoint')


def get_tsv_gene_names(wc):
    "Here be magic"
    ck_output = checkpoints.tsv_genes.get(**wc).output[0]
    df = pd.read_table(ck_output).set_index('gene_id')
    genes = list(df.index)
    return expand('figures/{label}/{gene}_coevo_with_boxes.png', 
                  gene=GENES, label=wc.label)


# rule covariance_gene_level:
#     input: 
#     shell: 
