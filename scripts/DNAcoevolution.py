import os
from prody import MSAFile, parseMSA, refineMSA, showMutinfoMatrix
from prody import buildMutinfoMatrix, applyMutinfoCorr, writeArray
import click
import pandas as pd
from pathlib import Path
import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
import seaborn as sns

REGIONS_OF_INTEREST = [
    'regularPromoterCoord',
    'regularutr5Region',
    'dominantPromoterCoord',
    'utr5Region',
    'coding_exonRegion',
    'utr3Region',
    'distalutr3Region',
]

COLORS = sns.color_palette('colorblind')[:len(REGIONS_OF_INTEREST)]

# Data
#file_path = "/data/hilgers/group/alfonso/projects/2021_LRS_paper/data/functional.analysis/fasta"
#staiFa = "/StaiGeneMSA.fasta"
# load
#msaFile = file_path+staiFa
#msafobj = MSAFile(msaFile)
#msa = parseMSA(msaFile)

# processing
# test params multiple resolutions

#msa_refine04 = refineMSA(msa, label='dm6', rowocc=0.4, seqid=0.98)
#writeMSA('staiTrimmed.fasta', msa_refine)
# calculate mutual information
#mutinfo04 = buildMutinfoMatrix(msa_refine04,ambiguity=False)
# mut cor
#mutinfo_corr04 = applyMutinfoCorr(mutinfo04, corr='apc')
#writeArray('stai.mutinfo_corr04.txt', mutinfo_corr04)
# plotting mut info corr
# plotting mut info corr
# plt.clf()
# showMutinfoMatrix(mutinfo_corr04, clim=[0, mutinfo_corr04.max()],vmin=0,vmax=0.5)
# plt.savefig('stai.mutinfo_corr04.other.vmax05.png', dpi=300, bbox_inches='tight')
# plt.show()
# plt.clf()
# showMutinfoMatrix(mutinfo_corr04, clim=[0, mutinfo_corr04.max()],vmin=0,vmax=0.4)
# plt.savefig('stai.mutinfo_corr04.othervmax04.png', dpi=300, bbox_inches='tight')

@click.option('-i', '--msa-file', type=click.Path(exists=True))
def build_mutinfo_from_msa(msa_file, out_array):
    msa = parseMSA(msa_file)
    msa_refine04 = refineMSA(msa, label='dm6', rowocc=0.4, seqid=0.98)
    mutinfo04 = buildMutinfoMatrix(msa_refine04, ambiguity=False)
    mutinfo_corr04 = applyMutinfoCorr(mutinfo04, corr='apc')
    writeArray(out_array, mutinfo_corr04)


@click.command()
@click.option('-i', '--msa-file', type=click.Path(exists=True))
@click.option('-o', '--outdir', type=click.Path(exists=False, path_type=Path))
@click.option('-g', '--gene-id', type=str)
@click.option('-t', '--tsv', type=click.File('r'), required=False)
@click.option('-a', '--from-array', type=click.Path(exists=True), 
              required=False)
def main(msa_file, outdir, tsv, gene_id, from_array):
    # initiate directories
    os.makedirs(outdir, exist_ok=True)

    # if we have already saved the array, just read it in
    if from_array is not None: 
        print('Reading in file from {}...'.format(from_array), file=sys.stderr)
        mutinfo_corr04 = np.loadtxt(from_array)
    else:
        # MSAFile(msa_file)
        msa = parseMSA(msa_file)
        msa_refine04 = refineMSA(msa, label='dm6', rowocc=0.4, seqid=0.98)
        mutinfo04 = buildMutinfoMatrix(msa_refine04, ambiguity=False)
        mutinfo_corr04 = applyMutinfoCorr(mutinfo04, corr='apc')

    array_path = outdir / 'array.txt'
    #writeArray(array_path, mutinfo_corr04)
    ga = sns.clustermap(mutinfo_corr04, row_cluster=False, 
                            col_cluster=False, vmin=0,vmax=0.5)
    ax = ga.ax_heatmap

    # dict of color, label for legend
    legend_handles = []

    if tsv is not None:
        df = pd.read_table(tsv).set_index('gene_id')
        if gene_id in df.index: 
            gene_row = df.loc[gene_id]
            for color, region in zip(COLORS, REGIONS_OF_INTEREST):
                if region not in gene_row:
                    continue
                start, end = (int(x) for x in gene_row[region].split('-'))
                gene_start = int(gene_row['geneStart'])
                start, end = (coord - gene_start for coord in (start, end))
                box_length = end - start
                ax.add_patch(Rectangle((start, start), box_length, box_length, 
                             fill=False, edgecolor=color, lw=3))
                legend_handles.append({'color': color, 'label': region})


    ax.legend(handles=[Patch(**kwargs) for kwargs in legend_handles], 
              loc='upper left', bbox_to_anchor=(0.7, 1.04))

    # save clustermap to file
    clustermap_path = outdir / 'clustermap.png'
    plt.savefig(clustermap_path, dpi=300, bbox_inches='tight')
    #import pdb; pdb.set_trace()
    plt.show()
    # plotting mut info corr
    # plt.clf()
    # showMutinfoMatrix(mutinfo_corr04, clim=[0, mutinfo_corr04.max()],vmin=0,vmax=0.5)
    # corr05_path = outdir / 'mutinfo_corr04.othervmax05.png'
    # plt.savefig(corr05_path, dpi=300, bbox_inches='tight')
    # plt.show()
    # plt.clf()
    # corr04_path = outdir / 'mutinfo_corr04.othervmax04.png'
    # showMutinfoMatrix(mutinfo_corr04, clim=[0, mutinfo_corr04.max()],vmin=0,vmax=0.4)
    # plt.savefig(corr04_path, dpi=300, bbox_inches='tight')


if __name__ == '__main__':
    main()