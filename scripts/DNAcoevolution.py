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
from collections import defaultdict
from itertools import combinations_with_replacement
from scipy.stats import ttest_1samp


REGIONS_OF_INTEREST = [
    'controlTSS',
    'control5utr',
    'targetsTSS',
    'target5utr',
    'cdsexons',
    'utr3Segment',
    'utr3segs.distalSegment',
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

def get_region(row, region, gene_start_col="geneStart"):
    "Get the region normalized to the geneStart"
    start, end = (int(x) for x in row[region].split('-'))
    gene_start = int(row[gene_start_col])
    start, end = (coord - gene_start for coord in (start, end))
    return (start, end)

def build_mutinfo_from_msa(msa_file):
    "Return raw and corrected mutinfo matrix"
    msa = parseMSA(msa_file)
    msa_refine04 = refineMSA(msa, label='dm6', rowocc=0.4, seqid=0.98)
    mutinfo04 = buildMutinfoMatrix(msa_refine04, ambiguity=False)
    return mutinfo04

@click.command()
@click.option('-i', '--msa-file', type=click.Path(exists=True), required=False)
@click.option('-o', '--outdir', type=click.Path(exists=False, path_type=Path), required=False)
@click.option('-g', '--gene-id', type=str, required=False)
@click.option('-c', '--csv', type=click.Path(exists=True), required=False)
@click.option('-a', '--from-array', type=click.Path(exists=True), 
              required=False)
@click.option('-s', '--save-array',  required=False,
              type=click.Path(exists=False, path_type=Path))
@click.option('-r', '--from-refined-array', required=False,
              type=click.Path(exists=True))
@click.option('-v', '--save-refined-array', required=False, 
              type=click.Path(exists=False, path_type=Path))
def main(msa_file, outdir, gene_id, csv, from_array, save_array, 
         from_refined_array, save_refined_array):
    # if we have already saved the array, just read it in
    if from_array is not None: 
        print('Reading in file from {}...'.format(from_array), file=sys.stderr)
        mutinfo04 = np.loadtxt(from_array)
    else:
        # MSAFile(msa_file)
        mutinfo04 = build_mutinfo_from_msa(msa_file)

    if save_array is not None:
        writeArray(str(save_array.resolve()), mutinfo04)

    if from_refined_array is not None: 
        mutinfo_corr04 = np.loadtxt(from_refined_array)

    else:
        # apply correction
        mutinfo_corr04 = applyMutinfoCorr(mutinfo04, corr='apc')
    
    if save_refined_array is not None:
        writeArray(str(save_refined_array.resolve()), mutinfo_corr04)

    # don't write to array for now as per discussion with Carlos
    #array_path = outdir / 'array.txt'
    #writeArray(array_path, mutinfo_corr04)
    if outdir is not None: 
        # initiate directories
        os.makedirs(outdir, exist_ok=True)


    del mutinfo04
    mat_label, mat = ('corrected_mutinfo', mutinfo_corr04)

    ga = sns.clustermap(mat, row_cluster=False, 
                        col_cluster=False, vmin=0,vmax=0.3)
    ax = ga.ax_heatmap
    
    # save clustermap to file
    clustermap_path = outdir / '{}_clustermap.png'.format(mat_label)
    plt.show()
    plt.savefig(clustermap_path, dpi=300, bbox_inches='tight')
        #import pdb; pdb.set_trace()

    showMutinfoMatrix(mat, clim=[0, mat.max()], 
                      vmin=0, vmax=0.3)
    corr03_path = outdir / '{}.othervmax03.png'.format(mat_label)
    plt.savefig(corr03_path, dpi=300, bbox_inches='tight')
    plt.clf()

    if csv is not None:
        df = pd.read_csv(csv, index_col=None, header=0)

        # subset on gene_id and tissue
        gene_rows = df[df['gene_id'] == gene_id]

        for __, gene_row in gene_rows.iterrows():
            tissue = gene_row['tissue']

            # dict of color, label for legend
            legend_handles = []
            avgs = []

            plt.clf()

            ga = sns.clustermap(mat, row_cluster=False, 
                                col_cluster=False, vmin=0,vmax=0.3)
            ax = ga.ax_heatmap
    
            # save clustermap to file
            clustermap_path = outdir / '{}_tissue-{}_clustermap.png'.format(mat_label, tissue)

            for color, region in zip(COLORS, REGIONS_OF_INTEREST):
                if region not in gene_row or pd.isna(gene_row[region]):
                    continue
                    
                # already normalized to gene start
                start, end = get_region(gene_row, region)
                box_length = end - start
                ax.add_patch(Rectangle((start, start), box_length, box_length, 
                            fill=False, edgecolor=color, lw=3))
                legend_handles.append({'color': color, 'label': region})

            # compare all regions to all regions
            for (region_x, region_y) in combinations_with_replacement(
                REGIONS_OF_INTEREST, 2
            ):

                # set start and end for both regions
                ((start_x, end_x), (start_y, end_y)) = (
                    get_region(gene_row, region) for region in (region_x, region_y)
                )

                region_of_interest = mat[start_x:end_x, start_y:end_y]
                row_of_interest = mat[:, start_y:end_y]
                mean_region = region_of_interest.mean()
                mean_row = row_of_interest.mean()
                pval = ttest_1samp(region_of_interest.flatten(), 
                                row_of_interest.mean(), 
                                alternative='greater').pvalue
                avgs.append([region_x, region_y, mean_region, mean_row, pval])

            if avgs: 
                avg_df = pd.DataFrame(avgs)
                avg_df.columns = ['REGION_X', 'REGION_Y', 'MEAN_X_VS_Y', 'MEAN_Y', 'P_VAL']
                avg_path = outdir / '{}_tissue-{}_regions.tsv'.format(mat_label, tissue)
                avg_df.to_csv(avg_path, sep='\t', header=True, index=False)

            ax.legend(handles=[Patch(**kwargs) for kwargs in legend_handles], 
                      loc='lower left', bbox_to_anchor=(0.7, 1.04))

            plt.savefig(clustermap_path, dpi=300, bbox_inches='tight')
            #import pdb; pdb.set_trace()
            plt.show()
            # plotting mut info corr

            plt.clf()

    # plt.show()
    # plt.clf()
    # corr04_path = outdir / 'mutinfo_corr04.othervmax04.png'
    # showMutinfoMatrix(mutinfo_corr04, clim=[0, mutinfo_corr04.max()],vmin=0,vmax=0.4)
    # plt.savefig(corr04_path, dpi=300, bbox_inches='tight')

    # TODO intersections between every feature compared toutr3Region and distalutr3Region
    # mean of backgrounds -> select 100 at random (or so) of same length (w/o replace)
    # compare to the mean of the column i.e. mutinfo_corr04


if __name__ == '__main__':
    main()