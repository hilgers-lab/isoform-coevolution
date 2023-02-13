import pandas as pd

dominant_tsv = "top50/top50.coEvoGenes.fullList.tsv"
nondominant_tsv = "top50/top50.control.coEvoGenes.tsv"
outfile = "config/all_box_genes.csv"

dominant_df = pd.read_csv(dominant_tsv, header=0)
nondominant_df = pd.read_csv(nondominant_tsv, header=0)

nondominant_df['tissue'] = 'controlGenes'

assert len(dominant_df.columns) == len(nondominant_df.columns)

nondominant_df.columns = list(dominant_df.columns)

geneType_rename = {
    'domPromoterGene': 'dominant',
    'controlGenes': 'non.dominant'
}

final_df = pd.concat([dominant_df, nondominant_df], axis=0)
final_df["geneType"] = final_df["geneType"].map(lambda x: geneType_rename[x])
final_df.to_csv(outfile, header=True, index=False)

