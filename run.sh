eval "$(conda shell.bash hook)"
conda activate snakemake
#snakemake --profile ~/.config/mpi-ie-slurm -j32
#snakemake --configfile config/config.yaml --profile ~/.config/mpi-ie-slurm -j32
snakemake --configfile config/config.yaml --profile ~/.config/mpi-ie-slurm-deep17 -j2 flags/all.done
#snakemake --configfile config/config.yaml --profile ~/.config/mpi-ie-slurm-deep17 -j2 output/coevo-box/label-dominant_gene-FBgn0014010/corrected_mutinfo_clustermap.png
