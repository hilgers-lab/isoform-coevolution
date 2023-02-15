eval "$(conda shell.bash hook)"
conda activate snakemake
#snakemake --profile ~/.config/mpi-ie-slurm -j32
#snakemake --configfile config/config.yaml --profile ~/.config/mpi-ie-slurm -j32
snakemake --configfile config/config.yaml --profile ~/.config/mpi-ie-slurm -j8 flags/all.done
