eval "$(conda shell.bash hook)"
conda activate snakemake
snakemake --profile ~/.config/mpi-ie-slurm -j16
