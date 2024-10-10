#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --account=che210028p
#SBATCH -N 1
#SBATCH -p GPU-shared
#SBATCH --ntasks-per-node=12
#SBATCH --mem=160G
#SBATCH --job-name=topoformer
#SBATCH --output=slurm_outputs/topoformer_%j.out
#SBATCH --error=slurm_outputs/topoformer_%j.err
#SBATCH --gres=gpu:v100-32:3
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=sabarik@colostate.edu

cd /ocean/projects/che210028p/skumar7/ProteinSol/topoformer 
APPTAINER_TMPDIR=/ocean/projects/che210028p/skumar7/singularity_temp APTAINERENV_WANDB_API_KEY=046f107c2d4757c1423c36e02a3b22da8ed18c71 APPTAINERENV_WANDB_PROJECT=topoformer apptainer exec --nv -B /ocean/projects/che210028p/skumar7/ProteinSol/ 20240927_topoformer_fixedgrads.sif /usr/bin/python /ocean/projects/che210028p/skumar7/ProteinSol/topoformer/scripts/sequential_hyperparams_bridges.py -6 -2 10 2 10 2 7 20 1
