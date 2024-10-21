#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --account=che210028p
#SBATCH -N 1
#SBATCH -p GPU-shared
#SBATCH --ntasks-per-node=12
#SBATCH --mem=160G
#SBATCH --job-name=topoformer
#SBATCH --output=$SLURM_OUTPUT_DIR_%j.out
#SBATCH --error=$SLURM_OUTPUT_DIR_%j.err
#SBATCH --gres=gpu:v100-32:3
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=sabarik@colostate.edu

cd /ocean/projects/che210028p/skumar7/ProteinSol/topoformer 
APPTAINER_TMPDIR=/ocean/projects/che210028p/skumar7/singularity_temp apptainer exec --nv -B /ocean/projects/che210028p/skumar7/ProteinSol/ /ocean/projects/che210028p/skumar7/ProteinSol/DuMD/20241015_dumd.sif /usr/bin/python /ocean/projects/che210028p/skumar7/ProteinSol/DuMD/run_simulation.py $PDB_FILE /ocean/projects/che210028p/skumar7/ProteinSol/DuMD/runtime/config/config.yaml
