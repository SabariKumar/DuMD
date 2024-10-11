### DuMD: DumbMD
![dumd](./src/DuMD.jpeg)


Submits SLURM jobs to mindlessly run OpenMM simulations for pdb files in a folder.
Nonstandard amino acids, heteroatoms, and missing atoms are treated per the pdbfixer defaults.
See `runtime/config/default_nvt_config.yaml` for other default run parameters. 


Developed for usage on bridges2 - modify the run_slurm_bridges.sh to your needs.


Usage:
1. scp over our OpenMM apptainer container sif file to your repo directory
2. Run the following bash command - replace the pdb_search_dir with the pdb directory and slurm_output_dir
   with the location for slurm output and error logs:


`find pdb_search_dir -name "*.pdb" -exec sbatch --export=PDB_FILE={} --export=SLURM_OUTPUT_DIR=/ocean/projects/che210028p/skumar7/slurm_dumps run_slurm_bridges.sh`
