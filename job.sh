#!/bin/bash
#SBATCH --cpus-per-task=64
#SBATCH --job-name=split-chains
#SBATCH --partition cpuq


python get_pdb_data.py --start_date now-1w --stop_date now --search_type polymer_type --search_term RNA &&
python split_chains.py &&
python get_xyz_data.py &&
python filter_by_interC1.py &&
python get_pub_dates_and_dedupe.py &&
python align_sequence.py &&
python generate_train.py
