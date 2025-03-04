#!/bin/bash

# Define SEEKR2 and SEEKRTools paths
SEEKR2_PATH="/mnt/home/aojha/seekr2"
SEEKRTOOLS_PATH="/mnt/home/aojha/seekrtools"

# Stop execution if any command fails
set -e

# Run the pipeline
python get_protein_ligand.py 
python parameterize.py 
python create_milestoning_input.py 
python prepare_milestoning.py --seekr2_path "$SEEKR2_PATH"
python run_initial_simulation.py --seekrtools_path "$SEEKRTOOLS_PATH"
python run_milestoning.py --seekr2_path "$SEEKR2_PATH"
python run_analysis.py --seekr2_path "$SEEKR2_PATH"

