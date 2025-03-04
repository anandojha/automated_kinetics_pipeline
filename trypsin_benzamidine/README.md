## Automated simulation pipeline for kinetic and thermodynamic predictions

The following series of scripts automates the process of running milestoning simulations to compute kinetic and thermodynamic properties for receptor-ligand interactions. Following this workflow, users can efficiently obtain binding and unbinding rates (kon, koff) for a given protein-ligand complex. To execute the pipeline, users must first install the necessary dependencies (refer to the README.md file for setup instructions) and ensure they work in the one_step_kinetics conda environment. Before running the scripts, the user must specify the correct paths to SEEKR2 and SEEKRTools in the relevant scripts. These paths need to be modified in:

 - run_initial_simulation.py (for SEEKRTools path)
 - run_milestoning.py (for SEEKR2 path)
 - run_analysis.py (for SEEKR2 path)


#### 1. get_protein_ligand.py
This script extracts the protein and ligand from a given PDB file containing a receptor-ligand complex. It separates the ligand based on a specified residue name and outputs two new PDB files: protein.pdb for the protein and ligand.pdb for the ligand. It also converts the ligand structure into SDF format (ligand.sdf) to prepare it for force field parameterization.

```sh
python get_protein_ligand.py 
```

#### 2. parameterize.py
This script applies force field parameters to the protein and ligand employing ESPALOMA (for the solute) and AMBER force fields (for the solvent). It creates a solvated and minimized structure that is ready for simulations. The system is prepared with the correct force field assignments, and necessary OpenMM-compatible system files are generated. These files are required for running molecular dynamics simulations. 

```sh
python parameterize.py 
```

#### 3. create_milestoning_input.py
This script determines the milestones (boundaries) for the SEEKR2 simulation by analyzing the center-of-mass (COM) distances between the protein and ligand. It processes the minimized PDB file, identifies key atoms, and generates the required model.xml file, which contains information about the reaction coordinates and milestones for the milestoning simulations.

```sh
python create_milestoning_input.py 
```

#### 4. prepare_milestoning.py
This script sets up the directory structure and necessary files for running SEEKR2 milestoning simulations. It ensures that all required files, including force field and system parameter files, are correctly formatted and placed in the appropriate locations. It processes the model.xml file and prepares SEEKR2-compatible input directories.

```sh
python prepare_milestoning.py --seekr2_path "/home/USERNAME/seekr2"
```

#### 5. run_initial_simulation.py
This script performs initial enhanced sampling using HIDR from SEEKRTools. It generates a set of starting structures at different milestones to initiate simulations. The initial sampling (Steered MD, RAMD, or metadynamics) helps populate all accessible regions of the system so that milestoning simulations can begin with well-distributed conformations.

```sh
python run_initial_simulation.py --seekrtools_path "/home/USERNAME/seekrtools"
```

#### 6. run_milestoning.py
This script runs the milestoning simulations. It executes multiple MD simulations for different milestones defined in model.xml. These simulations track the transitions between milestones, which are later used for computing kinetics and free energy differences.

```sh
python run_milestoning.py --seekr2_path "/home/USERNAME/seekr2"
```

#### 7. run_analysis.py
This script performs the final analysis of the milestoning simulation results. It calculates the transition rates between milestones, determines kinetic and thermodynamic properties, and generates the final results in analyze.out. The computed values help in understanding ligand-binding kinetics and other biophysical properties.

```sh
python run_analysis.py --seekr2_path "/home/USERNAME/seekr2"
```
After running the run_analysis.py script, an images_and_plots directory is created inside the SEEKR_SIMULATION folder. This directory contains visualization and analysis files that summarize the milestoning simulation results. 


Here is the summary of all the scripts to run in sequence:

```sh
python get_protein_ligand.py 
python parameterize.py 
python create_milestoning_input.py 
python prepare_milestoning.py --seekr2_path "/home/USERNAME/seekr2"
python run_initial_simulation.py --seekrtools_path "/home/USERNAME/seekrtools"
python run_milestoning.py --seekr2_path "/home/USERNAME/seekr2"
python run_analysis.py --seekr2_path "/home/USERNAME/seekr2"
```

The run_all.sh script automates the execution of the milestoning pipeline by defining two key variables:
 
```sh
SEEKR2_PATH="/home/USERNAME/seekr2"
SEEKRTOOLS_PATH="/home/USERNAME/seekrtools"
```
Replace "USERNAME" with the home directory name. The script then runs all required Python commands, passing these variables to the relevant scripts. To execute it, first make it executable:

```sh
chmod +x run_all.sh
```

```sh
./run_all.sh
```
