# Automated Simulation Pipeline for Kinetic and Thermodynamic Predictions

This repository provides an **end-to-end automated pipeline** for running **milestoning simulations** with machine-learned force fields to accelerate **drug-target kinetic and thermodynamic predictions**. The pipeline utilizes **SEEKR2**, **SEEKRTools**, and **espaloma=0.3.2** to automate simulations.

## **Getting Started**

To install and set up the necessary environment, run the following commands in sequence.

```sh
# Create and activate Conda environment
conda create -n one_step_kinetics python=3.10 --yes
conda activate one_step_kinetics

# Install Mamba for faster dependency resolution
conda install conda-forge::mamba --yes

# Install OpenMM and SEEKR2 Plugin
mamba install seekr2_openmm_plugin openmm=8.1 --yes

# Verify SEEKR2 OpenMM Plugin Installation
python -c "import seekr2plugin"

# Install Git (if not installed)
conda install conda-forge::git --yes

# Clone and install SEEKR2
git clone https://github.com/seekrcentral/seekr2.git
cd seekr2
python -m pip install .
pytest   # Run tests to verify installation
cd ..

# Clone and install SEEKRTools
git clone https://github.com/seekrcentral/seekrtools.git
cd seekrtools
python -m pip install .
pytest  # Run tests to verify installation
cd ..
