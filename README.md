## Automated Simulation Pipeline for Kinetic and Thermodynamic Predictions

This repository provides an **end-to-end automated pipeline** for running **milestoning simulations** with machine-learned force fields to accelerate **drug-target kinetic and thermodynamic predictions**. The pipeline utilizes **SEEKR2**, **SEEKRTools**, and **espaloma=0.3.2** to automate simulations.

### **Getting Started**

To install and set up the necessary environment, run the following commands in sequence.

#### 1. Create a new conda environment named 'one_step_kinetics' with Python 3.10
The --yes flag ensures it installs without prompting for confirmation.
```sh
conda create -n one_step_kinetics python=3.10 --yes
```

#### 2. Activate the conda environment
Once the environment is created, activate the conda environment, switching the current shell session to use the new environment.
```sh
conda activate one_step_kinetics
```

#### 3. Install mamba for faster dependency resolution
Mamba is a faster drop-in replacement for conda that speeds up package installations.
```sh
conda install conda-forge::mamba --yes
```

#### 4. Install SEEKR2 plugin
This step installs OpenMM plugin for SEEKR2 package.
```sh
mamba install seekr2_openmm_plugin openmm=8.1 --yes
```
#### 5. Verify SEEKR2 OpenMM Plugin Installation
Run the following command to check if the SEEKR2 OpenMM plugin is correctly installed. If no error message appears, the installation was successful.
```sh
python -c "import seekr2plugin"
```

#### 6. Install Git
Git is required to clone repositories. Ensure it is installed using:

```sh
conda install conda-forge::git --yes
```

#### 7. Install the SEEKR2 package.
SEEKR2 is the core package required for performing milestoning simulations. We will download it from GitHub, install it, and verify that it works correctly.

##### a. Before proceeding, go to the home directory (recommended).
```sh
cd ~
```
##### b. Clone the SEEKR2 repository from GitHub into the current directory.
```sh
git clone https://github.com/seekrcentral/seekr2.git
```
##### c. Navigate into the seekr2 directory where the cloned repository is located.
```sh
cd seekr2
```
##### d. Install SEEKR2 using pip, making it accessible as a Python package in the conda environment.
```sh
python -m pip install .
```
##### e. Run tests to verify that SEEKR2 has been installed correctly and is functioning as expected.
```sh
pytest
```
##### f. Return to the home directory (~), ensuring a clean workspace before proceeding to the next steps.
```sh
cd ~
```

#### 8. Install the SEEKRTools package
SEEKRTools is a companion package to SEEKR2 that provides utilities for preparing and facilitating multiscale milestoning simulations.

##### a. Before proceeding, go to the home directory (recommended).
```sh
cd ~
```
##### b. Clone the SEEKRTools repository from GitHub into the current directory.
```sh
git clone https://github.com/seekrcentral/seekrtools.git
```
##### c. Navigate into the seekrtools directory where the cloned repository is located.
```sh
cd seekrtools
```
##### d. Install SEEKRTools using pip, making it accessible as a Python package in the conda environment.
```sh
python -m pip install .
```
##### e. Run tests to verify that SEEKRTools has been installed correctly and is functioning as expected.
```sh
pytest
```
##### f. Return to the home directory (~), ensuring a clean workspace before proceeding to the next steps.
```sh
cd ~
```

#### 9. Install the Espaloma Force Field
Espaloma is a machine-learned molecular mechanics force field designed to improve the accuracy and flexibility of molecular simulations. To install Espaloma, we need to install its dependencies, including OpenFF Toolkit, OpenMM Force Fields, and OpenEye Toolkits.

##### a. Install the OpenFF Toolkit
The Open Force Field (OpenFF) Toolkit is required to assign and manipulate molecular mechanics parameters. 
```sh
conda install conda-forge::openff-toolkit --yes
```
##### b. Install OpenMM Force Fields
OpenMM Force Fields provide additional parameter sets for molecular simulations using OpenMM.

```sh
conda install conda-forge::openmmforcefields --yes
```
##### c. Install OpenEye Toolkits
The OpenEye Toolkits are used for quantum chemistry-based force field parameterization. This requires an OpenEye license, but it is free for academic users.

```sh
conda install openeye::openeye-toolkits --yes
```

##### d. Install espaloma mchine-learned force field
Install espaloma version 0.3.2, which includes the latest parameterization models.

```sh
conda install conda-forge::"espaloma=0.3.2" --yes
```

