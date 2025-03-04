## Automated simulation pipeline for kinetic and thermodynamic predictions

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
#### 5. Verify SEEKR2 OpenMM plugin installation
Run the following command to check if the SEEKR2 OpenMM plugin is correctly installed. If no error message appears, the installation was successful.
```sh
python -c "import seekr2plugin"
```

#### 6. Install Git
Git is required to clone repositories. Ensure it is installed using:

```sh
conda install conda-forge::git --yes
```

#### 7. Install the SEEKR2 package
SEEKR2 is the core package required for performing milestoning simulations. We will download it from GitHub, install it, and verify that it works correctly.

##### Before proceeding, go to the home directory
```sh
cd ~
```
##### Clone the SEEKR2 repository from GitHub into the current directory
```sh
git clone https://github.com/seekrcentral/seekr2.git
```
##### Navigate into the seekr2 directory, where the cloned repository is located
```sh
cd seekr2
```
##### Install SEEKR2 using pip, making it accessible as a Python package in the conda environment
```sh
python -m pip install .
```
##### Run tests to verify that SEEKR2 has been installed correctly and is functioning as expected
```sh
pytest
```
##### Return to the home directory (~), ensuring a clean workspace before proceeding to the next steps
```sh
cd ~
```

#### 8. Install the SEEKRTools package
SEEKRTools is a companion package to SEEKR2 that provides utilities for preparing and facilitating multiscale milestoning simulations.

##### Before proceeding, go to the home directory
```sh
cd ~
```
##### Clone the SEEKRTools repository from GitHub into the current directory
```sh
git clone https://github.com/seekrcentral/seekrtools.git
```
##### Navigate into the seekrtools directory where the cloned repository is located
```sh
cd seekrtools
```
##### Install SEEKRTools using pip, making it accessible as a Python package in the conda environment
```sh
python -m pip install .
```
##### Run tests to verify that SEEKRTools has been installed correctly and is functioning as expected
```sh
pytest
```
##### Return to the home directory (~), ensuring a clean workspace before proceeding to the next steps
```sh
cd ~
```

#### 9. Install the OpenFF toolkit
The Open Force Field (OpenFF) Toolkit is required to assign and manipulate molecular mechanics parameters. 

```sh
conda install conda-forge::openff-toolkit --yes
```

#### 10. Install OpenEye toolkits
The OpenEye Toolkits are used for quantum chemistry-based force field parameterization. The OpenEye toolkits require a valid OpenEye academic license, free for academic users but must be obtained directly from https://www.eyesopen.com/academic-licensing.

##### Install OpenEye toolkits
Run the following command to install OpenEye toolkits via conda:

```sh
conda install openeye::openeye-toolkits --yes
```
##### Obtain and place the license file
After obtaining an OpenEye academic license, save the provided oe_license.txt file in a secure location on your system.
For example, you may place it in:
```sh
/home/USERNAME/licenses/oe_license.txt
```
##### Add the license to your environment
To ensure that OpenEye toolkits can find the license file at runtime, export the license path by adding the following line to your ~/.bashrc. 

```sh
export OE_LICENSE="/home/USERNAME/licenses/oe_license.txt"
```
##### Source ~/.bashrc
To apply this change immediately in the current terminal session, run:

```sh
source ~/.bashrc
```

#### 11. Install OpenMM forcefields
OpenMM Force Fields provide additional parameter sets for molecular simulations using OpenMM.

```sh
conda install conda-forge::openmmforcefields --yes
```

#### 12. Install espaloma machine-learned force field
Install espaloma version 0.3.2, which includes the latest parameterization models.

```sh
conda install conda-forge::"espaloma=0.3.2" --yes
```

#### 13. Install BrownDye2
BrownDye2 is a package used for Brownian dynamics (BD) simulations, which are needed to compute association rate constants. If you plan to run BD simulations, follow these installation steps. Some of these steps require sudo privileges (administrator access). If you do not have sudo access, contact your system administrator.

##### Before proceeding, go to the home directory
```sh
cd ~
```
##### Install required dependencies
BrownDye2 requires several system libraries for compilation. Install them using:

```sh
sudo apt-get install libexpat1 make apbs liblapack-dev
```

```sh
sudo apt-get install ocaml ocaml-native-compilers
```

```sh
sudo apt-get install libexpat1-dev
```

##### Download the latest BrownDye2 source code
```sh
wget https://browndye.ucsd.edu/downloads/browndye2.tar.gz
```
##### Extract the downloaded archive
```sh
tar xvfz browndye2.tar.gz
```
##### Navigate into the browndye2 directory 
```sh
cd browndye2
```
##### Compile the software
```sh
make -j 4 all
```

##### Return to the home directory 
```sh
cd ~
```

##### Clean up unnecessary files
```sh
rm -rf browndye2.tar.gz
```
