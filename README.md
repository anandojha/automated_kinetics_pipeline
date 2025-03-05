## Automated milestoning pipeline for kinetic and thermodynamic predictions

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

## Instructions for setting up and running the automated milestoning pipeline

Once the conda environment is sett up with necessary package installations, please follow the step-by-step instructions on how to clone, navigate, and set up the milestoning simulation pipeline for kinetic and thermodynamic predictions.

#### 1. First, clone the repository from GitHub

```sh
git clone https://github.com/anandojha/automated_milestoning_pipeline.git
```

#### 2. Once the repository is cloned, change to the working directory to the project folder

```sh
cd automated_milestoning_pipeline
```

#### 3. Navigate to the trypsin-benzamidine folder

```sh
cd trypsin_benzamidine
```

#### 4. Activate the conda environment

```sh
conda activate one_step_kinetics
```

#### 5. Export the OpenEye license

Replace "USERNAME" with the actual home directory name and make sure the oe_license.txt file is stored in the specified location.
```sh
export OE_LICENSE="/home/USERNAME/licenses/oe_license.txt"
```
Once the above steps are successful, read the README.md file in the trypsin_benzamidine project folder for further instructions on running the series of scripts.

There is a separate folder named trypsin_benzamidine_simulation, where these scripts are executed for comparison purposes. Note that these simulations are run for a very short time, meaning that the computed kinetic, thermodynamic, and rate constants may not accurately represent absolute experimental values. The purpose of this execution is to demonstrate the workflow and methodology, but actual simulations should be run for extended durations to obtain scientifically meaningful results. The users can compare results across different simulation times by increasing sampling duration and milestone transitions for better accuracy.


### Relevant GitHub repositories
1. SEEKR2: https://github.com/seekrcentral/seekr2
2. SEEKR2 OpenMM Plugin: https://github.com/seekrcentral/seekr2_openmm_plugin
3. SEEKRTools: https://github.com/seekrcentral/seekrtools
4. QMrebind: https://github.com/seekrcentral/qmrebind

### Relevant milestoning papers
1. Ojha, Anupam Anand, Lane William Votapka, Gary Alexander Huber, Shang Gao, and Rommie Elizabeth Amaro. "An introductory tutorial to the SEEKR2 (Simulation enabled estimation of kinetic rates v. 2) multiscale milestoning software [Article v1. 0]." Living Journal of Computational Molecular Science 5, no. 1 (2023): 2359-2359.
2. Votapka, Lane W., Andrew M. Stokely, Anupam A. Ojha, and Rommie E. Amaro. "SEEKR2: Versatile multiscale milestoning utilizing the OpenMM molecular dynamics engine." Journal of chemical information and modeling 62, no. 13 (2022): 3253-3262.
3. Ojha, Anupam Anand, Lane William Votapka, and Rommie Elizabeth Amaro. "QMrebind: incorporating quantum mechanical force field reparameterization at the ligand binding site for improved drug-target kinetics through milestoning simulations." Chemical Science 14, no. 45 (2023): 13159-13175.
4. Ojha, Anupam Anand, Ambuj Srivastava, Lane William Votapka, and Rommie E. Amaro. "Selectivity and ranking of tight-binding JAK-STAT inhibitors using Markovian milestoning with Voronoi tessellations." Journal of chemical information and modeling 63, no. 8 (2023): 2469-2482.
5. Votapka, Lane W., Benjamin R. Jagger, Alexandra L. Heyneman, and Rommie E. Amaro. "SEEKR: simulation enabled estimation of kinetic rates, a computational tool to estimate molecular kinetics and its application to trypsin–benzamidine binding." The Journal of Physical Chemistry B 121, no. 15 (2017): 3597-3606.
6. Jagger, Benjamin R., Anupam A. Ojha, and Rommie E. Amaro. "Predicting ligand binding kinetics using a Markovian milestoning with voronoi tessellations multiscale approach." Journal of Chemical Theory and Computation 16, no. 8 (2020): 5348-5357.
7. Jagger, Benjamin R., Christopher T. Lee, and Rommie E. Amaro. "Quantitative ranking of ligand binding kinetics with a multiscale milestoning simulation approach." The journal of physical chemistry letters 9, no. 17 (2018): 4941-4948.
8. Votapka, Lane W., and Rommie E. Amaro. "Multiscale estimation of binding kinetics using Brownian dynamics, molecular dynamics and milestoning." PLoS computational biology 11, no. 10 (2015): e1004381.

### Authors and contributors 
The following people have contributed directly to the coding and validation efforts of automated milestoning pipeline for kinetic and thermodynamic predictions  (listed in alphabetical order of first name). The authors would like to thank everyone who has helped or will help improve this project by providing feedback, bug reports, or other comments.
1. Anupam A. Ojha, Flatiron Institute 
2. Lane W. Votapka, UC San Diego 
3. Rommie E. Amaro, UC San Diego 
4. Sonya M. Hanson, Flatiron Institute 

