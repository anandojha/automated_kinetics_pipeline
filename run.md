## Automated simulation pipeline for kinetic and thermodynamic predictions

The following series of scripts automates the process of running milestoning simulations to compute kinetic and thermodynamic properties for receptor-ligand interactions. Following this workflow, users can efficiently obtain binding and unbinding rates (kon, koff) for a given protein-ligand complex. To execute the pipeline, users must first install the necessary dependencies (refer to the README.md file for setup instructions) and ensure they work in the one_step_kinetics conda environment. Before running the scripts, the user must specify the correct paths to SEEKR2 and SEEKRTools in the relevant scripts. These paths need to be modified in:

 - run_initial_simulation.py (for SEEKRTools path)
 - run_milestoning.py (for SEEKR2 path)
 - run_analysis.py (for SEEKR2 path)


```
##### f. Compile the software
```sh
make -j 4 all
```

##### g. Return to the home directory 
```sh
cd ~
```

##### g. Clean up unnecessary files
```sh
rm -rf browndye2.tar.gz
```
