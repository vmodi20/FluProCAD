# FluProCAD
Computational modelling workflow to model fluorescent proteins and analyze the effect of mutations on thermodynamic properties.
FluProCad is a modelling workflow that helps a user to compute folding and dimerization free energies for fluorescent protien mutants using validated force field based atomistic simulations. Additionally, the workflow can also setup absorption and emission spectra calculations for the FP mutant model. 

To obtain the free energy associated with the mutation, the alchemical transformation is performed with the equilibrium approach, using 21 λ points. The free energy estimates if compute using the slow-growth Bennet’s acceptance ratio estimator.

## Installation
FluProCAD installation includes the following steps:

Download the code from gitlab with `git clone`

Install the dependencies with conda on your local machince to run the workflow from the FluProCAD folder
  	cd FluproCAD
	conda config --add channels conda-forge
	conda env create -f environment.yml
	source activate fpc

This will install the required dependencies. After it you should be able to run FluProCAD

The following sections will discuss the use of the workflow to model folding and dimerization free energy calculations as well as optical spectra for fluorescent protein mutants.

## Step-1: FPC_input.py
This script prepares the input structure for free energy calculations of fluorescent protein mutants. It performs several tasks including fixing missing atoms, introducing mutations, generating Gromacs topologies, and setting up the system for molecular dynamics simulations. A pre-requisite before running the this script is to create file named *mutations.txt* to specificy the mutations to be analysed through this workflow. An example of the format is available in the file *example_mutations.txt*, where the residue PHE at position 64 will be mutated to LEU. Runnig the script without a *.yaml* file will invoke user-prompt to enter some key input details. Optionally, the template *sample_FPC_input.yaml* can be edited an used as input.

### Example Usage:
``` python FPC_input.py ```

or

``` python FPC_input.py sample_FPC_input.yaml ```

The tasks performed by this script are described as follows:

**Prepare Input Structure**
1. **Fix Missing Residue and Atoms**: Fix missing atoms and residues in the input structure.
2. **Introduce Mutations**: Introduce user-provided mutations into the parent structure for structure prediction
3. **Add Hydrogens**: Add hydrogens to the input PDB file based on the specified pH.

**Free Energy MD setup**
1. **Prepare Structures**: Setup monomer, dimer, and unfolded state structure models for the mutant FP.
2. **Generate Hybrid Residues**: Generate hybrid residues and topologies for the mutated residues using the pmx module.
3. **Generate Hybrid Topologies**: Generate Gromacs topologies for the solvated system.
4. **Build Solvated Model**: Build a solvated model of the FP mutant with 0.15 M ionic (Na-Cl) concentration.
5. **Ion Charge Scaling in FE Setup**: Check net charge of the mutated state in free energy setup and scale Na ion charges to avoid having a non-neutral system in the mutated state.
6. **Generate Run Scripts**: Generate files 'runMD-eq.sh' and 'runMD-fe.sh' with Gromacs commands to run the equilibration MD and free-energy simulations. Run them one after another as the FE starting structures are take from 10ns equilibration MD results.
7. **Generate Unfolded State**: Generate a tripeptide GXG with the mutation residue as X to model the unfolded state.

**Structure Prediction MD Setup**
1.  **Prepare Structures**: Generate the FP structure and topologies with mutations from mutations.txt
2.  **Generate Run Scripts**: Generate the file 'runMD-eq.sh' with Gromacs commands to run a 100ns MD simulations.
3.  **Analysis Scripts**: Copy the trajectory analysis scripts. Run this script after the end of 100ns MD to extract predicted structures.

### Step-1 Output: 
- `<protien-name>/`: Generated directory with fixed protein structure
- `prepped.pdb`: Prepared PDB file ready generating free energy setup
- `prepped_structpred_input.pdb`: Prepared PDB file ready for structure prediction setup
- `unfolded/`: Directory containing the unfolded state peptide model for free energy calculations
- `monomer/`: Directory containing the monomer state of the FP mutant for free energy calculations
- `dimer/`: Directory containing the dimer state of the FP mutant for free energy calculations
- `structpred`: Directory containing monomer/dimer state of mutant FP for 100ns MD --> structure prediction analysis

## Step-2: Run MD Simulations
- Run the file 'runMD-eq.sh' under the directories monomer, dimer and unfolded to equilibrate the structures. If running locally with a non-MPI version of gromacs, replace the command 'gmx_mpi' by 
'gmx'. It wo
- Next, to run the free energy calculations, go the directory `slow-TI` where the 21 λ points have been setup to model the slow-growth alchemical transition. Run the file `runMD-fe.sh` present under each of the 21 directories of the λ points. The output files needed for the free energy estimate are the `dhdl.xvg` for each point.
**NOTE**: Running these simulations on the local computer will take more than few days, so it is advised to run them on a cluster which can allow running them in parallel in a much more effiecient time-scale.

## Step-3: Folding and Dimerization free energy estimates (analysis-feMD.py)
When each of the 21 points have compelted the 10ns simulation run (check the end of `md.log` file), we have all the data needed to estimate the foldign and dimerization free energy difference. The script `analysis-feMD.py` will evaluate the free energies using the **alchemlyb** python module and can be run under the main folder of the protien mutant which contains the directories `unfolded`, `monomer`, and `dimer`. For each protein state, the script will estimates the equilibrium free energy change *(ΔG)* and generate a plot of the free energy profile under the `slow-TI` directory. The key results (ΔΔG<sub>folding</sub> and ΔΔG<sub>dimerization</sub>) will be stored in the file **fe-results.dat** 

### Usage:
``` 
python analysis-feMD.py
```
The tasks performed by this script are described as follows:
1. **Estimate Free Energy**: For each state, estimate the free energy difference using the `estimate_equilb_fe` function.
2. **Calculate Free Energy Differences**: Calculate the free energy differences for folding and dimerization.
3. **Mutation Effects**: Determine the effect of mutations on the folding and dimerization free energy and write the results to the file `fe-results.dat`

### Step-2 Output:
- `ti_dhdl.png`: Plot of the free energy profile in the *slow-TI* directory for each state (monomer, dimer, unfolded) of the mutant.
- `fe-results.dat`: Summary of the free energy differences and the effect of mutations.
  
## Step-3: Structure Prediction
The structure and topology input files for this step are already prepared in step-1 and stored under the directory `structpred`. The directory contains both the monomer and dimer state model for the FP mutant, but you can choose the correct state based on the dimerization free energy results or known literature if available. Go the directory with correct oligomeric state of the mutant to run the file `runMD-eq.sh` to simulate the mutant FP for 100ns.

Check the *md.log* files to check if the MD simulations have run successfully. Next, run the analysis scripts to extract a solution structure(s) from the 100ns trajectory that represent the ensemble generated with MD. Here, the script will use a set of RMSD and clustering analysis (GROMOS algorithm) tools available in GROMACS.
**NOTE**: The clustering analysis part of the script can easily take several hours on a local computer depending on the size of the protein system. 

### Usage:
```
analyze-jobs.sh
```
The following tasks are performed at this step:
1. **Root Mean Square Deviation**:
2. **Root Mean Square Fluctuation**:
3. **Radius of Gyration**:
4. **Representative Solution State Structure**:


### Step-3 Output:
- `.xvg`: The *RMSD*, *RMSF*, and *radius of gyration* analysis are stored under a new sub-directory called `MD-analysis`. Examine the ‘.xvg’ plot files using *xmgrace* tool to evaluate the structural stability (RMSD), flexible residue regions (RMSF) and effect of the mutated residue on the compactness of the structure (radius of gyraiton).
- `clusters-FP.log` : The log file summarizes the number of clusters found in the trajectory and extracts a single representative structure for each cluster into the pdb file.
- `clusters-FP.pdb` : Representative structure for each cluster found in the trajectory. Clusters with reasonable population size (>200 frames) and distribution are considered as the representative solution structures for the selected mutation. The first structure is the representative model of largest cluster and the last structure of the smallest.

Visualize the structure with your choice of visualization program (PyMol, VMD,..).


