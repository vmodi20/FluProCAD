import os
import shutil
from openmm import app
import pdbfixer
from pmx import *
import pymol.cmd as pml
import argparse

# Parse user input info as arguments
parser = argparse.ArgumentParser(description='Prepare Gromacs format topology and structure files with aa mutations to setup the free energy MD simulations')
parser.add_argument("-proteindir", help='Directory with the prepped protein structure', required=True)
parser.add_argument("-mutfile", help='file with mutationd to introduce in the format: residue, residue ID, mutated_residue')
parser.add_argument("-ffmut", help='Force field for the mutated residues: amber99sbmut, amber03_gfp', required=False, default='amber99sb-star-ildn-mut_FP')
parser.add_argument("-water", help='Water model', required=False, default='tip3p')
parser.add_argument("-protoligmr", help="Please check if your protein complex is a monomer or a dimer", required=False, default="dimer")
parser.add_argument("-fe", help="Perform free energy calculations (yes/no)", required=False, default="yes" )
parser.add_argument("-structpred", help="Perform solution state structure prediction", required=False, default="yes")

args = parser.parse_args()

cwd = os.getcwd()
ff = cwd + '/data/forcefields'
mdppath = cwd + '/data/mdp-files/'
mutfile = cwd + '/' + args.mutfile

os.environ['GMXLIB'] = cwd + '/data/forcefields'

os.chdir(args.proteindir)

# Generate dimer structure for monomeric input structures
def gen_dimer(protfile):
    shutil.copy(cwd+'/data/1gfl.pdb', '1gfl.pdb')
    pml.load(protfile, "target")
    pml.load("1gfl.pdb", "ref")
    pml.select("refchainA", "ref and chain A")
    pml.select("refchainB", "ref and chain B")
    pml.create("targetchainA", "target")
    pml.create("targetchainB", "target")
    pml.super("targetchainA", "refchainA")
    pml.super("targetchainB", "refchainB")
    pml.copy("dimer", "targetchainA")
    pml.copy_to("dimer", "targetchainB")
    pml.select("water", "dimer and resname HOH")
    pml.create("solv", "water")
    pml.remove("water")
    pml.copy_to("dimer", "solv")
    pml.save("monomer.pdb", "targetchainA")
    pml.save("dimer.pdb", "dimer")
    os.remove("1gfl.pdb")

# Generate monomer structure for dimeric input structures
def gen_monomer(protfile):
    pml.load(protfile, "target")
    pml.select("chainA", "target and chain A")
    pml.create("monomer", "chainA")
    pml.save("monomer.pdb", "monomer")
    shutil.copy("prepped.pdb", "dimer.pdb")

#check if command gmx exists
def check_gmx():
    if shutil.which('gmx') is None:
        print("Gromacs is not installed or not in the PATH. Please install Gromacs and add it to the PATH.")
        exit()
    else:
        print("Gromacs is installed and in the PATH.")
        print("#"*50)

# Function to create bash file with gmx commands for equilibration MD step:
def runMDeq():
    bashfile = 'runMD-eq.sh'
    with open(bashfile, 'w') as file:
        file.write('#!/bin/bash\n')
        file.write('\n# Run minimization, NVT, and NPT simulations using Gromacs\n')
        file.write('gmx_mpi grompp -f min.mdp -c ions.pdb -p newtopol.top -o min.tpr -maxwarn 1\n')
        file.write('gmx_mpi mdrun -deffnm min\n')
        file.write('gmx_mpi grompp -f nvt.mdp -c min.gro -p newtopol.top -r min.gro -o nvt.tpr -maxwarn 1\n')
        file.write('gmx_mpi mdrun -deffnm nvt\n')
        file.write('gmx_mpi grompp -f npt.mdp -c nvt.gro -p newtopol.top -r nvt.gro -o npt.tpr -maxwarn 1\n')
        file.write('gmx_mpi mdrun -deffnm npt\n')
        file.write('gmx_mpi grompp -f md.mdp -c npt.gro -p newtopol.top -r npt.gro -o topol.tpr -maxwarn 1\n')
        file.write('gmx_mpi mdrun -s topol\n')
        file.write('echo "Protein System " | gmx_mpi trjconv -f traj.trr -s topol.tpr -pbc mol -center -o confout-MDeq.gro')
        file.close()

# Function to create bash file with gmx commands for free energy MD step:
def runMDfe():
    bashfile = 'runMD-fe.sh'
    shutil.copy(cwd+'/analysis-feMD.py', 'analysis-feMD.py')
    os.mkdir('slow-TI')
    os.chdir('slow-TI')
    for i in range(21):
        os.mkdir(str(i))
        for file in os.listdir(mdppath):
            if file.endswith('-fe.mdp'):
                # copy the file to the $i directory
                shutil.copy(mdppath+file, str(i))
                with open(str(i)+'/'+file, 'r') as f:
                    filedata = f.read()
                filedata = filedata.replace('@lambda@', str(i))
                with open(str(i)+'/'+file, 'w') as f:
                    f.write(filedata)
            
        with open(str(i)+'/'+bashfile, 'w') as file:
            file.write('#!/bin/bash\n')
            file.write('\n# For each lambda point, run a minimization, NVT, and NPT simulation\n')
            file.write('gmx_mpi grompp -f min-fe.mdp -c ../../confout-MDeq.gro -p ../../newtopol.top -o min-fe.tpr -maxwarn 1\n')
            file.write('gmx_mpi mdrun -deffnm min-fe\n')
            file.write('gmx_mpi grompp -f nvt-fe.mdp -c min-fe.gro -p ../../newtopol.top -r min-fe.gro -o nvt-fe.tpr -maxwarn 1\n')
            file.write('gmx_mpi mdrun -deffnm nvt-fe\n')
            file.write('gmx_mpi grompp -f npt-fe.mdp -c nvt-fe.gro -p ../../newtopol.top -r nvt-fe.gro -o npt-fe.tpr -maxwarn 1\n')
            file.write('gmx_mpi mdrun -deffnm npt-fe\n')
            file.write('gmx_mpi grompp -f md-fe.mdp -c npt-fe.gro -p ../../newtopol.top -r npt-fe.gro -o topol.tpr -maxwarn 1\n')
            file.write('gmx_mpi mdrun -s topol\n')
            file.close()
    os.chdir('..')

# Check the net charge of system in mutated state and scale the charge of ions to neutralize the system
def check_stateB_charge(pdbfile, topfile, olig):

    aa_charg_dict = {
        'ALA': 0,
        'ARG': 1,
        'ASN': 0,
        'ASP': -1,
        'CYS': 0,
        'GLN': 0,
        'GLU': -1,
        'GLH': 0,
        'GLY': 0,
        'HIS': 0,
        'HIE': 0,
        'HIP': 1,
        'HID': 0,
        'ILE': 0,
        'LEU': 0,
        'LYS': 1,
        'MET': 0,
        'PHE': 0,
        'PRO': 0,
        'SER': 0,
        'THR': 0,
        'TRP': 0,
        'TYR': 0,
        'VAL': 0,
        }

    # Read the mutations file and calculate the net charge of the mutated residues
    with open(mutfile, 'r') as file:
        mutline = [line.split() for line in file]
    net_charge_oresnm = 0
    net_charge_mresnm = 0
    for oresnm, mresid, mresnm in mutline:
        net_charge_mresnm += aa_charg_dict[mresnm]
        net_charge_oresnm += aa_charg_dict[oresnm]
    net_charge = net_charge_mresnm - net_charge_oresnm
    if olig == 'dimer':
        net_charge = 2*net_charge
    print("#"*50)
    print("Net change in charge with the introduced mutations: ", net_charge)

    NAscale = 0 
    if net_charge !=  0:
        print("Scaling the charge of ions to neutralize the system...")
        NAscale_tot = abs(net_charge*10) # Number of Na ions to scale

        # count the number of NA and CL ions in the system
        with open(pdbfile, 'r') as file:
            lines = file.readlines()
        NAtotal = 0
        for line in lines:
            if line.startswith('ATOM') and line[17:20] == ' NA':
                NAtotal += 1
        print("Number of NA ions in the system: ", NAtotal)
        if NAscale_tot > NAtotal:
            print("WARNING : The number of Na ions in the system is less than the required number to neutralize the system!")
            print("Please consider adding the charge changing mutations in two different set of simulations as the system")
            print("Running the simulations with net charge in the mutated state can affect the thermodynamics of the system.")
            exit()
        
        print("Scaling ", NAscale_tot, " ions to neutralize the system...")

        if net_charge < 0 :
            newNA = 'NAP'
        elif net_charge > 0 :
            newNA = 'NAN'

        pdblines = []
        for line in lines:
            if line.startswith('ATOM') and line[17:20] == ' NA':
                if NAscale < NAscale_tot:
                    line = line[:17] + newNA + line[20:]
                    NAscale += 1
            pdblines.append(line)
        
        shutil.copy(pdbfile, 'old_'+pdbfile)
        with open(pdbfile, 'w') as file:
            file.writelines(pdblines)
            file.close()

        updatedNA = NAtotal - NAscale_tot

        #Update the topology file to include scaled NA ions
        with open(topfile, 'r') as file:
            lines = file.readlines()

        new_lines = []
        for line in lines:
            if line.startswith('NA'):
                line = f'NA   {updatedNA}\n'
            new_lines.append(line)
            if line.startswith('NA'): 
                new_lines.append(f'{newNA} {NAscale_tot}\n')
        
        shutil.copy(topfile, 'old_'+topfile)
        with open(topfile, 'w') as file:
            file.writelines(new_lines)

        print("The system is now neutralized!")
    else:
        print("The system is already neutralized!")
    print("#"*50)

# Check the oligomeric state of the protein complex
if args.protoligmr not in ["monomer", "dimer"]:
    print("Please provide a valid protein complex type: monomer or dimer")
    exit()

if args.protoligmr == "monomer":
    gen_dimer('prepped.pdb')
elif args.protoligmr == "dimer":
    gen_monomer('prepped.pdb')

for olig in ['monomer', 'dimer']:
    if olig == 'monomer': num_chains = 1
    elif olig == 'dimer': num_chains = 2

    if os.path.exists(olig):
        if os.path.exists(olig+'_old'):
            shutil.rmtree(olig+'_old')
            shutil.move(olig, olig+'_old')
    os.makedirs(olig)
    shutil.move(olig+'.pdb', olig)
    os.chdir(olig)

    # Add N-terminal cap ACE and C-terminal cap NME to the protein complex
    p = Model(olig+'.pdb', renumber_residues=False)
    for chainindex in range(num_chains):
        chain = list(p.chains)[chainindex]
        chain.add_nterm_cap()
        chain.add_cterm_cap()
    p.write(olig+'_capped.pdb')

    # Generate temporary GMX topologies to match residue and atom names before introducing mutations with pmx
    gmx.pdb2gmx(f=olig+'_capped.pdb', o='conf.pdb', p='topol.top', ff='amber99sb-star-ildn-mut_FP', water=args.water, other_flags='-ignh')

    if args.fe == "yes":
        # Introduce mutations and generate hybrid residues
        with open(mutfile, 'r') as file:
            mutline = [line.split() for line in file]

        pmutref = Model('conf.pdb', renumber_residues=False)

        # mutate residues of chain A and B if protein is a dimer, else over chain A only
        for oresnm, mresid, mresnm in mutline:
            mchain = "A"
            pmut = mutate(m=pmutref, mut_resid=int(mresid), mut_resname=mresnm, mut_chain=mchain, ff='amber99sb-star-ildn-mut_FP')
            if olig == 'dimer':
                mchain = "B"
                pmut = mutate(m=pmut, mut_resid=int(mresid), mut_resname=mresnm, mut_chain=mchain, ff='amber99sb-star-ildn-mut_FP')
            pmutref = pmut
            
        pmutref.write(olig+'_mutated.pdb')
        # Generate GMX toplogies
        gmx.pdb2gmx(f=olig+'_mutated.pdb', o='conf.pdb', p='topol.top', ff='amber99sb-star-ildn-mut_FP', water=args.water)

        if olig == "monomer":
            gmxtop = Topology('topol.top', ff='amber99sb-star-ildn-mut_FP')
            pmxtop, _ = gen_hybrid_top(gmxtop)
            pmxtop.write('newtopol.top')
        elif olig == "dimer":
            shutil.copy('topol.top', 'newtopol.top')
            for itp in ['topol_Protein_chain_A.itp', 'topol_Protein_chain_B.itp']:
                gmxitp = Topology(itp, ff='amber99sb-star-ildn-mut_FP')
                pmxitp, _ = gen_hybrid_top(gmxitp)
                pmxitp.write('new'+itp)
                # find string itp in newtopol.top and replace it with newitp
                with open('newtopol.top', 'r') as file:
                    filedata = file.read()
                filedata = filedata.replace(itp, 'new'+itp)
                with open('newtopol.top', 'w') as file:
                    file.write(filedata)

        topol = 'newtopol.top'

        print("#"*50)
        print("Generated hybrid topology files for the "+olig+" structure!")
        print("Next, solvating the complex in a cubic box with 0.15M NaCl concentration...")

    # Setup the solvated system with NA+ and CL- ions at 0.15 M concentration
    gmx.editconf(f='conf.pdb', o='box.pdb', bt='cubic', d=1.0)
    gmx.solvate(cp='box.pdb', cs='spc216.gro', o='solv.pdb', p=topol)
    gmx.grompp(f=mdppath+'min.mdp', c='solv.pdb', p=topol, o='ions.tpr', maxwarn='2')
    gmx.genion(s='ions.tpr', o='ions.pdb', p=topol, neutral=True, conc=0.15)

    print("#"*50)
    print("Checking the total charge of system in B-state and scaling the charge of ions accordingly....")
    
    check_stateB_charge('ions.pdb', 'newtopol.top', olig)

    #copy all .mdp files from the mdpfiles path to current directory
    shutil.copytree(mdppath, os.getcwd(), dirs_exist_ok=True)
    print("Generating a bash files 'runMD-eq.sh' and 'runMD-fe.sh' with Gromacs commands to run the equilibration MD and Free-energy simulations")
    runMDeq()
    runMDfe()

    # remove all the temporary files generated during the mutation process
    for file in os.listdir('.'):
        if file.startswith("#"):
            os.remove(file)
    os.chdir('..')

print("#"*50)
print("Generated solvated system with ions!")
print("Proceeding to build the unfolded state peptide...")

# Generate a tripeptide GXG with the mutation residue as X to model the unfolded state
os.makedirs('unfolded')
os.chdir('unfolded')
for oresnm, mresid, mresnm in mutline:
    ocode = library._one_letter[oresnm]
    mcode = library._one_letter[mresnm]
    os.mkdir(ocode + '2' + mcode) 
    os.chdir(ocode + '2' + mcode)
    gxg_pdbpath = cwd + '/data/tripeptide_GGG.pdb' 
    shutil.copy(gxg_pdbpath, 'GGG.pdb') 
    xres = 'GLY-2-' + oresnm
    g = pdbfixer.PDBFixer(filename='GGG.pdb')
    if oresnm != 'GLY': g.applyMutations([xres], 'A')
    g.findMissingResidues()
    g.findMissingAtoms()
    g.addMissingAtoms()
    g.addMissingHydrogens(7)
    app.PDBFile.writeFile(g.topology, g.positions, open('unfolded.pdb', 'w'), keepIds=True) 

    gxg = Model('unfolded.pdb', rename_atoms=True)
    list(gxg.chains)[0].add_nterm_cap()
    list(gxg.chains)[0].add_cterm_cap()
    gxg.write('unfolded.pdb')
    gmx.pdb2gmx(f='unfolded.pdb', o='conf.pdb', p='topol.top', ff='amber99sb-star-ildn-mut_FP', water=args.water, other_flags='-ignh')
    gxgnew = Model('conf.pdb', renumber_residues=False)
    gxgmut = mutate(m=gxgnew, mut_resid=2, mut_resname=mresnm, ff='amber99sb-star-ildn-mut_FP')
    gxgmut.write('unfolded_'+ocode+'2'+mcode+'_mut.pdb')

    # Generate GMX toplogies for hybrid residue and solvate the system without ions 
    gmx.pdb2gmx(f='unfolded_'+ocode+'2'+mcode+'_mut.pdb', o='conf.pdb', p='topol.top', ff='amber99sb-star-ildn-mut_FP', water=args.water)
    unftop = Topology('topol.top', ff='amber99sb-star-ildn-mut_FP')
    unfpmxtop, _ = gen_hybrid_top(unftop)
    unfpmxtop.write('newtopol.top')
    topol = 'newtopol.top'
    gmx.editconf(f='conf.pdb', o='box.pdb', bt='cubic', other_flags='-box 5')
    gmx.solvate(cp='box.pdb', cs='spc216.gro', o='solv.pdb', p=topol)
    gmx.grompp(f=mdppath+'min.mdp', c='solv.pdb', p=topol, o='ions.tpr', maxwarn='1')
    gmx.genion(s='ions.tpr', o='ions.pdb', p=topol, neutral=True, conc=0.15)

    check_stateB_charge('ions.pdb', 'newtopol.top', 'unfolded')

    #copy all .mdp files from the mdpfiles path to current directory
    shutil.copytree(mdppath, os.getcwd(), dirs_exist_ok=True)
    print("Generating a bash files 'runMD-eq.sh' and 'runMD-fe.sh' with Gromacs commands to run the equilibration MD and Free-energy simulations")
    runMDeq()
    runMDfe()

    os.remove('GGG.pdb')
    # remove all the temporary files generated during the mutation process
    for file in os.listdir('.'):
        if file.startswith("#"):
            os.remove(file)
    os.chdir('..')

os.chdir('..')
print("#"*50)
print("Generated the unfolded state!")

def structpred_input():
    os.makedirs('structpred')
    os.chdir('structpred')
    shutil.move(cwd+'/'+args.proteindir+'/prepped_structpred_input.pdb', 'prepped_structpred_input.pdb')

    if args.protoligmr == "monomer":
        gen_dimer('prepped_structpred_input.pdb')
    elif args.protoligmr == "dimer":
        gen_monomer('prepped_structpred_input.pdb')

    for olig in ['monomer', 'dimer']:
        if olig == 'monomer': num_chains = 1
        elif olig == 'dimer': num_chains = 2

        if os.path.exists(olig):
            if os.path.exists(olig+'_old'): shutil.rmtree(olig+'_old')
            shutil.move(olig, olig+'_old')
        os.makedirs(olig)
        shutil.move(olig+'.pdb', olig)
        os.chdir(olig)

        p = Model(olig+'.pdb', renumber_residues=False)
        for chainindex in range(num_chains):
            chain = list(p.chains)[chainindex]
            chain.add_nterm_cap()
            chain.add_cterm_cap()
        p.write(olig+'_capped.pdb')

        # Generate GMX topologies for the mutated FP 
        print("Generating GMX topologies for mutated FP for the structure prediction MD simulations...")
        gmx.pdb2gmx(f=olig+'_capped.pdb', o='conf.pdb', p='topol.top', ff='amber99sb-star-ildn-mut_FP', water=args.water, other_flags='-ignh')
        print("Next, solvating the complex in a cubic box with 0.15M NaCl concentration...")
        gmx.editconf(f='conf.pdb', o='box.pdb', bt='cubic', d=1.0)
        gmx.solvate(cp='box.pdb', cs='spc216.gro', o='solv.pdb', p='topol.top')
        gmx.grompp(f=mdppath+'min.mdp', c='solv.pdb', p='topol.top', o='ions.tpr', maxwarn='2')
        gmx.genion(s='ions.tpr', o='ions.pdb', p='topol.top', neutral=True, conc=0.15)
        print(" DONE!")
        shutil.copytree(mdppath, os.getcwd(), dirs_exist_ok=True)
        print("Generating a bash files 'runMD-eq.sh' with Gromacs commands to run 100ns equilibration MD simulations")
        runMDeq()
        print(" DONE!")
        print(" Run the bash file 'runMD-eq.sh' locally or on a cluster to start equilibration MD simulations.")

        shutil.copy(cwd+'/analysis-MDeq.sh', 'analysis-MDeq.sh')
        shutil.copy(cwd+'/qmmm-inpsetup-firefly.sh', 'qmmm-inpsetup-firefly.sh')
        for file in os.listdir('.'):
            if file.startswith("#"): os.remove(file)
            if file.endswith('-fe.mdp'): os.remove(file)
        
        with open('md.mdp', 'r') as f:
            filedata = f.read()
        filedata = filedata.replace('= 5000000', '= 50000000')
        with open('md.mdp', 'w') as f:
            f.write(filedata)
        
        os.chdir('..')
    os.chdir('..')

if args.structpred == "yes":
    print("#"*50)
    print("# Proceeding to generate input topologies for structure prediction of the protein complex...")

    structpred_input()
    print("#"*50)    
    print("# Both monomer and dimer structures will be generated but run MD only for the correct oligomeric state.")
    print("# Determine the correct oligomeric state based on dimerization free energy results or know literature")
    print("#"*50)

os.chdir(cwd)
