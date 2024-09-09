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
parser.add_argument("-mut", help='Mutation to introduce in the format: residue, residue ID, mutated_residue, chainID')
parser.add_argument("-ffmut", help='Force field for the mutated residues: amber99sbmut, amber03_gfp', required=False, default='amber99sb-star-ildn-mut_FP')
parser.add_argument("-water", help='Water model', required=False, default='tip3p')
parser.add_argument("-protoligmr", help="Please check if your protein complex is a monomer or a dimer", required=False, default="dimer")
parser.add_argument("-fe", help="Perform free energy calculations (yes/no)", required=False, default="yes" )

args = parser.parse_args()

cwd = os.getcwd()
ff = cwd + '/data/forcefields'
mdppath = cwd + '/data/mdp-files/'
mutfile = cwd + '/' + args.mut

#os.environ['GMXLIB'] = '~/FPC/data/forcefields'

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

# Check the oligomeric state of the protein complex
if args.protoligmr not in ["monomer", "dimer"]:
    print("Please provide a valid protein complex type: monomer or dimer")
    exit()

if args.protoligmr == "monomer":
    #num_chains = 1
    gen_dimer('prepped.pdb')
elif args.protoligmr == "dimer":
    #num_chains = 2
    gen_monomer('prepped.pdb')

for olig in ['monomer', 'dimer']:
    if olig == 'monomer': num_chains = 1
    elif olig == 'dimer': num_chains = 2

    if os.path.exists(olig): shutil.move(olig, olig+'_old')
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

    # Generate temporary GMX toppologies to match residue and atom names before introducing mutations with pmx
    gmx.pdb2gmx(f=olig+'_capped.pdb', o='conf.pdb', p='topol.top', ff='amber99sb-star-ildn-mut_FP', water=args.water, other_flags='-ignh')

    if args.fe == "no":
        topol = 'topol.top'
        pass
    else:
        # Introduce mutations and generate hybrid residues
        with open(mutfile, 'r') as file:
            mutline = [line.split() for line in file]

        pmutref = Model('conf.pdb', renumber_residues=False)

        # mutate residues of chain A and B if protein is a dimer, else over chain A only
        for oresnm, mresid, mresnm, tmp in mutline:
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
    gmx.grompp(f=mdppath+'min.mdp', c='solv.pdb', p=topol, o='ions.tpr', maxwarn='1')
    gmx.genion(s='ions.tpr', o='ions.pdb', p=topol, neutral=True, conc=0.15)

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
for oresnm, mresid, mresnm, mchain in mutline:
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
    gmx.editconf(f='conf.pdb', o='box.pdb', bt='cubic', d=1.0)
    gmx.solvate(cp='box.pdb', cs='spc216.gro', o='solv.pdb', p=topol)
    gmx.grompp(f=mdppath+'min.mdp', c='solv.pdb', p=topol, o='ions.tpr', maxwarn='1')
    gmx.genion(s='ions.tpr', o='ions.pdb', p=topol, neutral=True, conc=0.15)

    os.remove('GGG.pdb')
    # remove all the temporary files generated during the mutation process
    for file in os.listdir('.'):
        if file.startswith("#"):
            os.remove(file)
    os.chdir('..')

print("#"*50)
print("Generated the unfolded state!")

os.chdir(cwd)