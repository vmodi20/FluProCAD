import os
import shutil
from openmm import app
import pdbfixer
from pmx import *
import argparse

# Parse command line arguments
parser = argparse.ArgumentParser(description='Prepare Gromacs format topology and structure files with aa mutations to setup the free energy MD simulations')
parser.add_argument("-proteindir", help='Directory with the prepped protein structure', required=True)
parser.add_argument("-mut", help='Mutation to introduce in the format: residue, residue ID, mutated_residue, chainID')
parser.add_argument("-ffmut", help='Force field for the mutated residues: amber99sbmut, amber03_gfp', required=False, default='amber99sb-star-ildn-mut_FP')
parser.add_argument("-water", help='Water model', required=False, default='tip3p')
parser.add_argument("-protoligmr", help="Please check if your protein complex is a monomer or a dimer", required=False, default="dimer")
parser.add_argument("--fe", help="Perform free energy calculations", action='store_true')

args = parser.parse_args()

cwd = os.getcwd()
ff = cwd + '/data/forcefields'
mdppath = cwd + '/data/mdp-files/'
mutfile = cwd + '/' + args.mut

#os.environ['GMXLIB'] = '~/FPC/data/forcefields'

os.chdir(args.proteindir)

# Check the oligomeric state of the protein complex
if args.protoligmr not in ["monomer", "dimer"]:
    print("Please provide a valid protein complex type: monomer or dimer")
    exit()
if args.protoligmr == "monomer": num_chains = 1
elif args.protoligmr == "dimer": num_chains = 2

# Load the PDB file
p = Model('prepped.pdb', renumber_residues=False)

# Add N-terminal cap ACE and C-terminal cap NME to the protein complex
for chainindex in range(num_chains):
    chain = list(p.chains)[chainindex]
    chain.add_nterm_cap()
    chain.add_cterm_cap()

p.write('capped.pdb')

# Generate temporary GMX toppologies to match residue and atom names before introducing mutations with pmx
gmx.pdb2gmx(f='capped.pdb', o='conf.pdb', p='topol.top', ff='amber99sb-star-ildn-mut_FP', water=args.water, other_flags='-ignh')

# Introduce mutations and generate hybrid residues
with open(mutfile, 'r') as file:
    mutline = [line.split() for line in file]

if args.fe is False:
    topol = 'topol.top'
    pass
else:
    pmutref = Model('conf.pdb', renumber_residues=False)
    for oresnm, mresid, mresnm, mchain in mutline:
        pmut = mutate(m=pmutref, mut_resid=int(mresid), mut_resname=mresnm, mut_chain=mchain, ff='amber99sb-star-ildn-mut_FP')
        pmutref = pmut

    pmutref.write('mutated.pdb')
    # Generate GMX toplogies
    gmx.pdb2gmx(f='mutated.pdb', o='conf.pdb', p='topol.top', ff='amber99sb-star-ildn-mut_FP', water=args.water)

    if args.protoligmr == "monomer":
        gmxtop = Topology('topol.top', ff='amber99sb-star-ildn-mut_FP')
        pmxtop, _ = gen_hybrid_top(gmxtop)
        pmxtop.write('newtopol.top')
    elif args.protoligmr == "dimer":
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
print("Generated hybrid topology files!")
print("Proceeding with the cubic box and solvation step...")

# Setup the solvated system with NA+ and CL- ions at 0.15 M concentration
gmx.editconf(f='conf.pdb', o='box.pdb', bt='cubic', d=1.0)
gmx.solvate(cp='box.pdb', cs='spc216.gro', o='solv.pdb', p=topol)
gmx.grompp(f=mdppath+'min.mdp', c='solv.pdb', p=topol, o='ions.tpr', maxwarn='1')
gmx.genion(s='ions.tpr', o='ions.pdb', p=topol, neutral=True, conc=0.15)

print("#"*50)
print("Generated solvated system with ions!")
print("Proceeding with building the unfolded state...")

# Generate a tripeptide GXG with the mutation residue as X to model the unfolded state
def gen_unfolded_peptide(mutations):
    os.mkdir('unfolded')
    os.chdir('unfolded')
    for oresnm, mresid, mresnm, mchain in mutations:
        ocode = library._one_letter[oresnm]
        mcode = library._one_letter[mresnm]
        os.mkdir(ocode + '2' + mcode) 
        os.chdir(ocode + '2' + mcode)
        gxg_pdbpath = cwd + '/data/tripeptide_GGG.pdb' 
        shutil.copy(gxg_pdbpath, 'GGG.pdb') 
        xres = 'GLY-2-' + oresnm
        g = pdbfixer.PDBFixer(filename='GGG.pdb')
        g.applyMutations([xres], 'A')
        g.findMissingResidues()
        g.findMissingAtoms()
        g.addMissingAtoms()
        g.addMissingHydrogens(7)
        app.PDBFile.writeFile(g.topology, g.positions, open('unfolded.pdb', 'w'), keepIds=True) 

        gxg = Model('unfolded.pdb', rename_atoms=True)
        list(gxg.chains)[0].add_nterm_cap()
        list(gxg.chains)[0].add_cterm_cap()
        gxg.write('unfolded.pdb')
        gxgmut = mutate(m='conf.pdb', mut_resid=3, mut_resname=mresnm, ff='amber99sb-star-ildn-mut_FP')
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
        os.chdir('..')

gen_unfolded_peptide(mutline)
print("#"*50)
print("Generated the unfolded state!")

os.chdir(cwd)