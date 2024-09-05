import os
import shutil
from openmm import app
import pdbfixer
import argparse

parser=argparse.ArgumentParser()
parser.add_argument("-pdbid" , required=False)
parser.add_argument("-pdbfile" , required=False)
parser.add_argument("-ph", type=int, help="pH value", required=False, default=7)
parser.add_argument("-mutfile", required=False, help="File with mutations to introduce in the format org_resname, org_resid, mut_resname, chainID")
parser.add_argument("-crotype", required=True, help="Type of chromophore: SYG or TYG")
parser.add_argument("-cresid", required=True, help="Index of the chromophore residue")
parser.add_argument("-protoligmr", help="Please check if your protein complex is a monomer or a dimer", required=False, default="dimer")
parser.add_argument("-fe", help="Perform free energy calculations", default="no")

args = parser.parse_args()

cwd = os.getcwd()
ff = cwd + '/data/forcefields'
mdppath = cwd + '/data/mdp-files/'

if args.crotype not in ["SYG", "TYG"]:
    print("Please provide a valid chromophore type: SYG or TYG")
    exit()

if args.protoligmr not in ["monomer", "dimer"]:
    print("Please provide a valid protein complex type: monomer or dimer")
    exit()
if args.protoligmr == "monomer": num_chains = 1
elif args.protoligmr == "dimer": num_chains = 2

def create_wrkdir(dirname):
    os.mkdir(dirname)
    os.chdir(dirname)

# Load the PDB file
if args.pdbid:
    protname = args.pdbid
    create_wrkdir(protname)
    fixer = pdbfixer.PDBFixer(pdbid=args.pdbid)
    app.PDBFile.writeFile(fixer.topology, fixer.positions, open(args.pdbid+'.pdb', 'w'))
elif args.pdbfile:
    protname = args.pdbfile.split('.')[0]
    create_wrkdir(protname)
    shutil.copy(cwd + '/' + args.pdbfile, args.pdbfile)
    fixer = pdbfixer.PDBFixer(filename=args.pdbfile)
else:
    print("Please provide a PDB file or a PDB ID!")
    exit()

# Special case for pdbid 1GFL: mixing the chromophore residues into a single residue with the name and index provided by the user
if args.pdbid == '1gfl':
    for chain in fixer.topology.chains():
        for residue in chain.residues():
            if residue.name == 'SER' and residue.id == '65':
                for atom in residue.atoms():
                    if atom.name == 'C': atom.name = 'C1'
                    if atom.name == 'CA': atom.name = 'CA1'
                    if atom.name == 'CB': atom.name = 'CB1'
                    if atom.name == 'OG': atom.name = 'OG1'
                residue.name = args.crotype
                residue.id = args.cresid

            if residue.name == 'TYR' and residue.id == '66':
                for atom in residue.atoms():
                    if atom.name == 'N': atom.name = 'N2'
                    if atom.name == 'C': atom.name = 'C2'
                    if atom.name == 'O': atom.name = 'O2'
                    if atom.name == 'CA': atom.name = 'CA2'
                    if atom.name == 'CB': atom.name = 'CB2'
                    if atom.name == 'CG': atom.name = 'CG2'
                residue.name = args.crotype
                residue.id = args.cresid

            if residue.name == 'GLY' and residue.id == '67':
                for atom in residue.atoms():
                    if atom.name == 'N': atom.name = 'N3'
                    if atom.name == 'CA': atom.name = 'CA3'
                residue.name = args.crotype
                residue.id = args.cresid

app.PDBFile.writeFile(fixer.topology, fixer.positions, open('tmp.pdb', 'w'), keepIds=True)

# Check if the pdb file includes a chromophore/flurophore residue
for chainindex in range(num_chains):
    chain = list(fixer.topology.chains())[chainindex]
    for residue in chain.residues():
        if residue.name == args.crotype: cro_res = residue
        elif residue.name == 'CRO': residue.name = args.crotype; cro_res = residue
    if cro_res is None:
        print("Chromophore residue not found in chain", chainindex)
        exit()
    else:
        for atom in cro_res.atoms():
            if atom.name == 'N1': atom.name = 'N'
            if atom.name == 'C3': atom.name = 'C'
            if atom.name == 'O3': atom.name = 'O'
        app.PDBFile.writeFile(fixer.topology, fixer.positions, open('tmp.pdb', 'w'), keepIds=True)
        
# Function to add missing residues and atoms
def add_missing_residues_and_atoms(fixer, ph=args.ph):
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    chains = list(fixer.topology.chains())
    keys = fixer.missingResidues.keys()
    for key in list(keys):
        chain = chains[key[0]]
        if key[1] == 0 or key[1] == len(list(chain.residues())):
            del fixer.missingResidues[key]

    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(ph)

    # Save the fixed PDB file
    app.PDBFile.writeFile(fixer.topology, fixer.positions, open('prepped.pdb', 'w'), keepIds=True)
    return fixer

fixer = add_missing_residues_and_atoms(fixer)
print("Missing residues and atoms added successfully!")

# Apply mutations to generate a variant structure without calling for free energy calculations
# Read mutations from mutations.txt file
def build_mutation(fixer, mut=args.mutfile):
    with open(mut, 'r') as file:
        mutline = [line.split() for line in file]

    # apply mutations
    for oresnm, mresid, mresnm, mchain in mutline:
        mutres = oresnm + '-' + mresid + '-' + mresnm
        fixer.applyMutations(mutres, mchain)
   
    return fixer

if args.fe == "no":
    fixer = build_mutation(fixer)
    # save the mutated structure to pdb file
    app.PDBFile.writeFile(fixer.topology, fixer.positions, open('prepped.pdb', 'w'), keepIds=True)

os.remove('tmp.pdb')
os.chdir(cwd)
print("Preparation of the protein structure is complete!")

