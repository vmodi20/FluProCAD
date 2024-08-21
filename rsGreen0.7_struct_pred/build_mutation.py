import sys
import os
import argparse

from modeller import *
from modeller.optimizers import molecular_dynamics, conjugate_gradients
from modeller.automodel import autosched
from modeller.scripts import complete_pdb

###########################################################################################################################
# DEFINE COMMAND LINE OPTIONS AND GENERATE --help AND ERROR HANDLING
parser=argparse.ArgumentParser()
parser.add_argument("-modelname")
parser.add_argument("-suffix")
parser.add_argument("-respos", nargs="+", type=int)
parser.add_argument("-restyp", nargs="+")
parser.add_argument("-chain", nargs="+")

# parse the command line 
args = parser.parse_args()

###########################################################################################################################

# RESIDUE DATABASE
residb=['ALA',
    'ARG',
    'ASN',
    'ASX',
    'ASP',
    'CYS',
    'CSS',
    'GLU',
    'GLX',
    'GLN',
    'GLY',
    'HIS',
    'HSE',
    'HSP',
    'ILE',
    'LEU',
    'LYS',
    'MET',
    'PHE',
    'PRO',
    'SER',
    'THR',
    'TRP',
    'TYR',
    'VAL']

###########################################################################################################################

#ENERGY MINIMIZATION
def optimize(atmsel, sched):
    #Conjugate gradient
    for step in sched:
        step.optimize(atmsel, max_iterations=200, min_atom_shift=0.001)
    #md
    refine(atmsel)
    cg = conjugate_gradients()
    cg.optimize(atmsel, max_iterations=200, min_atom_shift=0.001)


# MOLECULAR DYNAMICS
def refine(atmsel):
    # at T=1000, max_atom_shift for 4fs is cca 0.15 A.
    md = molecular_dynamics(cap_atom_shift=0.39, md_time_step=4.0,
                            md_return='FINAL')
    init_vel = True
    for (its, equil, temps) in ((200, 20, (150.0, 250.0, 400.0, 700.0, 1000.0)),
                                (200, 600,
                                 (1000.0, 800.0, 600.0, 500.0, 400.0, 300.0))):
        for temp in temps:
            md.optimize(atmsel, init_velocities=init_vel, temperature=temp,
                         max_iterations=its, equilibrate=equil)
            init_vel = False

#USE HOMOLOGS AND DIHEDRAL LIBRARY FOR DIHEDRAL ANGLE RESTRAINTS
def make_restraints(mdl1, aln):
   rsr = mdl1.restraints
   rsr.clear()
   s = selection(mdl1)
   for typ in ('stereo', 'phi-psi_binormal'):
       rsr.make(s, restraint_type=typ, aln=aln, spline_on_site=True)
   for typ in ('omega', 'chi1', 'chi2', 'chi3', 'chi4'):
       rsr.make(s, restraint_type=typ+'_dihedral', spline_range=4.0,
                spline_dx=0.3, spline_min_points = 5, aln=aln,
                spline_on_site=True)

###########################################################################################################################
log.verbose()

# SEED GENERATOR TO GET A DIFFERNT FINAL MODEL
env = environ(rand_seed=-49837)
env.io.water = True
env.io.hetatm = True
#soft sphere potential
env.edat.dynamic_sphere=False
#LENNARD_JONES POTENTIAL 
env.edat.dynamic_lennard=True
env.edat.contact_shell = 4.0
env.edat.update_dynamic = 0.39
# Read topology file 
env.libs.topology.read(file='$(LIB)/top_heav.lib')
# Read CHARMM parameter library 
env.libs.parameters.read(file='$(LIB)/par.lib')

###########################################################################################################################

# READ THE OROGINAL PDB FILE AND EXTRACT THE SEQUENCE FOR ALIGNMENT 
mdl1 = model(env, file=args.modelname)
ali = alignment(env)
ali.append_model(mdl1, atom_files=args.modelname, align_codes=args.modelname)

# SELECT RESIDUE ATOMS FOR MUTATION 
s=[]
for n in range(len(args.respos)):
	s.append(selection(mdl1.chains[args.chain[n]].residues[str(args.respos[n])]))
	print(s[n])
	if args.restyp[n] not in residb:
		print("\nUnknown residue type", args.restyp[n], "\nCheck and correct the 3-letter code for the mutant residue:", args.respos[n], args.restyp[n], args.chain[n])
	else:
		s[n].mutate(residue_type=args.restyp[n])

###########################################################################################################################

# GENERATE SEQUENCE, TOPOLOGY AND BUILD THE MUTANT STRUCTURE COORDINATES 
 
#Get sequence. 
ali.append_model(mdl1, align_codes=args.modelname)

#Generate molecular topology for mutant
mdl1.clear_topology()
mdl1.generate_topology(ali[-1])

#Transfer and build all the coordinates you can from the template native structure to the mutant
mdl1.transfer_xyz(ali)
mdl1.build(initialize_xyz=False, build_method='INTERNAL_COORDINATES')

#yes model2 is the same file as model1.  It's a modeller trick.
mdl2 = model(env, file=args.modelname)

#transfers from "model 2" to "model 1"
mdl1.res_num_from(mdl2,ali)

#Usually necessary to write the mutated sequence out and read it in before proceeding, because not all sequence related information about MODEL
#is changed by this command (e.g., internal coordinates, charges, and atom types and radii are not updated).
mdl1.write(file=args.modelname+'-'+args.suffix+'.tmp')
mdl1.read(file=args.modelname+'-'+args.suffix+'.tmp')

#Set up restraints before computing energy
make_restraints(mdl1, ali)

#a non-bonded pair has to have at least as many selected atoms
mdl1.env.edat.nonbonded_sel_atoms=1
sched = autosched.loop.make_for_model(mdl1)

###########################################################################################################################

#Update the selection to mutated residue atoms only 
r=[]
r_all = selection()

for m in range(len(args.respos)):
	r.append(selection(mdl1.chains[args.chain[m]].residues[str(args.respos[m])]))
	r_all = r_all|r[m]

#Restrain on the  mutated atoms
mdl1.restraints.unpick_all()
mdl1.restraints.pick(r_all)

#This will calculate the stereochemical energy (bonds, angles, dihedrals, impropers) for the model with defined restraints.
r_all.energy()

#This command will randomize the cartesion coordination of the selected atoms
r_all.randomize_xyz(deviation=4.0)

#Optimization-1: Turn off dynamic interactions between the selected and unselected regions by setting env.edat.nonbonded_sel_atoms to 2 (by default it is 1) 
mdl1.env.edat.nonbonded_sel_atoms=2
optimize(r_all, sched)

#Optimization-2: Feels environment (energy computed on pairs that have at least one member in the selected)
mdl1.env.edat.nonbonded_sel_atoms=1
optimize(r_all, sched)

#Re-calculate the energy of the model after optimization
r_all.energy()

#Write out the mutated protein PDB structure
mdl1.write(file=args.modelname+'-'+args.suffix+'.pdb')

#Delete the temporary file
os.remove(args.modelname+'-'+args.suffix+'.tmp')

