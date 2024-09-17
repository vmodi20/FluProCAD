import os
import importlib.util
import yaml
import sys

# Check if yaml file is provided, load the arguments from the yaml file
if len(sys.argv) > 1:
    if sys.argv[1].endswith('.yaml'):
        if os.path.exists(sys.argv[1]):  
            with open(sys.argv[1], 'r') as yaml_file:
                arguments = yaml.load(yaml_file, Loader=yaml.FullLoader)
                pdbid = arguments['pdbid']
                pdbfile = arguments['pdbfile']
                ph = arguments['ph']
                mutfile = arguments['mutfile']
                crotype = arguments['crotype']
                cresid = arguments['cresid']
                protoligmr = arguments['protoligmr']
                fe = arguments['fe']
        else:
            print("#"*50)
            print("ERROR: Please check that the yaml file exists in the current working directory!\n")
            exit()
    else:
        print("#"*50)
        print("ERROR: Please provide a valid yaml file!\n")
        exit()
elif len(sys.argv) > 2:
    print("#"*50)
    print("ERROR: Please provide yaml file as single argument or no arguments to generate a new yaml input file!\n")
    exit()
else:
    print("#"*50)
    # Convert the parsed arguments to user prompts
    pdbid = input("\nEnter value for pdbid (leave blank if not applicable): ") 
    pdbfile = input("\nEnter value for pdbfile (leave blank if not applicable): ") 
    ph = input("\nEnter pH value (default is 7): ")
    mutfile = input("\nEnter name of file with mutations: ") 
    crotype = input("\nEnter type of chromophore (SYG or TYG): ") 
    cresid = input("\nEnter residue position for the chromophore: ") 
    protoligmr = input("\nIs your protein complex a monomer or a dimer? (default is dimer): ") 
    fe = input("\nPerform free energy calculations (yes/no, default is yes): ") 

    arguments = {
        'pdbid': pdbid,
        'pdbfile': pdbfile,
        'ph': ph,
        'mutfile': mutfile,
        'crotype': crotype,
        'cresid': cresid,
        'protoligmr': protoligmr,
        'fe': fe
    }

# check if pdbid or pdbfilename is provided or if both are provided
if pdbid:
    protname = pdbid
elif pdbfile:
    protname = pdbfile.split('.')[0]
else: 
    print("#"*50)
    print("Please provide either a PDB file or a PDB ID!\n")
    exit()

if pdbid and pdbfile:
    print("#"*50)
    print("Please provide either a PDB file or a PDB ID!\n")
    exit()

# check if the chromophore type is valid
if crotype not in ["SYG", "TYG"]:
    print("#"*50)
    print("Please provide a valid chromophore type: SYG or TYG\n")
    exit()

# check if the protein complex type is valid
if protoligmr not in ["monomer", "dimer"]:
    print("#"*50)
    print("ERROR: Please provide a valid protein complex type: monomer or dimer\n")
    exit()

if ph:
    ph = round(ph)
elif not ph:
    ph = 7

# Save the arguments to a YAML file
yaml_file_path = os.path.join(os.path.dirname(__file__), protname+'_FPC_input.yaml')
with open(yaml_file_path, 'w') as yaml_file:
    yaml.dump(arguments, yaml_file)

print("#"*50)
print("The input has been saved to the file: ", yaml_file_path)    

# Run the python scripts "PreModellingPrep.py" and "MDsetup.py" with the user inputs
if pdbid:
    os.system('python PreModellingPrep.py -pdbid {} -ph {} -mutfile {} -crotype {} -cresid {} -protoligmr {} -fe {} > prep.log'.format(pdbid, ph, mutfile, crotype, cresid, protoligmr, fe))
elif pdbfile:
    os.system('python PreModellingPrep.py -pdbfile {} -ph {} -mutfile {} -crotype {} -cresid {} -protoligmr {} -fe {} > prep.log'.format(pdbid, pdbfile, ph, mutfile, crotype, cresid, protoligmr, fe))

os.system('python MDsetup.py -proteindir {} -mutfile {} -protoligmr {} -fe {} >> setup.log'.format(protname, mutfile, protoligmr, fe))

