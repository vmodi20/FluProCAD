import os
import shutil
import pmx
import alchemlyb
import math

if shutil.which('gmx') is None and shutil.which('gmx_mpi') is None:
    print("Gromacs is not installed or not in the PATH. Please install Gromacs and add it to the PATH.")
    exit()
else:
    print("Gromacs is installed and in the PATH.")
    print("#"*50)

if shutil.which('gmx') is not None:
    GMXCMD = "gmx"
elif shutil.which('gmx_mpi') is not None:
    GMXCMD = "gmx_mpi"
else:
    print("E> Neither 'gmx' nor 'gmx_mpi' command found. Please load the GROMACS module or install it.")
    exit()

# Estimate free energy usign the Bennet Acceptance Ratio method 
def estimate_fe_bar():
    os.chdir("slow-TI")
    for i in range(21):
        lambdapoint = str(i)
        src = os.path.join(lambdapoint, 'dhdl.xvg')
        if os.path.isfile(src):
            dest = f'dhdl-{lambdapoint}.xvg'
            shutil.copy(src, dest)

    dhdl_files = [f'dhdl-{i}.xvg' for i in range(21)]    
    dhdl_data = [alchemlyb.parsing.gmx.extract_dHdl(file) for file in dhdl_files]
    combined_data = alchemlyb.concat(dhdl_data)

    bar = alchemlyb.estimators.BAR()
    bar.fit(combined_data)

    return bar.delta_f_.loc[0.00,1.00], bar.d_delta_f_.loc[0.00,1.00]


df = {}
for prot_state in "monomer", "dimer", "unfolded": 
    os.chdir(prot_state)
    if prot_state == "unfolded" :
        for dir in os.listdir():
            if os.path.isdir(dir):
                os.chdir(dir)
                df[f'dF_{prot_state}_{dir}'] = [estimate_fe_bar()]
                unf_tot_df = unf_tot_df + df[f'dF_{prot_state}_{dir}'][0]
                err_unf_df = err_unf_df + df[f'dF_{prot_state}_{dir}'][1]**2
                os.chdir('..')
        df[f'dF_{prot_state}'] = [unf_tot_df, math.sqrt(err_unf_df)]
    else:
        df[f'dF_{prot_state}'] = [estimate_fe_bar()]
    os.chdir('..')

# Calculate the ddF values for folding and dimerization : 
dF_folding_monomer = df['dF_monomer'][0] - df['dF_unfolded'][0]
dF_folding_monomer_err = math.sqrt(df['dF_monomer'][1]**2 + df['dF_unfolded'][1]**2)

dF_folding_dimer = df['dF_dimer'][0] - df['dF_unfolded'][0]
dF_folding_dimer_err = math.sqrt(df['dF_dimer'][1]**2 + df['dF_unfolded'][1]**2)

df_dimerization = df['dF_monomer'][0] + df['dF_dimer'][0]
dF_dimerization_err = math.sqrt(df['dF_monomer'][1]**2 + df['dF_dimer'][1]**2)

if dF_folding_monomer > 0 :
    fold_monomer_out = "The ddF values suggest that the mutation(s) have a destabilising effect in comparison to the parent protein system"
else:
    fold_monomer_out = "The ddF values suggest that the mutation(s) have a stabilising effect in comparison to the parent protein system"

if dF_folding_dimer > 0 :
    fold_dimer_out = "The ddF values suggest that the mutation(s) have a destabilising effect in comparison to the parent protein system"
else:
    fold_dimer_out = "The ddF values suggest that the mutation has a stabilising effect in comparison to the parent protein system"

if df_dimerization < 0 :
    dimerize_out = "The ddF values suggest that the mutation(s) have a destabilising effect on the dimer state in comparison to the parent protien system"
else:
    dimerize_out = "The ddF values suggest that the mutation(s) do not destabilize the dimer state"



with open('fe-results.dat','w') as outfile:
    outfile.write(f"Effect of mutation on the folding and dimerization free energy ")
    outfile.write(f"\ndF_folding_monomer: {dF_folding_monomer} ± {dF_folding_monomer_err}\n")
    outfile.write(f"# {fold_monomer_out}\n")
    outfile.write(f"\ndF_folding_dimer: {dF_folding_dimer} ± {dF_folding_dimer_err}\n")
    outfile.write(f"# {fold_dimer_out}\n")
    outfile.write(f"dF_dimerization: {df_dimerization} ± {dF_dimerization_err}\n")
    outfile.write(f"# {dimerize_out}\n")
    outfile.close()




