import os
import shutil
import math
import alchemlyb
from alchemlyb.parsing.gmx import extract_dHdl
from alchemlyb.estimators import TI
from alchemlyb.visualisation import plot_ti_dhdl

# Estimate free energy usign the Bennet Acceptance Ratio method 
def estimate_equilb_fe():
    os.chdir("slow-TI")
    for i in range(21):
        lambdapoint = str(i)
        src = os.path.join(lambdapoint, 'dhdl.xvg')
        dest = os.path.join('dhdl-'+str(i)+'.xvg')
        if not os.path.isfile(dest):
            if os.path.isfile(src):
                shutil.copy(src, dest)
        
    dhdl_files = [f'dhdl-{i}.xvg' for i in range(21)]    
    dhdl_data = [extract_dHdl(file, T=300) for file in dhdl_files]
    combined_data = alchemlyb.concat(dhdl_data)

    ti_fe = TI()
    equilb_fe = ti_fe.fit(combined_data)

    # Generate the plot of the free energy profile
    plt = plot_ti_dhdl(equilb_fe)
    plt.figure.savefig('ti_dhdl.png')

    os.chdir('..')
    return equilb_fe.delta_f_.loc[0.00,1.00], equilb_fe.d_delta_f_.loc[0.00,1.00]


df = {}
for prot_state in "monomer", "dimer", "unfolded": 
    os.chdir(prot_state)
    if prot_state == "unfolded" :
        df_unf = 0
        df_unf_err = 0
        for dir in os.listdir():
            if os.path.isdir(dir):
                os.chdir(dir)
                df_tmp, df_tmp_err = estimate_equilb_fe()
                os.chdir('..')
            df_unf += df_tmp
            df_unf_err += df_tmp_err**2
        df[f'dF_{prot_state}'] = [df_unf, math.sqrt(df_unf_err)]
    else:
        df_tmp, df_tmp_err = estimate_equilb_fe()
        df[f'dF_{prot_state}'] = [df_tmp, df_tmp_err]
    os.chdir('..')

# Calculate the ddF values for folding and dimerization : 
df_mono = df['dF_monomer'][0]
df_unf = df['dF_unfolded'][0]
df_dimer = df['dF_dimer'][0]

df_mono_err = df['dF_monomer'][1]
df_unf_err = df['dF_unfolded'][1]
df_dimer_err = df['dF_dimer'][1]


ddf_fold_mono = df_mono - df_unf
ddf_fold_mono_err = math.sqrt(df_mono_err**2 + df_unf_err**2)

ddf_fold_dimer = df_dimer - df_unf*2 
ddf_fold_dimer_err = math.sqrt(df_dimer_err**2 + df_unf_err**2)

ddf_dimerization = df_mono*2 - df_dimer
ddf_dimerization_err = math.sqrt(df_mono_err**2 + df_dimer_err**2)

if ddf_fold_mono > 0 :
    fold_monomer_out = "The ddF values suggest that the mutation(s) have a destabilising effect in comparison to the parent protein system"
else:
    fold_monomer_out = "The ddF values suggest that the mutation(s) have a stabilising effect in comparison to the parent protein system"

if ddf_fold_dimer > 0 :
    fold_dimer_out = "The ddF values suggest that the mutation(s) have a destabilising effect in comparison to the parent protein system"
else:
    fold_dimer_out = "The ddF values suggest that the mutation has a stabilising effect in comparison to the parent protein system"

if ddf_dimerization < 0 :
    dimerize_out = "The ddF values suggest that the mutation(s) have a destabilising effect on the dimer state in comparison to the parent protien system"
else:
    dimerize_out = "The ddF values suggest that the mutation(s) do not destabilize the dimer state"


with open('fe-results.dat','w') as outfile:
    outfile.write(f"# Effect of mutation on the folding and dimerization free energy\n")
    outfile.write(f"\n (i) dF_folding_monomer: {ddf_fold_mono} ± {ddf_fold_mono_err}\n")
    outfile.write(f"> {fold_monomer_out}\n")
    outfile.write(f"\n (ii) dF_folding_dimer: {ddf_fold_dimer} ± {ddf_fold_dimer_err}\n")
    outfile.write(f"> {fold_dimer_out}\n")
    outfile.write(f"\n (iii) dF_dimerization: {ddf_dimerization} ± {ddf_dimerization_err}\n")
    outfile.write(f"> {dimerize_out}\n")
    outfile.close()
