#!/bin/bash

########################## READ AND CHECK FOR MISSING ARGUMENTS ######################
if [[ $# -ne 2 ]]; then
   echo "------------------------------------------------------------------------"
   echo
   echo " Missing arguments for either the PDB filename or the output suffix" 
   echo "$0  <pdb_filename> <output_suffix>"
   echo 
   exit 0
fi

modelname=$1
suffix=$2

##### Read and prepare structure for the mutations defined in the mutations.dat file ###
########################################################################################

if [[ ! -s mutations.dat ]] ; then
	printf '\n %s \n' "Missing file 'mutations.dat' with list of residue names and positions for mutation."
	exit 1
else
	resnr=($(awk '{print $1}' mutations.dat))
	resnm=($(awk '{print $2}' mutations.dat))
	chain=($(awk '{print $3}' mutations.dat))
fi

########################################################################################
pythonversion1="python2"
pythonversion2="python3"
DEFpythonversion="python"

echo "------------------------------------------------------------------------"
echo
echo "The currently supported architectures are as follows:"
echo
echo '   1) Python 2.3 to  2.7.'
echo '   2) Python 3.0 to 3.7.'

printf "Select the version of python installed on your system:"
read ans
if [ x$ans != x ] ; then
  pythonversion=$ans
else
  pythonversion=$DEFpythonversion
fi

if [ $pythonversion = 1 ] ; then PYTHON_VERSION=$pythonversion1 ;
elif [ $pythonversion = 2 ] ; then PYTHON_VERSION=$pythonversion2 ; fi

########################################################################################

printf '\n %s \n' "Generating your mutant............ "
$PYTHON_VERSION build_mutation.py -modelname ${modelname} -suffix ${suffix} -respos ${resnr[@]} -restyp ${resnm[@]} -chain ${chain[@]} > ${modelname}-${suffix}.log

if [[ ! -s ${modelname}-${suffix}.pdb ]]; then
        echo "E> Oops!! Something went wrong there. Re-check your input model for residue numbering/ missing residues and the log files"
        tail -n 2 ${modelname}-${suffix}.log
        exit 2
fi

################### CHECK GROMACS VERSION AND ADD PREFIX/SUFFIX ACCORDINGLY ############################
gmxversion1="4.x"
gmxversion2="5.x"
gmxversion3="2016.x"
gmxversion4="2018.x"
gmxversion5="2020.x"

echo "------------------------------------------------------------------------"
echo
echo " GROMACS versions compatible with this package:" 
echo
echo '   1) 4.x '
echo '   2) 5.x'
echo '   3) 2016.x' 
echo '   4) 2018.x'
echo '   5) 2020.x'

printf "Select the version of GROMACS loaded/installed  on your system:"
read ans
if [ x$ans != x ] ; then
  gmxversion=$ans
else
  gmxversion=1
fi

prefix1=""
prefix2="gmx"
prefix3="gmx_mpi"

echo "------------------------------------------------------------------------"
echo
echo " GROMACS compilations may differ based on the available hardware" 
echo " The prefix can be used to select from the available built support of the GROMACS available on your system"
echo
echo '   1) None (Local installation with MPI or GPU) '
echo '   2) gmx'
echo '   3) gmx_mpi'

printf "Select the prefix you intend you implement for the loaded GROMACS version:"
read ans
if [ x$ans != x ] ; then
  prefix=$ans
else
  prefix=""
fi

if [ $prefix = 1 ] ; then GMX_PREFIX=$prefix1;
elif [ $prefix = 2 ] ; then GMX_PREFIX=$prefix2; 
elif [ $prefix = 3 ] ; then GMX_PREFIX=$prefix3; fi

echo "$GMX_PREFIX"

shopt -s expand_aliases
if [ $gmxversion = 1 ] ; then alias genbox='${GMX_PREFIX} genbox' ; else alias genbox='${GMX_PREFIX} solvate'; fi
#alias genbox='${GMX_PREFIX} solvate'
alias pdb2gmx='${GMX_PREFIX} pdb2gmx'
alias editconf='${GMX_PREFIX} editconf'
alias grompp='${GMX_PREFIX} grompp'
alias genion='${GMX_PREFIX} genion'
alias mdrun='${GMX_PREFIX} mdrun'
#######################################################################################################


if [[ ! -s ${modelname}-${suffix}.pdb ]]; then 
	echo "E> Oops!! Something went wrong there. Re-check your input model for residue numbering/ missing residues and the log files"
	tail -n 2 ${modelname}-${suffix}.log
	exit 3
else
	sed -i 's/CLE /CLEU/;s/NLY /NLYP/' ${modelname}-${suffix}.pdb
	rm -rf ${modelname}-${suffix}
	mkdir ${modelname}-${suffix}
        cp -r amber03_gfp.ff residuetypes.dat specbond.dat xlateat.dat mdp_files/*.mdp ${modelname}-${suffix}.pdb ${modelname}-${suffix}/
	mv ${modelname}-${suffix}.pdb ${modelname}-${suffix}.log ${modelname}-${suffix}/
	cd ${modelname}-${suffix}

	for ((len=0; len<${#resnr[@]}; len++)) ; do
		if [[ ${resnm[len]} == "HIS" ]]; then
			read -p "Enter the protonation form of Histidine residue (HIE/HID/HIP):" his
			sed -i "s/${resnm[len]} ${chain[len]} ${resnr[len]}/$his ${chain[len]} ${resnr[len]}/" ${modelname}-${suffix}.pdb
		fi
	done

	echo " Generating GROMACS topologies for your mutant..."
	echo "1" | pdb2gmx -f ${modelname}-${suffix}.pdb -water tip3p -ignh >& pdb2gmx.log
	if [[ ! -s conf.gro ]]; then
		echo "E> (pdb2gmx) Failed .. check 'pdb2gmx.log' "
		exit 4
	else
		grep -A 1 "PLEASE NOTE" pdb2gmx.log | tail -n 1
		grep "Total charge in system" pdb2gmx.log
	fi

	editconf -f conf.gro -o ed.gro -d 1.0 -bt cubic >& editconf.log

	echo " Adding water to the box..."
	genbox -cp ed.gro -cs spc216.gro -p -o box.gro >& genbox.log
	if [[ ! -s box.gro ]]; then
               	echo "E> (genbox) Failed .. check 'genbox.log' "
		exit 5
	else
		echo $(grep 'SOL molecules' genbox.log | awk '{print "Added",$3,$4,$5}')
	fi

	echo " Adding ions to the box..."
	grompp -f ions.mdp -p topol.top -c box.gro -o genion >& grompp.log
	echo "SOL" | genion -s genion.tpr -p topol.top -o ions.gro -pname Na -nname Cl -neutral -conc 0.15 >& genion.log
	if [[ ! -s genion.tpr || ! -s ions.gro ]]; then
		echo "E> Adding ions failed .. check 'grompp.log and/or genion.log' "
		exit 6
	else
		tail -n 1 genion.log
	fi
	rm -f \#*
	sed -i "s/\.\/amber/\.\.\/amber/"  *.top
	rm -rf amber03_gfp.ff specbond.dat  xlateat.dat
	cd ..
	cp gromacs_jobs.sh analyze-jobs.sh analyze-traj.sh ${modelname}-${suffix}/
fi



