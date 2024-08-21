#!/bin/bash
#
#  - Extracts just the protein component from a given trajectory
#  - Calculates rmsds (global, individual chains)
#  - Extracts the last frame
#  - Extracts clusters and cluster center structures
#

function fdu()
{
local a=`du -sh $1 2> /dev/null | awk '{print $1}'`
echo $a
}

if [[ $# -ne 1 ]]; then
   echo "------------------------------------------------------------------------"
   echo 
   echo "$0 out-tag"
   echo "      Also needed are **protein** only top, tpr, ndx for fixing  trajectory"
   echo "               top: prot.top  tpr: prot.tpr  ndx: rmsd.ndx "
   echo 
   exit
fi

TPR=topol.tpr
XTC=traj_comp.xtc
TAG=$1

##########################################################################################################
################### SET  GROMACS VERSION AND ADD PREFIX/SUFFIX ACCORDINGLY ###############################
##########################################################################################################

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

shopt -s expand_aliases
if [ $gmxversion = 1 ] ; then alias genbox='${GMX_PREFIX} genbox' ; else alias genbox='${GMX_PREFIX} solvate'; fi
if [ $gmxversion = 1 ] ; then alias rms='g_rms' ; else alias rms='${GMX_PREFIX} rms'; fi
if [ $gmxversion = 1 ] ; then alias gyrate='g_gyrate' ; else alias gyrate='${GMX_PREFIX} gyrate'; fi
if [ $gmxversion = 1 ] ; then alias trjconv='g_trjconv' ; else alias trjconv='${GMX_PREFIX} trjconv'; fi
if [ $gmxversion = 1 ] ; then alias rmsf='g_rmsf' ; else alias rmsf='${GMX_PREFIX} rmsf'; fi
if [ $gmxversion = 1 ] ; then alias cluster='g_cluster' ; else alias cluster='${GMX_PREFIX} cluster'; fi
alias pdb2gmx='${GMX_PREFIX} pdb2gmx'
alias editconf='${GMX_PREFIX} editconf'
alias grompp='${GMX_PREFIX} grompp'
alias genion='${GMX_PREFIX} genion'
alias mdrun='${GMX_PREFIX} mdrun'
alias make_ndx='${GMX_PREFIX} make_ndx'

##########################################################################################################
######################## EXTRACT PROTEIN ONLY .top, .tpr., .ndx FILES FOR ANALYSIS #######################
##########################################################################################################

mkdir MD-analysis
cp topol.top prot.top
sed -i '/SOL /d; /Na  /d; /Cl  /d' prot.top
#echo "q " | make_ndx -f ions.gro -o rmsd.ndx
echo "------------------------------------------------------------------------"
echo
read -p "W> Enter 3 or 4 letter residue name of the CHromophore in the system (Ex: SYG or TYG or AYG):" cronm
if [ -z "$cronm" ]; then
   cronm=$(grep "YG" topol.top | tail -1 | awk '{print $4}')
fi

make_ndx -f ions.gro -o rmsd.ndx << EOF
q
EOF
sed -i 's/Protein/Solute/' rmsd.ndx
sed -i 's/Water/Solvent/' rmsd.ndx

rm -f dummy.mdp
cat >> dummy.mdp << endl
pbc=xyz
endl

echo "M> Extracting a protein only structure (for fitting)"
echo "Solute" | trjconv -f $XTC -s $TPR -n rmsd.ndx -dump 100 -o prot.gro -pbc nojump &> log
if [[ ! -s prot.gro ]]; then
   echo "E> (trjconv:0) Failed .. check 'log' "
   exit
fi

echo "M> Running grompp "
grompp -f dummy.mdp -p prot.top -c prot.gro -o prot.tpr >& log
if [[ ! -s prot.tpr ]]; then
   echo "E> (grompp) Failed .. check 'log' "
   exit
fi

echo "------------------------------------------------------------------------"
echo
echo " The sciprt will now compute the following properties ncecessary to perform the analysis on your MD trajectories using GROMACS tools:"
echo "        - Extracts just the protein component from a given trajectory"
echo "        - Generates a visualization trajectory file fitted to reference structure"
echo "        - Calculates RMSD, RMSF, radius of gyration "
echo "        - RMSD fitted to average structure"
echo "        - Extracts clusters and cluster center structures"
echo 

##########################################################################################################
######################## EXTRACT REFERENCE STRUCTURE FITTED VISUALIZATION TRAJECTORY ##################### 
##########################################################################################################

echo "M> Extracting Protein+chromophore trajectory"
echo "Solute" | trjconv -s $TPR -f $XTC -n rmsd.ndx -o ${TAG}-prot.xtc -pbc nojump -e 100000 >& log
if [[ ! -s ${TAG}-prot.xtc ]]; then
   echo "E> Failed .. check 'log' "
   exit
fi

echo "M> Extracting a structure at 100 ps (for fitting)"
echo "Solute Solute Solute" | trjconv -f ${TAG}-prot.xtc -n rmsd.ndx -dump 100 -o 100ps.gro -pbc cluster -s prot.tpr -center -ur compact &> log
if [[ ! -s 100ps.gro ]]; then
   echo "E> (trjconv:0) Failed .. check 'log' "
   exit
fi

echo "M> Running grompp "
grompp -f dummy.mdp -p prot.top -c 100ps.gro -o 100ps.tpr >& log
if [[ ! -s 100ps.tpr ]]; then
   echo "E> (grompp) Failed .. check 'log' "
   exit
fi

echo "M> Clustering trajectory"
echo "Solute Solute Solute" | trjconv -s 100ps.tpr -pbc cluster -n rmsd.ndx -f ${TAG}-prot.xtc -center -o ${TAG}-prot-clust.xtc -ur compact>& log
if [[ ! -s ${TAG}-prot-clust.xtc ]]; then
   echo "E> (clust) Failed ... check 'log' "
   exit
fi

echo "M> Fitting trajectory (visualization)"
echo "Solute Solute" | trjconv -fit progressive -s 100ps.tpr -n rmsd.ndx -f ${TAG}-prot-clust.xtc -o ${TAG}-prot-clust-fit.xtc >& log
if [[ ! -s ${TAG}-prot-clust-fit.xtc ]]; then
   echo  "E> (fit)  Failed .. check 'log' "
   exit
fi



##########################################################################################################
######################## EXTRACT REFERENCE STRUCTURE FITTED RMSD and  GYRATION ######################### 
##########################################################################################################

echo "M> Running g_rms "
echo "        Protein + chromophore"
echo "Solute Solute" | rms -f ${TAG}-prot-clust.xtc -o rmsd-prot-${TAG}.xvg -fit rot+trans -s prot.tpr -n rmsd.ndx >& log
if [[ ! -s rmsd-prot-${TAG}.xvg ]]; then
   echo "E> (g_rms:0) Failed ... check 'log' "
   exit
else
   mv rmsd-prot-${TAG}.xvg MD-analysis/
fi

echo "M> Running g_rms "
echo "        backbone"
echo "Backbone Backbone" | rms -f ${TAG}-prot-clust.xtc -o rmsd-bb-${TAG}.xvg -fit rot+trans -s prot.tpr -n rmsd.ndx >& log
if [[ ! -s rmsd-bb-${TAG}.xvg ]]; then
   echo "E> (g_rms:0) Failed ... check 'log' "
   exit
else
   mv rmsd-bb-${TAG}.xvg MD-analysis/
fi

echo "M> Running g_gyrate to extract Radius of Gyration "
echo "          Protein + chromophore"
echo "Solute" | gyrate -f ${TAG}-prot-clust-fit.xtc -s 100ps.tpr -o gyrate_${TAG}.xvg -n rmsd.ndx >& log
if [[ ! -s gyrate_${TAG}.xvg  ]]; then
   echo "E> (g_gyrate:0) Failed ... check 'log' "
   exit
else
   mv gyrate_${TAG}.xvg MD-analysis/
fi

##########################################################################################################
########################       EXTRACT AVERGAE STRUCTURE FITTED RMSD      ################################ 
##########################################################################################################

echo "M> Running g_rmsf to extract average structure "
echo "		Protein + chromophore"
echo "Solute" | rmsf -s $TPR -f $XTC -n rmsd.ndx -ox xaver_${TAG}.pdb -o rmsf_${TAG}.xvg -oq bfac_${TAG}.pdb >& log
if [[ ! -s rmsf_${TAG}.xvg || ! -s xaver_${TAG}.pdb  ]]; then
   echo "E> (g_rmsf:0) Failed ... check 'log' "
   exit
else
   mv rmsf_${TAG}.xvg xaver_${TAG}.pdb bfac_${TAG}.pdb MD-analysis/
fi

echo "M> Running g_rms against average structure from g_rmsf "
echo "          Protein + chromophore"
echo "Solute Solute"| rms -f ${TAG}-prot-clust.xtc -s MD-analysis/xaver_${TAG}.pdb -n rmsd.ndx -fit rot+trans -o rmsd-xaver-prot >& log
if [[ ! -s rmsd-xaver-prot.xvg ]]; then
   echo "E> (g_rms:0) Failed ... check 'log' "
   exit
else
   mv rmsd-xaver-prot.xvg MD-analysis/
fi

##########################################################################################################
########################   CLUSTER TRAJECTORIES TO EXTRACT SET OF REPRESENTATIVE FRAMES  #################
##########################################################################################################

echo "M> Running g_cluster to extract cluster size and center  structure "
echo "          Protein + chromophore"
read -p "W> Enter RMSD-cutoff for cluster MD-analysis (recommended values for GFP derivates is 0.30) :" csize 
echo "Solute Solute" | cluster -f $XTC -s $TPR -n rmsd.ndx -cutoff $csize -b 5000 -method gromos -o rmsd-clust-${csize} -g clusters-${csize}-${TAG} -cl clusters-${csize}-${TAG}.pdb -dist rmsd-dist-${csize}-${TAG} >& log
if [[ ! -s clusters-${csize}-${TAG}.log || ! -s clusters-${csize}-${TAG}.pdb  ]]; then
   echo "E> (g_cluster:0) Failed ... check 'log' "
   exit
else
   mv rmsd-clust-${csize}.xpm clusters-${csize}-${TAG}.log clusters-${csize}-${TAG}.pdb rmsd-dist-${csize}.xvg MD-analysis/
fi

##########################################################################################################

rm -f log dummy.mdp
rm -f ${TAG}-prot.xtc ${TAG}-prot-clust.xtc
echo  "M> $XTC: `fdu $XTC`  clust-fit:`fdu ${TAG}-prot-clust-fit.xtc` "
mv ${TAG}-prot-clust-fit.xtc 100ps.* prot.* rmsd.ndx MD-analysis/
echo "M> Done"
