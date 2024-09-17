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

#write a function to check if the command 'gmx' or 'gmx_mpi' exists
function check_gmx_command() {
   if command -v gmx &> /dev/null; then
      GMX_PREFIX="gmx"
   elif command -v gmx_mpi &> /dev/null; then
      GMX_PREFIX="gmx_mpi"
   else
      echo "E> Neither 'gmx' nor 'gmx_mpi' command found. Please install GROMACS."
      exit 1
   fi
}

check_gmx_command

TPR=topol.tpr
XTC=traj_comp.xtc
TAG=$1
cronm=$2


# Check if the file "newtopol.top" exists
if [[ ! -f "newtopol.top" ]]; then
   echo "E> The file 'newtopol.top' does not exist. Please provide the necessary file."
   TOP=topol.top
else
   TOP=newtopol.top
fi

if [ -z "$cronm" ]; then
   cronm=$(grep "YG" $TOP | tail -1 | awk '{print $4}')
fi

if [ -z "$TAG" ]; then
   TAG="FP"
shopt -s expand_aliases
alias genbox='${GMX_PREFIX} solvate'
alias rms='${GMX_PREFIX} rms'
alias gyrate='${GMX_PREFIX} gyrate'
alias trjconv='${GMX_PREFIX} trjconv'
alias rmsf='${GMX_PREFIX} rmsf'
alias cluster='${GMX_PREFIX} cluster'
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
cp $TOP prot.top
sed -i '/SOL /d; /NA  /d; /CL  /d; /NAN /d; /NAP /d; ' prot.top

make_ndx -f ions.pdb -o rmsd.ndx << EOF
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
echo " Running strutture analysis on the MD trajectories using GROMACS tools:"
echo "        - Extracts just the protein component from a given trajectory"
echo "        - Generates a visualization trajectory file fitted to reference structure"
echo "        - Calculates RMSD, RMSF, radius of gyration "
echo "        - RMSD fitted to average structure"
echo "        - Extracts clusters and cluster center structures"
echo " NOTE: This may take a few hours if running on a local computer "
echo
echo "------------------------------------------------------------------------"

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
echo "Solute Solute" | cluster -f $XTC -s $TPR -n rmsd.ndx -cutoff 0.13 -b 5000 -method gromos -o rmsd-clust -g clusters-${TAG} -cl clusters-${TAG}.pdb -dist rmsd-dist-${TAG} >& log
if [[ ! -s clusters-${TAG}.log || ! -s clusters-${TAG}.pdb  ]]; then
   echo "E> (g_cluster:0) Failed ... check 'log' "
   exit
else
   mv rmsd-clust.xpm clusters-${TAG}.log clusters-${TAG}.pdb rmsd-dist.xvg MD-analysis/
fi

##########################################################################################################

rm -f log dummy.mdp
rm -f ${TAG}-prot.xtc ${TAG}-prot-clust.xtc
echo  "M> $XTC: `fdu $XTC`  clust-fit:`fdu ${TAG}-prot-clust-fit.xtc` "
mv ${TAG}-prot-clust-fit.xtc 100ps.* prot.* rmsd.ndx MD-analysis/
echo "M> Done"
