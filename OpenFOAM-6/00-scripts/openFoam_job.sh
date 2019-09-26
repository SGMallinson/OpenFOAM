#!/bin/bash
#PBS -N rti_bm
#PBS -l nodes=1:ppn=24
#PBS -l mem=6GB
#PBS -l walltime=1:00:00
#PBS -j oe
#PBS -M s.mallinson@unsw.edu.au
#PBS -m ae

cd /home/z9104023/projects/huang

# Unload modules.
#module rm openmpi intel-cc intel-fc
 
# Load modules.
module load openfoam/1812
module load openmpi/4.0.1

export LD_LIBRARY_PATH=$FOAM_USER_LIBBIN:$LD_LIBRARY_PATH

. ${WM_PROJECT_DIR}/bin/tools/RunFunctions

runApplication blockMesh
cp -r 0/alpha.liquid.org 0/alpha.liquid

runApplication decomposePar
runParallel interFoam
runApplication reconstructPar

# post-process
#module load python/3.7.3
#runApplication -s parse python $PYFLAGS parse.py $APPLICATION
#runApplication -s plot python frontier.py ${APPLICATION}.csv
