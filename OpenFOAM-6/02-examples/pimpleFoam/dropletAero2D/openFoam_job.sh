#!/bin/bash
#PBS -N 2D2D
#PBS -l nodes=1:ppn=24
#PBS -l mem=6GB
#PBS -l walltime=1:00:00
#PBS -M s.mallinson@unsw.edu.au
#PBS -j oe
#PBS -m ae

cd /home/z9104023/projects/dropletAero/2d

# Load modules.
module load openfoam/1812
module load openmpi/4.0.1

export LD_LIBRARY_PATH=$FOAM_USER_LIBBIN:$LD_LIBRARY_PATH

. ${WM_PROJECT_DIR}/bin/tools/RunFunctions

runApplication blockMesh
runApplication topoSet

cp -r 0.org 0
runApplication decomposePar
runParallel $(getApplication)
runApplication reconstructPar
runApplication postProcess -func vorticity
runApplication foamToVTK

# post-process
module load python/3.7.3
cd plots
make clean
make
