#!/bin/bash

# Define path for job directories
BASE_PATH=/afs/cern.ch/work/e/ecarrera/lw/simulations/CMSSW_5_3_32/src/AODgeneration/joblaunch
#BASE_PATH=/path/to/job/directory
mkdir -p $BASE_PATH

# Set processes
PROCESSES=( \
    #LWSM200DnR\
    #LWSM300DnR\
    LWSM400DnR\
    LWSM500DnR\
    )

# Create JDL files and job directories
for PROCESS in ${PROCESSES[@]}
do
    python create_job.py $PROCESS $BASE_PATH
done

# Submit jobs
THIS_PWD=$PWD
for PROCESS in ${PROCESSES[@]}
do
    cd $BASE_PATH/$PROCESS
    condor_submit job.jdl
    cd $THIS_PWD
done
