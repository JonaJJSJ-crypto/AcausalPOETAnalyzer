#!/bin/bash

cd ${_CONDOR_SCRATCH_DIR}
cat <<'EndOfFile' > execute.sh
#!/bin/bash

# Exit on error
set -e

echo "### Begin of job"

ID=$1
echo "ID:" $ID

PROCESS=$2
echo "Process:" $PROCESS

FILE=$3
echo "File:" $FILE

EOS_HOME=/eos/user/e/ecarrera
#EOS_HOME=/eos/user/FIRST_LETTER/USERNAME
echo "EOS home:" $EOS_HOME

OUTPUT_DIR=${EOS_HOME}/lee-wick/signal/
#OUTPUT_DIR=${EOS_HOME}/opendata_files/
echo "Output directory:" $OUTPUT_DIR

CMSSW_BASE=/afs/cern.ch/work/e/ecarrera/lw/simulations/CMSSW_5_3_32
#CMSSW_BASE=/afs/cern.ch/work/FIRST_LETTER/USERNAME/CMSSW_5_3_32
echo "CMSSW base:" $CMSSW_BASE

CONFIG=${CMSSW_BASE}/src/AODgeneration/config/merge.py

echo "CMSSW config:" $CONFIG

echo "Hostname:" `hostname`

echo "How am I?" `id`

echo "Where am I?" `pwd`

echo "What is my system?" `uname -a`

echo "### Start working"

# Trigger auto mount of EOS
ls -la $EOS_HOME

# Make output directory
mkdir -p ${OUTPUT_DIR}/${PROCESS}/merged

# Setup CMSSW
THIS_DIR=$PWD
cd $CMSSW_BASE
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc481
eval `scramv1 runtime -sh`
cd $THIS_DIR

# Copy config file
mkdir -p configs/
CONFIG_COPY=configs/cfg_${ID}.py
cp $CONFIG $CONFIG_COPY



# Modify CMSSW config to run only a single file list
sed -i -e "s,data/CMSPrivate_MonteCarlo_LWSM200DnR_unmerged_file_index.txt,${CMSSW_BASE}/src/AODgeneration/${FILE},g" $CONFIG_COPY

# Print config
cat $CONFIG_COPY

# Run CMSSW config
cmsRun $CONFIG_COPY

# Copy output file
xrdcp -f mymerged.root root://eosuser.cern.ch/${OUTPUT_DIR}/${PROCESS}/merged/recoLW_merged_${PROCESS}_${ID}.root
rm mymerged.root

echo "### End of job"

EndOfFile

# Trigger auto mount of EOS
ls -la $EOS_HOME

# Make file executable
chmod +x execute.sh

export SINGULARITY_CACHEDIR="/tmp/$(whoami)/singularity"
singularity run -B /afs -B /eos -B /cvmfs -B /usr/libexec/condor -B /pool --no-home docker://cmssw/slc6:latest $(echo $(pwd)/execute.sh $1 $2 $3)
