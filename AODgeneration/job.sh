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

FILE=${PROCESS}.lhe
echo "File:" $FILE

MBIDXFILE=CMS_MonteCarlo2012_TuneZ2star_8TeV-pythia6_GEN-SIM_START50_V13-v3_file_index.txt
echo "Min Bias index file: " $MBIDXFILE

EVTS=$3
echo "Events: " $EVTS

SKIPEVTS=$4
echo "skipEvents: " $SKIPEVTS

FIRSTEVT=$5
echo "firstEvent: " $FIRSTEVT

FIRSTRUN=$6
echo "firstRun: " $FIRSTRUN

THEMASS=$7
echo "themass: " $THEMASS

THETAU=$8
echo "thetau: " $THETAU 


EOS_HOME=/eos/user/e/ecarrera
#EOS_HOME=/eos/user/FIRST_LETTER/USERNAME
echo "EOS home:" $EOS_HOME

OUTPUT_DIR=${EOS_HOME}/lee-wick/signal/
#OUTPUT_DIR=${EOS_HOME}/opendata_files/
echo "Output directory:" $OUTPUT_DIR

CMSSW_BASE=/afs/cern.ch/work/e/ecarrera/lw/simulations/CMSSW_5_3_32
#CMSSW_BASE=/afs/cern.ch/work/FIRST_LETTER/USERNAME/CMSSW_5_3_32
echo "CMSSW base:" $CMSSW_BASE

GENCONFIG=${CMSSW_BASE}/src/AODgeneration/config/gensimLW.py
HLTCONFIG=${CMSSW_BASE}/src/AODgeneration/config/hltLW.py
RECOCONFIG=${CMSSW_BASE}/src/AODgeneration/config/recoLW.py

LHEFILE=${CMSSW_BASE}/src/AODgeneration/data/${FILE}
MBFILES=${CMSSW_BASE}/src/AODgeneration/data/${MBIDXFILE}

echo "gen config:" $GENCONFIG
echo "hlt config:" $HLTCONFIG
echo "reco config:" $RECOCONFIG

echo "minbias index file:" $MBFILES

echo "Hostname:" `hostname`

echo "How am I?" `id`

echo "Where am I?" `pwd`

echo "What is my system?" `uname -a`

echo "### Start working"

# Trigger auto mount of EOS
ls -la $EOS_HOME

# Make output directory
mkdir -p ${OUTPUT_DIR}/${PROCESS}

# Setup CMSSW
THIS_DIR=$PWD
cd $CMSSW_BASE
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc481
eval `scramv1 runtime -sh`
cd $THIS_DIR

# Copy config files
mkdir -p configs/
GENCONFIG_COPY=configs/gencfg_${ID}.py
cp $GENCONFIG $GENCONFIG_COPY
HLTCONFIG_COPY=configs/hltcfg_${ID}.py
cp $HLTCONFIG $HLTCONFIG_COPY
RECOCONFIG_COPY=configs/recocfg_${ID}.py
cp $RECOCONFIG $RECOCONFIG_COPY



# Modify the gensim config 
sed -i -e "s,LWSM200DnR.lhe,${LHEFILE},g" $GENCONFIG_COPY
sed -i -e "s,input = cms.untracked.int32(20),input = cms.untracked.int32(${EVTS}),g" $GENCONFIG_COPY
sed -i -e "s,skipEvents = cms.untracked.uint32(0),skipEvents = cms.untracked.uint32(${SKIPEVTS}),g" $GENCONFIG_COPY
sed -i -e "s,firstEvent = cms.untracked.uint32(1),firstEvent = cms.untracked.uint32(${FIRSTEVT}),g" $GENCONFIG_COPY
sed -i -e "s,firstRun = cms.untracked.uint32(200000),firstRun = cms.untracked.uint32(${FIRSTRUN}),g" $GENCONFIG_COPY
sed -i -e "s,556:new = lwe- lwe+ 2 -3 0 200.0 0.0 200.0 200.0 2.70765e-02,556:new = lwe- lwe+ 2 -3 0 ${THEMASS} 0.0 ${THEMASS} ${THEMASS} ${THETAU},g" $GENCONFIG_COPY


# Print gen config
cat $GENCONFIG_COPY

# Run CMSSW gen config
cmsRun $GENCONFIG_COPY 

# Copy output file
xrdcp -f gensimLW.root root://eosuser.cern.ch/${OUTPUT_DIR}/${PROCESS}/gensimLW_${PROCESS}_${ID}.root

# Print hlt config
cat $HLTCONFIG_COPY

# Modify the hlt config to deal with pileup sim
THEMBFILE=`shuf -n 1 ${MBFILES}`
echo "Min bias file: " ${THEMBFILE} 
sed -i -e "s,root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2012/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0002/FEF2F4CC-0E6A-E111-96F6-0030487F1C57.root,${THEMBFILE},g" $HLTCONFIG_COPY

# Run CMSSW hlt config
cmsRun $HLTCONFIG_COPY

#Copy output file
xrdcp -f hltLW.root root://eosuser.cern.ch/${OUTPUT_DIR}/${PROCESS}/hltLW_${PROCESS}_${ID}.root

# Print reco config
cat $RECOCONFIG_COPY

# Run CMSSW reco config
cmsRun $RECOCONFIG_COPY

#Copy output file
xrdcp -f recoLW.root root://eosuser.cern.ch/${OUTPUT_DIR}/${PROCESS}/recoLW_${PROCESS}_${ID}.root

#Copy output file
xrdcp -f histProbFunction.root root://eosuser.cern.ch/${OUTPUT_DIR}/${PROCESS}/histProbFunction_${PROCESS}_${ID}.root

#Erase leftovers
rm gensimLW.root
rm hltLW.root
rm recoLW.root
rm recoLW_inDQM.root
rm histProbFunction.root

echo "### End of job"

EndOfFile

# Trigger auto mount of EOS
ls -la $EOS_HOME

# Make file executable
chmod +x execute.sh

export SINGULARITY_CACHEDIR="/tmp/$(whoami)/singularity"
singularity run -B /afs -B /eos -B /cvmfs -B /usr/libexec/condor -B /pool --no-home docker://cmssw/slc6:latest $(echo $(pwd)/execute.sh $1 $2 $3 $4 $5 $6 $7 $8)
