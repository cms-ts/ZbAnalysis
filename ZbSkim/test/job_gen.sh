#!/bin/sh

ulimit -c 0

unset PYTHIA8175DATA
unset G4NEUTRONXS

cd $CMS_PATH/slc6_amd64_gcc472/cms/cmssw/CMSSW_6_2_7
eval `scramv1 runtime -sh`
cd -

(time hadd rootTuple_$LSB_BATCH_JID.root $@ ) > job.log 2>&1

exit
