#!/bin/sh

ulimit -c 0

cd $CMS_PATH/slc5_amd64_gcc472/cms/cmssw/CMSSW_6_1_1
eval `scramv1 runtime -sh`
cd -

(time hadd $opts rootTuple_$LSB_BATCH_JID.root $@ ) > job.log 2>&1

exit
