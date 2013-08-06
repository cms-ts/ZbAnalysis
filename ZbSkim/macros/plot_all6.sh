#!/bin/sh

export ROOT_HIST=0

cd $CMS_PATH/slc5_amd64_gcc472/cms/cmssw/CMSSW_6_1_0
eval `scramv1 runtime -sh`
cd -

root -l -q -b DataMCComp6.C\(\"w_first_bjet_pt\",1,0\)
root -l -q -b DataMCComp6.C\(\"w_first_bjet_eta\",1,0\)
root -l -q -b DataMCComp6.C\(\"w_pt_Z\",1,0\)
root -l -q -b DataMCComp6.C\(\"w_Ht\",1,0\)
root -l -q -b DataMCComp6.C\(\"w_delta_phi\",1,0\)

root -l -q -b DataMCComp6.C\(\"w_first_bjet_pt\",1,1\)
root -l -q -b DataMCComp6.C\(\"w_first_bjet_eta\",1,1\)
root -l -q -b DataMCComp6.C\(\"w_pt_Z\",1,1\)
root -l -q -b DataMCComp6.C\(\"w_Ht\",1,1\)
root -l -q -b DataMCComp6.C\(\"w_delta_phi\",1,1\)

exit
