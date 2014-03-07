#!/bin/sh

export ROOT_HIST=0

unset PYTHIA8175DATA
unset G4NEUTRONXS

cd $CMS_PATH/slc6_amd64_gcc472/cms/cmssw/CMSSW_6_2_7
eval `scramv1 runtime -sh`
cd -

root -l -q -b DataMCComp8.C+\(\"w_first_bjet_pt\",1,0\)
root -l -q -b DataMCComp8.C+\(\"w_first_bjet_eta\",1,0\)
root -l -q -b DataMCComp8.C+\(\"w_pt_Z\",1,0\)
root -l -q -b DataMCComp8.C+\(\"w_y_Z\",1,0\)
root -l -q -b DataMCComp8.C+\(\"w_delta_phi\",1,0\)
root -l -q -b DataMCComp8.C+\(\"w_Ht\",1,0\)
root -l -q -b DataMCComp8.C+\(\"w_mass_Zj\",1,0\)

root -l -q -b DataMCComp8.C+\(\"w_first_bjet_pt\",1,1\)
root -l -q -b DataMCComp8.C+\(\"w_first_bjet_eta\",1,1\)
root -l -q -b DataMCComp8.C+\(\"w_pt_Z\",1,1\)
root -l -q -b DataMCComp8.C+\(\"w_y_Z\",1,1\)
root -l -q -b DataMCComp8.C+\(\"w_delta_phi\",1,1\)
root -l -q -b DataMCComp8.C+\(\"w_Ht\",1,1\)
root -l -q -b DataMCComp8.C+\(\"w_mass_Zj\",1,1\)

exit
