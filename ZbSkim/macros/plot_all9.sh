#!/bin/sh

export ROOT_HIST=0

unset PYTHIA8175DATA
unset G4NEUTRONXS

n=0
if [ ! -z "$1" ]; then
  n=$1
fi

cd $CMS_PATH/slc6_amd64_gcc472/cms/cmssw/CMSSW_6_2_7
eval `scramv1 runtime -sh`
cd -

root -l -q -b DataMCComp9.C+\(\"w_first_bjet_pt\",1,0,$n\)
root -l -q -b DataMCComp9.C+\(\"w_first_bjet_eta\",1,0,$n\)
if [ $n -eq 0 ]; then
  root -l -q -b DataMCComp9.C+\(\"w_first_bjet_eta_abs\",1,0,$n\)
fi

if [ $n -eq 0 ]; then
  root -l -q -b DataMCComp9.C+\(\"w_pt_Z\",1,0,$n\)
  root -l -q -b DataMCComp9.C+\(\"w_Ht\",1,0,$n\)
  root -l -q -b DataMCComp9.C+\(\"w_delta_phi\",1,0,$n\)
  root -l -q -b DataMCComp9.C+\(\"w_mass_Zj\",1,0,$n\)
fi

if [ $n -ne 0 ]; then
  root -l -q -b DataMCComp9.C+\(\"w_pt_Z_b\",1,0,$n\)
  root -l -q -b DataMCComp9.C+\(\"w_Ht_b\",1,0,$n\)
fi

if [ $n -eq 2 ]; then
  root -l -q -b DataMCComp9.C+\(\"w_DR_bb\",1,0,$n\)
  root -l -q -b DataMCComp9.C+\(\"w_DR_Zb_min\",1,0,$n\)
  root -l -q -b DataMCComp9.C+\(\"w_DR_Zb_max\",1,0,$n\)
  root -l -q -b DataMCComp9.C+\(\"w_A_Zb\",1,0,$n\)
  root -l -q -b DataMCComp9.C+\(\"w_Zbb_mass\",1,0,$n\)
  root -l -q -b DataMCComp9.C+\(\"w_bb_mass\",1,0,$n\)
  root -l -q -b DataMCComp9.C+\(\"w_delta_phi_2b\",1,0,$n\)
fi

if [ $n -eq 0 ]; then
  root -l -q -b DataMCComp9.C+\(\"w_first_bjet_pt\",1,1,$n\)
  root -l -q -b DataMCComp9.C+\(\"w_first_bjet_eta\",1,1,$n\)
  root -l -q -b DataMCComp9.C+\(\"w_first_bjet_eta_abs\",1,1,$n\)
fi

if [ $n -eq 0 ]; then
  root -l -q -b DataMCComp9.C+\(\"w_pt_Z\",1,1,$n\)
  root -l -q -b DataMCComp9.C+\(\"w_Ht\",1,1,$n\)
  root -l -q -b DataMCComp9.C+\(\"w_delta_phi\",1,1,$n\)
fi

exit
