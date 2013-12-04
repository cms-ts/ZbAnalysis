#!/bin/sh

d=0
if [ ! -z "$1" ]; then
  d=$1
fi

export ROOT_HIST=0

cd $CMS_PATH/slc5_amd64_gcc472/cms/cmssw/CMSSW_6_1_1
eval `scramv1 runtime -sh`
cd -

i=1
while [ $i -le 2 ]; do

  root -l -q -b DataMCComp3.C+\($d,\"w_first_jet_pt\",1,$i\)
  root -l -q -b DataMCComp3.C+\($d,\"w_first_jet_eta\",1,$i\)
  root -l -q -b DataMCComp3.C+\($d,\"w_first_jet_pt_b\",1,$i\)
  root -l -q -b DataMCComp3.C+\($d,\"w_first_jet_eta_b\",1,$i\)
  root -l -q -b DataMCComp3.C+\($d,\"w_first_bjet_pt\",1,$i\)
  root -l -q -b DataMCComp3.C+\($d,\"w_first_bjet_eta\",1,$i\)
  root -l -q -b DataMCComp3.C+\($d,\"w_pt_Z_ee\",1,$i\)
  root -l -q -b DataMCComp3.C+\($d,\"w_pt_Z_mm\",1,$i\)
  root -l -q -b DataMCComp3.C+\($d,\"w_pt_Z_ee_b\",1,$i\)
  root -l -q -b DataMCComp3.C+\($d,\"w_pt_Z_mm_b\",1,$i\)
  root -l -q -b DataMCComp3.C+\($d,\"w_Ht\",1,$i\)
  root -l -q -b DataMCComp3.C+\($d,\"w_Ht_b\",1,$i\)
  root -l -q -b DataMCComp3.C+\($d,\"w_delta_phi_mm\",1,$i\)
  root -l -q -b DataMCComp3.C+\($d,\"w_delta_phi_ee\",1,$i\)
  root -l -q -b DataMCComp3.C+\($d,\"w_delta_phi_mm_b\",1,$i\)
  root -l -q -b DataMCComp3.C+\($d,\"w_delta_phi_ee_b\",1,$i\)
  root -l -q -b DataMCComp3.C+\($d,\"w_mass_Zj_ee\",1,$i\)
  root -l -q -b DataMCComp3.C+\($d,\"w_mass_Zj_mm\",1,$i\)
  root -l -q -b DataMCComp3.C+\($d,\"w_mass_Zj_ee_b\",1,$i\)
  root -l -q -b DataMCComp3.C+\($d,\"w_mass_Zj_mm_b\",1,$i\)

  root -l -q -b DataMCComp3.C+\($d,\"w_single_bjet_pt\",1,$i\)
  root -l -q -b DataMCComp3.C+\($d,\"w_single_bjet_eta\",1,$i\)
  root -l -q -b DataMCComp3.C+\($d,\"w_single_pt_Z_ee_b\",1,$i\)
  root -l -q -b DataMCComp3.C+\($d,\"w_single_pt_Z_mm_b\",1,$i\)
  root -l -q -b DataMCComp3.C+\($d,\"w_single_delta_phi_ee_b\",1,$i\)
  root -l -q -b DataMCComp3.C+\($d,\"w_single_delta_phi_mm_b\",1,$i\)
  root -l -q -b DataMCComp3.C+\($d,\"w_single_Ht_b\",1,$i\)

  i=$((i+1))
done

exit
