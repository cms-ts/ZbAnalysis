#!/bin/sh

d=0
if [ ! -z "$1" ]; then
  d=$1
fi

n=0
if [ ! -z "$2" ]; then
  n=$2
fi

export ROOT_HIST=0

unset PYTHIA8175DATA
unset G4NEUTRONXS

cd $CMS_PATH/slc6_amd64_gcc472/cms/cmssw/CMSSW_6_2_7
eval `scramv1 runtime -sh`
cd -

i=1
while [ $i -le 2 ]; do

  if [ $n -eq 0 ]; then
    root -l -q -b DataMCComp3.C+\($d,\"w_first_jet_pt\",1,$i\)
    root -l -q -b DataMCComp3.C+\($d,\"w_first_jet_eta\",1,$i\)
    root -l -q -b DataMCComp3.C+\($d,\"w_first_jet_eta_abs\",1,$i\)

    root -l -q -b DataMCComp3.C+\($d,\"w_first_jet_pt_b\",1,$i\)
    root -l -q -b DataMCComp3.C+\($d,\"w_first_jet_eta_b\",1,$i\)
  fi

  root -l -q -b DataMCComp3.C+\($d,\"w_first_bjet_pt\",1,$i,$n\)
  root -l -q -b DataMCComp3.C+\($d,\"w_first_bjet_eta\",1,$i,$n\)
  if [ $n -eq 0 ]; then
    root -l -q -b DataMCComp3.C+\($d,\"w_first_bjet_eta_abs\",1,$i,$n\)
  fi

  if [ $n -eq 0 ]; then
    root -l -q -b DataMCComp3.C+\($d,\"w_pt_Z_ee\",1,$i\)
    root -l -q -b DataMCComp3.C+\($d,\"w_pt_Z_mm\",1,$i\)
    root -l -q -b DataMCComp3.C+\($d,\"w_Ht\",1,$i\)
    root -l -q -b DataMCComp3.C+\($d,\"w_delta_phi_ee\",1,$i\)
    root -l -q -b DataMCComp3.C+\($d,\"w_delta_phi_mm\",1,$i\)
    root -l -q -b DataMCComp3.C+\($d,\"w_mass_Zj_ee\",1,$i\)
    root -l -q -b DataMCComp3.C+\($d,\"w_mass_Zj_mm\",1,$i\)

    root -l -q -b DataMCComp3.C+\($d,\"w_y_Z_ee\",1,$i\)
    root -l -q -b DataMCComp3.C+\($d,\"w_y_Z_mm\",1,$i\)
    root -l -q -b DataMCComp3.C+\($d,\"w_Phi_star_ee\",1,$i,$n\)
    root -l -q -b DataMCComp3.C+\($d,\"w_Phi_star_mm\",1,$i,$n\)
  fi

  root -l -q -b DataMCComp3.C+\($d,\"w_Ht_b\",1,$i,$n\)
  root -l -q -b DataMCComp3.C+\($d,\"w_pt_Z_ee_b\",1,$i,$n\)
  root -l -q -b DataMCComp3.C+\($d,\"w_pt_Z_mm_b\",1,$i,$n\)

  if [ $n -eq 0 ]; then
    root -l -q -b DataMCComp3.C+\($d,\"w_delta_phi_mm_b\",1,$i,$n\)
    root -l -q -b DataMCComp3.C+\($d,\"w_delta_phi_ee_b\",1,$i,$n\)
    root -l -q -b DataMCComp3.C+\($d,\"w_mass_Zj_ee_b\",1,$i,$n\)
    root -l -q -b DataMCComp3.C+\($d,\"w_mass_Zj_mm_b\",1,$i,$n\)

    root -l -q -b DataMCComp3.C+\($d,\"w_y_Z_ee_b\",1,$i,$n\)
    root -l -q -b DataMCComp3.C+\($d,\"w_y_Z_mm_b\",1,$i,$n\)
    root -l -q -b DataMCComp3.C+\($d,\"w_Phi_star_ee_b\",1,$i,$n\)
    root -l -q -b DataMCComp3.C+\($d,\"w_Phi_star_mm_b\",1,$i,$n\)
  fi

  if [ $n -eq 2 ]; then
    root -l -q -b DataMCComp3.C+\($d,\"w_DR_bb\",1,$i,2\)
    root -l -q -b DataMCComp3.C+\($d,\"w_DR_eeb_min\",1,$i,2\)
    root -l -q -b DataMCComp3.C+\($d,\"w_DR_eeb_max\",1,$i,2\)
    root -l -q -b DataMCComp3.C+\($d,\"w_DR_mmb_min\",1,$i,2\)
    root -l -q -b DataMCComp3.C+\($d,\"w_DR_mmb_max\",1,$i,2\)
    root -l -q -b DataMCComp3.C+\($d,\"w_A_eeb\",1,$i,2\)
    root -l -q -b DataMCComp3.C+\($d,\"w_A_mmb\",1,$i,2\)
    root -l -q -b DataMCComp3.C+\($d,\"w_eebb_mass\",1,$i,2\)
    root -l -q -b DataMCComp3.C+\($d,\"w_mmbb_mass\",1,$i,2\)
    root -l -q -b DataMCComp3.C+\($d,\"w_bb_mass\",1,$i,2\)
    root -l -q -b DataMCComp3.C+\($d,\"w_delta_phi_2b\",1,$i,2\)
  fi

  i=$((i+1))
done

exit
