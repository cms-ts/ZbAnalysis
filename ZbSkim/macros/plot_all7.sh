#!/bin/sh

export ROOT_HIST=0

unset PYTHIA8175DATA
unset G4NEUTRONXS

cd $CMS_PATH/slc6_amd64_gcc472/cms/cmssw/CMSSW_6_2_7
eval `scramv1 runtime -sh`
cd -

n=0
if [ ! -z "$1" ]; then
  n=$1
fi

i=1
while [ $i -le 2 ]; do

  root -l -q -b DataMCComp7.C+\(\"w_first_bjet_pt\",1,$i,0,$n\)
  root -l -q -b DataMCComp7.C+\(\"w_first_bjet_eta\",1,$i,0,$n\)
  root -l -q -b DataMCComp7.C+\(\"w_first_bjet_eta_abs\",1,$i,0,$n\)

  if [ $n -eq 2 ]; then
    root -l -q -b DataMCComp7.C+\(\"w_second_bjet_pt\",1,$i,0,$n\)
    root -l -q -b DataMCComp7.C+\(\"w_second_bjet_eta\",1,$i,0,$n\)
    root -l -q -b DataMCComp7.C+\(\"w_second_bjet_eta_abs\",1,$i,0,$n\)
  fi

  if [ $n -eq 0 ]; then
    root -l -q -b DataMCComp7.C+\(\"w_pt_Z_ee\",1,$i,0,$n\)
    root -l -q -b DataMCComp7.C+\(\"w_pt_Z_mm\",1,$i,0,$n\)
    root -l -q -b DataMCComp7.C+\(\"w_Ht\",1,$i,0,$n\)
    root -l -q -b DataMCComp7.C+\(\"w_delta_phi_ee\",1,$i,0,$n\)
    root -l -q -b DataMCComp7.C+\(\"w_delta_phi_mm\",1,$i,0,$n\)
    root -l -q -b DataMCComp7.C+\(\"w_mass_Zj_ee\",1,$i,0,$n\)
    root -l -q -b DataMCComp7.C+\(\"w_mass_Zj_mm\",1,$i,0,$n\)
  fi

  if [ $n -ne 0 ]; then
    root -l -q -b DataMCComp7.C+\(\"w_pt_Z_ee_b\",1,$i,0,$n\)
    root -l -q -b DataMCComp7.C+\(\"w_pt_Z_mm_b\",1,$i,0,$n\)
    root -l -q -b DataMCComp7.C+\(\"w_Ht_b\",1,$i,0,$n\)
  fi

  if [ $n -eq 2 ]; then
    root -l -q -b DataMCComp7.C+\(\"w_DR_bb\",1,$i,0,$n\)
    root -l -q -b DataMCComp7.C+\(\"w_DR_eeb_min\",1,$i,0,$n\)
    root -l -q -b DataMCComp7.C+\(\"w_DR_eeb_max\",1,$i,0,$n\)
    root -l -q -b DataMCComp7.C+\(\"w_DR_mmb_min\",1,$i,0,$n\)
    root -l -q -b DataMCComp7.C+\(\"w_DR_mmb_max\",1,$i,0,$n\)
    root -l -q -b DataMCComp7.C+\(\"w_A_eeb\",1,$i,0,$n\)
    root -l -q -b DataMCComp7.C+\(\"w_A_mmb\",1,$i,0,$n\)
    root -l -q -b DataMCComp7.C+\(\"w_eebb_mass\",1,$i,0,$n\)
    root -l -q -b DataMCComp7.C+\(\"w_mmbb_mass\",1,$i,0,$n\)
    root -l -q -b DataMCComp7.C+\(\"w_bb_mass\",1,$i,0,$n\)
    root -l -q -b DataMCComp7.C+\(\"w_delta_phi_2b\",1,$i,0,$n\)
  fi

  i=$((i+1))
done

i=1
while [ $i -le 2 ]; do

  if [ $n -eq 0 ]; then
    root -l -q -b DataMCComp7.C+\(\"w_first_bjet_pt\",1,$i,1,$n\)
    root -l -q -b DataMCComp7.C+\(\"w_first_bjet_eta\",1,$i,1,$n\)
    root -l -q -b DataMCComp7.C+\(\"w_first_bjet_eta_abs\",1,$i,1,$n\)
  fi

  if [ $n -eq 0 ]; then
    root -l -q -b DataMCComp7.C+\(\"w_pt_Z_ee\",1,$i,1,$n\)
    root -l -q -b DataMCComp7.C+\(\"w_pt_Z_mm\",1,$i,1,$n\)
    root -l -q -b DataMCComp7.C+\(\"w_Ht\",1,$i,1,$n\)
    root -l -q -b DataMCComp7.C+\(\"w_delta_phi_ee\",1,$i,1,$n\)
    root -l -q -b DataMCComp7.C+\(\"w_delta_phi_mm\",1,$i,1,$n\)
  fi

  i=$((i+1))
done

exit
