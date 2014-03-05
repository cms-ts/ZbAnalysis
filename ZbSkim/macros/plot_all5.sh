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

  root -l -q -b DataMCComp5.C+\($d,\"w_mass_ee_wide\",1,$i,1,$n\)
  root -l -q -b DataMCComp5.C+\($d,\"w_mass_mm_wide\",1,$i,1,$n\)
  root -l -q -b DataMCComp5.C+\($d,\"w_mass_ee_b_wide\",1,$i,1,$n\)
  root -l -q -b DataMCComp5.C+\($d,\"w_mass_mm_b_wide\",1,$i,1,$n\)

  root -l -q -b DataMCComp5.C+\($d,\"w_MET\",1,$i,1,$n\) 
  root -l -q -b DataMCComp5.C+\($d,\"w_MET_b\",1,$i,1,$n\)
  root -l -q -b DataMCComp5.C+\($d,\"w_MET_sign\",1,$i,1,$n\)
  root -l -q -b DataMCComp5.C+\($d,\"w_MET_sign_b\",1,$i,1,$n\)

i=$((i+1))
done

i=1
while [ $i -le 2 ]; do

  if [ $n -eq 0 ]; then

    root -l -q -b DataMCComp5.C+\($d,\"h_pu_weights\",1,$i,0,$n\)
    root -l -q -b DataMCComp5.C+\($d,\"h_recoVTX\",1,$i,0,$n\)
    root -l -q -b DataMCComp5.C+\($d,\"w_recoVTX\",1,$i,0,$n\)
    root -l -q -b DataMCComp5.C+\($d,\"h_tracks\",1,$i,0,$n\)
    root -l -q -b DataMCComp5.C+\($d,\"w_tracks\",1,$i,0,$n\)

    root -l -q -b DataMCComp5.C+\($d,\"w_first_ele_pt\",1,$i,0,$n\)
    root -l -q -b DataMCComp5.C+\($d,\"w_first_ele_eta\",1,$i,0,$n\)
    root -l -q -b DataMCComp5.C+\($d,\"w_second_ele_pt\",1,$i,0,$n\)
    root -l -q -b DataMCComp5.C+\($d,\"w_second_ele_eta\",1,$i,0,$n\)

    root -l -q -b DataMCComp5.C+\($d,\"w_first_muon_pt\",1,$i,0,$n\)
    root -l -q -b DataMCComp5.C+\($d,\"w_first_muon_eta\",1,$i,0,$n\)
    root -l -q -b DataMCComp5.C+\($d,\"w_second_muon_pt\",1,$i,0,$n\)
    root -l -q -b DataMCComp5.C+\($d,\"w_second_muon_eta\",1,$i,0,$n\)

    root -l -q -b DataMCComp5.C+\($d,\"w_mass_Zj_ee\",1,$i,0,$n\)
    root -l -q -b DataMCComp5.C+\($d,\"w_mass_Zj_mm\",1,$i,0,$n\)

    root -l -q -b DataMCComp5.C+\($d,\"w_numberOfZ\",1,$i,0,$n\)

    root -l -q -b DataMCComp5.C+\($d,\"w_mass_ee\",1,$i,0,$n\)
    root -l -q -b DataMCComp5.C+\($d,\"w_mass_mm\",1,$i,0,$n\)
    root -l -q -b DataMCComp5.C+\($d,\"w_pt_Z_ee\",1,$i,0,$n\)
    root -l -q -b DataMCComp5.C+\($d,\"w_pt_Z_mm\",1,$i,0,$n\)
    root -l -q -b DataMCComp5.C+\($d,\"w_y_Z_ee_b\",1,$i,0,$n\)
    root -l -q -b DataMCComp5.C+\($d,\"w_y_Z_mm_b\",1,$i,0,$n\)

    root -l -q -b DataMCComp5.C+\($d,\"w_Phi_star_ee\",1,$i,0,$n\)
    root -l -q -b DataMCComp5.C+\($d,\"w_Phi_star_mm\",1,$i,0,$n\)

    root -l -q -b DataMCComp5.C+\($d,\"w_jetmultiplicity\",1,$i,0,$n\)
    root -l -q -b DataMCComp5.C+\($d,\"w_first_jet_pt\",1,$i,0,$n\)
    root -l -q -b DataMCComp5.C+\($d,\"w_first_jet_eta\",1,$i,0,$n\)

    root -l -q -b DataMCComp5.C+\($d,\"w_second_jet_pt\",1,$i,0,$n\)
    root -l -q -b DataMCComp5.C+\($d,\"w_second_jet_eta\",1,$i,0,$n\)
    root -l -q -b DataMCComp5.C+\($d,\"w_third_jet_pt\",1,$i,0,$n\)
    root -l -q -b DataMCComp5.C+\($d,\"w_third_jet_eta\",1,$i,0,$n\)

    root -l -q -b DataMCComp5.C+\($d,\"w_secondvtx_N\",1,$i,0,$n\)

    root -l -q -b DataMCComp5.C+\($d,\"w_Ht\",1,$i,0,$n\)

    root -l -q -b DataMCComp5.C+\($d,\"h_scaleFactor_first_ele\",1,$i,0,$n\)
    root -l -q -b DataMCComp5.C+\($d,\"h_scaleFactor_second_ele\",1,$i,0,$n\)
    root -l -q -b DataMCComp5.C+\($d,\"h_scaleFactor_first_muon\",1,$i,0,$n\)
    root -l -q -b DataMCComp5.C+\($d,\"h_scaleFactor_second_muon\",1,$i,0,$n\)
  fi

  root -l -q -b DataMCComp5.C+\($d,\"w_mass_ee_wide\",1,$i,0,$n\)
  root -l -q -b DataMCComp5.C+\($d,\"w_mass_mm_wide\",1,$i,0,$n\)

  root -l -q -b DataMCComp5.C+\($d,\"w_mass_ee_b_wide\",1,$i,0,$n\)
  root -l -q -b DataMCComp5.C+\($d,\"w_mass_mm_b_wide\",1,$i,0,$n\)

  root -l -q -b DataMCComp5.C+\($d,\"w_delta_phi_ee\",1,$i,0,$n\)
  root -l -q -b DataMCComp5.C+\($d,\"w_delta_phi_mm\",1,$i,0,$n\)

  root -l -q -b DataMCComp5.C+\($d,\"w_mass_mm_b\",1,$i,0,$n\)
  root -l -q -b DataMCComp5.C+\($d,\"w_mass_ee_b\",1,$i,0,$n\)
  root -l -q -b DataMCComp5.C+\($d,\"w_pt_Z_ee_b\",1,$i,0,$n\)
  root -l -q -b DataMCComp5.C+\($d,\"w_pt_Z_mm_b\",1,$i,0,$n\)
  root -l -q -b DataMCComp5.C+\($d,\"w_delta_phi_ee_b\",1,$i,0,$n\)
  root -l -q -b DataMCComp5.C+\($d,\"w_delta_phi_mm_b\",1,$i,0,$n\)

  root -l -q -b DataMCComp5.C+\($d,\"w_Phi_star_ee_b\",1,$i,0,$n\)
  root -l -q -b DataMCComp5.C+\($d,\"w_Phi_star_mm_b\",1,$i,0,$n\)

  root -l -q -b DataMCComp5.C+\($d,\"w_single_pt_Z_ee_b\",1,$i,0,$n\)
  root -l -q -b DataMCComp5.C+\($d,\"w_single_pt_Z_mm_b\",1,$i,0,$n\)
  root -l -q -b DataMCComp5.C+\($d,\"w_single_bjet_pt\",1,$i,0,$n\)
  root -l -q -b DataMCComp5.C+\($d,\"w_single_bjet_eta\",1,$i,0,$n\)
  root -l -q -b DataMCComp5.C+\($d,\"w_single_delta_phi_ee_b\",1,$i,0,$n\)
  root -l -q -b DataMCComp5.C+\($d,\"w_single_delta_phi_mm_b\",1,$i,0,$n\)
  root -l -q -b DataMCComp5.C+\($d,\"w_single_Ht_b\",1,$i,0,$n\)
  root -l -q -b DataMCComp5.C+\($d,\"w_mass_Zj_ee_b\",1,$i,0,$n\)
  root -l -q -b DataMCComp5.C+\($d,\"w_mass_Zj_mm_b\",1,$i,0,$n\)
  
  root -l -q -b DataMCComp5.C+\($d,\"w_MET\",1,$i,0,$n\)
  root -l -q -b DataMCComp5.C+\($d,\"w_MET_sign\",1,$i,0,$n\)

  root -l -q -b DataMCComp5.C+\($d,\"w_MET_b\",1,$i,0,$n\)
  root -l -q -b DataMCComp5.C+\($d,\"w_MET_sign_b\",1,$i,0,$n\)

  root -l -q -b DataMCComp5.C+\($d,\"w_Ht_b\",1,$i,0,$n\)

  root -l -q -b DataMCComp5.C+\($d,\"w_first_jet_pt_b\",1,$i,0,$n\)
  root -l -q -b DataMCComp5.C+\($d,\"w_first_jet_eta_b\",1,$i,0,$n\)
  root -l -q -b DataMCComp5.C+\($d,\"w_second_jet_pt_b\",1,$i,0,$n\)
  root -l -q -b DataMCComp5.C+\($d,\"w_second_jet_eta_b\",1,$i,0,$n\)
  root -l -q -b DataMCComp5.C+\($d,\"w_third_jet_pt_b\",1,$i,0,$n\)
  root -l -q -b DataMCComp5.C+\($d,\"w_third_jet_eta_b\",1,$i,0,$n\)

  root -l -q -b DataMCComp5.C+\($d,\"w_delta_phi_2b\",1,$i,0,$n\)
  root -l -q -b DataMCComp5.C+\($d,\"w_bjetmultiplicity_exc\",1,$i,0,$n\)
  root -l -q -b DataMCComp5.C+\($d,\"w_bjetmultiplicity\",1,$i,0,$n\)
  root -l -q -b DataMCComp5.C+\($d,\"w_first_bjet_pt\",1,$i,0,$n\)
  root -l -q -b DataMCComp5.C+\($d,\"w_first_bjet_eta\",1,$i,0,$n\)
  root -l -q -b DataMCComp5.C+\($d,\"w_second_bjet_pt\",1,$i,0,$n\)
  root -l -q -b DataMCComp5.C+\($d,\"w_second_bjet_eta\",1,$i,0,$n\)
  root -l -q -b DataMCComp5.C+\($d,\"w_third_bjet_pt\",1,$i,0,$n\)
  root -l -q -b DataMCComp5.C+\($d,\"w_third_bjet_eta\",1,$i,0,$n\)

  root -l -q -b DataMCComp5.C+\($d,\"w_DR_eeb_min\",1,$i,0,$n\)
  root -l -q -b DataMCComp5.C+\($d,\"w_DR_eeb_max\",1,$i,0,$n\)
  root -l -q -b DataMCComp5.C+\($d,\"w_DR_mmb_min\",1,$i,0,$n\)
  root -l -q -b DataMCComp5.C+\($d,\"w_DR_mmb_max\",1,$i,0,$n\)

  root -l -q -b DataMCComp5.C+\($d,\"w_A_eeb\",1,$i,0,$n\)
  root -l -q -b DataMCComp5.C+\($d,\"w_A_mmb\",1,$i,0,$n\)

  root -l -q -b DataMCComp5.C+\($d,\"w_bb_mass\",1,$i,0,$n\)
  root -l -q -b DataMCComp5.C+\($d,\"w_eebb_mass\",1,$i,0,$n\)
  root -l -q -b DataMCComp5.C+\($d,\"w_mmbb_mass\",1,$i,0,$n\)

  root -l -q -b DataMCComp5.C+\($d,\"w_SVTX_mass_jet\",1,$i,0,$n\)
  root -l -q -b DataMCComp5.C+\($d,\"w_SVTX_mass_trk\",1,$i,0,$n\)
  root -l -q -b DataMCComp5.C+\($d,\"w_SVTX_mass\",1,$i,0,$n\)
  root -l -q -b DataMCComp5.C+\($d,\"w_SVTX_mass_sub\",1,$i,0,$n\)

  root -l -q -b DataMCComp5.C+\($d,\"w_secondvtx_N_zoom\",1,$i,0,$n\)

  root -l -q -b DataMCComp5.C+\($d,\"w_BJP\",1,$i,0,$n\)
  root -l -q -b DataMCComp5.C+\($d,\"w_BJP_sub\",1,$i,0,$n\)
  root -l -q -b DataMCComp5.C+\($d,\"w_JBP\",1,$i,0,$n\)

  root -l -q -b DataMCComp5.C+\($d,\"w_secondvtx_N_mass\",1,$i,0,$n\)

  root -l -q -b DataMCComp5.C+\($d,\"w_secondvtx_N_nomass\",1,$i,0,$n\)

  root -l -q -b DataMCComp5.C+\($d,\"w_BJP_mass\",1,$i,0,$n\)
  root -l -q -b DataMCComp5.C+\($d,\"w_JBP_mass\",1,$i,0,$n\)

  root -l -q -b DataMCComp5.C+\($d,\"w_BJP_nomass\",1,$i,0,$n\)
  root -l -q -b DataMCComp5.C+\($d,\"w_JBP_nomass\",1,$i,0,$n\)

  root -l -q -b DataMCComp5.C+\($d,\"w_BJP0\",1,$i,0,$n\)
  root -l -q -b DataMCComp5.C+\($d,\"w_BJP1\",1,$i,0,$n\)
  root -l -q -b DataMCComp5.C+\($d,\"w_BJP2\",1,$i,0,$n\)

  i=$((i+1))
done

exit
