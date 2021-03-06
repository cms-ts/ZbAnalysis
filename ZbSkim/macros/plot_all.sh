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

  root -l -q -b DataMCComp.C+\($d,\"w_MET\",1,$i,0,1,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_MET_b\",1,$i,0,1,$n\)

  if [ $n -eq 2 ]; then
    root -l -q -b DataMCComp.C+\($d,\"w_BJP\",1,$i,1,4,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_BJP\",1,$i,0,4,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_SVTX_mass\",1,$i,1,4,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_SVTX_mass\",1,$i,0,4,$n\)
  fi

  i=$((i+1))
done

root -l -q -b DataMCComp.C+\($d,\"w_MET\",1,3,0,1,$n\)
root -l -q -b DataMCComp.C+\($d,\"w_MET_b\",1,3,0,1,$n\)

i=1
while [ $i -le 2 ]; do

  if [ $n -eq 0 ]; then
    root -l -q -b DataMCComp.C+\($d,\"w_Ht\",1,$i,0,0,$n\)

    root -l -q -b DataMCComp.C+\($d,\"h_pu_weights\",1,$i,0,0,$n\)
    root -l -q -b DataMCComp.C+\($d,\"h_recoVTX\",1,$i,0,0,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_recoVTX\",1,$i,0,0,$n\)
    root -l -q -b DataMCComp.C+\($d,\"h_tracks\",1,$i,0,0,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_tracks\",1,$i,0,0,$n\)

    root -l -q -b DataMCComp.C+\($d,\"w_first_ele_pt\",1,$i,0,0,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_first_ele_eta\",1,$i,0,0,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_second_ele_pt\",1,$i,0,0,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_second_ele_eta\",1,$i,0,0,$n\)

    root -l -q -b DataMCComp.C+\($d,\"w_first_muon_pt\",1,$i,0,0,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_first_muon_eta\",1,$i,0,0,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_second_muon_pt\",1,$i,0,0,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_second_muon_eta\",1,$i,0,0,$n\)

    root -l -q -b DataMCComp.C+\($d,\"w_numberOfZ\",1,$i,0,0,$n\)
    
    root -l -q -b DataMCComp.C+\($d,\"w_mass_ee_wide\",1,$i,0,0,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_mass_mm_wide\",1,$i,0,0,$n\)

    root -l -q -b DataMCComp.C+\($d,\"w_mass_ee\",1,$i,0,0,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_pt_Z_ee\",1,$i,0,0,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_y_Z_ee\",1,$i,0,0,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_y_Z_ee_abs\",1,$i,0,0,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_mass_mm\",1,$i,0,0,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_pt_Z_mm\",1,$i,0,0,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_y_Z_mm\",1,$i,0,0,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_y_Z_mm_abs\",1,$i,0,0,$n\)

    root -l -q -b DataMCComp.C+\($d,\"w_Phi_star_ee\",1,$i,0,0,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_Phi_star_mm\",1,$i,0,0,$n\)

    root -l -q -b DataMCComp.C+\($d,\"w_mass_Zj_ee\",1,$i,0,0,$n\)
  
    root -l -q -b DataMCComp.C+\($d,\"w_secondvtx_N\",1,$i,0,0,$n\)

    root -l -q -b DataMCComp.C+\($d,\"w_mass_Zj_mm\",1,$i,0,0,$n\)

    root -l -q -b DataMCComp.C+\($d,\"w_jetmultiplicity\",1,$i,0,0,$n\)

    root -l -q -b DataMCComp.C+\($d,\"w_first_jet_pt\",1,$i,0,0,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_first_jet_eta\",1,$i,0,0,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_first_jet_eta_abs\",1,$i,0,0,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_second_jet_pt\",1,$i,0,0,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_second_jet_eta\",1,$i,0,0,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_second_jet_eta_abs\",1,$i,0,0,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_third_jet_pt\",1,$i,0,0,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_third_jet_eta\",1,$i,0,0,$n\)

    root -l -q -b DataMCComp.C+\($d,\"h_scaleFactor_first_ele\",1,$i,0,0,$n\)
    root -l -q -b DataMCComp.C+\($d,\"h_scaleFactor_second_ele\",1,$i,0,0,$n\)
    root -l -q -b DataMCComp.C+\($d,\"h_scaleFactor_first_muon\",1,$i,0,0,$n\)
    root -l -q -b DataMCComp.C+\($d,\"h_scaleFactor_second_muon\",1,$i,0,0,$n\)
 
  fi

  root -l -q -b DataMCComp.C+\($d,\"w_mass_ee_b_wide\",1,$i,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_mass_mm_b_wide\",1,$i,0,0,$n\)

  root -l -q -b DataMCComp.C+\($d,\"w_delta_phi_ee\",1,$i,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_delta_phi_mm\",1,$i,0,0,$n\)

  root -l -q -b DataMCComp.C+\($d,\"w_mass_ee_b\",1,$i,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_mass_mm_b\",1,$i,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_pt_Z_ee_b\",1,$i,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_pt_Z_mm_b\",1,$i,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_y_Z_ee_b\",1,$i,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_y_Z_ee_b_abs\",1,$i,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_y_Z_mm_b\",1,$i,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_y_Z_mm_b_abs\",1,$i,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_delta_phi_ee_b\",1,$i,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_delta_phi_mm_b\",1,$i,0,0,$n\)

  root -l -q -b DataMCComp.C+\($d,\"w_Phi_star_ee_b\",1,$i,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_Phi_star_mm_b\",1,$i,0,0,$n\)
  
  if [ $n -ne 2 ]; then 
    root -l -q -b DataMCComp.C+\($d,\"w_single_pt_Z_ee_b\",1,$i,0,0,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_single_pt_Z_mm_b\",1,$i,0,0,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_single_bjet_pt\",1,$i,0,0,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_single_bjet_eta\",1,$i,0,0,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_single_delta_phi_ee_b\",1,$i,0,0,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_single_delta_phi_mm_b\",1,$i,0,0,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_single_Ht_b\",1,$i,0,0,$n\)
  fi  

  root -l -q -b DataMCComp.C+\($d,\"w_MET\",1,$i,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_MET_b\",1,$i,0,0,$n\)

  root -l -q -b DataMCComp.C+\($d,\"w_MET_sign\",1,$i,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_MET_sign_b\",1,$i,0,0,$n\)

  root -l -q -b DataMCComp.C+\($d,\"w_Ht_b\",1,$i,0,0,$n\)
 
  root -l -q -b DataMCComp.C+\($d,\"w_mass_Zj_ee_b\",1,$i,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_mass_Zj_mm_b\",1,$i,0,0,$n\)

  root -l -q -b DataMCComp.C+\($d,\"w_first_jet_pt_b\",1,$i,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_first_jet_eta_b\",1,$i,0,0,$n\)

  if [ $n -ne 1 ]; then 
    root -l -q -b DataMCComp.C+\($d,\"w_second_jet_pt_b\",1,$i,0,0,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_second_jet_eta_b\",1,$i,0,0,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_third_jet_pt_b\",1,$i,0,0,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_third_jet_eta_b\",1,$i,0,0,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_SVTX_mass_sub\",1,$i,0,0,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_BJP_sub\",1,$i,0,0,$n\)
  fi

  root -l -q -b DataMCComp.C+\($d,\"w_bjetmultiplicity_exc\",1,$i,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_bjetmultiplicity\",1,$i,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_first_bjet_pt\",1,$i,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_first_bjet_pt_SVTX\",1,$i,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_first_bjet_eta\",1,$i,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_first_bjet_eta_abs\",1,$i,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_SVTX_mass_jet\",1,$i,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_SVTX_mass_trk\",1,$i,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_SVTX_mass\",1,$i,0,0,$n\)
 
  root -l -q -b DataMCComp.C+\($d,\"w_flightd\",1,$i,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_flightd_sig\",1,$i,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_dxy\",1,$i,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_dxy_sig\",1,$i,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_SV_NTracks\",1,$i,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_N_SV\",1,$i,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_delta_j\",1,$i,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_delta_j_n\",1,$i,0,0,$n\)

  root -l -q -b DataMCComp.C+\($d,\"w_secondvtx_N_zoom\",1,$i,0,0,$n\)

  root -l -q -b DataMCComp.C+\($d,\"w_BJP\",1,$i,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_JBP\",1,$i,0,0,$n\)

  root -l -q -b DataMCComp.C+\($d,\"w_secondvtx_N_mass\",1,$i,0,0,$n\)

  root -l -q -b DataMCComp.C+\($d,\"w_BJP_mass\",1,$i,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_JBP_mass\",1,$i,0,0,$n\)

  if [ $n -ne 2 ]; then
    root -l -q -b DataMCComp.C+\($d,\"w_BJP0\",1,$i,0,0,$n\)
  fi

  if [ $n -ne 1 ]; then
    root -l -q -b DataMCComp.C+\($d,\"w_BJP1\",1,$i,0,0,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_BJP2\",1,$i,0,0,$n\)
  fi

  if [ $n -eq 2 ]; then
    root -l -q -b DataMCComp.C+\($d,\"w_delta_phi_2b\",1,$i,0,0,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_DR_bb\",1,$i,0,0,$n\)

    root -l -q -b DataMCComp.C+\($d,\"w_DR_eeb_min\",1,$i,0,0,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_DR_eeb_max\",1,$i,0,0,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_DR_mmb_min\",1,$i,0,0,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_DR_mmb_max\",1,$i,0,0,$n\)

    root -l -q -b DataMCComp.C+\($d,\"w_A_eeb\",1,$i,0,0,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_A_mmb\",1,$i,0,0,$n\)

    root -l -q -b DataMCComp.C+\($d,\"w_second_bjet_pt\",1,$i,0,0,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_second_bjet_eta\",1,$i,0,0,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_second_bjet_eta_abs\",1,$i,0,0,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_third_bjet_pt\",1,$i,0,0,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_third_bjet_eta\",1,$i,0,0,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_third_bjet_eta_abs\",1,$i,0,0,$n\)

    root -l -q -b DataMCComp.C+\($d,\"w_bb_mass\",1,$i,0,0,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_eebb_mass\",1,$i,0,0,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_mmbb_mass\",1,$i,0,0,$n\)
  fi

  i=$((i+1))
done

if [ $n -eq 0 ]; then
  root -l -q -b DataMCComp.C+\($d,\"w_secondvtx_N\",1,3,0,0,$n\)
  
  root -l -q -b DataMCComp.C+\($d,\"w_Ht\",1,3,0,0,$n\)

  root -l -q -b DataMCComp.C+\($d,\"w_mass_em_wide\",1,3,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_mass_em\",1,3,0,0,$n\)
 
  root -l -q -b DataMCComp.C+\($d,\"w_pt_Z_em\",1,3,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_y_Z_em\",1,3,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_y_Z_em_abs\",1,3,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_mass_Zj_em\",1,3,0,0,$n\)

  root -l -q -b DataMCComp.C+\($d,\"w_Phi_star_em\",1,3,0,0,$n\)

  root -l -q -b DataMCComp.C+\($d,\"w_jetmultiplicity\",1,3,0,0,$n\)

  root -l -q -b DataMCComp.C+\($d,\"w_first_jet_pt\",1,3,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_first_jet_eta\",1,3,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_first_jet_eta_abs\",1,3,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_second_jet_pt\",1,3,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_second_jet_eta\",1,3,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_second_jet_eta_abs\",1,3,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_third_jet_pt\",1,3,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_third_jet_eta\",1,3,0,0,$n\)
 
fi

root -l -q -b DataMCComp.C+\($d,\"w_mass_em_b_wide\",1,3,0,0,$n\)

root -l -q -b DataMCComp.C+\($d,\"w_delta_phi_em\",1,3,0,0,$n\)

root -l -q -b DataMCComp.C+\($d,\"w_mass_em_b\",1,3,0,0,$n\)
root -l -q -b DataMCComp.C+\($d,\"w_pt_Z_em_b\",1,3,0,0,$n\)
root -l -q -b DataMCComp.C+\($d,\"w_y_Z_em_b\",1,3,0,0,$n\)
root -l -q -b DataMCComp.C+\($d,\"w_y_Z_em_b_abs\",1,3,0,0,$n\)
root -l -q -b DataMCComp.C+\($d,\"w_delta_phi_em_b\",1,3,0,0,$n\)

if [ $n -ne 2 ]; then
  root -l -q -b DataMCComp.C+\($d,\"w_single_pt_Z_em_b\",1,3,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_single_bjet_pt\",1,3,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_single_bjet_eta\",1,3,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_single_delta_phi_em_b\",1,3,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_single_Ht_b\",1,3,0,0,$n\)
fi

root -l -q -b DataMCComp.C+\($d,\"w_Phi_star_em_b\",1,3,0,0,$n\)

root -l -q -b DataMCComp.C+\($d,\"w_MET\",1,3,0,0,$n\)
root -l -q -b DataMCComp.C+\($d,\"w_MET_sign\",1,3,0,0,$n\)

root -l -q -b DataMCComp.C+\($d,\"w_MET_b\",1,3,0,0,$n\)
root -l -q -b DataMCComp.C+\($d,\"w_MET_sign_b\",1,3,0,0,$n\)

root -l -q -b DataMCComp.C+\($d,\"w_Ht_b\",1,3,0,0,$n\)
  
root -l -q -b DataMCComp.C+\($d,\"w_mass_Zj_em_b\",1,3,0,0,$n\)

if [ $n -eq 2 ]; then
  root -l -q -b DataMCComp.C+\($d,\"w_delta_phi_2b\",1,3,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_DR_bb\",1,3,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_bjetmultiplicity_exc\",1,3,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_bjetmultiplicity\",1,3,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_first_bjet_pt\",1,3,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_first_bjet_pt_SVTX\",1,3,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_first_bjet_eta\",1,3,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_first_bjet_eta_abs\",1,3,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_second_bjet_pt\",1,3,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_second_bjet_eta\",1,3,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_second_bjet_eta_abs\",1,3,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_third_bjet_pt\",1,3,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_third_bjet_eta\",1,3,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_third_bjet_eta_abs\",1,3,0,0,$n\)
fi

root -l -q -b DataMCComp.C+\($d,\"w_flightd\",1,3,0,0,$n\)
root -l -q -b DataMCComp.C+\($d,\"w_flightd_sig\",1,3,0,0,$n\)
root -l -q -b DataMCComp.C+\($d,\"w_dxy\",1,3,0,0,$n\)
root -l -q -b DataMCComp.C+\($d,\"w_dxy_sig\",1,3,0,0,$n\)
root -l -q -b DataMCComp.C+\($d,\"w_SV_NTracks\",1,3,0,0,$n\)
root -l -q -b DataMCComp.C+\($d,\"w_N_SV\",1,3,0,0,$n\)
root -l -q -b DataMCComp.C+\($d,\"w_delta_j\",1,3,0,0,$n\)
root -l -q -b DataMCComp.C+\($d,\"w_delta_j_n\",1,3,0,0,$n\)

if [ $n -eq 2 ]; then
  root -l -q -b DataMCComp.C+\($d,\"w_DR_emb_min\",1,3,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_DR_emb_max\",1,3,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_A_emb\",1,3,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_bb_mass\",1,3,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_embb_mass\",1,3,0,0,$n\)
fi

root -l -q -b DataMCComp.C+\($d,\"w_SVTX_mass_jet\",1,3,0,0,$n\)
root -l -q -b DataMCComp.C+\($d,\"w_SVTX_mass_trk\",1,3,0,0,$n\)
root -l -q -b DataMCComp.C+\($d,\"w_SVTX_mass\",1,3,0,0,$n\)

root -l -q -b DataMCComp.C+\($d,\"w_secondvtx_N_zoom\",1,3,0,0,$n\)

root -l -q -b DataMCComp.C+\($d,\"w_BJP\",1,3,0,0,$n\)
root -l -q -b DataMCComp.C+\($d,\"w_JBP\",1,3,0,0,$n\)

root -l -q -b DataMCComp.C+\($d,\"w_secondvtx_N_mass\",1,3,0,0,$n\)

root -l -q -b DataMCComp.C+\($d,\"w_BJP_mass\",1,3,0,0,$n\)
root -l -q -b DataMCComp.C+\($d,\"w_JBP_mass\",1,3,0,0,$n\)

if [ $n -ne 2 ]; then
  root -l -q -b DataMCComp.C+\($d,\"w_BJP0\",1,3,0,0,$n\)
fi

if [ $n -ne 1 ]; then
  root -l -q -b DataMCComp.C+\($d,\"w_BJP1\",1,3,0,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_BJP2\",1,3,0,0,$n\)
fi

i=1
while [ $i -le 2 ]; do

  if [ $n -eq 0 ]; then
    root -l -q -b DataMCComp.C+\($d,\"w_mass_ee_wide\",1,$i,0,2,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_mass_mm_wide\",1,$i,0,2,$n\)
  fi

  root -l -q -b DataMCComp.C+\($d,\"w_MET_sign\",1,$i,0,1,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_MET_sign_b\",1,$i,0,1,$n\)

  root -l -q -b DataMCComp.C+\($d,\"w_mass_ee_b_wide\",1,$i,0,1,$n\)

  root -l -q -b DataMCComp.C+\($d,\"w_mass_mm_b_wide\",1,$i,0,1,$n\)

  if [ $n -ne 2 ]; then
    root -l -q -b DataMCComp.C+\($d,\"w_SVTX_mass\",1,$i,1,0,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_SVTX_mass\",1,$i,0,3,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_SVTX_mass\",1,$i,1,3,$n\)
  fi

  root -l -q -b DataMCComp.C+\($d,\"w_secondvtx_N_zoom\",1,$i,1,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_secondvtx_N_zoom\",1,$i,0,3,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_secondvtx_N_zoom\",1,$i,1,3,$n\)

  root -l -q -b DataMCComp.C+\($d,\"w_secondvtx_N_mass\",1,$i,1,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_secondvtx_N_mass\",1,$i,0,3,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_secondvtx_N_mass\",1,$i,1,3,$n\)

  root -l -q -b DataMCComp.C+\($d,\"w_BJP\",1,$i,1,0,$n\)

  if [ $n -ne 2 ]; then
    root -l -q -b DataMCComp.C+\($d,\"w_BJP\",1,$i,0,3,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_BJP\",1,$i,1,3,$n\)
  fi

  root -l -q -b DataMCComp.C+\($d,\"w_BJP_mass\",1,$i,1,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_BJP_mass\",1,$i,0,3,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_BJP_mass\",1,$i,1,3,$n\)

  root -l -q -b DataMCComp.C+\($d,\"w_JBP\",1,$i,1,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_JBP\",1,$i,0,3,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_JBP\",1,$i,1,3,$n\)

  root -l -q -b DataMCComp.C+\($d,\"w_JBP_mass\",1,$i,1,0,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_JBP_mass\",1,$i,0,3,$n\)
  root -l -q -b DataMCComp.C+\($d,\"w_JBP_mass\",1,$i,1,3,$n\)
 
  if [ $n -ne 2 ]; then 
    root -l -q -b DataMCComp.C+\($d,\"w_BJP0\",1,$i,1,0,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_BJP0\",1,$i,0,3,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_BJP0\",1,$i,1,3,$n\)
  fi

  if [ $n -ne 1 ]; then
    root -l -q -b DataMCComp.C+\($d,\"w_BJP1\",1,$i,0,3,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_BJP1\",1,$i,1,3,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_BJP1\",1,$i,1,0,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_BJP2\",1,$i,1,3,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_BJP2\",1,$i,0,3,$n\)
    root -l -q -b DataMCComp.C+\($d,\"w_BJP2\",1,$i,1,0,$n\) 
  fi
  

  i=$((i+1))
done

root -l -q -b DataMCComp.C+\($d,\"w_MET_sign\",1,3,0,1,$n\)
root -l -q -b DataMCComp.C+\($d,\"w_MET_sign_b\",1,3,0,1,$n\)

root -l -q -b DataMCComp.C+\($d,\"w_mass_em_wide\",1,3,0,1,$n\)
root -l -q -b DataMCComp.C+\($d,\"w_mass_em_b_wide\",1,3,0,1,$n\)

root -l -q -b DataMCComp.C+\($d,\"w_JBP\",1,3,0,1,$n\)
root -l -q -b DataMCComp.C+\($d,\"w_BJP\",1,3,0,1,$n\)
root -l -q -b DataMCComp.C+\($d,\"w_SVTX_mass\",1,3,0,1,$n\)

exit
