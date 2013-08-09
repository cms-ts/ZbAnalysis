#!/bin/sh

export ROOT_HIST=0

cd $CMS_PATH/slc5_amd64_gcc472/cms/cmssw/CMSSW_6_1_0
eval `scramv1 runtime -sh`
cd -

i=1
while [ $i -le 2 ]; do

  root -l -q -b DataMCComp.C\(\"h_pu_weights\",1,$i\)
  root -l -q -b DataMCComp.C\(\"h_recoVTX\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_recoVTX\",1,$i\)
  root -l -q -b DataMCComp.C\(\"h_tracks\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_tracks\",1,$i\)

  root -l -q -b DataMCComp.C\(\"w_first_ele_pt\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_first_ele_eta\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_second_ele_pt\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_second_ele_eta\",1,$i\)

  root -l -q -b DataMCComp.C\(\"w_first_muon_pt\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_first_muon_eta\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_second_muon_pt\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_second_muon_eta\",1,$i\)

  root -l -q -b DataMCComp.C\(\"w_numberOfZ\",1,$i\)

  root -l -q -b DataMCComp.C\(\"w_mass_ee_wide\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_mass_mm_wide\",1,$i\)

  root -l -q -b DataMCComp.C\(\"w_mass_ee_b_wide\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_mass_mm_b_wide\",1,$i\)

  root -l -q -b DataMCComp.C\(\"w_mass_ee\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_mass_mm\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_pt_Z_ee\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_pt_Z_mm\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_delta_phi_ee\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_delta_phi_mm\",1,$i\)

  root -l -q -b DataMCComp.C\(\"w_mass_mm_b\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_mass_ee_b\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_pt_Z_ee_b\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_pt_Z_mm_b\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_delta_phi_ee_b\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_delta_phi_mm_b\",1,$i\)

  root -l -q -b DataMCComp.C\(\"w_jetmultiplicity\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_first_jet_pt\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_first_jet_eta\",1,$i\)
  
  root -l -q -b DataMCComp.C\(\"w_single_pt_Z_ee_b\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_single_pt_Z_mm_b\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_single_bjet_pt\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_single_bjet_eta\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_single_delta_phi_ee_b\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_single_delta_phi_mm_b\",1,$i\)
  
  root -l -q -b DataMCComp.C\(\"w_second_jet_pt\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_second_jet_eta\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_third_jet_pt\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_third_jet_eta\",1,$i\)

  root -l -q -b DataMCComp.C\(\"w_MET\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_MET_sign\",1,$i\)

  root -l -q -b DataMCComp.C\(\"w_MET_b\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_MET_sign_b\",1,$i\)

  root -l -q -b DataMCComp.C\(\"w_Ht\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_Ht_b\",1,$i\)

  root -l -q -b DataMCComp.C\(\"w_first_jet_pt_b\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_first_jet_eta_b\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_second_jet_pt_b\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_second_jet_eta_b\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_third_jet_pt_b\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_third_jet_eta_b\",1,$i\)

  root -l -q -b DataMCComp.C\(\"w_bjetmultiplicity\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_first_bjet_pt\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_first_bjet_eta\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_second_bjet_pt\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_second_bjet_eta\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_third_bjet_pt\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_third_bjet_eta\",1,$i\)

  root -l -q -b DataMCComp.C\(\"h_scaleFactor_first_ele\",1,$i\)
  root -l -q -b DataMCComp.C\(\"h_scaleFactor_second_ele\",1,$i\)
  root -l -q -b DataMCComp.C\(\"h_scaleFactor_first_muon\",1,$i\)
  root -l -q -b DataMCComp.C\(\"h_scaleFactor_second_muon\",1,$i\)

  root -l -q -b DataMCComp.C\(\"w_SVTX_mass_jet\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_SVTX_mass_trk\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_SVTX_mass\",1,$i\)

  root -l -q -b DataMCComp.C\(\"w_secondvtx_N\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_secondvtx_N_zoom\",1,$i\)

  root -l -q -b DataMCComp.C\(\"w_BJP\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_JBP\",1,$i\)

  root -l -q -b DataMCComp.C\(\"w_secondvtx_N_mass\",1,$i\)

  root -l -q -b DataMCComp.C\(\"w_BJP_mass\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_JBP_mass\",1,$i\)

  i=$((i+1))
done

i=1
while [ $i -le 2 ]; do

  root -l -q -b DataMCComp.C\(\"w_MET\",1,$i,0,1\)
  root -l -q -b DataMCComp.C\(\"w_MET_sign\",1,$i,0,1\)

  root -l -q -b DataMCComp.C\(\"w_MET_b\",1,$i,0,1\)
  root -l -q -b DataMCComp.C\(\"w_MET_sign_b\",1,$i,0,1\)

  root -l -q -b DataMCComp.C\(\"w_mass_ee_wide\",1,$i,0,2\)
  root -l -q -b DataMCComp.C\(\"w_mass_mm_wide\",1,$i,0,2\)

  root -l -q -b DataMCComp.C\(\"w_mass_ee_b_wide\",1,$i,0,2\)
  root -l -q -b DataMCComp.C\(\"w_mass_mm_b_wide\",1,$i,0,2\)

  root -l -q -b DataMCComp.C\(\"w_SVTX_mass\",1,$i,1,0\)
  root -l -q -b DataMCComp.C\(\"w_SVTX_mass\",1,$i,0,3\)
  root -l -q -b DataMCComp.C\(\"w_SVTX_mass\",1,$i,1,3\)
  
  root -l -q -b DataMCComp.C\(\"w_secondvtx_N_zoom\",1,$i,1,0\)
  root -l -q -b DataMCComp.C\(\"w_secondvtx_N_zoom\",1,$i,0,3\)
  root -l -q -b DataMCComp.C\(\"w_secondvtx_N_zoom\",1,$i,1,3\)

  root -l -q -b DataMCComp.C\(\"w_secondvtx_N_mass\",1,$i,1,0\)
  root -l -q -b DataMCComp.C\(\"w_secondvtx_N_mass\",1,$i,0,3\)
  root -l -q -b DataMCComp.C\(\"w_secondvtx_N_mass\",1,$i,1,3\)

  root -l -q -b DataMCComp.C\(\"w_BJP\",1,$i,1,0\)
  root -l -q -b DataMCComp.C\(\"w_BJP\",1,$i,0,3\)
  root -l -q -b DataMCComp.C\(\"w_BJP\",1,$i,1,3\)

  root -l -q -b DataMCComp.C\(\"w_BJP_mass\",1,$i,1,0\)
  root -l -q -b DataMCComp.C\(\"w_BJP_mass\",1,$i,0,3\)
  root -l -q -b DataMCComp.C\(\"w_BJP_mass\",1,$i,1,3\)

  root -l -q -b DataMCComp.C\(\"w_JBP\",1,$i,1,0\)
  root -l -q -b DataMCComp.C\(\"w_JBP\",1,$i,0,3\)
  root -l -q -b DataMCComp.C\(\"w_JBP\",1,$i,1,3\)

  root -l -q -b DataMCComp.C\(\"w_JBP_mass\",1,$i,1,0\)
  root -l -q -b DataMCComp.C\(\"w_JBP_mass\",1,$i,0,3\)
  root -l -q -b DataMCComp.C\(\"w_JBP_mass\",1,$i,1,3\)

  i=$((i+1))
done

exit
