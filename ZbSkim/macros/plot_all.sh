#!/bin/sh

cd $CMS_PATH/slc5_amd64_gcc472/cms/cmssw/CMSSW_6_1_0
eval `scramv1 runtime -sh`
cd -

i=1
while [ $i -le 2 ]; do

  root -l -q -b DataMCComp.C\(\"h_pu_weights\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_secondvtx_N\",1,$i\)
  root -l -q -b DataMCComp.C\(\"recoVTX\",1,$i\)
  root -l -q -b DataMCComp.C\(\"recoVTXw\",1,$i\)
  root -l -q -b DataMCComp.C\(\"h_tracks\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_tracks\",1,$i\)

  root -l -q -b DataMCComp.C\(\"first_ele_pt\",1,$i\)
  root -l -q -b DataMCComp.C\(\"first_ele_eta\",1,$i\)
  root -l -q -b DataMCComp.C\(\"first_muon_pt\",1,$i\)
  root -l -q -b DataMCComp.C\(\"first_muon_eta\",1,$i\)

  root -l -q -b DataMCComp.C\(\"numberOfZ\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_mass_ee\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_mass_mm\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_mass_mm_b\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_mass_ee_b\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_pt_Z_ee\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_pt_Z_mm\",1,$i\)

  root -l -q -b DataMCComp.C\(\"w_jetmultiplicity\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_first_jet_pt\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_first_jet_eta\",1,$i\)

  root -l -q -b DataMCComp.C\(\"w_MET\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_MET_sign\",1,$i\)

  root -l -q -b DataMCComp.C\(\"w_Ht\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_Ht_b\",1,$i\)

  root -l -q -b DataMCComp.C\(\"w_jetmultiplicity_b\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_first_jet_pt_b\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_first_jet_eta_b\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_first_bjet_pt\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_first_bjet_eta\",1,$i\)

  root -l -q -b DataMCComp.C\(\"w_delta_phi_ee\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_delta_phi_mm\",1,$i\)

  root -l -q -b DataMCComp.C\(\"w_SVTX_mass_jet\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_SVTX_mass_trk\",1,$i\)
  root -l -q -b DataMCComp.C\(\"w_SVTX_mass\",1,$i\)

  i=$((i+1))
done

i=1
while [ $i -le 2 ]; do

  root -l -q -b DataMCComp.C\(\"w_MET\",1,$i,0,1\)

  root -l -q -b DataMCComp.C\(\"w_SVTX_mass\",1,$i,1,0\)
  root -l -q -b DataMCComp.C\(\"w_SVTX_mass\",1,$i,0,2\)
  root -l -q -b DataMCComp.C\(\"w_SVTX_mass\",1,$i,1,2\)
  
  root -l -q -b DataMCComp.C\(\"w_secondvtx_N\",1,$i,1,0\)
  root -l -q -b DataMCComp.C\(\"w_secondvtx_N\",1,$i,0,2\)
  root -l -q -b DataMCComp.C\(\"w_secondvtx_N\",1,$i,1,2\)

  i=$((i+1))
done

exit
