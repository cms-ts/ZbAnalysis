#!/bin/sh

export ROOT_HIST=0

cd $CMS_PATH/slc5_amd64_gcc472/cms/cmssw/CMSSW_6_1_1
eval `scramv1 runtime -sh`
cd -

i=1
while [ $i -le 2 ]; do

  root -l -q -b DataMCComp7.C\(\"w_first_bjet_pt\",1,$i,0\)
  root -l -q -b DataMCComp7.C\(\"w_first_bjet_eta\",1,$i,0\)
  root -l -q -b DataMCComp7.C\(\"w_pt_Z_ee\",1,$i,0\)
  root -l -q -b DataMCComp7.C\(\"w_pt_Z_mm\",1,$i,0\)
  root -l -q -b DataMCComp7.C\(\"w_Ht\",1,$i,0\)
  root -l -q -b DataMCComp7.C\(\"w_delta_phi_ee\",1,$i,0\)
  root -l -q -b DataMCComp7.C\(\"w_delta_phi_mm\",1,$i,0\)

  i=$((i+1))
done

i=1
while [ $i -le 2 ]; do

  root -l -q -b DataMCComp7.C\(\"w_first_bjet_pt\",1,$i,1\)
  root -l -q -b DataMCComp7.C\(\"w_first_bjet_eta\",1,$i,1\)
  root -l -q -b DataMCComp7.C\(\"w_pt_Z_ee\",1,$i,1\)
  root -l -q -b DataMCComp7.C\(\"w_pt_Z_mm\",1,$i,1\)
  root -l -q -b DataMCComp7.C\(\"w_Ht\",1,$i,1\)
  root -l -q -b DataMCComp7.C\(\"w_delta_phi_ee\",1,$i,1\)
  root -l -q -b DataMCComp7.C\(\"w_delta_phi_mm\",1,$i,1\)

  i=$((i+1))
done

exit
