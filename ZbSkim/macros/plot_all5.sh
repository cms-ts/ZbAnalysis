#!/bin/sh

export ROOT_HIST=0

cd $CMS_PATH/slc5_amd64_gcc472/cms/cmssw/CMSSW_6_1_0
eval `scramv1 runtime -sh`
cd -

i=1
while [ $i -le 2 ]; do

root -l -q -b DataMCComp5.C\(\"w_mass_ee_wide\",1,$i\)
root -l -q -b DataMCComp5.C\(\"w_mass_mm_wide\",1,$i\)
root -l -q -b DataMCComp5.C\(\"w_mass_ee_b_wide\",1,$i\)
root -l -q -b DataMCComp5.C\(\"w_mass_mm_b_wide\",1,$i\)

i=$((i+1))
done

exit
