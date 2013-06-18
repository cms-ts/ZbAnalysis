#!/bin/sh

ulimit -c 0

cd /gpfs/cms/users/candelis/CMSSW_5_3_9
eval `scramv1 runtime -sh`
cd -

#cp /gpfs/cms/users/candelis/work/ZbSkim/test/demoanalyzer_cfg.py job.py
cp /gpfs/cms/users/candelis/work/ZbSkim/test/demogenanalyzer_cfg.py job.py

pileup=$1
echo "process.demo_ee.pileup = cms.untracked.string('"$pileup"')" >> job.py
echo "process.demo_ee_up.pileup = cms.untracked.string('"$pileup"')" >> job.py
echo "process.demo_ee_down.pileup = cms.untracked.string('"$pileup"')" >> job.py
echo "process.demo_mm.pileup = cms.untracked.string('"$pileup"')" >> job.py
echo "process.demo_mm_up.pileup = cms.untracked.string('"$pileup"')" >> job.py
echo "process.demo_mm_down.pileup = cms.untracked.string('"$pileup"')" >> job.py
echo "process.demo_ee_btag.pileup = cms.untracked.string('"$pileup"')" >> job.py
echo "process.demo_mm_btag.pileup = cms.untracked.string('"$pileup"')" >> job.py
echo "process.demo2_ee.pileup = cms.untracked.string('"$pileup"')" >> job.py
echo "process.demo2_mm.pileup = cms.untracked.string('"$pileup"')" >> job.py

shift

echo "fileList = cms.untracked.vstring()" >> job.py
i=1
while [ $i -le $# ]; do
  file=`echo ${!i} | sed -e 's;/gpfs/grid/srm/cms;;'`
  echo "fileList.extend(['"$file"'])" >> job.py
  file=`basename ${!i} | sed -e 's/patTuple/rootTuple/'`
  echo "process.TFileService.fileName = cms.string('"$file"')" >> job.py
  i=$((i+1))
done
echo "process.source.fileNames = fileList" >> job.py

(time cmsRun -j job.xml job.py) > job.log 2>&1

exit
