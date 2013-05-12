#!/bin/sh

VERSION=v06

if [ $# -eq 0 ]; then
  echo 'Usage: merge.sh jobdir [version]'
  exit
fi

JOBDIR=$1

if [ ! -z "$2" ]; then
  VERSION=$2
fi

WORKDIR=/gpfs/cms/users/candelis/work/ZbSkim/test/data

if [ ! -e $WORKDIR/$VERSION ]; then
  echo 'ERROR: version "'$VERSION'" does not exist !'
  exit
fi

cd $CMS_PATH/slc5_amd64_gcc472/cms/cmssw/CMSSW_6_1_0
eval `scramv1 runtime -sh`
cd $WORKDIR

opts="-T -v 0"

if [ -d $WORKDIR/$VERSION/$JOBDIR ]; then
  echo 'Preparing '$WORKDIR/$VERSION/$JOBDIR.root
  rm -f $WORKDIR/$VERSION/$JOBDIR.root
  hadd $opts $WORKDIR/$VERSION/$JOBDIR.root $WORKDIR/$VERSION/$JOBDIR/LSFJOB_*/rootTuple_*.root
elif [ "$JOBDIR" == "data-all" ]; then
  echo 'Preparing '$WORKDIR/$VERSION/DoubleElectron_2012_merge.root
  rm -f $WORKDIR/$VERSION/DoubleElectron_2012_merge.root
  hadd $opts $WORKDIR/$VERSION/DoubleElectron_2012_merge.root $WORKDIR/$VERSION/DoubleElectron_2012*.root
  echo 'Preparing '$WORKDIR/$VERSION/DoubleMu_2012_merge.root
  rm -f $WORKDIR/$VERSION/DoubleMu_2012_merge.root
  hadd $opts $WORKDIR/$VERSION/DoubleMu_2012_merge.root $WORKDIR/$VERSION/DoubleMu*_2012*.root
else
  echo 'ERROR: jobdir "'$JOBDIR'" does not exist !'
fi

exit
