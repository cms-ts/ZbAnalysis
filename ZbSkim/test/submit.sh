#!/bin/sh

#QUEUE=ts_cms
#QUEUE=short
#QUEUE=normal
QUEUE=normal_io

VERSION=v10

DATADIR=/gpfs/grid/srm/cms/store/user/vieri/grid

WORKDIR=/gpfs/cms/users/candelis/work/ZbSkim/test

OUTDIR=/gpfs/cms/users/candelis/work/ZbSkim/test/data

if [ $# -eq 0 ]; then
  echo 'Usage: submit.sh jobdir [version]'
  exit
fi

JOBDIR=$1

if [ ! -z "$2" ]; then
  VERSION=$2
fi

if [ ! -e $DATADIR/$VERSION/$JOBDIR ]; then
  echo "ERROR: $DATADIR/$VERSION/$JOBDIR does not exist !"
  exit
fi

P=S10
if [ "$JOBDIR" == "QCD" ]; then
  P=S7
fi

rm -fr $OUTDIR/$VERSION/$JOBDIR
mkdir -p $OUTDIR/$VERSION/$JOBDIR
cd $OUTDIR/$VERSION/$JOBDIR

find $DATADIR/$VERSION/$JOBDIR -maxdepth 1 -name '*.root' | \
xargs -n 50 bsub -q $QUEUE -e /dev/null -o /dev/null $WORKDIR/job.sh $P

exit
