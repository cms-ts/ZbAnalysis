#!/bin/sh

#QUEUE=ts_cms
#QUEUE=short
#QUEUE=normal
QUEUE=normal_io

VERSION=v11
CUT=0

DATADIR=/gpfs/grid/srm/cms/store/user/vieri/grid

WORKDIR=/gpfs/cms/users/candelis/work/ZbSkim/test

OUTDIR=/gpfs/cms/users/candelis/work/ZbSkim/test/data

if [ $# -eq 0 ]; then
  echo 'Usage: submit.sh jobdir [version] [cut]'
  exit
fi

JOBDIR=$1

if [ ! -z "$2" ]; then
  VERSION=$2
fi

if [ ! -z "$3" ]; then
  CUT=$3
fi

if [ ! -e $DATADIR/$VERSION/$JOBDIR ]; then
  echo "ERROR: $DATADIR/$VERSION/$JOBDIR does not exist !"
  exit
fi

P=S10
if [ "$JOBDIR" == "QCD" ]; then
  P=S7
fi

if [ "$CUT" == "0" ]; then
  rm -fr $OUTDIR/$VERSION/$JOBDIR
  mkdir -p $OUTDIR/$VERSION/$JOBDIR
  cd $OUTDIR/$VERSION/$JOBDIR
else
  rm -fr $OUTDIR/$VERSION.$CUT/$JOBDIR
  mkdir -p $OUTDIR/$VERSION.$CUT/$JOBDIR
  cd $OUTDIR/$VERSION.$CUT/$JOBDIR
fi

find $DATADIR/$VERSION/$JOBDIR -maxdepth 1 -name '*.root' | \
xargs -n 50 bsub -q $QUEUE -e /dev/null -o /dev/null $WORKDIR/job.sh $P $CUT

exit
