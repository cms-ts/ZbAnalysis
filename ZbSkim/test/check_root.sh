#!/bin/sh

VERSION=v01
USER=$USER

do_size=0
if [ ! -z "$1" ]; then
  if [ "$1" == "-s" ]; then
    do_size=1
    shift
  fi
fi

if [ ! -z "$1" ]; then
  VERSION=$1
fi

if [ ! -z "$2" ]; then
  USER=$2
fi

SUBDIR=""
[ "$USER" == "dellaric" ] && SUBDIR="GDR"

DIR=/gpfs/cms/users/candelis/work/ZbSkim/test/$SUBDIR/data/$VERSION/

[ "$USER" == "lalicata" ] && DIR = /gpfs/cms/users/lalicata/CMSSW_5_3_9/src/ZbAnalysis/ZbSkim/test/data/$VERSION

if [ ! -e $DIR ]; then
  echo "ERROR: $DIR does not exist !"
  exit
fi

if [ ! -z "$3" ]; then
  DIRS=$3
else
  DIRS=`ls $DIR | grep -v root`
fi

NT=0

for D in $DIRS; do

  N=`find $DIR/$D/ -name '*.root' | wc -w`

  NT=$((NT+N))

  SIZE=""
  if [ $do_size -eq 1 ]; then
    SIZE=`du -sh $DIR/$D/ | awk '{print $1}'`
  fi

  printf "%30s\t%s%6i%10s\n" $D ":" $N $SIZE

done

SIZE=""
if [ $do_size -eq 1 ]; then
  if [ ! -z "$3" ]; then
    SIZE=`du -sh $DIR/$3 | awk '{print $1}'`
  else
    SIZE=`du -sh $DIR | awk '{print $1}'`
  fi
fi

printf "%30s\t%s%6i%10s\n" "TOTAL" ":" $NT $SIZE

exit
