#!/bin/sh

VERSION=v01
USER=$USER

do_check=1
if [ ! -z "$1" ]; then
  if [ "$1" == "-n" ]; then
    do_check=0
    shift
  fi
fi

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

DIR=/gpfs/grid/srm/cms/store/user/$USER/grid/$VERSION/

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

  N1=`ls $DIR/$D/ | grep root | awk -F_ '{print $2}' | sort | wc -w`
  N2=`ls $DIR/$D/ | grep root | awk -F_ '{print $2}' | sort | uniq | wc -w`

  NT=$((NT+N2))

  SIZE=""
  if [ $do_size -eq 1 ]; then
    SIZE=`du -sh $DIR/$D/ | awk '{print $1}'`
  fi

  printf "%30s\t%s%6i%10s\n" $D ":" $N2 $SIZE

  if [ $do_check -eq 1 ]; then

    if [ "$N1" -gt "$N2" ]; then
      echo "ERROR: possible duplicate in "$D" : "$N1" != "$N2
      ls $DIR/$D/ | grep root | awk -F_ '{print $2}' | sort > /tmp/l1
      ls $DIR/$D/ | grep root | awk -F_ '{print $2}' | sort | uniq > /tmp/l2
      ERRORS=`diff /tmp/l1 /tmp/l2 | grep -v d | sed -e 's/<//' | sed -e 's/ //g'`
      echo $ERRORS
      for E in $ERRORS; do
        ls -lt $DIR/$D/*Tuple_${E}_*_*.root
        ls -lt $DIR/$D/*Tuple_${E}_*_*.root | tail -1 | awk '{print $9}' | sed -e 's;/gpfs/grid/srm;lcg-del -l srm://gridsrm.ts.infn.it;'
      done
      rm /tmp/l1 /tmp/l2
    fi

  FILES=`ls $DIR/$D/ | grep root`

  for F in $FILES; do

    if [ -e $DIR/$D/$F ]; then
      if [ ! -s $DIR/$D/$F ]; then
        echo "ERROR: empty file "
        ls -lt $DIR/$D/$F
      fi
    fi

  done

  fi

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
