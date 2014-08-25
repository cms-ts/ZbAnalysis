#!/bin/sh

VERSION=v14
CUT=0

if [ ! -z "$1" ]; then
  VERSION=$1
fi

if [ ! -z "$2" ]; then
  CUT=$2
fi

./submit.sh DoubleElectron_2012A_22Jan13 $VERSION $CUT
./submit.sh DoubleElectron_2012B_22Jan13 $VERSION $CUT
./submit.sh DoubleElectron_2012C_22Jan13 $VERSION $CUT
./submit.sh DoubleElectron_2012D_22Jan13 $VERSION $CUT

./submit.sh DoubleMuParked_2012A_22Jan13 $VERSION $CUT
./submit.sh DoubleMuParked_2012B_22Jan13 $VERSION $CUT
./submit.sh DoubleMuParked_2012C_22Jan13 $VERSION $CUT
./submit.sh DoubleMuParked_2012D_22Jan13 $VERSION $CUT

./submit.sh MuEG_2012A_22Jan13 $VERSION $CUT
./submit.sh MuEG_2012B_22Jan13 $VERSION $CUT
./submit.sh MuEG_2012C_22Jan13 $VERSION $CUT
./submit.sh MuEG_2012D_22Jan13 $VERSION $CUT

./submit.sh DYJetsToLL $VERSION $CUT
./submit.sh DYJetsToLL2 $VERSION $CUT
./submit.sh QCD $VERSION $CUT
./submit.sh TTbar $VERSION $CUT
./submit.sh WW $VERSION $CUT
./submit.sh WZ $VERSION $CUT
./submit.sh Wj $VERSION $CUT
./submit.sh ZZ $VERSION $CUT

exit
