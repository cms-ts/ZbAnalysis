#!/bin/sh

VERSION=v09

if [ ! -z "$1" ]; then
  VERSION=$1
fi

if [ "${VERSION}" \< "v10" ]; then

  ./submit.sh DoubleElectron_2012A_22Jan13 $VERSION
  ./submit.sh DoubleElectron_2012B_22Jan13 $VERSION
  ./submit.sh DoubleElectron_2012C_22Jan13 $VERSION
  ./submit.sh DoubleElectron_2012D_22Jan13 $VERSION

  ./submit.sh DoubleMu_2012A_22Jan13 $VERSION
  ./submit.sh DoubleMuParked_2012B_22Jan13 $VERSION
  ./submit.sh DoubleMuParked_2012C_22Jan13 $VERSION
  ./submit.sh DoubleMuParked_2012D_22Jan13 $VERSION

  ./submit.sh MuEG_2012A_22Jan13 $VERSION
  ./submit.sh MuEG_2012B_22Jan13 $VERSION
  ./submit.sh MuEG_2012C_22Jan13 $VERSION
  ./submit.sh MuEG_2012D_22Jan13 $VERSION

else

  ./submit.sh DoubleElectron_2012A_22Jan13 $VERSION
  ./submit.sh DoubleElectron_2012B_22Jan13 $VERSION
  ./submit.sh DoubleElectron_2012C_22Jan13 $VERSION
  ./submit.sh DoubleElectron_2012D_22Jan13 $VERSION

  ./submit.sh DoubleMuParked_2012A_22Jan13 $VERSION
  ./submit.sh DoubleMuParked_2012B_22Jan13 $VERSION
  ./submit.sh DoubleMuParked_2012C_22Jan13 $VERSION
  ./submit.sh DoubleMuParked_2012D_22Jan13 $VERSION

  ./submit.sh MuEG_2012A_22Jan13 $VERSION
  ./submit.sh MuEG_2012B_22Jan13 $VERSION
  ./submit.sh MuEG_2012C_22Jan13 $VERSION
  ./submit.sh MuEG_2012D_22Jan13 $VERSION

fi

./submit.sh DYJetsToLL $VERSION
./submit.sh QCD $VERSION
./submit.sh TTbar $VERSION
./submit.sh WW $VERSION
./submit.sh WZ $VERSION
./submit.sh Wj $VERSION
./submit.sh ZZ $VERSION

exit
