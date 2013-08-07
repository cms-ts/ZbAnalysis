#!/bin/sh

VERSION=v09

if [ ! -z "$1" ]; then
  VERSION=$1
fi

if [ "${VERSION}" \< "v10" ]; then

  ./merge.sh DoubleElectron_2012A_22Jan13 $VERSION
  ./merge.sh DoubleElectron_2012B_22Jan13 $VERSION
  ./merge.sh DoubleElectron_2012C_22Jan13 $VERSION
  ./merge.sh DoubleElectron_2012D_22Jan13 $VERSION

  ./merge.sh DoubleMu_2012A_22Jan13 $VERSION
  ./merge.sh DoubleMuParked_2012B_22Jan13 $VERSION
  ./merge.sh DoubleMuParked_2012C_22Jan13 $VERSION
  ./merge.sh DoubleMuParked_2012D_22Jan13 $VERSION

  ./merge.sh MuEG_2012A_22Jan13 $VERSION
  ./merge.sh MuEG_2012B_22Jan13 $VERSION
  ./merge.sh MuEG_2012C_22Jan13 $VERSION
  ./merge.sh MuEG_2012D_22Jan13 $VERSION

else

  ./merge.sh DoubleElectron_2012A_22Jan13 $VERSION
  ./merge.sh DoubleElectron_2012B_22Jan13 $VERSION
  ./merge.sh DoubleElectron_2012C_22Jan13 $VERSION
  ./merge.sh DoubleElectron_2012D_22Jan13 $VERSION

  ./merge.sh DoubleMuParked_2012A_22Jan13 $VERSION
  ./merge.sh DoubleMuParked_2012B_22Jan13 $VERSION
  ./merge.sh DoubleMuParked_2012C_22Jan13 $VERSION
  ./merge.sh DoubleMuParked_2012D_22Jan13 $VERSION

  ./merge.sh MuEG_2012A_22Jan13 $VERSION
  ./merge.sh MuEG_2012B_22Jan13 $VERSION
  ./merge.sh MuEG_2012C_22Jan13 $VERSION
  ./merge.sh MuEG_2012D_22Jan13 $VERSION

fi

./merge.sh data-all $VERSION

./merge.sh DYJetsToLL $VERSION
./merge.sh QCD $VERSION
./merge.sh TTbar $VERSION
./merge.sh ZZ $VERSION
./merge.sh WZ $VERSION
./merge.sh Wj $VERSION
./merge.sh WW $VERSION

exit
