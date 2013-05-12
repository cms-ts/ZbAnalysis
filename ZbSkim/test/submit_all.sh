#!/bin/sh

VERSION=v06

if [ ! -z "$1" ]; then
  VERSION=$1
fi

if [ "${VERSION}" \< "v07" ]; then

  ./submit.sh DoubleElectron_2012A_13Jul12 $VERSION
  ./submit.sh DoubleElectron_2012A_06Aug12 $VERSION
  ./submit.sh DoubleElectron_2012B_13Jul12 $VERSION
  ./submit.sh DoubleElectron_2012C $VERSION
  ./submit.sh DoubleElectron_2012C_11Dec12 $VERSION
  ./submit.sh DoubleElectron_2012C_24Aug12 $VERSION
  ./submit.sh DoubleElectron_2012D $VERSION
  ./submit.sh DoubleElectron_2012D_16Jan13 $VERSION

  ./submit.sh DoubleMu_2012A_13Jul12 $VERSION
  ./submit.sh DoubleMu_2012A_06Aug12 $VERSION
  ./submit.sh DoubleMu_2012B_13Jul12 $VERSION
  ./submit.sh DoubleMu_2012C $VERSION
  ./submit.sh DoubleMu_2012C_11Dec12 $VERSION
  ./submit.sh DoubleMu_2012C_24Aug12 $VERSION
  ./submit.sh DoubleMu_2012D $VERSION
  ./submit.sh DoubleMu_2012D_16Jan13 $VERSION

else

  ./submit.sh DoubleElectron_2012A_22Jan13 $VERSION
  ./submit.sh DoubleElectron_2012B_22Jan13 $VERSION
  ./submit.sh DoubleElectron_2012C_22Jan13 $VERSION
  ./submit.sh DoubleElectron_2012D_22Jan13 $VERSION

  ./submit.sh DoubleMu_2012A_22Jan13 $VERSION
  ./submit.sh DoubleMuParked_2012B_22Jan13 $VERSION
  ./submit.sh DoubleMuParked_2012C_22Jan13 $VERSION
  ./submit.sh DoubleMuParked_2012D_22Jan13 $VERSION

fi

./submit.sh DYJetsToLL $VERSION
./submit.sh QCD $VERSION
./submit.sh TTbar $VERSION
./submit.sh WW $VERSION
./submit.sh WZ $VERSION
./submit.sh Wj $VERSION
./submit.sh ZZ $VERSION

exit
