#!/bin/sh

VERSION=v14

if [ ! -z "$1" ]; then
  VERSION=$1
fi

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

./merge.sh data-all $VERSION

./merge.sh DYJetsToLL $VERSION
./merge.sh DYJetsToLL2 $VERSION

./merge.sh DYJetsToLL_aMC $VERSION

./merge.sh DYJets_sherpa $VERSION

./merge.sh QCD $VERSION
./merge.sh TTbar $VERSION
./merge.sh ZZ $VERSION
./merge.sh WZ $VERSION
./merge.sh Wj $VERSION
./merge.sh WW $VERSION

./merge.sh TBar_s $VERSION
./merge.sh TBar_t $VERSION
./merge.sh TBar_tW $VERSION
./merge.sh T_s $VERSION
./merge.sh T_t $VERSION
./merge.sh T_tW $VERSION

exit
