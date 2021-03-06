#!/bin/sh

VERSION=v14

if [ ! -z "$1" ]; then
  VERSION=$1
fi

./submit_gen.sh DYJetsToLL_gen $VERSION
./submit_gen.sh DYJetsToLL_gen_weights $VERSION

./submit_gen.sh DYJetsToLL2_gen $VERSION

./submit_gen.sh DYToEE_powheg_gen $VERSION
./submit_gen.sh DYToMuMu_powheg_gen $VERSION

./submit_gen.sh DYJetsToLL_aMC_gen $VERSION
./submit_gen.sh DYJetsToLL_aMC_weights $VERSION

./submit_gen.sh DYJets_sherpa_gen $VERSION

exit
