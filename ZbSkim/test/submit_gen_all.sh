#!/bin/sh

VERSION=v09

if [ ! -z "$1" ]; then
  VERSION=$1
fi

./submit_gen.sh DYJetsToLL_gen $VERSION

./submit_gen.sh DYToEE_powheg_gen $VERSION
./submit_gen.sh DYToMuMu_powheg_gen $VERSION

#./submit_gen.sh DYJetsToLL_sherpa_gen $VERSION
./submit_gen.sh DYJets_sherpa_gen $VERSION

exit