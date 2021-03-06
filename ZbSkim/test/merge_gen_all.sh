#!/bin/sh

VERSION=v14

if [ ! -z "$1" ]; then
  VERSION=$1
fi

./merge.sh DYJetsToLL_gen $VERSION
./merge.sh DYJetsToLL_gen_weights $VERSION

./merge.sh DYJetsToLL2_gen $VERSION

./merge.sh DYToEE_powheg_gen $VERSION
./merge.sh DYToMuMu_powheg_gen $VERSION

./merge.sh DYJets_sherpa_gen $VERSION

exit
