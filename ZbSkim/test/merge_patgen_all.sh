#!/bin/sh

VERSION=v12

if [ ! -z "$1" ]; then
  VERSION=$1
fi

./merge.sh DYJetsToLL_patgen $VERSION

exit
