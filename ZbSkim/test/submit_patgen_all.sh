#!/bin/sh

VERSION=v11

if [ ! -z "$1" ]; then
  VERSION=$1
fi

./submit_patgen.sh DYJetsToLL $VERSION

exit
