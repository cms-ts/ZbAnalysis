#!/bin/sh

VERSION=v12

if [ ! -z "$1" ]; then
  VERSION=$1
fi

./submit_patgen.sh DYJetsToLL $VERSION

exit
