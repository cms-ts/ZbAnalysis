#!/bin/sh

./merge.sh DoubleElectron_2012A_13Jul12 $1
./merge.sh DoubleElectron_2012A_06Aug12 $1
./merge.sh DoubleElectron_2012B_13Jul12 $1
./merge.sh DoubleElectron_2012C $1
./merge.sh DoubleElectron_2012C_11Dec12 $1
./merge.sh DoubleElectron_2012C_24Aug12 $1
./merge.sh DoubleElectron_2012D $1
./merge.sh DoubleElectron_2012D_16Jan13 $1

./merge.sh DoubleMu_2012A_13Jul12 $1
./merge.sh DoubleMu_2012A_06Aug12 $1
./merge.sh DoubleMu_2012B_13Jul12 $1
./merge.sh DoubleMu_2012C $1
./merge.sh DoubleMu_2012C_11Dec12 $1
./merge.sh DoubleMu_2012C_24Aug12 $1
./merge.sh DoubleMu_2012D $1
./merge.sh DoubleMu_2012D_16Jan13 $1

./merge.sh data-all $1

./merge.sh DYJetsToLL $1
./merge.sh QCD $1
./merge.sh TTbar $1
./merge.sh ZZ $1
./merge.sh WZ $1
./merge.sh Wj $1
./merge.sh WW $1

