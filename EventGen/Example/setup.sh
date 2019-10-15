#!/bin/bash

setup_PYTHIA() {
    export PYTHIA8LOCATION=/home/jfilipek/Programs/pythia8235
    export PYTHIA8DATA=${PYTHIA8LOCATION}/share/Pythia8/xmldoc/
    export LD_LIBRARY_PATH=${PYTHIA8LOCATION}/lib/:$LD_LIBRARY_PATH
}

setup_ROOT() {
    source /home/jfilipek/Programs/root/bin/thisroot.sh
}

setup_fastjet() {
    export FASTJETLOCATION=/home/jfilipek/Programs/fastjet-install/
    export LD_LIBRARY_PATH=/home/jfilipek/Programs/fastjet-install/lib:$LD_LIBRARY_PATH
}

setup_boost() {
    export BOOSTINCDIR=/home/jfilipek/Programs/boost-install/include/
    export BOOSTLIBLOCATION=/home/jfilipek/Programs/boost-install/lib
    export LD_LIBRARY_PATH=${BOOSTLIBLOCATION}:$LD_LIBRARY_PATH
}

setup_ROOT
setup_PYTHIA
setup_fastjet
setup_boost
