#!/bin/bash

setup_PYTHIA() {
    export PYTHIA8LOCATION=/home/balbok/Programs/pythia8235/
    export PYTHIA8DATA=${PYTHIA8LOCATION}xmldoc/
    export LD_LIBRARY_PATH=${PYTHIA8LOCATION}lib/:$LD_LIBRARY_PATH
}

setup_ROOT() {
    source /home/balbok/Programs/root/bin/thisroot.sh
}

setup_fastjet() {
    export FASTJETLOCATION=/usr/local/
    export LD_LIBRARY_PATH=${FASTJETLOCATION}lib/:$LD_LIBRARY_PATH
}

setup_boost() {
    export BOOSTROOTLOCATION=/home/balbok/Programs/boost_1_68_0/
    export BOOSTINCDIR=${BOOSTROOTLOCATION}include/
    export BOOSTLIBLOCATION=${BOOSTROOTLOCATION}lib/
    export LD_LIBRARY_PATH=${BOOSTLIBLOCATION}:$LD_LIBRARY_PATH
}

setup_ROOT
setup_PYTHIA
setup_fastjet
setup_boost

