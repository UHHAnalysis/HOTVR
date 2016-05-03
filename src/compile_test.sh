#!/bin/bash

export HEPMCDIR=$HOME/HepMC
export PYTHIA8DIR=$HOME/pythia8
export FASTJETDIR=$HOME/fastjet

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:.libs:$HEPMCDIR/lib:$PYTHIA8DIR/lib/archive:$FASTJETDIR/lib

export PYTHIA8DATA=$PYTHIA8DIR/xmldoc

g++ -o TestAnalysis TestAnalysis.cxx ParseUtils.cxx -I$PYTHIA8DIR/include -I$HEPMCDIR/include `root-config --cflags --libs` `fastjet-config --cxxflags --libs` -L$PYTHIA8DIR/lib/archive -lpythia8 -llhapdfdummy -lhepmcinterface -L$HEPMCDIR/lib -lHepMC -L.libs -lDeconstruction


