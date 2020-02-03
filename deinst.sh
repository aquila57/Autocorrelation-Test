#!/bin/bash
make -f libcorr.mak clean
make -f corr.mak clean
make -f etauscorr.mak clean
make -f rucorr.mak clean
make -f lfsrcorr.mak clean
make -f fibocorr.mak clean
make -f sinecorr.mak clean
make -f gslcorr.mak clean
make -f diecorr.mak clean
rm -f alldie.out allgsl.out
