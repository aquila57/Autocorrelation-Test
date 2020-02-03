#!/bin/bash
make -f libcorr.mak
make -f corr.mak
make -f etauscorr.mak
make -f rucorr.mak
make -f lfsrcorr.mak
make -f fibocorr.mak
make -f sinecorr.mak
make -f gslcorr.mak
make -f diecorr.mak
