#  libcorr.mak - Compile Subroutine Library Version 0.1.0
#  Copyright (C) 2020 aquila57 at github.com

#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License as
#  published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.

#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU General Public License for more details.

#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to:

   #  Free Software Foundation, Inc.
   #  59 Temple Place - Suite 330
   #  Boston, MA 02111-1307, USA.

OBJ=diepfx.o \
	getdie.o \
	bld_maxint.o \
	bld_rngtbl.o \
	eegl.o \
	eeglabt.o \
	eeglbit.o \
	eeglcrc.o \
	eegldspl.o \
	eeglfrac.o \
	eeglinit.o \
	eeglint.o \
	eeglpwr.o \
	eeglsd.o \
	eeglstrt.o \
	eeglunif.o \
	etaus.o \
	etausabt.o \
	etausbit.o \
	etausdspl.o \
	etausfrac.o \
	etausinit.o \
	etausint.o \
	etauspwr.o \
	etausstrt.o \
	etausunif.o

CC=gcc

AR=ar

CFLAGS=-c -Wall -O2

LDFLAGS=

libcorr.a:			$(OBJ)
		rm -f libcorr.a
		$(AR) r libcorr.a $(OBJ)

diepfx.o:			diepfx.c
		$(CC) $(CFLAGS) diepfx.c

getdie.o:			getdie.c
		$(CC) $(CFLAGS) getdie.c

bld_maxint.o:			bld_maxint.c
		$(CC) $(CFLAGS) bld_maxint.c

bld_rngtbl.o:			bld_rngtbl.c
		$(CC) $(CFLAGS) bld_rngtbl.c

eegl.o:				eegl.c
		$(CC) $(CFLAGS) eegl.c

eeglabt.o:			eeglabt.c
		$(CC) $(CFLAGS) eeglabt.c

eeglbit.o:			eeglbit.c
		$(CC) $(CFLAGS) eeglbit.c

eeglcrc.o:			eeglcrc.c
		$(CC) $(CFLAGS) eeglcrc.c

eegldspl.o:			eegldspl.c
		$(CC) $(CFLAGS) eegldspl.c

eeglfrac.o:			eeglfrac.c
		$(CC) $(CFLAGS) eeglfrac.c

eeglinit.o:			eeglinit.c
		$(CC) $(CFLAGS) eeglinit.c

eeglint.o:			eeglint.c
		$(CC) $(CFLAGS) eeglint.c

eeglpwr.o:			eeglpwr.c
		$(CC) $(CFLAGS) eeglpwr.c

eeglsd.o:			eeglsd.c
		$(CC) $(CFLAGS) eeglsd.c

eeglstrt.o:			eeglstrt.c
		$(CC) $(CFLAGS) eeglstrt.c

eeglunif.o:			eeglunif.c
		$(CC) $(CFLAGS) eeglunif.c

etaus.o:			etaus.c
		$(CC) $(CFLAGS) etaus.c

etausabt.o:			etausabt.c
		$(CC) $(CFLAGS) etausabt.c

etausbit.o:			etausbit.c
		$(CC) $(CFLAGS) etausbit.c

etauscrc.o:			etauscrc.c
		$(CC) $(CFLAGS) etauscrc.c

etausdspl.o:			etausdspl.c
		$(CC) $(CFLAGS) etausdspl.c

etausfrac.o:			etausfrac.c
		$(CC) $(CFLAGS) etausfrac.c

etausinit.o:			etausinit.c
		$(CC) $(CFLAGS) etausinit.c

etausint.o:			etausint.c
		$(CC) $(CFLAGS) etausint.c

etauspwr.o:			etauspwr.c
		$(CC) $(CFLAGS) etauspwr.c

etaussd.o:			etaussd.c
		$(CC) $(CFLAGS) etaussd.c

etausstrt.o:			etausstrt.c
		$(CC) $(CFLAGS) etausstrt.c

etausunif.o:			etausunif.c
		$(CC) $(CFLAGS) etausunif.c

clean:
		rm -f $(OBJ) libcorr.a
