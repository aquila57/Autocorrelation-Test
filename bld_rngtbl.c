/* bld_rngtbl.c - Build dieharder RNG table  Version 0.1.0 */
/* Copyright (C) 2020 aquila57 at github.com */

/* This program is free software; you can redistribute it and/or     */
/* modify it under the terms of the GNU General Public License as    */
/* published by the Free Software Foundation; either version 2 of    */
/* the License, or (at your option) any later version.               */

/* This program is distributed in the hope that it will be useful,   */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of    */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */
/* GNU General Public License for more details.                      */

/* You should have received a copy of the GNU General Public License */
/* along with this program; if not, write to:                        */

   /* Free Software Foundation, Inc.                                 */
   /* 59 Temple Place - Suite 330                                    */
   /* Boston, MA 02111-1307, USA.                                    */

#include "corr.h"

/* Build dieharder random number generator table */

void bld_rngtbl(xxfmt *xx)
   {
   int i;
   for (i=0;i<1000;i++)
      {
      int j;
      for (j=0;j<64;j++) xx->rngtbl[i][j] = '\0';
      } /* for (i=0;i<1000;i++) */
   strcpy(xx->rngtbl[0],"borosh13");
   strcpy(xx->rngtbl[1],"cmrg");
   strcpy(xx->rngtbl[2],"coveyou");
   strcpy(xx->rngtbl[3],"fishman18");
   strcpy(xx->rngtbl[4],"fishman20");
   strcpy(xx->rngtbl[5],"fishman2x");
   strcpy(xx->rngtbl[6],"gfsr4");
   strcpy(xx->rngtbl[7],"knuthran");
   strcpy(xx->rngtbl[8],"knuthran2");
   strcpy(xx->rngtbl[9],"knuthran2002");
   strcpy(xx->rngtbl[10],"lecuyer21");
   strcpy(xx->rngtbl[11],"minstd");
   strcpy(xx->rngtbl[12],"mrg");
   strcpy(xx->rngtbl[13],"mt19937");
   strcpy(xx->rngtbl[14],"mt19937_1999");
   strcpy(xx->rngtbl[15],"mt19937_1998");
   strcpy(xx->rngtbl[16],"r250");
   strcpy(xx->rngtbl[17],"ran0");
   strcpy(xx->rngtbl[18],"ran1");
   strcpy(xx->rngtbl[19],"ran2");
   strcpy(xx->rngtbl[20],"ran3");
   strcpy(xx->rngtbl[21],"rand");
   strcpy(xx->rngtbl[22],"rand48");
   strcpy(xx->rngtbl[23],"random128-bsd");
   strcpy(xx->rngtbl[24],"random128-glibc2");
   strcpy(xx->rngtbl[25],"random128-libc5");
   strcpy(xx->rngtbl[26],"random256-bsd");
   strcpy(xx->rngtbl[27],"random256-glibc2");
   strcpy(xx->rngtbl[28],"random256-libc5");
   strcpy(xx->rngtbl[29],"random32-bsd");
   strcpy(xx->rngtbl[30],"random32-glibc2");
   strcpy(xx->rngtbl[31],"random32-libc5");
   strcpy(xx->rngtbl[32],"random64-bsd");
   strcpy(xx->rngtbl[33],"random64-glibc2");
   strcpy(xx->rngtbl[34],"random64-libc5");
   strcpy(xx->rngtbl[35],"random8-bsd");
   strcpy(xx->rngtbl[36],"random8-glibc2");
   strcpy(xx->rngtbl[37],"random8-libc5");
   strcpy(xx->rngtbl[38],"random-bsd");
   strcpy(xx->rngtbl[39],"random-glibc2");
   strcpy(xx->rngtbl[40],"random-libc5");
   strcpy(xx->rngtbl[41],"randu");
   strcpy(xx->rngtbl[42],"ranf");
   strcpy(xx->rngtbl[43],"ranlux");
   strcpy(xx->rngtbl[44],"ranlux389");
   strcpy(xx->rngtbl[45],"ranlxd1");
   strcpy(xx->rngtbl[46],"ranlxd2");
   strcpy(xx->rngtbl[47],"ranlxs0");
   strcpy(xx->rngtbl[48],"ranlxs1");
   strcpy(xx->rngtbl[49],"ranlxs2");
   strcpy(xx->rngtbl[50],"ranmar");
   strcpy(xx->rngtbl[51],"slatec");
   strcpy(xx->rngtbl[52],"taus");
   strcpy(xx->rngtbl[53],"taus2");
   strcpy(xx->rngtbl[54],"taus113");
   strcpy(xx->rngtbl[55],"transputer");
   strcpy(xx->rngtbl[56],"tt800");
   strcpy(xx->rngtbl[57],"uni");
   strcpy(xx->rngtbl[58],"uni32");
   strcpy(xx->rngtbl[59],"vax");
   strcpy(xx->rngtbl[200],"stdin_input_raw");
   strcpy(xx->rngtbl[201],"file_input_raw");
   strcpy(xx->rngtbl[202],"file_input");
   strcpy(xx->rngtbl[203],"ca");
   strcpy(xx->rngtbl[204],"uvag");
   strcpy(xx->rngtbl[205],"AES_OFB");
   strcpy(xx->rngtbl[206],"Threefish_OFB");
   strcpy(xx->rngtbl[207],"XOR");
   strcpy(xx->rngtbl[208],"kiss");
   strcpy(xx->rngtbl[400],"R_wichmann_hill");
   strcpy(xx->rngtbl[401],"R_marsaglia_multic.");
   strcpy(xx->rngtbl[402],"R_super_duper");
   strcpy(xx->rngtbl[403],"R_mersenne_twister");
   strcpy(xx->rngtbl[404],"R_knuth_taocp");
   strcpy(xx->rngtbl[405],"R_knuth_taocp2");
   } /* bld_rngtbl */
