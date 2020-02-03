/* corr.h - Autocorrelation Header File  Version 0.1.0 */
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

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <sys/times.h>
#include <math.h>
#include <errno.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_rng.h>
#include "eegl.h"
#include "etaus.h"

#define SIZE 1000000

#define EOFDIE (-999999999.0)

typedef struct xxstruct {
   char *rngname;
   char dieharder_rngname[128];
   char rngtbl[1000][64];
   int generator;
   int eofsw;
   unsigned int maxint;
   double dblsmpls;
   double count;
   double modulus;
   double maxint_tbl[2048];
   } xxfmt;

void calcchi(int *freqtbl);
void bld_rngtbl(xxfmt *xx);
void bld_maxint(xxfmt *xx);
void diepfx(xxfmt *xx);
double getdie(xxfmt *xx);
