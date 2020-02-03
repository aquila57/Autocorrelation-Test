/* corr.c - Autocorrelation, eegl64 generator  Version 0.1.0 */
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
#include <math.h>
#include <gsl/gsl_sf_erf.h>
#include "eegl.h"

#define SIZE 1000000

void putstx(char *pgm)
   {
   fprintf(stderr,"Usage: %s index lag\n", pgm);
   fprintf(stderr,"Where index is i, the starting offset\n");
   fprintf(stderr,"Where lag is m, the distance between samples\n");
   fprintf(stderr,"Example: %s 3 5\n", pgm);
   exit(1);
   } /* putstx */

int main(int argc, char **argv)
   {
   int i;
   int indx;
   int lag;
   int lmt;
   double M;
   double N;
   double dblindx;
   double dbllag;
   double sum;
   double rho;
   double sigma;
   double zedzero; 
   double abszedzero; 
   double pvalue; 
   double *p,*q;
   double *smpls;
   eefmt *ee;
   if (argc != 3) putstx(*argv);
   indx = atoi(*(argv+1));
   if (indx < 0 || indx > 1000)
      {
      fprintf(stderr,"main: parameter indx %s "
         "out of range\n", *(argv+1));
      putstx(*argv);
      } /* indx out of range */
   lag  = atoi(*(argv+2));
   if (lag < 1 || lag > 100)
      {
      fprintf(stderr,"main: parameter lag %s "
         "out of range\n", *(argv+2));
      putstx(*argv);
      } /* lag out of range */
   smpls = (double *) malloc(sizeof(double) * SIZE + 10);
   if (smpls == NULL)
      {
      fprintf(stderr,"main: out of memory "
         "allocating smpls\n");
      exit(1);
      } /* out of memory */
   ee = (eefmt *) eeglinit();
   N = (double) SIZE;
   dblindx = (double) indx;
   dbllag  = (double) lag;
   p = (double *) smpls;
   q = (double *) smpls + SIZE;
   while (p < q)
      {
      *p++ = eeglunif(ee);
      } /* for each sample */
   free(ee->state);
   free(ee);
   M = floor(((N - dblindx) / dbllag) - 1.0);
   sum = 0.0;
   lmt = SIZE - lag;
   for (i = indx+lag; i < lmt; i += lag)
      {
      sum += (smpls[i] * smpls[i+1]);
      } /* for each lag */
   rho = (sum / (M + 1.0)) - 0.25;
   sigma = sqrt((13.0 * M) + 7.0) / (12.0 * (M + 1.0));
   zedzero = rho / sigma;
   printf("Autocorrelation Test\n");
   printf("  eegl64 Generator\n");
   printf("N %.0f index %.0f lag %.0f M %.0f\n",
      N, dblindx, dbllag, M);
   abszedzero = fabs(zedzero);
   pvalue = gsl_sf_erf_Q(abszedzero);
   printf("Z0 %f  P-value %f\n", zedzero, pvalue + pvalue);
   if (abszedzero < 0.01 || abszedzero > 1.96)
      printf("*************************************"   
         "*************************************\n");
   free(smpls);
   return(0);
   } /* main */
