/* gslcorr.c - Autocorrelation, GSL generators  Version 0.1.0 */
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

/* Autocorrelation test for the GNU Scientific Library */
/* of random number generators */
/* Test with the script, tstgsl.sh */

#include "corr.h"

void putstx(char *pgm)
   {
   fprintf(stderr,"Autocorrelation test, GNU Scientific "
      "Library generators\n");
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
   gsl_rng *rng;
   gsl_rng_type *typ;
   unsigned int dttk;          /* combined date and #ticks */
   time_t now;                 /* current date and time */
   clock_t clk;                /* current number of ticks */
   struct tms t;               /* structure used by times() */
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
   /*************************************************************/
   /* obtain parameters                                         */
   /*************************************************************/
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
   /*************************************************************/
   smpls = (double *) malloc(sizeof(double) * SIZE + 10);
   if (smpls == NULL)
      {
      fprintf(stderr,"main: out of memory "
         "allocating smpls\n");
      exit(1);
      } /* out of memory */
   /***************************************************/
   /* Initialize the seed to date/time/ticks          */
   /***************************************************/
   /* get clock ticks since boot                       */
   clk = times(&t);
   /* get date & time                                  */
   time(&now);
   /* combine date, time, and ticks into a single UINT */
   dttk = (unsigned int) (now ^ clk);
   /****************************************************/
   /* Initialize the random number generator           */
   /* in the GNU Scientific Library                    */
   /****************************************************/
   gsl_rng_env_setup();
   typ = (gsl_rng_type *) gsl_rng_default;
   rng = (gsl_rng *) gsl_rng_alloc(typ);
   /* allocate GSL random number generator to set initial */
   /* values in the three fibonacci numbers               */
   gsl_rng_set(rng, dttk);
   /******************************************************************/
   N = (double) SIZE;
   dblindx = (double) indx;
   dbllag  = (double) lag;
   p = (double *) smpls;
   q = (double *) smpls + SIZE;
   while (p < q)
      {
      double num;
      num = gsl_rng_uniform(rng);
      *p++ = num;
      } /* for each sample */
   /******************************************************************/
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
   printf("      Autocorrelation Test\n");
   printf("GNU Scientific Library "
      "Generator: %s\n", gsl_rng_name(rng));
   printf("N %.0f index %.0f lag %.0f M %.0f\n",
      N, dblindx, dbllag, M);
   abszedzero = fabs(zedzero);
   pvalue = gsl_sf_erf_Q(abszedzero);
   printf("Z0 %f  P-value %f\n", zedzero, pvalue + pvalue);
   if (abszedzero < 0.01 || abszedzero > 1.96)
      printf("*************************************"   
         "*************************************\n");
   gsl_rng_free(rng);
   free(smpls);
   return(0);
   } /* main */
