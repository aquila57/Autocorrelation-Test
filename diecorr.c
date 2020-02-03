/* corr.c - Autocorrelation, dieharder generators  Version 0.1.0 */
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

void putstx(char *pgm)
   {
   fprintf(stderr,"Autocorrelation test, dieharder "
      "random number generators\n");
   fprintf(stderr,"Usage: %s index lag generator\n", pgm);
   fprintf(stderr,"Where index is i, the starting offset\n");
   fprintf(stderr,"Where lag is m, the distance between samples\n");
   fprintf(stderr,"Where generator is the dieharder ");
   fprintf(stderr,"random number generator\n");
   fprintf(stderr,"Example: %s 3 5 053\n", pgm);
   fprintf(stderr,"053 is the id for the taus2 generator\n");
   fprintf(stderr,"To see a list of dieharder generators, type:\n");
   fprintf(stderr,"dieharder -g -l\n");
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
   xxfmt *xx;
   /*************************************************************/
   /* Allocate memory for the global structure.                 */
   /* This is a re-entrant program.                             */
   /*************************************************************/
   xx = (xxfmt *) malloc(sizeof(xxfmt));
   if (xx == NULL)
      {
      fprintf(stderr,"main: out of memory "
         "allocating xx\n");
      exit(1);
      } /* out of memory */
   /*************************************************************/
   /* initialize global data                                    */
   /*************************************************************/
   xx->dieharder_rngname[0] = '\0';
   xx->rngname = (char *) NULL;
   xx->generator = -1;
   /*************************************************************/
   /* obtain dieharder generator number parameter               */
   /*************************************************************/
   if (argc != 4) putstx(*argv);
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
   xx->generator = atoi(*(argv+3));
   if (xx->generator < 0 || xx->generator > 405)
      {
      fprintf(stderr,"main: invalid random "
         "number generator %s\n", *(argv+3));
      putstx(*argv);
      } /* invalid generator */
   /*************************************************************/
   /* initialize numeric fields                                 */
   /*************************************************************/
   bld_maxint(xx);  /* build maxint table for dieharder RNGs */
   bld_rngtbl(xx);  /* build name table for dieharder RNGs */
   if (xx->maxint_tbl[xx->generator] == 0.0)
      {
      fprintf(stderr,"main: invalid random "
         "number generator %s\n", *(argv+1));
      putstx(*argv);
      } /* invalid generator */
   xx->modulus = xx->maxint_tbl[xx->generator];
   xx->maxint  = (unsigned int) xx->modulus;
   /*************************************************************/
   smpls = (double *) malloc(sizeof(double) * SIZE + 10);
   if (smpls == NULL)
      {
      fprintf(stderr,"main: out of memory "
         "allocating smpls\n");
      exit(1);
      } /* out of memory */
   /***************************************************/
   /* Initialize the dieharder generator              */
   /* Bypass the dieharder prefix                     */
   /***************************************************/
   diepfx(xx);
   /******************************************************************/
   N = (double) SIZE;
   dblindx = (double) indx;
   dbllag  = (double) lag;
   p = (double *) smpls;
   q = (double *) smpls + SIZE;
   while (p < q)
      {
      double num;
      num = getdie(xx);
      if (xx->eofsw)
         {
	 fprintf(stderr,"main: end of dieharder input\n");
	 exit(1);
	 } /* end of file */
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
   printf("  Autocorrelation Test\n");
   printf("Dieharder Generator: %s\n",
      xx->dieharder_rngname);
   printf("N %.0f index %.0f lag %.0f M %.0f\n",
      N, dblindx, dbllag, M);
   abszedzero = fabs(zedzero);
   pvalue = gsl_sf_erf_Q(abszedzero);
   printf("Z0 %f  P-value %f\n", zedzero, pvalue + pvalue);
   if (abszedzero < 0.01 || abszedzero > 1.96)
      printf("*************************************"   
         "*************************************\n");
   free(smpls);
   free(xx);
   return(0);
   } /* main */
