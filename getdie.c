/* getdie.c - Read one line from Dieharder  Version 0.1.0 */
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

/* Read one line of random output from Dieharder      */
/* The Dieharder command to produce random output is: */
/* dieharder -g nnn -t sss -o                         */
/* where nnn is the number of the Dieharder random    */
/* number generator and                               */
/* sss is the number of 32-bit words to generate      */
/* The output from Dieharder is in ASCII format.      */
/* Dieharder writes one line per 32-bit generation    */
/* of output.                                         */
/* The output of getdie() is a uniform random         */
/* number from zero to one.                           */

double getdie(xxfmt *xx)
   {
   int len;
   double num;
   double frac;
   char buf[256];
   len = scanf("%s", buf);
   if (len < 1)
      {
      if (feof(stdin))
	 {
	 xx->eofsw = 1;
         return(EOFDIE);
	 } /* if end of file */
      if (ferror(stdin))
         {
         perror("getdie: read error\n");
         exit(1);
         } /* read error */
      fprintf(stderr,"getdie: blank line\n");
      exit(1);
      } /* if eof or read error */
   if (len > 1)
      {
      fprintf(stderr,"getdie: read error\n");
      fprintf(stderr,"scanf read more than "
         "one token %d\n", len);
      exit(1);
      } /* invalid input */
   num = atof(buf);
   /* maxint is the maximum integer for the dieharder */
   /* random number generator in this run */
   frac = num / xx->modulus;
   return(frac);
   } /* getdie */
