/* diepfx.c - Bypass dieharder prefix  Version 0.1.0 */
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

void diepfx(xxfmt *xx)
   {
   int len;
   int rslt;
   int namesw;
   char buf[256];
   namesw = 0;
   xx->dieharder_rngname[0] = '\0';
   while (1)
      {
      len = scanf("%s", buf);
      if (len < 1)
         {
         if (feof(stdin))
	    {
	    fprintf(stderr,"diepfx: scanf end of file\n");
	    exit(1);
	    } /* if end of file */
         if (ferror(stdin))
            {
            perror("diepfx: read error\n");
            exit(1);
            } /* read error */
         fprintf(stderr,"diepfx: blank line, len %d\n", len);
         exit(1);
         } /* if eof or read error */
      if (len > 1)
         {
         fprintf(stderr,"diepfx: read error\n");
         fprintf(stderr,"scanf read more than "
            "one token %d\n", len);
         exit(1);
         } /* invalid input */
      rslt = strcmp(buf,"numbit:");
      if (!rslt) break;
      if (namesw)   /* buf = generator name? */
         {
	 strcpy(xx->dieharder_rngname,buf);
	 fprintf(stderr,"Dieharder RNG Name %s\n",
	    xx->dieharder_rngname);
	 namesw = 0;    /* turn off name switch */
	 } /* if generator name switch turned on */
      rslt = strcmp(buf,"generator");
      if (!rslt)
         {
	 namesw = 1;    /* turn on name switch */
	 } /* if generator statement */
      } /* read loop until numbit: */
   len = scanf("%s", buf);
   if (len < 1)
      {
      if (feof(stdin))
	 {
	 fprintf(stderr,"diepfx: scanf end of file\n");
	 exit(1);
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
   rslt = strcmp(buf,"32");
   if (rslt)
      {
      fprintf(stderr,"diepfx: numbit not = 32\n");
      fprintf(stderr,"buf = %s\n", buf);
      exit(1);
      } /* numbit: not equal to 32 */
   } /* diepfx */
