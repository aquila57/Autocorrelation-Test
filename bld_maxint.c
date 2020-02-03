/* bld_maxint.c - Build MAXINT Table  Version 0.1.0 */
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

/* This program builds the maxint_tbl  */

void bld_maxint(xxfmt *xx)
   {
   double *p,*q;
   /***********************************************/
   /* zero fill the maxint table                  */
   /***********************************************/
   p = (double *) xx->maxint_tbl;
   q = (double *) xx->maxint_tbl + 2048;
   while (p < q) *p++ = 0.0;
   p = (double *) xx->maxint_tbl;
   *p++ = (double) 0xffffffff;     /* 000 borosh13 */
   *p++ = (double) 0x7fffffff;     /* 001 cmrg */
   *p++ = (double) 0xffffffff;     /* 002 coveyou */
   *p++ = (double) 0x7fffffff;     /* 003 fishman18 */
   *p++ = (double) 0x7fffffff;     /* 004 fishman20 */
   *p++ = (double) 0x7fffffff;     /* 005 fishman2x */
   *p++ = (double) 0xffffffff;     /* 006 gfsr4 */
   *p++ = (double) 0x3fffffff;     /* 007 knuthran */
   *p++ = (double) 0x7fffffff;     /* 008 knuthran2 */
   *p++ = (double) 0x3fffffff;     /* 009 knuthran2002 */
   *p++ = (double) 0x7fffffff;     /* 010 lecuyer21 */
   *p++ = (double) 0x7fffffff;     /* 011 minstd */
   *p++ = (double) 0x7fffffff;     /* 012 mrg */
   *p++ = (double) 0xffffffff;     /* 013 mt19937 */
   *p++ = (double) 0xffffffff;     /* 014 mt19937_1999 */
   *p++ = (double) 0xffffffff;     /* 015 mt19937_1998 */
   *p++ = (double) 0xffffffff;     /* 016 r250 */
   *p++ = (double) 0x7fffffff;     /* 017 ran0 */
   *p++ = (double) 0x7fffffff;     /* 018 ran1 */
   *p++ = (double) 0x7fffffff;     /* 019 ran2 */
   *p++ = (double) 999999999;      /* 020 ran3 */
   *p++ = (double) 0x7fffffff;     /* 021 rand */
   *p++ = (double) 0xffffffff;     /* 022 rand48 */
   *p++ = (double) 0x7fffffff;     /* 023 random128-bsd    */
   *p++ = (double) 0x7fffffff;     /* 024 random128-glibc2 */
   *p++ = (double) 0x7fffffff;     /* 025 random128-libc5  */
   *p++ = (double) 0x7fffffff;     /* 026 random256-bsd    */
   *p++ = (double) 0x7fffffff;     /* 027 random256-glibc2 */
   *p++ = (double) 0x7fffffff;     /* 028 random256-libc5  */
   *p++ = (double) 0x7fffffff;     /* 029 random32-bsd     */
   *p++ = (double) 0x7fffffff;     /* 030 random32-glibc2  */
   *p++ = (double) 0x7fffffff;     /* 031 random32-libc5   */
   *p++ = (double) 0x7fffffff;     /* 032 random64-bsd     */
   *p++ = (double) 0x7fffffff;     /* 033 random64-glibc2  */
   *p++ = (double) 0x7fffffff;     /* 034 random64-libc5   */
   *p++ = (double) 0x7fffffff;     /* 035 random8-bsd      */
   *p++ = (double) 0x7fffffff;     /* 036 random8-glibc2   */
   *p++ = (double) 0x7fffffff;     /* 037 random8-libc5    */
   *p++ = (double) 0x7fffffff;     /* 038 random-bsd       */
   *p++ = (double) 0x7fffffff;     /* 039 random-glibc2    */
   *p++ = (double) 0x7fffffff;     /* 040 random-libc5     */
   *p++ = (double) 0x7fffffff;     /* 041 randu */
   *p++ = (double) 0xffffffff;     /* 042 ranf */
   *p++ = (double) 0xffffff;       /* 043 ranlux */
   *p++ = (double) 0xffffff;       /* 044 ranlux389 */
   *p++ = (double) 0xffffffff;     /* 045 ranlxd1 */
   *p++ = (double) 0xffffffff;     /* 046 ranlxd2 */
   *p++ = (double) 0xffffff;       /* 047 ranlxs0 */
   *p++ = (double) 0xffffff;       /* 048 ranlxs1 */
   *p++ = (double) 0xffffff;       /* 049 ranlxs2 */
   *p++ = (double) 0xffffff;       /* 050 ranmar */
   *p++ = (double) 0x3fffff;       /* 051 slatec */
   *p++ = (double) 0xffffffff;     /* 052 taus */
   *p++ = (double) 0xffffffff;     /* 053 taus2 */
   *p++ = (double) 0xffffffff;     /* 054 taus113 */
   *p++ = (double) 0xffffffff;     /* 055 transputer */
   *p++ = (double) 0xffffffff;     /* 056 tt800 */
   *p++ = (double) 0xffff;         /* 057 uni */
   *p++ = (double) 0x7fffffff;     /* 058 uni32 */
   *p++ = (double) 0xffffffff;     /* 059 vax */
   *p++ = (double) 0xffffffff;     /* 060 waterman14 */
   *p++ = (double) 0xffffff;       /* 061 zuf */
   p = (double *) xx->maxint_tbl + 203;
   *p++ = (double) 0xffffffff;     /* 203 ca */
   *p++ = (double) 0xffffffff;     /* 204 uvag */
   *p++ = (double) 0xffffffff;     /* 205 AES_OFB */
   *p++ = (double) 0xffffffff;     /* 206 Threefish_OFB */
   *p++ = (double) 0xffffffff;     /* 207 XOR (supergenerator) */
   *p++ = (double) 0xffffffff;     /* 208 kiss */
   *p++ = (double) 0xffffffff;     /* 209 superkiss */
   p = (double *) xx->maxint_tbl + 400;
   *p++ = (double) 0xffffffff;     /* 400 R_wichmann_hill */
   *p++ = (double) 0xffffffff;     /* 401 R_marsaglia_multic. */
   *p++ = (double) 0xffffffff;     /* 402 R_super_duper */
   *p++ = (double) 0xffffffff;     /* 403 R_mersenne_twister */
   *p++ = (double) 0x3fffffff;     /* 404 R_knuth_taocp */
   *p++ = (double) 0x3fffffff;     /* 405 R_knuth_taocp2 */
   } /* bld_maxint */
