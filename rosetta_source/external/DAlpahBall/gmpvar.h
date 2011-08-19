/* Author		Patrice Koehl
 * Revision #		1
 * Date			6/10/2005
 * 
 * This file defines all global variables for gmp.
 */

#include "defines.h"

#ifndef __GMPVAR__
#define __GMPVAR__


/******************************************************************
 * include gmp if necessary
 *
 **********************************************************************/
#include "gmp.h"

/**********************************************************************
 * Define gmp array for regular triangulation
 **********************************************************************/

mpz_t coord_gmp[MAX_COORD];
mpz_t radius_gmp[MAX_ATOM];
mpz_t weight_gmp[MAX_ATOM];
//mpz_t coord_gmp[30];
//mpz_t radius_gmp[10];
//mpz_t weight_gmp[10];

mpz_t a11_mp,a12_mp,a13_mp,a14_mp;
mpz_t a21_mp,a22_mp,a23_mp,a24_mp;
mpz_t a31_mp,a32_mp,a33_mp,a34_mp;
mpz_t a41_mp,a42_mp,a43_mp,a44_mp;
mpz_t a51_mp,a52_mp,a53_mp,a54_mp;

mpz_t temp1,temp2,temp3,temp4;
mpz_t val1,val2,val3;

mpz_t c11,c12,c13,c14,c21,c22,c23,c24,c31,c32,c33,c34,c41,c42,c43,c44;
mpz_t d1,d2,d3,e1,e2,e3,f1,f2,f3,g1,g2,g3;

/**********************************************************************
 * Define gmp array for dual complex
 **********************************************************************/

mpz_t ra2,rb2,dist2,dtest, num, den;
mpz_t r_11, r_22, r_33, r_14, r_313, r_212,diff, det0, det1, det2, det3, det4;
mpz_t Dabc, Dabd, Dacd, Dbcd, Dabcd;
mpz_t wa,wb,wc,wd;

mpz_t ra_mp,rb_mp, rc_mp, rd_mp;
mpz_t alp;

mpz_t res[4][5], res2_c[4][5];
mpz_t a_mp[5], b_mp[5],c_mp[5],d_mp[5];
mpz_t Tab[4], Sab[4], Dab[5];
mpz_t Sac[4], Sad[4], Sbc[4], Sbd[4], Scd[4];
mpz_t Sa[4], Sb[4], Sd[4];
mpz_t Sam1[4], Sbm1[4], Scm1[4], Sdm1[4];
mpz_t Deter[4];
mpz_t Tc[4],Sc[4];
mpz_t Mab[4][5], Mac[4][5], Mbc[4][5], S[4][5], T[3][4];


#endif /* __GMPVAR__ */
