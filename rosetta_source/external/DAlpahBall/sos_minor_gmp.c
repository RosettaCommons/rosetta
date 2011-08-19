/*	Sos_minor_gmp.c	Version 1 1/16/2002	Patrice Koehl             */
/*									  */
/*  This is the C version of sos_minor.f, which performs all operations   */
/*  with multi precision arithmetics, using the package GMP               */
/*									  */
/*------------------------------------------------------------------------*/
/*                                                                        */
/* Copyright (C) 2002 Patrice Koehl                                       */
/*                                                                        */
/* This library is free software; you can redistribute it and/or          */
/* modify it under the terms of the GNU Lesser General Public             */
/* License as published by the Free Software Foundation; either           */
/* version 2.1 of the License, or (at your option) any later version.     */
/*                                                                        */
/* This library is distributed in the hope that it will be useful,        */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of         */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU       */
/* Lesser General Public License for more details.                        */
/*                                                                        */
/* You should have received a copy of the GNU Lesser General Public       */
/* License along with this library; if not, write to the Free Software    */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA*/
/*                                                                        */
/*------------------------------------------------------------------------*/
/*                                                                        */
/* Includes :								  */

#include <stdio.h>
#include <math.h>
#include "gmpvar.h"

/*****************************************************************************************
 * Naming convention between C and Fortran
 *
 * Let us consider two Fortran subroutines : foo, and foo_with_underscore
 *
 * Case 1: the name of the Fortran subroutine (foo) does not contain an underscore (_)
 *         In the C program, we must add ONE underscore to the Fortran name:
 *         the C program refers to the subroutine as:
 *              foo_(...)
 *         while the Fortran code writes it as:
 *              foo(...)
 *         This is independent of the compiler pair (at least for gcc/f77 and
 *         Intel icc/ifort)
 *
 * Case 2: the name of the Fortran subroutine (foo_with_underscore) contains at least one
 * underscore (_)
 *         Treatment of case 2 is compiler dependent!
 *         - The Intel compiler treats this case as if it was case 1, i.e. requires
 *           ONE underscore at the end of the Fortran name in the C program
 *         - the gnu f77 however requires TWO underscores at the end of the Fortran name
 *           in the C program.
 *
 * I solve this by introducing two functions, FTName1 and FTName2, where FTName2 is
 * compiler dependent
 */

#define F77Name1(x) x##_   /* Need to add one underscore to Fortran program */

#if defined intel
#define F77Name2(x) x##_   /* Need to add one underscore to Fortran program for
                                Intel compiler */
#else
#define F77Name2(x) x##__   /* Need to add two underscores to Fortran program for GNU
                                f77 compiler */
#endif

#define round(x) ((x)>=0?(int)((x)+0.5):(int)((x)-0.5))

/*------------------------------------------------------------------------*/
/*Local procedures:							  */

void deter2_gmp(mpz_t deter, mpz_t a, mpz_t b);

void deter3_gmp(mpz_t deter, mpz_t a11, mpz_t a12, mpz_t a21, mpz_t a22,
		mpz_t a31, mpz_t a32);

void deter4_gmp(mpz_t deter, mpz_t a11, mpz_t a12, mpz_t a13, mpz_t a21,
		mpz_t a22, mpz_t a23, mpz_t a31, mpz_t a32,
		mpz_t a33, mpz_t a41, mpz_t a42, mpz_t a43);

void deter5_gmp(mpz_t deter, mpz_t a11, mpz_t a12, mpz_t a13, mpz_t a14,
		mpz_t a21, mpz_t a22, mpz_t a23, mpz_t a24,
		mpz_t a31, mpz_t a32, mpz_t a33, mpz_t a34,
		mpz_t a41, mpz_t a42, mpz_t a43, mpz_t a44,
		mpz_t a51, mpz_t a52, mpz_t a53, mpz_t a54);

#ifdef __cplusplus
extern "C"
#endif
void F77Name2(transfer_all_to_gmp)(double *coord, double *radius, double *scale,
	int *nvertex);

#ifdef __cplusplus
extern "C"
#endif
void F77Name2(transfer_one_to_gmp)(double *coord, double *radius, double *scale,
	int *nvertex);

#ifdef __cplusplus
extern "C"
#endif
void clear_all_gmp_arrays(int *nvertex);

#ifdef __cplusplus
extern "C"
#endif
void F77Name2(sos_minor2_gmp)(int *a, int *b, int *i1, int *res);
#ifdef __cplusplus
extern "C"
#endif
void F77Name2(sos_minor3_gmp)(int *a, int *b, int *c, int *i1, int *i2, int *res);

#ifdef __cplusplus
extern "C"
#endif
void F77Name2(sos_minor4_gmp)(int *a, int *b, int *c, int *d, int *res);
#ifdef __cplusplus
extern "C"
#endif
void F77Name2(sos_minor5_gmp)(int *a, int *b, int *c, int *d, int *e, int *res);

#ifdef __cplusplus
extern "C"
#endif
void F77Name2(minor4_gmp)(int *a, int *b, int *c, int *d, int *res);

/*									  */
/*------------------------------------------------------------------------*/
/* transfer_all_to_gmp:							  */
/* This subroutine converts all coordinates into mpz_t type
                                                                          */
#ifdef __cplusplus
extern "C" {
#endif
void F77Name2(transfer_all_to_gmp)(double *coord, double *radius, double *scale,
		int *nvertex)
{
	int i;
	double x,value;
	int ivalue;

	mpz_init(a11_mp);mpz_init(a12_mp); mpz_init(a13_mp); mpz_init(a14_mp);
	mpz_init(a21_mp);mpz_init(a22_mp); mpz_init(a23_mp); mpz_init(a24_mp);
	mpz_init(a31_mp);mpz_init(a32_mp); mpz_init(a33_mp); mpz_init(a34_mp);
	mpz_init(a41_mp);mpz_init(a42_mp); mpz_init(a43_mp); mpz_init(a44_mp);
	mpz_init(a51_mp);mpz_init(a52_mp); mpz_init(a53_mp); mpz_init(a54_mp);

	mpz_init(temp1); mpz_init(temp2); mpz_init(temp3); mpz_init(temp4);
	mpz_init(val1);mpz_init(val2); mpz_init(val3);

	mpz_init(c11); mpz_init(c12); mpz_init(c13); mpz_init(c14);
	mpz_init(c21); mpz_init(c22); mpz_init(c23); mpz_init(c24);
	mpz_init(c31); mpz_init(c32); mpz_init(c33); mpz_init(c34);
	mpz_init(c41); mpz_init(c42); mpz_init(c43); mpz_init(c44);
	mpz_init(d1); mpz_init(d2); mpz_init(d3);
	mpz_init(e1); mpz_init(e2); mpz_init(e3);
	mpz_init(f1); mpz_init(f2); mpz_init(f3);
	mpz_init(g1); mpz_init(g2); mpz_init(g3);

	for ( i=0; i< 3*(*nvertex); i++)
	{
		mpz_init(coord_gmp[i]);
	}
	for ( i=0; i<(*nvertex); i++)
	{
		mpz_init(radius_gmp[i]);
		mpz_init(weight_gmp[i]);
	}

	mpz_set_d(temp3, *scale);

	for (i=0; i<3*(*nvertex); i++)
	{
		ivalue= (int) coord[i];
		mpz_set_si(temp1,ivalue);
		mpz_mul(temp1,temp1,temp3);
		x = (coord[i]-ivalue)*(*scale);
		ivalue = (int) round(x);
		mpz_set_si(temp2,ivalue);
		mpz_add(coord_gmp[i],temp1,temp2);
	}

	for (i=0; i<(*nvertex); i++)
	{
		ivalue = (int) radius[i];
		mpz_set_si(temp1,ivalue);
		mpz_mul(temp1,temp1,temp3);
		x = (radius[i]-ivalue)*(*scale);
		ivalue = (int) round(x);
		mpz_set_si(temp2,ivalue);
		mpz_add(radius_gmp[i],temp1,temp2);

		mpz_mul(temp1,radius_gmp[i],radius_gmp[i]);
		mpz_mul(temp2,coord_gmp[3*i+2],coord_gmp[3*i+2]), mpz_sub(temp1,temp2,temp1);
		mpz_mul(temp2,coord_gmp[3*i+1],coord_gmp[3*i+1]), mpz_add(temp1,temp2,temp1);
		mpz_mul(temp2,coord_gmp[3*i],coord_gmp[3*i]), mpz_add(weight_gmp[i],temp2,temp1);
	}

}
/*									  */
/*------------------------------------------------------------------------*/
/* transfer_one_to_gmp:							  */
/* This subroutine converts the coordinate of one point into mpz_t type
   and add it at the end of the current mpz_t arrays
                                                                          */
#ifdef __cplusplus
extern "C" {
#endif
void F77Name2(transfer_one_to_gmp)(double *coord, double *radius, double *scale,
		int *nvertex)
{
	int i;
	double x,value;
	int ivalue;

	for (i = 0; i<3 ; i++)
	{
		mpz_init(coord_gmp[3*(*nvertex)+i]);
	}
	mpz_init(radius_gmp[(*nvertex)]);
	mpz_init(weight_gmp[(*nvertex)]);

	mpz_set_d(temp3, *scale);

	for (i=0; i<3; i++)
	{
		ivalue= (int) coord[3*(*nvertex)+i];
		mpz_set_si(temp1,ivalue);
		mpz_mul(temp1,temp1,temp3);
		x = (coord[3*(*nvertex)+i]-ivalue)*(*scale);
		ivalue = (int) round(x);
		mpz_set_si(temp2,ivalue);
		mpz_add(coord_gmp[3*(*nvertex)+i],temp1,temp2);
	}

	ivalue = (int) radius[(*nvertex)];
	mpz_set_si(temp1,ivalue);
	mpz_mul(temp1,temp1,temp3);
	x = (radius[(*nvertex)]-ivalue)*(*scale);
	ivalue = (int) round(x);
	mpz_set_si(temp2,ivalue);
	mpz_add(radius_gmp[(*nvertex)],temp1,temp2);

	mpz_mul(temp1,radius_gmp[(*nvertex)],radius_gmp[(*nvertex)]);
	mpz_mul(temp2,coord_gmp[3*(*nvertex)+2],coord_gmp[3*(*nvertex)+2]), mpz_sub(temp1,temp2,temp1);
	mpz_mul(temp2,coord_gmp[3*(*nvertex)+1],coord_gmp[3*(*nvertex)+1]), mpz_add(temp1,temp2,temp1);
	mpz_mul(temp2,coord_gmp[3*(*nvertex)],coord_gmp[3*(*nvertex)]), mpz_add(weight_gmp[(*nvertex)],temp2,temp1);

}
/*									  */
/*------------------------------------------------------------------------*/
/* clear_all_gmp_arrays:							  */
/* This subroutine clears all pending mpz_t arrays
                                                                          */

#ifdef __cplusplus
extern "C" {
#endif
void clear_all_gmp_arrays(int *nvertex)
{
	int i;

	for ( i=0; i< 3*(*nvertex); i++)
	{
		mpz_clear(coord_gmp[i]);
	}
	for ( i=0; i<(*nvertex); i++)
	{
		mpz_clear(radius_gmp[i]);
		mpz_clear(weight_gmp[i]);
	}

	mpz_clear(a11_mp);mpz_clear(a12_mp); mpz_clear(a13_mp); mpz_clear(a14_mp);
	mpz_clear(a21_mp);mpz_clear(a22_mp); mpz_clear(a23_mp); mpz_clear(a24_mp);
	mpz_clear(a31_mp);mpz_clear(a32_mp); mpz_clear(a33_mp); mpz_clear(a34_mp);
	mpz_clear(a41_mp);mpz_clear(a42_mp); mpz_clear(a43_mp); mpz_clear(a44_mp);
	mpz_clear(a51_mp);mpz_clear(a52_mp); mpz_clear(a53_mp); mpz_clear(a54_mp);
	mpz_clear(temp1); mpz_clear(temp2); mpz_clear(temp3); mpz_clear(temp4);
	mpz_clear(val1);mpz_clear(val2);mpz_clear(val3);

	mpz_clear(c11); mpz_clear(c12); mpz_clear(c13); mpz_clear(c14);
	mpz_clear(c21); mpz_clear(c22); mpz_clear(c23); mpz_clear(c24);
	mpz_clear(c31); mpz_clear(c32); mpz_clear(c33); mpz_clear(c34);
	mpz_clear(c41); mpz_clear(c42); mpz_clear(c43); mpz_clear(c44);
	mpz_clear(d1); mpz_clear(d2); mpz_clear(d3);
	mpz_clear(e1); mpz_clear(e2); mpz_clear(e3);
	mpz_clear(f1); mpz_clear(f2); mpz_clear(f3);
	mpz_clear(g1); mpz_clear(g2); mpz_clear(g3);
}
/*									  */
/*									  */
/*------------------------------------------------------------------------*/
/* deter2_gmp:								  */
/*									  */
/*	This subroutine evaluates the determinant:
		D = | b11 1 |
		    | b21 1 |

	Input:
		b11, b21
	Output:
		deter
									  */

void deter2_gmp(mpz_t deter, mpz_t b11, mpz_t b21)

{
	mpz_sub(deter,b11,b21);
}
/*									  */
/*------------------------------------------------------------------------*/
/* deter3_gmp:								  */
/*									  */
/*	This subroutine evaluates the determinant:
		D = | b11 b12 1 |
		    | b21 b22 1 |
		    | b31 b32 1 |

	Input:
		b11, b12, b21, b22, b31, b32
	Output:
		deter3_gmp
									  */

void deter3_gmp(mpz_t deter, mpz_t b11, mpz_t b12, mpz_t b21, 
		mpz_t b22, mpz_t b31, mpz_t b32)

{

	mpz_sub(temp1,b21,b11);
	mpz_sub(temp2,b22,b12);
	mpz_sub(temp3,b31,b11);
	mpz_sub(temp4,b32,b12);

	mpz_mul(val1,temp1,temp4);
	mpz_mul(val2,temp2,temp3);

	mpz_sub(deter,val1,val2);

}
/*									  */
/*------------------------------------------------------------------------*/
/* deter4_gmp:								  */
/*									  */
/*	This subroutine evaluates the determinant:
		D = | b11 b12 b13 1 |
		    | b21 b22 b23 1 |
		    | b31 b32 b33 1 |
		    | b41 b42 b43 1 |

	Input:
		b11, b12, b13, b21, b22, b23, b31, b32, b33
		b41, b42, b43
	Output:
		deter4_gmp
									  */

void deter4_gmp(mpz_t deter, mpz_t b11, mpz_t b12, mpz_t b13, mpz_t b21, 
		mpz_t b22, mpz_t b23, mpz_t b31, mpz_t b32, mpz_t b33, 
		mpz_t b41, mpz_t b42, mpz_t b43)

{
	mpz_sub(c11,b21,b11);mpz_sub(c12,b22,b12);mpz_sub(c13,b23,b13);
	mpz_sub(c21,b31,b11);mpz_sub(c22,b32,b12);mpz_sub(c23,b33,b13);
	mpz_sub(c31,b41,b11);mpz_sub(c32,b42,b12);mpz_sub(c33,b43,b13);

	mpz_mul(temp1,c22,c33);mpz_mul(temp2,c32,c23);mpz_sub(val1,temp1,temp2);
	mpz_mul(temp1,c12,c33);mpz_mul(temp2,c32,c13);mpz_sub(val2,temp1,temp2);
	mpz_mul(temp1,c12,c23);mpz_mul(temp2,c22,c13);mpz_sub(val3,temp1,temp2);

	mpz_mul(temp1,c21,val2);mpz_mul(temp2,c11,val1);mpz_mul(temp3,c31,val3);

	mpz_add(val1,temp2,temp3);
	mpz_sub(deter,temp1,val1);

}
/*									  */
/*------------------------------------------------------------------------*/
/* deter5_gmp:								  */
/*									  */
/*	This subroutine evaluates the determinant:
		D = | b11 b12 b13 b14 1 |
		    | b21 b22 b23 b24 1 |
		    | b31 b32 b33 b34 1 |
		    | b41 b42 b43 b44 1 |
		    | b51 b52 b53 b54 1 |

	Input:
		b11, b12, b13, b14, b21, b22, b23, b24, b31, b32, b33, b34
		b41, b42, b43, b44, b51, b52, b53, b54
	Output:
		deter5_gmp
									  */

void deter5_gmp(mpz_t deter, mpz_t b11, mpz_t b12, mpz_t b13, mpz_t b14, 
		 mpz_t b21, mpz_t b22, mpz_t b23, mpz_t b24,
		 mpz_t b31, mpz_t b32, mpz_t b33, mpz_t b34,
		 mpz_t b41, mpz_t b42, mpz_t b43, mpz_t b44,
		 mpz_t b51, mpz_t b52, mpz_t b53, mpz_t b54)
{

	mpz_sub(c11,b21,b11); mpz_sub(c12,b22,b12); mpz_sub(c13,b23,b13);
	mpz_sub(c14,b24,b14);
	mpz_sub(c21,b31,b11); mpz_sub(c22,b32,b12); mpz_sub(c23,b33,b13);
	mpz_sub(c24,b34,b14);
	mpz_sub(c31,b41,b11); mpz_sub(c32,b42,b12); mpz_sub(c33,b43,b13);
	mpz_sub(c34,b44,b14);
	mpz_sub(c41,b51,b11); mpz_sub(c42,b52,b12); mpz_sub(c43,b53,b13);
	mpz_sub(c44,b54,b14);

	mpz_mul(temp1,c32,c43); mpz_mul(temp2,c42,c33); mpz_sub(d1,temp1,temp2);
	mpz_mul(temp1,c32,c44); mpz_mul(temp2,c42,c34); mpz_sub(d2,temp1,temp2);
	mpz_mul(temp1,c33,c44); mpz_mul(temp2,c43,c34); mpz_sub(d3,temp1,temp2);

	mpz_mul(temp1,c12,c23); mpz_mul(temp2,c22,c13); mpz_sub(e1,temp1,temp2);
	mpz_mul(temp1,c12,c24); mpz_mul(temp2,c22,c14); mpz_sub(e2,temp1,temp2);
	mpz_mul(temp1,c13,c24); mpz_mul(temp2,c23,c14); mpz_sub(e3,temp1,temp2);

	mpz_mul(temp1,c11,c24); mpz_mul(temp2,c21,c14); mpz_sub(f1,temp1,temp2);
	mpz_mul(temp1,c11,c23); mpz_mul(temp2,c21,c13); mpz_sub(f2,temp1,temp2);
	mpz_mul(temp1,c11,c22); mpz_mul(temp2,c21,c12); mpz_sub(f3,temp1,temp2);

	mpz_mul(temp1,c31,c44); mpz_mul(temp2,c41,c34); mpz_sub(g1,temp1,temp2);
	mpz_mul(temp1,c31,c43); mpz_mul(temp2,c41,c33); mpz_sub(g2,temp1,temp2);
	mpz_mul(temp1,c31,c42); mpz_mul(temp2,c41,c32); mpz_sub(g3,temp1,temp2);
 
	mpz_mul(temp1,e3,g3); mpz_mul(temp2,e2,g2); mpz_sub(temp3,temp1,temp2);
	mpz_mul(temp1,e1,g1); mpz_add(temp3,temp3,temp1);
	mpz_mul(temp1,d3,f3); mpz_add(temp3,temp3,temp1);
	mpz_mul(temp1,d2,f2); mpz_sub(temp3,temp3,temp1);
	mpz_mul(temp1,d1,f1); mpz_add(deter,temp3,temp1);

}
/*									  */
/*------------------------------------------------------------------------*/
/* sos_minor2_gmp:							  */
/*									  
	This subroutine tests the sign of the determinant
		D = | a11 1 |
		    | a21 1 |
	If the determinant is found to be 0, then the SoS procedure is used:
	a development of the determinant with respect to a perturbation EPS
	applied to the coordinates in the determinant is computed, and
	the sign of the first non zero term defines the sign of the 
	determinant.
	In the case of a 2x2 determinant, the first term in the expansion
	is the coefficient 1 ...					*/
#ifdef __cplusplus
extern "C" {
#endif
void F77Name2(sos_minor2_gmp)(int *a, int *b, int *ia, int *res)

{
	int icomp,idx_a,idx_b;

/* Initialise local GMP variables */

/* Get coordinates */

	idx_a = (*a)*3 + (*ia) -4;
	idx_b = (*b)*3 + (*ia) -4;
	mpz_set(a11_mp,coord_gmp[idx_a]);
	mpz_set(a21_mp,coord_gmp[idx_b]);

/* Compute determinant */

	deter2_gmp(temp1,a11_mp,a21_mp);

	icomp = mpz_sgn(temp1);

	if (icomp != 0) {
		*res = icomp;
	}
	else {
		*res = 1;
	}

/* Clear local GMP variables */

}
#ifdef __cplusplus
}
#endif
/*									  */
/*------------------------------------------------------------------------*/
/* sos_minor3_gmp:							  */
/*									 
	This subroutine tests the sign of the determinant
		D = | a11 a12 1 |
		    | a21 a22 1 |
		    | a31 a32 1 |
	If the determinant is found to be 0, then the SoS procedure is used:
	a development of the determinant with respect to a perturbation EPS
	applied to the coordinates in the determinant is computed, and
	the sign of the first non zero term defines the sign of the 
	determinant.
	In the case of a 3x3 determinant, the maximum number of terms to be
	checked is 4 ...					*/
#ifdef __cplusplus
extern "C" {
#endif
void F77Name2(sos_minor3_gmp)(int *a, int *b, int *c, int *i1, int *i2, int *res)
{
	int icomp;
	int idx_a1,idx_a2,idx_b1,idx_b2,idx_c1,idx_c2;
	int ivalue;

	double value;

/* Initialise local GMP variables */


/* Transfer coordinates to GMP */

	idx_a1 = (*a)*3 + (*i1) -4; idx_a2 = (*a)*3 + (*i2) -4;
	idx_b1 = (*b)*3 + (*i1) -4; idx_b2 = (*b)*3 + (*i2) -4;
	idx_c1 = (*c)*3 + (*i1) -4; idx_c2 = (*c)*3 + (*i2) -4;

	mpz_set(a11_mp,coord_gmp[idx_a1]);
	mpz_set(a12_mp,coord_gmp[idx_a2]);
	mpz_set(a21_mp,coord_gmp[idx_b1]);
	mpz_set(a22_mp,coord_gmp[idx_b2]);
	mpz_set(a31_mp,coord_gmp[idx_c1]);
	mpz_set(a32_mp,coord_gmp[idx_c2]);

/* Compute determinant */

	deter3_gmp(temp1,a11_mp,a12_mp,a21_mp,a22_mp,a31_mp,a32_mp);

	icomp = mpz_sgn(temp1);

/* if major determinant is non 0, return its sign */

	if (icomp != 0) {
		*res = icomp;
		return;
	}

/* Look now at each term in the expansion of the determinant with
   respect to EPS                                               
	The initial determinant is:
		Minor3(i,j,k,1,2,0)				*/

/* Term 1: - Minor2(j,k,1,0) */

	deter2_gmp(temp1,a21_mp,a31_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res = - icomp;
		return;
	}

/* Term 2: Minor2(j,k,2,0)  */

	deter2_gmp(temp1,a22_mp,a32_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res = icomp;
		return;
	}

/* Term 3: Minor2(i,k,1,0)  */

	deter2_gmp(temp1,a11_mp,a31_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res = icomp;
		return;
	}

/* Term 4: 1 */

	*res = 1;

}
#ifdef __cplusplus
}
#endif
/*									  */
/*------------------------------------------------------------------------*/
/* sos_minor4_:								  */
/*									  
	This subroutine tests the sign of the determinant
		D = | a11 a12 a13 1 |
		    | a21 a22 a23 1 |
		    | a31 a32 a33 1 |
		    | a41 a42 a43 1 |
	If the determinant is found to be 0, then the SoS procedure is used:
	a development of the determinant with *respect to a perturbation EPS
	applied to the coordinates in the determinant is computed, and
	the sign of the first non zero term defines the sign of the 
	determinant.
	In the case of a 4x4 determinant, the maximum number of terms to be
	checked is 14 ...					*/
#ifdef __cplusplus
extern "C" {
#endif
void F77Name2(sos_minor4_gmp)(int *a, int *b, int *c, int *d, int *res)

{
	int icomp;
	int idx_a1,idx_a2,idx_a3,idx_b1,idx_b2,idx_b3;
	int idx_c1,idx_c2,idx_c3,idx_d1,idx_d2,idx_d3;
	int ivalue;
	int base=10;

	double value;


/* Transfer coordinates to gmp */

	idx_a1 = (*a)*3 -3; idx_a2 = idx_a1+1; idx_a3 = idx_a2 +1;
	idx_b1 = (*b)*3 -3; idx_b2 = idx_b1+1; idx_b3 = idx_b2 +1;
	idx_c1 = (*c)*3 -3; idx_c2 = idx_c1+1; idx_c3 = idx_c2 +1;
	idx_d1 = (*d)*3 -3; idx_d2 = idx_d1+1; idx_d3 = idx_d2 +1;

	mpz_set(a11_mp,coord_gmp[idx_a1]);
	mpz_set(a12_mp,coord_gmp[idx_a2]);
	mpz_set(a13_mp,coord_gmp[idx_a3]);
	mpz_set(a21_mp,coord_gmp[idx_b1]);
	mpz_set(a22_mp,coord_gmp[idx_b2]);
	mpz_set(a23_mp,coord_gmp[idx_b3]);
	mpz_set(a31_mp,coord_gmp[idx_c1]);
	mpz_set(a32_mp,coord_gmp[idx_c2]);
	mpz_set(a33_mp,coord_gmp[idx_c3]);
	mpz_set(a41_mp,coord_gmp[idx_d1]);
	mpz_set(a42_mp,coord_gmp[idx_d2]);
	mpz_set(a43_mp,coord_gmp[idx_d3]);

/*	printf("a11_mp = %s\n",mpz_get_str(NULL,base,a11_mp)); 
	printf("a12_mp = %s\n",mpz_get_str(NULL,base,a12_mp)); 
	printf("a13_mp = %s\n",mpz_get_str(NULL,base,a13_mp)); 
	printf("a21_mp = %s\n",mpz_get_str(NULL,base,a21_mp)); 
	printf("a22_mp = %s\n",mpz_get_str(NULL,base,a22_mp)); 
	printf("a23_mp = %s\n",mpz_get_str(NULL,base,a23_mp)); 
	printf("a31_mp = %s\n",mpz_get_str(NULL,base,a31_mp)); 
	printf("a32_mp = %s\n",mpz_get_str(NULL,base,a32_mp)); 
	printf("a33_mp = %s\n",mpz_get_str(NULL,base,a33_mp)); 
	printf("a41_mp = %s\n",mpz_get_str(NULL,base,a41_mp)); 
	printf("a42_mp = %s\n",mpz_get_str(NULL,base,a42_mp)); 
	printf("a43_mp = %s\n",mpz_get_str(NULL,base,a43_mp));  */

/* Compute determinant */

	deter4_gmp(temp1,a11_mp,a12_mp,a13_mp,a21_mp,a22_mp,a23_mp,
			   a31_mp,a32_mp,a33_mp,a41_mp,a42_mp,a43_mp);

//	printf("Deter4 = %s\n",mpz_get_str(NULL,base,temp1)); 

	icomp = mpz_sgn(temp1);

/* if major determinant is non 0, return its sign 
   (don't forget to clear GMP variables !!)	*/

	if (icomp != 0) {
		*res = icomp;
		return;
	}

/* Look now at each term in the expansion of the determinant with
   *respect to EPS                                                
	The initial determinant is:
		Minor4(i,j,k,l,1,2,3,0)				*/

/* Term 1:	Minor3(j,k,l,1,2,0)  */

	deter3_gmp(temp1,a21_mp,a22_mp,a31_mp,a32_mp,a41_mp,a42_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res = icomp;
		return;
	}

/* Term 2:	-Minor3(j,k,l,1,3,0) */

	deter3_gmp(temp1,a21_mp,a23_mp,a31_mp,a33_mp,a41_mp,a43_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res = - icomp;
		return;
	}

/* Term 3:	Minor3(j,k,l,2,3,0) */

	deter3_gmp(temp1,a22_mp,a23_mp,a32_mp,a33_mp,a42_mp,a43_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res = icomp;
		return;
	}

/* Term 4:	- Minor3(i,k,l,1,2,0) */

	deter3_gmp(temp1,a11_mp,a12_mp,a31_mp,a32_mp,a41_mp,a42_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res = - icomp;
		return;
	}

/* Term 5:	Minor2(k,l,1,0) */

	deter2_gmp(temp1,a31_mp,a41_mp);	
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res = icomp;
		return;
	}

/* Term 6:	-Minor2(k,l,2,0) */

	deter2_gmp(temp1,a32_mp,a42_mp);	
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res = - icomp;
		return;
	}

/* Term 7:	Minor3(i,k,l,1,3,0) */

	deter3_gmp(temp1,a11_mp,a13_mp,a31_mp,a33_mp,a41_mp,a43_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 8:	Minor2(k,l,3,0) */

	deter2_gmp(temp1,a33_mp,a43_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 9:	- Minor3(i,k,l,2,3,0) */

	deter3_gmp(temp1,a12_mp,a13_mp,a32_mp,a33_mp,a42_mp,a43_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  - icomp;
		return;
	}

/* Term 10:	Minor3(i,j,l,1,2,0) */

	deter3_gmp(temp1,a11_mp,a12_mp,a21_mp,a22_mp,a41_mp,a42_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 11: 	- Minor2(j,l,1,0) */

	deter2_gmp(temp1,a21_mp,a41_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  - icomp;
		return;
	}

/* Term 12:	Minor2(j,l,2,0) */

	deter2_gmp(temp1,a22_mp,a42_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 13:	Minor2(i,l,1,0) */

	deter2_gmp(temp1,a11_mp,a41_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 14:	1 */

	*res = 1;

}
#ifdef __cplusplus
}
#endif

/*									  */
/*------------------------------------------------------------------------*/
/* sos_minor5_gmp:							  */
/*									  
	This subroutine tests the sign of the determinant
		D = | a11 a12 a13 a14 1 |
		    | a21 a22 a23 a24 1 |
		    | a31 a32 a33 a34 1 |
		    | a41 a42 a43 a44 1 |
		    | a51 a52 a53 a54 1 |
	If the determinant is found to be 0, then the SoS procedure is used:
	a development of the determinant with *respect to a perturbation EPS
	applied to the coordinates in the determinant is computed, and
	the sign of the first non zero term defines the sign of the 
	determinant.
	In the case of a 5x5 determinant, the maximum number of terms to be
	checked is 49 ...					*/
#ifdef __cplusplus
extern "C" {
#endif
void F77Name2(sos_minor5_gmp)(int *a, int *b, int *c, int *d, int *e, int *res) 

{
	int icomp;

	int idx_a1,idx_a2,idx_a3,idx_a4,idx_b1,idx_b2,idx_b3,idx_b4;
	int idx_c1,idx_c2,idx_c3,idx_c4,idx_d1,idx_d2,idx_d3,idx_d4;
	int idx_e1,idx_e2,idx_e3,idx_e4;

	int ivalue;
	int base=10;

	double value;

/*	Initialise local GMP variables */

	idx_a1 = (*a)*3 -3; idx_a2 = idx_a1+1; idx_a3 = idx_a2 +1;
	idx_b1 = (*b)*3 -3; idx_b2 = idx_b1+1; idx_b3 = idx_b2 +1;
	idx_c1 = (*c)*3 -3; idx_c2 = idx_c1+1; idx_c3 = idx_c2 +1;
	idx_d1 = (*d)*3 -3; idx_d2 = idx_d1+1; idx_d3 = idx_d2 +1;
	idx_e1 = (*e)*3 -3; idx_e2 = idx_e1+1; idx_e3 = idx_e2 +1;

	idx_a4 = (*a) -1; idx_b4 = (*b)-1; idx_c4 = (*c)-1; 
	idx_d4 = (*d) -1; idx_e4 = (*e)-1;

	mpz_set(a11_mp,coord_gmp[idx_a1]);
	mpz_set(a12_mp,coord_gmp[idx_a2]);
	mpz_set(a13_mp,coord_gmp[idx_a3]);
	mpz_set(a21_mp,coord_gmp[idx_b1]);
	mpz_set(a22_mp,coord_gmp[idx_b2]);
	mpz_set(a23_mp,coord_gmp[idx_b3]);
	mpz_set(a31_mp,coord_gmp[idx_c1]);
	mpz_set(a32_mp,coord_gmp[idx_c2]);
	mpz_set(a33_mp,coord_gmp[idx_c3]);
	mpz_set(a41_mp,coord_gmp[idx_d1]);
	mpz_set(a42_mp,coord_gmp[idx_d2]);
	mpz_set(a43_mp,coord_gmp[idx_d3]);
	mpz_set(a51_mp,coord_gmp[idx_e1]);
	mpz_set(a52_mp,coord_gmp[idx_e2]);
	mpz_set(a53_mp,coord_gmp[idx_e3]);

	mpz_set(a14_mp,weight_gmp[idx_a4]);
	mpz_set(a24_mp,weight_gmp[idx_b4]);
	mpz_set(a34_mp,weight_gmp[idx_c4]);
	mpz_set(a44_mp,weight_gmp[idx_d4]);
	mpz_set(a54_mp,weight_gmp[idx_e4]);

/* Compute determinant */

	deter5_gmp(temp1,a11_mp,a12_mp,a13_mp,a14_mp,a21_mp,a22_mp,
	  a23_mp,a24_mp,a31_mp,a32_mp,a33_mp,a34_mp,a41_mp,a42_mp,
	   a43_mp,a44_mp,a51_mp,a52_mp,a53_mp,a54_mp);

//	printf("Deter5 = %s\n",mpz_get_str(NULL,base,temp1)); 

	icomp = mpz_sgn(temp1);

/* if major determinant is non 0, return its sign */

	if (icomp != 0) {
		*res = icomp;
		return;
	}

/* Look now at each term in the expansion of the determinant with
   *respect to EPS                                                
	The initial determinant is:
	Minor5(i,j,k,l,m,1,2,3,4,0)			*/

/* Term 1: 	-Minor4(j,k,l,m,1,2,3,0) */

	deter4_gmp(temp1,a21_mp,a22_mp,a23_mp,a31_mp,a32_mp,a33_mp,
		a41_mp,a42_mp,a43_mp,a51_mp,a52_mp,a53_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res = - icomp;
		return;
	}

/* Term 2:	Minor4(j,k,l,m,1,2,4,0) */
	
	deter4_gmp(temp1,a21_mp,a22_mp,a24_mp,a31_mp,a32_mp,a34_mp,
		a41_mp,a42_mp,a44_mp,a51_mp,a52_mp,a54_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 3:	- Minor4(j,k,l,m,1,3,4,0) */
	
	deter4_gmp(temp1,a21_mp,a23_mp,a24_mp,a31_mp,a33_mp,a34_mp,
		a41_mp,a43_mp,a44_mp,a51_mp,a53_mp,a54_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  - icomp;
		return;
	}

/* Term 4:	Minor4(j,k,l,m,2,3,4,0) */
	
	deter4_gmp(temp1,a22_mp,a23_mp,a24_mp,a32_mp,a33_mp,a34_mp,
		a42_mp,a43_mp,a44_mp,a52_mp,a53_mp,a54_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 5:	Minor4(i,k,l,m,1,2,3,0) */
	
	deter4_gmp(temp1,a11_mp,a12_mp,a13_mp,a31_mp,a32_mp,a33_mp,
		a41_mp,a42_mp,a43_mp,a51_mp,a52_mp,a53_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 6:	Minor3(k,l,m,1,2,0) */
	
	deter3_gmp(temp1,a31_mp,a32_mp,a41_mp,a42_mp,a51_mp,a52_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 7:	-Minor3(k,l,m,1,3,0) */
	
	deter3_gmp(temp1,a31_mp,a33_mp,a41_mp,a43_mp,a51_mp,a53_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  -icomp;
		return;
	}

/* Term 8:	Minor3(k,l,m,2,3,0) */
	
	deter3_gmp(temp1,a32_mp,a33_mp,a42_mp,a43_mp,a52_mp,a53_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 9:	-Minor4(i,k,l,m,1,2,4,0) */
	
	deter4_gmp(temp1,a11_mp,a12_mp,a14_mp,a31_mp,a32_mp,a34_mp,
		a41_mp,a42_mp,a44_mp,a51_mp,a52_mp,a54_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  -icomp;
		return;
	}

/* Term 10:	Minor3(k,l,m,1,4,0) */
	
	deter3_gmp(temp1,a31_mp,a34_mp,a41_mp,a44_mp,a51_mp,a54_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 11:	-Minor3(k,l,m,2,4,0) */
	
	deter3_gmp(temp1,a32_mp,a34_mp,a42_mp,a44_mp,a52_mp,a54_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  -icomp;
		return;
	}

/* Term 12:	Minor4(i,k,l,m,1,3,4,0) */
	
	deter4_gmp(temp1,a11_mp,a13_mp,a14_mp,a31_mp,a33_mp,a34_mp,
		a41_mp,a43_mp,a44_mp,a51_mp,a53_mp,a54_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 13:	Minor3(k,l,m,3,4,0) */
	
	deter3_gmp(temp1,a33_mp,a34_mp,a43_mp,a44_mp,a53_mp,a54_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 14:	-Minor4(i,k,l,m,2,3,4,0) */
	
	deter4_gmp(temp1,a12_mp,a13_mp,a14_mp,a32_mp,a33_mp,a34_mp,
		a42_mp,a43_mp,a44_mp,a52_mp,a53_mp,a54_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  -icomp;
		return;
	}

/* Term 15:	-Minor4(i,j,l,m,1,2,3,0) */
	
	deter4_gmp(temp1,a11_mp,a12_mp,a13_mp,a21_mp,a22_mp,a23_mp,
		a41_mp,a42_mp,a43_mp,a51_mp,a52_mp,a53_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  -icomp;
		return;
	}

/* Term 16:	-Minor3(j,l,m,1,2,0) */
	
	deter3_gmp(temp1,a21_mp,a22_mp,a41_mp,a42_mp,a51_mp,a52_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  -icomp;
		return;
	}

/* Term 17:	Minor3(j,l,m,1,3,0) */
	
	deter3_gmp(temp1,a21_mp,a23_mp,a41_mp,a43_mp,a51_mp,a53_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 18:	-Minor3(j,l,m,2,3,0) */
	
	deter3_gmp(temp1,a22_mp,a23_mp,a42_mp,a43_mp,a52_mp,a53_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  -icomp;
		return;
	}

/* Term 19:	Minor3(i,l,m,1,2,0) */
	
	deter3_gmp(temp1,a11_mp,a12_mp,a41_mp,a42_mp,a51_mp,a52_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 20:	-Minor2(l,m,1,0) */

	deter2_gmp(temp1,a41_mp,a51_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  -icomp;
		return;
	}

/* Term 21:	Minor2(l,m,2,0) */

	deter2_gmp(temp1,a42_mp,a52_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 22:	-Minor3(i,l,m,1,3,0) */
	
	deter3_gmp(temp1,a11_mp,a13_mp,a41_mp,a43_mp,a51_mp,a53_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  -icomp;
		return;
	}

/* Term 23:	-Minor2(l,m,3,0) */

	deter2_gmp(temp1,a43_mp,a53_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  -icomp;
		return;
	}

/* Term 24:	Minor3(i,l,m,2,3,0) */
	
	deter3_gmp(temp1,a12_mp,a13_mp,a42_mp,a43_mp,a52_mp,a53_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 25:	Minor4(i,j,l,m,1,2,4,0) */
	
	deter4_gmp(temp1,a11_mp,a12_mp,a14_mp,a21_mp,a22_mp,a24_mp,
		a41_mp,a42_mp,a44_mp,a51_mp,a52_mp,a54_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 26:	-Minor3(j,l,m,1,4,0) */
	
	deter3_gmp(temp1,a21_mp,a24_mp,a41_mp,a44_mp,a51_mp,a54_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  -icomp;
		return;
	}

/* Term 27:	Minor3(j,l,m,2,4,0) */
	
	deter3_gmp(temp1,a22_mp,a24_mp,a42_mp,a44_mp,a52_mp,a54_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 28:	Minor3(i,l,m,1,4,0) */
	
	deter3_gmp(temp1,a11_mp,a14_mp,a41_mp,a44_mp,a51_mp,a54_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 29:	Minor2(l,m,4,0) */

	deter2_gmp(temp1,a44_mp,a54_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 30:	-Minor3(i,l,m,2,4,0) */
	
	deter3_gmp(temp1,a12_mp,a14_mp,a42_mp,a44_mp,a52_mp,a54_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  -icomp;
		return;
	}

/* Term 31:	-Minor4(i,j,l,m,1,3,4,0) */
	
	deter4_gmp(temp1,a11_mp,a13_mp,a14_mp,a21_mp,a23_mp,a24_mp,
		a41_mp,a43_mp,a44_mp,a51_mp,a53_mp,a54_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  -icomp;
		return;
	}

/* Term 32:	-Minor3(j,l,m,3,4,0) */
	
	deter3_gmp(temp1,a23_mp,a24_mp,a43_mp,a44_mp,a53_mp,a54_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  -icomp;
		return;
	}

/* Term 33:	Minor3(i,l,m,3,4,0) */
	
	deter3_gmp(temp1,a13_mp,a14_mp,a43_mp,a44_mp,a53_mp,a54_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 34:	Minor4(i,j,l,m,2,3,4,0) */
	
	deter4_gmp(temp1,a12_mp,a13_mp,a14_mp,a22_mp,a23_mp,a24_mp,
		a42_mp,a43_mp,a44_mp,a52_mp,a53_mp,a54_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 35:	Minor4(i,j,k,m,1,2,3,0) */
	
	deter4_gmp(temp1,a11_mp,a12_mp,a13_mp,a21_mp,a22_mp,a23_mp,
		a31_mp,a32_mp,a33_mp,a51_mp,a52_mp,a53_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 36:	Minor3(j,k,m,1,2,0) */
	
	deter3_gmp(temp1,a21_mp,a22_mp,a31_mp,a32_mp,a51_mp,a52_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 37:	-Minor3(j,k,m,1,3,0) */
	
	deter3_gmp(temp1,a21_mp,a23_mp,a31_mp,a33_mp,a51_mp,a53_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  -icomp;
		return;
	}

/* Term 38:	Minor3(j,k,m,2,3,0) */
	
	deter3_gmp(temp1,a22_mp,a23_mp,a32_mp,a33_mp,a52_mp,a53_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 39:	-Minor3(i,k,m,1,2,0) */
	
	deter3_gmp(temp1,a11_mp,a12_mp,a31_mp,a32_mp,a51_mp,a52_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  -icomp;
		return;
	}

/* Term 40:	Minor2(k,m,1,0) */

	deter2_gmp(temp1,a31_mp,a51_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 41:	-Minor2(k,m,2,0) */

	deter2_gmp(temp1,a32_mp,a52_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  -icomp;
		return;
	}

/* Term 42:	Minor3(i,k,m,1,3,0) */
	
	deter3_gmp(temp1,a11_mp,a13_mp,a31_mp,a33_mp,a51_mp,a53_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 43:	Minor2(k,m,3,0) */

	deter2_gmp(temp1,a33_mp,a53_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 44:	-Minor3(i,k,m,2,3,0) */
	
	deter3_gmp(temp1,a12_mp,a13_mp,a32_mp,a33_mp,a52_mp,a53_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  -icomp;
		return;
	}

/* Term 45:	Minor3(i,j,m,1,2,0) */
	
	deter3_gmp(temp1,a11_mp,a12_mp,a21_mp,a22_mp,a51_mp,a52_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 46:	-Minor2(j,m,1,0) */

	deter2_gmp(temp1,a21_mp,a51_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  -icomp;
		return;
	}

/* Term 47:	Minor2(j,m,2,0) */

	deter2_gmp(temp1,a22_mp,a52_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 48:	Minor2(i,m,1,0) */

	deter2_gmp(temp1,a11_mp,a51_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 49:	1 */

	*res = 1;

}
#ifdef __cplusplus
}
#endif

/*------------------------------------------------------------------------*/
/* minor4_:								  */
/*									  
	This subroutine tests the sign of the determinant
		D = | a11 a12 a13 1 |
		    | a21 a22 a23 1 |
		    | a31 a32 a33 1 |
		    | a41 a42 a43 1 |
	and return 1 if positive, -1 if negative, 0 otherwise   */
#ifdef __cplusplus
extern "C" {
#endif
void F77Name2(minor4_gmp)(int *a, int *b, int *c, int *d, int *res)

{
	int icomp;
	int idx_a1,idx_a2,idx_a3,idx_b1,idx_b2,idx_b3;
	int idx_c1,idx_c2,idx_c3,idx_d1,idx_d2,idx_d3;
	int ivalue;

	double value;

/* Initialise local gmp variables */

/* Transfer coordinates to gmp */

	idx_a1 = (*a)*3 -3; idx_a2 = idx_a1+1; idx_a3 = idx_a2 +1;
	idx_b1 = (*b)*3 -3; idx_b2 = idx_b1+1; idx_b3 = idx_b2 +1;
	idx_c1 = (*c)*3 -3; idx_c2 = idx_c1+1; idx_c3 = idx_c2 +1;
	idx_d1 = (*d)*3 -3; idx_d2 = idx_d1+1; idx_d3 = idx_d2 +1;

	mpz_set(a11_mp,coord_gmp[idx_a1]);
	mpz_set(a12_mp,coord_gmp[idx_a2]);
	mpz_set(a13_mp,coord_gmp[idx_a3]);
	mpz_set(a21_mp,coord_gmp[idx_b1]);
	mpz_set(a22_mp,coord_gmp[idx_b2]);
	mpz_set(a23_mp,coord_gmp[idx_b3]);
	mpz_set(a31_mp,coord_gmp[idx_c1]);
	mpz_set(a32_mp,coord_gmp[idx_c2]);
	mpz_set(a33_mp,coord_gmp[idx_c3]);
	mpz_set(a41_mp,coord_gmp[idx_d1]);
	mpz_set(a42_mp,coord_gmp[idx_d2]);
	mpz_set(a43_mp,coord_gmp[idx_d3]);

/* Compute determinant */

	deter4_gmp(temp1,a11_mp,a12_mp,a13_mp,a21_mp,a22_mp,a23_mp,
			   a31_mp,a32_mp,a33_mp,a41_mp,a42_mp,a43_mp);

	icomp = mpz_sgn(temp1);

	*res = icomp;

}
#ifdef __cplusplus
}
#endif
