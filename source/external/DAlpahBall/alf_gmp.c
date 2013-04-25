/*	alf_gmp.c	Version 1 11/27/2000	Patrice Koehl             */
/*									  */
/*  This is the C version of alf.f, which performs all operations   */
/*  with multi precision arithmetics, using the package GMP               */
/*									  */
/*------------------------------------------------------------------------*/
/*									  */
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
/*									  */
/*------------------------------------------------------------------------*/
/*                                                                        */
/* Includes :								  */

#include <stdio.h>
#include <stdlib.h>
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

/*************************************************************************************/



/*------------------------------------------------------------------------*/
/*Local procedures:							  */

void vertex_attach_gmp(mpz_t *Dab, mpz_t ra_mp, mpz_t rb_mp, int *testa, int *testb);

void edge_attach_gmp(mpz_t *m_mp, mpz_t *n_mp, mpz_t *Dmn, mpz_t *Smn, 
		mpz_t *Sp, mpz_t *Tp, int *testa);

void edge_radius_gmp(mpz_t *m_mp, mpz_t *n_np, mpz_t wm, mpz_t wn,
		mpz_t *Dmn, mpz_t *Smn, int *testr, mpz_t alp);

void triangle_attach_gmp(mpz_t *S, mpz_t De1, mpz_t De2, mpz_t De3, mpz_t Dmnpq, 
		mpz_t Dmnp, int *testa);

void triangle_radius_gmp(mpz_t *S, mpz_t *T, mpz_t *U, mpz_t Dmnp,
		int *testr,  mpz_t alp);
#ifdef __cplusplus
extern "C"
#endif
void F77Name2(alf_gmp)( int *ia, int *ib, int *ic, int *id,
			int *testu, int *testr, int *testa, 
		int *trig_stat, int *edge_stat, double *scale, double *alpha);

void F77Name2(set_alf_gmp)();
void clear_alf_gmp();


/*------------------------------------------------------------------------*/
/* vertex_attach_gmp:
	This subroutine checks if a vertex is attached to another vertex 
	Input:
		a, b	: coordinates of the two points
		ra,rb	: radii of the two points
	Output:
		testa	: flag equal to 1 if a is attached to b
		testb	: flag equal to 1 if b is attached to a
*/
void	vertex_attach_gmp(mpz_t *Dab, mpz_t ra_mp, mpz_t rb_mp, int *testa, int *testb)

{

	(*testa = 0);
	(*testb = 0);

	mpz_mul(temp1, Dab[1],Dab[1]);
	mpz_mul(temp2, Dab[2],Dab[2]);
	mpz_mul(temp3, Dab[3],Dab[3]);
	mpz_add(temp1, temp1,temp2); mpz_add(dist2,temp1,temp3);

//	printf("Dist2 = %s\n",mpz_get_str(NULL,10,dist2));
	mpz_mul(ra2, ra_mp, ra_mp); mpz_mul(rb2, rb_mp, rb_mp);
//	printf("ra2 = %s\n",mpz_get_str(NULL,10,ra2));
//	printf("rb2 = %s\n",mpz_get_str(NULL,10,rb2));

	mpz_add(dtest,dist2,ra2);
	mpz_sub(dtest,dtest,rb2);
//	printf("Dtesta = %s\n",mpz_get_str(NULL,10,dtest));
	if(mpz_sgn(dtest) < 0) (*testa = 1);

	mpz_sub(dtest,dist2,ra2);
	mpz_add(dtest,dtest,rb2);
//	printf("Dtestb = %s\n",mpz_get_str(NULL,10,dtest));
	if(mpz_sgn(dtest) < 0) (*testb = 1);

}

/*------------------------------------------------------------------------*/
/* edge_attach_gmp:
	This subroutine checks if an edge of the regular triangulation is "attached"
	to another vertex (i.e.  if the vertex belongs to the smallest circumsphere
	of the edge). For that, it needs:

	Input:
		Dab	: minor(a,b,i,0) for all i=1,2,3
		Sab	: minor(a,b,i,j) for i = 1,2 and j =i+1,3
		Sc	: minor(a,b,c,i,j,0) for i=1,2 and j = i+1,3
			  c is the other vertex
		Tc	: minor(a,b,c,i,4,0) for i = 1,2,3
	Ouput:
		testa	: flag that defines if edge is attached or not
									  
For comments, see alfcx.f which contains the fortran equivalence of this
routine
									*/

void edge_attach_gmp(mpz_t *m_mp, mpz_t *n_mp, mpz_t *Dab, mpz_t *Sab, 
	mpz_t *Sc,mpz_t *Tc, int *testa)

{
	int i,j,coef;

/*	This is the "hidden1" part */

	(*testa) = 0;

	if( mpz_cmp(m_mp[1],n_mp[1]) != 0) 
	{
		for (i = 1; i < 4 ; i++)
		{
			mpz_set(res[0][i],Dab[i]);
			mpz_set(res2_c[i][4],Tc[i]);
		}
		mpz_set(res[1][2],Sab[1]); mpz_set(res[1][3],Sab[2]);
		mpz_set(res[2][3],Sab[3]);
		mpz_set(res2_c[1][2],Sc[1]); mpz_set(res2_c[1][3],Sc[2]);
		mpz_set(res2_c[2][3],Sc[3]);
	}
	else if ( mpz_cmp(m_mp[2],n_mp[2]) != 0)
	{
		mpz_set(res[0][1],Dab[2]); mpz_set(res[0][2],Dab[3]);
		mpz_set(res[0][3],Dab[1]);
		mpz_set(res[1][2],Sab[3]);
		mpz_neg(res[1][3],Sab[1]); mpz_neg(res[2][3],Sab[2]);
		mpz_set(res2_c[1][2],Sc[3]);
		mpz_neg(res2_c[1][3],Sc[1]); mpz_neg(res2_c[2][3],Sc[2]);
		mpz_set(res2_c[1][4],Tc[2]); mpz_set(res2_c[2][4],Tc[3]);
		mpz_set(res2_c[3][4],Tc[1]);
	}
	else if (  mpz_cmp(m_mp[3],n_mp[3]) != 0)
	{
		mpz_set(res[0][1],Dab[3]); mpz_set(res[0][2],Dab[1]);
		mpz_set(res[0][3],Dab[2]);
		mpz_neg(res[1][2],Sab[2]);
		mpz_neg(res[1][3],Sab[3]); mpz_set(res[2][3],Sab[1]);
		mpz_neg(res2_c[1][2],Sc[2]);
		mpz_neg(res2_c[1][3],Sc[3]); mpz_set(res2_c[2][3],Sc[1]);
		mpz_set(res2_c[1][4],Tc[3]); mpz_set(res2_c[2][4],Tc[1]);
		mpz_set(res2_c[3][4],Tc[2]);
	}
	else
	{
		mpz_clear(temp1); mpz_clear(temp2); mpz_clear(temp3);
		mpz_clear (r_11); mpz_clear (r_22); mpz_clear (r_33);
		mpz_clear (diff);
		mpz_clear (det0); mpz_clear (dtest);
		for (i= 0; i < 4; i++)
		{
			for (j=0; j < 5; j++)
			{
				mpz_clear(res[i][j]);
				mpz_clear(res2_c[i][j]);
			}
		}

		exit(1);
	}

	mpz_mul(r_11,res[0][1],res[0][1]);
	mpz_mul(r_22,res[0][2],res[0][2]);
	mpz_mul(r_33,res[0][3],res[0][3]);
	mpz_mul(temp1,res[0][3],res[1][2]); 
	mpz_mul(temp2,res[0][2],res[1][3]);
	mpz_sub(diff,temp1,temp2);

/* Compute det0 */

	mpz_add(temp1,r_22,r_33); mpz_add(temp1,temp1,r_11); 
	mpz_mul(temp1,temp1,res[0][1]); 
	coef = -2;
	mpz_mul_si(det0,temp1,coef);


/* Now check if edge (ab) is attached to c */

	mpz_mul(temp1,res[1][2],res2_c[1][2]);
	mpz_mul(temp2,res[1][3],res2_c[1][3]);
	mpz_add(temp1,temp1,temp2);
	coef = -2;
	mpz_mul_si(temp1,temp1,coef);

	mpz_set_si(temp2,0);
	for (i=1; i<4; i++)
	{
		mpz_mul(temp3,res[0][i],res2_c[i][4]);
		mpz_add(temp2,temp2,temp3);
	}
	mpz_add(temp1,temp2,temp1);
	mpz_mul(temp1,temp1,res[0][1]);

	mpz_mul(temp2,res2_c[2][3],diff);
	mpz_mul_si(temp2,temp2,coef);
	mpz_sub(temp3,temp1,temp2);
	mpz_mul(dtest,temp3,det0);
//	printf("Dtest (edge attach) = %s\n",mpz_get_str(NULL,10,dtest));

	if(mpz_sgn(dtest) < 0) (*testa = 1);

}

/*------------------------------------------------------------------------*/
/* edge_radius_gmp:
	This subroutine checks if the radius of the smallest circumsphere of an edge
	of the regular triangulation is smaller than the value of alpha
	For that, it needs:

	Input:
		a,b	: coordinates of the two vertices defining the edge
		wa,wb	: weights of these two vertices
		Dab	: minor(a,b,i,0) for all i=1,2,3
		Sab	: minor(a,b,i,j) for i = 1,2 and j =i+1,3
	Ouput:
		testr	: flag that defines if radius < alpha or not
									  
For comments, see alfcx.f which contains the fortran equivalence of this
routine                                                                           */

void edge_radius_gmp(mpz_t *m_mp, mpz_t *n_mp, mpz_t wm, mpz_t wn,
		mpz_t *Dab, mpz_t *Sab, int *testr, mpz_t alp)

{
	int i,j,coef;

/*	This is the "hidden1" part */

	(*testr) = 0;
	mpz_sub(res[0][4],wm,wn);


	if( mpz_cmp(m_mp[1],n_mp[1]) != 0) 
	{
		for (i = 1; i < 4 ; i++)
		{
			mpz_set(res[0][i],Dab[i]);
			mpz_mul(temp1,n_mp[i],wm);
			mpz_mul(temp2,m_mp[i],wn);
			mpz_sub(res[i][4],temp2,temp1);
		}
		mpz_set(res[1][2],Sab[1]); mpz_set(res[1][3],Sab[2]);
		mpz_set(res[2][3],Sab[3]);
	}
	else if ( mpz_cmp(m_mp[2],n_mp[2]) != 0)
	{
		mpz_set(res[0][1],Dab[2]); mpz_set(res[0][2],Dab[3]);
		mpz_set(res[0][3],Dab[1]);
		mpz_set(res[1][2],Sab[3]);
		mpz_neg(res[1][3],Sab[1]); mpz_neg(res[2][3],Sab[2]);
		mpz_mul(temp1,m_mp[2],wn); mpz_mul(temp2,n_mp[2],wm);
		mpz_sub(res[1][4],temp1,temp2);
		mpz_mul(temp1,m_mp[3],wn); mpz_mul(temp2,n_mp[3],wm);
		mpz_sub(res[2][4],temp1,temp2);
		mpz_mul(temp1,m_mp[1],wn); mpz_mul(temp1,n_mp[1],wm);
		mpz_sub(res[3][4],temp1,temp2);
	}
	else if (  mpz_cmp(m_mp[3],n_mp[3]) != 0)
	{
		mpz_set(res[0][1],Dab[3]); mpz_set(res[0][2],Dab[1]);
		mpz_set(res[0][3],Dab[2]);
		mpz_neg(res[1][2],Sab[2]);
		mpz_neg(res[1][3],Sab[3]); mpz_set(res[2][3],Sab[1]);
		mpz_mul(temp1,m_mp[3],wn); mpz_mul(temp2,n_mp[3],wm);
		mpz_sub(res[1][4],temp1,temp2);
		mpz_mul(temp1,m_mp[1],wn); mpz_mul(temp2,n_mp[1],wm);
		mpz_sub(res[2][4],temp1,temp2);
		mpz_mul(temp1,m_mp[2],wn); mpz_mul(temp1,n_mp[2],wm);
		mpz_sub(res[3][4],temp1,temp2);
	}
	else
	{
		mpz_clear(temp1); mpz_clear(temp2); mpz_clear(temp3);
		mpz_clear (r_11); mpz_clear (r_22); mpz_clear (r_33); mpz_clear (r_14);
		mpz_clear (r_313); mpz_clear (r_212); mpz_clear (diff);
		mpz_clear (det0); mpz_clear (det1); mpz_clear (det2); mpz_clear (det3);
		mpz_clear (det4); 
		mpz_clear (num); mpz_clear (den); mpz_clear (dtest);
		for (i= 0; i < 4; i++)
		{
			for (j=0; j < 5; j++)
			{
				mpz_clear(res[i][j]);
			}
		}

		exit(1);
	}

	mpz_mul(r_11,res[0][1],res[0][1]);
	mpz_mul(r_22,res[0][2],res[0][2]);
	mpz_mul(r_33,res[0][3],res[0][3]);
	mpz_mul(r_14,res[0][1],res[0][4]);
	mpz_mul(r_313,res[0][3],res[1][3]);
	mpz_mul(r_212,res[0][2],res[1][2]);
	mpz_mul(temp1,res[0][3],res[1][2]); 
	mpz_mul(temp2,res[0][2],res[1][3]);
	mpz_sub(diff,temp1,temp2);

/* Compute det0 */

	mpz_add(temp1,r_22,r_33); mpz_add(temp1,temp1,r_11); 
	mpz_mul(temp1,temp1,res[0][1]); 
	coef = -2;
	mpz_mul_si(det0,temp1,coef);

/* Compute det1 */

	mpz_add(temp1,r_313,r_212);
	coef = 2;
	mpz_mul_si(temp1,temp1,coef);
	mpz_sub(temp1,temp1,r_14);
	mpz_mul(det1,res[0][1],temp1);


/* Compute det2 */

	mpz_add(temp1,r_11,r_33);
	mpz_mul(temp1,temp1,res[1][2]);
	coef = -2;
	mpz_mul_si(temp1,temp1,coef);
	mpz_mul_si(temp2,r_313,coef);
	mpz_add(temp2,temp2,r_14);
	mpz_mul(temp2,temp2,res[0][2]);
	mpz_sub(det2,temp1,temp2);

/* Compute det3 */

	mpz_add(temp1,r_11,r_22);
	mpz_mul(temp1,temp1,res[1][3]);
	mpz_mul_si(temp1,temp1,coef);
	mpz_mul_si(temp2,r_212,coef);
	mpz_add(temp2,temp2,r_14);
	mpz_mul(temp2,temp2,res[0][3]);
	mpz_sub(det3,temp1,temp2);

/* Compute det4 */

	mpz_mul(temp1,res[0][3],res[3][4]);
	mpz_mul(temp2,res[0][2],res[2][4]);
	mpz_mul(temp3,res[0][1],res[1][4]);
	mpz_add(temp1,temp1,temp2); mpz_add(temp1,temp3,temp1);
	coef = 2;
	mpz_mul_si(temp1,temp1,coef);
	mpz_mul(temp1,temp1,res[0][1]);
	mpz_mul(temp2,res[1][3],res[1][3]);
	mpz_mul(temp3,res[1][2],res[1][2]);
	mpz_add(temp2,temp3,temp2);
	mpz_mul(temp2,temp2,res[0][1]);
	mpz_mul(temp3,res[2][3],diff);
	mpz_sub(temp2,temp3,temp2);
	coef = 4;
	mpz_mul_si(temp2,temp2,coef);
	mpz_add(det4,temp1,temp2);

/* Compute numerator of the radius of the smallest circumsphere of the edge */

	mpz_mul(temp1,det0,det4);
	mpz_mul(temp2,det3,det3);
	mpz_sub(temp2,temp2,temp1);
	mpz_mul(temp1,det2,det2);
	mpz_add(temp2,temp2,temp1);
	mpz_mul(temp1,det1,det1);
	mpz_add(num,temp1,temp2);
//	printf("Num (edge radius) = %s\n",mpz_get_str(NULL,10,num));

/* Compute denominator of the radius of the smallest circumsphere of the edge */

	mpz_mul(den,det0,det0);
//	printf("Den (edge radius) = %s\n",mpz_get_str(NULL,10,den));

/* check if radius is lower than ALPHA         */

	mpz_mul(temp1,den,alp);
	mpz_sub(temp2,num,temp1);

	if(mpz_sgn(temp2) < 0) (*testr)=1;

}

/*------------------------------------------------------------------------*/
/* triangle_radius_gmp:
   This program checks if the radius of the circumsphere of a facet of the
   regular triangulation is smaller than alpha

	Input:

	For the three points m,n,p that form the triangles, the program
	needs as input the following determinants:

	S(i+j-2) = Minor(m,n,p,i,j,0)= det | m(i)  m(j)  1 |
					   | n(i)  n(j)  1 |
					   | p(i)  p(j)  1 |

	T(i) = Minor(m,n,p,i,4,0) = det | m(i)  m(4)  1 |
					| n(i)  n(4)  1 |
					| p(i)  p(4)  1 |

	U(i) = Minor(m,n,p,i,j,4) = det | m(i) m(j) m(4) |
					| n(i) n(j) n(4) |
					| p(i) p(j) p(4) |

	Dmnp  = Minor(m,n,p,1,2,3)

	Output:

	testr	: flag set to 1 if ALPHA is larger than rho, the radius
		  of the circumsphere of the triangle

*/

void triangle_radius_gmp(mpz_t *S, mpz_t *T, mpz_t *U,
		 mpz_t Dmnp, int *testr, mpz_t alp)

{
	int i,coef;


	(*testr) = 0;

	mpz_set_si(temp1,0);
	for (i=1; i<4; i++)
	{
		mpz_mul(temp2,S[i],S[i]);
		mpz_add(temp1,temp1,temp2);
	}

/* Compute det0 */

	coef = 4;
	mpz_mul_si(det0,temp1,coef);

/* Compute det1 */

	mpz_mul(temp1,Dmnp,S[3]);
	coef = -2;
	mpz_mul_si(temp1,temp1,coef);
	mpz_mul(temp2,S[1],T[2]);
	mpz_add(temp1,temp2,temp1);
	mpz_mul(temp2,S[2],T[3]);
	mpz_add(temp1,temp2,temp1);
	mpz_mul_si(det1,temp1,coef);

/* Compute det2 */

	coef = 2;
	mpz_mul(temp1,Dmnp,S[2]);
	mpz_mul_si(temp1,temp1,coef);
	mpz_mul(temp2,S[3],T[3]);
	mpz_add(temp1,temp2,temp1);
	mpz_mul(temp2,S[1],T[1]);
	mpz_sub(temp1,temp2,temp1);
	mpz_mul_si(det2,temp1,coef);

/* Compute det3 */

	mpz_mul(temp1,Dmnp,S[1]);
	mpz_mul_si(temp1,temp1,coef);
	mpz_mul(temp2,S[2],T[1]);
	mpz_add(temp1,temp2,temp1);
	mpz_mul(temp2,S[3],T[2]);
	mpz_add(temp1,temp2,temp1);
	mpz_mul_si(det3,temp1,coef);

/* Compute det4 */

	mpz_mul(temp1,Dmnp,Dmnp);
	coef = -2;
	mpz_mul_si(temp1,temp1,coef);

	for (i=1; i<4; i++)
	{
		mpz_mul(temp2,S[i],U[i]);
		mpz_add(temp1,temp1,temp2);
	}
	coef = -4;
	mpz_mul_si(det4,temp1,coef);

/* Now compute numerator of the radius of the circumsphere of the triangle */

	mpz_mul(temp1,det0,det4);
	mpz_mul(temp2,det3,det3);
	mpz_sub(temp2,temp2,temp1);
	mpz_mul(temp1,det2,det2);
	mpz_add(temp2,temp2,temp1);
	mpz_mul(temp1,det1,det1);
	mpz_add(num,temp1,temp2);
//	printf("Num (triangle radius) = %s\n",mpz_get_str(NULL,10,num));

/* Now compute denominator of the radius of the circumsphere of the triangle */

	mpz_mul(den,det0,det0);
//	printf("Den (triangle radius) = %s\n",mpz_get_str(NULL,10,den));

/* Check if radius is lower than ALPHA */

	mpz_mul(temp1,den,alp);
	mpz_sub(temp2,num,temp1);

	if(mpz_sgn(temp2) < 0) (*testr) = 1;

}
/*------------------------------------------------------------------------*/
/* triangle_attach_gmp:
   This program checks if a facet is attached to a vertex

	Input:

	For the three points m,n,p that form the triangles, the program
	needs as input the following determinants:

	S(i+j-2) = Minor(m,n,p,i,j,0)= det | m(i)  m(j)  1 |
					   | n(i)  n(j)  1 |
					   | p(i)  p(j)  1 |

If q is the fourth point of the tetrahedron,

	De1 = Minor(m,n,p,q,2,3,4,0)
	De2 = Minor(m,n,p,q,1,3,4,0)
	De3 = Minor(m,n,p,q,1,2,4,0)
	Dmnpq= Minor(m,n,p,q,1,2,3,0)
	Dmnp  = Minor(m,n,p,1,2,3)

	Output:

	testa	: flag set to 1 if the fourth point d is inside the
		  circumsphere of (a,b,c)

*/

void triangle_attach_gmp(mpz_t *S, mpz_t De1, mpz_t De2, mpz_t De3, 
		mpz_t Dmnpq, mpz_t Dmnp, int *testa)

{
	int i,coef;


/* First check if triangle is attached                               */

	(*testa) = 0;

	mpz_set_si(temp1,0);
	for (i=1; i<4; i++)
	{
		mpz_mul(temp2,S[i],S[i]);
		mpz_add(temp1,temp1,temp2);
	}

	mpz_mul(temp2,Dmnpq,Dmnp);
	coef = -2;
	mpz_mul_si(temp2,temp2,coef);
	mpz_mul(temp3,De3,S[1]);
	mpz_add(temp2,temp3,temp2);
	mpz_mul(temp3,De2,S[2]);
	mpz_add(temp2,temp3,temp2);
	mpz_mul(temp3,De1,S[3]);
	mpz_add(temp2,temp3,temp2);
	mpz_mul(dtest,temp1,temp2);
//	printf("Dtest (triangle attach) = %s\n",mpz_get_str(NULL,10,dtest));

	if(mpz_sgn(dtest) > 0) (*testa)=1;

/* 	Clear local GMP variables */

}
/*------------------------------------------------------------------------*/
/* alf_gmp	Version 1 11/24/2000	Patrice Koehl
 
 	This subroutine computes the radius R of the circumsphere containing
 	a tetrahedron [A,B,C,D], as well as check if any fourth point L 
 	(in A, B, C, D) of the tetrahedron is "hidden" by its opposite 
 	face [I,J,K], i.e. is interior to the cicumsphere of [I,J,K]
 
 	Since we are only interested at how R compares to Alpha, we don't
 	output R, rather the result of the comparison
 
 	Computation extends to all four faces of the tetrehedron, as
 	well as to all six edges of the tetrahedron
 
	This procedure works with Multiple Precision Integer Arithmetics  (MPIA)
 	The package GMP is used for MPIA (with a C wrapper)
*/
#ifdef __cplusplus
extern "C" {
#endif
void F77Name2(alf_gmp)( int *a, int *b, int *c, int *d,
		int *testu, int *testr, int *testa, int *trig_stat, int *edge_stat, 
		double *scale, double *alpha)

{
	int i,j,k,coef,test_u,test_r,test_a,test_a2;
	int idx_a1,idx_a2,idx_b1,idx_b2;
	int idx_c1,idx_c2,idx_d1,idx_d2;
	int ivalue;
	double value;

/*	Transfer data in multiple precision */

	idx_a1 = (*a)*3 -3; 
	idx_b1 = (*b)*3 -3;
	idx_c1 = (*c)*3 -3;
	idx_d1 = (*d)*3 -3;

	for (i=0; i<4; i++)
	{
		mpz_set(a_mp[i+1],coord_gmp[idx_a1+i]);
		mpz_set(b_mp[i+1],coord_gmp[idx_b1+i]);
		mpz_set(c_mp[i+1],coord_gmp[idx_c1+i]);
		mpz_set(d_mp[i+1],coord_gmp[idx_d1+i]);
	}
	
	idx_a2 = (*a)-1; idx_b2 = (*b)-1; idx_c2 = (*c)-1; idx_d2 = (*d)-1;

	mpz_set(wa,weight_gmp[idx_a2]);
	mpz_set(wb,weight_gmp[idx_b2]);
	mpz_set(wc,weight_gmp[idx_c2]);
	mpz_set(wd,weight_gmp[idx_d2]);
		
	mpz_set(ra_mp,radius_gmp[idx_a2]);
	mpz_set(rb_mp,radius_gmp[idx_b2]);
	mpz_set(rc_mp,radius_gmp[idx_c2]);
	mpz_set(rd_mp,radius_gmp[idx_d2]);

	value = (*alpha)*(*scale); ivalue = (int) floor(value); 
	mpz_set_si(alp,ivalue);

/*	1. Computes all Minors Smn(i+j-2)= M(m,n,i,j) = Det | m(i)  m(j) |
						            | n(i)  n(j) |
	for all i in [1,2] and all j in [i+1,3]                         */

	for (i=1;  i<3; i++)
	{
		for (j=i+1; j<4 ; j++)
		{
			k=i+j-2;
			mpz_mul(temp1,a_mp[j],b_mp[i]); 
			mpz_mul(temp2,a_mp[i],b_mp[j]);
			mpz_sub(Sab[k],temp2,temp1);
			mpz_mul(temp1,a_mp[j],c_mp[i]); 
			mpz_mul(temp2,a_mp[i],c_mp[j]);
			mpz_sub(Sac[k],temp2,temp1);
			mpz_mul(temp1,a_mp[j],d_mp[i]); 
			mpz_mul(temp2,a_mp[i],d_mp[j]);
			mpz_sub(Sad[k],temp2,temp1);
			mpz_mul(temp1,b_mp[j],c_mp[i]); 
			mpz_mul(temp2,b_mp[i],c_mp[j]);
			mpz_sub(Sbc[k],temp2,temp1);
			mpz_mul(temp1,b_mp[j],d_mp[i]); 
			mpz_mul(temp2,b_mp[i],d_mp[j]);
			mpz_sub(Sbd[k],temp2,temp1);
			mpz_mul(temp1,c_mp[j],d_mp[i]); 
			mpz_mul(temp2,c_mp[i],d_mp[j]);
			mpz_sub(Scd[k],temp2,temp1);
		}
	}

/*	Now compute all Minors 
		Sq(i+j-2) = M(m,n,p,i,j,0) = Det | m(i) m(j) 1 |
		       			         | n(i) n(j) 1 |
						 | p(i) p(j) 1 |

	and all Minors
		Det(i+j-2) = M(m,n,p,q,i,j,4,0) = Det | m(i) m(j) m(4) 1 |
						      | n(i) n(j) n(4) 1 |
						      | p(i) p(j) p(4) 1 |
						      | q(i) q(j) q(4) 1 |

	m,n,p,q are the four vertices of the tetrahedron, i and j correspond
	to two of the coordinates of the vertices, and m(4) refers to the
	"weight" of vertices m                                           */
 
	for (i=1; i<4; i++)
	{
		mpz_sub(temp1,Scd[i],Sbd[i]); mpz_add(Sa[i],temp1,Sbc[i]);
		mpz_mul(temp2,Sa[i],wa);
		mpz_sub(temp1,Scd[i],Sad[i]); mpz_add(Sb[i],temp1,Sac[i]);
		mpz_mul(temp3,Sb[i],wb); mpz_sub(temp2,temp2,temp3);
		mpz_sub(temp1,Sbd[i],Sad[i]); mpz_add(Sc[i],temp1,Sab[i]);
		mpz_mul(temp3,Sc[i],wc); mpz_add(temp2,temp2,temp3);
		mpz_sub(temp1,Sbc[i],Sac[i]); mpz_add(Sd[i],temp1,Sab[i]);
		mpz_mul(temp3,Sd[i],wd); mpz_sub(Deter[i],temp2,temp3);
		mpz_neg(Sam1[i],Sa[i]); mpz_neg(Sbm1[i],Sb[i]);
		mpz_neg(Scm1[i],Sc[i]); mpz_neg(Sdm1[i],Sd[i]);
	}
 
/*
	Now compute the determinant needed to compute the radius of the
	circumsphere of the tetrahedron :

		Det1 = Minor(a,b,c,d,4,2,3,0)
		Det2 = Minor(a,b,c,d,1,3,4,0)
		Det3 = Minor(a,b,c,d,1,2,4,0)
		Det4 = Minor(a,b,c,d,1,2,3,0)
									*/

	mpz_set(Det1,Deter[3]);
	mpz_set(Det2,Deter[2]);
	mpz_set(Det3,Deter[1]);

	mpz_mul(temp1,a_mp[1],Sa[3]);mpz_mul(temp2,b_mp[1],Sb[3]);
	mpz_sub(temp3,temp1,temp2);
	mpz_mul(temp1,c_mp[1],Sc[3]);mpz_mul(temp2,d_mp[1],Sd[3]);
	mpz_sub(temp1,temp1,temp2);
	mpz_add(Det4,temp1,temp3);

/*
	Now compute all minors:
		Dmnp = Minor(m,n,p,1,2,3) = Det | m(1) m(2) m(3) |
						| n(1) n(2) n(3) |
						| p(1) p(2) p(3) |
									*/

	mpz_mul(temp1,a_mp[1],Sbc[3]); mpz_mul(temp2,b_mp[1],Sac[3]);
	mpz_sub(temp3,temp1,temp2);
	mpz_mul(temp1,c_mp[1],Sab[3]);mpz_add(Dabc,temp3,temp1);

	mpz_mul(temp1,a_mp[1],Sbd[3]); mpz_mul(temp2,b_mp[1],Sad[3]);
	mpz_sub(temp3,temp1,temp2);
	mpz_mul(temp1,d_mp[1],Sab[3]);mpz_add(Dabd,temp3,temp1);

	mpz_mul(temp1,a_mp[1],Scd[3]); mpz_mul(temp2,c_mp[1],Sad[3]);
	mpz_sub(temp3,temp1,temp2);
	mpz_mul(temp1,d_mp[1],Sac[3]);mpz_add(Dacd,temp3,temp1);

	mpz_mul(temp1,b_mp[1],Scd[3]); mpz_mul(temp2,c_mp[1],Sbd[3]);
	mpz_sub(temp3,temp1,temp2);
	mpz_mul(temp1,d_mp[1],Sbc[3]);mpz_add(Dbcd,temp3,temp1);


/*
	We also need :
		Det = Det | m(1) m(2) m(3) m(4) |
			  | n(1) n(2) n(3) n(4) |
			  | p(1) p(2) p(3) p(4) |
			  | q(1) q(2) q(3) q(4) |
								*/

	mpz_mul(temp1,wa,Dbcd); mpz_mul(temp2,wb,Dacd);
	mpz_sub(temp3,temp2,temp1);
	mpz_mul(temp1,wc,Dabd); mpz_mul(temp2,wd,Dabc);
	mpz_sub(temp1,temp2,temp1); mpz_add(Dabcd,temp3,temp1);

/*
	The radius of the circumsphere of the weighted tetrahedron is then:
	r_t = (Det1*Det1 + Det2*Det2 + Det3*Det3 + 4*Det4*Dabcd)/(4*Det4*Det4)
								*/

	mpz_mul(temp1,Det4,Det4); coef=4; mpz_mul_si(den,temp1,coef);

	mpz_mul(temp1,Det1,Det1); mpz_mul(temp2,Det2,Det2);
	mpz_add(temp1,temp1,temp2); mpz_mul(temp2,Det3,Det3);
	mpz_add(temp1,temp1,temp2); mpz_mul(temp2,Det4,Dabcd);
	mpz_mul_si(temp2,temp2,coef); mpz_add(num,temp1,temp2);
	
	mpz_mul(temp1,den,alp); mpz_sub(temp2,num,temp1);
//	printf("Dtest = %s\n",mpz_get_str(NULL,10,temp2));
	
/* 
 	If tetrahedron is part of the alpha shape, then the 4 triangles,
 	the 6 edges and the four vertices are also part of the alpha
 	complex
 								*/
	if(!(mpz_sgn(temp2) > 0)) 
//	if((mpz_sgn(temp2) < 0)) 
	{
		for (i=0; i<15;i++) 
		{
			testu[i]=1;
			testr[i]=1;
		}
	}
	else
	{
/* 
 	Now check all four faces of the tetrahedra:
 	We look both if the fourth vertex is "hidden" by the face of 
 	interest, and then we compute the radius of the circumsphere
 	of the face
 
 	We first need the minors:
 
 	Tm(i) = Minor(n,p,q,i,4,0) = det | n(i)  n(4)  1 |
 					 | p(i)  p(4)  1 |
 					 | q(i)  q(4)  1 |
 
 	Um(i) = Minor(n,p,q,i,j,4) = det | n(i) n(j) n(4) |
 					 | p(i) p(j) p(4) |
 					 | q(i) q(j) q(4) |
 									*/

		for (i=1; i<4; i++)
		{
			mpz_sub(Dab[i],a_mp[i],b_mp[i]);
			mpz_sub(Dac[i],a_mp[i],c_mp[i]);
			mpz_sub(Dad[i],a_mp[i],d_mp[i]);
			mpz_sub(Dbc[i],b_mp[i],c_mp[i]);
			mpz_sub(Dbd[i],b_mp[i],d_mp[i]);
			mpz_sub(Dcd[i],c_mp[i],d_mp[i]);
		}
		for (i=1; i<4; i++)
		{
			mpz_mul(temp1,wb,Dcd[i]); mpz_mul(temp2,wd,Dbc[i]);
			mpz_mul(temp3,wc,Dbd[i]); mpz_add(temp1,temp1,temp2);
			mpz_sub(Ta[i],temp3,temp1); 
			mpz_neg(Tam1[i],Ta[i]);
			mpz_mul(temp1,wa,Dcd[i]); mpz_mul(temp2,wd,Dac[i]);
			mpz_mul(temp3,wc,Dad[i]); mpz_add(temp1,temp1,temp2);
			mpz_sub(Tb[i],temp3,temp1); 
			mpz_neg(Tbm1[i],Tb[i]);
			mpz_mul(temp1,wa,Dbd[i]); mpz_mul(temp2,wd,Dab[i]);
			mpz_mul(temp3,wb,Dad[i]); mpz_add(temp1,temp1,temp2);
			mpz_sub(Tc[i],temp3,temp1); 
			mpz_neg(Tcm1[i],Tc[i]);
			mpz_mul(temp1,wa,Dbc[i]); mpz_mul(temp2,wc,Dab[i]);
			mpz_mul(temp3,wb,Dac[i]); mpz_add(temp1,temp1,temp2);
			mpz_sub(Td[i],temp3,temp1); 
			mpz_neg(Tdm1[i],Td[i]);
			mpz_mul(temp1,wb,Scd[i]); mpz_mul(temp2,wd,Sbc[i]);
			mpz_mul(temp3,wc,Sbd[i]); mpz_add(temp1,temp1,temp2);
			mpz_sub(Ua[i],temp1,temp3);
			mpz_mul(temp1,wa,Scd[i]); mpz_mul(temp2,wd,Sac[i]);
			mpz_mul(temp3,wc,Sad[i]); mpz_add(temp1,temp1,temp2);
			mpz_sub(Ub[i],temp1,temp3);
			mpz_mul(temp1,wa,Sbd[i]); mpz_mul(temp2,wd,Sab[i]);
			mpz_mul(temp3,wb,Sad[i]); mpz_add(temp1,temp1,temp2);
			mpz_sub(Uc[i],temp1,temp3);
			mpz_mul(temp1,wa,Sbc[i]); mpz_mul(temp2,wc,Sab[i]);
			mpz_mul(temp3,wb,Sac[i]); mpz_add(temp1,temp1,temp2);
			mpz_sub(Ud[i],temp1,temp3);
		}
		mpz_neg(val1,Det1); mpz_neg(val2,Det2);
		mpz_neg(val3,Det3); mpz_neg(val4,Det4);

/* First check face abc */

		if(!testu[1])
		{
//			printf("Test face abc \n");
			triangle_attach_gmp(Sd,Det1,Det2,Det3,Det4,Dabc,&test_a);
			if(!testa[1]) testa[1]=test_a;

			if(!trig_stat[0])
			{

				triangle_radius_gmp(Sd,Td,Ud,Dabc,&test_r,alp);
				testr[1]=test_r;

				edge_attach_gmp(a_mp, b_mp,Dab,Sab,Sd,Td,&test_a);
				if(!testa[5]) testa[5]=test_a;
				edge_attach_gmp(a_mp, c_mp, Dac,Sac,Sdm1,Tdm1,&test_a);
				if(!testa[6]) testa[6]=test_a;
				edge_attach_gmp(b_mp, c_mp, Dbc,Sbc,Sd,Td,&test_a);
				if(!testa[8]) testa[8]=test_a;

				vertex_attach_gmp(Dab,ra_mp,rb_mp,&test_a,&test_a2);
				if(!testa[11]) testa[11]=test_a;
				if(!testa[12]) testa[12]=test_a2;
				vertex_attach_gmp(Dac,ra_mp,rc_mp,&test_a,&test_a2);
				if(!testa[11]) testa[11]=test_a;
				if(!testa[13]) testa[13]=test_a2;
				vertex_attach_gmp(Dbc,rb_mp,rc_mp,&test_a,&test_a2);
				if(!testa[12]) testa[12]=test_a;
				if(!testa[13]) testa[13]=test_a2;

			}
			else
			{
				if((!testa[1]) && testr[1])
				{
					testu[5]=1;
					testu[6]=1;
					testu[8]=1;
					testr[5]=1;
					testr[6]=1;
					testr[8]=1;
					testu[11]=1;
					testu[12]=1;
					testu[13]=1;
				}
			}
		}

/* Now check face abd */
 
		if(!testu[2])
		{
//			printf("Test face abd \n");
			triangle_attach_gmp(Sc,val1,val2,val3,val4,Dabd,&test_a);
			if(!testa[2]) testa[2]=test_a;

			if(!trig_stat[1])
			{

				triangle_radius_gmp(Sc,Tc,Uc,Dabd,&test_r,alp);
				testr[2]=test_r;

				edge_attach_gmp(a_mp, b_mp, Dab,Sab,Sc,Tc,&test_a);
				if(!testa[5]) testa[5]=test_a;
				edge_attach_gmp(a_mp, d_mp,Dad,Sad,Scm1,Tcm1,&test_a);
				if(!testa[7]) testa[7]=test_a;
				edge_attach_gmp(b_mp, d_mp, Dbd,Sbd,Sc,Tc,&test_a);
				if(!testa[9]) testa[9]=test_a;

				vertex_attach_gmp(Dad,ra_mp,rd_mp,&test_a,&test_a2);
				if(!testa[11]) testa[11]=test_a;
				if(!testa[14]) testa[14]=test_a2;
				vertex_attach_gmp(Dbd,rb_mp,rd_mp,&test_a,&test_a2);
				if(!testa[12]) testa[12]=test_a;
				if(!testa[14]) testa[14]=test_a2;

			}
			else
			{
				if((!testa[2]) && testr[2])
				{
					testu[5]=1;
					testu[7]=1;
					testu[9]=1;
					testr[5]=1;
					testr[7]=1;
					testr[9]=1;
					testu[11]=1;
					testu[12]=1;
					testu[14]=1;
				}
			}
		}

/* Now check face acd */
 
		if(!testu[3])
		{
//			printf("Test face acd \n");
			triangle_attach_gmp(Sb,Det1,Det2,Det3,Det4,Dacd,&test_a);
			if(!testa[3]) testa[3]=test_a;

			if(!trig_stat[2])
			{

				triangle_radius_gmp(Sb,Tb,Ub,Dacd,&test_r,alp);
				testr[3]=test_r;

				edge_attach_gmp(a_mp, c_mp, Dac,Sac,Sb,Tb,&test_a);
				if(!testa[6]) testa[6]=test_a;
				edge_attach_gmp(a_mp, d_mp, Dad,Sad,Sbm1,Tbm1,&test_a);
				if(!testa[7]) testa[7]=test_a;
				edge_attach_gmp(c_mp, d_mp, Dcd,Scd,Sb,Tb,&test_a);
				if(!testa[10]) testa[10]=test_a;

				vertex_attach_gmp(Dcd,rc_mp,rd_mp,&test_a,&test_a2);
				if(!testa[13]) testa[13]=test_a;
				if(!testa[14]) testa[14]=test_a2;
			}
			else
			{
				if((!testa[3]) && testr[3])
				{
					testu[6]=1;
					testu[7]=1;
					testu[10]=1;
					testr[6]=1;
					testr[7]=1;
					testr[10]=1;
					testu[11]=1;
					testu[13]=1;
					testu[14]=1;
				}
			}
		}
 
/* Now check face bcd */
 
		if(!testu[4])
		{
//			printf("Test face bcd \n");
			triangle_attach_gmp(Sa,val1,val2,val3,val4,Dbcd,&test_a);
			if(!testa[4]) testa[4]=test_a;

			if(!trig_stat[3])
			{

				triangle_radius_gmp(Sa,Ta,Ua,Dbcd,&test_r,alp);
				testr[4]=test_r;

				edge_attach_gmp(b_mp, c_mp, Dbc,Sbc,Sa,Ta,&test_a);
				if(!testa[8]) testa[8]=test_a;
				edge_attach_gmp(b_mp, d_mp, Dbd,Sbd,Sam1,Tam1,&test_a);
				if(!testa[9]) testa[9]=test_a;
				edge_attach_gmp(c_mp, d_mp, Dcd,Scd,Sa,Ta,&test_a);
				if(!testa[10]) testa[10]=test_a;
			}
			else
			{
				if((!testa[4]) && testr[4])
				{
					testu[8]=1;
					testu[9]=1;
					testu[10]=1;
					testr[8]=1;
					testr[9]=1;
					testr[10]=1;
					testu[12]=1;
					testu[13]=1;
					testu[14]=1;
				}
			}
		}
/* 
 	Now check radii for each edge
									*/

		if((testu[5]==0) && (edge_stat[0]==0) ) {
			edge_radius_gmp(a_mp,b_mp, wa, wb, Dab, Sab,&test_r,alp);
			if(!testr[5]) testr[5]=test_r;
		}
		if((testu[6]==0) && (edge_stat[1]==0) ) {
			edge_radius_gmp(a_mp,c_mp, wa, wc, Dac, Sac,&test_r,alp);
			if(!testr[6]) testr[6]=test_r;
		}
		if((testu[7]==0) && (edge_stat[2]==0) ) {
			edge_radius_gmp(a_mp,d_mp, wa, wd, Dad, Sad,&test_r,alp);
			if(!testr[7]) testr[7]=test_r;
		}
		if((testu[8]==0) && (edge_stat[3]==0) ) {
			edge_radius_gmp(b_mp,c_mp, wb, wc, Dbc, Sbc,&test_r,alp);
			if(!testr[8]) testr[8]=test_r;
		}
		if((testu[9]==0) && (edge_stat[4]==0) ) {
			edge_radius_gmp(b_mp,d_mp, wb, wd, Dbd, Sbd,&test_r,alp);
			if(!testr[9]) testr[9]=test_r;
		}
		if((testu[10]==0) && (edge_stat[5]==0) ) {
			edge_radius_gmp(c_mp,d_mp, wc, wd, Dcd, Scd,&test_r,alp);
			if(!testr[10]) testr[10]=test_r;
		}
 

/*	There is no need to check radii of vertices, as they are equal to -radius**2, i.e. always
        negative                    */

/* End if */

	}

}
#ifdef __cplusplus
}
#endif

/*------------------------------------------------------------------------*/
/* set_alf_gmp	Version 1 11/24/2000	Patrice Koehl
 
	This procedure initialises all gmp variables that can be used for 
	computing the dual complex
*/

#ifdef __cplusplus
extern "C" {
#endif
void F77Name2(set_alf_gmp)()
{
/*	Initialise local GMP variables */

	int i,j;

	mpz_init(ra2); mpz_init(rb2); mpz_init(dist2);
	mpz_init(dtest);
	mpz_init (r_11); mpz_init (r_22); mpz_init (r_33);
	mpz_init (diff);
	mpz_init (det0);

	for (i= 0; i < 4; i++)
	{
		for (j=0; j < 5; j++)
		{
			mpz_init(res[i][j]);
			mpz_init(res2_c[i][j]);
		}
	}
	mpz_init (r_14); mpz_init (r_313); mpz_init (r_212);
	mpz_init (det1); mpz_init (det2); mpz_init (det3); mpz_init (det4); 

	for (i = 0; i < 4; i++) 
	{
		mpz_init(a_mp[i]);
		mpz_init(b_mp[i]);
		mpz_init(c_mp[i]);
		mpz_init(d_mp[i]);
	}
	mpz_init(wa);mpz_init(wb);mpz_init(wc);mpz_init(wd);
	mpz_init(ra_mp);mpz_init(rb_mp);mpz_init(rc_mp);mpz_init(rd_mp);

	mpz_init (num); mpz_init (den);

	mpz_init(val4);

	for (i=0; i < 4; i++)
	{
		mpz_init(Sab[i]); mpz_init(Sac[i]); mpz_init(Sad[i]);
		mpz_init(Sbc[i]); mpz_init(Sbd[i]); mpz_init(Scd[i]);
		mpz_init(Dab[i]); mpz_init(Dac[i]); mpz_init(Dad[i]);
		mpz_init(Dbc[i]); mpz_init(Dbd[i]); mpz_init(Dcd[i]);
		mpz_init(Sa[i]); mpz_init(Sb[i]); mpz_init(Sc[i]); 
		mpz_init(Sd[i]);
		mpz_init(Sam1[i]); mpz_init(Sbm1[i]); mpz_init(Scm1[i]); 
		mpz_init(Sdm1[i]);
		mpz_init(Ta[i]); mpz_init(Tb[i]); mpz_init(Tc[i]); 
		mpz_init(Td[i]);
		mpz_init(Tam1[i]); mpz_init(Tbm1[i]); mpz_init(Tcm1[i]); 
		mpz_init(Tdm1[i]);
		mpz_init(Ua[i]); mpz_init(Ub[i]); mpz_init(Uc[i]); 
		mpz_init(Ud[i]);
		mpz_init(Deter[i]);
	}

	mpz_init(Dabc); mpz_init(Dabd); mpz_init(Dacd); mpz_init(Dbcd);
	mpz_init(Det1); mpz_init(Det2); mpz_init(Det3); mpz_init(Dabcd);
	mpz_init(Det4);

	mpz_init(alp);

}
#ifdef __cplusplus
}
#endif

/*------------------------------------------------------------------------*/
/* clear_alf_gmp	Version 1 11/24/2000	Patrice Koehl
 
	This procedure clears all gmp variables that were used for 
	computing the dual complex
*/

#ifdef __cplusplus
extern "C" {
#endif
void clear_alf_gmp()
{
/*	Initialise local GMP variables */

	int i,j;

	mpz_clear(ra2); mpz_clear(rb2); mpz_clear(dist2);
	mpz_clear(dtest);
	mpz_clear (r_11); mpz_clear (r_22); mpz_clear (r_33);
	mpz_clear (diff);
	mpz_clear (det0);

	for (i= 0; i < 4; i++)
	{
		for (j=0; j < 5; j++)
		{
			mpz_clear(res[i][j]);
			mpz_clear(res2_c[i][j]);
		}
	}
	mpz_clear (r_14); mpz_clear (r_313); mpz_clear (r_212);
	mpz_clear (det1); mpz_clear (det2); mpz_clear (det3); mpz_clear (det4); 

	for (i = 0; i < 4; i++) 
	{
		mpz_clear(a_mp[i]);
		mpz_clear(b_mp[i]);
		mpz_clear(c_mp[i]);
		mpz_clear(d_mp[i]);
	}
	mpz_clear(wa);mpz_clear(wb);mpz_clear(wc);mpz_clear(wd);
	mpz_clear(ra_mp);mpz_clear(rb_mp);mpz_clear(rc_mp);mpz_clear(rd_mp);

	mpz_clear (num); mpz_clear (den);

	mpz_clear(val4);

	for (i=0; i < 4; i++)
	{
		mpz_clear(Sab[i]); mpz_clear(Sac[i]); mpz_clear(Sad[i]);
		mpz_clear(Sbc[i]); mpz_clear(Sbd[i]); mpz_clear(Scd[i]);
		mpz_clear(Dab[i]); mpz_clear(Dac[i]); mpz_clear(Dad[i]);
		mpz_clear(Dbc[i]); mpz_clear(Dbd[i]); mpz_clear(Dcd[i]);
		mpz_clear(Sa[i]); mpz_clear(Sb[i]); mpz_clear(Sc[i]); 
		mpz_clear(Sd[i]);
		mpz_clear(Sam1[i]); mpz_clear(Sbm1[i]); mpz_clear(Scm1[i]); 
		mpz_clear(Sdm1[i]);
		mpz_clear(Ta[i]); mpz_clear(Tb[i]); mpz_clear(Tc[i]); 
		mpz_clear(Td[i]);
		mpz_clear(Tam1[i]); mpz_clear(Tbm1[i]); mpz_clear(Tcm1[i]); 
		mpz_clear(Tdm1[i]);
		mpz_clear(Ua[i]); mpz_clear(Ub[i]); mpz_clear(Uc[i]); 
		mpz_clear(Ud[i]);
		mpz_clear(Deter[i]);
	}

	mpz_clear(Dabc); mpz_clear(Dabd); mpz_clear(Dacd); mpz_clear(Dbcd);
	mpz_clear(Det1); mpz_clear(Det2); mpz_clear(Det3); mpz_clear(Dabcd);
	mpz_clear(Det4);

	mpz_clear(alp);

}
#ifdef __cplusplus
}
#endif

