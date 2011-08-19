/*	alf_tools_gmp.c	Version 2 3/30/2007	Patrice Koehl             */
/*									  */
/*  This is the C version of alfcx_tools.f, which performs all operations */
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

void set_edge(int a, int b);

void set_triangle(int a, int b, int c);

#ifdef __cplusplus
extern "C"
#endif
void F77Name2(tetra_radius_gmp)(int *ia, int *ib, int *ic, int *id, int *test, double *scale, double *alpha);

#ifdef __cplusplus
extern "C"
#endif
void F77Name2(vertex_attach_gmp)(int *ia, int *ib, int *testa, int *testb);

#ifdef __cplusplus
extern "C"
#endif
void F77Name2(edge_attach_gmp)(int *ia, int *ib, int *ic, int *test, int *memory);

#ifdef __cplusplus
extern "C"
#endif
void F77Name2(edge_radius_gmp)(int *ia, int *ib, int *test, double *scale, double *alpha, int *memory);

#ifdef __cplusplus
extern "C"
#endif
void F77Name2(triangle_attach_gmp)(int *ia, int *ib, int *ic, int *id, int *test, int *memory);

#ifdef __cplusplus
extern "C"
#endif
void F77Name2(triangle_radius_gmp)(int *ia, int *ib, int *ic, int *test, double *scale, double *alpha, int *memory);

void F77Name2(set_alf_gmp)();
void clear_alf_gmp();

/*------------------------------------------------------------------------*/
/* set_triangle:
	This subroutine sets all common arrays for gmp calculation over edges
	Input:
		ia,ib: indices of the two points considered
*/
void set_triangle(int a, int b, int c)
{
	int idx_a1,idx_b1, idx_c1;
	int i,j;

	idx_a1 = (a)*3 -3; 
	idx_b1 = (b)*3 -3;
	idx_c1 = (c)*3 -3;


//       0. define coordinates


	for (i=0; i<3; i++)
	{
		mpz_set(a_mp[i+1],coord_gmp[idx_a1+i]);
		mpz_set(b_mp[i+1],coord_gmp[idx_b1+i]);
		mpz_set(c_mp[i+1],coord_gmp[idx_c1+i]);
	}
	
	idx_a1 = (a)-1; idx_b1 = (b)-1; idx_c1 = (c)-1;
	mpz_set(a_mp[4],weight_gmp[idx_a1]);
	mpz_set(b_mp[4],weight_gmp[idx_b1]);
	mpz_set(c_mp[4],weight_gmp[idx_c1]);

/*
       1. Computes all Minors Mab(i,j)= M(a,b,i,j)   = Det | a(i)  a(j) |
                                                           | b(i)  b(j) |
*/
	for (i=1;  i<4; i++)
	{
		for (j=i+1; j<5 ; j++)
		{
			mpz_mul(temp1,a_mp[j],b_mp[i]); 
			mpz_mul(temp2,a_mp[i],b_mp[j]);
			mpz_sub(Mab[i][j],temp2,temp1);
			mpz_mul(temp1,a_mp[j],c_mp[i]); 
			mpz_mul(temp2,a_mp[i],c_mp[j]);
			mpz_sub(Mac[i][j],temp2,temp1);
			mpz_mul(temp1,b_mp[j],c_mp[i]); 
			mpz_mul(temp2,b_mp[i],c_mp[j]);
			mpz_sub(Mbc[i][j],temp2,temp1);
		}
	}

/*
       Now compute all Minors
               S(i,j) = M(a,b,c,i,j,0)    = Det | a(i) a(j) 1 |
                                                | b(i) b(j) 1 |
                                                | c(i) c(j) 1 |

       a,b,c are the 3 vertices of the triangle, i and j correspond
       to two of the coordinates of the vertices

       for all i in [1,3] and all j in [i+1,4]
*/

	for (i=1;  i<4; i++)
	{
		for (j=i+1; j<5 ; j++)
		{
			mpz_sub(temp1,Mbc[i][j],Mac[i][j]);
			mpz_add(S[i][j],temp1,Mab[i][j]);
		}
	}

/*
       Now compute all Minors
               T(i,j) = M(a,b,c,i,j,4)    = Det | a(i) a(j) a(4) |
                                                | b(i) b(j) b(4) |
                                                | c(i) c(j) c(4) |

       for all i in [1,2] and all j in [i+1,3]
*/

	for (i=1;  i<3; i++)
	{
		for (j=i+1; j<4 ; j++)
		{
			mpz_mul(temp1,a_mp[4],Mbc[i][j]);
			mpz_mul(temp2,b_mp[4],Mac[i][j]);
			mpz_sub(temp1,temp1,temp2);
			mpz_mul(temp2,c_mp[4],Mab[i][j]);
			mpz_add(T[i][j],temp1,temp2);
		}
	}
/*
       Finally,  need Dabc = M(a,b,c,1,2,3) Det | a(1) a(2) a(3) |
                                                | b(2) b(2) b(3) |
                                                | c(3) c(2) c(3) |
*/
	mpz_mul(temp1,a_mp[1],Mbc[2][3]); mpz_mul(temp2,b_mp[1],Mac[2][3]); mpz_sub(temp1,temp1,temp2);
	mpz_mul(temp2,c_mp[1],Mab[2][3]); mpz_add(Dabc,temp1,temp2);

}
/*------------------------------------------------------------------------*/
/* set_edge:
	This subroutine sets all common arrays for gmp calculation over edges
	Input:
		ia,ib: indices of the two points considered
*/
void set_edge(int a, int b)
{
	int idx_a1,idx_b1;
	int i,j,k;

	idx_a1 = (a)*3 -3; 
	idx_b1 = (b)*3 -3;


//       0. define coordinates


	for (i=0; i<3; i++)
	{
		mpz_set(a_mp[i+1],coord_gmp[idx_a1+i]);
		mpz_set(b_mp[i+1],coord_gmp[idx_b1+i]);
	}
	
	idx_a1 = (a)-1; idx_b1 = (b)-1; 

	mpz_set(a_mp[4],weight_gmp[idx_a1]);
	mpz_set(b_mp[4],weight_gmp[idx_b1]);

/*
       1. Compute all Minors Dab(i) = M(a,b,i,0) = Det | a(i) 1 |
                                                       | b(i) 1 |
*/
	for (i=1; i<5; i++)
	{
		mpz_sub(Dab[i],a_mp[i],b_mp[i]);
	}

/*
       2. Computes all Minors Sab(i,j)= M(a,b,i,j)   = Det | a(i)  a(j) |
                                                           | b(i)  b(j) |
*/
	for (i=1;  i<3; i++)
	{
		for (j=i+1; j<4 ; j++)
		{
			k=i+j-2;
			mpz_mul(temp1,a_mp[j],b_mp[i]); 
			mpz_mul(temp2,a_mp[i],b_mp[j]);
			mpz_sub(Sab[k],temp2,temp1);
		}
	}

/*
      3. Computes all Minors Tab(i)= M(a,b,i,4)   = Det | a(i)  a(4) |
                                                         | b(i)  b(4) |
*/
	for (i=1; i<4; i++)
	{
		mpz_mul(temp1,a_mp[i],b_mp[4]); 
		mpz_mul(temp2,a_mp[4],b_mp[i]);
		mpz_sub(Tab[i],temp1,temp2);
	}
}

/*------------------------------------------------------------------------*/
/* tetra_radius_gmp	Version 1 11/24/2000	Patrice Koehl
 
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
void F77Name2(tetra_radius_gmp)( int *a, int *b, int *c, int *d, int *testr,
		double *scale, double *alpha)

{
	int i,j,k,coef;
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

	mpz_set(det1,Deter[3]);
	mpz_set(det2,Deter[2]);
	mpz_set(det3,Deter[1]);

	mpz_mul(temp1,a_mp[1],Sa[3]);mpz_mul(temp2,b_mp[1],Sb[3]);
	mpz_sub(temp3,temp1,temp2);
	mpz_mul(temp1,c_mp[1],Sc[3]);mpz_mul(temp2,d_mp[1],Sd[3]);
	mpz_sub(temp1,temp1,temp2);
	mpz_add(det4,temp1,temp3);

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

	mpz_mul(temp1,det4,det4); coef=4; mpz_mul_si(den,temp1,coef);

	mpz_mul(temp1,det1,det1); mpz_mul(temp2,det2,det2);
	mpz_add(temp1,temp1,temp2); mpz_mul(temp2,det3,det3);
	mpz_add(temp1,temp1,temp2); mpz_mul(temp2,det4,Dabcd);
	mpz_mul_si(temp2,temp2,coef); mpz_add(num,temp1,temp2);
	
	mpz_mul(temp1,den,alp); mpz_sub(temp2,num,temp1);
//	printf("Dtest = %s\n",mpz_get_str(NULL,10,temp2));
	
/* 
 	If tetrahedron is part of the alpha shape, then the 4 triangles,
 	the 6 edges and the four vertices are also part of the alpha
 	complex
 								*/
	if(!(mpz_sgn(temp2) > 0)) 
	{
		(*testr)=1;
	}
	else
	{
		(*testr)=0;
	}
}
#ifdef __cplusplus
}
#endif
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
#ifdef __cplusplus
extern "C" {
#endif
void	F77Name2(vertex_attach_gmp)(int *ia, int *ib, int *testa, int *testb)

{
	int idx_a1,idx_b1;
	int i;

	(*testa = 0);
	(*testb = 0);

	idx_a1 = (*ia)*3 -3; 
	idx_b1 = (*ib)*3 -3;

	for (i=0; i<4; i++)
	{
		mpz_set(temp1,coord_gmp[idx_a1+i]);
		mpz_set(temp2,coord_gmp[idx_b1+i]);
		mpz_sub(Dab[i],temp1,temp2);
	}

	mpz_set(ra_mp,radius_gmp[idx_a1]);
	mpz_set(rb_mp,radius_gmp[idx_b1]);

	mpz_mul(temp1, Dab[0],Dab[0]);
	mpz_mul(temp2, Dab[1],Dab[1]);
	mpz_mul(temp3, Dab[2],Dab[2]);
	mpz_add(temp1, temp1,temp2); mpz_add(dist2,temp1,temp3);

	mpz_mul(ra2, ra_mp, ra_mp); mpz_mul(rb2, rb_mp, rb_mp);

	mpz_add(dtest,dist2,ra2);
	mpz_sub(dtest,dtest,rb2);
	if(mpz_sgn(dtest) < 0) (*testa = 1);

	mpz_sub(dtest,dist2,ra2);
	mpz_add(dtest,dtest,rb2);
	if(mpz_sgn(dtest) < 0) (*testb = 1);

}
#ifdef __cplusplus
}
#endif

/*------------------------------------------------------------------------*/
/* edge_attach_gmp:
	This subroutine checks if an edge of the regular triangulation is "attached"
	to another vertex (i.e.  if the vertex belongs to the smallest circumsphere
	of the edge). 

									*/

#ifdef __cplusplus
extern "C" {
#endif
void F77Name2(edge_attach_gmp)(int *ia, int *ib, int *ic, int *testa, int *memory)

{
	int i,j,k,coef;
	int idx_c;

/* Set up calculation, if not already done) */

	if((*memory) != 1) set_edge(*ia,*ib);

	idx_c = (*ic)*3 -3; 

	for (i=0; i<3; i++)
	{
		mpz_set(c_mp[i+1],coord_gmp[idx_c+i]);
	}
	idx_c = (*ic) -1; 
	mpz_set(c_mp[4],weight_gmp[idx_c]);

/* 
       Need to compute:
       Sc      : minor(a,b,c,i,j,0) for i=1,2 and j = i+1,3
       Tc      : minor(a,b,c,i,4,0) for i = 1,2,3
*/

	for (i=1; i<3 ; i++)
	{
		for (j=i+1; j<4; j++)
		{
			k=i+j-2;
			mpz_mul(temp1,c_mp[i],Dab[j]);
			mpz_mul(temp2,c_mp[j],Dab[i]);
			mpz_sub(temp1,temp1,temp2);
			mpz_add(Sc[k],temp1,Sab[k]);
		}
	}

	for (i=1; i<4; i++)
	{
		mpz_mul(temp1,c_mp[i],Dab[4]);
		mpz_mul(temp2,c_mp[4],Dab[i]);
		mpz_sub(temp1,temp1,temp2);
		mpz_add(Tc[i],temp1,Tab[i]);
	}

/*	This is the "hidden1" part */

	(*testa) = 0;

	if( mpz_cmp(a_mp[1],b_mp[1]) != 0) 
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
	else if ( mpz_cmp(a_mp[2],b_mp[2]) != 0)
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
	else if (  mpz_cmp(a_mp[3],b_mp[3]) != 0)
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
//	printf("d0 (edge attach) = %s\n",mpz_get_str(NULL,10,det0));
//	printf("d5 (edge attach) = %s\n",mpz_get_str(NULL,10,temp3));

	if(mpz_sgn(dtest) < 0) (*testa = 1);

}
#ifdef __cplusplus
}
#endif

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

#ifdef __cplusplus
extern "C" {
#endif
void F77Name2(edge_radius_gmp)(int *ia, int *ib, int *testr, double *scale, double *alpha, int *memory)

{
	int i,j,coef, ivalue;
	double value;

	value = (*alpha)*(*scale); ivalue = (int) floor(value); 
	mpz_set_si(alp,ivalue);

	if((*memory) != 1) set_edge(*ia,*ib);

/*	This is the "hidden1" part */

	(*testr) = 0;
	mpz_set(res[0][4],Dab[4]);


	if( mpz_cmp(a_mp[1],b_mp[1]) != 0) 
	{
		for (i = 1; i < 4 ; i++)
		{
			mpz_set(res[0][i],Dab[i]);
			mpz_mul(temp1,b_mp[i],a_mp[4]);
			mpz_mul(temp2,a_mp[i],b_mp[4]);
			mpz_sub(res[i][4],temp2,temp1);
		}
		mpz_set(res[1][2],Sab[1]); mpz_set(res[1][3],Sab[2]);
		mpz_set(res[2][3],Sab[3]);
	}
	else if ( mpz_cmp(a_mp[2],b_mp[2]) != 0)
	{
		mpz_set(res[0][1],Dab[2]); mpz_set(res[0][2],Dab[3]);
		mpz_set(res[0][3],Dab[1]);
		mpz_set(res[1][2],Sab[3]);
		mpz_neg(res[1][3],Sab[1]); mpz_neg(res[2][3],Sab[2]);
		mpz_mul(temp1,a_mp[2],b_mp[4]); mpz_mul(temp2,b_mp[2],a_mp[4]);
		mpz_sub(res[1][4],temp1,temp2);
		mpz_mul(temp1,a_mp[3],b_mp[4]); mpz_mul(temp2,b_mp[3],a_mp[4]);
		mpz_sub(res[2][4],temp1,temp2);
		mpz_mul(temp1,a_mp[1],b_mp[4]); mpz_mul(temp1,b_mp[1],a_mp[4]);
		mpz_sub(res[3][4],temp1,temp2);
	}
	else if (  mpz_cmp(a_mp[3],b_mp[3]) != 0)
	{
		mpz_set(res[0][1],Dab[3]); mpz_set(res[0][2],Dab[1]);
		mpz_set(res[0][3],Dab[2]);
		mpz_neg(res[1][2],Sab[2]);
		mpz_neg(res[1][3],Sab[3]); mpz_set(res[2][3],Sab[1]);
		mpz_mul(temp1,a_mp[3],b_mp[4]); mpz_mul(temp2,b_mp[3],a_mp[4]);
		mpz_sub(res[1][4],temp1,temp2);
		mpz_mul(temp1,a_mp[1],b_mp[4]); mpz_mul(temp2,b_mp[1],a_mp[4]);
		mpz_sub(res[2][4],temp1,temp2);
		mpz_mul(temp1,a_mp[2],b_mp[4]); mpz_mul(temp1,b_mp[2],a_mp[4]);
		mpz_sub(res[3][4],temp1,temp2);
	}
	else
	{
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
#ifdef __cplusplus
}
#endif

/*------------------------------------------------------------------------*/
/* triangle_radius_gmp:
   This program checks if the radius of the circumsphere of a facet of the
   regular triangulation is smaller than alpha


       Input:

       For the three points a,b,c that form the triangles, the program
       needs as input the following determinants:

       S(i,j)   = Minor(a,b,c,i,j,0)= det | a(i)  a(j)  1 |
                                          | b(i)  b(j)  1 |
                                          | c(i)  c(j)  1 |

       for i in [1,3] and j in [i+1,4]

       and:

       T(i,j) = Minor(a,b,c,i,j,4)=det | a(i) a(j) a(4) |
                                       | b(i) b(j) b(4) |
                                       | c(i) c(j) c(4) |

       and

       Dabc  = Minor(a,b,c,1,2,3)

	Output:

	testr	: flag set to 1 if ALPHA is larger than rho, the radius
		  of the circumsphere of the triangle

*/

void F77Name2(triangle_radius_gmp)(int *ia, int *ib, int *ic, int *testr, double *scale, double *alpha, int *memory)

{
	int i,j,coef, ivalue;
	double value;

	value = (*alpha)*(*scale); ivalue = (int) floor(value); 
	mpz_set_si(alp,ivalue);

	if(*memory !=1) set_triangle(*ia, *ib, *ic);

	(*testr) = 0;

	mpz_set_si(temp1,0);
	for (i=1; i<3; i++)
	{
		for (j=i+1; j<4; j++)
		{
			mpz_mul(temp2,S[i][j],S[i][j]);
			mpz_add(temp1,temp1,temp2);
		}
	}

/* Compute det0 */

	coef = 4;
	mpz_mul_si(det0,temp1,coef);

/* Compute det1 */

	mpz_mul(temp1,Dabc,S[2][3]);
	coef = -2;
	mpz_mul_si(temp1,temp1,coef);
	mpz_mul(temp2,S[1][2],S[2][4]);
	mpz_add(temp1,temp2,temp1);
	mpz_mul(temp2,S[1][3],S[3][4]);
	mpz_add(temp1,temp2,temp1);
	mpz_mul_si(det1,temp1,coef);

/* Compute det2 */

	coef = 2;
	mpz_mul(temp1,Dabc,S[1][3]);
	mpz_mul_si(temp1,temp1,coef);
	mpz_mul(temp2,S[2][3],S[3][4]);
	mpz_add(temp1,temp2,temp1);
	mpz_mul(temp2,S[1][2],S[1][4]);
	mpz_sub(temp1,temp2,temp1);
	mpz_mul_si(det2,temp1,coef);

/* Compute det3 */

	mpz_mul(temp1,Dabc,S[1][2]);
	mpz_mul_si(temp1,temp1,coef);
	mpz_mul(temp2,S[1][3],S[1][4]);
	mpz_add(temp1,temp2,temp1);
	mpz_mul(temp2,S[2][3],S[2][4]);
	mpz_add(temp1,temp2,temp1);
	mpz_mul_si(det3,temp1,coef);

/* Compute det4 */

	mpz_mul(temp1,Dabc,Dabc);
	coef = -2;
	mpz_mul_si(temp1,temp1,coef);

	for (i=1; i<3; i++)
	{
		for (j=i+1; j<4; j++)
		{
			mpz_mul(temp2,S[i][j],T[i][j]);
			mpz_add(temp1,temp1,temp2);
		}
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
#ifdef __cplusplus
}
#endif
/*------------------------------------------------------------------------*/
/* triangle_attach_gmp:
   This program checks if a facet is attached to a vertex

	Input:

	For the three points a,b,c that form the triangles, the program
	needs as input the following determinants:

         S(i,j) = Minor(a,b,c,i,j,0)= det | a(i)  a(j)  1 |
                                          | b(i)  b(j)  1 |
                                          | c(i)  c(j)  1 |
       for all i in [1,3], j in [i+1,4]

       T(i,j) = M(a,b,c,i,j,4)    = Det | a(i) a(j) a(4) |
                                        | b(i) b(j) b(4) |
                                        | c(i) c(j) c(4) |

       for all i in [1,2] and all j in [i+1,3]

       Dabc = Det | a(1) a(2) a(3) |
                  | b(1) b(2) b(3) |
                  | c(1) c(2) c(3) |

       and the coordinates of the fourth vertex d


	testa	: flag set to 1 if the fourth point d is inside the
		  circumsphere of (a,b,c)

*/

#ifdef __cplusplus
extern "C" {
#endif
void F77Name2(triangle_attach_gmp)(int *ia, int *ib, int *ic, int *id, int *testa, int *memory)

{
	int i,coef, idx_d;

	if(*memory !=1) set_triangle(*ia, *ib, *ic);

	idx_d = (*id)*3 -3; 

	for (i=0; i<3; i++)
	{
		mpz_set(d_mp[i+1],coord_gmp[idx_d+i]);
	}
	idx_d = (*id) -1; 
	mpz_set(d_mp[4],weight_gmp[idx_d]);

/*
       We need to compute:

       det1 = Minor(a,b,c,d,2,3,4,0)
       det2 = Minor(a,b,c,d,1,3,4,0)
       det3 = Minor(a,b,c,d,1,2,4,0)
       det4 = Minor(a,b,c,d,1,2,3,0)
*/
	mpz_mul(temp1,d_mp[2],S[3][4]); mpz_mul(temp2,d_mp[3],S[2][4]); mpz_sub(temp1,temp2,temp1);
	mpz_mul(temp2,d_mp[4],S[2][3]); mpz_sub(temp2,T[2][3],temp2); mpz_add(det1,temp2,temp1);
	mpz_mul(temp1,d_mp[1],S[3][4]); mpz_mul(temp2,d_mp[3],S[1][4]); mpz_sub(temp1,temp2,temp1);
	mpz_mul(temp2,d_mp[4],S[1][3]); mpz_sub(temp2,T[1][3],temp2); mpz_add(det2,temp2,temp1);
	mpz_mul(temp1,d_mp[1],S[2][4]); mpz_mul(temp2,d_mp[2],S[1][4]); mpz_sub(temp1,temp2,temp1);
	mpz_mul(temp2,d_mp[4],S[1][2]); mpz_sub(temp2,T[1][2],temp2); mpz_add(det3,temp2,temp1);
	mpz_mul(temp1,d_mp[1],S[2][3]); mpz_mul(temp2,d_mp[2],S[1][3]); mpz_sub(temp1,temp2,temp1);
	mpz_mul(temp2,d_mp[3],S[1][2]); mpz_sub(temp2,Dabc,temp2); mpz_add(det4,temp2,temp1);

/* check if triangle is attached                               */

	(*testa) = 0;

	mpz_set_si(temp1,1);
/*	for (i=1; i<4; i++)
	{
		mpz_mul(temp2,S[i],S[i]);
		mpz_add(temp1,temp1,temp2);
	}
*/
	mpz_mul(temp2,det4,Dabc);
	coef = -2;
	mpz_mul_si(temp2,temp2,coef);
	mpz_mul(temp3,det3,S[1][2]);
	mpz_add(temp2,temp3,temp2);
	mpz_mul(temp3,det2,S[1][3]);
	mpz_add(temp2,temp3,temp2);
	mpz_mul(temp3,det1,S[2][3]);
	mpz_add(temp2,temp3,temp2);
	mpz_mul(dtest,temp1,temp2);
//	printf("Dtest (triangle attach) = %s\n",mpz_get_str(NULL,10,dtest));

	if(mpz_sgn(dtest) > 0) (*testa)=1;

/* 	Clear local GMP variables */

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
			mpz_init(Mab[i][j]);
			mpz_init(Mac[i][j]);
			mpz_init(Mbc[i][j]);
			mpz_init(S[i][j]);
			mpz_init(T[i][j]);
		}
	}
	mpz_init (r_14); mpz_init (r_313); mpz_init (r_212);
	mpz_init (det1); mpz_init (det2); mpz_init (det3); mpz_init (det4); 
	mpz_init (wa); mpz_init(wb); mpz_init(wc); mpz_init(wd);

	for (i = 0; i < 5; i++) 
	{
		mpz_init(a_mp[i]);
		mpz_init(b_mp[i]);
		mpz_init(c_mp[i]);
		mpz_init(d_mp[i]);
		mpz_init(Dab[i]);
	}
	mpz_init(ra_mp);mpz_init(rb_mp);
	mpz_init(rc_mp);mpz_init(rd_mp);

	mpz_init (num); mpz_init (den);

	for (i=0; i < 4; i++)
	{
		mpz_init(Sab[i]);
		mpz_init(Sac[i]);
		mpz_init(Sad[i]);
		mpz_init(Sbc[i]);
		mpz_init(Sbd[i]);
		mpz_init(Scd[i]);
		mpz_init(Tab[i]);
		mpz_init(Sa[i]);
		mpz_init(Sb[i]);
		mpz_init(Sc[i]);
		mpz_init(Sd[i]);
		mpz_init(Sam1[i]);
		mpz_init(Sbm1[i]);
		mpz_init(Scm1[i]);
		mpz_init(Sdm1[i]);
		mpz_init(Tc[i]);
		mpz_init(Deter[i]);
	}

	mpz_init(Dabc); mpz_init(Dabd); mpz_init(Dacd); mpz_init(Dbcd);
	mpz_init(alp);mpz_init(Dabcd);

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
			mpz_clear(Mab[i][j]);
			mpz_clear(Mac[i][j]);
			mpz_clear(Mbc[i][j]);
			mpz_clear(S[i][j]);
		}
	}
	for (i= 0; i < 3; i++)
	{
		for (j=0; j < 4; j++)
		{
			mpz_clear(T[i][j]);
		}
	}
	mpz_clear (r_14); mpz_clear (r_313); mpz_clear (r_212);
	mpz_clear (det1); mpz_clear (det2); mpz_clear (det3); mpz_clear (det4); 
	mpz_clear (wa); mpz_clear(wb); mpz_clear(wc); mpz_clear(wd);

	for (i = 0; i < 5; i++) 
	{
		mpz_clear(a_mp[i]);
		mpz_clear(b_mp[i]);
		mpz_clear(c_mp[i]);
		mpz_clear(d_mp[i]);
		mpz_clear(Dab[i]);
	}
	mpz_clear(ra_mp);mpz_clear(rb_mp);
	mpz_clear(rc_mp);mpz_clear(rd_mp);

	mpz_clear (num); mpz_clear (den);

	for (i=0; i < 4; i++)
	{
		mpz_clear(Sab[i]);
		mpz_clear(Sac[i]);
		mpz_clear(Sad[i]);
		mpz_clear(Sbc[i]);
		mpz_clear(Sbd[i]);
		mpz_clear(Scd[i]);
		mpz_clear(Tab[i]);
		mpz_clear(Sa[i]);
		mpz_clear(Sb[i]);
		mpz_clear(Sc[i]);
		mpz_clear(Sd[i]);
		mpz_clear(Sam1[i]);
		mpz_clear(Sbm1[i]);
		mpz_clear(Scm1[i]);
		mpz_clear(Sdm1[i]);
		mpz_clear(Tc[i]);
		mpz_clear(Deter[i]);
	}

	mpz_clear(Dabc);mpz_clear(Dabd); mpz_clear(Dacd); mpz_clear(Dbcd);
	mpz_clear(Dabcd);
	mpz_clear(alp);

}
#ifdef __cplusplus
}
#endif

