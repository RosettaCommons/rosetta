// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file external/calibur/rmsd.cc
/// @author YK Ng & SC Li (kalngyk@gmail.com)


#ifndef __WIN32__
#include <sys/resource.h>
#endif
#include <stdio.h>
#include <math.h>
#include <protocols/cluster/calibur/jacobi.hh>
#include <protocols/cluster/calibur/cubic.hh>
#include <protocols/cluster/calibur/rmsd.hh>

namespace protocols {
namespace cluster {
namespace calibur {

/******************************************************************************
* Brief description of the Kabsch algorithm implemented herein:
*
* (1) We start with the RMSD = min_T sqrt(sum_{i=1 to n} |p_i- T q_i|^2),
*  where T is a rigid transformation.
*
*  First of all, a basic result says that the centroids of P and Q must
*  coincide to minimize the RMSD. Hence the centroids are first aligned,
*  and the equation is changed to
*
*  RMSD = min_R sqrt(sum_{i=1 to n} |p_i- R q_i|^2), where R is a rotation.
*    = sqrt(sum_{i=1 to n} |p_i- Rmin q_i|^2), ie let Rmin be the min R
*    = sqrt(sum_{i=1 to n} p_i^2 + q_i^2 - 2p_i Rmin q_i)
*
*  Now (sum_{i=1 to n} p_i^2 + q_i^2) is invariant, so we only need to
*  find Rmin which maximizes (sum_{i=1 to n} pi Rmin q_i), or (Pt)RQ
*
*  However, (Pt)RQ is an nxn matrix. How can we compute it efficiently?
*
* (2) Consider the 3x3 matrix C=Q(Pt) instead
*
*  Algebraic manipulation shows that if we decompose C into its SVD,
*     C=UtSV
*  say, then the singular values, S, will give us the solution to the
*  problem (details omitted). More precisely,
*  (i)  [RMSD] Suppose the singular values of S are l1, l2, l3, then
*     RMSD = sum_{i=1 to n} p_i^2 + q_i^2 - 2(l1 + l2 + l3) / n
*  (ii) [Rotation] The optimal rotation is VUt.
*
*  But how do we compute the SVD of C efficiently?
*
* (3) Consider CCt instead (see "Parallel Algorithms for the Singular Value
*  Decomposition", Berry et al.)
*  Reasons:
*  (i)   CCt and C are related by CCt=(UtSV)(VtStU)=Ut S^2 U. Hence,
*     (a) The singular values of C are the square roots of the singular
*      values of CCt, and
*     (b) The right singular vectors of C and CCt are the same.
*  (ii)  Since CCt is real and symmetric,
*     (a) The jacobi method for eigenvalues & eigenvectors can be used
*     (b) Its eigenvalues are the same as its singular values
*     (c) Its eigenvectors give us U, the left singular vector of C
*     (d) The right singular vector V=(v1, v2, v3) of C can be computed
*      as v_i = Ct u_i/e_i, where (u1, u2, u3)=U, and e_i is the i-th
*      eigenvalue of C.
*
* (4) If the optimal rotation is not required, then the eigenvectors in
*  (3)(ii)(c-d) are not needed. In this case, instead of the jacobi method
*  in (3)(ii)(a), a more efficient method through solving cubic equations
*  can be used to obtain the eigenvalues of CtC (note that we used CCt
*  earlier -- both CCt and CtC will work with only slight differences).
*  This is what fast_rmsd() does.
*****************************************************************************/

#define cubic_roots(a,b,c,z) cubic_roots2(a,b,c,z)

double RMSD( std::vector<double> & coords1, std::vector<double> & coords2, int n, double R[3][3] )
{
	/*************************************************************************
	* Compute the centroids of coords1 and coords2 respectively    *
	*************************************************************************/

	double centroid1[3] = {0., 0., 0.};
	double centroid2[3] = {0., 0., 0.};
	for ( int j=0; j < n; j++ ) {
		int j3 = 3*j;
		centroid1[0] += coords1[j3  ];
		centroid2[0] += coords2[j3  ];
		centroid1[1] += coords1[j3+1];
		centroid2[1] += coords2[j3+1];
		centroid1[2] += coords1[j3+2];
		centroid2[2] += coords2[j3+2];
	}
	centroid1[0] /= n;
	centroid2[0] /= n;
	centroid1[1] /= n;
	centroid2[1] /= n;
	centroid1[2] /= n;
	centroid2[2] /= n;

	/*************************************************************************
	* Translate coords1 and coords2 such that their centroids are at the *
	* origin.                  *
	*************************************************************************/

	for ( int m=0; m < n; m++ ) {
		int m3 = 3*m;
		coords1[m3  ] -= centroid1[0];
		coords2[m3  ] -= centroid2[0];
		coords1[m3+1] -= centroid1[1];
		coords2[m3+1] -= centroid2[1];
		coords1[m3+2] -= centroid1[2];
		coords2[m3+2] -= centroid2[2];
	}

	/*************************************************************************
	* Compute C=Q(Pt) and the sum_{i=1 to n} p_i^2 + q_i^2 term (Eo)  *
	*************************************************************************/

	double C[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
	double Eo = 0.0;
	for ( int m=0; m < n; m++ ) {
		int m3 = 3*m;
		C[0][0] += coords2[m3  ] * coords1[m3  ];
		C[0][1] += coords2[m3  ] * coords1[m3+1];
		C[0][2] += coords2[m3  ] * coords1[m3+2];
		C[1][0] += coords2[m3+1] * coords1[m3  ];
		C[1][1] += coords2[m3+1] * coords1[m3+1];
		C[1][2] += coords2[m3+1] * coords1[m3+2];
		C[2][0] += coords2[m3+2] * coords1[m3  ];
		C[2][1] += coords2[m3+2] * coords1[m3+1];
		C[2][2] += coords2[m3+2] * coords1[m3+2];
		Eo += coords1[m3  ]*coords1[m3  ] + coords2[m3  ]*coords2[m3  ]
			+ coords1[m3+1]*coords1[m3+1] + coords2[m3+1]*coords2[m3+1]
			+ coords1[m3+2]*coords1[m3+2] + coords2[m3+2]*coords2[m3+2];
	}
	Eo *= 0.5;

#ifdef __DEBUG_RMSD__
	std::cout << "	[" <<C[0][0] << " " <<C[0][1] << " " <<C[0][2] << "]" << std::endl;
	std::cout << "C = [" <<C[1][0] << " " <<C[1][1] << " " <<C[1][2] << "]" << std::endl;
	std::cout << "	[" <<C[2][0] << " " <<C[2][1] << " " <<C[2][2] << "]" << std::endl;
	std::cout << std::endl;
#endif

	/*************************************************************************
	* Prepare CCt for eigenvector decomposition        *
	*************************************************************************/

	double Ct[3][3], CCt[3][3];
	for ( int i=0; i < 3; i++ ) {
		for ( int j=0; j < 3; j++ ) {
			Ct[i][j] = C[j][i];
		}
	}
	for ( int i=0; i<3; i++ ) {
		for ( int j=0; j < 3; j++ ) {
			CCt[i][j] = 0.;
			for ( int k = 0; k < 3; k++ ) {
				CCt[i][j] += C[i][k] * C[j][k]; // C[i][k]*Ct[k][j]
			}
		}
	}

	/*************************************************************************
	* Decompose CCt for its eigenvalues/eigenvectors, placed in    *
	* eigenval{0,1,2} and eigenvec{0,1,2} respectively.      *
	*************************************************************************/

	double eigenval[3], eigenvec[3][3];
	if ( !jacobi3(CCt, eigenval, eigenvec) ) { return -1.; }

	/* Sort the eigenvalues such that eigenval0 >= eigenval1 >= eigenval2 */
	double eigenval0, eigenval1, eigenval2;
	double eigenvec0[3], eigenvec1[3], eigenvec2[3];
	if ( eigenval[0] >= eigenval[1] ) {
		if ( eigenval[1] >= eigenval[2] ) {
			eigenval0 = eigenval[0];
			for ( int i=0; i < 3; i++ ) eigenvec0[i] = eigenvec[i][0];
			eigenval1 = eigenval[1];
			for ( int i=0; i < 3; i++ ) eigenvec1[i] = eigenvec[i][1];
			eigenval2 = eigenval[2];
			for ( int i=0; i < 3; i++ ) eigenvec2[i] = eigenvec[i][2];
		} else {
			eigenval2 = eigenval[1];
			for ( int i=0; i < 3; i++ ) eigenvec2[i] = eigenvec[i][1];
			if ( eigenval[0] >= eigenval[2] ) {
				eigenval0 = eigenval[0];
				for ( int i=0; i < 3; i++ ) eigenvec0[i] = eigenvec[i][0];
				eigenval1 = eigenval[2];
				for ( int i=0; i < 3; i++ ) eigenvec1[i] = eigenvec[i][2];
			} else {
				eigenval0 = eigenval[2];
				for ( int i=0; i < 3; i++ ) eigenvec0[i] = eigenvec[i][2];
				eigenval1 = eigenval[0];
				for ( int i=0; i < 3; i++ ) eigenvec1[i] = eigenvec[i][0];
			}
		}
	} else {
		if ( eigenval[2] >= eigenval[1] ) {
			eigenval0 = eigenval[2];
			for ( int i=0; i < 3; i++ ) eigenvec0[i] = eigenvec[i][2];
			eigenval1 = eigenval[1];
			for ( int i=0; i < 3; i++ ) eigenvec1[i] = eigenvec[i][1];
			eigenval2 = eigenval[0];
			for ( int i=0; i < 3; i++ ) eigenvec2[i] = eigenvec[i][0];
		} else {
			eigenval0 = eigenval[1];
			for ( int i=0; i < 3; i++ ) eigenvec0[i] = eigenvec[i][1];
			if ( eigenval[0] >= eigenval[2] ) {
				eigenval1 = eigenval[0];
				for ( int i=0; i < 3; i++ ) eigenvec1[i] = eigenvec[i][0];
				eigenval2 = eigenval[2];
				for ( int i=0; i < 3; i++ ) eigenvec2[i] = eigenvec[i][2];
			} else {
				eigenval1 = eigenval[2];
				for ( int i=0; i < 3; i++ ) eigenvec1[i] = eigenvec[i][2];
				eigenval2 = eigenval[0];
				for ( int i=0; i < 3; i++ ) eigenvec2[i] = eigenvec[i][0];
			}
		}
	}

	/*************************************************************************
	* Compute the singular vectors U and V:         *
	* The left singular vectors, U are the same as the eigenvectors, while  *
	* right singular vectors, V are computed from the left singular vectors.*
	*************************************************************************/

	double * U0 = eigenvec0;
	double * U1 = eigenvec1;
	double * U2 = eigenvec2;
	double V0[3], V1[3], V2[3];
	for ( int i=0; i < 3; i++ ) {
		V0[i] = Ct[i][0] * U0[0] + Ct[i][1] * U0[1] + Ct[i][2] * U0[2];
		V1[i] = Ct[i][0] * U1[0] + Ct[i][1] * U1[1] + Ct[i][2] * U1[2];
		V2[i] = Ct[i][0] * U2[0] + Ct[i][1] * U2[1] + Ct[i][2] * U2[2];
	}
	double mag; // get the magnitude of the singular vector to normalize it
	mag = sqrt(V0[0]*V0[0] + V0[1]*V0[1] + V0[2]*V0[2]);
	for ( int i=0; i < 3; i++ ) V0[i] /= mag;
	mag = sqrt(V1[0]*V1[0] + V1[1]*V1[1] + V1[2]*V1[2]);
	for ( int i=0; i < 3; i++ )  V1[i] /= mag;
	mag = sqrt(V2[0]*V2[0] + V2[1]*V2[1] + V2[2]*V2[2]);
	for ( int i=0; i < 3; i++ )  V2[i] /= mag;

#ifdef __DEBUG_RMSD__
	std::cout << "	[" << U0[0] << " " << U1[0] << " " << U2[0] << "]" << std::endl;
	std::cout << "U = [" << U0[1] << " " << U1[1] << " " << U2[1] << "]" << std::endl;
	std::cout << "	[" << U0[2] << " " << U1[2] << " " << U2[2] << "]" << std::endl;
	std::cout << std::endl;
	std::cout << "	[" << V0[0] << " " << V1[0] << " " << V2[0] << "]" << std::endl;
	std::cout << "V = [" << V0[1] << " " << V1[1] << " " << V2[1] << "]" << std::endl;
	std::cout << "	[" << V0[2] << " " << V1[2] << " " << V2[2] << "]" << std::endl;
	std::cout << std::endl;
#endif

	/*************************************************************************
	* Compute R=VUt, taking chirality into consideration     *
	*************************************************************************/

	double detV = V0[0]*(V1[1]*V2[2] - V2[1]*V1[2])
		- V1[0]*(V0[1]*V2[2] - V2[1]*V0[2])
		+ V2[0]*(V0[1]*V1[2] - V1[1]*V0[2]);
	double detU = U0[0]*(U1[1]*U2[2] - U2[1]*U1[2])
		- U1[0]*(U0[1]*U2[2] - U2[1]*U0[2])
		+ U2[0]*(U0[1]*U1[2] - U1[1]*U0[2]);
	int is_reflection = (detV * detU >= 0.)? 1: -1;
	if ( is_reflection == 1 ) {
		R[0][0] = V0[0]*U0[0] + V1[0]*U1[0] + V2[0]*U2[0];
		R[0][1] = V0[0]*U0[1] + V1[0]*U1[1] + V2[0]*U2[1];
		R[0][2] = V0[0]*U0[2] + V1[0]*U1[2] + V2[0]*U2[2];
		R[1][0] = V0[1]*U0[0] + V1[1]*U1[0] + V2[1]*U2[0];
		R[1][1] = V0[1]*U0[1] + V1[1]*U1[1] + V2[1]*U2[1];
		R[1][2] = V0[1]*U0[2] + V1[1]*U1[2] + V2[1]*U2[2];
		R[2][0] = V0[2]*U0[0] + V1[2]*U1[0] + V2[2]*U2[0];
		R[2][1] = V0[2]*U0[1] + V1[2]*U1[1] + V2[2]*U2[1];
		R[2][2] = V0[2]*U0[2] + V1[2]*U1[2] + V2[2]*U2[2];
	} else {
		R[0][0] = V0[0]*U0[0] + V1[0]*U1[0] - V2[0]*U2[0];
		R[0][1] = V0[0]*U0[1] + V1[0]*U1[1] - V2[0]*U2[1];
		R[0][2] = V0[0]*U0[2] + V1[0]*U1[2] - V2[0]*U2[2];
		R[1][0] = V0[1]*U0[0] + V1[1]*U1[0] - V2[1]*U2[0];
		R[1][1] = V0[1]*U0[1] + V1[1]*U1[1] - V2[1]*U2[1];
		R[1][2] = V0[1]*U0[2] + V1[1]*U1[2] - V2[1]*U2[2];
		R[2][0] = V0[2]*U0[0] + V1[2]*U1[0] - V2[2]*U2[0];
		R[2][1] = V0[2]*U0[1] + V1[2]*U1[1] - V2[2]*U2[1];
		R[2][2] = V0[2]*U0[2] + V1[2]*U1[2] - V2[2]*U2[2];
	}

#ifdef __DEBUG_RMSD__
	std::cout << "	[" <<R[0][0] << " " <<R[0][1] << " " <<R[0][2] << "]" << std::endl;
	std::cout << "R = [" <<R[1][0] << " " <<R[1][1] << " " <<R[1][2] << "]" << std::endl;
	std::cout << "	[" <<R[2][0] << " " <<R[2][1] << " " <<R[2][2] << "]" << std::endl;
	std::cout << std::endl;
#endif

	double residual = Eo - (double) sqrt(fabs(eigenval0))
		- (double) sqrt(fabs(eigenval1))
		- is_reflection * (double) sqrt(fabs(eigenval2));
	double rmsd = sqrt( fabs((residual)*2.0/n) );
	return rmsd;
}


// Same as above, but skips computation of R
double RMSD( std::vector<double> & coords1, std::vector<double> & coords2, int n )
{
	/*************************************************************************
	* Compute the centroids of coords1 and coords2 respectively    *
	*************************************************************************/

	double centroid1[3] = {0., 0., 0.};
	double centroid2[3] = {0., 0., 0.};
	for ( int j=0; j < n; j++ ) {
		int j3 = 3*j;
		centroid1[0] += coords1[j3  ];
		centroid2[0] += coords2[j3  ];
		centroid1[1] += coords1[j3+1];
		centroid2[1] += coords2[j3+1];
		centroid1[2] += coords1[j3+2];
		centroid2[2] += coords2[j3+2];
	}
	centroid1[0] /= n;
	centroid2[0] /= n;
	centroid1[1] /= n;
	centroid2[1] /= n;
	centroid1[2] /= n;
	centroid2[2] /= n;

	/*************************************************************************
	* Translate coords1 and coords2 such that their centroids are at the *
	* origin.                  *
	*************************************************************************/

	for ( int m=0; m < n; m++ ) {
		int m3 = 3*m;
		coords1[m3  ] -= centroid1[0];
		coords2[m3  ] -= centroid2[0];
		coords1[m3+1] -= centroid1[1];
		coords2[m3+1] -= centroid2[1];
		coords1[m3+2] -= centroid1[2];
		coords2[m3+2] -= centroid2[2];
	}

	/*************************************************************************
	* Compute C=Q(Pt) and the sum_{i=1 to n} p_i^2 + q_i^2 term (Eo)  *
	*************************************************************************/

	double C[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
	double Eo = 0.0;
	for ( int m=0; m < n; m++ ) {
		int m3 = 3*m;
		C[0][0] += coords2[m3  ] * coords1[m3  ];
		C[0][1] += coords2[m3  ] * coords1[m3+1];
		C[0][2] += coords2[m3  ] * coords1[m3+2];
		C[1][0] += coords2[m3+1] * coords1[m3  ];
		C[1][1] += coords2[m3+1] * coords1[m3+1];
		C[1][2] += coords2[m3+1] * coords1[m3+2];
		C[2][0] += coords2[m3+2] * coords1[m3  ];
		C[2][1] += coords2[m3+2] * coords1[m3+1];
		C[2][2] += coords2[m3+2] * coords1[m3+2];
		Eo += coords1[m3  ]*coords1[m3  ] + coords2[m3  ]*coords2[m3  ]
			+ coords1[m3+1]*coords1[m3+1] + coords2[m3+1]*coords2[m3+1]
			+ coords1[m3+2]*coords1[m3+2] + coords2[m3+2]*coords2[m3+2];
	}
	Eo *= 0.5;

#ifdef __DEBUG_RMSD__
	std::cout << "	[" <<C[0][0] << " " <<C[0][1] << " " <<C[0][2] << "]" << std::endl;
	std::cout << "C = [" <<C[1][0] << " " <<C[1][1] << " " <<C[1][2] << "]" << std::endl;
	std::cout << "	[" <<C[2][0] << " " <<C[2][1] << " " <<C[2][2] << "]" << std::endl;
	std::cout << std::endl;
#endif

	/*************************************************************************
	* Prepare CCt for eigenvector decomposition        *
	*************************************************************************/

	double Ct[3][3], CCt[3][3];
	for ( int i=0; i < 3; i++ ) {
		for ( int j=0; j < 3; j++ ) {
			Ct[i][j] = C[j][i];
		}
	}
	for ( int i=0; i<3; i++ ) {
		for ( int j=0; j < 3; j++ ) {
			CCt[i][j] = 0.;
			for ( int k = 0; k < 3; k++ ) {
				CCt[i][j] += C[i][k] * C[j][k]; // C[i][k]*Ct[k][j]
			}
		}
	}

	/*************************************************************************
	* Decompose CCt for its eigenvalues/eigenvectors, placed in    *
	* eigenval{0,1,2} and eigenvec{0,1,2} respectively.      *
	*************************************************************************/

	double eigenval[3], eigenvec[3][3];
	if ( !jacobi3(CCt, eigenval, eigenvec) ) return -1.;

	/* Sort the eigenvalues such that eigenval0 >= eigenval1 >= eigenval2 */
	double eigenval0, eigenval1, eigenval2;
	double eigenvec0[3], eigenvec1[3], eigenvec2[3];
	if ( eigenval[0] >= eigenval[1] ) {
		if ( eigenval[1] >= eigenval[2] ) {
			eigenval0 = eigenval[0];
			for ( int i=0; i < 3; i++ ) eigenvec0[i] = eigenvec[i][0];
			eigenval1 = eigenval[1];
			for ( int i=0; i < 3; i++ ) eigenvec1[i] = eigenvec[i][1];
			eigenval2 = eigenval[2];
			for ( int i=0; i < 3; i++ ) eigenvec2[i] = eigenvec[i][2];
		} else {
			eigenval2 = eigenval[1];
			for ( int i=0; i < 3; i++ ) eigenvec2[i] = eigenvec[i][1];
			if ( eigenval[0] >= eigenval[2] ) {
				eigenval0 = eigenval[0];
				for ( int i=0; i < 3; i++ ) eigenvec0[i] = eigenvec[i][0];
				eigenval1 = eigenval[2];
				for ( int i=0; i < 3; i++ ) eigenvec1[i] = eigenvec[i][2];
			} else {
				eigenval0 = eigenval[2];
				for ( int i=0; i < 3; i++ ) eigenvec0[i] = eigenvec[i][2];
				eigenval1 = eigenval[0];
				for ( int i=0; i < 3; i++ ) eigenvec1[i] = eigenvec[i][0];
			}
		}
	} else {
		if ( eigenval[2] >= eigenval[1] ) {
			eigenval0 = eigenval[2];
			for ( int i=0; i < 3; i++ ) eigenvec0[i] = eigenvec[i][2];
			eigenval1 = eigenval[1];
			for ( int i=0; i < 3; i++ ) eigenvec1[i] = eigenvec[i][1];
			eigenval2 = eigenval[0];
			for ( int i=0; i < 3; i++ ) eigenvec2[i] = eigenvec[i][0];
		} else {
			eigenval0 = eigenval[1];
			for ( int i=0; i < 3; i++ ) eigenvec0[i] = eigenvec[i][1];
			if ( eigenval[0] >= eigenval[2] ) {
				eigenval1 = eigenval[0];
				for ( int i=0; i < 3; i++ ) eigenvec1[i] = eigenvec[i][0];
				eigenval2 = eigenval[2];
				for ( int i=0; i < 3; i++ ) eigenvec2[i] = eigenvec[i][2];
			} else {
				eigenval1 = eigenval[2];
				for ( int i=0; i < 3; i++ ) eigenvec1[i] = eigenvec[i][2];
				eigenval2 = eigenval[0];
				for ( int i=0; i < 3; i++ ) eigenvec2[i] = eigenvec[i][0];
			}
		}
	}

	/*************************************************************************
	* Compute the singular vectors U and V:         *
	* The left singular vectors, U are the same as the eigenvectors, while  *
	* right singular vectors, V are computed from the left singular vectors.*
	*************************************************************************/

	double * U0 = eigenvec0;
	double * U1 = eigenvec1;
	double * U2 = eigenvec2;
	double V0[3], V1[3], V2[3];
	for ( int i=0; i < 3; i++ ) {
		V0[i] = Ct[i][0] * U0[0] + Ct[i][1] * U0[1] + Ct[i][2] * U0[2];
		V1[i] = Ct[i][0] * U1[0] + Ct[i][1] * U1[1] + Ct[i][2] * U1[2];
		V2[i] = Ct[i][0] * U2[0] + Ct[i][1] * U2[1] + Ct[i][2] * U2[2];
	}
	double mag; // get the magnitude of the singular vector to normalize it
	mag = sqrt(V0[0]*V0[0] + V0[1]*V0[1] + V0[2]*V0[2]);
	for ( int i=0; i < 3; i++ ) {
		V0[i] /= mag;
	}
	mag = sqrt(V1[0]*V1[0] + V1[1]*V1[1] + V1[2]*V1[2]);
	for ( int i=0; i < 3; i++ ) {
		V1[i] /= mag;
	}
	mag = sqrt(V2[0]*V2[0] + V2[1]*V2[1] + V2[2]*V2[2]);
	for ( int i=0; i < 3; i++ ) {
		V2[i] /= mag;
	}

#ifdef __DEBUG_RMSD__
	std::cout << "	[" << U0[0] << " " << U1[0] << " " << U2[0] << "]" << std::endl;
	std::cout << "U = [" << U0[1] << " " << U1[1] << " " << U2[1] << "]" << std::endl;
	std::cout << "	[" << U0[2] << " " << U1[2] << " " << U2[2] << "]" << std::endl;
	std::cout << std::endl;
	std::cout << "	[" << V0[0] << " " << V1[0] << " " << V2[0] << "]" << std::endl;
	std::cout << "V = [" << V0[1] << " " << V1[1] << " " << V2[1] << "]" << std::endl;
	std::cout << "	[" << V0[2] << " " << V1[2] << " " << V2[2] << "]" << std::endl;
	std::cout << std::endl;
#endif

	/*************************************************************************
	* Compute R=VUt, taking chirality into consideration     *
	*************************************************************************/

	double detV = V0[0]*(V1[1]*V2[2] - V2[1]*V1[2])
		- V1[0]*(V0[1]*V2[2] - V2[1]*V0[2])
		+ V2[0]*(V0[1]*V1[2] - V1[1]*V0[2]);
	double detU = U0[0]*(U1[1]*U2[2] - U2[1]*U1[2])
		- U1[0]*(U0[1]*U2[2] - U2[1]*U0[2])
		+ U2[0]*(U0[1]*U1[2] - U1[1]*U0[2]);
	int is_reflection = (detV * detU >= 0.)? 1: -1;

#ifdef __DEBUG_RMSD__
	std::cout << "	[" <<R[0][0] << " " <<R[0][1] << " " <<R[0][2] << "]" << std::endl;
	std::cout << "R = [" <<R[1][0] << " " <<R[1][1] << " " <<R[1][2] << "]" << std::endl;
	std::cout << "	[" <<R[2][0] << " " <<R[2][1] << " " <<R[2][2] << "]" << std::endl;
	std::cout << std::endl;
#endif

	double residual = Eo - (double) sqrt(fabs(eigenval0))
		- (double) sqrt(fabs(eigenval1))
		- is_reflection * (double) sqrt(fabs(eigenval2));
	double rmsd = sqrt( fabs((residual)*2.0/n) );
	return rmsd;
}



void rotate( std::vector<double> & coords, int n, double R[3][3], double * result )
{
	for ( int m=0; m < n; m++ ) {
		int m3 = 3*m;
		result[m3  ] = R[0][0]*coords[m3  ]
			+ R[0][1]*coords[m3+1]
			+ R[0][2]*coords[m3+2];
		result[m3+1] = R[1][0]*coords[m3  ]
			+ R[1][1]*coords[m3+1]
			+ R[1][2]*coords[m3+2];
		result[m3+2] = R[2][0]*coords[m3  ]
			+ R[2][1]*coords[m3+1]
			+ R[2][2]*coords[m3+2];
	}
}


//- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -

double fast_rmsd( std::vector<double> & coords1, std::vector<double> & coords2, int n)
{
	/*************************************************************************
	* Compute the centroids of coords1 and coords2 respectively    *
	*************************************************************************/

	double centroid1[3] = {0., 0., 0.};
	double centroid2[3] = {0., 0., 0.};
	for ( int j=0; j < n; j++ ) {
		int j3 = 3*j;
		centroid1[0] += coords1[j3  ];
		centroid2[0] += coords2[j3  ];
		centroid1[1] += coords1[j3+1];
		centroid2[1] += coords2[j3+1];
		centroid1[2] += coords1[j3+2];
		centroid2[2] += coords2[j3+2];
	}
	centroid1[0] /= n;
	centroid2[0] /= n;
	centroid1[1] /= n;
	centroid2[1] /= n;
	centroid1[2] /= n;
	centroid2[2] /= n;

	/*************************************************************************
	* Translate coords1 and coords2 such that their centroids are at the *
	* origin.                  *
	*************************************************************************/

	for ( int m=0; m < n; m++ ) {
		int m3 = 3*m;
		coords1[m3  ] -= centroid1[0];
		coords2[m3  ] -= centroid2[0];
		coords1[m3+1] -= centroid1[1];
		coords2[m3+1] -= centroid2[1];
		coords1[m3+2] -= centroid1[2];
		coords2[m3+2] -= centroid2[2];
	}

	/*************************************************************************
	* Compute C=Q(Pt) and the sum_{i=1 to n} p_i^2 + q_i^2 term (Eo)  *
	*************************************************************************/

	double C[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
	double Eo = 0.0;
	for ( int m=0; m < n; m++ ) {
		int m3 = 3*m;
		C[0][0] += coords2[m3  ] * coords1[m3  ];
		C[0][1] += coords2[m3  ] * coords1[m3+1];
		C[0][2] += coords2[m3  ] * coords1[m3+2];
		C[1][0] += coords2[m3+1] * coords1[m3  ];
		C[1][1] += coords2[m3+1] * coords1[m3+1];
		C[1][2] += coords2[m3+1] * coords1[m3+2];
		C[2][0] += coords2[m3+2] * coords1[m3  ];
		C[2][1] += coords2[m3+2] * coords1[m3+1];
		C[2][2] += coords2[m3+2] * coords1[m3+2];
		Eo += coords1[m3  ]*coords1[m3  ] + coords2[m3  ]*coords2[m3  ]
			+ coords1[m3+1]*coords1[m3+1] + coords2[m3+1]*coords2[m3+1]
			+ coords1[m3+2]*coords1[m3+2] + coords2[m3+2]*coords2[m3+2];
	}
	Eo *= 0.5;

	/*************************************************************************
	* Determine the chirality of the correlation matrix C.      *
	* That is, see if det(C) > 0.             *
	*************************************************************************/

	double detC = C[0][0] * (C[1][1]*C[2][2] - C[1][2]*C[2][1])
		- C[0][1] * (C[1][0]*C[2][2] - C[1][2]*C[2][0])
		+ C[0][2] * (C[1][0]*C[2][1] - C[1][1]*C[2][0]);
	int omega = (detC > 0.0)? 1: -1;

	/*************************************************************************
	*     [ 1  e01 e02 ]          *
	* Compute CtC/x = [ *  e11 e12 ]          *
	*     [ * * e22 ]          *
	*************************************************************************/

	double x   =  C[0][0]*C[0][0] + C[1][0]*C[1][0] + C[2][0]*C[2][0];
	double e01 = (C[0][1]*C[0][1] + C[1][1]*C[1][1] + C[2][1]*C[2][1]) / x;
	double e02 = (C[0][2]*C[0][2] + C[1][2]*C[1][2] + C[2][2]*C[2][2]) / x;

	double e11 = (C[0][0]*C[0][1] + C[1][0]*C[1][1] + C[2][0]*C[2][1]) / x;
	double e12 = (C[0][1]*C[0][2] + C[1][1]*C[1][2] + C[2][1]*C[2][2]) / x;

	double e22 = (C[0][0]*C[0][2] + C[1][0]*C[1][2] + C[2][0]*C[2][2]) / x;

	/*************************************************************************
	* Obtain the eigenvalues of CtC through solving the equation   *
	*                    *
	*   det( CtC-I*lambda ) = 0            *
	*                    *
	* This gives us               *
	*   lambda^3 + a2*lambda^2 + a1*lambda + a0 = 0       *
	* where...                 *
	*************************************************************************/

	double a2 = -1.0 - e01 - e02;
	double a1 = e01 + e02 + e01*e02 - e11*e11 - e12*e12 - e22*e22;
	double a0 = e11*e11*e02 + e12*e12 + e22*e22*e01 - e01*e02 - 2*e11*e22*e12;
	double roots[3];
	cubic_roots(a2, a1, a0, roots);
	roots[0] *= x;
	roots[1] *= x;
	roots[2] *= x;
	int min_index = 2; // place smallest eigenvalue in roots[2]
	if ( roots[0] < roots[2] ) {
		min_index = (roots[0] < roots[1])? 0: 1;
	} else {
		if ( roots[1] < roots[2] ) {
			min_index = 1;
		}
	}
	if ( min_index != 2 ) { // then, we swap roots[min_index] with roots[2]
		double aux = roots[min_index];
		roots[min_index] = roots[2];
		roots[2] = aux;
	}

	double residual = Eo - sqrt(roots[0]) - sqrt(roots[1])
		- omega * sqrt(roots[2]);
	return sqrt( residual*2.0 / n );
}

}
}
}

