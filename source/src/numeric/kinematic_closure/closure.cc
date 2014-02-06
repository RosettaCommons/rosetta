// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/* This file is contains all the same functions as `bridgeObjects.cc', but 
 * radians are assumed for all angles.  Although this is easier to use in most 
 * cases, it is not backwards compatible.  That's why a new file was created.  
 */

/// @file   bridge_objects.cc
/// @brief  Solves the triaxial loop closure problem for a system of atoms
/// @author Evangelos A. Coutsias
/// @author Daniel J. Mandell

// Unit Headers
#include <numeric/types.hh>
#include <numeric/constants.hh>
#include <numeric/wrap_angles.hh>
#include <numeric/conversions.hh>
#include <numeric/kinematic_closure/dixon.hh>
#include <numeric/kinematic_closure/vector.hh>
#include <numeric/kinematic_closure/kinematic_closure_helpers.hh>

// Utility headers
#include <utility/vector1.hh>

// C++ headers
#include <cmath>
#include <iostream>

namespace numeric {
namespace kinematic_closure {
namespace radians {

using namespace std;
using numeric::constants::r::pi;

void triangle ( // {{{1
		const utility::vector1<Real>& vbond,
		utility::vector1<Real>& calpha,
		utility::vector1<Real>& salpha) {

  Real d11, d12, d13, d22, d23, d33;
  d11 =   vbond[1]*vbond[1];
  d22 =   vbond[2]*vbond[2];
  d33 =   vbond[3]*vbond[3];
  d12 = 2*vbond[1]*vbond[2];
  d13 = 2*vbond[1]*vbond[3];
  d23 = 2*vbond[2]*vbond[3];
  calpha[1]=(d22-d33-d11)/d13; // exterior angle
  salpha[1]=std::sqrt(1-calpha[1]*calpha[1]);
  calpha[2]=(d33-d11-d22)/d12; // exterior angle
  salpha[2]=std::sqrt(1-calpha[2]*calpha[2]);
  calpha[3]=(d11-d22-d33)/d23; // exterior angle
  salpha[3]=std::sqrt(1-calpha[3]*calpha[3]);
}

void sincos ( // {{{1
		const utility::vector1<Real>& theta,
		const int flag,
		utility::vector1<Real>& cosine,
		utility::vector1<Real>& sine) {

  utility::vector1<Real> a, aa, a1;
  if (flag == 1) {
    for (int i=1; i<=3; i++) {
      cosine[i]=std::cos(theta[i]);
      sine[i]=std::sin(theta[i]);
    }
  }
  else {
    a.resize(3);
    aa.resize(3);
    a1.resize(3);
    for (int i=1; i<=3; i++) {
      a[i]=std::tan(theta[i]/2);
      aa[i]=a[i]*a[i];
      a1[i]=1+aa[i];
      sine[i]=2*a[i]/a1[i];
      cosine[i]=(1-aa[i])/a1[i];
    }
  }
}

// triaxialCoefficients {{{1
/*
 * triaxialCoefficients
 * --------------------
 * Sets up characteristic polynomial matrices A,B,C,D from tripeptide parameters
 */

void triaxialCoefficients(
		const utility::vector1<Real>& vb, 
		const utility::vector1<Real>& xi, 
		const utility::vector1<Real>& eta, 
		const utility::vector1<Real>& delta, 
		const utility::vector1<Real>& theta, 
		const utility::vector1<int>& order, 
		utility::vector1<utility::vector1<Real> >& A, 
		utility::vector1<utility::vector1<Real> >& B, 
		utility::vector1<utility::vector1<Real> >& C, 
		utility::vector1<utility::vector1<Real> >& D, 
		utility::vector1<Real>& cal, 
		utility::vector1<Real>& sal, 
		bool& feasible_triangle) {

  utility::vector1<Real> ctheta (4), calpha (4), salpha (4), cxi (3), sxi (3), ceta (4), seta (4), cdelta (3), sdelta (3), caceta (4), caseta (4), saceta (4), saseta (4), capeta (4), sapeta (4), cameta (4), sameta (4);
  Real am, ap, cdsx, sdsx, b, bm, bp, d, dp, dm;
  utility::vector1<utility::vector1<Real> > L, M, N;
  utility::vector1<utility::vector1<utility::vector1<Real> > > p (4);

  // allocate p
  for (int i=1; i<=4; i++) {
    p[i].resize(3);
    for (int j=1; j<=3; j++) {
      p[i][j].resize(3);
    }
  }

  // allocate polynomial matrices
  A.resize(3);
  B.resize(3);
  C.resize(3);
  D.resize(3);
  for (int i=1; i<=3; i++) {
    for (int j=1; j<=3; j++) {
      A[i].resize(3);
      B[i].resize(3);
      C[i].resize(3);
      D[i].resize(3);
    }
  }

  // if triangle is feasible, compute the polynomial matrices
  if ((std::abs(vb[1] - vb[2]) <= vb[3]) && (vb[1] + vb[2] >= vb[3])) {
    feasible_triangle = true; // we have a feasible triangle

    for (int i=1; i<=3; i++) {
      ctheta[i]=std::cos(theta[i]);
    }
    triangle(vb, calpha, salpha);
		cal.resize(3);
		sal.resize(3);
		for (int i=1; i<=3; i++) {
			cal[i]=calpha[i];
			sal[i]=salpha[i];
		}
    sincos(xi, 2, cxi, sxi);
    sincos(eta, 2, ceta, seta);
    sincos(delta, 2, cdelta, sdelta);
    calpha[4] = calpha[1];
    salpha[4] = salpha[1];
    ceta[4] = ceta[1];
    seta[4] = seta[1];
    ctheta[4] = ctheta[1];
    for (int i=1; i<=4; i++) {
      caceta[i] = calpha[i] * ceta[i];
      caseta[i] = calpha[i] * seta[i];
      saceta[i] = salpha[i] * ceta[i];
      saseta[i] = salpha[i] * seta[i];
      capeta[i] = caceta[i] - saseta[i];
      sapeta[i] = saceta[i] + caseta[i];
      cameta[i] = caceta[i] + saseta[i];
      sameta[i] = saceta[i] - caseta[i];
    }
    for (int i=2; i<=4; i++) {
      am    = cxi[i-1]*cameta[i] + ctheta[i];
      ap    = cxi[i-1]*capeta[i] + ctheta[i];
      cdsx  = cdelta[i-1]*sxi[i-1];
      sdsx  = 2*sdelta[i-1]*sxi[i-1];
      b     = 4*cdsx*seta[i];
      bm    = cdsx*sameta[i];
      bp    = cdsx*sapeta[i];
      d     = sdsx*seta[i];
      dp    = sdsx*sapeta[i];
      dm    = sdsx*sameta[i];
      p[i][3][3]= -am -bm;
      p[i][3][2]=     -d;
      p[i][3][1]= -ap -bp;
      p[i][2][3]=     -dm;
      p[i][2][2]=      b;
      p[i][2][1]=     -dp;
      p[i][1][3]= -am +bm;
      p[i][1][2]=      d;
      p[i][1][1]= -ap +bp;
    }
    for (int i=1; i<=3; i++) {
      for (int j=1; j<=3; j++) {
		p[1][i][j]=p[4][i][j];
      }
    }

    L.resize(3);
    M.resize(3);
    N.resize(3);
    for (int i=1; i<=3; i++) {
      L[i].resize(3);
      M[i].resize(3);
      N[i].resize(3);
    }
    for (int k=1; k<=3; k++) {
      for (int i=1; i<=3; i++) {
		L[k][i] = p[order[1]][i][k];
		M[k][i] = p[order[2]][k][i];
		N[k][i] = p[order[3]][k][i];
      }
    }

    for (int i=1; i<=3; i++) {
      for (int j=1; j<=3; j++) {
		// not using transpose indexing so we don't take transpose in dixon
		A[i][j] = M[i][2]*N[1][j] - M[i][1]*N[2][j];
        B[i][j] = M[i][3]*N[1][j] - M[i][1]*N[3][j];
        C[i][j] = M[i][3]*N[2][j] - M[i][2]*N[3][j];
        D[i][j] = L[i][j];
      }
    }
  }
  else {
    feasible_triangle = false; // we have an unfeasible triangle
    for (int i=1; i<=3; i++) {
      for (int j=1; j<=3; j++) {
		A[i][j]=0;
		B[i][j]=0;
		C[i][j]=0;
		D[i][j]=0;
      }
    }
  }
}

// void cross {{{1
/*
 * cross
 * -----
 * Output r is cross product of vectors L and r0
 */
void cross(const utility::vector1<Real>& L, const utility::vector1<Real>& r0, utility::vector1<Real>& r) {
	r.resize(3);
	r[1]=L[2]*r0[3]-L[3]*r0[2];
	r[2]=L[3]*r0[1]-L[1]*r0[3];
	r[3]=L[1]*r0[2]-L[2]*r0[1];
}

// void frame {{{1
/*
 * frame
 * -----
 * construct body frame for r1 origin
 *                          r2 r1-->r2: positive x axis
 *                          r3 on pos. y side of xy plane
 */
void frame(const utility::vector1<utility::vector1<Real> >& R, utility::vector1<utility::vector1<Real> >& U) {
	Real norm1, norm3;

// 	utility::vector1<Real> cross1 (3), cross3 (3), cross11 (3), cross12 (3), cross31 (3), cross32 (3);
	utility::vector1<Real> dR3R2 (3);
 	U.resize(3);
 	for (int i=1; i<=3; i++) {
 		U[i].resize(3);
 	}
 	for (int i=1; i<=3; i++) {
 		U[1][i]=R[2][i]-R[1][i];
 	}
	norm1=std::sqrt(U[1][1]*U[1][1] + U[1][2]*U[1][2] + U[1][3]*U[1][3]);
 	for (int i=1; i<=3; i++) {
 		U[1][i] /= norm1;
		dR3R2[i]=R[3][i]-R[2][i];
	}
	cross(U[1], dR3R2, U[3]);
 	norm3=std::sqrt(U[3][1]*U[3][1] + U[3][2]*U[3][2] + U[3][3]*U[3][3]);
 	for (int i=1; i<=3; i++) {
 		U[3][i] /= norm3;
	}
	cross(U[3], U[1], U[2]);
}

Real eucDistance( // {{{1
		const utility::vector1<Real>& a,
		const utility::vector1<Real>& b) {

	return std::sqrt(std::pow(a[1]-b[1],2) + std::pow(a[2]-b[2],2) + std::pow(a[3]-b[3],2));
}

Real scpn( // {{{1
		const utility::vector1<Real>& a,
		const utility::vector1<Real>& b,
		const utility::vector1<Real>& c) {

	Real d=0;
	for (int i=1; i<=3; i++) {
		d += (a[i]-b[i]) * (c[i]-b[i]);
	}
	return d;
}

// Real bondangle( // {{{1
/// @details @returns the bond angle between the given vectors in radians.
Real bondangle(
		const utility::vector1<Real>& a,
		const utility::vector1<Real>& b,
		const utility::vector1<Real>& c) {

	Real r = scpn(a,b,c) / (eucDistance(a,b) * eucDistance(b,c));
	Real ang = std::atan2(sqrt(1-r*r), (Real) r);
	return ang;
}

// Real torsion( // {{{1
/// @details @returns the torsion angle between the given vectors in radians.
Real torsion(
		const utility::vector1<Real>& a,
		const utility::vector1<Real>& b,
		const utility::vector1<Real>& c,
		const utility::vector1<Real>& d) {

	Real f, y, z, chi;
	utility::vector1<Real> r (3), sc1 (3), sc2 (3), sc3 (3), cs12, cs31;
	utility::vector1<utility::vector1<Real> > s (3);
	for (int i=1; i<=3; i++) {
		s[i].resize(3);
	}

	if (b==c) {
		return numeric::constants::r::pi;
	}

	for (int i=1; i<=3; i++) {
		r[i]    = d[i] - b[i];
		s[i][1] = c[i] - b[i];
		s[i][2] = a[i] - b[i];
		sc1[i]  = s[i][1];
		sc2[i]  = s[i][2];
	}
	cross(sc1, sc2, cs12);
	for (int i=1; i<=3; i++) {
		s[i][3] = cs12[i];
		sc3[i]  = s[i][3];
	}
	cross(sc3, sc1, cs31);
	for (int i=1; i<=3; i++) {
		s[i][2] = cs31[i];
	}
	for (int i=2; i<=3; i++) {
		f = std::sqrt( pow(s[1][i],2) + pow(s[2][i],2) + pow(s[3][i],2));
		s[1][i] /= f;
		s[2][i] /= f;
		s[3][i] /= f;
	}

	y = r[1]*s[1][2] + r[2]*s[2][2] + r[3]*s[3][2];
	z = r[1]*s[1][3] + r[2]*s[2][3] + r[3]*s[3][3];

	chi = atan2(z,y);
	return wrap_2pi(chi);
}

void chainParams( // {{{1
		const int& n,
		const utility::vector1<utility::vector1<Real> >& atoms,
		Real& vbond,
		Real& xi,
		Real& eta,
		Real& delta,
		utility::vector1<Real>& R0,
		utility::vector1<utility::vector1<Real> >& Q) {

	//utility::vector1<Real> ac1 (3), ac2 (3), acn (3), acneg1 (3);
	utility::vector1<utility::vector1<Real> > a1n2 (3);
// R0.resize(3);
// for (int i=1; i<=3; i++) {
// 		ac1[i] = atoms[i][1];
// 		ac2[i] = atoms[i][2];
// 		acn[i] = atoms[i][n];
// 		acneg1[i]= atoms[i][n-1];
// 		R0[i]  = ac1[i];
// 	}
 	a1n2[1]=atoms[1];
 	a1n2[2]=atoms[n];
 	a1n2[3]=atoms[2];

// 	printVector(acn);
// 	printVector(ac1);

	R0=atoms[1];
	frame(a1n2, Q);
	vbond=eucDistance(atoms[n], atoms[1]);
	xi=bondangle(atoms[n-1], atoms[n], atoms[1]);
	eta=bondangle(atoms[n], atoms[1], atoms[2]);
	delta=torsion(atoms[n-1], atoms[n], atoms[1], atoms[2]);

// 	frame(a1n2, Q);
// 	vbond=eucDistance(acn, ac1);
// 	xi=bondangle(acneg1, acn, ac1);
// 	eta=bondangle(acn, ac1, ac2);
// 	delta=torsion(acneg1, acn, ac1, ac2);
}


// void chainXYZ {{{1
/* Generates a set of cartesian coordinates (atoms) given the supplied bond angles, bond lengths, and torsions
 * of a peptide chain
 */
void chainXYZ(const int& n,
				const utility::vector1<Real>& b_len1,
				const utility::vector1<Real>& b_ang1,
				const utility::vector1<Real>& t_ang1,
				const bool space,
				const utility::vector1<Real>& R0,
				const utility::vector1<utility::vector1<Real> >& Q,
				utility::vector1<utility::vector1<Real> >& atoms
				)
{
	Real ca, sa, ct, st; // cos(b_ang[2]), sin(b_ang[2]), cos(t_ang[2]), sin(t_ang[2]),
	Real ba;
	int n1, n2, j1, j2;
	utility::vector1<Real> b_len, b_ang (n), t_ang (n);//, cos_a (n), cos_t (n), sin_a (n), sin_t (n);
	utility::vector1<Real> v(3), s(3);
	utility::vector1<utility::vector1<Real> > U (3); // rotation matrix

	atoms.resize(n);
	b_len = b_len1; // only defined for consistent naming
	// convert angles to radians and allocate atom vectors
	for (int i=1; i<=n; i++) {
		t_ang[i] = t_ang1[i];
		b_ang[i] = b_ang1[i];
		atoms[i].resize(3);
	}
	// place the first atom at the origin
	for (int i=1; i<=n; i++) {
		for (int j=1; j<=3; j++) {
			atoms[i][j]=0.0;
		}
	}
	// place the second atom at bond length along x-axis from first atom
	atoms[2][1] = b_len[1];
	// third atom shows that b_ang is interior
	ca = std::cos(b_ang[2]);
	sa = std::sin(b_ang[2]);
	atoms[3][1] = b_len[1] - b_len[2] * ca;
	atoms[3][2] =            b_len[2] * sa;

	// initialize rotation matrix
	for (int i=1; i<=3; i++) {
		U[i].resize(3);
	}
	U[1][1] =    -ca;
	U[2][1] =    -sa;
	U[1][2] =     sa;
	U[2][2] =    -ca;
	for (int i=1; i<=2; i++) {
		U[3][i] = 0;
		U[i][3] = 0;
	}
	U[3][3] = 1;

	// all other atoms
	n1 = n-1;
	n2 = n-2;
	j1 = 3;
	j2 = 2;
	for (int j=4; j<=n1; j++) {
		ca = std::cos(b_ang[j1]);
		sa = std::sin(b_ang[j1]);
		ct = std::cos(t_ang[j2]);
		st = std::sin(t_ang[j2]);
		for (int i=1; i<=3; i++) {
			s[i] = U[2][i];
		}
		for (int i=1; i<=3; i++) {
			U[2][i] = s[i]*ct + U[3][i]*st;
		}
		for (int i=1; i<=3; i++) {
			U[3][i] = -s[i]*st + U[3][i]*ct;
		}
		for (int i=1; i<=3; i++) {
			s[i] = U[2][i];
		}
		for (int i=1; i<=3; i++) {
			U[2][i] = -s[i]*ca - U[1][i]*sa;
		}
		for (int i=1; i<=3; i++) {
			U[1][i] = s[i]*sa - U[1][i]*ca;
		}
		for (int i=1; i<=3; i++) {
			atoms[j][i] = atoms[j1][i] + b_len[j1] * U[1][i];
		}
		j2 = j1;
		j1 = j;
	}
	ca = std::cos(b_ang[n1]);
	sa = std::sin(b_ang[n1]);
	ct = std::cos(t_ang[n2]);
	st = std::sin(t_ang[n2]);
	v[1] =  -b_len[n1] * ca;
	ba   =   b_len[n1] * sa;
	v[2] =   ba * ct;
	v[3] =   ba * st;
	// implementing s = U*v'
	for (int i=1; i<=3; i++) {
		s[i] = (U[1][i] * v[1]) + (U[2][i] * v[2]) + (U[3][i] * v[3]);
	}
	for (int i=1; i<=3; i++) {
		atoms[n][i] = atoms[n1][i] + s[i];
	}
	// check atoms here but still need to implement space and frame
	if (space == true) {
		// implement:
		// for i = 1:n
		//   atoms(:,i) = R0 + Q*atoms(:,i);
		// end
		utility::vector1<Real> spaceatoms (3); // holder for inner-loop coordinates
		for (int ind=1; ind<=n; ind++) {
			for (int i=1; i<=3; i++) {
				Real inner_sum=0.0;
				for (int j=1; j<=3; j++) {
					inner_sum += (Q[i][j] * atoms[ind][j]);
				}
				spaceatoms[i]=R0[i]+inner_sum;
			}
			for (int k=1; k<=3; k++) {
				atoms[ind][k]=spaceatoms[k];
			}
		}
	}
	/*
	else if (index[1] != 1 || index[2] != 2 || index[3] != 3) {
		utility::vector1<Real> Ra (3);
		utility::vector1<utility::vector1<Real> > RR (3), Qa (3), body_atoms(n);
		//printMatrix(atoms);
		for (int i=1; i<=3; i++) {
			RR[i].resize(3);
			Ra[i] = -atoms[index[1]][i];
		}
		//printVector(Ra);
		for (int i=1; i<=3; i++) {
			for (int j=1; j<=3; j++) {
				RR[j][i]=atoms[index[j]][i];
			}
		}
		//printMatrix(RR);
		frame(RR, Qa);
		//printMatrix(Qa);
		for (int ind=1; ind<=n; ind++) {
			// implement atoms(:,i) = Qa'*(atoms(:,i)+Ra);
			body_atoms[ind].resize(3);
			for (int i=1; i<=3; i++) {
				Real inner_sum=0.0;
				for (int j=1; j<=3; j++) {
					inner_sum += Qa[i][j] * (atoms[ind][j] + Ra[j]);
				}
				body_atoms[ind][i] = inner_sum;
			}
		}
		atoms = body_atoms;
		//exit(0);
	}
	*/
}


// void chainXYZ {{{1

void chainXYZ(
		const int& n,
		const utility::vector1<Real>& b_len,
		const utility::vector1<Real>& b_ang,
		const utility::vector1<Real>& t_ang,
		utility::vector1<utility::vector1<Real> >& atoms) {

	utility::vector1<Real> dummy_R0;
	utility::vector1<utility::vector1<Real> > dummy_Q;

	chainXYZ(n, b_len, b_ang, t_ang, false, dummy_R0, dummy_Q, atoms);
}

// void chainTORS {{{1
/*        chainTORS
//        ---------
//        calculates torsions, bond angles, bond lengths from atoms
//        for a linear chain
//        calculates angles in degrees
//
//        (j-1) --------  (j)
//                 b_ang(j) \
//                  t_ang(j) \b_len(j)
//                            \
//                            (j+1) ------  (j+2)
//
*/

void chainTORS (
		const int& n,
		const utility::vector1<utility::vector1<Real> >& atoms,
		utility::vector1<Real>& t_ang,
		utility::vector1<Real>& b_ang,
		utility::vector1<Real>& b_len,
		utility::vector1<Real>& R0,
		utility::vector1<utility::vector1<Real> >& Q) {

	t_ang.resize(n);
	b_ang.resize(n);
	b_len.resize(n);
	utility::vector1<utility::vector1<Real> > framein (3);
	for (int i=1; i<=3; i++) {
		framein[i]=atoms[i];
	}
	R0=atoms[1];
	frame(framein,Q);
	for (int j=1; j <= n-1; j++) {
		b_len[j] = eucDistance( atoms[j], atoms[j+1] );
	}
	b_len[n] = eucDistance( atoms[n], atoms[1] );
	for (int j=2; j <= n-1; j++) {
		b_ang[j] = bondangle( atoms[j-1], atoms[j], atoms[j+1] );
	}
	b_ang[n] = 0; //bondangle( atoms[n-1], atoms[n], atoms[1] );
	b_ang[1] = 0; //bondangle( atoms[n],   atoms[1], atoms[2] );

	if (n >= 4) {
		for (int j=2; j <= n-2; j++) {
			t_ang[j] = torsion( atoms[j-1], atoms[j], atoms[j+1], atoms[j+2] );
		}
		t_ang[n-1] = 0; //torsion( atoms[n-2], atoms[n-1], atoms[n], atoms[1] );
		t_ang[n]   = 0; //torsion( atoms[n-1], atoms[n],   atoms[1], atoms[2] );
		t_ang[1]   = 0; //torsion( atoms[n]  , atoms[1],   atoms[2], atoms[3] );
	}
	else {
		for (int j=1; j<=n; j++) {
			t_ang[j]=0.0;
		}
	}
}

void rotateX( // {{{1
		const utility::vector1<utility::vector1<Real> >& R, 
		const Real& c,
		const Real& s,
		utility::vector1<utility::vector1<Real> >& Rx) {

	int dim1=R.size();
	Rx.resize(dim1);
	for (int i=1; i<=dim1; i++) {
		Rx[i].resize(R[i].size());
		Rx[i][1] =    R[i][1];
		Rx[i][2] =  c*R[i][2]  - s*R[i][3];
		Rx[i][3] =  s*R[i][2]  + c*R[i][3];
	}
}

void rotateY( // {{{1
		const utility::vector1<utility::vector1<Real> >& R,
		const Real& c,
		const Real& s,
		utility::vector1<utility::vector1<Real> >& Ry) {

	int dim1=R.size();
	Ry.resize(dim1);
	for (int i=1; i<=dim1; i++) {
		Ry[i].resize(R[i].size());
		Ry[i][1] =   c*R[i][1]  + s*R[i][3];
		Ry[i][2] =     R[i][2];
		Ry[i][3] =  -s*R[i][1]  + c*R[i][3];
	}
}

void rotateZ( // {{{1
		const utility::vector1<utility::vector1<Real> >& R,
		const Real& c,
		const Real& s,
		utility::vector1<utility::vector1<Real> >& Rz) {

	int dim1=R.size();
	Rz.resize(dim1);
	for (int i=1; i<=dim1; i++) {
		Rz[i].resize(R[i].size());
		Rz[i][1] =   c*R[i][1]  - s*R[i][2];
		Rz[i][2] =   s*R[i][1]  + c*R[i][2];
		Rz[i][3] =     R[i][3];
	}
}

void to_radians ( // {{{1
		utility::vector1<Real> & degrees) {

	using numeric::conversions::to_radians;
	for (Size i = 1; i < degrees.size(); i++) {
		to_radians(degrees[i]);
	}
}
// }}}1
void to_degrees ( // {{{1
		utility::vector1<Real> & radians) {

	using numeric::conversions::to_degrees;
	for (Size i = 1; i < radians.size(); i++) {
		to_degrees(radians[i]);
	}
}
// }}}1

// void bridge_objects ( // {{{1
////////////////////////////////////////////////////////////////////////////////
/// @begin bridge_objects
///
/// @brief
/// Solve the triaxial loop closure problem for a system of atoms
///
/// @detailed
///
/// @param[in] atoms - matrix of cartesian coordinates of N-CA-C atoms indexed 
/// as atoms[atom][dimension]
/// @param[in] dt - desired torsions for each atom (radians)
/// @param[in] da - desired bond angle for each atom (radians)
/// @param[in] db - desired bond length for each atom (radians)
/// @param[in] pivots - 3 indices (base 1) of atoms to be used as loop closure pivots
/// @param[in] order - length 3 vector giving order to solve for the tau parameters (use 1,2,3 if unsure)
/// @param[out] t_ang - matrix giving torsion angles for each atom for each solution, indexed as t_ang[solution][atom] (radians)
/// @param[out] b_ang - matrix giving bond angles for each atom for each solution, indexed as b_ang[solution][atom] (radians)
/// @param[out] b_len - matrix giving bond lengths for each atom for each solution, indexed as b_len[solution][atom] (radians)
/// @param[out] nsol - number of solutions found
///
/// @global_read
///
/// @global_write
///
/// @return
///
/// @remarks
/// dt, da, and db are cast to Real precision when placed in t_ang1,2, b_ang1,2 and b_len1,2.
/// Solution is carried out in Real precision.
/// Solutions are cast back to Reals when placed into t_ang.
///
/// @references
///
/// @author Evangelos A. Coutsias
/// @author Daniel J. Mandell
///
/// @last_modified June 08, 2008
////////////////////////////////////////////////////////////////////////////////

void bridge_objects (
		const utility::vector1<utility::vector1<Real> >& atoms,
		const utility::vector1<Real> & dt,
		const utility::vector1<Real> & da,
		const utility::vector1<Real> & db,
		const utility::vector1<int>& pivots,
		const utility::vector1<int>& order,
		utility::vector1<utility::vector1<Real> >& t_ang,
		utility::vector1<utility::vector1<Real> >& b_ang,
		utility::vector1<utility::vector1<Real> >& b_len,
		int& nsol) {

	utility::vector1<Real> R0 (3) /* for chainXYZ interface */, R1 (3), R2 (3), R3 (3); // origins of triangle pieces
	utility::vector1<utility::vector1<Real> > Q0; // for chainXYZ interface
	utility::vector1<utility::vector1<Real> > chain1, chain2, chain3, fchain1, fchain2, fchain3; // 3 triangle pieces
	utility::vector1<utility::vector1<Real> > chain1a, chain1b, chain2a, chain2b, chain12, chain12rX; // reconstructions of the molecular chain
	utility::vector1<Real> t_ang1, b_ang1, b_len1, t_ang2, b_ang2, b_len2; // torsions, angles, lengths of chains 1,2
	utility::vector1<Real> vbond (3), xi (3), eta (3), delta (3), theta (3); // triaxial parameters for each pivot
	utility::vector1<Real> cal, sal;
	utility::vector1<utility::vector1<Real> > frame1, frame2, frame3;
	utility::vector1<utility::vector1<Real> > A, B, C, D; // dixon matrices
	utility::vector1<utility::vector1<Real> > cosines, sines, taus;
	utility::vector1<utility::vector1<utility::vector1<Real> > > loop; // the solutions
	int k1, k2, k3, l1, l2, l3, l3a, l3b;
	int ind;
	int N=atoms.size(); // number of atoms in chain
	bool feasible_triangle;

	k1 = pivots[1];
	k2 = pivots[2];
	k3 = pivots[3];
	l1 = k2-k1+1;
	l2 = k3-k2+1;
	l3a = N-k3+1;
	l3b = k1;
	l3  = l3a + l3b;
	chain1.resize(l1);
	chain2.resize(l2);
	chain3.resize(l3);
	t_ang1.resize(l1);
	t_ang2.resize(l2);
	b_ang1.resize(l1);
	b_ang2.resize(l2);
	b_len1.resize(l1);
	b_len2.resize(l2);

	// break the designed torsions, angles, and lengths into 3 chains (3rd 
	// chain's coords set directly)
	ind=1;
	for (int i=pivots[1]; i<=pivots[2]; i++) {
		t_ang1[ind++]=dt[i];
	}
	ind=1;
	for (int i=pivots[2]; i<=pivots[3]; i++) {
		t_ang2[ind++]=dt[i];
	}
	ind=1;
	for (int i=pivots[1]; i<=pivots[2]; i++) {
		b_ang1[ind++]=da[i];
	}
	ind=1;
	for (int i=pivots[2]; i<=pivots[3]; i++) {
		b_ang2[ind++]=da[i];
	}
	ind=1;
	for (int i=pivots[1]; i<=pivots[2]; i++) {
		b_len1[ind++] = db[i];
	}
	ind=1;
	for (int i=pivots[2]; i<=pivots[3]; i++) {
		b_len2[ind++] = db[i];
	}

	chainXYZ( l1, b_len1, b_ang1, t_ang1, false, R0, Q0, chain1 );
	chainXYZ( l2, b_len2, b_ang2, t_ang2, false, R0, Q0, chain2 );

	// set chain 3 directly from original atoms
	ind=1;
	for (int i=k3; i<=N; i++) {
		chain3[ind++] = atoms[i];
	}
	for (int i=1; i<=k1; i++) {
		chain3[ind++] = atoms[i];
	}

	// compute pivot triangle parameters
	chainParams( l1, chain1, vbond[1], xi[1], eta[1], delta[1], R1, frame1 );
	for (int i=1; i<=l1; i++) {
		for (int j=1; j<=3; j++) {
			chain1[i][j] -= R1[j];
		}
	}
	multTransMatrix(frame1, chain1, fchain1);
	chainParams( l2, chain2, vbond[2], xi[2], eta[2], delta[2], R2, frame2 );
	for (int i=1; i<=l2; i++) {
		for (int j=1; j<=3; j++) {
			chain2[i][j] -= R2[j];
		}
	}
	multTransMatrix(frame2, chain2, fchain2);
	chainParams( l3, chain3, vbond[3], xi[3], eta[3], delta[3], R3, frame3 );
	for (int i=1; i<=3; i++) {
		theta[i]  = da[pivots[i]]; // use the designed bond angles
	}
	fchain3=chain3; // for naming consistency

	// check feasibility of arrangement and define polynomial coefficients
	triaxialCoefficients( vbond, xi, eta, delta, theta, order, A, B, C, D, cal, sal, feasible_triangle);

	if (feasible_triangle == false) {
		nsol = 0;
		return;
	}

	// find the solutions
	dixon_eig(A, B, C, D, order, cosines, sines, taus, nsol);

	// reconstruct the molecular chain with new torsions
	loop.resize(nsol);
	t_ang.resize(nsol);
	b_ang.resize(nsol);
	b_len.resize(nsol);
	for (int j=1; j<=nsol; j++) {
		rotateX(fchain1, cosines[j][1], sines[j][1], chain1a);
		rotateZ(chain1a, cal[1], sal[1], chain1b);
		for (unsigned k=1; k<=chain1b.size(); k++) {
			chain1b[k][1] += vbond[3];
		}
		rotateX(fchain2, cosines[j][2], sines[j][2], chain2a);
		rotateZ(chain2a, cal[3], -sal[3], chain2b);
		for (unsigned k=1; k<=chain2b.size(); k++) {
			chain2b[k][1] -= chain2b[l2][1];
			chain2b[k][2] -= chain2b[l2][2];
		}
		// implementing chain12 = [chain1b(:,2:l1), chain2b(:,2:l2-1)];
		int chain12len = (l1-1) + (l2-2);
		chain12.resize(chain12len);
		int ind = 1;
		for (int k=2; k<=l1; k++) {
			chain12[ind].resize(3);
			for (int n=1; n<=3; n++) {
				chain12[ind][n] = chain1b[k][n];
			}
			ind++;
		}
		for (int k=2; k<=l2-1; k++) {
			chain12[ind].resize(3);
			for (int n=1; n<=3; n++) {
				chain12[ind][n] = chain2b[k][n];
			}
			ind++;
		}
		rotateX(chain12, cosines[j][3], -sines[j][3], chain12rX);
		multMatrix(frame3, chain12rX, chain12);
		for (int k=1; k<=chain12len; k++) {
			for (int n=1; n<=3; n++) {
				chain12[k][n] += R3[n];
			}
		}

		// populate the loop table with the cartesian solutions
		loop[j].resize(N);
		for (int k=1; k<=l3b; k++) {
			loop[j][k].resize(3);
			for (int kk=1; kk<=3; kk++) {
				loop[j][k][kk] = chain3[l3a+k][kk];
			}
		}
		for (int k=1; k<= l1 + l2 - 3; k++) {
			loop[j][l3b+k].resize(3);
			for (int kk=1; kk<=3; kk++) {
				loop[j][l3b+k][kk] = chain12[k][kk];
			}
		}
		for (int k=1; k<=l3a; k++) {
			loop[j][l3b+l1+l2-3+k].resize(3);
			for (int kk=1; kk<=3; kk++) {
				loop[j][l3b+l1+l2-3+k][kk] = chain3[k][kk];
			}
		}

		// calculate the 6 closure torsions and place them (with Real precision) 
		// into the output vector together with the designed torsions, angles, and 
		// lengths (which are copied directly from the input)
		int pivots2m1 = pivots[2]-1, pivots2p1 = pivots[2]+1,
			pivots3m1 = pivots[3]-1, pivots3p1 = pivots[3]+1;

		t_ang[j] = dt;
		b_ang[j] = da;
		b_len[j] = db;
		t_ang[j][4] = static_cast<Real> (torsion(loop[j][3],loop[j][4],loop[j][5],loop[j][6])); // pivot 1 always atom 5
		t_ang[j][5] = static_cast<Real> (torsion(loop[j][4],loop[j][5],loop[j][6],loop[j][7]));
		t_ang[j][pivots2m1] = static_cast<Real> (torsion(loop[j][pivots[2]-2],loop[j][pivots2m1],
														loop[j][pivots[2]],loop[j][pivots2p1]));
		t_ang[j][pivots[2]] = static_cast<Real> (torsion(loop[j][pivots2m1],loop[j][pivots[2]],
														  loop[j][pivots2p1],loop[j][pivots[2]+2]));
		t_ang[j][pivots3m1] = static_cast<Real> (torsion(loop[j][pivots[3]-2],loop[j][pivots3m1],
														loop[j][pivots[3]],loop[j][pivots3p1]));
		t_ang[j][pivots[3]] = static_cast<Real> (torsion(loop[j][pivots3m1],loop[j][pivots[3]],
														  loop[j][pivots3p1],loop[j][pivots[3]+2]));
	}
}
// }}}1

}
}
}
