// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   bridgeObjects.cc
/// @brief  Solves the triaxial loop closure problem for a system of atoms
/// @author Evangelos A. Coutsias
/// @author Daniel J. Mandell

// C++ headers
#include <cmath>
// Utility headers
#include <utility/vector1.hh>
// Rosetta Headers
#include <numeric/types.hh>
//#include <basic/Tracer.hh> // tracer output
#include <numeric/kinematic_closure/dixon.hh>
#include <numeric/kinematic_closure/kinematic_closure_helpers.hh>

#include <iostream>

// Constants  **DJM: replace with standard from mini or rosetta++**
#define rad2deg  57.2957795130823208767981548141051703
#define deg2rad 0.0174532925199432957692369076848861271
#define pi 3.14159265358979323846264338327950288
#define pidegs 180.0

using numeric::Real;
//static basic::Tracer std::cout( "numeric.kinematic_closure.bridgeObjects" );

namespace numeric {
namespace kinematic_closure {

void triangle (const utility::vector1<Real>& vbond, utility::vector1<Real>& calpha, utility::vector1<Real>& salpha) {
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

void sincos (const utility::vector1<Real>& theta,
	const int flag,
	utility::vector1<Real>& cosine,
	utility::vector1<Real>& sine)
{
	utility::vector1<Real> a, aa, a1;
	if ( flag == 1 ) {
		for ( int i=1; i<=3; i++ ) {
			cosine[i]=std::cos(theta[i]);
			sine[i]=std::sin(theta[i]);
		}
	} else {
		a.resize(3);
		aa.resize(3);
		a1.resize(3);
		for ( int i=1; i<=3; i++ ) {
			a[i]=std::tan(theta[i]/2);
			aa[i]=a[i]*a[i];
			a1[i]=1+aa[i];
			sine[i]=2*a[i]/a1[i];
			cosine[i]=(1-aa[i])/a1[i];
		}
	}
}

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
	int& f) {

	utility::vector1<Real> ctheta (4), calpha (4), salpha (4), cxi (3), sxi (3), ceta (4), seta (4), cdelta (3), sdelta (3), caceta (4), caseta (4), saceta (4), saseta (4), capeta (4), sapeta (4), cameta (4), sameta (4);
	utility::vector1<utility::vector1<Real> > L, M, N;
	utility::vector1<utility::vector1<utility::vector1<Real> > > p (4);

	// allocate p
	for ( int i=1; i<=4; i++ ) {
		p[i].resize(3);
		for ( int j=1; j<=3; j++ ) {
			p[i][j].resize(3);
		}
	}

	// allocate polynomial matrices
	A.resize(3);
	B.resize(3);
	C.resize(3);
	D.resize(3);
	for ( int i=1; i<=3; i++ ) {
		for ( int j=1; j<=3; j++ ) {
			A[i].resize(3);
			B[i].resize(3);
			C[i].resize(3);
			D[i].resize(3);
		}
	}

	// if triangle is feasible, compute the polynomial matrices
	if ( (std::abs(vb[1] - vb[2]) <= vb[3]) && (vb[1] + vb[2] >= vb[3]) ) {
		f=1; // we have a feasible triangle

		for ( int i=1; i<=3; i++ ) {
			ctheta[i]=std::cos(theta[i]);
		}
		triangle(vb, calpha, salpha);
		cal.resize(3);
		sal.resize(3);
		for ( int i=1; i<=3; i++ ) {
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
		for ( int i=1; i<=4; i++ ) {
			caceta[i] = calpha[i] * ceta[i];
			caseta[i] = calpha[i] * seta[i];
			saceta[i] = salpha[i] * ceta[i];
			saseta[i] = salpha[i] * seta[i];
			capeta[i] = caceta[i] - saseta[i];
			sapeta[i] = saceta[i] + caseta[i];
			cameta[i] = caceta[i] + saseta[i];
			sameta[i] = saceta[i] - caseta[i];
		}
		for ( int i=2; i<=4; i++ ) {
			Real am, ap, cdsx, sdsx, b, bm, bp, d, dp, dm;
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
		for ( int i=1; i<=3; i++ ) {
			for ( int j=1; j<=3; j++ ) {
				p[1][i][j]=p[4][i][j];
			}
		}

		L.resize(3);
		M.resize(3);
		N.resize(3);
		for ( int i=1; i<=3; i++ ) {
			L[i].resize(3);
			M[i].resize(3);
			N[i].resize(3);
		}
		for ( int k=1; k<=3; k++ ) {
			for ( int i=1; i<=3; i++ ) {
				L[k][i] = p[order[1]][i][k];
				M[k][i] = p[order[2]][k][i];
				N[k][i] = p[order[3]][k][i];
			}
		}

		for ( int i=1; i<=3; i++ ) {
			for ( int j=1; j<=3; j++ ) {
				// not using transpose indexing so we don't take transpose in dixon
				A[i][j] = M[i][2]*N[1][j] - M[i][1]*N[2][j];
				B[i][j] = M[i][3]*N[1][j] - M[i][1]*N[3][j];
				C[i][j] = M[i][3]*N[2][j] - M[i][2]*N[3][j];
				D[i][j] = L[i][j];
			}
		}
	} else {
		f=0; // we have an unfeasible triangle
		for ( int i=1; i<=3; i++ ) {
			for ( int j=1; j<=3; j++ ) {
				A[i][j]=0;
				B[i][j]=0;
				C[i][j]=0;
				D[i][j]=0;
			}
		}
	}
}

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

/*
* frame
* -----
* construct body frame for r1 origin
*                          r2 r1-->r2: positive x axis
*                          r3 on pos. y side of xy plane
*/
void frame(const utility::vector1<utility::vector1<Real> >& R, utility::vector1<utility::vector1<Real> >& U) {
	Real norm1, norm3;

	//  utility::vector1<Real> cross1 (3), cross3 (3), cross11 (3), cross12 (3), cross31 (3), cross32 (3);
	utility::vector1<Real> dR3R2 (3);
	U.resize(3);
	for ( int i=1; i<=3; i++ ) {
		U[i].resize(3);
	}
	for ( int i=1; i<=3; i++ ) {
		U[1][i]=R[2][i]-R[1][i];
	}
	norm1=std::sqrt(U[1][1]*U[1][1] + U[1][2]*U[1][2] + U[1][3]*U[1][3]);
	for ( int i=1; i<=3; i++ ) {
		U[1][i] /= norm1;
		dR3R2[i]=R[3][i]-R[2][i];
	}
	cross(U[1], dR3R2, U[3]);
	norm3=std::sqrt(U[3][1]*U[3][1] + U[3][2]*U[3][2] + U[3][3]*U[3][3]);
	for ( int i=1; i<=3; i++ ) {
		U[3][i] /= norm3;
	}
	cross(U[3], U[1], U[2]);
}

Real eucDistance(const utility::vector1<Real>& a, const utility::vector1<Real>& b) {
	return std::sqrt(std::pow(a[1]-b[1],2) + std::pow(a[2]-b[2],2) + std::pow(a[3]-b[3],2));
}

Real scpn(const utility::vector1<Real>& a, const utility::vector1<Real>& b, const utility::vector1<Real>& c) {
	Real d=0;
	for ( int i=1; i<=3; i++ ) {
		d += (a[i]-b[i]) * (c[i]-b[i]);
	}
	return d;
}

Real bondangle(const utility::vector1<Real>& a, const utility::vector1<Real>& b, const utility::vector1<Real>& c) {
	Real r = scpn(a,b,c) / (eucDistance(a,b) * eucDistance(b,c));
	Real ang = std::atan2(sqrt(1-r*r), (Real) r) * rad2deg;
	return ang;
}

Real torsion(const utility::vector1<Real>& a, const utility::vector1<Real>& b, const utility::vector1<Real>& c, const utility::vector1<Real>& d) {
	Real chi;
	utility::vector1<Real> r (3), sc1 (3), sc2 (3), sc3 (3), cs12, cs31;
	utility::vector1<utility::vector1<Real> > s (3);
	for ( int i=1; i<=3; i++ ) { // init s
		s[i].resize(3);
	}
	if ( b==c ) {
		chi=pidegs;
	} else {
		for ( int i=1; i<=3; i++ ) {
			r[i]    = d[i] - b[i];
			s[i][1] = c[i] - b[i];
			s[i][2] = a[i] - b[i];
			sc1[i]  = s[i][1];
			sc2[i]  = s[i][2];
		}
		cross(sc1, sc2, cs12);
		for ( int i=1; i<=3; i++ ) {
			s[i][3] = cs12[i];
			sc3[i]  = s[i][3];
		}
		cross(sc3, sc1, cs31);
		for ( int i=1; i<=3; i++ ) {
			s[i][2] = cs31[i];
		}
		for ( int i=2; i<=3; i++ ) {
			Real f = std::sqrt( pow(s[1][i],2) + pow(s[2][i],2) + pow(s[3][i],2));
			s[1][i] /= f;
			s[2][i] /= f;
			s[3][i] /= f;
		}
		Real y = r[1]*s[1][2] + r[2]*s[2][2] + r[3]*s[3][2];
		Real z = r[1]*s[1][3] + r[2]*s[2][3] + r[3]*s[3][3];
		chi = atan2(z,y) * rad2deg;
		if ( chi < 0 ) {
			chi += 360.0;
		}
	}
	return chi;
}

void chainParams(const int& n, const utility::vector1<utility::vector1<Real> >& atoms, Real& vbond, Real& xi, Real& eta, Real& delta, utility::vector1<Real>& R0, utility::vector1<utility::vector1<Real> >& Q) {
	//utility::vector1<Real> ac1 (3), ac2 (3), acn (3), acneg1 (3);
	utility::vector1<utility::vector1<Real> > a1n2 (3);
	// R0.resize(3);
	// for (int i=1; i<=3; i++) {
	//   ac1[i] = atoms[i][1];
	//   ac2[i] = atoms[i][2];
	//   acn[i] = atoms[i][n];
	//   acneg1[i]= atoms[i][n-1];
	//   R0[i]  = ac1[i];
	//  }
	a1n2[1]=atoms[1];
	a1n2[2]=atoms[n];
	a1n2[3]=atoms[2];

	//  printVector(acn);
	//  printVector(ac1);

	R0=atoms[1];
	frame(a1n2, Q);
	vbond=eucDistance(atoms[n], atoms[1]);
	xi=bondangle(atoms[n-1], atoms[n], atoms[1]);
	eta=bondangle(atoms[n], atoms[1], atoms[2]);
	delta=torsion(atoms[n-1], atoms[n], atoms[1], atoms[2]);

	//  frame(a1n2, Q);
	//  vbond=eucDistance(acn, ac1);
	//  xi=bondangle(acneg1, acn, ac1);
	//  eta=bondangle(acn, ac1, ac2);
	//  delta=torsion(acneg1, acn, ac1, ac2);
}


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

	//std::cout << "about to allocate atoms" << std::endl;
	atoms.resize(n);
	b_len = b_len1; // only defined for consistent naming
	// convert angles to radians and allocate atom vectors
	for ( int i=1; i<=n; i++ ) {
		t_ang[i] = deg2rad*t_ang1[i];
		b_ang[i] = deg2rad*b_ang1[i];
		atoms[i].resize(3);
	}
	//std::cout << "atoms allocated." << std::endl;
	// place the first atom at the origin
	for ( int i=1; i<=n; i++ ) {
		for ( int j=1; j<=3; j++ ) {
			atoms[i][j]=0.0;
		}
	}
	//std::cout << "placed first atom" << std::endl;
	//printMatrix(atoms);
	// place the second atom at bond length along x-axis from first atom
	atoms[2][1] = b_len[1];
	// third atom shows that b_ang is interior
	ca = std::cos(b_ang[2]);
	sa = std::sin(b_ang[2]);
	atoms[3][1] = b_len[1] - b_len[2] * ca;
	atoms[3][2] =            b_len[2] * sa;
	//std::cout << "placed second and third atoms" << std::endl;
	//printMatrix(atoms);

	// initialize rotation matrix
	for ( int i=1; i<=3; i++ ) {
		U[i].resize(3);
	}
	U[1][1] =    -ca;
	U[2][1] =    -sa;
	U[1][2] =     sa;
	U[2][2] =    -ca;
	for ( int i=1; i<=2; i++ ) {
		U[3][i] = 0;
		U[i][3] = 0;
	}
	U[3][3] = 1;

	// all other atoms
	n1 = n-1;
	n2 = n-2;
	j1 = 3;
	j2 = 2;
	//std::cout << "preparing to place all other atoms" << std::endl;
	//std::cout << "U before for loop: " << std::endl;
	//printMatrix(U);
	//exit(0);
	for ( int j=4; j<=n1; j++ ) {
		ca = std::cos(b_ang[j1]);
		sa = std::sin(b_ang[j1]);
		ct = std::cos(t_ang[j2]);
		st = std::sin(t_ang[j2]);
		for ( int i=1; i<=3; i++ ) {
			s[i] = U[2][i];
		}
		for ( int i=1; i<=3; i++ ) {
			U[2][i] = s[i]*ct + U[3][i]*st;
		}
		for ( int i=1; i<=3; i++ ) {
			U[3][i] = -s[i]*st + U[3][i]*ct;
		}
		for ( int i=1; i<=3; i++ ) {
			s[i] = U[2][i];
		}
		for ( int i=1; i<=3; i++ ) {
			U[2][i] = -s[i]*ca - U[1][i]*sa;
		}
		for ( int i=1; i<=3; i++ ) {
			U[1][i] = s[i]*sa - U[1][i]*ca;
		}
		for ( int i=1; i<=3; i++ ) {
			atoms[j][i] = atoms[j1][i] + b_len[j1] * U[1][i];
		}
		j2 = j1;
		j1 = j;
		//std::cout << "at end of j=" << j << std::endl;
		//printMatrix(atoms);
	}
	//std::cout << "ended atom placement for loop. printing U" << std::endl;
	//printMatrix(U);
	//exit(0);
	ca = std::cos(b_ang[n1]);
	sa = std::sin(b_ang[n1]);
	ct = std::cos(t_ang[n2]);
	st = std::sin(t_ang[n2]);
	v[1] =  -b_len[n1] * ca;
	ba   =   b_len[n1] * sa;
	v[2] =   ba * ct;
	v[3] =   ba * st;
	//std::cout << "about to calculate s = U*v'. first printing v" << std::endl;
	//printVector(v);
	// implementing s = U*v'
	for ( int i=1; i<=3; i++ ) {
		s[i] = (U[1][i] * v[1]) + (U[2][i] * v[2]) + (U[3][i] * v[3]);
		//std::cout << "check 1" << std::endl;
	}
	for ( int i=1; i<=3; i++ ) {
		atoms[n][i] = atoms[n1][i] + s[i];
		//std::cout << "check 2" << std::endl;
	}
	//printMatrix(atoms);
	// check atoms here but still need to implement space and frame
	if ( space == true ) {
		// implement:
		// for i = 1:n
		//   atoms(:,i) = R0 + Q*atoms(:,i);
		// end
		utility::vector1<Real> spaceatoms (3); // holder for inner-loop coordinates
		for ( int ind=1; ind<=n; ind++ ) {
			for ( int i=1; i<=3; i++ ) {
				Real inner_sum=0.0;
				for ( int j=1; j<=3; j++ ) {
					inner_sum += (Q[i][j] * atoms[ind][j]);
				}
				spaceatoms[i]=R0[i]+inner_sum;
			}
			for ( int k=1; k<=3; k++ ) {
				atoms[ind][k]=spaceatoms[k];
			}
		}
	}
	/*
	else if (index[1] != 1 || index[2] != 2 || index[3] != 3) {
	std::cout << "setting to body frame based on index" << std::endl;
	utility::vector1<Real> Ra (3);
	utility::vector1<utility::vector1<Real> > RR (3), Qa (3), body_atoms(n);
	//std::cout << "atoms is " << std::endl;
	//printMatrix(atoms);
	for (int i=1; i<=3; i++) {
	RR[i].resize(3);
	Ra[i] = -atoms[index[1]][i];
	}
	//std::cout << "Ra is " << std::endl;
	//printVector(Ra);
	for (int i=1; i<=3; i++) {
	for (int j=1; j<=3; j++) {
	RR[j][i]=atoms[index[j]][i];
	}
	}
	//std::cout << "RR is " << std::endl;
	//printMatrix(RR);
	frame(RR, Qa);
	//std::cout << "Qa is " << std::endl;
	//printMatrix(Qa);
	for (int ind=1; ind<=n; ind++) {
	// implement atoms(:,i) = Qa'*(atoms(:,i)+Ra);
	body_atoms[ind].resize(3);
	for (int i=1; i<=3; i++) {
	Real inner_sum=0.0;
	for (int j=1; j<=3; j++) {
	inner_sum += Qa[i][j] * (atoms[ind][j] + Ra[j]);
	//std::cout << "Qa[" << i << "][" << j << "]: " << Qa[i][j] << std::endl;
	//std::cout << "atoms[" << ind << "][" << j << "]: " << atoms[ind][j] << std::endl;
	}
	body_atoms[ind][i] = inner_sum;
	}
	}
	atoms = body_atoms;
	//exit(0);
	}
	*/
}

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

void chainTORS (const int& n, const utility::vector1<utility::vector1<Real> >& atoms, utility::vector1<Real>& t_ang, utility::vector1<Real>& b_ang, utility::vector1<Real>& b_len, utility::vector1<Real>& R0, utility::vector1<utility::vector1<Real> >& Q) {
	t_ang.resize(n);
	b_ang.resize(n);
	b_len.resize(n);
	utility::vector1<utility::vector1<Real> > framein (3);
	for ( int i=1; i<=3; i++ ) {
		framein[i]=atoms[i];
	}
	R0=atoms[1];
	frame(framein,Q);
	for ( int j=1; j <= n-1; j++ ) {
		b_len[j] = eucDistance( atoms[j], atoms[j+1] );
	}
	b_len[n] = eucDistance( atoms[n], atoms[1] );
	for ( int j=2; j <= n-1; j++ ) {
		b_ang[j] = bondangle( atoms[j-1], atoms[j], atoms[j+1] );
	}
	b_ang[n] = bondangle( atoms[n-1], atoms[n], atoms[1] );
	b_ang[1] = bondangle( atoms[n],   atoms[1], atoms[2] );

	if ( n >= 4 ) {
		for ( int j=2; j <= n-2; j++ ) {
			t_ang[j] = torsion( atoms[j-1], atoms[j], atoms[j+1], atoms[j+2] );
		}
		t_ang[n-1] = torsion( atoms[n-2], atoms[n-1], atoms[n], atoms[1] );
		t_ang[n]   = torsion( atoms[n-1], atoms[n],   atoms[1], atoms[2] );
		t_ang[1]   = torsion( atoms[n]  , atoms[1],   atoms[2], atoms[3] );
	} else {
		for ( int j=1; j<=n; j++ ) {
			t_ang[j]=0.0;
		}
	}
}

void rotateX(const utility::vector1<utility::vector1<Real> >& R, const Real& c, const Real& s, utility::vector1<utility::vector1<Real> >& Rx) {
	int dim1=R.size();
	Rx.resize(dim1);
	for ( int i=1; i<=dim1; i++ ) {
		Rx[i].resize(R[i].size());
		Rx[i][1] =    R[i][1];
		Rx[i][2] =  c*R[i][2]  - s*R[i][3];
		Rx[i][3] =  s*R[i][2]  + c*R[i][3];
	}
}

void rotateY(const utility::vector1<utility::vector1<Real> >& R, const Real& c, const Real& s, utility::vector1<utility::vector1<Real> >& Ry) {
	int dim1=R.size();
	Ry.resize(dim1);
	for ( int i=1; i<=dim1; i++ ) {
		Ry[i].resize(R[i].size());
		Ry[i][1] =   c*R[i][1]  + s*R[i][3];
		Ry[i][2] =     R[i][2];
		Ry[i][3] =  -s*R[i][1]  + c*R[i][3];
	}
}

void rotateZ(const utility::vector1<utility::vector1<Real> >& R, const Real& c, const Real& s, utility::vector1<utility::vector1<Real> >& Rz) {
	int dim1=R.size();
	Rz.resize(dim1);
	for ( int i=1; i<=dim1; i++ ) {
		Rz[i].resize(R[i].size());
		Rz[i][1] =   c*R[i][1]  - s*R[i][2];
		Rz[i][2] =   s*R[i][1]  + c*R[i][2];
		Rz[i][3] =     R[i][3];
	}
}


////////////////////////////////////////////////////////////////////////////////
///
/// @brief
/// Solve the triaxial loop closure problem for a system of atoms
///
/// @details
///
/// @param[in] atoms - matrix of cartesian coordiantes of N-CA-C atoms indexed as atoms[atom][dimension]
/// @param[in] dt - desired torsions for each atom
/// @param[in] da - desired bond angle for each atom
/// @param[in] db - desired bond length for each atom
/// @param[in] pivots - 3 indices (base 1) of atoms to be used as loop closure pivots
/// @param[in] order - length 3 vector giving order to solve for the tau parameters (use 1,2,3 if unsure)
/// @param[out] t_ang - matrix giving torsion angles for each atom for each solution, indexed as t_ang[solution][atom]
/// @param[out] b_ang - matrix giving bond angles for each atom for each solution, indexed as b_ang[solution][atom]
/// @param[out] b_len - matrix giving bond lengths for each atom for each solution, indexed as b_len[solution][atom]
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
////////////////////////////////////////////////////////////////////////////////

void bridgeObjects (const utility::vector1<utility::vector1<Real> >& atoms,
	const utility::vector1<Real> & dt,
	const utility::vector1<Real> & da,
	const utility::vector1<Real> & db,
	const utility::vector1<int>& pivots,
	const utility::vector1<int>& order,
	utility::vector1<utility::vector1<Real> >& t_ang,
	utility::vector1<utility::vector1<Real> >& b_ang,
	utility::vector1<utility::vector1<Real> >& b_len,
	int& nsol)
{
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
	int k1, k2, k3, l1, l2, l3, l3a, l3b, f;
	int ind;
	int N=atoms.size(); // number of atoms in chain

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

	// break the designed torsions, angles, and lengths into 3 chains (3rd chain's coords set directly)
	ind=1;
	for ( int i=pivots[1]; i<=pivots[2]; i++ ) {
		t_ang1[ind++]=dt[i];
	}
	ind=1;
	for ( int i=pivots[2]; i<=pivots[3]; i++ ) {
		t_ang2[ind++]=dt[i];
	}
	ind=1;
	for ( int i=pivots[1]; i<=pivots[2]; i++ ) {
		b_ang1[ind++]=da[i];
	}
	ind=1;
	for ( int i=pivots[2]; i<=pivots[3]; i++ ) {
		b_ang2[ind++]=da[i];
	}
	ind=1;
	for ( int i=pivots[1]; i<=pivots[2]; i++ ) {
		b_len1[ind++] = db[i];
	}
	ind=1;
	for ( int i=pivots[2]; i<=pivots[3]; i++ ) {
		b_len2[ind++] = db[i];
	}

	chainXYZ( l1, b_len1, b_ang1, t_ang1, false, R0, Q0, chain1 );
	chainXYZ( l2, b_len2, b_ang2, t_ang2, false, R0, Q0, chain2 );

	// set chain 3 directly from original atoms
	ind=1;
	for ( int i=k3; i<=N; i++ ) {
		chain3[ind++] = atoms[i];
	}
	for ( int i=1; i<=k1; i++ ) {
		chain3[ind++] = atoms[i];
	}

	/*
	std::cout << "chain 1: " << std::endl;
	printMatrix(chain1);
	std::cout << "chain 2: " << std::endl;
	printMatrix(chain2);
	std::cout << "chain 3: " << std::endl;
	printMatrix(chain3);
	*/

	// compute pivot traingle parameters
	chainParams( l1, chain1, vbond[1], xi[1], eta[1], delta[1], R1, frame1 );
	for ( int i=1; i<=l1; i++ ) {
		for ( int j=1; j<=3; j++ ) {
			chain1[i][j] -= R1[j];
		}
	}
	multTransMatrix(frame1, chain1, fchain1);
	chainParams( l2, chain2, vbond[2], xi[2], eta[2], delta[2], R2, frame2 );
	for ( int i=1; i<=l2; i++ ) {
		for ( int j=1; j<=3; j++ ) {
			chain2[i][j] -= R2[j];
		}
	}
	multTransMatrix(frame2, chain2, fchain2);
	chainParams( l3, chain3, vbond[3], xi[3], eta[3], delta[3], R3, frame3 );
	for ( int i=1; i<=3; i++ ) {
		theta[i]  = deg2rad * da[pivots[i]]; // use the designed bond angles
		eta[i]   *= deg2rad;
		xi[i]    *= deg2rad;
		delta[i] *= deg2rad;
	}
	fchain3=chain3; // for naming consistency

	// check feasibility of arrangement and define polynomial coefficients
	triaxialCoefficients( vbond, xi, eta, delta, theta, order, A, B, C, D, cal, sal, f );

	if ( f == 0 ) {
		nsol = 0;
		return;
	}

	// find the solutions
	dixon_eig( A, B, C, D, order, cosines, sines, taus, nsol );

	// reconstruct the molecular chain with new torsions
	loop.resize(nsol);
	t_ang.resize(nsol);
	b_ang.resize(nsol);
	b_len.resize(nsol);
	for ( int j=1; j<=nsol; j++ ) {
		rotateX(fchain1, cosines[j][1], sines[j][1], chain1a);
		rotateZ(chain1a, cal[1], sal[1], chain1b);
		for ( unsigned k=1; k<=chain1b.size(); k++ ) {
			chain1b[k][1] += vbond[3];
		}
		rotateX(fchain2, cosines[j][2], sines[j][2], chain2a);
		rotateZ(chain2a, cal[3], -sal[3], chain2b);
		for ( unsigned k=1; k<=chain2b.size(); k++ ) {
			chain2b[k][1] -= chain2b[l2][1];
			chain2b[k][2] -= chain2b[l2][2];
		}
		int chain12len = (l1-1) + (l2-2); // implementing chain12 = [chain1b(:,2:l1), chain2b(:,2:l2-1)];
		chain12.resize(chain12len);
		int ind = 1;
		for ( int k=2; k<=l1; k++ ) {
			chain12[ind].resize(3);
			for ( int n=1; n<=3; n++ ) {
				chain12[ind][n] = chain1b[k][n];
			}
			ind++;
		}
		for ( int k=2; k<=l2-1; k++ ) {
			chain12[ind].resize(3);
			for ( int n=1; n<=3; n++ ) {
				chain12[ind][n] = chain2b[k][n];
			}
			ind++;
		}
		rotateX(chain12, cosines[j][3], -sines[j][3], chain12rX);
		multMatrix(frame3, chain12rX, chain12);
		for ( int k=1; k<=chain12len; k++ ) {
			for ( int n=1; n<=3; n++ ) {
				chain12[k][n] += R3[n];
			}
		}

		// populate the loop table with the cartesian solutions
		loop[j].resize(N);
		for ( int k=1; k<=l3b; k++ ) {
			loop[j][k].resize(3);
			for ( int kk=1; kk<=3; kk++ ) {
				loop[j][k][kk] = chain3[l3a+k][kk];
			}
		}
		for ( int k=1; k<= l1 + l2 - 3; k++ ) {
			loop[j][l3b+k].resize(3);
			for ( int kk=1; kk<=3; kk++ ) {
				loop[j][l3b+k][kk] = chain12[k][kk];
			}
		}
		for ( int k=1; k<=l3a; k++ ) {
			loop[j][l3b+l1+l2-3+k].resize(3);
			for ( int kk=1; kk<=3; kk++ ) {
				loop[j][l3b+l1+l2-3+k][kk] = chain3[k][kk];
			}
		}

		// calculate the 6 closure torsions and place them (with Real precision) into the output vector
		// together with the designed torsions, angles, and lengths (which are copied directly from the input)

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


void test_bridgeObjects() {
	// input arguments
	int N=8;
	utility::vector1<Real> t_ang0 (N), b_ang0 (N), b_len0 (N), dt (N), db (N), da (N);
	utility::vector1<int> pivots (3), order (3);
	// input argument values
	// DJM [ if inputting torsions, uncomment next 4 lines ]
	Real t_ang0vals[] = {36.0213, 0.0, 253.8843, 52.8441, 36.0366, 359.8985, 253.9602, 52.8272};
	Real b_ang0vals[] = {114.9908, 114.9805, 115.0112, 114.9923, 114.9582, 115.0341, 114.9920, 115.0594};
	Real b_len0vals[] = {1.5200, 1.5202, 1.5197, 1.5205, 1.5198, 1.5201, 1.5197, 1.5196};
	Real dtvals[] = {180.0000, 180.0000, -5.3651, 180.0000, 180.0000, 180.0000, 180.0000, -5.0393};

	utility::vector1<utility::vector1<Real> > atoms (N);
	// output arguments
	int nsol;
	utility::vector1<utility::vector1<Real> > t_ang, b_ang, b_len;
	/*
	* DJM uncomment if inputting torsions
	*/
	int ind=0;
	for ( int i=1; i<=8; i++ ) {
		t_ang0[i]=t_ang0vals[ind++];
	}

	ind=0;
	for ( int i=1; i<=8; i++ ) {
		b_ang0[i]=b_ang0vals[ind++];
	}

	ind=0;
	for ( int i=1; i<=8; i++ ) {
		b_len0[i]=b_len0vals[ind++];
	}

	ind=0;
	for ( int i=1; i<=8; i++ ) {
		dt[i]=dtvals[ind++];
	}

	for ( int i=1; i<=8; i++ ) {
		db[i]=1.52;
	}

	for ( int i=1; i<=8; i++ ) {
		da[i]=115;
	}
	//*/
	pivots[1] = 2;
	pivots[2] = 5;
	pivots[3] = 7;
	order[1] = 1;
	order[2] = 2;
	order[3] = 3;

	//bridgeObjects(N, t_ang0, b_ang0, b_len0, dt, db, da, pivots, order, t_ang, b_ang, b_len, nsol);
	bridgeObjects(atoms, dt, db, da, pivots, order, t_ang, b_ang, b_len, nsol);
	//bridgeObjects(atoms, pivots, order, t_ang, b_ang, b_len, nsol);
	//printMatrix(t_ang);
	//printMatrix(b_ang);
	//printMatrix(b_len);
	////std::cout << nsol << std::endl;
}

/* // DJM: this was for benchmarking current version of bridgeObjects against earlier (slower) version.
// v2 is now current version. 6-8-08
void test_bridgeObjects_v2() {
// input arguments
int N=27;
utility::vector1<Real> dt (N), db (N), da (N);
utility::vector1<int> pivots (3), order (3);
// input argument values
// DJM [ if inputting torsions, uncomment next 4 lines ]
Real dt_ang_vals[] = {348.798, 337.654, 181.609, -170, 140, 179.743, -140, 80, 181.769, -130, 10, 180.419, -90, -30, 177.171, -90, -20, 184.154, -90, -40, 182.472, -80, 170, 184.776, 310.274, 150.675, 290.245};
Real db_ang_vals[] = {88.2267, 112.669, 115.978, 122.109, 110.442, 117.607, 121.606, 110.194, 115.415, 120.286, 110.888, 117.314, 125.907, 110.287, 117.773, 121.692, 109.691, 115.437, 122.169, 109.827, 119.212, 120.735, 111.479, 115.289, 121.921, 114.467, 85.1688};
Real db_len_vals[] = {1.41632, 1.48472, 1.32645, 1.46204, 1.54295, 1.31506, 1.45302, 1.51076, 1.28081, 1.46153, 1.51705, 1.33131, 1.46801, 1.5172, 1.34797, 1.45556, 1.49902, 1.31172, 1.41681, 1.49235, 1.36072, 1.4561, 1.52821, 1.31126, 1.44372, 1.48442, 9.90623};
// initialize the atoms vector
utility::vector1<utility::vector1<Real> > atoms (N);
Real atoms_vals[] = {13.92, 13.881, 12.607, 11.925, 10.671, 9.481, 9.657, 8.597, 8.771, 8.01, 8.116, 6.871, 5.803, 4.504, 3.703, 2.616, 1.752, 0.851, 0.206, -0.747, -0.251, 0.965, 1.458, 2.652, 2.912, 4.066, 4.232, 43.318, 44.382, 44.429, 43.293, 43.179, 43.781, 44.046, 44.619, 46.119, 46.785, 48.241, 48.954, 48.205, 48.696, 47.567, 47.894, 46.878, 46.229, 47.052, 46.598, 45.73, 45.989, 45.197, 44.279, 44.041, 43.277, 41.984, 25.637, 24.703, 23.942, 23.88, 23.137, 23.913, 25.189, 26.001, 26.047, 25.261, 25.191, 25.684, 25.95, 26.426, 27.047, 27.774, 28.357, 27.35, 26.558, 25.613, 24.505, 23.952, 22.834, 23.093, 24.356, 24.767, 24.057};
for (int i=1; i<=N; i++) {
atoms[i].resize(3);
}
int ind=0;
for (int j=1; j<=3; j++) {
for (int i=1; i<=N; i++) {
atoms[i][j] = atoms_vals[ind++];
}
}
ind=0;
// initialize the designed torsions, angles, and length
for (int i=1; i<=N; i++) {
dt[i]=dt_ang_vals[ind++];
}
ind=0;
for (int i=1; i<=N; i++) {
da[i]=db_ang_vals[ind++];
}
ind=0;
for (int i=1; i<=N; i++) {
db[i]=db_len_vals[ind++];
}
// these pivots specifically correspond to the designed inputs to produce the expected outputs
pivots[1] = 5;
pivots[2] = 14;
pivots[3] = 23;
order[1] = 1;
order[2] = 2;
order[3] = 3;

// DJM: debug
//std::cout << "dt: " << std::endl;
//printVector(dt);
//std::cout << "db: " << std::endl;
//printVector(db);
//std::cout << "da: " << std::endl;
//printVector(da);
// output arguments
int nits=100000; // number times to call bridgeObjects for performance benchmarking
int nsol;
utility::vector1<utility::vector1<Real> > t_ang, b_ang, b_len;

std::cout << "calling bridge2" << std::endl;
std::time_t const start_time_v2( std::time( NULL ) );
for (int d=1; d<=nits; d++) {
bridgeObjects_v2(atoms, dt, da, db, pivots, order, t_ang, b_ang, b_len, nsol);
}
std::time_t const end_time_v2( std::time( NULL ) );
std::cout << "start time: " << start_time_v2 << std::endl;
std::cout << "end time: " << end_time_v2 << std::endl;
std::cout << "start time minus end time v2: " << end_time_v2 - start_time_v2 << std::endl;
std::cout << "from bridgeObjects v2: " << std::endl;
printMatrix(t_ang);

nsol=0;

std::cout << "calling bridge1" << std::endl;
std::time_t const start_time_v1( std::time( NULL ) );
for (int d=1; d<=nits; d++) {
bridgeObjects(atoms, dt, da, db, pivots, order, t_ang, b_ang, b_len, nsol);
}
std::time_t const end_time_v1( std::time( NULL ) );
std::cout << "start time: " << start_time_v1 << std::endl;
std::cout << "end time: " << end_time_v1 << std::endl;
std::cout << "start time minus end time v1: " << end_time_v1 - start_time_v1 << std::endl;
std::cout << "from bridgeObjects v1: " << std::endl;
printMatrix(t_ang);

//std::cout << "output t_ang: " << std::endl;
//printMatrix(t_ang);
//std::cout << "output b_ang: " << std::endl;
//printMatrix(b_ang);
//std::cout << "output b_len: " << std::endl;
//printMatrix(b_len);
////std::cout << nsol << std::endl;
}
*/

void test_rotateX() {
	Real cos1 = -0.852235910816473;
	Real sin1 = -0.523157674429820;
	Real chain1vals[] = { 0.0, 0.0, 0.0, 0.644662961550365, 1.376520855637543, 0.0, 2.163059129530747, 1.374750785448438, -0.069784983442129, 2.807722091081112, 0.0, 0.0 };
	utility::vector1<utility::vector1<Real> > chain1 (4), chain1a;
	// populate chain1
	int ind=0;
	for ( int i=1; i<=4; i++ ) {
		chain1[i].resize(3);
		for ( int j=1; j<=3; j++ ) {
			chain1[i][j]=chain1vals[ind++];
		}
	}
	rotateX(chain1, cos1, sin1, chain1a);
	printMatrix(chain1a);
}

void test_chainTORS () {
	// inputs
	int n=8;
	utility::vector1<utility::vector1<Real> > atoms (n);
	Real atomVals[] = {2.1636, 1.3766, 0, 1.6816, 2.3109, 1.0978, 1.1202, 3.6372, 0.6119, -0.3753, 3.6409, 0.3402, -1.0556, 2.2881, 0.4729, -0.6424, 1.2567, -0.5643, 0.0, 0.0, 0.0, 1.5200, 0.0,  0.0}; // for populating the atoms input

	//outputs
	utility::vector1<Real> t_ang, b_ang, b_len, R0;
	utility::vector1<utility::vector1<Real> > Q;

	// allocate and populate the atoms
	int ind=0;
	for ( int i=1; i<=n; i++ ) {
		atoms[i].resize(3);
		for ( int j=1; j<=3; j++ ) {
			atoms[i][j]=atomVals[ind++];
		}
	}

	chainTORS(n, atoms, t_ang, b_ang, b_len, R0, Q);
	printVector(t_ang);
	printVector(b_ang);
	printVector(b_len);
	printVector(R0);
	printMatrix(Q);

	// output should be:
	//
	//   [t_ang,   b_ang,      b_len]
	//   236.3251  114.9908    1.5200
	//    90.4073  115.0000    1.5200
	//   354.6349  115.0000    1.5200
	//   293.3675  115.0000    1.5200
	//   119.1794  115.0000    1.5200
	//   263.4791  115.0000    1.5200
	//    24.1804  115.0000    1.5200
	//    52.8272  115.0594    1.5196
	//
	//    R0: [2.1636, 1.3766, 0]'
	//    Q:
	//    -0.3172   -0.2596   -0.9121
	//     0.6147    0.6762   -0.4062
	//     0.7222   -0.6895   -0.0549
}

void test_chainXYZ() {
	utility::vector1<utility::vector1<Real> > atoms, Q(3);
	int n = 4;
	utility::vector1<Real> b_len (n), b_ang(n), t_ang(n), R0 (3);

	// setup bond lengths, bond angles, and torsions
	for ( int i=1; i<=n; i++ ) {
		b_len[i]=1.52;
		b_ang[i]=115.0;
	}
	t_ang[1] = 334.5902397604552;
	t_ang[2] = 52.8271503523038;
	t_ang[3] = 334.6192090171498;
	t_ang[4] = 180.0;

	// Q is the identity matrix and R0 is the origin
	for ( int i=1; i<=3; i++ ) {
		R0[i]=0.0;
		Q[i].resize(3);
		for ( int j=1; j<=3; j++ ) {
			(i==j) ? Q[i][j]=1.0 : Q[i][j]=0.0;
		}
	}

	chainXYZ(n, b_len, b_ang, t_ang, false, R0, Q, atoms);
	printMatrix(atoms);

	// output should be:
	//                    0   1.520000000000000   2.163648967954409   1.681566958628807
	//                    0                   0   1.376581274771391   2.310875247341517
	//                    0                   0                   0   1.097766691562338
}

/* // DJM: this was for benchmarking current version of chainXYZ against earlier (slower) version.
// v2 is now current version. 6-8-08
void test_chainXYZv2() {
utility::vector1<utility::vector1<Real> > atoms_v1, atoms_v2, Q(3);
int n = 12;
utility::vector1<Real> b_len (n), b_ang(n), t_ang(n), R0 (3);
utility::vector1<int> index(3);
int nits = 1000000; // number of test iterations for time measurement

// setup bond lengths, bond angles, and torsions
for (int i=1; i<=n; i++) {
b_len[i]=1.52;
b_ang[i]=115.0;
}
t_ang[1] = 334.5902397604552;
t_ang[2] = 52.8271503523038;
t_ang[3] = 334.6192090171498;
t_ang[4] = 180.0;
t_ang[5] = 334.5902397604552;
t_ang[6] = 52.8271503523038;
t_ang[7] = 334.6192090171498;
t_ang[8] = 180.0;
t_ang[9] = 334.5902397604552;
t_ang[10] = 52.8271503523038;
t_ang[11] = 334.6192090171498;
t_ang[12] = 180.0;

// Q is the identity matrix and R0 is the origin
for (int i=1; i<=3; i++) {
R0[i]=0.0;
Q[i].resize(3);
for (int j=1; j<=3; j++) {
(i==j) ? Q[i][j]=1.0 : Q[i][j]=0.0;
}
}

// set up the index vector to select which atoms will establish the frame
index[1] = 1;
//index[2] = n;
//index[3] = 2;
index[2] = 2;
index[3] = 3;

std::time_t const start_time_v1( std::time( NULL ) );
for (int d=1; d<=nits; d++) {
chainXYZ(n, b_len, b_ang, t_ang, false, R0, Q, atoms_v1);
}
std::time_t const end_time_v1( std::time( NULL ) );
std::cout << "start time: " << start_time_v1 << std::endl;
std::cout << "end time: " << end_time_v1 << std::endl;
std::cout << "start time minus end time v1: " << end_time_v1 - start_time_v1 << std::endl;

//std::cout << "chainXYZ v1 says " << std::endl;
//printMatrix(atoms_v1);
std::time_t const start_time_v2( std::time( NULL ) );
for (int d=1; d<=nits; d++) {
chainXYZv2(n, b_len, b_ang, t_ang, false, R0, Q, atoms_v2);
}
std::time_t const end_time_v2( std::time( NULL ) );
std::cout << "start time: " << start_time_v2 << std::endl;
std::cout << "end time: " << end_time_v2 << std::endl;
std::cout << "start time minus end time v2: " << end_time_v2 - start_time_v2 << std::endl;
//std::cout << "chainXYZ v2 says " << std::endl;
//printMatrix(atoms_v2);

// output should be:
//                    0   1.520000000000000   2.163648967954409   1.681566958628807
//                    0                   0   1.376581274771391   2.310875247341517
//                    0                   0                   0   1.097766691562338
}
*/

void test_chainParams() {
	int n=4;
	Real vbond, xi, eta, delta;
	Real atomVals[] = { -0.160000000000000,  1.225000000000001, -2.247000000000002, \
		-0.642082009325601,  1.113981114034975, -0.809755503317661, \
		0.000939232598462, -0.000706370692062, -0.000808351874906, \
		-0.481142776727140, -0.111725256657088,  1.436436144807435 };
	utility::vector1<Real> R0;
	utility::vector1<utility::vector1<Real> > atoms (4), Q;
	int ind=0;
	for ( int i=1; i<=4; i++ ) {
		atoms[i].resize(3);
		for ( int j=1; j<=3; j++ ) {
			atoms[i][j]=atomVals[ind++];
		}
	}
	//printMatrix(atoms);
	chainParams(n, atoms, vbond, xi, eta, delta, R0, Q);
	//std::cout << "vbond: " << vbond << std::endl;
	//std::cout << "xi: " << xi << std::endl;
	//std::cout << "eta: " << eta << std::endl;
	//std::cout << "delta: " << delta << std::endl;
	//std::cout << "R0: " << std::endl;
	printVector(R0);
	//std::cout << "Q: " << std::endl;
	printMatrix(Q);
	// output should be:
	// vbond: 3.93162
	// xi: 20.511
	// eta: 20.511
	// delta: 180

}

void test_torsion() {
	utility::vector1<Real> a (3), b(3), c(3), d(3);
	a[1]=  0.000939232598461515;
	a[2]= -0.000706370692062031;
	a[3]= -0.000808351874906243;
	b[1]= -0.481142776727140;
	b[2]= -0.111725256657088;
	b[3]=  1.436436144807435;
	c[1]= -0.160000000000000;
	c[2]=  1.225000000000001;
	c[3]= -2.247000000000002;
	d[1]= -0.642082009325601;
	d[2]=  1.113981114034975;
	d[3]= -0.809755503317661;
	Real tors=torsion(a,b,c,d);
	std::cout << tors << std::endl;
	// output should be 180.0
}

void test_bondangle() {
	utility::vector1<Real> a (3), b (3), c (3);
	a[1]=  0.000939232598461515;
	a[2]= -0.000706370692062031;
	a[3]= -0.000808351874906243;
	b[1]= -0.481142776727140;
	b[2]= -0.111725256657088;
	b[3]=  1.436436144807435;
	c[1]= -0.160000000000000;
	c[2]=  1.225000000000001;
	c[3]= -2.247000000000002;
	Real angle=bondangle(a,b,c);
	std::cout << angle << std::endl;
	// output should be 20.510953731276413
}

void test_scpn() {
	utility::vector1<Real> a (3), b (3), c (3);
	a[1]=  0.000939232598461515;
	a[2]= -0.000706370692062031;
	a[3]= -0.000808351874906243;
	b[1]= -0.481142776727140;
	b[2]= -0.111725256657088;
	b[3]=  1.436436144807435;
	c[1]= -0.160000000000000;
	c[2]=  1.225000000000001;
	c[3]= -2.247000000000002;
	Real scpn_val=scpn(a,b,c);
	std::cout << scpn_val << std::endl;
	// output should be 5.597217231925712
}

void test_eucDistance() {
	utility::vector1<Real> a (3), b (3);
	a[1]= -0.481142776727140;
	a[2]= -0.111725256657088;
	a[3]=  1.436436144807435;
	b[1]= -0.160000000000000;
	b[2]=  1.225000000000001;
	b[3]= -2.247000000000002;
	Real dist=eucDistance(a,b);
	std::cout << dist << std::endl;
	// output should be 3.931624209878514
}

void test_frame() {
	utility::vector1<utility::vector1<Real> > R (3), U;
	Real Rvals[] = { -0.160000000000000,  1.225000000000001, -2.247000000000002, \
		-0.481142776727140, -0.111725256657088,  1.436436144807435, \
		-0.642082009325601,  1.113981114034975, -0.809755503317661 };

	// intialize R
	int ind=0;
	for ( int i=1; i<=3; i++ ) {
		R[i].resize(3);
		for ( int j=1; j<=3; j++ ) {
			R[i][j]=Rvals[ind++];
		}
	}
	frame(R, U);
	printMatrix(U);
}

void test_cross() {
	utility::vector1<Real> L (3), r0 (3), r;
	L[1]=  -0.081681961343163;
	L[2]=  -0.339993139043773;
	L[3]=   0.936873909656095;
	r0[1]= -0.160939232598461;
	r0[2]=  1.225706370692063;
	r0[3]= -2.246191648125096;
	cross(L, r0, r);
	printVector(r);
}

void test_triaxialCoefficients() {
	// inputs
	utility::vector1<Real> vb (3), xi (3), eta (3), delta (3), theta (3);
	utility::vector1<int> order (3);

	// outputs
	utility::vector1<Real> cal, sal;
	utility::vector1<utility::vector1<Real> > A, B, C, D;
	int f;

	// dixon outputs (if testing dixon as well)
	utility::vector1<utility::vector1<Real> > sines, cosines, tau;
	int nsol;

	// initialize inputs
	vb[1]   = 2.807722146231753;
	vb[2]   = 2.563909995271172;
	vb[3]   = 2.807373715195975;
	xi[1]   = 1.132805948744302;
	xi[2]   = 0.567232006898157;
	xi[3]   = 1.133000954486070;
	eta[1]  = 1.132805948744302;
	eta[2]  = 0.567232006898157;
	eta[3]  = 1.133000954486070;
	delta[1]= 6.232466452684100;
	delta[2]= 0.0;
	delta[3]= 6.235542719576548;
	theta[1]= 2.007128639793479;
	theta[2]= 2.007128639793479;
	theta[3]= 2.007128639793479;
	order[1]= 1;
	order[2]= 2;
	order[3]= 3;

	triaxialCoefficients(vb, xi, eta, delta, theta, order, A, B, C, D, cal, sal, f);
	printMatrix(B);
	dixon_sturm(A,B,C,D,order,sines,cosines,tau,nsol);
	//std::cout << nsol << std::endl;
	//printMatrix(sines);

	// cal and sal output should be
	//cal[1]= -0.583014272556966;
	//cal[2]= -0.456717749075032;
	//cal[3]= -0.456502620110421;
	//sal[1]=  0.812461911719480;
	//sal[2]=  0.889611655544056;
	//sal[3]=  0.889722067744934;
}

void test_sincos() {
	utility::vector1<Real> theta (3), cosine, sine;
	int flag=2;
	theta[1]=1.132805948744302;
	theta[2]=0.567232006898157;
	theta[3]=1.133000954486070;

	sincos(theta, flag, cosine, sine);
}

void test_triangle() {
	utility::vector1<Real> vbond (3), calpha, salpha;
	// initialize vbond
	vbond[1]=2.807722146231753;
	vbond[2]=2.563909995271172;
	vbond[3]=2.807373715195975;

	triangle(vbond, calpha, salpha);
	printVector(calpha);
	//std::cout << std::endl;
	printVector(salpha);
}

} // end namespace kinematic_closure
} // end namespace numeric

/*
int main(int argc, char** argv)
{
//test_triangle();
//test_sincos();
//test_triaxialCoefficients();
//test_cross();
//test_frame();
//test_eucDistance();
//test_scpn();
//test_bondangle();
//test_torsion();
//test_chainParams();
//test_chainXYZ();
//numeric::kinematic_closure::test_chainXYZv2();
//test_chainTORS();
//test_rotateX();
//test_bridgeObjects();
numeric::kinematic_closure::test_bridgeObjects_v2();
}
*/
