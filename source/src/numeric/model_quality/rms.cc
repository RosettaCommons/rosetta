// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/model_quality/rms.hh
/// @brief  RMS functions imported from rosetta++
/// @author James Thompson
/// @date   Wed Aug 22 12:10:37 2007


// Rosetta Headers

// numeric libraries
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/model_quality/rms.hh>
#include <numeric/model_quality/maxsub.hh>
#include <numeric/types.hh>

#include <utility/vector1.hh>

// ObjexxFCL libraries
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>

// C++ libraries
#include <complex>


namespace numeric {
namespace model_quality {

using ObjexxFCL::FArray2D;
using ObjexxFCL::FArray1D;
using ObjexxFCL::FArray2A;
using ObjexxFCL::FArray1A;
using ObjexxFCL::FArray2;
using ObjexxFCL::FArray1;

numeric::Real
calc_rms(
	utility::vector1< xyzVector< Real > > p1_coords,
	utility::vector1< xyzVector< Real > > p2_coords
) {
	assert( p1_coords.size() == p2_coords.size() );
	Size const natoms( p1_coords.size() );

	FArray2D< numeric::Real > p1a( 3, p1_coords.size() );
	FArray2D< numeric::Real > p2a( 3, p2_coords.size() );

	for ( Size i = 1; i <= natoms; ++i ) {
		for ( Size k = 1; k <= 3; ++k ) { // k = X, Y and Z
			p1a(k,i) = p1_coords[i](k);
			p2a(k,i) = p2_coords[i](k);
		}
	}
	return rms_wrapper( natoms, p1a, p2a );
}

numeric::Real
rms_wrapper_slow_and_correct(
	int natoms,
	FArray2D< numeric::Real > p1a,
	FArray2D< numeric::Real > p2a
) {
	// rotate and translate coordinates to minimize rmsd
	FArray1D< numeric::Real > ww(natoms,1.0);
	Real bogus = 0;
	rmsfitca2(natoms,p1a,p2a,ww,natoms,bogus);

	// manually calculate rmsd
	numeric::Real tot = 0;
	for ( int i = 1; i <= natoms; ++i ) {
		for ( int j = 1; j <= 3; ++j ) {
			tot += std::pow( p1a(i,j) - p2a(i,j), 2 );
		}
	}

	return std::sqrt(tot/natoms);
}

// Calculate an RMS based on aligned set of points in p1a and p2a composed
// each representing a list of natoms.
numeric::Real
rms_wrapper(
	int natoms,
	FArray2D< numeric::Real > p1a,
	FArray2D< numeric::Real > p2a
) {
	FArray1D< numeric::Real > ww( natoms, 1.0 );
	FArray2D< numeric::Real > uu( 3, 3, 0.0 );
	numeric::Real ctx;

	findUU( p1a, p2a, ww, natoms, uu, ctx );

	float fast_rms;
	calc_rms_fast( fast_rms, p1a, p2a, ww, natoms, ctx );
	return fast_rms;
} // rms_wrapper

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
/// companion function for findUU ( it is optional ) computes the minimum
/// RMS deviation beteen XX and YY as though it rotated the arrays, without
/// actually rotating them.
///
/// @details
///
/// @param  rms_out   [in/out]? - the real-valued output value of the rms deviation
/// @param  XX - [in/out]? - first set of points representing xyz coordinates
/// @param  YY - [in/out]? - second set of points representing xyz coordinates
/// @param  WW - [in/out]? - relative weight for each point
/// @param  npoints - [in/out]? - number of points
/// @param  ctx - [in/out]? - magic number computed during find UU that is needed
///               for this calculation
///
/// @global_read
///
/// @global_write
///
/// @remarks
/// the XX, YY, WW must be the same as call to findUU (remember that findUU
/// offsets the XX and YY weighted COM to the origin!)
///
/// @references
///
/// @author
///
/////////////////////////////////////////////////////////////////////////////////
void
calc_rms_fast(
	float & rms_out,
	FArray2A< numeric::Real > xx,
	FArray2A< numeric::Real > yy,
	FArray1A< numeric::Real > ww,
	int npoints,
	numeric::Real ctx
)
{
	xx.dimension( 3, npoints );
	yy.dimension( 3, npoints );
	ww.dimension( npoints );

	numeric::Real rms = 0.0;

	for ( int k = 1; k <= npoints; ++k ) {
		numeric::Real const ww_k = ww(k);
		for ( int j = 1; j <= 3; ++j ) {
			numeric::Real const xx_jk = xx(j,k);
			numeric::Real const yy_jk = yy(j,k);
			rms += ww_k * ( ( xx_jk * xx_jk ) + ( yy_jk * yy_jk ) );
		}
	}

	rms -= 2 * ctx;
	// abs on next line catches small difference of large numbers that turn out
	// to be accidental small but negative
	rms_out = std::sqrt(std::abs(rms/npoints)); // return a float
} // calc_rms_fast

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
/// intended to rotate one protein xyz array onto another one such that
/// the point-by-point rms is minimized.
///
/// @details
///   1) ORIGINAL PAPER HAD ERROR IN HANDEDNESS OF VECTORS, LEADING
///      TO INVERSION MATRICIES ON OCCASION. OOPS. NOW FIXED.
///       SEE ACTA CRYST(1978) A34 PAGE 827 FOR REVISED MATH
///   2) TRAP DIVIDE BY ZERO ERRORS WHEN NO ROTATIONS REQUIRED.
///   3) ADDED WEIGHTS (WEIGHTS NOW WORK)
///   4) ADDED FAST RMS CALC AUXILIRARY ROUTINE.
///   5) CHANGED TO numeric::Real TO DEAL WITH HIGHLY DISSIMILAR BUT LARGE PROTEINS.
///
/// switched order of array subscripts so that can use logical array sizes
/// XX and YY are lists of Npoints  XYZ vectors (3xNpoint matrix) to be co-aligned
/// these matrices are returned slightly modified: they are translated so their origins
/// are at the center of mass (see Weights below ).
/// WW is the weight or importance of each point (weighted RMS) vector of size Npoints
/// The center of mass is figured including this variable.
/// UU is a 3x3 symmetric orthonornal rotation matrix that will rotate YY onto XX such that
/// the weighted RMS distance is minimized.
///
/// @param[in]   XX - in - XX,YY  are 2D arrays of x,y,z position of each atom
/// @param[in]   YY - in -
/// @param[in]   WW - in -  a weight matrix for the points
/// @param[in]   Npoints - in - the number of XYZ points (need not be physical array size)
/// @param[out]   UU - out - 3x3 rotation matrix.
/// @param[out]   sigma3 - out - TO BE PASSED TO OPTIONAL FAST_RMS CALC ROUTINE.
///
/// @global_read
///
/// @global_write
///
/// @remarks
/// SIDEEFECTS: the Matrices XX, YY are Modified so that their weighted center
///         of mass is  moved to (0,0,0).
///
/// CAVEATS:
///      1) it is CRITICAL that the first physical dimension of XX and YY is 3
///
///      2) an iterative approx algorithm computes the diagonalization of
///         a 3x3 matrix.  if this program needs a speed up this could be
///         made into an analytic but uggggly diagonalization function.
///
/// @references
/// Mathethematical Basis from paper:
/// (Wolfgang Kabsch) acta Cryst (1976) A32 page 922
///
/// @author
///  Charlie Strauss 1999
///  Revised april 22
///
/////////////////////////////////////////////////////////////////////////////////
void
findUU(
	FArray2< numeric::Real > & XX,
	FArray2< numeric::Real > & YY,
	FArray1< numeric::Real > const & WW,
	int Npoints,
	FArray2< numeric::Real > & UU,
	numeric::Real & sigma3
)
{
	using numeric::xyzMatrix;
	using numeric::xyzVector;

	ObjexxFCL::FArray1D_int sort( 3 );
	FArray2D< numeric::Real > eVec( 3, 3 );
	FArray2D< numeric::Real > bb( 3, 3 );
	FArray1D< numeric::Real > w_w( 3 );
	FArray2D< numeric::Real > m_moment( 3, 3 );
	FArray2D< numeric::Real > rr_moment( 3, 3 );
	numeric::Real temp1;
	numeric::Real temp2;
	//numeric::Real temp3;
	FArray1D< numeric::Real > Ra( 3 );

	if ( Npoints < 1 ) {

		// return identity rotation matrix to moron
		for ( int i = 1; i <= 3; ++i ) {
			for ( int k = 1; k <= 3; ++k ) {
				UU(i,k) = 0.0;
				if ( i == k ) UU(i,i) = 1.0;
			}
		}
		sigma3 = 0.0;
		return;
	}

	// align center of mass to origin
	for ( int k = 1; k <= 3; ++k ) {
		temp1 = 0.0;
		temp2 = 0.0;
		numeric::Real temp3 = 0.0;
		for ( int j = 1; j <= Npoints; ++j ) {
			temp1 += XX(k,j) * WW(j);
			temp2 += YY(k,j) * WW(j);
			temp3 += WW(j);
		}
		if ( temp3 > 0.001 ) temp1 /= temp3;
		if ( temp3 > 0.001 ) temp2 /= temp3;

		for ( int j = 1; j <= Npoints; ++j ) {
			XX(k,j) -= temp1;
			YY(k,j) -= temp2;
		}
	}

	// Make cross moments matrix   INCLUDE THE WEIGHTS HERE
	for ( int k = 1; k <= 3; ++k ) {
		for ( int j = 1; j <= 3; ++j ) {
			temp1 = 0.0;
			for ( int i = 1; i <= Npoints; ++i ) {
				temp1 += WW(i) * YY(k,i) * XX(j,i);
			}
			m_moment(k,j) = temp1;
		}
	}

	// Multiply CROSS MOMENTS by transpose
	BlankMatrixMult(m_moment,3,3,1,m_moment,3,0,rr_moment);

	// Copy to/from xyzMatrix/xyzVector since rest of functions use FArrays
	xyzMatrix< numeric::Real > xyz_rr_moment( xyzMatrix< numeric::Real >::cols( &rr_moment( 1,1 ) ) );
	xyzVector< numeric::Real > xyz_w_w;
	xyzMatrix< numeric::Real > xyz_eVec;

	// Find eigenvalues, eigenvectors of symmetric matrix rr_moment
	xyz_w_w = eigenvector_jacobi( xyz_rr_moment, (numeric::Real) 1E-9, xyz_eVec );

	// Copy eigenvalues/vectors back to FArray
	for ( int i = 1; i <= 3; ++i ) {
		w_w( i ) = xyz_w_w( i );
		for ( int j = 1; j <= 3; ++j ) {
			eVec( i, j ) = xyz_eVec( i, j );
		}
	}

	// explicitly coded 3 level index sort using eigenvalues
	for ( int i = 1; i <= 3; ++i ) {
		sort(i) = i;
	}

	if ( w_w(1) < w_w(2) ) {
		sort(2) = 1;
		sort(1) = 2;
	}

	if ( w_w(sort(2)) < w_w(3) ) {
		sort(3) = sort(2);
		sort(2) = 3;

		if ( w_w(sort(1)) < w_w(3) ) {
			sort(2) = sort(1);
			sort(1) = 3;
		}
	}

	// sort is now an index to order of eigen values

	if ( w_w(sort(2)) == 0.0 ) { // holy smokes, two eigen values are zeros
		// return identity rotation matrix to moron
		for ( int i = 1; i <= 3; ++i ) {
			for ( int k = 1; k <= 3; ++k ) {
				UU(i,k) = 0.0;
			}
			UU(i,i) = 1.0;
		}
		if ( w_w(sort(1)) < 0.0 ) {
			w_w(sort(1)) = std::abs(w_w(sort(1)));
		}
		sigma3 = std::sqrt(w_w(sort(1)));

		return; // make like a prom dress and slip off
	}

	// sort eigen values
	temp1 = w_w(sort(1));
	temp2 = w_w(sort(2));
	w_w(3) = w_w(sort(3));
	w_w(2) = temp2;
	w_w(1) = temp1;
	// sort first two eigen vectors (dont care about third)
	for ( int i = 1; i <= 3; ++i ) {
		temp1 = eVec(i,sort(1));
		temp2 = eVec(i,sort(2));
		eVec(i,1) = temp1;
		eVec(i,2) = temp2;
	}


	// april 20: the fix not only fixes bad eigen vectors but solves a problem of
	// forcing a right-handed coordinate system

	fixEigenvector(eVec);
	// at this point we now have three good eigenvectors in a right hand
	// coordinate system.

	// make bb basis vectors   = moments*eVec

	BlankMatrixMult(m_moment,3,3,0,eVec,3,0,bb);
	//     std::cerr << "m_moment" << std::endl;
	// squirrel away a free copy of the third eigenvector before normalization/fix
	for ( int j = 1; j <= 3; ++j ) {
		Ra(j) = bb(j,3);
	}

	// normalize first two bb-basis vectors
	// dont care about third since were going to replace it with b1xb2
	// this also avoids problem of possible zero third eigen value
	for ( int j = 1; j <= 2; ++j ) {
		temp1 = 1.0/std::sqrt(w_w(j)); // zero checked for above
		for ( int k = 1; k <= 3; ++k ) { // x,y,z
			bb(k,j) *= temp1;
		}
	}

	//  fix things so that bb eigenvecs are right handed

	fixEigenvector(bb); // need to fix this one too
	// find  product of eVec and bb matrices

	BlankMatrixMult(eVec,3,3,0,bb,3,1,UU);
	// result is returned in UU.

	// and lastly determine a value used in another function to compute the rms
	sigma3 = 0.0;
	for ( int j = 1; j <= 3; ++j ) {
		sigma3 += bb(j,3)*Ra(j);
	}
	//cems the std::abs() fixes some round off error situations where the w_w values are
	//cems very small and accidentally negative.  (theoretically they are positive,
	//cems but in practice round off error makes them negative)
	if ( sigma3 < 0.0 ) {
		sigma3 = std::sqrt(std::abs(w_w(1))) + std::sqrt(std::abs(w_w(2))) -
			std::sqrt(std::abs(w_w(3)));
	} else {
		sigma3 = std::sqrt(std::abs(w_w(1))) + std::sqrt(std::abs(w_w(2))) +
			std::sqrt(std::abs(w_w(3)));
	}

} // findUU

/// @brief This is a helper function for using the above implementation of findUU.  There is some cost to the
/// conversion but everything else is probably slower and also you don't have to use FArrays everywhere
void
findUU(
	utility::vector1< numeric::xyzVector<numeric::Real> > & XX,
	utility::vector1< numeric::xyzVector<numeric::Real> > & YY,
	utility::vector1< numeric::Real > const & WW,
	int Npoints,
	numeric::xyzMatrix< numeric::Real > & UU,
	numeric::Real & sigma3
)
{
	FArray2D< numeric::Real > XX_Farray(numeric::vector_of_xyzvectors_to_FArray<numeric::Real>(XX));
	FArray2D< numeric::Real > YY_Farray(numeric::vector_of_xyzvectors_to_FArray<numeric::Real>(YY));
	FArray1D< numeric::Real > WW_Farray(WW.size());
	FArray2D< numeric::Real > UU_Farray(numeric::xyzmatrix_to_FArray<numeric::Real>(UU));

	for ( numeric::Size index = 1; index <= WW.size(); ++index ) {
		WW_Farray(index) = WW[index];
	}

	findUU(XX_Farray,YY_Farray,WW_Farray,Npoints,UU_Farray,sigma3);

	//XX, YY, and UU were updated, convert those back.  WW is const so don't bother

	XX = numeric::FArray_to_vector_of_xyzvectors(XX_Farray);
	YY = numeric::FArray_to_vector_of_xyzvectors(YY_Farray);
	UU = numeric::FArray_to_xyzmatrix(UU_Farray);


}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
///
/// @param  A - [in/out]? -
/// @param  n - [in/out]? -
/// @param  np - [in/out]? -
/// @param  transposeA - [in/out]? -
/// @param  B - [in/out]? -
/// @param  m - [in/out]? -
/// @param  transposeB - [in/out]? -
/// @param  AxB_out - [in/out]? -
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author
///
/////////////////////////////////////////////////////////////////////////////////
void
BlankMatrixMult(
	FArray2A< numeric::Real > A,
	int n,
	int np,
	int transposeA,
	FArray2A< numeric::Real > B,
	int m,
	int transposeB,
	FArray2A< numeric::Real > AxB_out
)
{
	A.dimension( np, n );
	B.dimension( np, m );
	AxB_out.dimension( m, n );

	// fills output matrix with zeros before calling matrix multiply
	AxB_out = 0.0;

	MatrixMult(A,n,np,transposeA,B,m,transposeB,AxB_out);
} // BlankMatrixMult

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
/// multiplys matrices A (npXn) and B (npXn). results in AxB.out
/// IF THE MATRICES are SQUARE.  you can also multiply the transposes of these matrices
/// to do so set the transposeA or transposeB flags to 1, otherwise they should be zero.
///
/// @param  A - [in/out]? -
/// @param  n - [in/out]? -
/// @param  np - [in/out]? -
/// @param  transposeA - [in/out]? -
/// @param  B - [in/out]? -
/// @param  m - [in/out]? -
/// @param  transposeB - [in/out]? -
/// @param  AxB_out - [in/out]? -
///
/// @global_read
///
/// @global_write
///
/// @remarks
/// transpose only works correctly for square matricies!
///
/// this function does SUMS to the old value of AxB_out...
/// you might want to call BlankMatrixMult instead (see above) 4/30/01 jjg
///
///this function works on numeric::Real values
///float version below. jjg
///
/// @references
///
/// @author
/// charlie strauss 1999
///
/////////////////////////////////////////////////////////////////////////////////
void
MatrixMult(
	FArray2A< numeric::Real > A,
	int n,
	int np,
	int transposeA,
	FArray2A< numeric::Real > B,
	int m,
	int transposeB,
	FArray2A< numeric::Real > AxB_out
)
{
	A.dimension( np, n );
	B.dimension( np, m );
	AxB_out.dimension( m, n );


	if ( transposeA == 0 ) {
		if ( transposeB == 0 ) {
			for ( int k = 1; k <= m; ++k ) {
				for ( int j = 1; j <= n; ++j ) {
					for ( int i = 1; i <= np; ++i ) {
						AxB_out(k,j) += A(k,i)*B(i,j);
					}
				}
			}
		} else {
			for ( int k = 1; k <= m; ++k ) {
				for ( int j = 1; j <= n; ++j ) {
					for ( int i = 1; i <= np; ++i ) {
						AxB_out(k,j) += A(k,i)*B(j,i);
					}
				}
			}
		}
	} else {
		if ( transposeB == 0 ) {
			for ( int k = 1; k <= m; ++k ) {
				for ( int j = 1; j <= n; ++j ) {
					for ( int i = 1; i <= np; ++i ) {
						AxB_out(k,j) += A(i,k)*B(i,j);
					}
				}
			}
		} else {
			for ( int k = 1; k <= m; ++k ) {
				for ( int j = 1; j <= n; ++j ) {
					for ( int i = 1; i <= np; ++i ) {
						AxB_out(k,j) += A(i,k)*B(j,i);
					}
				}
			}
		}
	}
} // MatrixMult

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
///       m_v is a 3x3  matrix of 3 eigen vectors
///       replaces the third  eigenvector by taking cross product of
///       of the first two eigenvectors
///
/// @param  m_v - [in/out]? -
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author
///
/////////////////////////////////////////////////////////////////////////////////
void
fixEigenvector( FArray2A< numeric::Real > m_v )
{
	m_v.dimension( 3, 3 );

	numeric::Real const m_v_13 = m_v(2,1)*m_v(3,2) - m_v(3,1)*m_v(2,2);
	numeric::Real const m_v_23 = m_v(3,1)*m_v(1,2) - m_v(1,1)*m_v(3,2);
	numeric::Real const m_v_33 = m_v(1,1)*m_v(2,2) - m_v(2,1)*m_v(1,2);
	//     normalize it to 1 (should already be one but lets be safe)
	numeric::Real const norm = std::sqrt( 1 /
		( ( m_v_13 * m_v_13 ) + ( m_v_23 * m_v_23 ) + ( m_v_33 * m_v_33 ) ) );

	m_v(1,3) = m_v_13 * norm;
	m_v(2,3) = m_v_23 * norm;
	m_v(3,3) = m_v_33 * norm;
} // fixEigenvector

void
rmsfitca2(
	int npoints,
	ObjexxFCL::FArray2A< double > xx,
	ObjexxFCL::FArray2A< double > yy,
	ObjexxFCL::FArray1A< double > ww,
	int natsel,
	double & esq,
	double const & offset_val,
	bool const realign
)
{
	xx.dimension( 3, npoints );
	yy.dimension( 3, npoints );
	ww.dimension( npoints );


	double det;
	int i,j,k;
	double temp1,temp3;
	FArray1D< double > ev( 3 );
	FArray2D< double > m_moment( 3, 3 );
	FArray2D< double > rr_moment( 3, 3 );
	double rms_ctx;
	double rms_sum;
	double handedness;
	FArray1D< double > t( 3 );

	FArray2D< double > R( 3, 3 );
	double XPC, YPC, ZPC, XEC, YEC, ZEC;
	//       //COMMON /TRANSFORM/ XPC,YPC,ZPC,XEC,YEC,ZEC,R

	// align center of mass to origin

	COMAS(xx,ww,npoints,XPC,YPC,ZPC);
	COMAS(yy,ww,npoints,XEC,YEC,ZEC);
	temp3 = 0.0;
	for ( i = 1; i <= npoints; ++i ) {
		temp3 += ww(i);
		// this is outrageous, but there are cases (e.g. in a single-residue pose) where
		// all the z's are at 0, and this makes the det zero.
		xx(3,i) -= offset_val;
		yy(3,i) += offset_val;
	}

	//       Make cross moments matrix   INCLUDE THE WEIGHTS HERE
	for ( k = 1; k <= 3; ++k ) {
		for ( j = 1; j <= 3; ++j ) {
			temp1 = 0.0;
			for ( i = 1; i <= npoints; ++i ) {
				temp1 += ww(i)*yy(k,i)*xx(j,i);
			}
			m_moment(k,j) = temp1    /(temp3); // rescale by temp3
		}
	}
	det = det3(m_moment); // will get handedness  of frame from determinant

	if ( std::abs(det) <= 1.0E-24 ) {
		//     //  std::cerr << "Warning:degenerate cross moments: det=" << det << std::endl;
		//     // might think about returning a zero rms, to avoid any chance of Floating Point Errors?

		esq = 0.0;
		return;

	}
	handedness = numeric::sign_transfered(det, 1.0);
	//  // weird but documented fortran "feature" of sign(a,b) (but not SIGN) is that if fails if a < 0

	//  //  multiply cross moments by itself

	for ( i = 1; i <= 3; ++i ) {
		for ( j = i; j <= 3; ++j ) {
			rr_moment(j,i) = rr_moment(i,j) = // well it is symmetric afterall
				m_moment(1,i)*m_moment(1,j) +
				m_moment(2,i)*m_moment(2,j) +
				m_moment(3,i)*m_moment(3,j);
		}
	}

	//            //  compute eigen values of cross-cross moments

	rsym_eigenval(rr_moment,ev);

	//               // reorder eigen values  so that ev(3) is the smallest eigenvalue

	if ( ev(2) > ev(3) ) {
		if ( ev(3) > ev(1) ) {
			temp1 = ev(3);
			ev(3) = ev(1);
			ev(1) = temp1;
		}
	} else {
		if ( ev(2) > ev(1) ) {
			temp1 = ev(3);
			ev(3) = ev(1);
			ev(1) = temp1;
		} else {
			temp1 = ev(3);
			ev(3) = ev(2);
			ev(2) = temp1;
		}
	}

	//                 // ev(3) is now the smallest eigen value.  the other two are not
	//                 //  sorted.  this is prefered order for rotation matrix


	rsym_rotation(m_moment,rr_moment,ev,R);

	//$$$             for ( i = 1; i <= npoints; ++i ) {
	//$$$               for ( j = 1; j <= 3; ++j ) {
	//$$$                 temp1 = 0.0;
	//$$$                for ( k = 1; k <= 3; ++k ) {
	//$$$                  temp1 += R(j,k)*yy(k,i);
	//$$$                }
	//$$$                t(j) = temp1;
	//$$$               }
	//$$$               yy(1,i) = t(1);
	//$$$               yy(2,i) = t(2);
	//$$$               yy(3,i) = t(3);
	//$$$             }

	for ( i = 1; i <= npoints; ++i ) {
		if ( realign ) { //Undo the offset:
			xx(3,i) += offset_val;
			yy(3,i) -= offset_val;
		}
		for ( j = 1; j <= 3; ++j ) { // compute rotation
			t(j) = R(j,1)*yy(1,i) + R(j,2)*yy(2,i) + R(j,3)*yy(3,i);
		}
		yy(1,i) = t(1);
		yy(2,i) = t(2);
		yy(3,i) = t(3);
	}
	//   // now we must catch the special case of the rotation with inversion.
	//   // we cannot allow inversion rotations.
	//   // fortunatley, and curiously, the optimal non-inverted rotation matrix
	//   // will have the similar eigen values.
	//   // we just have to make a slight change in how we handle things depending on determinant

	rms_ctx = std::sqrt(std::abs(ev(1))) + std::sqrt(std::abs(ev(2))) +
		handedness*std::sqrt(std::abs(ev(3)));

	rms_ctx *= temp3;

	//   // the std::abs() are theoretically unneccessary since the eigen values of a real symmetric
	//   // matrix are non-negative.  in practice sometimes small eigen vals end up just negative
	rms_sum = 0.0;
	for ( i = 1; i <= npoints; ++i ) {
		for ( j = 1; j <= 3; ++j ) {
			rms_sum += ww(i)*( ( yy(j,i) * yy(j,i) ) + ( xx(j,i) * xx(j,i) ) );
		}
	}
	// rms_sum = rms_sum; //   /temp3   (will use natsel instead)

	//  // and combine the outer and cross terms into the final calculation.
	//  //  (the std::abs() just saves us a headache when the roundoff error accidantally makes the sum negative)

	esq = std::sqrt( std::abs( rms_sum - ( 2.0 * rms_ctx ) ) / natsel );

} // rmsfitca2


void
rmsfitca3(
	int npoints, // number of points to fit
	ObjexxFCL::FArray2A< double > xx0,
	ObjexxFCL::FArray2A< double > xx,
	ObjexxFCL::FArray2A< double > yy0,
	ObjexxFCL::FArray2A< double > yy,
	double & esq
)
{

	xx0.dimension( 3, npoints );
	xx.dimension ( 3, npoints );
	yy0.dimension( 3, npoints );
	yy.dimension ( 3, npoints );

	// local
	double det;
	double temp1,mass;
	FArray1D< double > come( 3 );
	FArray1D< double > comp( 3 );
	FArray1D< double > ev( 3 );
	FArray2D< double > m_moment( 3, 3 );
	FArray2D< double > rr_moment( 3, 3 );
	double rms_ctx,rms2_sum;
	double handedness;

	FArray2D< double > r( 3, 3 );
	FArray1D< double > t( 3 );


	// compute center of mass
	numeric::model_quality::RmsData* rmsdata = RmsData::instance(); // get a pointer to the singleton class

	mass                   = rmsdata->count();
	double xre             = rmsdata->xre();
	double xrp             = rmsdata->xrp();
	FArray1D< double > xse = rmsdata->xse();
	FArray1D< double > xsp = rmsdata->xsp();
	FArray2D< double > xm = rmsdata->xm();

	// std::cerr << mass << " " << xre << " " << xrp << " " << xse(1) << " " << xse(2) << " " << xse(3) << " "
	//     << xsp(1) << " " << xsp(2) << " " << xsp(3) << " " << xm(1,1) << " " << xm(1, 2) << std::endl;


	come(1) = xse(1)/mass; // x_com
	come(2) = xse(2)/mass; // y_com
	come(3) = xse(3)/mass; // z_com

	comp(1) = xsp(1)/mass; // x_com
	comp(2) = xsp(2)/mass; // y_com
	comp(3) = xsp(3)/mass; // z_com


	//       Make cross moments matrix

	for ( int k = 1; k <= 3; ++k ) {
		for ( int j = 1; j <= 3; ++j ) {
			m_moment(k,j) = xm(k,j)/mass - come(k)*comp(j); // flopped com
		}
	}

	det = det3(m_moment); // get handedness  of frame from determinant

	if ( std::abs(det) <= 1.0E-24 ) {
		//   //std::cerr << "Warning:degenerate cross moments: det=" << det << std::endl;
		//   // might think about returning a zero rms, to avoid any chance of
		//   // Floating Point Errors?

		esq = 0.0;
		return;

	}
	handedness = numeric::sign_transfered(det, 1.0); // changed name of call from sign to sign_transfered in mini!
	/// OL and order of arguments -- damn

	// std::cerr << handedness << std::endl;
	//    // weird but documented "feature" of sign(a,b) (but not SIGN) is
	//    // that if fails if a < 0

	//    //  multiply cross moments by itself

	for ( int i = 1; i <= 3; ++i ) {
		for ( int j = i; j <= 3; ++j ) {
			rr_moment(j,i) = rr_moment(i,j) =
				m_moment(1,i)*m_moment(1,j) +
				m_moment(2,i)*m_moment(2,j) +
				m_moment(3,i)*m_moment(3,j);
		}
	}

	//    // compute eigen values of cross-cross moments

	rsym_eigenval(rr_moment,ev);


	//    // reorder eigen values  so that ev(3) is the smallest eigenvalue

	if ( ev(2) > ev(3) ) {
		if ( ev(3) > ev(1) ) {
			temp1 = ev(3);
			ev(3) = ev(1);
			ev(1) = temp1;
		}
	} else {
		if ( ev(2) > ev(1) ) {
			temp1 = ev(3);
			ev(3) = ev(1);
			ev(1) = temp1;
		} else {
			temp1 = ev(3);
			ev(3) = ev(2);
			ev(2) = temp1;
		}
	}

	//     // ev(3) is now the smallest eigen value.  the other two are not
	//     //  sorted.  This is prefered order for computing the rotation matrix

	rsym_rotation(m_moment,rr_moment,ev,r);

	//            // now we rotate and offset all npoints


	for ( int i = 1; i <= npoints; ++i ) {
		for ( int k = 1; k <= 3; ++k ) { // remove center of mass
			yy(k,i) = yy0(k,i)-come(k);
			xx(k,i) = xx0(k,i)-comp(k);
		}
		for ( int j = 1; j <= 3; ++j ) { // compute rotation
			//               // temp1 = 0.0;
			//               // for ( k = 1; k <= 3; ++k ) {
			//               //  temp1 += r(j,k)*yy(k,i);
			//               // }
			//               // t(j) = temp1;
			t(j) = r(j,1)*yy(1,i) + r(j,2)*yy(2,i) + r(j,3)*yy(3,i);

		}
		yy(1,i) = t(1);
		yy(2,i) = t(2);
		yy(3,i) = t(3);
	}

	//            // now we must catch the special case of the rotation with inversion.
	//            // fortunatley, and curiously, the optimal non-inverted rotation
	//            // matrix will have a similar relation between rmsd and the eigen values.
	//            // we just have to make a slight change in how we handle things
	//            // depending on determinant

	rms_ctx = std::sqrt(std::abs(ev(1))) + std::sqrt(std::abs(ev(2))) +
		handedness*std::sqrt(std::abs(ev(3)));
	// std::cerr << handedness << std::endl;
	//            // the std::abs() are theoretically unneccessary since the eigen values
	//            // of a real symmetric matrix are non-negative.
	//            // in practice sometimes small eigen vals end up as tiny negatives.

	temp1 = ( come(1) * come(1) ) + ( come(2) * come(2) ) + ( come(3) * come(3) ) +
		( comp(1) * comp(1) ) + ( comp(2) * comp(2) ) + ( comp(3) * comp(3) );

	rms2_sum = (xre + xrp)/mass - temp1;

	//            // and combine the outer and cross terms into the final calculation.
	//            //  (the std::abs() just saves us a headache when the roundoff error
	//            // accidantally makes the sum negative)

	esq = std::sqrt(std::abs(rms2_sum-2.0*rms_ctx));

} // rmsfitca3

double
det3( FArray2A< double > m )
{

	m.dimension( 3, 3 );

	return
		m(1,3)*( m(2,1)*m(3,2) - m(2,2)*m(3,1) ) -
		m(2,3)*( m(1,1)*m(3,2) - m(1,2)*m(3,1) ) +
		m(3,3)*( m(1,1)*m(2,2) - m(1,2)*m(2,1) );
}


void
rsym_eigenval(
	ObjexxFCL::FArray2A< double > m,
	ObjexxFCL::FArray1A< double > ev
)
{
	m.dimension( 3, 3 );
	ev.dimension( 3 );


	double xx,yy,zz,xy,xz,yz;
	double a,b,c,s0;
	std::complex< double > f1,f2,f3,f4,f5;

	static std::complex< double > const unity = std::complex< double >(1.0,0.0);
	static std::complex< double > const sqrt_3i =
		std::sqrt( 3.0 ) * std::complex< double >(0.0,1.0);

	// first, for lexical sanity, name some temporary variables
	xx = m(1,1);
	yy = m(2,2);
	zz = m(3,3);
	xy = m(1,2);
	xz = m(1,3);
	yz = m(2,3);

	// coefficients of characterisitic polynomial
	a = xx+yy+zz;
	b = -xx*zz-xx*yy-yy*zz+xy*xy+xz*xz+yz*yz;
	c = xx*yy*zz-xz*xz*yy-xy*xy*zz-yz*yz*xx+2*xy*xz*yz;

	// For numerical dynamic range we rescale the variables here
	// with rare exceptions this is unnessary but it doesn't add much to the calculation time.
	// this also allows use of double if desired. cems
	//  for complex< float > the rescaling trigger  should be  1e15  (or less)
	//  for complex< double > the rescaling trigger should be  1e150 (or less)
	// note we completely ignore the possiblity of needing rescaling to avoid
	// underflow due to too small numbers. left as an excercise to the reader.

	double norm = std::max(std::abs(a),std::max(std::abs(b),std::abs(c)));
	if ( norm > 1.0E50 ) {   // rescaling trigger
		a /= norm;
		b /= norm * norm;
		c /= norm * norm * norm;
	} else {
		norm = 1.0;
	}
	// we undo the scaling by de-scaling the eigen values at the end

	// Power constants
	double const a2 = a * a;
	double const a3 = a * a * a;

	// eigenvals are the roots of the characteristic polymonial  0 = c + b*e + a*e^2 - e^3
	// solving for the three roots now:
	// dont try to follow this in detail: its just a tricky
	// factorization of the formulas for cubic equation roots.

	s0 = ( -12.0 * ( b * b * b ) ) - ( 3.0 * ( b * b ) * a2 ) +
		( 54.0 * c * b * a ) + ( 81.0 * ( c * c ) ) + ( 12.0 * c * a3 );
	// butt ugly term

	f1 = b*a/6.0 + c/2.0 + a3/27.0 + std::sqrt(s0*unity)/18.0;

	f2 = std::pow( f1, (1.0/3.0) ); // note f1 is a complex number

	f3 = (-b/3.0 - a2/9.0)/f2;

	f4 = f2-f3;

	f5 = sqrt_3i * (f2+f3); // just our imaginary friend, mr. i

	s0 = a/3.0;
	ev(1) = f4.real();
	// note implicitly take real part, imag part "should" be zero
	ev(1) = norm*(ev(1)+s0 );
	// do addition after type conversion in previous line.
	ev(2) = (-f4+f5).real(); // note real part, imag part is zero
	ev(2) = norm*(ev(2)*0.5 + s0 );
	ev(3) = (-f4-f5).real(); // note real part, imag part is zero
	ev(3) = norm*(ev(3)*0.5 + s0);

}

void
rsym_rotation(
	ObjexxFCL::FArray2A< double > mm,
	ObjexxFCL::FArray2A< double > m,
	ObjexxFCL::FArray1A< double > ev,
	ObjexxFCL::FArray2A< double > rot
)
{
	mm.dimension( 3, 3 );
	m.dimension( 3, 3 );
	ev.dimension( 3 );
	rot.dimension( 3, 3 );


	FArray2D< double > temp( 3, 3 );
	FArray2D< double > mvec( 3, 3 );

	rsym_evector(m,ev,mvec);

	for ( int i = 1; i <= 2; ++i ) { // dont need no stinkin third component
		double norm = 1.0 / std::sqrt(std::abs(ev(i) )); // abs just fixes boo boos
		// if one was nervous here one could explicitly
		// compute the norm of temp.
		for ( int j = 1; j <= 3; ++j ) {
			temp(j,i) = 0.0;
			for ( int k = 1; k <= 3; ++k ) {
				temp(j,i) += mvec(k,i)*mm(j,k);
			}
			temp(j,i) *= norm;
		}
	}

	temp(1,3) =  temp(2,1)*temp(3,2) - temp(2,2)*temp(3,1);
	temp(2,3) = -temp(1,1)*temp(3,2) + temp(1,2)*temp(3,1);
	temp(3,3) =  temp(1,1)*temp(2,2) - temp(1,2)*temp(2,1);

	for ( int i = 1; i <= 3; ++i ) {
		for ( int j = 1; j <= 3; ++j ) {
			rot(j,i) = 0.0;
			for ( int k = 1; k <= 3; ++k ) {
				rot(j,i) += temp(i,k)*mvec(j,k);
			}
		}
	}
}


void
rsym_evector(
	ObjexxFCL::FArray2A< double > m,
	ObjexxFCL::FArray1A< double > ev,
	ObjexxFCL::FArray2A< double > mvec
)
{
	m.dimension( 3, 3 );
	ev.dimension( 3 );
	mvec.dimension( 3, 3 );

	// local
	double xx,yy,xy,zx,yz; //zz
	double e1,e2,e3,znorm;


	// first, for sanity only, name some temporary variables
	//zz = m(3,3);
	xy = m(1,2);
	zx = m(1,3);
	yz = m(2,3);

	if ( ev(1) != ev(2) ) {  // test for degenerate eigen values

		for ( int i = 1; i <= 2; ++i ) {
			// only computer first two eigen vectors using this method
			// note you could compute all three this way if you wanted to,
			// but you would run into problems with degenerate eigen values.

			xx = m(1,1)-ev(i);
			yy = m(2,2)-ev(i);
			// I marvel at how simple this is when you know the eigen values.
			e1 = xy*yz-zx*yy;
			e2 = xy*zx-yz*xx;
			e3 = xx*yy-xy*xy;

			znorm = std::sqrt( ( e1 * e1 ) + ( e2 * e2 ) + ( e3 * e3 ) );

			mvec(1,i) = e1/znorm;
			mvec(2,i) = e2/znorm;
			mvec(3,i) = e3/znorm;

		}

		// now compute the third eigenvector
		mvec(1,3) =  mvec(2,1)*mvec(3,2) - mvec(2,2)*mvec(3,1);
		mvec(2,3) = -mvec(1,1)*mvec(3,2) + mvec(1,2)*mvec(3,1);
		mvec(3,3) =  mvec(1,1)*mvec(2,2) - mvec(1,2)*mvec(2,1);

		// pathologically nervous people would explicitly normalize this vector too.

		return;

	} else {

		if ( ev(2) != ev(3) ) {
			std::cerr << " hey is this the right thing to be doing??? " << std::endl;

			for ( int i = 2; i <= 3; ++i ) {
				// Okay, since 1 and 2 are degenerate we will use 2 and 3 instead.

				xx = m(1,1)-ev(i);
				yy = m(2,2)-ev(i);
				// I marvel at how simple this is when you know the eigen values.
				e1 = xy*yz-zx*yy;
				e2 = xy*zx-yz*xx;
				e3 = xx*yy-xy*xy;
				// yes you sharp eyed person, its not quite symmetric here too.
				//                   life is odd.

				znorm = std::sqrt( ( e1 * e1 ) + ( e2 * e2 ) + ( e3 * e3 ) );

				mvec(1,i) = e1/znorm;
				mvec(2,i) = e2/znorm;
				mvec(3,i) = e3/znorm;

			}

			// now compute the third eigenvector
			mvec(1,1) =  mvec(2,2)*mvec(3,3) - mvec(2,3)*mvec(3,2);
			mvec(2,1) = -mvec(1,2)*mvec(3,3) + mvec(1,3)*mvec(3,2);
			mvec(3,1) =  mvec(1,2)*mvec(2,3) - mvec(1,3)*mvec(2,2);

			// pathologically nervous people would explicitly normalize this vector too.

			return;

		} else {

			std::cerr << "warning: all eigen values are equal" << std::endl;

			for ( int i = 1; i <= 3; ++i ) {
				mvec(1,i) = 0.0;
				mvec(2,i) = 0.0;
				mvec(3,i) = 0.0;
				mvec(i,i) = 1.0;
			}
			return;
		}
	}
} // rsym_evector


} // rms
} // numeric
