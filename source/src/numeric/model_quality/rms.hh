// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/model_quality/rms.hh
/// @brief  RMS stuff imported from rosetta++
/// @author James Thompson
/// @date   Wed Aug 22 12:10:37 2007


#ifndef INCLUDED_numeric_model_quality_rms_HH
#define INCLUDED_numeric_model_quality_rms_HH

#include <numeric/types.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>

#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray1A.hh>
#include <ObjexxFCL/FArray2A.hh>

#include <utility/vector1.hh>

namespace numeric {
namespace model_quality {

numeric::Real
calc_rms(
	utility::vector1< xyzVector< Real > > p1_coords,
	utility::vector1< xyzVector< Real > > p2_coords
);

numeric::Real
rms_wrapper (
	int natoms,
	ObjexxFCL::FArray2D< numeric::Real > p1a,
	ObjexxFCL::FArray2D< numeric::Real > p2a
);

numeric::Real
rms_wrapper_slow_and_correct(
	int natoms,
	ObjexxFCL::FArray2D< numeric::Real > p1a,
	ObjexxFCL::FArray2D< numeric::Real > p2a
);

void
BlankMatrixMult(
	ObjexxFCL::FArray2A< numeric::Real > A,
	int n,
	int np,
	int transposeA,
	ObjexxFCL::FArray2A< numeric::Real > B,
	int m,
	int transposeB,
	ObjexxFCL::FArray2A< numeric::Real > AxB_out
);

void
MatrixMult(
	ObjexxFCL::FArray2A< numeric::Real > A,
	int n,
	int np,
	int transposeA,
	ObjexxFCL::FArray2A< numeric::Real > B,
	int m,
	int transposeB,
	ObjexxFCL::FArray2A< numeric::Real > AxB_out
);

void
fixEigenvector( ObjexxFCL::FArray2A< numeric::Real > m_v );

void
calc_rms_fast(
	float & rms_out,
	ObjexxFCL::FArray2A< numeric::Real > xx,
	ObjexxFCL::FArray2A< numeric::Real > yy,
	ObjexxFCL::FArray1A< numeric::Real > ww,
	int npoints,
	numeric::Real ctx
);

void
findUU(
	utility::vector1< numeric::xyzVector<numeric::Real> > & XX,
	utility::vector1< numeric::xyzVector<numeric::Real> > & YY,
	utility::vector1< numeric::Real > const & WW,
	int Npoints,
	numeric::xyzMatrix< numeric::Real > & UU,
	numeric::Real & sigma3
);

void
findUU(
	ObjexxFCL::FArray2< numeric::Real > & XX,
	ObjexxFCL::FArray2< numeric::Real > & YY,
	ObjexxFCL::FArray1< numeric::Real > const & WW,
	int Npoints,
	ObjexxFCL::FArray2< numeric::Real > & UU,
	numeric::Real & sigma3
);

void
MatrixMult(
	ObjexxFCL::FArray2A< numeric::Real > A,
	int n,
	int np,
	int transposeA,
	ObjexxFCL::FArray2A< numeric::Real > B,
	int m,
	int transposeB,
	ObjexxFCL::FArray2A< numeric::Real > AxB_out
);

// functions required for maxsub to work properly


////////////////////////////////////////////////////////////////////////////////
///
/// @brief determinant of a 3x3 matrix
///
/// @details
/// cute factoid: det of a 3x3 is the dot product of one row with the cross
/// product of the other two. This explains why a right hand coordinate system
/// has a positive determinant. cute huh?
///
/// @param  m - [in/out]? -
///
/// @return
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author charlie strauss 2001
///
/////////////////////////////////////////////////////////////////////////////////
double
det3( ObjexxFCL::FArray2A< double > m );


////////////////////////////////////////////////////////////////////////////////
///
/// @brief computes the rms between two weighted point vectors.
///
/// @details
///   xx_0,yy_0 are the input vectors of of points and ww is their weights
///   xx,yy are by product output vectors of the same points offset to remove center of mass
///   det is an out value of the determinant of the cross moment matrix
///   returned value is the rms
///
///   most of this is double for good reasons.  first there are some large
///   differences of small numbers.  and second the rsymm_eignen() function can
///   internally have numbers larger than the largest double number.
///   (you could do some fancy foot work to rescale things if you really
///   had a problem with this.)
///
/// @param  npoints - [in/out]? -
/// @param  xx - [in/out]? -
/// @param  yy - [in/out]? -
/// @param  ww - [in/out]? -
/// @param  natsel - [in/out]? -
/// @param  esq - [in/out]? -
/// @param[in] offset_val - A small offset temporarily added or subtracted from z-coordinates to ensure that determinants are nonzero.  Defaults to 1e-7.
/// @param[in] relalign - If true, the small offset is subtracted off again before final output.  False by default (legacy behaviour).
///
/// @global_read
///
/// @global_write
///
/// @remarks
///   det is a double precision real
///  (xx,yy) can be same arrays as (xx_0,yy_0) if desired
///
/// @references
///
/// @author charlie strauss 2001
///
/////////////////////////////////////////////////////////////////////////////////
void
rmsfitca2(
	int npoints,
	ObjexxFCL::FArray2A< double > xx,
	ObjexxFCL::FArray2A< double > yy,
	ObjexxFCL::FArray1A< double > ww,
	int natsel,
	double & esq,
	double const & offset_val=1.0e-7,
	bool const realign = false
);

//////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
///   This function gets its alignment info via a namespace!
///   Alignment (rotation matrix) and rms(esq) are computed on the basis
///   of residues previously designated by calls to add_rms().
///   However, the rotation is applied to all Npoints of XX0,yy0 with the
///   results returned in xx,yy.
///
///   most of this is double for good reasons.
///   first there are some large differences of small numbers.
///   second the rsymm_eignen() function can internally have numbers
///   larger than the largest double number.  (you could do some fancy foot work
///   to rescale m_moment if you really had a problem with this.)
///
/// @param  npoints - [in/out]? -
/// @param  xx0 - [in/out]? -
/// @param  xx - [in/out]? -
/// @param  yy0 - [in/out]? -
/// @param  yy - [in/out]? -
/// @param  esq - [in/out]? -
///
/// @global_read
///
/// @global_write
///
/// @remarks
///   NOTE: det is a double precision real
///   NOTE: (xx,yy) can be same arrays as (xx_0,yy_0) if desired
///
/// @references
///
/// @author
///
////////////////////////////////////////////////////////////////////////////////
void
rmsfitca3(
	int npoints, // number of points to fit
	ObjexxFCL::FArray2A< double > xx0,
	ObjexxFCL::FArray2A< double > xx,
	ObjexxFCL::FArray2A< double > yy0,
	ObjexxFCL::FArray2A< double > yy,
	double & esq
);


////////////////////////////////////////////////////////////////////////////////
///
/// @brief computes the eigen values of a real symmetric 3x3 matrix
///
/// @details
/// the method used is a deterministic analytic result that I hand factored.(whew!)
/// Amusingly, while I suspect this factorization is not yet optimal in the
/// number of calcs required
/// I cannot find a signifcantly better one.
///  (if it were optimal I suspect I would not have to compute imaginary numbers that
///  I know must eventually cancel out to give a net real result.)
/// this method relys on the fact that an analytic factoring of an order 3 polynomial exists.
/// m(3,3) is a 3x3 real symmetric matrix: only the upper triangle is actually used
///  ev(3) is a real vector of eigen values, not neccesarily in sorted order.
///
/// @param  m - [in/out]? -
/// @param  ev - [in/out]? -
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author charlie strauss 2001
///
/////////////////////////////////////////////////////////////////////////////////
void
rsym_eigenval(
	ObjexxFCL::FArray2A< double > m,
	ObjexxFCL::FArray1A< double > ev
);

////////////////////////////////////////////////////////////////////////////////
///
/// @brief finds a (proper) rotation matrix that minimizes the rms.
///
/// @details this computes the rotation matrix based on the eigenvectors of m that gives the
/// the mimimum rms.  this is determined using mm the cross moments matrix.
///
/// for best results the third eigen value should be the
/// smallest eigen value!
///
///
/// @param  mm - [in/out]? -
/// @param  m - [in/out]? -
/// @param  ev - [in/out]? -
/// @param  rot - [in/out]? -
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author charlie strauss (cems) 2001
///
/////////////////////////////////////////////////////////////////////////////////
void
rsym_rotation(
	ObjexxFCL::FArray2A< double > mm,
	ObjexxFCL::FArray2A< double > m,
	ObjexxFCL::FArray1A< double > ev,
	ObjexxFCL::FArray2A< double > rot
);

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
/// Author: charlie strauss (cems) 2001
/// given original matrix plus its eigen values
/// compute the eigen vectors.
/// USAGE notice: for computing min rms rotations be sure to call this
/// with the lowest eigen value in Ev(3).
///
/// The minimal factorization of the eigenvector problem I have derived below has a puzzling
/// or, rather, interesting assymetry.  Namely, it doesn't matter what
///  either ZZ=M(3,3) is or what Ev(3) is!
///  One reason for this is that of course all we need to  know is contained in the
///  M matrix in the firstplace, so the eigen values overdetermine the problem
///  We are just exploiting the redundant info in the eigen values to hasten the solution.
///  What I dont know is if infact there exists any symmetic form using the eigen values.
///
/// we deliberately introduce another assymetry for numerical stability, and
/// to force proper rotations  (since eigen vectors are not unique within a sign change).
///   first, we explicitly numerically norm the vectors to 1 rather than assuming the
///   algebra will do it with enough accuracy. ( and its faster to boot!)
///   second, we explicitly compute the third vector as being the cross product of the
///   previous two rather than using the M matrix and third eigen value.
///   If you arrange the eigen values so that the third eigen value is the lowest
///   then this guarantees a stable result under the case of either very small
///   eigen value or degenerate eigen values.
///   this norm, and ignoring the third eigen value also gaurentee us the even if the
///   eigen vectors are not perfectly accurate that, at least, the matrix
///   that results is a pure orthonormal rotation matrix, which for most applications
///   is the most important form of accuracy.
///
/// @details
///
/// @param  m - [in/out]? -
/// @param  ev - [in/out]? -
/// @param  mvec - [in/out]? -
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
rsym_evector(
	ObjexxFCL::FArray2A< double > m,
	ObjexxFCL::FArray1A< double > ev,
	ObjexxFCL::FArray2A< double > mvec
);


} // end namespace model_quality
} // end namespace numeric

#endif
