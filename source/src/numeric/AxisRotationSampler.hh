// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/numeric/AxisRotationSampler.hh
/// @author Ryan Pavlovicz
// Uniformly sample rotation space about an axis
// modified from UniformRotationSampler

#ifndef INCLUDED_numeric_axisrotationsampler_hh
#define INCLUDED_numeric_axisrotationsampler_hh

#include <numeric/conversions.hh>
#include <numeric/constants.hh>
#include <numeric/Quaternion.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/angle.functions.hh>
#include <numeric/types.hh>

namespace numeric {

class AxisRotationSampler {
private:
	utility::vector1< xyzMatrix<Real> > rotlist_;
	Real theta_;

public:
	Size
	nrots() const { return rotlist_.size(); }

	void get(Size ii, xyzMatrix<Real> &R) const {
		if ( ii==0 ) {
			R = xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1);
		} else {
			R = rotlist_[ii];
		}
	}

	AxisRotationSampler( xyzVector<Real> axis, Real theta ) {
		theta_ = theta;
		if ( theta<=1e-6 ) {
			rotlist_.resize(1);
			return;
		}

		int nrot = (int)std::floor(360.0 / theta + 0.5);
		rotlist_.resize(nrot);

		Quaternion<Real> quat;
		xyzMatrix<Real> rot_matrix;
		for ( int j=0; j<nrot; ++j ) {
			Real rot_deg = j*theta;
			quat = Quaternion<Real>( axis.normalize(), conversions::radians(rot_deg), true );
			quat2R( quat, rot_matrix );
			rotlist_[j+1] = rot_matrix;
		}
	}

}; // class AxisRotationSampler

} // namespace numeric

#endif
