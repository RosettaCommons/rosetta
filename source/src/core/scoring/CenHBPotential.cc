// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/CenHBPotential.cc
/// @brief  Smooth, differentiable version of centroid hbond term
/// @author Frank DiMaio


#include <core/scoring/CenHBPotential.hh>


#include <basic/database/open.hh>
#include <utility/io/izstream.hh>

#include <utility/vector1.hh>
#include <boost/bind.hpp>


namespace core {
namespace scoring {

///////////////////////////////////////////////////////////////////////////////////////////////

Real CenHBPotential::func( Size seqsep, Real d, Real xd, Real xh ) const {
	using numeric::constants::f::pi;

	utility::vector1< Real > const & As       = (seqsep<=4) ? sr_As_ : lr_As_;
	utility::vector1< Vector > const & mus    = (seqsep<=4) ? sr_mus_ : lr_mus_;
	utility::vector1< Vector > const & sigmas = (seqsep<=4) ? sr_sigmas_ : lr_sigmas_;

	Real y = 0;
	for (Size i=1; i<=As.size(); ++i) {
		y += As[i] * exp( -(d-mus[i][0])*(d-mus[i][0]) / (2*sigmas[i][0]*sigmas[i][0]) )
		           * ( exp( sigmas[i][1]*cos( pi/180*(xd-mus[i][1]) ) ) + exp( sigmas[i][1]*cos( pi/180*(xd+mus[i][1]) ) ) )
		           * ( exp( sigmas[i][2]*cos( pi/180*(xh-mus[i][2]) ) ) + exp( sigmas[i][2]*cos( pi/180*(xh+mus[i][2]) ) ) );
	}
	return y;
}


Vector CenHBPotential::dfunc( Size seqsep, Real d, Real xd, Real xh ) const {
	using numeric::constants::f::pi;

	utility::vector1< Real > const & As       = (seqsep<=4) ? sr_As_ : lr_As_;
	utility::vector1< Vector > const & mus    = (seqsep<=4) ? sr_mus_ : lr_mus_;
	utility::vector1< Vector > const & sigmas = (seqsep<=4) ? sr_sigmas_ : lr_sigmas_;

	Vector dy = Vector(0,0,0);
	for (Size i=1; i<=As.size(); ++i) {
			Real s1 = exp( -(d-mus[i][0])*(d-mus[i][0]) / (2*sigmas[i][0]*sigmas[i][0]) );
			Real s2a = exp( sigmas[i][1]*cos( pi/180*(xd-mus[i][1]) ) );
			Real s2b = exp( sigmas[i][1]*cos( pi/180*(xd+mus[i][1]) ) );
			Real s3a = exp( sigmas[i][2]*cos( pi/180*(xh-mus[i][2]) ) );
			Real s3b = exp( sigmas[i][2]*cos( pi/180*(xh+mus[i][2]) ) );
			Real s2 = s2a+s2b;
			Real s3 = s3a+s3b;

			dy[0] += -As[i] * s1 * s2 * s3 * (d-mus[i][0]) / (sigmas[i][0]*sigmas[i][0]);
			dy[1] += -As[i] * s1 * s3 * pi/180 * sigmas[i][1] * ( sin( pi/180*(xd-mus[i][1])) * s2a + sin( pi/180*(xd+mus[i][1])) * s2b );
			dy[2] += -As[i] * s1 * s2 * pi/180 * sigmas[i][2] * ( sin( pi/180*(xh-mus[i][2])) * s3a + sin( pi/180*(xh+mus[i][2])) * s3b );
	}
	return dy;
}

///////////////////////////////////////////////////////////////////////////////////////////////

CenHBPotential::CenHBPotential() {
	// defaults
	cutoff_sr_ = 6.0;
	cutoff_lr_ = 6.0;

	// load the data
	std::string tag,line;
	utility::io::izstream stream;
	basic::database::open( stream, "scoring/score_functions/centroid_smooth/cen_hb_params.txt");

	while ( getline( stream, line ) ) {
		std::istringstream l(line);
		l >> tag;
		if (tag == "HBOND_BB_LR:") {
			l >> tag;
			if (tag == "CUTOFF") {
				l >> cutoff_lr_;
			} else if (tag == "GAUSSIAN3D") {
				Size ngauss; l >> ngauss;
				Real A;
				numeric::xyzVector<Real> mu, sigma;
				for (Size i=1; i<=ngauss; ++i) {
					l >> A >> mu[0] >> mu[1] >> mu[2] >> sigma[0] >> sigma[1] >> sigma[2];
					add_lr_gaussian(A,mu,sigma);
				}
			}
		} else if (tag == "HBOND_BB_SR:") {
			l >> tag;
			if (tag == "CUTOFF") {
				l >> cutoff_sr_;
			} else if (tag == "GAUSSIAN3D") {
				Size ngauss; l >> ngauss;
				Real A;
				numeric::xyzVector<Real> mu, sigma;
				for (Size i=1; i<=ngauss; ++i) {
					l >> A >> mu[0] >> mu[1] >> mu[2] >> sigma[0] >> sigma[1] >> sigma[2];
					add_sr_gaussian(A,mu,sigma);
				}
			}
		} else if (tag != "#") {
			utility_exit_with_message("bad format for cen_smooth_params.txt");
		}

		if ( l.fail() ) utility_exit_with_message("bad format for cen_smooth_params.txt");
	}
}


}
}
