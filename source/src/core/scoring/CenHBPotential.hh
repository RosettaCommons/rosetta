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

#ifndef INCLUDED_core_scoring_CenHBPotential_hh
#define INCLUDED_core_scoring_CenHBPotential_hh

#include <core/scoring/CenHBPotential.fwd.hh>

#include <core/types.hh>

#include <numeric/constants.hh>


#include <utility/vector1_bool.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>


namespace core {
namespace scoring {

////////////////////////
// fpd helper class holds some # of gaussians
class CenHBPotential : public utility::pointer::ReferenceCount {
public:
	CenHBPotential();

	Size nlr_gaussians() const { return lr_As_.size(); }
	Size nsr_gaussians() const { return lr_As_.size(); }

	void clear() {
		lr_As_.clear( );
		lr_mus_.clear( );
		lr_sigmas_.clear( );
		sr_As_.clear( );
		sr_mus_.clear( );
		sr_sigmas_.clear( );
	}

	void add_sr_gaussian(
		Real A_in,
		numeric::xyzVector< Real > mu_in,
		numeric::xyzVector< Real > sigma_in ) {
		using numeric::constants::f::pi;

		// normalize terms
		A_in /= sqrt(2*pi*sigma_in[0]*sigma_in[0]);
		A_in /= 2*pi*BesselI0(sigma_in[1]);
		A_in /= 2*pi*BesselI0(sigma_in[2]);

		sr_As_.push_back( A_in );
		sr_mus_.push_back( mu_in );
		sr_sigmas_.push_back( sigma_in );
	}

	void add_lr_gaussian(
		Real A_in,
		numeric::xyzVector< Real > mu_in,
		numeric::xyzVector< Real > sigma_in ) {
		using numeric::constants::f::pi;

		// normalize terms
		A_in /= sqrt(2*pi*sigma_in[0]*sigma_in[0]);
		A_in /= 2*pi*BesselI0(sigma_in[1]);
		A_in /= 2*pi*BesselI0(sigma_in[2]);

		lr_As_.push_back( A_in );
		lr_mus_.push_back( mu_in );
		lr_sigmas_.push_back( sigma_in );
	}

	void set_cutoff_sr(Real cut_in) { cutoff_sr_ = cut_in; }
	void set_cutoff_lr(Real cut_in) { cutoff_lr_ = cut_in; }

	// score
	Real func( Size seqsep, Real d, Real xd, Real xh ) const;

	// derivative wrt dist,angle,angle
	Vector dfunc( Size seqsep, Real d, Real xd, Real xh  ) const;

	// cutoff
	Real cutoff( Size seqsep ) const { return (seqsep<=4 ? cutoff_sr_:cutoff_lr_); }

private:
	// helper function computes fast approximation to besseli_0
	inline Real BesselI0( Real X ) {
		Real Y,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9;
		P1=1.0; P2=3.5156229; P3=3.0899424; P4=1.2067492;
		P5=0.2659732; P6=0.360768e-1; P7=0.45813e-2;
		Q1=0.39894228; Q2=0.1328592e-1; Q3=0.225319e-2;
		Q4=-0.157565e-2; Q5=0.916281e-2; Q6=-0.2057706e-1;
		Q7=0.2635537e-1; Q8=-0.1647633e-1; Q9=0.392377e-2;
		if (fabs(X) < 3.75) {
			Y=(X/3.75)*(X/3.75);
			return (P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))));
		} else {
			Real AX=fabs(X);
			Y=3.75/AX;
			Real BX=exp(AX)/sqrt(AX);
			AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))));
			return (AX*BX);
		}
	}

private:
	utility::vector1< Real > lr_As_, sr_As_; 
	utility::vector1< numeric::xyzVector< Real > > lr_mus_, lr_sigmas_;
	utility::vector1< numeric::xyzVector< Real > > sr_mus_, sr_sigmas_;
	Real cutoff_sr_, cutoff_lr_;
};


} // ns scoring
} // ns core

#endif
