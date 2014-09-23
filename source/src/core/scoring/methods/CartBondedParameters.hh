// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/CartBondedParameters.hh
/// @brief
/// @author

#ifndef INCLUDED_core_scoring_methods_CartBondedParameters_hh
#define INCLUDED_core_scoring_methods_CartBondedParameters_hh

#include <core/types.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>
#include <numeric/angle.functions.hh>
#include <numeric/constants.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <numeric/interpolation/spline/Bicubic_spline.hh>
#include <numeric/MathMatrix.hh>

#include <ObjexxFCL/FArray2D.hh>

namespace core {
namespace scoring {
namespace methods {

// Unit headers
class CartBondedParameters : public utility::pointer::ReferenceCount {
public:
	virtual
	Real mu(Real phi, Real psi) const = 0;

	virtual
	Real K(Real phi, Real psi) const = 0;

	virtual
	Size period() const = 0;

	virtual
	Real dmu_dphi(Real /*phi*/, Real /*psi*/) const {
		return 0;
	}

	virtual
	Real dK_dphi(Real /*phi*/, Real /*psi*/) const {
		return 0;
	}

	virtual
	Real dmu_dpsi(Real /*phi*/, Real /*psi*/) const {
		return 0;
	}

	virtual
	Real dK_dpsi(Real /*phi*/, Real /*psi*/) const {
		return 0;
	}

	virtual
	bool is_null() const {
		return false;
	}
};


class BBIndepCartBondedParameters : public CartBondedParameters {
public:
	BBIndepCartBondedParameters() {
		K0_ = 0;  // K==0 is a null constraint
		mu0_ = 0;
		period_ = 0;
	}

	BBIndepCartBondedParameters(core::Real mu0_in, core::Real K0_in, core::Size period_in=1) {
		mu0_ = mu0_in;
		K0_ = K0_in;
		period_ = period_in;
	}

	virtual
	Real mu(Real , Real ) const {
		return mu0_;
	}

	virtual
	Real K(Real , Real ) const {
		return K0_;
	}

	virtual
	Size period() const {
		return period_;
	}

	virtual
	bool is_null() const {
		return (K0_==0);
	}

private:
	core::Real mu0_, K0_;
	core::Size period_;
};

class BBDepCartBondedParameters : public CartBondedParameters {
public:
	BBDepCartBondedParameters( ) {}

	BBDepCartBondedParameters(
			ObjexxFCL::FArray2D< core::Real > const &mu,
			ObjexxFCL::FArray2D< core::Real > const &Ks,
			std::string tag_in="" ) {
		init( mu, Ks );
		tag_ = tag_in;
	}

	// init
	void init ( ObjexxFCL::FArray2D< core::Real > const &mu, ObjexxFCL::FArray2D< core::Real > const &Ks) {
		using namespace numeric;
		using namespace numeric::interpolation::spline;
		MathMatrix< Real > mu_copy( 36, 36 ), Ks_copy( 36, 36 );
		for ( Size jj = 0; jj < 36; ++jj ) {
			for ( Size kk = 0; kk < 36; ++kk ) {
				mu_copy( jj, kk ) = mu(jj+1,kk+1);
				Ks_copy( jj, kk ) = Ks(jj+1,kk+1);
			}
		}
		init( mu_copy, mus_spline_);
		init( Ks_copy, Ks_spline_);
	}

	void init ( numeric::MathMatrix< Real > const &x, numeric::interpolation::spline::BicubicSpline &x_spline) {
		using namespace numeric;
		using namespace numeric::interpolation::spline;

		BorderFlag periodic_boundary[2] = { e_Periodic, e_Periodic };
		Real start_vals[2] = {5.0, 5.0};
		Real deltas[2] = {10.0, 10.0};
		bool lincont[2] = {false,false};
		std::pair< Real, Real > unused[2];
		unused[0] = std::make_pair( 0.0, 0.0 );
		unused[1] = std::make_pair( 0.0, 0.0 );
		x_spline.train( periodic_boundary, start_vals, deltas, x, lincont, unused );
	}

	// interpolate mu
	virtual
	Real mu (Real phi, Real psi) const { return mus_spline_.F(phi,psi); }

	// interpolate Ks
	virtual
	Real K (Real phi, Real psi) const {
		core::Real retval = Ks_spline_.F(phi,psi);
		//if (retval<0) {
		//	std::cerr << tag_ << "   K(" << phi << "," << psi << ") = " << Ks_spline_.F(phi,psi) << std::endl;
		//}
		return retval;
	}

	// get dmu_dphi
	virtual
	Real dmu_dphi (Real phi, Real psi) const { return mus_spline_.dFdx(phi,psi); }

	// get dmu_dpsi
	virtual
	Real dmu_dpsi (Real phi, Real psi) const  { return mus_spline_.dFdy(phi,psi); }

	// get dmu_dphi
	virtual
	Real dK_dphi (Real phi, Real psi) const {
		//std::cerr << "   dKs_dphi(" << phi << "," << psi << ") = " << Ks_spline_.dFdx(phi,psi) << std::endl;
		return Ks_spline_.dFdx(phi,psi);
	}

	// get dmu_dpsi
	virtual
	Real dK_dpsi (Real phi, Real psi) const  {
		//std::cerr << "   dKs_dpsi(" << phi << "," << psi << ") = " << Ks_spline_.dFdy(phi,psi) << std::endl;
		return Ks_spline_.dFdy(phi,psi);
	}

	Size period() const {
		return 1;
	}

private:
	numeric::interpolation::spline::BicubicSpline mus_spline_;
	numeric::interpolation::spline::BicubicSpline Ks_spline_;
	std::string tag_;
};

// owning pointers
typedef  utility::pointer::weak_ptr< CartBondedParameters > CartBondedParametersAP;
typedef  utility::pointer::weak_ptr< CartBondedParameters const > CartBondedParametersCAP;
typedef  utility::pointer::shared_ptr< CartBondedParameters > CartBondedParametersOP;
typedef  utility::pointer::shared_ptr< CartBondedParameters const > CartBondedParametersCOP;



} // namespace methods
} // namespace scoring
} // namespace core


#endif // INCLUDED_core_scoring_methods_CartesianBondedEnergy_HH
