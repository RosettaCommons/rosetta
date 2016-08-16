// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

//////////////////////////////////////////////
///
/// @file protocols/scoring/methods/pcs2/TensorsOptimizer.hh
///
/// @brief
///
/// @details
///
/// @param
///
/// @return
///
/// @remarks
///
/// @references
///
/// @authorv Christophe Schmitz
///
////////////////////////////////////////////////

#ifndef INCLUDED_protocols_scoring_methods_pcs2_TensorsOptimizer_hh
#define INCLUDED_protocols_scoring_methods_pcs2_TensorsOptimizer_hh

// Package headers
#include <protocols/scoring/methods/pcs2/PcsDataCenter.fwd.hh>
// Project headers
#include <core/optimization/Multifunc.hh>

// Utility headers

// Numeric headers
#include <numeric/constants.hh>

#include <utility/vector1.hh>


// Objexx headers

// C++ headers


namespace protocols {
namespace scoring {
namespace methods {
namespace pcs2 {

class TensorsOptimizer : public core::optimization::Multifunc {

public:
	//PcsDataCenter  const & pcs_d_c_;
	PcsDataCenter  & pcs_d_c_;

	//TensorsOptimizer();

	TensorsOptimizer(PcsDataCenter /*const*/  & pcs_d_c);

	virtual
	~TensorsOptimizer();

	// @brief OptE func
	virtual
	core::Real
	operator ()( core::optimization::Multivec const & vars ) const;

	core::Real
	func( core::optimization::Multivec const & vars ) const;


	/// @brief OptE dfunc
	virtual
	void
	dfunc(core::optimization::Multivec const & vars,
		core::optimization::Multivec & dE_dvars
	) const;

	virtual
	bool
	abort_min(core::optimization::Multivec const & vars ) const;

	/// @brief exact derivative (fast)
	void
	dfunc_exact(core::optimization::Multivec const & vars,
		core::optimization::Multivec & dE_dvars
	) const;

	/// @brief numeric derivative (slow)
	void
	dfunc_numeric(core::optimization::Multivec const & vars,
		core::optimization::Multivec & dE_dvars
	) const;

private:


};

}//namespace pcs2
}//namespace methods
}//namespace scoring
}//namespace protocols

static const core::Real FACT_USI_PRECALC_FOR_A_3( (10000.0/12.0/ core::Real( numeric::constants::d::pi ) ) * 3.0 );

static const core::Real FACT_20_PI_OVER_10000(20 * numeric::constants::d::pi / 10000.0);
static const core::Real FACT_10000_OVER_4PI(10000.0/(4 * numeric::constants::d::pi));

#endif
