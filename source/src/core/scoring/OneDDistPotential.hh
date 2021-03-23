// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/OneDDistPotential.hh
/// @brief
/// @detailed
/// @author Andy Watkins (amw579@stanford.edu)


#ifndef INCLUDED_core_scoring_OneDDistPotential_HH
#define INCLUDED_core_scoring_OneDDistPotential_HH

#include <utility/VirtualBase.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/OneDDistPotential.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.hh>
#include <utility/json_spirit/json_spirit_value.h>
#include <numeric/MathNTensor.hh>
#include <numeric/interpolation/polycubic_catmull_rom.hh>

namespace core {
namespace scoring {

class OneDDistPotential: public utility::VirtualBase {

public:

	//constructor
	OneDDistPotential( std::string const & filename );

	OneDDistPotential( numeric::MathNTensor< core::Real, 1 > const & tensor,
		utility::json_spirit::mObject const & json );

	//destructor
	~OneDDistPotential() override;

public:

	core::Real
	evaluate( core::Real const dist ) const;

	core::Real
	get_derivative( core::Real const dist ) const;

	core::Real
	evaluate( core::Real const dist,
		bool const compute_deriv, core::Real & deriv ) const;

private:

	void
	initialize_from_json( utility::json_spirit::mObject const & json );

	// core::Real
	// evaluate_constraining_potential( Distance const & R ) const;

	// std::pair< Vector, Vector >
	// evaluate_constraining_potential_derivative( Vector const & t ) const;

private:

	numeric::MathNTensor< core::Real, 1 > tensor_;
	utility::fixedsizearray1< core::Real, 1 > minval_, binwidth_;
	// bounds past which a super-steep penalty is applied.
	// core::Real Rmax_, Emax_;

	utility::fixedsizearray1< numeric::interpolation::CatmullRomSplineBoundaryType, 1 > const boundary_;

};

} //scoring
} //core

#endif
