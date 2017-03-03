// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/loop_graph/evaluator/SixDTransRotPotential.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_core_scoring_loop_graph_SixDTransRotPotential_HH
#define INCLUDED_core_scoring_loop_graph_SixDTransRotPotential_HH

#include <utility/pointer/ReferenceCount.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/loop_graph/evaluator/SixDTransRotPotential.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.hh>
#include <utility/json_spirit/json_spirit_value.h>
#include <numeric/MathNTensor.hh>
#include <numeric/interpolation/polycubic_catmull_rom.hh>

namespace core {
namespace scoring {
namespace loop_graph {
namespace evaluator {

class SixDTransRotPotential: public utility::pointer::ReferenceCount {

public:

	//constructor
	SixDTransRotPotential( std::string const & filename );

	SixDTransRotPotential( numeric::MathNTensor< core::Real, 6 > const & tensor,
		utility::json_spirit::mObject const & json );

	//destructor
	~SixDTransRotPotential();

public:

	core::Real
	evaluate( core::kinematics::Jump const & j ) const;

	std::pair< Vector, Vector >
	get_derivative( core::kinematics::Jump const & j ) const;

	core::Real
	evaluate( Vector const & t, Vector const & rotation_vector,
		bool const compute_deriv, std::pair< Vector, Vector > & deriv ) const;

	void set_enforce_continuity_at_pi( bool const & setting ){ enforce_continuity_at_pi_ = setting; }
	bool enforce_continuity_at_pi() const { return enforce_continuity_at_pi_; }

private:

	void
	initialize_from_json( utility::json_spirit::mObject const & json );

	void
	adjust_near_pi( Vector const & rotation_vector,
		core::Real & value,
		utility::fixedsizearray1< core::Real, 6 > & outvals,
		bool const & compute_deriv,
		utility::fixedsizearray1< core::Real, 6 > & tensor_deriv ) const;

	core::Real
	evaluate_constraining_potential( Distance const & R ) const;

	std::pair< Vector, Vector >
	evaluate_constraining_potential_derivative( Vector const & t ) const;

private:

	numeric::MathNTensor< core::Real, 6 > tensor_;
	utility::fixedsizearray1< core::Real, 6 > minval_, binwidth_;
	// bounds past which a super-steep penalty is applied.
	core::Real Rmax_, Emax_;

	bool enforce_continuity_at_pi_;
	bool const turn_off_rotation_dependence_;
	bool const use_cubic_interp_;
	utility::fixedsizearray1< numeric::interpolation::CatmullRomSplineBoundaryType, 6 > const boundary_;

};

} //evaluator
} //loop_graph
} //scoring
} //core

#endif
