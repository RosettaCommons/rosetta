// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/helical_bundle_predict/HBP_FinalFullatomRefinementMoveGenerator.cc
/// @brief A class to generate ParsedProtocols for the final full-atom refinement step.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

// Project includes:
#include <protocols/helical_bundle_predict/HBP_FinalFullatomRefinementMoveGenerator.hh>

// Core includes:
#include <core/scoring/ScoreFunction.hh>

// Protocols includes:
#include <protocols/relax/FastRelax.hh>
#include <protocols/rosetta_scripts/ParsedProtocol.hh>
#include <protocols/cyclic_peptide/TryDisulfPermutations.hh>

// Basic includes:
#include <basic/Tracer.hh>

// Utility includes:
#include <utility/pointer/memory.hh>

static basic::Tracer TR( "protocols.helical_bundle_predict.HBP_FinalFullatomRefinementMoveGenerator" );


namespace protocols {
namespace helical_bundle_predict {

/// @brief Options constructor.
HBP_FinalFullatomRefinementMoveGenerator::HBP_FinalFullatomRefinementMoveGenerator(
	core::Size const fast_relax_rounds,
	bool const find_disulfides,
	core::scoring::ScoreFunctionCOP sfxn_in
):
	HBP_MoveGenerator(),
	fast_relax_rounds_(fast_relax_rounds < 1 ? 1 : fast_relax_rounds),
	find_disulfides_(find_disulfides),
	sfxn_() //Initialized below
{
	runtime_assert_string_msg( sfxn_in != nullptr, "Error in HBP_FinalFullatomRefinementMoveGenerator constructor: The scoring function provided must not be null!" );
	sfxn_ = sfxn_in->clone();
}

/// @brief Copy constructor.
HBP_FinalFullatomRefinementMoveGenerator::HBP_FinalFullatomRefinementMoveGenerator( HBP_FinalFullatomRefinementMoveGenerator const &src ) :
	HBP_MoveGenerator(src),
	fast_relax_rounds_(src.fast_relax_rounds_),
	find_disulfides_(src.find_disulfides_),
	sfxn_( src.sfxn_->clone() )
{
	runtime_assert( sfxn_ != nullptr );
	runtime_assert( fast_relax_rounds_ > 0 );
}

/// @brief Destructor.
HBP_FinalFullatomRefinementMoveGenerator::~HBP_FinalFullatomRefinementMoveGenerator() {}

/// @brief Clone this object: that is, make a copy and return an owning pointer to the copy.
HBP_MoveGeneratorOP
HBP_FinalFullatomRefinementMoveGenerator::clone() const {
	return utility::pointer::make_shared< HBP_FinalFullatomRefinementMoveGenerator >( *this );
}

//////////////////////////////////////////// PUBLIC METHODS ////////////////////////////////////////////

/// @brief Given the current step index, the total number of steps in the trajectory, and the pose for analysis,
/// construct a ParsedProtocol of things to do to this pose for this move in the Monte Carlo trajectory.
/// @details Pure virtual.  Must be implemented by derived classes.
protocols::rosetta_scripts::ParsedProtocolOP
HBP_FinalFullatomRefinementMoveGenerator::generate_monte_carlo_move(
	core::Size const ,//current_step,
	core::Size const ,//num_steps,
	core::pose::Pose const &//pose
) const {
	using namespace protocols::rosetta_scripts;
	using namespace protocols::relax;
	using namespace protocols::cyclic_peptide;

	ParsedProtocolOP pp( utility::pointer::make_shared< ParsedProtocol >() );

	if ( find_disulfides_ ) {
		TryDisulfPermutationsOP find_disulf( utility::pointer::make_shared< TryDisulfPermutations >() );
		pp->add_mover_filter_pair( find_disulf, "Find Disulfides", nullptr, false );
	}

	FastRelaxOP frlx( utility::pointer::make_shared< FastRelax >( sfxn_, fast_relax_rounds_ ) );
	pp->add_mover_filter_pair( frlx, "Fullatom FastRelax", nullptr, false );

	return pp;
}

} //protocols
} //helical_bundle_predict






