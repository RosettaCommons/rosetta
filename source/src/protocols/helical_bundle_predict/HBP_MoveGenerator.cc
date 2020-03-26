// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/helical_bundle_predict/HBP_MoveGenerator.cc
/// @brief A base class for a module to generate ParsedProtocols for the next move in a Monte Carlo trajectory, based on the current state of the pose.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#include <protocols/helical_bundle_predict/HBP_MoveGenerator.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.helical_bundle_predict.HBP_MoveGenerator" );


namespace protocols {
namespace helical_bundle_predict {

/// @brief Constructor.
HBP_MoveGenerator::HBP_MoveGenerator():
	utility::VirtualBase(),
	current_round_(0),
	max_rounds_(0)
{}

/// @brief Copy constructor.
HBP_MoveGenerator::HBP_MoveGenerator( HBP_MoveGenerator const &src ) :
	VirtualBase( src ),
	current_round_( src.current_round_ ),
	max_rounds_( src.max_rounds_ )
{}

/// @brief Destructor.
HBP_MoveGenerator::~HBP_MoveGenerator() {}

//////////////////////////////////////////// PUBLIC METHODS ////////////////////////////////////////////

/// @brief Does nothing in base class.  Can be overridden in derived classes to do something
/// more useful when a move is accepted.
void
HBP_MoveGenerator::mark_move_accepted() const { /*GNDN*/ }

/// @brief Does nothing in base class.  Can be overridden in derived classes to do something
/// more useful when a move is rejected.
void
HBP_MoveGenerator::mark_move_rejected() const { /*GNDN*/ }

/// @brief Set the current round.
void
HBP_MoveGenerator::set_current_round(
	core::Size const setting
) {
	if ( max_rounds_ > 0 ) runtime_assert( setting <= max_rounds_ );
	runtime_assert( setting > 0 );
	current_round_ = setting;
}

/// @brief Set the total number of rounds.
void
HBP_MoveGenerator::set_max_rounds(
	core::Size const setting
) {
	runtime_assert( setting > 0 );
	runtime_assert( current_round_ <= setting );
	max_rounds_ = setting;
}


} //protocols
} //helical_bundle_predict






