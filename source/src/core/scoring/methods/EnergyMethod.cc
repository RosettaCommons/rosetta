// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/method/EnergyMethod.cc
/// @brief  Base class for energy classes
/// @author Phil Bradley


// Unit headers
#include <core/scoring/methods/EnergyMethod.hh>

// Package headers
#include <core/scoring/methods/EnergyMethodCreator.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {

EnergyMethod::EnergyMethod( EnergyMethodCreatorOP creator ) :
	score_types_( creator->score_types_for_method() )
{
}

EnergyMethod::EnergyMethod( EnergyMethod const & src ) :
	parent(),
	score_types_( src.score_types_ )
{
}


EnergyMethod::~EnergyMethod() {}


bool
EnergyMethod::minimize_in_whole_structure_context( pose::Pose const & ) const
{
	/// APL -- Grandfathering mechanism.  As energy methods are updated to use
	/// the new derivative evaluation machinery, they can set this function to
	/// return "false", and therefore be used in rtmin.
	return true;
}

bool
EnergyMethod::defines_high_order_terms( pose::Pose const & ) const
{
	return false;
}

/// @details default implementation noop
void
EnergyMethod::setup_for_packing( pose::Pose &, utility::vector1< bool > const &, utility::vector1< bool > const & ) const {}


/// @details default implementation noop
void
EnergyMethod::prepare_rotamers_for_packing(
	pose::Pose const &,
	conformation::RotamerSetBase & ) const
{}


/// @details default implementation noop
void
EnergyMethod::update_residue_for_packing(
	pose::Pose &,
	Size ) const
{}

void
EnergyMethod::setup_for_scoring( pose::Pose &, ScoreFunction const & ) const {}

/// @details default implementation noop
void
EnergyMethod::setup_for_minimizing(
	pose::Pose & ,
	ScoreFunction const & ,
	kinematics::MinimizerMapBase const &
) const {}

/// @details default implementation noop
void
EnergyMethod::setup_for_derivatives( pose::Pose &, ScoreFunction const &  ) const {}

/// @details default implementation noop
void
EnergyMethod::finalize_after_derivatives( pose::Pose &, ScoreFunction const &  ) const
{}

/// @details default implementation does not alter either the F1 or F2 vectors.
void
EnergyMethod::eval_atom_derivative(
	id::AtomID const &,
	pose::Pose const &,
	kinematics::DomainMap const &,
	ScoreFunction const &,
	EnergyMap const &,
	Vector &,// F1,
	Vector & // F2
) const {}

/*
/// @details default behavior is to return 0
Real
EnergyMethod::eval_dof_derivative(
	id::DOF_ID const &,
	id::TorsionID const &,
	pose::Pose const &,
	ScoreFunction const &,
	EnergyMap const &
) const
{
	return 0.0;
}*/

/// called at the end of energy evaluation
void
EnergyMethod::finalize_total_energy(
	pose::Pose &,
	ScoreFunction const &,
	EnergyMap &
) const {}

/// @details Enforce that the derived class which has been
/// instantiated does not attempt to use an inactive score type.
/// If an inactive score type should be used, the ScoreType enumeration
/// in ScoreType.hh must be modified so that the desired type appears
/// before the n_score_types position in the enumeration.  Rosetta
/// must then be recompiled.

/* DEPRICATED!  ALL SCORE TYPE REGISTRATION SHOULD BE PERFORMED BY AN EnergyMethodCreator
void
EnergyMethod::add_score_type( ScoreType const & type )
{
	if ( type > n_score_types ) {
		std::cerr << "Critical error in EnergyMethod::add_score_type().\nAn EnergyMethod class has been instantiated ";
		std::cerr << "that will compute the inactive score type '" << type;
		std::cerr << "' defined at position " << (int) type << " in the ScoreType enumeration.\n";
		std::cerr << "This energy method would produce a segmentation fault down the road as it tried to write to ";
		std::cerr << "a position in an EnergyMap object beyond the end of the array contained in that object.\n";
		std::cerr << "Active score types must appear before the n_score_types element ";
		std::cerr << "(at position " << (int) n_score_types << ") as this element marks the end of the active score types.\n";
		std::cerr << "Rosetta must be recompiled after src/core/scoring/ScoreType.hh is modified to include " << type;
		std::cerr << " as an active score type." << std::endl;
		utility_exit_with_message( "ERROR: Attempted to use an inactive score type" );
	}
	score_types_.push_back( type );
}
*/

/// @brief Override the entirety of the score types list if they
/// were initialized incorrectly in a parent's constructor.
void
EnergyMethod::set_score_types( EnergyMethodCreatorOP creator ) {
	score_types_ = creator->score_types_for_method();
}


} // methods
} // scoring
} // core
