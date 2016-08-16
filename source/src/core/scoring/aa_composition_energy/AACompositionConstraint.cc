// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/scoring/aa_composition_energy/AACompositionConstraint.cc
/// @brief A constraint for constraining sequences to have a desired amino acid composition, analogous to a geometric constraint.
/// @details The corresponding energy term for this constraint is the AACompositionEnergy (aa_composition in wts files).
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#include <core/scoring/aa_composition_energy/AACompositionConstraint.hh>
#include <core/scoring/constraints/ConstraintIO.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/conformation/Conformation.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/scoring/aa_composition_energy/AACompositionEnergySetup.hh>
#include <basic/Tracer.hh>

#include <numeric/xyz.functions.hh>
#include <numeric/trig.functions.hh>
#include <numeric/deriv/dihedral_deriv.hh>

#include <utility/exit.hh>

#include <utility/vector1.hh>

namespace core {
namespace scoring {
namespace aa_composition_energy {

static THREAD_LOCAL basic::Tracer TR( "core.scoring.aa_composition_energy.AACompositionConstraint" );

/// @brief Constructor
///
AACompositionConstraint::AACompositionConstraint():
	core::scoring::aa_composition_energy::SequenceConstraint( core::scoring::aa_composition ),
	//TODO -- initialize variables here.
	selector_(),
	aa_comp_setup_( new AACompositionEnergySetup )
{}

/// @brief Copy constructor
///
AACompositionConstraint::AACompositionConstraint( AACompositionConstraint const &src ):
	core::scoring::aa_composition_energy::SequenceConstraint( core::scoring::aa_composition ),
	//TODO -- copy variables here.
	selector_(), //Cloned if present, below
	aa_comp_setup_() //Cloned below.
{
	if ( src.selector_ ) selector_ = src.selector_->clone();
	runtime_assert( src.aa_comp_setup_ );
	aa_comp_setup_ = src.aa_comp_setup_->clone();
}

/// @brief Destructor
///
AACompositionConstraint::~AACompositionConstraint() {}

/// @brief Clone operator
///
core::scoring::constraints::ConstraintOP
AACompositionConstraint::clone() const { return core::scoring::constraints::ConstraintOP( new AACompositionConstraint( *this ) ); }

bool AACompositionConstraint::operator == ( Constraint const & other ) const
{
	if ( ! other.same_type_as_me( *this ) ) return false;
	if ( !       same_type_as_me( other ) ) return false;

	// TO DO: implement ResidueSelector comparison operators.
	// TO DO: implement AACompositionEnergySetup comparison operator
	return false;
}

bool
AACompositionConstraint::same_type_as_me( Constraint const & other ) const
{
	return dynamic_cast< AACompositionConstraint const * > (&other);
}


/// @brief Set the selector to be used by this constraint.
/// @details Clones the input.
void
AACompositionConstraint::set_selector( core::select::residue_selector::ResidueSelectorCOP selector_in ) {
	selector_ = selector_in->clone();
	return;
}

select::residue_selector::ResidueSelectorCOP
AACompositionConstraint::selector() const {
	return selector_;
}

AACompositionEnergySetupCOP
AACompositionConstraint::aa_composition_energy_setup() const
{ return aa_comp_setup_; }

/// @brief Initialize the AACompositionEnergySetup object from a file.
///
void
AACompositionConstraint::initialize_from_file( std::string const &filename ) {
	runtime_assert( aa_comp_setup_ ); //The pointer should point at an actual object.
	aa_comp_setup_->initialize_from_file( filename );
	return;
}

/// @brief Initialize the AACompositionEnergySetup object from the contents of a file.
/// @details Allows external code to initialize a constriant object without having the
/// object read directly from disk.
void
AACompositionConstraint::initialize_from_file_contents( std::string const &filecontents ) {
	runtime_assert( aa_comp_setup_ ); //The pointer should point at an actual object.
	aa_comp_setup_->initialize_from_file_contents( filecontents );
	return;
}


} // aa_composition_energy
} // scoring
} // core
