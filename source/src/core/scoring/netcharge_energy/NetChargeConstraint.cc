// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/scoring/netcharge_energy/NetChargeConstraint.cc
/// @brief A constraint for constraining sequences to have a desired net charge, analogous to a geometric constraint.
/// @details The corresponding energy term for this constraint is the NetChargeEnergy (netcharge in wts files).
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#include <core/scoring/netcharge_energy/NetChargeConstraint.hh>
#include <core/scoring/constraints/ConstraintIO.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/conformation/Conformation.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/ResidueRanges.hh>
#include <core/scoring/netcharge_energy/NetChargeEnergySetup.hh>
#include <basic/Tracer.hh>

#include <numeric/xyz.functions.hh>
#include <numeric/trig.functions.hh>
#include <numeric/deriv/dihedral_deriv.hh>

#include <utility/exit.hh>

#include <utility/vector1.hh>

namespace core {
namespace scoring {
namespace netcharge_energy {

static basic::Tracer TR( "core.scoring.netcharge_energy.NetChargeConstraint" );

/// @brief Constructor
///
NetChargeConstraint::NetChargeConstraint():
	core::scoring::aa_composition_energy::SequenceConstraint( core::scoring::netcharge ),
	//TODO -- initialize variables here.
	selector_(),
	netcharge_setup_( new NetChargeEnergySetup )
{}

/// @brief Copy constructor
///
NetChargeConstraint::NetChargeConstraint( NetChargeConstraint const &src ):
	core::scoring::aa_composition_energy::SequenceConstraint( core::scoring::netcharge ),
	//TODO -- copy variables here.
	selector_(), //Cloned if present, below
	netcharge_setup_() //Cloned below.
{
	if ( src.selector_ ) selector_ = src.selector_->clone();
	runtime_assert( src.netcharge_setup_ );
	netcharge_setup_ = src.netcharge_setup_->clone();
}

/// @brief Destructor
///
NetChargeConstraint::~NetChargeConstraint() = default;

/// @brief Clone operator
///
core::scoring::constraints::ConstraintOP
NetChargeConstraint::clone() const { return core::scoring::constraints::ConstraintOP( new NetChargeConstraint( *this ) ); }

bool NetChargeConstraint::operator == ( Constraint const & other ) const
{
	if ( ! other.same_type_as_me( *this ) ) return false;
	if ( !       same_type_as_me( other ) ) return false;

	// TO DO: implement ResidueSelector comparison operators.
	// TO DO: implement NetChargeEnergySetup comparison operator
	return false;
}

bool
NetChargeConstraint::same_type_as_me( Constraint const & other ) const
{
	return dynamic_cast< NetChargeConstraint const * > (&other);
}


/// @brief Set the selector to be used by this constraint.
/// @details Clones the input.
void
NetChargeConstraint::set_selector( core::select::residue_selector::ResidueSelectorCOP selector_in ) {
	selector_ = selector_in->clone();
	return;
}

select::residue_selector::ResidueSelectorCOP
NetChargeConstraint::selector() const {
	return selector_;
}

NetChargeEnergySetupCOP
NetChargeConstraint::netcharge_energy_setup() const
{ return netcharge_setup_; }

/// @brief Initialize the NetChargeEnergySetup object from a file.
///
void
NetChargeConstraint::initialize_from_file( std::string const &filename ) {
	runtime_assert( netcharge_setup_ ); //The pointer should point at an actual object.
	netcharge_setup_->initialize_from_file( filename );
	return;
}

/// @brief Initialize the NetChargeEnergySetup object from the contents of a file.
/// @details Allows external code to initialize a constriant object without having the
/// object read directly from disk.
void
NetChargeConstraint::initialize_from_file_contents( std::string const &filecontents ) {
	runtime_assert( netcharge_setup_ ); //The pointer should point at an actual object.
	netcharge_setup_->initialize_from_file_contents( filecontents );
	return;
}

/// @brief Print info on the constraint
void
NetChargeConstraint::show_def (std::ostream &TO, pose::Pose const &pose) const {
	runtime_assert( netcharge_setup_ );
	select::residue_selector::ResidueRangesOP ranges( new select::residue_selector::ResidueRanges );
	ranges->from_subset( selector_->apply( pose ) );
	for ( auto const & range : *ranges ) {
		TO << "NetCharge Residue " << range.start() << " Residue " << range.stop() << std::endl;
		TO << netcharge_setup_->report() << std::endl;
	}
}

} // netcharge_energy
} // scoring
} // core
