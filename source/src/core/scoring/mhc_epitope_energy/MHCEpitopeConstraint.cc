// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/scoring/mhc_epitope_energy/MHCEpitopeConstraint.cc
/// @brief A constraint for MHC epitope scores
/// Follows analogous file for Vikram K. Mulligan's NetChargeEnergy
/// @author Chris Bailey-Kellogg, cbk@cs.dartmouth.edu; Brahm Yachnin, brahm.yachnin@rutgers.edu

#include <core/scoring/mhc_epitope_energy/MHCEpitopeConstraint.hh>
#include <core/scoring/constraints/ConstraintIO.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/conformation/Conformation.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/ResidueRanges.hh>
#include <core/scoring/mhc_epitope_energy/MHCEpitopeEnergySetup.hh>
#include <basic/Tracer.hh>

#include <numeric/xyz.functions.hh>
#include <numeric/trig.functions.hh>
#include <numeric/deriv/dihedral_deriv.hh>

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

#include <utility/exit.hh>

#include <utility/vector1.hh>
#include <utility/pointer/memory.hh>

namespace core {
namespace scoring {
namespace mhc_epitope_energy {

static basic::Tracer TR( "core.scoring.mhc_epitope_energy.MHCEpitopeConstraint" );

/// @brief Constructor
///
MHCEpitopeConstraint::MHCEpitopeConstraint():
	core::scoring::aa_composition_energy::SequenceConstraint( core::scoring::mhc_epitope ),
	selector_(),
	mhc_epitope_setup_( utility::pointer::make_shared<MHCEpitopeEnergySetup>() ),
	cst_weight_(1.0)
{}

/// @brief Copy constructor
///
MHCEpitopeConstraint::MHCEpitopeConstraint( MHCEpitopeConstraint const &src ):
	core::scoring::aa_composition_energy::SequenceConstraint( core::scoring::mhc_epitope ),
	selector_(), //Cloned if present, below
	mhc_epitope_setup_(), //Cloned below.
	cst_weight_()
{
	if ( src.selector_ ) selector_ = src.selector_->clone();
	runtime_assert( src.mhc_epitope_setup_ );
	mhc_epitope_setup_ = src.mhc_epitope_setup_->clone();
	cst_weight_ = src.get_cst_weight();
}

/// @brief Destructor
///
MHCEpitopeConstraint::~MHCEpitopeConstraint() = default;

/// @brief Clone operator
///
core::scoring::constraints::ConstraintOP
MHCEpitopeConstraint::clone() const { return utility::pointer::make_shared< MHCEpitopeConstraint >( *this ); }

bool MHCEpitopeConstraint::operator == ( Constraint const & other ) const
{
	if ( ! other.same_type_as_me( *this ) ) return false;
	if ( !       same_type_as_me( other ) ) return false;

	// TODO: implement ResidueSelector comparison operators.

	MHCEpitopeConstraint const *o = dynamic_cast< MHCEpitopeConstraint const * > (&other);
	if ( o->cst_weight_ != cst_weight_  ) return false;
	if ( o->cst_selector_name_ != cst_selector_name_  ) return false; // In the absence of a real ResidueSelector comparison, check that the name of the selector is the same.
	if ( o->selector_.get() != selector_.get()  ) return false; // Also check if the location in memory for both selectors is the same.
	if ( ! (*o->mhc_epitope_setup_ == *mhc_epitope_setup_) ) return false;
	return true;
}

bool
MHCEpitopeConstraint::same_type_as_me( Constraint const & other ) const
{
	return dynamic_cast< MHCEpitopeConstraint const * > (&other);
}


/// @brief Set the selector to be used by this constraint.
/// @details Clones the input.
void
MHCEpitopeConstraint::set_selector( core::select::residue_selector::ResidueSelectorCOP selector_in ) {
	selector_ = selector_in->clone();
	return;
}

select::residue_selector::ResidueSelectorCOP
MHCEpitopeConstraint::selector() const {
	return selector_;
}

MHCEpitopeEnergySetupCOP
MHCEpitopeConstraint::mhc_epitope_energy_setup() const
{ return mhc_epitope_setup_; }

/// @brief Initialize the MHCEpitopeEnergySetup object from a file.
///
void
MHCEpitopeConstraint::initialize_from_file( std::string const &filename ) {
	runtime_assert( mhc_epitope_setup_ ); //The pointer should point at an actual object.
	mhc_epitope_setup_->initialize_from_file( filename );
	return;
}

/// @brief Initialize the MHCEpitopeEnergySetup object from the contents of a file.
/// @details Allows external code to initialize a constriant object without having the
/// object read directly from disk.
void
MHCEpitopeConstraint::initialize_from_file_contents( std::string const &filecontents ) {
	runtime_assert( mhc_epitope_setup_ ); //The pointer should point at an actual object.
	mhc_epitope_setup_->initialize_from_file_contents( filecontents );
	return;
}

/// @brief Set the cst_weight_
void
MHCEpitopeConstraint::set_cst_weight( core::Real cst_weight ) {
	cst_weight_ = cst_weight;
	return;
}

/// @brief Get the cst_weight_
core::Real
MHCEpitopeConstraint::get_cst_weight() const {return cst_weight_;}

/// @brief Print info on the constraint
void
MHCEpitopeConstraint::show_def (std::ostream &TO, pose::Pose const &pose) const {
	runtime_assert( mhc_epitope_setup_ );
	select::residue_selector::ResidueRangesOP ranges( utility::pointer::make_shared<select::residue_selector::ResidueRanges>() );
	ranges->from_subset( selector_->apply( pose ) );
	for ( auto const & range : *ranges ) {
		TO << "MHCEpitope Residue " << range.start() << " Residue " << range.stop() << std::endl;
		TO << mhc_epitope_setup_->report() << std::endl;
	}
}

} // mhc_epitope_energy
} // scoring
} // core

#ifdef    SERIALIZATION
template< class Archive >
void
core::scoring::mhc_epitope_energy::MHCEpitopeConstraint::save( Archive & arc ) const {
	arc( cereal::base_class< core::scoring::aa_composition_energy::SequenceConstraint >( this ) );
	arc( CEREAL_NVP( selector_ ) );
	arc( CEREAL_NVP( mhc_epitope_setup_ ) );
	arc( CEREAL_NVP( cst_weight_ ) );
	arc( CEREAL_NVP( cst_selector_name_ ) );
}

template< class Archive >
void
core::scoring::mhc_epitope_energy::MHCEpitopeConstraint::load( Archive & arc ) {
	arc( cereal::base_class< core::scoring::aa_composition_energy::SequenceConstraint >( this ) );
	arc( selector_ );
	arc( mhc_epitope_setup_ );
	arc( cst_weight_ );
	arc( cst_selector_name_ );
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::mhc_epitope_energy::MHCEpitopeConstraint );
CEREAL_REGISTER_TYPE( core::scoring::mhc_epitope_energy::MHCEpitopeConstraint )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_mhc_epitope_energy_MHCEpitopeConstraint )
#endif // SERIALIZATION
