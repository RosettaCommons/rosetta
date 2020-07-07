// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/constraints/ResidueTypeLinkingConstraint.cc
///
/// @brief Implements linked constraints, where bonuses or penalties can be applied if pairs of residues match specific sequences.
/// @author Sarel Fleishman
/// @author Brahm Yachnin (byachnin@visterrainc.com)
/// @author Andrew Wollacott (awollacott@visterrainc.com)


#include <core/scoring/constraints/ResidueTypeLinkingConstraint.hh>

#include <core/conformation/Residue.hh>
#include <core/scoring/ScoreType.hh>
#include <basic/Tracer.hh>

#include <core/kinematics/Jump.hh>
#include <core/chemical/AA.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/func/XYZ_Func.hh>
#include <utility>
#include <utility/vector1.hh>
#include <utility/pointer/memory.hh>


#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/utility.hpp> // For std::pair
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace constraints {


static basic::Tracer TR( "core.scoring.constraints.ResidueTypeLinkingConstraint" );

ResidueTypeLinkingConstraint::ResidueTypeLinkingConstraint():
	Constraint( core::scoring::res_type_linking_constraint )

{}

/// @brief Constructor that sets up ResidueTypeLinkingConstraints, where all bonus for residue pairs are all defined in a std::map.
/// @details A negative number results in a bonus (lower score) and a positive number results in a penalty (higher
/// @details score) if the constraint requirements are fulfilled.
/// @details If the residue identities of the residues in the pose are not found, a default_bonus_ (defaults to 0) is scored.
/// @details Other ctors should construct the bonus_map_ from residue identities.
/// @author Brahm Yachnin (byachnin@visterrainc.com)
ResidueTypeLinkingConstraint::ResidueTypeLinkingConstraint(
	Size seqpos1,
	Size seqpos2,
	std::map<std::pair<core::chemical::AA, core::chemical::AA>, core::Real> const & bonus_map
):
	Constraint( core::scoring::res_type_linking_constraint ),

	seqpos1_( seqpos1 ),
	seqpos2_( seqpos2 )
{
	set_bonus_map(bonus_map);
}

/// @brief Apply a bonus whenever seqpos1 and seqpos2 are NOT the same residue.
/// @details bonus_map_ and default_bonus_ are constructed such that bonus is scored any time seqpos1 and seqpos2 aren't the same residue.
ResidueTypeLinkingConstraint::ResidueTypeLinkingConstraint(
	Size seqpos1,
	Size seqpos2,
	core::Real bonus
):
	Constraint( core::scoring::res_type_linking_constraint ),

	seqpos1_( seqpos1 ),
	seqpos2_( seqpos2 )
{
	// Set the default bonus to bonus
	set_default_bonus(bonus);

	// Iterate over canonical AAs.
	for ( auto aa( static_cast<core::Size>( core::chemical::first_l_aa )); aa<=static_cast<core::Size>( core::chemical::num_canonical_aas ); ++aa ) {
		// Generate a bonus_map_ key for identical residue pairs (e.g. "AA", "CC", etc.)
		std::pair<core::chemical::AA, core::chemical::AA> key(static_cast<core::chemical::AA>(aa), static_cast<core::chemical::AA>(aa));
		// Modify the bonus_map_ to give no bonus/penalty for this key
		bonus_map_[key] = 0;
	}
}


/// @brief Apply a bonus when seqpos1 == AA1_name AND seqpos2 == AA2_name.
/// @details We score the default_bonus_ (defaults to 0) whenever this isn't true.
/// @details bonus_map_ is constructed such that bonus is scored any time seqpos1 == AA1_name AND seqpos2 == AA2_name.
ResidueTypeLinkingConstraint::ResidueTypeLinkingConstraint(
	Size seqpos1,
	Size seqpos2,
	core::chemical::AA const & AA1,
	core::chemical::AA const & AA2,
	core::Real bonus
):
	Constraint( core::scoring::res_type_linking_constraint ),

	seqpos1_( seqpos1 ),
	seqpos2_( seqpos2 )
{
	// Generate the bonus_map_ key
	std::pair<core::chemical::AA, core::chemical::AA> key(AA1, AA2);
	// Add the key to bonus_map_, with bonus as value
	bonus_map_[key] = bonus;
}

ResidueTypeLinkingConstraint::~ResidueTypeLinkingConstraint() = default;

/// @brief Add a map entry to the current bonus_map_
/// @details By default, if the entry key already exists, this will throw an error.
/// @details Use a std::pair of core::chemical::AAs to specify the linked residue identities.
void
ResidueTypeLinkingConstraint::add_bonus_map_entry(
	std::pair<core::chemical::AA, core::chemical::AA> key,
	core::Real bonus, bool overwrite /*= false*/)
{
	// If we specify overwrite, just add/replace key entry with bonus.
	if ( overwrite ) {
		bonus_map_[key] = bonus;
	} else { // If we don't specify overwrite, check if the key exists first adding key.
		if ( bonus_map_.count(key) == 1 ) {
			utility_exit_with_message("Trying to add a bonus value to an existing bonus_map_ key, and overwrite is false.");
		} else {
			bonus_map_[key] = bonus;
		}
	}
}

/// @brief Add a map entry to the current bonus_map_
/// @details By default, if the entry key already exists, this will throw an error.
/// @details Use two core::chemical::AAs to specify the linked residue identities.
void
ResidueTypeLinkingConstraint::add_bonus_map_entry(core::chemical::AA aa1, core::chemical::AA aa2, core::Real bonus, bool overwrite /*= false*/)
{
	std::pair<core::chemical::AA, core::chemical::AA> key(aa1, aa2);
	add_bonus_map_entry(key, bonus, overwrite);
}

/// @brief Add a map entry to the current bonus_map_
/// @details By default, if the entry key already exists, this will throw an error.
/// @details Use a two-letter std::string to specify the linked residue identities.
void
ResidueTypeLinkingConstraint::add_bonus_map_entry(std::string key, core::Real bonus, bool overwrite /*= false*/)
{
	if ( key.length() != 2 ) {
		utility_exit_with_message( "Bonus map must have two-residue keys (e.g. 'AK'), not " + key );
	}

	std::pair<core::chemical::AA, core::chemical::AA>
		key_pair(core::chemical::aa_from_oneletter_code(key[0]), core::chemical::aa_from_oneletter_code(key[1]));

	add_bonus_map_entry(key_pair, bonus, overwrite);
}

ConstraintOP
ResidueTypeLinkingConstraint::clone() const
{
	return ConstraintOP( utility::pointer::make_shared< ResidueTypeLinkingConstraint >( *this ) );
}

utility::vector1< core::Size >
ResidueTypeLinkingConstraint::residues() const {
	utility::vector1< core::Size > pos_list;
	pos_list.push_back(seqpos1_);
	pos_list.push_back(seqpos2_);
	return pos_list;
}

void
ResidueTypeLinkingConstraint::show( std::ostream & out ) const {
	out << "ResidueTypeLinkingConstraint; ";
	out << "seqpos1: " << seqpos1_ << "; ";
	out << "seqpos2: " << seqpos2_ << "; ";
	out << "default_bonus: " << default_bonus_ << "; ";
	if ( bonus_map_.size() == 400 ) {
		out << "bonus map contains 400 elements (all AA pairs represented).";
	} else {
		out << "Bonus map:" << std::endl;

		for ( auto it = bonus_map_.cbegin(); it != bonus_map_.cend(); ++it ) {
			out << it->first << ": " << it->second << std::endl;
		}
	}
}

bool
ResidueTypeLinkingConstraint::operator == ( Constraint const & other_cst ) const
{
	if ( !           same_type_as_me( other_cst ) ) return false;
	if ( ! other_cst.same_type_as_me( *this ) ) return false;

	auto const & other( static_cast< ResidueTypeLinkingConstraint const & > (other_cst) );

	if ( seqpos1_ != other.seqpos1_ ) return false;
	if ( seqpos2_ != other.seqpos2_ ) return false;
	if ( default_bonus_ != other.default_bonus_ ) return false;
	if ( bonus_map_ != other.bonus_map_ ) return false;
	if ( score_type() != other.score_type() ) return false;

	return true;
}

bool ResidueTypeLinkingConstraint::same_type_as_me( Constraint const & other ) const {
	return dynamic_cast< ResidueTypeLinkingConstraint const * > (&other);
}

/// @brief Calculates a score for this constraint.
/// @details Calculates using XYZ_Func, and puts the UNWEIGHTED score intoemap. Although the current set of weights currently
/// @details is provided, Constraint objects should put unweighted scores into emap.
/// @author Sarel Fleishman
/// @author Brahm Yachnin (byachnin@visterrainc.com)
/// @author Andrew Wollacott (awollacott@visterrainc.com)
void
ResidueTypeLinkingConstraint::score( func::XYZ_Func const & xyz_func, EnergyMap const & weights, EnergyMap & emap ) const
{
	Real const weight(weights[ this->score_type() ] );
	// std::cout << "res_type_linking_cst weight " << weight << std::endl;

	if ( weight == 0 ) return; // what's the point?

	conformation::Residue const & rsd1( xyz_func.residue(seqpos1_) );
	conformation::Residue const & rsd2( xyz_func.residue(seqpos2_) );

	// pair is a string with the current residue names at seqpos1_ and seqpos2_
	// This is also the key to loop up in the std::map bonus_map_
	std::pair<core::chemical::AA, core::chemical::AA> pair(rsd1.aa(), rsd2.aa());

	// Look up pair in the std::map
	std::map<std::pair<core::chemical::AA, core::chemical::AA>, core::Real>::const_iterator it = bonus_map_.find(pair);

	// If we've reached the end of bonus_map_, the AA pair is not in the map.  Use the default.
	if ( it == bonus_map_.end() ) {
		emap[ this->score_type() ] += default_bonus_;
	} else {
		// pair is in bonus_map_.  Adjust the score in the emap by second (the bonus).
		emap[ this->score_type() ] += it->second;
	}
}

core::Real
ResidueTypeLinkingConstraint::dist( core::scoring::func::XYZ_Func const & xyz_func ) const {
	return xyz_func.residue(seqpos1_).aa() != xyz_func.residue(seqpos2_).aa();
}

void
ResidueTypeLinkingConstraint::fill_f1_f2(
	AtomID const & ,//atom,
	func::XYZ_Func const &,
	Vector & ,//F1,
	Vector & ,//F2,
	EnergyMap const & //weights
) const
{
	// Do nothing.
	// Derivative of this restraint is effectively zero
	// so we just "add zero" to F1 and F2.
}

// Getters
core::Size ResidueTypeLinkingConstraint::get_seqpos1() const {return seqpos1_;}
core::Size ResidueTypeLinkingConstraint::get_seqpos2() const {return seqpos2_;}
std::map<std::pair<core::chemical::AA, core::chemical::AA>, core::Real> const & ResidueTypeLinkingConstraint::get_bonus_map() {return bonus_map_;}
core::Real ResidueTypeLinkingConstraint::get_default_bonus() const {return default_bonus_;}

// Setters
void ResidueTypeLinkingConstraint::set_seqpos1(core::Size seqpos1) {seqpos1_ = seqpos1;}
void ResidueTypeLinkingConstraint::set_seqpos2(core::Size seqpos2) {seqpos2_ = seqpos2;}
void ResidueTypeLinkingConstraint::set_default_bonus(core::Real default_bonus) {default_bonus_ = default_bonus;}

void ResidueTypeLinkingConstraint::set_bonus_map(std::map<std::pair<core::chemical::AA, core::chemical::AA>, core::Real> const & bonus_map)
{
	// bonus_map is valid.  Set bonus_map_ to bonus_map
	bonus_map_ = bonus_map;
}

} // namespace constraints
} // namespace scoring
} // namespace core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::constraints::ResidueTypeLinkingConstraint::save( Archive & arc ) const {
	arc( cereal::base_class< Constraint >( this ) );
	arc( CEREAL_NVP( seqpos1_ ) ); // Size
	arc( CEREAL_NVP( seqpos2_ ) ); // Size
	arc( CEREAL_NVP( bonus_map_ ) ); // std::map
	arc( CEREAL_NVP( default_bonus_ ) ); // core::Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::constraints::ResidueTypeLinkingConstraint::load( Archive & arc ) {
	arc( cereal::base_class< Constraint >( this ) );
	arc( seqpos1_ ); // Size
	arc( seqpos2_ ); // Size
	arc( bonus_map_ ); // std::map
	arc( default_bonus_ ); // core::Real
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::constraints::ResidueTypeLinkingConstraint );
CEREAL_REGISTER_TYPE( core::scoring::constraints::ResidueTypeLinkingConstraint )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_constraints_ResidueTypeLinkingConstraint )
#endif // SERIALIZATION
