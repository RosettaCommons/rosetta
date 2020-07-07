// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/constraints/ResidueTypeLinkingConstraint.hh
///
/// @brief Implements linked constraints, where bonuses or penalties can be applied if pairs of residues match specific sequences.
/// @author Sarel Fleishman
/// @author Brahm Yachnin (byachnin@visterrainc.com)
/// @author Andrew Wollacott (awollacott@visterrainc.com)


#ifndef INCLUDED_core_scoring_constraints_ResidueTypeLinkingConstraint_hh
#define INCLUDED_core_scoring_constraints_ResidueTypeLinkingConstraint_hh

#include <core/scoring/constraints/ResidueTypeLinkingConstraint.fwd.hh>
#include <core/scoring/constraints/Constraint.hh>

#include <core/conformation/Conformation.fwd.hh>
#include <core/id/AtomID.hh>
#include <core/chemical/AA.hh>

#include <utility/vector1.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace constraints {

/// @brief This class favors a particular residue identity at a particular position by reducing its res_type energy.
class ResidueTypeLinkingConstraint : public Constraint
{
public:
	typedef core::Real Real;
public:
	ResidueTypeLinkingConstraint();

	/// @brief Constructor that sets up ResidueTypeLinkingConstraints, where all bonus for residue pairs are all defined in a std::map.
	/// @details A negative number results in a bonus (lower score) and a positive number results in a penalty (higher
	/// @details score) if the constraint requirements are fulfilled.
	/// @details If the residue identities of the residues in the pose are not found, a default_bonus_ (defaults to 0) is scored.
	/// @details Other ctors should construct the bonus_map_ from residue identities.
	/// @author Brahm Yachnin (byachnin@visterrainc.com)
	ResidueTypeLinkingConstraint(
		Size seqpos1,
		Size seqpos2,
		std::map<std::pair<core::chemical::AA, core::chemical::AA>, core::Real> const & bonus_map
	);

	/// @brief Apply a bonus whenever seqpos1 and seqpos2 are NOT the same residue.
	/// @details bonus_map_ and default_bonus_ are constructed such that bonus is scored any time seqpos1 and seqpos2 aren't the same residue.
	ResidueTypeLinkingConstraint(
		Size seqpos1,
		Size seqpos2,
		Real bonus
	);


	/// @brief Apply a bonus when seqpos1 == AA1_name AND seqpos2 == AA2_name.
	/// @details We score the default_bonus_ (defaults to 0) whenever this isn't true.
	/// @details bonus_map_ is constructed such that bonus is scored any time seqpos1 == AA1_name AND seqpos2 == AA2_name.
	ResidueTypeLinkingConstraint(
		Size seqpos1,
		Size seqpos2,
		core::chemical::AA const & AA1,
		core::chemical::AA const & AA2,
		Real bonus
	);

	~ResidueTypeLinkingConstraint() override;

	Size
	natoms() const override { return 0; }

	AtomID const &
	atom( Size const ) const override {
		utility_exit_with_message("ResidueTypeLinkingConstraint is not atom-based!.");
		return core::id::GLOBAL_BOGUS_ATOM_ID;  // required for compilation on Windows
	}

	utility::vector1< core::Size >
	residues() const override;

	void
	show( std::ostream & out ) const override;

	/// @brief possibility to compare constraint according to data
	/// and not just pointers
	bool operator == ( Constraint const & other ) const override;
	bool same_type_as_me( Constraint const & ) const override;

	// Needed to get the base class overloads
	using Constraint::score;
	using Constraint::dist;

	void
	score( func::XYZ_Func const & xyz_func, EnergyMap const & weights, EnergyMap & emap ) const override;

	/// @details Return 1.0 if constraint will get a bonus/penalty, 0.0 if not
	core::Real
	dist( core::scoring::func::XYZ_Func const & xyz ) const override;

	void
	fill_f1_f2(
		AtomID const & atom,
		func::XYZ_Func const & xyz,
		Vector & F1,
		Vector & F2,
		EnergyMap const & weights
	) const override;

	ConstraintOP
	clone() const override;

	/// @brief Add a map entry to the current bonus_map_
	/// @details By default, if the entry key already exists, this will throw an error.
	/// @details Use a std::pair of core::chemical::AAs to specify the linked residue identities.
	void add_bonus_map_entry(std::pair<core::chemical::AA, core::chemical::AA> key, core::Real bonus, bool overwrite=false);

	/// @brief Add a map entry to the current bonus_map_
	/// @details By default, if the entry key already exists, this will throw an error.
	/// @details Use two core::chemical::AAs to specify the linked residue identities.
	void add_bonus_map_entry(core::chemical::AA aa1, core::chemical::AA aa2, core::Real bonus, bool overwrite=false);

	/// @brief Add a map entry to the current bonus_map_
	/// @details By default, if the entry key already exists, this will throw an error.
	/// @details Use a two-letter std::string to specify the linked residue identities.
	void add_bonus_map_entry(std::string key, core::Real bonus, bool overwrite=false);

	// Getters
	core::Size get_seqpos1() const;
	core::Size get_seqpos2() const;
	std::map<std::pair<core::chemical::AA, core::chemical::AA>, core::Real> const & get_bonus_map();
	core::Real get_default_bonus() const;

	// Setters
	void set_seqpos1(core::Size seqpos1 );
	void set_seqpos2(core::Size seqpos2 );
	void set_bonus_map(std::map<std::pair<core::chemical::AA, core::chemical::AA>, core::Real> const & bonus_map );
	void set_default_bonus( core::Real );

private:
	Size seqpos1_;
	Size seqpos2_;
	std::map<std::pair<core::chemical::AA, core::chemical::AA>, core::Real> bonus_map_;
	core::Real default_bonus_ = 0;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; // RotamerConstraint


} // namespace constraints
} // namespace scoring
} // namespace core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_constraints_ResidueTypeLinkingConstraint )
#endif // SERIALIZATION


#endif // INCLUDED_core_scoring_constraints_ResidueTypeLinkingConstraint_HH
