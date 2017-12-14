// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/constraints/ResidueConstraint.hh
///
/// @brief implements favor_native_residue constraints, assigning an energy bonus
/// to a residue if it matches the identity within the constraint object
/// @author Sarel Fleishman


#ifndef INCLUDED_core_scoring_constraints_ResidueTypeLinkingConstraint_hh
#define INCLUDED_core_scoring_constraints_ResidueTypeLinkingConstraint_hh

#include <core/scoring/constraints/ResidueTypeLinkingConstraint.fwd.hh>
#include <core/scoring/constraints/Constraint.hh>

#include <core/conformation/Conformation.fwd.hh>
#include <core/id/AtomID.hh>

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

	ResidueTypeLinkingConstraint(
		pose::Pose const & pose,
		Size seqpos1,
		Size seqpos2,
		Real bonus
	);


	ResidueTypeLinkingConstraint(
		pose::Pose const & pose,
		Size seqpos1,
		Size seqpos2,
		std::string const & AA1name,
		std::string const & AA2name,
		Real bonus
	);
	/*
	ResidueTypeLinkingConstraint(
	Size seqpos1,
	Size seqpos2,
	std::string aa_in,
	std::string name3_in,
	core::Real bonus_in,
	utility::vector1< AtomID > const & atoms_in
	);
	*/
	virtual ~ResidueTypeLinkingConstraint();

	virtual
	Size
	natoms() const { return 0; }

	virtual
	AtomID const &
	atom( Size const ) const {
		utility_exit_with_message("ResidueTypeLinkingConstraint is not atom-based!.");
		return core::id::GLOBAL_BOGUS_ATOM_ID;  // required for compilation on Windows
	}

	virtual
	utility::vector1< core::Size >
	residues() const;

	void
	show( std::ostream & out ) const;
	/*
	virtual
	ConstraintOP
	remap_resid( core::id::SequenceMapping const &seqmap ) const;
	*/

	/// @brief possibility to compare constraint according to data
	/// and not just pointers
	virtual bool operator == ( Constraint const & other ) const;
	virtual bool same_type_as_me( Constraint const & ) const;

	/*
	virtual ConstraintOP remapped_clone(
	pose::Pose const & src,
	pose::Pose const & dest,
	id::SequenceMappingCOP map = NULL
	) const;
	*/

	// Needed to get the base class overloads
	using Constraint::score;
	using Constraint::dist;

	virtual
	void
	score( func::XYZ_Func const & xyz_func, EnergyMap const & weights, EnergyMap & emap ) const;

	/// @details Return 1.0 if constraint will get a bonus/penalty, 0.0 if not
	virtual
	core::Real
	dist( core::scoring::func::XYZ_Func const & xyz ) const;

	virtual
	void
	fill_f1_f2(
		AtomID const & atom,
		func::XYZ_Func const & xyz,
		Vector & F1,
		Vector & F2,
		EnergyMap const & weights
	) const;

	virtual ConstraintOP
	clone() const;
private:
	Size seqpos1_;
	Size seqpos2_;
	std::string AA1name;
	std::string AA2name;
	std::string rsd1_type_name3_;
	std::string rsd2_type_name3_;
	core::Real bonus_;
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
