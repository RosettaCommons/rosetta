// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/constraints/RotamerConstraint.hh
///
/// @brief
/// @author Ian W. Davis


#ifndef INCLUDED_core_pack_dunbrack_RotamerConstraint_hh
#define INCLUDED_core_pack_dunbrack_RotamerConstraint_hh

// Unit headers
#include <core/pack/dunbrack/RotamerConstraint.fwd.hh>

// Package headers
#include <core/scoring/constraints/Constraint.hh>

// Project headers
#include <core/chemical/ResidueType.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/id/AtomID.hh>
#include <core/pack/dunbrack/DunbrackRotamer.fwd.hh>
#include <core/pack/dunbrack/RotamerLibrary.fwd.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibrary.fwd.hh>
#include <utility/vector1_bool.hh>
#include <utility/vector1.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace pack {
namespace dunbrack {


/// @brief Convenience function adds constraints to the supplied pose based on
/// the -unboundrot command line flag.
void load_unboundrot(pose::Pose & pose);

/// @brief Convenience function adds constraints to the supplied pose based on
/// the list of provided poses.
void load_unboundrot(pose::Pose & pose, core::pose::PoseCOPs const & unboundrot_poses);


/// @brief This class favors a particular rotamer at a particular position by reducing its Dunbrack energy.
///
/// @details Specifically, the given rotamer well(s) receive a bonus equal to the
/// difference between their minimum energy (for the ideal rotamer in this well)
/// and the minimum energy of any rotamer (at the current phi,psi).
/// This class is used to implement the scoring component of the -unboundrot flag;
/// actually adding those rotamers to the library occurs in core/pack/rotamer_set/UnboundRotamersOperation.
class RotamerConstraint : public scoring::constraints::Constraint
{
public:
	virtual std::string type() const {
		return "Rotamer";
	}

	RotamerConstraint();

	RotamerConstraint(RotamerConstraint const & other);

	RotamerConstraint(
		pose::Pose const & pose,
		Size seqpos
	);

	virtual ~RotamerConstraint();

	virtual
	void
	add_residue( conformation::Residue const & rsd );

	virtual
	scoring::constraints::ConstraintOP
	clone() const;

	virtual bool operator == ( Constraint const & other ) const;
	virtual bool same_type_as_me( Constraint const & other ) const;

	virtual
	Size
	natoms() const;

	virtual
	AtomID const &
	atom( Size const index ) const;

	// Needed to get the base class overloads
	using Constraint::score;
	using Constraint::dist;

	virtual
	void
	score(
		scoring::func::XYZ_Func const & xyz_func,
		scoring::EnergyMap const & weights,
		scoring::EnergyMap & emap
	) const;

	/// @details Is ithere a "distance" that we can give for RotamerConstraints?
	virtual
	core::Real
	dist( core::scoring::func::XYZ_Func const & ) const { return 0; }

	virtual
	void
	fill_f1_f2(
		AtomID const & atom,
		scoring::func::XYZ_Func const & xyz,
		Vector & F1,
		Vector & F2,
		scoring::EnergyMap const & weights
	) const;

	/// @brief This gets used to compare one constraint to another, so it should
	/// uniquely reflect all the parameters.
	virtual void show( std::ostream & out ) const ;

private:
	Size seqpos_;
	std::string rsd_type_name_;
	utility::vector1< AtomID > atom_ids_;
	core::chemical::ResidueTypeCOP restype_; // for serialization purposes
	core::pack::rotamers::SingleResidueRotamerLibraryCOP rotlib_;
	utility::vector1< ChiVector > favored_rotamers_;
	utility::vector1< RotVector > favored_rotamer_numbers_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; // RotamerConstraint


} // namespace dunbrack
} // namespace pack
} // namespace core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_pack_dunbrack_RotamerConstraint )
#endif // SERIALIZATION


#endif // INCLUDED_core_pack_dunbrack_RotamerConstraint_HH
