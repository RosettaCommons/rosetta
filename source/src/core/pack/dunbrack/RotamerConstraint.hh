// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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

#include <core/conformation/Conformation.fwd.hh>

#include <core/id/AtomID.hh>
#include <core/pack/dunbrack/DunbrackRotamer.fwd.hh>
#include <core/pack/dunbrack/RotamerLibrary.fwd.hh>
#include <core/pack/dunbrack/SingleResidueRotamerLibrary.fwd.hh>
#include <utility/vector1_bool.hh>
#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace dunbrack {


/// @brief Convenience function adds constraints to the supplied pose based on
/// the -unboundrot command line flag.
void load_unboundrot(pose::Pose & pose);

/// @brief Convenience function adds constraints to the supplied pose based on
/// the list of provided poses.
void load_unboundrot(pose::Pose & pose, core::pose::PoseCOPs const & unboundrot_poses);


///@brief This class favors a particular rotamer at a particular position by reducing its Dunbrack energy.
///
///@details Specifically, the given rotamer well(s) receive a bonus equal to the
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
	clone() const { return new RotamerConstraint( *this ); }

	virtual
	Size
	natoms() const;

	virtual
	AtomID const &
	atom( Size const index ) const;

	virtual
	void
	score(
		scoring::func::XYZ_Func const & xyz_func,
		scoring::EnergyMap const & weights,
		scoring::EnergyMap & emap
	) const;

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
	SingleResidueRotamerLibraryCAP rotlib_;
	utility::vector1< ChiVector > favored_rotamers_;
	utility::vector1< RotVector > favored_rotamer_numbers_;

}; // RotamerConstraint


} // namespace dunbrack
} // namespace pack
} // namespace core

#endif // INCLUDED_core_pack_dunbrack_RotamerConstraint_HH
