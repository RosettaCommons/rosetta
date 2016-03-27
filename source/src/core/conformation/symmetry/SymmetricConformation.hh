// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/conformation/symmetry/SymmetricConformation.hh
/// @brief  symmetry conformation container.
//     Contains overloaded functions needed to
//     make changes in conformation symmetric
/// @author Phil Bradley, Ingemar Andre


#ifndef INCLUDED_core_conformation_symmetry_SymmetricConformation_hh
#define INCLUDED_core_conformation_symmetry_SymmetricConformation_hh


// Unit headers
#include <core/conformation/symmetry/SymmetricConformation.fwd.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>
#include <core/conformation/symmetry/SymmetryTransform.hh> // needs full decl

#include <core/conformation/Conformation.hh>

// Numeric
#include <numeric/HomogeneousTransform.hh>

#include <utility/vector1.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace conformation {
namespace symmetry {

/// @brief  A symmetric conformation: has an additional data member "SymmetryInfo" class
/// @details  Handles symmetrizing of *some* of the conformation-changing methods of Conformation

class SymmetricConformation : public core::conformation::Conformation {

public:

	/////////////////////////////////////////
	//// Construction / Copying
	/////////////////////////////////////////

	/// @brief  Default CTOR
	SymmetricConformation();

	/// @brief  Default CTOR
	SymmetricConformation( Conformation const & conf, SymmetryInfo const & symm_info );

	/// @brief copy constructor
	//SymmetricConformation( SymmetricConformation const & src );

	/// @brief virtual assignment operator
	virtual
	Conformation &
	operator=( Conformation const & src );

	virtual
	void
	detached_copy( Conformation const & src );

	ConformationOP
	clone() const;

	virtual
	bool
	same_type_as_me( Conformation const & other, bool recurse  /* = true */ ) const;

	virtual
	SymmetryInfoCOP Symmetry_Info() const;

	virtual
	SymmetryInfoOP Symmetry_Info();

	/////////////////////////////////////////
	//// Setters
	/////////////////////////////////////////

	/// DOF
	virtual
	void
	set_dof( DOF_ID const & id, Real setting );

	void
	set_secstruct( Size seqpos, char setting );

	/// BONDS/TORSIONS
	virtual
	void
	set_torsion( TorsionID const & id, Real setting );

	/// JUMPS
	/// @brief set a jump
	virtual
	void
	set_jump(
		int jump_number,
		Jump const & new_jump
	);

	/// @brief set a jump
	virtual
	void
	set_jump(
		AtomID const & id,
		Jump const & new_jump
	);


	virtual
	void
	set_bond_angle(
		AtomID const & atom1,
		AtomID const & atom2,
		AtomID const & atom3,
		Real setting
	);


	virtual
	void
	set_bond_length(
		AtomID const & atom1,
		AtomID const & atom2,
		Real setting
	);


	virtual
	void
	set_torsion_angle(
		AtomID const & atom1,
		AtomID const & atom2,
		AtomID const & atom3,
		AtomID const & atom4,
		Real setting,
		bool quiet=false
	);

	virtual
	utility::vector1<bool>
	get_residue_mask() const;

	virtual Real
	get_residue_weight(core::Size resid1, core::Size resid2) const;

	/// @brief replace residue
	virtual void
	replace_residue(
		Size seqpos,
		Residue const & new_rsd,
		bool orient_backbone
	);

	virtual void
	replace_residue(
		Size seqpos,
		Residue const & new_rsd,
		utility::vector1< std::pair< std::string, std::string > > const & atom_pairs
	);

	/// @brief set the fold_tree .. update symminfo if jump numbering changed
	virtual void
	fold_tree( FoldTree const & fold_tree_in );

	/// @brief FoldTree access
	virtual FoldTree const &
	fold_tree() const
	{
		return Conformation::fold_tree();
	}

	/// @brief Get the transformation controlling resid i
	numeric::HomogeneousTransform< core::Real >
	get_transformation( core::Size resid );

	/// @brief Remap coordinate X from resid i's frame to resid j's frame
	virtual
	PointPosition
	apply_transformation( PointPosition Xin, core::Size residfrom, core::Size residto, bool rotationonly=false );

	/// @brief Remap coordinate X from resid i's frame to resid j's frame
	///   assumes that the transformations are already computed (thus can be const)
	virtual
	PointPosition
	apply_transformation_norecompute( PointPosition Xin, core::Size residfrom, core::Size residto, bool rotationonly=false ) const;

	// @brief force recomputation of Tsymm_'s from the current conformation
	virtual
	void
	recalculate_transforms( );

	/// @brief Symmetric set_xyz
	virtual void
	set_xyz( AtomID const & id, PointPosition const & position );

	/// @brief Symmetric batch_set_xyz
	virtual void
	batch_set_xyz(
		utility::vector1<AtomID> const & ids,
		utility::vector1<PointPosition> const & positions
	);

	virtual
	void
	apply_transform_Rx_plus_v(
		numeric::xyzMatrix< Real > const & R,
		Vector const & v
	);

	virtual
	~SymmetricConformation();

	/// @brief Append a new residue by a jump; clones this append to all copies
	void
	append_residue_by_jump(
		conformation::Residue const & new_rsd,
		Size anchor_residue,
		std::string const& anchor_atom = "", // the atom in the anchor_residue
		std::string const& root_atom = "", // the atom in the new residue
		bool start_new_chain = false
	);

	/// @brief Append a new conformation by a jump; clones this append to all copies
	void
	insert_conformation_by_jump(
		Conformation const & conf,             // the conformation to be inserted
		Size insert_seqpos,              // rsd 1 in conf goes here
		Size insert_jumppos,             // jump#1 in conf goes here, see insert_fold_tree_by_jump
		Size anchor_pos,                 // in the current sequence numbering, ie before insertion of conf
		Size anchor_jump_number = 0,     // the desired jump number of the anchoring jump, default=0
		std::string const & anchor_atom = "",  // "" means take default anchor atom
		std::string const & root_atom   = ""   // "" means take default root   atom
	);

	//fpd eventually we should have symmetric implementations of all the insert/append/delete residue functions

	virtual
	void
	detect_disulfides( utility::vector1< Size > const & disulf_one = utility::vector1< core::Size >(), utility::vector1< Size > const & disulf_two = utility::vector1< core::Size >() );

	/// @brief Declare that a chemical bond exists between two residues
	/// @details This updates all symmetry copies, so that each one has a chemical
	/// bond between the residues in question.
	/// @author Frank DiMaio.
	/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
	virtual
	void
	declare_chemical_bond(
		Size seqpos1,
		std::string const & atom_name1,
		Size seqpos2,
		std::string const & atom_name2
	);

protected:

	// @brief invalidate current Tsymm settings
	void
	clear_Tsymm( );

	// @brief  invert one of the Tsymm transforms about Z
	void
	invert_Tsymm( char sub, core::Size subunit );

	// utility function gets the nearest upstream virtual
	core::Size
	get_upstream_vrt( core::Size seqpos ) const;


private:

	SymmetryInfoOP symm_info_;

	// stores the symmetric transformation between subunits
	//   computed when needed, invalidated when a jump changes
	// multicomp: store transforms for each component (indexed by character)
	//std::map< char, utility::vector1< numeric::HomogeneousTransform< core::Real > > > Tsymm_;
	std::map< char, utility::vector1<SymmetryTransform> > Tsymm_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} // symmetry
} // conformation
} // core


#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_conformation_symmetry_SymmetricConformation )
#endif // SERIALIZATION


#endif
