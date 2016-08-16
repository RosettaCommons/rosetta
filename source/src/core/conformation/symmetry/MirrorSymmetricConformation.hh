// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/conformation/symmetry/MirrorSymmetricConformation.hh
/// @brief  Symmetry conformation container.
/// @details Contains overloaded functions needed to make changes in conformation symmetric and
/// to handle mirror symmetry operations.
/// @author Frank DiMaio
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_core_conformation_symmetry_MirrorSymmetricConformation_hh
#define INCLUDED_core_conformation_symmetry_MirrorSymmetricConformation_hh

// Unit headers
#include <core/conformation/symmetry/MirrorSymmetricConformation.fwd.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>

#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryTransform.hh>

#include <utility/vector1.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace conformation {
namespace symmetry {

/// @brief Mirror symmetric conformation contains the same symminfo logic, but also
///        includes special logic for mirror symmetries
class MirrorSymmetricConformation : public SymmetricConformation {

public:
	/////////////////////////////////////////
	//// Construction / Copying
	/////////////////////////////////////////

	/// @brief  Default CTOR
	MirrorSymmetricConformation();

	/// @brief
	MirrorSymmetricConformation(
		Conformation const & conf,
		SymmetryInfo const & symm_info);

	/// @brief  Copy CTOR
	MirrorSymmetricConformation(
		MirrorSymmetricConformation const & conf);

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

	/////////////////////////////////////////
	//// Setters
	/////////////////////////////////////////

	/// DOF
	virtual
	void
	set_dof( DOF_ID const & id, Real setting );

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
	set_torsion_angle(
		AtomID const & atom1,
		AtomID const & atom2,
		AtomID const & atom3,
		AtomID const & atom4,
		Real setting,
		bool quiet=false
	);

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

	/// @brief set the fold_tree
	virtual void
	fold_tree( FoldTree const & fold_tree_in );

	// @brief force recomputation of Tsymm_'s from the current conformation
	virtual void
	recalculate_transforms( );

	virtual
	~MirrorSymmetricConformation();

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

	virtual
	void
	detect_disulfides( utility::vector1< Size > const & disulf_one = utility::vector1< core::Size >(), utility::vector1< Size > const & disulf_two = utility::vector1< core::Size >() );

	/// @brief Updates residue identities in symmetric subunits, ensuring that they are mirrored relative to the ASU in mirrored subunits
	/// and identical to the ASU in non-mirrored subunits.
	/// @details Assumes that the residue identities and variants (aside from D/L variants) already match.  That is, if I have ASN at position
	/// 5 in my asymmetric unit, I either have ASN or DASN at the equivalent position in each symmetry copy.  Safe to call repeatedly.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	void update_residue_identities( );

	/// @brief Is this residue mirrored relative to the asymmetric unit?
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	bool res_is_mirrored ( core::Size const seqpos ) const;

private:

	/// @brief Helper function to flip the chirality of a residue type.
	/// @details Assumes that the coordinates are already correct for the flipped type (i.e. this function does
	/// not alter atom positions -- only the ResidueType identity).
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void flip_chirality( ResidueOP & new_rsd );

	/// @brief helper functions called when #jumps or #residues changes, updating internal data
	void update_njumps_nres( );

	/// @brief helper function ensures pose jumps are in synch with 'jump_is_mirrored_'
	void synch_mirror_jumps_with_atomtree( );

	/// @brief keep track of which jumps and residues are mirrored w.r.t. master
	utility::vector1< bool > res_is_mirrored_;
	utility::vector1< std::pair< bool,bool > > jump_is_mirrored_;

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
CEREAL_FORCE_DYNAMIC_INIT( core_conformation_symmetry_MirrorSymmetricConformation )
#endif // SERIALIZATION


#endif
