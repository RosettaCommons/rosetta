// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   core/conformation/symmetry/SymmetricConformation.hh
/// @brief  symmetry conformation container.
//					Contains overloaded functions needed to
//					make changes in conformation symmetric
/// @author Phil Bradley, Ingemar Andre


#ifndef INCLUDED_core_conformation_symmetry_SymmetricConformation_hh
#define INCLUDED_core_conformation_symmetry_SymmetricConformation_hh


// Unit headers
#include <core/conformation/symmetry/SymmetricConformation.fwd.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>

#include <core/conformation/Conformation.hh>

// Numeric
#include <numeric/HomogeneousTransform.hh>

#include <utility/vector1.hh>


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
	SymmetricConformation( SymmetricConformation const & src );

	/// @brief operator
	Conformation &
	operator=( SymmetricConformation const & src );

  ConformationOP
  clone() const;

	virtual
	bool
	same_type_as_me( Conformation const & other, bool recurse  /* = true */ ) const;

	SymmetryInfoCOP Symmetry_Info() const;

	SymmetryInfoOP Symmetry_Info();

	/////////////////////////////////////////
	//// Setters
	/////////////////////////////////////////

	/// DOF
	virtual
	void
	set_dof( DOF_ID const & id, Real const setting );

	void
	set_secstruct( Size const seqpos, char const setting );

	/// BONDS/TORSIONS
	virtual
	void
	set_torsion( TorsionID const & id, Real const setting );

	/// JUMPS
	/// @brief set a jump
  virtual
  void
  set_jump(
    int const jump_number,
    Jump const & new_jump
  );

	/// @brief set a jump
  virtual
  void
  set_jump(
    AtomID const & id,
    Jump const & new_jump
  );

	///
	virtual
	void
	set_bond_angle(
		AtomID const & atom1,
		AtomID const & atom2,
		AtomID const & atom3,
		Real const setting
	);

	///
	virtual
	void
	set_bond_length(
		AtomID const & atom1,
		AtomID const & atom2,
		Real const setting
	);

	///
	virtual
	void
	set_torsion_angle(
		AtomID const & atom1,
		AtomID const & atom2,
		AtomID const & atom3,
		AtomID const & atom4,
		Real const setting,
		bool const quiet=false
	);

	virtual utility::vector1<bool>
	get_residue_mask() const;

	virtual Real
	get_residue_weight(core::Size resid1, core::Size resid2) const;

	/// @brief replace residue
	virtual void
	replace_residue(
		Size const seqpos,
		Residue const & new_rsd,
		bool const orient_backbone
	);

	virtual void
	replace_residue(
		Size const seqpos,
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
	PointPosition
	apply_transformation( PointPosition Xin, core::Size residfrom, core::Size residto );

	/// @brief Remap coordinate X from resid i's frame to resid j's frame
	PointPosition
	apply_transformation_norecompute( PointPosition Xin, core::Size residfrom, core::Size residto ) const;

	// @brief force recomputation of Tsymm_'s from the current conformation
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

	///
	virtual
	~SymmetricConformation();

	/// @brief Append a new residue by a jump; clones this append to all copies
	void
	append_residue_by_jump(
		conformation::Residue const & new_rsd,
		Size const anchor_residue,
		std::string const& anchor_atom = "", // the atom in the anchor_residue
		std::string const& root_atom = "", // the atom in the new residue
		bool const start_new_chain = false
	);

	/// @brief Append a new conformation by a jump; clones this append to all copies
  void
  insert_conformation_by_jump(
    Conformation const & conf,             // the conformation to be inserted
    Size const insert_seqpos,              // rsd 1 in conf goes here
    Size const insert_jumppos,             // jump#1 in conf goes here, see insert_fold_tree_by_jump
    Size const anchor_pos,                 // in the current sequence numbering, ie before insertion of conf
    Size const anchor_jump_number = 0,     // the desired jump number of the anchoring jump, default=0
    std::string const & anchor_atom = "",  // "" means take default anchor atom
    std::string const & root_atom   = ""   // "" means take default root   atom
  );

	//fpd eventually we should have symmetric implementations of all the insert/append/delete residue functions

	virtual
	void
	detect_disulfides();


private:

	/////////////////////////////////////////////////////////////////////////////
	// private methods
	/////////////////////////////////////////////////////////////////////////////

	// public: -- apl -- does this need to be public?

	// utility function gets the nearest upstream virtual
	core::Size
	get_upstream_vrt( core::Size seqpos ) const;

#ifdef USEBOOSTSERIALIZE
	friend class boost::serialization::access;

	template<class Archive>
	void serialize(Archive & ar, const unsigned int version){
		ar & boost::serialization::base_object<Conformation>(*this);
		ar & symm_info_;
	}

#endif


private:

	/////////////////////////////////////////////////////////////////////////////
	// data
	/////////////////////////////////////////////////////////////////////////////

	SymmetryInfoOP symm_info_;

	// stores the symmetric transformation between subunits
	// computed when needed, invalidated when a jump changes
	utility::vector1< numeric::HomogeneousTransform< core::Real > > Tsymm_;

};

} // symmetry
} // conformation
} // core



#endif
