// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//////////////////////////////////////////////////////////////////////
/// @begin ResidueType
///
/// @brief
/// A class for defining residue
///
/// @detailed
/// This class contains the "chemical" information for residues. This does not contain the actual
/// xyz coordinates of the class (xyz found in core/conformation/Residue.hh). A residue in Rosetta
/// can be a ligand, DNA, amino acid, or basically anything. A residue is read in through residue_io.cc
/// and read from parameter files, generally located in the database chemical/residuetypes. For ligands,
/// or anything that is not the natural 20 aa, a parameter has to be provided to rosetta through the -extra_res_fa
/// flag. Residue_io sets private member data in ResidueType. The primary data that are set are: atoms, mmatoms, orbitals,
/// properties of residues. These properties can be modified through patches, which is controlled through PatchOperations.cc. If
/// the residuetype is modified, the indices of atoms and mmatoms and everything associated with those indices must be redefined. This
/// redordering of indices is taken care of with the function reorder_primary_data().
///
/// Setting of primary data and then reordering is important. Primary data for the following are described:
///
/// atoms: setting of atoms includes indexing the atoms into vectors, saving their names into vectors/maps,
/// saving the associated mm_atom_type into a vector, saving bond connections into vectors, etc, etc. Since everything is
/// allocated into vectors, it is easy to reorder those vectors. On any given residue, the heavy atoms are put into the vector
/// first (their indices are first) and hydrogens are put in last.
///
/// properties: properties of a residue include things like DNA, PROTEIN, SC_ORBITALS, CHARGED, etc. These properties
/// indicate the type of residue it is and what properties that are asscociated with the residue. They
/// are set when read in. Several lines of code must be modified to get them to work, all found in residuetype.cc
///
/// orbitals: orbitals are indexed seperate from atoms. They function much the same way as atoms, except for some
/// key differences. To find atoms bonded to orbitals, you must provide the atom index, not the orbital index. I
/// havent figured out how to get the reverse to work because of the seperate indices. Orbital xyz coordinates are not updated when atom coordinates are.
/// This is to keep speed consistent with just having atoms. To output the orbitals, use the flag -output_orbitals
///
/// @authors
/// Phil Bradley
/// Steven Combs - these comments
///
///
////////////////////////////////////////////////////////////////////////



#ifndef INCLUDED_core_chemical_rna_RNA_ResidueType_hh
#define INCLUDED_core_chemical_rna_RNA_ResidueType_hh

// Package headers
#include <core/chemical/rna/RNA_ResidueType.fwd.hh>
#include <core/chemical/AtomType.fwd.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/Adduct.hh>
#include <core/chemical/AtomTypeSet.fwd.hh>
#include <core/chemical/ElementSet.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/MMAtomType.fwd.hh>
#include <core/chemical/MMAtomTypeSet.fwd.hh>
#include <core/chemical/ResidueConnection.hh>
// AUTO-REMOVED #include <core/chemical/VariantType.hh>


// Numeric headers
#include <numeric/xyzVector.hh>

// Utility headers
// AUTO-REMOVED #include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/keys/Key2Tuple.hh>
#include <utility/keys/Key4Tuple.hh>
#include <utility/keys/Key3Tuple.hh>

// C++ headers
#include <map>

#include <core/chemical/VariantType.fwd.hh>
#include <utility/vector1.hh>


namespace core {
namespace chemical {
namespace rna {


/////////////////////////////////////Parin S, Dec 23, 2011: RNA stuff//////////////////////////////////////////
class RNA_ResidueType : public utility::pointer::ReferenceCount {


	public:

	RNA_ResidueType();

	~RNA_ResidueType(){};

	///////////////////////////Implemented for fast lookup! Parin Sripakdeevong, June 25th, 2011//////////////////
	private:

	void
	rna_note_chi_controls_atom( core::Size const chi, core::Size const atomno, utility::vector1< core::Size > & last_controlling_chi);


	public:

	void
	update_derived_rna_data(ResidueTypeCOP const residue_type_in);

	void
	rna_update_last_controlling_chi(ResidueTypeCOP const residue_type_in,
																utility::vector1< core::Size > & last_controlling_chi,
																utility::vector1< AtomIndices >  &  atoms_last_controlled_by_chi);
	
	////////////Fast lookup functions///////

	utility::vector1< bool > const &
	Is_virtual_atom_list() const;

	bool
	atom_is_virtual( Size const atomno ) const;

	utility::vector1< bool > const &
	Is_phosphate_atom_list() const;

	/// @brief quick lookup: is the atom with the given index is part of the RNA phosphate or not?
	bool
	atom_is_phosphate( Size const atomno ) const;

	utility::vector1< bool > const &
	Is_RNA_base_atom_list() const;

	bool
	is_RNA_base_atom( Size const atomno ) const;

	AtomIndices const &
	RNA_base_atoms() const;

	Size
	ho2star_index() const;

	Size
	o2star_index() const;

	Size
	p_atom_index() const;

	Size
	o1p_atom_index() const;

	Size
	o2p_atom_index() const;

	Size
	o5star_atom_index() const;

	Size
	o3star_atom_index() const;

	Size
	o4star_atom_index() const;

	Size
	c1star_atom_index() const;

	Size
	c2star_atom_index() const;

	Size
	c4star_atom_index() const;

	///////////////////////////////////////////////////////////////////////////////////////

	public:

	//o2star atom
	core::Size o2star_index_;
	core::Size ho2star_index_;

	//Phosphate atoms
	core::Size p_atom_index_;
	core::Size o1p_atom_index_;
	core::Size o2p_atom_index_;
	core::Size o5star_index_;
	core::Size o3star_index_;

	//Ribose atoms
	core::Size o4star_index_;
	core::Size c1star_index_;
	core::Size c2star_index_;
	core::Size c4star_index_;

	AtomIndices base_atom_list_;

	utility::vector1< bool > Is_RNA_base_atom_list_;
	utility::vector1< bool > Is_phosphate_atom_list_;
	//For fast look whether atom is VIRTUAL type.
	utility::vector1< bool > Is_virtual_atom_list_;
	//AtomIndices base_heavy_atoms_;
	//AtomIndices base_hydrogen_atoms_;
	
	/// atom index lookup by atom name string

	ResidueTypeCOP residue_type_; //Pointer to the main ResidueType object that this RNA_ResidueType object belongs to.

}; //RNA_ResidueType



} // rna
} // chemical
} // core

#endif // INCLUDED_core_chemical_rna_RNA_ResidueType_hh
