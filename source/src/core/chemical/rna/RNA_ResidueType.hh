// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//////////////////////////////////////////////////////////////////////
/// @begin RNA_residueType
///
/// @brief
/// RNA specific properties
///
/// @author Parin Sripakdeevong (sripakpa@stanford.edu)
////////////////////////////////////////////////////////////////////////



#ifndef INCLUDED_core_chemical_rna_RNA_ResidueType_hh
#define INCLUDED_core_chemical_rna_RNA_ResidueType_hh


// Package headers
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/rna/RNA_ResidueType.fwd.hh>

// Numeric headers

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <utility/vector1.hh>


namespace core {
namespace chemical {
namespace rna {


/////////////////////////////////////Parin S, Dec 23, 2011: RNA stuff//////////////////////////////////////////
class RNA_ResidueType : public utility::pointer::ReferenceCount {


public:

	RNA_ResidueType();

	virtual ~RNA_ResidueType();

///////////////////////////Implemented for fast lookup! Parin Sripakdeevong, June 25th, 2011//////////////////
private:

	utility::vector1< Size > const
	figure_out_chi_order() const;

	void
	rna_note_chi_controls_atom( core::Size const chi, core::Size const atomno,
															utility::vector1< core::Size > & last_controlling_chi,
															utility::vector1< core::Size > const & chi_order );

	public:

	void
	update_derived_rna_data( ResidueTypeCOP const residue_type_in );

	void
	rna_update_last_controlling_chi( ResidueTypeCOP const residue_type_in,
																utility::vector1< core::Size > & last_controlling_chi,
																utility::vector1< AtomIndices >  &  atoms_last_controlled_by_chi );

	////////////Fast lookup functions///////

	utility::vector1< bool > const &
	is_virtual() const;

	bool
	atom_is_virtual( Size const atomno ) const;

	utility::vector1< bool > const &
	is_phosphate() const;

	/// @brief quick lookup: is the atom with the given index is part of the RNA phosphate or not?
	bool
	atom_is_phosphate( Size const atomno ) const;

	utility::vector1< bool > const &
	is_RNA_base() const;

	bool
	is_RNA_base_atom( Size const atomno ) const;

	AtomIndices const &
	RNA_base_atoms() const;

	Size
	ho2prime_index() const;

	Size
	o2prime_index() const;

	Size
	p_atom_index() const;

	Size
	op2_atom_index() const;

	Size
	op1_atom_index() const;

	Size
	o5prime_atom_index() const;

	Size
	o3prime_atom_index() const;

	Size
	o4prime_atom_index() const;

	Size
	c1prime_atom_index() const;

	Size
	c2prime_atom_index() const;

 	Size
	c4prime_atom_index() const;

	Size chi_number_pseudoalpha() const { return chi_number_pseudoalpha_;}
	Size chi_number_pseudobeta() const { return chi_number_pseudobeta_;}
	Size chi_number_pseudogamma() const { return chi_number_pseudogamma_;}
	Size chi_number_pseudoepsilon() const { return chi_number_pseudoepsilon_;}
	Size chi_number_pseudozeta() const { return chi_number_pseudozeta_;}

///////////////////////////////////////////////////////////////////////////////////////

public:

	//o2prime atom
	core::Size o2prime_index_;
	core::Size ho2prime_index_;

	//Phosphate atoms
	core::Size p_atom_index_;
	core::Size op2_atom_index_;
	core::Size op1_atom_index_;
	core::Size o5prime_index_;
	core::Size o3prime_index_;

	//Ribose atoms
	core::Size o4prime_index_;
	core::Size c1prime_index_;
	core::Size c2prime_index_;
	core::Size c4prime_index_;

	AtomIndices base_atom_list_;

	utility::vector1< bool > is_RNA_base_;
	utility::vector1< bool > is_phosphate_;
	//For fast look whether atom is VIRTUAL type.
	utility::vector1< bool > is_virtual_;

	ResidueTypeCOP residue_type_; //Pointer to the main ResidueType object that this RNA_ResidueType object belongs to.

	Size chi_number_pseudoalpha_;
	Size chi_number_pseudobeta_;
	Size chi_number_pseudogamma_;
	Size chi_number_pseudoepsilon_;
	Size chi_number_pseudozeta_;

}; //RNA_ResidueType



} // rna
} // chemical
} // core

#endif // INCLUDED_core_chemical_rna_RNA_ResidueType_hh
