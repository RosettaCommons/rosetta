// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/drug_design/substituion_support.hh
/// @brief use RDKit to substitute items based on matched templates
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_protocols_drug_design_substitution_support_hh
#define INCLUDED_protocols_drug_design_substitution_support_hh

#include <protocols/drug_design/substitution_support.fwd.hh>

#include <core/chemical/rdkit/RDKit.fwd.hh>
#include <core/chemical/rdkit/util.hh>

#include <core/chemical/AtomRefMapping.hh>

#include <utility/pointer/owning_ptr.hh>

#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/Geometry/Transform3D.h>
#include <rdkit/GraphMol/Substruct/SubstructMatch.h>

namespace protocols {
namespace drug_design {

/// @brief Get the RDKit bond type of the first bond to atom
::RDKit::Bond::BondType
get_first_bondtype( ::RDKit::ROMol const & mol, ::RDKit::Atom const * atom );

/// @brief Get the RDKit bond type of the first bond to the given atom
::RDKit::Bond::BondType
get_first_bondtype( ::RDKit::ROMol const & mol, unsigned int idx );

/// @brief Get the RDKit bond of the first bond to the given atom
::RDKit::Bond const &
get_first_bond( ::RDKit::ROMol const & mol, unsigned int idx );

/// @brief Get the atom number of the atom bonded to idx (it's the first atom)
unsigned int
get_bonded_atom( ::RDKit::ROMol const & mol, unsigned int idx );

/// @brief Is the given atom a dummy (atom number 0)?
bool
is_dummy( ::RDKit::ROMol const & mol, unsigned int idx );

/// @brief Populate the list with the index of all the dummies in the molecule
/// @details Existing contents will be cleared.
void
find_dummies( ::RDKit::ROMol const & mol, utility::vector1< unsigned int > & dummy_list );

///// @brief Is the bond in a ring?
//bool
//bond_is_in_ring( ::RDKit::Bond const & bond );

/// @brief Is the bond between the atoms in a ring?
bool
bond_is_in_ring( ::RDKit::ROMol const & rdmol, unsigned int atm1, unsigned int atm2 );

/// @brief Convert the RDKit MatcVect into a map (real molecule index to query index)
std::map< unsigned int, unsigned int >
convert_matchvect_to_map( ::RDKit::MatchVectType const & pairings );

enum MoleculeSelection {
	OriginalMol,
	TemplateMol,
	ReplaceMol,
	NewMol
};

/// @brief A class representing a substitution of an atom from an original molecule
/// through a pair of matched templates to a substituted molecule.
/// Intended to be used through MoleculeSubstitution
class AtomSubstitution {

public:
	AtomSubstitution(
		unsigned int mdx = invalid_index,
		unsigned int tdx = invalid_index,
		unsigned int rdx = invalid_index,
		unsigned int ndx = invalid_index ):
		mdx_(mdx),
		tdx_(tdx),
		rdx_(rdx),
		ndx_(ndx)
	{}

	unsigned int idx(MoleculeSelection sele) const;

	/// @brief The atom index in the original molecule
	unsigned int mdx() const { return mdx_; }
	/// @brief The atom index in the template which matches the original molecule
	unsigned int tdx() const { return tdx_; }
	/// @brief The atom index in the replacement template
	unsigned int rdx() const { return rdx_; }
	/// @brief The atom index in the new (post-replacement) molecule
	unsigned int ndx() const { return ndx_; }

	void set_idx(MoleculeSelection sele, unsigned int setting);

	void set_mdx(unsigned int setting) { mdx_ = setting; }
	void set_tdx(unsigned int setting) { tdx_ = setting; }
	void set_rdx(unsigned int setting) { rdx_ = setting; }
	void set_ndx(unsigned int setting) { ndx_ = setting; }

public:

	static unsigned int const invalid_index=65500; // Near the top end of the range for unsigned int

private:

	unsigned int mdx_;
	unsigned int tdx_;
	unsigned int rdx_;
	unsigned int ndx_;

};

/// @brief A class representing a substitution of a molecule from an original molecule
/// through a pair of matched templates to a substituted molecule.
/// @details The MoleculeSubstitution can be made progressively from molecule outwards.
/// It should contain an indexed entry for each item in all four molecules
/// (At least for those molecules which it's currently aware of.)
class MoleculeSubstitution {

public:

	MoleculeSubstitution()
	{}

	/// @brief Construct a MoleculeSubstitution from an RDKit molecule
	/// (The source molecule.)
	MoleculeSubstitution( ::RDKit::ROMolOP mol );

	/// @brief Add a template to this MoleculeSubstitution.
	void
	add_template(::RDKit::ROMolOP templt, ::RDKit::MatchVectType const & pairings );

	/// @brief Reconstitute the original pairings.
	/// This is a vector of pairs of (tdx, mdx)
	::RDKit::MatchVectType
	make_match_vector() const;

	/// @brief Create a mapping from mdx to ndx
	core::chemical::IndexIndexMapping
	find_mdx_to_ndx_mapping() const;

	/// @brief Make a *new* MoleculeSubstitution class with the replacement information
	/// @details Assumes that the new (post-replacement) molecule hasn't been set.
	MoleculeSubstitutionOP
	add_replacement(
		::RDKit::ROMolOP replacement,
		std::map< unsigned int, unsigned int > & r_to_t_mapping
	) const;

	void
	add_newmol( ::RDKit::RWMolOP newmol) { newmol_ = newmol; }

	::RDKit::ROMolOP get_romol( MoleculeSelection sele );

	::RDKit::ROMolOP mol() { return mol_; }
	::RDKit::ROMolOP templt() { return templt_; }
	::RDKit::ROMolOP replace() { return replace_; }
	::RDKit::RWMolOP newmol() { return newmol_; }

	AtomSubstitution & substitution_for_idx( MoleculeSelection sele, unsigned int idx );

	AtomSubstitution & substitution_for_mdx( unsigned int idx );
	AtomSubstitution & substitution_for_tdx( unsigned int idx );
	AtomSubstitution & substitution_for_rdx( unsigned int idx );
	//AtomSubstitution & substitution_for_ndx( unsigned int idx );

	bool tdx_is_dummy( unsigned int tdx ) const;
	bool rdx_is_dummy( unsigned int rdx ) const;

	/// @brief The indexes for dummy atoms on the matching template
	utility::vector1< unsigned int > const &
	template_dummies() const { return template_dummies_; }

	/// @brief The indexes for dummy atoms on the replacement template
	utility::vector1< unsigned int > const &
	replace_dummies() const { return replace_dummies_; }

private:

	::RDKit::ROMolOP mol_; // Original molecule
	::RDKit::ROMolOP templt_; // Template
	::RDKit::ROMolOP replace_; // Replacement
	::RDKit::RWMolOP newmol_; // Replacement

	std::map< unsigned int, AtomSubstitutionOP > by_mdx_;
	std::map< unsigned int, AtomSubstitutionOP > by_tdx_;
	std::map< unsigned int, AtomSubstitutionOP > by_rdx_;
	//std::map< unsigned int, AtomSubstitutionOP > by_ndx_;

	/// @brief What idx's correspond to dummy atoms?
	utility::vector1< unsigned int > template_dummies_;
	utility::vector1< unsigned int > replace_dummies_;

	/// @brief Storage for the individual atom substitutions
	utility::vector1< AtomSubstitutionOP > atom_substitutions_;

};

/// @brief Pick a template to use, return it and the atom pairing of template->rdmol
MoleculeSubstitutionOP
pick_template(
	::RDKit::ROMolOP rdmol,
	utility::vector1< ::RDKit::ROMolOP > & templates,
	bool dummy_only = true // Only match for templates where attached (non-hydrogen) atoms are connected through dummies.
);

/// @brief Is the replacement molecule template compatible with the molecule and template in template_frag_infos
/// If the replacement is appropriate, the return value will be a *new* MoleculeSubstitution which
/// has the appropriate replacement entries set. If the replacement doesn't work, it will return a null pointer
/// @details current_molsub should have the original and template molecules set.
MoleculeSubstitutionOP
test_replacement(
	MoleculeSubstitutionOP current_molsub,
	::RDKit::ROMolOP possible_replacement,
	core::Real dist_threshold
);

/// @brief Given the molecule and template info in template_frag_infos,
/// pick a possible replacement molecule from possible_replacements
/// @details current_molsub should have the original and template molecules set.
/// If multiple replacements are possible, one will be picked randomly,
/// based on the value of weighting_property in the replacement ROMol properties
MoleculeSubstitutionOP
pick_replacement(
	MoleculeSubstitutionOP current_molsub,
	utility::vector1< ::RDKit::ROMolOP > & possible_replacements,
	core::Real distance_threshold,
	std::string weighting_property = ""
);

/// @brief Copies all atom and bonds attached to start in the source molecule to the new molecule.
/// If skip is a valid index, it will be supressed in the source molecule, and will not be counted for copying or for attachment purposes
/// The transform is applied to the coordinates of the atom to get the coordinates of the new atom.
void
copy_attached_atoms( MoleculeSubstitution & molsub,
	MoleculeSelection source,
	::RDGeom::Transform3D const & transform,
	unsigned int start,
	unsigned int skip = AtomSubstitution::invalid_index );

}
}

#endif // Include guard
