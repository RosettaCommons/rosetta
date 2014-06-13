// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @author Oliver Lange

#ifndef INCLUDED_core_coarse_Translator_hh
#define INCLUDED_core_coarse_Translator_hh

// Unit headers
#include <core/coarse/Translator.fwd.hh>
// you cannot #include yourself #include <core/coarse/Translator.hh>

// Coarse headers
#include <core/coarse/Rules.hh>

// Project headers
//#include <core/chemical/ResidueTypeSet.hh>
//#include <core/chemical/ResidueTypeSet.fwd.hh>
//#include <core/chemical/AtomTypeSet.hh>
//#include <core/chemical/MMAtomTypeSet.fwd.hh>
#include <core/chemical/AA.hh>

// AUTO-REMOVED #include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/pose/Pose.fwd.hh>
#include <core/pack/dunbrack/RotamerLibrary.fwd.hh>
#ifdef WIN32
#include <core/pack/dunbrack/DunbrackRotamer.hh>
#endif
#include <core/pack/dunbrack/CoarseRotamer.fwd.hh>
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.hh>

// Utility headers
// AUTO-REMOVED #include <utility/vector1.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>



// std headers
#include <ostream>
#include <string>

//Auto Headers
#include <core/chemical/ResidueType.fwd.hh>
#include <core/pack/dunbrack/DunbrackRotamer.fwd.hh>



/* TODO:

make atoms in beads weighted
then a true centroid representation wouldn't need any change of the programming but just different rules
w=1 for CB and w=0 for all others
but all sidechain atoms in B1

*/

namespace core {
namespace coarse {

class Translator : public  utility::pointer::ReferenceCount {
	/* The responsibility of this class is to provide a
		 map of atom-names  between fine and coarse representations
		 for instance, bead1 contains fa_atoms { CB , HB1, HB2, HB3 etc}
		 FULL_ATOM contains fa_atoms{H CA O N}

		 organization:
		 beads[0] is AtomList of fa_atoms
		 beads[1..n] contains AtomLists for beads.
		 corresponding bead names are bead_names_[1..n]
		 (bead_names_[0] == "FULL_ATOM" )
	*/

public:

	//little helper struct. A BeadAtom is contained in the AtomList which describes
	// the ingredients of a bead. The name_ refers to the atom name in the fine-grain residue
	// the field weight_ controls how much the respective xyz coords, contribute to the bead-center
	struct BeadAtom {
		BeadAtom( std::string name ) : name_(name), weight_(1.0) {};
		BeadAtom( std::string name, Real weight ) : name_(name), weight_(weight) {};
		bool operator == ( BeadAtom const &other ) const { return other.name_==name_; };
		std::string name_;
		Real weight_;
	};

	//
	typedef std::vector<BeadAtom> AtomList;
	typedef std::vector<AtomList> BeadList;
	typedef std::vector<std::string> BeadNames;

public:

	/// @brief constructor
	Translator( RuleSet const &rules,
		chemical::ResidueTypeCOP fine_res,
		chemical::ResidueTypeAP coarse_res
	);

private:
	/// @brief private copy and assigment
	/// a translator is intimately connected to its coarse_res_type,
	/// (coarse_res_type points to its translator object and vice versa... better not to copy such things
	Translator( const Translator& t ) : ReferenceCount(t) {} ;
	Translator & operator= ( const Translator& ) { return *this; };

public:

	//
	void pretty_print( std::ostream &os ) const;

	/// @brief returns bead index for 'atom' in the fine-grained residue, e.g. 1 for CB, HB1, HB2, ...
	int map_atom_to_bead( std::string const atom ) const;

	/// @brief number of chi-angles in coarse residue
	int coarse_nchi() const;

	/// @brief  residue_type ID (they are the same for coarse and fine, return either)
	std::string const &name() const;

	/// @brief return a coarse residue with coordinates computed from the fine residue
	conformation::ResidueOP coarsify(const conformation::Residue &fine) const;


	/// @brief return a DunbrackRotamerLibrary, Rotamer's and their statistics are computed from the fine residues
	// DunbrackLibrary by calling coarsify (see below) (not inlined, to avoid circular includes)
	pack::dunbrack::SingleResidueRotamerLibraryCOP get_RotamerLibrary() const;

	/// @brief compute a coarse DunbrackLibrary from a fine RotamerLibrary  (not inlined, to avoid circular includes)
	pack::dunbrack::CoarseRotamerSetOP
	coarsify(
		utility::vector1< pack::dunbrack::DunbrackRotamer< pack::dunbrack::FOUR > > const & fine_rotamers
	) const;

protected:

	//helper routines for construction of Translator

	/// @brief add BeadAtom to AtomList
	void add_atom(AtomList &list, const chemical::ResidueType &res, Rule::AtomToken const &atom);
	void add_atom(AtomList &list, const chemical::ResidueType &res,  int pos);

	/// @brief add all non-assigned sidechain atoms to AtomList
	void add_remaining_sidechain(AtomList &list, const chemical::ResidueType &res);

	/// @brief add all non-assigned atoms to AtomList
	void add_all_remaining(AtomList &list, const chemical::ResidueType &res);

	/// @brief the PARAM files contain random geometries. fix them by coarsifying the fine-grained geometries
	/// automatically -- called by constructor
	void fix_coarsetype_geometry(chemical::ResidueTypeAP coarse_res_type);

private:
	/// @brief  Pointers to the connected Residue Sets... CHECK these should probably be access pointers CHECK
	chemical::ResidueTypeCOP coarse_res_type_;
	chemical::ResidueTypeCOP fine_res_type_;

	/// @brief list of beads, a bead is a AtomList --- which atoms of fine_rsd belong to bead of coarse_rsd
	BeadList beads_;

	/// @brief store atom and bead names of the coarse-grained residue,
	/// first bead is FULL_ATOM and contains the names of all atoms that remain in full.
	BeadNames bead_names_;

};


} //namespace coarse
} // namespace core

#endif
