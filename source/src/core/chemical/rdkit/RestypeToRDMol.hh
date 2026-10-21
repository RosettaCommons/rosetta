// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file RestypeToRDMol.chh
///
/// @brief
/// A class for creating an RDKit RWMol object from a Rosetta ResidueType
///
/// @author Rocco Moretti (rmorettiase@gmail.com)


#ifndef INCLUDED_core_chemical_rdkit_RestypeToRDMol_hh
#define INCLUDED_core_chemical_rdkit_RestypeToRDMol_hh


#include <core/chemical/MutableResidueType.fwd.hh>
#include <core/chemical/AtomRefMapping.hh>

#include <rdkit/GraphMol/RWMol.h>

namespace core {
namespace chemical {
namespace rdkit {

struct RestypeToRDMolOptions{
	bool neutralize = true;
	bool keep_hydro = false;
	bool sanitize = true;
	bool noImplicitHs = false;
	bool skipHs = false;
	bool aro2double = false;
};

class RestypeToRDMol  {
public:

	/// @brief Convert ResidueType to an RDKit RWMol object
	///
	/// If neutralize is true, attempt to reprotonate into a neutral (non-formally charged) molecule.
	/// (Charges due to heavy atoms, like quaternary amines, will still be present.)
	///
	/// If keep_hydro is true, represent hydrogens as physical atoms, rather than "explicit"
	/// annotations (see below).
	///
	/// Note that due to implementation details, "neutralize = true" will not play well with "keep_hydro = true".
	///
	/// @details RDKit is "aware" of three types of hydrogens.
	///
	/// 1) Physical hydrogens: actual atoms in the atom graph and can have coordinates.
	/// 2) "Explicit" hydrogens: don't have existence in the atom graph or coordinates,
	///    but are instead represented by a field on the heavy atom ("This atom has 3 hydrogens attached to it")
	/// 3) "Implicit" hydrogens: aren't annotated anywhere, but instead are implied by
	///    the difference in the expected valence of the heavy atom and the number of
	///    valences which are currently occupied by bonds and/or "explicit" hydrogens.
	///
	/// Most of RDKit is written around the assumption of non-physical hydrogens, the exception
	/// being things like energy minimization which needs coordinates.
	/// On the other hand, "implicit" hydrogens cause issues in certain cases with kekulization.
	/// Because of this, the default conversion here is to remove the physical hydrogens, and
	/// replace them with "explicit" hydrogen annotations.
	///
	/// Much of RDKit (especially the metric calculations) assumes neutral protonation. (e.g.
	/// what you'd see in an aprotic organic solvent), which is why we neutralize the residues
	/// by default.
	RestypeToRDMol(MutableResidueType const & res, bool neutralize = true, bool keep_hydro = false,
		bool sanitize = true, bool noImplicitHs = false, bool skipHs = false );

	RestypeToRDMol(MutableResidueType const & res, RestypeToRDMolOptions const & options );

	/// @brief Return the corespondance from the input ResidueType VD to the RDKit atom index
	VDIndexMapping const &
	vd_to_index() const { return vd_to_index_; }

	IndexVDMapping const &
	index_to_vd() const { return index_to_vd_; }

	std::map<core::Size, core::Size> const &
	get_atomIndexMap_Rd2Res() const { return atomIndexMap_Rd2Res_; }

	/// @brief Return an RDKit RWMol object which represents the residue type
	::RDKit::RWMOL_SPTR
	Mol();

private:
	RestypeToRDMol(); // unimplemented - must provide a residue type

	MutableResidueType const & res_;
	RestypeToRDMolOptions options_;
	bool neutralize_;
	bool keep_hydro_;
	bool sanitize_;
	bool noImplicitHs_;
	bool skipHs_;
	bool aro2double_ = false;

	/// @brief Mapping of restype vertex descriptors to indices of the rdkit object
	// Uses utility::get_undefined_size() for non-represented RDKit indices
	VDIndexMapping vd_to_index_;
	IndexVDMapping index_to_vd_;
	std::map<core::Size, core::Size> atomIndexMap_Rd2Res_;

};

}
}
}


#endif
