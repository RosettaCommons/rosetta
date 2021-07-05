// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/chemical/rdkit/RDMolToRestype.hh
///
/// @brief A class for creating a residuetype based on a RDKit molecule
///
/// @details This class takes a RDKit molecule and converts it into a ResidueType.
///
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_core_chemical_rdkit_RDMolToRestype_hh
#define INCLUDED_core_chemical_rdkit_RDMolToRestype_hh

#include <core/chemical/MutableResidueType.fwd.hh>
#include <core/chemical/AtomRefMapping.hh>


#include <rdkit/GraphMol/ROMol.h>

namespace core {
namespace chemical {
namespace rdkit {

class RDMolToRestype {
public:

	RDMolToRestype( ::RDKit::ROMol const & rdmol);

	/// @brief Convert the stored molecule into a ResidueType.
	MutableResidueTypeOP generate_restype(VDIndexMapping const & mapping = {});

	/// @brief Convert the stored molecule into a ResidueType,
	/// extracting additional information from the provided restype
	/// mapping is a mapping of the residue type atoms onto the rdmol atoms
	///
	/// Currently, it just pulls the following:
	/// * The type sets
	/// * name/name1/name3 and related
	/// * atom name information for atoms with valid mappings.
	MutableResidueTypeOP generate_restype(MutableResidueType const & orig_restype, VDIndexMapping const & mapping);

	MutableResidueTypeOP generate_restype(IndexNameMapping const & mapping);
	MutableResidueTypeOP generate_restype(MutableResidueType const & orig_restype, IndexNameMapping const & mapping);

	/// @brief Get how the most recently created ResidueType corresponds to the underlying fragment.
	IndexVDMapping const &
	index_to_vd() const { return index_to_vd_; }

	/// @brief Which atom in the fragment to use as the neighbor atom when the a restype is generated.
	void set_nbr( core::Size nbr ) { nbr_ = nbr; }

private:

	RDMolToRestype();

	/// @brief Which index in the fragment is used for the neighbor atom.
	/// utility::get_undefined_size() means autodetermine.
	core::Size nbr_;

	/// @brief How the fragments indicies map to the most recently created ResidueType
	IndexVDMapping index_to_vd_;

	::RDKit::ROMol rdmol_;

};


}
}
}


#endif
