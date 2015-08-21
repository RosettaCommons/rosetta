// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/core/chemical/sdf/MolFileIOReader.hh
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_core_chemical_sdf_MolFileIOReader_hh
#define INCLUDED_core_chemical_sdf_MolFileIOReader_hh

#include <core/chemical/sdf/MolFileIOReader.fwd.hh>
#include <core/chemical/sdf/MolFileIOData.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/AtomTypeSet.fwd.hh>
#include <core/chemical/ElementSet.fwd.hh>
#include <core/chemical/MMAtomTypeSet.fwd.hh>

#include <core/types.hh>

#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <numeric/xyzVector.hh>

#include <map>
#include <string>
#include <istream>

namespace core {
namespace chemical {
namespace sdf {

class MolFileIOReader : public utility::pointer::ReferenceCount
{
public:
	MolFileIOReader();
	virtual ~MolFileIOReader();

	/// @brief parse file, with the possibility of type autodetection.
	utility::vector1< MolFileIOMoleculeOP > parse_file( std::string const & filename, std::string type = "" );
	/// @brief parse file from stream, type must be specified.
	// type being passed by value is intentional - we need to lower case it.
	utility::vector1< MolFileIOMoleculeOP > parse_file( std::istream & file, std::string type );

private:
};

/// @brief Convert the vector of MolFileIOMolecules into a single residue type,
/// using multiple entries as rotamers
ResidueTypeOP convert_to_ResidueType( utility::vector1< MolFileIOMoleculeOP > molfile_data,
	std::string atom_type_tag = "fa_standard",
	std::string element_type_tag = "default",
	std::string mm_atom_type_tag = "fa_standard");

/// @brief Convert the vector of MolFileIOMolecules into a single residue type,
/// using multiple entries as rotamers
ResidueTypeOP convert_to_ResidueType( utility::vector1< MolFileIOMoleculeOP > molfile_data,
	AtomTypeSetCOP atom_types,
	ElementSetCOP element_types,
	MMAtomTypeSetCOP mm_atom_types);

/// @brief Convert the vector of MolFileIOMolecules into multiple residue types
/// If load_rotamers is false, each will be loaded as a single ResidueType
/// Otherwise, entries with the same name will be loaded as rotamers
utility::vector1< ResidueTypeOP > convert_to_ResidueTypes( utility::vector1< MolFileIOMoleculeOP > molfile_data,
	bool load_rotamers = true,
	std::string atom_type_tag = "fa_standard",
	std::string element_type_tag = "default",
	std::string mm_atom_type_tag = "fa_standard");

/// @brief Convert the vector of MolFileIOMolecules into multiple residue types
/// If load_rotamers is false, each will be loaded as a single ResidueType
/// Otherwise, entries with the same name will be loaded as rotamers
utility::vector1< ResidueTypeOP > convert_to_ResidueTypes( utility::vector1< MolFileIOMoleculeOP > molfile_data,
	bool load_rotamers,
	AtomTypeSetCOP atom_types,
	ElementSetCOP element_types,
	MMAtomTypeSetCOP mm_atom_types);


}
}
}

#endif
