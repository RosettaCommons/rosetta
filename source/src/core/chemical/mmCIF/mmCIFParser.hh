// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//////////////////////////////////////////////////////////////////////
///
/// @brief
/// Class to read an mmCIF file, which defines chemical compositions for
/// residues. Most notably used for small molecuels
///
/// @file   src/core/chemical/sdf/mmCIFParser.hh
///
/// @details
/// Parse an mmCIF file and create an MolFileIOMoleculeOP object from the parsed data
/// The reason why we create a MolFileIOMoleculeOP object is that there is
/// machinary available to take a MolFileIOMoleculeOP and convert it into
/// a residuetype. The function get_molfile_molecule() aggregates data
/// from tables in the mmCIF parser and adds them to the MolFileIOMoleculeOP.
/// Downstream process can convert the molfile into a residuetype.
///
/// @author
/// Steven Combs (steven.combs1@gmail.com)
///
///
/////////////////////////////////////////////////////////////////////////

/// @author Steven Combs

#ifndef INCLUDED_core_chemical_mmCIF_mmCIFParser_hh
#define INCLUDED_core_chemical_mmCIF_mmCIFParser_hh

#include <core/chemical/mmCIF/mmCIFParser.fwd.hh>

#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <map>
#include <core/types.hh>
#include <cifparse/CifFile.h>
#include <core/chemical/sdf/MolFileIOData.fwd.hh>

namespace core {
namespace chemical {
namespace mmCIF {

class mmCIFParser : public utility::pointer::ReferenceCount
{

public:

	mmCIFParser();
	virtual ~mmCIFParser(){};

	/// @brief parse the mmCIF file. This is generally used with the command line option
	/// -extra_res_mmCIF_fa. The mmCIF can contain multiple blocks that contain different
	/// residue paramaters. This will read each block and pushback to the MolFileIO object
	utility::vector1< sdf::MolFileIOMoleculeOP> parse( std::string const & filename);

	/// @brief parse a concatenated string (line) and create a MolFileIO object from the
	/// stream. Specificaly pull out the block from the pdb_id. This is used in conjection
	/// with streamed mmCIF files, ie, the lazy loader.
	sdf::MolFileIOMoleculeOP parse(  std::string const &lines, std::string const &pdb_id);

	/// @brief -the secret sauce. Pull tables form the block segment to construct parameters
	/// return a MolFileIOMolecule
	sdf::MolFileIOMoleculeOP
	get_molfile_molecule(Block& block);

private:
	/// @brief When you are lazy, create a map of strings and sizes to convert the string
	/// for bonds found in the mmCIF file to the sdf version
	std::map< std::string, core::Size> bond_string_to_sdf_size_;

};



}
}
}
#endif
