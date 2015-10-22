// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/chemical/rotamers/PDBRotamerLibrarySpecification.hh
/// @brief  The PDBRotamerLibrarySpecification class says to build a rotamer library from a PDB file
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_core_chemical_rotamers_PDBRotamerLibrarySpecification_HH
#define INCLUDED_core_chemical_rotamers_PDBRotamerLibrarySpecification_HH

// Unit headers
#include <core/chemical/rotamers/PDBRotamerLibrarySpecification.fwd.hh>
#include <core/chemical/rotamers/RotamerLibrarySpecification.hh>

// Package headers
#include <core/chemical/ResidueType.fwd.hh>

// C++ headers
#include <istream>

namespace core {
namespace chemical {
namespace rotamers {

class PDBRotamerLibrarySpecification : public RotamerLibrarySpecification {
public:
	PDBRotamerLibrarySpecification();
	PDBRotamerLibrarySpecification(std::string library_filename);
	PDBRotamerLibrarySpecification(std::istream & input);
	virtual ~PDBRotamerLibrarySpecification();

	std::string const &
	pdb_rotamers_file() const { return pdb_rotamers_file_; }

	void
	pdb_rotamers_file( std::string const & filename){ pdb_rotamers_file_ = filename; }

	virtual
	std::string
	keyname() const;

	virtual
	std::string
	cache_tag(core::chemical::ResidueType const &) const { return pdb_rotamers_file_; }

	static
	std::string
	library_name();

private:

	std::string pdb_rotamers_file_;

};


} //namespace rotamers
} //namespace chemical
} //namespace core


#endif
