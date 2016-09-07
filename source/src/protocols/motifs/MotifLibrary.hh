// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file MotifLibrary.hh
/// @brief class declaration for sets of interaction motifs between residues
/// @author havranek, sthyme (sthyme@gmail.com)

#ifndef INCLUDED_protocols_motifs_MotifLibrary_hh
#define INCLUDED_protocols_motifs_MotifLibrary_hh

// Unit Headers
#include <protocols/motifs/MotifLibrary.fwd.hh>

// Package Headers
#include <protocols/motifs/Motif.fwd.hh>

// Project Headers
#include <core/types.hh>

// Utility Headers
#include <utility/file/FileName.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// C++ Headers
#include <iosfwd>

//Auto Headers
namespace protocols {
namespace motifs {

class MotifLibrary : public utility::pointer::ReferenceCount
{

public:

	typedef utility::vector1< utility::file::FileName > const FileNames;

	// Constructor
	MotifLibrary();

	// Destructor
	~MotifLibrary() override;

	// Constructor for loading motifs from motif PDBs
	MotifLibrary(
		utility::vector1< utility::file::FileName > const & motif_filenames
	);

	// Constructor for loading motifs from a file with coordinates
	MotifLibrary(
		std::istream & motif_info
	);

	// Add motif to library
	void
	add_to_library(
		Motif const & new_motif
	);

	// Constructor for loading motifs from a file with coordinates for ligands
	MotifLibrary(
		std::istream & motif_info, core::Size ligand_marker
	);

	// Add ligand motif from a PDB file
	void
	add_ligand_from_file(
		std::string const & motif_filename
	);

	// Add motif from a PDB file
	void
	add_from_file(
		std::string const & motif_filename
	);

	// Number of motifs
	core::Size
	nmotifs();

	// Iterators
	MotifCOPs::const_iterator begin() { return library_.begin(); }
	MotifCOPs::const_iterator end() { return library_.end(); }

	MotifCOPs::const_iterator begin() const { return library_.begin(); }
	MotifCOPs::const_iterator end() const { return library_.end(); }

	// Accessors
	MotifCOPs const & library() const { return library_; }

	// Overloaded operator for output
	friend std::ostream & operator <<(
		std::ostream & os,
		MotifLibrary const & mot_lib
	);

private:
	MotifCOPs library_;

};

} // namespace motifs
} // namespace protocols

#endif // INCLUDED_protocols_motifs_MotifLibrary
