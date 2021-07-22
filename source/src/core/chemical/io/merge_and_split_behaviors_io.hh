// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/chemical/io/merge_and_split_behaviors_io.hh
/// @brief   Database input/output function declarations for residue-splitting behaviors.
/// @author  Labonte <JWLabonte@jhu.edu>


#ifndef INCLUDED_core_chemical_io_merge_and_split_behaviors_io_HH
#define INCLUDED_core_chemical_io_merge_and_split_behaviors_io_HH

// C++ headers
#include <map>

// Utility header
#include <utility/vector1.hh>


namespace core {
namespace chemical {
namespace io {

/// @brief  Per-residue setting for the behavior of each residue from an input file.
enum merge_residue_behavior
{
	mrb_do_not_merge,
	mrb_merge_w_prev,
	mrb_merge_w_next
};


typedef std::map<
	std::string,  // key: old PDB atom name
	std::string  // value: new Rosetta atom name
	> AtomRenamingMap;


typedef std::pair< merge_residue_behavior, io::AtomRenamingMap > ResidueMergeInstructions;

typedef std::map< std::string,  // The map key is a PDB 3-letter code (old residue).
	ResidueMergeInstructions > MergeBehaviorMap;


typedef std::pair<
	utility::vector1< std::pair<  // list of new base residue types
	std::string,  // Rosetta 3-letter code
	std::string > >,  // Rosetta base name
	utility::vector1< AtomRenamingMap >  // list of atom renaming instructions
	> SplitBehaviors;

typedef std::map< std::string,  // The map key is a PDB 3-letter code (old residue).
	SplitBehaviors > SplitBehaviorsMap;


/// @brief  Return a mapping of PDB 3-letter codes to a set of merging instructions.
MergeBehaviorMap read_merge_behaviors_from_database_file( std::string const & filename );

/// @brief  Return a mapping of PDB 3-letter codes to a set of splitting instructions.
SplitBehaviorsMap read_split_behaviors_from_database_file( std::string const & filename );

}  // namespace io
}  // namespace chemical
}  // namespace core

#endif  // INCLUDED_core_chemical_io_merge_and_split_behaviors_io_HH
