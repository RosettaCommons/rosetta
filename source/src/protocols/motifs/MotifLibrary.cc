// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/motifs/MotifLibrary.cc
/// @brief Implementation of interaction motifs
/// @author havranek, sthyme (sthyme@gmail.com)

// Unit Headers
#include <protocols/motifs/MotifLibrary.hh>

// Package Headers
#include <protocols/motifs/motif_utils.hh>
#include <protocols/motifs/Motif.hh>
#include <protocols/motifs/SingleMotif.hh>

// Project Headers

// Utility Headers
#include <utility/file/file_sys_util.hh>
#include <utility/file/FileName.hh>

// C++ Headers
#include <fstream>
#include <iostream>

#include <utility/vector1.hh>


namespace protocols {
namespace motifs {

MotifLibrary::MotifLibrary() = default;

MotifLibrary::~MotifLibrary() = default;

MotifLibrary::MotifLibrary(
	FileNames & motif_filenames
)
{
	for ( auto const & motif_filename : motif_filenames ) {
		if ( !utility::file::file_exists( motif_filename ) ) {
			continue;
		}
		SingleMotifOP new_motif = single_motif_from_filename( motif_filename );
		add_to_library( *new_motif );
	}
}

MotifLibrary::MotifLibrary(
	std::istream & motif_info
)
{
	std::string key_in;
	while ( motif_info >> key_in ) {
		if ( key_in == "SINGLE" ) {
			SingleMotifOP new_motif = single_motif_from_stream( motif_info );
			add_to_library( *new_motif );
		} else {
			std::cout << "ERROR - BAD MOTIF KEY " << key_in << "\n";
		}
	}
}

MotifLibrary::MotifLibrary(
	std::istream & motif_info, core::Size
)
{
	//std::cout << "In MotifLibrary.cc, in LigandMotifLibrary istream function" << std::endl;
	std::string key_in;
	while ( motif_info >> key_in ) {
		if ( key_in == "SINGLE" ) {
			//std::cout << "In MotifLibrary.cc, about to make single motif OP" << std::endl;

			SingleMotifOP new_motif = single_ligand_motif_from_stream( motif_info );
			add_to_library( *new_motif );
		} else {
			//std::cout << "ERROR - BAD MOTIF KEY " << key_in << "\n";
		}
	}
}

void
MotifLibrary::check_ctor_helper(std::istream & motif_info, bool check_for_bad_motifs, utility::vector1< std::string > const & ligand_atom_names)
{
	//std::cout << "In MotifLibrary.cc, in LigandMotifLibrary istream function" << std::endl;
	std::string key_in;
	while ( motif_info >> key_in ) {
		if ( key_in == "SINGLE" ) {
			//std::cout << "In MotifLibrary.cc, about to make single motif OP" << std::endl;

			//test to confirm that >> overload works
			//std::cout << "Trying >> overload" << std::endl;
			//SingleMotifOP new_motif;
			//motif_info >> new_motif;

			//original line for backup
			SingleMotifOP new_motif = single_ligand_motif_from_stream( motif_info, check_for_bad_motifs );

			//handling for if using an empty vector of ligand_atom_names
			if ( ligand_atom_names.size() == 0 ) {
				add_to_library( *new_motif );

				//continue, no reason to look down the rest of the function
				continue;
			}

			//run through the motif, and determine if it is relevant to the ligand
			//preliminary filter that will omit a motif if at least one atom in the motif ligand is not present in the ligand
			bool keep_motif = false;
			bool atom1_good = false;
			bool atom2_good = false;
			bool atom3_good = false;
			std::string motif_atom_1(new_motif->res2_atom1_name());
			std::string motif_atom_2(new_motif->res2_atom2_name());
			std::string motif_atom_3(new_motif->res2_atom3_name());

			//std::cout << "Motif atoms: " << motif_atom_1 << " " << motif_atom_2 << " " << motif_atom_3 << std::endl;

			for ( core::Size this_atom = 1; this_atom <= ligand_atom_names.size(); ++this_atom ) {
				//std::cout << "Comparing test: " << ligand_atom_names[this_atom] << " to " << motif_atom_1 << std::endl;
				if ( motif_atom_1 == ligand_atom_names[this_atom] ) {
					//std::cout << "Match!" << std::endl;
					atom1_good = true;
				}
				//std::cout << "Comparing test: " << ligand_atom_names[this_atom] << " to " << motif_atom_2 << std::endl;
				if ( motif_atom_2 == ligand_atom_names[this_atom] ) {
					//std::cout << "Match!" << std::endl;
					atom2_good = true;
				}
				//std::cout << "Comparing test: " << ligand_atom_names[this_atom] << " to " << motif_atom_3 << std::endl;
				if ( motif_atom_3 == ligand_atom_names[this_atom] ) {
					//std::cout << "Match!" << std::endl;
					atom3_good = true;
				}
			}

			if ( atom1_good && atom2_good && atom3_good ) {
				keep_motif = true;
			}

			if ( keep_motif ) {
				add_to_library( *new_motif );
			}
		}
	}
}

// @brief
MotifLibrary::MotifLibrary(
	std::istream & motif_info, bool check_for_bad_motifs, utility::vector1< std::string > const & ligand_atom_names
)
{
	check_ctor_helper(motif_info, check_for_bad_motifs, ligand_atom_names);
}

void
MotifLibrary::add_to_library( Motif const & add_me )
{
	library_.push_back( add_me.clone() );
}

core::Size
MotifLibrary::nmotifs()
{
	return library_.size();
}

void
MotifLibrary::add_from_file( std::string const & motif_filename )
{
	// Try to open the file
	std::ifstream motif_file;
	motif_file.open( motif_filename.c_str() );
	if ( !motif_file ) {
		std::cout << "ERROR:  No motif file " << motif_filename << " - FAILING!\n";
		return;
	}

	// Attempt to read in motifs until exhausted
	MotifLibrary new_library( motif_file );

	// Add to this library
	for ( auto const & pmot : new_library ) {
		add_to_library( *pmot );
	}
}

std::ostream & operator <<(
	std::ostream & os, MotifLibrary const & mlib
)
{
	for ( auto const & pmot : mlib ) {
		os << pmot->print();
	}
	return os;
}

void
MotifLibrary::add_ligand_from_file( std::string const & motif_filename )
{
	std::ifstream motif_file;
	motif_file.open( motif_filename.c_str() );
	if ( !motif_file ) {
		std::cout << "ERROR:  No motif file " << motif_filename << " - FAILING!\n";
		return;
	}

	// Attempt to read in motifs until exhausted
	core::Size ligand_marker = 1;
	MotifLibrary new_library( motif_file, ligand_marker );

	// Add to this library
	for ( auto const & pmot : new_library ) {
		add_to_library( *pmot );
	}
}

void
MotifLibrary::add_ligand_from_file( std::string const & motif_filename, bool check_for_bad_motifs, utility::vector1< std::string > const & ligand_atom_names )
{
	std::ifstream motif_file;
	motif_file.open( motif_filename.c_str() );
	if ( !motif_file ) {
		std::cout << "ERROR:  No motif file " << motif_filename << " - FAILING!\n";
		return;
	}

	// Attempt to read in motifs until exhausted
	MotifLibrary new_library( motif_file, check_for_bad_motifs, ligand_atom_names );

	// Add to this library
	for ( auto const & pmot : new_library ) {
		add_to_library( *pmot );
	}
}

} // namespace motifs
} // namespace protocols
