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
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <sstream>
#include <sys/stat.h>

#include <utility/vector1.hh>


namespace protocols {
namespace motifs {

namespace {

char const ligand_cache_magic[8] = { 'L', 'D', 'M', 'I', 'D', 'X', '2', '\0' };

template< class T >
bool read_binary(std::istream & input, T & value)
{
	input.read(reinterpret_cast< char * >(&value), sizeof(T));
	return static_cast< bool >(input);
}

bool read_cache_string(std::istream & input, std::string & value)
{
	std::uint16_t length = 0;
	if ( !read_binary(input, length) ) return false;
	value.resize(length);
	if ( length != 0 ) input.read(&value[0], length);
	return static_cast< bool >(input);
}

struct LigandCacheGroup {
	std::uint64_t offset;
	std::uint64_t size;
};

} // anonymous namespace

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
MotifLibrary::check_ctor_helper(
	std::istream & motif_info,
	bool check_for_bad_motifs,
	std::set< std::tuple< std::string, std::string, std::string > > const & ligand_atom_type_triplets
)
{
	std::string key_in;
	while ( motif_info >> key_in ) {
		if ( key_in != "SINGLE" ) continue;

		std::string res1, r1atom1, r1atom2, r1atom3;
		std::string res2, r2atom1, r2atom2, r2atom3;
		motif_info >> res1 >> r1atom1 >> r1atom2 >> r1atom3
			>> res2 >> r2atom1 >> r2atom2 >> r2atom3;

		bool const keep_motif(ligand_atom_type_triplets.empty() ||
			ligand_atom_type_triplets.count(std::make_tuple(r2atom1, r2atom2, r2atom3)) != 0);
		if ( !keep_motif ) {
			motif_info.ignore(std::numeric_limits< std::streamsize >::max(), '\n');
			motif_info.ignore(std::numeric_limits< std::streamsize >::max(), '\n');
			continue;
		}

		std::string rest_of_line;
		std::getline(motif_info, rest_of_line);
		std::string::size_type const remark_index(rest_of_line.find("REMARK "));

		core::kinematics::Jump motif_jump;
		motif_info >> motif_jump;
		bool const bad_motif(motif_info.rdstate() != 0);
		motif_info.clear();
		motif_info.setstate(std::ios_base::goodbit);
		if ( check_for_bad_motifs && bad_motif ) continue;

		SingleMotifOP new_motif(new SingleMotif(
			res1, r1atom1, r1atom2, r1atom3,
			r2atom1, r2atom2, r2atom3, motif_jump));
		if ( remark_index != std::string::npos ) {
			new_motif->store_remark(rest_of_line.substr(remark_index + 7));
		}
		library_.push_back(new_motif);
	}
}

MotifLibrary::MotifLibrary(
	std::istream & motif_info,
	bool check_for_bad_motifs,
	std::set< std::tuple< std::string, std::string, std::string > > const & ligand_atom_type_triplets
)
{
	check_ctor_helper(motif_info, check_for_bad_motifs, ligand_atom_type_triplets);
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

void
MotifLibrary::add_ligand_from_file(
	std::string const & motif_filename,
	bool check_for_bad_motifs,
	std::set< std::tuple< std::string, std::string, std::string > > const & ligand_atom_type_triplets
)
{
	char const * cache_environment = std::getenv("ROSETTA_LIGAND_MOTIF_CACHE");
	std::string const cache_filename = cache_environment && cache_environment[0] != '\0'
		? std::string(cache_environment) : motif_filename + ".ldcache";
	if ( utility::file::file_exists(cache_filename)
		&& load_ligand_index_cache(motif_filename, cache_filename, check_for_bad_motifs, ligand_atom_type_triplets) ) {
		return;
	}

	std::ifstream motif_file;
	motif_file.open( motif_filename.c_str() );
	if ( !motif_file ) {
		std::cout << "ERROR:  No motif file " << motif_filename << " - FAILING!\n";
		return;
	}

	check_ctor_helper( motif_file, check_for_bad_motifs, ligand_atom_type_triplets );
}

bool
MotifLibrary::load_ligand_index_cache(
	std::string const & motif_filename,
	std::string const & cache_filename,
	bool check_for_bad_motifs,
	std::set< std::tuple< std::string, std::string, std::string > > const & ligand_atom_type_triplets
)
{
	struct stat source_status;
	if ( ::stat(motif_filename.c_str(), &source_status) != 0 ) return false;

	std::ifstream cache(cache_filename.c_str(), std::ios::in | std::ios::binary);
	if ( !cache ) return false;

	char magic[sizeof(ligand_cache_magic)];
	cache.read(magic, sizeof(magic));
	std::uint64_t source_size = 0;
	std::int64_t source_mtime = 0;
	std::uint32_t group_count = 0;
	if ( !cache || std::memcmp(magic, ligand_cache_magic, sizeof(magic)) != 0
		|| !read_binary(cache, source_size) || !read_binary(cache, source_mtime)
		|| !read_binary(cache, group_count)
		|| source_size != static_cast< std::uint64_t >(source_status.st_size)
		|| source_mtime != static_cast< std::int64_t >(source_status.st_mtime) ) {
		return false;
	}

	utility::vector1< LigandCacheGroup > selected_groups;
	for ( std::uint32_t group_index = 0; group_index < group_count; ++group_index ) {
		std::string atom1, atom2, atom3;
		std::uint64_t offset = 0;
		std::uint64_t size = 0;
		if ( !read_cache_string(cache, atom1) || !read_cache_string(cache, atom2)
			|| !read_cache_string(cache, atom3) || !read_binary(cache, offset)
			|| !read_binary(cache, size) ) {
			return false;
		}
		auto const key = std::make_tuple(atom1, atom2, atom3);
		if ( ligand_atom_type_triplets.empty() || ligand_atom_type_triplets.count(key) != 0 ) {
			selected_groups.push_back({offset, size});
		}
	}

	std::set< std::tuple< std::string, std::string, std::string > > const accept_all;
	MotifLibrary cached_library;
	for ( auto const & group : selected_groups ) {
		cache.seekg(static_cast< std::streamoff >(group.offset));
		if ( !cache ) return false;
		std::string group_records(static_cast< std::size_t >(group.size), '\0');
		if ( group.size != 0 ) cache.read(&group_records[0], static_cast< std::streamsize >(group.size));
		if ( !cache ) return false;
		std::istringstream group_stream(group_records);
		cached_library.check_ctor_helper(group_stream, check_for_bad_motifs, accept_all);
	}

	for ( auto const & motif : cached_library.library() ) library_.push_back(motif);
	std::cout << "Loaded " << cached_library.nmotifs() << " motifs from indexed cache " << cache_filename << '\n';
	return true;
}

} // namespace motifs
} // namespace protocols
