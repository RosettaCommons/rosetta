// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/ligand_docking/scoring_grid/GridManager.cc
/// @author Sam DeLuca
/// @author Rocco Moretti (rmorettiase@gmail.com)

#include <protocols/qsar/scoring_grid/GridManager.hh>
#include <protocols/qsar/scoring_grid/SingleGrid.hh>
#include <protocols/qsar/scoring_grid/GridFactory.hh>
#include <protocols/qsar/scoring_grid/GridSet.hh>
#include <protocols/qsar/scoring_grid/ScoreNormalization.hh>
#include <protocols/qsar/qsarMap.hh>

#include <core/pose/util.hh>
#include <core/pose/chains_util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/UltraLightResidue.hh>
#include <basic/Tracer.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/qsar.OptionKeys.gen.hh>

// Utility headers
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/exit.hh>
#include <utility/tools/make_vector.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/json_spirit/json_spirit_writer.h>
#include <utility/json_spirit/json_spirit_reader.h>
#include <utility/file/file_sys_util.hh>

#include <numeric/xyzVector.string.hh>

//STL headers
#include <iostream>
#include <fstream>
#include <map>

// Boost headers

namespace protocols {
namespace qsar {
namespace scoring_grid {

static basic::Tracer TR( "protocols.qsar.scoring_grid.GridManager" );

GridManager::GridManager() = default;

GridSetCOP
GridManager::get_grids(
	GridSet const & prototype,
	core::pose::Pose const & pose,
	core::Vector const & center,
	std::string const & chain,
	bool exclude,
	core::Real fuzz_factor )
{
	utility::vector1< std::string > chains(1, chain );
	return get_grids( prototype, pose, center, chains, exclude, fuzz_factor );
}

GridSetCOP
GridManager::get_grids(
	GridSet const & prototype,
	core::pose::Pose const & pose,
	core::Vector const & center,
	char chain,
	bool exclude,
	core::Real fuzz_factor )
{
	utility::vector1< std::string > chains(1, utility::to_string( chain ) );
	return get_grids( prototype, pose, center, chains, exclude, fuzz_factor );
}

GridSetCOP
GridManager::get_grids(
	GridSet const & prototype,
	core::pose::Pose const & pose,
	core::Vector const & center,
	utility::vector1< std::string > const & chains,
	bool exclude /*=true*/,
	core::Real fuzz_factor /*=0.0*/ )
{
	std::string const chain_hash = compute_hash(prototype, pose, center, chains, exclude );

	GridSetCOP grid_set = find_grid_set( chain_hash );
	if ( grid_set ) { return grid_set; }

	if ( fuzz_factor > 0.0 ) {
		grid_set = fuzzy_find_grid_set( prototype, pose, center, chains, exclude, fuzz_factor );
		if ( grid_set ) { return grid_set; }
	}

	bool grid_directory_active = basic::options::option[basic::options::OptionKeys::qsar::grid_dir].user();

	if ( !grid_directory_active ) {
		TR.Warning << "option -qsar:grid_dir is not set.  Use this flag to specify a directory to store scoring grids.  This will save you a huge amount of time" <<std::endl;
	} else if ( !utility::file::file_exists(basic::options::option[basic::options::OptionKeys::qsar::grid_dir]()) ) {
		utility_exit_with_message(basic::options::option[basic::options::OptionKeys::qsar::grid_dir]()+" does not exist. specify a valid path with -qsar:grid_dir");
	}

	//Try to read it off the disk
	if ( grid_directory_active ) {
		//files are in the format grid_directory/hash.json.gz
		std::string directory_path(basic::options::option[basic::options::OptionKeys::qsar::grid_dir]());
		utility::io::izstream grid_file(directory_path+"/"+chain_hash+".json.gz");
		if ( grid_file ) {
			utility::json_spirit::mValue gridmap_data;
			utility::json_spirit::read(grid_file,gridmap_data);
			GridSetOP new_grid_set( new GridSet );
			new_grid_set->deserialize(gridmap_data.get_array());
			TR << "successfully read grids from the disk for conformation matching hash" << chain_hash <<std::endl;
			insert_into_cache( chain_hash, center, new_grid_set );
			return new_grid_set;
		}
	}

	// This is a new conformation
	TR << "No conformation matching hash: " << chain_hash << " Updating grid and adding it to the cache" <<std::endl;

	// Make a copy of the prototype, and then reinitialize for the current system.
	GridSetOP new_grid_set( prototype.clone() );
	new_grid_set->reinitialize( pose, center );

	insert_into_cache( chain_hash, center, new_grid_set );

	if ( grid_directory_active ) {
		//if we just made a grid, we should write it to the disk for safekeeping.
		std::string directory_path(basic::options::option[basic::options::OptionKeys::qsar::grid_dir]());
		std::string temp_path(directory_path+"/"+chain_hash+".inprogress");
		if ( !utility::file::file_exists(temp_path) ) {  //If the inprogress file is there something else is busy writing
			utility::io::ozstream progress_file(temp_path);
			progress_file << "temp" <<std::endl;
			progress_file.close();

			utility::io::ozstream grid_file(directory_path+"/"+chain_hash+".json.gz");

			grid_file << utility::json_spirit::write(new_grid_set->serialize()) << std::endl;
			grid_file.close();
			utility::file::file_delete(temp_path);

			TR << "wrote grid matching hash: " << chain_hash << " to disk" <<std::endl;
		}
	}

	return new_grid_set;
}

std::string
GridManager::compute_hash(
	GridSet const & prototype,
	core::pose::Pose const & pose,
	core::Vector const & center,
	utility::vector1< std::string > const & chains,
	bool exclude
) const {
	// Include the fingerprint of the grid setup and the center of the grid in with the pose hash.
	// We don't need to be all that accurate with the position of the center. If we're off by 0.1 Ang, it shouldn't make too much difference.
	if ( ! exclude ) {
		return core::pose::get_sha1_hash_from_chains(chains, pose, numeric::truncate_and_serialize_xyz_vector(center,1) + prototype.hash_fingerprint() );
	} else {
		return core::pose::get_sha1_hash_excluding_chains(chains, pose, numeric::truncate_and_serialize_xyz_vector(center,1) + prototype.hash_fingerprint() );
	}
}

GridSetCOP
GridManager::find_grid_set(std::string const & chain_hash ) const {
	std::map<std::string,GridSetCOP>::const_iterator grid_cache_entry(grid_set_cache_.find(chain_hash));

	if ( grid_cache_entry != grid_set_cache_.end() ) { //we've already seen this conformation, load the associated grid out of the map
		TR << "Found a conformation matching hash: " << chain_hash << " Loading from grid cache" <<std::endl;
		return grid_cache_entry->second;
	}
	return nullptr;
}

GridSetCOP
GridManager::fuzzy_find_grid_set(
	GridSet const & prototype,
	core::pose::Pose const & pose,
	core::Vector const & center,
	utility::vector1< std::string > const & chains,
	bool exclude,
	core::Real fuzz_factor ) const
{
	// See if any of the current centers is close enough.
	for ( core::Vector const & new_center: centers_ ) {
		if ( new_center.distance( center ) > fuzz_factor ) { continue; }

		// (We can't necessarily trust that the hash will be the same, even with the new center matching.)
		std::string const chain_hash = compute_hash(prototype, pose, new_center, chains, exclude );
		GridSetCOP grid_cache = find_grid_set( chain_hash );
		if ( grid_cache ) { return grid_cache; }
	}

	return nullptr; // didn't find.
}

/// @brief Insert the given GridSet into the grid_set_cache under the given index value.
void
GridManager::insert_into_cache( std::string const & hash_val, core::Vector const & center, GridSetOP const & grid_set ) {
	if ( grid_set_cache_.count( hash_val ) != 0 ) {
		utility_exit_with_message( "Error: Attempting to insert a new GridSet into the cache over an existing one.");
	}
	if ( basic::options::option[basic::options::OptionKeys::qsar::max_grid_cache_size].active() &&
			grid_set_cache_.size()  >= core::Size(basic::options::option[basic::options::OptionKeys::qsar::max_grid_cache_size]() ) ) {
		TR << "Grid cache exceeds max_cache_size of " << basic::options::option[basic::options::OptionKeys::qsar::max_grid_cache_size]() << ", clearing old scoring grids to save memory." <<std::endl;
		// Can/Should we be more intelligent about this than simply blowing away the cache?
		grid_set_cache_.clear();
		centers_.clear();
	}

	grid_set_cache_[ hash_val ] = grid_set;
	centers_.push_back( center );
	// We assume that there aren't any further modifications of the cached grid_set
}

}
}
}
