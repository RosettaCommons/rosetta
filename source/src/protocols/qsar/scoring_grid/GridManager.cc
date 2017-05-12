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

#include <protocols/qsar/scoring_grid/GridManager.hh>
#include <protocols/qsar/scoring_grid/SingleGrid.hh>
#include <protocols/qsar/scoring_grid/GridFactory.hh>
#include <protocols/qsar/scoring_grid/ScoreNormalization.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/qsar/qsarMap.hh>

#include <core/pose/util.hh>
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
#include <utility/tag/Tag.hh>
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
#include <boost/bind.hpp>
#include <boost/function.hpp>

namespace protocols {
namespace qsar {
namespace scoring_grid {

static THREAD_LOCAL basic::Tracer TR( "protocols.qsar.scoring_grid.GridManager" );

void GridManager::reset()
{
	grid_map_.clear();
	score_map_.clear();
	last_tag_ = "";
	width_ = 40;
	resolution_ = 0.25;
	qsar_map_.reset();
	initialized_ = false;
	chain_ = 'X';
	norm_function_.reset();
}

GridManager::GridManager() :
	last_tag_(""),
	width_(40),
	resolution_(0.25),
	qsar_map_(/* 0 */),
	initialized_(false),
	chain_('X'),
	norm_function_(/* 0 */)
{
	grid_map_.clear();
	score_map_.clear();
}

void GridManager::set_width(core::Real width)
{
	width_ = width;
}

void GridManager::set_resolution(core::Real resolution)
{
	resolution_=resolution;
}

void GridManager::set_normalization_function(std::string norm_function_name)
{
	norm_function_ = get_score_normalization_function(norm_function_name);
}

void GridManager::set_chain(char chain)
{
	chain_ = chain;
}

void GridManager::set_qsar_map(qsarMapOP qsar_map)
{
	qsar_map_ = qsar_map;
}

bool GridManager::is_qsar_map_attached()
{
	if ( qsar_map_ != 0 ) {
		return true;
	} else {
		return false;
	}
}

void GridManager::make_new_grid(utility::tag::TagCOP tag)
{

	std::string name= tag->getOption< std::string >( "grid_name" );

	//Sometimes creating grids can be time consuming
	//so check before we do something that calls a grid constructor.
	//Access the various options anyways so we don't die for unaccessed options.
	core::Real weight = tag->getOption<core::Real>("weight");
	if ( grid_map_.find(name) == grid_map_.end() ) {
		grid_weights_.insert(std::make_pair(name,weight));
		TR.Debug << "Making new grid " << name << std::endl;
		GridBaseOP new_grid(GridFactory::get_instance()->new_grid(tag));
		insert_grid(name, new_grid);
	}


}

void GridManager::insert_grid(std::string const & name, GridBaseOP const grid)
{
	grid_map_[name] = grid;
}

bool GridManager::is_normalization_enabled()
{
	if ( norm_function_ ) {
		return true;
	} else {
		return false;
	}
}

GridBaseCOP GridManager::get_grid(std::string const & grid_type)
{
	return grid_map_.find(grid_type)->second;
}

utility::vector1<std::string> GridManager::get_grid_names()
{
	utility::vector1<std::string> grid_names;

	for ( GridMap::value_type const & grid_map_entry: grid_map_ ) {
		grid_names.push_back(grid_map_entry.first);
	}
	return grid_names;
}

core::Real GridManager::ideal_score(utility::vector1<core::conformation::UltraLightResidue> & residues)
{

	core::Real score=0.0;

	//Does not use weighted average because the maximum score of each residue already scales with atom size
	//Hence, the contribution to total_score from larger ligands is already greater

	for ( core::conformation::UltraLightResidue const & residue: residues ) {
		score += ideal_score(residue);
	}

	score = (score)/(residues.size());
	return score;

}

core::Real GridManager::ideal_score(core::conformation::UltraLightResidue const & residue)
{

	core::Real score = 0;
	score= (core::Real)(-1) * residue.natoms();

	return score;

}

core::Real GridManager::total_score(utility::vector1<core::conformation::UltraLightResidue> & residues)
{
	core::Real score=0.0;

	//Does not use weighted average because the maximum score of each residue already scales with atom size
	//Hence, the contribution to total_score from larger ligands is already greater

	for ( core::conformation::UltraLightResidue const & residue: residues ) {
		score += total_score(residue);
	}

	score = (score)/(residues.size());
	return score;

}


core::Real GridManager::total_score(core::conformation::UltraLightResidue const & residue)
{
	score_map_.clear();

	core::Real total_score =0.0;
	const core::Real max_score = 9999.0;
	for ( GridMap::value_type const & grid_map_entry: grid_map_ ) {
		core::Real component_score =0;
		GridBaseCOP current_grid(grid_map_entry.second);

		//for(core::Size atom_index = 1; atom_index <= residue.nheavyatoms();++atom_index)
		//{
		core::Real current_score(current_grid->score(residue,max_score,qsar_map_));
		component_score += current_score;
		//}
		total_score += component_score;
		std::pair<std::string, core::Real> new_score(current_grid->get_type(),component_score);
		score_map_.insert(new_score);
	}

	if ( norm_function_ ) {
		core::Real normalized_score = (*norm_function_)(total_score,*residue.residue());
		TR.Trace << "Score normalized from " << total_score << " to "<< normalized_score << std::endl;
		return normalized_score;
	} else {
		return total_score;
	}
}

core::Real GridManager::total_score(core::conformation::Residue const & residue)
{
	score_map_.clear();

	core::Real total_score =0.0;
	const core::Real max_score = 9999.0;
	for ( GridMap::value_type const & grid_map_entry: grid_map_ ) {
		core::Real component_score =0;
		GridBaseCOP current_grid(grid_map_entry.second);

		//for(core::Size atom_index = 1; atom_index <= residue.nheavyatoms();++atom_index)
		//{
		core::Real current_score(current_grid->score(residue,max_score,qsar_map_));
		component_score += current_score;
		//}
		total_score += component_score;
		std::pair<std::string, core::Real> new_score(current_grid->get_type(),component_score);
		score_map_.insert(new_score);
	}

	if ( norm_function_ ) {
		core::Real normalized_score = (*norm_function_)(total_score,residue);
		TR.Trace << "Score normalized from " << total_score << " to "<< normalized_score << std::endl;
		return normalized_score;
	} else {
		return total_score;
	}
}

core::Real GridManager::total_score(core::pose::Pose const & pose, core::Size const chain_id)
{
	utility::vector1< core::Size > residues_in_chain( core::pose::get_resnums_for_chain_id( pose, chain_id ) );
	return total_score( pose, residues_in_chain );
}

core::Real GridManager::total_score(core::pose::Pose const & pose, utility::vector1< core::Size > const & residues)
{
	score_map_.clear();
	core::Real total_score = 0.0;
	const core::Real max_score = 9999.0;

	for ( GridMap::value_type const & grid_map_entry: grid_map_ ) {
		core::Real component_score = 0;
		GridBaseCOP current_grid(grid_map_entry.second);
		core::Real weight(grid_weights_[grid_map_entry.first]);
		for ( core::Size resi : residues ) {
			core::conformation::Residue const & residue( pose.residue( resi ) );
			core::Real current_score(current_grid->score(residue,max_score,qsar_map_));
			component_score += current_score;
		}
		total_score += component_score*weight;
		std::pair<std::string, core::Real> new_score(current_grid->get_type(),component_score);
		score_map_.insert(new_score);
	}

	if ( norm_function_ ) {
		core::conformation::ResidueCOPs residue_cops;
		for ( core::Size resi : residues ) {
			residue_cops.push_back( pose.residue( resi ).get_self_ptr() );
		}
		core::Real normalized_score = (*norm_function_)(total_score,residue_cops);
		TR.Trace << "Score normalized from " << total_score << " to "<< normalized_score << std::endl;
		return normalized_score;
	} else {
		return total_score;
	}
}

std::map<std::string,core::Real> GridManager::atom_score(core::pose::Pose const & /*pose*/, core::conformation::Residue const & residue, core::Size atomindex )
{
	std::map<std::string,core::Real> score_map;
	for ( GridMap::value_type const & grid_map_entry: grid_map_ ) {
		GridBaseCOP current_grid(grid_map_entry.second);
		core::Real weight(grid_weights_[grid_map_entry.first]);
		core::Real atom_score = current_grid->atom_score(residue,atomindex,qsar_map_);
		std::string grid_type = current_grid->get_type();
		score_map.insert(std::make_pair(grid_type,atom_score*weight));
	}
	return score_map;
}

void GridManager::update_grids(core::pose::Pose const & pose, core::Vector const & center,utility::vector1<core::Size> ligand_chain_ids_to_exclude)
{
	for ( GridMap::value_type const & grid_map_entry: grid_map_ ) {
		GridBaseOP current_grid(grid_map_entry.second);
		current_grid->initialize(center,width_,resolution_);
		current_grid->set_chain(chain_);
		current_grid->refresh(pose,center,ligand_chain_ids_to_exclude);
	}
}

//void GridManager::reset_grids()
//{
// std::map<std::string,GridBaseOP>::iterator map_iterator(grid_map_.begin());
// for(;map_iterator != grid_map_.end(); ++map_iterator)
// {
//  GridBaseOP current_grid(*map_iterator->second);
//  current_grid->reset();
// }
// grid_map_cache_.clear();
//}

void GridManager::update_grids(core::pose::Pose const & pose, core::Vector const & center, core::Size const & ligand_chain_id_to_exclude)
{
	for ( GridMap::value_type const & grid_map_entry: grid_map_ ) {
		GridBaseOP current_grid(grid_map_entry.second);
		current_grid->initialize(center,width_,resolution_);
		current_grid->set_chain(chain_);
		current_grid->refresh(pose,center,ligand_chain_id_to_exclude);
	}
}

void GridManager::update_grids(core::pose::Pose const & pose,  core::Vector const & center, bool multi) //def: multi=false
{

	std::string chain_hash;

	// We don't need to be all that accurate with the position of the center. If we're off by 0.1 Ang, it shouldn't make too much difference.
	if ( multi ) {
		chain_hash = core::pose::get_sha1_hash_from_chain(chain_, pose, numeric::truncate_and_serialize_xyz_vector(center,1) );
	} else {
		chain_hash = core::pose::get_sha1_hash_excluding_chain(chain_, pose, numeric::truncate_and_serialize_xyz_vector(center,1) );
	}

	std::map<std::string,ConstGridMap>::const_iterator grid_cache_entry(grid_map_cache_.find(chain_hash));

	if ( grid_cache_entry != grid_map_cache_.end() ) { //we've already seen this conformation, load the associated grid out of the map
		TR << "Found a conformation matching hash: " << chain_hash << " Loading from grid cache" <<std::endl;
		set_grid_map_from_cache( chain_hash );
		return;
	}

	bool grid_directory_active = basic::options::option[basic::options::OptionKeys::qsar::grid_dir].user();

	if ( !grid_directory_active ) {
		TR.Warning << "option -qsar:grid_dir is not set.  Use this flag to specify a directory to store scoring grids.  This will save you a huge amount of time" <<std::endl;
	} else if ( !utility::file::file_exists(basic::options::option[basic::options::OptionKeys::qsar::grid_dir]()) ) {
		utility_exit_with_message(basic::options::option[basic::options::OptionKeys::qsar::grid_dir]()+" does not exist. specify a valid path with -qsar:grid_dir");
	}

	if ( basic::options::option[basic::options::OptionKeys::qsar::max_grid_cache_size].user() &&
			grid_map_cache_.size()  >= core::Size(basic::options::option[basic::options::OptionKeys::qsar::max_grid_cache_size]() ) ) {
		TR << "Grid cache exceeds max_cache_size, clearing old scoring grids to save memory." <<std::endl;
		grid_map_cache_.clear();
	}

	//Try to read it off the disk
	if ( grid_directory_active ) {
		//files are in the format grid_directory/hash.json.gz
		std::string directory_path(basic::options::option[basic::options::OptionKeys::qsar::grid_dir]());
		utility::io::izstream grid_file(directory_path+"/"+chain_hash+".json.gz");
		if ( grid_file ) {
			utility::json_spirit::mValue gridmap_data;
			utility::json_spirit::read(grid_file,gridmap_data);
			deserialize(gridmap_data.get_array());
			//Now grid_map_ is whatever was in that file.  We never want to do this again, put it in the cache
			TR << "successfully read grids from the disk for conformation matching hash" << chain_hash <<std::endl;
			insert_into_cache( chain_hash, grid_map_ );
			return;
		}
	}

	// This is a new conformation

	TR << "No conformation matching hash: " << chain_hash << " Updating grid and adding it to the cache" <<std::endl;

	for ( GridMap::value_type const & grid_map_entry: grid_map_ ) {
		GridBaseOP current_grid(grid_map_entry.second);
		TR.Debug <<"updating grid " << grid_map_entry.first << std::endl;
		current_grid->initialize(center,width_,resolution_);
		current_grid->set_chain(chain_);
		current_grid->refresh(pose,center);
		TR.Debug <<"done updating grid" <<std::endl;
	}

	insert_into_cache( chain_hash, grid_map_ );

	if ( grid_directory_active ) {
		//if we just made a grid, we should write it to the disk for safekeeping.
		std::string directory_path(basic::options::option[basic::options::OptionKeys::qsar::grid_dir]());
		std::string temp_path(directory_path+"/"+chain_hash+".inprogress");
		if ( !utility::file::file_exists(temp_path) ) {  //If the inprogress file is there something else is busy writing
			utility::io::ozstream progress_file(temp_path);
			progress_file << "temp" <<std::endl;
			progress_file.close();

			utility::io::ozstream grid_file(directory_path+"/"+chain_hash+".json.gz");

			grid_file << utility::json_spirit::write(serialize()) << std::endl;
			grid_file.close();
			utility::file::file_delete(temp_path);

			TR << "wrote grid matching hash: " << chain_hash << " to disk" <<std::endl;
		}
	}

}

void GridManager::initialize_all_grids(core::Vector const & center)
{
	if ( grid_map_.size() == 0 ) {
		utility_exit_with_message("hmm, no grids in the grid manager. Are they defined in the XML script?");
	}
	if ( !initialized_ ) {
		for ( GridMap::value_type const & grid_map_entry: grid_map_ ) {
			GridBaseOP current_grid(grid_map_entry.second);
			current_grid->initialize(center,width_,resolution_);
		}
		initialized_ = true;
	}

}

core::Size GridManager::size()
{
	return grid_map_.size();
}

ScoreMap GridManager::get_cached_scores()
{
	return score_map_;
}

void GridManager::append_cached_scores(jd2::JobOP job)
{
	std::map<std::string,core::Real>::iterator map_iterator(score_map_.begin());
	int total_score =0;
	for ( ; map_iterator !=score_map_.end(); ++map_iterator ) {
		std::string type_name(map_iterator->first);
		int score(static_cast<int>(map_iterator->second));
		total_score +=score;

		//jd2::PDBJobOutputter::

		job->add_string_real_pair("grid_"+type_name,score);
	}
	job->add_string_real_pair("grid_total",total_score);
}

void GridManager::write_grids(std::string prefix)
{
	for ( GridMap::value_type const & grid_map_entry: grid_map_ ) {
		GridBaseCOP current_grid(grid_map_entry.second);
		current_grid->dump_BRIX(prefix);
	}
}

utility::json_spirit::Value GridManager::serialize()
{
	using utility::json_spirit::Value;
	std::vector<Value> gridmap_data;
	for ( GridMap::value_type const & grid_map_entry: grid_map_ ) {
		Value grid_name(grid_map_entry.first);
		Value grid(grid_map_entry.second->serialize());
		std::vector<Value> grid_pair_values;
		grid_pair_values.push_back(grid_name);
		grid_pair_values.push_back(grid);
		Value grid_pair(grid_pair_values);
		gridmap_data.push_back(grid_pair);
	}
	return Value(gridmap_data);
}

void GridManager::deserialize(utility::json_spirit::mArray data)
{
	grid_map_.clear();
	for ( utility::json_spirit::mArray::iterator it = data.begin(); it != data.end(); ++it ) {
		utility::json_spirit::mArray grid_data(it->get_array());
		std::string grid_name = grid_data[0].get_str();
		GridBaseOP grid(GridFactory::get_instance()->new_grid(grid_data[1].get_obj()));
		grid_map_[grid_name] = grid;
	}
}

bool GridManager::is_in_grid(utility::vector1<core::conformation::UltraLightResidue> const & residues)
{
	for ( core::conformation::UltraLightResidue const & residue: residues ) {
		if ( is_in_grid(residue) == false ) {
			return false;
		}
	}

	return true;

}


bool GridManager::is_in_grid(core::conformation::UltraLightResidue const & residue)
{
	for ( GridMap::value_type const & grid_map_entry: grid_map_ ) {
		GridBaseCOP current_grid(grid_map_entry.second);
		if ( !current_grid->is_in_grid(residue) ) {
			return false;
		}
	}
	return true;
}

bool GridManager::is_in_grid(core::conformation::Residue const & residue)
{
	for ( GridMap::value_type const & grid_map_entry: grid_map_ ) {
		GridBaseCOP current_grid(grid_map_entry.second);
		if ( !current_grid->is_in_grid(residue) ) {
			return false;
		}
	}
	return true;
}

/// @brief Insert the given GridMap into the grid_map_cache under the given index value.
void
GridManager::insert_into_cache( std::string const & hash_val, GridMap const & grid_map ) {
	if ( grid_map_cache_.count( hash_val ) != 0 ) {
		TR.Warning << "[ WARNING ] The Grid cache already has a value for setting! Reseting it anyway." << std::endl;
		// Probably should be a utility_exit, as this indicates a serious logic error.
	}
	ConstGridMap & cached_gridmap = grid_map_cache_[ hash_val ]; // Will default initialize the cache map.
	cached_gridmap.clear(); // Only needed if we're not utility_exit-ing above

	// We do a deep copy into the cache, so further modifications to the passed value don't change the cached value.
	for ( GridMap::value_type const & entry: grid_map ) {
		cached_gridmap[ entry.first ] = entry.second->clone();
	}
}

/// @brief Reset the grid_map_ variable to the cached state
/// @details As grid_map_ is modifiable, make a copy of cached value
void
GridManager::set_grid_map_from_cache( std::string const & hash_val ) {
	if ( grid_map_cache_.count( hash_val ) == 0 ) {
		utility_exit_with_message("Cannot get cached GridMap - it doesn't exist!");
	}
	grid_map_.clear(); // Dump the current values.

	ConstGridMap const & cached_gridmap = grid_map_cache_[ hash_val ];
	// We do a deep copy into grid_map_, so further modifications to it don't change the cached value.
	for ( ConstGridMap::value_type const & entry: cached_gridmap ) {
		grid_map_[ entry.first ] = entry.second->clone();
	}

}

}
}
}
