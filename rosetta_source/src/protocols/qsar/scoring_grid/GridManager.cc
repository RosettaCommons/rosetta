// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/ligand_docking/scoring_grid/GridManager.cc
/// @author Sam DeLuca

#include <protocols/qsar/scoring_grid/GridManager.hh>
#include <protocols/qsar/scoring_grid/SingleGrid.hh>
#include <protocols/qsar/scoring_grid/GridFactory.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/qsar/qsarMap.hh>

#include <core/pose/util.hh>
#include <core/conformation/Residue.hh>
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

//STL headers
#include <iostream>
#include <fstream>
#include <map>

//Auto Headers
#include <boost/bind.hpp>

namespace protocols {
namespace qsar {
namespace scoring_grid {

static basic::Tracer GridManagerTracer("protocols.qsar.scoring_grid.GridManager");

GridManager* GridManager::instance_(0);

GridManager* GridManager::get_instance()
{
	if(instance_ == 0)
	{
		instance_ = new GridManager();
	}
	return instance_;
}

void GridManager::reset()
{
	grid_map_.clear();
	score_map_.clear();
	last_tag_ = "";
	width_ = 40;
	resolution_ = 0.25;
	qsar_map_ = 0;
	initialized_ = false;
	chain_ = 'X';
	normalized_ = false;
}

GridManager::GridManager() :
	last_tag_(""),
	width_(40),
	resolution_(0.25),
	qsar_map_(0),
	initialized_(false),
	chain_('X'),
	normalized_(false)
{
	grid_map_.clear();
	score_map_.clear();
}

void GridManager::set_normalized(bool normalized)
{
	normalized_ = normalized;
}

void GridManager::set_width(core::Real width)
{
	width_ = width;
}

void GridManager::set_resolution(core::Real resolution)
{
	resolution_=resolution;
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
	if(qsar_map_ != 0)
	{
		return true;
	}else
	{
		return false;
	}
}

void GridManager::make_new_grid(utility::tag::TagPtr const tag)
{

	std::string name= tag->getName();
	GridManagerTracer.Debug << name <<std::endl;
	GridBaseOP new_grid(GridFactory::get_instance()->new_grid(tag));
	insert_grid(name, new_grid);

}

void GridManager::insert_grid(std::string const name, GridBaseOP const grid)
{
	//if(grid_map_.find(name) != grid_map_.end()){
	//	utility_exit_with_message("2 grids with the same name!");
	//}
	grid_map_[name] = grid;
}

GridBaseOP GridManager::get_grid(std::string const & grid_type)
{
	return grid_map_.find(grid_type)->second;
}

utility::vector1<std::string> GridManager::get_grid_names()
{
	utility::vector1<std::string> grid_names;

	std::map<std::string,GridBaseOP>::const_iterator it;
	for(it = grid_map_.begin();it != grid_map_.end();++it)
	{
		grid_names.push_back(it->first);
	}
	return grid_names;
}

core::Real GridManager::total_score(core::conformation::Residue const & residue)
{
	score_map_.clear();

	core::Real total_score =0.0;
	const core::Real max_score = 9999.0;
	std::map<std::string,GridBaseOP>::iterator map_iterator(grid_map_.begin());
	for(;map_iterator != grid_map_.end();++map_iterator)
	{
		core::Real component_score =0;
		GridBaseOP current_grid(*map_iterator->second);

		//for(core::Size atom_index = 1; atom_index <= residue.nheavyatoms();++atom_index)
		//{
		core::Real current_score(current_grid->score(residue,max_score,qsar_map_));
		component_score += current_score;
		//}
		total_score += component_score;
		std::pair<std::string, core::Real> new_score(current_grid->get_type(),component_score);
		score_map_.insert(new_score);
	}

	if(normalized_)
	{
		return total_score/static_cast<core::Real>(residue.natoms());
	}
	return total_score;
}

core::Real GridManager::total_score(core::pose::Pose const & pose, core::Size const chain_id)
{
	score_map_.clear();
	core::Real total_score = 0.0;
	const core::Real max_score = 9999.0;

	core::conformation::ResidueCOPs residue_vector = core::pose::get_chain_residues(pose,chain_id);

	std::map<std::string,GridBaseOP>::iterator map_iterator(grid_map_.begin());
	for(;map_iterator != grid_map_.end();++map_iterator)
	{
		core::Real component_score = 0;
		GridBaseOP current_grid(*map_iterator->second);
		for(core::Size residue_count = 1; residue_count <= residue_vector.size(); ++residue_count )
		{
			core::conformation::ResidueCOP residue(residue_vector[residue_count]);
			//for(core::Size atom_index = 1; atom_index <= residue->nheavyatoms(); ++atom_index)
			//{
			core::Real current_score(current_grid->score(*residue,max_score,qsar_map_));
			component_score += current_score;
			//}
		}
		total_score += component_score;
		std::pair<std::string, core::Real> new_score(current_grid->get_type(),component_score);
		score_map_.insert(new_score);
	}

	if(normalized_)
	{
		core::Size n_atoms = 0;
		for(core::Size i = 1; i < residue_vector.size();++i)
		{
			n_atoms += residue_vector[i]->natoms();
		}
		return total_score/static_cast<core::Real>(n_atoms);
	}

	return total_score;

}

void GridManager::update_grids(core::pose::Pose const & pose, core::Vector const & center,utility::vector1<core::Size> ligand_chain_ids_to_exclude)
{
	std::map<std::string,GridBaseOP>::iterator map_iterator(grid_map_.begin());
	for(;map_iterator != grid_map_.end(); ++map_iterator)
	{
		GridBaseOP current_grid(*map_iterator->second);
		current_grid->initialize(center,width_,resolution_);
		current_grid->set_chain(chain_);
		current_grid->refresh(pose,center,ligand_chain_ids_to_exclude);
	}
}


void GridManager::update_grids(core::pose::Pose const & pose, core::Vector const & center, core::Size const & ligand_chain_id_to_exclude)
{
	std::map<std::string,GridBaseOP>::iterator map_iterator(grid_map_.begin());
	for(;map_iterator != grid_map_.end();++map_iterator)
	{
		GridBaseOP current_grid(*map_iterator->second);
		current_grid->initialize(center,width_,resolution_);
		current_grid->set_chain(chain_);
		current_grid->refresh(pose,center,ligand_chain_id_to_exclude);
	}
}


void GridManager::update_grids(core::pose::Pose const & pose,  core::Vector const & center)
{

	core::Size chain_hash = core::pose::get_hash_excluding_chain(chain_,pose);
	std::string hash_string(utility::to_string(chain_hash));
	std::map<core::Size,GridMap>::const_iterator grid_cache_entry(grid_map_cache_.find(chain_hash));

	bool grid_directory_active = basic::options::option[basic::options::OptionKeys::qsar::grid_dir].user();

	if(!grid_directory_active)
	{
		GridManagerTracer << "WARNING: option -qsar:grid_dir is not set.  Use this flag to specify a directory to store scoring grids.  This will save you a huge amount of time" <<std::endl;
	}

	if(grid_cache_entry != grid_map_cache_.end()) //we've already seen this conformation, load the associated grid out of the map
	{
		GridManagerTracer << "Found a conformation matching hash: " << chain_hash << " Loading from grid cache" <<std::endl;
		grid_map_ = grid_cache_entry->second;
	}else // This is a new conformation
	{

		//Try to read it off the disk
		if(grid_directory_active)
		{
			//files are in the format grid_directory/hash.json.gz
			std::string directory_path(basic::options::option[basic::options::OptionKeys::qsar::grid_dir]());
			utility::io::izstream grid_file(directory_path+"/"+hash_string+".json.gz");
			if(grid_file)
			{
				utility::json_spirit::mValue gridmap_data;
				utility::json_spirit::read(grid_file,gridmap_data);
				deserialize(gridmap_data.get_array());
				//Now grid_map_ is whatever was in that file.  We never want to do this again, put it in the cache
				GridManagerTracer << "successfully read grids from the disk for conformation matching hash" << chain_hash <<std::endl;
				grid_map_cache_.insert(std::make_pair(chain_hash,grid_map_));
				return;

			}
		}

		GridManagerTracer << "No conformation matching hash: " << chain_hash << " Updating grid and adding it to the cache" <<std::endl;

		std::map<std::string,GridBaseOP>::iterator map_iterator(grid_map_.begin());

		for(;map_iterator != grid_map_.end();++map_iterator)
		{

			GridBaseOP current_grid(*map_iterator->second);
			GridManagerTracer.Debug <<"updating grid " << map_iterator->first << std::endl;
			current_grid->initialize(center,width_,resolution_);
			current_grid->set_chain(chain_);
			current_grid->refresh(pose,center);
			GridManagerTracer.Debug <<"done updating grid" <<std::endl;
		}
		grid_map_cache_.insert(std::make_pair(chain_hash,grid_map_));

		if(grid_directory_active)
		{
			//if we just made a grid, we should write it to the disk for safekeeping.
			std::string directory_path(basic::options::option[basic::options::OptionKeys::qsar::grid_dir]());
			std::string temp_path(directory_path+"/"+hash_string+".inprogress");
			if(!utility::file::file_exists(temp_path))  //If the inprogress file is there something else is busy writing
			{
				utility::io::ozstream progress_file(temp_path);
				progress_file << "temp" <<std::endl;
				progress_file.close();

				utility::io::ozstream grid_file(directory_path+"/"+hash_string+".json.gz");

				grid_file << utility::json_spirit::write(serialize()) << std::endl;
				grid_file.close();
				utility::file::file_delete(temp_path);

				GridManagerTracer << "wrote grid matching hash: " << chain_hash << " to disk" <<std::endl;

			}

		}

	}
}

void GridManager::initialize_all_grids(core::Vector const & center)
{
	if(grid_map_.size() == 0)
		utility_exit_with_message("hmm, no grids in the grid manager. Are they defined in the XML script?");
	if(!initialized_)
	{
		std::map<std::string,GridBaseOP>::iterator map_iterator(grid_map_.begin());
		for(;map_iterator != grid_map_.end();++map_iterator)
		{
			GridBaseOP current_grid(*map_iterator->second);
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
	for(;map_iterator !=score_map_.end(); ++map_iterator)
	{
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
	std::map<std::string, GridBaseOP>::iterator map_iterator(grid_map_.begin());
	for(;map_iterator != grid_map_.end(); ++map_iterator)
	{
		GridBaseOP current_grid(map_iterator->second);
		current_grid->dump_BRIX(prefix);
	}
}

utility::json_spirit::Value GridManager::serialize()
{
	using utility::json_spirit::Value;
	std::vector<Value> gridmap_data;
	for(GridMap::iterator it = grid_map_.begin(); it != grid_map_.end();++it)
	{
		Value grid_name(it->first);
		Value grid(it->second->serialize());
		gridmap_data.push_back(Value(utility::tools::make_vector(grid_name,grid)));
	}
	return Value(gridmap_data);
}

void GridManager::deserialize(utility::json_spirit::mArray data)
{
	grid_map_.clear();
	for(utility::json_spirit::mArray::iterator it = data.begin(); it != data.end();++it)
	{
		utility::json_spirit::mArray grid_data(it->get_array());
		std::string grid_name = grid_data[0].get_str();
		GridBaseOP grid(GridFactory::get_instance()->new_grid(grid_data[1].get_obj()));
		grid_map_[grid_name] = grid;
	}
}

}
}
}
