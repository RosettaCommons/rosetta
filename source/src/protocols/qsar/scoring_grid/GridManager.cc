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
#include <utility/thread/threadsafe_creation.hh>

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

static basic::Tracer GridManagerTracer("protocols.qsar.scoring_grid.GridManager");

GridManager* GridManager::instance_(0);

#ifdef MULTI_THREADED
#ifdef CXX11

std::mutex GridManager::singleton_mutex_;

std::mutex & GridManager::singleton_mutex() { return singleton_mutex_; }

#endif
#endif

/// @brief static function to get the instance of ( pointer to) this singleton class
GridManager * GridManager::get_instance()
{
	boost::function< GridManager * () > creator = boost::bind( &GridManager::create_singleton_instance );
	utility::thread::safely_create_singleton( creator, instance_ );
	return instance_;
}

GridManager *
GridManager::create_singleton_instance()
{
	return new GridManager;
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
	norm_function_ = 0;
}

GridManager::GridManager() :
	last_tag_(""),
	width_(40),
	resolution_(0.25),
	qsar_map_(0),
	initialized_(false),
	chain_('X'),
	norm_function_(0)
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
	if(qsar_map_ != 0)
	{
		return true;
	}else
	{
		return false;
	}
}

void GridManager::make_new_grid(utility::tag::TagCOP const tag)
{

	std::string name= tag->getName();
	
	//Sometimes creating grids can be time consuming
	//so check before we do something that calls a grid constructor.
	//Access the various options anyways so we don't die for unaccessed options. 
	core::Real weight = tag->getOption<core::Real>("weight");
	tag->getOption<std::string>("grid_type");
	if(grid_map_.find(name) == grid_map_.end())
	{
		grid_weights_.insert(std::make_pair(name,weight));
		GridManagerTracer.Debug << name <<std::endl;
		GridBaseOP new_grid(GridFactory::get_instance()->new_grid(tag));
		insert_grid(name, new_grid);
	}



}

void GridManager::insert_grid(std::string const name, GridBaseOP const grid)
{
	grid_map_[name] = grid;
}

bool GridManager::is_normalization_enabled()
{
	if(norm_function_)
	{
		return true;
	}
	else
	{
		return false;
	}
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

core::Real GridManager::total_score(core::conformation::UltraLightResidue const & residue)
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

	if(norm_function_)
	{
		core::Real normalized_score = (*norm_function_)(total_score,*residue.residue());
		GridManagerTracer.Trace << "Score normalized from " << total_score << " to "<< normalized_score << std::endl;
		return normalized_score;
	}else
	{
		return total_score;
	}
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

	if(norm_function_)
	{
		core::Real normalized_score = (*norm_function_)(total_score,residue);
		GridManagerTracer.Trace << "Score normalized from " << total_score << " to "<< normalized_score << std::endl;
		return normalized_score;
	}else
	{
		return total_score;
	}
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
		core::Real weight(grid_weights_[map_iterator->first]);
		for(core::Size residue_count = 1; residue_count <= residue_vector.size(); ++residue_count )
		{
			core::conformation::ResidueCOP residue(residue_vector[residue_count]);
			core::Real current_score(current_grid->score(*residue,max_score,qsar_map_));
			component_score += current_score;
		}
		total_score += component_score*weight;
		std::pair<std::string, core::Real> new_score(current_grid->get_type(),component_score);
		score_map_.insert(new_score);
	}

	if(norm_function_)
	{
		core::Real normalized_score = (*norm_function_)(total_score,residue_vector);
		GridManagerTracer.Trace << "Score normalized from " << total_score << " to "<< normalized_score << std::endl;
		return normalized_score;
	}else
	{
		return total_score;
	}
}


std::map<std::string,core::Real> GridManager::atom_score(core::pose::Pose const & /*pose*/, core::conformation::Residue const & residue, core::Size atomindex )
{
	std::map<std::string,core::Real> score_map;
	std::map<std::string,GridBaseOP>::iterator map_iterator(grid_map_.begin());
	for(;map_iterator != grid_map_.end();++map_iterator)
	{
		GridBaseOP current_grid(*map_iterator->second);
		core::Real weight(grid_weights_[map_iterator->first]);
		core::Real atom_score = current_grid->atom_score(residue,atomindex,qsar_map_);
		std::string grid_type = current_grid->get_type();
		score_map.insert(std::make_pair(grid_type,atom_score*weight));
	}
	return score_map;
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

	std::string chain_hash = core::pose::get_sha1_hash_excluding_chain(chain_,pose);
	std::map<std::string,GridMap>::const_iterator grid_cache_entry(grid_map_cache_.find(chain_hash));

	bool grid_directory_active = basic::options::option[basic::options::OptionKeys::qsar::grid_dir].user();

	if(!grid_directory_active)
	{
		GridManagerTracer << "WARNING: option -qsar:grid_dir is not set.  Use this flag to specify a directory to store scoring grids.  This will save you a huge amount of time" <<std::endl;
	}else if(!utility::file::file_exists(basic::options::option[basic::options::OptionKeys::qsar::grid_dir]()))
    {
        utility_exit_with_message(basic::options::option[basic::options::OptionKeys::qsar::grid_dir]()+" does not exist. specify a valid path with -qsar:grid_dir");
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
			utility::io::izstream grid_file(directory_path+"/"+chain_hash+".json.gz");
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

		if(basic::options::option[basic::options::OptionKeys::qsar::max_grid_cache_size].user() &&
			grid_map_cache_.size() >= core::Size(basic::options::option[basic::options::OptionKeys::qsar::max_grid_cache_size]() ) )
		{
			GridManagerTracer << "Grid cache exceeds max_cache_size, clearing old scoring grids to save memory." <<std::endl;
			grid_map_cache_.clear();
		}
		grid_map_cache_.insert(std::make_pair(chain_hash,grid_map_));

		if(grid_directory_active)
		{
			//if we just made a grid, we should write it to the disk for safekeeping.
			std::string directory_path(basic::options::option[basic::options::OptionKeys::qsar::grid_dir]());
			std::string temp_path(directory_path+"/"+chain_hash+".inprogress");
			if(!utility::file::file_exists(temp_path))  //If the inprogress file is there something else is busy writing
			{
				utility::io::ozstream progress_file(temp_path);
				progress_file << "temp" <<std::endl;
				progress_file.close();

				utility::io::ozstream grid_file(directory_path+"/"+chain_hash+".json.gz");

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
	for(utility::json_spirit::mArray::iterator it = data.begin(); it != data.end();++it)
	{
		utility::json_spirit::mArray grid_data(it->get_array());
		std::string grid_name = grid_data[0].get_str();
		GridBaseOP grid(GridFactory::get_instance()->new_grid(grid_data[1].get_obj()));
		grid_map_[grid_name] = grid;
	}
}

bool GridManager::is_in_grid(core::conformation::UltraLightResidue const & residue)
{
    std::map<std::string,GridBaseOP>::iterator map_iterator(grid_map_.begin());
    for(;map_iterator != grid_map_.end();++map_iterator)
    {
        GridBaseOP current_grid(*map_iterator->second);
        if(!current_grid->is_in_grid(residue))
        {
            return false;
        }
    }
    return true;
}

bool GridManager::is_in_grid(core::conformation::Residue const & residue)
{
    std::map<std::string,GridBaseOP>::iterator map_iterator(grid_map_.begin());
    for(;map_iterator != grid_map_.end();++map_iterator)
    {
        GridBaseOP current_grid(*map_iterator->second);
        if(!current_grid->is_in_grid(residue))
        {
            return false;
        }
    }
    return true;
}

}
}
}
