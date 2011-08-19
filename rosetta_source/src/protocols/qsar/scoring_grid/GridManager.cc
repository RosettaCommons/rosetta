// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   /git/src/protocols/ligand_docking/scoring_grid/GridManager.cc
/// @author Sam DeLuca

#include <protocols/qsar/scoring_grid/GridManager.hh>
#include <protocols/qsar/scoring_grid/GridBase.hh>
//#include <protocols/qsar/qsarTypeManager.hh>

#include <protocols/qsar/scoring_grid/GridFactory.hh>
#include <protocols/jd2/Job.hh>

#include <core/pose/util.hh>

//#include <protocols/jd2/PDBJobOutputter.hh>
//#include <utility/exit.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/qsar.OptionKeys.gen.hh>
#include <basic/database/open.hh>
#include <core/conformation/Residue.hh>
#include <basic/Tracer.hh>

#include <utility/io/izstream.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>

#include <iostream>
#include <fstream>
#include <map>

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

GridManager::GridManager() : last_tag_(""),width_(0), resolution_(0), qsar_map_(0),initialized_(false)
{
	std::string base_name = basic::options::option[basic::options::OptionKeys::qsar::weights]();

	utility::io::izstream infile;
	bool file_status = basic::database::open(infile,"qsar/weights/"+base_name+".wts");
	if(!file_status)
	{
		utility_exit_with_message("could not open "+ base_name+" weight file, check your paths");
	}
	std::string line;
	do
	{
		getline(infile,line);
		if(line.size()> 0 && line[0] != '#')
		{
			utility::vector1<std::string> split_line = static_cast<utility::vector1<std::string> >(utility::string_split(line,'\t'));
			if(split_line.size() == 2)
			{
				//std::cout << line <<std::endl;
				std::string qsar_name(split_line[1]);
				//qsarType type(qsarTypeManager::qsar_type_from_name(qsar_name));
				core::Real qsar_weight=utility::from_string(split_line[2],core::Real(0));
				weight_map_.insert(std::pair<std::string,core::Real>(qsar_name,qsar_weight));
			}
		}
	}while(!infile.eof());
	infile.close();
}


void GridManager::set_dimensions(core::Real width, core::Real resolution)
{
	width_ = width;
	resolution_=resolution;
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

void GridManager::make_new_grid(std::string grid_name)
{

	//qsar::qsarType grid_type(qsar::qsarTypeManager::qsar_type_from_name(grid_name));
	core::Real weight(weight_map_[grid_name]);
	//GridBaseOP new_grid;
	GridManagerTracer.Debug << grid_name <<std::endl;
	//GridFactoryOP grid_factory(GridFactory::get_instance());
	GridBaseOP new_grid(GridFactory::get_instance()->new_grid(grid_name));
	new_grid->set_weight(weight);
	insert_grid(new_grid);
}

void GridManager::insert_grid(GridBaseOP grid)
{
	std::string type(grid->get_type());
	std::pair<std::string, GridBaseOP> new_grid(type,grid);
	grid_map_.insert(new_grid);
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
	return total_score;

}

void GridManager::update_grids(core::pose::Pose const & pose, core::Vector const & center,utility::vector1<core::Size> ligand_chain_ids_to_exclude)
{
	std::map<std::string,GridBaseOP>::iterator map_iterator(grid_map_.begin());
	for(;map_iterator != grid_map_.end(); ++map_iterator)
	{
		GridBaseOP current_grid(*map_iterator->second);
		current_grid->refresh(pose,center,ligand_chain_ids_to_exclude);
	}
}


void GridManager::update_grids(core::pose::Pose const & pose, core::Vector const & center, core::Size const & ligand_chain_id_to_exclude)
{
	std::map<std::string,GridBaseOP>::iterator map_iterator(grid_map_.begin());
	for(;map_iterator != grid_map_.end();++map_iterator)
	{
		GridBaseOP current_grid(*map_iterator->second);
		current_grid->refresh(pose,center,ligand_chain_id_to_exclude);
	}
}


void GridManager::update_grids(core::pose::Pose const & pose,  core::Vector const & center)
{
	std::map<std::string,GridBaseOP>::iterator map_iterator(grid_map_.begin());

	for(;map_iterator != grid_map_.end();++map_iterator)
	{


		GridBaseOP current_grid(*map_iterator->second);
		GridManagerTracer.Debug <<"updating grid " << map_iterator->first << std::endl;
		current_grid->refresh(pose,center);
		GridManagerTracer.Debug <<"done updating grid" <<std::endl;
	}
}

void GridManager::update_grids(core::pose::Pose const & pose, core::Vector const & center,std::string const & tag)
{
	//std::cout << tag << " " <<last_tag_ <<std::endl;
	if(tag != last_tag_)
	{
		update_grids(pose,center);
	}
	last_tag_ = tag;
}

void GridManager::initialize_all_grids(core::Vector const & center)
{
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

std::map<std::string, core::Real> GridManager::get_cached_scores()
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

}
}
}
