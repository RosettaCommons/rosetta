// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/ligand_docking/qsar/qsarMover.cc
/// @author Sam DeLuca

#include <protocols/qsar/qsarMover.hh>
#include <protocols/qsar/qsarMap.hh>
#include <protocols/qsar/scoring_grid/GridManager.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/util.hh>
#include <protocols/geometry/RB_geometry.hh>
#include <basic/Tracer.hh>
#include <protocols/jd2/JobDistributor.hh>

namespace protocols {
namespace qsar {
static basic::Tracer TR("protocols.ligand_docking.qsar.qsarMover");

qsarMover::qsarMover(core::Real width, core::Real resolution):
		grid_manager_(scoring_grid::GridManager::get_instance()), qsar_map_(0), initialize_(false)
{
	grid_manager_->set_dimensions(width,resolution);
}
void qsarMover::set_chain(std::string chain_id)
{
	chain_id_ = chain_id;
}

//qsarMover::qsarMover(core::Real width, core::Real resolution) : grid_manager_(new scoring_grid::GridManager(width,resolution)), initialize_(false)
//{

//}

void qsarMover::add_grid(std::string grid_name)
{
	grids_to_make_.push_back(grid_name);
}

void qsarMover::write_all_grids(std::string prefix)
{
	grid_manager_->write_grids(prefix);
}

std::map<std::string,core::Real> qsarMover::get_cached_scores()
{
	return grid_manager_->get_cached_scores();
}

void qsarMover::apply(core::pose::Pose & pose)
{

	if(grids_to_make_.size()==0)
	{
		TR << "WARNING: no grids specified, QSAR scoring function will be empty!!" <<std::endl;
		return;
	}else if(!initialize_)
	{
		utility::vector1<std::string>::iterator grid_iterator(grids_to_make_.begin());
		for(; grid_iterator != grids_to_make_.end();++grid_iterator)
		{
			grid_manager_->make_new_grid(*grid_iterator);
			TR.Debug << "made new " << *grid_iterator << " grid" <<std::endl;
		}
		if(qsar_map_ == 0)
		{
			core::Size const chain_num(core::pose::get_chain_id_from_chain(chain_id_,pose));
			core::Size const begin(pose.conformation().chain_begin(chain_num));
			core::conformation::ResidueOP residue = new core::conformation::Residue(pose.residue(begin));
			qsar_map_ = new qsarMap("default",residue);

			qsar_map_->fill_with_value(1,grids_to_make_);
		}

		grid_manager_->set_qsar_map(qsar_map_);

	}

	if(!initialize_ || pose.conformation().structure_moved())
	{
		core::Size const chain_num(core::pose::get_chain_id_from_chain(chain_id_,pose));
		core::Size const jump_id(core::pose::get_jump_id_from_chain_id(chain_num,pose));
		core::Vector const center(geometry::downstream_centroid_by_jump(pose,jump_id));

		core::Size const begin(pose.conformation().chain_begin(chain_num));
		//core::Size const end(pose.conformation().chain_end(chain_num));

		core::conformation::Residue const & residue(pose.residue(begin));

		if(!initialize_)
		{
			TR.Debug <<"initializing grids" << std::endl;
			grid_manager_->initialize_all_grids(center);
			TR.Debug <<"grids initialized" <<std::endl;
		}
		grid_manager_->update_grids(pose,center);
		TR.Debug <<"grids updated, scoring.."<<std::endl;
		core::Real score(grid_manager_->total_score(residue));
		TR.Debug << "total score is " << score <<std::endl;
		std::map<std::string,core::Real> scores(grid_manager_->get_cached_scores());


		jd2::JobOP job(jd2::JobDistributor::get_instance()->current_job());
		grid_manager_->append_cached_scores(job);
		//grid_manager_->write_grids("test_");
	}
	if(!initialize_)
		initialize_=true;
}

std::string qsarMover::get_name() const
{
	return "qsarMover";
}

}
}
