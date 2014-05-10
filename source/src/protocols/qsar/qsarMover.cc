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
#include <protocols/qsar/qsarMoverCreator.hh>
#include <protocols/qsar/qsarMap.hh>
#include <protocols/qsar/scoring_grid/SingleGrid.hh>
#include <protocols/qsar/scoring_grid/GridManager.hh>
#include <protocols/rigid/RB_geometry.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/util.hh>
#include <core/pose/Pose.hh>

#include <basic/Tracer.hh>

#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <utility/vector0.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace qsar {
static basic::Tracer TR("protocols.ligand_docking.qsar.qsarMover");


std::string
qsarCreator::keyname() const
{
	return qsarCreator::mover_name();
}

protocols::moves::MoverOP
qsarCreator::create_mover() const {
	return new qsarMover;
}

std::string
qsarCreator::mover_name()
{
	return "qsar";
}

qsarMover::qsarMover():
		qsar_map_(0),
		chain_(),
		initialize_(false)
{}

//@brief parse XML (specifically in the context of the parser/scripting scheme)
void
qsarMover::parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & /*datamap*/,
		protocols::filters::Filters_map const & /*filters*/,
		protocols::moves::Movers_map const & /*movers*/,
		core::pose::Pose const & /*pose*/
)
{
	if ( tag->getName() != "qsarMover" ) throw utility::excn::EXCN_RosettaScriptsOption("This should be impossible");

	if ( ! tag->hasOption("chain") ) throw utility::excn::EXCN_RosettaScriptsOption("'qsar' mover requires chain tag");
	chain_= tag->getOption<std::string>("chain");

	if ( ! tag->hasOption("grids") ) throw utility::excn::EXCN_RosettaScriptsOption("'qsarMover' requires grids tag");

	std::string grids_string= tag->getOption<std::string>("grids");
	utility::vector1<std::string> grid_strings(utility::string_split(grids_string, ','));
	grids_to_use_ = grid_strings;
}

void qsarMover::apply(core::pose::Pose & pose)
{
	scoring_grid::GridManager* grid_manager(scoring_grid::GridManager::get_instance());

	core::Size const chain_id(core::pose::get_chain_id_from_chain(chain_,pose));
	core::Size const begin(pose.conformation().chain_begin(chain_id));
	/// TODO The next line assumes the chain is one residue. fix this.
	core::conformation::ResidueOP residue = new core::conformation::Residue(pose.residue(begin));

	if(grids_to_use_.size()==0)
	{
		TR << "WARNING: no grids specified, QSAR scoring function will be empty!!" <<std::endl;
		return;
	}else if(!initialize_)
	{
		utility::vector1<std::string>::iterator grid_iterator(grids_to_use_.begin());
		for(; grid_iterator != grids_to_use_.end();++grid_iterator)
		{
			//grid_manager->make_new_grid(*grid_iterator);
			grid_manager->get_grid(*grid_iterator);
			TR.Debug << "getting grid: " << *grid_iterator << std::endl;
		}
		if(qsar_map_ == 0)
		{

			qsar_map_ = new qsarMap("default",residue);

			qsar_map_->fill_with_value(1,grids_to_use_);
		}

		grid_manager->set_qsar_map(qsar_map_);

	}

	if(!initialize_ || pose.conformation().structure_moved())
	{
		core::Size const jump_id(core::pose::get_jump_id_from_chain_id(chain_id,pose));
		core::Vector const center(geometry::downstream_centroid_by_jump(pose,jump_id));

		if(!initialize_)
		{
			TR.Debug <<"initializing grids" << std::endl;
			grid_manager->initialize_all_grids(center);
			TR.Debug <<"grids initialized" <<std::endl;
		}
		grid_manager->update_grids(pose,center);
		TR.Debug <<"grids updated, scoring.."<<std::endl;
		core::Real score(grid_manager->total_score(*residue));
		TR.Debug << "total score is " << score <<std::endl;
		std::map<std::string,core::Real> scores(grid_manager->get_cached_scores());


		jd2::JobOP job(jd2::JobDistributor::get_instance()->current_job());
		grid_manager->append_cached_scores(job);
		//grid_manager->write_grids("test_");
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
