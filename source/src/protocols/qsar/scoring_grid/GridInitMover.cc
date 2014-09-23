// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/qsar/scoring_grid/GridInitMover.cc
/// @author Sam DeLuca

#include <protocols/qsar/scoring_grid/GridInitMover.hh>
#include <protocols/qsar/scoring_grid/GridInitMoverCreator.hh>
#include <protocols/qsar/scoring_grid/GridManager.hh>
#include <protocols/rigid/RB_geometry.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/kinematics/Jump.hh>

#include <utility/excn/Exceptions.hh>
#include <utility/tag/Tag.hh>

namespace protocols {
namespace qsar {
namespace scoring_grid {

std::string GridInitMoverCreator::keyname() const
{
	return GridInitMoverCreator::mover_name();
}

protocols::moves::MoverOP GridInitMoverCreator::create_mover() const
{
	return protocols::moves::MoverOP( new GridInitMover );
}

std::string GridInitMoverCreator::mover_name()
{
	return "GridInitMover";
}

GridInitMover::GridInitMover()
{

}

GridInitMover::~GridInitMover()
{

}

protocols::moves::MoverOP GridInitMover::clone() const
{
	return protocols::moves::MoverOP( new GridInitMover(*this) );
}

protocols::moves::MoverOP GridInitMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new GridInitMover );
}

std::string GridInitMover::get_name() const
{
	return "GridInitMover";
}

void GridInitMover::parse_my_tag
(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & /*data_map*/,
	protocols::filters::Filters_map const & /*filters*/,
	protocols::moves::Movers_map const & /*movers*/,
	core::pose::Pose const & /*pose*/
)
{
	if ( tag->getName() != "GridInitMover" )
	{
		throw utility::excn::EXCN_RosettaScriptsOption("This should be impossible");
	}

	if ( ! tag->hasOption("chain") ) throw utility::excn::EXCN_RosettaScriptsOption("'GridInitMover' mover requires chain tag");

	chain_ = tag->getOption<std::string>("chain");
}

void GridInitMover::apply(core::pose::Pose & pose)
{
	qsar::scoring_grid::GridManager* grid_manager(qsar::scoring_grid::GridManager::get_instance());

	assert(chain_.size() == 1);


	core::Size const chain_id = core::pose::get_chain_id_from_chain(chain_, pose);
	core::Size const jump_id = core::pose::get_jump_id_from_chain_id(chain_id,pose);
	core::Vector const center(protocols::geometry::downstream_centroid_by_jump(pose,jump_id));

	grid_manager->initialize_all_grids(center);
	grid_manager->update_grids(pose,center);

}

}
}
}



