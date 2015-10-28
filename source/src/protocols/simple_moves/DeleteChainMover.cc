// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file DeleteChainMover.cc
/// @brief

// Unit headers
#include <protocols/simple_moves/DeleteChainMover.hh>
#include <protocols/simple_moves/DeleteChainMoverCreator.hh>

// Core headers
#include <core/pose/util.hh>
#include <core/pose/Pose.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

namespace protocols {
namespace simple_moves {

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.DeleteChainMover" );

std::string
DeleteChainMoverCreator::keyname() const
{
	return DeleteChainMoverCreator::mover_name();
}

protocols::moves::MoverOP
DeleteChainMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new DeleteChainMover );
}

std::string
DeleteChainMoverCreator::mover_name()
{
	return "DeleteChain";
}

std::string
DeleteChainMover::get_name() const {
	return DeleteChainMoverCreator::mover_name();
}

moves::MoverOP
DeleteChainMover::clone() const
{
	return moves::MoverOP( new DeleteChainMover( *this ) );
}

moves::MoverOP
DeleteChainMover::fresh_instance() const
{
	return moves::MoverOP( new DeleteChainMover );
}

DeleteChainMover::DeleteChainMover():
	moves::Mover("DeleteChain"),
	chain_num_( 1 )
{}

///////////////END BOILER PLATE CODE////////////////

void
DeleteChainMover::chain_num(
	core::Size chain_num
) {
	chain_num_ = chain_num;
}

core::Size
DeleteChainMover::chain_num() {
	return chain_num_;
}

void
DeleteChainMover::apply( Pose & pose )
{
	utility::vector1<core::pose::PoseOP> chain_poses = pose.split_by_chain();
	TR << "Removing chain " << chain_num() << " from pose with " << chain_poses.size() << " chains" << std::endl;

	bool chain_added = false;
	for(core::Size i=1; i<=chain_poses.size(); ++i) {
		if(i != chain_num_){
			if(!chain_added) {
				pose = *chain_poses[i];
				chain_added = true;
			}
			else {
				pose.append_pose_by_jump(*chain_poses[i], 1);
			}
		}
	}
}


void
DeleteChainMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & ,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
	if(tag->hasOption("chain")) {
		chain_num_ = tag->getOption< core::Size >( "chain");
	}
}

} // simple_moves
} // protocols
