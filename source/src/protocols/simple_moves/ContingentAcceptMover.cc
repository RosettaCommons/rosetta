// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ContingentAcceptMover.cc
/// @brief

// Unit headers
#include <protocols/simple_moves/ContingentAcceptMover.hh>
#include <protocols/simple_moves/ContingentAcceptMoverCreator.hh>
#include <basic/datacache/DataMapObj.hh>
#include <basic/datacache/DataMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <basic/Tracer.hh>
#include <core/pose/util.hh>
using basic::T;
using basic::Error;
using basic::Warning;
static thread_local basic::Tracer TR( "protocols.simple_moves.ContingentAcceptMover" );
#include <utility/tag/Tag.hh>

#include <utility/vector1.hh>
#include <core/pose/Pose.hh>

namespace protocols {
namespace simple_moves {

std::string
ContingentAcceptMoverCreator::keyname() const
{
	return ContingentAcceptMoverCreator::mover_name();
}

protocols::moves::MoverOP
ContingentAcceptMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new ContingentAcceptMover );
}

std::string
ContingentAcceptMoverCreator::mover_name()
{
	return "ContingentAccept";
}

ContingentAcceptMover::ContingentAcceptMover()
    : moves::Mover("ContingentAccept"),
    filter_(/* NULL */),
    mover_(/* NULL */),
    delta_(5)
{
}
    
protocols::filters::FilterOP
ContingentAcceptMover::filter() const{
    return filter_;
}

void
ContingentAcceptMover::filter( protocols::filters::FilterOP f ){
    filter_ = f;
}
    
void
ContingentAcceptMover::apply( Pose & pose )
{
    core::pose::Pose const old_pose = pose;
    core::Real preMoveScore;
    preMoveScore = filter()->report_sm( pose );
    mover()->apply(pose);
    core::Real postMoveScore;
    postMoveScore = filter()->report_sm( pose );
	
    if (postMoveScore <= preMoveScore + delta()){
        TR<<"postMoveScore is "<<postMoveScore<<" preMoveScore is "<<preMoveScore<<" . Accepting new pose."<<std::endl;
	return;
    }
    else{
	TR<<"postMoveScore is "<<postMoveScore<<" preMoveScore is "<<preMoveScore<<" . Rejecting new pose."<<std::endl;
        pose = old_pose;
    }
}

std::string
ContingentAcceptMover::get_name() const {
	return ContingentAcceptMoverCreator::mover_name();
}

moves::MoverOP
ContingentAcceptMover::clone() const
{
	return moves::MoverOP( new ContingentAcceptMover( *this ) );
}

moves::MoverOP
ContingentAcceptMover::fresh_instance() const
{
	return moves::MoverOP( new ContingentAcceptMover );
}

/*
	returns the contains mover
*/
protocols::moves::MoverOP ContingentAcceptMover::mover() const {
	return mover_;
}
  
/*
	Setting the internal mover to point to m.
*/
void ContingentAcceptMover::mover( protocols::moves::MoverOP m ){
	mover_ = m;
}

void
ContingentAcceptMover::parse_my_tag(
                               TagCOP const tag,
                               basic::datacache::DataMap &,
                               protocols::filters::Filters_map const & filters,
                               protocols::moves::Movers_map const & movers,
                               Pose const & )
{
    filter( protocols::rosetta_scripts::parse_filter( tag->getOption< std::string >( "filter" ), filters ) );
    mover( protocols::rosetta_scripts::parse_mover( tag->getOption< std::string >( "mover" ), movers ) );
    delta( tag->getOption< core::Real >( "delta" , 5)); // do I need anything like ( "MaxDeltaFilterVal", 5 )?? What is the integer for?
}

    
    
} // simple_moves
} // protocols

