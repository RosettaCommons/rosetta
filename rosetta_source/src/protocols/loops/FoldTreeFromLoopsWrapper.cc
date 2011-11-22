// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/loops/FoldTreeFromLoops.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu)

// Unit headers
#include <protocols/loops/FoldTreeFromLoopsWrapper.hh>
#include <protocols/loops/FoldTreeFromLoopsWrapperCreator.hh>
// Package headers
#include <core/pose/Pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/loops/loops_main.hh>
#include <utility/vector1.hh>
#include <protocols/loops/util.hh>

//Auto Headers
#include <utility/vector0.hh>

namespace protocols {
namespace loops {

using namespace std;
using namespace core::kinematics;

static basic::Tracer TR( "protocols.loops.FoldTreeFromLoopsWrapper" );

std::string
FoldTreeFromLoopsCreator::keyname() const
{
	return FoldTreeFromLoopsCreator::mover_name();
}

protocols::moves::MoverOP
FoldTreeFromLoopsCreator::create_mover() const {
	return new FoldTreeFromLoops;
}

std::string
FoldTreeFromLoopsCreator::mover_name()
{
	return "FoldTreeFromLoops";
}

FoldTreeFromLoops::FoldTreeFromLoops() :
	Mover( FoldTreeFromLoopsCreator::mover_name() ), loop_str_( "" )
{
	loops_.clear();
}


FoldTreeFromLoops::~FoldTreeFromLoops() {}

void
FoldTreeFromLoops::apply( core::pose::Pose & pose )
{
	if( loops().empty() )
		loops( loops_from_string( loop_str(), pose ) );
	FoldTree f;
	fold_tree_from_loops( pose, loops(), f );
	TR<<"old foldtree "<<pose.fold_tree()<<"\nNew foldtree ";
	pose.fold_tree( f );
	TR<<pose.fold_tree()<<std::endl;
}

std::string
FoldTreeFromLoops::get_name() const {
	return FoldTreeFromLoopsCreator::mover_name();
}

void
FoldTreeFromLoops::parse_my_tag( TagPtr const tag, protocols::moves::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & )
{
	using namespace protocols::rosetta_scripts;

	loop_str( tag->getOption< string >( "loops" ) );

	TR<<"FoldTreeFromLoops with loops "<<loop_str()<<std::endl;
}

} //loops
} //protocols
