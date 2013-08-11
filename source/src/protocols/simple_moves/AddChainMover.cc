// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file AddChainMover.cc
/// @brief

// Unit headers
#include <protocols/simple_moves/AddChainMover.hh>
#include <protocols/simple_moves/AddChainMoverCreator.hh>
#include <protocols/moves/DataMapObj.hh>
#include <protocols/moves/DataMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/Tracer.hh>
#include <core/pose/util.hh>
using basic::T;
using basic::Error;
using basic::Warning;
static basic::Tracer TR("protocols.simple_moves.AddChainMover");
#include <utility/tag/Tag.hh>

#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>

namespace protocols {
namespace simple_moves {

std::string
AddChainMoverCreator::keyname() const
{
	return AddChainMoverCreator::mover_name();
}

protocols::moves::MoverOP
AddChainMoverCreator::create_mover() const {
	return new AddChainMover;
}

std::string
AddChainMoverCreator::mover_name()
{
	return "AddChain";
}

AddChainMover::AddChainMover()
	: moves::Mover("AddChain"),
	fname_( "" ),
	new_chain_( true )
{
}

void
AddChainMover::apply( Pose & pose )
{
	using namespace core::pose;

	Pose new_pose;

	TR<<"Before addchain, total residues: "<<pose.total_residue()<<std::endl;
	core::import_pose::pose_from_pdb( new_pose, fname() );
	new_pose.conformation().detect_disulfides();
	(*scorefxn()) ( new_pose );

	append_pose_to_pose( pose, new_pose, new_chain() );
	pose.conformation().detect_disulfides();
	pose.update_residue_neighbors();
	(*scorefxn())( pose );
	pose.pdb_info( new core::pose::PDBInfo( pose, true ) ); //reinitialize the PDBInfo
	TR<<"After addchain, total residues: "<<pose.total_residue()<<std::endl;
}

std::string
AddChainMover::get_name() const {
	return AddChainMoverCreator::mover_name();
}

moves::MoverOP
AddChainMover::clone() const
{
	return new AddChainMover( *this );
}

moves::MoverOP
AddChainMover::fresh_instance() const
{
	return new AddChainMover;
}

void
AddChainMover::parse_my_tag(
	utility::tag::TagPtr const tag,
	protocols::moves::DataMap & data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
	fname( tag->getOption< std::string >( "file_name" ) );
	new_chain( tag->getOption< bool >( "new_chain", 1 ) );
	scorefxn( protocols::rosetta_scripts::parse_score_function( tag, data ) );
	TR<<"AddChain sets fname: "<<fname()<<" new_chain: "<<new_chain()<<std::endl;
}

void AddChainMover::scorefxn( core::scoring::ScoreFunctionOP s ){ scorefxn_ = s; }
core::scoring::ScoreFunctionOP AddChainMover::scorefxn() const{ return scorefxn_; }
} // simple_moves
} // protocols
