// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SwitchChainOrderMover.cc
/// @brief

// Unit headers
#include <protocols/simple_moves/SwitchChainOrderMover.hh>
#include <protocols/simple_moves/SwitchChainOrderMoverCreator.hh>
#include <protocols/moves/DataMapObj.hh>
#include <protocols/moves/DataMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <basic/Tracer.hh>
#include <core/pose/util.hh>
using basic::T;
using basic::Error;
using basic::Warning;
static basic::Tracer TR("protocols.simple_moves.SwitchChainOrderMover");
#include <utility/tag/Tag.hh>

// AUTO-REMOVED #include <core/chemical/AtomType.hh>
#include <utility/vector1.hh>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>

namespace protocols {
namespace simple_moves {

std::string
SwitchChainOrderMoverCreator::keyname() const
{
	return SwitchChainOrderMoverCreator::mover_name();
}

protocols::moves::MoverOP
SwitchChainOrderMoverCreator::create_mover() const {
	return new SwitchChainOrderMover;
}

std::string
SwitchChainOrderMoverCreator::mover_name()
{
	return "SwitchChainOrder";
}

SwitchChainOrderMover::SwitchChainOrderMover()
	: moves::Mover("SwitchChainOrder"),
	residue_numbers_( NULL )
{
}

void
SwitchChainOrderMover::apply( Pose & pose )
{
	core::pose::Pose new_pose(pose);
	core::kinematics::FoldTree new_ft;
	new_ft.clear();
	core::conformation::Conformation const conf( pose.conformation() );
	core::Size chain_count( 1 );
	utility::vector1< core::Size > new_residue_numbers;
	new_residue_numbers.clear();
	utility::vector1< core::Size > positions_in_new_pose;
	positions_in_new_pose.clear();
	foreach( char const chaini, chain_order() ){
		core::Size const chain( chaini - '0' );
		core::Size const chain_begin( conf.chain_begin( chain ) );
		core::Size const chain_end(   conf.chain_end(   chain ) );

		core::Size const new_chain_begin( positions_in_new_pose.size() + 1 );
		for( core::Size i = chain_begin; i<=chain_end; ++i )
			positions_in_new_pose.push_back( i );
		core::Size const new_chain_end( positions_in_new_pose.size() );
		if( residue_numbers_() != NULL ){
			foreach( core::Size const residue_number, residue_numbers_->obj ){
				if( residue_number <= chain_begin && residue_number >= chain_end )
					new_residue_numbers.push_back( residue_number - ( chain_begin - new_chain_begin ) );
			}
		}
		new_ft.add_edge( new_chain_begin, new_chain_end, -1 );
		if( chain_count > 1 )
			new_ft.add_edge( 1, new_chain_begin, chain_count - 1 );
		chain_count++;
	}
	new_ft.reorder( 1 );
	core::pose::create_subpose( pose, positions_in_new_pose, new_ft, new_pose );
	new_pose.update_residue_neighbors();
	pose.clear();
	pose = new_pose;
	pose.conformation().detect_disulfides();
	pose.update_residue_neighbors();
	TR<<"New pose's foldtree "<<pose.fold_tree()<<std::endl;
	if( residue_numbers_() != NULL ){
		residue_numbers_->obj = new_residue_numbers;
		TR<<"new residue numbers: ";
		foreach( core::Size const res, residue_numbers_->obj )
			TR<<res<<", ";
		TR<<std::endl;
	}
}

std::string
SwitchChainOrderMover::get_name() const {
	return SwitchChainOrderMoverCreator::mover_name();
}

moves::MoverOP
SwitchChainOrderMover::clone() const
{
	return new SwitchChainOrderMover( *this );
}

moves::MoverOP
SwitchChainOrderMover::fresh_instance() const
{
	return new SwitchChainOrderMover;
}

void
SwitchChainOrderMover::parse_my_tag(
	utility::tag::TagPtr const tag,
	protocols::moves::DataMap & data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
	chain_order( tag->getOption< std::string >( "chain_order" ) );
	std::string const residue_numbers_setter( tag->getOption< std::string >( "residue_numbers_setter", "" ) );
	if( data.has( "residue_number_list", residue_numbers_setter ) ){
		residue_numbers_ = data.get< protocols::moves::DataMapObj< utility::vector1< core::Size > > * >( "residue_numbers", "" );
		TR<<"SwitchChainOrderMover will change the residue numbers listed under "<<residue_numbers_setter<<std::endl;
	}
}

void
SwitchChainOrderMover::chain_order( std::string const co ){
	chain_order_ = co;
}

std::string
SwitchChainOrderMover::chain_order() const
{
	return chain_order_;
}

} // simple_moves
} // protocols
