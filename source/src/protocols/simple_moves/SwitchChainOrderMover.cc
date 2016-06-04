// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
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
#include <basic/datacache/DataMapObj.hh>
#include <basic/datacache/DataMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <basic/Tracer.hh>
#include <core/pose/util.hh>
using basic::T;
using basic::Error;
using basic::Warning;
static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.SwitchChainOrderMover" );
#include <utility/tag/Tag.hh>

#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <boost/foreach.hpp>
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
//options Includes
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh> // for option[ out::file::silent  ] and etc.
#include <basic/options/keys/in.OptionKeys.gen.hh> // for option[ in::file::tags ] and etc.
#include <basic/options/keys/OptionKeys.hh>


namespace protocols {
namespace simple_moves {

std::string
SwitchChainOrderMoverCreator::keyname() const
{
	return SwitchChainOrderMoverCreator::mover_name();
}

protocols::moves::MoverOP
SwitchChainOrderMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new SwitchChainOrderMover );
}

std::string
SwitchChainOrderMoverCreator::mover_name()
{
	return "SwitchChainOrder";
}

SwitchChainOrderMover::SwitchChainOrderMover()
: moves::Mover("SwitchChainOrder"),
	residue_numbers_( /* NULL */ )
{
	scorefxn( core::scoring::get_score_function() );
}

void
SwitchChainOrderMover::apply( Pose & pose )
{
	core::pose::Pose new_pose(pose);
	core::kinematics::FoldTree new_ft;
	new_ft.clear();
	core::conformation::Conformation const & conf = pose.conformation();
	TR<<"Number of chains in pose: "<<conf.num_chains()<<std::endl;
	core::Size chain_count( 1 );
	utility::vector1< core::Size > new_residue_numbers;
	new_residue_numbers.clear();
	utility::vector1< core::Size > positions_in_new_pose;
	positions_in_new_pose.clear();
	BOOST_FOREACH ( char const chaini, chain_order() ) {
		core::Size const chain( chaini - '0' );
		TR<<"Now at chain: "<<chain<<std::endl;
		runtime_assert( chain > 0 && chain <= conf.num_chains() );
		core::Size const chain_begin( conf.chain_begin( chain ) );
		core::Size const chain_end(   conf.chain_end(   chain ) );

		core::Size const new_chain_begin( positions_in_new_pose.size() + 1 );
		for ( core::Size i = chain_begin; i<=chain_end; ++i ) {
			positions_in_new_pose.push_back( i );
		}
		core::Size const new_chain_end( positions_in_new_pose.size() );
		if ( residue_numbers_ != NULL ) {
			BOOST_FOREACH ( core::Size const residue_number, residue_numbers_->obj ) {
				if ( residue_number >= chain_begin && residue_number <= chain_end ) {
					new_residue_numbers.push_back( residue_number - ( chain_begin - new_chain_begin ) );
				}
			}
		}
		new_ft.add_edge( new_chain_begin, new_chain_end, -1 );
		if ( chain_count > 1 ) {
			new_ft.add_edge( 1, new_chain_begin, chain_count - 1 );
		}
		chain_count++;
	}
	new_ft.reorder( 1 );
	core::pose::create_subpose( pose, positions_in_new_pose, new_ft, new_pose );
	new_pose.update_residue_neighbors();
	new_pose.pdb_info( core::pose::PDBInfoOP( new core::pose::PDBInfo( new_pose, true ) ) ); //reinitialize the PDBInfo

	//When applying switch then comments are erased from the pose. adding condition that if -pdb comments true flag is turned on then copy comments to new pose. gideonla 1/5/13
	if ( basic::options::option[ basic::options::OptionKeys::out::file::pdb_comments ].value() ) {
		std::map< std::string, std::string > const comments = core::pose::get_all_comments( pose );
		for ( std::map< std::string, std::string >::const_iterator i = comments.begin(); i != comments.end(); ++i ) {
			core::pose::add_comment(new_pose,i->first,i->second);
		}
	}
	pose.clear();
	pose = new_pose;
	pose.conformation().detect_disulfides();
	( *scorefxn() ) ( pose );
	pose.update_residue_neighbors();
	TR<<"New pose's foldtree "<<pose.fold_tree()<<std::endl;
	if ( residue_numbers_ != NULL ) {
		residue_numbers_->obj = new_residue_numbers;
		TR<<"new residue numbers: ";
		BOOST_FOREACH ( core::Size const res, residue_numbers_->obj ) {
			TR<<res<<", ";
		}
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
	return moves::MoverOP( new SwitchChainOrderMover( *this ) );
}

moves::MoverOP
SwitchChainOrderMover::fresh_instance() const
{
	return moves::MoverOP( new SwitchChainOrderMover );
}

void
SwitchChainOrderMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & pose )
{
	if ( tag->hasOption("chain_order") && (tag->hasOption("chain_name") || tag->hasOption("chain_num")) )  {
		throw utility::excn::EXCN_RosettaScriptsOption("You can specify a chain_order string or the comma separated chain names or numbers, but not both");
	}

	if ( tag->hasOption("chain_order") ) {
		chain_order( tag->getOption< std::string >( "chain_order" ) );
	} else {
		if ( tag->hasOption("chain_num") ) {
			chain_ids_ = utility::string_split(tag->getOption<std::string>("chain_num"),',',core::Size());
		}

		if ( tag->hasOption("chain_name") ) {
			utility::vector1<std::string> chain_names = utility::string_split(tag->getOption<std::string>("chain_name"),',',std::string());
			for ( utility::vector1<std::string>::iterator chain_name_it = chain_names.begin(); chain_name_it != chain_names.end(); ++chain_name_it ) {
				chain_ids_.push_back(core::pose::get_chain_id_from_chain(*chain_name_it,pose));
			}
		}
	}

	if ( tag->getOption<bool>("invert_chains", 0) ) {
		// Invert the chain selection
		utility::vector1<core::Size> inverted_chain_ids_;
		for ( core::Size i=1 ; i <= pose.conformation().num_chains() ; i++ ) {
			if ( std::find(chain_ids_.begin(), chain_ids_.end(), i) == chain_ids_.end() ) {
				inverted_chain_ids_.push_back(i);
			}
		}
		chain_ids_ = inverted_chain_ids_;
	}

	std::string const residue_numbers_setter( tag->getOption< std::string >( "residue_numbers_setter", "" ) );
	scorefxn( protocols::rosetta_scripts::parse_score_function( tag, data ) );
	if ( residue_numbers_setter != "" ) {
		residue_numbers_ = basic::datacache::get_set_from_datamap< basic::datacache::DataMapObj< utility::vector1< core::Size > > >( "residue_numbers", residue_numbers_setter, data );
	}

}

void
SwitchChainOrderMover::chain_order( std::string const & co ){
	utility::vector1<core::Size> new_chain_ids;
	BOOST_FOREACH ( char const chaini, co ) {
		core::Size const chain( chaini - '0' );
		new_chain_ids.push_back( chain );
	}
	chain_ids_ = new_chain_ids;
}

std::string
SwitchChainOrderMover::chain_order() const
{
	std::string chain_order;
	BOOST_FOREACH ( core::Size const chain, chain_ids_ ) {
		chain_order += utility::to_string<core::Size>(chain);
	}
	return chain_order;
}

void
SwitchChainOrderMover::scorefxn( core::scoring::ScoreFunctionOP s ) { scorefxn_ = s; }

core::scoring::ScoreFunctionOP
SwitchChainOrderMover::scorefxn() const { return scorefxn_; }

} // simple_moves
} // protocols
