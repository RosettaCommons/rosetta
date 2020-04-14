// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file SwitchChainOrderMover.cc
/// @brief

// Unit headers
#include <protocols/simple_moves/SwitchChainOrderMover.hh>
#include <protocols/simple_moves/SwitchChainOrderMoverCreator.hh>
#include <basic/datacache/DataMapObj.hh>
#include <basic/datacache/DataMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <basic/Tracer.hh>
#include <core/pose/subpose_manipulation_util.hh>
#include <core/pose/util.hh>
#include <core/pose/chains_util.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <utility/tag/Tag.hh>

#include <utility/vector1.hh>
#include <utility/string_util.hh>
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
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


static basic::Tracer TR( "protocols.simple_moves.SwitchChainOrderMover" );

namespace protocols {
namespace simple_moves {




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
	utility::vector1< core::Size > positions_in_new_pose; // positions_in_new_pose[ new_pose_seqpos ] = old_pose_seqpos
	positions_in_new_pose.clear();
	for ( core::Size chain: chain_order(pose) ) {
		TR<<"Now at chain: "<<chain<<std::endl;
		runtime_assert( chain > 0 && chain <= conf.num_chains() );
		core::Size const chain_begin( conf.chain_begin( chain ) );
		core::Size const chain_end(   conf.chain_end(   chain ) );

		core::Size const new_chain_begin( positions_in_new_pose.size() + 1 );
		for ( core::Size i = chain_begin; i<=chain_end; ++i ) {
			positions_in_new_pose.push_back( i );
		}
		core::Size const new_chain_end( positions_in_new_pose.size() );
		if ( residue_numbers_ != nullptr ) {
			for ( core::Size const residue_number : residue_numbers_->obj ) {
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
	new_pose.pdb_info( utility::pointer::make_shared< core::pose::PDBInfo >( new_pose, true ) ); //reinitialize the PDBInfo

	// Salvage old PDBInfo information. Feel free to add anything else you need here.
	for ( core::Size new_seqpos = 1; new_seqpos <= new_pose.size(); new_seqpos++ ) {
		core::Size old_seqpos = positions_in_new_pose.at( new_seqpos );
		for ( auto reslabel : pose.pdb_info()->get_reslabels( old_seqpos ) ) {
			new_pose.pdb_info()->add_reslabel( new_seqpos, reslabel );
		}
	}

	//When applying switch then comments are erased from the pose. adding condition that if -pdb comments true flag is turned on then copy comments to new pose. gideonla 1/5/13
	if ( basic::options::option[ basic::options::OptionKeys::out::file::pdb_comments ].value() ) {
		std::map< std::string, std::string > const comments = core::pose::get_all_comments( pose );
		for ( auto const & comment : comments ) {
			core::pose::add_comment(new_pose,comment.first,comment.second);
		}
	}
	pose.clear();
	pose = new_pose;
	pose.conformation().detect_disulfides();
	( *scorefxn() ) ( pose );
	pose.update_residue_neighbors();
	TR<<"New pose's foldtree "<<pose.fold_tree()<<std::endl;
	if ( residue_numbers_ != nullptr ) {
		residue_numbers_->obj = new_residue_numbers;
		TR<<"new residue numbers: ";
		for ( core::Size const res : residue_numbers_->obj ) {
			TR<<res<<", ";
		}
		TR<<std::endl;
	}
}


moves::MoverOP
SwitchChainOrderMover::clone() const
{
	return utility::pointer::make_shared< SwitchChainOrderMover >( *this );
}

moves::MoverOP
SwitchChainOrderMover::fresh_instance() const
{
	return utility::pointer::make_shared< SwitchChainOrderMover >();
}

void
SwitchChainOrderMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data
)
{
	if ( tag->hasOption("chain_order") && (tag->hasOption("chain_name") || tag->hasOption("chain_num")) )  {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "You can specify a chain_order string or the comma separated chain names or numbers, but not both");
	}

	if ( tag->hasOption("chain_order") ) {
		chain_order( tag->getOption< std::string >( "chain_order" ) );
	} else {
		if ( tag->hasOption("chain_num") ) {
			chain_ids_ = utility::string_split(tag->getOption<std::string>("chain_num"),',');
		}

		if ( tag->hasOption("chain_name") ) {
			for ( auto & chain_name : utility::string_split(tag->getOption<std::string>("chain_name"),',') ) {
				chain_ids_.push_back( chain_name );
			}
		}
	}

	invert_chains_ = tag->getOption<bool>("invert_chains", false);

	std::string const residue_numbers_setter( tag->getOption< std::string >( "residue_numbers_setter", "" ) );
	scorefxn( protocols::rosetta_scripts::parse_score_function( tag, data ) );
	if ( residue_numbers_setter != "" ) {
		residue_numbers_ = basic::datacache::get_set_from_datamap< basic::datacache::DataMapObj< utility::vector1< core::Size > > >( "residue_numbers", residue_numbers_setter, data );
	}

}

void
SwitchChainOrderMover::chain_order( std::string const & co ) {
	chain_ids_.clear();

	for ( char const chaini : co ) {
		chain_ids_.push_back( std::string(1,chaini) );
	}
}

utility::vector1< core::Size >
SwitchChainOrderMover::chain_order( core::pose::Pose const & pose ) const
{
	utility::vector1< core::Size > chain_order;

	for ( std::string const & chaini: chain_ids_ ) {
		if ( chaini.size() == 1 && chaini[0] > '0' && chaini[0] <= '9' ) {
			chain_order.push_back( chaini[0] - '0' );
		} else {
			chain_order.push_back( core::pose::get_chain_id_from_chain( chaini, pose ) );
		}
	}

	if ( invert_chains_ ) {
		std::reverse( chain_order.begin(), chain_order.end() );
	}

	return chain_order;
}

void
SwitchChainOrderMover::scorefxn( core::scoring::ScoreFunctionOP s ) { scorefxn_ = s; }

core::scoring::ScoreFunctionOP
SwitchChainOrderMover::scorefxn() const { return scorefxn_; }

std::string SwitchChainOrderMover::get_name() const {
	return mover_name();
}

std::string SwitchChainOrderMover::mover_name() {
	return "SwitchChainOrder";
}

void SwitchChainOrderMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute( "chain_order", xs_string, "Order of final chains" )
		+ XMLSchemaAttribute( "chain_num", xsct_nnegative_int_cslist, "List of final chain positions, to be correlated to chain_name" )
		+ XMLSchemaAttribute( "chain_name", xsct_chain_cslist, "List of final chain names, to be correlated to chain_num" )
		+ XMLSchemaAttribute::attribute_w_default( "invert_chains", xsct_rosetta_bool, "If true, apply on the inverse of chain selection (as specified in name or number options)", "false" )
		+ XMLSchemaAttribute( "residue_numbers_setter", xs_string, "List of final chain names, to be correlated to chain_num" );
	rosetta_scripts::attributes_for_parse_score_function( attlist );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TO DO", attlist );
}

std::string SwitchChainOrderMoverCreator::keyname() const {
	return SwitchChainOrderMover::mover_name();
}

protocols::moves::MoverOP
SwitchChainOrderMoverCreator::create_mover() const {
	return utility::pointer::make_shared< SwitchChainOrderMover >();
}

void SwitchChainOrderMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SwitchChainOrderMover::provide_xml_schema( xsd );
}


} // simple_moves
} // protocols
