// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   SplitAndMixPoseMover.cc
/// @brief  Splits a Pose according to a ResidueSelector
/// @author jaumebonet (jaume.bonet@gmail.com), Correia's LPDI/EPFL


#include <protocols/fold_from_loops/SplitAndMixPoseMover.hh>
#include <protocols/fold_from_loops/SplitAndMixPoseMoverCreator.hh>

#include <utility/string_util.hh>
#include <utility/vector1.functions.hh>
#include <core/pose/Pose.hh>
#include <core/pose/subpose_manipulation_util.hh>
#include <core/pose/variant_util.hh>
#include <core/pose/chains_util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/select/residue_selector/util.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/NotResidueSelector.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <core/select/residue_selector/ResidueRanges.hh>

#include <protocols/hybridization/util.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/mover_schemas.hh>
#include <protocols/grafting/simple_movers/DeleteRegionMover.hh>
#include <protocols/fold_from_loops/SplitAndMixPoseMover.hh>

#include <utility/vector1.hh>

// XSD XRW Includes
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/Filter.hh>

namespace protocols {
namespace fold_from_loops {

static basic::Tracer TR( "protocols.fold_from_loops.SplitAndMixPoseMover", basic::t_trace );

SplitAndMixPoseMover::SplitAndMixPoseMover():
	selector_( nullptr ),
	ranges_( new core::select::residue_selector::ResidueRanges ),
	order_( utility::vector1< core::Size >() ),
	merge_chains_( default_merge_chains() ),
	pairs_( utility::vector1< SubPoseIDPair >() ),
	try_all_pairs_( default_try_all_pairs() ),
	exclude_consecutive_( default_exclude_consecutive() ),
	max_distance_( default_max_distance() ),
	global_i_( 1 )
{}

/// @brief Re-merge order Setter
void
SplitAndMixPoseMover::set_order( std::string order )
{
	reset_order();
	if ( order.size() > 0 ) {
		utility::vector1< std::string > const str_residues( utility::string_split( order, ',' ) );
		for ( auto const & res : str_residues ) {
			if ( utility::string2int( res ) != -1 ) {
				order_.push_back( utility::string2Size( res ) );
			} else {
				throw CREATE_EXCEPTION(utility::excn::Exception,  "Value " + res + " cannot be cast to integer." );
			}
		}
	}
}

/// @brief Calculate Number of Segments
core::Size
SplitAndMixPoseMover::count_segments( core::pose::Pose const & pose )
{
	TR << "count_segments" << std::endl;
	if ( count_segments() == 0 ) {
		core::select::residue_selector::ResidueRangesOP ranges( new core::select::residue_selector::ResidueRanges );
		TR << "apply selector to pose size " << pose.size() << std::endl;
		ranges->from_subset( selector_->apply( pose ) );
		TR << "Total ranges from selector: " << ranges->size() << std::endl;
		fix_ranges( pose, *ranges );
	}
	TR << "count_segments - end" << std::endl;
	return count_segments();
}

void
SplitAndMixPoseMover::apply( core::pose::Pose & pose )
{
	if ( selector_ != nullptr ) {
		count_segments( pose );
		utility::vector1< core::pose::PoseOP > pose_list = split_pose( pose );
		if ( pose_list.size() > 0 ) {

			fill_pair_list( pose_list.size() );
			adapt_order_to_pair( pose_list );

			if ( order_.size() == 0 ) set_order( pose_list.size() );
			if ( utility::max( order_ ) > pose_list.size() ) {
				throw CREATE_EXCEPTION(utility::excn::Exception,  "A subpose is requested identified with a number higher than the number of subposes available" );
			}
			if ( order_.size() > pose_list.size() ) { //  We could use less than all, never more.
				throw CREATE_EXCEPTION(utility::excn::Exception,  "Requested ordes defines a number of parts different that those provided" );
			}
			core::pose::PoseOP tmppose = merge_poses( pose_list );
			transfer_conformation( pose, *tmppose );
		} else {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "Selector does not match with any residue of the provided pose" );
		}
	}
	pose.pdb_info( core::pose::PDBInfoOP( new core::pose::PDBInfo( pose ) ) );
	TR << pose.annotated_sequence() << std::endl;
	TR << pose.fold_tree() << std::endl;
}

void
SplitAndMixPoseMover::fix_ranges( core::pose::Pose const & pose, core::select::residue_selector::ResidueRanges const & ranges )
{
	core::select::residue_selector::ResidueRangesOP tmpranges( new core::select::residue_selector::ResidueRanges );
	utility::vector1<core::Size> edges;
	for ( auto const & range : ranges ) {
		edges.push_back( range.start() );
		for ( core::Size i=range.start() + 1; i<range.stop(); ++i ) { // this will split both chain jumps and discontinued segments
			if ( hybridization::discontinued_lower( pose, i ) || hybridization::discontinued_upper( pose, i ) ) {
				edges.push_back( i );
			}
		}
		edges.push_back( range.stop() );
		if ( edges.size() % 2 != 0 ) {
			edges.push_back( range.stop() );
		}
	}
	TR << "Found edges:" << edges << std::endl;
	for ( core::Size i=1; i<=edges.size(); i+=2 ) {
		core::select::residue_selector::ResidueRange r( edges[i], edges[i+1] );
		tmpranges->push_back( r );
	}
	ranges_ = tmpranges;
}

utility::vector1< core::pose::PoseOP >
SplitAndMixPoseMover::split_pose( core::pose::Pose const & pose )
{
	using namespace core::select::residue_selector;

	utility::vector1< core::pose::PoseOP > pose_list;
	grafting::simple_movers::DeleteRegionMover deleter;
	TR << "Subposes requested: " << count_segments() << std::endl;
	for ( auto const & range : *ranges_ ) {
		core::pose::Pose copy_pose = pose;
		ResidueSelectorCOP index( new ResidueIndexSelector( range.to_string() ) );
		ResidueSelectorCOP notindex( new NotResidueSelector( index ) );
		deleter.set_residue_selector( notindex );
		TR << "From: " << range.start() << " To: " << range.stop() << std::endl;
		deleter.apply( copy_pose );
		TR << "Extracted segment of length " << copy_pose.size() << std::endl;
		copy_pose.fold_tree( core::kinematics::FoldTree( copy_pose.size() ) );
		TR << "Generated FoldTree " << copy_pose.fold_tree() << std::endl;
		pose_list.push_back( copy_pose.clone() );
	}
	if ( TR.visible() ) {
		for ( core::Size i=1; i<= order_.size(); ++i ) {
			TR << "Segment " << order_[i] << ": " << pose_list[order_[i]]->annotated_sequence()  << std::endl;
		}

	}
	return pose_list;
}

core::pose::PoseOP
SplitAndMixPoseMover::merge_poses( utility::vector1< core::pose::PoseOP > const & pose_list )
{
	if ( pose_list.size() == 1 ) return pose_list[1];

	core::pose::PoseOP pose( new core::pose::Pose );
	std::string chainID, newChainID;
	for ( core::Size i=1; i<= order_.size(); ++i ) {
		core::pose::PoseOP piece = pose_list[order_[i]];
		core::pose::remove_lower_terminus_type_from_pose_residue ( *piece, 1 );
		core::pose::remove_upper_terminus_type_from_pose_residue ( *piece, piece->size() );
		if ( pose->empty() ) {
			pose->detached_copy( *piece );
			chainID = core::pose::get_chain_from_chain_id( 1, *piece );
			pose->fold_tree( core::kinematics::FoldTree( pose->size() ) );
			core::pose::add_lower_terminus_type_to_pose_residue ( *pose, 1 );
		} else {
			newChainID = core::pose::get_chain_from_chain_id( 1, *piece );
			if ( merge_chains_ || ( chainID == newChainID ) ) {
				core::pose::append_pose_to_pose( *pose, *piece,  false );
			} else {
				core::pose::add_upper_terminus_type_to_pose_residue ( *pose , pose->size() );
				core::pose::add_lower_terminus_type_to_pose_residue ( *piece, 1 );
				core::pose::append_pose_to_pose( *pose, *piece,  true );
			}
			chainID = newChainID;
		}
	}
	core::pose::add_upper_terminus_type_to_pose_residue ( *pose , pose->size() );
	return pose;
}

void
SplitAndMixPoseMover::transfer_conformation( core::pose::Pose & pose, core::pose::Pose const & merged )
{
	pose.clear();
	pose.set_new_conformation( merged.conformation_ptr() );
	pose.fold_tree( merged.fold_tree() );
}

void
SplitAndMixPoseMover::fill_pair_list( core::Size total_subposes )
{
	if ( try_all_pairs_ && pairs_.size() == 0 ) {
		core::Size gap = exclude_consecutive_ ? 2 : 1;
		for ( core::Size i = 1; i <= total_subposes; ++i ) {
			for ( core::Size j = i + gap; j <= total_subposes; ++j ) {
				pairs_.push_back( boost::make_tuple( i, j ) );
				pairs_.push_back( boost::make_tuple( j, i ) );
			}
		}
		TR << "Created a total of " << pairs_.size() << " possible combinations" << std::endl;
	}
}

void
SplitAndMixPoseMover::adapt_order_to_pair( utility::vector1< core::pose::PoseOP > const & pose_list )
{
	if ( try_all_pairs_ ) {
		reset_order();
		bool new_case = false;
		while ( !new_case ) {
			utility::vector1< core::Size > wpair;
			wpair.push_back( boost::get<0>(pairs_[global_i_]) );
			wpair.push_back( boost::get<1>(pairs_[global_i_]) );
			TR << "Evaluating Pair " << wpair[1] << " - " << wpair[2] << std::endl;
			increase_global_count();
			if ( NC_distance_filter( *pose_list[wpair[1]], *pose_list[wpair[2]] ) ) {
				TR << "> Trying Pair " << wpair[1] << " - " << wpair[2] << std::endl;
				set_order( wpair );
				new_case = true;
			} else {
				TR << "> Distance limit excludes pair " << wpair[1] << " - " << wpair[2] << std::endl;
			}
		}
	}
}

bool
SplitAndMixPoseMover::NC_distance_filter( core::pose::Pose const & npose, core::pose::Pose const & cpose )
{
	if ( max_distance_ < 0 ) return true;

	core::conformation::Residue r1 = npose.residue( npose.size() );
	core::conformation::Residue r2 = cpose.residue( 1 );
	core::Real distance = r1.xyz( r1.atom_index("CA") ).distance( r2.xyz( r2.atom_index("CA") ) );
	return distance < max_distance_;
}

void
SplitAndMixPoseMover::parse_my_tag( utility::tag::TagCOP tag,
	basic::datacache::DataMap &data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
	set_residue_selector( core::select::residue_selector::parse_residue_selector( tag, data, "residue_selector" ) );
	set_order( tag->getOption< std::string >( "order", std::string() ) );
	set_merge_chains( tag->getOption< bool >( "merge_chains", default_merge_chains() ) );
	set_try_all_pairs( tag->getOption< bool >( "try_all_pairs", default_try_all_pairs() ) );
	set_exclude_consecutive( tag->getOption< bool >( "exclude_consecutive", default_exclude_consecutive() ) );
	set_max_distance( tag->getOption< core::Real >( "max_distance", default_max_distance() ) );
}

void
SplitAndMixPoseMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute( "order", xs_string, "Coma separated numbers defining the orders to rejoin the segments" )
		+ XMLSchemaAttribute::attribute_w_default( "merge_chains",  xsct_rosetta_bool,
		"When true, all the pieces are merged as a single chain. Otherwise the chain of a piece is assigned according to its original chain.", utility::to_string( default_merge_chains() ) )
		+ XMLSchemaAttribute::attribute_w_default( "try_all_pairs",  xsct_rosetta_bool,
		"When true, all the pieces are merged as a single chain", utility::to_string( default_try_all_pairs() ) )
		+ XMLSchemaAttribute::attribute_w_default( "exclude_consecutive",  xsct_rosetta_bool,
		"When true, and try_all_pairs true, it will avoit trying consecutive subposes. Only applies when try_all_pairs=True.", utility::to_string( default_exclude_consecutive() ) )
		+ XMLSchemaAttribute::attribute_w_default( "max_distance",  xsct_real,
		"When positive, and try_all_pairs true, limits the distance between N-term and C-term to generate a pair. Only applies when try_all_pairs=True.", utility::to_string( default_max_distance() ) );
	core::select::residue_selector::attributes_for_parse_residue_selector_when_required( attlist, "residue_selector", "ResidueSelector to define the regions of the pose to keep" );

	moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Splits a Pose into requested pieces and joins them together. Allows to swap order", attlist );
}

std::string
SplitAndMixPoseMoverCreator::keyname() const {
	return SplitAndMixPoseMover::mover_name();
}
protocols::moves::MoverOP
SplitAndMixPoseMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new SplitAndMixPoseMover );
}
void
SplitAndMixPoseMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	SplitAndMixPoseMover::provide_xml_schema( xsd );
}

}
}
