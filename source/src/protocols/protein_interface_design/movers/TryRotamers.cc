// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


/// @file protocols/protein_interface_design/movers/TryRotamers.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

// Unit headers
#include <protocols/protein_interface_design/movers/TryRotamers.hh>
#include <protocols/protein_interface_design/movers/TryRotamersCreator.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <utility/tag/Tag.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/kinematics/Jump.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/util.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/filters/BasicFilters.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/protein_interface_design/util.hh>
#include <utility/string_util.hh>

#include <utility/graph/Graph.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>

#include <core/pose/variant_util.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

using namespace core;
using namespace std;
using namespace core::scoring;
using namespace protocols::moves;

static THREAD_LOCAL basic::Tracer TR( "protocols.protein_interface_design.movers.TryRotamers" );

// XRW TEMP std::string
// XRW TEMP TryRotamersCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return TryRotamers::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP TryRotamersCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new TryRotamers );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP TryRotamers::mover_name()
// XRW TEMP {
// XRW TEMP  return "TryRotamers";
// XRW TEMP }


TryRotamers::TryRotamers() :
	protocols::moves::Mover( TryRotamers::mover_name() )
{ }

TryRotamers::TryRotamers( std::string const & resnum,
	ScoreFunction const & scorefxn,
	protocols::filters::Filter const & final_filter,
	Size explosion, // rotamer explosion
	Size jump_num,
	bool clash_check,
	bool solo_res,
	bool include_current
) :
	protocols::moves::Mover( TryRotamers::mover_name() ),
	scorefxn_(scorefxn.clone()),
	resnum_(resnum),
	jump_num_(jump_num),
	clash_check_(clash_check),
	solo_res_(solo_res),
	include_current_(include_current),
	explosion_(explosion),
	final_filter_(final_filter.clone())
{}

/// @note Pass everything through the final filter (True Filter)
TryRotamers::TryRotamers( std::string const & resnum,
	core::scoring::ScoreFunction const& scorefxn,
	Size explosion, // rotamer explosion
	Size jump_num,
	bool solo_res,
	bool clash_check,
	bool include_current
) :
	protocols::moves::Mover(),
	scorefxn_( scorefxn.clone() ),
	resnum_(resnum),
	jump_num_(jump_num),
	clash_check_(clash_check),
	solo_res_(solo_res),
	include_current_(include_current),
	explosion_(explosion),
	final_filter_(protocols::filters::FilterOP( new protocols::filters::TrueFilter ))
{}

TryRotamers::~TryRotamers() {}

void
TryRotamers::setup_rotamer_set( pose::Pose & pose, core::Size resnum )
{
	using namespace core::pack::task;
	using namespace core::pack::rotamer_set;

	PackerTaskOP ptask( TaskFactory::create_packer_task( pose ) );
	ptask->set_bump_check( clash_check_ );

	ResidueLevelTask & restask( ptask->nonconst_residue_task( resnum ) );
	utility::graph::GraphOP packer_graph( new utility::graph::Graph( pose.size() ) );
	restask.or_ex1( true );
	restask.or_ex2( true );
	restask.or_ex3( true );
	restask.or_ex4( true );
	if ( explosion_ > 0 ) restask.or_ex1_sample_level(core::pack::task::EX_FOUR_HALF_STEP_STDDEVS);
	if ( explosion_ > 1 ) restask.or_ex2_sample_level(core::pack::task::EX_FOUR_HALF_STEP_STDDEVS);
	if ( explosion_ > 2 ) restask.or_ex3_sample_level(core::pack::task::EX_FOUR_HALF_STEP_STDDEVS);
	if ( explosion_ > 3 ) restask.or_ex4_sample_level(core::pack::task::EX_FOUR_HALF_STEP_STDDEVS);
	restask.or_include_current( include_current_ );

	restask.restrict_to_repacking();
	RotamerSetFactory rsf;
	rotset_ = rsf.create_rotamer_set( pose.residue( resnum ) );

	rotset_->set_resid( resnum );
	rotset_->build_rotamers( pose, *scorefxn_, *ptask, packer_graph, false );
	rotamer_it_ = rotset_->begin();
	TR<<"building rotamer set of " <<rotset_->num_rotamers()<< " different rotamers ...\n";
}

void
TryRotamers::apply ( pose::Pose & pose )
{
	//using namespace rotamer_set;
	using namespace core::scoring;

	core::kinematics::FoldTree const saved_ft( pose.fold_tree() );

	TR << "current fold-tree:\n" << pose.fold_tree() << std::endl;

	if ( automatic_connection_ ) {
		core::kinematics::FoldTree const new_ft( make_hotspot_foldtree( pose ) );
		TR<<"New foldtree:\n"<<new_ft<<std::endl;
		pose.fold_tree( new_ft );
	}

	if ( shove_residues_ ) {
		for ( core::Size const resid : core::select::get_residues_from_subset( shove_residues_->apply(pose) ) ) {
			core::pose::add_variant_type_to_pose_residue( pose, core::chemical::SHOVE_BB, resid );
		}
	}

	if ( solo_res_ ) {
		core::pose::add_lower_terminus_type_to_pose_residue( pose, pose.size() ); //prolly critical so that the dunbrack library uses neutral phi
		core::pose::add_upper_terminus_type_to_pose_residue( pose, pose.size() );
	}

	pose.update_residue_neighbors();

	core::Size resnum( core::pose::parse_resnum( resnum_, pose ) );

	if ( !rotset_ || rotset_->num_rotamers() == 0 ) {
		setup_rotamer_set(pose, resnum);
	}//end building rotamer set

	// job distributor iterates ...
	if ( rotamer_it_ == rotset_->end() ) {
		set_last_move_status( protocols::moves::FAIL_DO_NOT_RETRY );
		TR<<"reached the end of the rotamer ensemble "  <<std::endl;
		pose.fold_tree( saved_ft );
		return;
	}

	if ( final_filter_->apply( pose ) ) {
		if ( jump_num_ > 0 ) { // Let 0 indicate no jumps
			core::kinematics::Jump const saved_jump( pose.jump( jump_num_ ) );
			pose.replace_residue ( resnum, **rotamer_it_, false/*orient bb*/ );
			pose.set_jump( jump_num_, saved_jump );
		} else { // no jumps to save
			pose.replace_residue ( resnum, **rotamer_it_, false/*orient bb*/ );
		}
		TR << "TryRotamers passed the final_filter" <<std::endl;
		set_last_move_status ( protocols::moves::MS_SUCCESS );
		++rotamer_it_;
		pose.fold_tree( saved_ft );
		return;
	} else {
		set_last_move_status( protocols::moves::FAIL_RETRY );
		++rotamer_it_;
		TR << std::endl;
	}
	pose.fold_tree( saved_ft );
}

// XRW TEMP std::string
// XRW TEMP TryRotamers::get_name() const {
// XRW TEMP  return TryRotamers::mover_name();
// XRW TEMP }

void
TryRotamers::parse_my_tag( TagCOP const tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const &filters,
	Movers_map const &,
	Pose const & )
{
	resnum_ = core::pose::get_resnum_string( tag );
	automatic_connection_ = tag->getOption< bool >( "automatic_connection", 1 );
	scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, data )->clone();
	jump_num_ = tag->getOption<core::Size>( "jump_num", 1);
	solo_res_ = tag->getOption<bool>( "solo_res", 0);
	std::string const final_filter_name( tag->getOption<std::string>( "final_filter", "true_filter" ) );
	protocols::filters::Filters_map::const_iterator find_filter( filters.find( final_filter_name ));

	clash_check_ = tag->getOption<bool>("clash_check", 0 );
	include_current_ = tag->getOption<bool>("include_current", 1 );

	explosion_ = tag->getOption<core::Size>( "explosion", 0);
	bool const filter_found( find_filter != filters.end() );
	if ( filter_found ) {
		final_filter_ = find_filter->second->clone();
	} else {
		if ( final_filter_name != "true_filter" ) {
			TR<<"***WARNING WARNING! Filter defined for TryRotamers not found in filter_list!!!! Defaulting to truefilter***"<<std::endl;
			runtime_assert( filter_found );
		} else {
			final_filter_ = protocols::filters::FilterOP( new protocols::filters::TrueFilter );
		}
	}

	if ( tag->hasOption( "shove" ) ) {
		shove_residues_ = core::pose::get_resnum_selector( tag, "shove" );
	}


	TR
		<< "TryRotamers was instantiated using scorefxn=" << rosetta_scripts::get_score_function_name(tag)
		<< ", jump_number=" << jump_num_
		<< ", solo_res=" << solo_res_
		<< ", clash_check=" << clash_check_
		<< ", include_current=" << include_current_
		<< ", and explosion=" << explosion_ << std::endl;
}

std::string TryRotamers::get_name() const {
	return mover_name();
}

std::string TryRotamers::mover_name() {
	return "TryRotamers";
}

void TryRotamers::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	core::pose::attributes_for_get_resnum_string( attlist );
	attlist + XMLSchemaAttribute::attribute_w_default( "automatic_connection", xsct_rosetta_bool, "Automatically set up a hotspot foldtree", "1" );
	rosetta_scripts::attributes_for_parse_score_function( attlist );
	attlist + XMLSchemaAttribute::attribute_w_default( "jump_num", xsct_non_negative_integer, "Jump across which to evaluate energies, numbered sequentially from one", "1" )
		+ XMLSchemaAttribute::attribute_w_default( "solo_res", xsct_rosetta_bool, "Add terminal types (upper and lower) to the hotspot residue", "1" )
		+ XMLSchemaAttribute::attribute_w_default( "final_filter", xs_string, "Add terminal types (upper and lower) to the hotspot residue", "true_filter" )
		+ XMLSchemaAttribute::attribute_w_default( "clash_check", xsct_rosetta_bool, "Perform a clash check and filter out clashing rotamers (helps with speed!)", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "include_current", xsct_rosetta_bool, "Include input conformation", "1" )
		+ XMLSchemaAttribute::attribute_w_default( "explosion", xsct_non_negative_integer, "Aggression with which to sample conformations, 0-4 (values more than four are equivalent thereto)", "0" )
		+ XMLSchemaAttribute( "shove", xsct_refpose_enabled_residue_number_cslist, "List of residues for which to use the shove type, given in refpose, pose, or PDB numbering" );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TO DO", attlist );
}

std::string TryRotamersCreator::keyname() const {
	return TryRotamers::mover_name();
}

protocols::moves::MoverOP
TryRotamersCreator::create_mover() const {
	return protocols::moves::MoverOP( new TryRotamers );
}

void TryRotamersCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	TryRotamers::provide_xml_schema( xsd );
}


} //movers
} //protein_interface_design
} //protocols

