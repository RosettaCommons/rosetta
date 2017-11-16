// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/protein_interface_design/movers/LoopOver.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

// Unit headers
#include <protocols/protein_interface_design/movers/LoopOver.hh>
#include <protocols/protein_interface_design/movers/LoopOverCreator.hh>
#include <protocols/moves/MoverStatus.hh>
#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>
#include <core/pose/Pose.hh>

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
using namespace protocols::moves;

static basic::Tracer TR( "protocols.protein_interface_design.movers.LoopOver" );

// XRW TEMP std::string
// XRW TEMP LoopOverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return LoopOver::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP LoopOverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new LoopOver );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP LoopOver::mover_name()
// XRW TEMP {
// XRW TEMP  return "LoopOver";
// XRW TEMP }

LoopOver::LoopOver():
	protocols::moves::Mover( LoopOver::mover_name() ),
	max_iterations_( 10 ),
	ms_whenfail_( protocols::moves::MS_SUCCESS )
{}

LoopOver::LoopOver(
	core::Size max_iterations,
	protocols::moves::MoverCOP mover,
	protocols::filters::FilterCOP condition,
	protocols::moves::MoverStatus ms_whenfail
):
	protocols::moves::Mover( LoopOver::mover_name() ),
	max_iterations_( max_iterations ),
	mover_( mover->clone() ),
	condition_( condition->clone() ),
	ms_whenfail_( ms_whenfail )
{}

LoopOver::~LoopOver() {}

void
LoopOver::apply( core::pose::Pose & pose )
{
	core::Size count( 0 );
	core::pose::Pose const saved_pose( pose );

	// clear out any existing info (from prior LoopOver::apply calls)
	info().clear();

	//fpd
	if ( mover_->get_additional_output() ) {
		utility_exit_with_message("Movers returning multiple poses are unsupported by LoopOver.");
	}

	bool filter_result( false );
	while ( !filter_result && count < max_iterations_ ) {
		if ( !drift_ ) {
			pose = saved_pose; // to prevent drift
		}
		TR<<"Loop iteration "<<count<<std::endl;
		// clear out any existing info (from prior mover_->apply calls)
		mover_->info().clear();
		mover_->apply( pose );
		// inherit Mover info: jd2 JobDistributor passes this info to Job,
		// and JobOutputters may then write this info to output files
		info().insert( info().end(), mover_->info().begin(), mover_->info().end() );

		// condition is evaluated only when mover has the status, MS_SUCCESS
		if ( mover_->get_last_move_status() == protocols::moves::MS_SUCCESS ) {
			filter_result = condition_->apply( pose );
		}
		++count;
	}
	if ( !filter_result ) {
		set_last_move_status( ms_whenfail_ );
		if ( !drift_ ) pose = saved_pose;
	} else {
		set_last_move_status( protocols::moves::MS_SUCCESS );
	}
}

// XRW TEMP std::string
// XRW TEMP LoopOver::get_name() const {
// XRW TEMP  return LoopOver::mover_name();
// XRW TEMP }

void
LoopOver::parse_my_tag( TagCOP const tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &filters, protocols::moves::Movers_map const &movers, core::pose::Pose const & )
{
	using namespace filters;

	TR<<"loop\n";
	drift_ = tag->getOption< bool >( "drift", 1 );
	std::string const mover_name( tag->getOption< std::string >( "mover_name" ));
	std::string const filter_name( tag->getOption< std::string >( "filter_name", "false_filter" ) );
	max_iterations_ = tag->getOption< core::Size >( "iterations", 10 );

	std::string const mover_status( tag->getOption< std::string >( "ms_whenfail", "MS_SUCCESS" ));
	ms_whenfail_ = protocols::moves::mstype_from_name( mover_status );

	Movers_map::const_iterator find_mover( movers.find( mover_name ) );
	Filters_map::const_iterator find_filter( filters.find( filter_name ));
	if ( find_mover == movers.end() ) {
		TR<<"WARNING WARNING!!! mover not found in map. skipping:\n"<<tag<<std::endl;
		runtime_assert( find_mover != movers.end() );
	}
	if ( find_filter == filters.end() ) {
		TR<<"WARNING WARNING!!! filter not found in map. skipping:\n"<<tag<<std::endl;
		runtime_assert( find_filter == filters.end() );
	}
	mover_ = find_mover->second; // no cloning to allow other movers to change this mover at their parse time
	condition_ = find_filter->second->clone();
	TR << "with mover \"" << mover_name << "\" and filter \"" << filter_name << "\" and  " << max_iterations_<< " max iterations\n";
	TR.flush();
}

std::string LoopOver::get_name() const {
	return mover_name();
}

std::string LoopOver::mover_name() {
	return "LoopOver";
}

void LoopOver::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	XMLSchemaRestriction mover_status;
	mover_status.name( "mover_status" );
	mover_status.base_type( xs_string );
	mover_status.add_restriction( xsr_enumeration, "MS_SUCCESS" );
	mover_status.add_restriction( xsr_enumeration, "FAIL_RETRY" );
	mover_status.add_restriction( xsr_enumeration, "FAIL_DO_NOT_RETRY" );
	mover_status.add_restriction( xsr_enumeration, "FAIL_BAD_INPUT" );
	mover_status.add_restriction( xsr_enumeration, "FAIL" );
	xsd.add_top_level_element( mover_status );

	attlist + XMLSchemaAttribute::attribute_w_default( "drift", xsct_rosetta_bool, "Avoid repeatedly restoring from a saved pose (which is intended to prevent drift)", "1" )
		+ XMLSchemaAttribute( "mover_name", xs_string, "Mover to use" )
		+ XMLSchemaAttribute( "filter_name", xs_string, "Filter to use" )
		+ XMLSchemaAttribute::attribute_w_default( "iterations", xsct_non_negative_integer, "Number of iterations", "10" )
		+ XMLSchemaAttribute::attribute_w_default( "ms_whenfail", "mover_status", "Mover status to emit upon failure", "MS_SUCCESS" );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TO DO", attlist );
}

std::string LoopOverCreator::keyname() const {
	return LoopOver::mover_name();
}

protocols::moves::MoverOP
LoopOverCreator::create_mover() const {
	return protocols::moves::MoverOP( new LoopOver );
}

void LoopOverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	LoopOver::provide_xml_schema( xsd );
}



} //movers
} //protein_interface_design
} //protocols

