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
#include <protocols/rosetta_scripts/util.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

// Basic headers
#include <basic/citation_manager/CitationCollection.hh>
#include <basic/citation_manager/UnpublishedModuleInfo.hh>
#include <basic/citation_manager/CitationManager.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

using namespace core;
using namespace std;
using namespace protocols::moves;

static basic::Tracer TR( "protocols.protein_interface_design.movers.LoopOver" );




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

LoopOver::~LoopOver() = default;

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


void
LoopOver::parse_my_tag( TagCOP const tag, basic::datacache::DataMap & data )
{
	using namespace filters;

	TR<<"loop\n";
	drift_ = tag->getOption< bool >( "drift", true );
	std::string const mover_name( tag->getOption< std::string >( "mover_name" ));
	std::string const filter_name( tag->getOption< std::string >( "filter_name", "false_filter" ) );
	max_iterations_ = tag->getOption< core::Size >( "iterations", 10 );

	std::string const mover_status( tag->getOption< std::string >( "ms_whenfail", "MS_SUCCESS" ));
	ms_whenfail_ = protocols::moves::mstype_from_name( mover_status );

	mover_ = protocols::rosetta_scripts::parse_mover( mover_name, data ); // no cloning to allow other movers to change this mover at their parse time
	condition_ = protocols::rosetta_scripts::parse_filter( filter_name, data )->clone();
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

/// @brief This mover provides citation info (returns true).
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
bool
LoopOver::mover_provides_citation_info() const {
	return true;
}

/// @brief Provide the citation (Fleishman 2011).
/// @returns A vector of citation collections for this mover and any published movers and filters that it invokes.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
utility::vector1< basic::citation_manager::CitationCollectionCOP >
LoopOver::provide_citation_info() const {
	using namespace basic::citation_manager;
	CitationCollectionOP mycitation( utility::pointer::make_shared<CitationCollection>( get_name(), CitedModuleType::Mover ));
	mycitation->add_citation( CitationManager::get_instance()->get_citation_by_doi( "10.1371/journal.pone.0020161" ) );
	utility::vector1< CitationCollectionCOP > returnvec { mycitation };

	//Add citations for modules that this mover invokes:
	if ( mover_ != nullptr ) {
		merge_into_citation_collection_vector( mover_->provide_citation_info(), returnvec );
	}
	if ( condition_ != nullptr ) {
		merge_into_citation_collection_vector( condition_->provide_citation_info(), returnvec );
	}

	return returnvec;
}

/// @brief Provides unpublished citation info for any unpublished movers and filters that this mover invokes.
/// @returns A list of pairs of (author, e-mail address).  Empty list if none are unpublished.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
utility::vector1< basic::citation_manager::UnpublishedModuleInfoCOP >
LoopOver::provide_authorship_info_for_unpublished() const {
	using namespace basic::citation_manager;
	utility::vector1< UnpublishedModuleInfoCOP > returnvec;

	//Add citations for modules that this mover invokes:
	if ( mover_ != nullptr ) {
		merge_into_unpublished_collection_vector( mover_->provide_authorship_info_for_unpublished(), returnvec );
	}
	if ( condition_ != nullptr ) {
		merge_into_unpublished_collection_vector( condition_->provide_authorship_info_for_unpublished(), returnvec );
	}

	return returnvec;
}

std::string LoopOverCreator::keyname() const {
	return LoopOver::mover_name();
}

protocols::moves::MoverOP
LoopOverCreator::create_mover() const {
	return utility::pointer::make_shared< LoopOver >();
}

void LoopOverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	LoopOver::provide_xml_schema( xsd );
}



} //movers
} //protein_interface_design
} //protocols

