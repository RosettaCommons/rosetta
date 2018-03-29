// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/filters/ConservedPosMutationFilter.cc
/// @author Florian Richter (floric@u.washington.edu), may 2011

// Unit Headers
#include <protocols/calc_taskop_filters/ConservedPosMutationFilter.hh>
#include <protocols/calc_taskop_filters/ConservedPosMutationFilterCreator.hh>

// Package Headers
#include <protocols/task_operations/SeqprofConsensusOperation.hh>

// Project Headers
#include <basic/Tracer.hh>

#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.hh>
#include <core/pose/Pose.hh>
#include <core/sequence/SequenceProfile.hh>
#include <core/types.hh>


// Utility headers
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>


//// C++ headers
static basic::Tracer tr( "protocols.simple_filters.ConservedPosMutationFilter" );

namespace protocols {
namespace calc_taskop_filters {


ConservedPosMutationFilter::ConservedPosMutationFilter() :
	parent("ConservedPosMutationFilter"),
	conserved_pos_taskop_(task_operations::RestrictConservedLowDdgOperationOP( new task_operations::RestrictConservedLowDdgOperation() )),
	max_allowed_conserved_pos_mutations_(0)
{
}

ConservedPosMutationFilter::~ConservedPosMutationFilter()= default;

filters::FilterOP
ConservedPosMutationFilter::clone() const {
	return filters::FilterOP( new ConservedPosMutationFilter( *this ) ); }

filters::FilterOP
ConservedPosMutationFilter::fresh_instance() const {
	return filters::FilterOP( new ConservedPosMutationFilter() ); }


bool
ConservedPosMutationFilter::apply( core::pose::Pose const & pose ) const {

	core::Size conserved_pos_mutations(0);
	bool verbose( conserved_pos_taskop_->verbose() );
	std::string mutstring;

	for ( core::Size i(1); i <= pose.size(); ++i ) {
		if ( !pose.residue_type(i).is_protein() ) continue;

		core::chemical::AA wt_aa( conserved_pos_taskop_->seqprof_wt_aa(i) );
		if ( conserved_pos_taskop_->position_untouchable(i, wt_aa) && (wt_aa != pose.residue_type(i).aa() ) ) {
			conserved_pos_mutations++;
			if ( !verbose && (conserved_pos_mutations > max_allowed_conserved_pos_mutations_) ) {
				tr << "Pose has at least " << conserved_pos_mutations << " mutations at conserved positions, but only " << max_allowed_conserved_pos_mutations_ << " are allowed, returnig false..." << std::endl;
				return false;
			} else if ( verbose ) {
				using namespace core::chemical;
				mutstring = mutstring + oneletter_code_from_aa( wt_aa ) + utility::to_string( i ) + oneletter_code_from_aa( pose.residue_type(i).aa() ) + ", ";
			}
		}
		//if the user is interested, we'll write out the conservation and ddg info for every mutation
		if ( verbose && (wt_aa != pose.residue_type(i).aa() ) ) {
			tr << "Mutation at pos " << i << " with ala_ddg of " << conserved_pos_taskop_->position_ala_ddG( i ) << " and wt conservation " << conserved_pos_taskop_->seqprof()->profile()[ i ][ conserved_pos_taskop_->seqprof_wt_aa(i) ]  << std::endl;
		}
	}
	if ( verbose ) tr << "Forbidden mutations detected: " << mutstring << std::endl;
	if ( conserved_pos_mutations > max_allowed_conserved_pos_mutations_ ) {
		tr << "Pose has " << conserved_pos_mutations << " mutations at conserved positions, but only " << max_allowed_conserved_pos_mutations_ << " are allowed, returnig false..." << std::endl;
		return false;
	}

	tr << "Pose has " << conserved_pos_mutations << " mutations at conserved positions, <= than the allowed value of " << max_allowed_conserved_pos_mutations_ << ", returnig true..." << std::endl;
	return true;
} // apply


void
ConservedPosMutationFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & datamap, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const &  )
{
	conserved_pos_taskop_->parse_tag( tag, datamap );

	if ( tag->hasOption("max_conserved_pos_mutations") ) {
		max_allowed_conserved_pos_mutations_ = tag->getOption<core::Size>("max_conserved_pos_mutations");
	}
}

std::string ConservedPosMutationFilter::name() const {
	return class_name();
}

std::string ConservedPosMutationFilter::class_name() {
	return "ConservedPosMutationFilter";
}

void ConservedPosMutationFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	auto ct_gen = protocols::task_operations::
		RestrictConservedLowDdgOperation::create_complex_type_generator(xsd);

	attlist + XMLSchemaAttribute(
		"max_conerved_pos_mutations", xsct_non_negative_integer,
		"Maximum number of allowed mutations of conserved positions.");

	ct_gen->add_attributes(attlist).complex_type_naming_func( filters::complex_type_name_for_filter )
		.description(
		"Returns true if the given pose passes the filter, false otherwise. "
		"In this case, a pose passes if it less than max_allowed_conserved_pos_mutations_ "
		"mutations at conserved position. The decision whether a given position "
		"counts as conserved is made by the RestrictConservedLowDdgOperation "
		"task operation to prevent duplication of code")
		.element_name(class_name())
		.write_complex_type_to_schema(xsd);
}

std::string ConservedPosMutationFilterCreator::keyname() const {
	return ConservedPosMutationFilter::class_name();
}

protocols::filters::FilterOP
ConservedPosMutationFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new ConservedPosMutationFilter );
}

void ConservedPosMutationFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ConservedPosMutationFilter::provide_xml_schema( xsd );
}



} // filters
} // protocols
