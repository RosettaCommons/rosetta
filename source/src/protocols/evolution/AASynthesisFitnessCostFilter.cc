// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/evolution/AASynthesisFitnessCostFilter.cc
/// @brief
/// @author Christoffer Norn


//Unit Headers
#include <protocols/evolution/AASynthesisFitnessCostFilter.hh>
#include <protocols/evolution/AASynthesisFitnessCostFilterCreator.hh>
#include <utility/tag/Tag.hh>
#include <map>

//Project Headers
#include <basic/Tracer.hh>
#include <core/pose/Pose.hh>
//#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/util.hh>

// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

namespace protocols {
namespace evolution {

static basic::Tracer TR( "protocols.evolution.AASynthesisFitnessCost" );

// XRW TEMP protocols::filters::FilterOP
// XRW TEMP AASynthesisFitnessCostFilterCreator::create_filter() const { return protocols::filters::FilterOP( new AASynthesisFitnessCost ); }

// XRW TEMP std::string
// XRW TEMP AASynthesisFitnessCostFilterCreator::keyname() const { return "AASynthesisFitnessCost"; }

//default ctor
AASynthesisFitnessCost::AASynthesisFitnessCost() :
	protocols::filters::Filter( "AASynthesisFitnessCost" ),
	threshold_( 1 ),
	fitnessCostPerATP_( 0.00001 )
{
}

AASynthesisFitnessCost::~AASynthesisFitnessCost() = default;

void
AASynthesisFitnessCost::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &, moves::Movers_map const &, core::pose::Pose const &)
{
	fitnessCostPerATP( tag->getOption< core::Real >( "fitnessCostPerATP", 0.00001 ) );
	threshold( tag->getOption< core::Size >( "threshold", 1 ) );
}

bool
AASynthesisFitnessCost::apply( core::pose::Pose const & pose ) const {
	core::Real val = compute( pose );
	TR << "The cost of synthesis yields a relative fitness reduction of " << val << " mutations." << std::endl;
	return( compute( pose ) <= threshold() );
}

void
AASynthesisFitnessCost::report( std::ostream & o, core::pose::Pose const & pose ) const {
	bool const val = ( compute( pose ) >= threshold() );
	o << "AASynthesisFitnessCost returns " << val << std::endl;
}

core::Real
AASynthesisFitnessCost::report_sm( core::pose::Pose const & pose ) const {
	return( core::Real (compute( pose )) );
}

core::Real
AASynthesisFitnessCost::compute(
	core::pose::Pose const & pose
) const {
	// Setup aa to aa cost map
	char aminoAcids[20] = {'V','L','I','M','F','D','N','E','Q','R','K','H','A','G','S','T','P','C','Y','W'};
	core::Real synCost[20] = {23.3,27.3,32.3,34.3,52,12.7,14.7,15.3,16.3,27.3,30.3,38.3,11.7,11.7,11.7,18.7,20.3,24.7,50,74.3};
	std::map< char, core::Real > aa2cost;
	for ( int i=0; i<20; i++ ) {
		aa2cost[ aminoAcids[i] ] = synCost[i];
	}

	// Calculate the sum cost of the sequence
	core::Real sum_cost = 0;

	std::string pose_sequence( "" );
	for ( core::Size chaini = 1 ; chaini <= pose.conformation().num_chains(); ++chaini ) {
		pose_sequence += pose.chain_sequence( chaini );
	}

	for ( char& c : pose_sequence ) {
		sum_cost += aa2cost[c];
	}

	// Convert this to a fitness penalty
	core::Real relativeFitnessLoss = 1 - ( fitnessCostPerATP() * sum_cost );
	if ( relativeFitnessLoss < 0 ) {
		relativeFitnessLoss = 0; // we don't want this to go negative
	}

	TR << "The relative fitness loss due to metabolic cost of synthesis is " << relativeFitnessLoss << std::endl;
	return relativeFitnessLoss;
}

std::string AASynthesisFitnessCost::name() const {
	return class_name();
}

std::string AASynthesisFitnessCost::class_name() {
	return "AASynthesisFitnessCost";
}

void AASynthesisFitnessCost::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default("fitnessCostPerATP", xsct_real, "Sets the fitness loss cost per metabolic cost (think ATP usage) for synthesis", "0.00001")
		+ XMLSchemaAttribute::attribute_w_default(
		"threshold", xsct_non_negative_integer,
		"threshold",
		"1000");

	protocols::filters::xsd_type_definition_w_attributes(
		xsd, class_name(),
		"Measures the fitness cost of an amino acid sequence",
		attlist );
}

std::string AASynthesisFitnessCostFilterCreator::keyname() const {
	return AASynthesisFitnessCost::class_name();
}

protocols::filters::FilterOP
AASynthesisFitnessCostFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new AASynthesisFitnessCost );
}

void AASynthesisFitnessCostFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	AASynthesisFitnessCost::provide_xml_schema( xsd );
}


}
}
