// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/protein_interface_design/movers/FavorNativeResiduePreCycle.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)
#include <protocols/protein_interface_design/movers/FavorNativeResiduePreCycle.hh>
#include <protocols/protein_interface_design/movers/FavorNativeResiduePreCycleCreator.hh>
#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>

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

static basic::Tracer TR( "protocols.protein_interface_design.movers.FavorNativeResiduePreCycle" );

// XRW TEMP std::string FavorNativeResiduePreCycleCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return FavorNativeResiduePreCycle::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP FavorNativeResiduePreCycleCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new FavorNativeResiduePreCycle );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP FavorNativeResiduePreCycle::mover_name() {
// XRW TEMP  return "FavorNativeResidue";
// XRW TEMP }


// XRW TEMP std::string
// XRW TEMP FavorNativeResiduePreCycle::get_name() const {
// XRW TEMP  return "FavorNativeResidue";
// XRW TEMP }

FavorNativeResiduePreCycle::~FavorNativeResiduePreCycle() = default;

void
FavorNativeResiduePreCycle::parse_my_tag( TagCOP const tag, basic::datacache::DataMap & data, protocols::filters::Filters_map const &, Movers_map const &, core::pose::Pose const & )
{
	using namespace utility::pointer;

	bonus_ = tag->getOption<core::Real>( "bonus", 1.5 );
	for ( std::map< std::string, ReferenceCountOP >::const_iterator it = (data)[ "scorefxns" ].begin(); it!=(data)[ "scorefxns" ].end(); ++it ) {
		ScoreFunctionOP scorefxn( data.get_ptr< ScoreFunction >( "scorefxns", it->first ) );
		if ( scorefxn->get_weight( res_type_constraint ) == 0.0 ) {
			scorefxn->set_weight( res_type_constraint, bonus_ );
			TR<<"Setting res_type_constraint weight in scorefxn "<<it->first<<" to "<<bonus_<<'\n';
		}
	}
	/*
	for( std::map< std::string, ReferenceCountOP >::const_iterator it=(data)[ "scorefxns_hshash" ].begin(); it!=(data)[ "scorefxns" ].end(); ++it ){ // scorefxns where the user defined hs_hash but not fnr
	ScoreFunctionOP scorefxn( *data.get< ScoreFunction * >( "scorefxns", it->first ) );
	scorefxn->set_weight( res_type_constraint, bonus_ );
	TR<<"Setting res_type_constraint weight in scorefxn "<<it->first<<" to "<<bonus_<<'\n';
	}
	*/
	TR<<"applying favor native residue to pose with weight: "<<bonus_<<std::endl;
}

std::string FavorNativeResiduePreCycle::get_name() const {
	return mover_name();
}

std::string FavorNativeResiduePreCycle::mover_name() {
	return "FavorNativeResidue";
}

void FavorNativeResiduePreCycle::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default( "bonus", xsct_real, "Bonus for the native residue", "1.5" );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TO DO", attlist );
}

std::string FavorNativeResiduePreCycleCreator::keyname() const {
	return FavorNativeResiduePreCycle::mover_name();
}

protocols::moves::MoverOP
FavorNativeResiduePreCycleCreator::create_mover() const {
	return protocols::moves::MoverOP( new FavorNativeResiduePreCycle );
}

void FavorNativeResiduePreCycleCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	FavorNativeResiduePreCycle::provide_xml_schema( xsd );
}


} //movers
} //protein_interface_design
} //protocols

