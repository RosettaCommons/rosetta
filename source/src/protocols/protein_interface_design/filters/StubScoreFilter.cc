// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @author Sarel Fleishman (sarelf@uw.edu)
#include <protocols/protein_interface_design/filters/StubScoreFilter.hh>
#include <protocols/protein_interface_design/filters/StubScoreFilterCreator.hh>
#include <core/pose/Pose.hh>
#include <protocols/hotspot_hashing/HotspotStub.hh>
#include <protocols/hotspot_hashing/HotspotStubSet.hh>
#include <utility/tag/Tag.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <numeric/xyzVector.hh>


#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <protocols/protein_interface_design/movers/PlacementMinimizationMover.hh>
#include <protocols/protein_interface_design/movers/PlaceUtils.hh>

#include <core/kinematics/Jump.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <protocols/score_filters/ScoreTypeFilter.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>


namespace protocols {
namespace protein_interface_design {
namespace filters {

static basic::Tracer TR( "protocols.protein_interface_design.filters.StubScoreFilter" );

/// @brief default ctor
StubScoreFilter::StubScoreFilter() :
	parent( "StubScore" ),
	host_chain_( 2 ),
	cb_force_( 0.5 )
{}

bool
StubScoreFilter::apply(core::pose::Pose const & pose ) const
{
	core::Real const stub_score( compute( pose ) );
	if ( stub_score >= -0.0001 ) {
		TR<<"bb_cst evalutes to 0. Failing"<<std::endl;
		return false;
	}
	return true;
}

core::Real
StubScoreFilter::compute( core::pose::Pose const & in_pose ) const{
	if ( !stub_sets_.size() ) {
		TR.Error<<"Stubsets not set in StubScoreFilter. Have I been parsed correctly?"<<std::endl;
		runtime_assert( stub_sets_.size() );
	}
	core::scoring::ScoreFunctionCOP stub_scorefxn( protocols::protein_interface_design::movers::make_stub_scorefxn() );
	core::pose::Pose pose( in_pose );
	(*stub_scorefxn)(pose);//for constraints to be active
	protocols::hotspot_hashing::remove_hotspot_constraints_from_pose( pose );
	protocols::protein_interface_design::movers::PlacementMinimizationMover dummy_min;
	dummy_min.stub_sets( stub_sets_ );
	dummy_min.host_chain( host_chain_ );
	dummy_min.cb_force( cb_force_ );
	dummy_min.refresh_bbstub_constraints( pose );

	protocols::score_filters::ScoreTypeFilter const stf( stub_scorefxn, core::scoring::backbone_stub_constraint, 1.0 );
	core::Real const stub_score( stf.compute( pose ) );
	return( stub_score );
}

core::Real
StubScoreFilter::report_sm( core::pose::Pose const & pose ) const
{
	return( compute( pose ) );
}

void
StubScoreFilter::report( std::ostream & out, core::pose::Pose const & pose ) const
{
	out<<"StubScoreFilter returns "<<compute( pose )<<std::endl;
}

void
StubScoreFilter::stub_sets( utility::vector1<  std::pair< protocols::hotspot_hashing::HotspotStubSetOP, std::pair< protocols::hotspot_hashing::HotspotStubOP, core::Size > > > const & stubsets ){
	stub_sets_ = stubsets;
}

void
StubScoreFilter::parse_my_tag( utility::tag::TagCOP tag,
	basic::datacache::DataMap &data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & pose )
{
	TR.Info << "StubScoreFilter"<<std::endl;
	host_chain_ = tag->getOption< core::Size >( "chain_to_design", 2 );
	cb_force_ = tag->getOption< core::Real >( "cb_force", 0.5 );
	runtime_assert( cb_force_ > -0.00001 );
	stub_sets_ = protocols::protein_interface_design::movers::parse_stub_sets( tag, pose, host_chain_, data );
	runtime_assert( stub_sets_.size() );
}

protocols::filters::FilterOP
StubScoreFilter::fresh_instance() const{
	return protocols::filters::FilterOP( new StubScoreFilter() );
}

StubScoreFilter::~StubScoreFilter()= default;


protocols::filters::FilterOP
StubScoreFilter::clone() const{
	return protocols::filters::FilterOP( new StubScoreFilter( *this ) );
}

// XRW TEMP protocols::filters::FilterOP
// XRW TEMP StubScoreFilterCreator::create_filter() const { return protocols::filters::FilterOP( new StubScoreFilter ); }

// XRW TEMP std::string
// XRW TEMP StubScoreFilterCreator::keyname() const { return "StubScore"; }

std::string StubScoreFilter::name() const {
	return class_name();
}

std::string StubScoreFilter::class_name() {
	return "StubScore";
}

void StubScoreFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default( "chain_to_design", xsct_non_negative_integer, "Chain that ought to be designed, numbered sequentially from 1", "2" )
		+ XMLSchemaAttribute::attribute_w_default( "cb_force", xsct_non_negative_integer, "Chain that ought to be designed, numbered sequentially from 1", "2" );

	XMLSchemaSimpleSubelementList ssl;
	protein_interface_design::movers::add_subelement_for_parse_stub_sets( ssl, xsd );

	protocols::filters::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, class_name(), "This filter checks whether in the current configuration the scaffold is 'feeling' any of the hotspot stub constraints. This is useful for quick triaging of hopeless configuration. If used in conjunction with the Mover PlaceSimultaneously, don't bother setting flags -- the Mover will take care of it.", attlist, ssl );
}

std::string StubScoreFilterCreator::keyname() const {
	return StubScoreFilter::class_name();
}

protocols::filters::FilterOP
StubScoreFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new StubScoreFilter );
}

void StubScoreFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	StubScoreFilter::provide_xml_schema( xsd );
}



} // filters
} // protein_interface_design
} // protocols
