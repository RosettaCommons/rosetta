// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @author Sarel Fleishman (sarelf@uw.edu)
#include <protocols/protein_interface_design/filters/StubScoreLoopsFilter.hh>
#include <protocols/protein_interface_design/filters/StubScoreLoopsFilterCreator.hh>
#include <core/pose/Pose.hh>
#include <protocols/hotspot_hashing/HotspotStub.hh>
#include <protocols/hotspot_hashing/HotspotStubSet.hh>
#include <utility/tag/Tag.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>

#include <protocols/protein_interface_design/movers/SetupHotspotConstraintsLoopsMover.hh>

#include <core/kinematics/Jump.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <protocols/jobdist/Jobs.hh>
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

static basic::Tracer tr( "protocols.protein_interface_design.filters.StubScoreLoopsFilter" );

/// @brief default ctor
StubScoreLoopsFilter::StubScoreLoopsFilter() :
	cb_force_( 0.5 )
{
	set_score_type( core::scoring::backbone_stub_constraint );
}

StubScoreLoopsFilter::~StubScoreLoopsFilter() = default;

void
StubScoreLoopsFilter::parse_my_tag( utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & pose )
{
	tr.Info << "StubScoreLoopsFilter"<<std::endl;
	cb_force_ = tag->getOption< core::Real >( "cb_force", 0.5 );
	runtime_assert( cb_force_ > 0.00001 );
	stub_set_ = hotspot_hashing::HotspotStubSetOP( new hotspot_hashing::HotspotStubSet );
	stub_set_->read_data( tag->getOption< std::string >("stubfile") );
	loop_start_ = tag->getOption<Size>("start", 0 );
	loop_stop_ = tag->getOption<Size>("stop", 0 );
	if ( loop_start_ <= 0 ) {
		utility_exit_with_message( "please provide loop-start with 'start' option" );
	}
	if ( loop_stop_ <= loop_start_ ) {
		utility_exit_with_message( "loop-stop has to be larger than loop-start. assign with 'stop'" );
	}
	protein_interface_design::movers::SetupHotspotConstraintsLoopsMover hspmover( stub_set_ );
	hspmover.set_loop_start( loop_start_ );
	hspmover.set_loop_stop( loop_stop_ );
	resfile_ = tag->getOption< std::string >("resfile","NONE");
	if ( resfile_ != "NONE" ) hspmover.set_resfile( resfile_ );
	core::scoring::constraints::ConstraintCOPs constraints;
	core::Size ncst = hspmover.generate_csts( pose, constraints );
	set_constraints( constraints );
	tr.Info << "Filter with " << ncst << " hotspots in " << constraints.size() << " constraints." << std::endl;
}

protocols::filters::FilterOP
StubScoreLoopsFilter::fresh_instance() const{
	return protocols::filters::FilterOP( new StubScoreLoopsFilter() );
}


protocols::filters::FilterOP
StubScoreLoopsFilter::clone() const{
	return protocols::filters::FilterOP( new StubScoreLoopsFilter( *this ) );
}

// XRW TEMP protocols::filters::FilterOP
// XRW TEMP StubScoreLoopsFilterCreator::create_filter() const { return protocols::filters::FilterOP( new StubScoreLoopsFilter ); }

// XRW TEMP std::string
// XRW TEMP StubScoreLoopsFilterCreator::keyname() const { return "StubScoreLoops"; }

std::string StubScoreLoopsFilter::name() const {
	return class_name();
}

std::string StubScoreLoopsFilter::class_name() {
	return "StubScoreLoops";
}

void StubScoreLoopsFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	attlist + XMLSchemaAttribute::attribute_w_default( "cb_force", xsct_non_negative_integer, "Chain that ought to be designed, numbered sequentially from 1", "2" )
		+ XMLSchemaAttribute::required_attribute( "start", xsct_positive_integer, "Starting residue of loop" )
		+ XMLSchemaAttribute::required_attribute( "stop", xsct_positive_integer, "Ending residue of loop" )
		+ XMLSchemaAttribute::required_attribute( "stubfile", xs_string, "Stub file filename" )
		+ XMLSchemaAttribute::attribute_w_default( "resfile", xs_string, "Resfile filename", "NONE" );

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "This filter checks whether in the current configuration the scaffold is 'feeling' any of the hotspot stub constraints. This is useful for quick triaging of hopeless configuration. If used in conjunction with the Mover PlaceSimultaneously, don't bother setting flags -- the Mover will take care of it.", attlist );
}

std::string StubScoreLoopsFilterCreator::keyname() const {
	return StubScoreLoopsFilter::class_name();
}

protocols::filters::FilterOP
StubScoreLoopsFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new StubScoreLoopsFilter );
}

void StubScoreLoopsFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	StubScoreLoopsFilter::provide_xml_schema( xsd );
}



} // filters
} // protein_interface_design
} // protocols
