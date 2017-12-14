// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/fldsgn/filters/StrandHelixGeometryFilter.cc
/// @brief Another filter used in Marcos & Basanta et al. 2017 that needs to be updated.
/// @author Benjamin Basanta (basantab@uw.edu)

#include <protocols/fldsgn/filters/StrandHelixGeometryFilter.hh>
#include <protocols/fldsgn/filters/StrandHelixGeometryFilterCreator.hh>

//XSD includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

// Package Headers
#include <utility/vector0.hh>
#include <utility/vector1.hh>

#include <protocols/fldsgn/topology/SheetFoldTypeManager.hh>
#include <protocols/fldsgn/topology/StrandPairing.hh>
#include <protocols/fldsgn/topology/util.hh>

#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <protocols/parser/BluePrint.hh>
#include <protocols/fldsgn/topology/SS_Info2.hh>



static basic::Tracer TR( "protocols.fldsgn.filters.StrandHelixGeometryFilter" );

namespace protocols {
namespace fldsgn {
namespace filters {

StrandHelixGeometryFilter::StrandHelixGeometryFilter():
	protocols::filters::Filter( "StrandHelixGeometryFilter" ),
	secstruct_( "" ),
	filter_min_orthoangle_(0.0),
	filter_max_orthoangle_(0.0),
	filter_min_planeangle_(0.0),
	filter_max_planeangle_(0.0),
	filter_min_dist_(0.0),
	filter_max_dist_(0.0),
	strand_id1_(0),
	strand_id2_(0),
	helix_id_(0),
	output_type_( "dist" ),
	output_value_( -990.0 )
{}

StrandHelixGeometryFilter::~StrandHelixGeometryFilter() = default;

// @brief set secondary structure
void
StrandHelixGeometryFilter::secstruct( std::string const & ss )
{
	secstruct_ = ss;
}

/// @brief minimum angle for filtering
void
StrandHelixGeometryFilter::filter_min_orthoangle( core::Real const r )
{
	filter_min_orthoangle_ = r;
}

/// @brief maximum angle for filtering
void
StrandHelixGeometryFilter::filter_max_orthoangle( core::Real const r )
{
	filter_max_orthoangle_ = r;
}

/// @brief minimum angle for filtering
void
StrandHelixGeometryFilter::filter_min_planeangle( core::Real const r )
{
	filter_min_planeangle_ = r;
}

/// @brief maximum angle for filtering
void
StrandHelixGeometryFilter::filter_max_planeangle( core::Real const r )
{
	filter_max_planeangle_ = r;
}

/// @brief minimum distance for filtering
void
StrandHelixGeometryFilter::filter_min_dist( core::Real const r )
{
	filter_min_dist_ = r;
}

/// @brief maximum distance for filtering
void
StrandHelixGeometryFilter::filter_max_dist( core::Real const r )
{
	filter_max_dist_ = r;
}
/// @brief Strand id number
void
StrandHelixGeometryFilter::strand_id1( core::Size const r )
{
	strand_id1_ = r;
}

/// @brief Strand id number
void
StrandHelixGeometryFilter::strand_id2( core::Size const r )
{
	strand_id2_ = r;
}

/// @brief Strand id number
void
StrandHelixGeometryFilter::helix_id( core::Size const r )
{
	helix_id_ = r;
}

/// @brief
void
StrandHelixGeometryFilter::output_type( std::string const & s )
{
	output_type_ = s;
}

void
StrandHelixGeometryFilter::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & ,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
	// Blueprint is for giving secondary structure information, otherwise dssp will run for ss definition
	// SSPAIR line is read for the topology of strand pairings
	filter_min_dist_ = tag->getOption<core::Real>( "min_dist", 8 );
	filter_max_dist_ = tag->getOption<core::Real>( "max_dist", 12.0 );
	filter_min_orthoangle_ = tag->getOption<core::Real>( "min_ortho_angle", -180.0 );
	filter_max_orthoangle_ = tag->getOption<core::Real>( "max_ortho_angle", 180.0 );
	filter_min_planeangle_ = tag->getOption<core::Real>( "min_plane_angle", -180.0 );
	filter_max_planeangle_ = tag->getOption<core::Real>( "max_plane_angle", 180.0 );
	strand_id1_ = tag->getOption<core::Size>( "StrandID1", 1 );
	strand_id2_ = tag->getOption<core::Size>( "StrandID2", 2 );
	helix_id_ = tag->getOption<core::Size>( "HelixID", 1 );
	output_type_ = tag->getOption<std::string>( "output_type", "dist" );

	//if( output_type_ != "dist" && output_type_ != "angle" ) {
	//  tr << "Invalid type of output_type, choose either dist or angle. " << std::endl;
	//  }


	std::string const blueprint = tag->getOption<std::string>( "blueprint", "" );
	if ( blueprint != "" ) {
		protocols::parser::BluePrint blue( blueprint );
		secstruct_ = blue.secstruct();
	}
}

protocols::filters::FilterOP
StrandHelixGeometryFilter::clone() const
{
	return protocols::filters::FilterOP( new StrandHelixGeometryFilter( *this ) );
}


protocols::filters::FilterOP
StrandHelixGeometryFilter::fresh_instance() const
{
	return protocols::filters::FilterOP( new StrandHelixGeometryFilter );
}

bool
StrandHelixGeometryFilter::apply( core::pose::Pose const & pose ) const
{
	using core::scoring::dssp::Dssp;
	using protocols::fldsgn::topology::SS_Info2_OP;
	using protocols::fldsgn::topology::NO_STRANDS;

	// set secondary structure
	if ( secstruct_ == "" ) {
		Dssp dssp( pose );
		secstruct_ = dssp.get_dssp_secstruct();
	}
	// set SS_Info
	//SS_Info2_OP  ssinfo = new SS_Info2( pose, secstruct_ ); // old format
	SS_Info2_OP  ssinfo( new topology::SS_Info2( pose, secstruct_ ) );

	TR << "ss of input pose=" << secstruct_ << std::endl;


	bool filter( true );
	core::Real distance = protocols::fldsgn::topology::calc_strand_helix_angle( pose, ssinfo,strand_id1_, strand_id2_, helix_id_ , "dist" );
	core::Real ortho_angle = protocols::fldsgn::topology::calc_strand_helix_angle( pose, ssinfo, strand_id1_, strand_id2_, helix_id_ , "ortho_angle" );
	core::Real plane_angle = protocols::fldsgn::topology::calc_strand_helix_angle( pose, ssinfo, strand_id1_, strand_id2_, helix_id_ , "plane_angle" );
	//core::Real distance = 2;
	//core::Real ortho_angle = 2;
	//core::Real plane_angle = 3;

	if ( ortho_angle < filter_min_orthoangle_ ||  ortho_angle > filter_max_orthoangle_ ) {
		filter = false;
	}
	if ( plane_angle < filter_min_planeangle_ ||  plane_angle > filter_max_planeangle_ ) {
		filter = false;
	}
	if ( distance < filter_min_dist_ || distance > filter_max_dist_ ) {
		filter = false;
	}

	if ( filter ) {
		TR << " Filter success ! " << std::endl;
	} else {
		TR << " Filter failed ! " << std::endl;
	}
	if ( output_type_ == "dist" ) {
		output_value_ = distance;
	} else if ( output_type_ == "plane_angle" ) {
		output_value_ = plane_angle;
	} else if ( output_type_ == "ortho_angle" ) {
		output_value_ = ortho_angle;
	}

	TR << " Ortho angle = " << ortho_angle << std::endl;
	TR << " Plane angle = " << plane_angle << std::endl;
	TR << " Distance = " << distance << std::endl;

	return filter;
}

core::Real
StrandHelixGeometryFilter::report_sm( core::pose::Pose const & ) const
{
	return output_value_;
}

void
StrandHelixGeometryFilter::report( std::ostream &, core::pose::Pose const & ) const
{

}

std::string StrandHelixGeometryFilter::name() const {
	return class_name();
}

std::string StrandHelixGeometryFilter::class_name() {
	return "StrandHelixGeometryFilter";
}

void StrandHelixGeometryFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	/*
	filter_min_dist_ = tag->getOption<core::Real>( "min_dist", 8 );
	filter_max_dist_ = tag->getOption<core::Real>( "max_dist", 12.0 );
	filter_min_orthoangle_ = tag->getOption<core::Real>( "min_ortho_angle", -180.0 );
	filter_max_orthoangle_ = tag->getOption<core::Real>( "max_ortho_angle", 180.0 );
	filter_min_planeangle_ = tag->getOption<core::Real>( "min_plane_angle", -180.0 );
	filter_max_planeangle_ = tag->getOption<core::Real>( "max_plane_angle", 180.0 );
	strand_id1_ = tag->getOption<core::Size>( "StrandID1", 1 );
	strand_id2_ = tag->getOption<core::Size>( "StrandID2", 2 );
	helix_id_ = tag->getOption<core::Size>( "HelixID", 1 );
	output_type_ = tag->getOption<std::string>( "output_type", "dist" );
	*/
	attlist
		+ XMLSchemaAttribute::attribute_w_default( "min_dist"  , xsct_real , "Minimum distance between the centers of the strand pair and the helix." , "8.0" )
		+ XMLSchemaAttribute::attribute_w_default( "max_dist"  , xsct_real , "Maximum distance between the centers of the strand pair and the helix." , "12.0" )
		+ XMLSchemaAttribute::attribute_w_default( "min_ortho_angle"  , xsct_real , "Minimum orthogonal angle between the strand pair and the helix." , "-180.0" )
		+ XMLSchemaAttribute::attribute_w_default( "max_ortho_angle"  , xsct_real , "Maximum orthogonal angle between the strand pair and the helix." , "180.0" )
		+ XMLSchemaAttribute::attribute_w_default( "min_plane_angle"  , xsct_real , "Minimum plane angle between the strand pair and the helix." , "-180.0" )
		+ XMLSchemaAttribute::attribute_w_default( "max_plane_angle"  , xsct_real , "Maximum plane angle between the strand pair and the helix." , "180.0" )
		+ XMLSchemaAttribute::attribute_w_default( "StrandID1"  , xsct_positive_integer , "Strand number of the first strand to be considered, according to blueprint" , "1" )
		+ XMLSchemaAttribute::attribute_w_default( "StrandID2"  , xsct_positive_integer , "Strand number of the second strand to be considered, according to blueprint" , "2" )
		+ XMLSchemaAttribute::attribute_w_default( "HelixID"  , xsct_positive_integer , "Helix number of helix to be considered, according to blueprint" , "1" )
		+ XMLSchemaAttribute::attribute_w_default( "output_type", xs_string , "What metric should be output? distance (dist), orthogonal angle (ortho_angle), or the plane angle (plane_angle)?" , "dist" )
		+ XMLSchemaAttribute::attribute_w_default( "blueprint", xs_string , "path to blueprint file from which to parse strands" , "" );

	protocols::filters::xsd_type_definition_w_attributes(
		xsd,
		class_name(),
		"This filter will check the distanca, plane and orthogonal angle between a helix and a pair of strands. The plane, used to calculate both angles, is the one described by the vector parallel to the first strand and one from the center of the first to the center of the sencond (perpendicular). The plane angle is the angle between the first of named vectors and the projection of the helix on the name plane. The orthogonal angle is the one between the plane and the helix. The distance is measured between the center of mass of both elements.",
		attlist );
}

/////////////// Creator ///////////////

protocols::filters::FilterOP
StrandHelixGeometryFilterCreator::create_filter() const
{
	return protocols::filters::FilterOP( new StrandHelixGeometryFilter );
}

std::string
StrandHelixGeometryFilterCreator::keyname() const
{
	return StrandHelixGeometryFilter::class_name();
}

void StrandHelixGeometryFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	StrandHelixGeometryFilter::provide_xml_schema( xsd );
}

} //protocols
} //fldsgn
} //filters
