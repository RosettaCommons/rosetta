// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/fldsgn/filters/HelixBendFilter.cc
/// @brief Filter used in 'Principles for designing proteins with cavities formed by curved b-sheets' to control helix geometry.
/// @author Enrique Marcos (emarcos82@gmail.com)
/// Update and integration test by Benjamin Basanta (basantab@uw.edu)

#include <protocols/fldsgn/filters/HelixBendFilter.hh>
#include <protocols/fldsgn/filters/HelixBendFilterCreator.hh>

//XSD includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/selection.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <protocols/parser/BluePrint.hh>
#include <protocols/fldsgn/topology/SS_Info2.hh>

#include <numeric/conversions.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/trig.functions.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/vector1.functions.hh>

static basic::Tracer TR( "protocols.fldsgn.filters.HelixBendFilter" );

namespace protocols {
namespace fldsgn {
namespace filters {

HelixBendFilter::HelixBendFilter():
	protocols::filters::Filter( "HelixBendFilter" ),
	secstruct_( "" ),
	threshold_(0.0),
	helix_id_(0),
	filter_status_(false)
{}

HelixBendFilter::~HelixBendFilter() = default;

/// @brief set secondary structure
void HelixBendFilter::secstruct( std::string const & ss )
{
	secstruct_ = ss;
}

/// @brief minimum angle for filtering
void
HelixBendFilter::threshold( core::Real const r )
{
	threshold_ = r;
}

/// @brief Strand id number
void
HelixBendFilter::helix_id( core::Size const r )
{
	helix_id_ = r;
}

/// @brief set secondary structure from blueprint
/// @details This triggers a read from disk; use with caution.
void
HelixBendFilter::set_secstruct_from_bp(std::string const & bp_file_path)
{
	protocols::parser::BluePrint blue( bp_file_path );
	HelixBendFilter::secstruct(blue.secstruct());
}

void
HelixBendFilter::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & ,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
	// Blueprint is for giving secondary structure information, otherwise dssp will run for ss definition
	// SSPAIR line is read for the topology of strand pairings
	threshold_ = tag->getOption<core::Real>( "threshold", 155.0 );
	helix_id_ = tag->getOption<core::Size>( "HelixID", 1 );
	std::string const blueprint = tag->getOption<std::string>( "blueprint", "" );
	std::string const secstruct_str = tag->getOption<std::string>( "secstruct", "" );
	if ( (not blueprint.empty()) and secstruct_str.empty() ) {
		set_secstruct_from_bp(blueprint);
	}
	// string has priority if blueprint and secondary structure string are defined.
	if ( not secstruct_str.empty() ) {
		secstruct(secstruct_str);
	}
	// If secstruct remains empty, DSSP is used (see compute).

}

protocols::filters::FilterOP
HelixBendFilter::clone() const
{
	return protocols::filters::FilterOP( new HelixBendFilter( *this ) );
}


protocols::filters::FilterOP
HelixBendFilter::fresh_instance() const
{
	return protocols::filters::FilterOP( new HelixBendFilter );
}

bool
HelixBendFilter::apply( core::pose::Pose const & pose ) const
{

	if ( pose.size() != secstruct_.length() ) {
		TR.Error << "Secstruct string length =! pose length" << std::endl;
		utility_exit_with_message("Unable to parse secstruct, string length =! pose length");
	}
	compute(pose);
	return filter_status_;
}

void
HelixBendFilter::compute( core::pose::Pose const & pose ) const
{
	using namespace protocols::fldsgn::topology;
	using core::scoring::dssp::Dssp;
	using protocols::fldsgn::topology::Helices;
	using protocols::fldsgn::topology::Strands;
	using protocols::fldsgn::topology::SS_Info2;
	using protocols::fldsgn::topology::SS_Info2_OP;

	// set secondary structure
	if ( secstruct_ == "" ) {
		Dssp dssp( pose );
		secstruct_ = dssp.get_dssp_secstruct();
	}
	// set SS_Info
	//SS_Info2_OP  ssinfo = new SS_Info2( pose, secstruct_ ); // old format
	SS_Info2_OP  ssinfo( new SS_Info2( pose, secstruct_ ) );

	//protocols::fldsgn::topology::Strands const & strands( ssinfo->strands() );
	protocols::fldsgn::topology::Helices const & helices( ssinfo->helices() );

	// Helix positions
	core::Size begin_res( helices[ helix_id_ ]->begin() );
	core::Size end_res( helices[ helix_id_ ]->end() );
	utility::vector1<core::Real> curves;
	core::Real min_curve;
	bool filter( true );
	numeric::xyzVector<core::Real> v1;
	numeric::xyzVector<core::Real> v2;
	core::Size h_len = end_res - begin_res + 1;
	if ( h_len < 8 ) {
		TR << "This helix is too short for measuring curvature, passing by default with max curve = 180.0" << std::endl;
		min_curve = 180.0;
	} else {
		for ( core::Size i=begin_res+4; i<=end_res-4; i++ ) {
			v1 = pose.residue(i-4).xyz("CA") - pose.residue(i).xyz("CA") ;
			v2 = pose.residue(i+4).xyz("CA") - pose.residue(i).xyz("CA") ;
			core::Real d1 = pose.residue( i-4 ).atom("O").xyz().distance( pose.residue(i).atom("N").xyz() ) ;
			core::Real d2 = pose.residue( i ).atom("O").xyz().distance( pose.residue(i+4).atom("N").xyz() ) ;
			core::Real angle = numeric::conversions::degrees( angle_of(v1,v2) ) ;
			if ( d1 > 3.6  ) {
				filter=false ;
				TR << "Hbond too large between residues "<< " " << i-4 << " " << i << std::endl;
			}

			if ( d2 > 3.6 ) {
				filter=false ;
				TR << "Hbond too large between residues "<< " " << i << " " << i+4 << std::endl;
			}


			if ( angle < threshold_ ) {
				TR << " Helix has local bend above threshold  = " << angle << std::endl;
				filter=false;
			}
			curves.push_back(angle);
			/// TR << " This angle  = " << angle << std::endl;
		}
		/// core::Size turns = ( end_res - 4 ) - ( begin_res + 4 ) + 1;
		min_curve = utility::min( curves );
	}

	filter_status_ = filter ;
	filter_val_ = min_curve ;
}

core::Real
HelixBendFilter::report_sm( core::pose::Pose const & pose) const
{
	compute(pose);
	if ( filter_status_ ) {
		TR << "Filter passed: all helix turns CA angle are above threshold "<< threshold_ << "." << std::endl;
		TR << "All H-bonding heavy atom pairs have a distance below 3.6A." << std::endl;
		TR << "Max bending (min angle) of helix observed: "<< filter_val_ << std::endl;
	} else {
		TR << "Filter did not pass." << std::endl;
	}

	return filter_val_;
}

void
HelixBendFilter::report( std::ostream &, core::pose::Pose const & ) const
{
	///compute(pose);
	if ( filter_status_ ) {
		TR << "Filter passed: all helix turns CA angle are above threshold "<< threshold_ << "." << std::endl;
		TR << "All H-bonding heavy atom pairs have a distance below 3.6A." << std::endl;
		TR << "Max bending (min angle) of helix observed: "<< filter_val_ << std::endl;
	} else {
		TR << "Filter did not pass." << std::endl;
	}
}

std::string HelixBendFilter::name() const {
	return class_name();
}

std::string HelixBendFilter::class_name() {
	return "HelixBendFilter";
}

void HelixBendFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default( "threshold" , xsct_real , "Minimum bending angle allowed." , "155.0" )
		+ XMLSchemaAttribute::attribute_w_default( "HelixID"  , xsct_positive_integer , "Helix number over which to calculate bend and Hbond pair distance according to blueprint numbers" , "1" )
		+ XMLSchemaAttribute::attribute_w_default( "blueprint", xs_string , "Path to blueprint file from which to parse strands" , "" )
		+ XMLSchemaAttribute::attribute_w_default( "secstruct", xs_string , "Secondary structure string to use for selecting helices" , "" );

	protocols::filters::xsd_type_definition_w_attributes(
		xsd,
		class_name(),
		"This filter will check each helix turn in target HelixID helix for further-than-optimal Hbond pairs and curvature below (smaller angles) than cutoff, it reports the minimum angle value. If no secondary structure string or blueprint are provided, then the pose secstruct is parsed using DSSP.",
		attlist );
}

/////////////// Creator ///////////////

protocols::filters::FilterOP
HelixBendFilterCreator::create_filter() const
{
	return protocols::filters::FilterOP( new HelixBendFilter );
}

std::string
HelixBendFilterCreator::keyname() const
{
	return HelixBendFilter::class_name();
}

void HelixBendFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	HelixBendFilter::provide_xml_schema( xsd );
}

} //protocols
} //fldsgn
} //filters
