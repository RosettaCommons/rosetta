// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/fldsgn/filters/StrandCurvatureByLevels.cc
/// @brief Newer version of filter used in Marcos & Basanta et al. 2017
/// @author Enrique Marcos (emarcos82@uw.edu)

#include <protocols/fldsgn/filters/StrandCurvatureByLevels.hh>
#include <protocols/fldsgn/filters/StrandCurvatureByLevelsCreator.hh>

//XSD includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

// Package Headers
#include <protocols/fldsgn/topology/SheetFoldTypeManager.hh>
#include <protocols/fldsgn/topology/StrandPairing.hh>
#include <protocols/fldsgn/topology/util.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <protocols/parser/BluePrint.hh>
#include <protocols/fldsgn/topology/SS_Info2.hh>

// Utility headers
#include <basic/Tracer.hh>

// Parser headers
#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

#include <numeric/conversions.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/trig.functions.hh>

static basic::Tracer TR( "protocols.fldsgn.filters.StrandCurvatureByLevels" );

namespace protocols {
namespace fldsgn {
namespace filters {

StrandCurvatureByLevels::StrandCurvatureByLevels():
	Filter( "StrandCurvatureByLevels" ),
	secstruct_( "" ),
	bend_level_(0),
	filter_min_bend_(0.0),
	filter_max_bend_(0.0),
	twist_level_(0),
	filter_min_twist_(0.0),
	filter_max_twist_(0.0),
	strand_id_(0),
	output_type_( "bend" ),
	output_value_( -990.0 ),
	concavity_reference_residue_("last"),
	concavity_direction_(1)
{}

StrandCurvatureByLevels::~StrandCurvatureByLevels()
{}
/*
// @brief copy constructor
StrandCurvatureByLevels::StrandCurvatureByLevels( StrandCurvatureByLevels const & rval ):
Super( rval ),
secstruct_( rval.secstruct_ ),
bend_level_( rval.bend_level_ ),
filter_min_bend_( rval.filter_min_bend_ ),
filter_max_bend_( rval.filter_max_bend_ ),
twist_level_( rval.twist_level_ ),
filter_min_twist_( rval.filter_min_twist_ ),
filter_max_twist_( rval.filter_max_twist_ ),
strand_id_( rval.strand_id_ ),
output_type_( rval.output_type_ ),
output_value_( rval.output_value_ ),
concavity_reference_residue_( rval.concavity_reference_residue_ ),
concavity_direction_( rval.concavity_direction_ )
{}
*/

// @brief set secondary structure
void StrandCurvatureByLevels::secstruct( String const & ss )
{
	secstruct_ = ss;
}

void
StrandCurvatureByLevels::bend_level( core::Size const r )
{
	bend_level_ = r;
}

/// @brief minimum angle for filtering
void
StrandCurvatureByLevels::filter_min_bend( core::Real const r )
{
	filter_min_bend_ = r;
}

/// @brief maximum angle for filtering
void
StrandCurvatureByLevels::filter_max_bend( core::Real const r )
{
	filter_max_bend_ = r;
}

void
StrandCurvatureByLevels::twist_level( core::Size const r )
{
	twist_level_ = r;
}

// @brief minimum twist for filtering
void
StrandCurvatureByLevels::filter_min_twist( core::Real const r )
{
	filter_min_twist_ = r;
}

/// @brief maximum twist for filtering
void
StrandCurvatureByLevels::filter_max_twist( core::Real const r )
{
	filter_max_twist_ = r;
}

/// @brief Strand id number
void
StrandCurvatureByLevels::strand_id( core::Size const r )
{
	strand_id_ = r;
}

/// @brief
void
StrandCurvatureByLevels::output_type( String const & s )
{
	output_type_ = s;
}

void
StrandCurvatureByLevels::concavity_reference_residue( String const & ss )
{
	concavity_reference_residue_ = ss;
}

void
StrandCurvatureByLevels::concavity_direction( bool const & r )
{
	concavity_direction_ = r;
}

// @brief returns true if the given pose passes the filter, false otherwise.
// In this case, the test is whether the give pose is the topology we want.
bool StrandCurvatureByLevels::apply( Pose const & pose ) const
{
	using core::scoring::dssp::Dssp;
	using protocols::fldsgn::topology::SS_Info2_OP;
	using protocols::fldsgn::topology::NO_STRANDS;

	// set secondary structure
	if ( secstruct_ == "" ) {
		Dssp dssp( pose );
		secstruct_ = dssp.get_dssp_secstruct();
	}
	//secstruct = secstruct_;
	// set SS_Info
	//SS_Info2_OP  ssinfo = new SS_Info2( pose, secstruct_ );  // old format now gives segfault.
	SS_Info2_OP  ssinfo( new SS_Info2( pose, secstruct_ ) ); // old format now gives segfault.
	TR << "ss of input pose=" << secstruct_ << std::endl;
	protocols::fldsgn::topology::Strands const & strands( ssinfo->strands() );

	core::Size begin_res( strands[ strand_id_ ]->begin() );
	core::Size end_res( strands[ strand_id_ ]->end() );
	core::Size st_len = end_res - begin_res + 1 ;
	core::Size center = begin_res + st_len / 2 ;

	core::Real mean_bend, mean_twist ;
	bool is_concave(false) ;

	// BENDING
	if ( bend_level_ != 0 ) {
		numeric::xyzVector<core::Real> v1;
		numeric::xyzVector<core::Real> v2;
		core::Size step ( 2*bend_level_) ;
		core::Real bend (0.0) ;
		core::Real bend_sum (0.0) ;
		core::Size n(0) ;
		core::Size p1(0), p2(0), p3(0) ;

		for ( core::Size k=begin_res; k<=end_res-step*2; k++ ) {
			p1 = k ;
			p2 = k + step*1 ;
			p3 = k + step*2 ;
			v1 = ( pose.residue( p2 ).xyz("CA") - pose.residue( p1 ).xyz("CA") ).normalized();
			v2 = ( pose.residue( p3 ).xyz("CA") - pose.residue( p2 ).xyz("CA") ).normalized();
			bend = numeric::conversions::degrees( angle_of(v1,v2) );
			bend_sum += bend ;
			n+=1 ;
		}
		mean_bend = bend_sum / n ;
		TR << " Average bend of strand " << strand_id_ << " at level " << bend_level_ << ": " << mean_bend << std::endl;
	} else {
		runtime_assert( bend_level_ != 0 );
		mean_bend = 0 ;
	}
	// TWIST
	if ( twist_level_ != 0 ) {
		core::Size step ( 2*twist_level_ );
		core::Size n(0) ;
		core::Real twist (0.0) ;
		core::Real twist_sum (0.0) ;
		core::Size p1(0), p2(0) ;
		for ( core::Size k=begin_res; k<=end_res-step; k++ ) {
			p1 = k ;
			p2 = k + step ;
			twist = numeric::dihedral_degrees( pose.residue(p1).xyz("CB"), pose.residue(p1).xyz("CA"), pose.residue(p2).xyz("CA"),pose.residue(p2).xyz("CB") ) ;
			twist_sum += twist ;
			n+=1 ;
		}
		mean_twist = twist_sum / n ;
		TR << " Average twist of strand " << strand_id_ << " at level " << twist_level_ << ": " << mean_twist << std::endl;
	} else {
		runtime_assert( twist_level_ != 0 );
		mean_twist = 0 ;
	}

	// Check concavity
	if ( concavity_reference_residue_ != "" ) {
		numeric::xyzVector<core::Real> v1;
		numeric::xyzVector<core::Real> v2;
		numeric::xyzVector<core::Real> w;
		numeric::xyzVector<core::Real> end_ba;
		numeric::xyzVector<core::Real> first_ba;
		numeric::xyzVector<core::Real> ref_ba;
		// vectors from edges to the center and summed give a vector pointing to the concave face
		v1 = pose.residue(begin_res).xyz("CA") - pose.residue(center).xyz("CA");
		v2 = pose.residue(end_res).xyz("CA") - pose.residue(center).xyz("CA") ;
		w = v1+v2 ;
		if ( concavity_reference_residue_ == "last" ) {
			if ( (center-end_res) % 2 == 0 ) { // same orientation in a regular strand
				//ref_ba = pose.residue(end_res).xyz("CB") - pose.residue(end_res).xyz("CA") ;
				ref_ba = pose.residue(center).xyz("CB") - pose.residue(center).xyz("CA") ;
			} else {
				ref_ba = pose.residue(center-1).xyz("CB") - pose.residue(center-1).xyz("CA") ;
			}
		} else if ( concavity_reference_residue_ == "first" ) {
			if ( (center-begin_res) % 2 == 0 ) { // same orientation in a regular strand
				//ref_ba = pose.residue(begin_res).xyz("CB") - pose.residue(begin_res).xyz("CA") ;
				ref_ba = pose.residue(center).xyz("CB") - pose.residue(center).xyz("CA") ;
			} else {
				ref_ba = pose.residue(center-1).xyz("CB") - pose.residue(center-1).xyz("CA") ;
			}
		} else {
			runtime_assert( concavity_reference_residue_ == "first" || concavity_reference_residue_ == "last"  );
			ref_ba = pose.residue(center-1).xyz("CB") - pose.residue(center-1).xyz("CA") ;
		}

		if ( concavity_direction_ == true  ) {
			if ( ref_ba.dot(w) > 0 ) {
				is_concave=true ;
			} else {
				is_concave=false;
			}
		} else {
			if ( ref_ba.dot(w) < 0 ) {
				is_concave=true ;
			} else {
				is_concave=false;
			}
		}
	}
	bool filter( true );
	if ( bend_level_ != 0 ) {
		if  ( mean_bend < filter_min_bend_ ||  mean_bend > filter_max_bend_ ) {
			filter = false;
			TR << " Average bend of strand " << strand_id_ << " at level " << bend_level_ << " has a value out of bounds "   << std::endl ;

		}
	}
	if ( twist_level_ != 0 ) {
		if ( mean_twist < filter_min_twist_ ||  mean_twist > filter_max_twist_ ) {
			filter = false;
			TR << " Average twist of strand " << strand_id_ << " at level " << twist_level_ << " has a value out of bounds "  << std::endl ;
		}
	}
	if ( concavity_reference_residue_ != "" ) {
		if ( is_concave == false ) {
			filter=false ;
			TR << " Strand " << strand_id_ << " is convex" << std::endl ;
		} else {
			TR << " Strand " << strand_id_ << " is concave" << std::endl ;
		}
	}

	if ( filter ) {
		TR << " Filter success ! " << std::endl;
	} else {
		TR << " Filter failed ! " << std::endl;
	}

	if ( output_type_ == "bend" ) {
		output_value_ = mean_bend;
	}
	if ( output_type_ == "twist" ) {
		output_value_ = mean_twist;
	}

	return filter;

} // apply_filter

void
StrandCurvatureByLevels::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & ,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
	// Blueprint is for giving secondary structure information, otherwise dssp will run for ss definition
	// SSPAIR line is read for the topology of strand pairings
	bend_level_ = tag->getOption<core::Size>( "bend_level", 1 );
	filter_min_bend_ = tag->getOption<Real>( "min_bend", 0.0 );
	filter_max_bend_ = tag->getOption<core::Real>( "max_bend", 180.0 );
	twist_level_ = tag->getOption<core::Size>( "twist_level", 1 );
	filter_min_twist_ = tag->getOption<core::Real>( "min_twist", 0.0 );
	filter_max_twist_ = tag->getOption<core::Real>( "max_twist", 90.0 );
	strand_id_ = tag->getOption<core::Size>( "StrandID", 1 );
	output_type_ = tag->getOption<String>( "output_type", "bend" );
	concavity_reference_residue_ = tag->getOption<String>( "concavity_reference_residue", "" );
	concavity_direction_ = tag->getOption<bool>( "concavity_direction", true );

	String const blueprint = tag->getOption<String>( "blueprint", "" );
	if ( blueprint != "" ) {
		protocols::parser::BluePrint blue( blueprint );
		secstruct_ = blue.secstruct();
	}
}

protocols::filters::FilterOP
StrandCurvatureByLevels::clone() const
{
	return protocols::filters::FilterOP( new StrandCurvatureByLevels( *this ) );
}


protocols::filters::FilterOP
StrandCurvatureByLevels::fresh_instance() const
{
	return protocols::filters::FilterOP( new StrandCurvatureByLevels );
}

core::Real
StrandCurvatureByLevels::report_sm( core::pose::Pose const & pose ) const
{
	return compute( pose );
}

// @brief returns computed value
StrandCurvatureByLevels::Real
StrandCurvatureByLevels::compute( Pose const & pose ) const
{
	apply( pose );
	return output_value_;
}

void
StrandCurvatureByLevels::report( std::ostream &, core::pose::Pose const & ) const
{

}

std::string StrandCurvatureByLevels::name() const {
	return class_name();
}

std::string StrandCurvatureByLevels::class_name() {
	return "StrandCurvatureByLevels";
}

void StrandCurvatureByLevels::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	attlist
		+ XMLSchemaAttribute::attribute_w_default( "bend_level"  , xsct_positive_integer , "Number of CA pairs left between the vertex and the end of the arm used to calculate the bend angle." , "1" )
		+ XMLSchemaAttribute::attribute_w_default( "min_bend" , xsct_real , "Minimum allowed for bend angle, the complementary to CA1-CA2-CA3, so that the higher this angle, the more curved the strand is. Value is in degrees." , "0.0" )
		+ XMLSchemaAttribute::attribute_w_default( "max_bend" , xsct_real , "Maximum allowed for bend angle, the complementary to CA1-CA2-CA3, so that the higher this angle, the more curved the strand is. Value is in degrees." , "180" )
		+ XMLSchemaAttribute::attribute_w_default( "twist_level"  , xsct_positive_integer , "Number of CA-CB vecgor pairs left between the vertex and the end of the arm used to calculate the twist angle (dihedral between two CA-CB vetors)." , "1" )
		+ XMLSchemaAttribute::attribute_w_default( "min_twist" , xsct_real , "Minimum allowed for strand twist angle, the complementary to a CA-CB - CA-CB dihedral, so that the higher this angle, the more twisted the strand is. Value is in degrees." , "0.0" )
		+ XMLSchemaAttribute::attribute_w_default( "max_twist" , xsct_real , "Maximum allowed for strand twist angle, the complementary to a CA-CB - CA-CB dihedral, so that the higher this angle, the more twisted the strand is. Value is in degrees." , "90" )
		+ XMLSchemaAttribute::attribute_w_default( "StrandID"  , xsct_positive_integer , "Strand number over which to calculate bend, twist and concavity, according to blueprint numbers" , "1" )
		+ XMLSchemaAttribute::attribute_w_default( "output_type", xs_string , "What magnitude to inform: \"bend\" or \"twist\" angle" , "bend" )
		+ XMLSchemaAttribute::attribute_w_default( "concavity_reference_residue", xs_string , "Use first or last residue as a reference for concavity: \"first\" or \"last\" are the only options. See \"concavity_direction\" for more info on how to use this option. Defaults should be fine for regular (not bulged) strands." , "first" )
		+ XMLSchemaAttribute::attribute_w_default( "concavity_direction", xs_boolean , "Does the reference residue CA-CB vector point in the same face as the residue in the center of the strand? If yes, then set this option to true, otherwise, false. Defaults should be fine for regular (not bulged) strands." , "true" )
		+ XMLSchemaAttribute::attribute_w_default( "blueprint", xs_string , "path to blueprint file from which to parse strands" , "" );

	protocols::filters::xsd_type_definition_w_attributes(
		xsd,
		class_name(),
		"Calculate the direction of the concave face of a strand, and its curvarure based on the angle formed by alternating CA atoms. The \"levels\" are the number of CA pais left between the vertex and the end of the arm.",
		attlist );
}

/////////////// Creator ///////////////

protocols::filters::FilterOP
StrandCurvatureByLevelsCreator::create_filter() const
{
	return protocols::filters::FilterOP( new StrandCurvatureByLevels );
}

std::string
StrandCurvatureByLevelsCreator::keyname() const
{
	return StrandCurvatureByLevels::class_name();
}

void StrandCurvatureByLevelsCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	StrandCurvatureByLevels::provide_xml_schema( xsd );
}

} //protocols
} //fldsgn
} //filters
