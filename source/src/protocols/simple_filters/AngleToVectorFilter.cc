// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/AngleToVectorFilter.cc
/// @brief
/// @author Ravit Netzer & Sarel Fleishman


//Unit Headers
#include <protocols/simple_filters/AngleToVectorFilter.hh>
#include <protocols/simple_filters/AngleToVectorFilterCreator.hh>
#include <utility/tag/Tag.hh>
//Project Headers
#include <basic/Tracer.hh>
#include <core/pose/Pose.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

namespace protocols {
namespace simple_filters {

using core::Size;
using core::Real;
using std::string;
using utility::vector1;

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_filters.AngleToVector" );

// XRW TEMP protocols::filters::FilterOP
// XRW TEMP AngleToVectorFilterCreator::create_filter() const { return protocols::filters::FilterOP( new AngleToVector ); }

// XRW TEMP std::string
// XRW TEMP AngleToVectorFilterCreator::keyname() const { return "AngleToVector"; }

//default ctor
AngleToVector::AngleToVector() :
	protocols::filters::Filter( "AngleToVector" ),
	min_angle_( 0.0 ),
	max_angle_( 90.0 ),
	refx_( 0.0 ), refy_( 0.0 ), refz_( 0.0 ),
	chain_( 2 ),
	atm1_( "" ), atm2_( "" )
{
}

AngleToVector::~AngleToVector() = default;

void
AngleToVector::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &, moves::Movers_map const &, core::pose::Pose const & )
{
	min_angle( tag->getOption< Real >( "min_angle", 0.0 ) );
	max_angle( tag->getOption< Real >( "max_angle", 90.0 ) );
	atm1( tag->getOption< string >( "atm1" ) );
	atm2( tag->getOption< string >( "atm2" ) );
	chain( tag->getOption< Size >( "chain", 2 ) );
	refx( tag->getOption< Real >( "refx" ) );
	refy( tag->getOption< Real >( "refy" ) );
	refz( tag->getOption< Real >( "refz" ) );

	TR<<"AngleToVector with options: chain "<<chain()<<" atm1: "<<atm1()<<" atm2: "<<atm2()<<" min_angle "<<min_angle()<<" max_angle "<<max_angle()<<std::endl;
}

core::Real
AngleToVector::compute( core::pose::Pose const & pose ) const {
	using namespace numeric;

	core::conformation::Residue const res = pose.conformation().residue( pose.conformation().chain_begin( chain() ) );
	xyzVector< Real > diff_vec = res.atom( atm1_ ).xyz() - res.atom( atm2_ ).xyz();
	diff_vec.normalize();
	TR<<"diff vec "<<diff_vec.x()<<' '<<diff_vec.y()<<' '<<diff_vec.z()<<std::endl;

	xyzVector< Real > ref_vec( refx(), refy(), refz() );
	ref_vec.normalize();

	Real const dot_prod = ref_vec.dot( diff_vec );
	TR<<"dot: "<<dot_prod<<std::endl;
	Real const angle = acos( dot_prod ) * 180.0 / 3.1415927;

	TR<<"angle: "<<angle<<std::endl;
	return angle;
}

bool
AngleToVector::apply( core::pose::Pose const & pose ) const{
	core::Real const angle( compute( pose ) );
	if ( angle >= min_angle() && angle <= max_angle() ) return true;
	return false;
}

void
AngleToVector::report( std::ostream &, core::pose::Pose const & ) const {
	// os<<"angle: "<<compute( pose )<<std::endl;
}

core::Real
AngleToVector::report_sm( core::pose::Pose const & pose ) const {
	return compute( pose );
}

std::string AngleToVector::name() const {
	return class_name();
}

std::string AngleToVector::class_name() {
	return "AngleToVector";
}

void AngleToVector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default( "min_angle", xsct_real, "Minimum angle to pass the filter", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default( "max_angle", xsct_real, "Maximum angle to pass the filter", "90.0" )
		+ XMLSchemaAttribute( "atm1", xs_string, "First atom in vector definition" )
		+ XMLSchemaAttribute( "atm2", xs_string, "Seocnd atom in vector definition" )
		+ XMLSchemaAttribute( "chain", xsct_char, "Chain whose first residue's two atoms are in question" )
		+ XMLSchemaAttribute( "refx", xsct_real, "x coordinate of vector" )
		+ XMLSchemaAttribute( "refy", xsct_real, "y coordinate of vector" )
		+ XMLSchemaAttribute( "refz", xsct_real, "z coordinate of vector" );

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "Filter on the value of the angle between a vector between two atoms in the first residue of a chain and a reference vector.", attlist );
}

std::string AngleToVectorFilterCreator::keyname() const {
	return AngleToVector::class_name();
}

protocols::filters::FilterOP
AngleToVectorFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new AngleToVector );
}

void AngleToVectorFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	AngleToVector::provide_xml_schema( xsd );
}


}
}
