// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/matdes/GetRBDOFValues.cc
/// @brief  Returns the translation and angle for movable jump
/// @author Jacob Bale (balej@u.washington.edu)

// Unit Headers
#include <protocols/matdes/GetRBDOFValues.hh>
#include <protocols/matdes/GetRBDOFValuesCreator.hh>
#include <protocols/matdes/SymDofMoverSampler.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/conformation/Conformation.hh>

// Utility headers
#include <utility>
#include <utility/vector1.fwd.hh>
#include <basic/Tracer.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <iostream>
#include <utility/exit.hh>

// Parser headers
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>

#include <utility/vector1.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/conversions.hh>
#include <numeric/xyzVector.string.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>


//// C++ headers
static THREAD_LOCAL basic::Tracer TR( "protocols.matdes.GetRBDOFValues" );

namespace protocols {
namespace matdes {

using core::pose::Pose;

// @brief default constructor
GetRBDOFValues::GetRBDOFValues():
	jump_id_( 1 ),
	sym_dof_name_( "" ),
	verbose_( 0 ),
	radial_disp_( 0 ),
	angle_( 0 ),
	get_init_value_( 0 ),
	axis_( 'x' ),
	init_disp_( 0 ),
	init_angle_( 0 )
{}


// @brief constructor with arguments
GetRBDOFValues::GetRBDOFValues( int jump, std::string dof_name, bool verb, char ax, bool disp, bool ang, core::Real init_d, core::Real init_a, bool get_init ):
	jump_id_( jump),
	sym_dof_name_(std::move( dof_name )),
	verbose_( verb ),
	radial_disp_( disp ),
	angle_( ang ),
	get_init_value_( get_init ),
	axis_( ax ),
	init_disp_( init_d ),
	init_angle_( init_a )
{}

// @brief copy constructor
GetRBDOFValues::GetRBDOFValues( GetRBDOFValues const & )= default;

// @brief destructor
GetRBDOFValues::~GetRBDOFValues() = default;

protocols::filters::FilterOP
GetRBDOFValues::fresh_instance() const{
	return protocols::filters::FilterOP( new GetRBDOFValues() );
}

protocols::filters::FilterOP
GetRBDOFValues::clone() const{
	return protocols::filters::FilterOP( new GetRBDOFValues( *this ) );
}

// @brief getters
core::Size GetRBDOFValues::jump_id() const { return jump_id_; }
std::string GetRBDOFValues::sym_dof_name() const { return sym_dof_name_; }
bool GetRBDOFValues::verbose() const { return verbose_; }
char GetRBDOFValues::axis() const { return axis_; }
bool GetRBDOFValues::radial_disp() const { return radial_disp_; }
bool GetRBDOFValues::angle() const { return angle_; }
core::Real GetRBDOFValues::init_disp() const { return init_disp_; }
core::Real GetRBDOFValues::init_angle() const { return init_angle_; }
bool GetRBDOFValues::get_init_value() const { return get_init_value_; }

// @brief setters
void GetRBDOFValues::jump_id( core::Size const jump ) { jump_id_ = jump; }
void GetRBDOFValues::sym_dof_name( std::string const & dof_name ) { sym_dof_name_ = dof_name; }
void GetRBDOFValues::verbose( bool const verb ) { verbose_ = verb; }
void GetRBDOFValues::axis( char const ax ) { axis_ = ax; }
void GetRBDOFValues::radial_disp( bool const disp ) { radial_disp_ = disp; }
void GetRBDOFValues::angle( bool const ang ) { angle_ = ang; }
void GetRBDOFValues::init_disp( core::Real const init_d ) { init_disp_ = init_d; }
void GetRBDOFValues::init_angle( core::Real const init_a ) { init_angle_ = init_a; }
void GetRBDOFValues::get_init_value( bool const get_init ) { get_init_value_ = get_init; }

/// @brief
core::Real GetRBDOFValues::compute(
	Pose const & pose, bool const & verb,
	std::string const & dof_name,
	int const & jump,
	char const & ax,
	bool const & disp,
	bool const & ang,
	core::Real const & init_d,
	core::Real const & init_a,
	bool const & get_init ) const
{
	typedef numeric::xyzVector<Real> Vec;
	typedef numeric::xyzMatrix<Real> Mat;

	// Get the jump_id corresponding to the user-specified jump or sym_dof_name.
	int sym_aware_jump_id;
	if ( dof_name != "" ) {
		sym_aware_jump_id = core::pose::symmetry::sym_dof_jump_num( pose, dof_name );
	} else {
		sym_aware_jump_id = core::pose::symmetry::get_sym_aware_jump_num(pose, jump);
	}

	std::string design_id = protocols::jd2::JobDistributor::get_instance()->current_output_name();
	std::ostringstream val_string;
	Real value;
	std::string fn;

	int index = -1;
	utility::vector1<std::string> sym_dof_names;
	if ( get_init ) {
		if ( dof_name == "" ) utility_exit_with_message("A sym_dof_name must be specified in order to access the initial values from the SymDofSampler.");
		/// WARNING WARNING WARNING THREAD UNSAFE!
		sym_dof_names = SymDofMoverSampler::get_instance()->get_sym_dof_names();
		for ( Size i = 1; i <= sym_dof_names.size(); i++ ) {
			if ( dof_name == sym_dof_names[i] ) index = i;
		}
	}

	// Get the radial displacement or angle for the user specified jump.
	if ( disp && !ang ) {
		Real starting_disp = 0;
		/// WARNING WARNING WARNING THREAD UNSAFE!
		if ( get_init ) starting_disp = SymDofMoverSampler::get_instance()->get_radial_disps()[index];
		TR.Debug << "starting disp " << sym_dof_names[index] << ": " << init_d + starting_disp << std::endl;
		if ( ax == 'x' ) value = -pose.jump(sym_aware_jump_id).get_translation().x() + init_d + starting_disp; // note: negative
		else if ( ax == 'y' ) value = -pose.jump(sym_aware_jump_id).get_translation().y() + init_d + starting_disp; // note: negative
		else if ( ax == 'z' ) value = -pose.jump(sym_aware_jump_id).get_translation().z() + init_d + starting_disp; // note: negative
		else utility_exit_with_message("The axis specified for GetRBDOFValues does not match with x, y, or z.");
		val_string << std::fixed << std::setprecision(1) << value;
		fn = "SymDofName: " + dof_name + " radial_disp along axis " + ax + ": " + val_string.str();
	} else if ( ang && !disp ) {
		Real theta;
		Mat const & R = pose.jump(sym_aware_jump_id).get_rotation();
		Vec axis = numeric::rotation_axis(R, theta );
		TR.Debug << "axis: " << numeric::truncate_and_serialize_xyz_vector(axis, 2) << std::endl; // Check the vector about which theta was calculated.
		TR.Debug << theta << std::endl;
		int sign;
		if ( ax == 'x' ) sign = (axis.x()>0?-1:1); // Determine if the axis vector was chosen to be parallel or antiparallel to the jump by the rotation_axis function.
		else if ( ax == 'y' ) sign = (axis.y()>0?-1:1); // Determine if the axis vector was chosen to be parallel or antiparallel to the jump by the rotation_axis function.
		else if ( ax == 'z' ) sign = (axis.z()>0?-1:1); // Determine if the axis vector was chosen to be parallel or antiparallel to the jump by the rotation_axis function.
		else utility_exit_with_message("The axis specified for GetRBDOFValues does not match with x, y, or z.");
		TR.Debug << "sign: " << sym_dof_names[index] << ": " << sign << std::endl;
		Real starting_angle = 0;
		/// WARNING WARNING WARNING THREAD UNSAFE!
		if ( get_init ) starting_angle = SymDofMoverSampler::get_instance()->get_angles()[index];
		TR.Debug << "starting angle " << sym_dof_names[index] << ": " << starting_angle + init_a << std::endl;
		value = numeric::conversions::degrees(sign*theta) + init_a + starting_angle;
		val_string << std::fixed << std::setprecision(1) << value;
		fn = "SymDofName: " + dof_name + " angle of rotation: " + val_string.str();
	} else utility_exit_with_message("GetRBDOFValues works for either radial_displacement along a specified axis OR an angle about a specified jump.");

	if ( verb ) {
		TR << "Design: " << design_id << fn << std::endl;
	}

	return( value );

} // compute

// @brief Dummy Filter apply function
bool GetRBDOFValues::apply( Pose const & pose ) const
{
	//core::Real value( // Unused variable causes warning.
	compute(pose, verbose(), sym_dof_name(), jump_id(), axis(), radial_disp(), angle(), init_disp(), init_angle(), get_init_value() );
	// );
	return( true );
}

/// @brief parse xml
void
GetRBDOFValues::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & /*data*/,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
	jump_id( tag->getOption< int >( "jump", 1 ) );
	sym_dof_name( tag->getOption< std::string >( "sym_dof_name" , "" ) );
	verbose( tag->getOption< bool >("verbose", 0) );
	axis( tag->getOption< char >("axis", 'x') );
	radial_disp( tag->getOption< bool >("get_disp", 0 ) );
	angle( tag->getOption< bool >("get_angle", 0 ) );
	init_angle( tag->getOption< core::Real >("init_angle", 0 ) );
	init_disp( tag->getOption< core::Real >("init_disp", 0 ) );
	get_init_value( tag->getOption< bool >("get_init_value", 0 ) );
}

core::Real
GetRBDOFValues::report_sm( Pose const & pose ) const
{
	return( compute( pose, false, sym_dof_name(), jump_id(), axis(), radial_disp(), angle(), init_disp(), init_angle(), get_init_value() ) );
}

void
GetRBDOFValues::report( std::ostream & out, Pose const & pose ) const
{
	out << "GetRBDOFValues returns " << compute( pose, false, sym_dof_name(), jump_id(), axis(), radial_disp(), angle(), init_disp(), init_angle(), get_init_value() ) << std::endl;
}

// XRW TEMP protocols::filters::FilterOP
// XRW TEMP GetRBDOFValuesCreator::create_filter() const { return protocols::filters::FilterOP( new GetRBDOFValues ); }

// XRW TEMP std::string
// XRW TEMP GetRBDOFValuesCreator::keyname() const { return "GetRBDOFValues"; }

std::string GetRBDOFValues::name() const {
	return class_name();
}

std::string GetRBDOFValues::class_name() {
	return "GetRBDOFValues";
}

void GetRBDOFValues::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	XMLSchemaRestriction xyz_char;
	xyz_char.name( "type_for_xyz_char" );
	xyz_char.base_type( xs_string );
	xyz_char.add_restriction( xsr_pattern, "x|y|z" );
	xsd.add_top_level_element( xyz_char );

	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default( "jump" , xs_integer , "Jump number of movable jump for which to calculate the translation or rotation" , "1" )
		+ XMLSchemaAttribute::required_attribute( "sym_dof_name" , xs_string, "Sym_dof_name for movable jump for which to calculate the translation or rotation." )
		+ XMLSchemaAttribute::attribute_w_default( "verbose" , xsct_rosetta_bool , "Output jump and corresponding displacement or angle to tracer." , "0" )
		+ XMLSchemaAttribute::attribute_w_default( "axis" , "type_for_xyz_char" , "Axis in local coordinate frame about which to calculate the translation or rotation (not currently set up to handle off axis values)." , "x" )
		+ XMLSchemaAttribute::attribute_w_default( "get_disp" , xsct_rosetta_bool , "If set to true (and get_disp is false), then will calculate the displacement across the specified jump." , "0" )
		+ XMLSchemaAttribute::attribute_w_default( "get_angle" , xsct_rosetta_bool , "If set to true (and get_angle is false), then will calculate the angle of rotation about the specified jump." , "0" )
		+ XMLSchemaAttribute::attribute_w_default( "init_angle" , xsct_real , "Initial angle value to add to each calculated value." , "0" )
		+ XMLSchemaAttribute::attribute_w_default( "init_disp" , xsct_real , "Initial displacement value to add to each calculated value." , "0" )
		+ XMLSchemaAttribute::attribute_w_default( "get_init_value" , xsct_rosetta_bool , "Get the initial displacement or angle for the specified jump from the SymDofMoverSampler" , "0" ) ;

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "Calculates either the current translation or rotation across a user specified jump (referenced by jump_id or sym_dof_name).", attlist );
}

std::string GetRBDOFValuesCreator::keyname() const {
	return GetRBDOFValues::class_name();
}

protocols::filters::FilterOP
GetRBDOFValuesCreator::create_filter() const {
	return protocols::filters::FilterOP( new GetRBDOFValues );
}

void GetRBDOFValuesCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	GetRBDOFValues::provide_xml_schema( xsd );
}



} // matdes
} // protocols
