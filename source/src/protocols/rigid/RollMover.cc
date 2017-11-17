// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/moves/RollMover.cc
/// @brief RollMover methods implemented
/// @author

// Unit Headers
#include <protocols/rigid/RollMover.hh>
#include <protocols/rigid/RollMoverCreator.hh>
// Package Headers
#include <protocols/rigid/RB_geometry.hh>

// Project Headers
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/id/AtomID.hh>
#include <numeric/xyz.functions.hh>
#include <protocols/rosetta_scripts/util.hh>
// Random number generator
#include <numeric/random/random.hh>
// Utility Headers
#include <basic/Tracer.hh>
#include <core/types.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/tag/Tag.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


// C++ Headers

using basic::Error;
using basic::Warning;
static basic::Tracer TR( "protocols.rigid.RollMover" );

using namespace core;

namespace protocols {
namespace rigid {

/// @details
void RollMover::apply( core::pose::Pose & pose ){

	utility::vector1< utility::vector1< numeric::xyzVector< core::Real> > >  coords; // some of these will change for sure

	// now we have a vector1 of vector1s with a numeric::xyzVector
	// access will look like coords[residue][atom] to get xyz

	core::Size const nres( pose.size() );
	coords.resize( nres );

	for ( Size i=1; i<= nres; ++i ) {
		core::conformation::Residue const & rsd( pose.residue(i) );
		core::Size const number_atoms_this_residue( rsd.natoms() );
		if ( number_atoms_this_residue ) {
			coords[i].resize( number_atoms_this_residue );
			for ( Size j=1; j <= number_atoms_this_residue; ++j ) {
				coords[i][j] = rsd.atom( j ).xyz();
			}
		}
	}

	if ( random_roll_ ) {
		translate_ = protocols::geometry::center_of_mass(pose, start_res_, stop_res_);
		translate_.x() += numeric::random::rg().gaussian() * random_roll_trans_;
		translate_.y() += numeric::random::rg().gaussian() * random_roll_trans_;
		translate_.z() += numeric::random::rg().gaussian() * random_roll_trans_;
		angle_ = numeric::random::rg().gaussian() * random_roll_angle_;
		axis_ = numeric::xyzVector< core::Real >( numeric::random::rg().gaussian(),numeric::random::rg().gaussian(),numeric::random::rg().gaussian() ).normalized();
	} else {
		angle_ = min_angle_ + ( max_angle_ - min_angle_ ) * numeric::random::rg().uniform();
	}

	numeric::xyzMatrix< core::Real > rotation_matrix( numeric::rotation_matrix_degrees(axis_, angle_ ) );
	//move to origin
	for ( core::Size i =start_res_; i <= stop_res_; ++i ) {
		for ( core::Size j = 1; j <= coords[i].size(); ++j ) {

			// this may look strange but in a global coordinate system
			// rotation about an axis is easily done by movement to the origin
			// rotation and then movement back

			coords[i][j] = coords[i][j] - translate_; // translate to origin
			coords[i][j] = rotation_matrix * coords[i][j]; // rotate atom
			coords[i][j] = coords[i][j] + translate_; // reverse translate

		}
	}

	// now update pose with new coordinates
	for ( core::Size i =start_res_; i <= stop_res_; ++i ) {
		for ( core::Size j = 1; j <= coords[i].size(); ++j ) {
			id::AtomID id( j, i );
			pose.set_xyz( id, coords[i][j]);
		}
	}



}//apply


void
RollMover::set_min_max_angles( core::Real min_angle, core::Real max_angle ) {
	min_angle_ = min_angle;
	max_angle_ = max_angle;
}

// XRW TEMP std::string
// XRW TEMP RollMover::get_name() const {
// XRW TEMP  return "RollMover";
// XRW TEMP }

void
RollMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & /*datamap*/,
	Filters_map const & /*filters*/,
	moves::Movers_map const & /*movers*/,
	Pose const & pose )
{

	start_res_ = ( tag->hasOption("start_res") ) ?  tag->getOption<core::Size>("start_res") : 1;
	stop_res_ = ( tag->hasOption("stop_res") ) ?  tag->getOption<core::Size>("stop_res") : pose.size();

	if ( tag->hasOption("chain") ) {
		if ( tag->hasOption("start_res") || tag->hasOption("stop_res") ) utility_exit_with_message("cannot specify start/stop res AND chain!");
		core::Size const chain = tag->getOption<core::Size>("chain");
		runtime_assert_msg(chain > 0 && chain <= pose.conformation().num_chains(),"RollMover bad chain");
		start_res_ = pose.conformation().chain_begin(chain);
		stop_res_ = pose.conformation().chain_end(chain);
	}

	random_roll_       = tag->getOption<bool>("random_roll",false);
	random_roll_angle_ = tag->getOption<core::Real>("random_roll_angle_mag",3.0);
	random_roll_trans_ = tag->getOption<core::Real>("random_roll_trans_mag",1.0);

	if ( ! random_roll_ ) {

		/*parse min_angle*/
		if ( tag->hasOption("min_angle") ) {
			min_angle_ = tag->getOption<core::Real>("min_angle");
		} else {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "RollMover requires min_angle option");
		}
		/*parse max_angle*/
		if ( tag->hasOption("max_angle") ) {
			max_angle_ = tag->getOption<core::Real>("max_angle");
		} else {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "RollMover requires max_angle option");
		}

		bool axis_option_parsed = false;
		bool translate_option_parsed = false;

		if ( tag->hasOption("axis") ) {
			switch(tag->getOption<char>("axis")) {
			case 'x' :
				axis_ = numeric::xyzVector< core::Real >( 1.0, 0.0, 0.0 );
				break;
			case 'y' :
				axis_ = numeric::xyzVector< core::Real >( 0.0, 1.0, 0.0 );
				break;
			case 'z' :
				axis_ = numeric::xyzVector< core::Real >( 0.0, 0.0, 1.0 );
				break;
			}
			axis_option_parsed = true;
		}

		for ( utility::tag::TagCOP child_tag : tag->getTags() ) {
			std::string name= child_tag->getName();

			if ( name == "axis" ) {
				/*parse axis x,y,z*/
				axis_ = protocols::rosetta_scripts::parse_xyz_vector(child_tag);
				axis_option_parsed = true;
			} else if ( name == "translate" ) {
				/*parse translate x,y,z*/
				translate_ = protocols::rosetta_scripts::parse_xyz_vector(child_tag);
				translate_option_parsed = true;
			}

		}

		if ( !axis_option_parsed ) {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "RollMover requires axis option");
		}
		if ( !translate_option_parsed ) {
			//throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "RollMover requires translate option");
			TR << "No translation given, using the pose's center of mass" << std::endl;
			translate_ = protocols::geometry::center_of_mass(pose, start_res_, stop_res_);
		}

	}
}

// XRW TEMP std::string
// XRW TEMP RollMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return RollMover::mover_name();
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP RollMover::mover_name()
// XRW TEMP {
// XRW TEMP  return "RollMover";
// XRW TEMP }


// XRW TEMP protocols::moves::MoverOP
// XRW TEMP RollMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new RollMover );
// XRW TEMP }

/// @brief required in the context of the parser/scripting scheme
moves::MoverOP
RollMover::fresh_instance() const
{
	return moves::MoverOP( new RollMover );
}

/// @brief required in the context of the parser/scripting scheme
moves::MoverOP
RollMover::clone() const
{
	return moves::MoverOP( new RollMover( *this ) );
}

/// @brief
RollMover::RollMover(
) : Mover()
{
	moves::Mover::type( "RollMover" );
}

RollMover::RollMover(
	core::Size start_res,
	core::Size stop_res,
	core::Real min_angle,
	core::Real max_angle,
	numeric::xyzVector< core::Real > axis,
	numeric::xyzVector< core::Real > translate
):
	Mover(),
	start_res_(start_res),
	stop_res_(stop_res),
	min_angle_(min_angle),
	max_angle_(max_angle),
	axis_(axis),
	translate_(translate),
	random_roll_(false)
{
	moves::Mover::type( "RollMover" );
}

RollMover::~RollMover()= default;

std::string RollMover::get_name() const {
	return mover_name();
}

std::string RollMover::mover_name() {
	return "RollMover";
}

void RollMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	XMLSchemaRestriction axis;
	axis.name( "axis_char");
	axis.base_type( xs_string );
	axis.add_restriction( xsr_pattern, "[xyz]" );
	xsd.add_top_level_element( axis );
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default("start_res", xsct_non_negative_integer, "First residue id of object to roll", "1")
		+ XMLSchemaAttribute("stop_res", xsct_non_negative_integer, "Last residue ID of object to roll" )
		+ XMLSchemaAttribute("chain", xsct_positive_integer, "Chain ID of object to roll")
		+ XMLSchemaAttribute::attribute_w_default("random_roll", xsct_rosetta_bool, "Rotate and/or translate pose over random axis/direction", "false")
		+ XMLSchemaAttribute::attribute_w_default("random_roll_angle_mag", xsct_real, "Standard deviation for a gaussium magnitude rotation about a random axis", "3.0")
		+ XMLSchemaAttribute::attribute_w_default("random_roll_trans_mag", xsct_real, "Standard deviation for a 3D gaussian random translation", "1.0")
		+ XMLSchemaAttribute("min_angle", xsct_real, "Minimum angle to roll about axis, requierd if random_roll is not set to true" )
		+ XMLSchemaAttribute("max_angle", xsct_real, "Maximum angle to roll about axis, required if random_roll is not set to true" )
		+ XMLSchemaAttribute("axis", xsct_char, "Which axis (x, y, or z) to roll around");

	AttributeList subelement_attributes;
	rosetta_scripts::attributes_for_parse_xyz_vector( subelement_attributes );
	XMLSchemaSimpleSubelementList subelements;
	subelements
		.add_simple_subelement( "axis", subelement_attributes, "Specifies XYZ vector defining axis" )
		.add_simple_subelement( "translate", subelement_attributes, "Specifies XYZ vector defining a 3D translation for the pose" );

	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, mover_name(), "Rigid body mover that rolls the pose about the specified axis", attlist, subelements );
}

std::string RollMoverCreator::keyname() const {
	return RollMover::mover_name();
}

protocols::moves::MoverOP
RollMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new RollMover );
}

void RollMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RollMover::provide_xml_schema( xsd );
}


} //rigid
} //protocols
