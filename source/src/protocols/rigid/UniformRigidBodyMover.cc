// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/rigid/UniformRigidBodyMover.cc
/// @author Justin Porter

// Unit Headers
#include <protocols/rigid/UniformRigidBodyMover.hh>
#include <protocols/rigid/UniformRigidBodyMoverCreator.hh>

// Package headers
#include <core/kinematics/Jump.hh>
// Project headers

//Utility Headers
#include <utility/excn/Exceptions.hh>
#include <utility/tag/Tag.hh>

#include <numeric/random/random.hh>
#include <numeric/random/random.functions.hh>

// tracer
#include <basic/Tracer.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

// C++ Headers

// ObjexxFCL Headers

namespace protocols {
namespace rigid {

int const NO_JUMP = 0;

static basic::Tracer tr( "protocols.rigid.UniformRigidBodyMover", basic::t_info );

using namespace core::environment;
using namespace protocols::environment;

// creator
// XRW TEMP std::string
// XRW TEMP UniformRigidBodyMoverCreator::keyname() const {
// XRW TEMP  return UniformRigidBodyMover::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP UniformRigidBodyMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new UniformRigidBodyMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP UniformRigidBodyMover::mover_name() {
// XRW TEMP  return "UniformRigidBodyMover";
// XRW TEMP }

UniformRigidBodyMover::UniformRigidBodyMover():
	ThermodynamicMover(),
	target_jump_( NO_JUMP ),
	rotation_mag_( 3.0 ),
	translation_mag_( 8.0 )
{}

UniformRigidBodyMover::UniformRigidBodyMover( JumpNumber target_jump,
	core::Real rotation_mag,
	core::Real translation_mag ):
	ThermodynamicMover(),
	target_jump_( target_jump ),
	rotation_mag_( rotation_mag ),
	translation_mag_( translation_mag )
{}

void UniformRigidBodyMover::apply( core::pose::Pose& pose ){
	using namespace numeric;
	using namespace numeric::random;

	if ( target_jump_ == NO_JUMP ) {
		std::ostringstream ss;
		ss << "The target jump number of " << this->get_name()
			<< " was " << NO_JUMP << ", which probably means this value wasn't set properly.  "
			<< "If you're using RosettaScripts, try using the 'target_jump' option." << std::endl;
		throw CREATE_EXCEPTION(utility::excn::BadInput,  ss.str() );
	}

	core::kinematics::Jump flexible_jump = pose.jump( target_jump_ );

	xyzMatrix< core::Real> const rot=flexible_jump.get_rotation();
	xyzVector< core::Real> const trans=flexible_jump.get_translation();

	xyzVector<core::Real> delta_trans = random_translation( translation_mag_, numeric::random::rg() );

	core::Real theta = random_rotation_angle<core::Real>( rotation_mag_, numeric::random::rg() );

	xyzVector<core::Real> axis = random_point_on_unit_sphere<core::Real>( numeric::random::rg() );

	xyzMatrix<core::Real> delta_rot = rotation_matrix_radians( axis, theta );

	flexible_jump.set_translation( delta_trans + trans );
	flexible_jump.set_rotation( delta_rot*rot );

	pose.set_jump( (int) target_jump_, flexible_jump );
}

void UniformRigidBodyMover::jump_number( JumpNumber jnum ) {
	target_jump_ = jnum;
}

UniformRigidBodyMover::JumpNumber UniformRigidBodyMover::jump_number() const {
	return target_jump_;
}

void UniformRigidBodyMover::parse_my_tag( utility::tag::TagCOP tag,
	basic::datacache::DataMap&,
	protocols::filters::Filters_map const&,
	protocols::moves::Movers_map const&,
	core::pose::Pose const& ) {
	rotation_mag_ = tag->getOption< core::Real >( "rotation_magnitude", 8.0 );
	translation_mag_ = tag->getOption< core::Real >( "translation_magnitude", 3.0 );
	target_jump_ = tag->getOption< JumpNumber >("target_jump", NO_JUMP );
}

// XRW TEMP std::string UniformRigidBodyMover::get_name() const {
// XRW TEMP  return "UniformRigidBodyMover";
// XRW TEMP }

utility::vector1<core::id::TorsionID_Range>
UniformRigidBodyMover::torsion_id_ranges( core::pose::Pose & ) {
	return utility::vector1<core::id::TorsionID_Range>();
}

moves::MoverOP UniformRigidBodyMover::fresh_instance() const {
	return moves::MoverOP( new UniformRigidBodyMover() );
}

moves::MoverOP UniformRigidBodyMover::clone() const{
	return moves::MoverOP( new UniformRigidBodyMover( *this ) );
}

std::string UniformRigidBodyMover::get_name() const {
	return mover_name();
}

std::string UniformRigidBodyMover::mover_name() {
	return "UniformRigidBodyMover";
}

void UniformRigidBodyMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	// TO DO!
	using namespace utility::tag;
	AttributeList attlist; // TO DO: add attributes to this list

	XMLSchemaComplexTypeGeneratorOP ct_gen = complex_type_gen();
	ct_gen->description( "XRW TO DO" );
	ct_gen->add_attributes( attlist );
	ct_gen->element_name( mover_name() );
	ct_gen->write_complex_type_to_schema( xsd );

	//protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TO DO", attlist );
}

utility::tag::XMLSchemaComplexTypeGeneratorOP
UniformRigidBodyMover::complex_type_gen()
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default("rotation_magnitude", xsct_real, "XRW TO DO", "8.0")
		+ XMLSchemaAttribute::attribute_w_default("translation_magnitude", xsct_real, "XRW TO DO", "3.0")
		+ XMLSchemaAttribute::attribute_w_default("target_jump", xs_integer, "XRW TO DO", utility::to_string(NO_JUMP));


	XMLSchemaComplexTypeGeneratorOP ct_gen ( new XMLSchemaComplexTypeGenerator );
	ct_gen->add_attributes( attlist );
	ct_gen->complex_type_naming_func( & protocols::moves::complex_type_name_for_mover );
	ct_gen->add_optional_name_attribute();
	return ct_gen;
}

/*
void UniformRigidBodyMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
using namespace utility::tag;
XMLSchemaComplexTypeGeneratorOP ct_gen = complex_type_gen();
ct_gen->description( "XRW TO DO" );
ct_gen->add_attributes( attlist );
ct_gen->element_name( mover_name() );
ct_gen->write_complex_type_to_schema( xsd );
}
*/

std::string UniformRigidBodyMoverCreator::keyname() const {
	return UniformRigidBodyMover::mover_name();
}

protocols::moves::MoverOP
UniformRigidBodyMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new UniformRigidBodyMover );
}

void UniformRigidBodyMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	UniformRigidBodyMover::provide_xml_schema( xsd );
}


} // rigid
} // protocols
