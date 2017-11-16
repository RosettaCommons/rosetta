// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file XYZClaim
/// @brief Claims access to a torsional angle.
/// @author Justin Porter

// Unit Headers
#include <protocols/environment/claims/XYZClaim.hh>

// Package Headers
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/tree/Atom.hh>

#include <core/environment/LocalPosition.hh>
#include <core/environment/SequenceAnnotation.hh>

#include <protocols/environment/claims/EnvClaim.hh>
#include <protocols/environment/claims/TorsionClaim.hh>
#include <protocols/environment/claims/ClaimStrength.hh>
#include <protocols/environment/claims/EnvLabelSelector.hh>

#include <protocols/environment/ClientMover.hh>
#include <protocols/environment/ProtectedConformation.hh>

// Project Headers
#include <core/id/TorsionID.hh>

// ObjexxFCL Headers

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <basic/Tracer.hh>

// C++ headers

// option key includes


static basic::Tracer tr( "protocols.environment.XYZClaim", basic::t_info );

namespace protocols {
namespace environment {
namespace claims {

using core::environment::LocalPosition;
using core::environment::LocalPositions;

XYZClaim::XYZClaim( ClientMoverOP owner,
	utility::tag::TagCOP tag,
	basic::datacache::DataMap const& datamap ):
	EnvClaim( owner ),
	c_str_( Parent::parse_ctrl_str( tag->getOption< std::string >( "control_strength" ) ) ),
	i_str_( Parent::parse_ctrl_str( tag->getOption< std::string >( "initialization_strength", "DOES_NOT_CONTROL" ) ) ),
	bRelative_( tag->getOption( "relative_only", false ) )
{
	std::string const& selection = tag->getOption< std::string >( "selection" );
	if ( datamap.has( "ResidueSelector", selection ) ) {
		tr.Debug << "Pulling ResidueSelector '" << selection << "' from datamap in " << *this << std::endl;
		selector_ = datamap.get_ptr< core::select::residue_selector::ResidueSelector const >( "ResidueSelector", tag->getOption<std::string>( "selection" ) );
	} else {
		tr.Debug << "Instantiating EnvLabelSelector for selection '" << selection << "' in " << *this << std::endl;
		selector_ = ResidueSelectorCOP( ResidueSelectorOP( new EnvLabelSelector( LocalPosition( selection ) ) ) );
	}
}

//TODO: Remove me
XYZClaim::XYZClaim(
	ClientMoverOP owner,
	LocalPosition const& local_pos
):
	EnvClaim( owner ),
	selector_( ResidueSelectorCOP( ResidueSelectorOP( new EnvLabelSelector( local_pos ) ) ) ),
	c_str_( DOES_NOT_CONTROL ),
	i_str_( DOES_NOT_CONTROL ),
	bRelative_( false )
{}

//TODO: remove me
XYZClaim::XYZClaim(
	ClientMoverOP owner,
	std::string const& label,
	std::pair< core::Size, core::Size > const& range
):
	EnvClaim( owner ),
	selector_( ResidueSelectorCOP( ResidueSelectorOP( new EnvLabelSelector( label, range ) ) ) ),
	c_str_( DOES_NOT_CONTROL ),
	i_str_( DOES_NOT_CONTROL ),
	bRelative_( false )
{}

XYZClaim::XYZClaim(
	ClientMoverOP owner,
	core::select::residue_selector::ResidueSelectorCOP selector
):
	EnvClaim( owner ),
	selector_( ResidueSelectorCOP( selector ) ),
	c_str_( DOES_NOT_CONTROL ),
	i_str_( DOES_NOT_CONTROL ),
	bRelative_( false )
{}


DOFElement XYZClaim::wrap_dof_id( core::id::DOF_ID const& id ) const {
	DOFElement e = Parent::wrap_dof_id( id );

	e.c_str = ctrl_strength();
	e.i_str = init_strength();

	return e;
}

void XYZClaim::yield_elements( core::pose::Pose const& pose, DOFElements& elements ) const {
	if ( ctrl_strength() == DOES_NOT_CONTROL &&
			init_strength() == DOES_NOT_CONTROL ) {
		tr.Warning << "XYZClaim owned by " << owner()->get_name()
			<< " has both initializaiton and sampling strength set to DOES_NOT_CONTROL."
			<< "  Did you forget to set these values?" << std::endl;
	}

	utility::vector1< bool > selection = selector()->apply( pose );

	if ( tr.Debug.visible() ) {
		tr.Debug << *this << " selection for DoF Claiming: ";
		for ( core::Size i = 1; i <= selection.size(); ++i ) {
			tr.Debug << ( selection[i] ? "T" : "F" );
		}
		tr.Debug << std::endl;
	}

	// Run through each position in the pose and claim produce DoFElements claiming it,
	// iff it's selected
	for ( Size seqpos = 1; seqpos <= selection.size(); ++seqpos ) {
		if ( selection[seqpos] ) {
			for ( Size i = 1; i <= pose.conformation().residue( seqpos ).atoms().size(); ++i ) {
				core::id::AtomID const atom_id( i, seqpos );

				// if the relative setting is activated, only DoFs that build *relative* positions are claimed.
				// in other words, if the parent is outside the selection, don't claim it.
				core::kinematics::tree::AtomCOP parent = pose.conformation().atom_tree().atom( atom_id ).parent();
				if ( !parent ) {
					tr.Debug << *this << " and strengths ctrl=" << ctrl_strength() << " and init=" << init_strength()
						<< " skipping " << atom_id << " because its parent is null." << std::endl;
				} else if ( !relative() || selection[ parent->id().rsd() ] ) {
					if ( pose.conformation().atom_tree().atom( atom_id ).is_jump() ) {
						for ( int j = core::id::RB1; j <= core::id::RB6; ++j ) {
							elements.push_back( wrap_dof_id( core::id::DOF_ID( atom_id,
								core::id::DOF_Type( j ) ) ) );
						}
					} else {
						for ( int j = core::id::PHI; j <= core::id::D; ++j ) {
							elements.push_back( wrap_dof_id( core::id::DOF_ID( atom_id,
								core::id::DOF_Type( j ) ) ) );
						}
					}
				} else {
					tr.Debug << *this << " and strengths ctrl=" << ctrl_strength() << " and init=" << init_strength()
						<< " skipping " << atom_id << " because its parent ("
						<< pose.conformation().atom_tree().atom( atom_id ).parent()->id()
						<< ") doesn't belong to the XYZ selection and it is configured in relative mode." << std::endl;
				}
			}
		}
	}
}

ControlStrength const& XYZClaim::ctrl_strength() const {
	return c_str_;
}

ControlStrength const& XYZClaim::init_strength() const {
	return i_str_;
}

void XYZClaim::strength( ControlStrength const& c_str, ControlStrength const& i_str ){
	if ( c_str > EXCLUSIVE || c_str < DOES_NOT_CONTROL ) {
		throw utility::excn::EXCN_RangeError( "Sampling ControlStrengths are limited to values between DOES_NOT_CONTROL and EXCLUSIVE" );
	} else {
		c_str_ = c_str;
	}

	if ( i_str > EXCLUSIVE || i_str < DOES_NOT_CONTROL ) {
		throw utility::excn::EXCN_RangeError( "Initialization ControlStrengths are limited to values between DOES_NOT_INITIALIZE and EXCLUSIVE" );
	} else {
		i_str_ = i_str;
	}
}

EnvClaimOP XYZClaim::clone() const {
	return EnvClaimOP( new XYZClaim( *this ) );
}


std::string
XYZClaim::class_name(){
	return "XYZClaim";
}

void
XYZClaim::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ){
	using namespace utility::tag;

	XMLSchemaRestriction restr;
	restr.name( "envclaim_ctrl_str" );
	restr.base_type( xs_string );
	restr.add_restriction( xsr_pattern, "((([Dd][Oo][Ee][Ss]_[Nn][Oo][Tt])|([Cc][Aa][Nn])|([Mm][Uu][Ss][Tt]))_[Cc][Oo][Nn][Tt][Rr][Oo][Ll])|([Ee][Xx][Cc][Ll][Uu][Ss][Ii][Vv][Ee])" );
	xsd.add_top_level_element( restr );

	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::required_attribute( "control_strength", "envclaim_ctrl_str", "XRW TO DO" )
		+ XMLSchemaAttribute::attribute_w_default( "initialization_strength", "envclaim_ctrl_str", "XRW TO DO", "DOES_NOT_CONTROL" )
		+ XMLSchemaAttribute::attribute_w_default( "relative_only", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute::required_attribute( "selection", xs_string, "XRW TO DO" );

	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.element_name( class_name() )
		.complex_type_naming_func ( & EnvClaim::envclaim_ct_namer )
		.description( "XRW TO DO" )
		.add_attributes( attlist )
		.write_complex_type_to_schema( xsd );
}


std::string XYZClaim::type() const {
	return "XYZClaim";
}

void XYZClaim::show( std::ostream& os ) const {
	os << type() << " owned by a " << owner()->get_name();
}

} //claims
} //environment
} //protocols
