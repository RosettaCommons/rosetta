// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/environment/EnvExcn.cc
/// @author Justin R Porter

// Unit Headers
#include <protocols/environment/EnvExcn.hh>

// Package headers
#include <protocols/environment/ProtectedConformation.hh>

// Project headers
#include <core/pose/Pose.hh>

// tracer
#include <basic/Tracer.hh>

// C++ Headers

// ObjexxFCL Headers

static thread_local basic::Tracer tr( "protocols.environment.EnvExcn", basic::t_info );

namespace protocols {
namespace environment {

std::string dof_id_to_string( core::id::DOF_ID const& id, ProtectedConformation const& conf ){
	std::ostringstream ss;
	if ( id.type() == core::id::PHI ) {
		ss << "Torsion Angle ";
	} else if ( id.type() == core::id::THETA ) {
		ss << "Bond Angle ";
	} else if ( id.type() == core::id::D ) {
		ss << "Bond Length ";
	} else if ( id.type() >= core::id::RB1 &&
			id.type() <= core::id::RB6 ) {
		ss << "RB" << id.type() - core::id::RB1 + 1;
	}
	ss << " on " << conf.residue( id.rsd() ).name3() << id.rsd() << " owned by atom " << id.atomno() << "("
		<< conf.residue( id.rsd() ).atom_name( (int) id.atomno() ) << ")";
	return ss.str();
}


EXCN_Env_Security_Exception::EXCN_Env_Security_Exception(
	core::id::DOF_ID const& dof,
	core::environment::DofPassportCOP pass,
	std::string const& mod_type,
	ProtectedConformation const& conf,
	std::string const& mover_name,
	EnvironmentCAP env
):
	Parent( env ),
	id_( dof ),
	pass_( pass )
{
	std::ostringstream msg;

	msg << "ProtectedConformation reported illegal access to '"
		<< dof_id_to_string( id(), conf )
		<< "' by mover '" << mover_name << "' via '" << mod_type << "'. ";

	if ( passport() ) {
		core::Size const seqpos = id().atom_id().rsd();
		msg << "The passport did not contain the DoF. Dofs of this type available for this residue were: [ ";
		for ( core::Size i = 1; i <= conf.residue( seqpos ).natoms(); ++i ) {
			core::id::DOF_ID const other_id( core::id::AtomID( i, seqpos ), id().type() );
			if ( passport()->dof_access( other_id ) ) {
				msg << i << "(" << conf.residue( seqpos ).atom_name( i ) << "), ";
			}
		}
		msg << ']';
	} else {
		msg << "No passport was found.";
	}

	add_msg( msg.str() );
}

EXCN_Env_Security_Exception::EXCN_Env_Security_Exception(
	std::string const& message,
	std::string const& mover_name,
	EnvironmentCAP env ) :
	Parent( env )
{
	std::ostringstream msg;
	msg << message << std::endl << "Failure due to a call by mover '" << mover_name << "'";
	add_msg( msg.str() );
}



} // environment
} // protocols
