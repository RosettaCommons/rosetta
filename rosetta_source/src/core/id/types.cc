// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/id/types.cc
/// @brief  core::id package type declarations
/// @author Matthew O'Meara (mattjomeara@gmail.com)

// Package headers
#include <core/id/types.hh>
#include <utility/exit.hh>

// C++ Headers
#include <sstream>
namespace core {
namespace id {


std::string
DOF_Type_to_string(
	DOF_Type dof_type
) {
	switch(dof_type){
	case PHI:
		return "PHI";
	case THETA:
		return "THETA";
	case D:
		return "D";
	case RB1:
		return "RB1";
	case RB2:
		return "RB2";
	case RB3:
		return "RB3";
	case RB4:
		return "RB4";
	case RB5:
		return "RB5";
	case RB6:
		return "RB6";
	default:
		std::stringstream err_msg;
		err_msg
			<< "Unrecognized DOF_Type '" << dof_type << "'";
		utility_exit_with_message(err_msg.str());
	}
}

DOF_Type
string_to_DOF_Type(
	std::string const & dof_type
) {
	if(dof_type == "PHI") {
		return PHI;
	}	else if(dof_type == "THETA") {
		return THETA;
	}	else if(dof_type == "D") {
		return D;
	}	else if(dof_type == "RB1") {
		return RB1;
	}	else if(dof_type == "RB2") {
		return RB2;
	}	else if(dof_type == "RB3") {
		return RB3;
	}	else if(dof_type == "RB4") {
		return RB4;
	}	else if(dof_type == "RB5") {
		return RB5;
	}	else if(dof_type == "RB6") {
		return RB6;
	} else {
		std::stringstream err_msg;
		err_msg
			<< "Unrecognized DOF_Type '" << dof_type << "'";
		utility_exit_with_message(err_msg.str());
	}
}

std::string
TorsionType_to_string(
	TorsionType torsion_type
) {
	switch(torsion_type){
	case BB:
		return "BB";
	case CHI:
		return "CHI";
	case JUMP:
		return "JUMP";
	default:
		std::stringstream err_msg;
		err_msg
			<< "Unrecognized TosionType '" << torsion_type << "'";
		utility_exit_with_message(err_msg.str());
	}
}

TorsionType
string_to_TorsionType(
	std::string const & torsion_type
) {
	if(torsion_type == "BB") {
		return BB;
	} else if(torsion_type == "CHI") {
		return CHI;
	} else if(torsion_type == "JUMP") {
		return JUMP;
	} else {
		std::stringstream err_msg;
		err_msg
			<< "Unrecognized TorsionType '" << torsion_type << "'";
		utility_exit_with_message(err_msg.str());
	}
}

} // namespace id
} // namespace core
