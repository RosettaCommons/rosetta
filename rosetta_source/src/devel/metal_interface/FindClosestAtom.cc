// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file devel/metal_interface/FindClosestAtom.hh
/// @brief This header contains helper functions common to metal interface pilot_apps
/// @author Steven Lewis

// Unit Headers
#include <devel/metal_interface/FindClosestAtom.hh>
// Project Headers
#include <core/conformation/Residue.hh>
#include <core/chemical/AtomType.hh>
// Numeric Headers
#include <numeric/xyzVector.hh>
// Utility Headers
#include <basic/Tracer.hh>
using basic::T;
using basic::Warning;

// C++ headers
#include <string>

#include <utility/vector1.hh>


//tracers
// using basic::Error;
// using basic::Warning;
// using basic::T;

namespace devel{
namespace metal_interface{

//local helper function - iterates over all sidechain non-carbon heavy atoms of res to find the one closest to xyz
std::string find_closest_atom( core::conformation::Residue const & res, core::Vector const & xyz ){
	core::Size index(0);
	core::Real dis(1000000000); //arbitrary large value
	core::Real temp_dis(dis);

	//atom ordering rules are in ResidueType.hh; briefly sidechain heavy atoms are between first_sidechain_heavyatom and nheavyatoms
	for( core::Size i(res.first_sidechain_atom()); i <= res.nheavyatoms(); ++i){
		temp_dis = res.xyz(i).distance(xyz);
		if (res.is_virtual(i)){
			continue;
		}
		if ( (temp_dis < dis) && (res.atom_type(i).element() != "C") ){ //Carbon is not a metal-coordinating atom
			dis = temp_dis;
			index = i; //accepting change
		}//if shorter and non-carbon
	}//for all sidechain heavyatoms
	if (!(res.atom_type(index).element() == "N") && !(res.atom_type(index).element() == "S") && !(res.atom_type(index).element() == "O")){
		Warning() << "irregular ligand atom - not N or S or O" << std::endl;
	}
	return res.atom_name(index);
}

}//namespace metal_interface
}//namespace devel
