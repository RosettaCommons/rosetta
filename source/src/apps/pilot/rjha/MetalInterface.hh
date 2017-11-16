// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/rjha/MetalInterface.hh
/// @brief This header contains helper functions common to metal interface pilot_apps
/// @author Steven Lewis

// Unit Headers

// Project Headers

#include <core/conformation/Residue.hh>

// Numeric Headers
#include <numeric/xyzVector.hh>

// Utility Headers
#include <basic/Tracer.hh>
using basic::Warning;

// C++ headers
#include <string>

//Auto Headers
#include <core/chemical/AtomType.hh>


//tracers

//local enum - this names the (hopefully 5) residues in the input from RosettaMatch
enum MatchPosition {
	p1_1 = 1, //partner 1, residue 1
	p1_2 = 2, //partner 1, residue 2
	p1_3 = 3, //partner 1, residue 3
	p2 = 4,   //residue on partner2 - replaces partner2_residue in that pose
	metal = 5
};

//local helper function - iterates over all sidechain non-carbon heavy atoms of res to find the one closest to xyz
std::string const find_closest_atom( core::conformation::Residue const & res, core::Vector const & xyz ){
	core::Size index(0);
	core::Real dis(1000000000); //arbitrary large value
	core::Real temp_dis(dis);

	//atom ordering rules are in ResidueType.hh; briefly sidechain heavy atoms are between first_sidechain_heavyatom and nheavyatoms
	for ( core::Size i(res.first_sidechain_atom()); i <= res.nheavyatoms(); ++i ) {
		temp_dis = res.xyz(i).distance(xyz);
		T("apps.pilot.rhja.find_closest_atom") << "testing atom name " << res.atom_name(i) << std::endl;
		//T("apps.pilot.rhja.find_closest_atom") << temp_dis << std::endl;
		if ( res.is_virtual(i) ) {
			T("apps.pilot.rhja.find_closest_atom") << "skipping virtual atom " << res.atom_name(i) << std::endl;
			continue;
		}
		if ( (temp_dis < dis) && (res.atom_type(i).element() != "C") ) {
			dis = temp_dis;
			index = i;
			//T("apps.pilot.rhja.find_closest_atom") << "accepting change" << std::endl;
		}//if shorter and non-carbon
	}//for all sidechain heavyatoms
	if ( !(res.atom_type(index).element() == "N") && !(res.atom_type(index).element() == "S") ) {
		Warning() << "irregular ligand atom - not N or S" << std::endl;
	}

	T("apps.pilot.rhja.find_closest_atom") << "returning " << res.atom_name(index) << std::endl;
	return res.atom_name(index);
}
