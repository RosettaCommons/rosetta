// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#include <core/chemical/orbitals/OrbitalTypeMapper.hh>
#include <map>

#include <string>

#include <utility/thread/threadsafe_creation.hh>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>

namespace core {
namespace chemical {
namespace orbitals {

OrbitalTypeMapper::OrbitalTypeMapper()
{
	map_orbital_name_to_enum();
}

OrbitalTypeMapper::~OrbitalTypeMapper() = default;

void OrbitalTypeMapper::map_orbital_name_to_enum()
{
	orbital_type_2_enum_["C.pi.sp2"]=C_pi_sp2;
	orbital_type_2_enum_["N.pi.sp2"]=N_pi_sp2;
	orbital_type_2_enum_["N.p.sp2"]=N_p_sp2;
	orbital_type_2_enum_["O.pi.sp2"]=O_pi_sp2;
	orbital_type_2_enum_["O.p.sp2"]=O_p_sp2;
	orbital_type_2_enum_["O.p.sp3"]=O_p_sp3;
	orbital_type_2_enum_["S.p.sp3"]=S_p_sp3;
	orbital_type_2_enum_["O.pi.sp2.bb"]= O_pi_sp2_bb;
	orbital_type_2_enum_["O.p.sp2.bb"] = O_p_sp2_bb;


}

orbital_type_enum OrbitalTypeMapper::get_orbital_enum(std::string & orbital_type_name)
{
	return orbital_type_2_enum_.find(orbital_type_name)->second;

}


}
}
}
