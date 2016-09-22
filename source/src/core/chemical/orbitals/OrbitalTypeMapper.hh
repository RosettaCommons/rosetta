// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


#ifndef INCLUDED_core_chemical_orbitals_OrbitalTypeMapper_hh
#define INCLUDED_core_chemical_orbitals_OrbitalTypeMapper_hh

#include <core/chemical/orbitals/OrbitalTypeMapper.fwd.hh>

// Utility Headers
#include <utility/SingletonBase.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ Headers
#include <map>


namespace core {
namespace chemical {
namespace orbitals {

class OrbitalTypeMapper : public utility::SingletonBase< OrbitalTypeMapper >
{
public:
	friend class utility::SingletonBase< OrbitalTypeMapper >;

	virtual ~OrbitalTypeMapper();

	orbital_type_enum get_orbital_enum(std::string & orbital_type_name);

	//orbital_type_enum get_orbital_type_enum(std::string const & orbital_type_name)


private:
	OrbitalTypeMapper();

	//unimplemented -- uncopyable
	OrbitalTypeMapper(OrbitalTypeMapper const & ) = delete;
	OrbitalTypeMapper const & operator = (OrbitalTypeMapper const &) = delete;

	void map_orbital_name_to_enum();

private:

	std::map<std::string, orbital_type_enum> orbital_type_2_enum_;

};

}
}
}


#endif /* ORBITALTYPEMAPPER_HH_ */
