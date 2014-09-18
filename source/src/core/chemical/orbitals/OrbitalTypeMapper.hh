// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


#ifndef INCLUDED_core_chemical_orbitals_OrbitalTypeMapper_hh
#define INCLUDED_core_chemical_orbitals_OrbitalTypeMapper_hh

#include <core/chemical/orbitals/OrbitalTypeMapper.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// C++ Headers
#include <map>

#ifdef MULTI_THREADED
#ifdef CXX11
// C++11 Headers
#include <atomic>
#include <mutex>
#endif
#endif

namespace core {
namespace chemical {
namespace orbitals {

class OrbitalTypeMapper
{
public:

	virtual ~OrbitalTypeMapper();

	static OrbitalTypeMapper * get_instance();

	orbital_type_enum get_orbital_enum(std::string & orbital_type_name);


	//orbital_type_enum get_orbital_type_enum(std::string const & orbital_type_name)

#ifdef MULTI_THREADED
#ifdef CXX11
public:

	/// @brief This public method is meant to be used only by the
	/// utility::thread::safely_create_singleton function and not meant
	/// for any other purpose.  Do not use.
	static std::mutex & singleton_mutex();

private:

#endif
#endif

private:
	OrbitalTypeMapper();
	//unimplemented -- uncopyable
	OrbitalTypeMapper(OrbitalTypeMapper const & );
	OrbitalTypeMapper const & operator = (OrbitalTypeMapper const &);

	/// @brief private singleton creation function to be used with
	/// utility::thread::threadsafe_singleton
	static OrbitalTypeMapper * create_singleton_instance();

	void map_orbital_name_to_enum();

private:
	/// @brief static data member holding pointer to the singleton class itself
#if defined MULTI_THREADED && defined CXX11
	static std::atomic< OrbitalTypeMapper * > instance_;
#else
	static OrbitalTypeMapper * instance_;
#endif


	std::map<std::string, orbital_type_enum> orbital_type_2_enum_;

};

}
}
}


#endif /* ORBITALTYPEMAPPER_HH_ */
