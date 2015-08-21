// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//////////////////////////////////////////////////////////////////////
///
/// @brief
/// A class for reading in the orbital type properties
///
/// @details
/// This class reads in the orbital_properties.txt file which contains the "chemical" information for orbitals.
/// This does not contain the actual properties, but sets the properties through the OrbitalType class.
/// This class is called by the ChemicalManager. Modeled off of atomtypeset.
///
///
///
/// @author
/// Steven Combs
///
///
/////////////////////////////////////////////////////////////////////////

#ifndef INCLUDED_core_chemical_orbitals_OrbitalTypeSet_hh
#define INCLUDED_core_chemical_orbitals_OrbitalTypeSet_hh


// Project headers
#include <core/chemical/orbitals/OrbitalType.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <map>
#include <core/types.hh>
#include <utility/vector1_bool.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace chemical {
namespace orbitals {


class OrbitalTypeSet : public utility::pointer::ReferenceCount {


public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~OrbitalTypeSet();
	OrbitalTypeSet(std::string const & directory);

	void read_file(std::string const & filename);

	/// @brief [ ] operator, simulating vector index behavior
	///
	/// @details look up an OrbitalTypeSet by 1-based indexing
	///
	OrbitalType const &
	operator[](core::Size const index) const
	{
		return *(orbitals_[index]);
	}


	/// @brief lookup the orbital type by the orbital type name string
	int
	orbital_type_index( std::string const & orbital_type_name ) const;

	/// @brief lookup the orbital type by the orbital type name string
	int
	orbital_type_index( std::string & orbital_type_name ) const;

private:
	/// lookup map: get orbital_type_index by orbital_type_name
	std::map< std::string, int > orbital_type_index_;


	/// @brief  Save the directory name for future use
	std::string directory_;

	/// @brief a collection of OrbitalTypes,
	///
	/// @details OrbitalType has data of atom properties, and it can be
	/// looked up by orbital_type_index.
	utility::vector1< OrbitalType* > orbitals_;


};


}
}
}


#endif /* ORBITALTYPESET_HH_ */
