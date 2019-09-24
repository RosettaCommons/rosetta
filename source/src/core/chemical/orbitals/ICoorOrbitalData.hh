// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//////////////////////////////////////////////////////////////////////
///
/// @brief
/// A class for reading in the orbital type properties
///
/// @details
/// This class contains the ORBITALS INTERNAL_ICOOR data that is read in from residue_io.cc. Actually,
/// the data is set when residue_io.cc calls the command from residuetype.cc set_orbital_icoor. The data
/// is set and chills out in memory until needed. The actual xyz coordinates are not stored at this point.
/// xyz coordinates are constructed when conformation/Residue.hh asks for the build function in this class.
/// At that point, the coordinates are built and returned.
///
/// But wait, you say, why do you store the names of the atoms instead of the index of the atoms!? Well, the
/// problem occurs when residuetype reorders the indices of the atoms. When this occurrs, the indices for the icoor
/// are not reordered here. Another problem ocurs because orbital indices are not reordered in this process because
/// orbital indices are seperate from the atom indices. Regardless, when you build the xyz coords, this step is transparent
/// because the function orbital_xyz() in residue.hh takes care of this conversion of indice to string.
///
/// @note NOTE!!!!!!!!!!! The internal coordinates cannot contain an orbital as the stub1, stub2, or stub3 atom.
/// This is because the xyz coordinates are not updated when the conformation changes. The stub1, stub2, stub2 atoms
/// must be actual atoms and not orbitals!!!!! (design feature or flaw? you decide)
///
///
/// @author
/// Steven Combs
///
///
/////////////////////////////////////////////////////////////////////////

#ifndef INCLUDED_core_chemical_orbitals_ICoorOrbitalData_hh
#define INCLUDED_core_chemical_orbitals_ICoorOrbitalData_hh

// Project headers

#include <core/chemical/orbitals/ICoorOrbitalData.fwd.hh>

// Utility headers

#include <core/types.hh>

#include <string>

namespace core {
namespace chemical {
namespace orbitals {

class ICoorOrbitalData{
public:

public:
	/// @brief default constructor
	ICoorOrbitalData();


	ICoorOrbitalData(
		Real phi,
		Real theata,
		Real distance,
		std::string const & stub1,
		std::string const & stub2,
		std::string const & stub3
	);

	Real phi() const;
	Real theta() const;
	Real distance() const;

	Vector
	build(
		Vector stub1_xyz,
		Vector stub2_xyz,
		Vector stub3_xyz
	) const;


	void replace_stub1(std::string const & atom1 )
	{
		stub1_=atom1;
	}

	void replace_stub2(std::string const & atom2)
	{
		stub2_=atom2;
	}
	void replace_stub3(std::string const & atom3)
	{
		stub3_=atom3;
	}

	std::string const & get_stub1() const
	{
		return stub1_;
	}
	std::string const & get_stub2() const
	{
		return stub2_;
	}
	std::string const & get_stub3() const
	{
		return stub3_;
	}


private:
	Real phi_;
	Real theta_;
	Real distance_;

	std::string stub1_;
	std::string stub2_;
	std::string stub3_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


}
}
}


#endif /* ICOORORBITALID_HH_ */
