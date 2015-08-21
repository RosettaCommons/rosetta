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

// Utility headers

#include <core/types.hh>
#include <string>
#include <core/chemical/ResidueGraphTypes.hh>

//Auto Headers
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
		Size stub1,
		Size stub2,
		Size stub3,
		VD v1,
		VD v2,
		VD v3
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


	void replace_stub1(const core::Size atom1 )
	{
		stub1_=atom1;
	}

	void replace_stub2(const core::Size atom2)
	{
		stub2_=atom2;
	}
	void replace_stub3(const core::Size atom3)
	{
		stub3_=atom3;
	}

	core::Size get_stub1() const
	{
		return stub1_;
	}
	core::Size get_stub2() const
	{
		return stub2_;
	}
	core::Size get_stub3() const
	{
		return stub3_;
	}

	VD vertex1() const
	{
		return vertex1_;
	}
	VD vertex2() const
	{
		return vertex2_;
	}
	VD vertex3() const
	{
		return vertex3_;
	}

	void vertex1(VD const vertex)
	{
		vertex1_=vertex;
	}
	void vertex2(VD const vertex)
	{
		vertex2_=vertex;
	}
	void vertex3(VD const vertex)
	{
		vertex3_=vertex;
	}


private:
	Real phi_;
	Real theta_;
	Real distance_;

	Size stub1_;
	Size stub2_;
	Size stub3_;

	VD vertex1_;
	VD vertex2_;
	VD vertex3_;


};


}
}
}


#endif /* ICOORORBITALID_HH_ */
