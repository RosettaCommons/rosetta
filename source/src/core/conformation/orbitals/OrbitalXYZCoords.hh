// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/conformation/orbitals/OrbitalXYZCoords.hh
/// @brief Created on: Jun 30, 2011
/// @author combss

#ifndef INCLUDED_core_conformation_orbitals_OrbitalXYZCoords_HH
#define INCLUDED_core_conformation_orbitals_OrbitalXYZCoords_HH

// Project headers
#include <core/types.hh>

// Utility headers
#include <numeric/xyzVector.hh>


namespace core {
namespace conformation {
namespace orbitals {

class OrbitalXYZCoords {
public:
	/// @brief default constructor and set atom type number to 0 and place the
	/// atom at the origin
	OrbitalXYZCoords():
		xyz_( 0.0 ),
		type_(0)
	{}

	/// @brief constructor with an atom type number
	// type is set at construction time -- atom is placed at the origin
	OrbitalXYZCoords( Vector const & xyz_in, ShortSize const & type):
		xyz_( xyz_in ),
		type_(type)
	{}


	/// @brief destructor
	~OrbitalXYZCoords() {}


	Vector const &
	xyz() const
	{
		return xyz_;
	}


	void
	xyz( Vector const & xyz_in )
	{
		xyz_ = xyz_in;
	}

	ShortSize const &
	type() const
	{
		return type_;
	}

	void
	type( ShortSize const & type_in )
	{
		type_ = type_in;
	}

#ifdef    SERIALIZATION
	/// @brief Serialization method
	template < class Archive >
	void
	save( Archive & arch ) const;

	/// @brief De-serialization method
	template < class Archive >
	void
	load( Archive & arch );
#endif // SERIALIZATION

private:
	/// xyz coordinates
	Vector xyz_;
	ShortSize type_;


};


}
}
}


#endif /* ORBITALXYZCOORDS_HH_ */
