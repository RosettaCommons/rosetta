// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/conformation/Atom.hh
/// @brief  Class definitions for conformation::Atom
/// @note   not to be confused with chemical::Atom
/// @author Phil Bradley
/// @author Steven Combs - some comments

#ifndef INCLUDED_core_conformation_Atom_HH
#define INCLUDED_core_conformation_Atom_HH


// Unit headers
#include <core/conformation/Atom.fwd.hh>

// Project headers
#include <core/types.hh>

// Numeric headers
#include <numeric/xyzVector.hh>

namespace core {
namespace conformation {

/// @details This atom class differs from the chameical::Atom, AtomType, and AtomTypeSet classes in that it only
/// contains the information about the xyz coordinates.  In the future, it might contain information about such things
/// as B-factor.  This information is generally initialized and set by conformation::Residue.
/// @note    chemical::Atoms are stored in chemical::ResidueType (within its ResidueGraph);
/// conformation::Atoms are stored in conformation::Residue
class Atom {

public:
	/// @brief   Default constructor
	/// @details Set atom type number to 0 and place the atom at the origin.
	Atom():
		xyz_( 0.0 ),
		type_( 0 ),
		mm_type_( 0 )
	{}

	/// @brief constructor with an atom type number
	// type is set at construction time -- atom is placed at the origin
	Atom( ShortSize const type_in, ShortSize const mm_type_in ):
		xyz_( 0.0 ),
		type_( type_in ),
		mm_type_( mm_type_in )
	{}

	/// @brief constructor with xyz and an atom type number
	// Atom( Vector const & xyz_in, int const type_in, Real temperature = 0.0 ):
	Atom( Vector const & xyz_in, ShortSize const type_in, ShortSize const mm_type_in ):
		xyz_( xyz_in ),
		type_( type_in ),
		mm_type_( mm_type_in )
	{}

	/// @brief destructor
	virtual
	~Atom() {}


	/// @brief  Generate string representation of conformation::Atom for debugging purposes.
	void show( std::ostream & output=std::cout ) const;


	/// @brief set the atom type number
	void type( ShortSize const type_in )
	{
		type_ = type_in;
	}

	/// @brief Returns the AtomType number
	ShortSize
	type() const
	{
		return type_;
	}

	/// @brief set the mm atom type number
	void mm_type( ShortSize const mm_type_in )
	{
		mm_type_ = mm_type_in;
	}

	/// @brief get the mm atom type number
	ShortSize
	mm_type() const
	{
		return mm_type_;
	}

	/// @brief Returns the atom coordinates as an xyzVector
	Vector const &
	xyz() const
	{
		return xyz_;
	}

	/// @brief Sets the atom coordinates using an xyzVector
	void
	xyz( Vector const & xyz_in )
	{
		xyz_ = xyz_in;
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
	// data

	// position
	Vector xyz_;

	// rosetta atom_type, stored in conformation::Residue for the sake of speed in the Lennard-Jones calculation
	ShortSize type_;

	// mm_atom_type, stored in conformation::Residue for the sake of speed in the the Charmm Lennard-Jones calculation
	ShortSize mm_type_;

	// b-factor from pdb file
	// Real temperature_;
};

std::ostream & operator << ( std::ostream & out, Atom const & atom );

}  // conformation
}  // core

#endif  // INCLUDED_core_conformation_Atom_HH
