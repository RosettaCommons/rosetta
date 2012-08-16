// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//////////////////////////////////////////////////////////////////////
/// @begin Atom.hh
///
/// @brief
/// The conformation class of the atom object
///
/// @detailed
/// This atom class differs from the AtomType/AtomTypeSet class in that it only contains the
/// information about the xyz coordinates, the index of the atom, and the index of the atomtype.
/// This information is generally initialized and set by conformation/Residue.hh.
/// Steven Combs - comments
///
/// @authors
/// Phil Bradley
///
///
///
/// @last_modified December 9 2010
/////////////////////////////////////////////////////////////////////////

#ifndef INCLUDED_core_conformation_Atom_hh
#define INCLUDED_core_conformation_Atom_hh


// Unit headers
#include <core/conformation/Atom.fwd.hh>

// Project headers
#include <core/types.hh>

// Utility headers
#include <numeric/xyzVector.hh>

// C++ headers

namespace core {
namespace conformation {

/// A simple object with atom's position and its chemical type
class Atom {

public:
	/// @brief default constructor and set atom type number to 0 and place the
	/// atom at the origin
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

	// /// @brief access temperature Real
	// Real temperature() const {
	// 	return temperature_;
	// }
	//
	// /// @brief set temperature Real
	// void
	// temperature( Real temp ) {
	// 	temperature_ = temp;
	// }

	// data
private:

	/// position
	Vector xyz_;

	/// rosetta atom_type
	ShortSize type_;

	/// mm_atom_type
	ShortSize mm_type_;

	/// b-factor from pdb file
	// Real temperature_;
};

std::ostream & operator << ( std::ostream & os, Atom const & atom );

} // pose
} // core



#endif // INCLUDED_core_pose_Residues_HH
