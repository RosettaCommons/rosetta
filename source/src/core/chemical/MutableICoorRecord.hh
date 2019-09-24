// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief A ICoord record object for a MutableResidueType
/// @author Rocco Moretti (rmorettiase@gmail.com)


#ifndef INCLUDED_core_chemical_MutableICoorRecord_hh
#define INCLUDED_core_chemical_MutableICoorRecord_hh


// Unit headers
#include <core/chemical/MutableICoorRecord.fwd.hh>
#include <core/chemical/MutableResidueType.fwd.hh>
#include <core/chemical/AtomICoor.fwd.hh>

//// Project headers
#include <core/types.hh>

//// Utility headers
#include <utility/exit.hh>

// C++ headers
#include <string>

namespace core {
namespace chemical {

/// @brief A basic class containing basic info of internal coordinates needed for building an atom within a ResidueType
/// This is a simplified representation, used for MutableResidueType. It contains all the information, but is intended to be somewhat easier to update for added/deleted atoms than the standard AtomICoor.
/**
In atom tree, each atom is defined by its internal coordinates, which include a bond distance,
a bond angle and a torsion angle. Of course, all these internal coordinates are only meaningful
in the context of three reference (stub) atoms. MutableICoorRecord information is stored in the residue param
files and some terms are defined as following:\n
- bond distance d_ is that between the atom to be built (child) and stub_atom1 (parent)
- bond angle theta_ is that defined by child-parent-stub2(angle)
- torsion angle phi_ is that defined by child-parent-stub2-stub3(torsion)
*/
class MutableICoorRecord {
public:
	/// @brief default constructor
	MutableICoorRecord();
	/// @brief constructor
	MutableICoorRecord(
		Real const phi_in,
		Real const theta_in,
		Real const d_in,
		std::string const & stub_atom1_name,
		std::string const & stub_atom2_name,
		std::string const & stub_atom3_name
	);

	bool
	operator==( MutableICoorRecord const & rhs ) const;

public:
	/// @brief accessor to stub_atom1 ICoorAtomID
	Real
	phi() const
	{
		return phi_;
	}


	Real
	theta() const
	{
		return theta_;
	}


	Real
	d() const
	{
		return d_;
	}


	std::string const &
	stub_atom1() const
	{
		return stub_atom1_;
	}

	std::string const &
	stub_atom2() const
	{
		return stub_atom2_;
	}

	std::string const &
	stub_atom3() const
	{
		return stub_atom3_;
	}

	/// @brief accessor to stub_atom ICoorAtomID
	std::string const &
	stub_atom( int const atm ) const
	{
		switch( atm ) {
		case 1 : return stub_atom1_;
		case 2 : return stub_atom2_;
		case 3 : return stub_atom3_;
		}
		utility_exit_with_message( "MutableICoorRecord::stub_atom(): passed value should be 1--3" );
		return stub_atom1_;
	}

	ICoordAtomIDType
	stub_type1() const
	{
		return stub_type1_;
	}

	ICoordAtomIDType
	stub_type2() const
	{
		return stub_type2_;
	}

	ICoordAtomIDType
	stub_type3() const
	{
		return stub_type3_;
	}

	/// @brief accessor to stub_type ICoorAtomID
	ICoordAtomIDType
	stub_type( int const atm ) const
	{
		switch( atm ) {
		case 1 : return stub_type1_;
		case 2 : return stub_type2_;
		case 3 : return stub_type3_;
		}
		utility_exit_with_message( "MutableICoorRecord::stub_type(): passed value should be 1--3" );
		return stub_type1_;
	}

	void show( std::ostream & out ) const;

	/// @brief Can valid coordinates be built for this MutableICoorRecord, given the residue type?
	bool
	buildable( MutableResidueType const & rsd_type, bool verbose = false ) const;

public:

	/// @brief Build the location of the built atom, given the other atoms in the residue type.
	Vector
	build( MutableResidueType const & rsd_type ) const;

	/// @brief Given a the stub designation (1/2/3) find the coordinate from the MutableResidueType
	Vector
	xyz( core::Size stubno, MutableResidueType const & restype ) const;

	/// @brief Given a stub designation (atom name, UPPER/LOWER/CONN, etc.) find the coordinate from the MutableResidueType
	static
	Vector
	build_xyz( std::string const & stub, MutableResidueType const & restype );

private:

	Real phi_;
	Real theta_;
	Real d_;
	std::string stub_atom1_;
	std::string stub_atom2_;
	std::string stub_atom3_;
	ICoordAtomIDType stub_type1_;
	ICoordAtomIDType stub_type2_;
	ICoordAtomIDType stub_type3_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

inline
std::ostream &
operator<<( std::ostream & out, MutableICoorRecord const & icoor ) {
	icoor.show( out );
	return out;
}

} // chemical
} // core


#endif
