// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/id/NamedStubID.hh
/// @author Phil Bradley


#ifndef INCLUDED_core_id_NamedStubID_hh
#define INCLUDED_core_id_NamedStubID_hh


// Unit headers
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/NamedAtomID.hh>
#include <core/types.hh>

// Utility Headers
#include <utility/vector1.hh>

// C++ headers


namespace core {
namespace id {


/////////////////////////////////////////////////////////////////////////////
// NamedStubID
/////////////////////////////////////////////////////////////////////////////

// data to uniquely specify an atom
// could change -- pointer?

///////////////////////////////////////////////////////////////////////////////

class NamedStubID {
public:
	typedef utility::vector1<std::string> AtomList;

	NamedStubID( NamedAtomID const & a1, NamedAtomID const & a2, NamedAtomID const & a3 ):
		center_(),
		atom1( a1 ),
		atom2( a2 ),
		atom3( a3 )
	{}

	NamedStubID( NamedAtomID const & c, NamedAtomID const & a1, NamedAtomID const & a2, NamedAtomID const & a3 ) :
		center_( c ),
		atom1( a1 ),
		atom2( a2 ),
		atom3( a3 )
	{}

	// convienience c'stor if the residue is the same for all atoms
	NamedStubID( std::string const& a1, std::string const& a2, std::string const& a3, core::Size rsd );

	// convienience c'stor if the residue is the same for all atoms
	NamedStubID( std::string const& c, std::string const& a1, std::string const& a2, std::string const& a3, core::Size rsd );

	// convienience c'stor
	NamedStubID( std::string const& a1, Size rsd1, std::string const& a2, Size rsd2, std::string a3, core::Size rsd3 );

	// convienience c'stor takes list of strings either size 3 or size 4. If size 4 first atom is center
	NamedStubID( AtomList const&, core::Size rsd );

	NamedStubID() :
		center_(),
		atom1(),
		atom2(),
		atom3()
	{}

	NamedAtomID const &
	atom( Size const index ) const;

	NamedAtomID const &
	center() const {
		return center_;
	}

	bool valid() const {
		return atom1.valid() && atom2.valid() && atom3.valid();
	}

	inline
	friend
	bool
	operator< ( NamedStubID const & a, NamedStubID const & b )
	{
		return ( ( a.atom1  < b.atom1 ) ||
			( a.atom1 == b.atom1 && a.atom2  < b.atom2 ) ||
			( a.atom1 == b.atom1 && a.atom2 == b.atom2 && a.atom3 < b.atom3 ) );
	}


	/// @brief input operator
	friend std::istream & operator >>(std::istream & is, NamedStubID& e);

	/// @brief output operator
	friend std::ostream & operator <<(std::ostream & os, NamedStubID const& e);


public: // tmp hack -- phil fix this
	NamedAtomID center_;
	NamedAtomID atom1;
	NamedAtomID atom2;
	NamedAtomID atom3;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} // namespace id
} // namespace core


#endif // INCLUDED_core_id_AtomID_HH
