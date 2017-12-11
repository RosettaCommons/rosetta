// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/id/AtomID.hh
/// @author Phil Bradley


#ifndef INCLUDED_core_id_AtomID_hh
#define INCLUDED_core_id_AtomID_hh

// Unit headers
#include <core/id/AtomID.fwd.hh>

#include <utility/exit.hh>
// C++ headers

#include <core/types.hh>


namespace core {
namespace id {


/////////////////////////////////////////////////////////////////////////////
// AtomID
/////////////////////////////////////////////////////////////////////////////

// data to uniquely specify an atom
// could change -- pointer?

/// @brief  Atom identifier class.  Defined by the atom number and the residue number.
class AtomID
{

public: // Creation

	/// @brief Default constructor
	inline
	constexpr
	AtomID() :
		atomno_( 0 ),
		rsd_( 0 )
	{}

	/// @brief Copy constructor
	inline
	AtomID( AtomID const & ) = default;

	/// @brief Property constructor
	inline
	AtomID(
		Size const atomno_in,
		Size const rsd_in
	) :
		atomno_( atomno_in ),
		rsd_( rsd_in )
	{}

	static constexpr AtomID BOGUS_ATOM_ID() { return AtomID(); }
	static constexpr AtomID CHAINBREAK_BOGUS_ATOM_ID() { return AtomID(); }

	///@brief Set the value of atom and residue.
	void
	set( core::Size atomno_in, core::Size rsd_in);

public: // Properties

	/// @brief Returns the AtomID residue number
	inline
	Size
	rsd() const { return rsd_; }

	inline
	Size &
	rsd() { return rsd_; }

	/// @brief Returns the AtomID atom number
	inline
	Size
	atomno() const { return atomno_; }

	inline
	Size &
	atomno() { return atomno_; }

	/// @brief Returns true if the AtomID is valid
	/// @note must return false for BOGUS_ATOM_ID
	inline
	bool
	valid() const { return ( atomno_ > 0 ) && ( rsd_ > 0 ); }


public: // Friends

	friend
	std::ostream &
	operator <<(
		std::ostream & os,
		AtomID const & a
	);

	/// @brief a and b are the same atom
	friend
	inline
	bool
	operator ==(
		AtomID const & a,
		AtomID const & b
	) { return a.atomno_ == b.atomno_ && a.rsd_ == b.rsd_; }

	/// @brief a and b are different atom
	friend
	inline
	bool
	operator !=(
		AtomID const & a,
		AtomID const & b
	) { return a.atomno_ != b.atomno_ || a.rsd_ != b.rsd_; }

	/// @brief a is LOWER than b (e.g., first by smaller residue index number then by smaller atom index number)
	friend
	inline
	bool
	operator <(
		AtomID const & a,
		AtomID const & b
	)
	{
		return ( a.rsd_ <  b.rsd_ ||
			( a.rsd_ == b.rsd_ && a.atomno_ < b.atomno_ ) );
	}

private: // Fields


	/// @brief Atom number within the Residue
	Size atomno_;

	/// @brief Residue number within the complex
	Size rsd_;


#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; // AtomID

/// Globals -- may not be used until after core::init is called.
extern AtomID const GLOBAL_BOGUS_ATOM_ID;
extern AtomID const GLOBAL_CHAINBREAK_BOGUS_ATOM_ID;

///////////////////////////////////////////////////////////////////////////////
/// two more classes, temporary for testing purposes
///////////////////////////////////////////////////////////////////////////////

class BondID {
public:

	BondID():
		atom1( AtomID::BOGUS_ATOM_ID() ),
		atom2( AtomID::BOGUS_ATOM_ID() )
	{}

	BondID( AtomID const & a1, AtomID const & a2 ):
		atom1( a1 ),
		atom2( a2 )
	{}

	bool
	has( AtomID const & id ) const
	{
		return ( id == atom1 || id == atom2 );
	}

	AtomID const &
	other_atom( AtomID const & id ) const;

	BondID
	reversed() const
	{
		return BondID( atom2, atom1 );
	}

	void
	reverse()
	{
		AtomID const tmp( atom1 );
		atom1 = atom2;
		atom2 = tmp;
	}

	friend
	bool
	operator==( BondID const & a, BondID const & b )
	{
		return ( a.atom1 == b.atom1 && a.atom2 == b.atom2 );
	}


public:
	AtomID atom1;
	AtomID atom2;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


///////////////////////////////////////////////////////////////////////////////

class StubID {
public:
	StubID( AtomID const & a1, AtomID const & a2, AtomID const & a3 ):
		atom1( a1 ),
		atom2( a2 ),
		atom3( a3 ),
		center_()
	{}

	StubID( AtomID const &c, AtomID const & a1, AtomID const & a2, AtomID const & a3 ):
		atom1( a1 ),
		atom2( a2 ),
		atom3( a3 ),
		center_( c )
	{}

	constexpr
	StubID():
		atom1(AtomID::BOGUS_ATOM_ID()),
		atom2(AtomID::BOGUS_ATOM_ID()),
		atom3(AtomID::BOGUS_ATOM_ID())
	{}

	static constexpr StubID const BOGUS_STUB_ID() { return StubID(); }

	AtomID const &
	atom( Size const index ) const {
		switch ( index ) {
		case 1 :
			return atom1;
		case 2 :
			return atom2;
		case 3 :
			return atom3;
		default :
			utility_exit_with_message("StubID's have exactly three atoms, 1-3");
		}
		return atom1; // won't get here
	}

	AtomID const &
	center() const {
		return center_;
	}

	bool valid() const {
		return atom1.valid() && atom2.valid() && atom3.valid() && center_.valid();
	}

	inline
	friend
	bool
	operator< ( StubID const & a, StubID const & b )
	{
		return ( ( a.atom1  < b.atom1 ) ||
			( a.atom1 == b.atom1 && a.atom2  < b.atom2 ) ||
			( a.atom1 == b.atom1 && a.atom2 == b.atom2 && a.atom3 < b.atom3 ) );
	}

	friend
	std::ostream &
	operator <<(
		std::ostream & os,
		StubID const &
	);

	bool operator == ( StubID const & rhs ) const;
	bool operator != ( StubID const & rhs ) const;

public: // tmp hack -- phil fix this
	AtomID atom1;
	AtomID atom2;
	AtomID atom3;
	AtomID center_;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


/// @brief Global -- may not be used until after core::init is called.
extern StubID const GLOBAL_BOGUS_STUB_ID;


} // namespace id
} // namespace core


#endif // INCLUDED_core_id_AtomID_HH
