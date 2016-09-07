// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief

#ifndef INCLUDED_devel_cartesian_frags_SafeID_hh
#define INCLUDED_devel_cartesian_frags_SafeID_hh

#include <devel/cartesian_frags/SafeID.fwd.hh>
#include <devel/cartesian_frags/Direction.hh>
#include <core/id/AtomID.hh>
#include <core/id/TorsionID.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>

#include <utility>
#include <utility/exit.hh>

#include <string>
#include <iosfwd>


namespace devel {
namespace cartesian_frags {


/// @brief  An AtomID which is "safe" in that it's defined with a string-atom-name rather than an index

class SafeAtomID {
public:
	///
	SafeAtomID():
		name(),
		seqpos(0)
	{}


	SafeAtomID( std::string const & name_in, int const seqpos_in );

	/// SUBTRACT SEQPOS_OFFSET
	SafeAtomID( core::id::AtomID const & atom_id, core::conformation::Conformation const & conf, int const seqpos_offset=0 );

	/// ADD SEQPOS_OFFSET
	core::id::AtomID
	id( core::conformation::Conformation const & conf, int const seqpos_offset = 0 ) const;


	friend
	bool
	operator< ( SafeAtomID const & a, SafeAtomID const & b )
	{
		return ( a.seqpos < b.seqpos || ( a.seqpos == b.seqpos && a.name < b.name ) );
	}


	friend
	bool
	operator== ( SafeAtomID const & a, SafeAtomID const & b )
	{
		return ( a.seqpos == b.seqpos && a.name == b.name );
	}


	friend
	std::ostream &
	operator << ( std::ostream & os, SafeAtomID const & a );


	/// tmp public
	std::string name;
	int seqpos;

};


///////////////////////////////////////////////////////////////////////////////
/// @brief  A "safe: BondID, analogous to SafeAtomID

class SafeBondID {
public:

	SafeBondID( SafeAtomID  a1, SafeAtomID  a2 ):
		atom1(std::move( a1 )),
		atom2(std::move( a2 ))
	{}


	core::id::BondID
	id( core::conformation::Conformation const & conf ) const
	{
		return core::id::BondID( atom1.id( conf ), atom2.id( conf ) );
	}

	// tmp hack
	SafeAtomID atom1;
	SafeAtomID atom2;


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief  An ID which defines a stub along a torsion bond. Uses TorsionID plus an additional direction.

class TorsionStubID {
public:

	TorsionStubID():
		id (),
		dir( Middle ) // nonsense
	{}

	TorsionStubID( core::id::TorsionID const & id_in, Direction const & dir_in ):
		id ( id_in ),
		dir( dir_in )
	{}

	TorsionStubID( TorsionStubID const & ) = default;

	TorsionStubID &
	operator=( TorsionStubID const & ) = default;


	/// tmp public
	core::id::TorsionID id;
	Direction dir;

};


core::id::StubID
stub_id_from_torsion_stub_id(
	TorsionStubID const & id,
	core::conformation::Conformation const & conf
);


core::kinematics::Stub
torsion_stub(
	core::id::TorsionID const & id,
	Direction const & dir,
	core::conformation::Conformation const & conf
);


core::kinematics::Stub
torsion_stub( TorsionStubID const & id, core::conformation::Conformation const & conf );

/// @brief orientation of atoms such that atom1 is the first stub atom, atom2 the second
core::id::BondID
torsion_bond( TorsionStubID const & id, core::conformation::Conformation const & conf );


///////////////////////////////////////////////////////////////////////////////
/// @brief  A "safe" stub ID.

class SafeStubID {

public:


	// ctor from torsionstubid, conf, seq-offset?
	SafeStubID():
		atom1(),
		atom2(),
		atom3()
	{}

	// ctor from torsionstubid, conf, seq-offset?
	SafeStubID(
		TorsionStubID const & tor,
		core::conformation::Conformation const & conf,
		int const seqpos_offset = 0
	)
	{
		// fix this code duplication
		core::id::StubID const stubid( stub_id_from_torsion_stub_id( tor, conf ) );
		atom1 = SafeAtomID( stubid.atom1, conf, seqpos_offset );
		atom2 = SafeAtomID( stubid.atom2, conf, seqpos_offset );
		atom3 = SafeAtomID( stubid.atom3, conf, seqpos_offset );
	}

	core::id::StubID
	id( core::conformation::Conformation const & conf, int const seqpos_offset = 0 ) const
	{
		return core::id::StubID( atom1.id( conf, seqpos_offset ),
			atom2.id( conf, seqpos_offset ),
			atom3.id( conf, seqpos_offset ) );
	}


	SafeAtomID const &
	atom( int const i ) const
	{
		switch( i ) {
		case 1 : return atom1;
		case 2 : return atom2;
		case 3 : return atom3;
		default :
			utility_exit_with_message( "SafeStubAtom::atom: bad atom!" );
		}
		return atom1; // wont get here
	}


	friend
	bool
	operator==( SafeStubID const & a, SafeStubID const & b )
	{
		return ( ( a.atom1 == b.atom1 ) && ( a.atom2 == b.atom2 ) && ( a.atom3 == b.atom3 ) );
	}


	friend
	bool
	operator!=( SafeStubID const & a, SafeStubID const & b )
	{
		return ( !( a == b ) );
	}


	friend
	std::ostream &
	operator << ( std::ostream & os, SafeStubID const & a );

	//tmp hack
	SafeAtomID atom1;
	SafeAtomID atom2;
	SafeAtomID atom3;
};


} // cartesian_frags
} // devel

#endif
