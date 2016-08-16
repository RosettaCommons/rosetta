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

#include <devel/cartesian_frags/SafeID.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/conformation/Conformation.hh>
#include <core/kinematics/Stub.hh>

#include <ObjexxFCL/string.functions.hh>

#include <utility/vector1.hh>

//Auto Headers
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/FoldTree.hh>


namespace devel {
namespace cartesian_frags {

using namespace core;

///////////////////////////////////////////////////////////////////////////////
/// must be a better place for this, probably already exists!

std::string
strip_whitespace( std::string const & name )
{
	using namespace ObjexxFCL;
	std::string trimmed_name( name );
	left_justify( trimmed_name ); trim( trimmed_name ); // simpler way to do this?
	return trimmed_name;
}


///////////////////////////////////////////////////////////////////////////////

SafeAtomID::SafeAtomID( std::string const & name_in, int const seqpos_in ):
	name( strip_whitespace( name_in ) ),
	seqpos( seqpos_in )
{}


///////////////////////////////////////////////////////////////////////////////
/// SUBTRACT SEQPOS_OFFSET
SafeAtomID::SafeAtomID( id::AtomID const & atom_id, conformation::Conformation const & conf, int const seqpos_offset ):
	name( strip_whitespace( conf.residue_type( atom_id.rsd() ).atom_name( atom_id.atomno() ) ) ),
	seqpos( atom_id.rsd() - seqpos_offset )
{}


///////////////////////////////////////////////////////////////////////////////
/// ADD SEQPOS_OFFSET
id::AtomID
SafeAtomID::id( conformation::Conformation const & conf, int const seqpos_offset ) const
{
	return id::AtomID( conf.residue_type( seqpos + seqpos_offset ).atom_index( name ), seqpos + seqpos_offset );
}

///////////////////////////////////////////////////////////////////////////////
//
std::ostream &
operator << ( std::ostream & os, SafeAtomID const & a )
{
	os << a.seqpos << ' ' << a.name;
	return os;
}


///////////////////////////////////////////////////////////////////////////////
id::StubID
stub_id_from_torsion_stub_id(
	TorsionStubID const & id,
	conformation::Conformation const & conf
)
{
	using namespace id;
	AtomID atom1, atom2, atom3, atom4;
	bool const fail( conf.get_torsion_angle_atom_ids( id.id, atom1, atom2, atom3, atom4 ) );
	if ( !fail ) {
		if ( id.dir == Backward ) {
			return StubID( atom2, atom3, atom4 );
		} else {
			assert( id.dir == Forward );
			return StubID( atom3, atom2, atom1 );
		}
	} else if ( id.id.type() == BB ) {
		// most likely spanning a cutpoint
		Size const seqpos ( id.id.rsd() );
		Size const torsion( id.id.torsion() );
		conformation::Residue const & rsd( conf.residue( seqpos ) );
		chemical::AtomIndices const & mainchain( rsd.mainchain_atoms() );
		Size const nbb( mainchain.size() );
		if ( torsion == 1 && conf.fold_tree().is_cutpoint( seqpos-1 ) && id.dir == Backward ) {
			return StubID( AtomID( mainchain[1], seqpos ), AtomID( mainchain[2], seqpos ), AtomID( mainchain[3], seqpos ));
		} else if ( torsion == nbb-1 && conf.fold_tree().is_cutpoint( seqpos ) && id.dir == Forward ) {
			return StubID( AtomID( mainchain[nbb  ], seqpos ),
				AtomID( mainchain[nbb-1], seqpos ),
				AtomID( mainchain[nbb-2], seqpos ));
		}
	}

	utility_exit_with_message( "stub_id_from_torsion_stub_id failed!");

	return StubID( atom1, atom2, atom3 ); // wont get here
}

///////////////////////////////////////////////////////////////////////////////
kinematics::Stub
torsion_stub(
	id::TorsionID const & id,
	Direction const & dir,
	conformation::Conformation const & conf
)
{
	return conf.atom_tree().stub_from_id( stub_id_from_torsion_stub_id( TorsionStubID( id, dir ), conf ) );
}


//  // the coord sys you'd get by building atom1 then atom2 then atom3, BUT centered at atom2 not atom3
//  //  NO LONGER CENTERED AT ATOM2
//  // ie: origin at atom2 -- no changed to atom3
//  //     x-axis from atom2->atom3
//  //     y-axis in atom1,atom2,atom3 plane w/ pos dotprod with atom2->atom1 vector
//  //
//  return kinematics::Stub( conf.xyz( atom3 ), conf.xyz( atom2 ), conf.xyz( atom1 ) );
// //  return kinematics::Stub( conf.xyz( atom2 ), conf.xyz( atom3 ), conf.xyz( atom2 ), conf.xyz( atom1 ) );


///////////////////////////////////////////////////////////////////////////////
kinematics::Stub
torsion_stub(
	TorsionStubID const & id,
	conformation::Conformation const & conf
)
{
	return torsion_stub( id.id, id.dir, conf );
}


///////////////////////////////////////////////////////////////////////////////
/// orientation of atoms such that atom1 is the first stub atom, atom2 the second


id::BondID
torsion_bond( TorsionStubID const & id, conformation::Conformation const & conf )
{
	id::StubID const stubid( stub_id_from_torsion_stub_id( id, conf ) );
	return id::BondID( stubid.atom1, stubid.atom2 );
}


std::ostream &
operator << ( std::ostream & os, SafeStubID const & a )
{
	os << '(' << a.atom1 << ' ' << a.atom2 << ' ' << a.atom3 << ')';
	return os;
}


} // cartesian_frags
} // devel
