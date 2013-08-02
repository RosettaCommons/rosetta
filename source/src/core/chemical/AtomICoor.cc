// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief declaration of implementation class for abstract class Residue
/// @author Phil Bradley
/// @author Modified by Sergey Lyskov


// Unit headers
#include <core/chemical/AtomICoor.hh>

// Project headers
#include <core/chemical/ResidueConnection.hh>
#include <core/conformation/Conformation.hh>
#include <basic/Tracer.hh>

// Numeric headers
// Commented by inclean daemon #include <numeric/xyz.functions.hh>

// Utility headers
// Commented by inclean daemon #include <utility/exit.hh>

// AUTO-REMOVED #include <ObjexxFCL/ObjexxFCL.hh>
#include <ObjexxFCL/string.functions.hh>

#include <utility/vector1.hh>
#include <numeric/xyz.functions.hh>

//Auto using namespaces
namespace ObjexxFCL { } using namespace ObjexxFCL; // AUTO USING NS
//Auto using namespaces end


// C++ headers
// Commented by inclean daemon #include <string>


namespace core {
namespace chemical {

static basic::Tracer tw( "core.chemical.AtomICoor", basic::t_warning );

//////////////////////////////////////////////////////////////////////////////////////////////////
/**
	 After atom name is read in from residue param file, ICoorAtomID type_ and atomno_ is defined as:

	 - Anything less than four character is considered as INTERNAL atom name and its atom index
	 number is assigned as atomno_;

	 - "LOWER" and "UPPER" for polymer lower and upper connections. Since they are unique, no atomno_
	 is given ( e.g., 0)

	 - Non-polymer connections are flagged by "CONN*" in which * represents the index number of
	 this connection in the ResidueType ( from 1 to ResidueType.n_connection() ). This number is assigned
	 as atomno_.
 */
ICoorAtomID::ICoorAtomID(
		std::string name,
		ResidueType const & rsd_type
)
{
	if ( name.size() <= 4 ) {
		type_ = INTERNAL;
		atomno_ = rsd_type.atom_index( name );
	} else if ( name.substr(0,4) == "CONN" ) {
		type_ = CONNECT;
		assert( is_int( name.substr(4) ) );
		atomno_ = int_of( name.substr(4) );
		assert( atomno_ > 0 && atomno_ <= rsd_type.n_residue_connections() );
	} else if ( name == "LOWER" ) {
		type_ = POLYMER_LOWER; atomno_ = 0;
		// atomno is unused
	} else if ( name == "UPPER" ) {
		type_ = POLYMER_UPPER;
		atomno_ = 0; // atomno is unused
	} else {
		utility_exit_with_message( "ICoorAtomID: unable to parse atom_name: "+name );
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////
Vector const &
ICoorAtomID::xyz(
		Residue const & rsd,
		Conformation const & conformation
) const
{
	static Vector NullVector( 0, 0, 0 );

	//std::cout << "ICoorAtomID::xyz " << rsd.name() << ' ' << type_ << ' ' << atomno_ << std::endl;
	if ( type_ == INTERNAL ) {
		//std::cout << "ICoorAtomID::xyz " << rsd.name() << ' ' << atomno_ << ' ' << rsd.atom_is_backbone(atomno_) <<
		//	' ' << rsd.atom_name( atomno_ ) << ' ' <<
		//	rsd.atom( atomno_ ).xyz()(1) << ' ' <<
		//	rsd.atom( atomno_ ).xyz()(2) << ' ' <<
		//	rsd.atom( atomno_ ).xyz()(3) << std::endl;
		return rsd.atom( atomno_ ).xyz();
	} else if ( type_ == POLYMER_LOWER ) {
		int const seqpos( rsd.seqpos() - 1 );
		int const atomno( conformation.residue_type( seqpos ).upper_connect_atom() );
		return conformation.xyz( id::AtomID( atomno, seqpos ) );
		// 		Residue const & rsd_lower( conformation.residue( rsd.seqpos() - 1 ) );
		// 		return rsd_lower.atom( rsd_lower.upper_connect_atom() ).xyz();
	} else if ( type_ == POLYMER_UPPER ) {
		int const seqpos( rsd.seqpos() + 1 );
		int const atomno( conformation.residue_type( seqpos ).lower_connect_atom() );
		return conformation.xyz( id::AtomID( atomno, seqpos ) );
		// 		Residue const & rsd_upper( conformation.residue( rsd.seqpos() + 1 ) );
		// 		return rsd_upper.atom( rsd_upper.lower_connect_atom() ).xyz();
	} else if ( type_ == CONNECT ) {
		// in this case atomno_ is the connection identifer (index)
		Size const connid( atomno_ );
		int const  partner_seqpos( rsd.residue_connection_partner( connid ) );
		if ( partner_seqpos < 1 || partner_seqpos > int( conformation.size() ) ) {
			tw << "ICoorAtomID xyz depends on invalid residue connection, returning BOGUS coords: this_rsd= " << rsd.name() <<
					' ' << rsd.seqpos() << " connid= " << connid << " partner_seqpos= " << partner_seqpos << '\n';
			return NullVector;
		}
		Size const partner_connid( rsd.residue_connection_conn_id( connid ) );
		Size const partner_atomno( conformation.residue_type( partner_seqpos ).residue_connect_atom_index( partner_connid));
		return conformation.xyz( id::AtomID( partner_atomno, partner_seqpos ) );
	} else {
		utility_exit_with_message( "unrecognized stub atom id type!" );
	}

	// to appease the compiler
	return NullVector;
}


//////////////////////////////////////////////////////////////////////////////////////////////////
Vector // const &
ICoorAtomID::xyz(
		ResidueType const & rsd_type
) const
{
	if ( type_ == INTERNAL ) {
		return rsd_type.atom( atomno_ ).ideal_xyz();
	} else if ( type_ == POLYMER_LOWER ) {
		return rsd_type.lower_connect().icoor().build( rsd_type );
	} else if ( type_ == POLYMER_UPPER ) {
		return rsd_type.upper_connect().icoor().build( rsd_type );
	} else if ( type_ == CONNECT ) {
		return rsd_type.residue_connection( atomno_ ).icoor().build( rsd_type );
	} else {
		utility_exit_with_message( "unrecognized stub atom id type!" );
	}

	// to appease the compiler
	static Vector NullVector( 0, 0, 0 );
	return NullVector;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief WARNING: Slightly dangerous function intended for black magic use only.
///    Only to be used for situations where you *know* the ICoorAtomID can't be anything but
///    a real atom on the given residue, and where a conformation is absolutely not availible.
///    If you /can/ use ICoorAtomID::xyz( Residue const &, Conformation const &), you /should/.
Vector
ICoorAtomID::xyz( conformation::Residue const & rsd ) const
{
	if ( type_ == INTERNAL ) {
		return rsd.atom( atomno_ ).xyz();
	} else if ( type_ == POLYMER_LOWER ) {
		return rsd.type().lower_connect().icoor().build( rsd );
	} else if ( type_ == POLYMER_UPPER ) {
		return rsd.type().upper_connect().icoor().build( rsd );
	} else if ( type_ == CONNECT ) {
		return rsd.type().residue_connection( atomno_ ).icoor().build( rsd );
		utility_exit_with_message( "unrecognized stub atom id type!" );
	}
}
	

//////////////////////////////////////////////////////////////////////////////////////////////////
id::AtomID
ICoorAtomID::atom_id( Size const seqpos, Conformation const & conformation ) const
{
	using id::AtomID;
	switch ( type_ ) {
	case INTERNAL:
		return AtomID( atomno_, seqpos );
	case POLYMER_LOWER:
		return AtomID( conformation.residue_type( seqpos-1 ).upper_connect_atom(), seqpos - 1 );
	case POLYMER_UPPER:
		return AtomID( conformation.residue_type( seqpos+1 ).lower_connect_atom(), seqpos + 1 );
	case CONNECT:
		// in this case atomno_ is the connection identifer (index)
		return conformation.inter_residue_connection_partner( seqpos, atomno_ );
	default:
		return id::BOGUS_ATOM_ID;
	}
	return id::BOGUS_ATOM_ID;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
Vector
AtomICoor::build(
		conformation::Residue const & rsd,
		conformation::Conformation const & conformation
) const
{
	assert( kinematics::Stub( stub_atom1_.xyz( rsd, conformation ),
			stub_atom2_.xyz( rsd, conformation ),
			stub_atom3_.xyz( rsd, conformation ) ).is_orthogonal( 0.001 ) );

	return kinematics::Stub( stub_atom1_.xyz( rsd, conformation ),
			stub_atom2_.xyz( rsd, conformation ),
			stub_atom3_.xyz( rsd, conformation ) ).spherical( phi_, theta_, d_ );
}

Vector
AtomICoor::build(
		ResidueType const & rsd_type
) const
{
	assert( kinematics::Stub( stub_atom1_.xyz( rsd_type ),
			stub_atom2_.xyz( rsd_type ),
			stub_atom3_.xyz( rsd_type ) ).is_orthogonal( 0.001 ) );

	return kinematics::Stub( stub_atom1_.xyz( rsd_type ),
			stub_atom2_.xyz( rsd_type ),
			stub_atom3_.xyz( rsd_type ) ).spherical( phi_, theta_, d_ );
}


/// @brief WARNING: Slightly dangerous function intended for black magic use only.
///    Only to be used for situations where you *know* the AtomICoor /and all it's stub atoms/ can't be
///    anything but real atoms on the given residue, and where a conformation is absolutely not availible.
///    If you /can/ use AtomICoor::build( Residue const &, Conformation const &), you /should/.
Vector
AtomICoor::build( 
		conformation::Residue const & rsd 
) const
{
	kinematics::Stub built_stub( stub_atom1_.xyz( rsd ),
		stub_atom2_.xyz( rsd ),
		stub_atom3_.xyz( rsd ));
			
	assert( built_stub.is_orthogonal( 0.001 ) );
	
	return built_stub.spherical( phi_, theta_, d_ );
}

} // chemical
} // core
