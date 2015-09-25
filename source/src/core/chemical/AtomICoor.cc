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

#include <basic/datacache/DataMapObj.hh> // To wrap std::set for prettyprinter

// Numeric headers

// Utility headers

#include <ObjexxFCL/string.functions.hh>

#include <utility/vector1.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/NumericTraits.hh>

//Auto using namespaces
namespace ObjexxFCL { } using namespace ObjexxFCL; // AUTO USING NS
//Auto using namespaces end


// C++ headers


namespace core {
namespace chemical {

static THREAD_LOCAL basic::Tracer tw( "core.chemical.AtomICoor", basic::t_warning );

ICoorAtomID::ICoorAtomID():
	type_( INTERNAL ),
	atomno_( 0 ),
	vd_(ResidueType::null_vertex)
{}

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
		vd_ = rsd_type.atom_vertex(atomno_);
	} else if ( name.substr(0,4) == "CONN" ) {
		type_ = CONNECT;
		debug_assert( is_int( name.substr(4) ) );
		atomno_ = int_of( name.substr(4) );
		vd_ = ResidueType::null_vertex;
		debug_assert( atomno_ > 0 && atomno_ <= rsd_type.n_residue_connections() );
		if ( atomno_ > rsd_type.n_residue_connections() ) { // using > is not a great check but it is better than !=
			tw.Warning << "The record for CONN" << atomno_ << " in the topology file for " << rsd_type.name() <<
				" either has an incorrect index or is listed out of order in the file." << std::endl;
		}
	} else if ( name == "LOWER" ) {
		type_ = POLYMER_LOWER; atomno_ = 0;
		atomno_ = 0; // atomno is unused
		vd_ = ResidueType::null_vertex;
		// atomno is unused
	} else if ( name == "UPPER" ) {
		type_ = POLYMER_UPPER;
		atomno_ = 0; // atomno is unused
		vd_ = ResidueType::null_vertex;
	} else {
		utility_exit_with_message( "ICoorAtomID: unable to parse atom_name: "+name );
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////
/**
A Icoord to an atom given via a vertex descriptor is by definition INTERNAL.
*/

ICoorAtomID::ICoorAtomID(
	VD vd,
	ResidueType const & rsd_type
)
{
	debug_assert( rsd_type.has(vd) );
	type_ = INTERNAL;
	atomno_ = rsd_type.atom_index( vd );
	vd_ = vd;
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
		//' ' << rsd.atom_name( atomno_ ) << ' ' <<
		//rsd.atom( atomno_ ).xyz()(1) << ' ' <<
		//rsd.atom( atomno_ ).xyz()(2) << ' ' <<
		//rsd.atom( atomno_ ).xyz()(3) << std::endl;
		return rsd.atom( atomno_ ).xyz();
	} else if ( type_ == POLYMER_LOWER ) {
		//Changed by VKM on 9 April 2014: we no longer assume that the residue connected at the POLYMER_LOWER connection is the i-1 residue in sequence; nor
		//do we assume that that residue is connected via its upper terminus.
		debug_assert(!rsd.is_lower_terminus());
		core::Size const lowerID = rsd.type().lower_connect_id(); //The index of the lower connection of THIS residue.
		//debug_assert(!rsd.connection_incomplete( lowerID ));
		core::Size const seqpos( rsd.connect_map( lowerID ).resid() ); //Get the index of the residue connected to this one at the lower connection.
		if ( seqpos == 0 ) {
			tw << "Warning from IcoorAtomID::xyz(): Cannot get xyz for POLYMER_LOWER of residue " << rsd.name() << " " << rsd.seqpos() << ".  Returning BOGUS coords (null vector)." << std::endl;
			return NullVector;
		}
		core::Size const atomno( conformation.residue_type(seqpos).residue_connect_atom_index( rsd.connect_map( lowerID ).connid()  ) ); //Get the index of the atom on the other residue that forms a connection to the lower terminus of this residue.
		//int const seqpos( rsd.seqpos() - 1 );
		//int const atomno( conformation.residue_type( seqpos ).upper_connect_atom() );
		return conformation.xyz( id::AtomID( atomno, seqpos ) );
		//   Residue const & rsd_lower( conformation.residue( rsd.seqpos() - 1 ) );
		//   return rsd_lower.atom( rsd_lower.upper_connect_atom() ).xyz();
	} else if ( type_ == POLYMER_UPPER ) {
		//Changed by VKM on 9 April 2014: we no longer assume that the residue connected at the POLYMER_UPPER connection is the i+1 residue in sequence; nor
		//do we assume that that residue is connected via its lower terminus.
		debug_assert(!rsd.is_upper_terminus());
		core::Size const upperID = rsd.type().upper_connect_id(); //The index of the upper connection of THIS residue.
		//debug_assert(!rsd.connection_incomplete( upperID ));
		core::Size const seqpos( rsd.connect_map( upperID ).resid() ); //Get the index of the residue connected to this one at the upper connection.
		if ( seqpos == 0 ) {
			tw << "Warning from IcoorAtomID::xyz(): Cannot get xyz for POLYMER_UPPER of residue " << rsd.name() << " " << rsd.seqpos() << ".  Returning BOGUS coords (null vector)." << std::endl;
			return NullVector;
		}
		core::Size const atomno( conformation.residue_type(seqpos).residue_connect_atom_index( rsd.connect_map( upperID ).connid()  ) ); //Get the index of the atom on the other residue that forms a connection to the upper terminus of this residue.
		//int const seqpos( rsd.seqpos() + 1 );
		//int const atomno( conformation.residue_type( seqpos ).lower_connect_atom() );
		return conformation.xyz( id::AtomID( atomno, seqpos ) );
		//   Residue const & rsd_upper( conformation.residue( rsd.seqpos() + 1 ) );
		//   return rsd_upper.atom( rsd_upper.lower_connect_atom() ).xyz();
	} else if ( type_ == CONNECT ) {
		// in this case atomno_ is the connection identifer (index)
		Size connid( atomno_ );
		if ( rsd.is_lower_terminus() ) --connid;
		if ( rsd.is_upper_terminus() ) --connid;
		//tw << "conid=" << connid << std::endl; tw.flush(); //DELETE ME
		//debug_assert(!rsd.connection_incomplete( connid ) );
		int const  partner_seqpos( rsd.residue_connection_partner( connid ) );
		if ( partner_seqpos < 1 || partner_seqpos > int( conformation.size() ) ) {
			tw << "Warning from IcoorAtomID::xyz(): ICoorAtomID xyz depends on invalid residue connection, returning BOGUS coords (null vector): this_rsd= " << rsd.name() <<
				' ' << rsd.seqpos() << " connid= " << connid << " partner_seqpos= " << partner_seqpos << std::endl;
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
		return rsd_type.atom( vd_ ).ideal_xyz();
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
	} else {
		utility_exit_with_message( "unrecognized stub atom id type!" );
	}
	return Vector( 0.0 ); //get rid of compiler warnings.
}


//////////////////////////////////////////////////////////////////////////////////////////////////
id::AtomID
ICoorAtomID::atom_id( Size const seqpos, Conformation const & conformation ) const
{
	using id::AtomID;
	switch ( type_ ) {
	case INTERNAL :
		return AtomID( atomno_, seqpos );
	case POLYMER_LOWER :
		return AtomID( conformation.residue_type( seqpos-1 ).upper_connect_atom(), seqpos - 1 ); //TODO: take out the seqpos-1 / seqpos+1 assumption, here.
	case POLYMER_UPPER :
		return AtomID( conformation.residue_type( seqpos+1 ).lower_connect_atom(), seqpos + 1 ); //TODO: take out the seqpos-1 / seqpos+1 assumption, here.
	case CONNECT :
		// in this case atomno_ is the connection identifer (index)
		return conformation.inter_residue_connection_partner( seqpos, atomno_ );
	default :
		return id::BOGUS_ATOM_ID;
	}
	return id::BOGUS_ATOM_ID;
}

//////////////////////////////////////////////////////////////////////////////////////////////////

AtomICoor::AtomICoor():
	//index_(0),
	phi_(0.0),
	theta_(0.0),
	d_(0.0),
	stub_atom1_(),
	stub_atom2_(),
	stub_atom3_()
{}

/// @brief constructor
AtomICoor::AtomICoor(
	std::string const & built_atom_name,
	Real const phi_in,
	Real const theta_in,
	Real const d_in,
	std::string const & stub_atom1_name,
	std::string const & stub_atom2_name,
	std::string const & stub_atom3_name,
	ResidueType const & rsd_type
):
	phi_( phi_in ),
	theta_( theta_in ),
	d_( d_in ),
	stub_atom1_( stub_atom1_name, rsd_type ),
	stub_atom2_( stub_atom2_name, rsd_type ),
	stub_atom3_( stub_atom3_name, rsd_type )
{
	if ( built_atom_name.size() <= 4 ) {
		built_vd_ = rsd_type.atom_vertex(built_atom_name);
	} else {
		built_vd_ = ResidueType::null_vertex;
	}
}

/// @brief Vertex descriptor version
AtomICoor::AtomICoor(
	VD const & built_atom_vd,
	Real const phi_in,
	Real const theta_in,
	Real const d_in,
	VD const & stub_atom1_vd,
	VD const & stub_atom2_vd,
	VD const & stub_atom3_vd,
	ResidueType const & rsd_type
):
	//index_(0),
	built_vd_( built_atom_vd ),
	phi_( phi_in ),
	theta_( theta_in ),
	d_( d_in ),
	stub_atom1_( stub_atom1_vd, rsd_type ),
	stub_atom2_( stub_atom2_vd, rsd_type ),
	stub_atom3_( stub_atom3_vd, rsd_type )
{}

Vector
AtomICoor::build(
	conformation::Residue const & rsd,
	conformation::Conformation const & conformation
) const
{

	debug_assert( kinematics::Stub( stub_atom1_.xyz( rsd, conformation ),
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
	debug_assert( kinematics::Stub( stub_atom1_.xyz( rsd_type ),
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

	debug_assert( built_stub.is_orthogonal( 0.001 ) );

	return built_stub.spherical( phi_, theta_, d_ );
}


typedef basic::datacache::DataMapObj< std::set< std::string > > AtomMemo;
typedef utility::pointer::shared_ptr< AtomMemo > AtomMemoOP;

//Memoized version, implementation of the general interface found below.

void pretty_print_atomicoor(std::ostream & out, AtomICoor const & start, ResidueType const & rsd_type, core::Size indent, AtomMemoOP memo ) {
	if ( ! memo ) {
		memo = AtomMemoOP( new AtomMemo );
	}
	for ( core::Size ii(1); ii <= indent; ++ii ) {
		out << "   ";
	}
	if ( indent > 0 ) {
		out << "* ";
	}
	VD built_vert(start.built_atom_vertex());
	if ( built_vert == ResidueType::null_vertex || built_vert == 0 ) {
		out << " -- UNNAMED POINT" << std::endl;
		return;
	}
	if ( ! rsd_type.has(built_vert) ) {
		out << " -- !!!! VERTEX NOT IN RESIUDE !!! " << built_vert << " for " << rsd_type.name() << std::endl;
		return;
	}
	std::string const & name( rsd_type.atom_name( built_vert ) );
	if ( memo->obj.count(name) ) {
		out << name << " -- SEEN BEFORE " << std::endl;
		return;
	}
	memo->obj.insert( name );
	out << name
		<< '\t' << start.phi() * numeric::NumericTraits< double >::rad2deg()
		<< '\t' << start.theta() * numeric::NumericTraits< double >::rad2deg()
		<< '\t' << start.d() << '\t';
	for ( core::Size ii(1); ii <= 3; ++ii ) {
		ICoorAtomID icoorid( start.stub_atom(ii) );
		switch( icoorid.type() ) {
		case ICoorAtomID::INTERNAL :
			if( rsd_type.has( icoorid.vertex() ) ) {
				out << rsd_type.atom_name( icoorid.vertex() ) << "  ";
			} else {
				out << "(missing vd: " << icoorid.vertex() <<")  ";
			}
			break;
		case ICoorAtomID::CONNECT :
			out << "CONNECT" << "  ";
			break;
		case ICoorAtomID::POLYMER_LOWER :
			out << "LOWER" << "  ";
			break;
		case ICoorAtomID::POLYMER_UPPER :
			out << "UPPER" << "  ";
			break;
		default :
			out << "UNKNOWN" << "  ";
			break;
		}
	}
	out << std::endl;
	for ( core::Size n(1); n <= rsd_type.natoms(); ++n ) {
		if ( rsd_type.icoor(n).stub_atom1().vertex() == built_vert ) {
			if ( rsd_type.icoor(n).built_atom_vertex() == built_vert ) {
				continue; // The root lists itself as it's stub_atom1
			}
			pretty_print_atomicoor(out, rsd_type.icoor(n), rsd_type, indent+1, memo );
		}
	}
}

void pretty_print_atomicoor(std::ostream & out, ResidueType const & rsd_type) {
	pretty_print_atomicoor(out, rsd_type.icoor( rsd_type.root_atom() ), rsd_type, 0, 0);
}

void pretty_print_atomicoor(std::ostream & out, AtomICoor const & start, ResidueType const & rsd_type, core::Size indent) {
	pretty_print_atomicoor(out, start, rsd_type, indent, 0);
}


} // chemical
} // core
