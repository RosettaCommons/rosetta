// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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

//#include <utility/vector1.hh>
//#include <numeric/xyz.functions.hh>
#include <numeric/NumericTraits.hh>

// C++ headers

#ifdef    SERIALIZATION
#include <core/chemical/ResidueGraphTypes.srlz.hh>

// Utility serialization headers
#include <utility/serialization/serialization.hh>
#endif // SERIALIZATION

namespace core {
namespace chemical {

static basic::Tracer TR( "core.chemical.AtomICoor" );

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
		debug_assert( ObjexxFCL::is_int( name.substr(4) ) );
		atomno_ = ObjexxFCL::int_of( name.substr(4) );
		vd_ = ResidueType::null_vertex;
		debug_assert( atomno_ > 0 && atomno_ <= rsd_type.n_possible_residue_connections() );
		if ( atomno_ > rsd_type.n_possible_residue_connections() ) { // using > is not a great check but it is better than !=
			TR.Warning << "The record for CONN" << atomno_ << " in the topology file for " << rsd_type.name() <<
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

ICoorAtomID::ICoorAtomID( ICoorAtomID const & ) = default;

void
ICoorAtomID::remap_atom_vds( std::map< VD, VD > const & old_to_new ) {
	if ( old_to_new.count( vd_ ) == 1 ) {
		vd_ = old_to_new.at( vd_ );
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////
Vector const &
ICoorAtomID::xyz(
	Residue const & rsd,
	Conformation const & conformation
) const
{
	static const Vector NullVector( 0, 0, 0 );

	switch ( type_ ) {
	case INTERNAL :
		return rsd.atom( atomno_ ).xyz();
	case POLYMER_LOWER : // These three cases are deliberately combined
	case POLYMER_UPPER :
	case CONNECT :
		{ // Variable scoping
		id::AtomID connected_id( this->atom_id( rsd, conformation ) );
		if ( connected_id == id::AtomID::BOGUS_ATOM_ID() ) {
			return NullVector;
		}
		return conformation.xyz( connected_id );
	}
	default :
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
	static const Vector NullVector( 0, 0, 0 );
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
ICoorAtomID::atom_id( Residue const & rsd, Conformation const & conformation ) const
{
	using id::AtomID;
	switch ( type_ ) {
	case INTERNAL :
		return AtomID( atomno_, rsd.seqpos() );
	case POLYMER_LOWER :
		// Don't assume that the POLYMER_LOWER connection is to the i-1 residue in sequence,
		// nor to the upper connect of the connection partner
		{ // variable scoping
		debug_assert( ! rsd.is_lower_terminus() );
		core::Size const lowerID = rsd.type().lower_connect_id(); //The index of the lower connection of THIS residue.
		// lowerID may be zero even if the residue is not lower terminus, in which case we can't
		// even USE rsd.residue_connection_partner
		if ( lowerID == 0 || lowerID > rsd.connect_map_size() ) {
			TR.Warning << "IcoorAtomID::atom_id(): Cannot get atom_id for POLYMER_LOWER of residue " << rsd.name() << " " << rsd.seqpos() << " because it has no lower_connect_id(). Returning BOGUS ID instead." << std::endl;
			return id::AtomID::BOGUS_ATOM_ID();
		}
		core::Size const partner_seqpos( rsd.residue_connection_partner( lowerID ) ); //Get the index of the residue connected to this one at the lower connection.
		if ( partner_seqpos == 0 || partner_seqpos > conformation.size() ) {
			TR.Warning << "IcoorAtomID::atom_id(): Cannot get atom_id for POLYMER_LOWER of residue " << rsd.name() << " " << rsd.seqpos() << ".  Returning BOGUS ID instead." << std::endl;
			return id::AtomID::BOGUS_ATOM_ID();
		}
		//Get the index of the atom on the other residue that forms a connection to the lower terminus of this residue.
		core::Size const atomno( conformation.residue_type(partner_seqpos).residue_connect_atom_index( rsd.connect_map( lowerID ).connid()  ) );
		return AtomID( atomno, partner_seqpos );
	}
	case POLYMER_UPPER :
		// Don't assume that the POLYMER_UPPER connection is to the i+1 residue in sequence,
		// nor to the lower connect of the connection partner
		{ // variable scoping
		debug_assert(!rsd.is_upper_terminus());
		core::Size const upperID = rsd.type().upper_connect_id(); //The index of the upper connection of THIS residue.
		// upperID may be zero even if the residue is not lower terminus, in which case we can't
		// even USE rsd.residue_connection_partner
		if ( upperID == 0 || upperID > rsd.connect_map_size() ) {
			TR.Warning << "IcoorAtomID::atom_id(): Cannot get atom_id for POLYMER_UPPER of residue " << rsd.name() << " " << rsd.seqpos() << " because it has no upper_connect_id(). Returning BOGUS ID instead." << std::endl;
			return id::AtomID::BOGUS_ATOM_ID();
		}
		core::Size const partner_seqpos( rsd.residue_connection_partner( upperID ) ); //Get the index of the residue connected to this one at the upper connection.
		if ( partner_seqpos == 0 || partner_seqpos > conformation.size() ) {
			TR.Warning << "IcoorAtomID::atom_id(): Cannot get atom_id for POLYMER_UPPER of residue " << rsd.name() << " " << rsd.seqpos() << ". BOGUS ID instead." << std::endl;
			return id::AtomID::BOGUS_ATOM_ID();
		}
		//Get the index of the atom on the other residue that forms a connection to the upper terminus of this residue.
		core::Size const atomno( conformation.residue_type(partner_seqpos).residue_connect_atom_index( rsd.connect_map( upperID ).connid()  ) );
		return id::AtomID( atomno, partner_seqpos );
	}
	case CONNECT :
		return rsd.inter_residue_connection_partner( atomno_, conformation );
	default :
		utility_exit_with_message( "unrecognized stub atom id type!" );
	}
	return id::AtomID::BOGUS_ATOM_ID(); // To appease the compiler
}

//////////////////////////////////////////////////////////////////////////////////////////////////
bool
ICoorAtomID::buildable( Residue const & rsd, Conformation const & conformation ) const
{
	switch ( type_ ) {
	case INTERNAL :
		return (atomno_ != 0) && ( atomno_ <= rsd.natoms() );
	case POLYMER_LOWER :
	case POLYMER_UPPER :
	case CONNECT :
		return atom_id(rsd, conformation) != id::GLOBAL_BOGUS_ATOM_ID;
	default :
		return false;
	}
	return false;
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

/// @brief AtomICoorID version
AtomICoor::AtomICoor(
	std::string const & built_atom_name,
	Real const phi_in,
	Real const theta_in,
	Real const d_in,
	ICoorAtomID const & stub_atom1,
	ICoorAtomID const & stub_atom2,
	ICoorAtomID const & stub_atom3,
	ResidueType const & rsd_type
) :
	phi_( phi_in ),
	theta_( theta_in ),
	d_( d_in ),
	stub_atom1_( stub_atom1 ),
	stub_atom2_( stub_atom2 ),
	stub_atom3_( stub_atom3 )
{
	if ( built_atom_name.size() <= 4 ) {
		built_vd_ = rsd_type.atom_vertex(built_atom_name);
	} else {
		built_vd_ = ResidueType::null_vertex;
	}
}

void
AtomICoor::remap_atom_vds( std::map< VD, VD > const & old_to_new ) {
	if ( old_to_new.count( built_vd_) == 1 ) {
		built_vd_ = old_to_new.at( built_vd_ );
	}
	stub_atom1_.remap_atom_vds(old_to_new);
	stub_atom2_.remap_atom_vds(old_to_new);
	stub_atom3_.remap_atom_vds(old_to_new);
}


Vector
AtomICoor::build(
	conformation::Residue const & rsd,
	conformation::Conformation const & conformation
) const
{
	return build( stub_atom1_.xyz( rsd, conformation ),
		stub_atom2_.xyz( rsd, conformation ),
		stub_atom3_.xyz( rsd, conformation ) );
}

Vector
AtomICoor::build(
	ResidueType const & rsd_type
) const
{
	return build( stub_atom1_.xyz( rsd_type ),
		stub_atom2_.xyz( rsd_type ),
		stub_atom3_.xyz( rsd_type ) );
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
	return build( stub_atom1_.xyz( rsd ),
		stub_atom2_.xyz( rsd ),
		stub_atom3_.xyz( rsd ));
}

Vector
AtomICoor::build(
	Vector v1, Vector v2, Vector v3
) const
{
	kinematics::Stub built_stub( v1, v2, v3 );

	if ( ! built_stub.is_orthogonal( 0.001 ) ) {
		// Throw in a tiny shift and try again.
		core::Real delta = 1e-6;
		kinematics::Stub delta_stub(
			v1 + Vector( delta, 0, 0 ),
			v2 + Vector( 0, delta, 0 ),
			v3 + Vector( 0, 0, delta ) );
		built_stub = delta_stub;
	}

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
	if ( built_vert == ResidueType::null_vertex || built_vert == nullptr ) {
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
			if ( rsd_type.has( icoorid.vertex() ) ) {
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
	pretty_print_atomicoor(out, rsd_type.icoor( rsd_type.root_atom() ), rsd_type, 0, nullptr);
}

void pretty_print_atomicoor(std::ostream & out, AtomICoor const & start, ResidueType const & rsd_type, core::Size indent) {
	pretty_print_atomicoor(out, start, rsd_type, indent, nullptr);
}


} // chemical
} // core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::ICoorAtomID::save( Archive & arc ) const {
	arc( CEREAL_NVP( type_ ) ); // enum core::chemical::ICoorAtomID::Type
	arc( CEREAL_NVP( atomno_ ) ); // Size
	SERIALIZE_VD( arc, vd_, "vd_" ); // VD;
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::ICoorAtomID::load( Archive & arc ) {
	arc( type_ ); // enum core::chemical::ICoorAtomID::Type
	arc( atomno_ ); // Size
	DESERIALIZE_VD( arc, vd_ ); // VD;
}
SAVE_AND_LOAD_SERIALIZABLE( core::chemical::ICoorAtomID );

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::AtomICoor::save( Archive & arc ) const {
	SERIALIZE_VD( arc, built_vd_, "built_vd_" ); // VD;
	arc( CEREAL_NVP( phi_ ) ); // Real
	arc( CEREAL_NVP( theta_ ) ); // Real
	arc( CEREAL_NVP( d_ ) ); // Real
	arc( CEREAL_NVP( stub_atom1_ ) ); // class core::chemical::ICoorAtomID
	arc( CEREAL_NVP( stub_atom2_ ) ); // class core::chemical::ICoorAtomID
	arc( CEREAL_NVP( stub_atom3_ ) ); // class core::chemical::ICoorAtomID
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::AtomICoor::load( Archive & arc ) {
	DESERIALIZE_VD( arc, built_vd_ ); // VD;
	arc( phi_ ); // Real
	arc( theta_ ); // Real
	arc( d_ ); // Real
	arc( stub_atom1_ ); // class core::chemical::ICoorAtomID
	arc( stub_atom2_ ); // class core::chemical::ICoorAtomID
	arc( stub_atom3_ ); // class core::chemical::ICoorAtomID
}
SAVE_AND_LOAD_SERIALIZABLE( core::chemical::AtomICoor );
#endif // SERIALIZATION
