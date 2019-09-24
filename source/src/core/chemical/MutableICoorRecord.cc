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

// Unit headers
#include <core/chemical/MutableICoorRecord.hh>
#include <core/chemical/MutableResidueType.hh>

// Project headers
#include <core/chemical/AtomICoor.fwd.hh>
#include <core/chemical/MutableResidueConnection.hh>
#include <core/kinematics/Stub.hh>

// Numeric headers

// Utility headers

//#include <ObjexxFCL/string.functions.hh>
#include <basic/Tracer.hh>

//#include <utility/vector1.hh>
//#include <numeric/xyz.functions.hh>
//#include <numeric/NumericTraits.hh>

// C++ headers

#ifdef    SERIALIZATION
//#include <core/chemical/ResidueGraphTypes.srlz.hh>

// Utility serialization headers
#include <utility/serialization/serialization.hh>
#include <cereal/types/string.hpp>
#endif // SERIALIZATION

namespace core {
namespace chemical {

static basic::Tracer TR( "core.chemical.MutableICoorRecord" );

//////////////////////////////////////////////////////////////////////////////////////////////////

MutableICoorRecord::MutableICoorRecord():
	phi_(0.0),
	theta_(0.0),
	d_(0.0),
	stub_atom1_(),
	stub_atom2_(),
	stub_atom3_(),
	stub_type1_( ICoordAtomIDType::INTERNAL ),
	stub_type2_( ICoordAtomIDType::INTERNAL ),
	stub_type3_( ICoordAtomIDType::INTERNAL )
{}

/// @brief constructor
MutableICoorRecord::MutableICoorRecord(
	Real const phi_in,
	Real const theta_in,
	Real const d_in,
	std::string const & stub_atom1_name,
	std::string const & stub_atom2_name,
	std::string const & stub_atom3_name
):
	phi_( phi_in ),
	theta_( theta_in ),
	d_( d_in ),
	stub_atom1_( stub_atom1_name ),
	stub_atom2_( stub_atom2_name ),
	stub_atom3_( stub_atom3_name ),
	stub_type1_( string_to_icoord_type( stub_atom1_name ) ),
	stub_type2_( string_to_icoord_type( stub_atom2_name ) ),
	stub_type3_( string_to_icoord_type( stub_atom3_name ) )
{}

bool
MutableICoorRecord::operator==( MutableICoorRecord const & rhs ) const {
	return phi_ == rhs.phi_ &&
		theta_ == rhs.theta_ &&
		d_ == rhs.d_ &&
		stub_atom1_ == rhs.stub_atom1_ &&
		stub_atom2_ == rhs.stub_atom2_ &&
		stub_atom3_ == rhs.stub_atom3_;
}

/// @brief Can valid coordinates be built for this MutableICoorRecord, given the residue type?
bool
MutableICoorRecord::buildable( MutableResidueType const & rsd_type, bool verbose ) const {
	for ( int ii(1); ii <= 3; ++ii ) {
		std::string const & stub = stub_atom( ii );
		switch ( stub_type(ii) ) {
		case ICoordAtomIDType::INTERNAL :
			if ( ! rsd_type.has( stub ) ) {
				if ( verbose ) {
					TR.Error << "Atom `" << stub << "` missing for building ICoor for residue type " << rsd_type.name() << std::endl;
				}
				return false;
			}
			break;
		case ICoordAtomIDType::POLYMER_LOWER :
			if ( rsd_type.lower_connect_id() == 0 ) {
				if ( verbose ) {
					TR.Error << "Residue type " << rsd_type.name() << " does not have a lower connection, needed for building ICoor." << std::endl;
				}
				return false;
			}
			break;
		case ICoordAtomIDType::POLYMER_UPPER :
			if ( rsd_type.upper_connect_id() == 0 ) {
				if ( verbose ) {
					TR.Error << "Residue type " << rsd_type.name() << " does not have an upper connection, needed for building ICoor." << std::endl;
				}
				return false;
			}
			break;
		case ICoordAtomIDType::CONNECT :
			{ // scope for conn_num
			auto conn_num = get_connection_number(stub);
			if ( conn_num == 0 || conn_num > rsd_type.n_possible_residue_connections() ) {
				if ( verbose ) {
					TR.Error << "Residue type " << rsd_type.name() << " does not have an connection number " << conn_num << " needed for building ICoor." << std::endl;
				}
				return false;
			}
			break;
		}
		default :
			if ( verbose ) {
				TR.Error << "Stub type " << stub << " not understood for building ICoor." << std::endl;
			}
			return false;
		}
	}
	return true;
}

Vector
MutableICoorRecord::build( MutableResidueType const & rsd_type ) const {
	Vector v1( xyz( 1, rsd_type ) );
	Vector v2( xyz( 2, rsd_type ) );
	Vector v3( xyz( 3, rsd_type ) );
	return kinematics::Stub::create_orthogonal( v1, v2, v3 ).spherical( phi_, theta_, d_ );
}

Vector
MutableICoorRecord::xyz( core::Size stubno, MutableResidueType const & rsd_type ) const {
	switch ( stub_type(stubno) ) {
	case ICoordAtomIDType::INTERNAL :
		return rsd_type.atom( stub_atom(stubno) ).ideal_xyz();
	case ICoordAtomIDType::POLYMER_LOWER :
		return rsd_type.lower_connect().icoor().build( rsd_type );
	case ICoordAtomIDType::POLYMER_UPPER :
		return rsd_type.upper_connect().icoor().build( rsd_type );
	case ICoordAtomIDType::CONNECT :
		return rsd_type.residue_connection( get_connection_number(stub_atom(stubno)) ).icoor().build( rsd_type );
	default :
		utility_exit_with_message( "unrecognized stub atom id type!" );
	}
}

Vector
MutableICoorRecord::build_xyz( std::string const & stub, MutableResidueType const & rsd_type ) {
	switch ( string_to_icoord_type( stub ) ) {
	case ICoordAtomIDType::INTERNAL :
		return rsd_type.atom( stub ).ideal_xyz();
	case ICoordAtomIDType::POLYMER_LOWER :
		return rsd_type.lower_connect().icoor().build( rsd_type );
	case ICoordAtomIDType::POLYMER_UPPER :
		return rsd_type.upper_connect().icoor().build( rsd_type );
	case ICoordAtomIDType::CONNECT :
		return rsd_type.residue_connection( get_connection_number(stub) ).icoor().build( rsd_type );
	default :
		utility_exit_with_message( "unrecognized stub atom id type!" );
	}
}

void
MutableICoorRecord::show( std::ostream & out ) const {
	out << "MutableICoorRecord " << phi_ << "  " << theta_ << "  " << d_ << "  "
		<< stub_atom1_ << " " << stub_atom2_ << stub_atom3_ << std::endl;
}

} // chemical
} // core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::MutableICoorRecord::save( Archive & arc ) const {
	arc( CEREAL_NVP( phi_ ) ); // Real
	arc( CEREAL_NVP( theta_ ) ); // Real
	arc( CEREAL_NVP( d_ ) ); // Real
	arc( CEREAL_NVP( stub_atom1_ ) ); // std::string
	arc( CEREAL_NVP( stub_atom2_ ) ); // std::string
	arc( CEREAL_NVP( stub_atom3_ ) ); // std::string
	arc( CEREAL_NVP( stub_type1_ ) ); // ICoordAtomIDType
	arc( CEREAL_NVP( stub_type2_ ) ); // ICoordAtomIDType
	arc( CEREAL_NVP( stub_type3_ ) ); // ICoordAtomIDType
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::MutableICoorRecord::load( Archive & arc ) {
	arc( phi_ ); // Real
	arc( theta_ ); // Real
	arc( d_ ); // Real
	arc( stub_atom1_ ); // std::string
	arc( stub_atom2_ ); // std::string
	arc( stub_atom3_ ); // std::string
	arc( stub_type1_ ); // ICoordAtomIDType
	arc( stub_type2_ ); // ICoordAtomIDType
	arc( stub_type3_ ); // ICoordAtomIDType
}
SAVE_AND_LOAD_SERIALIZABLE( core::chemical::MutableICoorRecord );
#endif // SERIALIZATION
