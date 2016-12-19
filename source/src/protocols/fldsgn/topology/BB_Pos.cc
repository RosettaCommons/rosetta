// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ./src/protocols/fldsgn/topology/BB_Pos.cc
/// @brief
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu), Nobuyasu Koga ( nobuyasu@u.washington.edu )

/// Unit headers
#include <protocols/fldsgn/topology/BB_Pos.hh>

/// Package headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>

#include <utility/vector1.hh>


#ifdef    SERIALIZATION
// Project serialization headers
#include <core/chemical/ResidueType.srlz.hh>

// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Numeric serialization headers
#include <numeric/xyz.serialization.hh>
#endif // SERIALIZATION

namespace protocols {
namespace fldsgn {
namespace topology {

/// @details.  After a change in size, the residue types and the integer indices are all wrong.
/// Erase the old information.
void
BB_Pos::resize( Size const nres )
{
	if ( N_.size() == (Size) nres ) return;
	Vector NullVector( 0.0 );
	N_.resize( nres, NullVector );
	CA_.resize( nres, NullVector );
	C_.resize( nres, NullVector );
	O_.resize( nres, NullVector );
	CB_.resize( nres, NullVector );

	residue_types_.resize( nres ); std::fill( residue_types_.begin(), residue_types_.end(), core::chemical::ResidueTypeOP() );
	N_index_.resize( nres );       std::fill( N_index_.begin(), N_index_.end(), 0 );
	CA_index_.resize( nres );      std::fill( CA_index_.begin(), CA_index_.end(), 0 );
	CB_index_.resize( nres );      std::fill( CB_index_.begin(), CB_index_.end(), 0 );
	C_index_.resize( nres );       std::fill( C_index_.begin(), C_index_.end(), 0 );
	O_index_.resize( nres );       std::fill( O_index_.begin(), O_index_.end(), 0 );

}

/// @details.  After a change in size, the residue types and the integer indices are all wrong.
/// Erase the old information.
void
BB_Pos::clear()
{
	N_.clear();
	CA_.clear();
	C_.clear();
	O_.clear();
	CB_.clear();
	residue_types_.clear();
	N_index_.clear();
	CA_index_.clear();
	CB_index_.clear();
	C_index_.clear();
	O_index_.clear();
}

/// @details: Optimize the common case where the sequence of the pose is not changing from
/// score function evaluation to evaluation (e.g. abinitio!)
void
BB_Pos::take_coordinates_from_pose( Pose const & pose )
{
	if ( ! bbindices_up_to_date( pose ) ) {
		update_indices( pose );
	}

	for ( Size i=1; i<= pose.size(); ++i ) {
		core::conformation::Residue const & rsd( pose.residue(i) );
		if ( rsd.is_protein() ) {
			assert( N_index_[ i ]  );
			assert( CA_index_[ i ] );
			assert( C_index_[ i ]  );
			assert( O_index_[ i ]  );

			N_ [i] = rsd.xyz( N_index_[ i ]  );
			CA_[i] = rsd.xyz( CA_index_[ i ] );
			C_ [i] = rsd.xyz( C_index_[ i ]  );
			O_ [i] = rsd.xyz( O_index_[ i ]  );
			if ( rsd.aa() == core::chemical::aa_gly ) {
				CB_[i] = 0.0;
			} else {
				CB_[i] = rsd.xyz( CB_index_[ i ] );
			}
		} else {
			N_ [i] = 0.0;
			CA_[i] = 0.0;
			C_ [i] = 0.0;
			O_ [i] = 0.0;
			CB_[i] = 0.0;
		}
	}

}

bool
BB_Pos::bbindices_up_to_date( Pose const & pose ) const
{
	if ( N_.size() != pose.size() ) return false;

	for ( Size ii = 1; ii <= pose.size(); ++ii ) {
		if ( residue_types_[ ii ].get() != & ( pose.residue_type( ii ) ) ) return false;
	}
	return true;
}

void
BB_Pos::update_indices( Pose const & pose )
{
	resize( pose.size() );

	static std::string const bbN("N");
	static std::string const bbCA("CA");
	static std::string const bbC("C");
	static std::string const bbO("O");
	static std::string const scCB("CB");

	for ( Size i = 1; i <= pose.size(); ++i ) {
		core::conformation::Residue const & rsd( pose.residue(i) );
		residue_types_[ i ] = rsd.type().get_self_ptr();
		if ( rsd.is_protein() ) {

			N_index_[  i ] = rsd.atom_index( bbN );
			CA_index_[ i ] = rsd.atom_index( bbCA );
			C_index_[  i ] = rsd.atom_index( bbC );
			O_index_[  i ] = rsd.atom_index( bbO );

			if ( rsd.aa() == core::chemical::aa_gly ) {
				CB_index_[ i ] = 0;
			} else {
				CB_index_[ i ] = rsd.atom_index( scCB );
			}
		} else {
			N_index_[  i ] = 0;
			CA_index_[ i ] = 0;
			C_index_[  i ] = 0;
			O_index_[  i ] = 0;
			CB_index_[ i ] = 0;
		}
	}

}


} // ns topology
} // ns fldsgn
} // ns protocols


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::fldsgn::topology::BB_Pos::save( Archive & arc ) const {
	arc( CEREAL_NVP( N_ ) ); // utility::vector1<Vector>
	arc( CEREAL_NVP( CA_ ) ); // utility::vector1<Vector>
	arc( CEREAL_NVP( CB_ ) ); // utility::vector1<Vector>
	arc( CEREAL_NVP( C_ ) ); // utility::vector1<Vector>
	arc( CEREAL_NVP( O_ ) ); // utility::vector1<Vector>
	arc( CEREAL_NVP( residue_types_ ) );
	arc( CEREAL_NVP( N_index_ ) ); // utility::vector1<Size>
	arc( CEREAL_NVP( CA_index_ ) ); // utility::vector1<Size>
	arc( CEREAL_NVP( CB_index_ ) ); // utility::vector1<Size>
	arc( CEREAL_NVP( C_index_ ) ); // utility::vector1<Size>
	arc( CEREAL_NVP( O_index_ ) ); // utility::vector1<Size>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::fldsgn::topology::BB_Pos::load( Archive & arc )
{
	arc( N_ ); // utility::vector1<Vector>
	arc( CA_ ); // utility::vector1<Vector>
	arc( CB_ ); // utility::vector1<Vector>
	arc( C_ ); // utility::vector1<Vector>
	arc( O_ ); // utility::vector1<Vector>
	arc( CEREAL_NVP( residue_types_ ) );
	arc( N_index_ ); // utility::vector1<Size>
	arc( CA_index_ ); // utility::vector1<Size>
	arc( CB_index_ ); // utility::vector1<Size>
	arc( C_index_ ); // utility::vector1<Size>
	arc( O_index_ ); // utility::vector1<Size>
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::fldsgn::topology::BB_Pos );
#endif // SERIALIZATION
