// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/ResidueNeighborList.hh
/// @brief  A container class for use by the Etable and FA_Elec classes for storing
///         lists of atom neighbors
/// @author Andrew Leaver-Fay

// Unit headers
#include <core/scoring/ResidueNeighborList.hh>

// Package headers
#include <core/scoring/etable/count_pair/CountPairFunction.hh>

// Project headers
#include <core/conformation/Atom.hh>
#include <core/conformation/Residue.hh>

/// Utility headers
#include <utility/vector0.hh>

#include <utility/vector1.hh>


/// #define APL_TEMP_DEBUG

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {

ResidueNblistData::ResidueNblistData() {}
ResidueNblistData::~ResidueNblistData() {}

ResiduePairNeighborList::CacheableDataOP
ResidueNblistData::clone() const
{
	return ResiduePairNeighborList::CacheableDataOP( new ResidueNblistData( *this ) );
}

void ResidueNblistData::initialize(
	conformation::Residue const & res,
	CountPairFunctionCOP cpfxn,
	Real vvd2,
	Real hvd2,
	Real hhd2
) {
	atom_neighbors_.clear();

	if ( ! cpfxn ) return;

	utility::vector0< Real > cutoffs( 3 );
	cutoffs[ 0 ] = vvd2;
	cutoffs[ 1 ] = hvd2;
	cutoffs[ 2 ] = hhd2;

	for ( Size ii = 1; ii < res.natoms(); ++ii ) {
		int ii_isH = res.atom_is_hydrogen( ii );
		for ( Size jj = ii+1; jj <= res.natoms(); ++jj ) {
			int jj_isH = res.atom_is_hydrogen( jj );
			Real weight( 1.0 );
			Size path_dist( 0 );
			if ( cpfxn->count( ii, jj, weight, path_dist ) ) {
				Real d2 = res.xyz( ii ).distance_squared( res.xyz( jj ) );
				if ( d2 < cutoffs[ ii_isH + jj_isH ] ) {
					atom_neighbors_.push_back( SmallAtNb( ii, jj, path_dist, weight ) );
				}
			}
		}
	}
}

ResiduePairNeighborList::ResiduePairNeighborList() {}
ResiduePairNeighborList::~ResiduePairNeighborList() {}

ResiduePairNeighborList::CacheableDataOP ResiduePairNeighborList::clone() const { return ResiduePairNeighborList::CacheableDataOP( new ResiduePairNeighborList( *this ) ); }

void
ResiduePairNeighborList::initialize_from_residues(
	Real vvd2,
	Real hvd2,
	Real hhd2,
	conformation::Residue const & r1,
	conformation::Residue const & r2,
	etable::count_pair::CountPairFunctionCOP cpfxn
) {
	atom_neighbors_.clear();

	//std::cout << "ResiduePairNeighborList::initialize_from_residues " << r1.seqpos() << " " << r2.seqpos() << std::endl;

	utility::vector0< Real > cutoffs( 3 );
	cutoffs[ 0 ] = vvd2;
	cutoffs[ 1 ] = hvd2;
	cutoffs[ 2 ] = hhd2;

	for ( Size ii = 1; ii <= r1.natoms(); ++ii ) {
		int ii_isH = r1.atom_is_hydrogen( ii );
		for ( Size jj = 1; jj <= r2.natoms(); ++jj ) {
			int jj_isH = r2.atom_is_hydrogen( jj );
			Real weight( 1.0 );
			Size path_dist( 0 );
			if ( cpfxn->count( ii, jj, weight, path_dist ) ) {
				Real d2 = r1.xyz( ii ).distance_squared( r2.xyz( jj ) );
				///std::cout << "  atoms "  << ii << " " << jj << " " << r1.atom_name( ii ) << " " << r2.atom_name( jj )
				/// << " d2: " << d2 << " d " << std::sqrt( d2 ) <<  " cutoff2:" <<  cutoffs[ ii_isH + jj_isH ]
				/// << " cutoff " << std::sqrt( cutoffs[ ii_isH + jj_isH ] ) << std::endl;
				if ( d2 < cutoffs[ ii_isH + jj_isH ] ) {
					atom_neighbors_.push_back( SmallAtNb( ii, jj, path_dist, weight ) );
				}
			}
		}
	}
}

void
ResiduePairNeighborList::initialize_from_residues(
	Real vvd2,
	Real hvd2,
	Real hhd2,
	conformation::Residue const & r1,
	conformation::Residue const & r2,
	etable::count_pair::CountPairFunctionCOP cpfxn,
	std::map<core::Size,core::Size> const &r1_map,
	std::map<core::Size,core::Size> const &r2_map
) {
	atom_neighbors_.clear();

	//std::cout << "ResiduePairNeighborList::initialize_from_residues " << r1.seqpos() << " " << r2.seqpos() << std::endl;

	utility::vector0< Real > cutoffs( 3 );
	cutoffs[ 0 ] = vvd2;
	cutoffs[ 1 ] = hvd2;
	cutoffs[ 2 ] = hhd2;

	for ( Size ii = 1; ii <= r1.natoms(); ++ii ) {
		std::map<core::Size,core::Size>::const_iterator iter_i = r1_map.find(ii);
		int ii_rep = iter_i == r1_map.end() ? ii : iter_i->second;
		int ii_isH = r1.atom_is_hydrogen( ii );
		for ( Size jj = 1; jj <= r2.natoms(); ++jj ) {
			std::map<core::Size,core::Size>::const_iterator iter_j = r2_map.find(jj);
			int jj_rep = iter_j == r2_map.end() ? jj : iter_j->second;
			int jj_isH = r2.atom_is_hydrogen( jj );
			Real weight( 1.0 );
			Size path_dist( 0 );
			if ( cpfxn->count( ii_rep, jj_rep, weight, path_dist ) ) {
				Real d2 = r1.xyz( ii ).distance_squared( r2.xyz( jj ) );
				///std::cout << "  atoms "  << ii << " " << jj << " " << r1.atom_name( ii ) << " " << r2.atom_name( jj )
				/// << " d2: " << d2 << " d " << std::sqrt( d2 ) <<  " cutoff2:" <<  cutoffs[ ii_isH + jj_isH ]
				/// << " cutoff " << std::sqrt( cutoffs[ ii_isH + jj_isH ] ) << std::endl;
				if ( d2 < cutoffs[ ii_isH + jj_isH ] ) {
					atom_neighbors_.push_back( SmallAtNb( ii, jj, path_dist, weight ) );
				}
			}
		}
	}
}


}
}


#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
core::scoring::SmallAtNb::SmallAtNb() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::SmallAtNb::save( Archive & arc ) const {
	arc( CEREAL_NVP( atomno1_ ) ); // Size
	arc( CEREAL_NVP( atomno2_ ) ); // Size
	arc( CEREAL_NVP( path_dist_ ) ); // Size
	arc( CEREAL_NVP( weight_ ) ); // Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::SmallAtNb::load( Archive & arc ) {
	arc( atomno1_ ); // Size
	arc( atomno2_ ); // Size
	arc( path_dist_ ); // Size
	arc( weight_ ); // Real
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::SmallAtNb );

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::ResidueNblistData::save( Archive & arc ) const {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( CEREAL_NVP( atom_neighbors_ ) ); // utility::vector1<SmallAtNb>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::ResidueNblistData::load( Archive & arc ) {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( atom_neighbors_ ); // utility::vector1<SmallAtNb>
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::ResidueNblistData );
CEREAL_REGISTER_TYPE( core::scoring::ResidueNblistData )


/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::ResiduePairNeighborList::save( Archive & arc ) const {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( CEREAL_NVP( atom_neighbors_ ) ); // utility::vector1<SmallAtNb>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::ResiduePairNeighborList::load( Archive & arc ) {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( atom_neighbors_ ); // utility::vector1<SmallAtNb>
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::ResiduePairNeighborList );
CEREAL_REGISTER_TYPE( core::scoring::ResiduePairNeighborList )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_ResidueNeighborList )
#endif // SERIALIZATION
