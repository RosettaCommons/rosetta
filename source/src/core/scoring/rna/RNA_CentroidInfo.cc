// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/RNA_BaseBasePotential.cc
/// @brief  Statistically derived rotamer pair potential class implementation
/// @author Rhiju Das

// Unit headers
#include <core/scoring/rna/RNA_CentroidInfo.hh>

// Package headers

// Project headers
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>

// Utility headers
#include <numeric/xyzMatrix.hh>

#include <utility/vector1.hh>
#include <core/chemical/rna/util.hh>

using namespace core::chemical::rna;

// C++

///////////////////////////////////////////////////////
// Keep track of some base geometry that is
// useful for RNA scoring.
///////////////////////////////////////////////////////

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Numeric serialization headers
#include <numeric/xyz.serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {
namespace rna {

/// @details Copy constructors must copy all data, not just some...
RNA_CentroidInfo::RNA_CentroidInfo( RNA_CentroidInfo const & src ) :
	CacheableData()
{
	base_centroids_ = src.base_centroids_;
	base_stubs_ = src.base_stubs_;
	calculated_ = src.calculated_;
}

using Matrix = numeric::xyzMatrix<Real>;


///////////////////////////////////////////////////////////////////////////////
Vector
RNA_CentroidInfo::get_base_centroid( conformation::Residue const & rsd ) const
{
	bool const verbose = false;

	return get_rna_base_centroid( rsd, verbose );
}

//////////////////////////////////////////////////////////////////////////////////////
kinematics::Stub
RNA_CentroidInfo::get_base_coordinate_system( conformation::Residue const & rsd ) const
{
	Vector const centroid = get_base_centroid( rsd );
	return kinematics::Stub( get_rna_base_coordinate_system( rsd, centroid ), centroid );
}

//////////////////////////////////////////////////////////////////////////////////////
kinematics::Stub
RNA_CentroidInfo::get_base_coordinate_system( conformation::Residue const & rsd, Vector const & centroid ) const
{
	return kinematics::Stub( get_rna_base_coordinate_system( rsd, centroid ), centroid );
}

//////////////////////////////////////////////////////////////////////////////////////
void
RNA_CentroidInfo::initialize_base_centroids_and_stubs( pose::Pose const & pose )
{
	base_centroids_.clear();
	base_stubs_.clear();

	for ( Size i = 1; i <= pose.size(); i++ ) {
		conformation::Residue const & res_i  ( pose.residue( i ) );

		Vector centroid_i( 0.0  );
		kinematics::Stub stub_i;

		if ( res_i.is_NA() ) {
			centroid_i = get_base_centroid( res_i );
			stub_i     = get_base_coordinate_system( res_i, centroid_i );
		}

		base_centroids_.push_back( centroid_i );
		base_stubs_.push_back( stub_i );
	}
}

//////////////////////////////////////////////////////////////////////////////////////
void
RNA_CentroidInfo::update( pose::Pose const & pose )
{
	initialize_base_centroids_and_stubs( pose );
}


} //rna
} //scoring
} //core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::rna::RNA_CentroidInfo::save( Archive & arc ) const {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( CEREAL_NVP( base_centroids_ ) ); // utility::vector1<Vector>
	arc( CEREAL_NVP( base_stubs_ ) ); // utility::vector1<kinematics::Stub>
	arc( CEREAL_NVP( calculated_ ) ); // _Bool
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::rna::RNA_CentroidInfo::load( Archive & arc ) {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( base_centroids_ ); // utility::vector1<Vector>
	arc( base_stubs_ ); // utility::vector1<kinematics::Stub>
	arc( calculated_ ); // _Bool
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::rna::RNA_CentroidInfo );
CEREAL_REGISTER_TYPE( core::scoring::rna::RNA_CentroidInfo )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_rna_RNA_CentroidInfo )
#endif // SERIALIZATION
