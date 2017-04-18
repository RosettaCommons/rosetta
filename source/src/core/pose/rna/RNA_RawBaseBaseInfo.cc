// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pose/rna/RNA_RawBaseBaseInfo.cc
/// @brief  Statistically derived rotamer pair potential class implementation
/// @author Rhiju Das


// Unit headers
#include <core/pose/rna/RNA_RawBaseBaseInfo.hh>

// Package headers
#include <core/chemical/rna/util.hh>

// Project headers
#include <core/chemical/AA.hh>

// Utility headers

#include <utility/vector1.hh>


// C++

using namespace core::chemical::rna;


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// ObjexxFCL serialization headers
#include <utility/serialization/ObjexxFCL/FArray2D.srlz.hh>
#include <utility/serialization/ObjexxFCL/FArray3D.srlz.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace pose {
namespace rna {


/// @details Copy constructors must copy all data, not just some...
RNA_RawBaseBaseInfo::RNA_RawBaseBaseInfo( RNA_RawBaseBaseInfo const & src ) :
	CacheableData(),
	base_pair_array_(src.base_pair_array_),
	base_axis_array_(src.base_axis_array_),
	base_stagger_array_(src.base_stagger_array_),
	base_stack_array_(src.base_stack_array_),
	base_stack_axis_array_(src.base_stack_axis_array_),
	base_geometry_orientation_array_(src.base_geometry_orientation_array_),
	base_geometry_height_array_(src.base_geometry_height_array_),
	calculated_(src.calculated_)
{}

void
RNA_RawBaseBaseInfo::resize( Size const & total_residue )
{

	if ( base_pair_array_.size1() == total_residue ) return;

	base_pair_array_.dimension( total_residue, total_residue, 3 );
	base_stagger_array_.dimension( total_residue, total_residue, 3 );
	base_axis_array_.dimension( total_residue, total_residue, 3 );
	base_stack_array_.dimension( total_residue, total_residue );
	base_stack_axis_array_.dimension( total_residue, total_residue );
	base_geometry_orientation_array_.dimension( total_residue, total_residue );
	base_geometry_height_array_.dimension( total_residue, total_residue );
	zero();
}

////////////////////////////////////////////////////////
void
RNA_RawBaseBaseInfo::zero()
{
	base_pair_array_ = 0.0;
	base_stagger_array_ = 0.0;
	base_axis_array_ = 0.0;
	base_stack_array_ = 0.0;
	base_stack_axis_array_ = 0.0;
	base_geometry_orientation_array_ = 0.0;
	base_geometry_height_array_ = 0.0;
}

////////////////////////////////////////////////////////
void
RNA_RawBaseBaseInfo::copy_values( RNA_RawBaseBaseInfo const & src, Size const & i, Size const & j )
{
	for ( Size k = 1; k <= NUM_EDGES; k++ ) {
		base_pair_array_( i, j, k ) = src.base_pair_array_( i, j, k );
		base_stagger_array_( i, j, k ) = src.base_stagger_array_( i, j, k );
		base_axis_array_( i, j, k ) = src.base_axis_array_( i, j, k );
	}
	base_stack_array_( i, j ) = src.base_stack_array_( i, j );
	base_stack_axis_array_( i, j ) = src.base_stack_axis_array_( i, j );
	base_geometry_orientation_array_( i, j ) = src.base_geometry_orientation_array_( i, j );
	base_geometry_height_array_( i, j ) = src.base_geometry_height_array_( i, j );
}

} //rna
} //pose
} //core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::pose::rna::RNA_RawBaseBaseInfo::save( Archive & arc ) const {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( CEREAL_NVP( base_pair_array_ ) ); // ObjexxFCL::FArray3D<Real>; ObjexxFCL: ObjexxFCL::FArray3D<Real>
	arc( CEREAL_NVP( base_axis_array_ ) ); // ObjexxFCL::FArray3D<Real>; ObjexxFCL: ObjexxFCL::FArray3D<Real>
	arc( CEREAL_NVP( base_stagger_array_ ) ); // ObjexxFCL::FArray3D<Real>; ObjexxFCL: ObjexxFCL::FArray3D<Real>
	arc( CEREAL_NVP( base_stack_array_ ) ); // ObjexxFCL::FArray2D<Real>
	arc( CEREAL_NVP( base_stack_axis_array_ ) ); // ObjexxFCL::FArray2D<Real>
	arc( CEREAL_NVP( base_geometry_orientation_array_ ) ); // ObjexxFCL::FArray2D<Real>
	arc( CEREAL_NVP( base_geometry_height_array_ ) ); // ObjexxFCL::FArray2D<Real>
	arc( CEREAL_NVP( calculated_ ) ); // _Bool
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::pose::rna::RNA_RawBaseBaseInfo::load( Archive & arc ) {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( base_pair_array_ ); // ObjexxFCL::FArray3D<Real>; ObjexxFCL: ObjexxFCL::FArray3D<Real>
	arc( base_axis_array_ ); // ObjexxFCL::FArray3D<Real>; ObjexxFCL: ObjexxFCL::FArray3D<Real>
	arc( base_stagger_array_ ); // ObjexxFCL::FArray3D<Real>; ObjexxFCL: ObjexxFCL::FArray3D<Real>
	arc( base_stack_array_ ); // ObjexxFCL::FArray2D<Real>
	arc( base_stack_axis_array_ ); // ObjexxFCL::FArray2D<Real>
	arc( base_geometry_orientation_array_ ); // ObjexxFCL::FArray2D<Real>
	arc( base_geometry_height_array_ ); // ObjexxFCL::FArray2D<Real>
	arc( calculated_ ); // _Bool
}

SAVE_AND_LOAD_SERIALIZABLE( core::pose::rna::RNA_RawBaseBaseInfo );
CEREAL_REGISTER_TYPE( core::pose::rna::RNA_RawBaseBaseInfo )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_rna_RNA_RawBaseBaseInfo )
#endif // SERIALIZATION
