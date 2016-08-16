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
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <core/scoring/rna/RNA_ScoringInfo.fwd.hh>
#include <core/scoring/rna/data/RNA_DataInfo.hh>

// Package headers

// Project headers
#include <core/chemical/AA.hh>
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>

// Utility headers
#include <numeric/xyzMatrix.hh>

#include <utility/vector1.hh>

// C++

///////////////////////////////////////////////////////
// Keep track of some base geometry that is
// useful for RNA scoring.
///////////////////////////////////////////////////////

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {
namespace rna {

RNA_ScoringInfo::RNA_ScoringInfo()
{}


/// @details Copy constructors must copy all data, not just some...
RNA_ScoringInfo::RNA_ScoringInfo( RNA_ScoringInfo const & src ) :
	CacheableData(),
	rna_centroid_info_( src.rna_centroid_info_ ),
	rna_raw_base_base_info_( src.rna_raw_base_base_info_ ),
	rna_filtered_base_base_info_( src.rna_filtered_base_base_info_ ),
	rna_data_info_( src.rna_data_info_ ),
	atom_numbers_for_vdw_calculation_( src.atom_numbers_for_vdw_calculation_ ),
	atom_numbers_for_mg_calculation_( src.atom_numbers_for_mg_calculation_ ),
	is_magnesium_( src.is_magnesium_ ),
	vdw_calculation_annotated_sequence_( src.vdw_calculation_annotated_sequence_ ),
	mg_calculation_annotated_sequence_( src.mg_calculation_annotated_sequence_ )
{
}


/// @details Pose must already contain a rna_scoring_info object or this method will fail.
RNA_ScoringInfo const &
rna_scoring_info_from_pose( pose::Pose const & pose )
{
	//using core::pose::datacache::CacheableDataType::RNA_SCORING_INFO;

	debug_assert( pose.data().has( core::pose::datacache::CacheableDataType::RNA_SCORING_INFO ) );
	return *( utility::pointer::static_pointer_cast< core::scoring::rna::RNA_ScoringInfo const > ( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::RNA_SCORING_INFO ) ) );
}

/// @details Either returns a non-const reference to the rna_scoring object already stored
/// in the pose, or creates a new rna scoring info object, places it in the pose, and returns
/// a non-const reference to it.
RNA_ScoringInfo &
nonconst_rna_scoring_info_from_pose( pose::Pose & pose )
{
	if ( pose.data().has( core::pose::datacache::CacheableDataType::RNA_SCORING_INFO ) ) {
		return *( utility::pointer::static_pointer_cast< core::scoring::rna::RNA_ScoringInfo > ( pose.data().get_ptr( core::pose::datacache::CacheableDataType::RNA_SCORING_INFO ) ) );
	}
	// else
	RNA_ScoringInfoOP rna_scoring_info( new RNA_ScoringInfo() );

	// This should get initialized and dimensioned later -- when it is filled!
	rna_scoring_info->rna_raw_base_base_info().resize( pose.total_residue() );

	pose.data().set( core::pose::datacache::CacheableDataType::RNA_SCORING_INFO, rna_scoring_info );
	return *rna_scoring_info;
}


/// @details remove RNA scoring info object so that it will be re-initialized if necessary
void
clear_rna_scoring_info( pose::Pose & pose )
{
	pose.data().clear( core::pose::datacache::CacheableDataType::RNA_SCORING_INFO );
}


} //rna
} //scoring
} //core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::rna::RNA_ScoringInfo::save( Archive & arc ) const {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( CEREAL_NVP( rna_centroid_info_ ) ); // rna::RNA_CentroidInfo
	arc( CEREAL_NVP( rna_raw_base_base_info_ ) ); // rna::RNA_RawBaseBaseInfo
	arc( CEREAL_NVP( rna_filtered_base_base_info_ ) ); // rna::RNA_FilteredBaseBaseInfo
	arc( CEREAL_NVP( rna_data_info_ ) ); // rna::data::RNA_DataInfo
	arc( CEREAL_NVP( atom_numbers_for_vdw_calculation_ ) ); // utility::vector1<utility::vector1<Size> >
	arc( CEREAL_NVP( atom_numbers_for_mg_calculation_ ) ); // utility::vector1<utility::vector1<Size> >
	arc( CEREAL_NVP( is_magnesium_ ) ); // utility::vector1<_Bool>
	arc( CEREAL_NVP( vdw_calculation_annotated_sequence_ ) ); // std::string
	arc( CEREAL_NVP( mg_calculation_annotated_sequence_ ) ); // std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::rna::RNA_ScoringInfo::load( Archive & arc ) {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( rna_centroid_info_ ); // rna::RNA_CentroidInfo
	arc( rna_raw_base_base_info_ ); // rna::RNA_RawBaseBaseInfo
	arc( rna_filtered_base_base_info_ ); // rna::RNA_FilteredBaseBaseInfo
	arc( rna_data_info_ ); // rna::data::RNA_DataInfo
	arc( atom_numbers_for_vdw_calculation_ ); // utility::vector1<utility::vector1<Size> >
	arc( atom_numbers_for_mg_calculation_ ); // utility::vector1<utility::vector1<Size> >
	arc( is_magnesium_ ); // utility::vector1<_Bool>
	arc( vdw_calculation_annotated_sequence_ ); // std::string
	arc( mg_calculation_annotated_sequence_ ); // std::string
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::rna::RNA_ScoringInfo );
CEREAL_REGISTER_TYPE( core::scoring::rna::RNA_ScoringInfo )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_rna_RNA_ScoringInfo )
#endif // SERIALIZATION
