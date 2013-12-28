// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/RNA_BaseBasePotential.cc
/// @brief  Statistically derived rotamer pair potential class implementation
/// @author Rhiju Das

// Unit headers
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <core/scoring/rna/RNA_ScoringInfo.fwd.hh>
#include <core/scoring/rna/RNA_DataInfo.hh>

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

	assert( pose.data().has( core::pose::datacache::CacheableDataType::RNA_SCORING_INFO ) );
	return *( static_cast< RNA_ScoringInfo const * >( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::RNA_SCORING_INFO)() ) );
}

/// @details Either returns a non-const reference to the rna_scoring object already stored
/// in the pose, or creates a new rna scoring info object, places it in the pose, and returns
/// a non-const reference to it.
RNA_ScoringInfo &
nonconst_rna_scoring_info_from_pose( pose::Pose & pose )
{
	//using core::pose::datacache::CacheableDataType::RNA_SCORING_INFO;

	if ( pose.data().has( core::pose::datacache::CacheableDataType::RNA_SCORING_INFO ) ) {
		return *( static_cast< RNA_ScoringInfo * >( pose.data().get_ptr( core::pose::datacache::CacheableDataType::RNA_SCORING_INFO )() ));
	}
	// else
	RNA_ScoringInfoOP rna_scoring_info = new RNA_ScoringInfo();

	// This should get initialized and dimensioned later -- when it is filled!
	rna_scoring_info->rna_raw_base_base_info().resize( pose.total_residue() );

	pose.data().set( core::pose::datacache::CacheableDataType::RNA_SCORING_INFO, rna_scoring_info );
	return *rna_scoring_info;
}


}
}
}
