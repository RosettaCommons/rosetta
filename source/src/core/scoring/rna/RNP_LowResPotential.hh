// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/RNP_LowResPotential.hh
/// @brief  Statistically derived RNP low resolution protential
/// @author Kalli Kappel

#ifndef INCLUDED_core_scoring_RNP_LowResPotential_hh
#define INCLUDED_core_scoring_RNP_LowResPotential_hh

#include <core/types.hh>

// Unit headers
#include <core/scoring/rna/RNP_LowResPotential.fwd.hh>

// Package headers
#include <core/conformation/Residue.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

#include <basic/datacache/CacheableData.hh>

// Utility headers

#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/FArray4D.hh>

#include <utility/vector1.hh>


// C++


namespace core {
namespace scoring {
namespace rna {

////////////////////////////////////////////////////////////////////////////////////////////////////
class RNP_LowResPotential : public utility::pointer::ReferenceCount {

public:
	RNP_LowResPotential();

	void
	initialize_rnp_base_pair();

	void
	initialize_rnp_pair();

	void
	initialize_rnp_aa_rna_backbone();

	//void
	//initialize_rnp_stack_xy();

	void
	evaluate_rnp_base_pair_score(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		Real const & x,
		Real const & y,
		Real & rnp_bp_score
	) const;

	//void
	//evaluate_rnp_stack_xy_score(
	// conformation::Residue const & rsd1,
	// conformation::Residue const & rsd2,
	// Real const & x,
	// Real const & y,
	// Real & rnp_stack_score
	//) const;

	Real
	evaluate_rnp_aa_rna_backbone_score(
		conformation::Residue const & protein_rsd,
		Real const & dist_to_backbone ) const;

	void
	evaluate_rnp_pair_score(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		bool const & rsd1_is_interface,
		bool const & rsd1_is_buried,
		bool const & rsd2_is_interface,
		bool const & rsd2_is_buried,
		Real const & d,
		Real & rnp_pair_score
	) const;

	// void
	// compute_centroid_environment(
	//  pose::Pose & pose
	// ) const;
	//
	// void
	// finalize( pose::Pose & pose ) const;
	//
	//
	// void
	// evaluate_env_and_cbeta_scores(
	//  pose::Pose const & pose,
	//  conformation::Residue const & rsd,
	//  Real & env_score,
	//  Real & cb_score6,
	//  Real & cb_score12
	// ) const;
	//
	//
	// void
	// evaluate_pair_and_cenpack_score(
	//  conformation::Residue const & rsd1,
	//  conformation::Residue const & rsd2,
	//  Real const cendist,
	//  Real & pair_contribution,
	//  Real & cenpack_contribution
	// ) const;

	//private: // functions
	//
	// void
	// fill_cenlist(
	//  CenListInfo & cenlist,
	//  Size const res1,
	//  Size const res2,
	//  Real const cendist
	// ) const;
	//
	// void
	// truncate_cenlist_values( CenListInfo & cenlist ) const;

private: // data

	Size const max_aa_;
	Size const max_base_;
	Size const num_xbins_;
	Size const num_ybins_;
	Size const num_dbins_;
	Size const num_backbone_dbins_;

	ObjexxFCL::FArray4D < Real > rnp_basepair_xy_;
	//ObjexxFCL::FArray4D < Real > rnp_stack_xy_;
	//ObjexxFCL::FArray3D < Real > rnp_pair_;
	ObjexxFCL::FArray3D< Real > rnp_pair_base_interface_protein_buried_;
	ObjexxFCL::FArray3D< Real > rnp_pair_base_interface_protein_notburied_;
	ObjexxFCL::FArray2D< Real > rnp_aa_rna_backbone_;

};

} // rna
} // scoring
} // core

#endif
