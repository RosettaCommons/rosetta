// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/rna/data/RNA_DMS_LowResolutionPotential.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_core_scoring_rna_data_RNA_DMS_LowResolutionPotential_HH
#define INCLUDED_core_scoring_rna_data_RNA_DMS_LowResolutionPotential_HH

#include <utility/pointer/ReferenceCount.hh>
#include <core/scoring/rna/data/RNA_DMS_LowResolutionPotential.fwd.hh>
#include <core/pose/rna/RNA_DataInfo.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/hbonds/HBondSet.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <numeric/MathVector.hh>
#include <numeric/MathMatrix.hh>
#include <set>

namespace core {
namespace scoring {
namespace rna {
namespace data {

class RNA_DMS_LowResolutionPotential: public utility::pointer::ReferenceCount {

public:

	//constructor
	RNA_DMS_LowResolutionPotential();

	//destructor
	~RNA_DMS_LowResolutionPotential();

public:

	void
	initialize( core::pose::Pose & pose,
		bool const rna_base_pair_computed = false );

	core::Real
	evaluate( core::pose::Pose const & pose,
		pose::rna::RNA_Reactivity const & rna_reactivity );


	void
	update_edge_paired(  core::Size const i, core::Size const k,
		utility::vector1< bool > & wc_edge_paired,
		utility::vector1< bool > & hoogsteen_edge_paired,
		utility::vector1< bool > & sugar_edge_paired );

	void
	get_rna_base_pairing_status( core::pose::Pose & pose,
		utility::vector1< bool > & wc_edge_paired,
		utility::vector1< bool > & hoogsteen_edge_paired,
		utility::vector1< bool > & sugar_edge_paired,
		utility::vector1< bool > & is_bulged,
		bool const already_scored = false );

	bool
	get_wc_near_o2prime( core::pose::Pose const & pose, core::Size const i );

	utility::vector1< bool > const & wc_edge_paired() const { return wc_edge_paired_; }
	utility::vector1< bool > const & hoog_edge_paired() const { return hoog_edge_paired_; }
	utility::vector1< bool > const & sugar_edge_paired() const { return sugar_edge_paired_; }
	utility::vector1< bool > const & is_bulged() const { return is_bulged_; }

	void set_careful_base_pair_classifier( bool const & setting ){ careful_base_pair_classifier_ = setting; }
	bool careful_base_pair_classifier() const { return careful_base_pair_classifier_; }

private:

	void
	initialize_DMS_low_resolution_potential();

	numeric::MathVector< core::Real >
	read_DMS_low_resolution_stats_file( std::string const & potential_file );

	void
	figure_out_low_resolution_potential( numeric::MathMatrix< Real > & all_DMS_stats );

private:

	bool careful_base_pair_classifier_;
	utility::vector1< bool > wc_edge_paired_, hoog_edge_paired_, sugar_edge_paired_, is_bulged_;
	utility::vector1< bool > is_protected_values_;
	std::set< Real > DMS_values_;
	numeric::MathMatrix< Real > DMS_low_resolution_potential_;

};

} //data
} //rna
} //scoring
} //core

#endif
