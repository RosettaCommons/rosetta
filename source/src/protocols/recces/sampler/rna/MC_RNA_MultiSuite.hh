// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/recces/sampler/rna/MC_RNA_MultiSuite.hh
/// @brief Markov chain sampler for mulitple RNA suite.
/// @author Fang-Chieh Chou

#ifndef INCLUDED_protocols_sampler_rna_MC_RNA_MultiSuite_HH
#define INCLUDED_protocols_sampler_rna_MC_RNA_MultiSuite_HH

// Unit headers
#include <protocols/recces/sampler/rna/MC_RNA_MultiSuite.fwd.hh>

// Package headers
#include <protocols/recces/sampler/MC_Comb.hh>
#include <protocols/recces/sampler/rna/MC_RNA_Suite.fwd.hh>

// Project headers
#include <core/id/TorsionID.hh>

namespace protocols {
namespace recces {
namespace sampler {
namespace rna {

class MC_RNA_MultiSuite : public MC_Comb {
public:

	MC_RNA_MultiSuite();

	/// @brief Set the standard deviation of Gaussian sampler
	void set_gaussian_stdev( core::Real const setting );

	/// @brief Set the standard deviation of a sampler
	void set_gaussian_stdev(
		core::Real const setting,
		core::Size const rotamer_id
	);

	/// @brief Set the stored angle from a pose
	void set_angle( core::pose::Pose const & pose );

	/// @brief Set the flip rate of pucker
	void set_pucker_flip_rate( core::Real const setting );

	/// @brief Set the flip rate of pucker for a sampler
	void set_pucker_flip_rate(
		core::Real const setting,
		core::Size const rotamer_id
	);

	/// @brief Set starting torsions and pucker states from pose
	void set_init_from_pose( core::pose::Pose const & pose );

	/// @brief Clear all rotamer samplers stored in this sampler
	void clear_rotamer();

	/// @brief Add one more rotamer sampler to this sampler
	void add_external_loop_rotamer( MC_RNA_SuiteOP rotamer );

	/// @brief Name of the class
	std::string get_name() const { return "MC_RNA_MultiSuite"; }

	/// @brief Type of class (see enum in toolbox::SamplerPlusPlusTypes.hh)
	virtual toolbox::SamplerPlusPlusType type() const { return toolbox::MC_RNA_MULTI_SUITE; }

private:
	using MC_Comb::add_external_loop_rotamer;

	// same as MC_Comb::samplers_, but cast as MC_RNA_Suite
	utility::vector1<MC_RNA_SuiteOP> suite_samplers_;
};

} //rna
} //sampler
} //recces
} //protocols

#endif
