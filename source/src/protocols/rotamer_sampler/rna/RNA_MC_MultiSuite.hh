// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/rna/RNA_MC_MultiSuite.hh
/// @brief Markov chain sampler for mulitple RNA suite.
/// @author Fang-Chieh Chou

#ifndef INCLUDED_protocols_rotamer_sampler_rna_RNA_MC_MultiSuite_HH
#define INCLUDED_protocols_rotamer_sampler_rna_RNA_MC_MultiSuite_HH

// Unit headers
#include <protocols/rotamer_sampler/rna/RNA_MC_MultiSuite.fwd.hh>

// Package headers
#include <protocols/rotamer_sampler/MC_Comb.hh>
#include <protocols/rotamer_sampler/rna/RNA_MC_Suite.fwd.hh>

// Project headers
#include <core/id/TorsionID.hh>

namespace protocols {
namespace rotamer_sampler {
namespace rna {

class RNA_MC_MultiSuite : public MC_Comb {
public:

	RNA_MC_MultiSuite();

	/// @brief Set the standard deviation of Gaussian sampler
	void set_gaussian_stdev( core::Real const setting );

	/// @brief Set the standard deviation of a sampler
	void set_gaussian_stdev(
		core::Real const setting,
		core::Size const rotamer_id
	);

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
	void add_external_loop_rotamer( MC_RotamerSamplerOP const & rotamer );

	/// @brief Name of the class
	std::string get_name() const { return "RNA_MC_MultiSuite"; }

	/// @brief Type of class (see enum in RotamerSamplerTypes.hh)
	virtual RotamerSamplerType type() const { return RNA_MC_MULTI_SUITE; }

private:
	using MC_Comb::add_external_loop_rotamer;
	utility::vector1<RNA_MC_SuiteOP> suite_samplers_;
};

} //rna
} //rotamer_sampler
} //protocols

#endif
