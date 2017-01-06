// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/recces/sampler/rna/MC_RNA_KIC_Sampler.hh
/// @brief Sample and torsions and close an RNA loop.
/// @author Fang-Chieh Chou

// NOTE: The sampler makes the assumption that the sugar pucker of
// moving_suite + 1 is moving (with pucker and chi rotamer being sampled).
// This corresponds to the "Append" case in Parin's SWA code.
// I think the "Prepend" case is not closable using the current
// RNA_KinematicClosure class. But if we stick with "Append" then it does not
// matter.

#ifndef INCLUDED_protocols_sampler_rna_MC_RNA_KIC_Sampler_HH
#define INCLUDED_protocols_sampler_rna_MC_RNA_KIC_Sampler_HH

// Unit headers
#include <protocols/recces/sampler/rna/MC_RNA_KIC_Sampler.fwd.hh>

// Package headers
#include <protocols/toolbox/SamplerPlusPlus.hh>
#include <protocols/recces/sampler/MC_Sampler.hh>
#include <protocols/recces/sampler/MC_OneTorsion.fwd.hh>
#include <protocols/recces/sampler/MC_Comb.fwd.hh>
#include <protocols/recces/sampler/rna/MC_RNA_Sugar.fwd.hh>
#include <protocols/stepwise/monte_carlo/mover/TransientCutpointHandler.fwd.hh>
#include <protocols/stepwise/sampler/rna/RNA_KinematicCloser.fwd.hh>
#include <protocols/stepwise/sampler/rna/RNA_ChiStepWiseSampler.fwd.hh>
#include <protocols/stepwise/sampler/screener/RNA_TorsionScreener.fwd.hh>
#include <core/chemical/rna/util.hh>
#include <core/pose/Pose.hh>
#include <numeric/random/WeightedSampler.hh>
#include <numeric/random/random.fwd.hh>

namespace protocols {
namespace recces {
namespace sampler {
namespace rna {

class MC_RNA_KIC_Sampler : public MC_Sampler {
public:
	MC_RNA_KIC_Sampler(
		core::pose::PoseCOP mc_pose,
		core::Size const moving_suite,
		core::Size const chainbreak_suite,
		bool const change_ft = true
	);

	~MC_RNA_KIC_Sampler()
	{}

	/// @brief Initialization
	virtual void init();

	/// @brief Reset to the first (or random if is_random()) rotamer.
	virtual void reset() {return;}

	/// @brief Move to next rotamer
	virtual void operator++();

	/// @brief Update rotamer
	virtual void update();

	void
	get_next_solutions( core::pose::Pose const & pose );

	void
	choose_solution();

	void next( core::pose::Pose const & pose );

	/// @brief Did the pose move at all?
	bool check_moved() const;

	/// @brief Apply the current rotamer to pose
	void apply( core::pose::Pose & pose );

	/// @brief Name of the class
	virtual
	std::string get_name() const { return "MC_RNA_KIC_Sampler"; }

	/// @brief Type of class (see enum in toolbox::SamplerPlusPlusTypes.hh)
	virtual toolbox::SamplerPlusPlusType type() const { return toolbox::MC_RNA_KIC; }

	/// @brief Set the standard deviation of Gaussian sampler
	void set_gaussian_stdev( core::Real const setting );

	/// @brief Max # of step for trying in random samping
	void set_max_tries( core::Size const setting ) {
		set_and_reinit( max_tries_, setting );
	}

	void set_angle_range_from_init_torsions( core::Real const range );

	/// @brief output summary of class
	virtual
	void show( std::ostream & out, Size const indent = 0 ) const;

	/// @brief return OP to the subsampler that controls exactly this torsion_id (assume only one).
	virtual
	MC_SamplerOP
	find( core::id::TorsionID const & torsion_id );

private:

	core::Real get_jacobian( core::pose::Pose & pose );

	bool check_angles_in_range( const core::pose::Pose & pose );

	//Disable copy constructor and assignment
	MC_RNA_KIC_Sampler( const MC_RNA_KIC_Sampler & );
	void operator=( const MC_RNA_KIC_Sampler & );

	core::Real vector_sum( utility::vector1< core::Real > const & vector );

private:

	core::pose::Pose stored_pose_;

	core::Size const moving_suite_, chainbreak_suite_;
	core::Real gaussian_stdev_, sum_;
	core::Size max_tries_;
	bool verbose_, did_close_, used_current_solution_;

	stepwise::sampler::rna::RNA_KinematicCloserOP loop_closer_, stored_loop_closer_;
	stepwise::monte_carlo::mover::TransientCutpointHandlerOP cutpoint_handler_;
	utility::vector1<MC_OneTorsionOP> bb_samplers_, chi_samplers_;
	MC_CombOP chi_rotamer_;
	utility::vector1<MC_RNA_SugarOP> sugar_samplers_;
	utility::vector1< core::Real > stored_jacobians_, current_jacobians_;
	Size solution_;
	numeric::random::WeightedSampler jacobian_sampler_;
	utility::vector1<core::Real> all_jacobians_, copy_stored_jacobians_, initial_torsions_,
		angle_min_, angle_max_;
	utility::vector1< core::id::TorsionID > TorsionIDs;
};

} //rna
} //sampler
} //recces
} //protocols

#endif
