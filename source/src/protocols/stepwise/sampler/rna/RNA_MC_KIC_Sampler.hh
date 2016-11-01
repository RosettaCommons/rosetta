// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampler/rna/RNA_MC_KIC_Sampler.hh
/// @brief Sample and torsions and close an RNA loop.
/// @author Fang-Chieh Chou

// NOTE: The sampler makes the assumption that the sugar pucker of
// moving_suite + 1 is moving (with pucker and chi rotamer being sampled).
// This corresponds to the "Append" case in Parin's SWA code.
// I think the "Prepend" case is not closable using the current
// RNA_KinematicClosure class. But if we stick with "Append" then it does not
// matter.

#ifndef INCLUDED_protocols_sampler_rna_RNA_MC_KIC_Sampler_HH
#define INCLUDED_protocols_sampler_rna_RNA_MC_KIC_Sampler_HH

// Unit headers
#include <protocols/stepwise/sampler/rna/RNA_MC_KIC_Sampler.fwd.hh>

// Package headers
#include <protocols/stepwise/sampler/StepWiseSamplerBase.hh>
#include <protocols/stepwise/sampler/StepWiseSamplerSizedComb.fwd.hh>
#include <protocols/stepwise/sampler/rna/RNA_KinematicCloser_DB.fwd.hh>
#include <protocols/stepwise/sampler/rna/RNA_ChiStepWiseSampler.fwd.hh>
#include <protocols/stepwise/sampler/screener/RNA_TorsionScreener.fwd.hh>
#include <core/chemical/rna/util.hh>
#include <protocols/stepwise/sampler/MC_StepWiseSampler.hh>
#include <protocols/stepwise/sampler/MC_OneTorsion.fwd.hh>
#include <protocols/stepwise/sampler/MC_Comb.fwd.hh>
#include <protocols/stepwise/sampler/rna/RNA_MC_Sugar.fwd.hh>
#include <protocols/stepwise/monte_carlo/mover/TransientCutpointHandler.fwd.hh>
#include <numeric/random/WeightedSampler.hh>
#include <numeric/random/WeightedSampler.fwd.hh>
#include <numeric/random/random.fwd.hh>
#include <core/pose/Pose.hh>

namespace protocols {
namespace stepwise {
namespace sampler {
namespace rna {

class RNA_MC_KIC_Sampler : public MC_StepWiseSampler {
public:
	RNA_MC_KIC_Sampler(
		//  core::pose::Pose const & ref_pose,
		core::pose::PoseOP const & ref_pose,
		core::Size const moving_suite,
		core::Size const chainbreak_suite
	);
	
	RNA_MC_KIC_Sampler(
		//  core::pose::Pose const & ref_pose,
		core::pose::PoseOP const & ref_pose,
		core::Size const moving_suite,
		core::Size const chainbreak_suite,
		bool const change_foldtree
	);

	// ~RNA_MC_KIC_Sampler();

	/// @brief Initialization
	virtual void init();

	/// @brief Reset to the first (or random if is_random()) rotamer.
	virtual void reset() {return;}
	void reset( core::pose::Pose & pose );

	/// @brief Move to next rotamer
	virtual void operator++();

	/// @brief Update rotamer
	void update( const core::pose::Pose & pose );
	virtual void update();

	void next( const core::pose::Pose & pose );

	/// @brief Was the chain able to close?
	bool check_closed() const;

	/// @brief Did the pose move at all?
	bool check_moved() const;

	/// @brief Check if reach the end of rotamer list
	virtual bool not_end() const;

	/// @brief Apply the current rotamer to pose
	void apply( core::pose::Pose & pose );

	/// @brief Apply the current rotamer to ref_pose_
	void apply();

	/// @brief If the chain is closable (random modeler only)
	bool closable() const { return random_chain_closed_; }

	/// @brief Name of the class
	std::string get_name() const;

	/// @brief Type of class (see enum in StepWiseSamplerTypes.hh)
	virtual StepWiseSamplerType type() const { return RNA_KIC; }

	// Set functions
	void set_verbose( bool const setting ) {
		set_and_reinit( verbose_, setting );
	}

	/// @brief Set the standard deviation of Gaussian sampler
	void set_gaussian_stdev( core::Real const setting );

	void set_extra_epsilon( bool const setting ) {
		set_and_reinit( extra_epsilon_, setting );
	}

	void set_extra_chi( bool const setting ) {
		set_and_reinit( extra_chi_, setting );
	}

	void set_skip_same_pucker( bool const setting ) {
		set_and_reinit( skip_same_pucker_, setting );
	}

	void set_idealize_coord( bool const setting ) {
		set_and_reinit( idealize_coord_, setting );
	}

	void set_torsion_screen( bool const setting ) {
		set_and_reinit( torsion_screen_, setting );
	}

	void set_pucker_state( core::chemical::rna::PuckerState const setting ) {
		set_and_reinit( init_pucker_, setting );
		//set_and_reinit( pucker_state_, setting );
	}

	void set_base_state( core::chemical::rna::ChiState const setting ) {
		set_and_reinit( base_state_, setting );
	}

	void set_sample_nucleoside( core::Size const setting ) { sample_nucleoside_ = setting; }

	void set_fast( core::Real const setting ) {
		if ( setting ) {
			extra_chi_ = false;
			extra_epsilon_ = false;
		}
	}

	/// @brief Max # of step for trying in random samping
	void set_max_tries( core::Size const setting ) {
		set_and_reinit( max_tries_, setting );
	}

	core::Real get_jacobian( core::pose::Pose & pose );

	void set_angle_range_from_init_torsions( core::Real const range );
	bool check_angles_in_range( const core::pose::Pose & pose );

	core::Size moving_suite() const { return moving_suite_; }
	core::Size chainbreak_suite() const { return chainbreak_suite_; }

private:
	//Disable copy constructor and assignment
	RNA_MC_KIC_Sampler( const RNA_MC_KIC_Sampler & );
	void operator=( const RNA_MC_KIC_Sampler & );

	core::Real vector_sum( utility::vector1< core::Real > const & vector );

	core::Real jacobian_, angle_;
	core::pose::PoseOP ref_pose_;
	core::pose::Pose stored_pose_;
	// I think this is what I really want core::pose::Pose const & ref_pose_;
	// core::pose::PoseOP test_pose;
	//core::pose::PoseOP const test_pose;
	core::Size const moving_suite_, chainbreak_suite_;
	core::chemical::rna::PuckerState init_pucker_;
	core::Real pucker_flip_rate_;
	core::Real gaussian_stdev_, sum_;
	core::chemical::rna::ChiState base_state_;
	core::Size sample_nucleoside_;
	core::Size max_tries_;
	bool verbose_, extra_epsilon_, extra_chi_, skip_same_pucker_,
		idealize_coord_, torsion_screen_, random_chain_closed_, did_close,
	used_current_solution_, did_move;
	bool const change_foldtree_;

	// StepWiseSamplerSizedCombOP bb_rotamer_;
	MC_CombOP bb_rotamer_;
	// MC_Comb bb_rotamer_;
	RNA_KinematicCloser_DBOP loop_closer_, stored_loop_closer_, copy_loop_closer_;
	// RNA_ChiStepWiseSamplerOP chi_rotamer_;
	screener::RNA_TorsionScreenerOP screener_;
	monte_carlo::mover::TransientCutpointHandlerOP cutpoint_handler_, cutpoint_end_loop_;
	utility::vector1<MC_OneTorsionOP> bb_samplers_, chi_samplers_;
	// utility::vector1<MC_OneTorsionOP> chi_rotamer_;
	// MC_OneTorsionOP chi_rotamer_;
	MC_CombOP chi_rotamer_;
	utility::vector1<RNA_MC_SugarOP> sugar_samplers_;
	utility::vector1< core::Real > stored_jacobians_, current_jacobians_;
	Size solution_;
	numeric::random::WeightedSampler jacobian_sampler_;
	utility::vector1<core::Real> all_jacobians_, copy_stored_jacobians_, initial_torsions_,
		angle_min_, angle_max_;
	utility::vector1< core::id::TorsionID > TorsionIDs;
};

} //rna
} //sampler
} //stepwise
} //protocols

#endif
