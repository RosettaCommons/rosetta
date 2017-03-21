// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


/// @author Oliver Lange


#ifndef INCLUDED_protocols_abinitio_IterativeBase_hh
#define INCLUDED_protocols_abinitio_IterativeBase_hh

// Unit Headers
//#include <protocols/abinitio/IterativeAbrelax.fwd.hh>

// Package Headers
#include <protocols/jd2/archive/NormalizedEvaluatedArchive.hh>
#include <protocols/jd2/archive/ArchiveManager.fwd.hh>
#include <protocols/abinitio/HedgeArchive.hh>

// Project Headers
#include <protocols/abinitio/PairingStatistics.fwd.hh>
#include <protocols/loops/Loops.hh>

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/fragment/FragSet.fwd.hh>


#include <utility/options/OptionCollection.hh>


// Utility headers
//for dynamic patching
#include <utility/options/keys/FileVectorOptionKey.hh>
#include <utility/io/ozstream.fwd.hh>

// Third-party Headers
#include <boost/functional/hash.hpp>


//// C++ headers
#include <string>

#include <protocols/noesy_assign/NoesyModule.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace abinitio {

class IterativeBase : public jd2::archive::NormalizedEvaluatedArchive {
private:
	typedef jd2::archive::NormalizedEvaluatedArchive Parent;
protected:
	typedef utility::vector1< core::io::silent::SilentStructOP > SilentStructVector;
public:
	enum IterationStage {
		ENUMERATION = 1,
		TOPO_RESAMPLING,
		PURE_TOPO_RESAMPLING,
		STAGE2_RESAMPLING,
		NOESY_PHASEII_TOPO,
		NOESY_PHASEII_S2_RESAMPLING,
		CEN2FULLATOM,
		//  CEN2FULLATOM_NON_POOL_DECOYS,
		LAST_CENTROID_START = CEN2FULLATOM,
		//  FLEX_CORE_RESAMPLING,
		RIGID_CORE_RESAMPLING,
		FINISHED //keep last
	};

	IterativeBase(std::string name );
	~IterativeBase() override;

	/// @brief archive is finished when at last stage
	bool finished() const override { return stage_ >= finish_stage_; };

	/// @brief do initializing work that requires fully setup object here
	void initialize() override;

	/// @brief where to stop ?
	void set_finish_stage( IterationStage setting ) {
		finish_stage_ = setting;
	}

	/// @brief  calls increment_stage() if appropriate
	void test_for_stage_end();

	/// @brief overloaded to make input decoys appear the same as decoys coming from batches
	void init_from_decoy_set( core::io::silent::SilentFileData const& sfd ) override;

	/// @brief we are always ready to generate a new batch
	virtual bool ready_for_batch() const { return true; };

	/// @brief we are not interested in batches that were generated in old stages
	bool still_interested( jd2::archive::Batch const& batch ) const override;

	/// @brief generate a new batch, use different recipe according to current stage

	/// @brief generate a new batch, use different recipe according to current stage
	void generate_batch() override;
	core::Size generate_batch( jd2::archive::Batch&, core::Size repeat_id ) override;

	/// @brief while waiting for jobs to finish
	void idle() override;
	void rescore() override;
	void save_status( std::ostream& ) const override;
	void restore_status( std::istream& ) override;

	/// @brief overloaded to handel special convergence check 'pool_converged_rmsd'
	/// @brief add structure to Archive.. return false if structure is rejected.
	bool add_structure (
		core::io::silent::SilentStructOP new_decoy,
		core::io::silent::SilentStructOP alternative_decoy,
		jd2::archive::Batch const&
	) override;

	/// @brief setup JumpNrEvaluator
	void setup_default_evaluators();

	/// @brief overloaded so we can test for end of IterationStage after reading
	void read_structures(
		core::io::silent::SilentFileData& sfd,
		core::io::silent::SilentFileData& alternative_decoys,
		jd2::archive::Batch const& batch
	) override;


	/// @brief generate flags and stuff for the out-sourced evaluation ---> such that score_final column is returned for each decoy
	/// note needs to be public, since IterativeCentroid calls this from IterativeFullatom to prepare evaluation for soon to be full-atom decoys
	// cen2fullatom-stage ( stage5 )
	virtual void gen_evaluation_output( jd2::archive::Batch& batch, bool fullatom = false );

	// overload by IterativeCentroid for rerouting to fullatom archive
	virtual void gen_diversity_pool( jd2::archive::Batch&, bool fullatom = false );

	// cen2fullatom-stage ( stage5 )
	virtual void gen_dynamic_patches( jd2::archive::Batch& batch );

	virtual void update_noesy_filter_files(
		std::string const& current,
		bool fullatom
	);

	// /// @brief need to get these from the IterativeCentroid to IterativeFullatom at end of stage5 ;
	//  std::string const& first_noesy_fa_cst_file() const { return first_noesy_fa_cst_file_; }
protected:
	//void set_first_noesy_fa_cst_file( std::string setting ) { first_noesy_fa_cst_file_ = setting; }

	//  core::Real noesy_assign_float_cycle() const { return noesy_assign_float_cycle_; }
	void set_noesy_assign_float_cycle( core::Real setting ) { noesy_assign_float_cycle_ = setting; }
	bool never_switched_noe_filter_;
	loops::Loops scored_core_;

	bool super_quick_relax_of_centroids() const { return super_quick_relax_of_centroids_; }
	/// ------------- helper functions to be used from generate_batch() --------------------

	void gen_resample_topologies( jd2::archive::Batch& batch );
	void gen_start_structures( jd2::archive::Batch& batch );
	void gen_enumerate_pairings( jd2::archive::Batch& batch );
	void gen_resample_stage2( jd2::archive::Batch& batch );
	void gen_resample_fragments( jd2::archive::Batch& batch );
	void gen_cen2fullatom( jd2::archive::Batch& batch );
	void gen_cen2fullatom_non_pool_decoys( jd2::archive::Batch& batch );
	void collect_hedgeing_decoys_from_batches(
		jd2::archive::Batch const& batch,
		core::io::silent::SilentStructOPs& start_decoys,
		core::Real score_cut_per_batch
	);
	void add_fullatom_flags( jd2::archive::Batch& batch );

	/// actually run the assignment machinery (only after batch is started to keep archive from hogging the queue... )
	void reassign_noesy_data( jd2::archive::Batch& batch );

	/// generate cst-input from current assigned noesy data
	void gen_noe_assignments( jd2::archive::Batch& batch );

	/// some helpers for the helpers
	PairingStatisticsOP compute_beta_topology();
	void guess_pairings_from_secondary_structure(
		core::fragment::FragSet const& frags,
		std::string const& out_pairings_file,
		std::string const& out_frag_ss_file
	) const;
	void compute_cores();

	///these are set by the cmd-line options iterative::fa_score and iterative::fa_score_patch
	std::string const& fa_score() const {
		return fa_score_;
	}
	std::string const& fa_score_patch() const {
		return fa_score_patch_;
	}

	///these are set by the cmd-line options iterative::cen_score and iterative::cen_score_patch
	std::string const& cen_score() const {
		return cen_score_;
	}
	std::string const& cen_score_patch() const {
		return cen_score_patch_;
	}

	/// @brief this is set from score::atom_pair_constraint of the pool-scorefunction
	core::Real overall_cstfilter_weight() const {
		return overall_cstfilter_weight_;
	}

	/// @brief this is set from score::atom_pair_constraint of the pool-scorefunction
	void set_overall_cstfilter_weight( core::Real setting ) {
		overall_cstfilter_weight_ = setting;
	}

	///OBSOLET cores are computed by compute_cores() in idle()
	loops::Loops const& core( core::Size i ) {
		if ( i == 1 ) { return core15_; };
		if ( i == 2 ) { return core2_; };
		if ( i == 3 ) { return core3_; };
		if ( i == 4 ) { return core4_; };
		std::cerr << "Cannot handle a value of " << i << " in IterativeBase::core(). Must be 1-4." << std::endl;
		utility_exit_with_message("Improper value passed to IterativeBase::core().");
		return core2_; //happy compiler
	}

	/// @brief current stage?
	IterationStage stage() const {
		return stage_;
	}

	/// @brief needed for writing of psi-pred fiels (guess_pairings_from_secondary_structure)
	std::string const& target_sequence() const {
		return target_sequence_;
	}

	void set_stage( IterationStage setting ) {
		stage_ = setting;
	}


	/// @brief cluster structures with min_diversity_list_[ stage_ ] as cluster:radius
	void cluster();

	std::string const& chemshift_column() const {
		return chemshift_column_;
	}

	void test_broker_settings( jd2::archive::Batch const& batch );
	void setup_filter_cst( core::Real weight );

	virtual void collect_alternative_decoys(
		SilentStructs /*primary_decoys*/,
		std::string /*alternative_decoy_file*/,
		SilentStructVector& /*output_decoys*/
	) {};

private:

	void collect_hedge_structures( core::io::silent::SilentStructOP evaluated_decoy, jd2::archive::Batch const& batch );
	/// @brief score a pose with Pool-Scoring function (adds necessary data to pose (RDC, constraints,  etc ) )
	void score( core::pose::Pose& pose ) const override;

	/// @brief what is the expected lowest acceptance ratio at the current stage ?
	core::Real target_accept_ratio() const { return target_accept_ratio_[ stage_ ]; }

	/// @brief [OBSOLET] add score_coreX and rms_coreX evaluators (and columns) with 0.0 weight
	void add_core_evaluator( loops::Loops const& core, std::string const& core_tag );

	/// @brief necessary steps to go to next stage... e.g., saving snapshot of archive
	void increment_stage();

	// void read_noisy_assign_data_from_last_batch();

	void replace_noesy_filter_constraints();
	void rescore_nonlocal_archive();

	void setup_autoNOE();
	void do_dynamic_patching(
		jd2::archive::Batch& batch,
		utility::io::ozstream& flags,
		std::string score,
		utility::options::FileVectorOptionKey const& key
	) const;
private:
	///  ----------------- -- private data members -- --------------------

	/// ------------------------------ stage - control --------------------------
	/// @brief current stage
	IterationStage stage_;

	/// @brief end-condition
	IterationStage finish_stage_;

	/// @brief indices of prominent batches ( STATUS file )
	core::Size first_batch_this_stage_;
	core::Size first_fullatom_batch_;

	/// --------------------------------- other ------------------------------------
	/// @brief toggle to keep track of the enumerate-pairings mode .. want to run this only every 2nd batch
	bool bEnumeratedLastTime_;

	/// @brief [OBSOLET?] keep track when idle() has been run ...
	core::Size last_accepted_decoys_in_idle_;

	///core-regions --- used in IterativeFullatom for the "rigid-core" sampling step...
	loops::Loops core15_;
	loops::Loops core2_;
	loops::Loops core3_;
	loops::Loops core4_;

	/// ----------- some cmd-line controlled settings -----------------

	/// @brief how many structures are maximally produced in stage X
	utility::vector1< int > max_nstruct_list_; //-1 --> skip stage, 0 infinite, N>0 make a maximum of N structures.

	/// @brief cluster:radius for minimum diversity in stage X
	utility::vector1< core::Real > min_diversity_list_;

	/// @brief minimum acceptance ratio .. current_accept_ratio < target_accept_ratio[ stage_ ] --> increment_stage
	utility::vector1< core::Real > target_accept_ratio_;

	/// @brief for RMSD ... e.g., in:file:native
	core::pose::PoseCOP reference_pose_;

	/// @brief from cmd-line cen-score and patch names
	std::string cen_score_;
	std::string cen_score_patch_;

	/// @brief from cmd-line fa-score and patch names
	std::string fa_score_;
	std::string fa_score_patch_;

	/// @brief this is set from score::atom_pair_constraint of the pool-scorefunction
	core::Real overall_cstfilter_weight_;

	std::string target_sequence_; //read from in:file:fasta in c'stor

	protocols::noesy_assign::NoesyModuleOP noesy_module_;
	core::Real noesy_assign_float_cycle_;

	std::string first_noesy_cst_file_;
	std::string first_noesy_fa_cst_file_;

	std::string current_noesy_sampling_file_;
	bool bCombineNoesyCst_;


	//hash value to see if decoys for noesy assign have changed
	size_t noesy_assign_hash_;
	boost::hash<std::string> hasher;

	/// @brief even in centroid mode the end of abinitio will have a fast relax... enables cs-score and noe-assignment
	// KAB - below line commented out by warnings removal script (-Wunused-private-field) on 2014-09-11
	// bool recover_centroid_structure_for_pool_;

	std::string chemshift_column_;
	bool bDoBetaJumping_;

	/// @brief even in centroid mode the end of abinitio will have a fast relax... enables cs-score and noe-assignment
	bool super_quick_relax_of_centroids_;

	/// @brief use the score_variations from EvaluatedArchive to determine new sampling weights
	bool use_dynamic_weights_for_sampling_;

	core::Size delay_noesy_reassign_;

	HedgeArchiveOP hedge_archive_;
	/// @brief cache some of the experimental data so we don't reload from file for each evaluation
	mutable core::scoring::ResidualDipolarCouplingOP rdc_data_; //need to cache this to avoid reading RDC file each time...
	mutable core::scoring::constraints::ConstraintSetOP cst_data_;
	mutable core::scoring::constraints::ConstraintSetOP cst_fa_data_;
	/// ------------------ register cmdline options ---------------------------

	utility::options::OptionCollection const vanilla_options_; //options before stage-dpd auto-noe-options were added


private:
	static bool options_registered_;
public:
	static void register_options();

};


}
}

#endif
