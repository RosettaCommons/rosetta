// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software legacy and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampling/rna/legacy/suite/StepWiseRNA_ClassicResidueSampler.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_rna_StepWiseRNA_ClassicResidueSampler_HH
#define INCLUDED_protocols_stepwise_rna_StepWiseRNA_ClassicResidueSampler_HH

#include <protocols/stepwise/sampling/rna/legacy/suite/StepWiseRNA_ClassicResidueSampler.fwd.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_ModelerOptions.fwd.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_JobParameters.fwd.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_Classes.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_PoseSelection.hh>
#include <protocols/stepwise/sampling/rna/o2prime/O2PrimePacker.fwd.hh>
#include <protocols/stepwise/sampling/rna/phosphate/MultiPhosphateSampler.fwd.hh>
#include <protocols/stepwise/sampling/rna/checker/RNA_BaseCentroidChecker.fwd.hh>
#include <protocols/stepwise/sampling/rna/checker/RNA_VDW_BinChecker.fwd.hh>
#include <protocols/stepwise/sampling/rna/checker/AtrRepChecker.fwd.hh>
#include <protocols/stepwise/sampling/rna/checker/ChainClosureChecker.fwd.hh>
#include <protocols/stepwise/sampling/rna/checker/ChainClosableGeometryChecker.fwd.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <protocols/rotamer_sampler/RotamerBase.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/types.hh>
#include <protocols/moves/MoverForPoseList.hh>

using namespace core;
using namespace core::pose;
using namespace protocols::stepwise::sampling::rna;

namespace protocols {
namespace stepwise {
namespace sampling {
namespace rna {
namespace legacy {
namespace suite {

	class StepWiseRNA_ClassicResidueSampler: public protocols::moves::MoverForPoseList {

	public:

		//constructor
		StepWiseRNA_ClassicResidueSampler( StepWiseRNA_JobParametersCOP & job_parameters_ );

		//destructor
		~StepWiseRNA_ClassicResidueSampler();

		virtual void apply( core::pose::Pose & pose_to_visualize );

		virtual std::string get_name() const;

		using MoverForPoseList::apply;

	public:

		void
		set_extra_tag( std::string const setting ){ extra_tag_ = setting; }

		void set_silent_file( std::string const & setting ){ silent_file_ = setting; }

		void set_scorefxn( core::scoring::ScoreFunctionOP const & scorefxn );

		void set_pose_list( utility::vector1< PoseOP > &	pose_list );

		utility::vector1< PoseOP > & pose_list();

		void
		set_base_centroid_checker( checker::RNA_BaseCentroidCheckerOP & checker );

		void
		set_user_input_VDW_bin_checker( checker::RNA_VDW_BinCheckerOP const & user_input_VDW_bin_checker );

		void
		set_rebuild_bulge_mode( bool const & setting ){ rebuild_bulge_mode_ = setting; }

		void
		set_options( StepWiseRNA_ModelerOptionsCOP options );

	private:

		void
		initialize_poses_and_checkers( core::pose::Pose & pose );

		bool
		break_early_for_integration_tests();

		void
		output_count_data();

		Size
		get_num_pose_kept();

		bool
		bulge_variant_decision( core::pose::Pose & pose, Real const & delta_atr_score );

		void
		apply_bulge_variant( core::pose::Pose & pose ) const;

	private:

		StepWiseRNA_JobParametersCOP job_parameters_; //need to use the full_to_sub map...should convert to const style.. Parin Feb 28, 2010

		Size const moving_res_;
		Size const moving_suite_;
		bool const is_prepend_;
		bool const is_internal_;
		Size const actually_moving_res_;
		Size const gap_size_;
		utility::vector1 < core::Size > const working_moving_partition_pos_;
		Size const num_nucleotides_;
		bool const is_dinucleotide_;
		bool const close_chain_to_distal_;
		Size const five_prime_chain_break_res_;

		Size const last_append_res_;
		Size const last_prepend_res_;
		Real const atom_atom_overlap_dist_cutoff_;
		std::string extra_tag_;

		bool const build_pose_from_scratch_;
		bool kic_sampling_;
		bool rebuild_bulge_mode_;

		utility::vector1< PoseOP > pose_list_;
		core::scoring::ScoreFunctionOP scorefxn_;
		StepWiseRNA_CountStruct count_data_;

		std::string silent_file_;
		bool native_rmsd_screen_;
		core::Real const bin_size_; /*ALWAYS 20!!!*/
		checker::AtrRepCheckerOP atr_rep_checker_;
		checker::RNA_VDW_BinCheckerOP user_input_VDW_bin_checker_;
		checker::ChainClosableGeometryCheckerOP chain_closable_geometry_checker_;
		checker::ChainClosureCheckerOP chain_closure_checker_;
		checker::RNA_BaseCentroidCheckerOP base_centroid_checker_;
		pose::PoseOP screening_pose_;
		StepWiseRNA_PoseSelectionOP pose_selection_;
		o2prime::O2PrimePackerOP o2prime_packer_;
		phosphate::MultiPhosphateSamplerOP phosphate_sampler_;

		StepWiseRNA_ModelerOptionsCOP options_;

		rotamer_sampler::RotamerBaseOP rotamer_sampler_;

		};

} //suite
} //legacy
} //rna
} //sampling
} //stepwise
} //protocols

#endif
