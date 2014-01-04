// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/enumerate/rna/screener/AtrRepScreener.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_rna_screener_AtrRepScreener_HH
#define INCLUDED_protocols_stepwise_rna_screener_AtrRepScreener_HH

#include <utility/pointer/ReferenceCount.hh>
#include <protocols/stepwise/enumerate/rna/screener/AtrRepScreener.fwd.hh>
#include <protocols/stepwise/enumerate/rna/StepWiseRNA_Classes.hh>
#include <protocols/stepwise/enumerate/rna/StepWiseRNA_JobParameters.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pose/Pose.hh>

using namespace core;

namespace protocols {
namespace stepwise {
namespace enumerate {
namespace rna {
namespace screener {

	class AtrRepScreener: public utility::pointer::ReferenceCount {

	public:

		//Constructor
		AtrRepScreener( pose::Pose const & pose, StepWiseRNA_JobParametersCOP & job_parameters );

		AtrRepScreener( pose::Pose const & pose,
										Size const moving_res,
										Size const reference_res,
										Size const gap_size,
										bool const is_internal = false,
										bool const separate_moving_residue_to_estimate_baseline = true,
										bool const sample_both_sugar_base_rotamer = false );

		//destructor
		~AtrRepScreener();

		Real delta_atr_score() const{ return delta_atr_score_; }
		Real delta_rep_score() const{ return delta_rep_score_; }
		Real base_atr_score() const{ return base_atr_score_; }
		Real base_rep_score() const{ return base_rep_score_; }

	public:

		bool
		check_screen( pose::Pose & current_pose_screen );

		void
		set_kic_sampling( bool const & setting ) { kic_sampling_ = setting;	}

	private:

		void
		get_base_atr_rep_score( core::pose::Pose const & pose );

		void
		initialize_scorefxn();

		void
		initialize_parameters();

	private:

		Size const working_moving_res_;
		Size const working_reference_res_;
		Size const gap_size_;
		bool const is_prepend_;
		bool const is_internal_;
		bool const sample_both_sugar_base_rotamer_;
		bool const separate_moving_residue_to_estimate_baseline_;

		Real rep_cutoff_;
		Real base_atr_score_;
		Real base_rep_score_;
		Real delta_atr_score_;
		Real delta_rep_score_;
		bool output_pdb_;
		bool verbose_;
		bool kic_sampling_;

		core::scoring::ScoreFunctionOP atr_rep_screening_scorefxn_;

		StepWiseRNA_CountStruct count_data_;

	};

} //screener
} //rna
} //enumerate
} //stepwise
} //protocols

#endif
