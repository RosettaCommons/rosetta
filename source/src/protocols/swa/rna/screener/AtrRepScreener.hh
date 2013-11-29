// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/swa/rna/screener/AtrRepScreener.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_swa_rna_screener_AtrRepScreener_HH
#define INCLUDED_protocols_swa_rna_screener_AtrRepScreener_HH

#include <utility/pointer/ReferenceCount.hh>
#include <protocols/swa/rna/screener/AtrRepScreener.fwd.hh>
#include <protocols/swa/rna/StepWiseRNA_Classes.hh>
#include <protocols/swa/rna/StepWiseRNA_JobParameters.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pose/Pose.hh>

using namespace core;

namespace protocols {
namespace swa {
namespace rna {
namespace screener {

	class AtrRepScreener: public utility::pointer::ReferenceCount {

	public:

		//Constructor
		AtrRepScreener( pose::Pose const & pose, StepWiseRNA_JobParametersCOP & job_parameters );

		//destructor
		~AtrRepScreener();

		void
		set_sample_both_sugar_base_rotamer( bool const & setting ){ sample_both_sugar_base_rotamer_ = setting; }

		Real delta_atr_score() const{ return delta_atr_score_; }
		Real delta_rep_score() const{ return delta_rep_score_; }

	public:

		void
		get_base_atr_rep_score( core::pose::Pose const & pose );

		bool
		check_screen( pose::Pose & current_pose_screen,
									Size const & gap_size,
									bool const & is_internal,
									bool const & kic_sampling );

	private:

		StepWiseRNA_JobParametersCOP job_parameters_; //need to use the full_to_sub map...should convert to const style.. Parin Feb 28, 2010
		Real const rep_cutoff_;
		Real base_atr_score_;
		Real base_rep_score_;
		Real delta_atr_score_;
		Real delta_rep_score_;
		bool sample_both_sugar_base_rotamer_;
		bool output_pdb_;
		bool verbose_;

		core::scoring::ScoreFunctionOP atr_rep_screening_scorefxn_;

		StepWiseRNA_CountStruct count_data_;

	};

} //screener
} //rna
} //swa
} //protocols

#endif
