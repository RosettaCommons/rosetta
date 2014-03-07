// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampling/rna/sugar/StepWiseRNA_VirtualSugarJustInTimeInstantiatorOLD.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_rna_StepWiseRNA_VirtualSugarJustInTimeInstantiatorOLD_HH
#define INCLUDED_protocols_stepwise_rna_StepWiseRNA_VirtualSugarJustInTimeInstantiatorOLD_HH

#include <protocols/moves/Mover.hh>
#include <protocols/stepwise/sampling/rna/sugar/legacy/StepWiseRNA_VirtualSugarJustInTimeInstantiatorOLD.fwd.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_ModelerOptions.fwd.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_JobParameters.fwd.hh>
#include <protocols/stepwise/sampling/rna/sugar/SugarModeling.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

using namespace core;

namespace protocols {
namespace stepwise {
namespace sampling {
namespace rna {
namespace sugar {
namespace legacy {

	class StepWiseRNA_VirtualSugarJustInTimeInstantiatorOLD: public protocols::moves::Mover {

	public:

		//constructor
		StepWiseRNA_VirtualSugarJustInTimeInstantiatorOLD( StepWiseRNA_JobParametersCOP & job_parameters_ );

		//destructor
		~StepWiseRNA_VirtualSugarJustInTimeInstantiatorOLD();

		virtual void apply( pose::Pose & pose_to_visualize );

		virtual std::string get_name() const;

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		bool
		do_the_sampling( core::pose::Pose & pose );

		bool
		sampling_sugar() const;

		bool
		sampling_sugar_at_chain_break() const;

		void
		prepare_from_prior_sampled_sugar_jobs( pose::Pose const & pose,
																					 utility::vector1< pose::PoseOP > & starting_pose_data_list,
																					 bool const pose_explosion_legacy = true );

		void
		prepare_from_prior_sampled_sugar_jobs_for_chain_break( pose::Pose const & pose,
																													utility::vector1< pose::PoseOP > & starting_pose_data_list );

		SugarModeling const &
		anchor_sugar_modeling();

		void
		set_scorefxn( core::scoring::ScoreFunctionOP const & scorefxn );

		void
		set_keep_base_fixed( bool const & setting ) { keep_base_fixed_ = setting;	}

		bool const & success() const{ return success_; }

		void
		set_options( StepWiseRNA_ModelerOptionsCOP options );

	private:

		bool
		initialize_parameters( pose::Pose const & pose );

		bool
		do_consistency_checks();

		bool
		do_sugar_sampling( pose::Pose & viewer_pose, SugarModeling & sugar_modeling, std::string const name );

		bool
		setup_sugar_modeling( pose::Pose const & pose, Size const moving_res, SugarModeling & sugar_modeling );

		bool
		is_anchor_sugar_virtual( pose::Pose const & pose ) ;

		bool
		is_current_sugar_virtual( pose::Pose const & pose ) ;

		bool
		is_five_prime_chain_break_sugar_virtual( pose::Pose const & pose ) ;

		bool
		is_three_prime_chain_break_sugar_virtual( pose::Pose const & pose ) ;

		void
		instantiate_sugar( pose::Pose & pose, SugarModeling const & sugar_modeling, Size const sugar_ID );

	private:

		StepWiseRNA_JobParametersCOP job_parameters_;
		bool const is_prepend_;
		Size const moving_res_;
		Size const five_prime_chain_break_res_;
		Size const three_prime_chain_break_res_;
		Size const num_nucleotides_;
		Size const gap_size_;
		bool const rebuild_bulge_mode_;

		bool is_anchor_sugar_virt_;
		bool is_current_sugar_virt_;
		bool is_five_prime_chain_break_sugar_virt_;
		bool is_three_prime_chain_break_sugar_virt_;
		Size num_virtual_sugar_;

		SugarModeling anchor_sugar_modeling_;
		SugarModeling current_sugar_modeling_;
		SugarModeling five_prime_chain_break_sugar_modeling_;
		SugarModeling three_prime_chain_break_sugar_modeling_;

		bool keep_base_fixed_;
		StepWiseRNA_ModelerOptionsCOP options_;
		core::scoring::ScoreFunctionOP scorefxn_;

		std::map< Size, Size > reference_res_for_each_virtual_sugar_;

		bool success_;
	};

} //legacy
} //sugar
} //rna
} //sampling
} //stepwise
} //protocols

#endif
