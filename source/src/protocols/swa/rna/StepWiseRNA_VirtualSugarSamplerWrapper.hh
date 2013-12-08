// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/swa/rna/StepWiseRNA_VirtualSugarSamplerWrapper.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_swa_rna_StepWiseRNA_VirtualSugarSamplerWrapper_HH
#define INCLUDED_protocols_swa_rna_StepWiseRNA_VirtualSugarSamplerWrapper_HH

#include <protocols/moves/Mover.hh>
#include <protocols/swa/rna/StepWiseRNA_VirtualSugarSamplerWrapper.fwd.hh>
#include <protocols/swa/rna/StepWiseRNA_JobParameters.fwd.hh>
#include <protocols/swa/rna/SugarModeling.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

using namespace core;

namespace protocols {
namespace swa {
namespace rna {

	class StepWiseRNA_VirtualSugarSamplerWrapper: public protocols::moves::Mover {

	public:

		//constructor
		StepWiseRNA_VirtualSugarSamplerWrapper( StepWiseRNA_JobParametersCOP & job_parameters_ );

		//destructor
		~StepWiseRNA_VirtualSugarSamplerWrapper();

		virtual void apply( pose::Pose & pose_to_visualize );

		virtual std::string get_name() const;

		bool
		sampling_sugar() const;

		void
		prepare_from_prior_sampled_sugar_jobs( pose::Pose const & pose,
																					 utility::vector1< pose::PoseOP > & starting_pose_data_list );

		SugarModeling const &
		anchor_sugar_modeling();

		void
		set_scorefxn( core::scoring::ScoreFunctionOP const & scorefxn );

		void set_do_not_sample_multiple_virtual_sugar( bool const & setting ){ do_not_sample_multiple_virtual_sugar_ = setting; }

		void set_sample_ONLY_multiple_virtual_sugar( bool const & setting ){ sample_ONLY_multiple_virtual_sugar_ = setting; }

		void set_assert_no_virt_sugar_sampling( bool const & setting ){ assert_no_virt_sugar_sampling_ = setting; }

		void
		set_integration_test_mode( bool const & setting ){ integration_test_mode_ = setting; }

		void
		set_use_phenix_geo( bool const & setting ) { use_phenix_geo_ = setting;	}

		void
		set_legacy_mode( bool const & setting ) { legacy_mode_ = setting;	}

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
		bool integration_test_mode_;
		bool use_phenix_geo_;
		bool legacy_mode_;

		bool do_not_sample_multiple_virtual_sugar_;
		bool sample_ONLY_multiple_virtual_sugar_;
		bool assert_no_virt_sugar_sampling_;

		core::scoring::ScoreFunctionOP scorefxn_;

		std::map< Size, Size > reference_res_for_each_virtual_sugar_;
	};

} //rna
} //swa
} //protocols

#endif
