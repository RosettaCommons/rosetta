// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampling/protein/StepWiseProteinModeler.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_sampling_protein_StepWiseProteinModeler_HH
#define INCLUDED_protocols_stepwise_sampling_protein_StepWiseProteinModeler_HH

#include <protocols/moves/Mover.hh>
#include <protocols/stepwise/sampling/protein/StepWiseProteinModeler.fwd.hh>
#include <protocols/stepwise/sampling/protein/StepWiseProteinModelerOptions.fwd.hh>
#include <protocols/stepwise/sampling/protein/StepWiseProteinJobParameters.fwd.hh>
#include <protocols/stepwise/sampling/protein/InputStreamWithResidueInfo.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pose/Pose.fwd.hh>

////////////////////////////////////////////////////////////////////////
// This copies code in StepWiseRNA_Modeler, which is not good.
// The plan is eventually to unify the two into StepWiseModeler (which
//  is currently a wrapper).
// For now, writing this separately to pilot use of protein SWA functions
//  into a StepWise MonteCarlo; after this, we can identify commonalities with
//  RNA and integrate.
//

using namespace core;

namespace protocols {
namespace stepwise {
namespace sampling {
namespace protein {

	class StepWiseProteinModeler: public protocols::moves::Mover  {

	public:

		//constructor
		StepWiseProteinModeler( scoring::ScoreFunctionOP scorefxn );

		StepWiseProteinModeler( core::scoring::ScoreFunctionOP scorefxn,
														utility::vector1< Size > const & moving_res_list );

		StepWiseProteinModeler( StepWiseProteinModeler const & src );

		//destructor
		~StepWiseProteinModeler();

	public:

		StepWiseProteinModelerOP clone_modeler() const;

		StepWiseProteinModeler & operator=( StepWiseProteinModeler const & src );

		virtual void apply( pose::Pose & pose );

		virtual std::string get_name() const { return "StepWiseProteinModeler"; }

		void set_moving_res_and_reset( Size const moving_res );

		void set_working_minimize_res( utility::vector1< Size > const & setting ){ working_minimize_res_ = setting; }

		void set_skip_sampling( bool const & setting );

		StepWiseProteinModelerOptionsCOP options();

		void
		set_options( StepWiseProteinModelerOptionsCOP options );

		void
		set_frag_files( utility::vector1< std::string > const & frag_files );

		void
		set_input_streams( utility::vector1< InputStreamWithResidueInfoOP > const & input_streams );

		void
		set_job_parameters( StepWiseProteinJobParametersCOP job_parameters );

		void set_moving_res_list( utility::vector1< Size > const & setting ){ moving_res_list_ = setting; }
		void set_bridge_res( utility::vector1< Size > const & setting ){ bridge_res_ = setting; }

	private:

		void
		reinitialize();

		void
		initialize_job_parameters_and_root( pose::Pose & pose );


		void
		do_residue_sampling( core::pose::Pose & pose,
												 utility::vector1< pose::PoseOP > & pose_list );

		void
		do_minimizing( core::pose::Pose & pose,
									 utility::vector1< pose::PoseOP > & pose_list );

		StepWiseProteinJobParametersOP
		setup_job_parameters_for_stepwise_with_full_model_info( core::pose::Pose & pose );

		void figure_out_moving_res_list( pose::Pose const & pose );

		void
		figure_out_protein_moving_res_list_from_most_distal_res( pose::Pose const & pose, Size const moving_res );

		void
		align_pose_and_add_rmsd_constraints( pose::Pose & pose ) const;

		bool
		sampling_successful( utility::vector1< pose::PoseOP > & pose_list );

	private:

		core::scoring::ScoreFunctionOP scorefxn_;
		core::scoring::ScoreFunctionOP pack_scorefxn_;
		utility::vector1< Size > moving_res_list_;

		Size moving_res_;
		bool full_optimize_;
		bool do_ccd_;

		StepWiseProteinModelerOptionsCOP options_;
		StepWiseProteinJobParametersCOP job_parameters_;

		// setup in swa_protein_main, but not implemented/necessary for stepwise monte carlo.
		utility::vector1< std::string > frag_files_;
		utility::vector1< InputStreamWithResidueInfoOP > input_streams_;

		utility::vector1< Size > bridge_res_;
		utility::vector1< Size > working_minimize_res_;
	};

} //protein
} //sampling
} //stepwise
} //protocols

#endif
