// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/modeler/rna/sugar/VirtualSugarJustInTimeInstantiator.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_rna_VirtualSugarJustInTimeInstantiator_HH
#define INCLUDED_protocols_stepwise_rna_VirtualSugarJustInTimeInstantiator_HH

#include <protocols/moves/Mover.hh>
#include <protocols/stepwise/modeler/rna/sugar/VirtualSugarJustInTimeInstantiator.fwd.hh>
#include <protocols/stepwise/modeler/options/StepWiseModelerOptions.fwd.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.fwd.hh>
#include <protocols/stepwise/modeler/rna/sugar/SugarModeling.fwd.hh>
#include <protocols/stepwise/sampler/copy_dofs/ResidueAlternativeSet.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/types.hh>

using namespace core;

namespace protocols {
namespace stepwise {
namespace modeler {
namespace rna {
namespace sugar {

	class VirtualSugarJustInTimeInstantiator: public protocols::moves::Mover {

	public:

		//constructor
		VirtualSugarJustInTimeInstantiator( working_parameters::StepWiseWorkingParametersCOP & working_parameters_ );

		//destructor
		~VirtualSugarJustInTimeInstantiator();

		virtual void apply( pose::Pose & pose_to_visualize );

		virtual std::string get_name() const;

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		bool
		do_the_modeler( core::pose::Pose & pose );

		bool
		modeler_sugar() const;

		bool
		modeler_sugar_at_chain_break() const;

		void
		prepare_from_prior_sampled_sugar_jobs( pose::Pose const & pose,
																					 utility::vector1< pose::PoseOP > & starting_pose_data_list,
																					 bool const pose_explosion_legacy /* = false */ ) const;

		void
		prepare_from_prior_sampled_sugar_jobs_for_chain_break( pose::Pose const & pose,
																													 utility::vector1< pose::PoseOP > & starting_pose_data_list ) const;

		void
		set_scorefxn( core::scoring::ScoreFunctionCOP scorefxn );

		void
		set_keep_base_fixed( bool const & setting ) { keep_base_fixed_ = setting;	}

		bool const & success() const{ return success_; }

		void
		set_options( options::StepWiseModelerOptionsCOP options );

		Size
		num_sets() const { return sugar_modeling_sets_.size(); }

		sampler::copy_dofs::ResidueAlternativeSet const &
		residue_alternative_set( Size const n );

		//legacy:
		SugarModeling const &
		anchor_sugar_modeling();

		void
		instantiate_sugars_at_cutpoint_closed( pose::Pose & pose ) const;

	private:

		Size
		sampled_sugar_index( Size const i );

		bool
		do_sugar_sampling( pose::Pose & viewer_pose, SugarModeling & sugar_modeling, std::string const name );

		bool
		setup_sugar_modeling( pose::Pose const & pose, Size const moving_res, SugarModeling & sugar_modeling );

		bool
		get_sugar_modeling_set( pose::Pose & viewer_pose, Size const i );

		utility::vector1< SugarModelingOP >	get_sugar_modeling_sets_for_chainbreak() const;

		void
		instantiate_sugar( pose::Pose & pose, SugarModeling const & sugar_modeling, Size const sugar_ID ) const;

		void
		instantiate_sugars_recursively(  pose::Pose const & pose,
																		 utility::vector1< pose::PoseOP > & pose_data_list,
																		 utility::vector1< SugarModelingOP > const & sugar_modeling_sets,
																		 utility::vector1< Size > const & sugar_modeling_set_indices ) const;

		void
		minimize_sugar_sets_legacy( pose::Pose const & pose,
																utility::vector1< pose::PoseOP > & pose_data_list ) const;


	private:

		working_parameters::StepWiseWorkingParametersCOP working_parameters_;
		Size const moving_res_;
		bool const rebuild_bulge_mode_;

		utility::vector1< SugarModelingOP > sugar_modeling_sets_;
		utility::vector1< sampler::copy_dofs::ResidueAlternativeSetOP> residue_alternative_sets_;
		bool keep_base_fixed_;
		bool const moving_res_legacy_;

		options::StepWiseModelerOptionsCOP options_;
		core::scoring::ScoreFunctionCOP scorefxn_;

		bool success_;
	};

} //sugar
} //rna
} //modeler
} //stepwise
} //protocols

#endif
