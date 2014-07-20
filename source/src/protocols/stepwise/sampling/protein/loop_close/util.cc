// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampling/protein/loop_close/util.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/sampling/protein/loop_close/util.hh>
#include <protocols/stepwise/sampling/working_parameters/StepWiseWorkingParameters.hh>
#include <protocols/stepwise/sampling/protein/loop_close/StepWiseProteinKIC_LoopBridger.hh>
#include <protocols/stepwise/sampling/modeler_options/StepWiseModelerOptions.hh>
#include <protocols/stepwise/sampling/protein/StepWiseProteinBackboneSampler.hh>
#include <protocols/rotamer_sampler/RotamerSamplerSizedComb.hh>
#include <protocols/rotamer_sampler/protein/ProteinMainChainRotamerSampler.hh>
#include <core/pose/Pose.hh>
#include <core/id/TorsionID.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.stepwise.sampling.protein.loop_close.util" );

namespace protocols {
namespace stepwise {
namespace sampling {
namespace protein {
namespace loop_close {

	///////////////////////////////////////////////////////////////////////////////////////////////////
	// this is largely untested now -- if we want to reconstitute, we should
	// probably fold KIC into rotamer_sampler, as is done for RNA.
	void
	kic_close_loops_in_samples( rotamer_sampler::RotamerSamplerSizedOP & sampler,
															core::pose::Pose const & pose,
															working_parameters::StepWiseWorkingParametersCOP working_parameters_,
															modeler_options::StepWiseModelerOptionsCOP options_ ){

		using namespace protocols::rotamer_sampler::protein;
		using namespace loop_close;

		// Kinematic Inversion (a.k.a., 'analytical' or 'triaxial') loop closure
		// Involves 5 residues: takeoff - bridge_res1 - bridge_res2 -bridge_res3 - landing
		runtime_assert( options_->kic_sampling_if_relevant() );
		runtime_assert( working_parameters_->working_bridge_res().size() == 3 );

		if ( options_->dump() ) pose.dump_pdb("before_loop_close.pdb");

		// close loops here, figuring out what torsions are compatible with input poses.
		// Will go through all combinations of poses in the sample generator, and close loops.
		utility::vector1< id::TorsionID > which_torsions;  // will include entire loop.
		utility::vector1<  utility::vector1< Real > > main_chain_torsion_set_lists; // torsions that correspond to closed loops.

		// sample N-terminal-phi of takeoff and C-terminal-psi of landing.
		core::pose::Pose kic_sampler_pose = pose;
		if  ( !options_->disable_sampling_of_loop_takeoff() ) enable_sampling_of_loop_takeoff( sampler, kic_sampler_pose, working_parameters_, options_ );

		// stepwiseproteinloopbridger figures out loop residues as those positions that are not 'fixed'.
		StepWiseProteinKIC_LoopBridger stepwise_loop_bridger( sampler, working_parameters_ );
		stepwise_loop_bridger.apply( kic_sampler_pose );

		which_torsions = stepwise_loop_bridger.which_torsions();
		main_chain_torsion_set_lists = stepwise_loop_bridger.main_chain_torsion_set_lists();

		TR << "Using ProteinMainChainRotamerSampler with loops" << std::endl;
		ProteinMainChainRotamerSamplerOP new_sampler = new ProteinMainChainRotamerSampler( which_torsions,
																																				 main_chain_torsion_set_lists,
																																				 options_->choose_random() );
		new_sampler->init();
		sampler = new_sampler;
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	// This is meant for KIC-loop sampling cases in which we get a pose and just want to do closure across
	// three residues with minimal sampling. It still makes sense to sample 'takeoff' phi and psi into the
	// three-residue loop.
	void
	enable_sampling_of_loop_takeoff( rotamer_sampler::RotamerSamplerSizedOP & sampler,
																	 pose::Pose & pose,
																	 working_parameters::StepWiseWorkingParametersCOP working_parameters_,
																	 modeler_options::StepWiseModelerOptionsCOP options_ ) {

		using namespace protocols::stepwise;
		using namespace protocols::stepwise::sampling::protein;
		using namespace protocols::rotamer_sampler::protein;

		// this is to prevent confusion with sampling of loop takeoff that would occur above if moving_res_list
		// is non-empty and StepWiseProteinBackboneSampler is in use.
		runtime_assert( working_parameters_->working_moving_res_list().size() == 0 ); // may cause fail.
		runtime_assert( working_parameters_->working_bridge_res().size() == 3 );

		StepWiseProteinBackboneSampler backbone_sampler( working_parameters_ );
		backbone_sampler.set_n_sample( options_->n_sample() );

		utility::vector1< Size > takeoff_res;
		Size const pre_loop_res  =  working_parameters_->working_bridge_res()[1] - 1;
		Size const post_loop_res =  working_parameters_->working_bridge_res()[3] + 1;
		takeoff_res.push_back( pre_loop_res  );
		takeoff_res.push_back( post_loop_res );
		backbone_sampler.set_moving_residues( takeoff_res );

		utility::vector1< Size > fixed_res_for_backbone_sampler = takeoff_res;
		if ( pre_loop_res  > 1                    ) fixed_res_for_backbone_sampler.push_back( pre_loop_res  - 1 );
		if ( post_loop_res < pose.total_residue() ) fixed_res_for_backbone_sampler.push_back( post_loop_res + 1 );
		backbone_sampler.set_fixed_residues( fixed_res_for_backbone_sampler );
		backbone_sampler.apply( pose );

		TR << "Going to sample this many takeoff psi/phi combinations: " <<
			backbone_sampler.main_chain_torsion_set_lists_real().size() << std::endl;

		ProteinMainChainRotamerSamplerOP sampler_for_takeoff_res =
          new ProteinMainChainRotamerSampler( backbone_sampler.which_torsions(),
																			 backbone_sampler.main_chain_torsion_set_lists_real(),
																			 false /*choose_random*/ );
		sampler = new rotamer_sampler::RotamerSamplerSizedComb( sampler /*input pose sample generator from above*/, sampler_for_takeoff_res /*inner loop*/);

	}

} //loop_close
} //protein
} //sampling
} //stepwise
} //protocols
