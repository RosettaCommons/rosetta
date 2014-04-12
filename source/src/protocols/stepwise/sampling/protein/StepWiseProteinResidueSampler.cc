// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampling/protein/StepWiseProteinResidueSampler.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/sampling/protein/StepWiseProteinResidueSampler.hh>
#include <protocols/stepwise/sampling/protein/StepWiseProteinModelerOptions.hh>
#include <protocols/stepwise/sampling/protein/StepWiseProteinJobParameters.hh>
#include <protocols/stepwise/sampling/protein/StepWiseProteinBackboneSampler.hh>
#include <protocols/stepwise/sampling/protein/StepWiseProteinPacker.hh>
#include <protocols/stepwise/sampling/protein/MainChainTorsionSet.hh>
#include <protocols/stepwise/sampling/protein/StepWiseProteinUtil.hh>
#include <protocols/stepwise/sampling/protein/checker/ProteinAtrRepChecker.hh>
#include <protocols/stepwise/sampling/general/StepWiseClusterer.hh>
#include <protocols/stepwise/sampling/protein/loop_close/StepWiseProteinCCD_Closer.hh>
#include <protocols/stepwise/sampling/protein/loop_close/StepWiseProteinKIC_LoopBridger.hh>
#include <protocols/stepwise/screener/SampleApplier.hh>
#include <protocols/stepwise/screener/SimpleRMSD_Screener.hh>
#include <protocols/stepwise/screener/ProteinAtrRepScreener.hh>
#include <protocols/stepwise/screener/ProteinPackScreener.hh>
#include <protocols/stepwise/screener/SimplePoseSelection.hh>
#include <protocols/stepwise/StepWiseSampleAndScreen.hh>
#include <protocols/rotamer_sampler/NoOpRotamer.hh>
#include <protocols/rotamer_sampler/RotamerSizedComb.hh>
#include <protocols/rotamer_sampler/input_streams/InputStreamRotamer.hh>
#include <protocols/rotamer_sampler/protein/ProteinMainChainRotamer.hh>
#include <protocols/rotamer_sampler/protein/ProteinBetaAntiParallelRotamer.hh>
#include <protocols/rotamer_sampler/protein/ProteinFragmentRotamer.hh>
#include <protocols/stepwise/StepWiseUtil.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pose/Pose.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.stepwise.sampling.protein.StepWiseProteinResidueSampler" );

namespace protocols {
namespace stepwise {
namespace sampling {
namespace protein {

	//Constructor
	StepWiseProteinResidueSampler::StepWiseProteinResidueSampler( StepWiseProteinJobParametersCOP & job_parameters ):
		job_parameters_( job_parameters ),
		close_loops_( false )
	{}

	//Destructor
	StepWiseProteinResidueSampler::~StepWiseProteinResidueSampler()
	{}

	///////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinResidueSampler::apply( pose::Pose & pose ){

		initialize_sampler( pose );
		initialize_screeners( pose );

		StepWiseSampleAndScreen sample_and_screen( sampler_, screeners_ );
		sample_and_screen.set_max_ntries( 500 /*get_max_ntries()*/ );
		sample_and_screen.set_num_random_samples( options_->num_random_samples() );

		sample_and_screen.run();
		//sample_and_screen.output_counts();

		pose_selection_->finalize(); /* does clustering of all collected poses */
		pose_list_ = pose_selection_->pose_list();
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinResidueSampler::initialize_sampler( pose::Pose & pose ){
		sampler_ = get_basic_sampler( pose );
		if ( moving_res_list_.size() == 0 ) return;

		// loop closure stuff. Later will actually move KIC into sampler, and CCD into screener.
		close_loops_ =  (options_->ccd_close() || do_ccd_ || job_parameters_->working_bridge_res().size() > 0 );
		if ( close_loops_ ) sampler_ = close_loops_in_samples( pose, sampler_ );
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinResidueSampler::initialize_screeners( pose::Pose & pose ){

		using namespace screener;

		utility::vector1< Size > const moving_partition_res = figure_out_moving_partition_res( pose, moving_res_list_ );

		screeners_.push_back( new SampleApplier( pose ) );

		// later can put CCD closure in here as a screener.

		if ( options_->rmsd_screen() > 0.0 ){
			runtime_assert( job_parameters_->working_native_pose() );
			screeners_.push_back( new SimpleRMSD_Screener( pose, moving_partition_res, job_parameters_->working_native_pose(),
																										 options_->rmsd_screen() ) ); // [ replacement for PoseFilter_RMSD_Screen. ]
		}

		StepWiseProteinPackerOP stepwise_packer = get_packer( pose );
	 	if ( !full_optimize_ && options_->prepack() ) {
			stepwise_packer->set_moving_res_list( moving_res_list_ );
			stepwise_packer->do_prepack( pose );
		}
		stepwise_packer->initialize( pose );

  	if ( !full_optimize_ && !close_loops_ && options_->atr_rep_screen() ) {
			screening_pose_ = pose.clone();
			screeners_.push_back( new ProteinAtrRepScreener( *screening_pose_,
																											 new checker::ProteinAtrRepChecker( pose, moving_res_list_ ) ) );
		}

		screeners_.push_back( new ProteinPackScreener( pose, stepwise_packer ) );

		pose_selection_ = new SimplePoseSelection( pose, moving_res_list_, options_, full_optimize_ );
		screeners_.push_back( pose_selection_ );
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////
	rotamer_sampler::RotamerSizedOP
	StepWiseProteinResidueSampler::get_basic_sampler( pose::Pose & pose ){

		using namespace protocols::stepwise::sampling::protein;
		using namespace protocols::rotamer_sampler;
		using namespace protocols::rotamer_sampler::input_streams;
		using namespace protocols::rotamer_sampler::protein;

		RotamerSizedOP sampler;
		full_optimize_ = options_->global_optimize();

		if ( frag_files_.size() > 0 ){
			std::string const frag_file  = frag_files_[ 1 ];
			utility::vector1< Size > const & slice_res = job_parameters_->working_res_list();
			sampler = new ProteinFragmentRotamer( frag_file, slice_res, moving_res_list_ );
			if ( input_streams_.size() == 1 ) {
				RotamerSizedOP sampler_identity = new InputStreamRotamer( input_streams_[1] );
				sampler = new RotamerSizedComb( sampler_identity /*outer*/, sampler /*inner, fragment generator above*/);
			}
		} else if ( input_streams_.size() == 2 ){
			// assume that we want to "combine" two streams of poses...
			// This would be the mode if we have a bunch of templates from which we will graft chunks.
			// Or if we have SWA-based little fragments that we want to paste in.
			//			runtime_assert( stepwise_pose_setup_ != 0 );
			RotamerSizedOP input_stream_sampler1 = new InputStreamRotamer( input_streams_[1] );
			RotamerSizedOP input_stream_sampler2 = new InputStreamRotamer( input_streams_[2] );
			sampler= new RotamerSizedComb( input_stream_sampler1 /*outer*/, input_stream_sampler2 /*inner*/);
		} else	if ( options_->sample_beta() ) {

			if ( moving_res_list_.size() !=  1 ) utility_exit_with_message( "Sample beta only works for adding one residue to a beta sheet...");
			sampler =	new ProteinBetaAntiParallelRotamer( pose, moving_res_list_[1] );
		} else if ( moving_res_list_.size() > 0 && !options_->disallow_backbone_sampling() ){

			//////////////////////////////////////////////////////////////////////
			//  DEFAULT -- enumeration of conformations for moving residues.
			// Screen that predefines (phi, psi, omega)  for moving residues
			// --> input for later mover carries out all the green packer moves.
			//////////////////////////////////////////////////////////////////////
			StepWiseProteinBackboneSampler backbone_sampler( job_parameters_ );
			backbone_sampler.set_n_sample( options_->n_sample() );
			backbone_sampler.set_native_pose( job_parameters_->working_native_pose() );

			//when called in StepWiseMonteCarlo -- had a different convention in original StepWiseAssembly,
			// perhaps should be deprecated. If moving_res_list is 3-4, following samples psi/omega at 2 and phi at 5:
			backbone_sampler.set_expand_loop_takeoff( options_->expand_loop_takeoff() );

			///////////////////////////////////////////////////////////////
			// Following has not been tested in a while, may not work:
			backbone_sampler.set_filter_native_big_bins( options_->filter_native_big_bins()  );
			std::string const silent_file_centroid = get_file_name( options_->silent_file(), "_centroid" );
			if ( options_->centroid_output() ) backbone_sampler.set_silent_file( silent_file_centroid );
			if ( options_->centroid_screen() ) {
				backbone_sampler.setup_centroid_screen( options_->centroid_score_diff_cut(), options_->centroid_weights(),
																								 options_->nstruct_centroid(), options_->ghost_loops() );
			}

			// Could also put loose chainbreak closure check here.
			backbone_sampler.apply( pose );
			sampler =  new ProteinMainChainRotamer( backbone_sampler.which_torsions(),
																							 backbone_sampler.main_chain_torsion_set_lists_real(),
																							 options_->choose_random() );
			TR << "Using ProteinMainChainRotamer. Num poses: " << backbone_sampler.main_chain_torsion_set_lists_real().size() << std::endl;
	} else if ( input_streams_.size() > 0 ){ // Just give the poses one at a time... useful in prepacks.
			sampler = new InputStreamRotamer( input_streams_[1] );
			full_optimize_ = true;
		} else {
			sampler = new NoOpRotamer();
			full_optimize_ = true;
		}

		if ( full_optimize_ ) {
			// This should force clustering, minimizing based on RMSD over all residues.
			moving_res_list_.clear();
			for ( Size i = 1; i <= pose.total_residue(); i++ ) moving_res_list_.push_back( i );
		}

		sampler->init();
		return sampler;
	}


	///////////////////////////////////////////////////////////////////////////////////////////////////
	rotamer_sampler::RotamerSizedOP
	StepWiseProteinResidueSampler::close_loops_in_samples( core::pose::Pose & pose,
																								  rotamer_sampler::RotamerSizedOP sampler ){

		using namespace protocols::rotamer_sampler::protein;

		runtime_assert( close_loops_ );
		// close loops here, figuring out what torsions are compatible with input poses.
		// Will go through all combinations of poses in the sample generator, and close loops.
		utility::vector1< id::TorsionID > which_torsions;  // will include entire loop.
		utility::vector1<  utility::vector1< Real > > main_chain_torsion_set_lists; // torsions that correspond to closed loops.

		if ( options_->dump() ) pose.dump_pdb("before_loop_close.pdb");
		if ( do_ccd_ || options_->ccd_close() ) { // CCD, as opposed to KIC

			// CCD closure -- heuristic closer but will accept fewer than 6 torsions (and does not require
			//  the torsions to be triaxial like kinematic closer).
			// Involves up to 5 residues: takeoff - bridge_res1 - bridge_res2 -bridge_res3 - landing
			if ( job_parameters_->working_bridge_res().size() >= 3 ) utility_exit_with_message( "cannot specify more than 2 bridge_res for CCD loop closure" );
			StepWiseProteinCCD_Closer stepwise_ccd_closer( sampler, job_parameters_ );
			stepwise_ccd_closer.set_choose_random( options_->choose_random() );
			stepwise_ccd_closer.set_num_random_samples( options_->num_random_samples() );
			stepwise_ccd_closer.set_working_moving_res_list( moving_res_list_ );
			stepwise_ccd_closer.apply( pose );

			which_torsions = stepwise_ccd_closer.which_torsions();
			main_chain_torsion_set_lists = stepwise_ccd_closer.main_chain_torsion_set_lists();

		} else {
			// Kinematic Inversion (a.k.a., 'analytical' or 'triaxial') loop closure
			// Involves 5 residues: takeoff - bridge_res1 - bridge_res2 -bridge_res3 - landing
			if ( job_parameters_->working_bridge_res().size() != 3 ) utility_exit_with_message( "must specify exactly 3 bridge_res for kinematic loop closure" );

			// sample N-terminal-phi of takeoff and C-terminal-psi of landing.
			if  ( !options_->disable_sampling_of_loop_takeoff() ) enable_sampling_of_loop_takeoff( pose, sampler );

			// stepwiseproteinloopbridger figures out loop residues as those positions that are not 'fixed'.
			StepWiseProteinLoopBridger stepwise_loop_bridger( sampler, job_parameters_ );
			stepwise_loop_bridger.apply( pose );

			which_torsions = stepwise_loop_bridger.which_torsions();
			main_chain_torsion_set_lists = stepwise_loop_bridger.main_chain_torsion_set_lists();

		}

		moving_res_list_ = merge_vectors( moving_res_list_, job_parameters_->working_bridge_res() );
		TR << "Using ProteinMainChainRotamer with loops" << std::endl;
		ProteinMainChainRotamerOP new_sampler = new ProteinMainChainRotamer( which_torsions,
																																				 main_chain_torsion_set_lists,
																																				 options_->choose_random() );
		new_sampler->init();
		return new_sampler;
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	// This is meant for KIC-loop sampling cases in which we get a pose and just want to do closure across
	// three residues with minimal sampling. It still makes sense to sample 'takeoff' phi and psi into the
	// three-residue loop.
	void
	StepWiseProteinResidueSampler::enable_sampling_of_loop_takeoff( pose::Pose & pose,
																													 rotamer_sampler::RotamerSizedOP sampler ) {

		using namespace protocols::stepwise;
		using namespace protocols::stepwise::sampling::protein;
		using namespace protocols::rotamer_sampler::protein;

		// this is to prevent confusion with sampling of loop takeoff that would occur above if moving_res_list
		// is non-empty and StepWiseProteinBackboneSampler is in use.
		runtime_assert( moving_res_list_.size() == 0 );
		runtime_assert( job_parameters_->working_bridge_res().size() == 3 );

		StepWiseProteinBackboneSampler backbone_sampler( job_parameters_ );
		backbone_sampler.set_n_sample( options_->n_sample() );

		utility::vector1< Size > takeoff_res;
		Size const pre_loop_res  =  job_parameters_->working_bridge_res()[1] - 1;
		Size const post_loop_res =  job_parameters_->working_bridge_res()[3] + 1;
		takeoff_res.push_back( pre_loop_res  );
		takeoff_res.push_back( post_loop_res );
		backbone_sampler.set_moving_residues( takeoff_res );

		utility::vector1< Size > fixed_res_for_backbone_sampler = takeoff_res;
		if ( pre_loop_res  > 1                    ) fixed_res_for_backbone_sampler.push_back( pre_loop_res  - 1 );
		if ( post_loop_res < pose.total_residue() ) fixed_res_for_backbone_sampler.push_back( post_loop_res + 1 );
		backbone_sampler.set_fixed_residues( fixed_res_for_backbone_sampler ); //
		backbone_sampler.apply( pose );

		TR << "Going to sample this many takeoff psi/phi combinations: " <<
			backbone_sampler.main_chain_torsion_set_lists_real().size() << std::endl;

		ProteinMainChainRotamerOP sampler_for_takeoff_res =
          new ProteinMainChainRotamer( backbone_sampler.which_torsions(),
																			 backbone_sampler.main_chain_torsion_set_lists_real(),
																			 false /*choose_random*/ );
		sampler = new rotamer_sampler::RotamerSizedComb( sampler /*input pose sample generator from above*/, sampler_for_takeoff_res /*inner loop*/);

	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	StepWiseProteinPackerOP
	StepWiseProteinResidueSampler::get_packer( pose::Pose & pose ) {

		using namespace core::scoring;

		if ( !scorefxn_ || scorefxn_->get_name() != options_->pack_weights() ) scorefxn_ = ScoreFunctionFactory::create_score_function( options_->pack_weights() );
		check_scorefxn_has_constraint_terms_if_pose_has_constraints( pose, scorefxn_ );
		scorefxn_->set_weight( linear_chainbreak, 0.2 /*arbitrary*/ );
		if ( options_->mapfile_activated() && scorefxn_->get_weight( elec_dens_atomwise ) == 0.0 ) scorefxn_->set_weight( elec_dens_atomwise, 10.0 );
		// may want to put in proper handling of moving_partition_res -- i.e., pack any residues that make new contacts:
		StepWiseProteinPackerOP stepwise_packer = new StepWiseProteinPacker( moving_res_list_ );
		stepwise_packer->set_native_pose( job_parameters_->working_native_pose() );
		stepwise_packer->set_scorefxn( scorefxn_ );
		stepwise_packer->set_use_green_packer( options_->use_green_packer() );
		stepwise_packer->set_use_packer_instead_of_rotamer_trials( options_->use_packer_instead_of_rotamer_trials() );
		stepwise_packer->set_rescore_only( options_->rescore_only() );
		stepwise_packer->set_allow_virtual_side_chains( options_->allow_virtual_side_chains() );

		return stepwise_packer;
	}

} //protein
} //sampling
} //stepwise
} //protocols
