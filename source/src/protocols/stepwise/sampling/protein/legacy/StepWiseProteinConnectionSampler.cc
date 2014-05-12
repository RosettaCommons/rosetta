// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampling/protein/legacy/StepWiseProteinConnectionSampler.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/sampling/protein/legacy/StepWiseProteinConnectionSampler.hh>
#include <protocols/stepwise/sampling/protein/MainChainTorsionSet.hh>
#include <protocols/stepwise/sampling/protein/loop_close/StepWiseProteinCCD_Closer.hh>
#include <protocols/stepwise/sampling/protein/loop_close/util.hh>
#include <protocols/stepwise/sampling/packer/StepWisePacker.hh>
#include <protocols/stepwise/sampling/packer/util.hh>
#include <protocols/stepwise/sampling/protein/util.hh>
#include <protocols/stepwise/sampling/protein/checker/ProteinAtrRepChecker.hh>
#include <protocols/stepwise/sampling/protein/StepWiseProteinBackboneSampler.hh>
#include <protocols/stepwise/sampling/protein/InputStreamWithResidueInfo.hh>
#include <protocols/stepwise/sampling/modeler_options/StepWiseModelerOptions.hh>
#include <protocols/stepwise/sampling/working_parameters/StepWiseWorkingParameters.hh>
#include <protocols/stepwise/sampling/align/StepWiseLegacyClusterer.hh>
#include <protocols/stepwise/sampling/util.hh>
#include <protocols/stepwise/screener/SampleApplier.hh>
#include <protocols/stepwise/screener/SimpleRMSD_Screener.hh>
#include <protocols/stepwise/screener/ProteinCCD_ClosureScreener.hh>
#include <protocols/stepwise/screener/legacy/ProteinAtrRepScreener.hh>
#include <protocols/stepwise/screener/PackScreener.hh>
#include <protocols/stepwise/screener/SimplePoseSelection.hh>
#include <protocols/stepwise/StepWiseSampleAndScreen.hh>
#include <protocols/rotamer_sampler/RotamerSizedComb.hh>
#include <protocols/rotamer_sampler/protein/ProteinMainChainRotamer.hh>
#include <protocols/rotamer_sampler/protein/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/Pose.hh>

#include <basic/Tracer.hh>

using namespace core;

static basic::Tracer TR( "protocols.stepwise.sampling.protein.StepWiseProteinConnectionSampler" );

namespace protocols {
namespace stepwise {
namespace sampling {
namespace protein {

	//Constructor
	StepWiseProteinConnectionSampler::StepWiseProteinConnectionSampler( working_parameters::StepWiseWorkingParametersCOP & working_parameters ):
		working_parameters_( working_parameters ),
		moving_res_list_( working_parameters->working_moving_res_list() ),
		working_obligate_pack_res_( initialize_working_obligate_pack_res( working_parameters ) ),
		pack_all_side_chains_( false )
	{}

	//Destructor
	StepWiseProteinConnectionSampler::~StepWiseProteinConnectionSampler()
	{}

	///////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinConnectionSampler::apply( pose::Pose & pose ){

		initialize_poses_and_checkers( pose );
		initialize_sampler( pose );
		initialize_screeners( pose );

		StepWiseSampleAndScreen sample_and_screen( sampler_, screeners_ );
		sample_and_screen.set_max_ntries( 500 /*get_max_ntries()*/ );
		sample_and_screen.set_num_random_samples( options_->num_random_samples() );

		sample_and_screen.run();
		sample_and_screen.output_counts();

		pose_selection_->finalize(); /* does clustering of all collected poses */
		pose_list_ = pose_selection_->pose_list();

	}


	///////////////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinConnectionSampler::initialize_poses_and_checkers( pose::Pose & pose ){
		using namespace packer;

		// loop closure stuff. Later will actually move KIC into sampler. CCD has been moved into screener.
		utility::vector1< Size > all_moving_res = get_moving_res_including_domain_boundaries( pose, moving_res_list_ );
		protein_cutpoints_closed_ = protein::just_protein( figure_out_moving_cutpoints_closed_from_moving_res( pose, all_moving_res ), pose );

		if ( !options_->kic_sampling_if_relevant() ){
			// CCD closure -- heuristic closer but will accept fewer than 6 torsions (and does not require
			//  the torsions to be triaxial like kinematic closer).
			// Involves up to 5 residues: takeoff - bridge_res1 - bridge_res2 -bridge_res3 - landing
			for ( Size n = 1; n <= protein_cutpoints_closed_.size(); n++ ){
				loop_close::StepWiseProteinCCD_CloserOP stepwise_ccd_closer = new loop_close::StepWiseProteinCCD_Closer( working_parameters_ );
				stepwise_ccd_closer->set_ccd_close_res( protein_cutpoints_closed_[ n ] );
				stepwise_ccd_closer->set_working_moving_res_list( moving_res_list_ );
				stepwise_ccd_closers_.push_back( stepwise_ccd_closer );
				ccd_poses_.push_back( pose.clone() ); // does this really need to be an independent pose?
			}
		}

		pack_all_side_chains_ = ( options_->global_optimize() || ( moving_res_list_.size() == 0 ) );
		if ( pack_all_side_chains_ ){
			working_obligate_pack_res_.clear();
			for ( Size n = 1; n <= pose.total_residue(); n++ ) working_obligate_pack_res_.push_back( n );
		}

		//		TR << "WORKING OBLIGATE PACK RES " << working_obligate_pack_res_ << std::endl;
		stepwise_packer_ = get_packer( scorefxn_, working_obligate_pack_res_, options_ );
	 	if ( !pack_all_side_chains_ && options_->prepack() ) stepwise_packer_->do_prepack( pose, moving_res_list_, working_obligate_pack_res_ );
		stepwise_packer_->initialize( pose );

	}

	///////////////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinConnectionSampler::initialize_sampler( pose::Pose & pose ){
		using namespace protocols::rotamer_sampler::protein;
		sampler_ = get_basic_protein_sampler( pose, moving_res_list_, skip_sampling_,
																					working_parameters_, options_, input_streams_ );

		if ( protein_cutpoints_closed_.size() > 0 && options_->kic_sampling_if_relevant() && !skip_sampling_ ) {
			runtime_assert( protein_cutpoints_closed_.size() == 1 );
			loop_close::kic_close_loops_in_samples( sampler_, pose, working_parameters_, options_ );
		}
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinConnectionSampler::initialize_screeners( pose::Pose & pose ){
		using namespace screener;

		for ( Size n = 1; n <= stepwise_ccd_closers_.size(); n++ ) 	screeners_.push_back( new ProteinCCD_ClosureScreener( stepwise_ccd_closers_[n], *ccd_poses_[n] ) );

		screeners_.push_back( new SampleApplier( pose ) );

		// shouldn't following alread be set up in working_parameters_->moving_partition_res?
		utility::vector1< Size > const moving_partition_res = figure_out_moving_partition_res( pose, moving_res_list_ );
		if ( options_->rmsd_screen() > 0.0 ){
			runtime_assert( working_parameters_->working_native_pose() );
			screeners_.push_back( new SimpleRMSD_Screener( pose, moving_partition_res, working_parameters_->working_native_pose(),
																										 options_->rmsd_screen() ) );
		}

  	if ( !pack_all_side_chains_ && (protein_cutpoints_closed_.size() == 0)  && options_->atr_rep_screen() ) {
			atr_rep_screening_pose_ = pose.clone(); // no difference in variants actually. Perhaps don't even need this.
			screeners_.push_back( new ProteinAtrRepScreener( *atr_rep_screening_pose_,
																											 new checker::ProteinAtrRepChecker( pose, moving_res_list_ ) ) );
		}
		screeners_.push_back( new PackScreener( pose, stepwise_packer_ ) );

		pose_selection_ = new SimplePoseSelection( pose, moving_res_list_, options_, pack_all_side_chains_ );
		screeners_.push_back( pose_selection_ );
	}

} //protein
} //sampling
} //stepwise
} //protocols
