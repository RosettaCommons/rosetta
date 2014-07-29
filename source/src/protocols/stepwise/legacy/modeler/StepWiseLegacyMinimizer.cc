// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/modeler/StepWiseMinimizer.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu

#include <protocols/stepwise/legacy/modeler/StepWiseLegacyMinimizer.hh>
#include <protocols/stepwise/legacy/modeler/protein/StepWiseProteinMinimizer.hh>
#include <protocols/stepwise/modeler/protein/util.hh>
#include <protocols/stepwise/legacy/modeler/rna/StepWiseRNA_Minimizer.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_BaseCentroidChecker.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_VDW_BinChecker.hh>
#include <protocols/stepwise/modeler/align/StepWiseLegacyClusterer.hh>
#include <protocols/stepwise/modeler/options/StepWiseModelerOptions.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.hh>
#include <protocols/stepwise/modeler/working_parameters/util.hh>
#include <protocols/stepwise/modeler/movemap/util.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>
#include <utility/stream_util.hh>
#include <utility/tools/make_vector1.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.stepwise.legacy.modeler.StepWiseLegacyMinimizer" );

using namespace core;
using namespace protocols::stepwise::modeler::rna;
using namespace protocols::stepwise::modeler::protein;
using utility::tools::make_vector1;

////////////////////////////////////////////////////////////////////////////////
//
// Deprecated wrapper around StepWiseRNA_Minimizer and StepWiseProteinMinimizer.
//
////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace stepwise {
namespace legacy {
namespace modeler {

	//Constructor
	StepWiseLegacyMinimizer::StepWiseLegacyMinimizer( utility::vector1< pose::PoseOP > & pose_list,
																				stepwise::modeler::working_parameters::StepWiseWorkingParametersCOP working_parameters,
																				options::StepWiseModelerOptionsCOP options,
																				core::scoring::ScoreFunctionCOP scorefxn):
		pose_list_( pose_list ), // where work will be saved.
		working_parameters_( working_parameters ),
		options_( options ),
		scorefxn_( scorefxn ),
		moving_res_list_( working_parameters->working_moving_res() ),
		protein_full_optimize_( false )
	{}

	//Destructor
	StepWiseLegacyMinimizer::~StepWiseLegacyMinimizer()
	{}

	//////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseLegacyMinimizer::apply( pose::Pose & pose ){
		if ( is_protein( pose, moving_res_list_ ) ){
			do_protein_minimizing( pose );
		} else {
			do_rna_minimizing( pose );
		}
	}

	////////////////////////////////////////////////////////////////////
	// deprecate soon...
	////////////////////////////////////////////////////////////////////
	void
	StepWiseLegacyMinimizer::do_protein_minimizing( pose::Pose & pose ) {

		using namespace protocols::stepwise::modeler::align;

		if ( pose_list_.size() == 0 ) pose_list_ = make_vector1( pose.clone() );
		StepWiseProteinMinimizerOP stepwise_protein_minimizer = new StepWiseProteinMinimizer( pose_list_, moving_res_list_ );
		stepwise_protein_minimizer->set_scorefxn( get_minimize_scorefxn( pose, scorefxn_, options_ ) );

		std::string const silent_file_minimize = get_file_name( options_->silent_file(), "_minimize" );
		stepwise_protein_minimizer->set_move_jumps_between_chains( options_->move_jumps_between_chains() );
		stepwise_protein_minimizer->set_native_pose( working_parameters_->working_native_pose() );
		stepwise_protein_minimizer->set_calc_rms_res( working_parameters_->working_calc_rms_res() ); // used for calculating rmsds to native.
		stepwise_protein_minimizer->set_fixed_res( working_parameters_->working_fixed_res() );
		stepwise_protein_minimizer->set_move_takeoff_torsions( !options_->disable_sampling_of_loop_takeoff() );
		stepwise_protein_minimizer->set_rescore_only( options_->rescore_only() );
		stepwise_protein_minimizer->set_cartesian( options_->cart_min() );
		stepwise_protein_minimizer->set_min_type( options_->min_type() );
		stepwise_protein_minimizer->set_min_tolerance( options_->min_tolerance() );
		stepwise_protein_minimizer->set_use_coordinate_constraints( !options_->skip_coord_constraints() );
		stepwise_protein_minimizer->set_num_pose_minimize( options_->num_pose_minimize() );

		if ( !options_->skip_minimize() ){
			stepwise_protein_minimizer->apply( pose );
			StepWiseLegacyClustererOP stepwise_clusterer = new StepWiseLegacyClusterer( stepwise_protein_minimizer->pose_list(),
																																			moving_res_list_, options_,
																																			protein_full_optimize_ /*force_align*/ );
			stepwise_clusterer->cluster();
			pose_list_ = stepwise_clusterer->get_pose_list();
			pose = *pose_list_[ 1 ];
		}

		if ( options_->output_minimized_pose_list() ) output_pose_list( pose_list_, working_parameters_->working_native_pose(),
																																						options_->silent_file(), working_parameters_->working_calc_rms_res() );
	}

	////////////////////////////////////////////////////////////////////
	// deprecate soon...
	////////////////////////////////////////////////////////////////////
	void
	StepWiseLegacyMinimizer::do_rna_minimizing( pose::Pose & pose ) {

		if ( options_->minimize_and_score_native_pose() ) {
			runtime_assert ( get_native_pose() );
			add_to_pose_list( pose_list_, *get_native_pose(), "working_native_pose" );
		}

		rna::StepWiseRNA_MinimizerOP stepwise_rna_minimizer_ = new StepWiseRNA_Minimizer( pose_list_, working_parameters_ );
		stepwise_rna_minimizer_->set_scorefxn( scorefxn_ );
		stepwise_rna_minimizer_->set_silent_file( options_->silent_file() );
		stepwise_rna_minimizer_->set_base_centroid_checker( base_centroid_checker_ );
		stepwise_rna_minimizer_->set_user_input_VDW_bin_checker( user_input_VDW_bin_checker_ );

		options::StepWiseModelerOptionsOP minimizer_options = options_->clone();
		if ( options_->rmsd_screen() > 0.0 )	minimizer_options->set_rmsd_screen( options_->rmsd_screen() + 1.0 );
		stepwise_rna_minimizer_->set_options( minimizer_options );
		stepwise_rna_minimizer_->set_working_extra_minimize_res( full_to_sub( const_full_model_info( pose ).extra_minimize_res(), pose ) );
		stepwise_rna_minimizer_->apply ( pose );

		pose_list_ = stepwise_rna_minimizer_->minimized_pose_list();
	}




} //modeler
} //legacy
} //stepwise
} //protocols
