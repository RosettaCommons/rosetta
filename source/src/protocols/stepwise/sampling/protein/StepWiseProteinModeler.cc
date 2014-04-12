// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampling/protein/StepWiseProteinModeler.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/sampling/protein/StepWiseProteinModeler.hh>
#include <protocols/stepwise/sampling/protein/StepWiseProteinModelerOptions.hh>
#include <protocols/stepwise/sampling/protein/StepWiseProteinJobParameters.hh>
#include <protocols/stepwise/sampling/protein/StepWiseProteinJobParametersUtil.hh>
#include <protocols/stepwise/sampling/protein/StepWiseProteinPoseMinimizer.hh>
#include <protocols/stepwise/sampling/protein/StepWiseProteinResidueSampler.hh>
#include <protocols/stepwise/sampling/protein/StepWiseProteinUtil.hh>
#include <protocols/stepwise/sampling/protein/InputStreamWithResidueInfo.hh>
#include <protocols/stepwise/sampling/general/StepWiseClusterer.hh>
#include <protocols/stepwise/sampling/general/StepWisePoseAligner.hh>
#include <protocols/stepwise/monte_carlo/StepWiseMonteCarloUtil.hh>
#include <protocols/stepwise/StepWiseUtil.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/full_model_info/FullModelInfoUtil.hh>
#include <utility/stream_util.hh>
#include <utility/tools/make_vector1.hh>

#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.stepwise.sampling.protein.StepWiseProteinModeler" );

using namespace core;
using utility::tools::make_vector1;
using utility::operator<<;

/////////////////////////////////////////////////////////////////////////////////////////////
//
// Mimics StepWiseRNA_Modeler -- intermediate wrapper between StepWiseModeler and
//   ResidueSampler and Minimizer. Will soon be unified with StepWiseRNA_Modeler into
//   StepWiseModeler.
//
//                     -- rhiju, 2014
//
/////////////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace stepwise {
namespace sampling {
namespace protein {

	//Constructor
	StepWiseProteinModeler::StepWiseProteinModeler( core::scoring::ScoreFunctionOP scorefxn ):
		scorefxn_( scorefxn ),
		moving_res_( 0 ), // will need to be supplied later, and moving_res_list inferred from pose.
		full_optimize_( false ),
		do_ccd_( false )
	{
	}

	//Constructor from swa_protein_main
	StepWiseProteinModeler::StepWiseProteinModeler( core::scoring::ScoreFunctionOP scorefxn,
																									utility::vector1< Size > const & moving_res_list ):
		scorefxn_( scorefxn ),
		moving_res_list_( moving_res_list ),
		moving_res_( 0 ),
		full_optimize_( false ),
		do_ccd_( false )
	{
	}

	StepWiseProteinModeler::StepWiseProteinModeler( StepWiseProteinModeler const & src ):
		Mover( src )
	{
		*this = src;
	}

	//Destructor
	StepWiseProteinModeler::~StepWiseProteinModeler()
	{}

	/// @brief clone the conformation
	StepWiseProteinModelerOP
	StepWiseProteinModeler::clone_modeler() const
	{
		return new StepWiseProteinModeler( *this );
	}

	StepWiseProteinModeler &
	StepWiseProteinModeler::operator=( StepWiseProteinModeler const & src )
	{
		scorefxn_ = src.scorefxn_;
		pack_scorefxn_ = src.pack_scorefxn_;
		moving_res_list_ = src.moving_res_list_;
		moving_res_ = src.moving_res_;
		full_optimize_ = src.full_optimize_;
		do_ccd_ = src.do_ccd_;
		options_ = src.options_;
		job_parameters_ = src.job_parameters_;
		frag_files_ = src.frag_files_;
		input_streams_ = src.input_streams_;
		bridge_res_ = src.bridge_res_;
		working_minimize_res_ = src.working_minimize_res_;

		return *this;
	}

	//////////////////////////////////////////////////////////////////////////////
	// make this closer to StepWiseRNA_Modeler in anticipation of unification.
	void
	StepWiseProteinModeler::apply( pose::Pose & pose ){
		initialize_job_parameters_and_root( pose );

		utility::vector1< PoseOP > pose_list;
		do_residue_sampling( pose, pose_list );
		if ( sampling_successful( pose_list ) ) do_minimizing( pose, pose_list );

		pose.remove_constraints();
		reinitialize();
	}

//////////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinModeler::reinitialize(){
		// Important: make sure that the next time this is used, job parameters is set explicitly -- or it will be reset.
		job_parameters_ = 0;
		moving_res_list_.clear();
		bridge_res_.clear();
		working_minimize_res_.clear();
		do_ccd_ = false;
	}

	//////////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinModeler::initialize_job_parameters_and_root( pose::Pose & pose ){
		pose::full_model_info::make_sure_full_model_info_is_setup( pose );
		if ( job_parameters_ != 0 ) return;
		figure_out_moving_res_list( pose );
		revise_root_and_moving_res_list( pose, moving_res_list_ ); // specify reference_res_? [i.e. anchor_res?]
		align_pose_and_add_rmsd_constraints( pose );
		job_parameters_ = setup_job_parameters_for_stepwise_with_full_model_info( pose );
	}

	////////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinModeler::do_residue_sampling( core::pose::Pose & pose,
																							 utility::vector1< PoseOP > & pose_list ){

		StepWiseProteinResidueSampler stepwise_protein_residue_sampler( job_parameters_ );
		stepwise_protein_residue_sampler.set_options( options_ );
		stepwise_protein_residue_sampler.set_do_ccd( do_ccd_ );
		stepwise_protein_residue_sampler.set_frag_files( frag_files_ );
		stepwise_protein_residue_sampler.set_input_streams( input_streams_ );
		stepwise_protein_residue_sampler.set_scorefxn( pack_scorefxn_ );
		stepwise_protein_residue_sampler.set_moving_res_list( moving_res_list_ );

		stepwise_protein_residue_sampler.apply( pose );

		pose_list = stepwise_protein_residue_sampler.get_pose_list();
		moving_res_list_ = stepwise_protein_residue_sampler.moving_res_list();
		full_optimize_ = stepwise_protein_residue_sampler.full_optimize();
	}


	///////////////////////////////////////////////////////////////////////////////////////////////////
	// Wrapper for routine below (which requires a const pose). If you call this,
	// note that the pose's full_model_info object will be initialized based on its
	// current fold_tree, cutpoint_variants, and any chain/residue-numbering information in
	// PDBInfo.
	StepWiseProteinJobParametersOP
	StepWiseProteinModeler::setup_job_parameters_for_stepwise_with_full_model_info( core::pose::Pose & pose )
	{
		pose::full_model_info::make_sure_full_model_info_is_setup( pose );
		return setup_job_parameters_for_protein_swa( moving_res_list_, pose,
																								 get_native_pose(),
																								 bridge_res_,
																								 working_minimize_res_ /* specify in only minimizing subset of minimizable residues */);
	}

	////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinModeler::do_minimizing( pose::Pose & pose,
																				 utility::vector1< PoseOP > & pose_list ) {

		using namespace core::scoring;
		using namespace protocols::stepwise::sampling::general;

		StepWiseProteinPoseMinimizerOP stepwise_pose_minimizer;
		if ( pose_list.size() == 0 ) pose_list = make_vector1( pose.clone() );
		stepwise_pose_minimizer = new StepWiseProteinPoseMinimizer( pose_list, moving_res_list_ );

		ScoreFunctionOP minimize_scorefxn = scorefxn_->clone();
		if (minimize_scorefxn->get_weight( atom_pair_constraint ) == 0.0) minimize_scorefxn->set_weight( atom_pair_constraint, 1.0 ); //go ahead and turn these on
		if (minimize_scorefxn->get_weight( coordinate_constraint) == 0.0) minimize_scorefxn->set_weight( coordinate_constraint, 1.0 ); // go ahead and turn these on
		check_scorefxn_has_constraint_terms_if_pose_has_constraints( pose, minimize_scorefxn );
		minimize_scorefxn->set_weight( linear_chainbreak, 150.0 );
		if ( options_->cart_min() && ( minimize_scorefxn->get_weight( cart_bonded ) == 0.0 ) ) minimize_scorefxn->set_weight( cart_bonded, 1.0 );
		if ( options_->mapfile_activated() && minimize_scorefxn->get_weight( elec_dens_atomwise ) == 0.0 ) minimize_scorefxn->set_weight( elec_dens_atomwise, 10.0 );

		stepwise_pose_minimizer->set_scorefxn( minimize_scorefxn );

		std::string const silent_file_minimize = get_file_name( options_->silent_file(), "_minimize" );
		stepwise_pose_minimizer->set_move_jumps_between_chains( options_->move_jumps_between_chains() );
		stepwise_pose_minimizer->set_native_pose( job_parameters_->working_native_pose() );
		stepwise_pose_minimizer->set_calc_rms_res( job_parameters_->working_calc_rms_res() ); // used for calculating rmsds to native.
		stepwise_pose_minimizer->set_fixed_res( job_parameters_->working_fixed_res() );
		stepwise_pose_minimizer->set_move_takeoff_torsions( !options_->disable_sampling_of_loop_takeoff() );
		stepwise_pose_minimizer->set_rescore_only( options_->rescore_only() );
		stepwise_pose_minimizer->set_cartesian( options_->cart_min() );
		stepwise_pose_minimizer->set_min_type( options_->min_type() );
		stepwise_pose_minimizer->set_min_tolerance( options_->min_tolerance() );
		stepwise_pose_minimizer->set_use_coordinate_constraints( !options_->skip_coord_constraints() );
		stepwise_pose_minimizer->set_num_pose_minimize( options_->num_pose_minimize() );

		if ( !options_->skip_minimize() ){
			stepwise_pose_minimizer->apply( pose );
			StepWiseClustererOP stepwise_clusterer = new StepWiseClusterer( stepwise_pose_minimizer->pose_list(),
																																			moving_res_list_, options_, full_optimize_ /*force_align*/ );
			stepwise_clusterer->cluster();
			pose_list = stepwise_clusterer->get_pose_list();
			pose = *pose_list[ 1 ];
		}

		if ( options_->output_minimized_pose_list() ) {
			core::io::silent::SilentFileDataOP sfd = new core::io::silent::SilentFileData; // silly
			for ( Size n = 1; n <= pose_list.size(); n++ ){
				Pose & pose = *pose_list[ n ];
				std::string const tag = "S_"+ ObjexxFCL::string_of( n-1 /* start with zero */);
				protocols::stepwise::sampling::protein::output_silent_struct( pose, get_native_pose(), options_->silent_file(),
																																			tag, sfd,
																																			job_parameters_->working_calc_rms_res() );
			}
		}
	}

	////////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinModeler::set_moving_res_and_reset( Size const moving_res ){
		moving_res_ = moving_res;
		moving_res_list_.clear(); // signal to recompute when pose is supplied.
	}

	/////////////////////////////////////////////////////
	StepWiseProteinModelerOptionsCOP
	StepWiseProteinModeler::options(){
		return options_;
	}

	//////////////////////////////////////////////////////////////////
	void
	StepWiseProteinModeler::set_options( StepWiseProteinModelerOptionsCOP options ){
		options_ = options;
	}

////////////////////////////////////////////////////////////////////
void
StepWiseProteinModeler::figure_out_moving_res_list( pose::Pose const & pose ){
	if ( moving_res_list_.size() == 0 ){ // in principle, SWM can explicitly give moving_res_list
		// otherwise, figure it out from moving_res_ (single residue!)
		figure_out_protein_moving_res_list_from_most_distal_res( pose, moving_res_ );
	}
}

////////////////////////////////////////////////////////////////////
void
StepWiseProteinModeler::figure_out_protein_moving_res_list_from_most_distal_res( pose::Pose const & pose, Size const moving_res ) {

	moving_res_list_.clear();
	if ( moving_res == 0 ) return;

	moving_res_list_.push_back( moving_res );

	// go back another residue -- this was the default choice in protein SWA.
	utility::vector1< Size > const & fixed_domain_map = get_fixed_domain_from_full_model_info_const( pose );
	Size const upstream_res = pose.fold_tree().get_parent_residue( moving_res );
	if ( fixed_domain_map[ upstream_res ] == 0 ) moving_res_list_.push_back( upstream_res );

	// check if CCD-closure  is necessary.
	utility::vector1< Size > const cutpoints_closed = figure_out_moving_cutpoints_closed_from_moving_res( pose, moving_res );
	do_ccd_ = ( cutpoints_closed.size() > 0 );
	if ( do_ccd_ ){
		runtime_assert( cutpoints_closed.size() == 1 ); // for now. can generalize soon.
		int cutpoint_closed = static_cast<int>( cutpoints_closed[ 1 ] );

		// as in protein SWA, choose two bridge residues for CCD closure on the 'other side' of sampled residue.
		utility::vector1< int > offsets;
		if ( moving_res_list_.has_value( cutpoint_closed ) ){
			offsets = make_vector1( +1, +2 );
		} else if ( moving_res_list_.has_value( cutpoint_closed+1 ) ){
			offsets = make_vector1( -1, 0 );
		} else {
			offsets = make_vector1( 0, +1 ); // bracket CCD closure point, since moving_res does not.
		}
		utility::vector1< Size > working_bridge_res;
		for ( Size n = 1; n <= offsets.size(); n++ ) {
			int const bridge_res = cutpoint_closed + offsets[n];
			if( bridge_res < 1 ) continue;
			runtime_assert( !moving_res_list_.has_value( bridge_res ) );
			if ( fixed_domain_map[ bridge_res ] == 0 ) working_bridge_res.push_back( bridge_res );
		}
		bridge_res_ = const_full_model_info( pose ).sub_to_full( working_bridge_res );
	} else {
		bridge_res_.clear();
	}
}


//////////////////////////////////////////////////////////////////
void
StepWiseProteinModeler::align_pose_and_add_rmsd_constraints( pose::Pose & pose ) const {
	if ( !get_native_pose() ) return;

	utility::vector1< Size > root_partition_res = figure_out_root_partition_res( pose, moving_res_list_ );
	if ( root_partition_res.size() == 0 ) root_partition_res.push_back( pose.fold_tree().root() );

	// can later generalize to use 'reference_pose', not necessarily native_pose.
	general::StepWisePoseAligner pose_aligner( *get_native_pose() );
	pose_aligner.set_root_partition_res( root_partition_res );
	pose_aligner.apply( pose );
	if ( options_->rmsd_screen() > 0 ) pose_aligner.create_coordinate_constraints( pose, options_->rmsd_screen() );
}

//////////////////////////////////////////////////////////////////
void
StepWiseProteinModeler::set_skip_sampling( bool const & setting ){
	if( options_->disallow_backbone_sampling() != setting ){
		// needed to get aroud constant status of options -- which is desirable in most contexts.
		sampling::protein::StepWiseProteinModelerOptionsOP new_options = options_->clone();
		new_options->set_disallow_backbone_sampling( setting );
		set_options( new_options );
	}
}

//////////////////////////////////////////////////////////////////
void
StepWiseProteinModeler::set_frag_files( utility::vector1< std::string > const & frag_files ){
	frag_files_ = frag_files;
}

//////////////////////////////////////////////////////////////////
void
StepWiseProteinModeler::set_input_streams( utility::vector1< InputStreamWithResidueInfoOP > const & input_streams ){
	input_streams_ = input_streams;
}

//////////////////////////////////////////////////////////////////
void
StepWiseProteinModeler::set_job_parameters( StepWiseProteinJobParametersCOP job_parameters ){
	job_parameters_ = job_parameters;
}


//////////////////////////////////////////////////////////////////////////////
bool
StepWiseProteinModeler::sampling_successful( utility::vector1< PoseOP > & pose_list ){
	Size const num_sampled = pose_list.size();
	if ( num_sampled == 0 ){
		TR << "WARNING! WARNING! WARNING! pose_list_.size() == 0! " << std::endl;
		if ( !options_->output_minimized_pose_list() ) return false; // don't do a minimize...
	}
	return true;
}


} //protein
} //sampling
} //stepwise
} //protocols
