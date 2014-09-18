// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/modeler/protein/StepWiseProteinModeler.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/modeler/protein/StepWiseProteinModeler.hh>
#include <protocols/stepwise/modeler/options/StepWiseModelerOptions.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.hh>
#include <protocols/stepwise/modeler/working_parameters/util.hh>
#include <protocols/stepwise/legacy/modeler/protein/StepWiseProteinMinimizer.hh>
#include <protocols/stepwise/legacy/modeler/protein/StepWiseProteinConnectionSampler.hh>
#include <protocols/stepwise/modeler/protein/util.hh>
#include <protocols/stepwise/modeler/protein/InputStreamWithResidueInfo.hh>
#include <protocols/stepwise/modeler/align/StepWiseLegacyClusterer.hh>
#include <protocols/stepwise/modeler/align/StepWisePoseAligner.hh>
#include <protocols/stepwise/monte_carlo/util.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/full_model_info/util.hh>
#include <utility/stream_util.hh>
#include <utility/tools/make_vector1.hh>

#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>

#include <basic/Tracer.hh>

static thread_local basic::Tracer TR( "protocols.stepwise.legacy.modeler.protein.StepWiseProteinModeler" );

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
namespace legacy {
namespace modeler {
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
		working_parameters_ = src.working_parameters_;
		input_streams_ = src.input_streams_;
		bridge_res_ = src.bridge_res_;
		working_minimize_res_ = src.working_minimize_res_;

		return *this;
	}

	//////////////////////////////////////////////////////////////////////////////
	// make this closer to StepWiseRNA_Modeler in anticipation of unification.
	void
	StepWiseProteinModeler::apply( pose::Pose & pose ){
		initialize_working_parameters_and_root( pose );

		utility::vector1< PoseOP > pose_list;
		do_residue_modeler( pose, pose_list );
		if ( modeler_successful( pose_list ) ) do_minimizing( pose, pose_list );

		pose.remove_constraints();
		reinitialize();
	}

//////////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinModeler::reinitialize(){
		// Important: make sure that the next time this is used, job parameters is set explicitly -- or it will be reset.
		working_parameters_ = 0;
		moving_res_list_.clear();
		bridge_res_.clear();
		working_minimize_res_.clear();
		do_ccd_ = false;
	}

	//////////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinModeler::initialize_working_parameters_and_root( pose::Pose & pose ){
		pose::full_model_info::make_sure_full_model_info_is_setup( pose );
		if ( working_parameters_ != 0 ) return;
		figure_out_moving_res_list( pose );
		revise_root_and_moving_res_list( pose, moving_res_list_ ); // specify reference_res_? [i.e. anchor_res?]
		align_pose_and_add_rmsd_constraints( pose, get_native_pose(), moving_res_list_, options_->rmsd_screen() );
		working_parameters_ = setup_working_parameters_for_stepwise_with_full_model_info( pose );
	}

	////////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinModeler::do_residue_modeler( core::pose::Pose & pose,
																							 utility::vector1< PoseOP > & pose_list ){

		StepWiseProteinConnectionSampler stepwise_protein_residue_sampler( working_parameters_ );
		stepwise_protein_residue_sampler.set_options( options_ );
		stepwise_protein_residue_sampler.set_do_ccd( do_ccd_ );
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
	stepwise::modeler::working_parameters::StepWiseWorkingParametersOP
	StepWiseProteinModeler::setup_working_parameters_for_stepwise_with_full_model_info( core::pose::Pose & pose )
	{
		pose::full_model_info::make_sure_full_model_info_is_setup( pose );
		return setup_protein_working_parameters_for_swa( moving_res_list_, pose,
																										 get_native_pose(),
																										 bridge_res_,
																										 working_minimize_res_ /* specify in only minimizing subset of minimizable residues */);
	}

	////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinModeler::do_minimizing( pose::Pose & pose,
																				 utility::vector1< PoseOP > & pose_list ) {

		using namespace core::scoring;
		using namespace protocols::stepwise::modeler::align;

		StepWiseProteinMinimizerOP stepwise_pose_minimizer;
		if ( pose_list.size() == 0 ) pose_list = make_vector1( pose.clone() );
		stepwise_pose_minimizer = new StepWiseProteinMinimizer( pose_list, moving_res_list_ );

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
		stepwise_pose_minimizer->set_native_pose( working_parameters_->working_native_pose() );
		stepwise_pose_minimizer->set_calc_rms_res( working_parameters_->working_calc_rms_res() ); // used for calculating rmsds to native.
		stepwise_pose_minimizer->set_fixed_res( working_parameters_->working_fixed_res() );
		stepwise_pose_minimizer->set_move_takeoff_torsions( !options_->disable_sampling_of_loop_takeoff() );
		stepwise_pose_minimizer->set_rescore_only( options_->rescore_only() );
		stepwise_pose_minimizer->set_cartesian( options_->cart_min() );
		stepwise_pose_minimizer->set_min_type( options_->min_type() );
		stepwise_pose_minimizer->set_min_tolerance( options_->min_tolerance() );
		stepwise_pose_minimizer->set_use_coordinate_constraints( !options_->skip_coord_constraints() );
		stepwise_pose_minimizer->set_num_pose_minimize( options_->num_pose_minimize() );

		if ( !options_->skip_minimize() ){
			stepwise_pose_minimizer->apply( pose );
			StepWiseLegacyClustererOP stepwise_clusterer = new StepWiseLegacyClusterer( stepwise_pose_minimizer->pose_list(),
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
				protocols::stepwise::modeler::protein::output_silent_struct( pose, get_native_pose(), options_->silent_file(),
																																			tag, sfd,
																																			working_parameters_->working_calc_rms_res() );
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
	options::StepWiseModelerOptionsCOP
	StepWiseProteinModeler::options(){
		return options_;
	}

	//////////////////////////////////////////////////////////////////
	void
	StepWiseProteinModeler::set_options( options::StepWiseModelerOptionsCOP options ){
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

	figure_out_protein_modeling_info( pose, moving_res, moving_res_list_, do_ccd_, bridge_res_ );
}


//////////////////////////////////////////////////////////////////
void
StepWiseProteinModeler::set_skip_modeler( bool const & setting ){
	if( options_->disallow_backbone_modeler() != setting ){
		// needed to get aroud constant status of options -- which is desirable in most contexts.
		modeler::options::StepWiseModelerOptionsOP new_options = options_->clone();
		new_options->set_disallow_backbone_modeler( setting );
		set_options( new_options );
	}
}

//////////////////////////////////////////////////////////////////
void
StepWiseProteinModeler::set_input_streams( utility::vector1< InputStreamWithResidueInfoOP > const & input_streams ){
	input_streams_ = input_streams;
}

//////////////////////////////////////////////////////////////////
void
StepWiseProteinModeler::set_working_parameters( stepwise::modeler::working_parameters::StepWiseWorkingParametersCOP working_parameters ){
	working_parameters_ = working_parameters;
}


//////////////////////////////////////////////////////////////////////////////
bool
StepWiseProteinModeler::modeler_successful( utility::vector1< PoseOP > & pose_list ){
	Size const num_sampled = pose_list.size();
	if ( num_sampled == 0 ){
		TR << "WARNING! WARNING! WARNING! pose_list_.size() == 0! " << std::endl;
		if ( !options_->output_minimized_pose_list() ) return false; // don't do a minimize...
	}
	return true;
}

StepWiseProteinWorkingParametersOP
setup_protein_working_parameters_for_swa( utility::vector1< Size > const & moving_res_list,
																					pose::Pose const & pose,
																					pose::PoseCOP native_pose,
																					utility::vector1< Size > const & bridge_res,
																					utility::vector1< Size > const & working_minimize_res ){

		using namespace core::pose::full_model_info;

		FullModelInfo const & full_model_info = const_full_model_info( pose );
		std::string const & full_sequence = full_model_info.full_sequence();
		utility::vector1< Size > const & res_list = full_model_info.res_list();
		utility::vector1< Size > const & fixed_domain_map = full_model_info.fixed_domain_map();

		utility::vector1< Size > const & working_moving_res_list = moving_res_list;
		utility::vector1< Size > working_moving_suite_list;
		for ( Size n = 1; n <= moving_res_list.size(); n++ ){
			Size const & moving_res = moving_res_list[n];
			Size const parent_res = pose.fold_tree().get_parent_residue( moving_res );
			if ( pose.fold_tree().jump_nr( moving_res, parent_res ) == 0 ){
				working_moving_suite_list.push_back( std::min( moving_res, parent_res ) );
			}
		}

		utility::vector1< Size > is_working_res( full_sequence.size(), Size(0) );
		std::string working_sequence;
		for ( Size n = 1; n <= res_list.size(); n++ ) {
			is_working_res[ res_list[ n ] ] = 1;
			// following magic numbers match what is in StepWiseProteinPoseSetup:
			if ( moving_res_list.has_value( n ) ) is_working_res[ res_list[ n ] ] = 999 /*horrible HACK!*/;
			if ( bridge_res.has_value( res_list[ n ] ) ) is_working_res[ res_list[ n ] ] = 123 /*horrible HACK!*/;
			working_sequence += full_sequence[ n-1 ];
		}

		utility::vector1< Size > calc_rms_res_ = full_model_info.full_model_parameters()->get_res_list( CALC_RMS );
		if ( calc_rms_res_.size() == 0 ) calc_rms_res_ = full_model_info.sub_to_full( moving_res_list );

		utility::vector1< Size > fixed_res;
		utility::vector1< Size > const & extra_minimize_res = const_full_model_info( pose ).extra_minimize_res();
		for ( Size n = 1; n <= full_sequence.size(); n++ )	{
			if ( extra_minimize_res.has_value( n ) ) continue;
			if ( fixed_domain_map[ n ] > 0 ) {
				fixed_res.push_back( n );
				continue;
			}
			if ( working_minimize_res.size() > 0 ){ // finer specification  of what to minimize
				if( res_list.has_value( n ) && !working_minimize_res.has_value( res_list.index( n ) ) ) fixed_res.push_back( n );
			}
		}

		pose::PoseOP working_native_pose;
		if ( native_pose ) {
			working_native_pose = new pose::Pose;
			protocols::stepwise::modeler::pdbslice( *working_native_pose, *native_pose, res_list );
		}
		// ALSO WILL NEED A DEFAULT BRIDGE RES SETTER!
		StepWiseProteinWorkingParametersOP working_parameters = new StepWiseProteinWorkingParameters;
		working_parameters->set_full_to_sub( full_model_info.full_to_sub() ); // will update sub_to_full
		working_parameters->set_is_working_res( is_working_res );
		working_parameters->set_sequence( full_model_info.full_sequence() );
		//		working_parameters->set_working_sequence( working_sequence ); // auto-updated with set_sequence
		//		working_parameters->set_working_res_list( res_list ); // figured out from full_to_sub
		working_parameters->set_working_moving_res_list( working_moving_res_list );
		//working_parameters->set_working_moving_suite_list( working_moving_suite_list );  // set at the same time as working_moving_res_list
		//working_parameters->set_chain_boundaries(); // Not in use later. legacy of PoseSetup
		//working_parameters->set_which_chain_has_moving_res(); // Not in use later. legacy of PoseSetup
		//working_parameters->gap_size(); // Not in use later. legacy of PoseSetup
		//working_parameters->first_chain_break_res(); // Not in use later. legacy of PoseSetup
		//working_parameters->is_prepend(); // Not in use later. legacy of PoseSetup
		//working_parameters->is_internal(); // Not in use later. legacy of PoseSetup
		//working_parameters->partition_definition(); // Not in use later. legacy of PoseSetup
		working_parameters->set_working_native_pose( working_native_pose );

		working_parameters->set_working_fixed_res( full_model_info.full_to_sub( fixed_res ) );
		//working_parameters->set_working_terminal_res( full_model_info.full_to_sub( terminal_res ) ); // not in use for proteins.
		//		working_parameters->set_working_superimpose_res( full_model_info.full_to_sub( superimpose_res ) ); // not in use! assumes pose is already superimposed.
		working_parameters->set_working_calc_rms_res( full_model_info.full_to_sub( calc_rms_res_ ) );
		working_parameters->set_working_bridge_res( full_model_info.full_to_sub( bridge_res ) );
		//		working_parameters->set_working_moving_pos( full_model_info.full_to_sub( moving_pos ) ); // not in use! legacy of PoseSetup

		return working_parameters;

	}




} //protein
} //modeler
} //legacy
} //stepwise
} //protocols
