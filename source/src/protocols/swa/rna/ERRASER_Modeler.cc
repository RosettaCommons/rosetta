// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/swa/rna/ERRASER_Modeler.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu

#include <protocols/swa/rna/ERRASER_Modeler.hh>
#include <protocols/swa/rna/StepWiseRNA_AnalyticalLoopCloseSampler.hh>
#include <protocols/swa/rna/StepWiseRNA_BaseCentroidScreener.hh>
#include <protocols/swa/rna/StepWiseRNA_JobParameters.hh>
#include <protocols/swa/rna/StepWiseRNA_JobParametersSetup.hh>
#include <protocols/swa/rna/StepWiseRNA_Minimizer.hh>
#include <protocols/swa/rna/StepWiseRNA_BaseCentroidScreener.hh>

#include <core/types.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/Pose.hh>

#include <utility/tools/make_vector1.hh>
#include <basic/Tracer.hh>

#include <ObjexxFCL/string.functions.hh>

static basic::Tracer TR( "protocols.swa.rna.ERRASER_Modeler" );

////////////////////////////////////////////////////////////////////////////////////////
// This needs to be unified with StepWiseRNA_Modeler!
// Indeed ERRASER_Modeler should probably be removed entirely, and can be called
// with appropriate flags (e.g., -analytical_loop_close or something) form swa_rna_main.
////////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace swa {
namespace rna {

	//Constructor
	ERRASER_Modeler::ERRASER_Modeler( core::Size const sample_res,
																		core::scoring::ScoreFunctionOP scorefxn ) :
		moving_res_( utility::tools::make_vector1( sample_res ) ),
		scorefxn_( scorefxn ),
		silent_file_( "" ),
		sampler_num_pose_kept_( 108 ),
		num_pose_minimize_( 99999 ),
		nstruct_( 1 ),
		num_sampled_( 0 ),
		sampler_native_screen_rmsd_cutoff_( 2.0 ),
		cluster_rmsd_( 0.5 ),
		native_edensity_score_cutoff_( -1 ),
		fast_( false ),
		medium_fast_( false ),
		sampler_native_rmsd_screen_( false ),
		o2star_screen_( true ),
		verbose_( false ),
		distinguish_pucker_( true ),
		finer_sampling_at_chain_closure_( false ),
		PBP_clustering_at_chain_closure_( false ),
		allow_syn_pyrimidine_( false ),
		extra_syn_chi_rotamer_( false ),
		extra_anti_chi_rotamer_( false ),
		use_phenix_geo_( false ),
		centroid_screen_( true ),
		VDW_atr_rep_screen_( true ),
		force_centroid_interaction_( false ),
		choose_random_( false ),
		skip_sampling_( false ),
		perform_minimize_( true ),
		minimize_and_score_sugar_( true ),
		minimize_and_score_native_pose_( false ),
		rm_virt_phosphate_( false ),
		VDW_rep_alignment_RMSD_CUTOFF_( 0.001 ),
		output_pdb_( false ),
		output_minimized_pose_data_list_( false )
	{}

	//Destructor
	ERRASER_Modeler::~ERRASER_Modeler()
	{}

	/////////////////////
	std::string
	ERRASER_Modeler::get_name() const {
		return "ERRASER_Modeler";
	}

	//////////////////////////////////////////////////////////////////////////////
	void
	ERRASER_Modeler::apply( core::pose::Pose & pose ){

		using namespace core::pose;
		using namespace core::chemical;
		using namespace core::kinematics;
		using namespace core::scoring;
		using namespace protocols::swa::rna;

		if ( !job_parameters_ ) job_parameters_ = setup_job_parameters_for_erraser( moving_res_, pose );

		StepWiseRNA_AnalyticalLoopCloseSampler stepwise_rna_residue_sampler ( job_parameters_ );
		stepwise_rna_residue_sampler.set_silent_file ( silent_file_ + "_sampling" );
		stepwise_rna_residue_sampler.set_scorefxn ( scorefxn_ );
		stepwise_rna_residue_sampler.set_num_pose_kept ( sampler_num_pose_kept_ );
		stepwise_rna_residue_sampler.set_fast ( fast_ );
		stepwise_rna_residue_sampler.set_medium_fast ( medium_fast_ );
		stepwise_rna_residue_sampler.set_native_rmsd_screen ( sampler_native_rmsd_screen_ );
		stepwise_rna_residue_sampler.set_native_screen_rmsd_cutoff ( sampler_native_screen_rmsd_cutoff_ );
		stepwise_rna_residue_sampler.set_o2star_screen ( o2star_screen_ );
		stepwise_rna_residue_sampler.set_verbose ( verbose_ );
		stepwise_rna_residue_sampler.set_cluster_rmsd (	cluster_rmsd_	);
		stepwise_rna_residue_sampler.set_distinguish_pucker ( distinguish_pucker_ );
		stepwise_rna_residue_sampler.set_finer_sampling_at_chain_closure ( finer_sampling_at_chain_closure_ );
		stepwise_rna_residue_sampler.set_PBP_clustering_at_chain_closure ( PBP_clustering_at_chain_closure_ );
		stepwise_rna_residue_sampler.set_allow_syn_pyrimidine( allow_syn_pyrimidine_ );
		stepwise_rna_residue_sampler.set_extra_syn_chi_rotamer ( extra_syn_chi_rotamer_ );
		stepwise_rna_residue_sampler.set_extra_anti_chi_rotamer ( extra_anti_chi_rotamer_ );
		stepwise_rna_residue_sampler.set_use_phenix_geo ( use_phenix_geo_  );
		stepwise_rna_residue_sampler.set_centroid_screen ( centroid_screen_ );
		stepwise_rna_residue_sampler.set_VDW_atr_rep_screen ( VDW_atr_rep_screen_ );

		StepWiseRNA_BaseCentroidScreenerOP base_centroid_screener = new StepWiseRNA_BaseCentroidScreener ( pose, job_parameters_ );
		stepwise_rna_residue_sampler.set_base_centroid_screener ( base_centroid_screener );

		StepWiseRNA_VDW_BinScreenerOP user_input_VDW_bin_screener = new StepWiseRNA_VDW_BinScreener();
		if ( VDW_rep_screen_info_.size() > 0 ) {
			user_input_VDW_bin_screener->set_VDW_rep_alignment_RMSD_CUTOFF ( VDW_rep_alignment_RMSD_CUTOFF_ );
			user_input_VDW_bin_screener->setup_using_user_input_VDW_pose( VDW_rep_screen_info_, pose, job_parameters_ );
			user_input_VDW_bin_screener->set_output_pdb( output_pdb_ );
		}
		stepwise_rna_residue_sampler.set_user_input_VDW_bin_screener ( user_input_VDW_bin_screener );
		stepwise_rna_residue_sampler.set_force_centroid_interaction ( force_centroid_interaction_ );

		if ( choose_random_ ){
			stepwise_rna_residue_sampler.set_choose_random( true );
			stepwise_rna_residue_sampler.set_cluster_rmsd( 0.0 ); // don't cluster.
		}

		// nstruct is usually 1, unless choose_random is on -- then we need to call several times to
		// get a bunch of poses back.
		if ( ! skip_sampling_ ) {
			for ( Size n = 1; n <= nstruct_; n++ ){
				stepwise_rna_residue_sampler.apply ( pose );
			}
		}

		utility::vector1< pose_data_struct2 > & pose_data_list = stepwise_rna_residue_sampler.get_pose_data_list();
		num_sampled_ = pose_data_list.size();
		if ( num_sampled_ == 0 ){

			TR << "WARNING! WARNING! WARNING! pose_data_list.size() == 0! " << std::endl;
			//			if ( ! skip_sampling_ ) utility_exit_with_message( "No op in ERRASER_modeler!" );
			pose_data_struct2 data_struct;
			data_struct.pose_OP = new Pose( pose );
			data_struct.score = 0.0;
			data_struct.tag = "input_pose";
			pose_data_list.push_back( data_struct );
		}

		if ( minimize_and_score_native_pose_ ) {
			if ( !get_native_pose() ) utility_exit_with_message ( "minimize_and_score_native_pose == True but user did not pass in native pose" );
			pose_data_struct2 native_data_struct;
			native_data_struct.pose_OP = new Pose( *get_native_pose() ); //hopefully this clones...
			native_data_struct.score = 0.0;
			native_data_struct.tag = "working_native_pose";
			pose_data_list.push_back ( native_data_struct );
		}

		// let's output the final pose_data_list, just to have a look
		if ( verbose_ ) stepwise_rna_residue_sampler.output_pose_data_list ( silent_file_ + "_final_sample" );

		////////////////////////////////////////////////////////////////
		StepWiseRNA_Minimizer stepwise_rna_minimizer ( stepwise_rna_residue_sampler.get_pose_data_list(), job_parameters_ );
		stepwise_rna_minimizer.set_silent_file ( silent_file_ );
		stepwise_rna_minimizer.set_verbose (  verbose_ );
		stepwise_rna_minimizer.set_scorefxn ( scorefxn_ );
		stepwise_rna_minimizer.set_centroid_screen (  centroid_screen_ );
		stepwise_rna_minimizer.set_base_centroid_screener ( base_centroid_screener );
		stepwise_rna_minimizer.set_perform_minimize(  perform_minimize_ );
		stepwise_rna_minimizer.set_native_rmsd_screen (  sampler_native_rmsd_screen_ );
		stepwise_rna_minimizer.set_native_edensity_score_cutoff ( native_edensity_score_cutoff_ );
		stepwise_rna_minimizer.set_rm_virt_phosphate (  rm_virt_phosphate_ );
		stepwise_rna_minimizer.set_native_screen_rmsd_cutoff (  sampler_native_screen_rmsd_cutoff_ + 1 ); //+1 for leniency Sept 20, 2010

		if ( num_pose_minimize_ > 0 ) stepwise_rna_minimizer.set_num_pose_minimize ( num_pose_minimize_ );
		stepwise_rna_minimizer.set_minimize_and_score_sugar ( minimize_and_score_sugar_ );
		stepwise_rna_minimizer.set_user_input_VDW_bin_screener ( user_input_VDW_bin_screener );
		stepwise_rna_minimizer.set_output_minimized_pose_data_list( output_minimized_pose_data_list_ );
		stepwise_rna_minimizer.apply ( pose );

		// Need to make sure that final pose output is the lowest scoring one. Is it?


	}


	///////////////////////////////////////////////////////////////////////////////////////////////////
	// This could go into a Util.hh if it ends up being more useful.
	// Briefly, we need to make a StepWiseRNA_JobParameters object that will
	// be fed into various StepWiseRNA movers. Setting this up can be quite complicated...
	// its become a grab bag of residue lists referring to the global pose, the working pose,
	// sequence mappings, "Is_Prepend_map", etc.
	//
	StepWiseRNA_JobParametersOP
	ERRASER_Modeler::setup_job_parameters_for_erraser( utility::vector1< Size > moving_res, core::pose::Pose const & pose ){

		using namespace core::pose;

		if ( moving_res.size() != 1 ) utility_exit_with_message( "ERRASER requires exactly 1 number in -moving_res" );
		Size const rebuild_res =  moving_res[1];

		std::string full_sequence = pose.sequence();
		Size nres = pose.total_residue();

		// what if there is a virtual residue? need to remove it, actually, before running stepwise_rna_job_parameters_setup.
		if ( full_sequence[nres - 1] == 'X' ){
			full_sequence = full_sequence.substr( 0, nres - 1 );
			nres -= 1;
		}
		utility::vector1< Size > not_rebuild_res;
		for ( Size n = 1; n <= nres; n++ ) if ( n != rebuild_res ) not_rebuild_res.push_back( n );

		utility::vector1< Size > input_res1, input_res2 /*blank*/, cutpoint_open /*blank*/;
		input_res1 = not_rebuild_res;

		Size cutpoint_closed = rebuild_res;
		if ( !pose.fold_tree().is_cutpoint( cutpoint_closed ) ) utility_exit_with_message( "ERRASER requires a chainbreak right at sampled residue" );

		StepWiseRNA_JobParametersSetup stepwise_rna_job_parameters_setup( moving_res,
																																			 full_sequence,
																																			 input_res1,
																																			 input_res2,
																																			 cutpoint_open,
																																			 cutpoint_closed );

		utility::vector1< Size > fixed_res = not_rebuild_res;
		stepwise_rna_job_parameters_setup.set_fixed_res( fixed_res );

		utility::vector1< Size > rmsd_res_list;
		rmsd_res_list.push_back( rebuild_res );
		stepwise_rna_job_parameters_setup.set_rmsd_res_list( rmsd_res_list );

		if ( rebuild_res == 1 || rebuild_res == pose.total_residue() ) utility_exit_with_message( "ERRASER requires that residue is not at terminus!" );

		// not sure about the following -- instead, how about reading jump residues from within pose itself?
		utility::vector1< std::string > jump_point_pair_list;
		jump_point_pair_list.push_back( ObjexxFCL::string_of( rebuild_res - 1 ) + "-" + ObjexxFCL::string_of( rebuild_res + 1 ) );
		stepwise_rna_job_parameters_setup.set_jump_point_pair_list( jump_point_pair_list ); //Important!: Needs to be called after set_fixed_res

		utility::vector1< std::string > alignment_res; //why is this a string vector?????
		for ( Size n = 1; n <= fixed_res.size(); n++ ) alignment_res.push_back( ObjexxFCL::string_of( fixed_res[n] ) );
		stepwise_rna_job_parameters_setup.set_alignment_res( alignment_res );
		stepwise_rna_job_parameters_setup.set_native_alignment_res( fixed_res );

		// could use this later to minimize more residues...
		//stepwise_rna_job_parameters_setup.set_global_sample_res_list( option[ global_sample_res_list ]() ); //March 20, 2011

		// NOT SURE ABOUT THIS. false by deafult, but shows up as true in 'normal' erraser runs.
		stepwise_rna_job_parameters_setup.set_allow_chain_boundary_jump_partner_right_at_fixed_BP ( true );

		stepwise_rna_job_parameters_setup.set_add_virt_res_as_root( true );

		stepwise_rna_job_parameters_setup.apply();

		StepWiseRNA_JobParametersOP job_parameters = stepwise_rna_job_parameters_setup.job_parameters();

		job_parameters->set_working_native_pose( get_native_pose() );

		// should we also set the fold_tree here -- just take the pose's actual fold tree?
		// that fold_tree is only used in PoseSetup, and in JobParametersSetup, and not downstream
		// in any modelers...  -- rhiju

		// user input fixed_res...
		if ( fixed_res_.size() > 0 ) 	job_parameters->set_working_fixed_res( fixed_res_ );

		return job_parameters;
	}


	///////////////////////////////////////////////////////////////////////////////
	void
	ERRASER_Modeler::set_job_parameters( StepWiseRNA_JobParametersCOP job_parameters ){ job_parameters_ = job_parameters;	}

	///////////////////////////////////////////////////////////////////////////////
	void
	ERRASER_Modeler::set_native_pose( core::pose::PoseCOP native_pose ){ native_pose_ = native_pose; }

	///////////////////////////////////////////////////////////////////////////////
	core::pose::PoseCOP
	ERRASER_Modeler::get_native_pose(){ return native_pose_; }


} //rna
} //swa
} //protocols
