// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @file swa_monte_Carlo.cc
/// @author Rhiju Das (rhiju@stanford.edu)

// libRosetta headers
#include <core/types.hh>

// do we need all these?
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueSelector.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Conformation.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rna/RNA_FittedTorsionInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/rms_util.hh>
#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/id/types.hh>
#include <core/id/TorsionID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.hh>
#include <core/init/init.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>

//////////////////////////////////////////////////
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh> // for option[ out::file::silent  ] and etc.
#include <basic/options/keys/in.OptionKeys.gen.hh> // for option[ in::file::tags ] and etc.
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/option_macros.hh>
///////////////////////////////////////////////////
#include <protocols/idealize/idealize.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

#include <utility/vector1.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <numeric/xyzVector.hh>
#include <numeric/conversions.hh>

#include <core/pose/MiniPose.hh>
#include <core/pose/util.hh>

#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/BinaryRNASilentStruct.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>
#include <core/import_pose/import_pose.hh>

//////////////////////////////////////////////////////////
#include <protocols/viewer/viewers.hh>
#include <protocols/rna/RNA_ProtocolUtil.hh>
#include <protocols/rna/RNA_LoopCloser.hh>
#include <protocols/swa/rna/StepWiseRNA_OutputData.hh>
#include <protocols/swa/rna/StepWiseRNA_ResidueSampler.hh>
#include <protocols/swa/rna/StepWiseRNA_PoseSetup.fwd.hh>
#include <protocols/swa/rna/StepWiseRNA_PoseSetup.hh>
#include <protocols/swa/rna/StepWiseRNA_JobParametersSetup.hh>
#include <protocols/swa/rna/StepWiseRNA_JobParameters.hh>

#include <core/pose/full_model_info/FullModelInfo.hh>
#include <protocols/swa/monte_carlo/RNA_AddMover.hh>
#include <protocols/swa/monte_carlo/RNA_DeleteMover.hh>
#include <protocols/swa/monte_carlo/RNA_AddOrDeleteMover.hh>
#include <protocols/swa/monte_carlo/RNA_O2StarMover.hh>
#include <protocols/swa/monte_carlo/RNA_TorsionMover.hh>
#include <protocols/swa/monte_carlo/RNA_SWA_MonteCarloMover.hh>
#include <protocols/swa/monte_carlo/RNA_SWA_MonteCarloUtil.hh>
#include <protocols/swa/monte_carlo/types.hh>

#include <numeric/random/random.hh>
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>


// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>  //Test
#include <cctype>
#include <iomanip>
#include <map>
#include <cstdlib>
#include <ctime>
#include <unistd.h>


#include <list>
#include <stdio.h>
#include <math.h>

using namespace core;
using namespace protocols;
using namespace ObjexxFCL;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using utility::vector1;
using io::pdb::dump_pdb;

typedef  numeric::xyzMatrix< Real > Matrix;

// SWA Monte Carlo -- July 12, 2012 -- Rhiju Das
//
// TO DO
//
//  Clean up pose setup, inherited from SWA code -- all these options and setup functions
//    should go into their own .cc file.
//
//  Add screen for building new base [base atr/rep] -- like SWA
//
//  Set up GAAA tetraloop run [should be faster test case], with buildup
//   from either end.
//
//  Set up chainbreak variants when all residues are formed.
//
//  Set up loop closer move (perhaps just analytical loop close move?)
//
//  Set up constraints when gap is 1,2, etc.
//
//  encapsulate -- move into a namespace? lay out plans for others?
//

static numeric::random::RandomGenerator RG(2391121);  // <- Magic number, do not change it!


// A lot of these options should be placed into an 'official' namespace
// and the SWA RNA pose setup should go to its own function
OPT_KEY( Boolean, do_not_sample_multiple_virtual_sugar)
OPT_KEY( Boolean, sample_ONLY_multiple_virtual_sugar)
OPT_KEY( Boolean, skip_sampling)
OPT_KEY( Boolean, skip_clustering)
OPT_KEY( Boolean, minimize_and_score_native_pose)
OPT_KEY( Integer, num_pose_minimize)
OPT_KEY( Boolean, combine_long_loop_mode)
OPT_KEY( Integer, job_queue_ID)
OPT_KEY( String, filter_output_filename)
OPT_KEY( Boolean, filter_for_previous_contact)
OPT_KEY( Boolean, filter_for_previous_clash)
OPT_KEY( Boolean, filterer_undercount_ribose_rotamers)
OPT_KEY( Boolean, combine_helical_silent_file)
OPT_KEY( Boolean, exclude_alpha_beta_gamma_sampling)
OPT_KEY( Boolean, debug_eplison_south_sugar_mode)
OPT_KEY( Boolean, rebuild_bulge_mode)
OPT_KEY( Boolean, sampler_include_torsion_value_in_tag)
OPT_KEY( Boolean, sampler_extra_anti_chi_rotamer)
OPT_KEY( Boolean, sampler_extra_syn_chi_rotamer)
OPT_KEY( Boolean, sampler_extra_beta_rotamer)
OPT_KEY( Boolean, sampler_extra_epsilon_rotamer)
OPT_KEY( Boolean, sample_both_sugar_base_rotamer)
OPT_KEY( Boolean, reinitialize_CCD_torsions)
OPT_KEY( Boolean, PBP_clustering_at_chain_closure)
OPT_KEY( Boolean, finer_sampling_at_chain_closure)
OPT_KEY( StringVector, 	VDW_rep_screen_info)
OPT_KEY( Real, 	VDW_rep_alignment_RMSD_CUTOFF)
OPT_KEY( Boolean, graphic )
OPT_KEY( Real, Real_parameter_one )
OPT_KEY( Boolean, add_lead_zero_to_tag )
OPT_KEY( Boolean, distinguish_pucker )
OPT_KEY( Boolean, include_syn_chi  )
OPT_KEY( Boolean, sampler_allow_syn_pyrimidine )
OPT_KEY( IntegerVector, native_virtual_res  )
OPT_KEY( Real, whole_struct_cluster_radius )
OPT_KEY( Real, suite_cluster_radius )
OPT_KEY( Real, loop_cluster_radius )
OPT_KEY( StringVector, alignment_res )
OPT_KEY( Boolean, filter_user_alignment_res )
OPT_KEY( IntegerVector, native_alignment_res )
OPT_KEY( StringVector, jump_point_pairs )
OPT_KEY( Boolean, floating_base )
OPT_KEY( Boolean, parin_favorite_output )
OPT_KEY( IntegerVector, sample_res )
OPT_KEY( IntegerVector, input_res )
OPT_KEY( IntegerVector, input_res2 )
OPT_KEY( IntegerVector, missing_res )
OPT_KEY( IntegerVector, missing_res2 )
OPT_KEY( IntegerVector, cutpoint_open )
OPT_KEY( Integer, cutpoint_closed )
OPT_KEY( IntegerVector, fixed_res )
OPT_KEY( IntegerVector, minimize_res )
OPT_KEY( IntegerVector, virtual_res )
OPT_KEY( IntegerVector, bulge_res )
OPT_KEY( IntegerVector, terminal_res )
OPT_KEY( IntegerVector, rmsd_res )
OPT_KEY( Boolean, centroid_screen )
OPT_KEY( Boolean, allow_base_pair_only_centroid_screen )
OPT_KEY( Boolean, VDW_atr_rep_screen )
OPT_KEY( Boolean, sampler_perform_o2star_pack )
OPT_KEY( Boolean, fast )
OPT_KEY( Boolean, medium_fast )
OPT_KEY( Boolean, allow_bulge_at_chainbreak )
OPT_KEY( Boolean, VERBOSE )
OPT_KEY( Boolean, sampler_native_rmsd_screen )
OPT_KEY( Real, sampler_native_screen_rmsd_cutoff )
OPT_KEY( Real, score_diff_cut )
OPT_KEY( Boolean, clusterer_perform_score_diff_cut )
OPT_KEY( Integer, sampler_num_pose_kept)
OPT_KEY( Integer, clusterer_num_pose_kept)
OPT_KEY( Boolean, recreate_silent_struct )
OPT_KEY( Boolean, clusterer_use_best_neighboring_shift_RMSD )
OPT_KEY( Boolean, allow_chain_boundary_jump_partner_right_at_fixed_BP )
OPT_KEY( Boolean, allow_fixed_res_at_moving_res )
OPT_KEY( Boolean, clusterer_rename_tags )
OPT_KEY( Boolean, simple_append_map )
OPT_KEY( IntegerVector, global_sample_res_list )
OPT_KEY( IntegerVector, force_syn_chi_res_list )
OPT_KEY( IntegerVector, force_north_ribose_list )
OPT_KEY( IntegerVector, force_south_ribose_list )
OPT_KEY( IntegerVector, protonated_H1_adenosine_list )
OPT_KEY( Boolean,  output_pdb )
OPT_KEY( String, 	start_silent)
OPT_KEY( String, 	start_tag)
OPT_KEY( Boolean,  simple_full_length_job_params )
OPT_KEY( Real, sampler_cluster_rmsd )
OPT_KEY( Boolean, 	output_extra_RMSDs)
OPT_KEY( Boolean, 	integration_test)
OPT_KEY( Boolean, 	add_virt_root ) //For Fang's electron density code.

OPT_KEY( Real, stddev_small )
OPT_KEY( Real, stddev_large )
OPT_KEY( Real, output_score_cutoff )
OPT_KEY( Real, kT )
OPT_KEY( Real, unfolded_weight )
OPT_KEY( Integer, n_sample )
OPT_KEY( Integer, output_period )
OPT_KEY( Boolean, skip_randomize )
OPT_KEY( Boolean, sample_all_o2star )
OPT_KEY( Boolean, do_add_delete )
OPT_KEY( Boolean, presample_added_residue )
OPT_KEY( Integer, presample_internal_cycles )
OPT_KEY( Boolean, start_added_residue_in_aform )
OPT_KEY( Boolean, skip_delete )
OPT_KEY( Boolean, allow_deletion_of_last_residue )

using namespace protocols::swa::monte_carlo;


//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
// Copied from swa_rna_main
utility::vector1< core::Size >
get_fixed_res(core::Size const nres){

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace protocols::swa::rna;

	utility::vector1< Size > actual_fixed_res_list;
	actual_fixed_res_list.clear();

	utility::vector1< core::Size > const fixed_res_list = option[ fixed_res  ]();
	utility::vector1< core::Size > const minimize_res_list= option[ minimize_res ]();

	if(fixed_res_list.size()!=0 && minimize_res_list.size()!=0 ){
		utility_exit_with_message( "User Cannot specify both  fixed_res and minimize_res!" );
	}


	if( fixed_res_list.size()!=0  ){
		actual_fixed_res_list=fixed_res_list;

	}else if( minimize_res_list.size()!=0){

		for(Size seq_num=1; seq_num<=nres; seq_num++){
			if( minimize_res_list.has_value( seq_num) ) continue;
			actual_fixed_res_list.push_back(seq_num);
		}

	}else{ //here I am being a little stringent and require user specify one of these option. Could just return empty list...
		utility_exit_with_message( "User did not specify both fixed res and minimize_res!" );
	}

	return actual_fixed_res_list;
}


//////////////////////////////////////////////////////////////////////////////////////
// Copied from swa_rna_main
utility::vector1< core::Size >
get_input_res(core::Size const nres , std::string const pose_num){


	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace protocols::swa::rna;

	utility::vector1< core::Size > input_res_list;
	utility::vector1< core::Size > missing_res_list;

	if(pose_num=="1"){
		input_res_list= option[ input_res ]();
		missing_res_list= option[ missing_res ]();
	}else if(pose_num=="2"){
		input_res_list= option[ input_res2 ]();
		missing_res_list= option[ missing_res2 ]();
	}else{
		utility_exit_with_message( "Invalid pose_num " + pose_num + ", must by either 1 or 2 !" );
	}


	if( input_res_list.size()!=0 && missing_res_list.size()!=0 ){
		utility_exit_with_message( "User Cannot specify both input_res" + pose_num + " and missing_res" + pose_num + "!" );
	}

	utility::vector1< core::Size > actual_input_res_list;
	actual_input_res_list.clear();

	if( input_res_list.size()!=0){
		actual_input_res_list=input_res_list;

	}else if( missing_res_list.size()!=0){

		for(Size seq_num=1; seq_num<=nres; seq_num++){
			if( missing_res_list.has_value( seq_num) ) continue;
			actual_input_res_list.push_back(seq_num);
		}

	}else{ //did not specify both input_res and missing_res, return empty list
		std::cout << "user did not specify both input_res" << pose_num << " and missing_res" << pose_num << std::endl;
	}

	return actual_input_res_list;

}


///////////////////////////////////////////////////////////////////////////
// Copied from swa_rna_main
core::scoring::ScoreFunctionOP
create_scorefxn(){

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::scoring;


	std::string score_weight_file;

	Size num_score_weight_file=0;

	if ( option[ basic::options::OptionKeys::score::weights ].user() ) {
		score_weight_file= option[ basic::options::OptionKeys::score::weights ]();
		std::cout << "User passed in score:weight option: " << score_weight_file << std::endl;
		num_score_weight_file++;
	}


	if(num_score_weight_file==0){
		//rna_loop_hires_04092010.wts is same as 5X_linear_quarter_fa_stack_and_adjust_bulge_ss_benchmark.wts
		//change default from single_strand_benchmark to 5X_linear_chainbreak_single_strand_benchmark on May 24, 2010
		//change default to 5X_linear_quarter_fa_stack_and_adjust_bulge_ss_benchmark.wts" on April 9th, 2011
		//score_weight_file="rna_loop_hires_04092010.wts";
		utility_exit_with_message("User to need to pass in score:weights"); //Remove the default weight on Sept 28, 2011 Parin S.
	}

	if(num_score_weight_file>1){
		std::cout << "num_score_weight_file (inputted by user)=" << num_score_weight_file << std::endl;
		utility_exit_with_message("num_score_weight_file>1");
	}

	core::scoring::ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( score_weight_file );


	std::cout << "---------score function weights----------" << std::endl;
	scorefxn->show(std::cout);
	std::cout << "-----------------------------------------" << std::endl;


	return scorefxn;
}


void
setup_copy_DOF_input(protocols::swa::rna::StepWiseRNA_PoseSetupOP & stepwise_rna_pose_setup){

	/////////////////////////////////////////////////////////////////////////////////////////
	// StepWisePoseSetup should create the starting pose.
	// This class might eventually be united with the protein StepWisePoseSetup.
	utility::vector1< std::string > input_tags;
	utility::vector1< std::string > silent_files_in;

	if ( option[ in::file::s ].user() ) {
		// Then any pdbs that need to be read in from disk.
		utility::vector1< std::string > const	pdb_tags_from_disk( option[ in::file::s ]() );
		for ( Size n = 1; n <= pdb_tags_from_disk.size(); n++ ) {
			input_tags.push_back( pdb_tags_from_disk[ n ] );
		}
	}

	if(input_tags.size() > 2 ){
		utility_exit_with_message( "input_tags.size() > 2!!" );
	}

	std::cout << "Input structures for COPY DOF" << std::endl;
	for(Size n=1; n<=input_tags.size(); n++){
		if(n<=silent_files_in.size()){
			std::cout << "silent_file tag= " << input_tags[n] << " silent_file= " << silent_files_in[n] << std::endl;
		}else{
			std::cout << "input_tag= " << input_tags[n] << std::endl;
		}
	}

	//////////////////////////////////////////////////////////////////////////////////////////
	stepwise_rna_pose_setup->set_input_tags( input_tags);
	stepwise_rna_pose_setup->set_silent_files_in( silent_files_in);


}


//////////////////////////////////////////////////////////////////////////////////////
// Copied from swa_rna_main -- removed stuff to check silent files.
// probably don't need to set all these options... anyway.
//
protocols::swa::rna::StepWiseRNA_JobParametersOP
setup_rna_job_parameters(){

	using namespace protocols::swa::rna;
	using namespace ObjexxFCL;
	///////////////////////////////
	// Read in sequence.
	if ( !option[ in::file::fasta ].user() ) utility_exit_with_message( "Must supply in::file::fasta!" );
	std::string const fasta_file = option[ in::file::fasta ]()[1];
	core::sequence::SequenceOP fasta_sequence = core::sequence::read_fasta_file( fasta_file )[1];
	std::string const full_sequence = fasta_sequence->sequence();
	core::Size const nres=full_sequence.length();

	if ( !option[ sample_res ].user() ) utility_exit_with_message( "Must supply sample_res!" );


	/////////////////////////////////////////////////////

	StepWiseRNA_JobParametersSetup stepwise_rna_job_parameters_setup( option[ sample_res ](), /*the first element of moving_res_list is the sampling_res*/
																								 										 full_sequence,
																								 										 get_input_res(nres, "1" ),
																								 										 get_input_res(nres, "2" ),
																								 										 option[ cutpoint_open ](),
																								 										 option[ cutpoint_closed ]() );
	stepwise_rna_job_parameters_setup.set_simple_append_map( option[ simple_append_map]() );
	stepwise_rna_job_parameters_setup.set_allow_fixed_res_at_moving_res( option[ allow_fixed_res_at_moving_res ]() ); //Hacky just to get Hermann Duplex working. Need to called before set_fixed_res

	utility::vector1< Size > fixed_res_ = get_fixed_res(nres);
	stepwise_rna_job_parameters_setup.set_fixed_res( fixed_res_ );
	stepwise_rna_job_parameters_setup.set_terminal_res( option[ terminal_res ]() );
	stepwise_rna_job_parameters_setup.set_rmsd_res_list( option[ rmsd_res ]() );
	stepwise_rna_job_parameters_setup.set_jump_point_pair_list( option[ jump_point_pairs ]() ); //Important!: Need to be called after set_fixed_res
	stepwise_rna_job_parameters_setup.set_filter_user_alignment_res( option[ filter_user_alignment_res ]() );

	utility::vector1< std::string > alignment_res_; //why is this a string vector?????
	if ( option[ alignment_res ].user() ) {
		alignment_res_ = option[ alignment_res ]();
	} else {
		for ( Size n = 1; n <= fixed_res_.size(); n++ ) alignment_res_.push_back( string_of( fixed_res_[n] ) );
	}
	stepwise_rna_job_parameters_setup.set_alignment_res( alignment_res_ );

	if ( option[ native_alignment_res ].user() ) {
		stepwise_rna_job_parameters_setup.set_native_alignment_res( option[ native_alignment_res ]() );
	} else {
		stepwise_rna_job_parameters_setup.set_native_alignment_res( fixed_res_ );
	}

	stepwise_rna_job_parameters_setup.set_global_sample_res_list( option[ global_sample_res_list ]() ); //March 20, 2011

	stepwise_rna_job_parameters_setup.set_force_syn_chi_res_list( option[ force_syn_chi_res_list]() ); //April 29, 2011
	stepwise_rna_job_parameters_setup.set_force_north_ribose_list( option[ force_north_ribose_list ]() ); //April 29, 2011
	stepwise_rna_job_parameters_setup.set_force_south_ribose_list( option[ force_south_ribose_list ]() ); //April 29, 2011
	stepwise_rna_job_parameters_setup.set_protonated_H1_adenosine_list( option[ protonated_H1_adenosine_list ]() ); //May 02, 2011

	stepwise_rna_job_parameters_setup.set_allow_chain_boundary_jump_partner_right_at_fixed_BP( option[ allow_chain_boundary_jump_partner_right_at_fixed_BP ]() ); //Hacky just to get Square RNA working.

	stepwise_rna_job_parameters_setup.set_output_extra_RMSDs( option[ output_extra_RMSDs ]() );
	stepwise_rna_job_parameters_setup.set_add_virt_res_as_root( option[ add_virt_root]() );


	stepwise_rna_job_parameters_setup.set_skip_complicated_stuff( true ); // new by Rhiju.

	stepwise_rna_job_parameters_setup.apply();

	return stepwise_rna_job_parameters_setup.job_parameters();

}



//////////////////////////////////////////////////////////////////////////////////////////
// Copied from swa_rna_main
protocols::swa::rna::StepWiseRNA_PoseSetupOP
setup_pose_setup_class(protocols::swa::rna::StepWiseRNA_JobParametersOP & job_parameters, bool const copy_DOF=true){

  using namespace core::pose;
  using namespace core::chemical;
  using namespace core::kinematics;
  using namespace core::scoring;
	using namespace protocols::swa::rna;

	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( RNA );

	// Read in native_pose.
	PoseOP native_pose;
	if (option[ in::file::native ].user() ) {
		native_pose = PoseOP( new Pose );
		import_pose::pose_from_pdb( *native_pose, *rsd_set, option[ in::file::native ]() );
		std::cout << "native_pose->fold_tree(): " << native_pose->fold_tree();
		std::cout << "native_pose->annotated_sequence(true): " << native_pose->annotated_sequence( true ) << std::endl;
		protocols::rna::make_phosphate_nomenclature_matches_mini( *native_pose);
	}

	StepWiseRNA_PoseSetupOP stepwise_rna_pose_setup = new StepWiseRNA_PoseSetup( job_parameters);
	stepwise_rna_pose_setup->set_copy_DOF(copy_DOF);

	if(copy_DOF==true){
		setup_copy_DOF_input(stepwise_rna_pose_setup);
	}


	stepwise_rna_pose_setup->set_virtual_res( option[ virtual_res ]() );
	stepwise_rna_pose_setup->set_bulge_res( option[ bulge_res ]() );
	stepwise_rna_pose_setup->set_native_pose( native_pose );
	stepwise_rna_pose_setup->set_native_virtual_res( option[ native_virtual_res]() );
	stepwise_rna_pose_setup->set_rebuild_bulge_mode( option[rebuild_bulge_mode]() );
	stepwise_rna_pose_setup->set_output_pdb( option[ output_pdb ]() );
	stepwise_rna_pose_setup->set_apply_virtual_res_variant_at_dinucleotide( false );
	stepwise_rna_pose_setup->set_align_to_native( true );

	return stepwise_rna_pose_setup;
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
swa_rna_sample()
{

  using namespace core::pose;
  using namespace core::chemical;
  using namespace core::kinematics;
  using namespace core::scoring;
  using namespace core::io::silent;
  using namespace core::pose::full_model_info;
	using namespace protocols::swa::rna;
	using namespace protocols::moves;

	clock_t const time_start( clock() );

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// This is for RNA, for now.
	ResidueTypeSetCAP rsd_set  = core::chemical::ChemicalManager::get_instance()->residue_type_set( RNA );

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Scorefunction -- choose a default. Put into its own setup function?
	core::scoring::ScoreFunctionOP scorefxn;
	if ( option[ score::weights ].user() ) scorefxn = getScoreFunction();
	else scorefxn = ScoreFunctionFactory::create_score_function( "rna/rna_hires_07232011_with_intra_base_phosphate.wts" ); // Parin's latest weights.

	if ( !scorefxn->has_nonzero_weight( unfolded ) ) { 	// must have unfolded term!
		std::cout << "Putting 'unfolded' term into scorefunction! Use -unfolded_weight to reduce weight or turn off." << std::endl;
		scorefxn->set_weight( unfolded, 1.0 );
	}
	if ( option[ unfolded_weight ].user() ) scorefxn->set_weight( unfolded,  option[ unfolded_weight ]());

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	//  Pose setup -- shared with SWA stuff, for now. Gets native_pose, sample_res, etc. -- put into its own little function?
	StepWiseRNA_JobParametersOP	job_parameters = setup_rna_job_parameters(); // note -- hacked this to include option skip_complicated_stuff
	StepWiseRNA_PoseSetupOP stepwise_rna_pose_setup = setup_pose_setup_class(job_parameters);
	stepwise_rna_pose_setup->set_align_to_native( true );

  Pose pose;
	stepwise_rna_pose_setup->apply( pose );
	stepwise_rna_pose_setup->setup_native_pose( pose ); //NEED pose to align native_pose to pose.

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// graphics!
	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 800, 800 );

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Monte Carlo machinery
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// put this into its own little function?
	std::string const & full_sequence = job_parameters->full_sequence();
	utility::vector1< Size > const start_moving_res_list = job_parameters->working_moving_res_list();
	FullModelInfoOP full_model_info_op =
		new FullModelInfo(  job_parameters->working_res_list(),
												start_moving_res_list,
												full_sequence,
												option[ cutpoint_open ]() );
	pose.data().set( core::pose::datacache::CacheableDataType::FULL_MODEL_INFO, full_model_info_op );

	// put this into its own little function?
	// Set up Movers that go into Main Loop (RNA_SWA_MonteCarloMover). This could also go into RNA_SWA_MonteCarloMover. Hmm.
	Real const kT_ = option[ kT ]();
	Real const sample_range_small = option[ stddev_small ]();
	Real const sample_range_large = option[ stddev_large ]();

	RNA_TorsionMoverOP rna_torsion_mover = new RNA_TorsionMover;

	RNA_DeleteMoverOP rna_delete_mover = new RNA_DeleteMover;

	RNA_AddMoverOP rna_add_mover = new RNA_AddMover( rsd_set, scorefxn );
	rna_add_mover->set_start_added_residue_in_aform( option[ start_added_residue_in_aform ]() );
	rna_add_mover->set_presample_added_residue(  option[ presample_added_residue ]() );
	rna_add_mover->set_internal_cycles( option[ presample_internal_cycles ]() );
	rna_add_mover->set_sample_range_small( sample_range_small );
	rna_add_mover->set_sample_range_large( sample_range_large );
	rna_add_mover->set_kT( kT_ );

	RNA_AddOrDeleteMoverOP rna_add_or_delete_mover = new RNA_AddOrDeleteMover( rna_add_mover, rna_delete_mover );
	rna_add_or_delete_mover->set_allow_deletion_of_last_residue( option[ allow_deletion_of_last_residue ]() );

	RNA_O2StarMoverOP rna_o2star_mover = new RNA_O2StarMover( scorefxn, option[ sample_all_o2star ](), sample_range_small, sample_range_large );

	RNA_SWA_MonteCarloMoverOP rna_swa_montecarlo_mover = new RNA_SWA_MonteCarloMover(  rna_add_or_delete_mover, rna_torsion_mover, rna_o2star_mover, scorefxn );
	rna_swa_montecarlo_mover->set_native_pose( stepwise_rna_pose_setup->get_native_pose() );
	rna_swa_montecarlo_mover->set_silent_file( option[ out::file::silent ]() );
	rna_swa_montecarlo_mover->set_output_period( option[ output_period ]() );
	rna_swa_montecarlo_mover->set_num_cycles( option[ n_sample ]() );
	rna_swa_montecarlo_mover->set_do_add_delete( option[ do_add_delete ]() );

	// put into its own function?
	if ( ! option[ skip_delete ]() && option[ do_add_delete ]() ) rna_delete_mover->wipe_out_moving_residues( pose );

	rna_swa_montecarlo_mover->apply( pose );

	if ( option[ out::file::o].user()	) pose.dump_pdb( option[ out::file::o ]() );

	std::cout << "Total time for monte carlo: " << static_cast<Real>( clock() - time_start ) / CLOCKS_PER_SEC << " seconds." << std::endl;


}



///////////////////////////////////////////////////////////////
void*
my_main( void* )
{

	clock_t const my_main_time_start( clock() );

	swa_rna_sample();

	protocols::viewer::clear_conformation_viewers();

	std::cout << "Total time to run " << static_cast<Real>( clock() - my_main_time_start ) / CLOCKS_PER_SEC << " seconds." << std::endl;

  exit( 0 );

}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{

	try {

  using namespace basic::options;

	utility::vector1< Size > blank_size_vector;
	utility::vector1< std::string > blank_string_vector;

	////////////////////////////////////////////////////
	// should be able to get rid of all the following,
	// once job setup is in its own .cc file, shared by
	// swa, swa_monte_carlo.
	////////////////////////////////////////////////////

	//////////////General/////////////////////////////
	NEW_OPT( graphic, "Turn graphic on/off", true);
	NEW_OPT( Real_parameter_one, "free_variable for testing purposes ", 0.0);
	NEW_OPT( distinguish_pucker, "distinguish pucker when cluster:both in sampler and clusterer", true);
	NEW_OPT( output_pdb, "output_pdb: If true, then will dump the pose into a PDB file at different stages of the stepwise assembly process.", false); //Sept 24, 2011

	//////////////Job_Parameters///////////
	NEW_OPT( sample_res, "residues to build, the first element is the actual sample res while the other are the bulge residues", blank_size_vector );
	NEW_OPT( input_res, "Residues already present in starting pose_1", blank_size_vector );
	NEW_OPT( input_res2, "Residues already present in starting  pose_2", blank_size_vector );
	NEW_OPT( missing_res, "Residues missing in starting pose_1, alternative to input_res", blank_size_vector );
	NEW_OPT( missing_res2, "Residues missing in starting pose_2, alternative to input_res2", blank_size_vector );
	NEW_OPT( rmsd_res, "residues that will be use to calculate rmsd (for clustering as well as RMSD to native_pdb if specified)", blank_size_vector );
	NEW_OPT( alignment_res , "align_res_list", blank_string_vector ); // can this be a size vector now?
	NEW_OPT( global_sample_res_list, "A list of all the nucleotide to be build/sample over the entire dag.", blank_size_vector); //March 20, 2011

	NEW_OPT( cutpoint_open, "optional: chainbreak in full sequence", blank_size_vector );
	NEW_OPT( cutpoint_closed, "optional: cutpoint at which to apply chain closure", 0 );
	NEW_OPT( jump_point_pairs , "optional: extra jump_points specified by the user for setting up the fold_tree ", blank_string_vector );

	NEW_OPT( native_virtual_res , " optional: native_virtual_res ", blank_size_vector );
	NEW_OPT( native_alignment_res , "optional: native_alignment_res ", blank_size_vector );
	NEW_OPT( fixed_res, "optional: residues to be held fixed in minimizer", blank_size_vector );
	NEW_OPT( minimize_res, "optional: residues to be minimize in minimizer, alternative to fixed_res", blank_size_vector );
	NEW_OPT( virtual_res, "optional: residues to be made virtual", blank_size_vector );
	NEW_OPT( terminal_res, "optional: residues that are not allowed to stack during sampling", blank_size_vector );
	NEW_OPT( bulge_res, "optional: residues to be turned into a bulge variant", blank_size_vector );
	NEW_OPT( force_syn_chi_res_list, "optional: sample only syn chi for the res in sampler.", blank_size_vector); //April 29, 2011
	NEW_OPT( force_north_ribose_list, "optional: sample only north ribose for the res in sampler.", blank_size_vector); //April 29, 2011
	NEW_OPT( force_south_ribose_list, "optional: sample only south ribose for the res in sampler.", blank_size_vector); //April 29, 2011
	NEW_OPT( protonated_H1_adenosine_list, "optional: protonate_H1_adenosine_list", blank_size_vector); //May 02, 2011

	//////////////Pose setup///////
	NEW_OPT( job_queue_ID, " swa_rna_sample()/combine_long_loop mode: Specify the tag pair in filter_output_filename to be read in and imported (start from 0!)" , 0);

	///////////////Sampler////////////
	NEW_OPT( rebuild_bulge_mode, "rebuild_bulge_mode", false);
	NEW_OPT( floating_base , " floating_base ", false ); //DO NOT CHANGE TO TRUE, since single-nucleotide sampling need this to be false! April 9th, 2011

	//////////////CombineLongLoopFilterer/////////////
	NEW_OPT( filter_output_filename, "CombineLongLoopFilterer: filter_output_filename", "filter_struct.txt"); //Sept 12, 2010
	NEW_OPT( filter_for_previous_contact, "CombineLongLoopFilterer: filter_for_previous_contact", false); //Sept 12, 2010
	NEW_OPT( filter_for_previous_clash, "CombineLongLoopFilterer: filter_for_previous_clash", false); //Sept 12, 2010
	NEW_OPT( combine_helical_silent_file, "CombineLongLoopFilterer: combine_helical_silent_file", false); //Nov 27, 2010

	//////////////post_rebuild_bulge_assembly//////
	NEW_OPT( start_silent, "start_silent", ""); //Oct 22, 2011
	NEW_OPT( start_tag, "start_tag", ""); //Oct 22, 2011

	///////The options below are for testing purposes. Please do not make any changes without first consulting/////////////
	///////Parin Sripakdeevong (sripakpa@stanford.edu) or Rhiju Das (rhiju@stanford.edu) //////////////////////////////////
	//////////////General/////////////////////////////
	NEW_OPT( VERBOSE, "VERBOSE", false );
	NEW_OPT( parin_favorite_output , " parin_favorite_output ", true ); //Change to true on Oct 10, 2010
	NEW_OPT( integration_test , " integration_test ", false ); //March 16, 2012


	//////////////Job_Parameters///////////
	NEW_OPT( filter_user_alignment_res, " filter_user_alignment_res ", true ); //General want this to be true except for special cases! June 13, 2011
	NEW_OPT( simple_append_map , "simple_append_map", false);
	NEW_OPT( add_virt_root, "add_virt_root", false); //For Fang's electron density code.
	NEW_OPT( allow_chain_boundary_jump_partner_right_at_fixed_BP, "mainly just to get Hermann nano-square RNA modeling to work", false);
	NEW_OPT( allow_fixed_res_at_moving_res, "mainly just to get Hermann Duplex modeling to work", false); //Nov 15, 2010
	NEW_OPT( simple_full_length_job_params, "simple_full_length_job_params", false); //Oct 31, 2011
	NEW_OPT( output_extra_RMSDs, "output_extra_RMSDs", false); //March 16, 2012

	///////////////Sampler////////////
	NEW_OPT( sampler_cluster_rmsd, " Clustering rmsd of conformations in the sampler", 0.5); //DO NOT CHANGE THIS!
	NEW_OPT( skip_sampling, "no sampling step in rna_swa residue sampling", false );
	NEW_OPT( do_not_sample_multiple_virtual_sugar, " Samplerer: do_not_sample_multiple_virtual_sugar " , false);
	NEW_OPT( sample_ONLY_multiple_virtual_sugar, " Samplerer: sample_ONLY_multiple_virtual_sugar " , false);
	NEW_OPT( filterer_undercount_ribose_rotamers, "Undercount all ribose_rotamers as 1 count", false); //July 29, 2011
	NEW_OPT( exclude_alpha_beta_gamma_sampling, "Speed up the debug eplison south sugar mode", false);
	NEW_OPT( debug_eplison_south_sugar_mode, "Check why when eplison is roughly -160 and pucker is south, energy is not favorable", false);
	NEW_OPT( sampler_extra_anti_chi_rotamer, "Samplerer: extra_anti_chi_rotamer", false);
	NEW_OPT( sampler_extra_syn_chi_rotamer, "Samplerer: extra_syn_chi_rotamer", false);
	NEW_OPT( sampler_extra_beta_rotamer, "Samplerer: extra_beta_rotamer", false);
	NEW_OPT( sampler_extra_epsilon_rotamer, "Samplerer: extra_epsilon_rotamer", true); //Change this to true on April 9, 2011
	NEW_OPT( sample_both_sugar_base_rotamer, "Samplerer: Super hacky for SQAURE_RNA", false);
	NEW_OPT( reinitialize_CCD_torsions, "Samplerer: reinitialize_CCD_torsions: Reinitialize_CCD_torsion to zero before every CCD chain closure", false);
	NEW_OPT( PBP_clustering_at_chain_closure, "Samplerer: PBP_clustering_at_chain_closure", false);
	NEW_OPT( finer_sampling_at_chain_closure, "Samplerer: finer_sampling_at_chain_closure", false); //Jun 9, 2010
	NEW_OPT( sampler_include_torsion_value_in_tag, "Samplerer:include_torsion_value_in_tag", true);
	NEW_OPT( include_syn_chi, "include_syn_chi", true); //Change to true on Oct 10, 2010
	NEW_OPT( sampler_allow_syn_pyrimidine, "sampler_allow_syn_pyrimidine", false); //Nov 15, 2010
	NEW_OPT( fast, "quick runthrough for debugging", false );
	NEW_OPT( medium_fast, "quick runthrough for debugging (keep more poses and not as fast as fast option)", false );
	NEW_OPT( centroid_screen, "centroid_screen", true);
	NEW_OPT( allow_base_pair_only_centroid_screen, "allow_base_pair_only_centroid_screen", false); //This only effect floating base sampling + dinucleotide.. deprecate option
	NEW_OPT( sampler_perform_o2star_pack, "perform O2' hydrogen packing inside StepWiseRNA_ResidueSampler", true );
	NEW_OPT( allow_bulge_at_chainbreak, "Allow sampler to replace chainbreak res with virtual_rna_variant if it looks have bad fa_atr score.", true );

	NEW_OPT( add_lead_zero_to_tag, "Add lead zero to clusterer output tag ", false);
	NEW_OPT( n_sample, "Sample number for Random sampling", 0 );
	NEW_OPT( output_period, "How often to output structure", 5000 );
	NEW_OPT( stddev_small, "Sampling standard deviation in degree", 5.0 );
	NEW_OPT( stddev_large, "Sampling standard deviation in degree", 40.0 );
	NEW_OPT( output_score_cutoff, "Score cutoff for output to disk", 0.0 );
	NEW_OPT( kT, "kT of simulation in RU", 2.0 );
	NEW_OPT( unfolded_weight, "weight on unfolded term", 1.0 );
	NEW_OPT( skip_randomize, "do not randomize...", false );
	NEW_OPT( sample_all_o2star, "do not focus o2star sampling at residue of interest...", false );
	NEW_OPT( start_added_residue_in_aform, "on add move, take starting configuration to be A-form, not random", false );
	NEW_OPT( do_add_delete, "try add & delete moves...", false );
	NEW_OPT( presample_added_residue, "when adding a residue, do a little monte carlo to try to get it in place", false );
	NEW_OPT( presample_internal_cycles, "when adding a residue, number of monte carlo cycles", 100 );
	NEW_OPT( skip_delete, "normally wipe out all residues before building", false );
	NEW_OPT( allow_deletion_of_last_residue, "in add/delete allow complete erasure of moving residues during sampling", false );


	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////
  // setup
  ////////////////////////////////////////////////////////////////////////////
  core::init::init(argc, argv);


  ////////////////////////////////////////////////////////////////////////////
  // end of setup
  ////////////////////////////////////////////////////////////////////////////

  protocols::viewer::viewer_main( my_main );


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

}



