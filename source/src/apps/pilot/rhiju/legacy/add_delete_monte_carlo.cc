// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
#include <core/chemical/rna/RNA_FittedTorsionInfo.hh>
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
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>
#include <core/import_pose/import_pose.hh>

//////////////////////////////////////////////////////////
#include <protocols/viewer/viewers.hh>
#include <protocols/farna/util.hh>
#include <protocols/farna/RNA_LoopCloser.hh>
#include <protocols/stepwise/legacy/modeler/rna/StepWiseRNA_PoseSetupFromCommandLine.hh>
#include <protocols/stepwise/modeler/rna/StepWiseRNA_OutputData.hh>
#include <protocols/stepwise/modeler/rna/StepWiseRNA_ResidueSampler.hh>
#include <protocols/stepwise/legacy/modeler/rna/StepWiseRNA_PoseSetup.fwd.hh>
#include <protocols/stepwise/legacy/modeler/rna/StepWiseRNA_PoseSetup.hh>
#include <protocols/stepwise/legacy/modeler/rna/StepWiseRNA_WorkingParametersSetup.hh>
#include <protocols/stepwise/modeler/rna/StepWiseWorkingParameters.hh>
#include <protocols/stepwise/modeler/rna/StepWiseWorkingParametersUtil.hh>

#include <core/pose/full_model_info/FullModelInfo.hh>
#include <protocols/stepwise/monte_carlo/rna/RNA_AddMover.hh>
#include <protocols/stepwise/monte_carlo/rna/RNA_DeleteMover.hh>
#include <protocols/stepwise/monte_carlo/rna/RNA_AddOrDeleteMover.hh>
#include <protocols/stepwise/monte_carlo/rna/RNA_FromScratchMover.hh>
#include <protocols/stepwise/monte_carlo/rna/RNA_O2PrimeMover.hh>
#include <protocols/stepwise/monte_carlo/rna/RNA_TorsionMover.hh>
#include <protocols/stepwise/monte_carlo/rna/RNA_AddDeleteMonteCarlo.hh>
#include <protocols/stepwise/monte_carlo/SWA_Move.hh>

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


// A lot of these options should be placed into an 'official' namespace
// and the SWA RNA pose setup should go to its own function
OPT_KEY( Boolean, do_not_sample_multiple_virtual_sugar)
OPT_KEY( Boolean, sample_ONLY_multiple_virtual_sugar)
OPT_KEY( Boolean, skip_modeler)
OPT_KEY( Boolean, skip_clustering)
OPT_KEY( Boolean, minimize_and_score_native_pose)
OPT_KEY( Integer, num_pose_minimize)
OPT_KEY( Boolean, combine_long_loop_mode)
OPT_KEY( Integer, job_queue_ID)
OPT_KEY( String, filter_output_filename)
OPT_KEY( Boolean, filter_for_previous_contact)
OPT_KEY( Boolean, filter_for_previous_clash)
OPT_KEY( Boolean, filterer_undercount_sugar_rotamers)
OPT_KEY( Boolean, combine_helical_silent_file)
OPT_KEY( Boolean, exclude_alpha_beta_gamma_modeler)
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
OPT_KEY( Boolean, finer_modeler_at_chain_closure)
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
OPT_KEY( Boolean, sampler_perform_o2prime_pack )
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
OPT_KEY( IntegerVector, force_north_sugar_list )
OPT_KEY( IntegerVector, force_south_sugar_list )
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
OPT_KEY( Boolean, sample_all_o2prime )
OPT_KEY( Boolean, do_add_delete )
OPT_KEY( Boolean, presample_added_residue )
OPT_KEY( Integer, presample_internal_cycles )
OPT_KEY( Boolean, start_added_residue_in_aform )
OPT_KEY( Boolean, skip_delete )
OPT_KEY( Boolean, disallow_deletion_of_last_residue )

using namespace protocols::stepwise::monte_carlo;

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
	using namespace protocols::stepwise::modeler::rna;
	using namespace protocols::stepwise::modeler::rna::legacy;
	using namespace protocols::stepwise::monte_carlo::rna;
	using namespace protocols::moves;

	clock_t const time_start( clock() );

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// This is for RNA, for now.
	ResidueTypeSetCAP rsd_set  = core::chemical::ChemicalManager::get_instance()->residue_type_set( RNA );

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Scorefunction -- choose a default. Put into its own setup function?
	core::scoring::ScoreFunctionOP scorefxn;
	if ( option[ score::weights ].user() ) scorefxn = get_score_function();
	else scorefxn = ScoreFunctionFactory::create_score_function( "stepwise/rna/farna/rna_hires_07232011_with_intra_base_phosphate.wts" ); // Parin's latest weights.

	if ( !scorefxn->has_nonzero_weight( unfolded ) ) { 	// must have unfolded term!
		std::cout << "Putting 'unfolded' term into scorefunction! Use -unfolded_weight to reduce weight or turn off." << std::endl;
		scorefxn->set_weight( unfolded, 1.0 );
	}
	if ( option[ unfolded_weight ].user() ) scorefxn->set_weight( unfolded,  option[ unfolded_weight ]());

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	//  Pose setup -- shared with SWA stuff, for now. Gets native_pose, sample_res, etc. -- put into its own little function?
	working_parameters::StepWiseWorkingParametersOP	working_parameters = setup_rna_working_parameters(); // note -- hacked this to include option skip_complicated_stuff
	StepWiseRNA_PoseSetupOP stepwise_rna_pose_setup = setup_pose_setup_class(working_parameters);
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
	std::string const & full_sequence = working_parameters->full_sequence();
	utility::vector1< Size > const start_moving_res_list = working_parameters->working_moving_res_list();

	FullModelInfoOP full_model_info_op =	new FullModelInfo( pose, full_sequence, option[ cutpoint_open](), option[ input_res ]()  );
	pose.data().set( core::pose::datacache::CacheableDataType::FULL_MODEL_INFO, full_model_info_op );

	// put this into its own little function?
	// Set up Movers that go into Main Loop (RNA_AddDeleteMonteCarlo). This could also go into RNA_AddDeleteMonteCarlo. Hmm.
	Real const kT_ = option[ kT ]();
	Real const sample_range_small = option[ stddev_small ]();
	Real const sample_range_large = option[ stddev_large ]();

	RNA_TorsionMoverOP rna_torsion_mover = new RNA_TorsionMover;

	RNA_DeleteMoverOP rna_delete_mover = new RNA_DeleteMover;

	RNA_AddMoverOP rna_add_mover = new RNA_AddMover( scorefxn );
	rna_add_mover->set_start_added_residue_in_aform( option[ start_added_residue_in_aform ]() );
	rna_add_mover->set_presample_added_residue(  option[ presample_added_residue ]() );
	rna_add_mover->set_internal_cycles( option[ presample_internal_cycles ]() );
	rna_add_mover->set_sample_range_small( sample_range_small );
	rna_add_mover->set_sample_range_large( sample_range_large );
	rna_add_mover->set_kT( kT_ );

	RNA_AddOrDeleteMoverOP rna_add_or_delete_mover = new RNA_AddOrDeleteMover( rna_add_mover, rna_delete_mover, 0 );
	rna_add_or_delete_mover->set_disallow_deletion_of_last_residue( option[ disallow_deletion_of_last_residue ]() );

	RNA_O2PrimeMoverOP rna_o2prime_mover = new RNA_O2PrimeMover( scorefxn, option[ sample_all_o2prime ](), sample_range_small, sample_range_large );

	RNA_AddDeleteMonteCarloOP rna_swa_montecarlo_mover = new RNA_AddDeleteMonteCarlo(  rna_add_or_delete_mover, rna_torsion_mover, rna_o2prime_mover, scorefxn );
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
	NEW_OPT( terminal_res, "optional: residues that are not allowed to stack during modeler", blank_size_vector );
	NEW_OPT( bulge_res, "optional: residues to be turned into a bulge variant", blank_size_vector );
	NEW_OPT( force_syn_chi_res_list, "optional: sample only syn chi for the res in sampler.", blank_size_vector); //April 29, 2011
	NEW_OPT( force_north_sugar_list, "optional: sample only north sugar for the res in sampler.", blank_size_vector); //April 29, 2011
	NEW_OPT( force_south_sugar_list, "optional: sample only south sugar for the res in sampler.", blank_size_vector); //April 29, 2011
	NEW_OPT( protonated_H1_adenosine_list, "optional: protonate_H1_adenosine_list", blank_size_vector); //May 02, 2011

	//////////////Pose setup///////
	NEW_OPT( job_queue_ID, " swa_rna_sample()/combine_long_loop mode: Specify the tag pair in filter_output_filename to be read in and imported (start from 0!)" , 0);

	///////////////Sampler////////////
	NEW_OPT( rebuild_bulge_mode, "rebuild_bulge_mode", false);
	NEW_OPT( floating_base , " floating_base ", false ); //DO NOT CHANGE TO TRUE, since single-nucleotide modeler need this to be false! April 9th, 2011

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
	NEW_OPT( skip_modeler, "no modeler step in rna_swa residue modeler", false );
	NEW_OPT( do_not_sample_multiple_virtual_sugar, " Samplerer: do_not_sample_multiple_virtual_sugar " , false);
	NEW_OPT( sample_ONLY_multiple_virtual_sugar, " Samplerer: sample_ONLY_multiple_virtual_sugar " , false);
	NEW_OPT( filterer_undercount_sugar_rotamers, "Undercount all sugar_rotamers as 1 count", false); //July 29, 2011
	NEW_OPT( exclude_alpha_beta_gamma_modeler, "Speed up the debug eplison south sugar mode", false);
	NEW_OPT( debug_eplison_south_sugar_mode, "Check why when eplison is roughly -160 and pucker is south, energy is not favorable", false);
	NEW_OPT( sampler_extra_anti_chi_rotamer, "Samplerer: extra_anti_chi_rotamer", false);
	NEW_OPT( sampler_extra_syn_chi_rotamer, "Samplerer: extra_syn_chi_rotamer", false);
	NEW_OPT( sampler_extra_beta_rotamer, "Samplerer: extra_beta_rotamer", false);
	NEW_OPT( sampler_extra_epsilon_rotamer, "Samplerer: extra_epsilon_rotamer", true); //Change this to true on April 9, 2011
	NEW_OPT( sample_both_sugar_base_rotamer, "Samplerer: Super hacky for SQAURE_RNA", false);
	NEW_OPT( reinitialize_CCD_torsions, "Samplerer: reinitialize_CCD_torsions: Reinitialize_CCD_torsion to zero before every CCD chain closure", false);
	NEW_OPT( PBP_clustering_at_chain_closure, "Samplerer: PBP_clustering_at_chain_closure", false);
	NEW_OPT( finer_modeler_at_chain_closure, "Samplerer: finer_modeler_at_chain_closure", false); //Jun 9, 2010
	NEW_OPT( sampler_include_torsion_value_in_tag, "Samplerer:include_torsion_value_in_tag", true);
	NEW_OPT( include_syn_chi, "include_syn_chi", true); //Change to true on Oct 10, 2010
	NEW_OPT( sampler_allow_syn_pyrimidine, "sampler_allow_syn_pyrimidine", false); //Nov 15, 2010
	NEW_OPT( fast, "quick runthrough for debugging", false );
	NEW_OPT( medium_fast, "quick runthrough for debugging (keep more poses and not as fast as fast option)", false );
	NEW_OPT( centroid_screen, "centroid_screen", true);
	NEW_OPT( allow_base_pair_only_centroid_screen, "allow_base_pair_only_centroid_screen", false); //This only effect floating base modeler + dinucleotide.. deprecate option
	NEW_OPT( sampler_perform_o2prime_pack, "perform O2' hydrogen packing inside StepWiseRNA_ResidueSampler", true );
	NEW_OPT( allow_bulge_at_chainbreak, "Allow sampler to replace chainbreak res with virtual_rna_variant if it looks have bad fa_atr score.", true );

	NEW_OPT( add_lead_zero_to_tag, "Add lead zero to clusterer output tag ", false);
	NEW_OPT( n_sample, "Sample number for Random modeler", 0 );
	NEW_OPT( output_period, "How often to output structure", 5000 );
	NEW_OPT( stddev_small, "Sampling standard deviation in degree", 5.0 );
	NEW_OPT( stddev_large, "Sampling standard deviation in degree", 40.0 );
	NEW_OPT( output_score_cutoff, "Score cutoff for output to disk", 0.0 );
	NEW_OPT( kT, "kT of simulation in RU", 2.0 );
	NEW_OPT( unfolded_weight, "weight on unfolded term", 1.0 );
	NEW_OPT( skip_randomize, "do not randomize...", false );
	NEW_OPT( sample_all_o2prime, "do not focus o2prime modeler at residue of interest...", false );
	NEW_OPT( start_added_residue_in_aform, "on add move, take starting configuration to be A-form, not random", false );
	NEW_OPT( do_add_delete, "try add & delete moves...", false );
	NEW_OPT( presample_added_residue, "when adding a residue, do a little monte carlo to try to get it in place", false );
	NEW_OPT( presample_internal_cycles, "when adding a residue, number of monte carlo cycles", 100 );
	NEW_OPT( skip_delete, "normally wipe out all residues before building", false );
	NEW_OPT( disallow_deletion_of_last_residue, "in add/delete allow complete erasure of moving residues during modeler", false );


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
		return -1;
	}

}


