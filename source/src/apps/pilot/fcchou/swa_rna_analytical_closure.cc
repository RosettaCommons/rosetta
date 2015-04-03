// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief


// libRosetta headers
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueMatcher.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueSelector.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.hh>
#include <core/init/init.hh>
#include <core/io/pdb/pose_io.hh>

//////////////////////////////////////////////////
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/option_macros.hh>
#include <protocols/idealize/idealize.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <protocols/viewer/viewers.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/chemical/rna/RNA_FittedTorsionInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/import_pose/import_pose.hh>
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>

#include <numeric/xyzVector.hh>
#include <numeric/conversions.hh>


#include <string>
#include <map>

#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>
#include <core/import_pose/import_pose.hh>


//////////////////////////////////////////////////////////
#include <protocols/stepwise/sampling/rna/StepWiseRNA_CombineLongLoopFilterer.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_CombineLongLoopFilterer.fwd.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_BaseCentroidScreener.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_BaseCentroidScreener.fwd.hh>
#include <protocols/stepwise/legacy/sampling/rna/StepWiseRNA_Minimizer.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_AnalyticalLoopCloseSampler.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_PoseSetup.fwd.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_PoseSetup.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_Clusterer.hh>
#include <protocols/stepwise/sampling/rna/StepWiseWorkingParametersSetup.hh>
#include <protocols/stepwise/sampling/rna/StepWiseWorkingParameters.hh>
#include <protocols/stepwise/StepWiseLegacyClusterer.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_VDW_BinScreener.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_VDW_BinScreener.fwd.hh>
#include <protocols/farna/util.hh>

#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>

#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>

#include <utility/excn/Exceptions.hh>


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
#define GetCurrentDir getcwd

//Added by Parin
//#include <core/scoring/ScoreType.hh>
#include <list>
#include <stdio.h>
#include <math.h>

using namespace core;
using namespace protocols;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace ObjexxFCL;
using utility::vector1;
using io::pdb::dump_pdb;

typedef  numeric::xyzMatrix< Real > Matrix;


OPT_KEY ( Boolean, add_virt_root )
OPT_KEY ( Boolean, skip_sampling )
OPT_KEY ( Boolean, skip_clustering )
OPT_KEY ( Boolean, minimize_and_score_native_pose )
OPT_KEY ( Integer, num_pose_minimize )
OPT_KEY ( Integer, job_queue_ID )
OPT_KEY ( String, filter_output_filename )
OPT_KEY ( Boolean, filter_for_previous_contact )
OPT_KEY ( Boolean, filter_for_previous_clash )
OPT_KEY ( Boolean, minimize_and_score_sugar )
OPT_KEY ( Boolean, sampler_extra_anti_chi_rotamer )
OPT_KEY ( Boolean, sampler_extra_syn_chi_rotamer )
OPT_KEY ( Boolean, PBP_clustering_at_chain_closure )
OPT_KEY ( Boolean, clusterer_two_stage_clustering )
OPT_KEY ( Boolean, clusterer_keep_pose_in_memory )
OPT_KEY ( Boolean, finer_sampling_at_chain_closure )
OPT_KEY ( StringVector, 	VDW_rep_screen_info )
OPT_KEY ( Real, 	VDW_rep_alignment_RMSD_CUTOFF )
OPT_KEY ( Boolean, graphic )
OPT_KEY ( Real, Real_parameter_one )
OPT_KEY ( Boolean, clusterer_quick_alignment )
OPT_KEY ( Boolean, clusterer_optimize_memory_usage )
OPT_KEY ( Integer, clusterer_min_struct )
OPT_KEY ( Boolean, clusterer_write_score_only )
OPT_KEY ( Boolean, add_lead_zero_to_tag )
OPT_KEY ( Boolean, distinguish_pucker )
OPT_KEY ( IntegerVector, native_virtual_res )
OPT_KEY ( Real, suite_cluster_radius )
OPT_KEY ( Real, loop_cluster_radius )
OPT_KEY ( StringVector, alignment_res )
OPT_KEY ( IntegerVector, native_alignment_res )
OPT_KEY ( StringVector, jump_point_pairs )
OPT_KEY ( IntegerVector, sample_res )
OPT_KEY ( IntegerVector, input_res )
OPT_KEY ( IntegerVector, input_res2 )
OPT_KEY ( IntegerVector, missing_res )
OPT_KEY ( IntegerVector, missing_res2 )
OPT_KEY ( IntegerVector, cutpoint_open )
OPT_KEY ( Integer, cutpoint_closed )
OPT_KEY ( IntegerVector, fixed_res )
OPT_KEY ( IntegerVector, minimize_res )
OPT_KEY ( IntegerVector, virtual_res )
OPT_KEY ( IntegerVector, bulge_res )
OPT_KEY ( IntegerVector, terminal_res )
OPT_KEY ( IntegerVector, rmsd_res )
OPT_KEY ( Boolean, prepend )
OPT_KEY ( Boolean, centroid_screen )
OPT_KEY ( Boolean, VDW_atr_rep_screen )
OPT_KEY ( Boolean, sampler_perform_o2prime_pack )
OPT_KEY ( Boolean, fast )
OPT_KEY ( Boolean, medium_fast )
OPT_KEY ( Boolean, VERBOSE )
OPT_KEY ( Boolean, sampler_native_rmsd_screen )
OPT_KEY ( Real, sampler_native_screen_rmsd_cutoff )
OPT_KEY ( Real, native_edensity_score_cutoff )
OPT_KEY ( Boolean, auto_tune )
OPT_KEY ( Boolean, skip_minimize )
OPT_KEY ( Real, score_diff_min )
OPT_KEY ( Real, score_diff_cut )
OPT_KEY ( Real, score_diff_cut_tier_two )
OPT_KEY ( Real, score_diff_cut_tier_three )
OPT_KEY ( String, 	algorithm )
OPT_KEY ( String, 	cluster_type )
OPT_KEY ( Integer, sampler_num_pose_kept )
OPT_KEY ( StringVector, input_tag_list )
OPT_KEY ( Boolean, recreate_silent_struct )
OPT_KEY ( String, FARFAR_start_pdb )
OPT_KEY ( Boolean, allow_chain_boundary_jump_partner_right_at_fixed_BP )
OPT_KEY ( Boolean, allow_fixed_res_at_moving_res )
OPT_KEY( Real, sampler_cluster_rmsd )
OPT_KEY( Boolean,  output_pdb )

//////////////////////////////////////////////////////////////////////////////////////
std::string
get_working_directory() {
	char cCurrentPath[FILENAME_MAX];
	std::string current_directory_string;

	if ( !GetCurrentDir ( cCurrentPath, sizeof ( cCurrentPath ) ) ) {
		utility_exit_with_message ( "!GetCurrentDir(cCurrentPath, sizeof(cCurrentPath))" );
	}

	//cCurrentPath[sizeof(cCurrentPath) - 1] = '/0'; /* not really required */
	//std::cout << "current_directory= " << cCurrentPath << std::endl;
	std::stringstream ss;
	ss << cCurrentPath;
	ss >> current_directory_string;
	std::cout << "current_directory= " << current_directory_string << std::endl;
	return current_directory_string;
}


utility::vector1< core::Size >
get_fixed_res ( core::Size const nres ) {
	using namespace protocols::stepwise::sampling::rna;
	utility::vector1< Size > actual_fixed_res_list;
	actual_fixed_res_list.clear();
	utility::vector1< core::Size > const fixed_res_list = option[ fixed_res  ]();
	utility::vector1< core::Size > const minimize_res_list = option[ minimize_res ]();

	if ( fixed_res_list.size() != 0 && minimize_res_list.size() != 0 ) {
		utility_exit_with_message ( "User Cannot specify both  fixed_res and minimize_res!" );
	}

	if ( fixed_res_list.size() != 0 ) {
		actual_fixed_res_list = fixed_res_list;
	} else if ( minimize_res_list.size() != 0 ) {
		for ( Size seq_num = 1; seq_num <= nres; seq_num++ ) {
			if ( Contain_seq_num ( seq_num, minimize_res_list ) ) continue;

			actual_fixed_res_list.push_back ( seq_num );
		}
	} else { //here I am being a little stringent and require user specify one of these option. Could just return empty list...
		utility_exit_with_message ( "User did not specify both fixed res and minimize_res!" );
	}

	return actual_fixed_res_list;
}

//////////////////////////////////////////////////////////////////////////////////////
bool
is_nonempty_input_silent_file ( std::string const input_silent_file, std::string const exit_key_string ) {
	std::cout << "Checking that input_silent_file " << input_silent_file << " contain actual silent_structs or the correct exit_key_string" << std::endl;
	std::ifstream infile;
	infile.open ( input_silent_file.c_str() );

	if ( infile.fail() ) {
		utility_exit_with_message ( "Error! \"" + input_silent_file + "\" could not be opened!" );
	} else {
		std::cout << "Open \"" << input_silent_file << "\" successful!" << std::endl;
	}

	std::string line;
	bool found_queue_ID = false;
	bool found_line = getline ( infile, line );

	if ( found_line == false ) utility_exit_with_message ( "No line exist in input_silent_file= " + input_silent_file );

	size_t found_substring = line.find ( exit_key_string );

	if ( found_substring != std::string::npos ) {
		std::cout << "input_silent_file: " << input_silent_file << " contain no silent struct" << std::endl;
		std::cout << line << std::endl;
		//consistency check:////////////////////////////////////////////////////////////////////////////////////////////////////
		std::string next_line;
		bool found_next_line = getline ( infile, next_line );

		if ( found_next_line ) std::cout << "input silent_file contain more than one line! next_line= " << next_line << std::endl;

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		return false;
	} else {
		return true;
	}
}
//////////////////////////////////////////////////////////////////////////////////////

utility::vector1< core::Size >
get_input_res ( core::Size const nres , std::string const pose_num ) {
	using namespace protocols::stepwise::sampling::rna;
	utility::vector1< core::Size > input_res_list;
	utility::vector1< core::Size > missing_res_list;

	if ( pose_num == "1" ) {
		input_res_list = option[ input_res ]();
		missing_res_list = option[ missing_res ]();
	} else if ( pose_num == "2" ) {
		input_res_list = option[ input_res2 ]();
		missing_res_list = option[ missing_res2 ]();
	} else {
		utility_exit_with_message ( "Invalid pose_num " + pose_num + ", must by either 1 or 2 !" );
	}

	if ( input_res_list.size() != 0 && missing_res_list.size() != 0 ) {
		utility_exit_with_message ( "User Cannot specify both input_res" + pose_num + " and missing_res" + pose_num + "!" );
	}

	utility::vector1< core::Size > actual_input_res_list;
	actual_input_res_list.clear();

	if ( input_res_list.size() != 0 ) {
		actual_input_res_list = input_res_list;
	} else if ( missing_res_list.size() != 0 ) {
		for ( Size seq_num = 1; seq_num <= nres; seq_num++ ) {
			if ( Contain_seq_num ( seq_num, missing_res_list ) ) continue;

			actual_input_res_list.push_back ( seq_num );
		}
	} else { //did not specify both input_res and missing_res, return empty list
		std::cout << "user did not specify both input_res" << pose_num << " and missing_res" << pose_num << std::endl;
	}

	return actual_input_res_list;
}

//////////////////////////////////////////////////////////////////////////////////////


utility::vector1< std::string >
get_silent_file_tags() {
	using namespace protocols::stepwise::sampling::rna;
	bool tags_from_command_line = false;
	bool tags_from_filterer_outfile = false;
	utility::vector1< std::string > input_silent_file_tags;

	if ( option[ in::file::tags ].user() ) {
		tags_from_command_line = true;
		input_silent_file_tags = option[ in::file::tags ]();
	}

	if ( option[ job_queue_ID ].user() && option[ filter_output_filename ].user() ) {
		output_title_text( "importing tag from filter_outfile", TR );
		tags_from_filterer_outfile = true;
		std::string const filtered_tag_file = option[ filter_output_filename ]();
		std::ifstream infile;
		infile.open ( filtered_tag_file.c_str() );

		if ( infile.fail() ) {
			utility_exit_with_message ( "Error! \"" + filtered_tag_file + "\" could not be opened!" );
		} else {
			std::cout << "Open \"" << filtered_tag_file << "\" successful!" << std::endl;
		}

		//Becareful here... job_queue_ID start from ZERO!
		int const queue_ID = option[ job_queue_ID ]();
		int ID = 0;
		std::cout << "queue_ID= " << queue_ID << std::endl;
		std::string tag_pair_string;
		bool found_queue_ID = false;

		while ( getline ( infile, tag_pair_string ) ) {
			if ( queue_ID == ID ) {
				found_queue_ID = true;
				break;
			}

			ID++;
		}

		//Warning queue_ID start at ZERO!
		if ( found_queue_ID == false ) utility_exit_with_message ( "found_queue_ID==false, queue_ID= " + string_of ( queue_ID ) + " num_tag_string_in_file= " + string_of ( ID ) );

		std::cout << "import silent_file_tags: " << tag_pair_string << " from filter_output_filename= " << filtered_tag_file << std::endl;
		infile.close();
		utility::vector1< std::string > const line_list = Tokenize ( tag_pair_string, " \t\n\f\v" ); //Oct 19, 2010..now filterer_outfile contain other terms.
		input_silent_file_tags.clear();
		input_silent_file_tags.push_back ( line_list[1] );
		input_silent_file_tags.push_back ( line_list[2] );
		output_title_text( "", TR );
	}

	if ( ( tags_from_command_line == false ) && ( tags_from_filterer_outfile == false ) ) {
		utility_exit_with_message ( "(tags_from_command_line==false) && (tags_from_filterer_outfile==false)" );
	}

	if ( ( tags_from_command_line == true ) && ( tags_from_filterer_outfile == true ) ) {
		utility_exit_with_message ( "(tags_from_command_line==true) && (tags_from_filterer_outfile==true)" );
	}

	return input_silent_file_tags;
}


//////////////////////////////////////////////////////////////////////////////////////

core::scoring::ScoreFunctionOP
create_scorefxn() {
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


	if(option[minimize_and_score_sugar]()==false){
		std::cout << "WARNING minimize_and_score_sugar is false, SET rna_sugar_close weight to 0.0 " << std::endl;
    scorefxn->set_weight( rna_sugar_close, 0.000000000000 );
 	}

	std::cout << "---------score function weights----------" << std::endl;
	scorefxn->show(std::cout);
	std::cout << "-----------------------------------------" << std::endl;


	return scorefxn;
}

//////////////////////////////////////////////////////////////////////////////////////

protocols::stepwise::sampling::rna::working_parameters::StepWiseWorkingParametersOP
setup_rna_working_parameters ( bool check_for_previously_closed_cutpoint_with_input_pose = false ) {
	using namespace protocols::stepwise::sampling::rna;
	using namespace ObjexxFCL;
	///////////////////////////////
	// Read in sequence.
	std::string const fasta_file = option[ in::file::fasta ]() [1];
	core::sequence::SequenceOP fasta_sequence = core::sequence::read_fasta_file ( fasta_file ) [1];
	std::string const full_sequence = fasta_sequence->sequence();
	core::Size const nres = full_sequence.length();

	if ( !option[ sample_res ].user() ) utility_exit_with_message ( "Must supply sample_res!" );

	/////////////////////////////////////////////////////
	StepWiseWorkingParametersSetup stepwise_rna_working_parameters_setup ( option[ sample_res ](), /*the first element of moving_res_list is the sampling_res*/
	    full_sequence,
	    get_input_res ( nres, "1" ),
	    get_input_res ( nres, "2" ),
	    option[ cutpoint_open ](),
	    option[ cutpoint_closed ]() );
	stepwise_rna_working_parameters_setup.set_add_virt_res_as_root ( option[ add_virt_root]() );
	stepwise_rna_working_parameters_setup.set_allow_fixed_res_at_moving_res ( option[ allow_fixed_res_at_moving_res ]() ); //Hacky just to get Hermann Duplex working. Need to called before set_fixed_res
	stepwise_rna_working_parameters_setup.set_fixed_res ( get_fixed_res ( nres ) );
	stepwise_rna_working_parameters_setup.set_terminal_res ( option[ terminal_res ]() );
	stepwise_rna_working_parameters_setup.set_rmsd_res_list ( option[ rmsd_res ]() );
	stepwise_rna_working_parameters_setup.set_jump_point_pair_list ( option[ jump_point_pairs ]() ); //Important!: Need to be called after set_fixed_res
	stepwise_rna_working_parameters_setup.set_alignment_res ( option[ alignment_res ]() ); //Important!: Need to be called after set_fixed_res
	stepwise_rna_working_parameters_setup.set_native_alignment_res ( option[ native_alignment_res ]() );
	stepwise_rna_working_parameters_setup.set_FARFAR_start_pdb ( option[ FARFAR_start_pdb ]() );
	stepwise_rna_working_parameters_setup.set_allow_chain_boundary_jump_partner_right_at_fixed_BP ( option[ allow_chain_boundary_jump_partner_right_at_fixed_BP ]() ); //Hacky just to get Square RNA working.
	/////////////////////////////Sept 1, 2010////////////
	if ( check_for_previously_closed_cutpoint_with_input_pose ) {
		utility::vector1< std::string > input_tags;
		utility::vector1< std::string > silent_files_in;

		// First read in any information on pdb read in from silent files.
		// Assume one to one correspondence between number of tags and number of silent_file
		if ( option[ in::file::silent ].user() ) {
			silent_files_in = option[ in::file::silent ]();
			input_tags = get_silent_file_tags();

			if ( silent_files_in.size() != input_tags.size() ) {
				utility_exit_with_message ( "silent_files_in.size(" + string_of ( silent_files_in.size() ) + ")!=input_tags.size(" + string_of ( input_tags.size() ) + ")" );
			}
		}

		stepwise_rna_working_parameters_setup.set_input_tags ( input_tags );
		stepwise_rna_working_parameters_setup.set_silent_files_in ( silent_files_in );
	}

	///////////////////////////////////////////////////////
	stepwise_rna_working_parameters_setup.apply();
	return stepwise_rna_working_parameters_setup.working_parameters();
}


void
setup_copy_DOF_input ( protocols::stepwise::sampling::rna::StepWiseRNA_PoseSetupOP & stepwise_rna_pose_setup ) {
	/////////////////////////////////////////////////////////////////////////////////////////
	// StepWiseProteinPoseSetup should create the starting pose.
	// This class might eventually be united with the protein StepWiseProteinPoseSetup.
	utility::vector1< std::string > input_tags;
	utility::vector1< std::string > silent_files_in;

	// First read in any information on pdb read in from silent files.
	// Assume one to one correspondence between number of tags and number of silent_file
	if ( option[ in::file::silent ].user() ) {
		silent_files_in = option[ in::file::silent ]();
		input_tags = get_silent_file_tags();

		if ( silent_files_in.size() != input_tags.size() ) {
			utility_exit_with_message ( "silent_files_in.size()!=input_tags.size()" );
		}
	}

	if ( option[ in::file::s ].user() ) {
		// Then any pdbs that need to be read in from disk.
		utility::vector1< std::string > const	pdb_tags_from_disk ( option[ in::file::s ]() );

		for ( Size n = 1; n <= pdb_tags_from_disk.size(); n++ ) {
			input_tags.push_back ( pdb_tags_from_disk[ n ] );
		}
	}

	if ( input_tags.size() > 2 ) {
		utility_exit_with_message ( "input_tags.size() > 2!!" );
	}

	std::cout << "Input structures for COPY DOF" << std::endl;

	for ( Size n = 1; n <= input_tags.size(); n++ ) {
		if ( n <= silent_files_in.size() ) {
			std::cout << "silent_file tag= " << input_tags[n] << " silent_file= " << silent_files_in[n] << std::endl;
		} else {
			std::cout << "input_tag= " << input_tags[n] << std::endl;
		}
	}

	//////////////////////////////////////////////////////////////////////////////////////////
	stepwise_rna_pose_setup->set_input_tags ( input_tags );
	stepwise_rna_pose_setup->set_silent_files_in ( silent_files_in );
}

protocols::stepwise::sampling::rna::StepWiseRNA_PoseSetupOP
setup_pose_setup_class(protocols::stepwise::sampling::rna::working_parameters::StepWiseWorkingParametersOP & working_parameters, bool const copy_DOF=true){

  using namespace core::pose;
  using namespace core::chemical;
  using namespace core::kinematics;
  using namespace core::scoring;
	using namespace protocols::stepwise::sampling::rna;

	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( RNA );

	// Read in native_pose.
	PoseOP native_pose;
	if (option[ in::file::native ].user() ) {
		native_pose = PoseOP( new Pose );
		import_pose::pose_from_pdb( *native_pose, *rsd_set, option[ in::file::native ]() );
		std::cout << "native_pose->fold_tree(): " << native_pose->fold_tree();
		std::cout << "native_pose->annotated_sequence(true): " << native_pose->annotated_sequence( true ) << std::endl;
		protocols::farna::make_phosphate_nomenclature_matches_mini( *native_pose);
	}


//	StepWiseRNA_PoseSetup stepwise_rna_pose_setup( pdb_tags, silent_files_in, working_parameters);

	StepWiseRNA_PoseSetupOP stepwise_rna_pose_setup = new StepWiseRNA_PoseSetup( working_parameters);
	stepwise_rna_pose_setup->set_copy_DOF(copy_DOF);

	if(copy_DOF==true){
		setup_copy_DOF_input(stepwise_rna_pose_setup);
	}

	stepwise_rna_pose_setup->set_virtual_res( option[ virtual_res ]() );
	stepwise_rna_pose_setup->set_bulge_res( option[ bulge_res ]() );
	stepwise_rna_pose_setup->set_native_pose( native_pose );
	stepwise_rna_pose_setup->set_native_virtual_res( option[ native_virtual_res]() );
	stepwise_rna_pose_setup->set_FARFAR_start_pdb( option[FARFAR_start_pdb]() );
	stepwise_rna_pose_setup->set_output_pdb( option[ output_pdb ]() );

	return stepwise_rna_pose_setup;
}


void
rna_resample_test() {
	using namespace core::pose;
	using namespace core::chemical;
	using namespace core::kinematics;
	using namespace core::scoring;
	using namespace protocols::stepwise::sampling::rna;
	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set ( RNA );

	if ( option[ minimize_and_score_native_pose]() == true ) {
		if ( option[ in::file::native ].user() == false ) utility_exit_with_message ( "minimize_and_score_native_pose==True but user did not pass in native pose" );
	}

	///////////////////////////////
	working_parameters::StepWiseWorkingParametersOP	working_parameters = setup_rna_working_parameters ( true  /*check_for_previously_closed_cutpoint_with_input_pose */ );
	working_parameters::StepWiseWorkingParametersCOP working_parameters_COP ( working_parameters );
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	StepWiseRNA_PoseSetupOP stepwise_rna_pose_setup = setup_pose_setup_class ( working_parameters );
	Pose pose;
	stepwise_rna_pose_setup->apply ( pose );
	stepwise_rna_pose_setup->setup_native_pose ( pose ); //NEED pose to align native_pose to pose.
	StepWiseRNA_VDW_BinScreenerOP user_input_VDW_bin_screener = new StepWiseRNA_VDW_BinScreener();

	if ( option[ VDW_rep_screen_info].user() ) {
		user_input_VDW_bin_screener->set_VDW_rep_alignment_RMSD_CUTOFF ( option[ VDW_rep_alignment_RMSD_CUTOFF]() );
		user_input_VDW_bin_screener->setup_using_user_input_VDW_pose ( option[ VDW_rep_screen_info](), pose, working_parameters_COP );
		user_input_VDW_bin_screener->set_output_pdb( option[ output_pdb ]() );
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	StepWiseRNA_BaseCentroidScreenerOP base_centroid_screener = new StepWiseRNA_BaseCentroidScreener ( pose, working_parameters_COP );
	//////////////////////////////////////////////////////////////////
	core::scoring::ScoreFunctionOP scorefxn = create_scorefxn();

	if ( option[ graphic ]() ) {
		std::string const current_directory_string = get_working_directory();
		protocols::viewer::add_conformation_viewer ( pose.conformation(), current_directory_string, 400, 400 );
	}

	StepWiseRNA_AnalyticalLoopCloseSampler stepwise_rna_residue_sampler ( working_parameters_COP );
	std::string const silent_file = option[ out::file::silent  ]();
	stepwise_rna_residue_sampler.set_silent_file ( silent_file + "_sampling" );
	stepwise_rna_residue_sampler.set_scorefxn ( scorefxn );
	stepwise_rna_residue_sampler.set_num_pose_kept ( option[ sampler_num_pose_kept ]() );
	stepwise_rna_residue_sampler.set_fast ( option[ fast ]() );
	stepwise_rna_residue_sampler.set_medium_fast ( option[ medium_fast ]() );
	stepwise_rna_residue_sampler.set_native_rmsd_screen ( option[ sampler_native_rmsd_screen ]() );
	stepwise_rna_residue_sampler.set_native_screen_rmsd_cutoff ( option[ sampler_native_screen_rmsd_cutoff ]() );
	stepwise_rna_residue_sampler.set_o2prime_screen ( option[ sampler_perform_o2prime_pack ]() );
	stepwise_rna_residue_sampler.set_verbose ( option[ VERBOSE ]() );
	stepwise_rna_residue_sampler.set_cluster_rmsd (	option[ sampler_cluster_rmsd ]()	);
	stepwise_rna_residue_sampler.set_distinguish_pucker ( option[ distinguish_pucker]() );
	stepwise_rna_residue_sampler.set_finer_sampling_at_chain_closure ( option[ finer_sampling_at_chain_closure]() );
	stepwise_rna_residue_sampler.set_PBP_clustering_at_chain_closure ( option[ PBP_clustering_at_chain_closure]() );
	stepwise_rna_residue_sampler.set_extra_syn_chi_rotamer ( option[ sampler_extra_syn_chi_rotamer]() );
	stepwise_rna_residue_sampler.set_extra_anti_chi_rotamer ( option[ sampler_extra_anti_chi_rotamer]() );
	stepwise_rna_residue_sampler.set_centroid_screen ( option[ centroid_screen ]() );
	stepwise_rna_residue_sampler.set_VDW_atr_rep_screen ( option[ VDW_atr_rep_screen ]() );
	stepwise_rna_residue_sampler.set_base_centroid_screener ( base_centroid_screener );
	stepwise_rna_residue_sampler.set_user_input_VDW_bin_screener ( user_input_VDW_bin_screener );

	if ( option[skip_sampling]() == false ) stepwise_rna_residue_sampler.apply ( pose );

	utility::vector1< pose_data_struct2 > & pose_data_list = stepwise_rna_residue_sampler.get_pose_data_list();

	if ( option[minimize_and_score_native_pose]() == true ) {
		pose::Pose native_pose = ( *working_parameters_COP->working_native_pose() );
		pose_data_struct2 native_data_struct;
		native_data_struct.pose_OP = new pose::Pose;
		( *native_data_struct.pose_OP ) = native_pose;
		native_data_struct.score = 0.0;
		native_data_struct.tag = "working_native_pose";
		pose_data_list.push_back ( native_data_struct );
	}

	// let's output the final pose_data_list, just to have a look
	if ( option[ VERBOSE ]() ) {
		stepwise_rna_residue_sampler.output_pose_data_list ( silent_file + "_final_sample" );
	}

	////////////////////////////////////////////////////////////////
	StepWiseRNA_Minimizer stepwise_rna_minimizer ( stepwise_rna_residue_sampler.get_pose_data_list() , working_parameters_COP );
	stepwise_rna_minimizer.set_silent_file ( silent_file );
	stepwise_rna_minimizer.set_verbose ( option[ VERBOSE ]() );
	stepwise_rna_minimizer.set_scorefxn ( scorefxn );
	stepwise_rna_minimizer.set_centroid_screen ( option[ centroid_screen ]() );
	stepwise_rna_minimizer.set_base_centroid_screener ( base_centroid_screener );
	stepwise_rna_minimizer.set_skip_minimize ( option[ skip_minimize ]() );
	stepwise_rna_minimizer.set_native_rmsd_screen ( option[ sampler_native_rmsd_screen ]() );
	stepwise_rna_minimizer.set_native_edensity_score_cutoff ( option[native_edensity_score_cutoff]() );
	stepwise_rna_minimizer.set_native_screen_rmsd_cutoff ( option[ sampler_native_screen_rmsd_cutoff ]() + 1 ); //+1 for leniency Sept 20, 2010

	if ( option[num_pose_minimize].user() ) 	stepwise_rna_minimizer.set_num_pose_minimize ( option[ num_pose_minimize ]() );

	stepwise_rna_minimizer.set_minimize_and_score_sugar ( option[  minimize_and_score_sugar ]() );
	stepwise_rna_minimizer.set_user_input_VDW_bin_screener ( user_input_VDW_bin_screener );
	stepwise_rna_minimizer.apply ( pose );
}

///////////////////////////////////////////////////////////////
void*
my_main ( void* ) {
	rna_resample_test();
	protocols::viewer::clear_conformation_viewers();
	std::cout << "JOB_SUCCESSFULLY_COMPLETED" << std::endl;
	exit ( 0 );
}


///////////////////////////////////////////////////////////////////////////////
int
main ( int argc, char * argv [] ) {
    try {
    utility::vector1< Size > blank_size_vector;
	utility::vector1< std::string > blank_string_vector;
	NEW_OPT ( add_virt_root, "add_virt_root", false );
	NEW_OPT ( job_queue_ID, " rna_resample_test()/combine_long_loop mode: Specify the tag pair in filter_output_filename to be read in and imported (start from 0!)" , 0 );
	NEW_OPT ( filter_output_filename, "CombineLongLoopFilterer: filter_output_filename", "filter_struct.txt" ); //Sept 12, 2010
	NEW_OPT ( filter_for_previous_contact, "CombineLongLoopFilterer: filter_for_previous_contact", false ); //Sept 12, 2010
	NEW_OPT ( filter_for_previous_clash, "CombineLongLoopFilterer: filter_for_previous_clash", false ); //Sept 12, 2010
//New option Aug 15 2010 //Reinitialize_CCD_torsion to zero before every CCD chain closure
	NEW_OPT ( minimize_and_score_native_pose, "minimize_and_score_native_pose ", false ); //Sept 15, 2010
	NEW_OPT ( minimize_and_score_sugar, "minimize and sugar torsion+angle? and include the rna_sugar_close_score_term ", true ); //Sept 15, 2010
	NEW_OPT ( sampler_extra_syn_chi_rotamer, "Samplerer: extra_chi_rotamer", false );
	NEW_OPT ( sampler_extra_anti_chi_rotamer, "Samplerer: extra_chi_rotamer", false );
	NEW_OPT ( PBP_clustering_at_chain_closure, "Samplerer: PBP_clustering_at_chain_closure", false );
	NEW_OPT ( clusterer_two_stage_clustering, "Cluster is two stage..using triangle inequaility to speed up clustering", true ); //Change to true on Oct 10, 2010
	NEW_OPT ( clusterer_keep_pose_in_memory, "reduce memory usage for the clusterer", true ); //Aug 6, 2010
	NEW_OPT ( finer_sampling_at_chain_closure, "Samplerer: finer_sampling_at_chain_closure", false ); //Jun 9, 2010
	NEW_OPT ( VDW_rep_screen_info, "VDW_rep_screen_info to create VDW_rep_screen_bin (useful when building loop from large poses)", blank_string_vector ); //Jun 9, 2010
	NEW_OPT ( VDW_rep_alignment_RMSD_CUTOFF, "use with VDW_rep_screen_info", 0.001 ); //Nov 12, 2010
	NEW_OPT ( graphic, "Turn graphic on/off", true ); //May 5, 2010
	NEW_OPT ( recreate_silent_struct, "Special mode to recreate_silent_struct for clusterer output...for analysis purposes", false ); //May 5, 2010
	NEW_OPT ( Real_parameter_one, "free_variable for testing purposes ", 0.0 );
	NEW_OPT ( clusterer_quick_alignment, "quick alignment during clusterer...only work if the alignment residues are fixed ", false );
	NEW_OPT ( clusterer_optimize_memory_usage, "clusterer_optimize_memory_usage ", false );
	NEW_OPT ( clusterer_min_struct, "clusterer_min_struct ", 400 ); //Oct 13, 2010
	NEW_OPT ( clusterer_write_score_only, "clusterer_write_score_only/ only effect recreate_silent_struct mode  ", false ); //Oct 20, 2010
	NEW_OPT ( add_lead_zero_to_tag, "Add lead zero to clusterer output tag ", false );
	NEW_OPT ( distinguish_pucker, "distinguish pucker when cluster:both in sampler and clusterer", true );
	NEW_OPT ( input_tag_list, "input_tag_list", blank_string_vector );
	NEW_OPT ( native_virtual_res , " native_virtual_res ", blank_size_vector );
	NEW_OPT ( suite_cluster_radius , " individual_suite_cluster_radius ", 999.99 ); //IMPORTANT DO NOT CHANGE
	NEW_OPT ( loop_cluster_radius , " loop_cluster_radius ", 999.99 ); //IMPORTANT DO NOT CHANGE
	NEW_OPT ( alignment_res , " align_res_list ", blank_string_vector );
	NEW_OPT ( native_alignment_res , " native_alignment_res ", blank_size_vector );
	NEW_OPT ( jump_point_pairs , " jump_point_pairs ", blank_string_vector );
	NEW_OPT ( sample_res, "residues to build, the first element is the actual sample res while the other are the bulge residues", blank_size_vector );
	NEW_OPT ( prepend , "prepend ", true );
	// Note that this could be specified in a PDB INFO file -- but how about for silent files?
	NEW_OPT ( cluster_type, "cluster_type", "all_atom" );
	NEW_OPT ( input_res, "Residues already present in starting pose_1", blank_size_vector );
	NEW_OPT ( input_res2, "Residues already present in starting  pose_2", blank_size_vector );
	NEW_OPT ( missing_res, "Residues missing in starting pose_1, alternative to input_res", blank_size_vector );
	NEW_OPT ( missing_res2, "Residues missing in starting pose_2, alternative to input_res2", blank_size_vector );
	NEW_OPT ( cutpoint_open, "optional: chainbreak in full sequence", blank_size_vector );
	NEW_OPT ( cutpoint_closed, "optional: cutpoint at which to apply chain closure", 0 );
	NEW_OPT ( fixed_res, "optional: residues to be held fixed in minimizer", blank_size_vector );
	NEW_OPT ( minimize_res, "optional: residues to be minimize in minimizer, alternative to fixed_res", blank_size_vector );
	NEW_OPT ( virtual_res, "optional: residues to be made virtual", blank_size_vector );
	NEW_OPT ( num_pose_minimize, "optional: set_num_pose_minimize by Minimizer", 99999 );
	NEW_OPT ( sampler_num_pose_kept, "optional: set_num_pose_kept by ResidueSampler", 108 );
	NEW_OPT ( terminal_res, "optional: residues that are not allowed to stack during sampling", blank_size_vector );
	NEW_OPT ( rmsd_res, "optional: residues that will be use to calculate rmsd", blank_size_vector );
	NEW_OPT ( bulge_res, "optional: residues to be turned into a bulge variant", blank_size_vector );
	NEW_OPT ( centroid_screen, "centroid_screen", true );
	NEW_OPT ( VDW_atr_rep_screen, "classic VDW_atr_rep_screen", true );
	NEW_OPT ( sampler_perform_o2prime_pack, "perform O2' hydrogen packing inside StepWiseRNA_ResidueSampler", true );
	NEW_OPT ( fast, "quick runthrough for debugging", false );
	NEW_OPT ( medium_fast, "quick runthrough for debugging (keep more poses and not as fast as fast option)", false );
	NEW_OPT ( VERBOSE, "VERBOSE", false );
	NEW_OPT ( sampler_native_rmsd_screen, "native_rmsd_screen ResidueSampler", false );
	NEW_OPT ( sampler_native_screen_rmsd_cutoff, "sampler_native_screen_rmsd_cutoff", 2.0 );
	NEW_OPT ( native_edensity_score_cutoff, "native_edensity_score_cutoff", -1 );
	NEW_OPT ( auto_tune, "autotune rmsd for clustering between 0.1A up to 2.0A", false );
	NEW_OPT ( skip_minimize, "no minimize step in rna_swa residue sampling", false );
	NEW_OPT ( skip_sampling, "no sampling step in rna_swa residue sampling", false );
	NEW_OPT ( skip_clustering, "keep every pose, no clustering", false );
	NEW_OPT ( score_diff_min, "minimum score_diff before max_decoy_ condition applies", 0.0 ); //Oct 3, 2010
	NEW_OPT ( score_diff_cut, "score difference cut for clustering", 1000000.0 );
	NEW_OPT ( score_diff_cut_tier_two, "score_tier_two difference cut for clustering", 0.0 ); //Sept 24, 2010
	NEW_OPT ( score_diff_cut_tier_three, "score_tier_three difference cut for clustering", 0.0 ); //Sept 24, 2010
	NEW_OPT ( algorithm, "Specify algorithm to execute", "" );
	NEW_OPT ( FARFAR_start_pdb, "if pass in name, then create pdb with input chunks and the less with idealize bond length and angle residue", "" ); //Oct 16, 2010
	NEW_OPT ( allow_chain_boundary_jump_partner_right_at_fixed_BP, "allow_chain_boundary_jump_partner_right_at_fixed_BP, mainly just to get SQUARE RNA working", false ); //Nov 6, 2010
	NEW_OPT ( allow_fixed_res_at_moving_res, "allow_fixed_res_at_moving_res, mainly just to get Hermann Duplex working", false ); //Nov 15, 2010
	NEW_OPT( sampler_cluster_rmsd, " Clustering rmsd of conformations in the sampler", 0.5); //DO NOT CHANGE THIS!
	NEW_OPT( output_pdb, "output_pdb: If true, then will dump the pose into a PDB file at different stages of the stepwise assembly process.", false); //Sept 24, 2011
	////////////////////////////////////////////////////////////////////////////
	// setup
	////////////////////////////////////////////////////////////////////////////
	core::init::init ( argc, argv );
	////////////////////////////////////////////////////////////////////////////
	// end of setup
	////////////////////////////////////////////////////////////////////////////
	protocols::viewer::viewer_main ( my_main );
    } catch ( utility::excn::EXCN_Base const & e ) {
                              std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
                                  }
        return 0;
}


