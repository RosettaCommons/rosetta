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

/// @file swa_rna_main.cc
/// @author Parin Sripakdeevong (sripakpa@stanford.edu), Rhiju Das (rhiju@stanford.edu)

// libRosetta headers
#include <core/types.hh>

#include <core/chemical/AA.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueSelector.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/ResidueMatcher.hh>
#include <core/conformation/Conformation.hh>

#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rna/RNA_FittedTorsionInfo.hh>
#include <core/scoring/rna/RNA_Util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/constraints/CharmmPeriodicFunc.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>

#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.hh>
#include <core/init/init.hh>


#include <utility/file/file_sys_util.hh> //Add by Parin on May 04, 2011.

//////////////////////////////////////////////////
#include <basic/options/keys/rna.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh> // for option[ out::file::silent  ] and etc.
#include <basic/options/keys/in.OptionKeys.gen.hh> // for option[ in::file::tags ] and etc.
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>
///////////////////////////////////////////////////
#include <protocols/idealize/idealize.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <utility/excn/Exceptions.hh>
#include <numeric/xyzVector.hh>
#include <numeric/conversions.hh>

#include <core/pose/MiniPose.hh>
#include <core/pose/util.hh>

#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>
#include <core/import_pose/import_pose.hh>

//////////////////////////////////////////////////////////
#include <protocols/viewer/viewers.hh>
#include <protocols/rna/RNA_ProtocolUtil.hh>
#include <protocols/rna/RNA_LoopCloser.hh>
#include <protocols/swa/rna/StepWiseRNA_OutputData.hh>
#include <protocols/swa/rna/StepWiseRNA_CombineLongLoopFilterer.hh>
#include <protocols/swa/rna/StepWiseRNA_CombineLongLoopFilterer.fwd.hh>
#include <protocols/swa/rna/StepWiseRNA_VirtualRiboseSampler.hh>
#include <protocols/swa/rna/StepWiseRNA_BaseCentroidScreener.hh>
#include <protocols/swa/rna/StepWiseRNA_BaseCentroidScreener.fwd.hh>
#include <protocols/swa/rna/StepWiseRNA_Minimizer.hh>
#include <protocols/swa/rna/StepWiseRNA_ResidueSampler.hh>
#include <protocols/swa/rna/StepWiseRNA_Modeler.hh>
#include <protocols/swa/rna/StepWiseRNA_PoseSetup.fwd.hh>
#include <protocols/swa/rna/StepWiseRNA_PoseSetup.hh>
#include <protocols/swa/rna/StepWiseRNA_JobParametersSetup.hh>
#include <protocols/swa/rna/StepWiseRNA_JobParameters.hh>
#include <protocols/swa/rna/StepWiseRNA_Clusterer.hh>
#include <protocols/swa/rna/StepWiseRNA_VDW_BinScreener.hh>
#include <protocols/swa/rna/StepWiseRNA_VDW_BinScreener.fwd.hh>

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
#include <libgen.h>
#define GetCurrentDir getcwd


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

static basic::Tracer TR( "swa_rna_main" );


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
OPT_KEY( Boolean, minimize_and_score_sugar)
OPT_KEY( Boolean, debug_epsilon_south_sugar_mode)
OPT_KEY( Boolean, rebuild_bulge_mode)
OPT_KEY( Boolean, sampler_include_torsion_value_in_tag)
OPT_KEY( Boolean, sampler_extra_anti_chi_rotamer)
OPT_KEY( Boolean, sampler_extra_syn_chi_rotamer)
OPT_KEY( Boolean, sampler_extra_beta_rotamer)
OPT_KEY( Boolean, sampler_extra_epsilon_rotamer)
OPT_KEY( Boolean, sample_both_sugar_base_rotamer)
OPT_KEY( Boolean, reinitialize_CCD_torsions)
OPT_KEY( Boolean, PBP_clustering_at_chain_closure)
OPT_KEY( Boolean, clusterer_two_stage_clustering)
OPT_KEY( Boolean, clusterer_keep_pose_in_memory)
OPT_KEY( Boolean, finer_sampling_at_chain_closure)
OPT_KEY( StringVector, 	VDW_rep_screen_info)
OPT_KEY( Real, 	VDW_rep_alignment_RMSD_CUTOFF)
OPT_KEY( Boolean, graphic )
OPT_KEY( Real, Real_parameter_one )
OPT_KEY( Boolean, clusterer_quick_alignment )
OPT_KEY( Boolean, clusterer_align_only_over_base_atoms )
OPT_KEY( Boolean, clusterer_optimize_memory_usage )
OPT_KEY( Integer, clusterer_min_struct )
OPT_KEY( Boolean, clusterer_write_score_only )
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
OPT_KEY( Real, native_edensity_score_cutoff )
OPT_KEY( Real, score_diff_cut )
OPT_KEY( Boolean, clusterer_perform_score_diff_cut )
OPT_KEY( String, 	algorithm)
OPT_KEY( Integer, sampler_num_pose_kept)
OPT_KEY( Integer, clusterer_num_pose_kept)
OPT_KEY( Boolean, recreate_silent_struct )
OPT_KEY( Boolean, clusterer_use_best_neighboring_shift_RMSD )
OPT_KEY( Boolean, allow_chain_boundary_jump_partner_right_at_fixed_BP )
OPT_KEY( Boolean, allow_fixed_res_at_moving_res )
OPT_KEY( Boolean, clusterer_rename_tags )
OPT_KEY( Boolean, simple_append_map )
OPT_KEY( Boolean, minimizer_perform_o2star_pack )
OPT_KEY( Boolean, minimizer_output_before_o2star_pack )
OPT_KEY( Boolean, minimizer_rename_tag )
OPT_KEY( Boolean, minimizer_perform_minimize )
OPT_KEY( StringVector, 	VDW_rep_delete_matching_res)
OPT_KEY( IntegerVector, global_sample_res_list )
OPT_KEY( Boolean, clusterer_perform_VDW_rep_screen )
OPT_KEY( Boolean, clusterer_perform_filters )
OPT_KEY( Integer, clusterer_min_num_south_ribose_filter )
OPT_KEY( Real,  VDW_rep_screen_physical_pose_clash_dist_cutoff )
OPT_KEY( Boolean,  clusterer_full_length_loop_rmsd_clustering )
OPT_KEY( IntegerVector, force_syn_chi_res_list )
OPT_KEY( IntegerVector, force_north_ribose_list )
OPT_KEY( IntegerVector, force_south_ribose_list )
OPT_KEY( IntegerVector, protonated_H1_adenosine_list )
OPT_KEY( StringVector, 	sample_virtual_ribose_list)
OPT_KEY( Boolean,  sampler_assert_no_virt_ribose_sampling )
OPT_KEY( Boolean,  clusterer_ignore_FARFAR_no_auto_bulge_tag )
OPT_KEY( Boolean,  clusterer_ignore_FARFAR_no_auto_bulge_parent_tag )
OPT_KEY( Boolean,  clusterer_ignore_unmatched_virtual_res )
OPT_KEY( Boolean,  output_pdb )
OPT_KEY( String, 	start_silent)
OPT_KEY( String, 	start_tag)
OPT_KEY( Boolean,  simple_full_length_job_params )
OPT_KEY( Real, sampler_cluster_rmsd )
OPT_KEY( Boolean, 	output_extra_RMSDs)
OPT_KEY( Boolean, 	integration_test)
OPT_KEY( Boolean, 	add_virt_root ) //For Fang's electron density code.
OPT_KEY ( Boolean, constraint_chi )
OPT_KEY ( Boolean, rm_virt_phosphate )
OPT_KEY ( Boolean, choose_random )

//////////////////////////////////////////////////////////////////////////////////////
//Apply chi angle constraint to the purines --
// WHY IS THIS CODE COPIED IN SWA_RNA_ANALYTICAL_CLOSURE?
void apply_chi_cst(core::pose::Pose & pose, core::pose::Pose const & ref_pose) {
	using namespace core::conformation;
	using namespace core::id;
	using namespace core::scoring;
	using namespace core::scoring::rna;
	using namespace core::scoring::constraints;
	using namespace core::chemical;

	Size const nres = pose.total_residue();
	ConstraintSetOP cst_set = new ConstraintSet;
	for (Size i = 1; i <= nres; ++i) {
		Residue const & res = pose.residue(i);
		if ( res.is_RNA() && (res.aa() == na_rad || res.aa() == na_rgu)) {
			Real const chi = numeric::conversions::radians( ref_pose.torsion( TorsionID( i, id::CHI, 1 ) ) );
			FuncOP chi_cst_func ( new CharmmPeriodicFunc( chi, 1.0, 1.0 ) );
			AtomID const atom1 (res.atom_index("C2'"), i);
			AtomID const atom2 (res.atom_index("C1'"), i);
			AtomID const atom3 ( is_purine(res) ? res.atom_index("N9") : res.atom_index("N1"), i);
			AtomID const atom4 (is_purine(res) ? res.atom_index("C4") : res.atom_index("C2"), i);
			cst_set->add_constraint( new DihedralConstraint( atom1, atom2, atom3, atom4, chi_cst_func ) );
		}
	}
	pose.constraint_set ( cst_set );
}


//////////////////////////////////////////////////////////////////////////////////////
std::string
get_working_directory(){

	char cCurrentPath[FILENAME_MAX];
	std::string current_directory_string;

 	if (!GetCurrentDir(cCurrentPath, sizeof(cCurrentPath))){
 	 utility_exit_with_message( "!GetCurrentDir(cCurrentPath, sizeof(cCurrentPath))" );
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
bool
Is_nonempty_input_silent_file(std::string const input_silent_file, std::string const exit_key_string){

    std::cout << "Checking that input_silent_file " << input_silent_file << " contain actual silent_structs or the correct exit_key_string" << std::endl;

		std::ifstream infile;
	 	infile.open(input_silent_file.c_str());

		if (infile.fail()){
		 utility_exit_with_message("Error! \"" + input_silent_file + "\" could not be opened!");
		}else{
			std::cout << "Open \"" << input_silent_file << "\" successful!" << std::endl;
		}

		std::string line;

		bool found_line=getline(infile, line);

		if(found_line==false) utility_exit_with_message("No line exist in input_silent_file= " + input_silent_file);

		size_t found_substring=line.find(exit_key_string);

		if(found_substring!=std::string::npos){
			std::cout << "input_silent_file: " << input_silent_file << " contain no silent struct" << std::endl;
			std::cout << line << std::endl;

			//consistency check:////////////////////////////////////////////////////////////////////////////////////////////////////
			std::string next_line;
			bool found_next_line=getline(infile, next_line);
			if(found_next_line) std::cout << "input silent_file contain more than one line! next_line= " << next_line << std::endl;
			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			return false;
		}else{
			return true;
		}
}
//////////////////////////////////////////////////////////////////////////////////////

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

//////////////////////////////////////////////////////////////////////////////////////


utility::vector1< std::string >
get_silent_file_tags(){

	using namespace protocols::swa::rna;

	bool tags_from_command_line=false;
	bool tags_from_filterer_outfile=false;

//	bool const verbose=option[ VERBOSE ]();

	utility::vector1< std::string > input_silent_file_tags;

	if( option[ in::file::tags ].user()){
		tags_from_command_line=true;

		input_silent_file_tags=option[ in::file::tags ]();

	}


	if(option[ job_queue_ID ].user() && option[ filter_output_filename ].user()){

		Output_title_text("importing tag from filter_outfile", TR );

		bool combine_long_loop=option[ combine_long_loop_mode ]();

		bool combine_helical=option[ combine_helical_silent_file ]();

		if(combine_long_loop==false && combine_helical==false) utility_exit_with_message("combine_long_loop ==false && combine_helical==false");

		tags_from_filterer_outfile=true;

		std::string const filtered_tag_file=option[ filter_output_filename ]();

		std::ifstream infile;
	 	infile.open(filtered_tag_file.c_str());

		if (infile.fail()){
		 utility_exit_with_message("Error! \"" + filtered_tag_file + "\" could not be opened!");
		}else{
			std::cout << "Open \"" << filtered_tag_file << "\" successful!" << std::endl;
		}



		//Becareful here... job_queue_ID start from ZERO!
		int const queue_ID=option[ job_queue_ID ]();
		int ID=0;

		std::cout << "queue_ID= " << queue_ID << std::endl;

		std::string tag_pair_string;

		bool found_queue_ID=false;

		while(getline(infile, tag_pair_string) ){

			if(queue_ID==ID){
				found_queue_ID=true;
				break;
			}

			ID++;
		}

		//Warning queue_ID start at ZERO!
		if(found_queue_ID==false) utility_exit_with_message("found_queue_ID==false, queue_ID= " + string_of(queue_ID) + " num_tag_string_in_file= " + string_of(ID) );


		std::cout << "import silent_file_tags: " << tag_pair_string << " from filter_output_filename= " << filtered_tag_file << std::endl;

		infile.close();

		utility::vector1< std::string > const line_list=Tokenize(tag_pair_string," \t\n\f\v"); //Oct 19, 2010..now filterer_outfile contain other terms.

		input_silent_file_tags.clear();
		input_silent_file_tags.push_back(line_list[1]);
		input_silent_file_tags.push_back(line_list[2]);

		Output_title_text("", TR );

	}


	if( (tags_from_command_line==false) && (tags_from_filterer_outfile==false) ){
		utility_exit_with_message("(tags_from_command_line==false) && (tags_from_filterer_outfile==false)");
	}

	if( (tags_from_command_line==true) && (tags_from_filterer_outfile==true) ){
		utility_exit_with_message("(tags_from_command_line==true) && (tags_from_filterer_outfile==true)");
	}

	return input_silent_file_tags;

}


//////////////////////////////////////////////////////////////////////////////////////

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

	core::scoring::ScoreFunctionOP scorefxn = getScoreFunction();


	if(option[minimize_and_score_sugar]()==false){
		std::cout << "WARNING minimize_and_score_sugar is false, SET rna_sugar_close weight to 0.0 " << std::endl;
    scorefxn->set_weight( rna_sugar_close, 0.000000000000 );

		//Sept 16, 2010. Thought about include a very small weight for rna_sugar_close so that column # will not change. HOWEVER this significant change the minimization results!
		//scorefxn->set_weight( rna_sugar_close, 0.000000000001 );
 	}

	std::cout << "---------score function weights----------" << std::endl;
	scorefxn->show(std::cout);
	std::cout << "-----------------------------------------" << std::endl;


	return scorefxn;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///Oct 28, 2011:This function is an alternative way to setup the job_parameters////////////////////////////////////////////////////
///The original job_parameters was designed primarily for use during the SWA building step/////////////////////////////////////////
///However since then, many other other SWA and FARFAR class make use of it especially the full pose length version////////////////
///Hence the original job_parameters setup has become clumbersome for most usage (containing many options and independencies)//////
///The new job_parameter setup function is a work in progress and will be updated as needed.///////////////////////////////////////
///For example right now (Oct 28, 2011), it only contains parameters needed by Output_data()///////////////////////////////////////

protocols::swa::rna::StepWiseRNA_JobParametersOP
setup_simple_full_length_rna_job_parameters(){

	using namespace protocols::swa::rna;
	using namespace ObjexxFCL;
  using namespace core::pose;
  using namespace core::chemical;

	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( RNA );

	/////////////////////////////////////////////////////
	if ( !option[ in::file::fasta].user() ) utility_exit_with_message( "Must supply in::file::fasta!" );
	if ( !option[ rmsd_res ].user() ) utility_exit_with_message( "Must supply rmsd_res!" );
	if ( !option[ alignment_res ].user() ) utility_exit_with_message( "Must supply alignment_res!" );
	if ( !option[ global_sample_res_list ].user() ) utility_exit_with_message( "Must supply global_sample_res_list!" );

	/////////////Read in sequence.///////////////////////

	std::string const fasta_file = option[ in::file::fasta ]()[1];
	core::sequence::SequenceOP fasta_sequence = core::sequence::read_fasta_file( fasta_file )[1];
	std::string const full_sequence = fasta_sequence->sequence();
	core::Size const nres=full_sequence.length();

	/////////////////////////////////////////////////////

	StepWiseRNA_JobParametersOP job_parameters = new StepWiseRNA_JobParameters;

	job_parameters->set_Is_simple_full_length_job_params(true); //FORGOT THIS! ONLY INCLUDED ON MARCH 03, 2012!
	//LUCKILY, BEFORE MARCH 03, 2012, only called Is_simple_full_length_job_params() for utility_exit_with_message() check in the following three functions of StepWiseRNA_Clusterer.cc:
	//i. initialize_VDW_rep_screener(),
	//ii. recalculate_rmsd_and_output_silent_file()
	//iii. get_best_neighboring_shift_RMSD_and_output_silent_file()
	/////////////////////////////////////////////////////

	job_parameters->set_output_extra_RMSDs( option[ output_extra_RMSDs ]() );

	utility::vector1< core::Size > is_working_res(nres, 1); //All res belong to 'mock pose' 1
	job_parameters->set_is_working_res( is_working_res );
	job_parameters->set_full_sequence( full_sequence ); //working_sequence is automatical init after BOTH is_working_res and full_sequence is initialized


	job_parameters->set_rmsd_res_list( option[ rmsd_res ]() );
 	job_parameters->set_gap_size( 0 );

	if ( option[ native_alignment_res ].user() ){
		//If working_native_align_res is not specified that most function will internally default to using working_best_alignment_list in its place
		job_parameters->set_native_alignment( option[ native_alignment_res ]() );
		job_parameters->set_working_native_alignment( option[ native_alignment_res ]() );
	}
	/////////////////////////////////////////////////////
	//Oct 31, 2011: Better to let code raise error if the there is a statement asking for the partition_definition from the job_params later in the run!
	//ObjexxFCL::FArray1D_bool partition_definition( nres, false ); //All res belong in partition 0!
	//job_parameters->set_partition_definition( partition_definition ); //this is a useful decomposition.
	/////////////////////////////////////////////////////

	std::map< core::Size, core::Size > full_to_sub;
	std::map< core::Size, bool > Is_prepend_map;

	utility::vector1< Size > input_res;
	utility::vector1< Size > input_res2;

	for( Size seq_num = 1; seq_num <= full_sequence.size(); seq_num++ ){
		input_res.push_back(seq_num);
		full_to_sub[ seq_num ] = seq_num;
		Is_prepend_map[seq_num]=false; //all append..arbitrary choice
	}

	job_parameters->set_full_to_sub( full_to_sub );  //res_map
	job_parameters->set_Is_prepend_map( Is_prepend_map );
	job_parameters->set_Is_prepend( Is_prepend_map[job_parameters->actually_moving_res()] );

	utility::vector1< utility::vector1< Size > > input_res_vectors;
	input_res_vectors.push_back( input_res );
	input_res_vectors.push_back( input_res2 );
	job_parameters->set_input_res_vectors(input_res_vectors);

	/////////////////////////////////////////////////////
	utility::vector1< Size > working_moving_res_list;
	working_moving_res_list.push_back( full_sequence.size() ); //arbitrary choice, choose residue of pose

	job_parameters->set_working_moving_res_list( working_moving_res_list); //This sets actually_moving_res()
	//////////////////////////////////////////////////////////////////////////////////////////////////////////

	utility::vector1< std::string > alignment_res_string_list=option[ alignment_res ]();
	utility::vector1< core::Size > best_alignment_list;

	for(Size n=1; n<=alignment_res_string_list.size(); n++){

		utility::vector1< std::string > alignments_res_string=Tokenize(alignment_res_string_list[n], "-");
		utility::vector1< core::Size >  alignment_list;

		for(Size ii=1; ii<=alignments_res_string.size(); ii++){
			alignment_list.push_back( string_to_int( alignments_res_string[ii] ) );
		}

		if( alignment_list.size()>best_alignment_list.size() ) best_alignment_list=alignment_list;
	}

	job_parameters->set_working_best_alignment(best_alignment_list);

	//////////////////////////////////////////////////////////////////////////////////////////////////////////


	job_parameters->set_working_fixed_res(  get_fixed_res(nres) );
	job_parameters->set_global_sample_res_list(  option[ global_sample_res_list ]() );
	job_parameters->set_force_syn_chi_res_list(  option[ force_syn_chi_res_list]() );
	job_parameters->set_force_north_ribose_list(  option[ force_north_ribose_list ]() );
	job_parameters->set_force_south_ribose_list(  option[ force_south_ribose_list ]() );
	job_parameters->set_protonated_H1_adenosine_list( option[ protonated_H1_adenosine_list ]() );


	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	PoseOP native_pose_OP;
	if (option[ in::file::native ].user() ){
		native_pose_OP = PoseOP( new Pose );
		import_pose::pose_from_pdb( *native_pose_OP, *rsd_set, option[ in::file::native ]() );
		protocols::rna::make_phosphate_nomenclature_matches_mini(*native_pose_OP);

		utility::vector1< core::Size > const native_virtual_res_list = option[ native_virtual_res]();

		for(Size n=1; n<=native_virtual_res_list.size(); n++){
			apply_virtual_rna_residue_variant_type( (*native_pose_OP), native_virtual_res_list[n] , false /*apply_check*/);
		}

		job_parameters->set_working_native_pose( native_pose_OP );

	}
	//////////////////////////////////////////////////////////////////////////////////////////////////////////

	print_JobParameters_info(job_parameters, "simple_full_length_job_parameters", true /*Is_simple_full_length_JP*/);

	return job_parameters;

}


//////////////////////////////////////////////////////////////////////////////////////
protocols::swa::rna::StepWiseRNA_JobParametersOP
setup_rna_job_parameters(bool check_for_previously_closed_cutpoint_with_input_pose=false){
	// Describe use_case of flag check_for_previously_closed_cutpoint_with_input_pose? -- rhiju
	//
	//   i. Consider tetraloop-receptor motif:
	//
	//
	//        Tetraloop
	//            |      G14-C15
	//            |      G13-C16
	//           A4 A5   U12 U17  <--Receptor
	//          G3   A6  U11 U18
	//           C2-G7   C10-C19
	//           C1-G8   C09-C20
	//
	//
	//    Assumes a fold-tree with jump-points C1-C8, G8-C09, C09-C20 and G14-C15
	//
	//   ii. Suppose tetraloop was built first:
	//
	//           A4*A5
	//          G3   A6
	//           C2-G7   C10-C19
	//           C1-G8   C09-C20
	//
	//       Becuase C1-C8 is a jump-point, there must be a corresponding closed-
	//       cutpoint somwhere in the loop (e.g. between A4 and A5; see *). This
	//       closed-cutpoint position can also be different for the various poses
	//       in the silent_file (assuming that multiple built-up paths were allowed).
	//
	//
	//   iii. Now we start building the receptor:
	//
	//
	//           A4*A5
	//          G3   A6  U11
	//           C2-G7   C10-C19
	//           C1-G8   C09-C20
	//
	//       Rhiju: Why is this so important that it is required as input for
	//              setup_rna_job_parameters?
	//
	//       The setup_rna_job_parameters class is responsible for figuring out
	//       the fold_tree. For simple motifs, such as a single hairpin or a
	//       single double-stranded motif, setup_rna_job_parameters can figure
	//       what the fold_tree is using just information specified from command-line
	//       (jump_point_pairs, cutpoint_open, cutpoint_closed and etcs).
	//
	//       However, for more complex motif like tetraloop-receptor, it cannot
	//       figure out the previously_closed_cupoint position of the tetraloop (i.e
	//       where the * position is located) by using just the command line
	//       information.
	//
	//       Setting check_for_previously_closed_cutpoint_with_input_pose to true
	//       tells the class to look inside the input_pose to figure out where this
	//       closed_cutpoint is located (see get_previously_closed_cutpoint_from_
	//       imported_silent_file() function insetup_rna_job_parameters class for
	//       futher details).
	//
	//       Lastly, note that the setup_rna_job_parameters class will still be able
	//       to come up with a fold_tree even when
	//       check_for_previously_closed_cutpoint_with_input_pose is false (default)
	//       In this case, however, the fold_tree might not have the correct
	//       closed-cutpoint position.
	//
	//       This is actually the default behavior since in most cases
	//       (e.g. clustering, filtering), the specific details of the fold_tree
	//       is not used. Additionally clustering and filtering deals with multiple
	//       poses and not all of them have the same previously_closed_cutpoint
	//       anyways.
	//
	// For "check_for_previously_closed_cutpoint_with_input_pose", can we just make
	//    it true by default? -- rhiju
	//
	//        I would say no.
	//
	//        If "check_for_previously_closed_cutpoint_with_input_pose" is true, then
	//        user will be required to pass in --in:file:silent and -tags options
	//        from command line. These flags are appropriate for the sampler but
	//        might not be appropriate for filterer and clusterer where -tags option is
	//        not required. Essentially the clusterer and filterer deals with multiple
	//        poses and the specific details each pose's fold_tree (e.g. the position
	//        of the previously_closed_cutpoint) is not used.
	//
	//        To simplify the code, I would suggest removing the
	//        "check_for_previously_closed_cutpoint_with_input_pose" option all together.
	//        Perhaps, the benefit of this option (as outlined in the last email) isn't
	//        worth the extra complexity it brings to the code.
	//             -- parin (2013)

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
	stepwise_rna_job_parameters_setup.set_fixed_res( get_fixed_res(nres) );
	stepwise_rna_job_parameters_setup.set_terminal_res( option[ terminal_res ]() );
	stepwise_rna_job_parameters_setup.set_rmsd_res_list( option[ rmsd_res ]() );

	// jump_point_pairs is string of pairs  "1-16 8-9", assumed for now to be connected by dashes.  See note in StepWiseRNA_JobParametersSetup.cc
	stepwise_rna_job_parameters_setup.set_jump_point_pair_list( option[ jump_point_pairs ]() ); //Important!: Need to be called after set_fixed_res

	//  Alignment_res is a string vector to allow user to specify multiple possible alignments. Note that 1-6 7-12 is different from 1-6-7-12. See note in StepWiseRNA_JobParametersSetup.cc
	stepwise_rna_job_parameters_setup.set_alignment_res( option[ alignment_res ]() );
	stepwise_rna_job_parameters_setup.set_filter_user_alignment_res( option[ filter_user_alignment_res ]() );
	stepwise_rna_job_parameters_setup.set_native_alignment_res( option[ native_alignment_res ]() );

	stepwise_rna_job_parameters_setup.set_global_sample_res_list( option[ global_sample_res_list ]() ); //March 20, 2011

	stepwise_rna_job_parameters_setup.set_force_syn_chi_res_list( option[ force_syn_chi_res_list]() ); //April 29, 2011
	stepwise_rna_job_parameters_setup.set_force_north_ribose_list( option[ force_north_ribose_list ]() ); //April 29, 2011
	stepwise_rna_job_parameters_setup.set_force_south_ribose_list( option[ force_south_ribose_list ]() ); //April 29, 2011
	stepwise_rna_job_parameters_setup.set_protonated_H1_adenosine_list( option[ protonated_H1_adenosine_list ]() ); //May 02, 2011

	stepwise_rna_job_parameters_setup.set_allow_chain_boundary_jump_partner_right_at_fixed_BP( option[ allow_chain_boundary_jump_partner_right_at_fixed_BP ]() ); //Hacky just to get Square RNA working.

	stepwise_rna_job_parameters_setup.set_output_extra_RMSDs( option[ output_extra_RMSDs ]() );
	stepwise_rna_job_parameters_setup.set_add_virt_res_as_root( option[ add_virt_root]() );


	/////////////////////////////Sept 1, 2010////////////
	if(check_for_previously_closed_cutpoint_with_input_pose){
		utility::vector1< std::string > input_tags;
		utility::vector1< std::string > silent_files_in;

		// First read in any information on pdb read in from silent files.
		// Assume one to one correspondence between number of tags and number of silent_file
		if ( option[ in::file::silent ].user() ) {
			silent_files_in = option[ in::file::silent ]();
			input_tags = get_silent_file_tags();

			if(silent_files_in.size()!=input_tags.size()){
			 	utility_exit_with_message("silent_files_in.size(" + string_of(silent_files_in.size()) + ")!=input_tags.size(" + string_of(input_tags.size())+ ")");
			}
		}
		stepwise_rna_job_parameters_setup.set_input_tags( input_tags);
		stepwise_rna_job_parameters_setup.set_silent_files_in( silent_files_in);
	}
	///////////////////////////////////////////////////////

	stepwise_rna_job_parameters_setup.apply();

	return stepwise_rna_job_parameters_setup.job_parameters();

}



void
setup_copy_DOF_input(protocols::swa::rna::StepWiseRNA_PoseSetupOP & stepwise_rna_pose_setup){

	/////////////////////////////////////////////////////////////////////////////////////////
	// StepWisePoseSetup should create the starting pose.
	// This class might eventually be united with the protein StepWisePoseSetup.
	utility::vector1< std::string > input_tags;
	utility::vector1< std::string > silent_files_in;

	// First read in any information on pdb read in from silent files.
	// Assume one to one correspondence between number of tags and number of silent_file
	if ( option[ in::file::silent ].user() ) {
		silent_files_in = option[ in::file::silent ]();
		input_tags = get_silent_file_tags();

		if(silent_files_in.size()!=input_tags.size()){
			utility_exit_with_message("silent_files_in.size()!=input_tags.size()");
		}
	}

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

//////////////////////////////////////////////////////////////////////////////////////////

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


//	StepWiseRNA_PoseSetup stepwise_rna_pose_setup( pdb_tags, silent_files_in, job_parameters);

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

	return stepwise_rna_pose_setup;
}


void
filter_combine_long_loop()
{

  using namespace core::pose;
  using namespace core::chemical;
  using namespace core::kinematics;
  using namespace core::scoring;
	using namespace protocols::swa::rna;

	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( RNA );

	///////////////////////////////
	StepWiseRNA_JobParametersOP	job_parameters = setup_rna_job_parameters(false);
	StepWiseRNA_JobParametersCOP job_parameters_COP( job_parameters );

	if ( option[ in::file::silent ].user() ==false) utility_exit_with_message(" option[ in::file::silent ].user() ==false");

	utility::vector1< std::string > const silent_files_in=option[ in::file::silent ]();
	std::string const output_filename = option[ filter_output_filename ]();

	if(silent_files_in.size()!=2) utility_exit_with_message("silent_files_in.size()!=2");

	bool at_least_one_empty_silent_file=false;
	for(Size n=1; n<=2; n++){
		std::string const input_silent_file=silent_files_in[n];

		if(Is_nonempty_input_silent_file(input_silent_file, "empty cluster silent_file since all input_silent_file are empty.")==false){
			std::cout << input_silent_file << " is empty" << std::endl;
			at_least_one_empty_silent_file=true;
		}

	}

	if ( utility::file::file_exists( output_filename ) ) { //Feb 08, 2012
		std::cout << "WARNING: output_filename " << output_filename << " already exist! removing..." << std::endl;

		int remove_file_return_value=std::remove( output_filename.c_str() );
		std::cout << "remove_file_return_value= " <<  remove_file_return_value << " for std::remove(" << output_filename << ")" << std::endl;

		if(remove_file_return_value!=0){
			utility_exit_with_message("remove_file_return_value=" +ObjexxFCL::string_of(remove_file_return_value )+ "!=0 for std::remove(" + output_filename + ")" );
		}
	}


	if(at_least_one_empty_silent_file){
		std::cout << "Early Exit: since at_least_one_empty_silent_file, outputting empty filterer outfile " << std::endl;

		std::ofstream outfile;
		outfile.open(output_filename.c_str()); //Opening the file with this command removes all prior content..


		outfile << "empty cluster silent_file since at least one of the two input_silent_file is empty.";
		outfile << " input_silent_file:";
		for(Size n=1; n<=silent_files_in.size(); n++){
			outfile << " " << silent_files_in[n];
		}
		outfile << "\n";

		outfile.flush();
		outfile.close();
		return; //Early return;
	}


	StepWiseRNA_CombineLongLoopFiltererOP stepwise_combine_long_loop_filterer = new StepWiseRNA_CombineLongLoopFilterer( job_parameters_COP, option[combine_helical_silent_file]);

	stepwise_combine_long_loop_filterer->set_max_decoys( option[clusterer_num_pose_kept]() ); //Updated on Jan 12, 2012

	stepwise_combine_long_loop_filterer->set_silent_files_in( silent_files_in );
	stepwise_combine_long_loop_filterer->set_output_filename( output_filename );

	//Remove score filtering on Jan 12, 2012

	stepwise_combine_long_loop_filterer->set_filter_for_previous_contact( option[filter_for_previous_contact] );
	stepwise_combine_long_loop_filterer->set_filter_for_previous_clash( option[filter_for_previous_clash] );
	stepwise_combine_long_loop_filterer->set_undercount_ribose_rotamers( option[filterer_undercount_ribose_rotamers] );


	stepwise_combine_long_loop_filterer->set_parin_favorite_output( option[ parin_favorite_output ]());

	stepwise_combine_long_loop_filterer->filter();

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
rna_sample_virtual_ribose(){ //July 19th, 2011...rebuild the bulge nucleotides after floating base step to properly close the chain.

  using namespace core::pose;
  using namespace core::chemical;
  using namespace core::kinematics;
  using namespace core::scoring;
	using namespace protocols::swa::rna;

	Output_title_text("Enter rna_sample_virtual_ribose()", TR );

	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( RNA );

	core::scoring::ScoreFunctionOP const scorefxn=create_scorefxn();

	StepWiseRNA_JobParametersOP	job_parameters = setup_rna_job_parameters(false);
	StepWiseRNA_JobParametersCOP job_parameters_COP( job_parameters );

	StepWiseRNA_PoseSetupOP stepwise_rna_pose_setup = setup_pose_setup_class(job_parameters, false /*COPY DOF*/);

	utility::vector1< std::string > const sample_virtual_ribose_string_list= option[ sample_virtual_ribose_list ]();

	////////////////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< std::string > input_tags;
	utility::vector1< std::string > silent_files_in;

	if ( option[ in::file::silent ].user()==false ) utility_exit_with_message( "option[ in::file::silent ].user()==false");

	silent_files_in = option[ in::file::silent ]();
	input_tags = get_silent_file_tags();

	if(silent_files_in.size()!=1) utility_exit_with_message("silent_files_in.size()!=1");

	if(silent_files_in.size()!=input_tags.size()) utility_exit_with_message("silent_files_in.size()!=input_tags.size()");

  pose::Pose pose;
	import_pose_from_silent_file(pose, silent_files_in[ 1 ], input_tags[1] );
	protocols::rna::assert_phosphate_nomenclature_matches_mini(pose);

  std::string const silent_file_out = option[ out::file::silent  ]();
	////////////////////////////////////////////////////////////////////////////////////////////////////

	if(option[ graphic ]()){
		std::string const current_directory_string=get_working_directory();
		protocols::viewer::add_conformation_viewer( pose.conformation(), current_directory_string, 400, 400 );
	}

	stepwise_rna_pose_setup->setup_native_pose( pose ); //NEED pose to align native_pose to pose.

	// Hey! This should be in a class! -- rhiju, 2013.
	sample_user_specified_virtual_riboses(pose, sample_virtual_ribose_string_list, job_parameters_COP, scorefxn, silent_file_out, input_tags[1], option[ integration_test ]()  );

	Output_title_text("Exit rna_sample_virtual_ribose()", TR );

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void check_if_silent_file_exists(){

	////////////////////////////////////////////////////////////////////////////////
	//May 04, 2011. Make sure that the "FINAL" output silent_file doesn't exist at the beginning of the run.
	//Mainly need this becuase the BIOX2-cluster (Stanford) is not robust...The node containing a slave_job can suddenly become unaviable and slave_job will need to be resubmitted to a new node.
	//For example:
	//		Your job <90107> has been killed because the execution host <node-9-30> is no longer available.
	//		The job will be re-queued and re-run with the same jobId.

	std::string const silent_file = option[ out::file::silent ]();

	std::cout << "Does following file exist? " <<  silent_file << " " << utility::file::file_exists( silent_file )  << std::endl;

	if ( utility::file::file_exists( silent_file ) ) {

		std::cout << "WARNING: silent_file " << silent_file << " already exists! removing..." << std::endl;

		int remove_file_return_value=std::remove( silent_file.c_str() );
		std::cout << "remove_file_return_value= " <<  remove_file_return_value << " for std::remove(" << silent_file << ")" << std::endl;
		if( !remove_file_return_value )	utility_exit_with_message("remove_file_return_value=" +ObjexxFCL::string_of(remove_file_return_value )+ "!=0 for std::remove(" + silent_file + ")" );
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
swa_rna_sample()
{

  using namespace core::pose;
  using namespace core::chemical;
  using namespace core::kinematics;
  using namespace core::scoring;
	using namespace protocols::swa::rna;

	Output_title_text("Enter swa_rna_sample()", TR );

	StepWiseRNA_JobParametersOP	job_parameters = setup_rna_job_parameters(true  /*check_for_previously_closed_cutpoint_with_input_pose */);
	StepWiseRNA_JobParametersCOP job_parameters_COP( job_parameters );
	StepWiseRNA_PoseSetupOP stepwise_rna_pose_setup = setup_pose_setup_class(job_parameters);

  Pose pose;
	stepwise_rna_pose_setup->apply( pose );
	stepwise_rna_pose_setup->setup_native_pose( pose ); //NEED pose to align native_pose to pose.
	PoseCOP native_pose = job_parameters_COP->working_native_pose();

	if ( option[ graphic ]() ) protocols::viewer::add_conformation_viewer ( pose.conformation(), get_working_directory(), 400, 400 );

	core::scoring::ScoreFunctionOP scorefxn=create_scorefxn();
	if ( option[ constraint_chi ]() )  apply_chi_cst( pose, *job_parameters_COP->working_native_pose() );

	// following is temporarily turned off... should be redundant anyway with ensure_directory_for_out_silent_file_exists.
	//	check_if_silent_file_exists();

	// Most of the following exactly matches ERRASER_Modeler setup in swa_rna_analytical_closure. Get rid of that other file!!
	Size const working_moving_res = job_parameters_COP->working_moving_res();
	StepWiseRNA_ModelerOP stepwise_rna_modeler = new StepWiseRNA_Modeler( working_moving_res, scorefxn );

	stepwise_rna_modeler->set_job_parameters( job_parameters );

	stepwise_rna_modeler->set_native_pose( native_pose );
	stepwise_rna_modeler->set_silent_file( option[ out::file::silent  ] );
	stepwise_rna_modeler->set_sampler_num_pose_kept ( option[ sampler_num_pose_kept ]() );
	stepwise_rna_modeler->set_fast ( option[ fast ]() );
	stepwise_rna_modeler->set_medium_fast ( option[ medium_fast ]() );
	// following should probably be 'reference pose' --> does not have to be native.
	stepwise_rna_modeler->set_sampler_native_rmsd_screen ( option[ sampler_native_rmsd_screen ]() );
	stepwise_rna_modeler->set_sampler_native_screen_rmsd_cutoff ( option[ sampler_native_screen_rmsd_cutoff ]() );
	stepwise_rna_modeler->set_o2star_screen ( option[ sampler_perform_o2star_pack ]() );
	stepwise_rna_modeler->set_verbose ( option[ VERBOSE ]() );
	stepwise_rna_modeler->set_cluster_rmsd (	option[ sampler_cluster_rmsd ]()	);
	stepwise_rna_modeler->set_distinguish_pucker ( option[ distinguish_pucker]() );
	stepwise_rna_modeler->set_finer_sampling_at_chain_closure ( option[ finer_sampling_at_chain_closure]() );
	stepwise_rna_modeler->set_PBP_clustering_at_chain_closure ( option[ PBP_clustering_at_chain_closure]() );
	stepwise_rna_modeler->set_allow_syn_pyrimidine( option[ sampler_allow_syn_pyrimidine ]() );
	stepwise_rna_modeler->set_extra_syn_chi_rotamer ( option[ sampler_extra_syn_chi_rotamer]() );
	stepwise_rna_modeler->set_extra_anti_chi_rotamer ( option[ sampler_extra_anti_chi_rotamer]() );
	stepwise_rna_modeler->set_use_phenix_geo ( option[ basic::options::OptionKeys::rna::corrected_geo ]() );
	stepwise_rna_modeler->set_centroid_screen ( option[ centroid_screen ]() );
	stepwise_rna_modeler->set_VDW_atr_rep_screen ( option[ VDW_atr_rep_screen ]() );
	stepwise_rna_modeler->set_VDW_rep_screen_info ( option[ VDW_rep_screen_info ]() );
	stepwise_rna_modeler->set_VDW_rep_alignment_RMSD_CUTOFF ( option[ VDW_rep_alignment_RMSD_CUTOFF ]() );
	// not yet implemented.
	//	stepwise_rna_modeler->set_force_centroid_interaction ( option[force_centroid_interaction]() );
	stepwise_rna_modeler->set_choose_random( option[ choose_random ]()  );
	stepwise_rna_modeler->set_nstruct( option[ out::nstruct ]() );
	stepwise_rna_modeler->set_skip_sampling( option[ skip_sampling ]() );
	stepwise_rna_modeler->set_perform_minimize( option[ minimizer_perform_minimize ]() );
	stepwise_rna_modeler->set_native_edensity_score_cutoff ( option[native_edensity_score_cutoff]() );
	stepwise_rna_modeler->set_rm_virt_phosphate ( option[rm_virt_phosphate]() );
	stepwise_rna_modeler->set_minimize_and_score_sugar ( option[ minimize_and_score_sugar ]() );
	stepwise_rna_modeler->set_minimize_and_score_native_pose ( option[ minimize_and_score_native_pose ]() );
	if ( option[ num_pose_minimize ].user() ) stepwise_rna_modeler->set_num_pose_minimize( option[ num_pose_minimize ]() );
	stepwise_rna_modeler->set_output_minimized_pose_data_list( true );

	// newer options, not yet shared with ERRASER modeler. I think.
	stepwise_rna_modeler->set_VDW_rep_delete_matching_res ( option[ VDW_rep_delete_matching_res ]() );
	stepwise_rna_modeler->set_VDW_rep_screen_physical_pose_clash_dist_cutoff ( option[ VDW_rep_screen_physical_pose_clash_dist_cutoff ]() );
	stepwise_rna_modeler->set_integration_test_mode( option[ integration_test ]() ); //Should set after setting sampler_native_screen_rmsd_cutoff, fast, medium_fast options.
	stepwise_rna_modeler->set_allow_bulge_at_chainbreak( option[ allow_bulge_at_chainbreak ]() );
	stepwise_rna_modeler->set_parin_favorite_output( option[ parin_favorite_output ]());
	stepwise_rna_modeler->set_floating_base( option[ floating_base ]() );
	stepwise_rna_modeler->set_include_syn_chi( option[ include_syn_chi ]() );
	stepwise_rna_modeler->set_reinitialize_CCD_torsions( option[ reinitialize_CCD_torsions]() );
	stepwise_rna_modeler->set_sampler_extra_epsilon_rotamer( option[ sampler_extra_epsilon_rotamer]() );
	stepwise_rna_modeler->set_sampler_extra_beta_rotamer( option[ sampler_extra_beta_rotamer]() );
	stepwise_rna_modeler->set_sample_both_sugar_base_rotamer( option[ sample_both_sugar_base_rotamer]() ); //Nov 12, 2010
	stepwise_rna_modeler->set_sampler_include_torsion_value_in_tag( option[ sampler_include_torsion_value_in_tag]() );
	stepwise_rna_modeler->set_rebuild_bulge_mode( option[ rebuild_bulge_mode ]() );
	stepwise_rna_modeler->set_debug_epsilon_south_sugar_mode( option[ debug_epsilon_south_sugar_mode ]() );
	stepwise_rna_modeler->set_exclude_alpha_beta_gamma_sampling( option[ exclude_alpha_beta_gamma_sampling ]() );
	stepwise_rna_modeler->set_combine_long_loop_mode( option[ combine_long_loop_mode]() );
	stepwise_rna_modeler->set_do_not_sample_multiple_virtual_sugar( option[ do_not_sample_multiple_virtual_sugar]() );
	stepwise_rna_modeler->set_sample_ONLY_multiple_virtual_sugar( option[ sample_ONLY_multiple_virtual_sugar]() );
	stepwise_rna_modeler->set_sampler_assert_no_virt_ribose_sampling( option[ sampler_assert_no_virt_ribose_sampling ]() );
	stepwise_rna_modeler->set_allow_base_pair_only_centroid_screen( option[ allow_base_pair_only_centroid_screen ]() );

	// this is new, not in ERRASER (swa_rna_analytical_closure)
	stepwise_rna_modeler->set_minimizer_perform_o2star_pack( option[ minimizer_perform_o2star_pack ]() );
	stepwise_rna_modeler->set_minimizer_output_before_o2star_pack( option[ minimizer_output_before_o2star_pack ]() );
	stepwise_rna_modeler->set_minimizer_rename_tag( option[ minimizer_rename_tag ]() );

	// currently creates silent file -- instead we should be able to output those silent structs if we want them.
	// probably should output best scoring pose, not whatever comes out randomly.
	stepwise_rna_modeler->apply( pose );


}


///////////////////////////////////////////////////////////////
void
swa_rna_cluster(){

	using namespace protocols::swa::rna;

	StepWiseRNA_JobParametersOP job_parameters;

	bool job_parameters_exist=false;

	if(option[ simple_full_length_job_params ]()){ //Oct 31, 2011.
		std::cout << "USING simple_full_length_job_params!" << std::endl;
		job_parameters_exist=true;
		job_parameters=setup_simple_full_length_rna_job_parameters();

	}else if(option[ rmsd_res ].user()){
		job_parameters_exist=true;
		job_parameters=setup_rna_job_parameters();
		print_JobParameters_info(job_parameters, "standard_clusterer_job_params", false /*Is_simple_full_length_JP*/);
	}

	StepWiseRNA_JobParametersCOP job_parameters_COP=job_parameters;

	//////////////////////////////////////////////////////////////

	StepWiseRNA_VDW_BinScreenerOP user_input_VDW_bin_screener = new StepWiseRNA_VDW_BinScreener();
	if(option[ VDW_rep_screen_info].user() ){ //This is used for post_processing only. Main VDW_rep_screener should be in the samplerer.

		user_input_VDW_bin_screener->set_VDW_rep_alignment_RMSD_CUTOFF(option[ VDW_rep_alignment_RMSD_CUTOFF]());

		user_input_VDW_bin_screener->set_VDW_rep_delete_matching_res( option[ VDW_rep_delete_matching_res ]() );

		user_input_VDW_bin_screener->set_physical_pose_clash_dist_cutoff(option[ VDW_rep_screen_physical_pose_clash_dist_cutoff ]() );

	}

	//////////////////////////////////////////////////////////////

	utility::vector1< std::string > const silent_files_in( option[ in::file::silent ]() );

	utility::vector1< std::string > non_empty_silent_files_in;
	non_empty_silent_files_in.clear();

	/////////////Check for empty silent_files/////////////////
	for(Size n=1; n<=silent_files_in.size(); n++){

		std::string const input_silent_file=silent_files_in[n];

		if(Is_nonempty_input_silent_file(input_silent_file, "empty filtered silent_file since no non-empty sampler silent_file.")){
			std::cout << "adding input_silent_file " << input_silent_file << " to non_empty_silent_files_in " << std::endl;
			non_empty_silent_files_in.push_back(input_silent_file);
		}

	}
	//////////////////////////////////////////////////
	std::string const silent_file_out= option[ out::file::silent  ]();


	////////////////////////////////////////////////////////////////////////////////
	//May 04, 2011. Make sure that silent_file_out doesn't exist before the clustering process.
	//Mainly need this becuase the BIOX2-cluster (Stanford) is not robust...The node containing a slave_job can suddenly become unaviable and slave_job will need to be resubmitted to a new node.
	//For example:
	//		Your job <90107> has been killed because the execution host <node-9-30> is no longer available.
	//		The job will be re-queued and re-run with the same jobId.

	if ( utility::file::file_exists( silent_file_out ) ) {
		std::cout << "WARNING: silent_file_out " << silent_file_out << " already exist! removing..." << std::endl;

		int remove_file_return_value=std::remove( silent_file_out.c_str() );
		std::cout << "remove_file_return_value= " <<  remove_file_return_value << " for std::remove(" << silent_file_out << ")" << std::endl;

		if(remove_file_return_value!=0){
			utility_exit_with_message("remove_file_return_value=" +ObjexxFCL::string_of(remove_file_return_value )+ "!=0 for std::remove(" + silent_file_out + ")" );
		}
	}
	////////////////////////////////////////////////////////////////////////////////


	if(non_empty_silent_files_in.size()==0){
		std::cout << "Early Exit: non_empty_silent_files_in.size()==0, outputting empty clustered outfile " << std::endl;

		std::ofstream outfile;
		outfile.open(silent_file_out.c_str()); //Opening the file with this command removes all prior content..


		outfile << "empty cluster silent_file since all input_silent_file are empty.";
		outfile << " input_silent_file:";
		for(Size n=1; n<=silent_files_in.size(); n++){
			outfile << " " << silent_files_in[n];
		}
		outfile << "\n";

		outfile.flush();
		outfile.close();
		return; //Early return;
	}


	//////////////////////////////////////////////////
	protocols::swa::rna::StepWiseRNA_Clusterer stepwise_rna_clusterer( non_empty_silent_files_in );

	stepwise_rna_clusterer.set_max_decoys( option[ clusterer_num_pose_kept ]() );

	stepwise_rna_clusterer.set_score_diff_cut( option[ score_diff_cut ]() );
	stepwise_rna_clusterer.set_perform_score_diff_cut( option[ clusterer_perform_score_diff_cut ] );

	stepwise_rna_clusterer.set_cluster_radius(	option[ whole_struct_cluster_radius ]()	);
	stepwise_rna_clusterer.set_rename_tags( option[ clusterer_rename_tags ]() );
	stepwise_rna_clusterer.set_job_parameters( job_parameters_COP );
	stepwise_rna_clusterer.set_job_parameters_exist( job_parameters_exist );
	stepwise_rna_clusterer.set_suite_cluster_radius( option[ suite_cluster_radius]() );
	stepwise_rna_clusterer.set_loop_cluster_radius( option[ loop_cluster_radius]() );
	stepwise_rna_clusterer.set_distinguish_pucker( option[ distinguish_pucker]() );
	stepwise_rna_clusterer.set_add_lead_zero_to_tag( option[ add_lead_zero_to_tag ]() );
	stepwise_rna_clusterer.set_quick_alignment( option[ clusterer_quick_alignment ]() );
	stepwise_rna_clusterer.set_align_only_over_base_atoms( option[clusterer_align_only_over_base_atoms]() );
	stepwise_rna_clusterer.set_optimize_memory_usage( option [clusterer_optimize_memory_usage]() );
	stepwise_rna_clusterer.set_keep_pose_in_memory( option [clusterer_keep_pose_in_memory]() );
	stepwise_rna_clusterer.set_two_stage_clustering( option [clusterer_two_stage_clustering]() );
	stepwise_rna_clusterer.set_PBP_clustering_at_chain_closure( option[ PBP_clustering_at_chain_closure]() );
	stepwise_rna_clusterer.set_verbose( option[ VERBOSE ]() );
	stepwise_rna_clusterer.set_skip_clustering( option [skip_clustering]() );
	stepwise_rna_clusterer.set_full_length_loop_rmsd_clustering( option[clusterer_full_length_loop_rmsd_clustering]() );

	stepwise_rna_clusterer.set_filter_virtual_res_list( option[ virtual_res ]() );
	stepwise_rna_clusterer.set_perform_VDW_rep_screen( option[ clusterer_perform_VDW_rep_screen ]() );
	stepwise_rna_clusterer.set_perform_filters( option[ clusterer_perform_filters ]() );
	stepwise_rna_clusterer.set_min_num_south_ribose_filter( option[ clusterer_min_num_south_ribose_filter ]() );
	stepwise_rna_clusterer.set_VDW_rep_screen_info( option[ VDW_rep_screen_info]() );
	stepwise_rna_clusterer.set_user_input_VDW_bin_screener( user_input_VDW_bin_screener );

	stepwise_rna_clusterer.set_ignore_FARFAR_no_auto_bulge_tag( option[ clusterer_ignore_FARFAR_no_auto_bulge_tag ]() );
	stepwise_rna_clusterer.set_ignore_FARFAR_no_auto_bulge_parent_tag( option[ clusterer_ignore_FARFAR_no_auto_bulge_parent_tag ]() );
	stepwise_rna_clusterer.set_ignore_unmatched_virtual_res( option[ clusterer_ignore_unmatched_virtual_res ]() );

	stepwise_rna_clusterer.set_output_pdb( option[ output_pdb ]() );


	stepwise_rna_clusterer.cluster();


	bool const recreate_silent_struct_for_output=( option[ recreate_silent_struct ]() );

	bool const use_best_neighboring_shift_RMSD_for_output=( option[ clusterer_use_best_neighboring_shift_RMSD ]() );

	Size num_special_mode=0;

	if(recreate_silent_struct_for_output) num_special_mode++;

	if(use_best_neighboring_shift_RMSD_for_output) num_special_mode++;

	if(num_special_mode>1) utility_exit_with_message("num_special_mode(" + ObjexxFCL::string_of(num_special_mode) + ")>1");

	if( recreate_silent_struct_for_output){
		//For analysis purposes....for example rescore with a different force-field...change native_pose and etc..

		if(job_parameters_exist==false) utility_exit_with_message("need job_parameters!");

		//core::scoring::ScoreFunctionOP scorefxn=create_scorefxn();
		StepWiseRNA_PoseSetupOP stepwise_rna_pose_setup = setup_pose_setup_class(job_parameters, false /*COPY DOF*/); //This contains the native_pose

		stepwise_rna_clusterer.recalculate_rmsd_and_output_silent_file(silent_file_out,
																														stepwise_rna_pose_setup,
																														option[clusterer_write_score_only]());

	}else if(use_best_neighboring_shift_RMSD_for_output){

		if(job_parameters_exist==false) utility_exit_with_message("need job_parameters!");

		stepwise_rna_clusterer.get_best_neighboring_shift_RMSD_and_output_silent_file(silent_file_out);


	}else{ //default, just output existing silent_struct

		stepwise_rna_clusterer.output_silent_file( silent_file_out );

	}
}


///////////////////////////////////////////////////////////////

void
post_rebuild_bulge_assembly() ///Oct 22, 2011
{

  using namespace core::pose;
  using namespace core::chemical;
  using namespace core::kinematics;
  using namespace core::scoring;
	using namespace protocols::swa::rna;
	using namespace core::io::silent;
	using namespace core::conformation;
	using namespace ObjexxFCL;
	using namespace core::id;
	using namespace core::optimization;

	Output_title_text("Enter post_rebuild_bulge_assembly()", TR );

	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( RNA );

	core::scoring::ScoreFunctionOP scorefxn=create_scorefxn();

	////////////////////////////////////////////////////////////////////////////////////

	if(option[ start_silent ].user()==false ) utility_exit_with_message("User need to pass in start_silent!");
	if(option[ start_tag ].user()==false ) utility_exit_with_message("User need to pass in start_tag!");

	if(option[ in::file::silent ].user()==false ) utility_exit_with_message("User need to pass in in::file::silent!");
	if(option[ in::file::tags ].user()==false ) utility_exit_with_message("User need to pass in in::file::tags!");

	if(option[ out::file::silent ].user()==false ) utility_exit_with_message("User need to pass in out::file::silent!");

	bool const OUTPUT_PDB= option[ output_pdb ]();

	std::string const start_silent_file=option[start_silent ]();
	std::string const start_tag_name   =option[start_tag ]();

	utility::vector1< std::string > const rebuild_silent_file_list=option[in::file::silent ]();
	utility::vector1< std::string > const rebuild_tag_name_list=option[in::file::tags]();

	std::string const rebuild_silent_file=rebuild_silent_file_list[1];
	std::string const rebuild_tag_name=rebuild_tag_name_list[1];

	std::string const out_silent_file = option[ out::file::silent ]();

	std::string const out_tag_name= start_tag_name;
	//std::string const out_tag_name= "R_" + start_tag_name.substr(2, start_tag_name.size()-2)

	///////////////////Create a 'mock' job_parameters only with the parameters called by Output_data()////////////////////////
	StepWiseRNA_JobParametersOP	job_parameters = setup_simple_full_length_rna_job_parameters();
	StepWiseRNA_JobParametersCOP job_parameters_COP( job_parameters );

	Size const total_res=(job_parameters->full_sequence()).size();
	////////////////////////////////////////////////////////////////////////////////////

	//StepWiseRNA_JobParametersOP	job_parameters = setup_rna_job_parameters(false);
	//StepWiseRNA_JobParametersCOP job_parameters_COP( job_parameters );
	//StepWiseRNA_PoseSetupOP stepwise_rna_pose_setup = setup_pose_setup_class(job_parameters, false /*COPY DOF*/); //This contains the native_pose

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  Pose start_pose;
  Pose rebuild_pose;

	import_pose_from_silent_file(start_pose, start_silent_file, start_tag_name );

	import_pose_from_silent_file(rebuild_pose, rebuild_silent_file, rebuild_tag_name );

	///I think that should able to extract the old_score line directly from the pose!

	if(start_pose.total_residue()!=total_res) utility_exit_with_message("start_pose.total_residue()!=total_res");

	if(rebuild_pose.total_residue()!=total_res) utility_exit_with_message("rebuild_pose.total_residue()!=total_res");

	if(OUTPUT_PDB) dump_pdb( start_pose, "start_pose_" + start_tag_name + ".pdb" );
	if(OUTPUT_PDB) dump_pdb( rebuild_pose, "rebuild_pose_" + rebuild_tag_name + ".pdb" );

	/////////Dec 18, 2011: Deal with protonated adenosine. Actually code works fine without this..but still a good consistency test!/////////
	utility::vector1< core::Size > const protonated_H1_adenosine_list = job_parameters->protonated_H1_adenosine_list();


	for(Size seq_num=1; seq_num<=total_res; seq_num++){

		if(has_virtual_rna_residue_variant_type(rebuild_pose, seq_num)) utility_exit_with_message("rebuild_pose has virtaul_residue at seq_num("+string_of(seq_num)+")!");

		if(protonated_H1_adenosine_list.has_value(seq_num)){

			if(rebuild_pose.residue(seq_num).has_variant_type("PROTONATED_H1_ADENOSINE")==false){
				utility_exit_with_message("seq_num(" + string_of(seq_num)+") is in protonated_H1_adenosine_list but rebuild_pose does not have PROTONATED_H1_ADENOSINE variant type!");
			}

			if(has_virtual_rna_residue_variant_type(start_pose, seq_num)){

				if(start_pose.residue(seq_num).has_variant_type("PROTONATED_H1_ADENOSINE")){
					utility_exit_with_message("seq_num(" + string_of(seq_num)+") of start_pose has PROTONATED_H1_ADENOSINE variant type but is a virtual_residue!");
				}

				//This ensures that the Adenosine base have the same atom_list in the start and rebuild pose!
				pose::remove_variant_type_from_pose_residue( rebuild_pose, "PROTONATED_H1_ADENOSINE", seq_num );
				std::cout << "removing PROTONATED_H1_ADENOSINE from seq_num " << seq_num << " of rebuild_pose since this seq_num is a virtual_residue in start_pose!" << std::endl;

			}else{

				if(start_pose.residue(seq_num).has_variant_type("PROTONATED_H1_ADENOSINE")==false){
					utility_exit_with_message("seq_num(" + string_of(seq_num)+") is in protonated_H1_adenosine_list but start_pose does not have PROTONATED_H1_ADENOSINE variant type!");
				}

			}

		}else{

			if(rebuild_pose.residue(seq_num).has_variant_type("PROTONATED_H1_ADENOSINE")){
				utility_exit_with_message("seq_num(" + string_of(seq_num)+") is not in protonated_H1_adenosine_list but rebuild_pose has PROTONATED_H1_ADENOSINE variant type!");
			}

		}

	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	core::pose::MiniPose const mini_rebuild_pose= *( core::pose::MiniPoseOP( new core::pose::MiniPose( rebuild_pose ) ) );


	if(mini_rebuild_pose.total_residue()!=total_res){
		std::cout << "mini_rebuild_pose.total_residue()= " << mini_rebuild_pose.total_residue() << std::endl;
		std::cout << "total_res= " << total_res << std::endl;
		utility_exit_with_message("mini_rebuild_pose.total_residue()!=total_res");
	}

	utility::vector1< utility::vector1< std::string > > const & chunk_atom_names_list = mini_rebuild_pose.atom_names_list();


	std::map< core::Size, core::Size > res_map;

	for ( Size n = 1; n <= start_pose.total_residue(); n++ ) {
		res_map[ n ] = n;

		std::cout << "chunk_seq_num= " << n << "| chunk_atom_names_list[chunk_seq_num].size()= " << chunk_atom_names_list[n].size() << std::endl;

	}


	/////Copy the conformation but nothing else. No energy and no cache data (having cache data can cause problem with column_name order in output silent_file!)//////
	/////OK...this should also copy the virtual_types and the fold_tree?//////////////////////////////////////////////////////////////////////////////////////////////
	//Pose output_pose=start_pose;

	ConformationOP copy_conformation = new Conformation();

	(*copy_conformation)=start_pose.conformation();

	pose::Pose output_pose;
	output_pose.set_new_conformation( copy_conformation );

	if(output_pose.total_residue()!=total_res) utility_exit_with_message("output_pose.total_residue()!=total_res");

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	//NEW_copy_dofs( output_pose, mini_rebuild_pose, res_map );
	copy_dofs_match_atom_names( output_pose, mini_rebuild_pose, res_map); //Dec 28, 2011..STILL NEED TO VERIFY THAT THIS WORKS PROPERLY!


	//OK THE copy_dofs seem to correctly position every atom except for the OVU1, OVL1 and OVL2 at the chainbreak(s). This is becuase the rebuild_pose does not necessaringly have the same cutpoint position as the start_pose

	protocols::rna::assert_phosphate_nomenclature_matches_mini( output_pose );//Just to be safe

	utility::vector1< Size > virtual_res_list;

	////////////////////////////////////////////////////////////////////////////////////
	for(Size seq_num=1; seq_num<=output_pose.total_residue(); seq_num++){
		///Oct 24, 2011 Should remove the virtual_variant type before minimizing to prevent artificial clashes!
		///Also currently have plan to not score torsional_potential if OVU1, OVL1 and OVL2 if corresponding atoms (O3', O5' and P1) are virtual (not implemented yet!)
		if(has_virtual_rna_residue_variant_type(output_pose, seq_num)){
			virtual_res_list.push_back(seq_num);
			remove_virtual_rna_residue_variant_type(output_pose, seq_num);
		}
	}

	Output_seq_num_list("virtual_res_list= ", virtual_res_list, TR, 30 );
	////////////////////////////////////////////////////////////////////////////////////

	for(Size seq_num=1; seq_num<=output_pose.total_residue(); seq_num++){

		if(output_pose.residue(seq_num).has_variant_type("CUTPOINT_LOWER")){

			if(seq_num==output_pose.total_residue()) utility_exit_with_message("seq_num==output_pose.total_residue() has CUTPOINT_LOWER variant_type!");

			if(output_pose.fold_tree().is_cutpoint( seq_num )==false) utility_exit_with_message("seq_num (" + string_of(seq_num) +") is not a cutpoint!");

			if(output_pose.residue(seq_num+1).has_variant_type("CUTPOINT_UPPER")==false){
				utility_exit_with_message("seq_num+1 (" + string_of(seq_num+1) +") does not have CUTPOINT_UPPER variant_type!");
			}


			/* //Cutpoint_open is not necessarily a virtual_res!
			if(output_pose.residue(seq_num).has_variant_type("VIRTUAL_RNA_RESIDUE")==false){
				utility_exit_with_message("seq_num (" + string_of(seq_num) +") does not have VIRTUAL_RNA_RESIDUE variant_type!");
			}

			if(output_pose.residue(seq_num+1).has_variant_type("VIRTUAL_RNA_RESIDUE_UPPER")==false){
				 utility_exit_with_message("seq_num+1 (" + string_of(seq_num+1) +") does not have VIRTUAL_RNA_RESIDUE_UPPER variant_type!");
			}
			*/

			/*///////////NOTE:This part didn't change result compare to copy_torsion alone////////////////////////
			pose::remove_variant_type_from_pose_residue( output_pose, chemical::CUTPOINT_LOWER, seq_num   );
			pose::remove_variant_type_from_pose_residue( output_pose, chemical::CUTPOINT_UPPER, seq_num+1 );

			Correctly_position_cutpoint_phosphate_torsions( output_pose, seq_num, false );

			pose::add_variant_type_to_pose_residue( output_pose, chemical::CUTPOINT_LOWER, seq_num   );
			pose::add_variant_type_to_pose_residue( output_pose, chemical::CUTPOINT_UPPER, seq_num+1 );
			///////////NOTE:This part didn't change result compare to copy_torsion alone////////////////////////*/

			////////////////////Copy the chain_break torsions///////////////////////////////////
			/////Important to copy torsion before calling rna_loop_closer///////////////////////
			/////Since copy_torsion with position OVU1, OVL1, OVL2 near the correct solution////
			/////This minimize changes in the other torsions (i.e. beta and gamma of 3'-res)////

			//alpha 3'-res aligns OVU1 of 3'-res to O3' of 5'-res
			//This doesn't seem to change the torsion? Possibly becuase OVU1 is already repositioned when COPY_DOF set the OP1 and OP2 pos?
			//output_pose.set_torsion( TorsionID( seq_num+1, id::BB,  1 ), rebuild_pose.residue(seq_num+1).mainchain_torsion(1)  );
			Real const output_alpha_torsion=output_pose.residue(seq_num+1).mainchain_torsion(1);
			Real const rebuild_alpha_torsion=output_pose.residue(seq_num+1).mainchain_torsion(1);

			Real const abs_diff=std::abs(output_alpha_torsion-rebuild_alpha_torsion);
			std::cout << "std::abs(output_alpha_torsion-rebuild_alpha_torsion)= " << abs_diff << std::endl;
			if(abs_diff>0.000001) utility_exit_with_message("std::abs(output_alpha_torsion-rebuild_alpha_torsion)>0.000001");

			for(Size n=5; n<=6; n++){
				//epsilon 5'-res aligns OVL1 of 5'-res ot P1  of 3'-res.
				//zeta of 5'-res aligns OVL2 of 5'-res to O5' of 3'-res
				output_pose.set_torsion( TorsionID( seq_num, id::BB,  n ), rebuild_pose.residue(seq_num).mainchain_torsion(n) );
			}
			////////////////////Copy the chain_break torsions////////////////////////////////////

			//This does slightly improve the CCD torsion (for TEST_MODE, no_minimize case!)
			//protocols::rna::RNA_LoopCloser rna_loop_closer;
			//rna_loop_closer.set_three_prime_alpha_only(true);
			//rna_loop_closer.apply( output_pose, seq_num );
			//////////////////////////////////////////////////////////////////////////////////////

			AtomTreeMinimizer minimizer;
	    float const dummy_tol( 0.00000025);
	    bool const use_nblist( true );
	    MinimizerOptions options( "dfpmin", dummy_tol, use_nblist, false, false );
	    options.nblist_auto_update( true );

			core::kinematics::MoveMap mm;

			mm.set_bb( false );
			mm.set_chi( false );
			mm.set_jump( false );

			mm.set( TorsionID( seq_num  , id::BB,  5 ) , true );	//5'-res epsilon
			mm.set( TorsionID( seq_num  , id::BB,  6 ) , true );	//5'-res zeta
			mm.set( TorsionID( seq_num+1, id::BB,  1 ) , true );	//3'-res alpha

			minimizer.run( output_pose, mm, *(scorefxn), options );



			//////////////////////////////////////////////////////////////////////////////////////

		}
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////
	for(Size seq_num=1; seq_num<=output_pose.total_residue(); seq_num++){
		if( virtual_res_list.has_value(seq_num) ){
			apply_virtual_rna_residue_variant_type(output_pose, seq_num, true /*apply_check*/);
		}
	}


	//////////////OK setup_native_pose does two things, slice_pose to match working length and align native to working_pose////
	//////////////However, for full pose length, don't need to slices. ////////////////////////////////////////////////////////
	//////////////Furthermore, if write_score_only==False then automatically inside Output_data()//////////////////////////////
	//stepwise_rna_pose_setup->set_verbose(true); //New OPTION, Mar 22
	//stepwise_rna_pose_setup->setup_native_pose( output_pose ); //Setup native_pose;


	/////////////////////////////////////////////////////////////////////////////////////////////////////////////

	SilentFileData silent_file_data;

	(*scorefxn)(output_pose);

	Output_data(silent_file_data, out_silent_file, out_tag_name, false /*write_score_only*/, output_pose, job_parameters->working_native_pose(),  job_parameters);

}


///////////////////////////////////////////////////////////////
void
ensure_directory_for_out_silent_file_exists(){

	if ( option[ out::file::silent ].user() ){

		std::string outfile =  option[ out::file::silent]();

		std::ofstream outstream;
	 	outstream.open( outfile.c_str() ); // for writing

		if (outstream.fail()){
			// wow, this is tortuous -- libgen.h has dirname, but requires and output C-style char.
			std::cout <<  "Could not create silent file output " << outfile << " so making the directory!" << std::endl;
			char * outfile_char = strdup( outfile.c_str() );
			char * outdir =  dirname( outfile_char );
			std::stringstream mkdir_command;
			mkdir_command << "mkdir -p " << outdir;
			system( mkdir_command.str().c_str() );
		} else {
			outstream.close();
			std::remove( outfile.c_str() ); // note that this removes the prior outfile if it exists...
		}
	}

}


///////////////////////////////////////////////////////////////
void*
my_main( void* )
{

	clock_t const my_main_time_start( clock() );

  using namespace basic::options;

	std::string algorithm_input = option[algorithm];

	ensure_directory_for_out_silent_file_exists();

	if (algorithm_input=="cluster_old" or algorithm_input=="rna_cluster"){
		swa_rna_cluster();
	}	else if (algorithm_input=="rna_resample_test" or algorithm_input=="rna_sample"){
	  swa_rna_sample();
	}	else if (algorithm_input=="rna_sample_virtual_ribose"){
	  rna_sample_virtual_ribose();
	} else if (algorithm_input=="post_rebuild_bulge_assembly"){
		post_rebuild_bulge_assembly();
	}	else if (algorithm_input=="filter_combine_long_loop"){
		filter_combine_long_loop();
	} else {
		utility_exit_with_message("Invalid User-specified algorithm ("+ algorithm_input +")!");
	}

	protocols::viewer::clear_conformation_viewers();

	std::cout << "Total time took to run algorithm (" << algorithm_input << "): " << static_cast<Real>( clock() - my_main_time_start ) / CLOCKS_PER_SEC << " seconds." << std::endl;

	std::cout << "JOB_SUCCESSFULLY_COMPLETED" << std::endl;

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

	//////////////General/////////////////////////////
	NEW_OPT( graphic, "Turn graphic on/off", true);
	NEW_OPT( Real_parameter_one, "free_variable for testing purposes ", 0.0);
	NEW_OPT( distinguish_pucker, "distinguish pucker when cluster:both in sampler and clusterer", true);
	NEW_OPT( output_pdb, "output_pdb: If true, then will dump the pose into a PDB file at different stages of the stepwise assembly process.", false); //Sept 24, 2011
	NEW_OPT( algorithm, "Specify algorithm to execute", "");

	//////////////Job_Parameters///////////
	NEW_OPT( sample_res, "residues to build, the first element is the actual sample res while the other are the bulge residues", blank_size_vector );
	NEW_OPT( input_res, "Residues already present in starting pose_1", blank_size_vector );
	NEW_OPT( input_res2, "Residues already present in starting  pose_2", blank_size_vector );
	NEW_OPT( missing_res, "Residues missing in starting pose_1, alternative to input_res", blank_size_vector );
	NEW_OPT( missing_res2, "Residues missing in starting pose_2, alternative to input_res2", blank_size_vector );
	NEW_OPT( rmsd_res, "residues that will be use to calculate rmsd (for clustering as well as RMSD to native_pdb if specified)", blank_size_vector );
	NEW_OPT( alignment_res , "align_res_list", blank_string_vector );
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
	NEW_OPT( sample_virtual_ribose_list, "optional: sample_virtual_ribose_list", blank_string_vector); //July 20, 2011
	NEW_OPT( force_syn_chi_res_list, "optional: sample only syn chi for the res in sampler.", blank_size_vector); //April 29, 2011
	NEW_OPT( force_north_ribose_list, "optional: sample only north ribose for the res in sampler.", blank_size_vector); //April 29, 2011
	NEW_OPT( force_south_ribose_list, "optional: sample only south ribose for the res in sampler.", blank_size_vector); //April 29, 2011
	NEW_OPT( protonated_H1_adenosine_list, "optional: protonate_H1_adenosine_list", blank_size_vector); //May 02, 2011

	//////////////Pose setup///////
	NEW_OPT( job_queue_ID, " swa_rna_sample()/combine_long_loop mode: Specify the tag pair in filter_output_filename to be read in and imported (start from 0!)" , 0);

	///////////////Sampler////////////
	NEW_OPT( combine_long_loop_mode, " Sampler: combine_long_loop_mode " , false);
	NEW_OPT( rebuild_bulge_mode, "rebuild_bulge_mode", false);
	NEW_OPT( floating_base , " floating_base ", false ); //DO NOT CHANGE TO TRUE, since single-nucleotide sampling need this to be false! April 9th, 2011
	NEW_OPT( sampler_native_rmsd_screen, "native_rmsd_screen ResidueSampler", false );
	NEW_OPT( sampler_native_screen_rmsd_cutoff, "sampler_native_screen_rmsd_cutoff", 2.0 );
	NEW_OPT( sampler_num_pose_kept, "optional: set_num_pose_kept by ResidueSampler", 108 );

	//////////////Minimizer////////////
	NEW_OPT( num_pose_minimize, "optional: set_num_pose_minimize by Minimizer", 999999 );

	//////////////Clusterer////////////
	NEW_OPT( clusterer_num_pose_kept, "optional: Num_pose_kept by the clusterer", 1000 );
	NEW_OPT( suite_cluster_radius , " individual_suite_cluster_radius ", 999.99); 							//IMPORTANT, DO NOT CHANGE DEFAULT VALUE!
	NEW_OPT( loop_cluster_radius , " loop_cluster_radius ", 999.99); 													//IMPORTANT, DO NOT CHANGE DEFAULT VALUE!
	NEW_OPT( clusterer_full_length_loop_rmsd_clustering, "use the full_length_rmsd function to calculate loop_rmsd", false); //April 06, 2011: Should switch to true for all length_full clustering steps.

	//////////VDW_bin_screener//////////
	NEW_OPT( VDW_rep_screen_info, "VDW_rep_screen_info to create VDW_rep_screen_bin (useful when building loop from large poses)", blank_string_vector); //Jun 9, 2010
	NEW_OPT( VDW_rep_alignment_RMSD_CUTOFF, "use with VDW_rep_screen_info", 0.001); //Nov 12, 2010
	NEW_OPT( VDW_rep_delete_matching_res, "delete residues in VDW_rep_pose that exist in the working_pose", blank_string_vector); //Feb 20, 2011
	NEW_OPT( VDW_rep_screen_physical_pose_clash_dist_cutoff, "The distance cutoff for VDW_rep_screen_with_physical_pose", 1.2); //March 23, 2011

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
	NEW_OPT( exclude_alpha_beta_gamma_sampling, "Speed up the debug epsilon south sugar mode", false);
	NEW_OPT( debug_epsilon_south_sugar_mode, "Check why when epsilon is roughly -160 and pucker is south, energy is not favorable", false);
	// FCC: Not doing anythin now... Just for consistency with swa_analytical_closure
	/////////
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
	NEW_OPT( sampler_assert_no_virt_ribose_sampling, "sampler_assert_no_virt_ribose_sampling", false); //July 28, 2011
	NEW_OPT( centroid_screen, "centroid_screen", true);
	NEW_OPT( allow_base_pair_only_centroid_screen, "allow_base_pair_only_centroid_screen", false); //This only effect floating base sampling + dinucleotide.. deprecate option
	NEW_OPT( VDW_atr_rep_screen, "classic VDW_atr_rep_screen", true);
	NEW_OPT( sampler_perform_o2star_pack, "perform O2' hydrogen packing inside StepWiseRNA_ResidueSampler", true );
	NEW_OPT( allow_bulge_at_chainbreak, "Allow sampler to replace chainbreak res with virtual_rna_variant if it looks have bad fa_atr score.", true );

	//////////////Minimizer////////////
	NEW_OPT( minimize_and_score_native_pose, "minimize_and_score_native_pose ", false); //Sept 15, 2010
	NEW_OPT( minimize_and_score_sugar, "minimize and sugar torsion+angle? and include the rna_sugar_close_score_term ", true); //Sept 15, 2010
	NEW_OPT( minimizer_perform_minimize, "minimizer_perform_minimize", true );
	NEW_OPT( minimizer_output_before_o2star_pack, "minimizer_output_before_o2star_pack", false);
	NEW_OPT( minimizer_perform_o2star_pack, "perform O2' hydrogen packing inside StepWiseRNA_Minimizer", true); //Jan 19, 2012
	NEW_OPT( minimizer_rename_tag, "Reorder and rename the tag by the energy_score", true); //March 15, 2012
	NEW_OPT( native_edensity_score_cutoff, "native_edensity_score_cutoff", -1.0 ); //Fang's electron density code

	//////////////Clusterer ///////////////////////
	NEW_OPT( clusterer_two_stage_clustering, "Cluster is two stage..using triangle inequaility to speed up clustering", false); //Change to false on April 9th 2011
	NEW_OPT( clusterer_keep_pose_in_memory, "reduce memory usage for the clusterer", true); //Aug 6, 2010
	NEW_OPT( clusterer_quick_alignment, "quick alignment during clusterer...only work if the alignment residues are fixed ", false);
	NEW_OPT( clusterer_align_only_over_base_atoms, "align_only_over_base_atoms in clusterer alignment ", true); //Add option in Aug 20, 2011
	NEW_OPT( clusterer_optimize_memory_usage, "clusterer_optimize_memory_usage ", false);
	NEW_OPT( score_diff_cut, "score difference cut for clustering", 1000000.0 );							//IMPORTANT, DO NOT CHANGE DEFAULT VALUE!
	NEW_OPT( clusterer_perform_score_diff_cut, "score difference cut for clustering", false ); //IMPORTANT, LEAVE THE DEFAULT AS FALSE!

	//////////////Clusterer Post-Analyis////////////
	NEW_OPT( recreate_silent_struct, "Special mode to recreate_silent_struct for clusterer output...for analysis purposes", false);
	NEW_OPT( clusterer_write_score_only, "clusterer_write_score_only/ only effect recreate_silent_struct mode  ", false);
	NEW_OPT( clusterer_perform_filters, "Other filters such as for specific puckers and chi conformations", false); //June 13, 2011
	NEW_OPT( clusterer_perform_VDW_rep_screen, "filter for VDW clash with the VDW_rep_pose in clusterer ", false); //March 20, 2011
	NEW_OPT( clusterer_use_best_neighboring_shift_RMSD, "Special mode for clusterer output...for analysis purposes", false); //Dec 10, 2011
	NEW_OPT( clusterer_ignore_FARFAR_no_auto_bulge_tag, "clusterer_ignore_FARFAR_no_auto_bulge_tag", false); //Sept 06, 2011
	NEW_OPT( clusterer_ignore_FARFAR_no_auto_bulge_parent_tag, "clusterer_ignore_FARFAR_no_auto_bulge_parent_tag", false); //Sept 06, 2011
	NEW_OPT( clusterer_ignore_unmatched_virtual_res, "clusterer_ignore_unmatched_virtual_res", false); //Sept 06, 2011
	NEW_OPT( clusterer_min_num_south_ribose_filter, "clusterer_min_num_south_ribose_filter", 0); //Oct 02, 2011
	NEW_OPT( add_lead_zero_to_tag, "Add lead zero to clusterer output tag ", false);
	NEW_OPT( skip_clustering, "keep every pose, no clustering", false );
	NEW_OPT( clusterer_rename_tags , "clusterer_rename_tags", true);
	NEW_OPT( whole_struct_cluster_radius , " whole_struct_cluster_radius ", 0.5); //IMPORTANT DO NOT CHANGE
	NEW_OPT ( constraint_chi, "Constrain the chi angles", false );
	NEW_OPT ( rm_virt_phosphate, "Remove virtual phosphate patches during minimization", false );
	NEW_OPT ( choose_random, "ask swa residue sampler for a random solution", false );


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



