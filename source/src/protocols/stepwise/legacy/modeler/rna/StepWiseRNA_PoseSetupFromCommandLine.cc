// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/legacy/modeler/rna/StepWiseRNA_PoseSetupFromCommandLine.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu, based on Parin Sripakdeevong's work

// libRosetta headers
#include <core/types.hh>

#include <core/chemical/AA.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/ResidueMatcher.hh>
#include <core/conformation/Conformation.hh>

#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/chemical/rna/RNA_FittedTorsionInfo.hh>
#include <core/chemical/rna/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/func/CharmmPeriodicFunc.hh>
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
#include <basic/options/option.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>
#include <basic/options/keys/stepwise.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/full_model.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
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
#include <core/pose/rna/util.hh>


#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/util.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>
#include <core/import_pose/import_pose.hh>

//////////////////////////////////////////////////////////
#include <protocols/viewer/viewers.hh>
#include <protocols/farna/util.hh>
#include <protocols/farna/movers/RNA_LoopCloser.hh>
#include <protocols/stepwise/modeler/rna/StepWiseRNA_OutputData.hh>
#include <protocols/stepwise/modeler/rna/StepWiseRNA_CombineLongLoopFilterer.hh>
#include <protocols/stepwise/modeler/rna/StepWiseRNA_CombineLongLoopFilterer.fwd.hh>
#include <protocols/stepwise/modeler/rna/sugar/VirtualSugarSampler.hh>
#include <protocols/stepwise/legacy/modeler/rna/StepWiseRNA_Minimizer.hh>
#include <protocols/stepwise/legacy/modeler/rna/StepWiseRNA_PoseSetup.fwd.hh>
#include <protocols/stepwise/legacy/modeler/rna/StepWiseRNA_PoseSetup.hh>
#include <protocols/stepwise/legacy/modeler/rna/StepWiseRNA_WorkingParametersSetup.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.hh>
#include <protocols/stepwise/legacy/modeler/rna/StepWiseRNA_Clusterer.hh>
#include <protocols/stepwise/legacy/modeler/rna/util.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_VDW_BinChecker.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_BaseCentroidChecker.hh>
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

#ifdef WIN32
#include <direct.h>
#define GetCurrentDir _getcwd

#else
#define GetCurrentDir getcwd
#include <libgen.h>
#include <unistd.h>
#endif

#include <list>
#include <stdio.h>
#include <math.h>

using namespace core;
using namespace protocols;
using namespace ObjexxFCL;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace protocols::stepwise::modeler;
using namespace protocols::stepwise::modeler::rna;
using namespace protocols::stepwise::legacy::modeler::rna;
using utility::vector1;

typedef  numeric::xyzMatrix< Real > Matrix;


static THREAD_LOCAL basic::Tracer TR( "protocols.stepwise.legacy.modeler.rna.StepWiseRNA_PoseSetupFromCommandLine" );

///////////////////////////////////////////////////////////////////////////
//
// All utility functions for setting up job parameters, poses, etc. for SWA
// and ERRASER. This used to all be in stepwise.rna_main, the main application!
//
///////////////////////////////////////////////////////////////////////////
namespace protocols {
namespace stepwise {
namespace legacy {
namespace modeler {
namespace rna {

//////////////////////////////////////////////////////////////////////////////////////
//Apply chi angle constraint to the purines --
void apply_chi_cst( core::pose::Pose & pose, core::pose::Pose const & ref_pose ) {

	using namespace core::conformation;
	using namespace core::id;
	using namespace core::scoring::constraints;
	using namespace core::chemical;
	using namespace core::chemical::rna;

	Size const nres = pose.size();
	ConstraintSetOP cst_set( new ConstraintSet );
	for ( Size i = 1; i <= nres; ++i ) {
		ResidueType const & res = pose.residue_type( i );
		if ( res.is_RNA() && ( res.aa() == na_rad || res.aa() == na_rgu ) ) {
			Real const chi = numeric::conversions::radians( ref_pose.torsion( TorsionID( i, id::CHI, 1 ) ) );
			core::scoring::func::FuncOP chi_cst_func( new core::scoring::func::CharmmPeriodicFunc( chi, 1.0, 1.0 ) );
			AtomID const atom1 ( res.atom_index( "C2'" ), i );
			AtomID const atom2 ( res.atom_index( "C1'" ), i );
			AtomID const atom3 ( res.is_purine() ? res.atom_index( "N9" ) : res.atom_index( "N1" ), i );
			AtomID const atom4 ( res.is_purine() ? res.atom_index( "C4" ) : res.atom_index( "C2" ), i );
			cst_set->add_constraint( ConstraintCOP( ConstraintOP( new DihedralConstraint( atom1, atom2, atom3, atom4, chi_cst_func ) ) ) );
		}
	}
	pose.constraint_set( cst_set );
}


//////////////////////////////////////////////////////////////////////////////////////
std::string
get_working_directory(){

	char cCurrentPath[FILENAME_MAX];
	std::string current_directory_string;

	if ( !GetCurrentDir( cCurrentPath, sizeof( cCurrentPath ) ) ) {
		utility_exit_with_message( "!GetCurrentDir( cCurrentPath, sizeof( cCurrentPath ) )" );
	}

	std::stringstream ss;
	ss << cCurrentPath;
	ss >> current_directory_string;

	return current_directory_string;
}

///////////////////////////////////////////////////////////////////////////////////////
utility::vector1< core::Size >
get_fixed_res( core::Size const nres ){

	utility::vector1< Size > blank_size_vector;

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace protocols::stepwise::modeler::rna;

	utility::vector1< Size > actual_fixed_res_list;
	actual_fixed_res_list.clear();

	utility::vector1< core::Size > const fixed_res_list    = option[ OptionKeys::stepwise::fixed_res ]();
	utility::vector1< core::Size > const minimize_res_list = option[ OptionKeys::stepwise::rna::minimize_res ]();

	if ( fixed_res_list.size() != 0 && minimize_res_list.size() != 0 ) {
		utility_exit_with_message( "User Cannot specify both  fixed_res and minimize_res!" );
	}

	if ( fixed_res_list.size() != 0  ) {
		actual_fixed_res_list = fixed_res_list;

	} else if ( minimize_res_list.size() != 0 ) {

		for ( Size seq_num = 1; seq_num <= nres; seq_num++ ) {
			if ( minimize_res_list.has_value( seq_num ) ) continue;
			actual_fixed_res_list.push_back( seq_num );
		}

	} else { //here I am being a little stringent and require user specify one of these option. Could just return empty list...
		//utility_exit_with_message( "User did not specify both fixed res and minimize_res!" );
		TR << " WARNING! User did not specify either fixed res and minimize_res!" << std::endl;
		TR << " Assuming -input_res gives fixed_res" << std::endl;
		actual_fixed_res_list = option[ in::file::input_res ]();
	}

	return actual_fixed_res_list;
}

//////////////////////////////////////////////////////////////////////////////////////
bool
is_nonempty_input_silent_file( std::string const & input_silent_file, std::string const & exit_key_string ){

	TR << "Checking that input_silent_file " << input_silent_file << " contain actual silent_structs or the correct exit_key_string" << std::endl;

	std::ifstream infile;
	infile.open( input_silent_file.c_str() );

	if ( infile.fail() ) {
		utility_exit_with_message( "Error! \"" + input_silent_file + "\" could not be opened!" );
	} else {
		std::cout << "Open \"" << input_silent_file << "\" successful!" << std::endl;
	}

	std::string line;

	getline( infile, line );
	bool found_line = infile.good();
	runtime_assert( found_line );

	size_t found_substring = line.find( exit_key_string );

	if ( found_substring != std::string::npos ) {
		std::cout << "input_silent_file: " << input_silent_file << " contain no silent struct" << std::endl;
		std::cout << line << std::endl;

		//consistency check:////////////////////////////////////////////////////////////////////////////////////////////////////
		std::string next_line;
		getline( infile, next_line );
		bool found_next_line = infile.good();
		if ( found_next_line ) std::cout << "input silent_file contain more than one line! next_line = " << next_line << std::endl;
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		return false;
	} else {
		return true;
	}
}
//////////////////////////////////////////////////////////////////////////////////////

utility::vector1< core::Size >
get_input_res( core::Size const nres, std::string const & pose_num ){


	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace protocols::stepwise::modeler::rna;

	utility::vector1< core::Size > input_res_list;
	utility::vector1< core::Size > missing_res_list;

	if ( pose_num == "1" ) {
		input_res_list = option[ in::file::input_res ]();
		missing_res_list = option[ OptionKeys::stepwise::rna::missing_res ]();
	} else if ( pose_num == "2" ) {
		input_res_list = option[ OptionKeys::stepwise::input_res2 ]();
		missing_res_list = option[ OptionKeys::stepwise::rna::missing_res2 ]();
	} else {
		utility_exit_with_message( "Invalid pose_num " + pose_num + ", must by either 1 or 2 !" );
	}


	if ( input_res_list.size() != 0 && missing_res_list.size() != 0 ) {
		utility_exit_with_message( "User Cannot specify both input_res" + pose_num + " and missing_res" + pose_num + "!" );
	}

	utility::vector1< core::Size > actual_input_res_list;
	actual_input_res_list.clear();

	if ( input_res_list.size() != 0 ) {
		actual_input_res_list = input_res_list;

	} else if ( missing_res_list.size() != 0 ) {

		for ( Size seq_num = 1; seq_num <= nres; seq_num++ ) {
			if ( missing_res_list.has_value( seq_num ) ) continue;
			actual_input_res_list.push_back( seq_num );
		}

	} else { //did not specify both input_res and missing_res, return empty list
		TR.Debug << "user did not specify either input_res" << pose_num << " and missing_res" << pose_num << std::endl;
	}

	return actual_input_res_list;

}

//////////////////////////////////////////////////////////////////////////////////////

utility::vector1< std::string >
get_silent_file_tags(){

	using namespace protocols::stepwise::modeler::rna;

	bool tags_from_command_line = false;
	bool tags_from_filterer_outfile = false;

	utility::vector1< std::string > input_silent_file_tags;

	if ( option[ in::file::tags ].user() ) {
		tags_from_command_line = true;
		input_silent_file_tags = option[ in::file::tags ]();
	}

	if ( option[ OptionKeys::stepwise::rna::job_queue_ID ].user() && option[ OptionKeys::stepwise::rna::filter_output_filename ].user() ) {

		stepwise::modeler::rna::output_title_text( "importing tag from filter_outfile", TR );
		bool combine_long_loop = option[ OptionKeys::stepwise::rna::combine_long_loop_mode ]();
		bool combine_helical = option[ OptionKeys::stepwise::rna::combine_helical_silent_file ]();
		runtime_assert( combine_long_loop || combine_helical );

		tags_from_filterer_outfile = true;

		std::string const filtered_tag_file = option[ OptionKeys::stepwise::rna::filter_output_filename ]();

		std::ifstream infile;
		infile.open( filtered_tag_file.c_str() );

		if ( infile.fail() ) {
			utility_exit_with_message( "Error! \"" + filtered_tag_file + "\" could not be opened!" );
		} else {
			TR << "Open \"" << filtered_tag_file << "\" successful!" << std::endl;
		}


		//Be careful here... job_queue_ID start from ZERO!
		int const queue_ID = option[ OptionKeys::stepwise::rna::job_queue_ID ]();
		int ID = 0;

		TR << "queue_ID = " << queue_ID << std::endl;

		std::string tag_pair_string;

		bool found_queue_ID = false;

		while ( getline( infile, tag_pair_string ) ) {

			if ( queue_ID == ID ) {
				found_queue_ID = true;
				break;
			}

			ID++;
		}

		//Warning queue_ID start at ZERO!
		if ( found_queue_ID == false ) utility_exit_with_message( "found_queue_ID == false, queue_ID = " + string_of( queue_ID ) + " num_tag_string_in_file = " + string_of( ID ) );


		TR << "import silent_file_tags: " << tag_pair_string << " from filter_output_filename = " << filtered_tag_file << std::endl;

		infile.close();

		utility::vector1< std::string > const line_list = tokenize( tag_pair_string, " \t\n\f\v" ); //Oct 19, 2010..now filterer_outfile contain other terms.

		input_silent_file_tags.clear();
		input_silent_file_tags.push_back( line_list[1] );
		input_silent_file_tags.push_back( line_list[2] );

		stepwise::modeler::rna::output_title_text( "", TR );

	}


	if ( ( tags_from_command_line == false ) && ( tags_from_filterer_outfile == false ) ) {
		utility_exit_with_message( "( tags_from_command_line == false ) && ( tags_from_filterer_outfile == false )" );
	}

	if ( ( tags_from_command_line == true ) && ( tags_from_filterer_outfile == true ) ) {
		utility_exit_with_message( "( tags_from_command_line == true ) && ( tags_from_filterer_outfile == true )" );
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

	Size num_score_weight_file = 0;

	if ( option[ basic::options::OptionKeys::score::weights ].user() ) {
		score_weight_file = option[ basic::options::OptionKeys::score::weights ]();
		TR << "User passed in score:weight option: " << score_weight_file << std::endl;
		num_score_weight_file++;
	}


	if ( num_score_weight_file == 0 ) {
		//rna_loop_hires_04092010.wts is same as 5X_linear_quarter_fa_stack_and_adjust_bulge_ss_benchmark.wts
		//change default from single_strand_benchmark to 5X_linear_chainbreak_single_strand_benchmark on May 24, 2010
		//change default to 5X_linear_quarter_fa_stack_and_adjust_bulge_ss_benchmark.wts" on April 9th, 2011
		//score_weight_file="rna_loop_hires_04092010.wts";
		utility_exit_with_message( "User to need to pass in score:weights" ); //Remove the default weight on Sept 28, 2011 Parin S.
	}

	if ( num_score_weight_file > 1 ) {
		TR << "num_score_weight_file ( inputted by user ) = " << num_score_weight_file << std::endl;
		utility_exit_with_message( "num_score_weight_file > 1" );
	}

	core::scoring::ScoreFunctionOP scorefxn = get_score_function();


	if ( ! option[ OptionKeys::stepwise::rna::minimize_and_score_sugar]() ) {
		TR << "WARNING minimize_and_score_sugar is false, SET rna_sugar_close weight to 0.0 " << std::endl;
		scorefxn->set_weight( rna_sugar_close, 0 );

		//Sept 16, 2010. Thought about include a very small weight for rna_sugar_close so that column # will not change. HOWEVER this significant change the minimization results!
		//scorefxn->set_weight( rna_sugar_close, 0.000000000001 );
	}

	TR.Debug << "---------score function weights----------" << std::endl;
	scorefxn->show( TR.Debug );
	TR.Debug << "-----------------------------------------" << std::endl;


	return scorefxn;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///Oct 28, 2011:This function is an alternative way to setup the working_parameters////////////////////////////////////////////////////
///The original working_parameters was designed primarily for use during the SWA building step/////////////////////////////////////////
///However since then, many other other SWA and FARFAR class make use of it especially the full pose length version////////////////
///Hence the original working_parameters setup has become clumbersome for most usage (containing many options and independencies)//////
///The new job_parameter setup function is a work in progress and will be updated as needed.///////////////////////////////////////
///For example right now (Oct 28, 2011), it only contains parameters needed by output_data()///////////////////////////////////////

stepwise::modeler::working_parameters::StepWiseWorkingParametersOP
setup_simple_full_length_rna_working_parameters(){

	using namespace protocols::stepwise::modeler::rna;
	using namespace ObjexxFCL;
	using namespace core::pose;
	using namespace core::chemical;

	ResidueTypeSetCOP rsd_set( core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );

	/////////////////////////////////////////////////////
	if ( !option[ in::file::fasta].user() ) utility_exit_with_message( "Must supply in::file::fasta!" );
	if ( !option[ OptionKeys::stepwise::rna::rmsd_res ].user() ) utility_exit_with_message( "Must supply rmsd_res!" );
	if ( !option[ OptionKeys::stepwise::rna::alignment_res ].user() ) utility_exit_with_message( "Must supply alignment_res!" );
	if ( !option[ OptionKeys::stepwise::rna::global_sample_res_list ].user() ) utility_exit_with_message( "Must supply global_sample_res_list!" );

	/////////////Read in sequence.///////////////////////

	std::string const fasta_file = option[ in::file::fasta ]()[1];
	std::string const full_sequence = core::sequence::read_fasta_file_and_concatenate( fasta_file );
	core::Size const nres = full_sequence.length();

	/////////////////////////////////////////////////////

	stepwise::modeler::working_parameters::StepWiseWorkingParametersOP working_parameters( new stepwise::modeler::working_parameters::StepWiseWorkingParameters );

	working_parameters->set_is_simple_full_length_job_params( true ); //FORGOT THIS! ONLY INCLUDED ON MARCH 03, 2012!
	//LUCKILY, BEFORE MARCH 03, 2012, only called is_simple_full_length_job_params() for utility_exit_with_message() check in the following three functions of StepWiseRNA_Clusterer.cc:
	//i. initialize_VDW_rep_checker(),
	//ii. recalculate_rmsd_and_output_silent_file()
	//iii. get_best_neighboring_shift_RMSD_and_output_silent_file()
	/////////////////////////////////////////////////////

	working_parameters->set_output_extra_RMSDs( option[ OptionKeys::stepwise::rna::output_extra_RMSDs ]() );

	utility::vector1< core::Size > is_working_res( nres, 1 ); //All res belong to 'mock pose' 1
	working_parameters->set_is_working_res( is_working_res );
	working_parameters->set_full_sequence( full_sequence ); //working_sequence is automatical init after BOTH is_working_res and full_sequence is initialized

	if ( option[ OptionKeys::stepwise::rna::rmsd_res ].user() ) {
		working_parameters->set_calc_rms_res( option[ OptionKeys::stepwise::rna::rmsd_res ]() );
	} else {
		TR << "Warning! rmsd_res not specified. Assuming sample_res is rmsd_res" << std::endl;
		working_parameters->set_calc_rms_res( option[ OptionKeys::full_model::sample_res ]() );
	}
	working_parameters->set_gap_size( 0 );

	if ( option[ OptionKeys::stepwise::rna::native_alignment_res ].user() ) {
		//If working_native_align_res is not specified that most function will internally default to using working_best_alignment_list in its place
		working_parameters->set_native_alignment( option[ OptionKeys::stepwise::rna::native_alignment_res ]() );
		working_parameters->set_working_native_alignment( option[ OptionKeys::stepwise::rna::native_alignment_res ]() );
	}
	/////////////////////////////////////////////////////
	//Oct 31, 2011: Better to let code raise error if the there is a statement asking for the partition_definition from the job_params later in the run!
	//ObjexxFCL::FArray1D_bool partition_definition( nres, false ); //All res belong in partition 0!
	//working_parameters->set_partition_definition( partition_definition ); //this is a useful decomposition.
	/////////////////////////////////////////////////////

	std::map< core::Size, core::Size > full_to_sub;
	std::map< core::Size, bool > is_prepend_map;

	utility::vector1< Size > input_res;
	utility::vector1< Size > input_res2;

	for ( Size seq_num = 1; seq_num <= full_sequence.size(); seq_num++ ) {
		input_res.push_back( seq_num );
		full_to_sub[ seq_num ] = seq_num;
		is_prepend_map[seq_num] = false; //all append..arbitrary choice
	}

	working_parameters->set_full_to_sub( full_to_sub );  //res_map
	working_parameters->set_is_prepend_map( is_prepend_map );
	working_parameters->set_is_prepend( is_prepend_map[working_parameters->actually_moving_res()] );

	utility::vector1< utility::vector1< Size > > input_res_vectors;
	input_res_vectors.push_back( input_res );
	input_res_vectors.push_back( input_res2 );
	working_parameters->set_input_res_vectors( input_res_vectors );

	/////////////////////////////////////////////////////
	utility::vector1< Size > working_moving_res_list;
	working_moving_res_list.push_back( full_sequence.size() ); //arbitrary choice, choose residue of pose

	working_parameters->set_working_moving_res_list( working_moving_res_list ); //This sets actually_moving_res()
	//////////////////////////////////////////////////////////////////////////////////////////////////////////

	utility::vector1< std::string > alignment_res_string_list = option[ OptionKeys::stepwise::rna::alignment_res ]();
	utility::vector1< core::Size > best_alignment_list;

	for ( Size n = 1; n <= alignment_res_string_list.size(); n++ ) {

		utility::vector1< std::string > alignments_res_string = tokenize( alignment_res_string_list[n], "-" );
		utility::vector1< core::Size >  alignment_list;

		for ( Size ii = 1; ii <= alignments_res_string.size(); ii++ ) {
			alignment_list.push_back( string_to_int( alignments_res_string[ii] ) );
		}

		if ( alignment_list.size() > best_alignment_list.size() ) best_alignment_list = alignment_list;
	}

	working_parameters->set_working_best_alignment( best_alignment_list );

	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	working_parameters->set_working_fixed_res(  get_fixed_res( nres ) );
	working_parameters->set_global_sample_res_list(  option[ OptionKeys::stepwise::rna::global_sample_res_list ]() );
	working_parameters->set_force_syn_chi_res_list(  option[ OptionKeys::full_model::rna::force_syn_chi_res_list]() );
	working_parameters->set_force_anti_chi_res_list(  option[ OptionKeys::full_model::rna::force_anti_chi_res_list]() );
	working_parameters->set_force_north_sugar_list(  option[ OptionKeys::full_model::rna::force_north_sugar_list ]() );
	working_parameters->set_force_south_sugar_list(  option[ OptionKeys::full_model::rna::force_south_sugar_list ]() );
	working_parameters->set_protonated_H1_adenosine_list( option[ OptionKeys::stepwise::rna::protonated_H1_adenosine_list ]() );

	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	PoseOP native_pose_OP;
	if ( option[ in::file::native ].user() ) {
		native_pose_OP = PoseOP( new Pose );
		import_pose::pose_from_file( *native_pose_OP, *rsd_set, option[ in::file::native ](), core::import_pose::PDB_file );
		core::pose::rna::make_phosphate_nomenclature_matches_mini( *native_pose_OP );

		utility::vector1< core::Size > const native_virtual_res_list = option[ OptionKeys::stepwise::rna::native_virtual_res]();

		for ( Size n = 1; n <= native_virtual_res_list.size(); n++ ) {
			core::pose::rna::apply_virtual_rna_residue_variant_type( ( *native_pose_OP ), native_virtual_res_list[n], false /*apply_check*/ );
		}

		working_parameters->set_working_native_pose( native_pose_OP );

	}
	//////////////////////////////////////////////////////////////////////////////////////////////////////////

	print_WorkingParameters_info( working_parameters, "simple_full_length_working_parameters", TR, true /*is_simple_full_length_WP*/ );

	return working_parameters;

}


//////////////////////////////////////////////////////////////////////////////////////
stepwise::modeler::working_parameters::StepWiseWorkingParametersOP
setup_rna_working_parameters( bool check_for_previously_closed_cutpoint_with_input_pose /* = false */ ){
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
	//              setup_rna_working_parameters?
	//
	//       The setup_rna_working_parameters class is responsible for figuring out
	//       the fold_tree. For simple motifs, such as a single hairpin or a
	//       single double-stranded motif, setup_rna_working_parameters can figure
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
	//       imported_silent_file() function insetup_rna_working_parameters class for
	//       futher details).
	//
	//       Lastly, note that the setup_rna_working_parameters class will still be able
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

	using namespace protocols::stepwise::modeler::rna;
	using namespace ObjexxFCL;
	///////////////////////////////
	// Read in sequence.
	if ( !option[ in::file::fasta ].user() ) utility_exit_with_message( "Must supply in::file::fasta!" );
	std::string const fasta_file = option[ in::file::fasta ]()[1];
	std::string const full_sequence = core::sequence::read_fasta_file_and_concatenate( fasta_file );
	core::Size const nres = full_sequence.length();

	if ( !option[ OptionKeys::full_model::sample_res ].user() ) utility_exit_with_message( "Must supply sample_res!" );

	/////////////////////////////////////////////////////
	Size const cutpoint_closed_ = option[ OptionKeys::full_model::cutpoint_closed ].user() ? option[ OptionKeys::full_model::cutpoint_closed ]()[1] : 0;
	StepWiseWorkingParametersSetup stepwise_rna_working_parameters_setup(
		option[ OptionKeys::full_model::sample_res ](), /*the first element of moving_res_list is the modeler_res*/
		full_sequence,
		get_input_res( nres, "1" ),
		get_input_res( nres, "2" ),
		option[ OptionKeys::full_model::cutpoint_open ](),
		cutpoint_closed_ );
	stepwise_rna_working_parameters_setup.set_simple_append_map( option[ OptionKeys::stepwise::rna::simple_append_map]() );
	stepwise_rna_working_parameters_setup.set_add_virt_res_as_root( option[ OptionKeys::stepwise::rna::add_virt_root]() );
	stepwise_rna_working_parameters_setup.set_allow_fixed_res_at_moving_res( option[ OptionKeys::stepwise::rna::allow_fixed_res_at_moving_res ]() ); //Hacky just to get Hermann Duplex working. Need to called before set_fixed_res
	stepwise_rna_working_parameters_setup.set_fixed_res( get_fixed_res( nres ) );
	stepwise_rna_working_parameters_setup.set_terminal_res( option[ OptionKeys::full_model::rna::terminal_res ]() );

	if ( option[ OptionKeys::stepwise::rna::rmsd_res ].user() ) {
		stepwise_rna_working_parameters_setup.set_calc_rms_res( option[ OptionKeys::stepwise::rna::rmsd_res ]() );
	} else {
		TR << "Warning! rmsd_res not specified. Assuming sample_res is rmsd_res" << std::endl;
		stepwise_rna_working_parameters_setup.set_calc_rms_res( option[ OptionKeys::full_model::sample_res ]() );
	}

	// jump_point_pairs is string of pairs  "1-16 8-9", assumed for now to be connected by dashes.  See note in StepWiseWorkingParametersSetup.cc
	stepwise_rna_working_parameters_setup.set_jump_point_pair_list( option[ OptionKeys::stepwise::rna::jump_point_pairs ]() ); //Important!: Need to be called after set_fixed_res
	stepwise_rna_working_parameters_setup.set_force_user_defined_jumps( option[ OptionKeys::stepwise::rna::force_user_defined_jumps ]() );

	//  Alignment_res is a string vector to allow user to specify multiple possible alignments. Note that 1-6 7-12 is different from 1-6-7-12. See note in StepWiseWorkingParametersSetup.cc
	stepwise_rna_working_parameters_setup.set_alignment_res( option[ OptionKeys::stepwise::rna::alignment_res ]() );
	stepwise_rna_working_parameters_setup.set_filter_user_alignment_res( option[ OptionKeys::stepwise::rna::filter_user_alignment_res ]() );
	stepwise_rna_working_parameters_setup.set_native_alignment_res( option[ OptionKeys::stepwise::rna::native_alignment_res ]() );

	stepwise_rna_working_parameters_setup.set_global_sample_res_list( option[ OptionKeys::stepwise::rna::global_sample_res_list ]() ); //March 20, 2011

	stepwise_rna_working_parameters_setup.set_force_syn_chi_res_list( option[ OptionKeys::full_model::rna::force_syn_chi_res_list]() ); //April 29, 2011
	stepwise_rna_working_parameters_setup.set_force_anti_chi_res_list( option[ OptionKeys::full_model::rna::force_anti_chi_res_list]() ); //April 29, 2011
	stepwise_rna_working_parameters_setup.set_force_north_sugar_list( option[ OptionKeys::full_model::rna::force_north_sugar_list ]() ); //April 29, 2011
	stepwise_rna_working_parameters_setup.set_force_south_sugar_list( option[ OptionKeys::full_model::rna::force_south_sugar_list ]() ); //April 29, 2011
	stepwise_rna_working_parameters_setup.set_protonated_H1_adenosine_list( option[ OptionKeys::stepwise::rna::protonated_H1_adenosine_list ]() ); //May 02, 2011

	stepwise_rna_working_parameters_setup.set_allow_chain_boundary_jump_partner_right_at_fixed_BP( option[ OptionKeys::stepwise::rna::allow_chain_boundary_jump_partner_right_at_fixed_BP ]() ); //Hacky just to get Square RNA working.

	stepwise_rna_working_parameters_setup.set_output_extra_RMSDs( option[ OptionKeys::stepwise::rna::output_extra_RMSDs ]() );
	stepwise_rna_working_parameters_setup.set_floating_base( option[ OptionKeys::stepwise::rna::floating_base ]() ||
		option[ OptionKeys::stepwise::rna::floating_base_anchor_res ].user() );
	stepwise_rna_working_parameters_setup.set_floating_base_anchor_res( option[ OptionKeys::stepwise::rna::floating_base_anchor_res ]() );
	if ( option[ OptionKeys::stepwise::rna::floating_base_anchor_res ]() ) runtime_assert( option[ OptionKeys::stepwise::rna::force_user_defined_jumps ]() );
	stepwise_rna_working_parameters_setup.set_rebuild_bulge_mode( option[ basic::options::OptionKeys::stepwise::rna::rebuild_bulge_mode]() );
	stepwise_rna_working_parameters_setup.set_sample_both_sugar_base_rotamer( option[ basic::options::OptionKeys::stepwise::rna::sample_both_sugar_base_rotamer]() );

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
				utility_exit_with_message( "silent_files_in.size( " + string_of( silent_files_in.size() ) + " ) != input_tags.size( " + string_of( input_tags.size() ) + " )" );
			}
		}
		stepwise_rna_working_parameters_setup.set_input_tags( input_tags );
		stepwise_rna_working_parameters_setup.set_silent_files_in( silent_files_in );
	}
	///////////////////////////////////////////////////////

	stepwise_rna_working_parameters_setup.apply();

	return stepwise_rna_working_parameters_setup.working_parameters();

}


////////////////////////////////////////////////////////////////////////////////////////////////////
void
setup_copy_DOF_input( StepWiseRNA_PoseSetupOP & stepwise_rna_pose_setup ){

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
			utility_exit_with_message( "silent_files_in.size() != input_tags.size()" );
		}
	}

	if ( option[ in::file::s ].user() ) {
		// Then any pdbs that need to be read in from disk.
		utility::vector1< std::string > const pdb_tags_from_disk( option[ in::file::s ]() );
		for ( Size n = 1; n <= pdb_tags_from_disk.size(); n++ ) {
			input_tags.push_back( pdb_tags_from_disk[ n ] );
		}
	}

	if ( input_tags.size() > 2 ) {
		utility_exit_with_message( "input_tags.size() > 2!!" );
	}

	TR << "Input structures for COPY DOF" << std::endl;
	for ( Size n = 1; n <= input_tags.size(); n++ ) {
		if ( n <= silent_files_in.size() ) {
			TR << "silent_file tag = " << input_tags[n] << " silent_file = " << silent_files_in[n] << std::endl;
		} else {
			TR << "input_tag = " << input_tags[n] << std::endl;
		}
	}

	//////////////////////////////////////////////////////////////////////////////////////////
	stepwise_rna_pose_setup->set_input_tags( input_tags );
	stepwise_rna_pose_setup->set_silent_files_in( silent_files_in );


}

//////////////////////////////////////////////////////////////////////////////////////////
StepWiseRNA_PoseSetupOP
setup_pose_setup_class( stepwise::modeler::working_parameters::StepWiseWorkingParametersOP & working_parameters, bool const copy_DOF /*= true*/ ){

	using namespace core::pose;
	using namespace core::chemical;
	using namespace core::kinematics;
	using namespace core::scoring;
	using namespace protocols::stepwise::modeler::rna;

	ResidueTypeSetCOP rsd_set( core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );

	// Read in native_pose.
	PoseOP native_pose;
	if ( option[ in::file::native ].user() ) {
		native_pose = PoseOP( new Pose );
		import_pose::pose_from_file( *native_pose, *rsd_set, option[ in::file::native ]() , core::import_pose::PDB_file);
		TR.Debug << "native_pose->fold_tree(): " << native_pose->fold_tree();
		TR.Debug << "native_pose->annotated_sequence( true ): " << native_pose->annotated_sequence( true ) << std::endl;
		core::pose::rna::make_phosphate_nomenclature_matches_mini( *native_pose );
	}


	// StepWiseRNA_PoseSetup stepwise_rna_pose_setup( pdb_tags, silent_files_in, working_parameters);

	StepWiseRNA_PoseSetupOP stepwise_rna_pose_setup( new StepWiseRNA_PoseSetup( working_parameters ) );
	stepwise_rna_pose_setup->set_copy_DOF( copy_DOF );

	if ( copy_DOF == true ) {
		setup_copy_DOF_input( stepwise_rna_pose_setup );
	}

	stepwise_rna_pose_setup->set_virtual_res( option[ basic::options::OptionKeys::full_model::virtual_res ]() );
	stepwise_rna_pose_setup->set_bulge_res( option[ basic::options::OptionKeys::full_model::rna::bulge_res ]() );
	stepwise_rna_pose_setup->set_native_pose( native_pose );
	stepwise_rna_pose_setup->set_native_virtual_res( option[ basic::options::OptionKeys::stepwise::rna::native_virtual_res]() );
	stepwise_rna_pose_setup->set_output_pdb( option[ OptionKeys::stepwise::dump ]() );
	stepwise_rna_pose_setup->set_use_phenix_geo ( option[ basic::options::OptionKeys::rna::corrected_geo ]() );

	return stepwise_rna_pose_setup;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Conventionally, SWA only outputted structures from a full enumeration, and
//  its own classes ('StepWiseRNA_Minimizer') handled the final output.
// As we start to explore Monte Carlo variants, its useful to run the whole mode
//  nstruct times. In that case, we only output the final lowest energy pose into
//  the silent file. We can keep other 'temporary' silent files along the way (which
//  are now named in the swa_silent_file object).
//
bool
get_tag_and_silent_file_for_struct( std::string & swa_silent_file,
	std::string & out_tag,
	Size const & n,
	bool const & multiple_shots,
	std::string const & silent_file ){

	swa_silent_file = silent_file; // default

	if ( multiple_shots ) {
		runtime_assert( option[ OptionKeys::stepwise::choose_random ]() );

		out_tag = "S_"+lead_zero_string_of( n, 6 );
		//  swa_silent_file = out_tag + "_" + swa_silent_file;
		swa_silent_file = "";

		std::map< std::string, bool > tag_is_done;
		static bool init_tag_is_done( false );
		if ( !init_tag_is_done ) {
			tag_is_done = core::io::silent::initialize_tag_is_done( silent_file );
			init_tag_is_done = true;
		}

		if ( tag_is_done[ out_tag ] ) {
			TR << "Already done: " << out_tag << std::endl;
			return false;
		}
	}
	return true;
}

///////////////////////////////////////////////////////////////
void
ensure_directory_for_out_silent_file_exists(){

	if ( option[ out::file::silent ].user() ) {

		std::string outfile =  option[ out::file::silent]();

		std::ofstream outstream;
		outstream.open( outfile.c_str() ); // for writing

		if ( outstream.fail() ) {
			// wow, this is tortuous -- libgen.h has dirname, but requires and output C-style char.
			TR <<  "Could not create silent file output " << outfile << " so making the directory!" << std::endl;
#ifdef WIN32
			//char * outdir; memory leak if we're just gonna exit right after?
			utility_exit_with_message( "protocols/stepwise/legacy/modeler/rna/StepWiseRNA_PoseSetupFromCommandLine.cc dirname is not implemented under Windows!" );
#else
			char * outfile_char = strdup( outfile.c_str() );
			char * outdir =  dirname( outfile_char );

			std::stringstream mkdir_command;
			mkdir_command << "mkdir -p " << outdir;
			int return_code = system( mkdir_command.str().c_str() );
			if ( return_code != 0 ) {
				TR.Error << "Could not make directory! Error code: " << return_code << std::endl;
			}
			// AMW cppcheck: deleting outdir if we're not on win32
			delete[] outdir;
#endif
		} else {
			outstream.close();
			std::remove( outfile.c_str() ); // note that this removes the prior outfile if it exists...
		}
	}

}


} //rna
} //modeler
} //legacy
} //stepwise
} //protocols
