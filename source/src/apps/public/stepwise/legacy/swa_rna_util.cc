// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file swa_rna_main.cc
/// @author Parin Sripakdeevong (sripakpa@stanford.edu), Rhiju Das (rhiju@stanford.edu)


// libRosetta headers
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeSelector.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueMatcher.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
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

#include <core/pose/util.hh>
#include <core/pose/util.tmpl.hh>
#include <core/pose/annotated_sequence.hh>

#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.fwd.hh>
#include <core/import_pose/import_pose.hh>
//////////////////////////////////////////////////
#include <basic/options/option.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh> // for option[ out::file::silent  ] and etc.
#include <basic/options/keys/in.OptionKeys.gen.hh> // for option[ in::file::tags ] and etc.
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/rna/util.hh>
#include <core/chemical/rna/RNA_FittedTorsionInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh> //for EnergyMap
#include <core/scoring/EnergyMap.fwd.hh> //for EnergyMap
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/rotamer_trials.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <utility/file/file_sys_util.hh> //Parin April 15, 2012
#include <utility/excn/Exceptions.hh>

#include <numeric/xyzVector.hh>
#include <numeric/conversions.hh>
//////////////////////////////////////////////////////////
#include <protocols/viewer/viewers.hh>
#include <protocols/farna/RNA_Minimizer.hh>
#include <protocols/farna/util.hh>
#include <protocols/stepwise/modeler/output_util.hh>
#include <protocols/stepwise/modeler/rna/util.hh>
#include <protocols/stepwise/modeler/rna/StepWiseRNA_ResidueInfo.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.hh>
#include <protocols/stepwise/legacy/modeler/rna/util.hh>
#include <protocols/farna/RNA_LoopCloser.hh>
#include <protocols/farna/RNA_LoopCloser.fwd.hh>

#include <core/pose/rna/RNA_BaseDoubletClasses.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/BinarySilentStruct.hh>
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

//Added by Parin
//#include <core/scoring/ScoreType.hh>
#include <list>
#include <stdio.h>
#include <math.h>


using namespace core;
using namespace core::pose::rna;
using namespace protocols;
using namespace ObjexxFCL;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using utility::vector1;
using io::pdb::dump_pdb;
using namespace protocols::stepwise::modeler;
using namespace protocols::stepwise::modeler::rna;
using namespace protocols::stepwise::legacy::modeler::rna;

typedef  numeric::xyzMatrix< Real > Matrix;

static thread_local basic::Tracer TR( "swa_rna_util" );

OPT_KEY( String,  algorithm )
OPT_KEY( Real, surrounding_radius )
OPT_KEY( IntegerVector, sample_res )
OPT_KEY( String, rebuild_sequence )
OPT_KEY( Boolean, reset_o2prime_torsion )
OPT_KEY( StringVector, rmsd_res_pairs )
OPT_KEY( StringVector, alignment_res_pairs )
OPT_KEY( Real, alignment_RMSD_CUTOFF )
OPT_KEY( String,  tag_name )
OPT_KEY( String,  output_silent_file )
OPT_KEY( StringVector, list_of_virtual_res )
OPT_KEY( RealVector, list_of_energy )
OPT_KEY( IntegerVector, virtual_res )
OPT_KEY( IntegerVector, virtual_sugar )
OPT_KEY( Boolean, align_only_over_base_atoms )
OPT_KEY( IntegerVector, additional_slice_res )
OPT_KEY( IntegerVector, native_virtual_res )
OPT_KEY( String, native_tag_name )
OPT_KEY( StringVector, decoy_tag_name )
OPT_KEY( Boolean, dump )
OPT_KEY( StringVector, input_tag_list )
OPT_KEY( Boolean, graphic )
OPT_KEY( Boolean, minimizer_deriv_check )
OPT_KEY( String,  minimizer_min_type )
OPT_KEY( Boolean, minimizer_skip_o2prime_trials )
OPT_KEY( Boolean, minimizer_perform_minimizer_run )

//////////////////////////////////////////////////////////////////////////////////////
core::scoring::ScoreFunctionOP
create_scorefxn(){ //Copy from rna_swa_test.cc on Oct 11, 2011

	using namespace core::scoring;


	std::string score_weight_file;

	Size num_score_weight_file = 0;

	if ( option[ basic::options::OptionKeys::score::weights ].user() ) {
		score_weight_file = option[ basic::options::OptionKeys::score::weights ]();
		std::cout << "User passed in score:weight option: " << score_weight_file << std::endl;
		num_score_weight_file++;
	}


	if ( num_score_weight_file == 0 ) {
		score_weight_file = "stepwise/rna/rna_loop_hires_04092010.wts";
		std::cout << "Using default score_weight_file = " << score_weight_file << std::endl;
	}

	if ( num_score_weight_file > 1 ) {
		std::cout << "num_score_weight_file ( inputted by user ) = " << num_score_weight_file << std::endl;
		utility_exit_with_message( "num_score_weight_file > 1" );
	}

	core::scoring::ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( score_weight_file );

	std::cout << "---------score function weights----------" << std::endl;
	scorefxn->show( std::cout );
	std::cout << "-----------------------------------------" << std::endl;


	return scorefxn;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///Jan 01, 2012: SHOULD INTEGRATE THIS WITH THE VERSION in protocols/stepwise/modeler/rna/util.hh
void
align_pose_general( core::pose::Pose const & static_pose, std::string const & static_tag, core::pose::Pose & moving_pose, std::string const & moving_tag, utility::vector1< std::pair< Size, Size > > const & alignment_res_pair_list, bool const base_only ){


	bool found_non_virtual_base = false;
	for ( Size n = 1; n <= alignment_res_pair_list.size(); n++ ) {

		Size const static_seq_num = alignment_res_pair_list[n].first;
		Size const moving_seq_num = alignment_res_pair_list[n].second;

		if ( is_virtual_base( static_pose.residue( static_seq_num ) ) || is_virtual_base( moving_pose.residue( moving_seq_num ) ) ) continue;

		found_non_virtual_base = true; //ok found a non-virtual base nucleotide that can be used for alignment
		break;
	}

	if ( found_non_virtual_base == false ) {
		for ( Size n = 1; n <= alignment_res_pair_list.size(); n++ ) {

			Size const static_seq_num = alignment_res_pair_list[n].first;
			Size const moving_seq_num = alignment_res_pair_list[n].second;

			std::cout << "static_seq_num = " << static_seq_num << " moving_seq_num = " << moving_seq_num;
			output_boolean( "  is_virtual_base( " + static_tag + " ):", is_virtual_base( static_pose.residue( static_seq_num ) ), TR );
			output_boolean( "  is_virtual_base( " + moving_tag + " ):", is_virtual_base( moving_pose.residue( moving_seq_num ) ), TR );
			std::cout << std::endl;
		}
		std::string error_message = "Error in aligning " + moving_tag + " to " + static_tag + ". No non - virtual_base in working_best_alignment to align the poses!";
		std::cout << error_message << std::endl;
		utility_exit_with_message( error_message );
	}

	////////////////////////////////////create the alignment map////////////////////////////////////////////////////////////////////
	id::AtomID_Map < id::AtomID > atom_ID_map;
	pose::initialize_atomid_map( atom_ID_map, moving_pose, id::BOGUS_ATOM_ID );

	std::string const static_sequence = static_pose.sequence();
	std::string const moving_sequence = moving_pose.sequence();


	for ( Size n = 1; n <= alignment_res_pair_list.size(); n++ ) {

		Size const static_seq_num = alignment_res_pair_list[n].first;
		Size const moving_seq_num = alignment_res_pair_list[n].second;

		if ( ( static_sequence[static_seq_num - 1] ) != ( moving_sequence[moving_seq_num - 1] ) ) {
			std::cout << "static_seq_num = " << static_seq_num << " static_sequence = " << static_sequence << " static_sequence[static_seq_num - 1] = " << static_sequence[static_seq_num - 1] << std::endl;
			std::cout << "moving_sequence = " << moving_seq_num << " moving_sequence = " << moving_sequence << " moving_sequence[moving_seq_num - 1] = " << moving_sequence[moving_seq_num - 1] << std::endl;
			utility_exit_with_message( "( static_sequence[static_seq_num - 1] ) != ( moving_sequence[moving_seq_num - 1] )" );
		}

		setup_suite_atom_id_map( moving_pose.residue( moving_seq_num ), static_pose.residue( static_seq_num ),  atom_ID_map, base_only );

	}

	core::scoring::superimpose_pose( moving_pose, static_pose, atom_ID_map );

	if ( check_for_messed_up_structure( moving_pose, moving_tag ) == true ) {
		std::string error_message = "Error in aligning " + moving_tag + " to " + static_tag + "!";
		std::cout << error_message << std::endl;
		utility_exit_with_message( moving_tag + " is messed up ...this is probably an alignment problem" );
	};


}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
core::Real
check_alignment_RMSD_cutoff( core::pose::Pose const & static_pose, std::string const & static_tag, core::pose::Pose & moving_pose, std::string const & moving_tag, utility::vector1< std::pair< Size, Size > > const & alignment_res_pair_list, bool const base_only, core::Real const alignment_RMSD_cutoff ){


	Size total_atom_count = 0;
	Real total_sum_sd = 0.0;

	for ( Size n = 1; n <= alignment_res_pair_list.size(); n++ ) {

		Size const static_seq_num = alignment_res_pair_list[n].first;
		Size const moving_seq_num = alignment_res_pair_list[n].second;


		if ( is_virtual_base( static_pose.residue( static_seq_num ) ) || is_virtual_base( static_pose.residue( static_seq_num ) ) ) continue;

		Size atom_count = 0;
		Real sum_sd = 0.0;

		if ( base_only ) {
			base_atoms_square_deviation( static_pose, moving_pose, static_seq_num, moving_seq_num, atom_count, sum_sd, false /*verbose*/, false /*ignore_virtual_atom*/ );
		} else {
			suite_square_deviation( static_pose, moving_pose, false /*is_prepend*/, static_seq_num, moving_seq_num, atom_count, sum_sd, false /*verbose*/, false /*ignore_virtual_atom*/ );
		}

		sum_sd = sum_sd/( atom_count );
		Real rmsd = sqrt( sum_sd );

		if ( atom_count == 0 ) rmsd = 99.99; //This is different from suite_rmsd function..!!


		if ( rmsd > alignment_RMSD_cutoff ) { //change on Sept 26, 2010..problem arise when use this in non-long-loop mode...
			std::cout << "rmsd = " << rmsd  << " is greater than " << alignment_RMSD_cutoff << " Angstrom between res " << moving_seq_num << " of moving_pose and res " << static_seq_num << " of static_pose" << std::endl;
			std::cout << "static_tag = " << static_tag << " moving_tag = " << moving_tag << std::endl;
			utility_exit_with_message( "rmsd > alignment_RMSD_cutoff!" );
		}

		total_atom_count += atom_count;
		total_sum_sd += sum_sd;

	}

	total_sum_sd = total_sum_sd/( total_atom_count );
	Real const all_rmsd = sqrt( total_sum_sd );

	return all_rmsd;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///Jan 01, 2012: SHOULD INTEGRATE THIS WITH THE VERSION in protocols/stepwise/modeler/rna/util.hh
Real
full_length_rmsd_over_reside_list_general( pose::Pose const & pose_one, pose::Pose const & pose_two, utility::vector1 < std::pair< Size, Size > > const & rmsd_res_pair_list, bool const verbose, bool const ignore_virtual_atom ){

	using namespace ObjexxFCL;

	if ( verbose ) {
		output_title_text( "Enter full_length_rmsd_over_residue_list_general function", TR );
		output_boolean( "ignore_virtual_atom = ", ignore_virtual_atom, TR ); std::cout << std::endl;
		output_pair_size( rmsd_res_pair_list, "rmsd_res_pair_list = ", TR );
	}

	utility::vector1< Size > calc_rms_res_one;
	utility::vector1< Size > calc_rms_res_two;

	for ( Size i = 1; i <= rmsd_res_pair_list.size(); i++ ) {

		Size const seq_num_one = rmsd_res_pair_list[i].first;
		Size const seq_num_two = rmsd_res_pair_list[i].second;

		calc_rms_res_one.push_back( seq_num_one );
		calc_rms_res_two.push_back( seq_num_two );

	}

	Size atom_count = 0;
	Real sum_sd = 0;

	std::string const sequence_one = pose_one.sequence();
	std::string const sequence_two = pose_two.sequence();

	for ( Size i = 1; i <= rmsd_res_pair_list.size(); i++ ) {

		Size const seq_num_one = rmsd_res_pair_list[i].first;
		Size const seq_num_two = rmsd_res_pair_list[i].second;

		if ( ( sequence_one[seq_num_one - 1] ) != ( sequence_two[seq_num_two - 1] ) ) {
			std::cout << "seq_num_one = " << seq_num_one << " sequence_one = " << sequence_one << " sequence[seq_num_one - 1] = " << sequence_one[seq_num_one - 1] << std::endl;
			std::cout << "seq_num_two = " << seq_num_two << " sequence_two = " << sequence_two << " sequence[seq_num_two - 1] = " << sequence_two[seq_num_two - 1] << std::endl;
			utility_exit_with_message( "( sequence_one[seq_num_one - 1] ) != ( sequence_two[seq_num_two - 1] )" );
		}


		bool is_prepend = false;
		bool both_pose_res_is_virtual = false;

		if ( pose_one.residue( seq_num_one ).has_variant_type( core::chemical::VIRTUAL_RNA_RESIDUE ) &&
				pose_two.residue( seq_num_two ).has_variant_type( core::chemical::VIRTUAL_RNA_RESIDUE ) ) {
			both_pose_res_is_virtual = true;
		}

		if ( ( seq_num_one + 1 ) <= pose_one.total_residue() ) {
			if ( pose_one.residue( seq_num_one ).has_variant_type( core::chemical::VIRTUAL_RNA_RESIDUE ) ) {
				if ( ! pose_one.residue( seq_num_one + 1 ).has_variant_type( core::chemical::VIRTUAL_PHOSPHATE ) ) { //consistency_check
					utility_exit_with_message( "pose_one's seq_num_one = " + string_of( seq_num_one ) +
						" is a virtual res but seq_num_one + 1 is not a virtual_res_upper!" );
				}
			}

			if ( pose_two.residue( seq_num_two ).has_variant_type( core::chemical::VIRTUAL_RNA_RESIDUE ) ) {
				if ( ! pose_two.residue( seq_num_two + 1 ).has_variant_type( core::chemical::VIRTUAL_PHOSPHATE ) ) { //consistency_check
					utility_exit_with_message( "pose_two's seq_num_two = " + string_of( seq_num_two ) +
						" is a virtual res but seq_num_two + 1 is not a virtual_res_upper!" );
				}
			}
		}

		if ( verbose ) {
			std::cout << "seq_num_one = " << seq_num_one << " seq_num_two = " << seq_num_two;
			output_boolean( " is_prepend = ", is_prepend, TR );
			output_boolean( " both_pose_res_is_virtual = ", both_pose_res_is_virtual, TR ); std::cout << std::endl;
		}

		if ( both_pose_res_is_virtual ) continue;

		//add atom in the suites to atom_count
		//add sd of each atom to sum_sd
		suite_square_deviation( pose_one, pose_two, is_prepend, seq_num_one, seq_num_two, atom_count, sum_sd, false, ignore_virtual_atom );


		if ( ( ( seq_num_one + 1 ) <= pose_one.total_residue() ) != ( ( seq_num_two + 1 ) <= pose_two.total_residue() ) ) {
			std::cout << "seq_num_one = " << seq_num_one << " pose_one.total_residue() = " <<  pose_one.total_residue() << std::endl;
			std::cout << "seq_num_two = " << seq_num_two << " pose_two.total_residue() = " <<  pose_two.total_residue() << std::endl;
			utility_exit_with_message( "( ( seq_num_one + 1 ) <= pose_one.total_residue() ) != ( ( seq_num_two + 1 ) <= pose_two.total_residue() )" );
		}

		if ( ( seq_num_one + 1 ) <= pose_one.total_residue() ) {

			if ( ( sequence_one[( seq_num_one + 1 ) - 1] ) != ( sequence_two[( seq_num_two + 1 ) - 1] ) ) {
				std::cout << "( seq_num_one + 1 ) - 1 = " << seq_num_one << " sequence_one = " << sequence_one << " sequence[( seq_num_one + 1 ) - 1] = " << sequence_one[( seq_num_one + 1 ) - 1] << std::endl;
				std::cout << "( seq_num_two + 1 ) - 1 = " << seq_num_two << " sequence_two = " << sequence_two << " sequence[( seq_num_two + 1 ) - 1] = " << sequence_two[( seq_num_two + 1 ) - 1] << std::endl;
				utility_exit_with_message( "( sequence_one[( seq_num_one + 1 ) - 1] ) != ( sequence_two[( seq_num_two + 1 ) - 1] )" );
			}

			bool is_phosphate_edge_res_one = false;
			bool is_phosphate_edge_res_two = false;

			if ( calc_rms_res_one.has_value( seq_num_one + 1 ) == false ) is_phosphate_edge_res_one = true;
			if ( calc_rms_res_two.has_value( seq_num_two + 1 ) == false ) is_phosphate_edge_res_two = true;

			if ( is_phosphate_edge_res_one != is_phosphate_edge_res_two ) {
				std::cout << "seq_num_one + 1 = " << seq_num_one + 1 << " seq_num_two + 1 = " << seq_num_two + 1 << std::endl;
				utility_exit_with_message( "is_phosphate_edge_res_one != is_phosphate_edge_res_two" );
			}

			if ( is_phosphate_edge_res_one ) {

				if ( verbose ) std::cout << "is_phosphate_edge_res = true! for seq_num_one = " << seq_num_one << " seq_num_two = " << seq_num_two << "" << std::endl;

				phosphate_square_deviation( pose_one, pose_two, seq_num_one + 1, seq_num_two + 1, atom_count, sum_sd, false, ignore_virtual_atom );

			}

		}


	}

	Real const rmsd_square = sum_sd/( atom_count );
	Real rmsd = sqrt( rmsd_square );

	if ( atom_count == 0 ) rmsd = 0.0; //special case...implement this on May 5, 2010

	if ( verbose ) {
		std::cout << "sum_sd = " << sum_sd << " atom_count = " << atom_count << " rmsd = " << rmsd << std::endl;
		output_title_text( "Exit In full_length_rmsd_over_residue_list function_general", TR );
	}

	return ( std::max( 0.01, rmsd ) );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Update on Aug 28, 2011 Parin S.
void
align_pdbs(){

	using namespace core::io::silent;
	using namespace core::scoring;
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::pose;
	using namespace protocols::farna;
	using namespace protocols::stepwise::modeler::rna;

	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_RNA );

	if ( !option[ in::file::native ].user() ) utility_exit_with_message( "User must supply in::file::native!" );
	if ( !option[ in::file::s ].user() ) utility_exit_with_message( "User must supply in::file::s!" );
	if ( !option[ alignment_res_pairs ].user() ) utility_exit_with_message( "User must supply alignment_res_pairs!" );
	if ( !option[ alignment_RMSD_CUTOFF ].user() ) utility_exit_with_message( "User must supply alignment_RMSD_CUTOFF!" );

	bool const align_base_only = option[align_only_over_base_atoms ]();

	output_boolean( "align_base_only = ", align_base_only, TR ); std::cout << std::endl;


	std::string const static_pdb_tag = option[ in::file::native ]();

	utility::vector1< std::string > const moving_pdb_tag_list = option[ in::file::s ]();

	utility::vector1< std::string > const alignment_res_string_pair_list = option[ alignment_res_pairs ]();

	utility::vector1< core::Size > const native_virtual_res_list = option[ native_virtual_res ]();

	output_seq_num_list( "native_virtual_res_list = ", native_virtual_res_list, TR );
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////

	if ( alignment_res_string_pair_list.size() == 0 ) utility_exit_with_message( "alignment_res_string_pair_list.size() == 0" );

	utility::vector1< std::pair < Size, Size > > alignment_res_pair_list; //first:static_res to second:moving_res!

	for ( Size n = 1; n <= alignment_res_string_pair_list.size(); n++ ) {

		utility::vector1< std::string > const alignment_res_string_pair = tokenize( alignment_res_string_pair_list[n], "-" );
		if ( alignment_res_string_pair.size() != 2 ) {
			utility_exit_with_message( "alignment_res_string_pair.size() != 2, alignment_res_string_pair_list[n] = " + alignment_res_string_pair_list[n] );
		}

		alignment_res_pair_list.push_back( std::make_pair( string_to_int( alignment_res_string_pair[1] ), string_to_int( alignment_res_string_pair[2] ) ) );

		//PREVIOUSLY:
		//static_pdb_align_res.push_back(string_to_int(alignment_res_pair[1]));
		//moving_pdb_align_res.push_back(string_to_int(alignment_res_pair[2]));

	}

	std::cout << "native_alignment_res to decoy_alignment_res:" << std::endl;
	for ( Size ii = 1; ii <= alignment_res_pair_list.size(); ii++ ) {
		std::cout << alignment_res_pair_list[ii].first << " ---> " << alignment_res_pair_list[ii].second << std::endl;
	}
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	Real const alignment_RMSD_cutoff = option[ alignment_RMSD_CUTOFF ]();

	std::cout << "alignment_RMSD_cutoff = " << alignment_RMSD_cutoff << std::endl;

	SilentFileData silent_file_data;
	std::string silent_file = "aligned_to_" + get_tag_from_pdb_filename( static_pdb_tag ) + ".out";

	pose::Pose static_pose;
	import_pose::pose_from_pdb( static_pose, *rsd_set, static_pdb_tag );

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///By applying the virtual_res, align_pose_general()/setup_suite_atom_id_map() will skip over the virtual_atoms///
	///Only need to apply to static/native pdb pose since only need virtual atoms in one of the two structures being aligned.
	for ( Size virtual_res_ID = 1; virtual_res_ID <= native_virtual_res_list.size(); virtual_res_ID++ ) {
		apply_virtual_rna_residue_variant_type( static_pose, native_virtual_res_list[virtual_res_ID], true /*apply_check*/ );
	}
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	std::string const static_tag_name = get_tag_from_pdb_filename( static_pdb_tag );

	for ( Size n = 1; n <= moving_pdb_tag_list.size(); n++ ) {

		std::string const moving_pdb_tag = moving_pdb_tag_list[n];

		pose::Pose moving_pose;
		import_pose::pose_from_pdb( moving_pose, *rsd_set, moving_pdb_tag );

		align_pose_general( static_pose, static_tag_name, moving_pose, get_tag_from_pdb_filename( moving_pdb_tag ), alignment_res_pair_list, align_base_only );

		std::string output_moving_pdb_tag = "";

		if ( option[ tag_name ].user() ) {
			std::string const user_tag = option[ tag_name ]();
			if ( moving_pdb_tag_list.size() != 1 ) {
				utility_exit_with_message( "User input tag_name but moving_pdb_tag_list.size() != 1" );
			}

			output_moving_pdb_tag = "aligned_" + user_tag + "_to_" + get_tag_from_pdb_filename( static_pdb_tag ) ;
		} else {
			output_moving_pdb_tag = "aligned_" + get_tag_from_pdb_filename( moving_pdb_tag ) + "_to_" + get_tag_from_pdb_filename( static_pdb_tag ) ;
		}

		dump_pdb( moving_pose, output_moving_pdb_tag + ".pdb" );

		Real const align_rmsd = check_alignment_RMSD_cutoff( static_pose, static_tag_name, moving_pose, get_tag_from_pdb_filename( moving_pdb_tag ), alignment_res_pair_list, align_base_only, alignment_RMSD_cutoff );

		BinarySilentStruct s( moving_pose, output_moving_pdb_tag );
		if ( align_base_only ) {
			s.add_energy( "align_base_rmsd", align_rmsd );
		} else {
			s.add_energy( "align_rmsd", align_rmsd );
		}

		silent_file_data.write_silent_struct( s, silent_file, true );

	}

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Aug 28, 2011
void
calculate_pairwise_RMSD(){

	using namespace core::io::silent;
	using namespace core::scoring;
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::pose;
	using namespace protocols::farna;
	using namespace protocols::stepwise::modeler::rna;
	using namespace ObjexxFCL;


	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_RNA );

	if ( !option[ in::file::native ].user() ) utility_exit_with_message( "User must supply in::file::native!" );
	if ( !option[ in::file::s ].user() ) utility_exit_with_message( "User must supply in::file::s!" );
	if ( !option[ alignment_res_pairs ].user() ) utility_exit_with_message( "User must supply alignment_res_pairs!" );
	if ( !option[ rmsd_res_pairs ].user() ) utility_exit_with_message( "User must supply rmsd_res_pairs!" );
	if ( !option[ alignment_RMSD_CUTOFF ].user() ) utility_exit_with_message( "User must supply alignment_RMSD_CUTOFF!" );

	bool const align_base_only = option[align_only_over_base_atoms ]();
	bool const do_dump_pdb = option[ dump ]();
	Real const alignment_RMSD_cutoff = option[ alignment_RMSD_CUTOFF ]();

	output_boolean( "align_base_only = ", align_base_only, TR ); std::cout << std::endl;
	output_boolean( "do_dump_pdb = ", do_dump_pdb, TR ); std::cout << std::endl;
	std::cout << "alignment_RMSD_cutoff = " << alignment_RMSD_cutoff << std::endl;


	std::string const native_silent_file = option[ in::file::native ]();

	std::string const native_tag = option[ native_tag_name ]();

	std::string const decoy_silent_file = option[ in::file::s ]()[1];

	std::string const decoy_tag = option[ decoy_tag_name ]()[1];

	pose::Pose native_pose;
	pose::Pose decoy_pose;

	import_pose_from_silent_file( native_pose, native_silent_file, native_tag );

	import_pose_from_silent_file( decoy_pose, decoy_silent_file, decoy_tag );

	utility::vector1< std::string > const alignment_res_string_pair_list = option[ alignment_res_pairs ]();

	utility::vector1< std::string > const rmsd_res_string_pair_list = option[ rmsd_res_pairs ]();

	//////////////////////////////////////////////////////////////////////////////////////////////////////
	if ( alignment_res_string_pair_list.size() == 0 ) utility_exit_with_message( "alignment_res_string_pair_list.size() == 0" );

	utility::vector1< std::pair < Size, Size > > alignment_res_pair_list; //first:native_res to second:decoy_res!

	for ( Size n = 1; n <= alignment_res_string_pair_list.size(); n++ ) {

		utility::vector1< std::string > const alignment_res_string_pair = tokenize( alignment_res_string_pair_list[n], "-" );
		if ( alignment_res_string_pair.size() != 2 ) {
			utility_exit_with_message( "alignment_res_string_pair.size() != 2, alignment_res_string_pair_list[n] = " + alignment_res_string_pair_list[n] );
		}

		alignment_res_pair_list.push_back( std::make_pair( string_to_int( alignment_res_string_pair[1] ), string_to_int( alignment_res_string_pair[2] ) )  );
	}


	std::cout << "native_alignment_res to decoy_alignment_res:" << std::endl;
	for ( Size ii = 1; ii <= alignment_res_pair_list.size(); ii++ ) {
		std::cout << alignment_res_pair_list[ii].first << " ---> " << alignment_res_pair_list[ii].second << std::endl;
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////
	if ( rmsd_res_string_pair_list.size() == 0 ) utility_exit_with_message( "rmsd_res_string_pair_list.size() == 0" );


	utility::vector1< std::pair < Size, Size > > rmsd_res_pair_list; //first:native_res to second:rmsd_res!

	for ( Size n = 1; n <= rmsd_res_string_pair_list.size(); n++ ) {

		utility::vector1< std::string > const rmsd_res_string_pair = tokenize( rmsd_res_string_pair_list[n], "-" );
		if ( rmsd_res_string_pair.size() != 2 ) {
			utility_exit_with_message( "rmsd_res_string_pair.size() != 2, rmsd_res_string_pair_list[n] = " + rmsd_res_string_pair_list[n] );
		}
		rmsd_res_pair_list.push_back( std::make_pair( string_to_int( rmsd_res_string_pair[1] ), string_to_int( rmsd_res_string_pair[2] ) ) );

	}


	std::cout << "native_calc_rms_res to decoy_calc_rms_res:" << std::endl;
	for ( Size ii = 1; ii <= rmsd_res_pair_list.size(); ii++ ) {
		std::cout << rmsd_res_pair_list[ii].first << " ---> " << rmsd_res_pair_list[ii].second << std::endl;
	}

	for ( Size n = 1; n <= rmsd_res_pair_list.size(); n++ ) {
		Size native_rmsd_res = rmsd_res_pair_list[n].first;
		Size decoy_rmsd_res = rmsd_res_pair_list[n].second;
		output_boolean( "is_native_virtual_res( " + string_of( native_rmsd_res ) + " ) = ",
			native_pose.residue( native_rmsd_res ).has_variant_type( core::chemical::VIRTUAL_RNA_RESIDUE ), TR );
		output_boolean( " | is_decoy_virtual_res( " + string_of( decoy_rmsd_res ) + " ) = ",
			decoy_pose.residue( decoy_rmsd_res ).has_variant_type( core::chemical::VIRTUAL_RNA_RESIDUE ), TR );
		std::cout << std::endl;
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////

	align_pose_general( native_pose, native_tag, decoy_pose, decoy_tag, alignment_res_pair_list, align_base_only );

	check_alignment_RMSD_cutoff( native_pose, native_tag, decoy_pose, decoy_tag, alignment_res_pair_list, align_base_only, alignment_RMSD_cutoff );

	Real const full_rmsd = full_length_rmsd_over_reside_list_general( native_pose, decoy_pose, rmsd_res_pair_list, true /*verbose*/,  false /*ignore_virtual_atom*/ );


	if ( do_dump_pdb ) {
		dump_pdb( native_pose, "ALIGNED_native_" + native_tag + ".pdb" );
		dump_pdb( decoy_pose, "ALIGNED_decoy_" + decoy_tag + ".pdb" );
	}

	std::cout << "RMSD between native_tag( " << native_tag << " ) and decoy_tag( " << decoy_tag << " ) IS " << full_rmsd << std::endl;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
import_and_dump_pdb(){

	using namespace core::chemical;
	using namespace protocols::stepwise::modeler::rna;
	using namespace core::id;

	if ( option[ in::file::s ].user() == false ) utility_exit_with_message( "User must supply in::file::s!" );

	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_RNA );

	utility::vector1< std::string > const pdb_file_list = option[ in::file::s ]();

	for ( Size pdb_ID = 1; pdb_ID <= pdb_file_list.size(); pdb_ID++ ) {

		std::string const pdb_file = pdb_file_list[pdb_ID];

		pose::Pose pose;
		import_pose::pose_from_pdb( pose, *rsd_set, pdb_file );

		std::string const output_pdb_file = "rosetta_" + pdb_file;

		if ( utility::file::file_exists( output_pdb_file ) ) {
			utility_exit_with_message( "output_pdb_file ( " + output_pdb_file + " ) already exist!" );
		}

		dump_pdb( pose, output_pdb_file );
	}

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
o2prime_packer(){

	// using namespace core::io::silent;
	// using namespace core::scoring;
	using namespace core::chemical;
	// using namespace core::conformation;
	// using namespace core::pose;
	// using namespace protocols::farna;
	using namespace protocols::stepwise::modeler::rna;
	using namespace core::id;

	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_RNA );

	// core::scoring::ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( "single_strand_benchmark" );

	core::scoring::ScoreFunctionOP scorefxn = create_scorefxn(); //replace this on Jun 11, 2010

	utility::vector1< std::string > const pdb_tags_from_disk( option[ in::file::s ]() );


	for ( Size n = 1; n <= pdb_tags_from_disk.size(); n++ ) {

		std::string pose_name = pdb_tags_from_disk[n];
		std::cout << pose_name << std::endl;
		pose::Pose pose;
		import_pose::pose_from_pdb( pose, *rsd_set, pose_name );

		if ( option[reset_o2prime_torsion]() ) {
			for ( Size seq_num = 1; seq_num <= pose.total_residue(); seq_num++ ) {
				pose.set_torsion( TorsionID( seq_num, id::CHI, 4 ), 0.0 );
			}
			dump_pdb( pose, "RESETTED_BEFORE_o2prime_pack_" + pose_name );
		}


		o2prime_trials( pose, scorefxn ); //replace this on Jun 11, 2010

		/*
		pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ) );
		task->initialize_from_command_line();

		for ( Size i = 1; i <= pose.total_residue(); i++ ) {
		if ( !pose.residue( i ).is_RNA() ) continue;
		task->nonconst_residue_task( i ).and_extrachi_cutoff( 0 );
		//  task->nonconst_residue_task(i).or_ex4( true );
		task->nonconst_residue_task( i ).or_include_current( true );
		// How about bump check?
		}

		pack::rotamer_trials( pose, *scorefxn, task );
		*/
		dump_pdb( pose, "o2prime_pack_" + pose_name );
	}


}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
mutate_residue( pose::Pose & pose, Size const seq_num, std::string const & res_name )
{
	using namespace core::id;
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::kinematics;
	// using namespace core::optimization;
	// using namespace core::io::silent;
	//
	// using namespace protocols::farna;

	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_RNA );

	pose::Pose start_pose = pose;

	std::cout << "before_mutation: " << std::endl; print_torsion_info( pose, seq_num, 1, "side_chain" );

	core::chemical::AA res_aa = aa_from_name( res_name );
	ResidueOP new_rsd = conformation::ResidueFactory::create_residue( *( rsd_set->get_representative_type_aa( res_aa ) ), start_pose.residue( seq_num ), start_pose.conformation(), true );
	// choose true to preserve the angle connecting bb and sidechain*/

	pose.replace_residue( seq_num, *new_rsd, true /*orient_backbone*/ ); //false doesn't work when input residue is not identical to replace residue???

	pose.set_torsion( TorsionID( seq_num, id::CHI, 1 ), start_pose.residue( seq_num ).chi( 1 ) );

	std::cout << "after_mutation: " << std::endl; print_torsion_info( pose, seq_num, 1, "side_chain" );

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
mutate_residues_wrapper()
{

	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace core::kinematics;
	using namespace core::optimization;
	using namespace core::io::silent;
	using namespace core::id;
	using namespace protocols::farna;

	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_RNA );

	pose::Pose pose;
	std::string pdb_file  = option[ in::file::s ][1];
	import_pose::pose_from_pdb( pose, *rsd_set, pdb_file );
	protocols::farna::make_phosphate_nomenclature_matches_mini( pose );

	// pose::Pose start_pose= pose; //Hard copy
	// dump_pdb( start_pose, "start.pdb");

	std::string rebuild_residue_string = option[rebuild_sequence];
	utility::vector1 < Residue_info > strand_residue_list = Convert_rebuild_residue_string_to_list( rebuild_residue_string );

	for ( Size n = 1; n <= strand_residue_list.size(); n++ ) {
		std::cout << "mutate_residue: "; output_residue_struct( strand_residue_list[n] );
		mutate_residue( pose, strand_residue_list[n].seq_num, strand_residue_list[n].name );
	}


	dump_pdb( pose, "mutate_" + rebuild_residue_string + "_" + pdb_file );

}

//////////////////////////////////////////////////////////////////////////////////////


void
slice_ellipsoid_envelope(){


	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::kinematics;
	using namespace core::optimization;
	using namespace core::io::silent;
	using namespace protocols::stepwise::modeler::rna;
	using namespace core::conformation;
	using namespace ObjexxFCL;
	using namespace core::id;

	clock_t const time_start( clock() );
	/////////////////////////
	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_RNA );

	core::scoring::ScoreFunctionOP scorefxn = create_scorefxn(); //replace this on Jun 11, 2010


	//////////////////import from -s (need to setup pose)/////////////////////////////////////////////

	if ( !option[ in::file::s ].user() ) utility_exit_with_message( "Must supply in::file::s!" );

	if ( !option[ sample_res ].user() ) utility_exit_with_message( "Must supply sample_res!" );

	std::string const pdb_file =   option[ in::file::s ]()[1];

	std::cout << "importing " << pdb_file << std::endl;
	pose::Pose pose;

	import_pose::pose_from_pdb( pose, *rsd_set, pdb_file );


	setup_simple_fold_tree( pose );

	utility::vector1< core::Size > input_sample_res_list = option[ sample_res ]();
	sort_seq_num_list( input_sample_res_list );
	utility::vector1< core::Size > const sample_res_list = input_sample_res_list;

	utility::vector1< core::Size > additional_slice_res_list = option[ additional_slice_res]();

	output_seq_num_list( "additional_slice_res_list = ", additional_slice_res_list, TR );

	for ( Size ii = 1; ii <= additional_slice_res_list.size(); ii++ ) {
		if ( additional_slice_res_list[ii] < 1 )  utility_exit_with_message( "additional_slice_res_list[" + string_of( ii ) + "] < 1" );
		if ( additional_slice_res_list[ii] > pose.total_residue() ) utility_exit_with_message( "additional_slice_res_list[" + string_of( ii ) + "] > pose.total_residue()" );
	}


	Size const five_prime_boundary = sample_res_list[1] - 1;

	Size const three_prime_boundary = sample_res_list[sample_res_list.size()] + 1;

	Size const num_loop_res = ( three_prime_boundary - five_prime_boundary ) - 1;

	std::cout << "five_prime_boundary = " << five_prime_boundary << " three_prime_boundary = " << three_prime_boundary << " num_loop_res = " << num_loop_res << std::endl;

	utility::vector1< core::Size > sample_res_final_seq_num = sample_res_list;

	utility::vector1< core::Size > keep_res_list;

	numeric::xyzVector< Real > const five_prime_foci = pose.residue( five_prime_boundary ).xyz( " O3'" );

	numeric::xyzVector< Real > const three_prime_foci = pose.residue( three_prime_boundary ).xyz( " C5'" );

	Real const expand_radius = option[ surrounding_radius ]();
	Real const int_expand_radius = int( 10*expand_radius );

	Real const foci_sep_dist = ( three_prime_foci - five_prime_foci ).length();

	Real const max_loop_length = ( num_loop_res*O3I_O3I_PLUS_ONE_MAX_DIST ) + O3I_C5I_PLUS_ONE_MAX_DIST;

	Real const major_diameter = max_loop_length + expand_radius;

	std::cout << "foci_sep_dist = " << foci_sep_dist << " major_diameter = " << major_diameter << " max_loop_length = " << max_loop_length << " expand_radius = " << expand_radius << std::endl;

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	std::string pose_name;

	size_t found = pdb_file.rfind( '/' );

	if ( found != std::string::npos ) {
		pose_name = pdb_file.substr( found + 1 );
	} else {
		pose_name = pdb_file;
	}


	//////////////Create a rotated pose which have the major axis as the z-axis and the center of the ellipsoid as the origin!/////////

	pose::Pose rotated_pose = pose;

	numeric::xyzVector< Real > const new_origin = ( five_prime_foci + three_prime_foci )/2;

	numeric::xyzVector< Real > internal_z_axis = numeric::xyzVector< Real > ( 0.0, 0.0, 1.0 );

	numeric::xyzVector< Real > z_axis = ( three_prime_foci - five_prime_foci );
	z_axis.normalize();

	numeric::xyzVector< Real > y_axis = cross( internal_z_axis, z_axis );
	y_axis.normalize();

	numeric::xyzVector< Real > x_axis = cross( y_axis, z_axis );
	x_axis.normalize();


	numeric::xyzMatrix< core::Real > const & coordinate_matrix = numeric::xyzMatrix< core::Real > ::cols( x_axis, y_axis, z_axis );

	numeric::xyzMatrix< core::Real > rotation_matrix = inverse( coordinate_matrix );

	for ( Size i = 1; i <= rotated_pose.total_residue(); ++i ) {
		for ( Size j = 1; j <= rotated_pose.residue_type( i ).natoms(); ++j ) { // use residue_type to prevent internal coord update

			id::AtomID const id( j, i );

			rotated_pose.set_xyz( id, rotated_pose.xyz( id ) - new_origin );
			rotated_pose.set_xyz( id, rotation_matrix * rotated_pose.xyz( id ) );

		}
	}

	dump_pdb( rotated_pose, "rotated_" + pose_name );

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	for ( Size seq_num = 1; seq_num <= pose.total_residue(); seq_num++ ) {

		if ( sample_res_list.has_value( seq_num ) ) {
			std::cout << "res " << seq_num << " is a sample_res" << std::endl;
			keep_res_list.push_back( seq_num );
			continue;
		}

		if ( additional_slice_res_list.has_value( seq_num ) ) {
			std::cout << "res " << seq_num << " is in additional_slice_res_list" << std::endl;
			keep_res_list.push_back( seq_num );
			continue;
		}


		core::conformation::Residue const & surrounding_rsd = pose.residue( seq_num );

		core::conformation::Residue const & rotated_surr_rsd = rotated_pose.residue( seq_num );

		if ( surrounding_rsd.natoms() != rotated_surr_rsd.natoms() ) {
			utility_exit_with_message( "surrounding_rsd.natoms() !=\trotated_surr_rsd.natoms()" );
		}

		for ( Size atomno = 1; atomno <= surrounding_rsd.natoms(); atomno++ ) {

			/////////////////////////////////////////////////////////////////////////////

			numeric::xyzVector< Real > const atom_xyz = surrounding_rsd.xyz( atomno );

			Real const sum_distances_to_foci = ( atom_xyz - five_prime_foci ).length() + ( atom_xyz - three_prime_foci ).length();

			bool const method_1 = ( sum_distances_to_foci <= major_diameter ) ? true : false;

			//////////////////////////////////////////////////////////////////////////////
			numeric::xyzVector< Real > const rotated_atom_xyz = rotated_surr_rsd.xyz( atomno );

			Real const r_sq = rotated_atom_xyz.x()*rotated_atom_xyz.x() + rotated_atom_xyz.y()*rotated_atom_xyz.y();
			Real const z_sq = rotated_atom_xyz.z()*rotated_atom_xyz.z();

			Real const denominator_one = ( ( major_diameter*major_diameter ) - ( foci_sep_dist*foci_sep_dist ) )/4.0;

			Real const denominator_two = ( major_diameter*major_diameter )/4.0;

			Real const ellipse_equation = ( r_sq/denominator_one ) + ( z_sq/denominator_two );

			bool const method_2 = ( ellipse_equation <= 1.0 ) ? true: false;

			//////////////////////////////////////////////////////////////////////////////


			if ( method_1 != method_2 ) {
				std::cout << "method_1 != method_2" << std::endl;
				output_boolean( "method_1 = ", method_1, TR ); std::cout << std::endl;
				output_boolean( "method_2 = ", method_2, TR ); std::cout << std::endl;
				utility_exit_with_message( "method_1 != method_2" );
			}


			if ( method_1 ) {
				if ( keep_res_list.has_value( seq_num ) == false ) { //Not already in the list!
					std::cout << "seq_num ( " << seq_num << " ) is a surrounding_res, sum_distances_to_foci = " << sum_distances_to_foci << " major_diameter = " << major_diameter << std::endl;
					keep_res_list.push_back( seq_num );
					//break;
				}
			}
		}
	}

	pose::Pose output_pose = pose;

	for ( Size seq_num = pose.total_residue(); seq_num >= 1; seq_num-- ) {

		if ( keep_res_list.has_value( seq_num ) == false ) {

			for ( Size n = 1; n <= sample_res_final_seq_num.size(); n++ ) {
				if ( sample_res_list[n] > seq_num ) {
					sample_res_final_seq_num[n]--;
				}
			}

			output_pose.conformation().delete_residue_slow( seq_num );

		}
	}


	for ( Size n = 1; n <= sample_res_final_seq_num.size(); n++ ) {
		std::cout << sample_res_list[n] << " --> " << sample_res_final_seq_num[n] << std::endl;
	}

	std::cout << "pose_name = " << pose_name << std::endl;
	dump_pdb( pose, "input_" + pose_name );


	pose::Pose no_loop_output_pose = output_pose; //copy before perform o2prime minimize with loop as part of struct.

	o2prime_trials( output_pose, scorefxn );
	dump_pdb( output_pose, "ellipsoid_expand_radius_" + string_of( int_expand_radius ) + "_" + pose_name );


	//make sure that sample_res_list is sorted
	utility::vector1< core::Size > sorted_sample_res_final_seq_num_list = sample_res_final_seq_num;
	sort_seq_num_list( sorted_sample_res_final_seq_num_list );


	for ( Size ii = sorted_sample_res_final_seq_num_list.size(); ii >= 1; ii-- ) {
		Size const seq_num = sorted_sample_res_final_seq_num_list[ii];
		no_loop_output_pose.conformation().delete_residue_slow( seq_num );
	}


	//Reset o2prime torsion....Dec 7, 2010///////////////////////////////////
	for ( Size seq_num = 1; seq_num <= no_loop_output_pose.total_residue(); seq_num++ ) {
		no_loop_output_pose.set_torsion( TorsionID( seq_num, id::CHI, 4 ), 0.0 );
	}
	dump_pdb( no_loop_output_pose, "RESETTED_o2prime_no_loop_ellipsoid_expand_radius_" + string_of( int_expand_radius )  + "_" + pose_name );
	///////////////////////////////////////////////////////////////////////


	o2prime_trials( no_loop_output_pose, scorefxn );
	dump_pdb( no_loop_output_pose, "no_loop_ellipsoid_expand_radius_" + string_of( int_expand_radius )  + "_" + pose_name );

	std::cout << "Sliced out " << output_pose.total_residue() << " out of " << pose.total_residue()  << " nucleotides" << std::endl;

	std::cout << "Total time in slice_ellipsoid_envelope: " << static_cast< Real > ( clock() - time_start ) / CLOCKS_PER_SEC << std::endl;


}
//////////////////////////////////////////////////////////////////////////////////////

void
slice_sample_res_and_surrounding(){

	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::kinematics;
	using namespace core::optimization;
	using namespace core::io::silent;
	using namespace protocols::stepwise::modeler::rna;
	using namespace core::conformation;
	using namespace ObjexxFCL;
	using namespace core::id;

	clock_t const time_start( clock() );
	/////////////////////////
	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_RNA );

	core::scoring::ScoreFunctionOP scorefxn = create_scorefxn(); //replace this on Jun 11, 2010


	//////////////////import from -s (need to setup pose)/////////////////////////////////////////////

	if ( !option[ in::file::s ].user() ) utility_exit_with_message( "Must supply in::file::s!" );

	if ( !option[ sample_res ].user() ) utility_exit_with_message( "Must supply sample_res!" );

	std::string const pdb_file =   option[ in::file::s ]()[1];

	std::cout << "importing " << pdb_file << std::endl;
	pose::Pose pose;

	import_pose::pose_from_pdb( pose, *rsd_set, pdb_file );


	setup_simple_fold_tree( pose );

	utility::vector1< core::Size > input_sample_res_list = option[ sample_res ]();
	sort_seq_num_list( input_sample_res_list );
	utility::vector1< core::Size > const sample_res_list = input_sample_res_list;

	// get_surrounding_O2prime_hydrogen(pose, sample_res_list, true);


	utility::vector1< core::Size > sample_res_final_seq_num = sample_res_list;

	utility::vector1< core::Size > keep_res_list;

	Real const expand_radius = option[ surrounding_radius ]();
	Real const int_expand_radius = int( 10*expand_radius );

	std::cout << "expand_radius = " << expand_radius << std::endl;

	utility::vector1< core::Size > additional_slice_res_list = option[ additional_slice_res]();

	output_seq_num_list( "additional_slice_res_list = ", additional_slice_res_list, TR );


	for ( Size ii = 1; ii <= additional_slice_res_list.size(); ii++ ) {
		if ( additional_slice_res_list[ii] < 1 )  utility_exit_with_message( "additional_slice_res_list[" + string_of( ii ) + "] < 1" );
		if ( additional_slice_res_list[ii] > pose.total_residue() ) utility_exit_with_message( "additional_slice_res_list[" + string_of( ii ) + "] > pose.total_residue()" );
	}


	pose::Pose output_pose = pose;

	for ( Size seq_num = pose.total_residue(); seq_num >= 1; seq_num-- ) {

		if ( sample_res_list.has_value( seq_num ) ) {
			std::cout << "res " << seq_num << " is a sample_res" << std::endl;
			keep_res_list.push_back( seq_num );
			continue;
		}

		if ( additional_slice_res_list.has_value( seq_num ) ) {
			std::cout << "res " << seq_num << " is in additional_slice_res_list" << std::endl;
			keep_res_list.push_back( seq_num );
			continue;
		}

		core::conformation::Residue const & surrounding_rsd = pose.residue( seq_num );

		bool is_surrounding_res = false;

		for ( Size ii = 1; ii <= sample_res_list.size(); ii++ ) {
			if ( is_surrounding_res == true ) break;

			bool is_very_far_from_each_other = false;

			Size const sample_res = sample_res_list[ii];

			core::conformation::Residue const & sample_rsd = pose.residue( sample_res );

			for ( Size surr_at = 1; surr_at <= surrounding_rsd.natoms(); surr_at++ ) {
				if ( is_surrounding_res == true ) break;
				if ( is_very_far_from_each_other == true ) break;

				for ( Size sample_at = 1; sample_at <= sample_rsd.natoms(); sample_at++ ) {
					if ( is_surrounding_res == true ) break;
					if ( is_very_far_from_each_other == true ) break;

					//     if( (surrounding_rsd.xyz(surr_at)-sample_rsd.xyz(sample_at) ).length_squared() > 35*35) {
					//      std::cout << "res " << seq_num << " is very_far_away from sample_res " << sample_res << ", length()= " << ( surrounding_rsd.xyz(surr_at)-sample_rsd.xyz(sample_at) ).length()  << std::endl;
					//      is_very_far_from_each_other=true;
					//      break;
					//     }


					if ( ( surrounding_rsd.xyz( surr_at ) - sample_rsd.xyz( sample_at ) ).length_squared() < expand_radius*expand_radius ) {
						std::cout << "res " << seq_num << " is a surrounding res, length() = " << ( surrounding_rsd.xyz( surr_at ) - sample_rsd.xyz( sample_at ) ).length()  << std::endl;
						keep_res_list.push_back( seq_num );
						is_surrounding_res = true;
						break;
					}
				}
			}
		}

		if ( is_surrounding_res == false ) {

			for ( Size n = 1; n <= sample_res_final_seq_num.size(); n++ ) {
				if ( sample_res_list[n] > seq_num ) {
					sample_res_final_seq_num[n]--;
				}
			}

			output_pose.conformation().delete_residue_slow( seq_num );
		}
	}


	std::string pose_name;

	size_t found = pdb_file.rfind( '/' );

	if ( found != std::string::npos ) {
		pose_name = pdb_file.substr( found + 1 );
	} else {
		pose_name = pdb_file;
	}


	for ( Size n = 1; n <= sample_res_final_seq_num.size(); n++ ) {
		std::cout << sample_res_list[n] << " --> " << sample_res_final_seq_num[n] << std::endl;
	}

	std::cout << "pose_name = " << pose_name << std::endl;
	dump_pdb( pose, "input_" + pose_name );


	pose::Pose no_loop_output_pose = output_pose; //copy before perform o2prime minimize with loop as part of struct.

	o2prime_trials( output_pose, scorefxn );
	dump_pdb( output_pose, "expand_radius_" + string_of( int_expand_radius ) + "_" + pose_name );


	//make sure that sample_res_list is sorted
	utility::vector1< core::Size > sorted_sample_res_final_seq_num_list = sample_res_final_seq_num;
	sort_seq_num_list( sorted_sample_res_final_seq_num_list );


	for ( Size ii = sorted_sample_res_final_seq_num_list.size(); ii >= 1; ii-- ) {
		Size const seq_num = sorted_sample_res_final_seq_num_list[ii];
		no_loop_output_pose.conformation().delete_residue_slow( seq_num );
	}


	//Reset o2prime torsion....Dec 7, 2010///////////////////////////////////
	for ( Size seq_num = 1; seq_num <= no_loop_output_pose.total_residue(); seq_num++ ) {
		no_loop_output_pose.set_torsion( TorsionID( seq_num, id::CHI, 4 ), 0.0 );
	}
	dump_pdb( no_loop_output_pose, "RESETTED_o2prime_no_loop_expand_radius_" + string_of( int_expand_radius )  + "_" + pose_name );
	///////////////////////////////////////////////////////////////////////


	o2prime_trials( no_loop_output_pose, scorefxn );
	dump_pdb( no_loop_output_pose, "no_loop_expand_radius_" + string_of( int_expand_radius )  + "_" + pose_name );


	std::cout << "Total time in slice_sample_res_and_surrounding: " << static_cast< Real > ( clock() - time_start ) / CLOCKS_PER_SEC << std::endl;

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
pdb_to_silent_file(){

	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::kinematics;
	using namespace core::optimization;
	using namespace core::io::silent;
	using namespace protocols::stepwise::modeler::rna;
	using namespace core::conformation;
	using namespace ObjexxFCL;
	using namespace core::io::silent;
	using namespace core::scoring;
	using namespace core::pose;

	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_RNA );
	SilentFileData silent_file_data;

	pose::Pose viewer_pose;

	protocols::viewer::add_conformation_viewer( viewer_pose.conformation(), "test", 400, 400 );


	if ( !option[ in::file::s ].user() ) utility_exit_with_message( "User must supply in::file::s!" );

	if ( !option[ output_silent_file].user() ) utility_exit_with_message( "User must supply output_silent_file!" );

	std::string const silent_outfile = option[ output_silent_file]();

	utility::vector1< std::string > const pdb_file_list = option[ in::file::s ]();

	utility::vector1< Size > const virtual_res_list = option[ virtual_res]();

	utility::vector1< Size > const virtual_sugar_list = option[ virtual_sugar]();

	utility::vector1< Real > const list_of_pose_energy = option[ list_of_energy ]();

	utility::vector1< std::string > const list_of_pose_virtual_res = option[ list_of_virtual_res]();

	if ( virtual_res_list.size() > 0 && list_of_pose_virtual_res.size() > 0 ) utility_exit_with_message( "virtual_pose_res_list.size() > 0 && list_of_virtual_res > 0" );

	if ( list_of_pose_virtual_res.size() > 0 ) {
		if ( list_of_pose_virtual_res.size() != pdb_file_list.size() ) utility_exit_with_message( "list_of_pose_virtual_res.size() != pdb_file_list.size()" );
	}

	for ( Size n = 1; n <= pdb_file_list.size(); n++ ) {

		std::string const pdb_file = pdb_file_list[n];

		std::string tag = get_tag_from_pdb_filename( pdb_file );

		if ( option[ tag_name ].user() ) {
			if ( pdb_file_list.size() != 1 ) utility_exit_with_message( "User passed in a single tag_name but pdb_file_list.size() != 1 " );
			tag = option[ tag_name ]();
		}

		std::cout << "importing pdb_file: " << pdb_file << std::endl;


		pose::Pose pose;

		import_pose::pose_from_pdb( pose, *rsd_set, pdb_file );

		protocols::farna::make_phosphate_nomenclature_matches_mini( pose );

		utility::vector1< Size > act_virtual_res_list;

		if ( virtual_res_list.size() > 0 ) {
			act_virtual_res_list = virtual_res_list;
		} else if ( option[ list_of_virtual_res].user() ) {
			utility::vector1< std::string > virtual_res_string_list = tokenize( list_of_pose_virtual_res[n], "-" );

			for ( Size ii = 1; ii <= virtual_res_string_list.size(); ii++ ) {
				Size const virtual_seq_num = string_to_int( virtual_res_string_list[ii] );
				if ( virtual_seq_num == 0 ) continue;
				act_virtual_res_list.push_back( virtual_seq_num );
			}
		}


		for ( Size ii = 1; ii <= act_virtual_res_list.size(); ii++ ) {
			apply_virtual_rna_residue_variant_type( pose, act_virtual_res_list[ii], false /*apply_check*/ ) ;
		}

		for ( Size ii = 1; ii <= virtual_sugar_list.size(); ii++ ) {
			add_variant_type_to_pose_residue( pose, core::chemical::VIRTUAL_RIBOSE, virtual_sugar_list[ii] );
		}

		viewer_pose = pose;

		std::cout << "converting pdb_file: " << pdb_file << " to silent_struct " << tag << std::endl;

		BinarySilentStruct s( pose, tag );

		if ( option[ list_of_energy ].user() ) {
			s.add_energy( "score", list_of_pose_energy[n] );
		}

		silent_file_data.write_silent_struct( s, silent_outfile, false );

	}

}


void
rna_fullatom_minimize_test()
{

	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::kinematics;
	using namespace core::optimization;
	using namespace core::io::silent;
	using namespace protocols::stepwise::modeler::rna;

	/////////////////////////
	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_RNA );


	if ( option[ in::file::silent ].user() &&  option[ input_tag_list ].user() ) {
		utility_exit_with_message( "user specify both silent_file and input_tag as pose source! ONLY one pose source allow!" );
	}

	//////////////////import from -s (need to setup pose)/////////////////////////////////////////////
	pose::Pose pose;

	std::string output_pose_name;

	if ( option[ input_tag_list ].user() ) {

		std::string const pdb_file =  option[ input_tag_list ]()[1];

		std::cout << "importing " << pdb_file << std::endl;

		output_pose_name =  "minimize_" + path_basename( pdb_file );

		import_pose::pose_from_pdb( pose, *rsd_set, pdb_file );


	} else if ( option[ in::file::silent ].user() ) {

		core::import_pose::pose_stream::SilentFilePoseInputStreamOP input( new core::import_pose::pose_stream::SilentFilePoseInputStream() );
		input->set_order_by_energy( true );

		utility::vector1< std::string > silent_files = option[in::file::silent ]();

		input->filenames( silent_files ); //triggers read in of files, too.

		bool found_silent_struct = false;
		if ( input->has_another_pose() ) { //Just get the first structure
			found_silent_struct = true;

			core::io::silent::SilentStructOP silent_struct( input->next_struct() );
			silent_struct->fill_pose( pose );
			std::cout << "importing silent_ = " << silent_struct->decoy_tag() << std::endl;

			output_pose_name = "minimize_silent_" + silent_struct->decoy_tag() + ".pdb";

		}

		if ( found_silent_struct == false ) utility_exit_with_message( "found_silent_struct == false!" );

	} else {

		utility_exit_with_message( "user need to specify input pose source!" );
	}

	////////////////////////

	if ( option[ graphic ]() ) protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 400, 400 );

	protocols::farna::RNA_Minimizer rna_minimizer;

	rna_minimizer.deriv_check( option[ minimizer_deriv_check ]() );

	rna_minimizer.set_include_default_linear_chainbreak( false );

	rna_minimizer.use_coordinate_constraints( false );
	rna_minimizer.set_verbose( true );
	rna_minimizer.vary_bond_geometry( false );
	rna_minimizer.skip_o2prime_trials( option[ minimizer_skip_o2prime_trials ] );
	rna_minimizer.set_perform_minimizer_run( option[ minimizer_perform_minimizer_run ] );

	rna_minimizer.set_do_dump_pdb( true );

	if ( option[ minimizer_min_type ].user() ) {
		rna_minimizer.set_min_type( option[ minimizer_min_type ]() );
	}

	rna_minimizer.apply( pose );


	dump_pdb( pose, output_pose_name );

	/////////////////////////////////////////////
	std::string silent_file = "output_silent.out";

	SilentFileData silent_file_data;

	BinarySilentStruct s( pose, output_pose_name );

	std::cout << "Outputting " << output_pose_name << " to silent file: " << silent_file << std::endl;
	silent_file_data.write_silent_struct( s, silent_file, false /*write score only*/ );


}


///////////////////////////////////////////////////////////////
void*
my_main( void* )
{

	std::string algorithm_input = option[algorithm];

	// if(option[ parin_favorite_output ]()){
	//  system(std::string("mkdir pose/").c_str()); //Output all the poses generated by the code in here. Parin S Jan 28, 2010
	// }

	if ( algorithm_input == "align_pdbs" ) {
		align_pdbs();
	} else if ( algorithm_input == "calculate_pairwise_RMSD" ) {
		calculate_pairwise_RMSD();
	} else if ( algorithm_input == "o2prime_packer" ) {
		o2prime_packer();
	} else if ( algorithm_input == "import_and_dump_pdb" ) {
		import_and_dump_pdb();
	} else if ( algorithm_input == "mutate_residues" ) {
		mutate_residues_wrapper();
	} else if ( algorithm_input == "slice_ellipsoid_envelope" ) {
		slice_ellipsoid_envelope();
	} else if ( algorithm_input == "slice_sample_res_and_surrounding" ) {
		slice_sample_res_and_surrounding();
	} else if ( algorithm_input == "pdb_to_silent_file" ) {
		pdb_to_silent_file();
	} else if ( algorithm_input == "rna_fullatom_minimize_test" ) {
		rna_fullatom_minimize_test();
	} else {
		std::cout << "Error no algorithm selected" << std::endl;
	}

	protocols::viewer::clear_conformation_viewers();
	exit( 0 );

}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {
		utility::vector1< Size > blank_size_vector;
		utility::vector1< std::string > blank_string_vector;
		utility::vector1< Real > blank_real_vector;

		NEW_OPT( graphic, "Turn graphic on/off", true );
		NEW_OPT( algorithm, "Specify algorithm to execute", "" );
		NEW_OPT( surrounding_radius, "expand_radius for slice_sample_res_and_surrounding function", 10.0 );
		NEW_OPT( sample_res, "sample_res", blank_size_vector );
		NEW_OPT( rebuild_sequence, "Specify the rebuild nucleotides in the order you want to rebuild the", "" );
		NEW_OPT( reset_o2prime_torsion, "use in the o2prime_packer function ", false );
		NEW_OPT( rmsd_res_pairs, "rmsd_res_pairs", blank_string_vector );  //1-3 4-5,res 1 of static to res 3 of moving...res 4 of static to res 5 of moving.
		NEW_OPT( alignment_res_pairs, "alignment_res_pairs", blank_string_vector );  //1-3 4-5,res 1 of static to res 3 of moving...res 4 of static to res 5 of moving.
		NEW_OPT( alignment_RMSD_CUTOFF, "alignment_RMSD_CUTOFF", 0.0 );
		NEW_OPT( tag_name, "tag_name", "BLAH_TAG" );
		NEW_OPT( output_silent_file, "output_silent_file", "" );
		NEW_OPT( list_of_virtual_res, " list of virtual_res of each corresponding imported pdb", blank_string_vector );
		NEW_OPT( list_of_energy, " list of energy of each corresponding imported pdb", blank_real_vector );
		NEW_OPT( virtual_res, " virtual_res ", blank_size_vector );
		NEW_OPT( virtual_sugar, " virtual_sugar ", blank_size_vector );
		NEW_OPT( align_only_over_base_atoms, "align_only_over_base_atoms", true );
		NEW_OPT( additional_slice_res, "additional_slice_res", blank_size_vector );
		NEW_OPT( native_virtual_res, " native_virtual_res ( use in align_pdbs() function )", blank_size_vector );
		NEW_OPT( native_tag_name, "native tag from a silent_file", "" );
		NEW_OPT( decoy_tag_name, "decoy tag from a silent_file", blank_string_vector );
		NEW_OPT( dump, "dump pdb", false );
		NEW_OPT( input_tag_list, "input_tag_list", blank_string_vector );
		NEW_OPT( minimizer_deriv_check, "deriv_check", true );
		NEW_OPT( minimizer_min_type, "minimizer_min_type", "" );
		NEW_OPT( minimizer_skip_o2prime_trials, "minimizer_skip_o2prime_trials", true );        //Parin Jan 08, 2012 (Avoid randomness)
		NEW_OPT( minimizer_perform_minimizer_run, "minimizer_perform_minimizer_run", true );  //Parin Jan 20, 2012 (for testing purposes)


		////////////////////////////////////////////////////////////////////////////
		// setup
		////////////////////////////////////////////////////////////////////////////
		core::init::init( argc, argv );

		option[ OptionKeys::chemical::include_patches ].push_back( "VIRTUAL_RIBOSE" );
		option[ OptionKeys::chemical::patch_selectors ].push_back( "TERMINAL_PHOSPHATE" ); // 5prime_phosphate and 3prime_phosphate
		option[ OptionKeys::chemical::include_patches ].push_back( "patches/nucleic/rna/Virtual_RNA_Residue.txt" );
		option[ OptionKeys::chemical::include_patches ].push_back( "patches/nucleic/rna/Virtual_Phosphate.txt" );
		option[ OptionKeys::chemical::include_patches ].push_back( "patches/nucleic/rna/Protonated_H1_Adenosine.txt" );
		option[ OptionKeys::chemical::include_patches ].push_back( "patches/nucleic/rna/Virtual_Backbone_Except_C1prime.txt" );

		////////////////////////////////////////////////////////////////////////////
		// end of setup
		////////////////////////////////////////////////////////////////////////////

		protocols::viewer::viewer_main( my_main );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}


