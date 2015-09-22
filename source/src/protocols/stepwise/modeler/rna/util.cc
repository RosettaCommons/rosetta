// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/stepwise/modeler/rna/util.cc
/// @brief  Util functions for Stepwise Assembly RNA.
/// @author Parin Sripakdeevong


//////////////////////////////////
#include <protocols/stepwise/modeler/rna/util.hh>
#include <protocols/stepwise/modeler/rna/StepWiseRNA_Classes.hh>
#include <protocols/stepwise/modeler/rna/StepWiseRNA_ResidueInfo.hh>
#include <protocols/stepwise/modeler/rna/phosphate/util.hh>
#include <protocols/stepwise/modeler/rna/sugar/util.hh>
#include <protocols/stepwise/modeler/rna/bulge/util.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <protocols/stepwise/modeler/output_util.hh>
#include <core/pose/rna/RNA_BasePairClassifier.hh>
#include <protocols/farna/util.hh>
#include <protocols/toolbox/rigid_body/util.hh>
#include <core/pose/rna/RNA_BaseDoubletClasses.hh>
#include <core/scoring/ScoreType.hh> //Parin Sept 20, 2011.

#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ResidueProperties.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/util.hh>
#include <core/pose/util.hh>
#include <core/pose/util.tmpl.hh>
#include <core/pose/rna/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/rotamer_trials.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/chemical/rna/util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/func/FadeFunc.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <protocols/toolbox/AllowInsert.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/TorsionID.hh>
#include <numeric/conversions.hh>
#include <basic/Tracer.hh>
#include <numeric/random/random.hh>
#include <utility/tools/make_vector1.hh>

#include <iostream>
#include <fstream>
#include <sstream>
#include <ObjexxFCL/format.hh>
#include <set>
#include <time.h>
#include <map>

#include <stdio.h> //Sept 26, 2011

//for process_mem_usage:
#include <ios>

using namespace core;
using namespace core::chemical::rna;


static THREAD_LOCAL basic::Tracer TR( "protocols.stepwise.modeler.rna.util" );

namespace protocols {
namespace stepwise {
namespace modeler {
namespace rna {

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
minimize_with_constraints( core::pose::Pose & pose, core::kinematics::MoveMap const & mm, core::scoring::ScoreFunctionOP const & scorefxn, core::optimization::MinimizerOptions const & options ){

	using namespace core::scoring;
	using namespace core::optimization;

	AtomTreeMinimizer minimizer;
	scoring::constraints::ConstraintSetOP save_pose_constraints = pose.constraint_set()->clone();
	core::scoring::constraints::add_coordinate_constraints( pose );
	minimizer.run( pose, mm, *( scorefxn ), options );
	pose.constraint_set( save_pose_constraints );

}

////////////////////////////////////////////////May 04, 2011////////////////////////////////////////////////////
void
apply_protonated_H1_adenosine_variant_type( core::pose::Pose & pose, core::Size const & seq_num, bool const apply_check ){

	bool verbose = true;

	if ( verbose ) TR << "Applying PROTONATED_H1_ADENOSINE variant_type to seq_num " << seq_num << std::endl;

	if ( apply_check ) {
		//Basically the two variant type are not compatible, VIRTUAL_RNA_RESIDUE variant type currently does not virtualize the protonated H1 atom.
		if ( pose.residue( seq_num ).has_variant_type( core::chemical::VIRTUAL_RNA_RESIDUE ) ) {
			utility_exit_with_message( "Cannot apply PROTONATED_H1_ADENOSINE variant_type to seq_num: " +
				ObjexxFCL::string_of( seq_num ) + ". This residue have a incompatible VIRTUAL_RNA_RESIDUE variant type."  );
		}
	}

	if ( pose.residue( seq_num ).has_variant_type( core::chemical::PROTONATED_H1_ADENOSINE ) ) {
		TR << "WARNING pose already have PROTONATED_H1_ADENOSINE variant_type at seq_num = " << seq_num <<
			", early RETURN!" << std::endl;
		return;
		//utility_exit_with_message("pose already have PROTONATED_H1_ADENOSINE variant_type at seq_num= " + ObjexxFCL::string_of(seq_num));
	}

	if ( pose.total_residue() < seq_num ) {
		utility_exit_with_message(  "Cannot apply PROTONATED_H1_ADENOSINE variant_type to seq_num: " + ObjexxFCL::string_of( seq_num ) + ". pose.total_residue() < seq_num"  );
	}

	if ( pose.residue( seq_num ).aa() != core::chemical::na_rad ) {
		utility_exit_with_message( "working_seq_num = " + ObjexxFCL::string_of( seq_num ) + " cannot have PROTONATED_H1_ADENOSINE variant type since it is not a adenosine!" );
	}

	pose::add_variant_type_to_pose_residue( pose, core::chemical::PROTONATED_H1_ADENOSINE, seq_num );

}

//////////////////////////////////////////////////////////////////////////////////////
void
remove_all_variant_types( pose::Pose & pose ){

	for ( Size seq_num = 1; seq_num <= pose.total_residue(); seq_num++ ) {

		if ( pose.residue( seq_num ).aa() == core::chemical::aa_vrt ) continue; //Fang's electron density code
		utility::vector1< std::string > target_variants( pose.residue( seq_num ).type().properties().get_list_of_variants() );
		runtime_assert ( target_variants.size() == pose.residue( seq_num ).type().properties().get_list_of_variants().size() );
		Size skip_variant_count = 0;
		for ( Size i = 1; i <= target_variants.size(); i++ ) {
			runtime_assert ( pose.residue( seq_num ).type().has_variant_type( target_variants[i] ) );
			bool skip_this_variant = false;
			if ( target_variants[i] == "LOWER_TERMINUS_VARIANT" ) skip_this_variant = true;
			if ( target_variants[i] == "UPPER_TERMINUS_VARIANT" ) skip_this_variant = true;
			if ( skip_this_variant ) {
				skip_variant_count++;
				continue;
			}
			pose::remove_variant_type_from_pose_residue( pose,
				core::chemical::ResidueProperties::get_variant_from_string( target_variants[i] ), seq_num );
		}
		runtime_assert( pose.residue( seq_num ).type().properties().get_list_of_variants().size() == skip_variant_count );
	}

}

//////////////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
apply_full_to_sub_mapping( utility::vector1< Size > const & res_vector, utility::vector1< core::Size > const & is_working_res, std::map< core::Size, core::Size > const & full_to_sub ){

	using namespace ObjexxFCL;
	runtime_assert( is_working_res.size() > 0 );
	runtime_assert ( !full_to_sub.empty() );

	Size const total_res = is_working_res.size();
	utility::vector1< core::Size > working_res_vector;
	for ( Size n = 1; n <= res_vector.size(); n++ ) {
		runtime_assert ( res_vector[ n ] <= total_res );
		if ( !is_working_res[ res_vector[ n ] ] ) continue;
		working_res_vector.push_back( full_to_sub.find( res_vector[ n ] )->second );
	}

	return working_res_vector;

}

//////////////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
apply_full_to_sub_mapping( utility::vector1< Size > const & res_vector, working_parameters::StepWiseWorkingParametersCOP working_parameters ){

	utility::vector1< core::Size > const & is_working_res = working_parameters->is_working_res();
	std::map< core::Size, core::Size > const & full_to_sub = working_parameters->const_full_to_sub();

	return apply_full_to_sub_mapping( res_vector, is_working_res, full_to_sub );

}


//////////////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
apply_sub_to_full_mapping( utility::vector1< Size > const & working_res_vector, working_parameters::StepWiseWorkingParametersCOP working_parameters ){
	std::map< core::Size, core::Size > const & sub_to_full( working_parameters->const_sub_to_full() );
	utility::vector1< core::Size > full_res_vector;
	for ( Size n = 1; n <= working_res_vector.size(); n++ ) {
		runtime_assert( sub_to_full.find( working_res_vector[ n ] ) != sub_to_full.end() );
		full_res_vector.push_back( sub_to_full.find( working_res_vector[ n ] )->second );
	}
	return full_res_vector;
}


///////////////////////////This should be a function of the working_parameters class///////////////////////
void
ensure_valid_full_seq_num( Size const full_seq_num, working_parameters::StepWiseWorkingParametersCOP const & working_parameters ){
	using namespace ObjexxFCL;
	utility::vector1< core::Size > const & is_working_res = working_parameters->is_working_res();
	runtime_assert(  full_seq_num >= 1 );
	runtime_assert(  full_seq_num <= is_working_res.size() );
}

//////////////////////////This should be a function of the working_parameters class/////////////////////////
bool
check_is_working_res( Size const full_seq_num, working_parameters::StepWiseWorkingParametersCOP const & working_parameters ){
	using namespace ObjexxFCL;
	utility::vector1< core::Size > const & is_working_res = working_parameters->is_working_res();
	ensure_valid_full_seq_num( full_seq_num, working_parameters );
	return is_working_res[full_seq_num];
}

//////////////////////////This should be a function of the working_parameters class/////////////////////////
core::Size
check_validity_and_get_working_res( Size const full_seq_num, working_parameters::StepWiseWorkingParametersCOP const & working_parameters ){

	using namespace ObjexxFCL;

	std::map< core::Size, core::Size > const & full_to_sub = working_parameters->const_full_to_sub();

	std::string const & working_sequence = working_parameters->working_sequence();

	if ( check_is_working_res( full_seq_num, working_parameters ) == false ) {
		utility_exit_with_message( "full_seq_num ( " + string_of( full_seq_num ) + " ) is NOT working_res" );
	}

	Size const working_seq_num = full_to_sub.find( full_seq_num )->second;

	if ( working_seq_num < 1 ) utility_exit_with_message( "working_seq_num ( " + string_of( working_seq_num ) + " ) is lesser then 1" );

	if ( working_seq_num > working_sequence.size() ) utility_exit_with_message( "working_seq_num ( " + string_of( working_seq_num ) + " ) is greater than working_sequence.size() ( " + string_of( working_sequence.size() ) + " )" );

	return working_seq_num;

}

//////////////////////////////////////////////////////////////////////////////////////
std::map< core::Size, core::Size >
create_full_to_input_res_map( utility::vector1< core::Size > const & input_res_vector ){

	std::map< core::Size, core::Size > full_to_input_res_map;

	for ( Size n = 1; n <= input_res_vector.size(); n++ ) {
		full_to_input_res_map[input_res_vector[n]] = n;
	}

	return full_to_input_res_map;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
core::Size
string_to_int( std::string const & input_string ){

	Size int_of_string; //misnomer
	std::stringstream ss ( std::stringstream::in | std::stringstream::out );

	ss << input_string;

	if ( ss.fail() ) utility_exit_with_message( "In string_to_real(): ss.fail() for ss << input_string | string ( " + input_string + " )" );

	ss >> int_of_string;

	if ( ss.fail() ) utility_exit_with_message( "In string_to_real(): ss.fail() for ss >> int_of_string | string ( " + input_string + " )" );

	return int_of_string;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
core::Real
string_to_real( std::string const & input_string ){

	Real real_of_string;
	std::stringstream ss ( std::stringstream::in | std::stringstream::out );

	ss << input_string;

	if ( ss.fail() ) utility_exit_with_message( "In string_to_real(): ss.fail() for ss << input_string | string ( " + input_string + " )" );

	ss >> real_of_string;

	if ( ss.fail() ) utility_exit_with_message( "In string_to_real(): ss.fail() for ss >> real_of_string | string ( " + input_string + " )" );

	return real_of_string;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// AMW: cppcheck wants you to pass delimiters by reference, but don't try--it'll cause more problems than it's worth
utility::vector1< std::string >
tokenize( std::string const str, std::string delimiters ){
	using namespace std;

	utility::vector1< std::string > tokens;

	// Skip delimiters at beginning.
	string::size_type lastPos = str.find_first_not_of( delimiters, 0 );
	// Find first "non-delimiter".
	string::size_type pos     = str.find_first_of( delimiters, lastPos );

	while ( string::npos != pos || string::npos != lastPos ) {
		// Found a token, add it to the vector.
		tokens.push_back( str.substr( lastPos, pos - lastPos ) );
		// Skip delimiters.  Note the "not_of"
		lastPos = str.find_first_not_of( delimiters, pos );
		// Find next "non-delimiter"
		pos = str.find_first_of( delimiters, lastPos );
	}
	return tokens;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
is_virtual_base( conformation::Residue const & rsd ){

	using namespace chemical;

	if ( !rsd.is_RNA() ) return false;

	//Cytosine and Uracil contain 8 heavy base atoms. Adenine contains 10 heavy base atoms. Guanine contains 11 heavy base atoms.
	if ( ( rsd.nheavyatoms() - rsd.first_sidechain_atom() + 2 ) < 8 ) { //plus 2 since need to count both start and end atom.
		utility_exit_with_message( "The rna base " + name_from_aa( rsd.aa() ) + " contain lesser than 8 heavy atoms" );
	}

	Size non_virtual_atom_count = 0;
	for ( Size atomno = rsd.first_sidechain_atom() + 1; atomno <= rsd.nheavyatoms(); ++atomno ) { //iterate over base atoms....+1 to exclude the O2prime oxygen
		if ( !rsd.is_virtual( atomno )  ) {
			runtime_assert( rsd.atom_type( atomno ).name() != "VIRT" );
			non_virtual_atom_count++;
		}
	}

	bool const method_1 = ( non_virtual_atom_count < 8 );
	bool const method_2 = ( rsd.has_variant_type( core::chemical::VIRTUAL_RNA_RESIDUE ) ||
		rsd.has_variant_type( core::chemical::BULGE ) );

	if ( method_1 != method_2 ) {
		TR << "residue " << rsd.seqpos() << " non virtual atoms: " << non_virtual_atom_count << std::endl;
		output_boolean( "is_virtual_base determination by method_1: ", method_1, TR );
		output_boolean( "is_virtual_base determination by method_2: ", method_2, TR );
		TR << std::endl;
		utility_exit_with_message( "is_virtual_base determination by the two methods are not the same!" );
	}

	return method_1;

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Sept 10, 2010.OK actually don't need to past in the seq_pos...since that info is avialable in the residue object. Left the old version of the code below for backward compatibility.
void
setup_suite_atom_id_map( conformation::Residue const & rsd_1, conformation::Residue const & rsd_2, id::AtomID_Map < id::AtomID > & atom_ID_map, bool const base_only ){

	setup_suite_atom_id_map( rsd_1, rsd_2, rsd_1.seqpos(), rsd_2.seqpos(), atom_ID_map, base_only );

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
setup_suite_atom_id_map( conformation::Residue const & rsd_1,
	conformation::Residue const & rsd_2,
	Size const res_num_1,
	Size const res_num_2, //allow for the possibility that two poses have different sizes Jun 9, 2010
	id::AtomID_Map < id::AtomID > & atom_ID_map,
	bool const base_only ){

	using namespace ObjexxFCL;

	if ( name_from_aa( rsd_1.aa() ) != name_from_aa( rsd_2.aa() ) ) {
		utility_exit_with_message( "rsd_1.aa() != rsd_2.aa(). res_num_1 = " + string_of( res_num_1 ) + " name_from_aa( rsd_1.aa() ) = " + name_from_aa( rsd_1.aa() ) + " res_num_2 = " +  string_of( res_num_2 ) + " name_from_aa( rsd_2.aa() ) = " + name_from_aa( rsd_2.aa() ) );
	}

	Size const first_atom = ( base_only ) ? rsd_1.first_sidechain_atom() + 1 : 1; //+1 to exclude the O2prime oxygen

	for ( Size atomno_1 = first_atom; atomno_1 <= rsd_1.nheavyatoms(); ++atomno_1 ) {

		std::string const atom_name_1 = rsd_1.type().atom_name( atomno_1 );

		if ( !rsd_2.has( atom_name_1 ) ) continue;

		Size const atomno_2 = rsd_2.atom_index( atom_name_1 );

		//Check
		std::string const atom_name_2 = rsd_2.type().atom_name( atomno_2 );
		if ( atom_name_1 != atom_name_2 ) {
			utility_exit_with_message( "atom_name_1 != atom_name_2, atom_name_1 = " + atom_name_1 + " atom_name_2 = " + atom_name_2 );
		}

		if ( rsd_1.is_virtual( atomno_1 )  ) continue; //Check for virtual atoms
		if ( rsd_2.is_virtual( atomno_2 )  ) continue; //Check for virtual atoms

		id::AtomID const id1( atomno_1, res_num_1 );
		id::AtomID const id2( atomno_2, res_num_2 );
		atom_ID_map.set( id1, id2 );

	}

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Easiest if just align over the base atoms
//Virtual types mess up numbering...This will make sure that numbering is correct.
void
setup_suite_atom_id_map( pose::Pose const & pose_1, pose::Pose const & pose_2, Size const base_res, id::AtomID_Map < id::AtomID > & atom_ID_map, bool const base_only ){

	conformation::Residue const & base_rsd_1 = pose_1.residue( base_res );
	conformation::Residue const & base_rsd_2 = pose_2.residue( base_res );

	setup_suite_atom_id_map( base_rsd_1, base_rsd_2, base_res, base_res, atom_ID_map, base_only );

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////Dec 23, 2011.
void
setup_suite_atom_id_map( pose::Pose const & pose_1, pose::Pose const & pose_2, Size const base_res_1, Size const base_res_2,  id::AtomID_Map < id::AtomID > & atom_ID_map, bool const base_only ){

	conformation::Residue const & base_rsd_1 = pose_1.residue( base_res_1 );
	conformation::Residue const & base_rsd_2 = pose_2.residue( base_res_2 );

	setup_suite_atom_id_map( base_rsd_1, base_rsd_2, base_res_1, base_res_2, atom_ID_map, base_only );

}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
id::AtomID_Map < id::AtomID >
create_alignment_id_map_legacy( pose::Pose & mod_pose, pose::Pose const & ref_pose, utility::vector1< core::Size > const & rmsd_residue_list, bool const base_only ){
	using namespace chemical;

	id::AtomID_Map < id::AtomID > atom_ID_map;

	pose::initialize_atomid_map( atom_ID_map, mod_pose, id::BOGUS_ATOM_ID );

	if ( ref_pose.sequence() != mod_pose.sequence() ) {
		TR << "ref_pose.sequence() = " << ref_pose.sequence() << std::endl;
		TR << "mod_pose.sequence() = " << mod_pose.sequence() << std::endl;
		utility_exit_with_message( "ref_pose.sequence() != mod_pose.sequence()" );
	}

	for ( Size seq_num = 1; seq_num <= mod_pose.total_residue(); ++seq_num ) {
		if ( mod_pose.residue( seq_num ).aa() == core::chemical::aa_vrt ) continue; //Fang's electron density code.
		if ( !rmsd_residue_list.has_value( seq_num ) ) continue;
		setup_suite_atom_id_map( mod_pose, ref_pose, seq_num, atom_ID_map, base_only );
	}

	return atom_ID_map;

}

///////////////////////////////////////////////////////////////////////////////////////////////////
void
align_poses( core::pose::Pose & moving_pose, std::string const & moving_tag, core::pose::Pose const & static_pose, std::string const & static_tag, utility::vector1< core::Size > const & working_best_alignment, bool const base_only ){

	bool found_non_virtual_base = false;
	for ( Size n = 1; n <= working_best_alignment.size(); n++ ) {
		Size const seq_num = working_best_alignment[n];
		if ( is_virtual_base( moving_pose.residue( seq_num ) ) || is_virtual_base( static_pose.residue( seq_num ) ) ) continue;
		found_non_virtual_base = true; //ok found a non-virtual base nucleotide that can be used for alignment
		break;
	}

	if ( !found_non_virtual_base ) {
		for ( Size n = 1; n <= working_best_alignment.size(); n++ ) {
			Size const seq_num = working_best_alignment[n];
			TR.Debug << "seq_num = " << seq_num;
			output_boolean( "  is_virtual_base( " + moving_tag + " ):", is_virtual_base( moving_pose.residue( seq_num ) ), TR.Debug );
			output_boolean( "  is_virtual_base( " + static_tag + " ):", is_virtual_base( static_pose.residue( seq_num ) ), TR.Debug );
			TR.Debug << std::endl;
		}
		std::string error_message = "Error in aligning " + moving_tag + " to " + static_tag + ". No non - virtual_base in working_best_alignment to align the poses!";
		TR << error_message << std::endl;
		utility_exit_with_message( error_message );
	}

	//align current_pose to pose_output_list.
	id::AtomID_Map < id::AtomID > const & alignment_atom_id_map = create_alignment_id_map_legacy( moving_pose, static_pose, working_best_alignment, base_only );
	core::scoring::superimpose_pose( moving_pose, static_pose, alignment_atom_id_map );

	if ( check_for_messed_up_structure( moving_pose, moving_tag ) ) {
		std::string error_message = "Error in aligning " + moving_tag + " to " + static_tag + "!";
		TR << error_message << std::endl;
		utility_exit_with_message( moving_tag + " is messed up ...this is probably an alignment problem" );
	};
}

///////////////////////////////////////////////////////////////////////////////////////////////////
bool
seq_num_sort_criterion( core::Size seq_num_1, core::Size seq_num_2 ){
	return ( seq_num_1 < seq_num_2 );
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void
sort_seq_num_list( utility::vector1< core::Size > & seq_num_list ) {  //Low number on the top of the list
	sort( seq_num_list.begin(), seq_num_list.end(), seq_num_sort_criterion );
}


///////////////////////////////////////////////////////////////////////////////////////////////////
void
output_seq_num_list( std::string const & tag, utility::vector1< core::Size > const & seq_num_list, std::ostream & outstream /* = std::cout */, core::Size const spacing ){

	using namespace ObjexxFCL;
	using namespace ObjexxFCL::format;

	outstream <<  std::setw( spacing ) << tag;

	utility::vector1< core::Size > sorted_seq_num_list = seq_num_list;
	sort_seq_num_list( sorted_seq_num_list );

	Size seq_num = 1;
	for ( Size n = 1; n <= sorted_seq_num_list.size(); n++ ) {

		while ( seq_num < sorted_seq_num_list[n] ) {
			outstream << A( 4, " " );
			seq_num++;
		}
		outstream << I( 4, sorted_seq_num_list[n] );
		seq_num++;
	}

	outstream << std::endl;

}

///////////////////////////////////////////////////////////////////////////////////////////////////
void
output_title_text( std::string const & title, std::ostream & outstream /* = std::cout */ ){

	outstream << std::endl;

	Size title_length = title.size();
	Size char_per_line = 80;
	Size dash_length( 0 );
	if ( title_length < char_per_line ) dash_length = char_per_line - title_length;

	for ( Size i = 1; i <= dash_length/2; i++ )  outstream << "-";

	outstream << title;

	for ( Size i = 1; i <= dash_length/2; i++ ) outstream << "-";

	outstream << std::endl;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//April 15, 2012: Umm should switch to utility::file::file_exists( full_filename )!

bool
file_exists( std::string const & file_name ){

	std::ifstream my_file;

	my_file.open( file_name.c_str() );

	bool const file_exist = my_file.is_open();

	my_file.close();

	return file_exist;

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
remove_file( std::string const & file_name ){

	//Note that remove is part of #include <stdio.h>
	//perror: interprets the value of the global variable errno into a string and prints that string to stderr

	using namespace ObjexxFCL;

	if ( file_exists( file_name ) == false ) {
		utility_exit_with_message( "file_name ( " + file_name + " ) doesn't exist!" );
	}

	int const retcode = std::remove( file_name.c_str() );

	if ( retcode != 0 ) {
		std::string const error_message = "The following error occurs when attempting to remove file_name ( " + file_name + " )";

		std::perror( error_message.c_str() );
		utility_exit_with_message( error_message + ", retcode = " + string_of( retcode ) );
	}

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
output_rotamer( utility::vector1 < Real > & rotamer ){

	Size const spacing = 10;

	TR <<  std::setw( 18 ) << "Torsions = ";
	TR << std::setw( spacing ) << rotamer[1] << " ";
	TR << std::setw( spacing ) << rotamer[2] << " ";
	TR << std::setw( spacing ) << rotamer[3] << " ";
	TR << std::setw( spacing ) << rotamer[4] << " ";
	TR << std::setw( spacing ) << rotamer[5] << " ";
	TR << std::setw( spacing ) << rotamer[6] << " ";
	TR << std::setw( spacing ) << rotamer[7] << " ";
	TR << std::setw( spacing ) << rotamer[8] << " ";
	TR << std::setw( spacing ) << rotamer[9] << " ";
	TR << std::setw( spacing ) << rotamer[10] << " ";
	TR << std::setw( spacing ) << rotamer[11] << " ";
	TR << std::setw( spacing ) << rotamer[12] << " ";
	TR << std::setw( spacing ) << rotamer[13] << " ";
	TR << std::endl;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
add_virtual_O2Prime_hydrogen( core::pose::Pose & pose ){
	for ( core::Size i = 1; i <= pose.total_residue(); i++ ) {
		if ( !pose.residue( i ).is_RNA() ) continue;
		pose::add_variant_type_to_pose_residue( pose, core::chemical::VIRTUAL_O2PRIME_HYDROGEN, i );
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
check_instantiated_O2Prime_hydrogen( core::pose::Pose const & pose ){
	for ( core::Size i = 1; i <= pose.total_residue(); i++ ) {
		if ( !pose.residue_type( i ).is_RNA() ) continue;
		runtime_assert( !pose.residue_type( i ).has_variant_type( chemical::VIRTUAL_O2PRIME_HYDROGEN ) );
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
remove_virtual_O2Prime_hydrogen( pose::Pose & pose ){

	for ( Size i = 1; i <= pose.total_residue(); i++ ) {
		if ( !pose.residue( i ).is_RNA() ) continue;
		if ( pose.residue_type( i ).has_variant_type( core::chemical::VIRTUAL_O2PRIME_HYDROGEN ) ) {
			pose::remove_variant_type_from_pose_residue( pose, core::chemical::VIRTUAL_O2PRIME_HYDROGEN, i );
		}
	}
	return true;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//This now works for any rebuild residue.
Real
suite_rmsd( pose::Pose const & pose1, pose::Pose const & pose2, Size const & moving_res_num, bool const prepend_res, bool const ignore_virtual_atom ){

	Size atom_count = 0;
	Real sum_sd = 0;

	suite_square_deviation( pose1, pose2, prepend_res, moving_res_num, moving_res_num, atom_count, sum_sd, false, ignore_virtual_atom );

	sum_sd = sum_sd/( atom_count );
	Real rmsd = sqrt( sum_sd );

	// AMW: we just divided by atom_count, so it shouldn't be 0!
	if ( atom_count == 0 ) rmsd = 0.0; //special case...implement this on June_11, 2010, took me a whole day to debug this since buggy only on Biox compiler!

	return ( std::max( 0.01, rmsd ) );

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
Real
rmsd_over_residue_list(
	pose::Pose const & pose1,
	pose::Pose const & pose2,
	utility::vector1 < Size > const & residue_list, std::map< core::Size, core::Size > const & full_to_sub,
	std::map< core::Size, bool > const & is_prepend_map,
	bool const verbose,
	bool const ignore_virtual_atom )
{
	if ( verbose ) {
		output_title_text( "Enter rmsd_over_residue_list function", TR );
		output_boolean( "ignore_virtual_atom = ", ignore_virtual_atom, TR ); TR << std::endl;
		output_seq_num_list( "residue_list = ", residue_list, TR, 30 );
	}

	Size atom_count = 0;
	Real sum_sd = 0;

	for ( Size i = 1; i <= residue_list.size(); i++ ) {

		Size const full_seq_num = residue_list[i];
		Size const seq_num = full_to_sub.find( full_seq_num )->second;
		bool is_prepend = is_prepend_map.find( full_seq_num )->second;
		bool both_pose_res_is_virtual = false;
		if ( pose1.residue( seq_num ).has_variant_type( core::chemical::VIRTUAL_RNA_RESIDUE ) &&
				pose2.residue( seq_num ).has_variant_type( core::chemical::VIRTUAL_RNA_RESIDUE ) ) {
			both_pose_res_is_virtual = true;
		}

		if ( verbose ) {
			TR << "Full_seq_num = " << full_seq_num << " partial_seq_num = " << seq_num;
			output_boolean( " is_prepend = ", is_prepend, TR );
			output_boolean( " both_pose_res_is_virtual = ", both_pose_res_is_virtual, TR ); TR << std::endl;
		}

		if ( both_pose_res_is_virtual ) continue;

		//add atom in the suites to atom_count
		//add sd of each atom to sum_sd
		suite_square_deviation( pose1, pose2, is_prepend, seq_num, seq_num, atom_count, sum_sd, verbose, ignore_virtual_atom );

	}

	sum_sd = sum_sd/( atom_count );
	Real rmsd = sqrt( sum_sd );

	// AMW: we just divided by atom_count, so it shouldn't be 0!
	if ( atom_count == 0 ) rmsd = 0.0; //special case...implement this on May 5, 2010

	if ( verbose ) {
		TR << "sum_sd = " << sum_sd << " atom_count = " << atom_count << " rmsd = " << rmsd << std::endl;
		output_title_text( "Exit rmsd_over_residue_list function", TR );
	}

	return ( std::max( 0.01, rmsd ) );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////


Real
rmsd_over_residue_list( pose::Pose const & pose1, pose::Pose const & pose2, working_parameters::StepWiseWorkingParametersCOP working_parameters_, bool const ignore_virtual_atom ){

	utility::vector1 < core::Size > const & calc_rms_res = working_parameters_->calc_rms_res();
	std::map< core::Size, core::Size > const & full_to_sub = working_parameters_->const_full_to_sub();
	std::map< core::Size, bool > const & is_prepend_map = working_parameters_->is_prepend_map();

	return rmsd_over_residue_list( pose1, pose2, calc_rms_res, full_to_sub, is_prepend_map, false /*verbose*/, ignore_virtual_atom );

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
print_heavy_atoms( Size const & suite_num_1, Size const & suite_num_2, pose::Pose const & pose1, pose::Pose const & pose2 ){

	using namespace conformation;
	Size num_atoms;

	num_atoms = std::max( pose1.residue( suite_num_1 ).nheavyatoms(), pose2.residue( suite_num_2 ).nheavyatoms() );

	TR << "num_atoms: " << num_atoms << std::endl;

	for ( Size n = 1; n <= num_atoms;  n++ ) {

		TR << " atom num = " <<  n;
		TR << "  atom_name of the pose1 " <<  pose1.residue( suite_num_1 ).atom_name( n );
		TR << "  atom_name of the pose2 " <<  pose2.residue( suite_num_2 ).atom_name( n ) << std::endl;

	}
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
Size
get_num_side_chain_atom_from_res_name( chemical::AA const & res_aa, bool const verbose ){

	Size num_side_chain_atom;

	if ( name_from_aa( res_aa ) == "RAD" ) {
		if ( verbose ) TR << "name_from_aa: RAD" << std::endl;
		num_side_chain_atom = 11;
	} else if ( name_from_aa( res_aa ) == "RCY" ) {
		if ( verbose ) TR << "name_from_aa: RCY" << std::endl;
		num_side_chain_atom = 9;
	} else if ( name_from_aa( res_aa ) == "RGU" ) {
		if ( verbose ) TR << "name_from_aa: RGU" << std::endl;
		num_side_chain_atom = 12;
	} else if ( name_from_aa( res_aa ) == "URA" ) {
		if ( verbose ) TR << "name_from_aa: URA" << std::endl;
		num_side_chain_atom = 9;
	} else {
		TR << "Error, cannot identify residue type" << std::endl;
		num_side_chain_atom = 0;
		exit ( 1 );
	}

	return num_side_chain_atom;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
Real
atom_square_deviation( conformation::Residue const & rsd_1, conformation::Residue const & rsd_2, Size const & atomno_1, Size const & atomno_2, bool verbose ){

	//////////////Might take this out later if it significantly slow down the code///////////////
	std::string const & atom_name_1 = rsd_1.type().atom_name( atomno_1 );
	std::string const & atom_name_2 = rsd_2.type().atom_name( atomno_2 );

	//This can be turned on for debugging, but by default is mod out since string comparison might be slow Parin Jan 28, 2009
	if ( atom_name_1 != atom_name_2 ) {
		utility_exit_with_message( "atom_name_1 != atom_name_2, atom_name_1 = " + atom_name_1 + " atom_name_2 = " + atom_name_2 );
	}

	//  if(rsd_1.is_virtual(atomno_1)   || rsd_2.is_virtual(atomno_2) ) {
	//   TR << "atom_name_1= " << atom_name_1 << " atom_name_2= " << atom_name_2 << std::endl;
	//   utility_exit_with_message( "rsd_1.atom_type(n).name()==\"VIRT\"  || rsd_2.atom_type(n).name()==\"VIRT\" =TRUE!");
	//  }
	////////////////////////////////////////////////////////////////////////////////////////////////

	Distance const dist_squared = ( rsd_1.xyz( atomno_1 ) - rsd_2.xyz( atomno_2 ) ).length_squared();


	if ( verbose ) {
		TR << " atom_name of the atom1 = " << atom_name_1 << " " << rsd_1.seqpos();
		TR << " atom_name of the atom2 = " << atom_name_2 << " " << rsd_2.seqpos();
		TR << " Dist_squared = " << dist_squared << std::endl;
	}

	return dist_squared;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//create on Sept 24, 2010...should really integrate these rmsd square deviation functions....cause now it is just copy and paste copy
void
base_atoms_square_deviation( pose::Pose const & pose1, pose::Pose const & pose2, Size const & moving_res_1, Size const & moving_res_2, Size& atom_count, Real& sum_sd, bool verbose, bool const ignore_virtual_atom ){


	chemical::AA const & res_aa =  pose1.residue( moving_res_1 ).aa();
	chemical::AA const & res_aa2 =  pose2.residue( moving_res_2 ).aa();

	if ( res_aa != res_aa2 ) utility_exit_with_message( "res_aa ( " + name_from_aa( res_aa ) + " ) != res_aa2 ( " + name_from_aa( res_aa2 ) + " ) " );

	Size const first_sidechain_atom1 = pose1.residue( moving_res_1 ).first_sidechain_atom();
	Size const first_sidechain_atom2 = pose2.residue( moving_res_2 ).first_sidechain_atom();

	Size const num_side_chain_atom = get_num_side_chain_atom_from_res_name( res_aa, verbose );

	if ( verbose ) TR << " MOVING_RES_1: " << moving_res_1 << " MOVING_RES_2: " << moving_res_2 << std::endl;

	//Need to use num_side_chain_atom from pose1 since a silly bug in Rosetta miscalculate num_heavy_atom by considering
	//the virtaul O2prime hydrogen to be heavy_atom when it is set to virtual in the current_pose_screen
	for ( Size n = 1; n <= num_side_chain_atom; n++ ) { //This INCLUDE the O2prime oxygen

		Size const atomno_1 = ( n - 1 ) + first_sidechain_atom1;
		Size const atomno_2 = ( n - 1 ) + first_sidechain_atom2;

		conformation::Residue const & rsd_1 = pose1.residue( moving_res_1 );
		conformation::Residue const & rsd_2 = pose2.residue( moving_res_2 );

		if ( rsd_1.type().atom_name( atomno_1 ) == " O2'" ) continue; //Exclude the O2prime oxygen
		if ( rsd_2.type().atom_name( atomno_2 ) == " O2'" ) continue; //Exclude the O2prime oxygen

		if ( ignore_virtual_atom ) {
			if ( rsd_1.is_virtual( atomno_1 )   || rsd_2.is_virtual( atomno_2 )  ) continue;
		}

		if ( rsd_1.is_virtual( atomno_1 )   && rsd_2.is_virtual( atomno_2 )  ) { //Change this to "AND" on Apr 5
			//     TR << "Both atoms are VIRTUAL! moving_res_1= " << moving_res_1 << " moving_res_2= " << moving_res_2 << " atomno_1= " << atomno_1 << " atomno_2= " << atomno_2 << std::endl;
			continue;
		}

		atom_count++;
		sum_sd = sum_sd + atom_square_deviation( rsd_1, rsd_2, atomno_1, atomno_2, verbose );
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
phosphate_square_deviation( pose::Pose const & pose1, pose::Pose const & pose2, Size const & moving_res_1, Size const & moving_res_2, Size& atom_count, Real& sum_sd, bool verbose, bool const ignore_virtual_atom ){


	chemical::AA const & res_aa =  pose1.residue( moving_res_1 ).aa();
	chemical::AA const & res_aa2 =  pose2.residue( moving_res_2 ).aa();

	if ( res_aa != res_aa2 ) utility_exit_with_message( "res_aa ( " + name_from_aa( res_aa ) + " ) != res_aa2 ( " + name_from_aa( res_aa2 ) + " ) " );

	if ( verbose ) TR << " MOVING_RES_1: " << moving_res_1 << " MOVING_RES_2: " << moving_res_2 << std::endl;

	for ( Size atomno = 1; atomno <= 4; atomno++ ) {

		conformation::Residue const & rsd_1 = pose1.residue( moving_res_1 );
		conformation::Residue const & rsd_2 = pose2.residue( moving_res_2 );

		if ( ignore_virtual_atom ) {
			if ( rsd_1.is_virtual( atomno )   || rsd_2.is_virtual( atomno )  ) continue;
		}

		if ( rsd_1.is_virtual( atomno )   && rsd_2.is_virtual( atomno )  ) continue;

		atom_count++;
		sum_sd = sum_sd + atom_square_deviation( rsd_1, rsd_2, atomno, atomno, verbose );

	}

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// unify with suite_rmsd (see below).
core::Real
phosphate_base_phosphate_rmsd( pose::Pose const & pose1, pose::Pose const & pose2, Size const & moving_res_num, bool const ignore_virtual_atom ){

	Size atom_count = 0;
	Real sum_sd = 0;

	phosphate_base_phosphate_square_deviation( pose1, pose2, moving_res_num, moving_res_num, atom_count, sum_sd, false, ignore_virtual_atom );

	sum_sd = sum_sd/( atom_count );
	Real rmsd = sqrt( sum_sd );

	// AMW: we just divided by atom_count, so it shouldn't be 0!
	if ( atom_count == 0 ) rmsd = 0.0; //special case...implement this on June_11, 2010, took me a whole day to debug this since buggy only on Biox compiler!

	return ( std::max( 0.01, rmsd ) );

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// should be unified with suite_square_deviation (see next).
void
phosphate_base_phosphate_square_deviation( pose::Pose const & pose1, pose::Pose const & pose2, Size const & moving_res_1, Size const & moving_res_2, Size& atom_count, Real& sum_sd, bool verbose, bool const ignore_virtual_atom ){

	chemical::AA const & res_aa  =  pose1.residue( moving_res_1 ).aa();
	chemical::AA const & res_aa2 =  pose2.residue( moving_res_2 ).aa();
	runtime_assert( res_aa == res_aa2 );

	Size const first_sidechain_atom1 = pose1.residue( moving_res_1 ).first_sidechain_atom();
	Size const first_sidechain_atom2 = pose2.residue( moving_res_2 ).first_sidechain_atom();

	Size const num_side_chain_atom = get_num_side_chain_atom_from_res_name( res_aa, verbose );

	Size const num_heavy_backbone_atoms = 11; //RNA contain 11 heavy backbone atoms.

	for ( Size atomno = 1; atomno <= num_heavy_backbone_atoms; atomno++ ) {

		Size const res_count = ( atomno <= 4 ) ? 2 : 1; //atom 1-4 are " P  ", " OP2", " OP1" and " O5'"

		// check out this extremely old bug! fixed in rd2014, after deprecation of the function.
		//   for ( Size ii = 1; ii < res_count; ii++ ){
		for ( Size ii = 1; ii <= res_count; ii++ ) {

			Size const res_num_1 = moving_res_1 + ( ii - 1 );
			Size const res_num_2 = moving_res_2 + ( ii - 1 );

			conformation::Residue const & rsd_1 = pose1.residue( res_num_1 );
			conformation::Residue const & rsd_2 = pose2.residue( res_num_2 );

			if ( ignore_virtual_atom && ( rsd_1.is_virtual( atomno )   || rsd_2.is_virtual( atomno ) ) ) continue;
			if ( rsd_1.is_virtual( atomno )   && rsd_2.is_virtual( atomno )  ) continue;

			atom_count++;
			sum_sd = sum_sd + atom_square_deviation( rsd_1, rsd_2, atomno, atomno, verbose );
		}
	}

	//Need to use num_side_chain_atom from pose1 since a silly bug in Rosetta miscalculate num_heavy_atom by considering
	//the virtaul O2prime hydrogen to be heavy_atom when it is set to virtual in the current_pose_screen
	for ( Size n = 1; n <= num_side_chain_atom; n++ ) { //INCLUDE the O2prime oxygen

		Size const atomno_1 = ( n - 1 ) + first_sidechain_atom1;
		Size const atomno_2 = ( n - 1 ) + first_sidechain_atom2;

		conformation::Residue const & rsd_1 = pose1.residue( moving_res_1 );
		conformation::Residue const & rsd_2 = pose2.residue( moving_res_2 );

		if ( ignore_virtual_atom && ( rsd_1.is_virtual( atomno_1 ) || rsd_2.is_virtual( atomno_2 ) ) ) continue;
		if ( rsd_1.is_virtual( atomno_1 )   && rsd_2.is_virtual( atomno_2 )  ) continue;

		atom_count++;
		sum_sd = sum_sd + atom_square_deviation( rsd_1, rsd_2, atomno_1, atomno_2, verbose );
	}

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Could get the atom_num using the atom_name but I don't think this will be fast Size atom_P = rsd.atom_index( " P  " ); Jan 28, 2010 Parin S///////////
//   (This could/should be made more robust in order to handle noncanonical bases.
//    Look at code in StepWisePoseAligner -- rhiju )
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
suite_square_deviation( pose::Pose const & pose1, pose::Pose const & pose2, bool const & prepend_res,
	Size const & moving_res_1, Size const & moving_res_2,
	Size & atom_count, Real & sum_sd, bool verbose, bool const ignore_virtual_atom ){

	chemical::AA const & res_aa  =  pose1.residue( moving_res_1 ).aa();
	chemical::AA const & res_aa2 =  pose2.residue( moving_res_2 ).aa();
	runtime_assert( res_aa == res_aa2 );

	Size const first_sidechain_atom1 = pose1.residue( moving_res_1 ).first_sidechain_atom();
	Size const first_sidechain_atom2 = pose2.residue( moving_res_2 ).first_sidechain_atom();

	Size const num_side_chain_atom = get_num_side_chain_atom_from_res_name( res_aa, verbose );

	Size const num_heavy_backbone_atoms = 11; //RNA contains 11 heavy backbone atoms.

	for ( Size atomno = 1; atomno <= num_heavy_backbone_atoms; atomno++ ) {

		//atom 1-4 are " P  ", " OP2", " OP1" and " O5'"
		Size const res_num_1 = ( prepend_res && atomno <= 4 ) ? moving_res_1 + 1: moving_res_1;
		Size const res_num_2 = ( prepend_res && atomno <= 4 ) ? moving_res_2 + 1: moving_res_2;

		conformation::Residue const & rsd_1 = pose1.residue( res_num_1 );
		conformation::Residue const & rsd_2 = pose2.residue( res_num_2 );

		if ( ignore_virtual_atom && ( rsd_1.is_virtual( atomno )   || rsd_2.is_virtual( atomno ) ) ) continue;
		if ( rsd_1.is_virtual( atomno )   && rsd_2.is_virtual( atomno )  ) continue;

		atom_count++;
		sum_sd = sum_sd + atom_square_deviation( rsd_1, rsd_2, atomno, atomno, verbose );
	}

	//Need to use num_side_chain_atom from pose1 since a silly bug in Rosetta miscalculates num_heavy_atom by considering
	//the virtual O2prime hydrogen to be heavy_atom when it is set to virtual in the current_pose_screen
	for ( Size n = 1; n <= num_side_chain_atom; n++ ) { //INCLUDE the O2prime oxygen

		Size const atomno_1 = ( n - 1 ) + first_sidechain_atom1;
		Size const atomno_2 = ( n - 1 ) + first_sidechain_atom2;

		conformation::Residue const & rsd_1 = pose1.residue( moving_res_1 );
		conformation::Residue const & rsd_2 = pose2.residue( moving_res_2 );

		if ( ignore_virtual_atom && ( rsd_1.is_virtual( atomno_1 ) || rsd_2.is_virtual( atomno_2 ) ) ) continue;
		if ( rsd_1.is_virtual( atomno_1 )   && rsd_2.is_virtual( atomno_2 )  ) continue;

		atom_count++;
		sum_sd = sum_sd + atom_square_deviation( rsd_1, rsd_2, atomno_1, atomno_2, verbose );
	}

}

//////////////////////////////////////////////////////////////////////////////////////////////////////
void
freeze_sugar_torsions( core::kinematics::MoveMap & mm, Size const total_residue ){

	using namespace core::id;

	TR << "Freeze pose sugar torsions, total_residue = " << total_residue << std::endl;
	for ( Size i = 1; i <= total_residue; i++ ) {
		mm.set( TorsionID( i , id::BB,  4 ), false ); //delta_i
		mm.set( TorsionID( i , id::CHI, 2 ), false ); //nu2_i
		mm.set( TorsionID( i , id::CHI, 3 ), false ); //nu1_i
	}

}

/////////////////////////////////////////////////////////////////////////////////////
utility::vector1< core::Size >
get_surrounding_O2prime_hydrogen( pose::Pose const & pose, utility::vector1< core::Size > const & moving_res, bool verbose /*= false */ ){

	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::kinematics;
	using namespace ObjexxFCL;

	//Consistency_check
	for ( Size seq_num = 1; seq_num <= pose.total_residue(); seq_num++ ) {

		if ( !pose.residue_type( seq_num ).is_RNA() ) continue;

		// these are actually redundant with virtual 2'-OH atom-wise check below.
		if ( pose.residue_type( seq_num ).has_variant_type( chemical::VIRTUAL_RNA_RESIDUE ) ) continue;
		if ( pose.residue_type( seq_num ).has_variant_type( chemical::VIRTUAL_RIBOSE ) ) continue;

		// allow for 2'-OH virtualization during packing. Decide during packer task setup, not in here.
		// if ( pose.residue_type( seq_num ).has_variant_type( chemical::VIRTUAL_O2PRIME_HYDROGEN ) ) continue;

		core::conformation::Residue const & rsd = pose.residue( seq_num );
		Size const at = rsd.first_sidechain_atom();
		runtime_assert ( rsd.type().atom_name( at ) == " O2'" );
	}

	utility::vector1< bool > is_O2prime_hydrogen_virtual_list;

	for ( Size seq_num = 1; seq_num <= pose.total_residue(); seq_num++ ) {

		if ( !pose.residue( seq_num ).is_RNA() ) {
			if ( verbose ) TR << "res " << seq_num << " is not RNA " << std::endl;
			is_O2prime_hydrogen_virtual_list.push_back( false ); //false since not virtual O2prime_hydrogen
			continue;
		}

		core::conformation::Residue const & rsd = pose.residue( seq_num );
		Size at = rsd.atom_index( "HO2'" );

		if ( rsd.is_virtual( at )  ) {
			if ( verbose ) TR << "res " << seq_num << " has a virtual o2prime hydrogen! " << std::endl;
			is_O2prime_hydrogen_virtual_list.push_back( true );
		} else {
			is_O2prime_hydrogen_virtual_list.push_back( false );
		}
	}

	//March 17, 2012 extra precaution.
	runtime_assert ( is_O2prime_hydrogen_virtual_list.size() == pose.total_residue() );

	utility::vector1< core::Size > surrounding_O2prime_hydrogen;
	Size num_o2prime_moving_res = 0;
	for ( Size n = 1; n <= moving_res.size(); n++ ) {
		Size const seq_num = moving_res[n];

		if ( !pose.residue_type( seq_num ).is_RNA() ) continue;
		if ( is_O2prime_hydrogen_virtual_list[seq_num] ) continue;

		surrounding_O2prime_hydrogen.push_back( seq_num );
		num_o2prime_moving_res++;
	}


	//1st layer, interaction between surrounding O2prime and moving_res
	for ( Size seq_num = 1; seq_num <= pose.total_residue(); seq_num++ ) {

		if ( !pose.residue_type( seq_num ).is_RNA() ) continue;
		if ( surrounding_O2prime_hydrogen.has_value( seq_num ) ) continue;

		bool is_surrounding_res = false;

		core::conformation::Residue const & surrounding_rsd = pose.residue( seq_num );
		Size const surr_at = surrounding_rsd.first_sidechain_atom();

		if ( is_O2prime_hydrogen_virtual_list[seq_num] ) continue;

		//3.5 Angstrom for O2prime-normal

		for ( Size ii = 1; ii <= moving_res.size(); ii++ ) {
			if ( is_surrounding_res == true ) break;

			core::conformation::Residue const & moving_rsd = pose.residue( moving_res[ii] );

			for ( Size moving_at = 1; moving_at <= moving_rsd.natoms(); moving_at++ ) {

				if ( moving_rsd.is_virtual( moving_at )  ) continue;

				Real const cutoff_dist = ( moving_at == moving_rsd.first_sidechain_atom() ) ? 4.5: 3.5 ;

				//4.5 Angstrom interaction dist_cutoff for surrounding O2prime- moving_res O2prime
				//3.5 Angstrom interaction dist_cutoff for surround O2prime and other atoms of moving_res

				Real const dist_squared = ( surrounding_rsd.xyz( surr_at ) - moving_rsd.xyz( moving_at ) ).length_squared();

				if ( dist_squared  < cutoff_dist*cutoff_dist ) {
					if ( verbose ) TR << "res " << seq_num << " is a 1st layer surrounding O2prime_hydrogen res, dist_squared = " << dist_squared << std::endl;
					surrounding_O2prime_hydrogen.push_back( seq_num );
					is_surrounding_res = true;
					break;
				}
			}
		}
	}


	//layer 2+, interaction between surrounding O2prime hydrogen themselves
	Size layer_num = 2;
	while ( true ) {
		bool add_new_O2prime_hydrogen = false;

		for ( Size seq_num = 1; seq_num <= pose.total_residue(); seq_num++ ) {

			if ( !pose.residue_type( seq_num ).is_RNA() ) continue;
			if ( surrounding_O2prime_hydrogen.has_value( seq_num ) ) continue;

			core::conformation::Residue const & rsd_1 = pose.residue( seq_num );
			Size const at_1 = rsd_1.first_sidechain_atom();

			if ( is_O2prime_hydrogen_virtual_list[seq_num] ) continue;

			for ( Size ii = 1; ii <= surrounding_O2prime_hydrogen.size(); ii++ ) {

				core::conformation::Residue const & rsd_2 = pose.residue( surrounding_O2prime_hydrogen[ii] );
				Size const at_2 = rsd_2.first_sidechain_atom();

				Real const cutoff_dist = 4.5;      //4.5 Angstrom interaction dist_cutoff for surrounding O2prime themselves.

				Real const dist_squared = ( rsd_1.xyz( at_1 ) - rsd_2.xyz( at_2 ) ).length_squared();

				if ( dist_squared  < cutoff_dist*cutoff_dist ) {
					if ( verbose ) TR << "res " << seq_num << " is layer " << layer_num << " surrounding O2prime_hydrogen res, dist_squared = " << dist_squared << std::endl;
					surrounding_O2prime_hydrogen.push_back( seq_num );
					add_new_O2prime_hydrogen = true;
					break;
				}
			}
		}

		layer_num++;
		if ( !add_new_O2prime_hydrogen ) break;
	}

	//consistency_check..
	for ( Size ii = 1; ii <= surrounding_O2prime_hydrogen.size(); ii++ ) {

		Size const seq_num = surrounding_O2prime_hydrogen[ii];

		if ( is_O2prime_hydrogen_virtual_list[seq_num] ) {
			std::string const exit_message = "surrounding_O2prime_hydrogen res " + string_of( seq_num ) + " has a virtual o2prime hydrogen!! ";
			utility_exit_with_message( exit_message );
		}
	}

	if ( verbose ) {
		TR << "num_o2prime_moving_res = " << num_o2prime_moving_res << "  surrounding_O2prime_hydrogen.size() = " << surrounding_O2prime_hydrogen.size() << std::endl;
	}

	return surrounding_O2prime_hydrogen;

}

/////////////////////////////////////////////////////////////////////////////////////
void
o2prime_trials( pose::Pose & pose, core::scoring::ScoreFunctionCOP const & packer_scorefxn,
	bool const pack_virtual_o2prime_hydrogen /* = false */ ){ //O2prime pack every position..

	utility::vector1< core::Size > O2prime_pack_seq_num;

	for ( Size seq_num = 1; seq_num <= pose.total_residue(); seq_num++ ) {
		if ( !pose.residue_type( seq_num ).is_RNA() ) continue;
		O2prime_pack_seq_num.push_back( seq_num );
	}

	o2prime_trials( pose, packer_scorefxn, O2prime_pack_seq_num, pack_virtual_o2prime_hydrogen );

}

/////////////////////////////////////////////////////////////////////////////////////
void
o2prime_trials( pose::Pose& pose, core::scoring::ScoreFunctionCOP const & packer_scorefxn,
	utility::vector1< core::Size > const & O2prime_pack_seq_num,
	bool const pack_virtual_o2prime_hydrogen /* = false */ ){

	output_seq_num_list( "O2prime_pack_seq_num = ", O2prime_pack_seq_num, TR.Debug );

	pack::task::PackerTaskOP task = create_standard_o2prime_pack_task( pose, O2prime_pack_seq_num, pack_virtual_o2prime_hydrogen );

	pack::rotamer_trials( pose, *packer_scorefxn, task );

}

/////////////////////////////////////////////////////////////////////////////////////
//created on Jan 20, 2012: Meant as the standard function to setup o2prime_pack_task to be called by both StepWiseRNA_ResidueSampler and StepWiseRNA_Minimizer...and in the future also FARFAR_minimizer.
pack::task::PackerTaskOP
create_standard_o2prime_pack_task( pose::Pose const & pose, utility::vector1< core::Size > const & O2prime_pack_seq_num,
	bool const pack_virtual_o2prime_hydrogen /* = false */ ){

	pack::task::PackerTaskOP o2prime_pack_task = pack::task::TaskFactory::create_packer_task( pose );

	for ( Size seq_num = 1; seq_num <= pose.total_residue(); seq_num++ ) {

		if ( O2prime_pack_seq_num.has_value( seq_num ) &&
				pose.residue( seq_num ).is_RNA() &&
				( pack_virtual_o2prime_hydrogen || !pose.residue_type( seq_num ).has_variant_type( chemical::VIRTUAL_O2PRIME_HYDROGEN ) ) ) {
			//pack this residue!

			/*
			NOTE: the ex4 option allows modeler of extra torsion angles ( "ex" ) at
			the 2' - OH ( which is torsion "4" ). To pack the bases, you would have to
			add "ex1"; it should not happen! If you see bases getting repacked,
			that is a bug, and please let me know. ( From Rhiju email response )

			ex4 ( HO2prime ) is "on" by default | ex1 ( base chi ) is "off" by default
			*/
			o2prime_pack_task->nonconst_residue_task( seq_num ).and_extrachi_cutoff( 0 );
			o2prime_pack_task->nonconst_residue_task( seq_num ).or_ex4( true ); //extra O2prime modeler
			o2prime_pack_task->nonconst_residue_task( seq_num ).or_include_current( true );
			o2prime_pack_task->nonconst_residue_task( seq_num ).or_include_virtual_side_chain( pack_virtual_o2prime_hydrogen );

		} else {
			o2prime_pack_task->nonconst_residue_task( seq_num ).prevent_repacking();
		}
	}

	return o2prime_pack_task;
}

//////////////////////////////////////////////////////////////////////////
void
setup_chain_break_variants( core::pose::Pose & pose,  Size const cutpoint ){
	pose::correctly_add_cutpoint_variants( pose, cutpoint );
}

//////////////////////////////////////////////////////////////////////////
void
remove_chain_break_variants( core::pose::Pose & pose,  Size const & cutpoint ){
	runtime_assert( is_cutpoint_closed( pose, cutpoint ) );
	pose::remove_variant_type_from_pose_residue( pose, chemical::CUTPOINT_LOWER, cutpoint );
	pose::remove_variant_type_from_pose_residue( pose, chemical::CUTPOINT_UPPER, cutpoint + 1 );
}

///////////////////////////////////////////////////////////////////////
utility::vector1< bool >
get_partition_definition_floating_base( pose::Pose const & pose, Size const & moving_res ){
	Size const jump_nr = look_for_unique_jump_to_moving_res( pose.fold_tree(), moving_res );
	return get_partition_definition_by_jump( pose, jump_nr);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Size
get_anchor_res( Size const rebuild_res, pose::Pose const & pose ){
	kinematics::FoldTree const & f = pose.fold_tree();
	Size const jump_nr = look_for_unique_jump_to_moving_res( f, rebuild_res );
	return ( f.upstream_jump_residue( jump_nr ) == static_cast<int>(rebuild_res) ) ?
		f.downstream_jump_residue( jump_nr ) : f.upstream_jump_residue( jump_nr );
}
////////////////////////////////////////////////////////////////////////
bool
check_for_messed_up_structure( core::pose::Pose const & pose, std::string const & tag ){
	using namespace core::scoring;
	using namespace core::chemical::rna;

	for ( Size seq_num = 1; seq_num <= pose.total_residue(); seq_num++ ) {

		if ( !pose.residue( seq_num ).is_RNA() ) continue;

		conformation::Residue const & rsd( pose.residue( seq_num ) );
		Real delta = numeric::principal_angle_degrees( rsd.mainchain_torsion( DELTA ) );
		Real chi = numeric::principal_angle_degrees( rsd.chi( 1 ) );
		Real nu_2 = numeric::principal_angle_degrees( rsd.chi( 2 ) );
		Real nu_1 = numeric::principal_angle_degrees( rsd.chi( 3 ) );

		//   TR.Debug << " tag= " << tag << " seq_num= " << seq_num << " delta= " << delta << " chi= " << chi << " nu_2= " << nu_2 << " nu_1= " << nu_1 << std::endl;

		if ( ( delta >  - 0.001 && delta < 0.001 ) || ( nu_2 >  - 0.001 && nu_2 < 0.001 ) || ( nu_1 >  - 0.001 && nu_1 < 0.001 ) ) { //observation is that messed up structure will have delta value of zero
			TR.Debug << "Warning: " << tag << " is probably a messed up pose, will be ignored" << std::endl;
			TR.Debug << " seq_num = " << seq_num << " delta = " << delta << " chi = " << chi << " nu_2 = " << nu_2 << " nu_1 = " << nu_1 << std::endl;
			if ( ( rsd.has_variant_type( core::chemical::VIRTUAL_RNA_RESIDUE ) ) ||
					( rsd.has_variant_type( core::chemical::VIRTUAL_RIBOSE ) ) ) { //Implement on Oct 28,2010
				TR.Debug << "OK lets NOT ignore yet since this rsd has virtual_res or virtual_sugar variant type..continue checking.. " << std::endl;
			} else {
				return true;
			}
		}
	}
	return false;
}


void
sleep( core::Size mseconds )
{
	clock_t endwait;
	endwait = clock () + mseconds * CLOCKS_PER_SEC/1000 ;
	while ( clock() < endwait ) {}
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool
is_residues_in_contact( core::Size const & res_ONE, core::pose::Pose const & pose_ONE, core::Size const & res_TWO, core::pose::Pose const & pose_TWO, core::Real const atom_atom_overlap_dist_cutoff, core::Size const num_atom_contacts_cutoff, bool const verbose ){

	using namespace ObjexxFCL;

	core::conformation::Residue const & rsd_ONE = pose_ONE.residue( res_ONE );
	core::conformation::Residue const & rsd_TWO = pose_TWO.residue( res_TWO );

	Size num_atom_contacts_so_far = 0;

	if ( num_atom_contacts_cutoff < 1 ) {
		utility_exit_with_message( "num_atom_contacts_cutoff( " + string_of( num_atom_contacts_cutoff ) + " ) < 1!" );
	}

	for ( Size at_ONE = 1; at_ONE <= rsd_ONE.natoms(); at_ONE++ ) { //include hydrogen atoms
		for ( Size at_TWO = 1; at_TWO <= rsd_TWO.natoms(); at_TWO++ ) { //include hydrogen atoms

			if ( rsd_ONE.is_virtual( at_ONE )  ) continue;
			if ( rsd_TWO.is_virtual( at_TWO )  ) continue;

			Real const VDW_radius_ONE = rsd_ONE.atom_type( at_ONE ).lj_radius();
			Real const VDW_radius_TWO = rsd_TWO.atom_type( at_TWO ).lj_radius();

			Real const cutoff_sum_VDW_radius = VDW_radius_ONE + VDW_radius_TWO - atom_atom_overlap_dist_cutoff;

			if ( cutoff_sum_VDW_radius < 0 ) utility_exit_with_message( "( VDW_radius_ONE + VDW_radius_TWO - atom_atom_overlap_dist_cutoff ) < 0!!" );

			Real const atom_atom_dist_squared = ( rsd_ONE.xyz( at_ONE ) - rsd_TWO.xyz( at_TWO ) ).length_squared();

			if ( atom_atom_dist_squared < ( cutoff_sum_VDW_radius*cutoff_sum_VDW_radius ) ) {
				num_atom_contacts_so_far++;
				if ( verbose ) {
					TR << "res_ONE = " << res_ONE << " res_TWO = " << res_TWO << " num_atom_contacts_so_far = " << num_atom_contacts_so_far << "|";
					TR << " VDW_radius_ONE = " << VDW_radius_ONE << " VDW_radius_TWO = " << VDW_radius_TWO;
					TR << " atom_atom_overlap_dist_cutoff = " << atom_atom_overlap_dist_cutoff << " cutoff_sum_VDW_radius = " << cutoff_sum_VDW_radius;
					TR << " atom_atom_dist = " << ( rsd_ONE.xyz( at_ONE ) - rsd_TWO.xyz( at_TWO ) ).length();
					TR << " " << rsd_ONE.atom_name( at_ONE ) << " of res_ONE and " << rsd_TWO.atom_name( at_TWO ) << " of res_TWO are in contact! " << std::endl;
				}
			}

			if ( num_atom_contacts_so_far == num_atom_contacts_cutoff ) {
				return true;
			}

			if ( num_atom_contacts_so_far > num_atom_contacts_cutoff ) { //consistency_check
				utility_exit_with_message( "num_atom_contacts_so_far( " + string_of( num_atom_contacts_so_far ) + " ) > num_atom_contacts_cutoff( " + string_of( num_atom_contacts_cutoff ) + " )" );
			}


		}
	}

	return false;

}

//////////////This function was originally part of SWA_RNA_Sampler and used for the Richardson code. Move to Util, so that FloatingBaseSamplerUtil will have access to it////////////
void
set_CCD_torsions_to_zero( core::pose::Pose & pose, Size const five_prime_res ){

	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::id;
	//TR << "set_CCD_torsions_to_zero" << std::endl;

	//Size const five_prime_res = working_parameters_->five_prime_chain_break_res();
	Size const three_prime_res = five_prime_res + 1;

	//Even through there is the chain_break, alpha of 3' and epl and gamma of 5' should be defined due to the existence of the upper and lower variant type atoms.

	for ( Size n = 1; n <= 3; n++ ) { //alpha, beta, gamma of 3' res
		pose.set_torsion( TorsionID( three_prime_res, id::BB,  n ), 0.0 );
	}

	for ( Size n = 5; n <= 6; n++ ) { //epsilon and zeta of 5' res
		pose.set_torsion( TorsionID( five_prime_res, id::BB,  n ), 0.0 );
	}

}


/////////////////////////////////////////////////////////////////////////////////////////////
void
print_atom_info( pose::Pose const & pose, Size const seq_num, std::string const & pose_name ){
	TR << "print_atom_info for pose: " << pose_name << " seq_num = " << seq_num << std::endl;

	conformation::Residue const & rsd = pose.residue( seq_num ); //static_pose

	for ( Size at = 1; at <= rsd.natoms(); at++ ) { //I wonder if we should just consider heavy atom? (rsd_1.nheavyatoms())

		TR << "atom = " << at  << "|name = " << rsd.type().atom_name( at ) << "|type = " << rsd.atom_type( at ).name();
		TR << "|element() = " << rsd.atom_type( at ).element() << "|" << std::endl;

	}

}
/////////////////////////////////////////////////////////////////////////////////////////////
void
print_individual_atom_info( core::conformation::Residue const & rsd, Size const atomno, std::string const & rsd_name ){
	TR << "individual_atom_info: rsd_name " << rsd_name;

	TR << " atom = " << atomno  << "|name = " << rsd.type().atom_name( atomno ) << "|type = " << rsd.atom_type( atomno ).name();
	TR << "|element() = " << rsd.atom_type( atomno ).element() << "|" << std::endl;

}

/////////////////////////////////////////////////////////////////////////////////////////////

void
print_base_state( std::string const & tag, core::Size const base_state, std::ostream & outstream /* = std::cout */ ){

	std::string base_state_string = "";

	if ( base_state == NO_CHI ) {
		base_state_string = "NO_CHI";
	} else if ( base_state == ANTI ) {
		base_state_string = "ANTI";
	} else if ( base_state == SYN ) {
		base_state_string = "SYN" ;
	} else if ( base_state == ANY_CHI ) {
		base_state_string = "ANY_CHI";
	} else {
		outstream << "Invalid base state = " << base_state << std::endl;
		utility_exit_with_message( "Invalid base state!" );
	}

	outstream << tag << std::setw( 4 ) << std::left <<  base_state_string << " ";
}

/////////////////////////////////////////////////////////////////////////////////////////////

void
print_sugar_pucker_state( std::string const & tag, core::Size const pucker_state, std::ostream & outstream /* = std::cout */ ){

	std::string pucker_state_string = "";

	if ( pucker_state == NORTH ) {
		pucker_state_string = "NORTH";
	} else if ( pucker_state == SOUTH ) {
		pucker_state_string = "SOUTH";
	} else if ( pucker_state == ANY_PUCKER ) {
		pucker_state_string = "ANY_PUCKER";
	} else if ( pucker_state == NO_PUCKER ) {
		pucker_state_string = "NO_PUCKER";
	} else {
		outstream << "Invalid pucker state = " << pucker_state << std::endl;
		utility_exit_with_message( "Invalid pucker state!" );
	}

	outstream << tag << std::setw( 5 ) << std::left << pucker_state_string << " ";

}

/////////////////////////////////////////////////////////////////////////////
scoring::ScoreFunctionOP
get_modeler_scorefxn( scoring::ScoreFunctionCOP scorefxn ){

	using namespace core::scoring;

	ScoreFunctionOP modeler_scorefxn_ = scorefxn->clone();
	//  modeler_scorefxn_->set_weight( rna_sugar_close, 0.0 ); (TURN IT BACK ON: RD 01/31/2010)

	modeler_scorefxn_->set_weight( fa_rep, 0.12 );
	//Only important only if fa_rep score in weight file is not 0.12..want to make sure that changing fa_rep in weight file doesn't effect modeler process. May 23 2010 Parin S.

	modeler_scorefxn_->set_weight( linear_chainbreak, 0.0 );
	modeler_scorefxn_->set_weight( chainbreak, 0.0 );
	modeler_scorefxn_->set_weight( angle_constraint, 0.0 );
	modeler_scorefxn_->set_weight( atom_pair_constraint, 0.0 );
	//This makes sure that there are no chain_break score involve in the full_score screening.
	//This works since by the time a pose reach full_score screening, it must already pass chain_break screening, May 23, 2010 Parin S.

	return modeler_scorefxn_;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
void
copy_all_o2prime_torsions( core::pose::Pose & mod_pose, core::pose::Pose const & template_pose ){
	using namespace core::id;
	runtime_assert ( template_pose.total_residue() == mod_pose.total_residue() );

	for ( Size seq_num = 1; seq_num <= template_pose.total_residue(); seq_num++ ) {

		if ( !template_pose.residue( seq_num ).is_RNA() ) continue;

		if ( std::abs( template_pose.torsion( TorsionID( seq_num, id::CHI, 4 ) ) - mod_pose.torsion( TorsionID( seq_num, id::CHI, 4 ) ) ) < 0.001 ) continue;

		mod_pose.set_torsion( TorsionID( seq_num, id::CHI, 4 ), template_pose.torsion( TorsionID( seq_num, id::CHI, 4 ) ) );
	}
}


/////////////////////////////////////////////////////////////////////////////////////////////
///According to Kyle B., DFPMIN should be start enough to determine what is the ideal step size. However, the exception is the first minimizing step which could lead to "blow up" error!
///To prevent can create new scorefxn with  scaling_factor=0.1 and minimize with this new score function Sept 20, 2011. Parin S.
core::scoring::ScoreFunctionOP
rescale_scorefxn( core::scoring::ScoreFunctionOP const & starting_scorefxn, Real const scaling_factor ){

	using namespace core::scoring;

	core::scoring::ScoreFunctionOP rescaled_scorefxn = starting_scorefxn->clone();

	core::Size non_zero_weight = 0;

	show_scorefxn_weight_lines( rescaled_scorefxn, "BEFORE REWEIGHT" );

	for ( Size n = 1; n <= n_score_types; n++ ) {

		core::Real const old_weight = rescaled_scorefxn->get_weight( ScoreType( n ) );

		if ( old_weight != 0.0 ) {
			non_zero_weight += 1;
			rescaled_scorefxn->set_weight( ScoreType( n ), old_weight*scaling_factor );
		}
	}

	TR.Debug << std::endl;

	TR.Debug << "n_score_types = " << int( n_score_types ) << " non_zero_weight = " << non_zero_weight << " scaling_factor = " << scaling_factor << std::endl;

	show_scorefxn_weight_lines( rescaled_scorefxn, "AFTER REWEIGHT" );

	return rescaled_scorefxn;

}

/////////////////////////////////////////////////////////////////////////////////////////////

void
show_scorefxn_weight_lines( core::scoring::ScoreFunctionOP const & scorefxn, std::string const & title ){

	using namespace core::scoring;
	using namespace ObjexxFCL::format;

	TR.Debug << "----------------" << title << "----------------" << std::endl;

	TR.Debug << "----------------------------------------------\n";
	TR.Debug << " Scores                             Weight\n";
	TR.Debug << "----------------------------------------------\n";

	pose::Pose empty_pose = *( new pose::Pose );

	for ( Size n = 1; n <= n_score_types; n++ ) {

		core::Real const weight = scorefxn->get_weight( ScoreType( n ) );

		if ( weight != 0.0 ) {
			TR.Debug << ' ' << LJ( 30, ScoreType( n ) ) << ' ' << F( 11, 5, weight ) << std::endl;
		}
	}

	std::string dash_string = "";
	for ( Size n = 1; n <= title.size(); n++ ) {
		dash_string += "-";
	}

	TR.Debug << "----------------" << dash_string << "----------------" << std::endl;

}


////////////////////////////////////////////////////////////////////
void
choose_random_if_unspecified_nucleotide( char & newrestype ) {

	std::string const rna_chars = "acgu";
	if ( newrestype == 'n' ) {
		newrestype = rna_chars[ numeric::random::rg().random_range( 1, rna_chars.size() ) - 1 ];
		TR << "Choosing random nucleotide: " << newrestype << std::endl;
	}

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
mutate_res_if_allowed( pose::Pose & pose, Size const mutate_res, Real const mutation_frequency /* = 0.5 */ ){

	using namespace core::pose::full_model_info;
	using namespace protocols::farna;

	// first need to slice up native_pose to match residues in actual pose.
	// define atoms over which to compute RMSD, using rmsd_res.
	FullModelInfo const & full_model_info = const_full_model_info( pose );
	utility::vector1< Size > const sub_to_full = get_res_list_from_full_model_info( pose );
	std::string const full_sequence = full_model_info.full_sequence();

	if ( numeric::random::rg().uniform() < mutation_frequency ) {

		char nt = full_sequence[ sub_to_full[ mutate_res ] - 1 ];
		choose_random_if_unspecified_nucleotide( nt );

		char const nt_orig = pose.sequence()[ mutate_res - 1 ];
		if ( nt != nt_orig ) {
			pose::rna::mutate_position( pose, mutate_res, nt );
			return true;
		}
	}
	return false;

}

/////////////////////////////////////////////////////////////////////
std::string
create_tag( std::string const & prestring, Size const i ) {
	using namespace ObjexxFCL;
	std::string tag = prestring;
	tag.append( "_" + lead_zero_string_of( i, 9 ) );
	return tag;
}


////////////////////////////////////////////////////////////////////////
std::string //silly function to convert to real to string
create_torsion_value_string( core::Real const & torsion_value ) {

	using namespace ObjexxFCL;

	std::string torsion_string = "";

	core::Real const principal_torsion = numeric::principal_angle_degrees( torsion_value );

	Size const principal_torsion_SIZE = Size( std::abs( principal_torsion + 0.00001 ) ); //0.00001 is to prevent random ambiguity if the torsion decimal value is exactly .0000 Oct 12, 2010


	if ( principal_torsion > 0 ) {
		torsion_string = "p" + lead_zero_string_of( principal_torsion_SIZE, 3 );
	} else {
		torsion_string = "n" + lead_zero_string_of( principal_torsion_SIZE, 3 );
	}

	return torsion_string;
}

////////////////////////////////////////////////////////////////////////
std::string //silly function used for appending the rotamer value to the tag
create_rotamer_string( core::pose::Pose const & pose, Size const moving_res, bool const is_prepend ) {
	Size const reference_res  = ( is_prepend ) ? ( moving_res + 1 ) : (moving_res - 1);
	return create_rotamer_string( pose, moving_res, reference_res );
}

////////////////////////////////////////////////////////////////////////
std::string //silly function used for appending the rotamer value to the tag
create_rotamer_string( core::pose::Pose const & pose, Size const moving_res, Size const reference_res ) {
	std::string rotamer_tag = "";

	Size const five_prime_res  = ( moving_res < reference_res ) ? moving_res : reference_res;
	Size const three_prime_res = ( moving_res < reference_res ) ? reference_res : moving_res;
	runtime_assert( five_prime_res == three_prime_res - 1 );
	conformation::Residue const & five_prime_rsd  = pose.residue( five_prime_res );
	conformation::Residue const & three_prime_rsd = pose.residue( three_prime_res );

	rotamer_tag.append( "_E" + create_torsion_value_string( five_prime_rsd.mainchain_torsion( 5  ) ) );
	rotamer_tag.append( "_Z" + create_torsion_value_string( five_prime_rsd.mainchain_torsion( 6  ) ) );
	rotamer_tag.append( "_A" + create_torsion_value_string( three_prime_rsd.mainchain_torsion( 1 ) ) );
	rotamer_tag.append( "_B" + create_torsion_value_string( three_prime_rsd.mainchain_torsion( 2 ) ) );
	rotamer_tag.append( "_G" + create_torsion_value_string( three_prime_rsd.mainchain_torsion( 3 ) ) );


	if ( moving_res < reference_res ) {
		rotamer_tag.append( "_D" + create_torsion_value_string( five_prime_rsd.mainchain_torsion( 4 ) ) );
		rotamer_tag.append( "_C" + create_torsion_value_string( five_prime_rsd.chi(  1 ) ) );
	} else {
		rotamer_tag.append( "_D" + create_torsion_value_string( three_prime_rsd.mainchain_torsion( 4 ) ) );
		rotamer_tag.append( "_C" + create_torsion_value_string( three_prime_rsd.chi( 1 ) ) );
	}

	return rotamer_tag;

}

/////////////////////////////////////////////////////////////////////
void
add_fade_chain_break_constraint_across_gap( pose::Pose & pose,
	Size const five_prime_res,
	Size const three_prime_res,
	Size const gap_size ){

	using namespace core::conformation;
	using namespace core::scoring::constraints;
	using namespace core::id;

	runtime_assert( gap_size > 0 );
	if ( gap_size == GAP_SIZE_DUMMY ) return; // this signifies that the residues are actually on different strands.

	Distance min_dist( 0.0 ), max_dist( 0.0 );
	get_possible_O3prime_C5prime_distance_range( gap_size, min_dist, max_dist );

	ConstraintSetOP cst_set( pose.constraint_set()->clone() );
	assert( cst_set ); //if ( !cst_set ) cst_set = new ConstraintSet();

	Distance fade_zone( 2.0 );
	Real const well_depth( -10.0 ), well_offset( +10.0 );
	core::scoring::func::FuncOP const distance_func( new core::scoring::func::FadeFunc( min_dist, max_dist, fade_zone,
		well_depth, well_offset ) );

	Residue const & rsd1( pose.residue( five_prime_res ) );
	Residue const & rsd2( pose.residue( three_prime_res ) );
	AtomID const O3_id( rsd1.atom_index( "O3'" ), five_prime_res );
	AtomID const C5_id( rsd2.atom_index( "C5'" ), three_prime_res );

	// distance from O3' to C5'
	cst_set->add_constraint( ConstraintCOP( ConstraintOP( new AtomPairConstraint( O3_id, C5_id, distance_func ) ) ) );

	pose.constraint_set( cst_set );

}

/////////////////////////////////////////////////////////////////////
void
add_harmonic_chain_break_constraint( pose::Pose & pose, Size const five_prime_res ){

	using namespace core::conformation;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace core::id;

	using numeric::conversions::radians;

	Size three_prime_res = five_prime_res + 1;
	//    From RAD.param file
	//   ICOOR_INTERNAL  UPPER -175.907669   60.206192    1.607146   O3'   C3'   C4'  , Upper is P1
	//    ICOOR_INTERNAL  LOWER  -64.027359   71.027062    1.593103   P     O5'   C5'  , Lower is O3'
	//    Bug that bond distance is not the same. Rhiju suggest using 1.593103

	//Original (Amber?) parameter 1.608, 119.8, 103.4
	Real const O3_P_distance( 1.593 ); //amber=1.608
	Real const O3_angle( 119.8 ); // 180-60.206192
	Real const  P_angle( 108.97 ); // Quite off from original (original==Amber's parameter??),180-71.027062

	Real const distance_stddev( 0.0659 ); // amber is 0.0659
	Real const angle_stddev_degrees_P( 8.54 ); // amber is 8.54 (P angle), 5.73 (O3 angle) //did I made a mistake with the variable naming or is that a actual error? May 25, 2011.
	Real const angle_stddev_degrees_O3( 5.73 ); //did I made a mistake with the variable naming or is that a actual error? May 25, 2011.

	ConstraintSetOP cst_set( pose.constraint_set()->clone() );
	assert( cst_set ); //if ( !cst_set ) cst_set = new ConstraintSet();

	core::scoring::func::FuncOP const distance_func( new core::scoring::func::HarmonicFunc( O3_P_distance, distance_stddev ) );
	core::scoring::func::FuncOP const O3_angle_func( new core::scoring::func::HarmonicFunc( radians( O3_angle ), radians( angle_stddev_degrees_P ) ) ); //did I made a mistake with the variable naming or is that a actual error? May 25, 2011.
	core::scoring::func::FuncOP const  P_angle_func( new core::scoring::func::HarmonicFunc( radians(  P_angle ), radians( angle_stddev_degrees_O3 ) ) ); //did I made a mistake with the variable naming or is that a actual error? May 25, 2011.

	Residue const & rsd1( pose.residue( five_prime_res ) );
	Residue const & rsd2( pose.residue( three_prime_res ) );

	AtomID const C3_id( rsd1.atom_index( "C3'" ), five_prime_res );
	AtomID const O3_id( rsd1.atom_index( "O3'" ), five_prime_res );
	AtomID const  P_id( rsd2.atom_index( "P"   ), three_prime_res );
	AtomID const O5_id( rsd2.atom_index( "O5'" ), three_prime_res );

	// actually sometimes these constraints are turned on in an effort to get 'reasonable'
	//  geometries at virtual residue chainbreaks that don't have explicit cutpoint variants.
	//  runtime_assert( !pose.residue( C3_id.rsd() ).is_virtual( C3_id.atomno() ) );
	//  runtime_assert( !pose.residue( O3_id.rsd() ).is_virtual( O3_id.atomno() ) );
	//  runtime_assert( !pose.residue(  P_id.rsd() ).is_virtual(  P_id.atomno() ) );
	//  runtime_assert( !pose.residue( O5_id.rsd() ).is_virtual( O5_id.atomno() ) );

	// distance from O3' to P
	cst_set->add_constraint( ConstraintCOP( ConstraintOP( new AtomPairConstraint( O3_id, P_id, distance_func ) ) ) );

	// angle at O3'
	cst_set->add_constraint( ConstraintCOP( ConstraintOP( new AngleConstraint( C3_id, O3_id, P_id, O3_angle_func ) ) ) );

	// angle at P
	cst_set->add_constraint( ConstraintCOP( ConstraintOP( new AngleConstraint( O3_id, P_id, O5_id,  P_angle_func ) ) ) );

	pose.constraint_set( cst_set );
}


/////////////////////////////////////////////////////////////////////
void
get_possible_O3prime_C5prime_distance_range( Size const gap_size_, Distance & min_dist, Distance & max_dist ){
	min_dist = 0.0;
	if ( gap_size_ > 1 ) { //new option Sept 18, 2010, most useful for long loop mode...
		max_dist = O3I_C5I_PLUS_TWO_MAX_DIST + ( ( gap_size_ - 1 )*( C5I_C5I_PLUS_ONE_MAX_DIST ) );
	} else if ( gap_size_ == 1 ) {
		//previously used 11.0138 as O3I_C5I_PLUS_TWO_MAX_DIS which is slight underestimate...TOO strict;
		max_dist = O3I_C5I_PLUS_TWO_MAX_DIST;
	} else if ( gap_size_ == 0 ) {
		min_dist = O3I_C5I_MIN_DIST;
		max_dist = O3I_C5I_MAX_DIST;
	}
}

// used in ERRASER.
void
remove_all_virtual_phosphates( core::pose::Pose & pose ){
	for ( Size n = 1; n <= pose.total_residue(); n++ ) {
		pose::remove_variant_type_from_pose_residue( pose, core::chemical::VIRTUAL_PHOSPHATE, n );
	}
}

//////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
just_rna( utility::vector1< Size > const & res_list, pose::Pose const & pose ){
	utility::vector1< Size > rna_res_list;
	for ( Size n = 1; n <= res_list.size(); n++ ) { if ( pose.residue_type( res_list[n] ).is_RNA() ) rna_res_list.push_back( res_list[n] ); }
	return rna_res_list;
}

//////////////////////////////////////////////////////////////////////////////////
void
figure_out_moving_rna_chain_breaks( pose::Pose const & pose,
	utility::vector1< Size > const & moving_partition_res,
	utility::vector1< Size > & rna_cutpoints_closed,
	utility::vector1< Size > & rna_five_prime_chain_breaks,
	utility::vector1< Size > & rna_three_prime_chain_breaks,
	utility::vector1< Size > & rna_chain_break_gap_sizes ){

	utility::vector1< Size > cutpoints_closed;
	utility::vector1< Size > five_prime_chain_breaks;
	utility::vector1< Size > three_prime_chain_breaks;
	utility::vector1< Size > chain_break_gap_sizes;
	figure_out_moving_chain_breaks( pose, moving_partition_res, cutpoints_closed, five_prime_chain_breaks, three_prime_chain_breaks, chain_break_gap_sizes );

	rna_cutpoints_closed = just_rna( cutpoints_closed, pose );

	rna_five_prime_chain_breaks.clear();
	rna_three_prime_chain_breaks.clear();
	rna_chain_break_gap_sizes.clear();

	for ( Size n = 1; n <= chain_break_gap_sizes.size(); n++ ) {
		if ( !pose.residue_type( five_prime_chain_breaks[n] ).is_RNA() ) continue;
		runtime_assert( pose.residue_type( three_prime_chain_breaks[n] ).is_RNA() );
		rna_five_prime_chain_breaks.push_back( five_prime_chain_breaks[n] );
		rna_three_prime_chain_breaks.push_back( three_prime_chain_breaks[n] );
		rna_chain_break_gap_sizes.push_back( chain_break_gap_sizes[n] );
	}
}

//////////////////////////////////////////////////////////////////////
// this is reasonably sophisticated...
//
// * looks at base contacts -- at least 5 atoms in base contact some other residue.
// * also checks if 2'-OH is H-bonded, and in that case, keeps backbone instantiated (just virtualize bases).
// * checks if phosphate makes hydrogen bonds, and then keeps those instantiated.
//
// Good test case is second U in UUCG in 2KOC.pdb.
//
//  -- rhiju, 2014
//////////////////////////////////////////////////////////////////////
void
virtualize_free_rna_moieties( pose::Pose & pose ){

	utility::vector1< bool > base_makes_contact      = rna::bulge::detect_base_contacts( pose );
	utility::vector1< bool > sugar_makes_contact     = rna::sugar::detect_sugar_contacts( pose );
	utility::vector1< bool > phosphate_makes_contact = rna::phosphate::detect_phosphate_contacts( pose );

	for ( Size i = 1; i <= pose.total_residue(); i++ ) {
		if ( !pose.residue_type( i ).is_RNA() ) continue;
		if ( base_makes_contact[ i ] ) continue;

		// base is virtual.
		if ( sugar_makes_contact[ i ] ) {
			add_variant_type_to_pose_residue( pose, chemical::VIRTUAL_BASE, i );
			continue;
		}

		// base and sugar are virtual.
		add_variant_type_to_pose_residue( pose, chemical::VIRTUAL_RNA_RESIDUE, i );

		// phosphate 5' to base
		if ( i > 1 && phosphate_makes_contact[ i ] &&
				pose.residue_type( i-1 ).is_RNA() &&
				!pose.fold_tree().is_cutpoint(i-1) ) {
			rna::phosphate::setup_three_prime_phosphate_based_on_next_residue( pose, i-1 );
		}

		// phosphate 3' to base
		if ( i < pose.total_residue() && !phosphate_makes_contact[ i+1 ] &&
				pose.residue_type( i+1 ).is_RNA() &&
				!pose.fold_tree().is_cutpoint( i ) ) {
			add_variant_type_to_pose_residue( pose, chemical::VIRTUAL_PHOSPHATE, i+1 );
		}

	}

	// check if all virtual
	bool all_virtual( true );
	for ( Size i = 1; i <= pose.total_residue(); i++ ) {
		if ( !pose.residue( i ).is_virtual_residue() && !pose.residue( i ).has_variant_type( chemical::VIRTUAL_RNA_RESIDUE ) ) all_virtual = false;
	}
	if ( all_virtual ) utility_exit_with_message( "Turned native into an all virtual pose. May want to fix native or rerun with -virtualize_free_moieties_in_native false." );
	TR.Debug << pose.annotated_sequence() << std::endl;

}


/////////////////////////////////////////////////////////////////////////////
bool
just_modeling_RNA( std::string const & sequence ) {
	std::string const rna_letters( "acgunZ" );
	for ( Size k = 1; k <= sequence.size(); k++ ) {
		if ( rna_letters.find( sequence[k-1] ) == std::string::npos ) return false;
	}
	return true;
}

} //rna
} //modeler
} //stepwise
} //protocols
