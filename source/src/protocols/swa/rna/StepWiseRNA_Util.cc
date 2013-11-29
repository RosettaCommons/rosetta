// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseRNA_Util
/// @brief Util functions for Stepwise Assembly RNA.
/// @detailed
/// @author Parin Sripakdeevong


//////////////////////////////////
#include <protocols/swa/rna/StepWiseRNA_Classes.hh>
#include <protocols/swa/rna/StepWiseRNA_Util.hh>
#include <protocols/swa/rna/StepWiseRNA_ResidueInfo.hh>
#include <protocols/swa/rna/StepWiseRNA_JobParameters.hh> //Oct 22, 2011..for some reason Util.cc have JP.hh BUT JP.cc also have Util.hh!!! SHOULD RESOLVE THIS!

#include <protocols/rna/RNA_BasePairClassifier.hh>
#include <protocols/rna/RNA_ProtocolUtil.hh>
#include <core/scoring/rna/RNA_BaseDoubletClasses.hh>
#include <core/scoring/ScoreType.hh> //Parin Sept 20, 2011.
//////////////////////////////////

#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/util.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/FullModelInfoUtil.hh>
#include <core/pose/util.hh>
#include <core/pose/util.tmpl.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/BinaryRNASilentStruct.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/rotamer_trials.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/chemical/rna/RNA_FittedTorsionInfo.hh>
#include <core/chemical/rna/RNA_Util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/util.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/optimization/AtomTreeMinimizer.hh>

#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/TorsionID.hh>

#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>

#include <numeric/conversions.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>
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

static numeric::random::RandomGenerator RG(257572);  // <- Magic number, do not change it!

static basic::Tracer TR( "protocols.swa.rna.StepWiseRNA_Util" );

namespace protocols {
namespace swa {
namespace rna {


	bool is_OP2_atom( std::string const & atom_name ){			return ( atom_name == " OP2" ); }

	bool is_OP1_atom( std::string const & atom_name ){			return ( atom_name == " OP1" ); }

	bool is_P_atom( std::string const & atom_name ){				return ( atom_name == " P  " ); }

	bool is_O2prime_atom( std::string const & atom_name ){  return ( atom_name == " O2'" ); }

	bool is_O3prime_atom( std::string const & atom_name ){  return ( atom_name == " O3'" ); }

	bool is_O4prime_atom( std::string const & atom_name ){  return ( atom_name == " O4'" ); }

	bool is_O5prime_atom( std::string const & atom_name ){	return ( atom_name == " O5'" ); }

	bool is_C2prime_atom( std::string const & atom_name ){  return ( atom_name == " C2'" ); }

	bool is_C3prime_atom( std::string const & atom_name ){  return ( atom_name == " C3'" ); }

	bool is_C4prime_atom( std::string const & atom_name ){  return ( atom_name == " C4'" ); }

	bool is_C5prime_atom( std::string const & atom_name ){  return ( atom_name == " C5'" ); }

	bool is_1H5prime_atom( std::string const & atom_name ){ return ( atom_name == " H5'" ); }

	bool is_2H5prime_atom( std::string const & atom_name ){ return ( atom_name == "H5''" ); }

	bool is_H3prime_atom( std::string const & atom_name ){  return ( atom_name == " H3'" ); }

	bool is_H4prime_atom( std::string const & atom_name ){  return ( atom_name == " H4'" ); }

	bool is_three_prime_phosphate_atom( std::string const & atom_name ){ return ( is_O3prime_atom( atom_name ) ); }

	bool is_five_prime_phosphate_atom( std::string const & atom_name ){ return ( is_O5prime_atom( atom_name ) || is_OP2_atom( atom_name ) || is_OP1_atom( atom_name ) || is_P_atom( atom_name ) ) ; }

	bool is_phosphate_atom( std::string const & atom_name ){ return ( is_three_prime_phosphate_atom( atom_name ) || is_five_prime_phosphate_atom( atom_name ) ); }


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

	bool
	check_can_prepend( utility::vector1< core::Size > const & seq_num_list ){

		for ( Size n = 1; n <= seq_num_list.size() - 1; n++ ){ //[11, 12, 13]
			if ( ( seq_num_list[n] + 1 ) != seq_num_list[n + 1] ) return false;
		}
		return true;
	}

	bool
	check_can_append( utility::vector1< core::Size > const & seq_num_list ){

		for ( Size n = 1; n <= seq_num_list.size() - 1; n++ ){ //[14, 13, 12]
			if ( ( seq_num_list[n] - 1 ) != seq_num_list[n + 1] ) return false;
		}
		return true;
	}

	////////////////////////////////////////////////May 04, 2011////////////////////////////////////////////////////
	void
	apply_protonated_H1_adenosine_variant_type( core::pose::Pose & pose, core::Size const & seq_num, bool const apply_check ){

		bool verbose = true;

		if ( verbose ) TR << "Applying PROTONATED_H1_ADENOSINE variant_type to seq_num " << seq_num << std::endl;

		if ( apply_check ){
			//Basically the two variant type are not compatible, VIRTUAL_RNA_RESIDUE variant type currently does not virtualize the protonated H1 atom.
			if ( pose.residue( seq_num ).has_variant_type( "VIRTUAL_RNA_RESIDUE" ) ){
				utility_exit_with_message( "Cannot apply PROTONATED_H1_ADENOSINE variant_type to seq_num: " + ObjexxFCL::string_of( seq_num ) + ". This residue have a incompatible VIRTUAL_RNA_RESIDUE variant type."  );
			}
		}


		if ( pose.residue( seq_num ).has_variant_type( "PROTONATED_H1_ADENOSINE" ) ){
			TR << "WARNING pose already have PROTONATED_H1_ADENOSINE variant_type at seq_num = " << seq_num << ", early RETURN!" << std::endl;
			return;
			//utility_exit_with_message("pose already have PROTONATED_H1_ADENOSINE variant_type at seq_num= " + ObjexxFCL::string_of(seq_num));
		}

		if ( pose.total_residue() < seq_num ){
			utility_exit_with_message(  "Cannot apply PROTONATED_H1_ADENOSINE variant_type to seq_num: " + ObjexxFCL::string_of( seq_num ) + ". pose.total_residue() < seq_num"  );
		}

		if ( pose.residue( seq_num ).aa() != core::chemical::na_rad ){
			utility_exit_with_message( "working_seq_num = " + ObjexxFCL::string_of( seq_num ) + " cannot have PROTONATED_H1_ADENOSINE variant type since it is not a adenosine!" );
		}

		pose::add_variant_type_to_pose_residue( pose, "PROTONATED_H1_ADENOSINE", seq_num );

	}
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	void
	apply_virtual_rna_residue_variant_type( core::pose::Pose & pose, core::Size const & seq_num, bool const apply_check ){

		utility::vector1< Size > working_cutpoint_closed_list;
		working_cutpoint_closed_list.clear(); //empty list

		apply_virtual_rna_residue_variant_type( pose, seq_num, working_cutpoint_closed_list, apply_check );
	}

	void
	apply_virtual_rna_residue_variant_type( core::pose::Pose & pose, core::Size const & seq_num, utility::vector1< core::Size > const & working_cutpoint_closed_list, bool const apply_check ){

		using namespace core::chemical;
		using namespace ObjexxFCL;

		if ( pose.total_residue() < seq_num ){
			utility_exit_with_message(  "Cannot apply VIRTUAL_RNA_RESIDUE VARIANT TYPE to seq_num: " + string_of( seq_num ) + ". pose.total_residue() < seq_num"  );
		}

		if ( pose.residue_type( seq_num ).has_variant_type( "VIRTUAL_RNA_RESIDUE" ) ) return;

		//Basically the two variant type are not compatible, VIRTUAL_RNA_RESIDUE variant type currently does not virtualize the protonated H1 atom.
		if ( pose.residue( seq_num ).has_variant_type( "PROTONATED_H1_ADENOSINE" ) ){
			TR << "Removing PROTONATED_H1_ADENOSINE variant_type from seq_num = " << seq_num << " before adding VIRTUAL_RNA_RESIDUE variant_type since the two variant_types are not compatible!" << std::endl;
			pose::remove_variant_type_from_pose_residue( pose, "PROTONATED_H1_ADENOSINE", seq_num );
		}
		//OK PROTONATED_H1_ADENOSINE variant type should also be removed when adding VIRTUAL_RNA_RESIDUE_EXCLUDE_PHOSPHATE variant type or BULGE variant type.
		//However these two variant type are not currently used in standard SWA run (May 04, 2011)

		bool is_cutpoint_closed = false;

		if ( pose.residue( seq_num ).has_variant_type( chemical::CUTPOINT_LOWER ) ){
			if ( pose.residue( seq_num + 1 ).has_variant_type( chemical::CUTPOINT_UPPER ) == false ){
				utility_exit_with_message( "seq_num " + string_of( seq_num ) + " is a CUTPOINT_LOWER but seq_num " + string_of( seq_num + 1 ) + " is not a cutpoint CUTPOINT_UPPER??" );
			}
			is_cutpoint_closed = true;
		}

		//Ok another possibility is that the CUTPOINT_LOWER and CUTPOINT_UPPER variant type had not been applied yet..so check the working_cutpoint_closed_list
		for ( Size n = 1; n <= working_cutpoint_closed_list.size(); n++ ){
			if ( seq_num == working_cutpoint_closed_list[n] ) {
				is_cutpoint_closed = true;
				break;
			}
		}


		bool const is_cutpoint_open = ( pose.fold_tree().is_cutpoint( seq_num ) && is_cutpoint_closed == false );

		if ( apply_check ){
			if ( is_cutpoint_open ) {
				utility_exit_with_message( "Cannot apply VIRTUAL_RNA_RESIDUE VARIANT TYPE to seq_num: " + string_of( seq_num ) + ". The residue is 5' of a OPEN cutpoint" );
			}

			if ( pose.total_residue() == seq_num ) {
				utility_exit_with_message( "Cannot apply VIRTUAL_RNA_RESIDUE VARIANT TYPE to seq_num: " + string_of( seq_num ) + ". pose.total_residue() == seq_num" );
			}

		}

		pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_RNA_RESIDUE", seq_num );

		//if( !(pose.total_residue()==seq_num) &&  (is_cutpoint_open==false) ) {
		if ( ( pose.total_residue() != seq_num ) &&  ( is_cutpoint_open == false ) ) { //April 6, 2011
			pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_RNA_RESIDUE_UPPER", seq_num + 1 );
		}
	}

	//////////////////////////////////////////////////////////////////////////////////////
	void
	remove_virtual_rna_residue_variant_type( pose::Pose & pose, Size const & seq_num ){

		using namespace core::chemical;
		using namespace ObjexxFCL;

		if ( pose.total_residue() < seq_num ){
			utility_exit_with_message(  "Cannot remove VIRTUAL_RNA_RESIDUE VARIANT TYPE to seq_num: " + string_of( seq_num ) + ". pose.total_residue() < seq_num"  );
		}

		if ( pose.total_residue() == seq_num ){ //April 6, 2011
			utility_exit_with_message(  "Cannot remove VIRTUAL_RNA_RESIDUE VARIANT TYPE to seq_num: " + string_of( seq_num ) + ". pose.total_residue() == seq_num"  );
		}

		pose::remove_variant_type_from_pose_residue( pose, "VIRTUAL_RNA_RESIDUE", seq_num );
		pose::remove_variant_type_from_pose_residue( pose, "VIRTUAL_RNA_RESIDUE_UPPER", seq_num + 1 );

	}
	//////////////////////////////////////////////////////////////////////////////////////
	bool
	has_virtual_rna_residue_variant_type( pose::Pose & pose, Size const & seq_num ){

		using namespace ObjexxFCL;

		if ( pose.residue( seq_num ).has_variant_type( "VIRTUAL_RNA_RESIDUE" ) == false ) return false;

		if ( ( seq_num + 1 ) > pose.total_residue() ){ //Check in range
			TR << "( seq_num + 1 ) = " << ( seq_num + 1 )  << std::endl;
			utility_exit_with_message( "( seq_num + 1 ) > pose.total_residue()!" );
		}

		if ( pose.residue( seq_num + 1 ).has_variant_type( "VIRTUAL_RNA_RESIDUE_UPPER" ) == false ){
			TR << "Problem seq_num = " << seq_num << std::endl;
			utility_exit_with_message( "res ( " + string_of( seq_num ) + " ) has_variant_type VIRTUAL_RNA_RESIDUE but res seq_num + 1 ( " + string_of( seq_num + 1 ) + " )does not have variant_type VIRTUAL_RNA_RESIDUE_UPPER" );
		}

		return true;

	}

	//////////////////////////////////////////////////////////////////////////////////////
	void
	remove_all_variant_types( pose::Pose & pose ){


		for ( Size seq_num = 1; seq_num <= pose.total_residue(); seq_num++ ) {

			if ( pose.residue( seq_num ).aa() == core::chemical::aa_vrt ) continue; //Fang's electron density code

			utility::vector1< core::chemical::VariantType > target_variants( pose.residue( seq_num ).type().variant_types() );

			if ( target_variants.size() != pose.residue( seq_num ).type().variant_types().size() ){
				utility_exit_with_message( "target_variants.size() != pose.residue( seq_num ).type().variant_types().size()" );
			}

			Size skip_variant_count = 0;

			for ( Size i = 1; i <= target_variants.size(); i++ ) {
				//TR << "seq_num=" << seq_num << " variant_type[" << i << "]=" << target_variants[i] << std::endl;

				if ( pose.residue( seq_num ).type().has_variant_type( target_variants[i] ) == false ) utility_exit_with_message( "pose.residue( seq_num ).type().has_variant_type( target_variants[i] ) == false!" );

				bool skip_this_variant = false;

				if ( target_variants[i] == "LOWER_TERMINUS" ) skip_this_variant = true;

				if ( target_variants[i] == "UPPER_TERMINUS" ) skip_this_variant = true;

				if ( skip_this_variant ){
					skip_variant_count++;
					continue;
				}

				pose::remove_variant_type_from_pose_residue( pose, target_variants[i], seq_num );

			}

			if ( pose.residue( seq_num ).type().variant_types().size() != skip_variant_count ){
				TR << "pose.residue( seq_num ).type().variant_types().size() = " << pose.residue( seq_num ).type().variant_types().size() << std::endl;
				TR << "skip_variant_count = " << skip_variant_count << std::endl;
				utility_exit_with_message( "pose.residue( seq_num ).type().variant_types().size() != skip_variant_count" );
			}

		}

	}

	//////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< Size >
	apply_full_to_sub_mapping( utility::vector1< Size > const & res_vector, utility::vector1< core::Size > const & is_working_res, std::map< core::Size, core::Size > const & full_to_sub ){

		using namespace ObjexxFCL;

		if ( is_working_res.size() == 0 ){
			utility_exit_with_message( "is_working_res.size() == 0" );
		}

		if ( full_to_sub.empty() == true ){
			utility_exit_with_message( "full_to_sub.empty() == true" );
		}

		Size const total_res = is_working_res.size();

		utility::vector1< core::Size > working_res_vector;
		for ( Size n = 1; n <= res_vector.size(); n++ ) {

			if ( res_vector[ n ] > total_res ) utility_exit_with_message( "res_vector[ n ] ( " + string_of( res_vector[ n ] ) + " ) > total_res ( " + string_of( total_res ) + " )!" );

			if ( !is_working_res[ res_vector[ n ] ] ) continue;
			working_res_vector.push_back( full_to_sub.find( res_vector[ n ] )->second );
		}

		return working_res_vector;

	}

	//////////////////////////////////////////////////////////////////////////////////////

	utility::vector1< Size >
	apply_full_to_sub_mapping( utility::vector1< Size > const & res_vector, StepWiseRNA_JobParametersCOP job_parameters ){

		utility::vector1< core::Size > const & is_working_res = job_parameters->is_working_res();
		std::map< core::Size, core::Size > const & full_to_sub = job_parameters->const_full_to_sub();

		return apply_full_to_sub_mapping( res_vector, is_working_res, full_to_sub );

	}


	//////////////////////////////////////////////////////////////////////////////////////

	utility::vector1< Size >
	apply_sub_to_full_mapping( utility::vector1< Size > const & working_res_vector, StepWiseRNA_JobParametersCOP job_parameters ){

		std::map< core::Size, core::Size > const & sub_to_full( job_parameters->const_sub_to_full() );

		utility::vector1< core::Size > full_res_vector;
		for ( Size n = 1; n <= working_res_vector.size(); n++ ){

			if ( sub_to_full.find( working_res_vector[ n ] ) == sub_to_full.end() ){
				utility_exit_with_message( "sub_to_full.find( working_res_vector[ n ] ).end() == sub_to_full.end()!" );
			}

			full_res_vector.push_back( sub_to_full.find( working_res_vector[ n ] )->second );
		}

		return full_res_vector;
	}


	///////////////////////////This should be a function of the job_parameters class///////////////////////
	void
	ensure_valid_full_seq_num( Size const full_seq_num, StepWiseRNA_JobParametersCOP const & job_parameters ){

		using namespace ObjexxFCL;

		utility::vector1< core::Size > const & is_working_res = job_parameters->is_working_res();

		if ( full_seq_num < 1 ) utility_exit_with_message( "full_seq_num ( " + string_of( full_seq_num ) + " ) is lesser then 1" );

		if ( full_seq_num > is_working_res.size() ) utility_exit_with_message( "full_seq_num ( " + string_of( full_seq_num ) + " ) is greater than is_working_res.size() ( " + string_of( is_working_res.size() ) + " )" );

	}

	//////////////////////////This should be a function of the job_parameters class/////////////////////////
	bool
	check_is_working_res( Size const full_seq_num, StepWiseRNA_JobParametersCOP const & job_parameters ){

		using namespace ObjexxFCL;

		utility::vector1< core::Size > const & is_working_res = job_parameters->is_working_res();

		ensure_valid_full_seq_num( full_seq_num, job_parameters );

		return is_working_res[full_seq_num];

	}
	//////////////////////////This should be a function of the job_parameters class/////////////////////////
	core::Size
	check_validity_and_get_working_res( Size const full_seq_num, StepWiseRNA_JobParametersCOP const & job_parameters ){

		using namespace ObjexxFCL;

		std::map< core::Size, core::Size > const & full_to_sub = job_parameters->const_full_to_sub();

		std::string const & working_sequence = job_parameters->working_sequence();

		if ( check_is_working_res( full_seq_num, job_parameters ) == false ){
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

		for ( Size n = 1; n <= input_res_vector.size(); n++ ){
			full_to_input_res_map[input_res_vector[n]] = n;
		}

		return full_to_input_res_map;
	}
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/*
	core::Size
	string_to_int( std::string const string ){
		Size int_of_string; //misnomer
		std::stringstream ss ( std::stringstream::in | std::stringstream::out );
		ss << string;
		ss >> int_of_string;

	//	TR << "The string  " <<  string << " have the corresponding int value " << int_of_string << std::endl;

		return int_of_string;
	} */

	core::Size
	string_to_int( std::string const input_string ){

		Size int_of_string; //misnomer
		std::stringstream ss ( std::stringstream::in | std::stringstream::out );

		ss << input_string;

		if ( ss.fail() ) utility_exit_with_message( "In string_to_real(): ss.fail() for ss << input_string | string ( " + input_string + " )" );

		ss >> int_of_string;

		if ( ss.fail() ) utility_exit_with_message( "In string_to_real(): ss.fail() for ss >> int_of_string | string ( " + input_string + " )" );

		return int_of_string;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/*
	core::Real
	string_to_real( std::string const string ){
		Real real_of_string;
		std::stringstream ss ( std::stringstream::in | std::stringstream::out );
		ss << string;
		ss >> real_of_string;
		return real_of_string;

		//TR << "The string  " <<  string << " have the corresponding real value " << real_of_string << std::endl;
	}
	*/

	core::Real
	string_to_real( std::string const input_string ){

		Real real_of_string;
		std::stringstream ss ( std::stringstream::in | std::stringstream::out );

		ss << input_string;

		if ( ss.fail() ) utility_exit_with_message( "In string_to_real(): ss.fail() for ss << input_string | string ( " + input_string + " )" );

		ss >> real_of_string;

		if ( ss.fail() ) utility_exit_with_message( "In string_to_real(): ss.fail() for ss >> real_of_string | string ( " + input_string + " )" );

		return real_of_string;

	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	utility::vector1< std::string >
	tokenize( std::string const str, std::string delimiters ){
	  using namespace std;

		utility::vector1< std::string > tokens;

    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of( delimiters, 0 );
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of( delimiters, lastPos );

    while ( string::npos != pos || string::npos != lastPos ){
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

		//SML PHENIX conference
		if ( basic::options::option[basic::options::OptionKeys::rna::rna_prot_erraser].value() ){
			if ( !rsd.is_RNA() ) return false;
		}

		//Cytosine and Uracil contain 8 heavy base atoms. Adenine contains 10 heavy base atoms. Guanine contains 11 heavy base atoms.
		if ( ( rsd.nheavyatoms() - rsd.first_sidechain_atom() + 2 ) < 8 ){ //plus 2 since need to count both start and end atom.
			utility_exit_with_message( "The rna base " + name_from_aa( rsd.aa() ) + " contain lesser than 8 heavy atoms" );
		}

		Size non_virtual_atom_count = 0;
  	for ( Size atomno = rsd.first_sidechain_atom() + 1; atomno <= rsd.nheavyatoms(); ++atomno ) { //iterate over base atoms....+1 to exclude the O2prime oxygen
			if ( rsd.atom_type( atomno ).name() != "VIRT" ){
				non_virtual_atom_count++;
			}
		}

		bool const method_1 = ( non_virtual_atom_count < 8 ) ? true : false;
		bool const method_2 = ( rsd.has_variant_type( "VIRTUAL_RNA_RESIDUE" ) || rsd.has_variant_type( "BULGE" ) ) ? true : false;

		if ( method_1 != method_2 ){
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

 		if ( name_from_aa( rsd_1.aa() ) != name_from_aa( rsd_2.aa() ) ){
 			utility_exit_with_message( "rsd_1.aa() != rsd_2.aa(). res_num_1 = " + string_of( res_num_1 ) + " name_from_aa( rsd_1.aa() ) = " + name_from_aa( rsd_1.aa() ) + " res_num_2 = " +  string_of( res_num_2 ) + " name_from_aa( rsd_2.aa() ) = " + name_from_aa( rsd_2.aa() ) );
 		}

		Size const first_atom = ( base_only ) ? rsd_1.first_sidechain_atom() + 1 : 1; //+1 to exclude the O2prime oxygen

	  	for ( Size atomno_1 = first_atom; atomno_1 <= rsd_1.nheavyatoms(); ++atomno_1 ) {

 			std::string const atom_name_1 = rsd_1.type().atom_name( atomno_1 );

			if ( !rsd_2.has( atom_name_1 ) ) continue;

			Size const atomno_2 = rsd_2.atom_index( atom_name_1 );

 			//Check
 			std::string const atom_name_2 = rsd_2.type().atom_name( atomno_2 );
 			if ( atom_name_1 != atom_name_2 ){
 				utility_exit_with_message( "atom_name_1 != atom_name_2, atom_name_1 = " + atom_name_1 + " atom_name_2 = " + atom_name_2 );
 			}

 			if ( rsd_1.atom_type( atomno_1 ).name() == "VIRT" ) continue; //Check for virtual atoms
 			if ( rsd_2.atom_type( atomno_2 ).name() == "VIRT" ) continue; //Check for virtual atoms

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
 	create_alignment_id_map(	pose::Pose & mod_pose, pose::Pose const & ref_pose, utility::vector1< core::Size > const & rmsd_residue_list, bool const base_only ){
 		using namespace chemical;

 		id::AtomID_Map < id::AtomID > atom_ID_map;

		pose::initialize_atomid_map( atom_ID_map, mod_pose, id::BOGUS_ATOM_ID );

 		if ( ref_pose.sequence() != mod_pose.sequence() ){
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

	void
	align_poses( core::pose::Pose & moving_pose, std::string const moving_tag, core::pose::Pose const & static_pose, std::string const static_tag, utility::vector1< core::Size > const & working_best_alignment, bool const base_only ){

//		using namespace core::conformation;

		bool found_non_virtual_base = false;
		for ( Size n = 1; n <= working_best_alignment.size(); n++ ){
			Size const seq_num = working_best_alignment[n];
			if ( is_virtual_base( moving_pose.residue( seq_num ) ) == true || is_virtual_base( static_pose.residue( seq_num ) ) == true ) continue;

			found_non_virtual_base = true; //ok found a non-virtual base nucleotide that can be used for alignment
			break;
		}

		if ( found_non_virtual_base == false ){
			for ( Size n = 1; n <= working_best_alignment.size(); n++ ){
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
		id::AtomID_Map < id::AtomID > const & alignment_atom_id_map = create_alignment_id_map( moving_pose, static_pose, working_best_alignment, base_only );
		core::scoring::superimpose_pose( moving_pose, static_pose, alignment_atom_id_map );

//				current_pose.dump_pdb( tag+ "_current_pose_after_alignment");
		if ( check_for_messed_up_structure( moving_pose, moving_tag ) == true ){
			std::string error_message = "Error in aligning " + moving_tag + " to " + static_tag + "!";
			TR << error_message << std::endl;
			utility_exit_with_message( moving_tag + " is messed up ...this is probably an alignment problem" );
		};
	}

	void
	output_pair_size( std::pair < Size, Size > const & pair_size, std::ostream & outstream /* = std::cout */ ){
		outstream << "( " << pair_size.first << ", " << pair_size.second << " ) ";
	}

	void
	output_pair_size( utility::vector1 < std::pair < Size, Size > > const & pair_size_vector, std::string const & output_string, std::ostream & outstream /* = std::cout */, core::Size const spacing ){
		outstream << std::setw( spacing ) << std::left << output_string << " :";
		for ( Size n = 1; n <= pair_size_vector.size(); n++ ){
			output_pair_size( pair_size_vector[n], outstream );
		}
		outstream << std::endl;
	}

	//Sort by the first element. Low number on the top of the list
	bool
	pair_sort_citeria( std::pair < Size, Size > pair_one, std::pair < Size, Size > pair_two ){
		return ( pair_one.first < pair_two.first );
	}

	void
	sort_pair_list( utility::vector1< std::pair < Size, Size > > pair_list ){
		sort( pair_list.begin(), pair_list.end(), pair_sort_citeria );
	}


	bool
	seq_num_sort_citeria( core::Size seq_num_1, core::Size seq_num_2 ){
		return ( seq_num_1 < seq_num_2 );
	}

	void
	sort_seq_num_list( utility::vector1< core::Size > & seq_num_list ) { 	//Low number on the top of the list
		sort( seq_num_list.begin(), seq_num_list.end(), seq_num_sort_citeria );
	}
	void
	output_seq_num_list( std::string const tag, utility::vector1< core::Size > const & seq_num_list, std::ostream & outstream /* = std::cout */, core::Size const spacing ){

		using namespace ObjexxFCL;
		using namespace ObjexxFCL::format;

		outstream <<  std::setw( spacing ) << tag;

		utility::vector1< core::Size > sorted_seq_num_list = seq_num_list;
		sort_seq_num_list( sorted_seq_num_list );

		Size seq_num = 1;
		for ( Size n = 1; n <= sorted_seq_num_list.size(); n++ ){

			while ( seq_num < sorted_seq_num_list[n] ){
				outstream << A( 4, " " );
				seq_num++;
			}
			outstream << I( 4, sorted_seq_num_list[n] );
			seq_num++;
		}

		outstream << std::endl;

	}

	bool
	is_equivalent_vector( utility::vector1< core::Size > const & seq_num_list_1, utility::vector1< core::Size > const & seq_num_list_2 ){

		utility::vector1< core::Size > sorted_seq_num_list_1 = seq_num_list_1;
		utility::vector1< core::Size > sorted_seq_num_list_2 = seq_num_list_2;

		sort_seq_num_list( sorted_seq_num_list_1 );
		sort_seq_num_list( sorted_seq_num_list_2 );

		if ( sorted_seq_num_list_1.size() != sorted_seq_num_list_2.size() ) return false;

		for ( Size n = 1; n <= sorted_seq_num_list_1.size(); n++ ){
			if ( sorted_seq_num_list_1[n] != sorted_seq_num_list_2[n] ) return false;
		}

		return true;

	}

	void output_is_prepend_map( std::string const tag, std::map< core::Size, bool > const & my_map, core::Size const max_seq_num, std::ostream & outstream /* = std::cout */, core::Size const tag_spacing ){

		using namespace ObjexxFCL;
		using namespace ObjexxFCL::format;

		outstream << std::setw( tag_spacing ) << tag;

		Size spacing = 4;
//		outstream << std::setw(30) << "is_residue_prepend:";
		for ( Size seq_num = 1; seq_num <= max_seq_num; seq_num++ ){
			char prepend_char;
			if ( my_map.find( seq_num ) != my_map.end() ){
				prepend_char = ( my_map.find( seq_num )->second ) ? 'P' : 'A';
			} else{
				prepend_char = '-';
			}
			outstream << std::setw( spacing ) << prepend_char;
		}
		outstream << std::endl;

	}

	void
	output_bool_list( std::string const tag, utility::vector1< Size > const & size_list, std::ostream & outstream /* = std::cout */, core::Size const spacing ){
		utility::vector1< bool > bool_list;

		for ( Size n = 1; n <= size_list.size(); n++ ){
			bool_list.push_back( size_list[n] );
		}
		output_bool_list( tag, bool_list, outstream, spacing );
	}

	void
	output_bool_list( std::string const tag, utility::vector1< bool > const & bool_list, std::ostream & outstream /* = std::cout */, core::Size const spacing ){

		using namespace ObjexxFCL;
		using namespace ObjexxFCL::format;

		outstream <<  std::setw( spacing ) << tag;

		for ( Size seq_num = 1; seq_num <= bool_list.size(); seq_num++ ){
			output_boolean( bool_list[seq_num], outstream );
		}
		outstream << std::endl;
	}

	void
	output_size_list( std::string const tag, utility::vector1< Size > const & size_list, std::ostream & outstream /* = std::cout */, core::Size const spacing ){

		using namespace ObjexxFCL;
		using namespace ObjexxFCL::format;

		outstream <<  std::setw( spacing ) << tag;

		for ( Size seq_num = 1; seq_num <= size_list.size(); seq_num++ ){
			outstream << I( 4, size_list[seq_num] );
		}
		outstream << std::endl;
	}


	bool
	is_close_chain_break( pose::Pose const & pose ){

		for ( Size seq_num = 1; seq_num <= pose.total_residue(); seq_num++ ) {
			if ( pose.residue( seq_num ).aa() == core::chemical::aa_vrt ) continue; //Fang's electron density code.
			if ( !pose.residue( seq_num  ).has_variant_type( chemical::CUTPOINT_LOWER )  ) continue;
			if ( !pose.residue( seq_num + 1 ).has_variant_type( chemical::CUTPOINT_UPPER )  ) continue;

			return true;
		}
		return false;
	}

/*
	Size
	get_five_prime_chain_break( pose::Pose const & pose ){

		if ( !is_close_chain_break( pose ) ){
			outstream << "In StepWiseRNA_ResidueSampler::get_five_prime_chain_break, the function should not be called! " << std::endl;
			exit( 1 );
		};

		for ( Size seq_num = 1; seq_num <= pose.total_residue(); seq_num++ ){

			if ( !pose.residue( seq_num   ).has_variant_type( chemical::CUTPOINT_LOWER )  ) continue;
			if ( !pose.residue( seq_num + 1 ).has_variant_type( chemical::CUTPOINT_UPPER )  ) continue;

			return seq_num;
		}

		outstream << "In StepWiseRNA_ResidueSampler::get_five_prime_chain_break, cannot find five_prime_chain_break" << std::endl;
		exit( 1 );
	}
*/


	///////////////////////////////////////////////////////////////////////////////////////////////////
	void
	output_title_text( std::string const title, std::ostream & outstream /* = std::cout */ ){

		outstream << std::endl;

		Size title_length = title.size();
		Size char_per_line = 80;
		Size dash_length( 0 );
		if ( title_length < char_per_line ) dash_length = char_per_line - title_length;

		for ( Size i = 1; i <= dash_length/2; i++ ) 	outstream << "-";

		outstream << title;

		for ( Size i = 1; i <= dash_length/2; i++ ) outstream << "-";

		outstream << std::endl;

	}

	//////////////////////////////////////////////////////////////////////////////////////////
	void
	output_fold_tree_info( kinematics::FoldTree const & fold_tree, std::string const pose_name, std::ostream & outstream /* = std::cout */ ){

		outstream << "fold tree of " << pose_name << ": " << std::endl;
		for ( int i = 1; i <= fold_tree.num_cutpoint(); i++ ){
			outstream << std::setw( 30 ) << "jump_point_num = " << i;
			outstream << "   cutpoint = " << fold_tree.cutpoint( i );
			outstream << "   5' jump_point = " << fold_tree.jump_point( 1, i ) << ", " << fold_tree.upstream_atom( i );
			outstream << "   3' jump_point = " << fold_tree.jump_point( 2, i ) << ", " << fold_tree.downstream_atom( i ) << std::endl;
		}
	}

	void
	output_fold_tree_info( pose::Pose const & pose, std::string pose_name, std::ostream & outstream /* = std::cout */ ){
		output_fold_tree_info( pose.fold_tree(), pose_name, outstream );
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

		if ( file_exists( file_name ) == false ){
			utility_exit_with_message( "file_name ( " + file_name + " ) doesn't exist!" );
		}

		int const retcode = std::remove( file_name.c_str() );

		if ( retcode != 0 ){
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
	add_virtual_O2Star_hydrogen( core::pose::Pose & pose ){

	  for ( core::Size i = 1; i <= pose.total_residue(); i++ ){
	    if ( pose.residue( i ).aa() == core::chemical::aa_vrt ) continue; //Fang's electron density code
			if ( !pose.residue( i ).is_RNA() ) continue;
			pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_O2PRIME_HYDROGEN", i );
	  }
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	remove_virtual_O2Star_hydrogen( pose::Pose & pose ){

		for ( Size i = 1; i <= pose.total_residue(); i++ ){
			if ( pose.residue( i ).aa() == core::chemical::aa_vrt ) continue; //Fang's electron density code
			if ( !pose.residue( i ).is_RNA() ) continue;
			if ( pose.residue_type( i ).has_variant_type( "VIRTUAL_O2PRIME_HYDROGEN" ) ){
				pose::remove_variant_type_from_pose_residue( pose, "VIRTUAL_O2PRIME_HYDROGEN", i );
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

		if ( atom_count == 0 ) rmsd = 0.0; //special case...implement this on June_11, 2010, took me a whole day to debug this since buggy only on Biox compiler!

		return ( std::max( 0.01, rmsd ) );

	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Ok this function should only be called if pose contain full sequence.
	//This one include edge phosphates
	core::Real
	full_length_rmsd_over_residue_list( pose::Pose const & pose1, pose::Pose const & pose2, utility::vector1 < Size > const & residue_list, std::string const & full_sequence, bool const verbose, bool const ignore_virtual_atom ){

		using namespace ObjexxFCL;


		if ( pose1.sequence() != full_sequence ){
			TR << "pose1.sequence() = " << pose1.sequence() << std::endl;
			TR << "pose2.sequence() = " << pose2.sequence() << std::endl;
			TR << "full_sequence = " << full_sequence << std::endl;
			utility_exit_with_message( "pose1.sequence() != full_sequence" );
		}

		if ( pose2.sequence() != full_sequence ){
			TR << "pose1.sequence() = " << pose1.sequence() << std::endl;
			TR << "pose2.sequence() = " << pose2.sequence() << std::endl;
			TR << "full_sequence = " << full_sequence << std::endl;
			utility_exit_with_message( "pose2.sequence() != full_sequence" );
		}

		Size const total_res = pose1.total_residue();

		if ( verbose ){
			output_title_text( "Enter full_length_rmsd_over_residue_list function", TR );
			output_boolean( "ignore_virtual_atom = ", ignore_virtual_atom, TR ); TR << std::endl;
			output_seq_num_list( "residue_list = ", residue_list, TR, 30 );
		}

		Size atom_count = 0;
		Real sum_sd = 0;

		for ( Size i = 1; i <= residue_list.size(); i++ ){

			Size const full_seq_num = residue_list[i];

			bool is_prepend = false;
			bool both_pose_res_is_virtual = false;

			if ( pose1.residue( full_seq_num ).has_variant_type( "VIRTUAL_RNA_RESIDUE" ) && pose2.residue( full_seq_num ).has_variant_type( "VIRTUAL_RNA_RESIDUE" ) ){
				both_pose_res_is_virtual = true;
			}

			if ( ( full_seq_num + 1 ) <= total_res ){
				if ( pose1.residue( full_seq_num ).has_variant_type( "VIRTUAL_RNA_RESIDUE" ) ){
					if ( pose1.residue( full_seq_num + 1 ).has_variant_type( "VIRTUAL_RNA_RESIDUE_UPPER" ) == false ){ //consistency_check
						utility_exit_with_message( "pose1's full_seq_num = " + string_of( full_seq_num ) + "  is a virtual res but seq_num + 1 is not a virtual_res_upper!" );
					}
				}

				if ( pose2.residue( full_seq_num ).has_variant_type( "VIRTUAL_RNA_RESIDUE" ) ){
					if ( pose2.residue( full_seq_num + 1 ).has_variant_type( "VIRTUAL_RNA_RESIDUE_UPPER" ) == false ){ //consistency_check
						utility_exit_with_message( "pose2's full_seq_num = " + string_of( full_seq_num ) + "  is a virtual res but seq_num + 1 is not a virtual_res_upper!" );
					}
				}
			}

			if ( verbose ){
				TR << "full_seq_num = " << full_seq_num;
				output_boolean( " is_prepend = ", is_prepend, TR );
				output_boolean( " both_pose_res_is_virtual = ", both_pose_res_is_virtual, TR ); TR << std::endl;
			}

			if ( both_pose_res_is_virtual ) continue;

			//add atom in the suites to atom_count
			//add sd of each atom to sum_sd
			suite_square_deviation( pose1, pose2, is_prepend, full_seq_num, full_seq_num, atom_count, sum_sd, verbose, ignore_virtual_atom );

			if ( ( ( full_seq_num + 1 ) <= total_res ) && residue_list.has_value( full_seq_num + 1 ) == false ){

				if ( verbose ) TR << "Phosphate_edge_res_( full_seq_num + 1 ) = " << full_seq_num + 1 << std::endl;

				phosphate_square_deviation( pose1, pose2, full_seq_num + 1, full_seq_num + 1, atom_count, sum_sd, verbose, ignore_virtual_atom );
			}

		}


		sum_sd = sum_sd/( atom_count );
		Real rmsd = sqrt( sum_sd );

		if ( atom_count == 0 ) rmsd = 0.0; //special case...implement this on May 5, 2010

		if ( verbose ){
			TR << "sum_sd = " << sum_sd << " atom_count = " << atom_count << " rmsd = " << rmsd << std::endl;
			output_title_text( "Exit In full_length_rmsd_over_residue_list function", TR );
		}

		return ( std::max( 0.01, rmsd ) );

	}
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	Real
	rmsd_over_residue_list( pose::Pose const & pose1, pose::Pose const & pose2, utility::vector1 < Size > const & residue_list, std::map< core::Size, core::Size > const & full_to_sub, std::map< core::Size, bool > const & is_prepend_map, bool const verbose, bool const ignore_virtual_atom ){

		if ( verbose ){
			output_title_text( "Enter rmsd_over_residue_list function", TR );
			output_boolean( "ignore_virtual_atom = ", ignore_virtual_atom, TR ); TR << std::endl;
			output_seq_num_list( "residue_list = ", residue_list, TR, 30 );
		}


		Size atom_count = 0;
		Real sum_sd = 0;

		for ( Size i = 1; i <= residue_list.size(); i++ ){

			Size const full_seq_num = residue_list[i];
 			Size const seq_num = full_to_sub.find( full_seq_num )->second;
			bool is_prepend = is_prepend_map.find( full_seq_num )->second;
			bool both_pose_res_is_virtual = false;
			if ( pose1.residue( seq_num ).has_variant_type( "VIRTUAL_RNA_RESIDUE" ) && pose2.residue( seq_num ).has_variant_type( "VIRTUAL_RNA_RESIDUE" ) ){
				both_pose_res_is_virtual = true;
			}


			if ( verbose ){
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

		if ( atom_count == 0 ) rmsd = 0.0; //special case...implement this on May 5, 2010

		if ( verbose ){
			TR << "sum_sd = " << sum_sd << " atom_count = " << atom_count << " rmsd = " << rmsd << std::endl;
			output_title_text( "Exit rmsd_over_residue_list function", TR );
		}

	  	return ( std::max( 0.01, rmsd ) );
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////


	Real
	rmsd_over_residue_list( pose::Pose const & pose1, pose::Pose const & pose2, StepWiseRNA_JobParametersCOP job_parameters_, bool const ignore_virtual_atom ){

		utility::vector1 < core::Size > const & rmsd_res_list = job_parameters_->rmsd_res_list();
		std::map< core::Size, core::Size > const & full_to_sub = job_parameters_->const_full_to_sub();
		std::map< core::Size, bool > const & is_prepend_map = job_parameters_->is_prepend_map();

		return rmsd_over_residue_list( pose1, pose2, rmsd_res_list, full_to_sub, is_prepend_map, false /*verbose*/, ignore_virtual_atom );

	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	print_heavy_atoms( Size const & suite_num_1, Size const & suite_num_2, pose::Pose const & pose1, pose::Pose const & pose2 ){

		using namespace conformation;
		Size num_atoms;

		num_atoms = std::max( pose1.residue( suite_num_1 ).nheavyatoms(), pose2.residue( suite_num_2 ).nheavyatoms() );

		TR << "num_atoms: " << num_atoms << std::endl;

		for ( Size n = 1; n <= num_atoms;  n++ ){

			TR << " atom num = " <<  n;
			TR << "  atom_name of the pose1 " <<  pose1.residue( suite_num_1 ).atom_name( n );
			TR << "  atom_name of the pose2 " <<  pose2.residue( suite_num_2 ).atom_name( n ) << std::endl;

		}
	}
//

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
		if ( atom_name_1 != atom_name_2 ){
			utility_exit_with_message( "atom_name_1 != atom_name_2, atom_name_1 = " + atom_name_1 + " atom_name_2 = " + atom_name_2 );
		}

//		if(rsd_1.atom_type(atomno_1).name()=="VIRT"  || rsd_2.atom_type(atomno_2).name()=="VIRT") {
//			TR << "atom_name_1= " << atom_name_1 << " atom_name_2= " << atom_name_2 << std::endl;
//			utility_exit_with_message( "rsd_1.atom_type(n).name()==\"VIRT\"  || rsd_2.atom_type(n).name()==\"VIRT\" =TRUE!");
//		}
		////////////////////////////////////////////////////////////////////////////////////////////////

  		Distance const dist_squared = ( rsd_1.xyz( atomno_1 ) - rsd_2.xyz( atomno_2 ) ).length_squared();


		if ( verbose ){
			TR << " atom_name of the atom1 = " << atom_name_1 << " " << rsd_1.seqpos();
			TR << " atom_name of the atom2 = " << atom_name_2 << " " << rsd_2.seqpos();
			TR << " Dist_squared = " << dist_squared << std::endl;
		}

		return dist_squared;

	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Create on Sept 24, 2010...should really integrate these rmsd square deviation functions....cause now it is just copy and paste copy

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
  		for ( Size n = 1; n <= num_side_chain_atom; n++ ){ //This INCLUDE the O2prime oxygen

			Size const atomno_1 = ( n - 1 ) + first_sidechain_atom1;
			Size const atomno_2 = ( n - 1 ) + first_sidechain_atom2;

			conformation::Residue const & rsd_1 = pose1.residue( moving_res_1 );
			conformation::Residue const & rsd_2 = pose2.residue( moving_res_2 );

			if ( rsd_1.type().atom_name( atomno_1 ) == " O2'" ) continue; //Exclude the O2prime oxygen
			if ( rsd_2.type().atom_name( atomno_2 ) == " O2'" ) continue; //Exclude the O2prime oxygen

			if ( ignore_virtual_atom ){
				if ( rsd_1.atom_type( atomno_1 ).name() == "VIRT"  || rsd_2.atom_type( atomno_2 ).name() == "VIRT" ) continue;
			}

			if ( rsd_1.atom_type( atomno_1 ).name() == "VIRT"  && rsd_2.atom_type( atomno_2 ).name() == "VIRT" ) { //Change this to "AND" on Apr 5
//					TR << "Both atoms are VIRTUAL! moving_res_1= " << moving_res_1 << " moving_res_2= " << moving_res_2 << " atomno_1= " << atomno_1 << " atomno_2= " << atomno_2 << std::endl;
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

		for ( Size atomno = 1; atomno <= 4; atomno++ ){

			conformation::Residue const & rsd_1 = pose1.residue( moving_res_1 );
			conformation::Residue const & rsd_2 = pose2.residue( moving_res_2 );

			if ( ignore_virtual_atom ){
				if ( rsd_1.atom_type( atomno ).name() == "VIRT"  || rsd_2.atom_type( atomno ).name() == "VIRT" ) continue;
			}

			if ( rsd_1.atom_type( atomno ).name() == "VIRT"  && rsd_2.atom_type( atomno ).name() == "VIRT" ) continue;

			atom_count++;
			sum_sd = sum_sd + atom_square_deviation( rsd_1, rsd_2, atomno, atomno, verbose );

		}

	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	core::Real
	phosphate_base_phosphate_rmsd( pose::Pose const & pose1, pose::Pose const & pose2, Size const & moving_res_num, bool const ignore_virtual_atom ){

		Size atom_count = 0;
  		Real sum_sd = 0;

		phosphate_base_phosphate_square_deviation( pose1, pose2, moving_res_num, moving_res_num, atom_count, sum_sd, false, ignore_virtual_atom );

		sum_sd = sum_sd/( atom_count );
  		Real rmsd = sqrt( sum_sd );

		if ( atom_count == 0 ) rmsd = 0.0; //special case...implement this on June_11, 2010, took me a whole day to debug this since buggy only on Biox compiler!

		return ( std::max( 0.01, rmsd ) );

	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	phosphate_base_phosphate_square_deviation( pose::Pose const & pose1, pose::Pose const & pose2, Size const & moving_res_1, Size const & moving_res_2, Size& atom_count, Real& sum_sd, bool verbose, bool const ignore_virtual_atom ){

  		chemical::AA const & res_aa =  pose1.residue( moving_res_1 ).aa();
  		chemical::AA const & res_aa2 =  pose2.residue( moving_res_2 ).aa();

		if ( res_aa != res_aa2 ) utility_exit_with_message( "res_aa ( " + name_from_aa( res_aa ) + " ) != res_aa2 ( " + name_from_aa( res_aa2 ) + " ) " );

  		Size const first_sidechain_atom1 = pose1.residue( moving_res_1 ).first_sidechain_atom();
  		Size const first_sidechain_atom2 = pose2.residue( moving_res_2 ).first_sidechain_atom();

		Size const num_side_chain_atom = get_num_side_chain_atom_from_res_name( res_aa, verbose );

		if ( verbose ) TR << " MOVING_RES_1: " << moving_res_1 << " MOVING_RES_2: " << moving_res_2 << std::endl;

		Size const num_heavy_backbone_atoms = 11; //RNA contain 11 heavy backbone atoms.

		for ( Size atomno = 1; atomno <= num_heavy_backbone_atoms; atomno++ ){

			Size const res_count = ( atomno <= 4 ) ? 2 : 1; //atom 1-4 are " P  ", " OP2", " OP1" and " O5'"

			for ( Size ii = 1; ii < res_count; ii++ ){

				Size const res_num_1 = moving_res_1 + ( ii - 1 );
				Size const res_num_2 = moving_res_2 + ( ii - 1 );

				conformation::Residue const & rsd_1 = pose1.residue( res_num_1 );
				conformation::Residue const & rsd_2 = pose2.residue( res_num_2 );

				if ( ignore_virtual_atom ){
					if ( rsd_1.atom_type( atomno ).name() == "VIRT"  || rsd_2.atom_type( atomno ).name() == "VIRT" ) continue;
				}

				if ( rsd_1.atom_type( atomno ).name() == "VIRT"  && rsd_2.atom_type( atomno ).name() == "VIRT" ) continue;

				atom_count++;
				sum_sd = sum_sd + atom_square_deviation( rsd_1, rsd_2, atomno, atomno, verbose );
			}
		}

  		//Need to use num_side_chain_atom from pose1 since a silly bug in Rosetta miscalculate num_heavy_atom by considering
		//the virtaul O2prime hydrogen to be heavy_atom when it is set to virtual in the current_pose_screen
  		for ( Size n = 1; n <= num_side_chain_atom; n++ ){ //INCLUDE the O2prime oxygen

			Size const atomno_1 = ( n - 1 ) + first_sidechain_atom1;
			Size const atomno_2 = ( n - 1 ) + first_sidechain_atom2;

			conformation::Residue const & rsd_1 = pose1.residue( moving_res_1 );
			conformation::Residue const & rsd_2 = pose2.residue( moving_res_2 );

			if ( ignore_virtual_atom ){
				if ( rsd_1.atom_type( atomno_1 ).name() == "VIRT"  || rsd_2.atom_type( atomno_2 ).name() == "VIRT" ) continue;
			}

			if ( rsd_1.atom_type( atomno_1 ).name() == "VIRT"  && rsd_2.atom_type( atomno_2 ).name() == "VIRT" ) { //Change this to "AND" on Apr 5
//					TR << "Both atoms are VIRTUAL! moving_res_1= " << moving_res_1 << " moving_res_2= " << moving_res_2 << " atomno_1= " << atomno_1 << " atomno_2= " << atomno_2 << std::endl;
				continue;
			}

			atom_count++;
			sum_sd = sum_sd + atom_square_deviation( rsd_1, rsd_2, atomno_1, atomno_2, verbose );
	 	}

	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	//	Could get the atom_num using the atom_name but I don't think this will be fast Size atom_P = rsd.atom_index( " P  " ); Jan 28, 2010 Parin S////////////////////////////////////////////

	void
	suite_square_deviation( pose::Pose const & pose1, pose::Pose const & pose2, bool const & prepend_res, Size const & moving_res_1, Size const & moving_res_2, Size& atom_count, Real& sum_sd, bool verbose, bool const ignore_virtual_atom ){

  		chemical::AA const & res_aa =  pose1.residue( moving_res_1 ).aa();
  		chemical::AA const & res_aa2 =  pose2.residue( moving_res_2 ).aa();

		if ( res_aa != res_aa2 ) utility_exit_with_message( "res_aa ( " + name_from_aa( res_aa ) + " ) != res_aa2 ( " + name_from_aa( res_aa2 ) +  " ) " );

  		Size const first_sidechain_atom1 = pose1.residue( moving_res_1 ).first_sidechain_atom();
  		Size const first_sidechain_atom2 = pose2.residue( moving_res_2 ).first_sidechain_atom();

		Size const num_side_chain_atom = get_num_side_chain_atom_from_res_name( res_aa, verbose );

		if ( false && verbose ){
			TR << " residue type1 = " <<  res_aa << " residue type2 = " <<  res_aa2;
			TR << " 1st_side_atom1 = " <<  first_sidechain_atom1;
			TR << " 1st_side_atom2 = " <<  first_sidechain_atom2;
			TR << " nheavyatoms1 = " <<  pose1.residue( moving_res_1 ).nheavyatoms();
			TR << " nheavyatoms2 = " <<  pose2.residue( moving_res_2 ).nheavyatoms() << std::endl;
			TR << " num_side_chain_atom = " <<  num_side_chain_atom << std::endl;
			print_heavy_atoms( moving_res_1, moving_res_2, pose1, pose2 );
		}

		if ( verbose ) TR << " MOVING_RES_1: " << moving_res_1 << " MOVING_RES_2: " << moving_res_2 << "    PREPEND? " << prepend_res << std::endl;

		Size const num_heavy_backbone_atoms = 11; //RNA contain 11 heavy backbone atoms.

		for ( Size atomno = 1; atomno <= num_heavy_backbone_atoms; atomno++ ){

			//atom 1-4 are " P  ", " OP2", " OP1" and " O5'"
			Size const res_num_1 = ( prepend_res && atomno <= 4 ) ? moving_res_1 + 1: moving_res_1;
			Size const res_num_2 = ( prepend_res && atomno <= 4 ) ? moving_res_2 + 1: moving_res_2;


			conformation::Residue const & rsd_1 = pose1.residue( res_num_1 );
			conformation::Residue const & rsd_2 = pose2.residue( res_num_2 );

			if ( ignore_virtual_atom ){
				if ( rsd_1.atom_type( atomno ).name() == "VIRT"  || rsd_2.atom_type( atomno ).name() == "VIRT" ) continue;
			}

			if ( rsd_1.atom_type( atomno ).name() == "VIRT"  && rsd_2.atom_type( atomno ).name() == "VIRT" ) {
//					TR << "Both atoms are VIRTUAL! res_num_1= " << res_num_1 << " res_num_2= " << res_num_2 << " atomno= " << atomno << std::endl;
				continue;
			}

			atom_count++;
			sum_sd = sum_sd + atom_square_deviation( rsd_1, rsd_2, atomno, atomno, verbose );
		}

  		//Need to use num_side_chain_atom from pose1 since a silly bug in Rosetta miscalculate num_heavy_atom by considering
			//the virtaul O2prime hydrogen to be heavy_atom when it is set to virtual in the current_pose_screen
  		for ( Size n = 1; n <= num_side_chain_atom; n++ ){ //INCLUDE the O2prime oxygen


			Size const atomno_1 = ( n - 1 ) + first_sidechain_atom1;
			Size const atomno_2 = ( n - 1 ) + first_sidechain_atom2;

			conformation::Residue const & rsd_1 = pose1.residue( moving_res_1 );
			conformation::Residue const & rsd_2 = pose2.residue( moving_res_2 );

			if ( ignore_virtual_atom ){
				if ( rsd_1.atom_type( atomno_1 ).name() == "VIRT"  || rsd_2.atom_type( atomno_2 ).name() == "VIRT" ) continue;
			}

			if ( rsd_1.atom_type( atomno_1 ).name() == "VIRT"  && rsd_2.atom_type( atomno_2 ).name() == "VIRT" ) { //Change this to "AND" on Apr 5
	//					TR << "Both atoms are VIRTUAL! moving_res_1= " << moving_res_1 << " moving_res_2= " << moving_res_2 << " atomno_1= " << atomno_1 << " atomno_2= " << atomno_2 << std::endl;
				continue;
			}

			atom_count++;
			sum_sd = sum_sd + atom_square_deviation( rsd_1, rsd_2, atomno_1, atomno_2, verbose );
	 	}


			// OP2<-->OP1 check on phosphate positions closest to fixed side
		if ( verbose && false ){
			Size const res_num_1 = ( prepend_res ) ? moving_res_1 + 1: moving_res_1;
			Size const res_num_2 = ( prepend_res ) ? moving_res_2 + 1: moving_res_2;
			Distance dist_squared = ( pose1.residue( res_num_1 ).xyz( 2 ) - pose2.residue( res_num_2 ).xyz( 3 ) ).length_squared();
			TR << " atom_name of the pose1 = " << pose1.residue( res_num_1 ).atom_name( 2 )  << " " << res_num_1;
			TR << " atom_name of the pose2 = " << pose2.residue( res_num_2 ).atom_name( 3 )  << " " << res_num_2;
			TR << " Switch Phosphate1 = " << " dist_squared = " << dist_squared << std::endl;
			dist_squared = ( pose1.residue( res_num_1 ).xyz( 3 ) - pose2.residue( res_num_2 ).xyz( 2 ) ).length_squared();
			TR << " atom_name of the pose1 = " << pose1.residue( res_num_1 ).atom_name( 3 )  << " " << res_num_1;
			TR << " atom_name of the pose2 = " << pose2.residue( res_num_2 ).atom_name( 2 )  << " " << res_num_2;
			TR << " Switch Phosphate2 = " << " dist_squared = " << dist_squared << std::endl;
		}


	}
	/////////////////////////////////////////////////////////////////////////////////////////////////
	//WARNING THIS FUNCTION ASSUME THAT THE GAP RESIDUES DOESN'T EXIST!!! COULD LEAD TO SERIOUS ERROR! Sept 17, 2010
	bool
	check_chain_closable( pose::Pose const & pose, Size const five_prime_chain_break_res, Size const gap_size ){
		//Use to call Calculate_dist_to_close_chain function, but that requires rebuild_unit_struct
		Size const three_prime_chain_break_res = five_prime_chain_break_res + 1; //THIS ASSUME THAT GAP RESIDUES DOESN'T EXIST

		return check_chain_closable( pose.residue( three_prime_chain_break_res ).xyz( "C5'" ), pose.residue( five_prime_chain_break_res ).xyz( "O3'" ), gap_size );
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	check_chain_closable( numeric::xyzVector< core::Real >  const & xyz_1, numeric::xyzVector< core::Real > const & xyz_2, core::Size const gap_size ){

		//Two possibilities
		//1. xyz_1-> five_prime_O3_xyz &&  xyz_2->three_prime_C5_xyz
		//2. xyz_2-> five_prime_O3_xyz &&  xyz_1->three_prime_C5_xyz

		using namespace ObjexxFCL;

		Distance const dist_squared = ( xyz_1 - xyz_2 ).length_squared();

		if ( gap_size > 1 ){ //new option Sept 18, 2010, most useful for long loop mode...

			Real const dist_cutoff = O3I_C5IPLUS2_MAX_DIST + ( ( gap_size - 1 )*( C5I_C5I_PLUS_ONE_MAX_DIST ) );

			//C5I_C5I_PLUS_ONE_MAX_DIST is the right choice and NOT O3I_O3I_PLUS_ONE_MAX_DIST, since we want to go from O3I to C5(I+gapsize+1)
			//Also C5I_C5I_PLUS_ONE_MAX_DIST is slightly larger than O3I_O3I_PLUS_ONE_MAX_DIST and hence the safe choice

			return( dist_squared  < ( dist_cutoff*dist_cutoff ) );

		} else if ( gap_size == 1 ) {

			//previously used 11.0138 as O3I_C5IPLUS2_MAX_DIS which is slight underestimate...TOO strict;
			return ( dist_squared < O3I_C5IPLUS2_MAX_DIST*O3I_C5IPLUS2_MAX_DIST );

		} else if ( gap_size == 0 ){

			static Distance const cutoff_distance_min_squared( 2.0 * 2.0 );
			static Distance const cutoff_distance_max_squared( 4.627 * 4.627 );

			//basically cannot close chain if the C5_O3_distance is either too short or too long.
			if ( ( dist_squared > cutoff_distance_max_squared ) || dist_squared < cutoff_distance_min_squared ) return false;

			return true;

		} else{
			std::string const exit_message = "Invalid gap_size = " + string_of( gap_size ) + " !!";
			utility_exit_with_message( exit_message );
		}

		utility_exit_with_message( "Should not reach this point of the function!" );
		exit( 1 ); //prevent compiler warning!
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////
	//Optimization for floating_base_chain_closure
	//Need to integrate this with standard check_chain_closable()
	bool
	check_chain_closable_floating_base( pose::Pose const & five_prime_pose, pose::Pose const & three_prime_pose,
																		 Size const five_prime_chain_break_res, Size const gap_size ){

		runtime_assert ( gap_size == 0 );

		//This is potentially dangerous since even in the case where gap_size != 0 ... it assumes that there is not residues between 3' and 5' res....
		Size const three_prime_chain_break_res = five_prime_chain_break_res + 1;

		Distance C5_O3_dist = ( three_prime_pose.residue( three_prime_chain_break_res ).xyz( "C5'" ) - five_prime_pose.residue( five_prime_chain_break_res ).xyz( "O3'" ) ).length();

		static Distance const C5_O3_min( 2.866000 );
		static Distance const C5_O3_max( 3.968000 );
		static Distance const leniency_dist( 0.0 );

		//basically cannot close chain if the C5_O3_distance is either too short or too long.
		if ( ( C5_O3_dist > ( C5_O3_max + leniency_dist ) ) || ( C5_O3_dist < ( C5_O3_min - leniency_dist ) ) ) return false;

		conformation::Residue const & five_prime_rsd = five_prime_pose.residue( five_prime_chain_break_res );
		conformation::Residue const & three_prime_rsd = three_prime_pose.residue( three_prime_chain_break_res );

		Distance C4_C3_min = 0.0;
		Distance C4_C3_max = 0.0;

		get_C4_C3_distance_range( five_prime_rsd, three_prime_rsd, C4_C3_min, C4_C3_max );

		Distance C4_C3_dist = ( three_prime_pose.residue( three_prime_chain_break_res ).xyz( " C4'" ) - five_prime_pose.residue( five_prime_chain_break_res ).xyz( " C3'" ) ).length();

		if ( ( C4_C3_dist > ( C4_C3_max + leniency_dist ) ) || ( C4_C3_dist < ( C4_C3_min - leniency_dist ) ) ) return false;

		return true;
	}

	void
	get_C4_C3_distance_range( conformation::Residue const & five_prime_rsd,
													 conformation::Residue const & three_prime_rsd,
												 	 Distance & C4_C3_dist_min,
												 	 Distance & C4_C3_dist_max ){


			numeric::xyzVector< Real > start_vector = five_prime_rsd.xyz( " O3'" ) - five_prime_rsd.xyz( " C3'" );
			numeric::xyzVector< Real > end_vector = three_prime_rsd.xyz( " C4'" ) - three_prime_rsd.xyz( " C5'" );

			start_vector.normalize();
			end_vector.normalize();

			Real dot_product = dot( start_vector, end_vector );


			if ( dot_product > -1.00 && dot_product < -0.95  ) { C4_C3_dist_min = 2.428;  C4_C3_dist_max = 4.337; }
			else if ( dot_product > -0.95 && dot_product < -0.90  ) { C4_C3_dist_min = 2.238;  C4_C3_dist_max = 4.582; }
			else if ( dot_product > -0.90 && dot_product < -0.85  ) { C4_C3_dist_min = 2.064;  C4_C3_dist_max = 4.743; }
			else if ( dot_product > -0.85 && dot_product < -0.80  ) { C4_C3_dist_min = 1.979;  C4_C3_dist_max = 4.882; }
			else if ( dot_product > -0.80 && dot_product < -0.75  ) { C4_C3_dist_min = 1.833;  C4_C3_dist_max = 4.995; }
			else if ( dot_product > -0.75 && dot_product < -0.70  ) { C4_C3_dist_min = 1.735;  C4_C3_dist_max = 5.099; }
			else if ( dot_product > -0.70 && dot_product < -0.65  ) { C4_C3_dist_min = 1.659;  C4_C3_dist_max = 5.195; }
			else if ( dot_product > -0.65 && dot_product < -0.60  ) { C4_C3_dist_min = 1.590;  C4_C3_dist_max = 5.273; }
			else if ( dot_product > -0.60 && dot_product < -0.55  ) { C4_C3_dist_min = 1.500;  C4_C3_dist_max = 5.347; }
			else if ( dot_product > -0.55 && dot_product < -0.50  ) { C4_C3_dist_min = 1.418;  C4_C3_dist_max = 5.417; }
			else if ( dot_product > -0.50 && dot_product < -0.45  ) { C4_C3_dist_min = 1.337;  C4_C3_dist_max = 5.488; }
			else if ( dot_product > -0.45 && dot_product < -0.40  ) { C4_C3_dist_min = 1.282;  C4_C3_dist_max = 5.552; }
			else if ( dot_product > -0.40 && dot_product < -0.35  ) { C4_C3_dist_min = 1.223;  C4_C3_dist_max = 5.611; }
			else if ( dot_product > -0.35 && dot_product < -0.30  ) { C4_C3_dist_min = 1.145;  C4_C3_dist_max = 5.659; }
			else if ( dot_product > -0.30 && dot_product < -0.25  ) { C4_C3_dist_min = 1.075;  C4_C3_dist_max = 5.713; }
			else if ( dot_product > -0.25 && dot_product < -0.20  ) { C4_C3_dist_min = 1.022;  C4_C3_dist_max = 5.769; }
			else if ( dot_product > -0.20 && dot_product < -0.15  ) { C4_C3_dist_min = 0.963;  C4_C3_dist_max = 5.812; }
			else if ( dot_product > -0.15 && dot_product < -0.10  ) { C4_C3_dist_min = 1.019;  C4_C3_dist_max = 5.861; }
			else if ( dot_product > -0.10 && dot_product < -0.05  ) { C4_C3_dist_min = 1.331;  C4_C3_dist_max = 5.904; }
			else if ( dot_product > -0.05 && dot_product < 0.00  ) { C4_C3_dist_min = 1.532;  C4_C3_dist_max = 5.942; }
			else if ( dot_product > 0.00 && dot_product < 0.05  ) { C4_C3_dist_min = 1.768;  C4_C3_dist_max = 5.979; }
			else if ( dot_product > 0.05 && dot_product < 0.10  ) { C4_C3_dist_min = 1.953;  C4_C3_dist_max = 6.017; }
			else if ( dot_product > 0.10 && dot_product < 0.15  ) { C4_C3_dist_min = 2.121;  C4_C3_dist_max = 6.046; }
			else if ( dot_product > 0.15 && dot_product < 0.20  ) { C4_C3_dist_min = 2.292;  C4_C3_dist_max = 6.083; }
			else if ( dot_product > 0.20 && dot_product < 0.25  ) { C4_C3_dist_min = 2.424;  C4_C3_dist_max = 6.118; }
			else if ( dot_product > 0.25 && dot_product < 0.30  ) { C4_C3_dist_min = 2.563;  C4_C3_dist_max = 6.140; }
			else if ( dot_product > 0.30 && dot_product < 0.35  ) { C4_C3_dist_min = 2.726;  C4_C3_dist_max = 6.171; }
			else if ( dot_product > 0.35 && dot_product < 0.40  ) { C4_C3_dist_min = 2.849;  C4_C3_dist_max = 6.200; }
			else if ( dot_product > 0.40 && dot_product < 0.45  ) { C4_C3_dist_min = 2.998;  C4_C3_dist_max = 6.219; }
			else if ( dot_product > 0.45 && dot_product < 0.50  ) { C4_C3_dist_min = 3.128;  C4_C3_dist_max = 6.245; }
			else if ( dot_product > 0.50 && dot_product < 0.55  ) { C4_C3_dist_min = 3.261;  C4_C3_dist_max = 6.261; }
			else if ( dot_product > 0.55 && dot_product < 0.60  ) { C4_C3_dist_min = 3.380;  C4_C3_dist_max = 6.284; }
			else if ( dot_product > 0.60 && dot_product < 0.65  ) { C4_C3_dist_min = 3.523;  C4_C3_dist_max = 6.298; }
			else if ( dot_product > 0.65 && dot_product < 0.70  ) { C4_C3_dist_min = 3.658;  C4_C3_dist_max = 6.315; }
			else if ( dot_product > 0.70 && dot_product < 0.75  ) { C4_C3_dist_min = 3.785;  C4_C3_dist_max = 6.329; }
			else if ( dot_product > 0.75 && dot_product < 0.80  ) { C4_C3_dist_min = 3.914;  C4_C3_dist_max = 6.340; }
			else if ( dot_product > 0.80 && dot_product < 0.85  ) { C4_C3_dist_min = 4.065;  C4_C3_dist_max = 6.350; }
			else if ( dot_product > 0.85 && dot_product < 0.90  ) { C4_C3_dist_min = 4.209;  C4_C3_dist_max = 6.356; }
			else if ( dot_product > 0.90 && dot_product < 0.95  ) { C4_C3_dist_min = 4.374;  C4_C3_dist_max = 6.357; }
			else if ( dot_product > 0.95 && dot_product < 1.00  ) { C4_C3_dist_min = 4.570;  C4_C3_dist_max = 6.349; }
			else{
				TR << "dot_product = " << dot_product << std::endl;
				utility_exit_with_message( "Invalid dot_product!" );
			}

	}

/*
dot_min = -1.000000  dot_max = -0.950000  C4_C3_dist_min = 2.428000  C4_C3_dist_max 4.337000 C5_O3_dist_min = 2.866000  C5_O3_dist_max 3.881000
dot_min = -0.950000  dot_max = -0.900000  C4_C3_dist_min = 2.238000  C4_C3_dist_max 4.582000 C5_O3_dist_min = 2.866000  C5_O3_dist_max 3.923000
dot_min = -0.900000  dot_max = -0.850000  C4_C3_dist_min = 2.064000  C4_C3_dist_max 4.743000 C5_O3_dist_min = 2.866000  C5_O3_dist_max 3.951000
dot_min = -0.850000  dot_max = -0.800000  C4_C3_dist_min = 1.979000  C4_C3_dist_max 4.882000 C5_O3_dist_min = 2.866000  C5_O3_dist_max 3.966000
dot_min = -0.800000  dot_max = -0.750000  C4_C3_dist_min = 1.833000  C4_C3_dist_max 4.995000 C5_O3_dist_min = 2.866000  C5_O3_dist_max 3.968000
dot_min = -0.750000  dot_max = -0.700000  C4_C3_dist_min = 1.735000  C4_C3_dist_max 5.099000 C5_O3_dist_min = 2.866000  C5_O3_dist_max 3.968000
dot_min = -0.700000  dot_max = -0.650000  C4_C3_dist_min = 1.659000  C4_C3_dist_max 5.195000 C5_O3_dist_min = 2.866000  C5_O3_dist_max 3.968000
dot_min = -0.650000  dot_max = -0.600000  C4_C3_dist_min = 1.590000  C4_C3_dist_max 5.273000 C5_O3_dist_min = 2.866000  C5_O3_dist_max 3.968000
dot_min = -0.600000  dot_max = -0.550000  C4_C3_dist_min = 1.500000  C4_C3_dist_max 5.347000 C5_O3_dist_min = 2.866000  C5_O3_dist_max 3.968000
dot_min = -0.550000  dot_max = -0.500000  C4_C3_dist_min = 1.418000  C4_C3_dist_max 5.417000 C5_O3_dist_min = 2.866000  C5_O3_dist_max 3.968000
dot_min = -0.500000  dot_max = -0.450000  C4_C3_dist_min = 1.337000  C4_C3_dist_max 5.488000 C5_O3_dist_min = 2.866000  C5_O3_dist_max 3.968000
dot_min = -0.450000  dot_max = -0.400000  C4_C3_dist_min = 1.282000  C4_C3_dist_max 5.552000 C5_O3_dist_min = 2.866000  C5_O3_dist_max 3.968000
dot_min = -0.400000  dot_max = -0.350000  C4_C3_dist_min = 1.223000  C4_C3_dist_max 5.611000 C5_O3_dist_min = 2.866000  C5_O3_dist_max 3.968000
dot_min = -0.350000  dot_max = -0.300000  C4_C3_dist_min = 1.145000  C4_C3_dist_max 5.659000 C5_O3_dist_min = 2.866000  C5_O3_dist_max 3.968000
dot_min = -0.300000  dot_max = -0.250000  C4_C3_dist_min = 1.075000  C4_C3_dist_max 5.713000 C5_O3_dist_min = 2.866000  C5_O3_dist_max 3.968000
dot_min = -0.250000  dot_max = -0.200000  C4_C3_dist_min = 1.022000  C4_C3_dist_max 5.769000 C5_O3_dist_min = 2.866000  C5_O3_dist_max 3.968000
dot_min = -0.200000  dot_max = -0.150000  C4_C3_dist_min = 0.963000  C4_C3_dist_max 5.812000 C5_O3_dist_min = 2.866000  C5_O3_dist_max 3.968000
dot_min = -0.150000  dot_max = -0.100000  C4_C3_dist_min = 1.019000  C4_C3_dist_max 5.861000 C5_O3_dist_min = 2.866000  C5_O3_dist_max 3.968000
dot_min = -0.100000  dot_max = -0.050000  C4_C3_dist_min = 1.331000  C4_C3_dist_max 5.904000 C5_O3_dist_min = 2.866000  C5_O3_dist_max 3.968000
dot_min = -0.050000  dot_max = 0.000000  C4_C3_dist_min = 1.532000  C4_C3_dist_max 5.942000 C5_O3_dist_min = 2.866000  C5_O3_dist_max 3.968000
dot_min = 0.000000  dot_max = 0.050000  C4_C3_dist_min = 1.768000  C4_C3_dist_max 5.979000 C5_O3_dist_min = 2.866000  C5_O3_dist_max 3.968000
dot_min = 0.050000  dot_max = 0.100000  C4_C3_dist_min = 1.953000  C4_C3_dist_max 6.017000 C5_O3_dist_min = 2.866000  C5_O3_dist_max 3.968000
dot_min = 0.100000  dot_max = 0.150000  C4_C3_dist_min = 2.121000  C4_C3_dist_max 6.046000 C5_O3_dist_min = 2.866000  C5_O3_dist_max 3.968000
dot_min = 0.150000  dot_max = 0.200000  C4_C3_dist_min = 2.292000  C4_C3_dist_max 6.083000 C5_O3_dist_min = 2.866000  C5_O3_dist_max 3.968000
dot_min = 0.200000  dot_max = 0.250000  C4_C3_dist_min = 2.424000  C4_C3_dist_max 6.118000 C5_O3_dist_min = 2.866000  C5_O3_dist_max 3.968000
dot_min = 0.250000  dot_max = 0.300000  C4_C3_dist_min = 2.563000  C4_C3_dist_max 6.140000 C5_O3_dist_min = 2.866000  C5_O3_dist_max 3.968000
dot_min = 0.300000  dot_max = 0.350000  C4_C3_dist_min = 2.726000  C4_C3_dist_max 6.171000 C5_O3_dist_min = 2.866000  C5_O3_dist_max 3.968000
dot_min = 0.350000  dot_max = 0.400000  C4_C3_dist_min = 2.849000  C4_C3_dist_max 6.200000 C5_O3_dist_min = 2.866000  C5_O3_dist_max 3.968000
dot_min = 0.400000  dot_max = 0.450000  C4_C3_dist_min = 2.998000  C4_C3_dist_max 6.219000 C5_O3_dist_min = 2.866000  C5_O3_dist_max 3.968000
dot_min = 0.450000  dot_max = 0.500000  C4_C3_dist_min = 3.128000  C4_C3_dist_max 6.245000 C5_O3_dist_min = 2.866000  C5_O3_dist_max 3.968000
dot_min = 0.500000  dot_max = 0.550000  C4_C3_dist_min = 3.261000  C4_C3_dist_max 6.261000 C5_O3_dist_min = 2.866000  C5_O3_dist_max 3.968000
dot_min = 0.550000  dot_max = 0.600000  C4_C3_dist_min = 3.380000  C4_C3_dist_max 6.284000 C5_O3_dist_min = 2.866000  C5_O3_dist_max 3.968000
dot_min = 0.600000  dot_max = 0.650000  C4_C3_dist_min = 3.523000  C4_C3_dist_max 6.298000 C5_O3_dist_min = 2.866000  C5_O3_dist_max 3.968000
dot_min = 0.650000  dot_max = 0.700000  C4_C3_dist_min = 3.658000  C4_C3_dist_max 6.315000 C5_O3_dist_min = 2.866000  C5_O3_dist_max 3.968000
dot_min = 0.700000  dot_max = 0.750000  C4_C3_dist_min = 3.785000  C4_C3_dist_max 6.329000 C5_O3_dist_min = 2.866000  C5_O3_dist_max 3.968000
dot_min = 0.750000  dot_max = 0.800000  C4_C3_dist_min = 3.914000  C4_C3_dist_max 6.340000 C5_O3_dist_min = 2.866000  C5_O3_dist_max 3.968000
dot_min = 0.800000  dot_max = 0.850000  C4_C3_dist_min = 4.065000  C4_C3_dist_max 6.350000 C5_O3_dist_min = 2.866000  C5_O3_dist_max 3.968000
dot_min = 0.850000  dot_max = 0.900000  C4_C3_dist_min = 4.209000  C4_C3_dist_max 6.356000 C5_O3_dist_min = 2.866000  C5_O3_dist_max 3.968000
dot_min = 0.900000  dot_max = 0.950000  C4_C3_dist_min = 4.374000  C4_C3_dist_max 6.357000 C5_O3_dist_min = 2.866000  C5_O3_dist_max 3.968000
dot_min = 0.950000  dot_max = 1.000000  C4_C3_dist_min = 4.570000  C4_C3_dist_max 6.349000 C5_O3_dist_min = 2.866000  C5_O3_dist_max 3.968000

*/

	//////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	Freeze_sugar_torsions( core::kinematics::MoveMap& mm, Size const total_residue ){

	 using namespace core::id;

	 TR << "Freeze pose sugar torsions, total_residue = " << total_residue << std::endl;

	 for ( Size i = 1; i <= total_residue; i++ ){

			mm.set( TorsionID( i , id::BB,  4 ), false ); //delta_i
			mm.set( TorsionID( i , id::CHI, 2 ), false ); //nu2_i
			mm.set( TorsionID( i , id::CHI, 3 ), false );	//nu1_i

	 }
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	output_boolean( std::string const & tag, bool boolean, std::ostream & outstream /* = std::cout */ ){

		using namespace ObjexxFCL;
		using namespace ObjexxFCL::format;
		outstream << tag;

		if ( boolean == true ){
			outstream << A( 4, "T" );
		} else {
			outstream << A( 4, "F" );
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////
	void
	output_boolean( bool boolean, std::ostream & outstream /* = std::cout */ ){

		using namespace ObjexxFCL;
		using namespace ObjexxFCL::format;

		if ( boolean == true ){
			outstream << A( 4, "T" );
		} else {
			outstream << A( 4, "F" );
		}
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	output_movemap( kinematics::MoveMap const & mm, core::pose::Pose const & pose, std::ostream & outstream /* = std::cout */ ){

		using namespace ObjexxFCL;
		using namespace ObjexxFCL::format;
		using namespace core::id;
		using namespace core::kinematics;

		Size const total_residue = pose.total_residue();

		Size spacing = 10;

		outstream << "--------------------------------------------------------------------------------------" << std::endl;
		outstream << "Movemap ( in term of partial_pose seq_num ): " << std::endl;
		outstream << A( spacing, "res_num" ) << A( spacing, " alpha  " ) << A( spacing, "  beta  " ) << A( spacing, " gamma  " );
		outstream << A( spacing, " delta  " ) << A( spacing, "eplison " ) << A( spacing, "  zeta  " ) << A( spacing, " chi_1  " );
		outstream << A( spacing, "  nu_2  " ) << A( spacing, "  nu_1  " ) << A( spacing, "chi_O2' " ) << std::endl;

		for ( Size n = 1; n <= total_residue; n++ ){

			outstream << I( spacing, 3, n );
			outstream << A(  - 1 + ( spacing - 4 )/2, "" ); output_boolean( mm.get( TorsionID( n, id::BB,  1 ) ), outstream ); outstream << A( 1 + ( spacing - 4 )/2, "" );
			outstream << A(  - 1 + ( spacing - 4 )/2, "" ); output_boolean( mm.get( TorsionID( n, id::BB,  2 ) ), outstream ); outstream << A( 1 + ( spacing - 4 )/2, "" );
			outstream << A(  - 1 + ( spacing - 4 )/2, "" ); output_boolean( mm.get( TorsionID( n, id::BB,  3 ) ), outstream ); outstream << A( 1 + ( spacing - 4 )/2, "" );
			outstream << A(  - 1 + ( spacing - 4 )/2, "" ); output_boolean( mm.get( TorsionID( n, id::BB,  4 ) ), outstream ); outstream << A( 1 + ( spacing - 4 )/2, "" );
			outstream << A(  - 1 + ( spacing - 4 )/2, "" ); output_boolean( mm.get( TorsionID( n, id::BB,  5 ) ), outstream ); outstream << A( 1 + ( spacing - 4 )/2, "" );
			outstream << A(  - 1 + ( spacing - 4 )/2, "" ); output_boolean( mm.get( TorsionID( n, id::BB,  6 ) ), outstream ); outstream << A( 1 + ( spacing - 4 )/2, "" );
			outstream << A(  - 1 + ( spacing - 4 )/2, "" ); output_boolean( mm.get( TorsionID( n, id::CHI, 1 ) ), outstream ); outstream << A( 1 + ( spacing - 4 )/2, "" );
			outstream << A(  - 1 + ( spacing - 4 )/2, "" ); output_boolean( mm.get( TorsionID( n, id::CHI, 2 ) ), outstream ); outstream << A( 1 + ( spacing - 4 )/2, "" );
			outstream << A(  - 1 + ( spacing - 4 )/2, "" ); output_boolean( mm.get( TorsionID( n, id::CHI, 3 ) ), outstream ); outstream << A( 1 + ( spacing - 4 )/2, "" );
			outstream << A(  - 1 + ( spacing - 4 )/2, "" ); output_boolean( mm.get( TorsionID( n, id::CHI, 4 ) ), outstream ); outstream << A( 1 + ( spacing - 4 )/2, "" );
			outstream << std::endl;
		}
		outstream << "--------------------------------------------------------------------------------------" << std::endl;



		outstream << "print movemap jump_points [explicit method]: " << std::endl;
		for ( Size n = 1; n <= pose.fold_tree().num_jump(); n++ ){
			Size const jump_pos1( pose.fold_tree().upstream_jump_residue( n ) );
			Size const jump_pos2( pose.fold_tree().downstream_jump_residue( n ) );

			outstream << "n = " << n << " | jump_pos1 = " << jump_pos1 << " | jump_pos2 = " << jump_pos2;
			outstream << " | mm.get_jump( n ) = "; output_boolean( mm.get_jump( n ), outstream );
			outstream << std::endl;

		}
		outstream << "--------------------------------------------------------------------------------------" << std::endl;

		//From core/kinematic/MoveMap.hh
		//typedef std::map< id::JumpID, bool > JumpID_Map
		outstream << "print movemap jump_points [iterator method]: " << std::endl;
	  for ( std::map< id::JumpID, bool > ::const_iterator it = mm.jump_id_begin(); it != mm.jump_id_end(); it++ ){
	    outstream << "movemap jump == true for jump_pos1 = " << it->first << " | jump_pos2 = " << it->second << std::endl;
		}
		outstream << "--------------------------------------------------------------------------------------" << std::endl;

	}

	/////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< core::Size >
	get_surrounding_O2prime_hydrogen( pose::Pose const & pose, utility::vector1< core::Size > const & moving_res, bool verbose ){

		using namespace core::chemical;
		using namespace core::scoring;
		using namespace core::kinematics;
		using namespace ObjexxFCL;


		//Consistency_check
		for ( Size seq_num = 1; seq_num <= pose.total_residue(); seq_num++ ){

			if ( pose.residue( seq_num ).aa() == core::chemical::aa_vrt ) continue; //Fang's electron density code
			if ( !pose.residue( seq_num ).is_RNA() ) continue;

			core::conformation::Residue const & rsd = pose.residue( seq_num );
			Size const at = rsd.first_sidechain_atom();

			if ( rsd.type().atom_name( at ) != " O2'" ) {
				std::string const exit_message = "seq_num = " + string_of( seq_num ) + ", rsd.type().atom_name( at ) != \" O2'\" ";
				utility_exit_with_message( exit_message );
			}
		}

		utility::vector1< bool > is_O2prime_hydrogen_virtual_list;

		for ( Size seq_num = 1; seq_num <= pose.total_residue(); seq_num++ ){

			if ( pose.residue( seq_num ).aa() == core::chemical::aa_vrt ){
				if ( verbose ) TR << "res " << seq_num << " is core::chemical::aa_vrt! " << std::endl;
				is_O2prime_hydrogen_virtual_list.push_back( false ); //false since not virtual O2prime_hydrogen
				continue;
			}

			if ( !pose.residue( seq_num ).is_RNA() ){
				if ( verbose ) TR << "res " << seq_num << " is not RNA " << std::endl;
				is_O2prime_hydrogen_virtual_list.push_back( false ); //false since not virtual O2prime_hydrogen
				continue;
			}

			core::conformation::Residue const & rsd = pose.residue( seq_num );
			Size at = rsd.atom_index( "HO2'" );

			if ( rsd.atom_type( at ).name() == "VIRT" ){
				if ( verbose ) TR << "res " << seq_num << " has a virtual o2prime hydrogen! " << std::endl;
				is_O2prime_hydrogen_virtual_list.push_back( true );
			} else{
				is_O2prime_hydrogen_virtual_list.push_back( false );
			}
		}

		//March 17, 2012 extra precaution.
		if ( is_O2prime_hydrogen_virtual_list.size() != pose.total_residue() ) utility_exit_with_message( "is_O2prime_hydrogen_virtual_list.size() != pose.total_residue()" );

		utility::vector1< core::Size > surrounding_O2prime_hydrogen;

		Size num_o2prime_moving_res = 0;
		for ( Size n = 1; n <= moving_res.size(); n++ ){
			Size const seq_num = moving_res[n];

			if ( pose.residue( seq_num ).aa() == core::chemical::aa_vrt ) continue; //Fang's electron density code
			if ( is_O2prime_hydrogen_virtual_list[seq_num] ) continue;

			//if(verbose) TR << "res " << seq_num << " is a working_moving_res" << std::endl;
			surrounding_O2prime_hydrogen.push_back( seq_num );
			num_o2prime_moving_res++;
		}


		//1st layer, interaction between surrounding O2prime and moving_res
		for ( Size seq_num = 1; seq_num <= pose.total_residue(); seq_num++ ){

			if ( pose.residue( seq_num ).aa() == core::chemical::aa_vrt ) continue; //Fang's electron density code
			if ( surrounding_O2prime_hydrogen.has_value( seq_num ) ) continue;

			bool is_surrounding_res = false;

			core::conformation::Residue const & surrounding_rsd = pose.residue( seq_num );
			Size const surr_at = surrounding_rsd.first_sidechain_atom();

			if ( is_O2prime_hydrogen_virtual_list[seq_num] ) continue;

			//3.5 Angstrom for O2prime-normal

			for ( Size ii = 1; ii <= moving_res.size(); ii++ ){
				if ( is_surrounding_res == true ) break;

				core::conformation::Residue const & moving_rsd = pose.residue( moving_res[ii] );

				for ( Size moving_at = 1; moving_at <= moving_rsd.natoms(); moving_at++ ){

					if ( moving_rsd.atom_type( moving_at ).name() == "VIRT" ) continue;

					Real const cutoff_dist = ( moving_at == moving_rsd.first_sidechain_atom() ) ? 4.5: 3.5 ;

					//4.5 Angstrom interaction dist_cutoff for surrounding O2prime- moving_res O2prime
					//3.5 Angstrom interaction dist_cutoff for surround O2prime and other atoms of moving_res

					Real const dist_squared = ( surrounding_rsd.xyz( surr_at ) - moving_rsd.xyz( moving_at ) ).length_squared();

					if ( dist_squared  < cutoff_dist*cutoff_dist ){
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
		while ( true ){
			bool add_new_O2prime_hydrogen = false;

			for ( Size seq_num = 1; seq_num <= pose.total_residue(); seq_num++ ){

				if ( pose.residue( seq_num ).aa() == core::chemical::aa_vrt ) continue; //Fang's electron density code
				if ( surrounding_O2prime_hydrogen.has_value( seq_num ) ) continue;

				core::conformation::Residue const & rsd_1 = pose.residue( seq_num );
				Size const at_1 = rsd_1.first_sidechain_atom();

				if ( is_O2prime_hydrogen_virtual_list[seq_num] ) continue;

				for ( Size ii = 1; ii <= surrounding_O2prime_hydrogen.size(); ii++ ){

					core::conformation::Residue const & rsd_2 = pose.residue( surrounding_O2prime_hydrogen[ii] );
					Size const at_2 = rsd_2.first_sidechain_atom();

					Real const cutoff_dist = 4.5; 					//4.5 Angstrom interaction dist_cutoff for surrounding O2prime themselves.

					Real const dist_squared = ( rsd_1.xyz( at_1 ) - rsd_2.xyz( at_2 ) ).length_squared();

					if ( dist_squared  < cutoff_dist*cutoff_dist ){
						if ( verbose ) TR << "res " << seq_num << " is layer " << layer_num << " surrounding O2prime_hydrogen res, dist_squared = " << dist_squared << std::endl;
						surrounding_O2prime_hydrogen.push_back( seq_num );
						add_new_O2prime_hydrogen = true;
						break;
					}
				}
			}

			layer_num++;
			if ( add_new_O2prime_hydrogen == false ) break;
		}

		//consistency_check..
		for ( Size ii = 1; ii <= surrounding_O2prime_hydrogen.size(); ii++ ){

			Size const seq_num = surrounding_O2prime_hydrogen[ii];

			if ( is_O2prime_hydrogen_virtual_list[seq_num] ){
				std::string const exit_message = "surrounding_O2prime_hydrogen res " + string_of( seq_num ) + " has a virtual o2prime hydrogen!! ";
				utility_exit_with_message( exit_message );
			}
		}

		if ( verbose ){
			TR << "num_o2prime_moving_res = " << num_o2prime_moving_res << "  surrounding_O2prime_hydrogen.size() = " << surrounding_O2prime_hydrogen.size() << std::endl;
		}

		return surrounding_O2prime_hydrogen;

	}

	/////////////////////////////////////////////////////////////////////////////////////
	void
	o2prime_minimize( pose::Pose& pose, core::scoring::ScoreFunctionOP const & packer_scorefxn ){ //O2prime pack every position..

		 utility::vector1< core::Size > O2prime_pack_seq_num;

		for ( Size seq_num = 1; seq_num <= pose.total_residue(); seq_num++ ){
			if ( pose.residue( seq_num ).aa() == core::chemical::aa_vrt ) continue; //Fang's electron density code
			O2prime_pack_seq_num.push_back( seq_num );
		}

		o2prime_minimize( pose, packer_scorefxn, O2prime_pack_seq_num );

	}

	/////////////////////////////////////////////////////////////////////////////////////
	void
	o2prime_minimize( pose::Pose& pose, core::scoring::ScoreFunctionOP const & packer_scorefxn, utility::vector1< core::Size > const & O2prime_pack_seq_num ){

		output_seq_num_list( "O2prime_pack_seq_num = ", O2prime_pack_seq_num, TR );

		pack::task::PackerTaskOP task = create_standard_o2prime_pack_task( pose, O2prime_pack_seq_num );
		task->initialize_from_command_line();

		pack::rotamer_trials( pose, *packer_scorefxn, task );

	}

	/////////////////////////////////////////////////////////////////////////////////////
	//Created on Jan 20, 2012: Meant as the standard function to setup o2prime_pack_task to be called by both StepWiseRNA_ResidueSampler and StepWiseRNA_Minimizer...and in the future also FARFAR_minimizer.
	pack::task::PackerTaskOP
	create_standard_o2prime_pack_task( pose::Pose const & pose, utility::vector1< core::Size > const & O2prime_pack_seq_num ){

		pack::task::PackerTaskOP o2prime_pack_task = pack::task::TaskFactory::create_packer_task( pose );

		//Commented out the initialize_from_command_line() call below  on Jan 20, 2012: REASONS:
		/////1. The line is actually missing in the StepWiseRNA_ResidueSampler version (see below) which now calls this function.
		/////2. This is potentially dangerous since behavior can be influence by command_line flags (for example from command_line could turn on ex1 (base sampling) and etc which is not what we want!
		//o2prime_pack_task->initialize_from_command_line(); //Found only in original version packer_setup in StepWiseRNA_Util.cc BUT not the one in StepWiseRNA_ResidueSampler version.

		for ( Size seq_num = 1; seq_num <= pose.total_residue(); seq_num++ ){

			if ( pose.residue( seq_num ).aa() == core::chemical::aa_vrt ) continue; //Fang's electron density code
			if ( O2prime_pack_seq_num.has_value( seq_num ) && pose.residue( seq_num ).is_RNA() ){ //pack this residue!

				/*
				NOTE: the ex4 option allows sampling of extra torsion angles ( "ex" ) at
				the 2' - OH ( which is torsion "4" ). To pack the bases, you would have to
				add "ex1"; it should not happen! If you see bases getting repacked,
				that is a bug, and please let me know. ( From Rhiju email response )

				ex4 ( HO2prime ) is "on" by default | ex1 ( base chi ) is "off" by default
				*/

				o2prime_pack_task->nonconst_residue_task( seq_num ).and_extrachi_cutoff( 0 );
				o2prime_pack_task->nonconst_residue_task( seq_num ).or_ex4( true ); //extra O2prime sampling
				o2prime_pack_task->nonconst_residue_task( seq_num ).or_include_current( true ); //Jan 21, 2012: umm this might be set to true inside call to rotamer_trials() BUT even so, doesn't hurt so have it here as well!
				// How about bump check?
			} else{
				o2prime_pack_task->nonconst_residue_task( seq_num ).prevent_repacking();
			}
		}

		return o2prime_pack_task;
	}

	/////////////////////////////////////////////////////////////////////////////////////
	void
	print_backbone_torsions( pose::Pose const & pose, Size const five_prime_chainbreak ){

		using namespace core::id;

		conformation::Residue const & suite_lower_res = pose.residue( five_prime_chainbreak );
		TR << std::setw( 5 ) << " ep = " << std::setw( 15 ) << suite_lower_res.mainchain_torsion( 5 );
		TR << std::setw( 5 ) << " z = "  << std::setw( 15 ) << suite_lower_res.mainchain_torsion( 6 );


		Size const three_prime_chainbreak = five_prime_chainbreak + 1;

		if ( three_prime_chainbreak <= pose.total_residue() ) {
			conformation::Residue const & suite_upper_res = pose.residue( three_prime_chainbreak );
			TR << std::setw( 5 ) << " a = "  << std::setw( 15 ) << suite_upper_res.mainchain_torsion( 1 );
			TR << std::setw( 5 ) << " b = "  << std::setw( 15 ) << suite_upper_res.mainchain_torsion( 2 );
			TR << std::setw( 5 ) << " g = "  << std::setw( 15 ) << suite_upper_res.mainchain_torsion( 3 );
		}

		TR << std::endl;

	}
	/////////////////////////////////////////////////////////////////////////////////////


	//When a CUTPOINT_UPPER is added to 3' chain_break residue, the EXISTENCE of the CUTPOINT_UPPER atoms means that the alpha torsion which previously DOES NOT exist due to the chain_break now exist. The alpha value is automatically defined to the A-form value by Rosetta. However Rosetta does not automatically adjust the OP2 and OP1 atom position to account for this fact. So it is important that the OP2 and OP1 atoms position are correctly set to be consistent with A-form alpha torsion before the CUTPOINT_UPPER IS ADDED Parin Jan 2, 2009
	//Uncomment print_backbone_torsion on Dec 26, 2010
	void
	correctly_position_cutpoint_phosphate_torsions( pose::Pose & current_pose, Size const five_prime_chainbreak,  bool verbose /* = false*/ ){

		using namespace core::chemical;
		using namespace core::conformation;
		using namespace core::id;
		using namespace core::io::pdb;

		static const ResidueTypeSetCAP rsd_set = core::chemical::ChemicalManager::get_instance() ->
			residue_type_set(	core::chemical::RNA );

		if ( verbose ) output_title_text( "ENTER correctly_position_cutpoint_phosphate_torsions function", TR );

		chemical::AA res_aa = aa_from_name( "RAD" );
		ResidueOP new_rsd = conformation::ResidueFactory::create_residue( *( rsd_set->aa_map( res_aa )[1] ) ) ;
		if ( verbose ) TR << "  res_aa: " << res_aa << std::endl;

		Size three_prime_chainbreak = five_prime_chainbreak + 1;

		if ( verbose ){
			//dump_pdb(current_pose, "Before_prepending_dummy_nucleotide");
			TR <<  std::setw( 50 ) << "Before_prepending_dummy_nucleotide";
			print_backbone_torsions( current_pose, five_prime_chainbreak );
		}

		current_pose.prepend_polymer_residue_before_seqpos( *new_rsd, three_prime_chainbreak, true );
		chemical::rna::RNA_FittedTorsionInfo const rna_fitted_torsion_info;

		if ( verbose ){
			//dump_pdb(current_pose, "Before_setting_torsion_to_A_form.pdb");
			TR << std::setw( 50 ) << "Before_setting_torsion_to_A_form";
			print_backbone_torsions( current_pose, five_prime_chainbreak + 1 );
		}
		//Actually just by prepending the residue causes the alpha torsion to automatically be set to -64.0274,
		//so the manual setting below is actually not needed, May 24, 2010.. Parin S.
		//These are the initial value of virtual upper and lower cutpoint atom.
		//Actaully only the alpha (id::BB, 1) is important here since it set the position of O3' (LOWER) atom which in turn determines  OP2 and OP1 atom
		current_pose.set_torsion( TorsionID( three_prime_chainbreak + 1, id::BB, 1 ), -64.027359 );

		/* BEFORE AUG 24, 2011
		//Where the hell did I get these numbers from value...by appending with ideal geometry and look at the initalized value? Oct 13, 2009
		current_pose.set_torsion( TorsionID( five_prime_chainbreak + 1, id::BB, 5 ), -151.943 ); //Not Important?
		current_pose.set_torsion( TorsionID( five_prime_chainbreak + 1, id::BB, 6 ), -76.4185 ); //Not Important?
		current_pose.set_torsion( TorsionID( three_prime_chainbreak + 1, id::BB, 1 ), -64.0274 );
		*/

		//RAD.params
		//ICOOR_INTERNAL  LOWER  -64.027359   71.027062    1.593103   P     O5'   C5'
		//ICOOR_INTERNAL    OP2 -111.509000   71.937134    1.485206   P     O5' LOWER
		//ICOOR_INTERNAL    OP1 -130.894000   71.712189    1.485010   P     O5'   OP2

		//RCY.params
		//ICOOR_INTERNAL  LOWER  -64.027359   71.027062    1.593103   P     O5'   C5'
		//ICOOR_INTERNAL    OP2 -111.509000   71.937134    1.485206   P     O5' LOWER
		//ICOOR_INTERNAL    OP1 -130.894000   71.712189    1.485010   P     O5'   OP2

		//RGU.params
		//ICOOR_INTERNAL  LOWER  -64.027359   71.027062    1.593103   P     O5'   C5'
		//ICOOR_INTERNAL    OP2 -111.509000   71.937134    1.485206   P     O5' LOWER
		//ICOOR_INTERNAL    OP1 -130.894000   71.712189    1.485010   P     O5'   OP2

		//URA.parms
		//ICOOR_INTERNAL  LOWER  -64.027359   71.027062    1.593103   P     O5'   C5'
		//ICOOR_INTERNAL    OP2 -111.509000   71.937134    1.485206   P     O5' LOWER
		//ICOOR_INTERNAL    OP1 -130.894000   71.712189    1.485010   P     O5'   OP2

		if ( verbose ){
			//dump_pdb(current_pose, "After_setting_torsion_to_A_form.pdb");
			TR << std::setw( 50 ) << "After_setting_torsion_to_A_form"; print_backbone_torsions( current_pose, five_prime_chainbreak + 1 );
		}

		current_pose.delete_polymer_residue( five_prime_chainbreak + 1 );

		if ( verbose ){
			//dump_pdb(current_pose, "After_deleting_dummy_nucleotide");
			TR << std::setw( 50 ) << "After_deleting_dummy_nucleotide"; print_backbone_torsions( current_pose, five_prime_chainbreak );
		}

		if ( verbose ) output_title_text( "EXIT correctly_position_cutpoint_phosphate_torsions function", TR );

	}
	//////////////////////////////////////////////////////////////////////////

	void
	copy_torsions_FROM_TO( core::id::TorsionID const start_torsion_ID, core::id::TorsionID const end_torsion_ID, core::pose::Pose const & template_pose, core::pose::Pose & pose ){

 		using namespace core::chemical;
		using namespace core::conformation;
		using namespace core::id;
		using namespace core::chemical::rna;
		using namespace core::kinematics;

		if ( template_pose.total_residue() != pose.total_residue() ){
			TR << "template_pose.total_residue() = " << template_pose.total_residue() << " pose.total_residue() = " << pose.total_residue() << std::endl;
			utility_exit_with_message( "template_pose.total_residue() != pose.total_residue()" );
		}

		if ( ( template_pose.fold_tree() == pose.fold_tree() ) == false ){
			output_fold_tree_info( template_pose.fold_tree(), "template_pose", TR );
			output_fold_tree_info( pose.fold_tree(), "pose", TR );
			utility_exit_with_message( "( template_pose.fold_tree() == pose.fold_tree() ) == false" );
		}

		if ( start_torsion_ID.type() != id::BB ){
			TR << "start_torsion_ID: " << start_torsion_ID;
			utility_exit_with_message( "start_torsion_ID.type() != id::BB" );
		}

		if ( end_torsion_ID.type() != id::BB ){
			TR << "end_torsion_ID: " << start_torsion_ID;
			utility_exit_with_message( "end_torsion_ID.type() != id::BB" );
		}


		for ( Size seq_num = start_torsion_ID.rsd(); seq_num <= end_torsion_ID.rsd(); seq_num++ ){

			//Backbone
			for ( Size n = 1; n <= NUM_RNA_MAINCHAIN_TORSIONS; n++ ){

				if ( seq_num == start_torsion_ID.rsd() && n < start_torsion_ID.torsion() ) continue;

				if ( seq_num == end_torsion_ID.rsd() && n > end_torsion_ID.torsion() ) continue;

				pose.set_torsion( TorsionID( seq_num, id::BB, n ), template_pose.residue( seq_num ).mainchain_torsion( n ) );

			}

			//Side_chain
			bool copy_side_chain = true;

			if ( seq_num == start_torsion_ID.rsd() && start_torsion_ID.torsion() > 4 /*Delta*/ ) copy_side_chain = false;
			if ( seq_num == end_torsion_ID.rsd() && end_torsion_ID.torsion() < 4 /*Delta*/ ) copy_side_chain = false;

			if ( copy_side_chain ){
				pose.set_torsion( TorsionID( seq_num, id::CHI, 1 ), template_pose.residue( seq_num ).chi( 1 ) ); /*CHI*/
				pose.set_torsion( TorsionID( seq_num, id::CHI, 2 ), template_pose.residue( seq_num ).chi( 2 ) ); /*NU2*/
				pose.set_torsion( TorsionID( seq_num, id::CHI, 3 ), template_pose.residue( seq_num ).chi( 3 ) ); /*NU1*/
				pose.set_torsion( TorsionID( seq_num, id::CHI, 4 ), template_pose.residue( seq_num ).chi( 4 ) ); /*chi_O2prime*/

			}

		}
	}


	//////////////////////////////////////////////////////////////////////////
	core::Size
	setup_chain_break_jump_point( core::pose::Pose & pose, core::Size const jump_point_one, core::Size const jump_point_two, core::Size const five_prime_cutpoint, bool const verbose ){

		Size cutpoint = setup_bulge_jump_point( pose, jump_point_one, jump_point_two, verbose );

		correctly_position_cutpoint_phosphate_torsions( pose, five_prime_cutpoint, false /*verbose*/ );

		pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_LOWER, five_prime_cutpoint   );
		pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_UPPER, five_prime_cutpoint + 1 );

		return cutpoint;
	}

	//////////////////////////////////////////////////////////////////////////

	void
	remove_chain_break_jump_point( core::pose::Pose & pose, core::Size const five_prime_cutpoint, core::kinematics::FoldTree const fold_tree_without_cutpoint ){

		TR << "remove_chain_break_jump_point " << std::endl;

		pose::remove_variant_type_from_pose_residue( pose, chemical::CUTPOINT_LOWER, five_prime_cutpoint );
		pose::remove_variant_type_from_pose_residue( pose, chemical::CUTPOINT_UPPER, five_prime_cutpoint + 1 );

		pose.fold_tree( fold_tree_without_cutpoint );


	}

	//////////////////////////////////////////////////////////////////////////

	core::Size
	setup_bulge_jump_point( pose::Pose & pose, Size const & moving_base, Size const & reference_base, bool const verbose ){


		using namespace core::conformation;

		if ( moving_base == reference_base ){
			utility_exit_with_message( "moving_base == reference_base!" );
		}

		int i, j;

		Size cutpoint;
		if ( moving_base > reference_base ){
			i = reference_base;
			j = moving_base;
			cutpoint = moving_base - 1;
		} else{
			i = moving_base;
			j = reference_base;
			cutpoint = moving_base;
		}

		core::kinematics::FoldTree fold_tree = pose.fold_tree();	//HARD COPY?

		if ( verbose ) output_fold_tree_info( fold_tree, "Before add bulge jump point", TR );

		//		fold_tree.new_jump( five_prime_seq_num, three_prime_seq_num, cut_point );
		fold_tree.new_jump( reference_base, moving_base, cutpoint ); //Choose the residue five_prime of the actual cutpoint position

		if ( verbose ) TR << "after add new jump point" << std::endl;

		Residue const & rsd1( pose.residue( i ) );
		Residue const & rsd2( pose.residue( j ) );

		Size jump_num = 9999;

		for ( Size n = 1; n <= fold_tree.num_jump(); n++ ) {
			if ( ( fold_tree.upstream_jump_residue( n ) == i && fold_tree.downstream_jump_residue( n ) == j ) ||
				  ( fold_tree.upstream_jump_residue( n ) == j && fold_tree.downstream_jump_residue( n ) == i ) ){
				jump_num = n;
				break;
			}
		}

		//	fold_tree.set_jump_atoms( n, five_prime_seq_num, five_prime_atom, three_prime_seq_num, three_prime_atom);
		fold_tree.set_jump_atoms( jump_num, rsd1.atom_name( rsd1.chi_atoms( 1 )[4] ), rsd2.atom_name( rsd2.chi_atoms( 1 )[4] ) ); //Base atoms...

		if ( verbose ) output_fold_tree_info( fold_tree, "New fold_tree with bulge jump point", TR );
		pose.fold_tree( fold_tree );

		return cutpoint;
	}


	///////////////////////////////////////////////////////////////////////
	utility::vector1< bool >
	get_partition_definition( pose::Pose const & pose, Size const & moving_suite ){

		ObjexxFCL::FArray1D<bool> partition_definition( pose.total_residue(), false );
		pose.fold_tree().partition_by_residue( moving_suite, partition_definition );

		//silly conversion. There may be a faster way to do this actually.
		utility::vector1< bool > partition_definition_vector1;
		for ( Size n = 1; n <= pose.total_residue(); n++ )	partition_definition_vector1.push_back( partition_definition(n) );

		return partition_definition_vector1;

	}


	////////////////////////////////////////////////////////////////////////
	void
	apply_rotamer( pose::Pose & pose, utility::vector1< Torsion_Info > const & rotamer_list ){
		for ( Size i = 1; i <= rotamer_list.size(); i++ ) {
			pose.set_torsion( rotamer_list[ i ].id, rotamer_list[i].value );
		}
	}

	//////////////////////////////////////////////////////////////////////////////////////
/*
/// @brief Principal value of angle in degrees on ( -180, 180 ]
template < typename T >
inline
T
principal_angle_degrees( T const & angle )
{
	return remainder( angle, T( 360.0 ) );
}
*/
	bool
	check_for_messed_up_structure( core::pose::Pose const & pose, std::string const & tag ){
		using namespace core::scoring;
		using namespace core::chemical::rna;

		for ( Size seq_num = 1; seq_num <= pose.total_residue(); seq_num++ ){

			if ( pose.residue( seq_num ).aa() == core::chemical::aa_vrt ) continue; //Fang's electron density code

			//SML PHENIX conference
			if ( basic::options::option[basic::options::OptionKeys::rna::rna_prot_erraser].value() ){
				if ( !pose.residue( seq_num ).is_RNA() ) continue;
			}

			conformation::Residue const & rsd( pose.residue( seq_num ) );
			Real delta = numeric::principal_angle_degrees( rsd.mainchain_torsion( DELTA ) );
			Real chi = numeric::principal_angle_degrees( rsd.chi( 1 ) );
			Real nu_2 = numeric::principal_angle_degrees( rsd.chi( 2 ) );
			Real nu_1 = numeric::principal_angle_degrees( rsd.chi( 3 ) );


//			TR << " tag= " << tag << " seq_num= " << seq_num << " delta= " << delta << " chi= " << chi << " nu_2= " << nu_2 << " nu_1= " << nu_1 << std::endl;

			if ( ( delta >  - 0.01 && delta < 0.01 ) || ( nu_2 >  - 0.01 && nu_2 < 0.01 ) || ( nu_1 >  - 0.01 && nu_1 < 0.01 ) ){ //observation is that messed up structure will have delta value of zero
				TR << "Warning: " << tag << " is probably a messed up pose, will be ignored" << std::endl;
				TR << " seq_num = " << seq_num << " delta = " << delta << " chi = " << chi << " nu_2 = " << nu_2 << " nu_1 = " << nu_1 << std::endl;
				if ( ( rsd.has_variant_type( "VIRTUAL_RNA_RESIDUE" ) == true ) || ( rsd.has_variant_type( "VIRTUAL_RIBOSE" ) == true ) ){ //Implement on Oct 28,2010
					TR << "OK lets NOT ignore yet since this rsd has virtual_res or virtual_sugar variant type..continue checking.. " << std::endl;
				} else{
					return true;
				}
			}
		}
		return false;
	}

	//////////////////////////////////////////////////////////////////////////////////////
	core::Size
	get_residue_base_state( core::pose::Pose const & pose, Size const seq_num ){

		using namespace core::scoring;
		using namespace core::chemical::rna;

		Real const CHI_CUTOFF = 15.0; //Kinda RANDOM..ROUGH average between north chi_anti (~79) and north chi_syn (~-50)

		conformation::Residue const & rsd = pose.residue( seq_num );
		Real const chi = numeric::principal_angle_degrees( rsd.chi( CHI - NUM_RNA_MAINCHAIN_TORSIONS ) );

		if ( chi <= CHI_CUTOFF ){
			return SYN;
		} else{
			return ANTI;
		}

	}

	//////////////////////////////////////////////////////////////////////////////////////

	core::Size
	get_residue_pucker_state( core::pose::Pose const & pose, Size const seq_num, bool const verbose ){

		using namespace core::scoring;
		using namespace core::chemical::rna;

		static RNA_FittedTorsionInfo const rna_fitted_torsion_info;
		Real const DELTA_CUTOFF( rna_fitted_torsion_info.delta_cutoff() );

		if ( verbose ) TR.Debug << "  DELTA_CUTOFF angle = " << DELTA_CUTOFF;

		conformation::Residue const & rsd( pose.residue( seq_num ) );
		Real delta = numeric::principal_angle_degrees( rsd.mainchain_torsion( DELTA ) );

		if ( ( delta > 1.0 && delta < 179.00 ) == false ){
			TR.Debug << " seq_num = " << seq_num << " delta angle = " << delta << std::endl;

/////////////////////////
			if ( pose.residue( seq_num ).has_variant_type( "VIRTUAL_RNA_RESIDUE" ) ){
				TR.Debug << "Warning: delta angle is out of range for virtual_residue at seq_num " << seq_num << "!" << std::endl;
			} else{
				//This part is now obsolete .... Apr 30, 2010...
				//A possibility for a out of range delta is for imported pdbs (upper element of 1q93 for example).
				Real principal_delta = numeric::principal_angle_degrees( delta );

				//Consistency check
				if ( delta <  - 180 ){
					if ( ( delta + 359 < principal_delta ) == false || ( delta + 361 > principal_delta ) == false ){
						utility_exit_with_message( "delta <  - 180 but delta + 359 < principal_delta ) == false || ( delta + 361 > principal_delta ) == false!" );
					}
				}

				if ( delta > 180 ){
					if ( ( delta - 359 > principal_delta ) == false || ( delta - 361 < principal_delta ) == false ){
						utility_exit_with_message( "delta > 180 but delta - 359 > principal_delta ) == false || ( delta - 361 < principal_delta ) == false!" );
					}
				}

				delta = principal_delta;

				//Check again
				if ( ( delta > 1.0 && delta < 179.00 ) == false ) utility_exit_with_message( "principal delta angle out of range!" );
//////////////////////////
			}
		}

		if ( verbose ) TR.Debug << "  delta angle = " << delta << std::endl;

		if ( delta <= DELTA_CUTOFF ) {
			return NORTH;
		} else {
			return SOUTH;
		}
	}

	bool
	is_same_sugar_pucker( core::pose::Pose const & current_pose, core::pose::Pose const & cluster_center_pose, Size const seq_num ){
		if ( get_residue_pucker_state( current_pose, seq_num ) == get_residue_pucker_state( cluster_center_pose, seq_num ) ){
			return true;
		} else{
			return false;
		}
	}


	void
	sleep( core::Size mseconds )
	{
 	 clock_t endwait;
 	 endwait = clock () + mseconds * CLOCKS_PER_SEC/1000 ;
 	 while ( clock() < endwait ) {}
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	void
	setup_simple_fold_tree( core::pose::Pose & pose ){

//		using namespace core::chemical;

		Size const nres = pose.total_residue();

		kinematics::FoldTree simple_fold_tree( nres ); //Create a simple fold tree

		simple_fold_tree.simple_tree( nres ); //Just to make sure.

		pose.fold_tree( simple_fold_tree );
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////


	void
	get_atom_coordinates( utility::vector1< std::pair < id::AtomID, numeric::xyzVector< core::Real > > > & xyz_list, Size const & seq_num, core::conformation::Residue const & rsd_at_origin, core::kinematics::Stub const & moving_res_base_stub ){

		xyz_list.clear();

		numeric::xyzVector< core::Real > const & new_centroid = moving_res_base_stub.v;
		numeric::xyzMatrix< core::Real > const & new_coordinate_matrix = moving_res_base_stub.M;

		for ( Size at = 1; at <= rsd_at_origin.natoms(); at++ ){

			id::AtomID const id( at, seq_num );

			numeric::xyzVector< core::Real > atom_pos;

			atom_pos = new_coordinate_matrix * rsd_at_origin.xyz( at ); //I think the order here does matter.
			atom_pos = atom_pos + new_centroid; //I think the order here does matter.

			xyz_list.push_back( std::make_pair( id, atom_pos ) );
		}

	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	void
	import_pose_from_silent_file( core::pose::Pose & import_pose, std::string const & silent_file , std::string const & input_tag ){

		using namespace core::chemical;
		using namespace core::conformation;

		static const ResidueTypeSetCAP rsd_set = core::chemical::ChemicalManager::get_instance() ->
			residue_type_set( core::chemical::RNA );

		core::io::silent::SilentFileData silent_file_data;
		silent_file_data.read_file( silent_file );

		Size num_matching_tag = 0;

		for ( core::io::silent::SilentFileData::iterator iter = silent_file_data.begin(), end = silent_file_data.end(); iter != end; ++iter ){
			if ( iter->decoy_tag() != input_tag ) continue;
			num_matching_tag += 1;
			iter->fill_pose( import_pose, *rsd_set );
		}

		if ( num_matching_tag != 1 ){
			utility_exit_with_message( "num_matching_tag = ( " + ObjexxFCL::string_of( num_matching_tag ) + " ) != 1 for tag " + input_tag + " in silent file ( " + silent_file + " )!" );
		}

		if ( check_for_messed_up_structure( import_pose, input_tag ) == true ){
		 	utility_exit_with_message( "import_pose " + input_tag + " from silent_file " + silent_file + " is a messed up pose!" );
		}

	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	std::string
	path_basename( std::string const full_path ){

		size_t found = full_path.rfind( '/' );

		std::string basename;

		if ( found != std::string::npos ){
			basename = full_path.substr( found + 1 );
		} else {
			basename = full_path;
		}

		return basename;
 	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////


	bool
	is_residues_in_contact( core::Size const & res_ONE, core::pose::Pose const & pose_ONE, core::Size const & res_TWO, core::pose::Pose const & pose_TWO, core::Real const atom_atom_overlap_dist_cutoff, core::Size const num_atom_contacts_cutoff, bool const verbose ){

		using namespace ObjexxFCL;

		core::conformation::Residue const & rsd_ONE = pose_ONE.residue( res_ONE );
		core::conformation::Residue const & rsd_TWO = pose_TWO.residue( res_TWO );

		Size num_atom_contacts_so_far = 0;

		if ( num_atom_contacts_cutoff < 1 ){
			utility_exit_with_message( "num_atom_contacts_cutoff( " + string_of( num_atom_contacts_cutoff ) + " ) < 1!" );
		}

		for ( Size at_ONE = 1; at_ONE <= rsd_ONE.natoms(); at_ONE++ ){ //include hydrogen atoms
			for ( Size at_TWO = 1; at_TWO <= rsd_TWO.natoms(); at_TWO++ ){ //include hydrogen atoms

				if ( rsd_ONE.atom_type( at_ONE ).name() == "VIRT" ) continue;
				if ( rsd_TWO.atom_type( at_TWO ).name() == "VIRT" ) continue;

				Real const VDW_radius_ONE = rsd_ONE.atom_type( at_ONE ).lj_radius();
				Real const VDW_radius_TWO = rsd_TWO.atom_type( at_TWO ).lj_radius();

				Real const cutoff_sum_VDW_radius = VDW_radius_ONE + VDW_radius_TWO - atom_atom_overlap_dist_cutoff;

				if ( cutoff_sum_VDW_radius < 0 ) utility_exit_with_message( "( VDW_radius_ONE + VDW_radius_TWO - atom_atom_overlap_dist_cutoff ) < 0!!" );

				Real const atom_atom_dist_squared = ( rsd_ONE.xyz( at_ONE ) - rsd_TWO.xyz( at_TWO ) ).length_squared();

				if ( atom_atom_dist_squared < ( cutoff_sum_VDW_radius*cutoff_sum_VDW_radius ) ){
					num_atom_contacts_so_far++;
					if ( verbose ){
						TR << "res_ONE = " << res_ONE << " res_TWO = " << res_TWO << " num_atom_contacts_so_far = " << num_atom_contacts_so_far << "|";
						TR << " VDW_radius_ONE = " << VDW_radius_ONE << " VDW_radius_TWO = " << VDW_radius_TWO;
						TR << " atom_atom_overlap_dist_cutoff = " << atom_atom_overlap_dist_cutoff << " cutoff_sum_VDW_radius = " << cutoff_sum_VDW_radius;
						TR << " atom_atom_dist = " << ( rsd_ONE.xyz( at_ONE ) - rsd_TWO.xyz( at_TWO ) ).length();
						TR << " " << rsd_ONE.atom_name( at_ONE ) << " of res_ONE and " << rsd_TWO.atom_name( at_TWO ) << " of res_TWO are in contact! " << std::endl;
					}
				}

				if ( num_atom_contacts_so_far == num_atom_contacts_cutoff ){
					return true;
				}

				if ( num_atom_contacts_so_far > num_atom_contacts_cutoff ){ //consistency_check
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

		//Size const five_prime_res = job_parameters_->five_prime_chain_break_res();
		Size const three_prime_res = five_prime_res + 1;

		//Even through there is the chain_break, alpha of 3' and epl and gamma of 5' should be defined due to the existence of the upper and lower variant type atoms.

		for ( Size n = 1; n <= 3; n++ ){ //alpha, beta, gamma of 3' res
			pose.set_torsion( TorsionID( three_prime_res, id::BB,  n ), 0.0 );
		}

		for ( Size n = 5; n <= 6; n++ ){ //epsilon and zeta of 5' res
			pose.set_torsion( TorsionID( five_prime_res, id::BB,  n ), 0.0 );
		}

	}

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void
    get_default_allowed_bulge_res(
        utility::vector1< core::Size > & allow_bulge_res_list,
        core::pose::Pose const & pose,
        bool const verbose ){

        if ( allow_bulge_res_list.size() != 0 ){
            utility_exit_with_message( "allow_bulge_res_list.size() != 0" );
        }

        if ( verbose ){
            TR << "allow_bulge_res_list.size() == 0, ";
            TR << "Getting default_allowed_bulge_res!" << std::endl;
        }



    for ( Size seq_num = 1; seq_num <= pose.total_residue(); seq_num++ ){

        //exclude edge residues:
        if ( seq_num == 1 ) continue;

        if ( seq_num == pose.total_residue() ) continue;

        //bool is_cutpoint_closed=false;

        bool is_cutpoint_lower = pose.residue( seq_num ).has_variant_type(
                                 chemical::CUTPOINT_LOWER );

        bool is_cutpoint_upper = pose.residue( seq_num ).has_variant_type(
                                 chemical::CUTPOINT_UPPER );

        bool near_cutpoint_closed = is_cutpoint_lower || is_cutpoint_upper;

        bool near_cutpoint = pose.fold_tree().is_cutpoint( seq_num ) ||
                             pose.fold_tree().is_cutpoint( seq_num - 1 );


        bool near_cutpoint_open = near_cutpoint && !near_cutpoint_closed;

        if ( near_cutpoint_open ) continue;

        allow_bulge_res_list.push_back( seq_num );

    }

    }
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////


    core::Size
    virtualize_bulges( core::pose::Pose & input_pose,
        utility::vector1< core::Size > const & in_allow_bulge_res_list,
        core::scoring::ScoreFunctionOP const & scorefxn,
        std::string const & tag,
        bool const allow_pre_virtualize,
        bool const allow_consecutive_bulges,
        bool const verbose ){

        using namespace core::pose;
        using namespace core::scoring;
        using namespace ObjexxFCL;


        Size const total_res = input_pose.total_residue();

        Real const rna_bulge_bonus = ( scorefxn->get_weight( rna_bulge ) )*10;

        utility::vector1< core::Size > allow_bulge_res_list = in_allow_bulge_res_list;

        if ( allow_bulge_res_list.size() == 0 ){
            get_default_allowed_bulge_res( allow_bulge_res_list, input_pose, verbose );
        }


		if ( verbose ){
			TR << "Enter virtualize_bulges() " << std::endl;
			TR << "rna_bulge_bonus = " << rna_bulge_bonus << std::endl;
			output_boolean( "allow_pre_virtualize = ", allow_pre_virtualize, TR ); TR << std::endl;
			output_boolean( "allow_consecutive_bulges = ", allow_consecutive_bulges, TR ); TR << std::endl;
			output_seq_num_list( "allow_bulge_res_list = ", allow_bulge_res_list, TR );

			//Testing to see if checking in apply_virtual_rna_residue_variant_type can be violated!//////////
			pose::Pose testing_pose = input_pose;

			for ( Size seq_num = 1; seq_num <= total_res; seq_num++ ){
				if ( testing_pose.residue( seq_num ).aa() == core::chemical::aa_vrt ) continue; //Fang's electron density code
				if ( allow_bulge_res_list.has_value( seq_num ) == false ) continue;
				apply_virtual_rna_residue_variant_type( testing_pose, seq_num, true /*apply_check*/ );
			}
			//////////////////////////////////////////////////////////////////////////////////////////////////
		}

		//////////////////////////////////////////////////////////////////////////////////////////////////////
		if ( allow_pre_virtualize == false ){
			for ( Size seq_num = 1; seq_num <= total_res; seq_num++ ){
				if ( input_pose.residue( seq_num ).aa() == core::chemical::aa_vrt ) continue; //Fang's electron density code
				if ( input_pose.residue( seq_num ).has_variant_type( "VIRTUAL_RNA_RESIDUE" ) ){
					utility_exit_with_message( "allow_pre_virtualize == false but seq_num = " + string_of( seq_num ) + "  is already virtualized!!" );
				}
			}
		}
		//////////////////////////////////////////////////////////////////////////////////////////////////////
		pose::Pose working_pose = input_pose;
		Real const start_score = ( *scorefxn )( working_pose );


		Size num_res_virtualized = 0;
		Size round_num = 0;

		while ( true ){	//Need multiple round, since virtualizing a particular res will reduce the energy score of neighoring res and might effect wether neighoring res should be virtualized.
			pose::Pose base_pose = input_pose;
			Real const base_score = ( *scorefxn )( base_pose );

			round_num++;

			Size num_res_virtualized_in_this_round = 0;

			for ( Size seq_num = 1; seq_num <= total_res; seq_num++ ){
				if ( input_pose.residue( seq_num ).aa() == core::chemical::aa_vrt ) continue; //Fang's electron density code
				if ( allow_bulge_res_list.has_value( seq_num ) == false ) continue;

				if ( input_pose.residue( seq_num ).has_variant_type( "VIRTUAL_RNA_RESIDUE" ) ){
					if ( input_pose.residue( seq_num + 1 ).has_variant_type( "VIRTUAL_RNA_RESIDUE_UPPER" ) == false ){ //consistency_check
						utility_exit_with_message( "seq_num = " + string_of( seq_num ) + "  is a virtual res but seq_num + 1 is not a virtual_res_upper!" );
					}

					if ( base_pose.residue( seq_num ).has_variant_type( "VIRTUAL_RNA_RESIDUE" ) == false ){ //consistency check
						utility_exit_with_message( "input_pose have virtual at seq_num = " + string_of( seq_num ) + "  but input_pose doesn't!" );
					}

					continue;
				}

				if ( allow_consecutive_bulges == false ){
					if ( ( seq_num + 1 ) <= total_res ){
						if ( input_pose.residue( seq_num + 1 ).has_variant_type( "VIRTUAL_RNA_RESIDUE" ) ) continue;
					}

					if ( ( seq_num - 1 ) >= 1 ){
						if ( input_pose.residue( seq_num - 1 ).has_variant_type( "VIRTUAL_RNA_RESIDUE" ) ) continue;
					}
				}


				working_pose = base_pose; //reset working_pose to base_pose
				apply_virtual_rna_residue_variant_type( working_pose, seq_num, true )	;
				Real const new_score = ( *scorefxn )( working_pose );

				if ( new_score < base_score ){
					num_res_virtualized++;
					num_res_virtualized_in_this_round++;

					TR << "tag = " << tag << " round_num = " << round_num << " seq_num = " << seq_num << ". new_score ( " << new_score << " ) is lesser than base_score ( " << base_score << " ). " << std::endl;

					apply_virtual_rna_residue_variant_type( input_pose, seq_num, true )	;

				}
			}

			if ( num_res_virtualized_in_this_round == 0 ) break;

		}


		working_pose = input_pose;
		Real const final_score = ( *scorefxn )( working_pose );

		if ( num_res_virtualized > 0 ){
			TR << "----------------------------------------------------------" << std::endl;
			TR << "Inside virtualize_bulges() " << std::endl;
			TR << "TOTAL_NUM_ROUND = " << round_num << std::endl;
			TR << "tag = " << tag << std::endl;
			TR << "num_res_virtualized = " << num_res_virtualized << std::endl;
			TR << "start_score = " << start_score << std::endl;
			TR << "final_score = " << final_score << std::endl;
			TR << "----------------------------------------------------------" << std::endl;

		}

		return num_res_virtualized;

	}

	/////////////////New function on Nov 11, 2010///////////////

	std::string
	get_tag_from_pdb_filename( std::string const pdb_filename ){

		std::string tag;

		size_t found = pdb_filename.rfind( '/' );

		if ( found != std::string::npos ){
			tag = pdb_filename.substr( found + 1 );
		} else {
			tag = pdb_filename;
		}

		size_t found_2 = tag.rfind( ".pdb" );

		if ( found_2 != std::string::npos ){
			tag = tag.substr( 0, tag.size() - 4 );
		}

		return tag;
	}

	//DUPLICATE OF CODE IN StepWiseRNA_JobParametersSetup.cc
	void
	move_jump_atom_to_base( core::kinematics::FoldTree & fold_tree, std::string const & working_sequence ){

		Size const num_cutpoint = fold_tree.num_cutpoint();

		for ( Size i = 1; i <= num_cutpoint; i++ ) {
			Size const k = fold_tree.upstream_jump_residue( i );
			Size const m = fold_tree.downstream_jump_residue( i );

			char upstream_res = working_sequence[k - 1];
			char downstream_res = working_sequence[m - 1];

			//Base atoms...
			//chi_atoms(1)[4] )=  C2 if URA or RCY
			//chi_atoms(1)[4] )=  C4 if RGU or RAD
			std::string upstream_jump_atom;
			std::string downstream_jump_atom;

			if ( upstream_res == 'u' || upstream_res == 'c' ){
				upstream_jump_atom = " C2 ";
			} else if ( upstream_res == 'a' || upstream_res == 'g' ){
				upstream_jump_atom = " C4 ";
			} else{
				utility_exit_with_message( "Invalid upstream_res!!" );
			}

			if ( downstream_res == 'u' || downstream_res == 'c' ){
				downstream_jump_atom = " C2 ";
			} else if ( downstream_res == 'a' || downstream_res == 'g' ){
				downstream_jump_atom = " C4 ";
			} else{
				utility_exit_with_message( "Invalid downstream_res!!" );
			}

			TR << "upstream_res = " << k << upstream_res << " upstream_jump_atom = " << upstream_jump_atom;
			TR << " downstream_res = " << k << downstream_res << " downstream_jump_atom = " << downstream_jump_atom << std::endl;

			fold_tree.set_jump_atoms( i, downstream_jump_atom, upstream_jump_atom );

		}
	}


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	print_JobParameters_info( StepWiseRNA_JobParametersCOP const & const_JP, std::string const JP_name, std::ostream & outstream /* = std::cout */, bool const is_simple_full_length_JP  ){

		StepWiseRNA_JobParametersOP JP = new StepWiseRNA_JobParameters;

		( *JP ) = ( *const_JP );

		print_JobParameters_info( JP, JP_name, outstream, is_simple_full_length_JP );

	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	print_JobParameters_info( StepWiseRNA_JobParametersOP const & JP, std::string const JP_name, std::ostream & outstream /* = std::cout */, bool const is_simple_full_length_JP ){

		using namespace ObjexxFCL;

		output_title_text( "printing JobParameters Information for " + JP_name, outstream );

		utility::vector1< Size > empty_seq_num_list;
		empty_seq_num_list.clear();

		outstream << "full_sequence = " <<  JP->full_sequence();
		outstream << " working_sequence = " <<  JP->working_sequence();
		outstream << " moving_res = " <<  JP->moving_res();
		outstream << " working_moving_res = " << JP->working_moving_res();
		outstream << " working_moving_suite = " <<  JP->working_moving_suite() << std::endl;


		outstream << "gap_size = " << JP->gap_size();
		outstream << " five_prime_chain_break_res = " << JP->five_prime_chain_break_res();
		output_boolean( " is_prepend = ",  JP->is_prepend(), outstream ) ;
		output_boolean( " is_internal = ", JP->is_internal(), outstream );
		output_boolean( " output_extra_RMSDs = ", JP->	output_extra_RMSDs(), outstream ); outstream << std::endl;


		//std::map< core::Size, core::Size > full_to_sub_;
		//std::map< core::Size, core::Size > sub_to_full_;

		//utility::vector1< std::pair< core::Size, core::Size > > chain_boundaries_;

		//	ObjexxFCL::FArray1D< bool > partition_definition_;

		//core::pose::PoseOP working_native_pose_;

		outstream << "------------full_stuff------------" << std::endl;
		output_bool_list( "is_working_res = ", JP->is_working_res(), outstream );

		for ( Size n = 1; n <= JP->input_res_vectors().size(); n++ ){
			output_seq_num_list( "input_res_vectors[" + string_of( n ) + "]", JP->input_res_vectors()[n], outstream );
		}
		output_seq_num_list( "global_sample_res_list = ", JP->global_sample_res_list(), outstream );

		if ( JP->force_syn_chi_res_list().size() > 0 )  output_seq_num_list( "force_syn_chi_res_list = ", JP->force_syn_chi_res_list(), outstream );
		if ( JP->force_north_sugar_list().size() > 0 ) output_seq_num_list( "force_north_sugar_list = ", JP->force_north_sugar_list(), outstream );
		if ( JP->force_south_sugar_list().size() > 0 ) output_seq_num_list( "force_south_sugar_list = ", JP->force_south_sugar_list(), outstream );
		if ( JP->protonated_H1_adenosine_list().size() > 0 ) output_seq_num_list( "protonated_H1_adenosine_list ", JP->protonated_H1_adenosine_list(), outstream );





		output_is_prepend_map( "is_prepend_map = ", JP->is_prepend_map(), JP->full_sequence().size(), outstream );
		output_seq_num_list( "rmsd_res_list = ", 								JP->rmsd_res_list() 								, outstream );

		output_seq_num_list( "native_alignment = ",  							JP->native_alignment(), outstream );
		output_seq_num_list( "cutpoint_closed_list = ", 					JP->cutpoint_closed_list(), outstream );


		outstream << "------------working_stuff------------" << std::endl;

		output_seq_num_list( "working_global_sample_res_list = ", JP->working_global_sample_res_list(), outstream );
		if ( JP->force_syn_chi_res_list().size() > 0 )  output_seq_num_list( "working_force_syn_chi_res_list = ", JP->working_force_syn_chi_res_list(), outstream );
		if ( JP->force_north_sugar_list().size() > 0 ) output_seq_num_list( "working_force_north_sugar_list = ", JP->working_force_north_sugar_list(), outstream );
		if ( JP->force_south_sugar_list().size() > 0 ) output_seq_num_list( "working_force_south_sugar_list = ", JP->working_force_south_sugar_list(), outstream );
		if ( JP->protonated_H1_adenosine_list().size() > 0 ) output_seq_num_list( "working_protonated_H1_adenosine_list = ", JP->working_protonated_H1_adenosine_list(), outstream );


		output_seq_num_list( "working_fixed_res = ",						JP->working_fixed_res() 						, outstream );

		if ( is_simple_full_length_JP == false ){
			output_seq_num_list( "working_moving_res_list = ", JP->working_moving_res_list(), outstream );
			output_seq_num_list( "working_moving_suite_list = ", JP->working_moving_suite_list(), outstream );
		} else{
			output_seq_num_list( "line_filler ", empty_seq_num_list, outstream );
			output_seq_num_list( "line_filler ", empty_seq_num_list, outstream );
		}

		output_seq_num_list( "working_terminal_res = ", 					JP->working_terminal_res() 					, outstream );
		output_seq_num_list( "working_moving_partition_pos = ",  JP->working_moving_partition_pos(), outstream );

		output_seq_num_list( "working_best_alignment = ", 				JP->working_best_alignment(), outstream );
		output_seq_num_list( "working_native_alignment = ", 			JP->working_native_alignment(), outstream );

		if ( is_simple_full_length_JP == false ){
			utility::vector1< bool > vector1_partition_definition;

			for ( Size n = 1; n <= JP->partition_definition().size(); n++ ){
				vector1_partition_definition.push_back( JP->partition_definition()( n ) );
			}

			output_bool_list( "partition_definition = ", 	vector1_partition_definition, outstream );

			outstream << "root_res = " << JP->fold_tree().root() << std::endl;
			output_fold_tree_info( JP->fold_tree(), "fold_tree", outstream );
		}

		output_title_text( "", outstream );


	}


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
	void
	process_mem_usage( double& vm_usage, double& resident_set, core::Size count ){
		 using std::ios_base;
		 using std::ifstream;
		 using std::string;

		 vm_usage     = 0.0;
		 resident_set = 0.0;

		 // 'file' stat seems to give the most reliable results
		 //
		 ifstream stat_stream( "/proc/self/stat", ios_base::in );

		 // dummy vars for leading entries in stat that we don't care about
		 //
		 string pid, comm, state, ppid, pgrp, session, tty_nr;
		 string tpgid, flags, minflt, cminflt, majflt, cmajflt;
		 string utime, stime, cutime, cstime, priority, nice;
		 string O, itrealvalue, starttime;

		 // the two fields we want
		 //
		 unsigned long vsize;
		 long rss;

		 stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
		             >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
		             >> utime >> stime >> cutime >> cstime >> priority >> nice
		             >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

		 long page_size_kb = sysconf( _SC_PAGE_SIZE ) / 1024; // in case x86-64 is configured to use 2MB pages
		 vm_usage     = vsize / 1024.0;
		 resident_set = rss * page_size_kb;


	 TR << "memory_usage = " << count << " " << pid << " " << comm << " " << state << " " << ppid << " " << pgrp << " " << session << " " << tty_nr << " " << tpgid << " " << flags << " " << minflt <<  " " <<  cminflt << " " << majflt << " " << cmajflt << " " << utime << " " << stime << " " << cutime << " " << cstime << " " << priority << " " << nice << " " <<  O << " " <<  itrealvalue << " " <<  starttime << " " <<  vsize << " " << rss << std::endl;

	}
*/

/*
	void
	add_native_base_pair_stats( core::io::silent::SilentStructOP & s, pose::Pose const & native_pose, utility::vector1 < core::Size > const & input_rmsd_res_list ){

		using namespace protocols::rna;
		using namespace conformation;
		using namespace scoring::rna;
		//using namespace core::pose;

		utility::vector1 < core::Size > act_rmsd_res_list;
		act_rmsd_res_list.clear();

		for ( Size n = 1; n <= input_rmsd_res_list.size(); n++ ){

			Size const seq_num = input_rmsd_res_list[n];

			if ( is_virtual_base( native_pose.residue( seq_num ) ) ) continue;

			act_rmsd_res_list.push_back( seq_num );

		}

		utility::vector1< core::scoring::rna::Base_pair > const base_pair_list = classify_base_pairs_parin( native_pose, act_rmsd_res_list );

		Size native_WC( 0 ), native_NWC( 0 );

		for ( Size n = 1; n <= base_pair_list.size(); n++ ) {

			core::scoring::rna::Base_pair const base_pair = base_pair_list[ n ];

			Residue const & rsd_1 = native_pose.residue( base_pair.res1 );
			Residue const & rsd_2 = native_pose.residue( base_pair.res2 );

			if ( ( base_pair.edge1 == WATSON_CRICK && base_pair.edge2 == WATSON_CRICK && base_pair.orientation == 1 )  && core::scoring::rna::possibly_canonical( rsd_1.aa(), rsd_2.aa() ) ){
				native_WC++;
			} else {
				native_NWC++;
			}
		}
		s->add_energy( "NAT_WC", native_WC );
		s->add_energy( "NAT_NWC", native_NWC );

	}
*/

/* This requires my version of the RNA_BasePairClassifier.cc which I decide not to commit to trunk ( Parin Sripakdeevong 27/12/2011 )
	void
	add_base_pair_stats( core::io::silent::SilentStructOP & s, pose::Pose const & pose, pose::Pose const & native_pose, utility::vector1 < core::Size > const & input_rmsd_res_list ){


		using namespace scoring::rna;
		using namespace core::conformation;
		using namespace protocols::rna;

		bool verbose = false;

		utility::vector1 < core::Size > act_rmsd_res_list;
		act_rmsd_res_list.clear();

		for ( Size n = 1; n <= input_rmsd_res_list.size(); n++ ){

			Size const seq_num = input_rmsd_res_list[n];

			if ( is_virtual_base( native_pose.residue( seq_num ) ) ) continue;

			act_rmsd_res_list.push_back( seq_num );

		}

		utility::vector1< core::scoring::rna::Base_pair > const model_base_pair_list = classify_base_pairs_strict( pose, act_rmsd_res_list );

		utility::vector1< core::scoring::rna::Base_pair > const native_base_pair_list = classify_base_pairs_strict( native_pose, act_rmsd_res_list );

		//TR << "BLAH" << std::endl;

		Size native_WC = 0;
		Size native_NWC = 0;
		Size recovered_WC = 0;
		Size recovered_NWC = 0;


		for ( Size n = 1; n <= native_base_pair_list.size(); n++ ) {

			core::scoring::rna::Base_pair const native_base_pair = native_base_pair_list[ n ];

			Residue const & rsd_1( native_pose.residue( native_base_pair.res1 ) );
			Residue const & rsd_2( native_pose.residue( native_base_pair.res2 ) );

			if ( ( native_base_pair.edge1 == WATSON_CRICK && native_base_pair.edge2 == WATSON_CRICK && native_base_pair.orientation == 1 )  &&  possibly_canonical( rsd_1.aa(), rsd_2.aa() ) )		{
				native_WC++;

				if ( check_in_base_pair_list( native_base_pair, model_base_pair_list ) ) recovered_WC++;

			} else {
				native_NWC++;

				if ( check_in_base_pair_list( native_base_pair, model_base_pair_list ) ){
					recovered_NWC++;
				} else {
					if ( verbose ){
						TR << "Missing native base pair " << pose.residue( native_base_pair.res1 ).name1() << " " << pose.residue( native_base_pair.res2 ).name1() << "  ";
						native_base_pair.print_info(); TR << std::endl;
					}
				}
			}
		}

		s->add_energy( "NAT_WC", native_WC );
		s->add_energy( "NAT_NWC", native_NWC );

		s->add_energy( "REC_WC", recovered_WC );
		s->add_energy( "REC_NWC", recovered_NWC );

	}
*/

	void
	set_nucleotide_to_A_form( pose::Pose & pose, Size const seq_num ){
		//Torsion value extracted from 3DNA (webiste) (A-U BP repeating) idealized A-form helix. Note that bond angle and bond length of idealized Rosetta doesn't exactly match the values in 3DNA

		using namespace core::id;

		pose.set_torsion( TorsionID( seq_num, id::BB,  1 ), -68.9 ); //alpha
		pose.set_torsion( TorsionID( seq_num, id::BB,  2 ), 179.5 ); //beta
		pose.set_torsion( TorsionID( seq_num, id::BB,  3 ), 54.5 ); //gamma
		pose.set_torsion( TorsionID( seq_num, id::BB,  5 ), -154.0 ); //epsilon
		pose.set_torsion( TorsionID( seq_num, id::BB,  6 ), -70.8 ); //zeta

		pose.set_torsion( TorsionID( seq_num, id::BB,  4 ), 82.2 ); //delta
		pose.set_torsion( TorsionID( seq_num, id::CHI, 1 ), 79.2 ); //chi
		pose.set_torsion( TorsionID( seq_num, id::CHI, 2 ), 36.9 ); //nu2
		pose.set_torsion( TorsionID( seq_num, id::CHI, 3 ), 94.7 ); //nu1


	}
/////////////////////////////////////////////////////////////////////////////////////////////
	void
	print_atom_info( pose::Pose const & pose, Size const seq_num, std::string const pose_name ){
		TR << "print_atom_info for pose: " << pose_name << " seq_num = " << seq_num << std::endl;

		conformation::Residue const & rsd = pose.residue( seq_num ); //static_pose

		for ( Size at = 1; at <= rsd.natoms(); at++ ){ //I wonder if we should just consider heavy atom? (rsd_1.nheavyatoms())

			TR << "atom = " << at  << "|name = " << rsd.type().atom_name( at ) << "|type = " << rsd.atom_type( at ).name();
			TR << "|element() = " << rsd.atom_type( at ).element() << "|" << std::endl;

		}

	}
/////////////////////////////////////////////////////////////////////////////////////////////
	void
	print_individual_atom_info( core::conformation::Residue const & rsd, Size const atomno, std::string const rsd_name ){
		TR << "individual_atom_info: rsd_name " << rsd_name;

		TR << " atom = " << atomno  << "|name = " << rsd.type().atom_name( atomno ) << "|type = " << rsd.atom_type( atomno ).name();
		TR << "|element() = " << rsd.atom_type( atomno ).element() << "|" << std::endl;

	}

/////////////////////////////////////////////////////////////////////////////////////////////

	void
	print_base_state( std::string const tag, core::Size const base_state, std::ostream & outstream /* = std::cout */ ){

		std::string base_state_string = "";

		if ( base_state == WHATEVER ){
			base_state_string = "WHATEVER";
		} else if ( base_state == ANTI ){
			base_state_string = "ANTI";
		} else if ( base_state == SYN ){
			base_state_string = "SYN" ;
		} else if ( base_state == NONE ){
			base_state_string = "NONE";
		} else{
			outstream << "Invalid base state = " << base_state << std::endl;
			utility_exit_with_message( "Invalid base state!" );
		}

		outstream << tag << std::setw( 4 ) << std::left <<  base_state_string << " ";
	}

/////////////////////////////////////////////////////////////////////////////////////////////

	void
	print_sugar_pucker_state( std::string const tag, core::Size const pucker_state, std::ostream & outstream /* = std::cout */ ){

		std::string pucker_state_string = "";

		if ( pucker_state == NORTH ){
			pucker_state_string = "NORTH";
		} else if ( pucker_state == SOUTH ){
			pucker_state_string = "SOUTH";
		} else if ( pucker_state == WHATEVER ){
			pucker_state_string = "WHATEVER";
		} else{
			outstream << "Invalid pucker state = " << pucker_state << std::endl;
			utility_exit_with_message( "Invalid pucker state!" );
		}

		outstream << tag << std::setw( 5 ) << std::left << pucker_state_string << " ";

	}

	//////////July 20, 2011, common scorefunctions..Used to be part of StepWiseRNA_ResidueSampler.cc
	void
	initialize_common_scorefxns( core::scoring::ScoreFunctionOP const & scorefxn_,
													 	core::scoring::ScoreFunctionOP & sampling_scorefxn_,
													 	core::scoring::ScoreFunctionOP & atr_rep_screening_scorefxn_,
													 	core::scoring::ScoreFunctionOP & chainbreak_scorefxn_,
													 	core::scoring::ScoreFunctionOP & o2prime_pack_scorefxn_ ){


		using namespace core::scoring;

		///////////////////////////////////////////////////////////////////
		// Bare minimum to check for contact (fa_atr) but not clash (fa_rep)
		atr_rep_screening_scorefxn_ =  new ScoreFunction;
		atr_rep_screening_scorefxn_->set_weight( fa_atr , 0.23 );
		atr_rep_screening_scorefxn_->set_weight( fa_rep , 0.12 );

		///////////////////////////////////////////////////////////////////
		chainbreak_scorefxn_ =  new ScoreFunction;
		chainbreak_scorefxn_->set_weight( angle_constraint, 1.0 );
		chainbreak_scorefxn_->set_weight( atom_pair_constraint, 1.0 );

		////////////////////Setup sampling scoring//////////////////////////////////////////////////////////////////////////////
    //1. Want to increase fa_rep during the minimization phase but want to keep it at 0.12 during the sample phase
	  //2. Sugar scoring is always turned off during sampling stage since it screw up pose selection. (TURN IT BACK ON: RD 01/31/2010)
		//3. Harmonic and Linear Chain_break scoring is always turned off during sampling stage
		sampling_scorefxn_ = scorefxn_->clone();

		//		sampling_scorefxn_->set_weight( rna_sugar_close, 0.0 ); (TURN IT BACK ON: RD 01/31/2010)

		sampling_scorefxn_->set_weight( fa_rep, 0.12 );
		//Only important only if fa_rep score in weight file is not 0.12..want to make sure that changing fa_rep in weight file doesn't effect sampling process. May 23 2010 Parin S.

		sampling_scorefxn_->set_weight( linear_chainbreak, 0.0 );
		sampling_scorefxn_->set_weight( chainbreak, 0.0 );
		sampling_scorefxn_->set_weight( angle_constraint, 0.0 );
		sampling_scorefxn_->set_weight( atom_pair_constraint, 0.0 );
		//This makes sure that there are no chain_break score involve in the full_score screening.
		//This works since by the time a pose reach full_score screening, it must already pass chain_break screening, May 23, 2010 Parin S.



		//////////////////////////////////////////////TESTING/////////////////////////////////////////////////
		//////////////////////////////////////////////TESTING/////////////////////////////////////////////////
		// Note that: rna_torsion, rna_sugar_close, fa_stack not optimized -- also irrelevant for 2'-OH sampling.
		//just a comparison. This is extremely slow. Would need to implement trie
		// for geom_sol, lk_nonpolar, and elec... Not too hard, but I don't feel like doing it now. (Rhiju)
		//o2prime_pack_scorefxn_ = sampling_scorefxn_->clone();
		/*
		//Parin July 21, 2011
		o2prime_pack_scorefxn_ = new ScoreFunction;
		// Each of the following terms have been pretty optimized for the packer (trie, etc.)
		o2prime_pack_scorefxn_->set_weight( fa_atr, sampling_scorefxn_->get_weight( fa_atr ) );
		o2prime_pack_scorefxn_->set_weight( fa_rep, sampling_scorefxn_->get_weight( fa_rep ) );
		o2prime_pack_scorefxn_->set_weight( hbond_lr_bb_sc, sampling_scorefxn_->get_weight( hbond_lr_bb_sc ) );
		o2prime_pack_scorefxn_->set_weight( hbond_sr_bb_sc, sampling_scorefxn_->get_weight( hbond_sr_bb_sc ) );
		o2prime_pack_scorefxn_->set_weight( hbond_sc, sampling_scorefxn_->get_weight( hbond_sc ) );
		o2prime_pack_scorefxn_->set_energy_method_options( sampling_scorefxn_->energy_method_options() );
		o2prime_pack_scorefxn_->set_weight( lk_nonpolar, sampling_scorefxn_->get_weight( lk_nonpolar ) ); //FIX THIS ON JULY 21, 2011. TOO EXPENSIVE!
		o2prime_pack_scorefxn_->set_weight( lk_polar, sampling_scorefxn_->get_weight( geom_sol ) ); //FIX THIS ON JULY 21, 2011. TOO EXPENSIVE!
		// note that geom_sol is not optimized well --> replace with lk_sol for now.
		*/
		//////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////


		o2prime_pack_scorefxn_ = new ScoreFunction;
		// Each of the following terms have been pretty optimized for the packer (trie, etc.)
		o2prime_pack_scorefxn_->set_weight( fa_atr, sampling_scorefxn_->get_weight( fa_atr ) );
		o2prime_pack_scorefxn_->set_weight( fa_rep, sampling_scorefxn_->get_weight( fa_rep ) );
		o2prime_pack_scorefxn_->set_weight( hbond_lr_bb_sc, sampling_scorefxn_->get_weight( hbond_lr_bb_sc ) );
		o2prime_pack_scorefxn_->set_weight( hbond_sr_bb_sc, sampling_scorefxn_->get_weight( hbond_sr_bb_sc ) );
		o2prime_pack_scorefxn_->set_weight( hbond_sc, sampling_scorefxn_->get_weight( hbond_sc ) );
		///Warning, don't include hbond_intra, since hbond_intra HAS NOT been been optimized for packing!

		o2prime_pack_scorefxn_->set_weight( fa_sol, sampling_scorefxn_->get_weight( lk_nonpolar ) ); //// note that geom_sol is not optimized well --> replace with lk_sol for now. //IS THIS A MISTAKE???

		o2prime_pack_scorefxn_->set_energy_method_options( sampling_scorefxn_->energy_method_options() ); //This set NO_HB_ENV_DEP, INCLUDE_INTRA_RES_RNA_HB and etcs.

	}

	void
	copy_all_o2prime_torsions( core::pose::Pose & mod_pose, core::pose::Pose const & template_pose ){

		using namespace core::id;
		using namespace core::conformation;

		if ( template_pose.total_residue() != mod_pose.total_residue() ){
			utility_exit_with_message( "template_pose.total_residue() != mod_pose.total_residue()" );
		}

		utility::vector1< bool > do_update_list( template_pose.total_residue(), false );

		for ( Size seq_num = 1; seq_num <= template_pose.total_residue(); seq_num++ ){
			if ( template_pose.residue( seq_num ).aa() == core::chemical::aa_vrt ) continue;	//Fang's electron density code
			if ( std::abs( template_pose.torsion( TorsionID( seq_num, id::CHI, 4 ) ) - mod_pose.torsion( TorsionID( seq_num, id::CHI, 4 ) ) ) > 0.001 ){
				//the two o2prime torsions are not the same! Basically don't want to trigger refolding if need not to.
				do_update_list[seq_num] = true;
			}
		}

		for ( Size seq_num = 1; seq_num <= template_pose.total_residue(); seq_num++ ){
			if ( mod_pose.residue( seq_num ).aa() == core::chemical::aa_vrt ) continue;	//Fang's electron density code
			if ( do_update_list[seq_num] == true ){
				mod_pose.set_torsion( TorsionID( seq_num, id::CHI, 4 ), template_pose.torsion( TorsionID( seq_num, id::CHI, 4 ) ) );
			}
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

	for ( Size n = 1; n <= n_score_types; n++ ){

		core::Real const old_weight = rescaled_scorefxn->get_weight( ScoreType( n ) );

		if ( old_weight != 0.0 ) {
			non_zero_weight += 1;
			rescaled_scorefxn->set_weight( ScoreType( n ), old_weight*scaling_factor );
		}
	}

	TR << std::endl;

	TR << "n_score_types = " << int( n_score_types ) << " non_zero_weight = " << non_zero_weight << " scaling_factor = " << scaling_factor << std::endl;

	show_scorefxn_weight_lines( rescaled_scorefxn, "AFTER REWEIGHT" );

	return rescaled_scorefxn;

}

/////////////////////////////////////////////////////////////////////////////////////////////

void
show_scorefxn_weight_lines( core::scoring::ScoreFunctionOP const & scorefxn, std::string const title ){

	using namespace core::scoring;
	using namespace ObjexxFCL::format;

	TR << "----------------" << title << "----------------" << std::endl;

	TR << "----------------------------------------------\n";
	TR << " Scores                             Weight\n";
	TR << "----------------------------------------------\n";

	pose::Pose empty_pose = *( new pose::Pose );

	for ( Size n = 1; n <= n_score_types; n++ ){

		core::Real const weight = scorefxn->get_weight( ScoreType( n ) );

		if ( weight != 0.0 ){
			TR << ' ' << LJ( 30, ScoreType( n ) ) << ' ' << F( 11, 5, weight ) << std::endl;
		}
	}

	std::string dash_string = "";
	for ( Size n = 1; n <= title.size(); n++ ){
		dash_string += "-";
	}

	TR << "----------------" << dash_string << "----------------" << std::endl;

}

	/////////////////////////////////////////////////////////////////////////////////////////////
	void
	figure_out_swa_rna_movemap( core::kinematics::MoveMap & mm, core::pose::Pose const & pose, utility::vector1< Size > const & minimize_res ){

		Size const nres( pose.total_residue() );
		ObjexxFCL::FArray1D < bool > allow_insert( nres, false );
		for ( Size i = 1; i <= minimize_res.size(); i++ ) allow_insert( minimize_res[ i ] ) = true;
		figure_out_swa_rna_movemap( mm, pose, allow_insert );

	}


	/////////////////////////////////////////////////////////////////////////////////////////////
	void
	figure_out_swa_rna_movemap( core::kinematics::MoveMap & mm, core::pose::Pose const & pose, ObjexxFCL::FArray1D < bool > const & allow_insert ){

		using namespace core::id;
		using namespace core::scoring::rna;
		using namespace core::chemical;

		Size const nres = pose.total_residue();
		runtime_assert( nres == allow_insert.size() );

		mm.set_bb( false );
		mm.set_chi( false );
		mm.set_jump( false );

		// New -- rhiju, june 2013
		for ( Size i = 1; i <= nres; i++ ){
			if ( allow_insert( i ) ) {
				mm.set_bb( i, true );
				mm.set_chi( i, true );
			}
		}

		for ( Size i = 1; i <= nres; i++ ){

			if ( pose.residue( i ).aa() == core::chemical::aa_vrt ) continue; //Fang's electron density code.

			//SML PHENIX conference
			if ( basic::options::option[basic::options::OptionKeys::rna::rna_prot_erraser].value() ){
				if ( !pose.residue( i ).is_RNA() ) continue;
			}

			utility::vector1< TorsionID > torsion_ids;

			for ( Size rna_torsion_number = 1; rna_torsion_number <= NUM_RNA_MAINCHAIN_TORSIONS; rna_torsion_number++ ) {
				torsion_ids.push_back( TorsionID( i, id::BB, rna_torsion_number ) );
			}
			for ( Size rna_torsion_number = 1; rna_torsion_number <= NUM_RNA_CHI_TORSIONS; rna_torsion_number++ ) {
				torsion_ids.push_back( TorsionID( i, id::CHI, rna_torsion_number ) );
			}


			for ( Size n = 1; n <= torsion_ids.size(); n++ ) {

				TorsionID const & torsion_id  = torsion_ids[ n ];

				id::AtomID id1, id2, id3, id4;
				bool fail = pose.conformation().get_torsion_angle_atom_ids( torsion_id, id1, id2, id3, id4 );
				if ( fail ) continue; //This part is risky, should also rewrite...

				// Dec 19, 2010..Crap there is a mistake here..should have realize this earlier...
				//Should allow torsions at the edges to minimize...will have to rewrite this. This effect the gamma and beta angle of the 3' fix res.

				// If there's any atom that is in a moving residue by this torsion, let the torsion move.
				//  should we handle a special case for cutpoint atoms? I kind of want those all to move.
				utility::vector1< AtomID > torsion_atom_ids = utility::tools::make_vector1( id1, id2, id3, id4 );

				for ( Size k = 1; k <= torsion_atom_ids.size(); k++ ){
					if ( allow_insert( torsion_atom_ids[k].rsd() ) ) {
						mm.set(  torsion_id, true );
						break;
					}
				}
				if ( mm.get( torsion_id ) ) continue;


				//
				// there is a note above from parin 'Should allow torsions at the edges to minimize...'
				// it appears fixed, except for cutpoints. Has this caused a problem in SWA & ERRASER at closed cutpoints?
				// Following code introduced by rhiju on aug. 2013:
				//
				if ( pose.residue(i).has_variant_type( CUTPOINT_LOWER ) && allow_insert( i+1 ) ){
					for ( Size k = 1; k <= torsion_atom_ids.size(); k++ ){
						if ( pose.residue_type( torsion_atom_ids[k].rsd() ).atom_name( torsion_atom_ids[k].atomno() ) == "OVL1" ||
								 pose.residue_type( torsion_atom_ids[k].rsd() ).atom_name( torsion_atom_ids[k].atomno() ) == "OVL2" )  {
							mm.set(  torsion_id, true );
							break;
						}
					}
				}
				if ( mm.get( torsion_id ) ) continue;

				if ( pose.residue(i).has_variant_type( CUTPOINT_UPPER ) && allow_insert( i-1 ) ){
					for ( Size k = 1; k <= torsion_atom_ids.size(); k++ ){
						if ( pose.residue_type( torsion_atom_ids[k].rsd() ).atom_name( torsion_atom_ids[k].atomno() ) == "OVU1" ){
							mm.set(  torsion_id, true );
							break;
						}
					}
				}
				if ( mm.get( torsion_id ) ) continue;

			}
		}


		// why is this in the internal loop? -- rhiju
		TR.Debug << "pose.fold_tree().num_jump() = " << pose.fold_tree().num_jump() << std::endl;

		for ( Size n = 1; n <= pose.fold_tree().num_jump(); n++ ){
			Size const jump_pos1( pose.fold_tree().upstream_jump_residue( n ) );
			Size const jump_pos2( pose.fold_tree().downstream_jump_residue( n ) );

			if ( pose.residue( jump_pos1 ).aa() == core::chemical::aa_vrt ) continue; //Fang's electron density code
			if ( pose.residue( jump_pos2 ).aa() == core::chemical::aa_vrt ) continue; //Fang's electron density code

			//SML PHENIX conference
			if ( basic::options::option[basic::options::OptionKeys::rna::rna_prot_erraser].value() ){
				if ( !pose.residue( jump_pos1 ).is_RNA() ) continue;
				if ( !pose.residue( jump_pos2 ).is_RNA() ) continue;
			}

			bool const move_jump = allow_insert( jump_pos1 ) || allow_insert( jump_pos2 );
			if ( move_jump )	mm.set_jump( n, true );
			TR.Debug << "jump_pos1 = " << jump_pos1 << " jump_pos2 = " << jump_pos2 << " mm.jump = "; output_boolean( move_jump, TR.Debug );  TR.Debug << std::endl;

		}

	}


	////////////////////////////////////////////////////////////////////
	void
	choose_random_if_unspecified_nucleotide( char & newrestype ) {

		std::string const rna_chars = "acgu";

		if ( newrestype == 'n' ){
			newrestype = rna_chars[ RG.random_range( 1, rna_chars.size() ) - 1 ];
			TR << "Choosing random nucleotide: " << newrestype << std::endl;
		}

	}



	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	mutate_res_if_allowed( pose::Pose & pose, Size const mutate_res, Real const mutation_frequency /* = 0.5 */ ){

		using namespace core::pose::full_model_info;
		using namespace protocols::rna;

		// first need to slice up native_pose to match residues in actual pose.
		// define atoms over which to compute RMSD, using rmsd_res.
		FullModelInfo const & full_model_info = const_full_model_info( pose );
		utility::vector1< Size > const sub_to_full = get_res_list_from_full_model_info( pose );
		std::string const full_sequence = full_model_info.full_sequence();

		if ( RG.uniform() < mutation_frequency) {

			char nt = full_sequence[ sub_to_full[ mutate_res ] - 1 ];
			choose_random_if_unspecified_nucleotide( nt );

			char const nt_orig = pose.sequence()[ mutate_res - 1 ];
			if ( nt != nt_orig ) {
				protocols::rna::mutate_position( pose, mutate_res, nt );
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

		core::Real	const principal_torsion = numeric::principal_angle_degrees( torsion_value );

		Size const principal_torsion_SIZE = Size( std::abs( principal_torsion + 0.00001 ) ); //0.00001 is to prevent random ambiguity if the torsion decimal value is exactly .0000 Oct 12, 2010


		if ( principal_torsion > 0 ){
			torsion_string = "p" + lead_zero_string_of( principal_torsion_SIZE, 3 );
		} else{
			torsion_string = "n" + lead_zero_string_of( principal_torsion_SIZE, 3 );
		}

		return torsion_string;
	}

	////////////////////////////////////////////////////////////////////////
	std::string //silly function used for appending the rotamer value to the tag
	create_rotamer_string( core::pose::Pose const & pose, Size const moving_res, bool const is_prepend ) {

		std::string rotamer_tag = "";

		conformation::Residue const & five_prime_rsd = ( is_prepend ) ? pose.residue( moving_res ): pose.residue( moving_res - 1 );
		conformation::Residue const & three_prime_rsd = ( is_prepend ) ?  pose.residue( moving_res + 1 ) : pose.residue( moving_res );


		rotamer_tag.append( "_E" + create_torsion_value_string( five_prime_rsd.mainchain_torsion( 5  ) ) );
		rotamer_tag.append( "_Z" + create_torsion_value_string( five_prime_rsd.mainchain_torsion( 6  ) ) );
		rotamer_tag.append( "_A" + create_torsion_value_string( three_prime_rsd.mainchain_torsion( 1 ) ) );
		rotamer_tag.append( "_B" + create_torsion_value_string( three_prime_rsd.mainchain_torsion( 2 ) ) );
		rotamer_tag.append( "_G" + create_torsion_value_string( three_prime_rsd.mainchain_torsion( 3 ) ) );


		if ( is_prepend ){
			rotamer_tag.append( "_D" + create_torsion_value_string( five_prime_rsd.mainchain_torsion( 4 ) ) );
			rotamer_tag.append( "_C" + create_torsion_value_string( five_prime_rsd.chi(  1 ) ) );

		} else{
			rotamer_tag.append( "_D" + create_torsion_value_string( three_prime_rsd.mainchain_torsion( 4 ) ) );
			rotamer_tag.append( "_C" + create_torsion_value_string( three_prime_rsd.chi( 1 ) ) );
		}

		return rotamer_tag;

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
		//	  ICOOR_INTERNAL  UPPER -175.907669   60.206192    1.607146   O3'   C3'   C4'  , Upper is P1
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

		FuncOP const distance_func( new HarmonicFunc( O3_P_distance, distance_stddev ) );
		FuncOP const O3_angle_func( new HarmonicFunc( radians( O3_angle ), radians( angle_stddev_degrees_P ) ) ); //did I made a mistake with the variable naming or is that a actual error? May 25, 2011.
		FuncOP const  P_angle_func( new HarmonicFunc( radians(  P_angle ), radians( angle_stddev_degrees_O3 ) ) ); //did I made a mistake with the variable naming or is that a actual error? May 25, 2011.

		Residue const & rsd1( pose.residue( five_prime_res ) );
		Residue const & rsd2( pose.residue( three_prime_res ) );

		AtomID const C3_id( rsd1.atom_index( "C3'" ), five_prime_res );
		AtomID const O3_id( rsd1.atom_index( "O3'" ), five_prime_res );
		AtomID const  P_id( rsd2.atom_index( "P"   ), three_prime_res );
		AtomID const O5_id( rsd2.atom_index( "O5'" ), three_prime_res );

		// distance from O3' to P
		cst_set->add_constraint( new AtomPairConstraint( O3_id, P_id, distance_func ) );

		// angle at O3'
		cst_set->add_constraint( new AngleConstraint( C3_id, O3_id, P_id, O3_angle_func ) );

		// angle at P
		cst_set->add_constraint( new AngleConstraint( O3_id, P_id, O5_id,  P_angle_func ) );

		pose.constraint_set( cst_set );
	}

}
}
}
