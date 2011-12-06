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

/// @file
/// @brief


// libRosetta headers
#include <core/scoring/rms_util.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueMatcher.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueSelector.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/VariantType.hh>

#include <core/chemical/ChemicalManager.hh>

//#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rna/RNA_Util.hh>
#include <core/scoring/rna/RNA_CentroidInfo.hh>
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <core/scoring/rna/RNA_TorsionPotential.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/sasa.hh>

#include <core/sequence/util.hh>

//Mmmm.. constraints.
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/constraints/HarmonicFunc.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/io/silent/RNA_SilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>

#include <core/pack/pack_rotamers.hh>
#include <core/pack/rotamer_trials.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

//////////////////////////////////////////////////
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

///////////////////////////////////////////////////

#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/util.hh>

#include <basic/options/option_macros.hh>
#include <protocols/idealize/idealize.hh>

#include <protocols/viewer/viewers.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

#include <basic/basic.hh>

#include <basic/database/open.hh>

#include <devel/init.hh>

#include <core/io/pdb/pose_io.hh>

#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>

#include <numeric/xyzVector.hh>
#include <numeric/conversions.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

//RNA stuff.
#include <protocols/rna/RNA_FragmentsClasses.hh>
#include <protocols/rna/RNA_DeNovoProtocol.hh>
#include <protocols/rna/RNA_Minimizer.hh>
#include <protocols/rna/RNA_LoopCloser.hh>
#include <protocols/rna/RNA_StructureParameters.hh>
#include <protocols/rna/RNA_ProtocolUtil.hh>

//For RNA jumps.
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <core/scoring/rna/RNA_FilteredBaseBaseInfo.hh>
#include <core/scoring/rna/RNA_BaseDoubletClasses.hh>

#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/EnergyMap.hh> //for EnergyMap
#include <core/scoring/EnergyMap.fwd.hh> //for EnergyMap


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



//silly using/typedef
#include <basic/Tracer.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/util.hh>


//#include <basic/tracer.hh>
using basic::T;
using basic::Error;
using basic::Warning;

using namespace core;
using namespace protocols;
using namespace basic::options::OptionKeys;

using utility::vector1;

using io::pdb::dump_pdb;

typedef  numeric::xyzMatrix< Real > Matrix;

//Definition of new OptionKeys
// these will be available in the top-level OptionKey namespace:
// i.e., OPT_KEY( Type, key ) -->  OptionKey::key
// to have them in a namespace use OPT_1GRP_KEY( Type, grp, key ) --> OptionKey::grp::key

OPT_KEY( Boolean, RNA_torsion_bug_function )
OPT_KEY( String, input_pose )
OPT_KEY( Boolean, old_file_format )
OPT_KEY( Boolean, exclude_bulge_residue_rmsd )
OPT_KEY( Boolean, vary_geometry )
OPT_KEY( String, pucker)
OPT_KEY( String, base)
OPT_KEY( String, outfile_name)
OPT_KEY( Boolean, native_sequence_mode)
OPT_KEY( Boolean, final_output)
OPT_KEY( Boolean, convert_star_to_dash)
OPT_KEY( Boolean, split_pdb_file)
OPT_KEY( Boolean, cut_bulge)
OPT_KEY( Boolean, print_hbonds )
OPT_KEY( Boolean, more_rotamers )
OPT_KEY( Boolean, quick_test )
OPT_KEY( Boolean, main_function ) //Create this one
OPT_KEY( Boolean, implement_rotamers ) //Create this one
OPT_KEY( Boolean, Full_minimize_function ) //Create this one
//OPT_KEY( Boolean, Minimize_chi_angle_function) //Create this one
OPT_KEY( Boolean, Chi_angle_function) //Create this one
//OPT_KEY( Boolean, Rescore_pose_list_function) //Create this one
OPT_KEY( Boolean, Torsion_analysis_function) //Create this one
OPT_KEY( Boolean, central_computer )
OPT_KEY( Boolean, read_job_file)
OPT_KEY( String, starting_pose)
OPT_KEY( Boolean, twoside_bound_original_poses)
OPT_KEY( String, HO2star_sampling)
//OPT_KEY( Boolean, o2star_trials_during_sampling )
//OPT_KEY( Boolean, Virtual_H2star)
OPT_KEY( Boolean, print_score_only )
OPT_KEY( Boolean, build_pose)
OPT_KEY( Boolean, full_base_pose_minimize)
OPT_KEY( Boolean, full_rebuild_pose_minimize)
OPT_KEY( Boolean, build_loop)
OPT_KEY( Boolean, test_mode)
OPT_KEY( Boolean, Finer_sampling)
OPT_KEY( Boolean, all_atom_cluster)
OPT_KEY( Boolean, minimize_and_score_sugar)
OPT_KEY( Integer, build_one_res_mode)
OPT_KEY( Boolean, Verbose)
OPT_KEY( String, score_function)
OPT_KEY( String, rebuild_sequence)
OPT_KEY( String, rotamer_set)
OPT_KEY( String, rebuild_algorithm)
OPT_KEY( String, rebuilt_pose_name )
OPT_KEY( String, apply_bound_function )
OPT_KEY( String, select_base_pose)
OPT_KEY( String, native_screen)
OPT_KEY( Real, torsion_bin)
OPT_KEY( Real, fine_torsion_bin)
OPT_KEY( Real, torsion_cutoff)
OPT_KEY( Real, bound_cutoff_value)
OPT_KEY( Boolean, apply_clustering )
OPT_KEY( Boolean, central_clustering )
OPT_KEY( Integer, pose_kept)
//OPT_KEY( Integer, first_res)
//OPT_KEY( Integer, last_res)
OPT_KEY( Integer, num_branch_kept)
OPT_KEY( Integer, hack_first_rebuild_num)
OPT_KEY( Integer, num_residue_rebuilt)
OPT_KEY( Integer, total_nodes_input)
OPT_KEY( Integer, node_number_input)
OPT_KEY( Integer, score_choice_option)
//OPT_KEY( Boolean, prepend_residue )
OPT_KEY( Real, fa_stack_weight )
OPT_KEY( Real, temperature )
OPT_KEY( Real, jump_change_frequency )
OPT_KEY( Real, suite_rmsd_cutoff )
OPT_KEY( Real, cluster_rmsd )
OPT_KEY( Real, fa_rep_cutoff)
OPT_KEY( String, bonus_stacking_score)
OPT_KEY( Integer, cycles )
OPT_KEY( String,  vall_torsions )
OPT_KEY( String,  assemble_file )
OPT_KEY( String,  basepair_file )
OPT_KEY( String,  jump_library_file )
OPT_KEY( String,  params_file )
OPT_KEY( String,  cst_file )
OPT_KEY( Real, atom_pair_constraint_weight )
OPT_KEY( Real, coordinate_constraint_weight )
OPT_KEY( Boolean, color_by_geom_sol )
OPT_KEY( Boolean, graphics )



////////////////////////////////////////////////////////////////////////////////////////////////
Real
convert_string_to_real(std::string& string){
	Real real_of_string;
	std::stringstream ss (std::stringstream::in | std::stringstream::out);
	ss << string;
	ss >> real_of_string;
	return real_of_string;

	std::cout << "The string  " <<  string << " have the corresponding real value " << real_of_string << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////
Size
convert_string_to_int(std::string& string){
	Size int_of_string; //misnomer
	std::stringstream ss (std::stringstream::in | std::stringstream::out);
	ss << string;
	ss >> int_of_string;


//	std::cout << "The string  " <<  string << " have the corresponding int value " << int_of_string << std::endl;

	return int_of_string;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


struct Residue_info_struct{

	std::string name;
	Size seq_num; //Full_pose_seq_number

};

//GLOBAL VARIABLE!!!
//std::map< Size_Pair, Size > Full_pose_to_partial_pose_seq_num_map;


namespace RunTimeParameters
{
	std::map< std::pair<Size, Size>, Size > Full_pose_to_partial_pose_seq_num_map;
	utility::vector1 <Residue_info_struct> rebuild_residue_list;
	utility::vector1 <bool> append_prepend_choice_list;
	utility::vector1 <Residue_info_struct> original_pose_residue_list;
}


Size
get_max_res(utility::vector1 <Residue_info_struct>& residue_list){

	Size current_max_res=0;
	for(Size i=1; i<=residue_list.size(); i++){
		 if(residue_list[i].seq_num>current_max_res) current_max_res=residue_list[i].seq_num;
	}
	return current_max_res;

}

std::string
Get_one_letter_name(std::string three_letter_name){
	if(three_letter_name=="RAD") return "A";
	if(three_letter_name=="RCY") return "C";
	if(three_letter_name=="URA") return "U";
	if(three_letter_name=="RGU") return "G";
	std::cout << "In get_one_letter_name_function, an invalid three_letter_name was passed into the function" << std::endl;
	exit (1);
}

std::string
Get_three_letter_name(std::string one_letter_name){
	if(one_letter_name=="A") return "RAD";
	if(one_letter_name=="C") return "RCY";
	if(one_letter_name=="U") return "URA";
	if(one_letter_name=="G") return "RGU";
	std::cout << "In get_three_letter_name_function, an invalid one_letter_name was passed into the function" << std::endl;
	exit (1);
}

bool
Contain_residue_at_seq_num(Size seq_num, utility::vector1 <Residue_info_struct>& residue_list){

		for(Size j=1; j<=residue_list.size(); j++){
			if(seq_num==residue_list[j].seq_num) {
				std::cout << Get_one_letter_name(residue_list[j].name);
				std::cout << lead_zero_string_of(residue_list[j].seq_num, 2);
				std::cout << " ";
				return true;
			}
		}
		return false;
}

void
Output_residue_list(utility::vector1 <Residue_info_struct>& residue_list){

	Size max_res=get_max_res(residue_list);
//	std::cout<< "max_res= " << max_res << std::endl;

	for(Size i=1; i<= max_res; i++){

		if(Contain_residue_at_seq_num(i, residue_list)){
		} else {
				std::cout << "    ";
		}

	}
	std::cout << std::endl;

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void Tokenize_dash_delimiter(std::string& str, utility::vector1<std::string>& tokens)
{
	  using namespace std;
	  std::string delimiters= "-";

    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
get_missing_residue_list(utility::vector1 <Residue_info_struct>& missing_residue_list, Size current_rebuild_num){

	using namespace RunTimeParameters;
	missing_residue_list.clear();

		for(Size j=rebuild_residue_list.size(); j> current_rebuild_num; j--){
				missing_residue_list.push_back(rebuild_residue_list[j]);
		}
}

//The union of partial_pose_residue_list and missing_residue_list is the original_pose_residue_list

void
get_partial_pose_residue_list(utility::vector1 <Residue_info_struct>& partial_pose_residue_list, Size current_rebuild_num){

	partial_pose_residue_list=RunTimeParameters::original_pose_residue_list;

	utility::vector1 <Residue_info_struct> missing_residue_list;
	get_missing_residue_list(missing_residue_list, current_rebuild_num);

	for(Size i=partial_pose_residue_list.size(); i>=1 ; i--){
		for(Size j= missing_residue_list.size(); j>0; j--){
				if(partial_pose_residue_list[i].seq_num==missing_residue_list[j].seq_num){
					partial_pose_residue_list.erase(partial_pose_residue_list.begin()+(i-1));
				}
		}
	}
}



//////////The three functions below and Initialize_append_prepend_choice_list() create the four objects in RunTimeParameter namespace and should be called once//////////////////////////////

void
Create_rebuild_residue_list(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace RunTimeParameters;

	 std::string rebuild_sequence_string = option[ rebuild_sequence];

	 utility::vector1<std::string> nucleotides_token;

 	 Tokenize_dash_delimiter(rebuild_sequence_string, nucleotides_token);

		for(Size i=1; i<=nucleotides_token.size(); i++){
			std::cout<< nucleotides_token[i] << " ";

			Residue_info_struct res_info;
			res_info.name=Get_three_letter_name(nucleotides_token[i].substr(0,1));
			std::string seq_num_string=nucleotides_token[i].substr(1,2);
			res_info.seq_num= convert_string_to_int(seq_num_string);
			RunTimeParameters::rebuild_residue_list.push_back( res_info);

		}
		std::cout << std::endl;

	std::cout << "Create rebuild_residue_list: "; Output_residue_list(RunTimeParameters::rebuild_residue_list);
}

//This function should be only call once in the whole code

void
Create_original_pose_residue_list(pose::Pose& pose){

	Residue_info_struct res_info;

	for(Size i=1; i<=pose.total_residue(); i++){

		chemical::AA res_aa =  pose.residue( i ).aa();
		res_info.name = name_from_aa( res_aa );
		res_info.seq_num=i;
		RunTimeParameters::original_pose_residue_list.push_back(res_info );
	}

	std::cout << "Create original_pose_residue_list: "; Output_residue_list(RunTimeParameters::original_pose_residue_list);

}


void
Create_rebuild_num_to_seq_num_map(pose::Pose& original_pose){

	using namespace RunTimeParameters;

	for(Size current_rebuild_num=0; current_rebuild_num<=RunTimeParameters::rebuild_residue_list.size(); current_rebuild_num++){
			utility::vector1 <Residue_info_struct> partial_pose_residue_list;
			get_partial_pose_residue_list(partial_pose_residue_list, current_rebuild_num);

			utility::vector1 <Residue_info_struct> missing_residue_list;
			get_missing_residue_list(missing_residue_list, current_rebuild_num);

			std::cout<< "Current_rebuild_num: " << current_rebuild_num << std::endl;
			Output_residue_list(partial_pose_residue_list);
			Output_residue_list(missing_residue_list);
	}


  ////////////////////Initialize the Full_pose_to_partial_pose_seq_num_map///////////////////////////////////////////////////////

	for(Size current_rebuild_num=0; current_rebuild_num<=rebuild_residue_list.size(); current_rebuild_num++){
		for(Size full_pose_seq_num=1; full_pose_seq_num<=original_pose.total_residue(); full_pose_seq_num++){

			utility::vector1 <Residue_info_struct> missing_residue_list;
			get_missing_residue_list(missing_residue_list, current_rebuild_num);

			Size num_missing_residues=0;

			for(Size i=1; i<=missing_residue_list.size(); i++){
				if(missing_residue_list[i].seq_num<full_pose_seq_num){
					num_missing_residues++;
				}
			}

			Size partial_pose_seq_num=full_pose_seq_num-num_missing_residues;

			//If residue corresponding full_pose_seq_num doesn't exist in rebuild pose, overwrite and set partial_pose_seq_num to zero.
			for(Size i=1; i<=missing_residue_list.size(); i++){
				if(missing_residue_list[i].seq_num==full_pose_seq_num){
					partial_pose_seq_num=0;
					break;
				}
			}

			 Full_pose_to_partial_pose_seq_num_map.insert ( std::pair<std::pair<Size, Size>, Size>(std::make_pair(full_pose_seq_num, current_rebuild_num), partial_pose_seq_num ));
		}
	}

  /////Debug:Output map values//////////////////////////////////////////////////////
	Size spacing=4;

	std::cout << "    ";
	for(Size full_pose_seq_num=1; full_pose_seq_num<=original_pose.total_residue(); full_pose_seq_num++){

		std::cout << std::setw(spacing) << full_pose_seq_num;
	}
	std::cout << std::endl;

	for(Size current_rebuild_num=0; current_rebuild_num<=rebuild_residue_list.size(); current_rebuild_num++){
		std::cout << std::setw(spacing) << current_rebuild_num;

		for(Size full_pose_seq_num=1; full_pose_seq_num<=original_pose.total_residue(); full_pose_seq_num++){

			Size partial_pose_seq_num=Full_pose_to_partial_pose_seq_num_map[std::make_pair(full_pose_seq_num, current_rebuild_num)];

			std::cout << std::setw(spacing) << partial_pose_seq_num;

		}
		std::cout << std::endl;
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool
residue_list_sort_citeria(Residue_info_struct residue_info_1, Residue_info_struct residue_info_2){


	//Sort by seq_number, lowest sequence number at the top of the vector list.
	return (residue_info_1.seq_num < residue_info_2.seq_num);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
sort_residue_list(utility::vector1<Residue_info_struct>& residue_list) {

		//Need to check if this work with vector1, if not switch to std::vector
		sort(residue_list.begin(), residue_list.end(), residue_list_sort_citeria);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//More general version of create_rebuild_residue_groups, allow for any residue_list.....will make create_rebuild_residue_groups call this function after I test that this function work
void
create_residue_group_list(utility::vector1 < utility::vector1 <Residue_info_struct> >& residue_group_list, utility::vector1 <Residue_info_struct> residue_list){

	utility::vector1 <Residue_info_struct> Sorted_residue_list = residue_list;
	sort_residue_list(Sorted_residue_list);

	Size j=1;
	while(j<=Sorted_residue_list.size()){

			Size first_element=j;

			//Test if Sorted_residue_list[j] contain an adjacent residue at the three_prime_end
			while((j<Sorted_residue_list.size()) && ((Sorted_residue_list[j].seq_num+1)==Sorted_residue_list[j+1].seq_num)){
				j++;
			}
			Size last_element=j;

			utility::vector1 <Residue_info_struct> residue_group;

			for(Size element=first_element; element<=last_element; element++){
				residue_group.push_back(Sorted_residue_list[element]);
			}
			residue_group_list.push_back(residue_group);

			j++;
	}

}


void
create_rebuild_residue_groups(utility::vector1 < utility::vector1 <Residue_info_struct> >& rebuild_residue_group_list){

	utility::vector1 <Residue_info_struct> Sorted_rebuild_residue_list = RunTimeParameters::rebuild_residue_list;
	sort_residue_list(Sorted_rebuild_residue_list);

	Size j=1;
	while(j<=Sorted_rebuild_residue_list.size()){

			Size first_element=j;

	//		Size five_prime=Sorted_rebuild_residue_list[j].seq_num-1;

			//Test if Sorted_rebuild_residue_list[j] contain an adjacent rebuild residue at the three_prime_end
			while((j<Sorted_rebuild_residue_list.size()) && ((Sorted_rebuild_residue_list[j].seq_num+1)==Sorted_rebuild_residue_list[j+1].seq_num)){
				j++;
			}
			Size last_element=j;

//			Size three_prime=Sorted_rebuild_residue_list[j].seq_num+1;

//			jump_points_list.push_back(std::make_pair(five_prime, three_prime));

			utility::vector1 <Residue_info_struct> rebuild_residue_group;

			for(Size element=first_element; element<=last_element; element++){
				rebuild_residue_group.push_back(Sorted_rebuild_residue_list[element]);
			}
			rebuild_residue_group_list.push_back(rebuild_residue_group);

			j++;
	}

//	std::cout << "In create_rebuild_residue_groups function:" << std::endl;
//	for(Size i=1; i<= rebuild_residue_group_list.size(); i++){
//		Output_residue_list(rebuild_residue_group_list[i]);
//	}


}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool
Is_five_prime_residue_build(utility::vector1 <Residue_info_struct>& missing_residue_list, Size & rebuild_res_seq_num){

	for(Size i=1; i<= missing_residue_list.size(); i++){
		if(missing_residue_list[i].seq_num==(rebuild_res_seq_num-1)) return false; //The adjacent five_prime residue to the rebuild_residue is not build yet...cannot append!
	}
	return true; //Note, this assume that there is a five_prime resiude to the current_rebuild_residue in the full pose, that is, the residue is not the first residue in the sequence.
}

bool
Is_three_prime_residue_build(utility::vector1 <Residue_info_struct>& missing_residue_list, Size & rebuild_res_seq_num){

	for(Size i=1; i<= missing_residue_list.size(); i++){
		if(missing_residue_list[i].seq_num==(rebuild_res_seq_num+1)) return false; //The adjacent three_prime residue to the rebuild_residue is not build yet...cannot prepend!
	}
	return true; //Note, this assume that there is a three_prime resiude to the current_rebuild_residue in the full pose, that is, the residue is not the last residue in the sequence.
}

bool
Initialize_Is_prepend(Size& current_rebuild_num){

	//If the seq gap is 1 (such as during last residue rebuild), then could either prepend or append. Will choose prepend for this situation.

	if(current_rebuild_num>RunTimeParameters::rebuild_residue_list.size()){
		std::cout << "Error: The value of current_rebuild_num passed into the function is greater than the number of residue rebuilt!" << std::endl;
		exit (1);
	}

	utility::vector1 <Residue_info_struct> missing_residue_list;
	get_missing_residue_list(missing_residue_list, current_rebuild_num);

	Size rebuild_res_seq_num=RunTimeParameters::rebuild_residue_list[current_rebuild_num].seq_num;

	if(Is_three_prime_residue_build(missing_residue_list, rebuild_res_seq_num)) return true; // Choose prepend

	if(Is_five_prime_residue_build(missing_residue_list, rebuild_res_seq_num)) return false; // Choose append

	std::cout << "Error: In Initialize_Is_prepend function, cannot rebuild residue: ";
	std::cout << Get_one_letter_name(RunTimeParameters::rebuild_residue_list[current_rebuild_num].name);
	std::cout << RunTimeParameters::rebuild_residue_list[current_rebuild_num].seq_num;
	std::cout << " (residue_rebuild_num= " << current_rebuild_num << ") by either prepending and appending!" << std::endl;
	exit (1);

}


/////////////////////////////////////////////////////////////////////////////////////////
void
Initialize_append_prepend_choice_list(){

	Size num_residue_reb = RunTimeParameters::rebuild_residue_list.size();

	utility::vector1 <bool> temp_append_prepend_choice_list;

	for(Size current_rebuild_num=1; current_rebuild_num<=num_residue_reb; current_rebuild_num++){
		temp_append_prepend_choice_list.push_back(Initialize_Is_prepend(current_rebuild_num));
	}

	std::cout << "prepend/append choice:" << std::endl;
	for(Size current_rebuild_num=1; current_rebuild_num<=num_residue_reb; current_rebuild_num++){
		std::cout << "current_rebuild_num: " << current_rebuild_num;
		std::cout << "  choice:  ";
		if(temp_append_prepend_choice_list[current_rebuild_num]==true) std::cout << "prepend" << std::endl;
		if(temp_append_prepend_choice_list[current_rebuild_num]==false) std::cout << "append" << std::endl;

	}

	RunTimeParameters::append_prepend_choice_list=temp_append_prepend_choice_list;


}


void
Initialize_RunTimeParameters(){

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace core::kinematics;
	using namespace core::optimization;
	using namespace core::io::silent;
	using namespace core::id;
	using namespace protocols::rna;

	ResidueTypeSetCAP rsd_set; //This line is repeated
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "rna" ); //This line is repeated

	//The sole purpose of this original_pose is to initialize the RunTimeParameters
	pose::Pose original_pose;

	std::string infile = option[ in ::file::s ][1];

	core::import_pose::pose_from_pdb( original_pose, *rsd_set, infile );

	std::cout << "Initializing RunTimeParameters" << std::endl;
	std::cout << "Native pose : " << infile  << std::endl;

	Create_rebuild_residue_list(); //This function create the vector rebuild_residue_list which is stored in the namespace RunTimeParameters
	Initialize_append_prepend_choice_list(); //This function create the constant vector append_prepend_choice_list which is stored in the namespace RunTimeParameters

	Create_original_pose_residue_list(original_pose); //This function create the vector original_pose_residue_list which is stored in the namespace RunTimeParameters
	Create_rebuild_num_to_seq_num_map(original_pose); //This function create the map full_pose_to_partial_pose_seq_num_map which is stored in the namespace RunTimeParameters


	//Check that base type of rebuild residue match that of the original pose
	if( option[ native_sequence_mode]) {
		for(Size i=1; i<=RunTimeParameters::rebuild_residue_list.size(); i++){
				for(Size j=1; j<RunTimeParameters::original_pose_residue_list.size(); j++){

					if(RunTimeParameters::rebuild_residue_list[i].seq_num==RunTimeParameters::original_pose_residue_list[j].seq_num){
						if(RunTimeParameters::rebuild_residue_list[i].name!=RunTimeParameters::original_pose_residue_list[j].name){
							std::cout << "In Initialize_RunTimeParameters function:" << std::endl;
							std::cout << "	Non-native base is not allowed in native_sequence_mode" << std::endl;
							std::cout << "	For seq_num =" << RunTimeParameters::rebuild_residue_list[i].seq_num << std::endl;
							std::cout << "		Original base is :" << RunTimeParameters::original_pose_residue_list[j].name << std::endl;
							std::cout << "  	Input rebuild base is :" << RunTimeParameters::rebuild_residue_list[i].name << std::endl;
							exit (1);
						}
					}

				}
		}
	}

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void
Output_fold_tree_info(pose::Pose& pose, std::string pose_name){

		std::cout << "fold tree of " << pose_name << ": " << std::endl;
		std::cout << "	number of chainbreak= " << pose.fold_tree().num_cutpoint() << std::endl;
		for(Size i=1; i<=pose.fold_tree().num_cutpoint(); i++){
			std::cout << "	Chain break #" << i << " :";;
			std::cout << "  cutpoint= " << pose.fold_tree().cutpoint(i);
			std::cout << "  5 prime jump point= " << pose.fold_tree().jump_point( 1, i );
			std::cout << "  3 prime jump point= " << pose.fold_tree().jump_point( 2, i ) << std::endl;
		}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
Is_prepend(Size const & current_rebuild_num){

		if(current_rebuild_num < 1 || current_rebuild_num > RunTimeParameters::rebuild_residue_list.size()){
			std::cout << "Error: In Is_prepend function:" << std::endl;
			std::cout << "The value of current_rebuild_num passed into the function is " <<current_rebuild_num;
			std::cout << ". This value is out of range! " << std::endl;
   		exit (1);
		}

		return RunTimeParameters::append_prepend_choice_list[current_rebuild_num];
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Size
get_full_pose_reb_res(Size const & current_rebuild_num){
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

//		Size full_pose_reb_res_old = get_full_pose_reb_res_old(current_rebuild_num);

		//Check that current_rebuild_num is in range;
		if(current_rebuild_num < 1 || current_rebuild_num > RunTimeParameters::rebuild_residue_list.size()){
			std::cout << "Error: In get_full_pose_reb_res function:" << std::endl;
			std::cout << "The value of current_rebuild_num passed into the function is " <<current_rebuild_num;
			std::cout << ". This value is out of range! " << std::endl;
   		exit (1);
		}

		Size full_pose_reb_res = RunTimeParameters::rebuild_residue_list[current_rebuild_num].seq_num;

//		std::cout << "full_pose_reb_res: " << full_pose_reb_res;
//		std::cout << " full_pose_reb_res_old: " << full_pose_reb_res_old << std::endl;

//		if(full_pose_reb_res!=full_pose_reb_res_old){
//			std::cout << "Error: full_pose_reb_res does not equal full_pose_reb_res_old" << std::endl;
//			exit (1);
//		}

		return full_pose_reb_res;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Size
Convert_to_partial_pose_seq_num(Size const & full_pose_seq_num, Size  const & current_rebuild_num){

//	std::pair<Size, Size> full_pose_seq_num_and_current_rebuild_num (full_pose_seq_num, current_rebuild_num);

	if(current_rebuild_num < 0 || current_rebuild_num > RunTimeParameters::rebuild_residue_list.size()){
			std::cout << "Error: In Convert_to_partial_pose_seq_num function:" << std::endl;
			std::cout << "The value of current_rebuild_num passed into the function is " << current_rebuild_num;
			std::cout << ". This value is out of range! " << std::endl;
   		exit (1);
	}


	Size partial_pose_seq_num=RunTimeParameters::Full_pose_to_partial_pose_seq_num_map[std::make_pair(full_pose_seq_num, current_rebuild_num)];

	if(partial_pose_seq_num==0){
		std::cout << "Error: full_pose_seq_num does not exist in partial_pose:" << std::endl;
		std::cout << "  Value passed into the function:" << std::endl;
		std::cout << "    full_pose_seq_num= " << full_pose_seq_num;
		std::cout << "	current_rebuild_num= " << current_rebuild_num << std::endl;
		exit (1);
	}


//	Size partial_pose_seq_num_old=Convert_to_partial_pose_seq_num_old(full_pose_seq_num, current_rebuild_num);

//	std::cout<< "partial_pose_seq_num: " << partial_pose_seq_num;
//	std::cout<< " partial_pose_seq_num_old: " << partial_pose_seq_num_old << std::endl;

//	if(partial_pose_seq_num!=partial_pose_seq_num_old){
//		std::cout << "Error: partial_pose_seq_num does not equal partial_pose_seq_num_old" << std::endl;
//		exit (1);
//	}

	return partial_pose_seq_num;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Size
get_rebuild_pose_res_num(Size const & current_rebuild_num){

	Size full_pose_reb_res=get_full_pose_reb_res(current_rebuild_num);
	Size rebuild_pose_reb_res=Convert_to_partial_pose_seq_num(full_pose_reb_res, current_rebuild_num);
	return rebuild_pose_reb_res;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
get_jump_points(utility::vector1 <std::pair<Size, Size> >& jump_points_list, utility::vector1 < utility::vector1 <Residue_info_struct> >& rebuild_residue_group_list){

	jump_points_list.clear();

	for(Size i=1; i<=rebuild_residue_group_list.size(); i++){

		Size first_element=1;
		Size last_element=rebuild_residue_group_list[i].size();

		Size five_prime=rebuild_residue_group_list[i][first_element].seq_num-1;
		Size three_prime=rebuild_residue_group_list[i][last_element].seq_num+1;

		jump_points_list.push_back(std::make_pair(five_prime, three_prime));
	}

}

//Get cut_points of full_pose
//If the seq gap is 1 (just as during last residue prepend), then could either prepend or append. Will choose prepend for this situation.
//Then to be consistent need to choose five_prime cut_point
void
get_cut_points(utility::vector1 <Size>& cut_points_list, utility::vector1 <std::pair<Size, Size> >& jump_points_list){
	cut_points_list.clear();


	//find seq_num of residue that will be

	for(Size i=1; i<= jump_points_list.size(); i++){

		for(Size j=RunTimeParameters::rebuild_residue_list.size(); j>=1; j--){
			if((RunTimeParameters::rebuild_residue_list[j].seq_num>jump_points_list[i].first) && (RunTimeParameters::rebuild_residue_list[j].seq_num<jump_points_list[i].second)){
				//This work since rebuild_residue_list is ordered according the order in which the residue is rebuild
				cut_points_list.push_back(RunTimeParameters::rebuild_residue_list[j].seq_num-1); //To be consistent, need to choose five_prime cut_point
				break;
			}
		}
//		std::cout << "Error: no cut_point found for chain_break #" << i << std::endl;
//		exit (1);
	}
}

//Stole from RNA_Until.cc
bool
is_rna_chainbreak_parin( core::pose::Pose & pose, Size const & i ) {

	static Real const CHAINBREAK_CUTOFF2 ( 2.5 * 2.5 ); //03'--P bond length is roughly 1.6 Angstrom

	if ( i >= pose.total_residue() ) return true;
	if ( i < 1 ) return true;

	conformation::Residue const & current_rsd( pose.residue( i   ) ) ;
	conformation::Residue const &    next_rsd( pose.residue( i+1 ) ) ;

	//A little inefficient, since atom indices for these backbone
	// atoms should be the same for all RNA residue types. I think.
	Size atom_O3star = current_rsd.atom_index( " O3*" );
	Size atom_P      =    next_rsd.atom_index( " P  " );
	Real const dist2 =
		( current_rsd.atom( atom_O3star ).xyz() - next_rsd.atom( atom_P ).xyz() ).length_squared();

	if ( dist2 > CHAINBREAK_CUTOFF2 ) {
		//std::cout << "Found chainbreak at residue "<< i << " .  O3*-P distance: " << sqrt( dist2 ) << std::endl;
		return true;
	}

	return false;

}


//Stole from RNA_Protocol_Util.cc
void
figure_out_reasonable_rna_fold_tree_parin( pose::Pose & pose )
{
	using namespace core::conformation;

	//Look for chainbreaks in PDB.
	Size const nres = pose.total_residue();
	kinematics::FoldTree f( nres );

	Size m( 0 );

	for (Size i=1; i < nres; ++i) {

		if ( is_rna_chainbreak_parin( pose, i ) ){

			f.new_jump( i, i+1, i );
			m++;

			Residue const & current_rsd( pose.residue( i   ) ) ;
			Residue const &    next_rsd( pose.residue( i+1 ) ) ;
			//			Size dummy( 0 ), jump_atom1( 0 ), jump_atom2( 0 );
			//rna_basepair_jump_atoms( current_rsd.aa(), jump_atom1, dummy, dummy );
			//rna_basepair_jump_atoms( next_rsd.aa(), jump_atom2, dummy, dummy );
			//f.set_jump_atoms( m, current_rsd.atom_name( jump_atom1 ), next_rsd.atom_name( jump_atom2 ) );

//			f.set_jump_atoms( m,core::scoring::rna::chi1_torsion_atom( current_rsd ), core::scoring::rna::chi1_torsion_atom( next_rsd )   );
			//Between last atom of five_prime residue and first atom of three residue. Note that this is completely opposite compare rebuild residue loop chainbreak
			f.set_jump_atoms( m, i , "O3*" , i+1 , "P" );
		}

	}
	pose.fold_tree( f );
	std::cout << "  "; Output_fold_tree_info(pose, 	"Natural_cutpoint_of_original_pose ");
}


void
Setup_jump_and_cut_point(pose::Pose& original_pose){

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	//Setup natural cut_point of original pose (exist when original pose is composed of multiple strands)
	figure_out_reasonable_rna_fold_tree_parin(original_pose);

	//Cut points for loop rebuild (below)
	utility::vector1 <std::pair<Size, Size> > jump_points_list;
	utility::vector1 <Size> cut_points_list;


	utility::vector1 < utility::vector1 <Residue_info_struct> > rebuild_residue_group_list;
	create_rebuild_residue_groups(rebuild_residue_group_list); //A rebuild_residue_group is a list of Residue_info_struct sorted by seq_num. There is one group for each chain break.

	get_jump_points(jump_points_list, rebuild_residue_group_list); // One jump_point_pair per chain break, ordered by seq_num
	get_cut_points(cut_points_list, jump_points_list); // One cut_point per chain break, ordered by seq_num

	kinematics::FoldTree f=original_pose.fold_tree();
	Size num_natural_cut_points=f.num_cutpoint();
	std::cout << "num_natural_cut_points= " << num_natural_cut_points << std::endl;

//	kinematics::FoldTree f(original_pose.total_residue());

	for(Size i=1; i<=jump_points_list.size(); i++){

		Size jump_num=f.new_jump( jump_points_list[i].first, jump_points_list[i].second, cut_points_list[i]);

//		std::cout << "check jump_number: " << jump_num << std::endl;

		f.set_jump_atoms( i+num_natural_cut_points, jump_points_list[i].first , "P", jump_points_list[i].second , "O3*" ); //Weird that will have to restate five_prime and three_prime jump_point...
		//This atom choice is necessary to allow the five_prime and three_prime residue itself to be movable.
		//!!Problem then is that there must be at least two residue between adjacent chain_break (last (three_prime_most) residue of one loop/chainbreak
		//and first (five_prime_most) residue of the next loop/chainbreak;!!

		//There is alternative version of set_jump_atoms which require as out AtomID..The AtomID should be specified as follow (haven't test this yet).
		//AtomID const P_id( rsd1.atom_index( "P*" ), jump_points_list[i].first);
		//AtomID const O3_id( rsd1.atom_index( "O3*" ), jump_points_list[i].second);

	}

	original_pose.fold_tree( f );
  kinematics::FoldTree fold_tree=original_pose.fold_tree();
	Output_fold_tree_info(original_pose, "original_pose");
}




///////////////////////////////////////////////////////////////////////////////
void
color_by_geom_sol_RNA_test()
{
/*
	// Read in pdb.
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::pose;
	using namespace core::id;
	using namespace core::scoring;
	using namespace core::kinematics;
	using namespace protocols::rna;

	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "rna" );
	utility::vector1 < std::string> pdb_files( option[ in::file::s ]() );

	ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( RNA_HIRES_WTS );
	hbonds::GeometricSolEnergy geometric_sol_energy_method( scorefxn->energy_method_options() );

	Pose pose;

	for (Size n = 1; n <= pdb_files.size(); n++) {
		std::string const pdb_file = pdb_files[n];
		core::import_pose::pose_from_pdb( pose, *rsd_set, pdb_file );
		protocols::rna::ensure_phosphate_nomenclature_matches_mini( pose );

		PDBInfoOP pdb_info(  new PDBInfo( pose, true ) );

		(*scorefxn)( pose );

		for ( Size i = 1; i <= pose.total_residue(); i++ ) {
			for (Size j = 1; j <= pose.residue( i ) .natoms(); j++ ) {
				pdb_info->temperature( i, j, geometric_sol_energy_method.eval_atom_energy( AtomID( j, i ), pose ) );
			}
		}

		pose.pdb_info( pdb_info );

		Size pos( pdb_file.find( ".pdb" ) );
		std::string outfile( pdb_file );
		outfile.replace( pos, 4, "_color_geomsol.pdb" );
		pose.dump_pdb( outfile );
	}
*/
}


///////////////////////////////////////////////////////////////////////////////
void
print_hbonds_test()
{
	// Read in pdb.
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::kinematics;
	using namespace core::optimization;
	using namespace core::io::silent;
	using namespace protocols::rna;

	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "rna" );
	utility::vector1 < std::string> pdb_files( option[ in::file::s ]() );

	ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( option[score::weights] );

	for (Size n = 1; n <= pdb_files.size(); n++) {
		std::string const pdb_file = pdb_files[n];

		pose::Pose pose;
		core::import_pose::pose_from_pdb( pose, *rsd_set, pdb_file );

		(*scorefxn)(pose);
		hbonds::HBondSetOP hbond_set( new hbonds::HBondSet() );
		hbond_set->use_hb_env_dep( false );

		hbonds::fill_hbond_set( pose, false /*calc deriv*/, *hbond_set );

		for (Size i = 1; i <= hbond_set->nhbonds(); i++ ) {
			hbonds::HBond const & hbond( hbond_set->hbond( i ) );

			Size const don_res_num = hbond.don_res();
			Size const don_hatm = hbond.don_hatm();

			Size const acc_res_num = hbond.acc_res();
			Size const acc_atm = hbond.acc_atm();


			std::cout << "HBOND: " << pose.residue( don_res_num ).name1() << don_res_num <<
				" " << pose.residue( don_res_num ).atom_name( don_hatm ) << " --- " <<
				pose.residue( acc_res_num).name1() << acc_res_num << " " << pose.residue( acc_res_num ).atom_name( acc_atm ) << " ==> " << hbond.energy()
								<< std::endl;

		}

	}
}

///////////////////////////////////////////////////////////////////////////////
void
create_dihedral_constraint( core::scoring::constraints::ConstraintSetOP & cst_set,
														core::pose::Pose & pose,
														core::id::TorsionID torsion_id )
{
	using namespace core::id;
	using namespace core::scoring;
	using namespace core::scoring::constraints;

	id::AtomID id1,id2,id3,id4;
	bool fail = pose.conformation().get_torsion_angle_atom_ids( torsion_id, id1, id2, id3, id4 );

	Real torsion_value( pose.torsion( torsion_id ) );

	HarmonicFuncOP harm_func  (new HarmonicFunc( numeric::conversions::radians( torsion_value ),
																							 numeric::conversions::radians( 20.0 ) ) );

	ConstraintOP dihedral1 =
		new DihedralConstraint( id1, id2, id3, id4,
																				 harm_func,
																				 rna_torsion );
	cst_set->add_constraint( dihedral1 );
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
get_correct_rotamers( utility::vector1< utility::vector1 <utility::vector1 <Real > > >& backbone_rotamers_groups, pose::Pose template_pose){

	using namespace core::scoring;
	using namespace core::scoring::rna;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	Real error= option[torsion_bin];

	std::cout << "check point 1" << std::endl;


	Size const nres( template_pose.total_residue() );


//	utility::vector1 <Real > rotamer_torsions(13 /*# torsion angles*/, 9999);

//	utility::vector1 <utility::vector1 <Real > > residue_rotamers(1 /*for now */, rotamer_torsions);


//	backbone_rotamers_groups.assign(nres, residue_rotamers);



	//backbone_rotamers_groups_north is used differently than in standard.
	//backbone_rotamers_groups_north[res_num][rotamers][torsions] instead of [group_rotamer][subgroup_rotamer][torsions];


	//conformation::Residue const &  current_res(current_pose.residue(last_res_append));

//	utility::vector1 < conformation::ResidueOP> residues;

//	for(Size res_num=1; res_num <= nres; res_num++){
//			residues.push_back(template_pose.residue(res_num));
//	}




	for(Size res_num=1; res_num< nres; res_num++){ //can't the nres residue since nres+1 doesn't exist

			utility::vector1 <utility::vector1 <Real > > residue_rotamers;

			conformation::Residue const & current_res=template_pose.residue(res_num);
			conformation::Residue const & next_res=template_pose.residue(res_num+1);

			utility::vector1 <Real > mean_rotamer_torsions(13 /*# torsion angles*/, 9999);


			for(Size n=1; n<=1; n++){ //delta1
				Size i=n+3;
				mean_rotamer_torsions[n]=current_res.mainchain_torsion(i);
			}

			for(Size n=2; n<=4; n++){ // chi_1, nu2_1, nu1_1 (recall that chi_1 2OH is determined seperately)
				Size i=n-1;
					mean_rotamer_torsions[n]=current_res.chi(i);
			}

		  for(Size n=5; n<=6; n++){ //epsilon1, zeta1
				Size i=n;
					mean_rotamer_torsions[n]=current_res.mainchain_torsion(i);
			}


			for(Size n=7; n<=10; n++){ //alpha2, beta2, gamma2, delta2
				Size i=n-6;
					mean_rotamer_torsions[n]=next_res.mainchain_torsion(i);
			}

			for(Size n=11; n<=13; n++){ // chi_2, nu2_2, nu1_2
				Size i=n-10;
					mean_rotamer_torsions[n]=next_res.chi(i);
			}

		Size spacing=8;

		std::cout << "correct torsions of residue, " << res_num <<  std::endl;
		for(Size n=1; n<=13; n++){
		std::cout << std::setw(spacing)<< mean_rotamer_torsions[n] << " ";
		}


		std::cout<< std::endl;

		utility::vector1 <Real > variation_rotamer_torsions=mean_rotamer_torsions;

		std::cout << "error " << error <<std::endl;

		//Why not -1 to 1??
		for(int z1=0; z1<=2; z1++){
			variation_rotamer_torsions[6]=mean_rotamer_torsions[6]+error*z1;
			for(int a2=0; a2<=2; a2++){
				variation_rotamer_torsions[7]=mean_rotamer_torsions[7]+error*a2;
				for(int b2=0; b2<=2; b2++){
					variation_rotamer_torsions[8]=mean_rotamer_torsions[8]+error*b2;
					for(int g2=0; g2<=2; g2++){
						variation_rotamer_torsions[9]=mean_rotamer_torsions[9]+error*g2;

						residue_rotamers.push_back(variation_rotamer_torsions);


					}
				}
			}
		}


		backbone_rotamers_groups.push_back(residue_rotamers);

	}

	std::cout << "End of get_correct_rotamers function " << std::endl;

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
Rebuild_is_missing_in_pose(Size seq_num, utility::vector1 <Residue_info_struct>& partial_pose_residue_list){

	for(Size i=1; i<= partial_pose_residue_list.size(); i++){
		if(partial_pose_residue_list[i].seq_num==seq_num) return false;
	}
	return true;
}


utility::vector1 <Real>
get_O2star_torsion_list(pose::Pose& template_pose, Size current_rebuild_num){



	utility::vector1 <Residue_info_struct> partial_pose_residue_list;
	get_partial_pose_residue_list(partial_pose_residue_list, current_rebuild_num);
	Output_residue_list(partial_pose_residue_list);

	Size total_residue_num=RunTimeParameters::original_pose_residue_list.size();

	utility::vector1 <Real > O2_star_torsion_list(total_residue_num, 9999);


	for(Size n=1; n<= total_residue_num; n++){

		if(Rebuild_is_missing_in_pose(n, partial_pose_residue_list)==false){
			Size res_num=Convert_to_partial_pose_seq_num(n, current_rebuild_num);
			conformation::Residue const & current_residue=template_pose.residue(res_num);
			O2_star_torsion_list[n]=current_residue.chi(4); //chi_O2star torsion
		}

	}

	std::cout << "total_residue_num= " << total_residue_num << std::endl;

	for(Size n=1; n<= total_residue_num; n++){

			std::cout << "O2_star_torsion of " << Get_one_letter_name(RunTimeParameters::original_pose_residue_list[n].name);
			std::cout << RunTimeParameters::original_pose_residue_list[n].seq_num << "= " << O2_star_torsion_list[n] << std::endl;

	}

	return O2_star_torsion_list;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1 <Real>
get_native_rotamer(pose::Pose& template_pose, Size const & current_rebuild_num){


	bool const prepend_res = 	Is_prepend(current_rebuild_num);
	////////////////////Correct rotamer////////////////////////////////////////////////////////////////
	Size full_pose_res_reb = get_full_pose_reb_res(current_rebuild_num);

	Size suite_upper_res_num, suite_lower_res_num;
	if(prepend_res){
		suite_upper_res_num=full_pose_res_reb+1;
		suite_lower_res_num=full_pose_res_reb;
	} else {
		suite_upper_res_num=full_pose_res_reb;
		suite_lower_res_num=full_pose_res_reb-1;
	}

	conformation::Residue const & suite_lower_res=template_pose.residue(suite_lower_res_num);
	conformation::Residue const & suite_upper_res=template_pose.residue(suite_upper_res_num);
	utility::vector1 <Real > native_rotamer_torsions(13 /*# torsion angles*/, 9999);

			for(Size n=1; n<=1; n++){ //delta1
				Size i=n+3;
				native_rotamer_torsions[n]=suite_lower_res.mainchain_torsion(i);
			}

			for(Size n=2; n<=4; n++){ // chi_1, nu2_1, nu1_1 (recall that chi_1 2OH is determined seperately)
				Size i=n-1;
				native_rotamer_torsions[n]=suite_lower_res.chi(i);
			}

		  for(Size n=5; n<=6; n++){ //epsilon1, zeta1
				Size i=n;
				native_rotamer_torsions[n]=suite_lower_res.mainchain_torsion(i);
			}


			for(Size n=7; n<=10; n++){ //alpha2, beta2, gamma2, delta2
				Size i=n-6;
				native_rotamer_torsions[n]=suite_upper_res.mainchain_torsion(i);
			}

			for(Size n=11; n<=13; n++){ // chi_2, nu2_2, nu1_2
				Size i=n-10;
				native_rotamer_torsions[n]=suite_upper_res.chi(i);
			}

			return  native_rotamer_torsions;

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool
Is_NearNative_Rotamer(utility::vector1 <Real >& backbone_rotamer, utility::vector1 <Real >& native_rotamer, bool prepend_res){

			using namespace basic::options;
			using namespace basic::options::OptionKeys;
			Real cut_off=option[torsion_cutoff];

//			std::cout<< std::endl;



			for(Size i=1; i<=13 /*prepend*/; i++){//d1, c1, v1_1, v2_1, e1, g1, a2, b2, g2, d2, c2, v1_2, v2_2,

    		if((prepend_res==true) && (i==10 || i== 11)) continue; //don't evaluate d2 and c2 in this case
    		if((prepend_res==false) && (i==1 || i== 2)) continue; //don't evaluate d1 and c1 in this case
    		if(i==3 || i==4 || i==12 || i==13) continue; //v1_1, v2_1, v1_2 or v2_2;

				Real diff= std::min(std::abs(backbone_rotamer[i]-native_rotamer[i]+360), std::abs(backbone_rotamer[i]-native_rotamer[i]-360));
				diff = std::min(std::abs(backbone_rotamer[i]-native_rotamer[i]), diff);

//				std::cout << std::setw(spacing)<< diff << " ";

				if(diff > cut_off) return false;
			}
/*
			Size spacing=8;
				std::cout << std::setw(18)<<"native_torsion= ";
			for(Size i=1; i<=9; i++){
				if(i==3 || i==4) continue; //v1_1 and v2_1;
				std::cout << std::setw(spacing)<< native_rotamer[i] << " ";
			}
			std::cout<< std::endl;

			std::cout << std::setw(18)<<"current_torsion= ";
			for(Size i=1; i<=9; i++){
				if(i==3 || i==4) continue; //v1_1 and v2_1;
				std::cout << std::setw(spacing)<< backbone_rotamer[i] << " ";
			}
			std::cout<< std::endl;

			std::cout << std::setw(18)<<"diff_torsions= ";
*/

	return true;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

enum PuckerState{ ALL, NORTH, SOUTH };


///////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct	rotamer_ID_struct {
		Size d1;
		Size c1;
		Size e1;
		Size a2;
		Size z1;
		Size b2;
		Size g2;
		Size e1_std;
		Size a2_std;
		Size z1_std;
		Size b2_std;
		Size g2_std;
		Size group_count;
		Size subgroup_count;
};

struct	fine_rotamer_ID_struct {
//		Size d1;
		Size chi;
		Size e1;
		Size a2;
		Size z1;
		Size b2;
		Size g2;
		bool finish;
		Size count;
};



/////////////////////////////////////////////////////////////////////////////////////////////////////////
//Updated June 14, still need to verify if function is working properly
void
get_fine_rotamer(utility::vector1 <Real >& base_rotamer,utility::vector1 <Real >& current_rotamer, fine_rotamer_ID_struct& fine_rotamer_ID, Size current_rebuild_num){

			using namespace basic::options;
			using namespace basic::options::OptionKeys;

			bool const prepend_res = Is_prepend(current_rebuild_num);

			Real bin_size= option[fine_torsion_bin];
			Size bins=20/bin_size; //1,2 or 4



			if(fine_rotamer_ID.g2 == (bins)) {fine_rotamer_ID.g2=0; fine_rotamer_ID.b2++;}
			if(fine_rotamer_ID.b2 == (bins)) {fine_rotamer_ID.b2=0; fine_rotamer_ID.z1++;}
			if(fine_rotamer_ID.z1 == (bins)) {fine_rotamer_ID.z1=0; fine_rotamer_ID.a2++;}
			if(fine_rotamer_ID.a2 == (bins)) {fine_rotamer_ID.a2=0; fine_rotamer_ID.e1++;}
			if(fine_rotamer_ID.e1 == (bins)) {fine_rotamer_ID.e1=0; fine_rotamer_ID.chi++;}
			if(fine_rotamer_ID.chi == (bins)) {fine_rotamer_ID.finish=true;}

		//	return;

			//Sugar,sidechain i
			current_rotamer[1]= base_rotamer[1]; //delta i
			if(prepend_res==true){
				current_rotamer[2]= base_rotamer[2] + (fine_rotamer_ID.chi*bin_size);
			}else {
				current_rotamer[2]= base_rotamer[2];
			}
			current_rotamer[3]= base_rotamer[3]; //nu2 i
			current_rotamer[4]= base_rotamer[4]; // nu1 i

			//Backbone
			current_rotamer[5]= base_rotamer[5] +(fine_rotamer_ID.e1*bin_size);
			current_rotamer[6]= base_rotamer[6] +(fine_rotamer_ID.z1*bin_size);
			current_rotamer[7]= base_rotamer[7] +(fine_rotamer_ID.a2*bin_size);
			current_rotamer[8]= base_rotamer[8] +(fine_rotamer_ID.b2*bin_size);
			current_rotamer[9]= base_rotamer[9] +(fine_rotamer_ID.g2*bin_size);

			//Sugar,sidechain i+1
			current_rotamer[10]= base_rotamer[10]; //delta i+1
			if(prepend_res==false){
				current_rotamer[11]= base_rotamer[11] + (fine_rotamer_ID.chi*bin_size);
			}else {
				current_rotamer[11]= base_rotamer[11];
			}
			current_rotamer[12]= base_rotamer[12]; //nu2 i+1
			current_rotamer[13]= base_rotamer[13]; // nu1 i+1

			fine_rotamer_ID.count++;

/*
			std::cout << "bins= " << bins;
			std::cout << " g2= " << fine_rotamer_ID.g2;
			std::cout << " b2= " << fine_rotamer_ID.b2;
			std::cout << " z1= " << fine_rotamer_ID.z1;
			std::cout << " a2= " << fine_rotamer_ID.a2;
			std::cout << " e1= " << fine_rotamer_ID.e1;
			std::cout << " c1= " << fine_rotamer_ID.c1;
			std::cout << " count= " << fine_rotamer_ID.count;
			std::cout << " finish= " << fine_rotamer_ID.finish <<std::endl;
*/
			fine_rotamer_ID.g2++;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////

void
update_rotamer_ID(rotamer_ID_struct& rotamer_ID, Size& bins1, Size& bins2, Size& bins3, Size& bins4){



			if(rotamer_ID.g2_std == (bins3+1)) {rotamer_ID.g2_std=1; rotamer_ID.b2_std++;} else  return;
			if(rotamer_ID.b2_std == (bins1+1)) {rotamer_ID.b2_std=1; rotamer_ID.z1_std++;} else return;
			if(rotamer_ID.z1_std == (bins2+1)) {rotamer_ID.z1_std=1; rotamer_ID.a2_std++;} else return;
			if(rotamer_ID.a2_std == (bins3+1)) {rotamer_ID.a2_std=1; rotamer_ID.e1_std++;} else return;
			if(rotamer_ID.e1_std == (bins4+1)) {rotamer_ID.e1_std=1; rotamer_ID.g2++;} else return;
			rotamer_ID.group_count++; rotamer_ID.subgroup_count=0;
//			if(rotamer_ID.g2 == (3+1)) rotamer_ID.g2=1; rotamer_ID.b2++; else return;
//			if(rotamer_ID.b2 == (1+1)) rotamer_ID.b2=1; rotamer_ID.z1++; else return;
//			if(rotamer_ID.z1 == (2+1)) rotamer_ID.z1=1; rotamer_ID.a2++; else return;
//			if(rotamer_ID.a2 == (3+1)) rotamer_ID.a2=1; rotamer_ID.e1++; else return;
//			if(rotamer_ID.e1 == (3+1)) rotamer_ID.e1=1; rotamer_ID.c1++; else return;
			if(rotamer_ID.c1 == (bins4+1)) {rotamer_ID.c1=1; rotamer_ID.d1++;} else return;
			//d1 should not reach 3.


}


void
get_next_rotamer(rotamer_ID_struct& rotamer_ID, utility::vector1 <Real >& current_rotamer){

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	scoring::rna::Gaussian_parameter gp( 0.0, 0.0, 0.0);
	Real bin_size= option[torsion_bin];
	Size bins1=360/bin_size; //This is total bins
	Size bins2=180/bin_size; //This is total bins divided by 2;
	Size bins3=120/bin_size; //This is total bins divided by 3;
	Size bins4=40/bin_size+1; //This is the bin for chi and episilon, these two angle vary from -20+mean to 20+mean

	Real delta1, nu2_1, nu1_1, chi_1; //Sugar 1
	Real alpha2, beta2, gamma2, epsilon1, zeta1; //BB

	//Sugar 2, these are not initialized for now.
	Real delta2=0.0;
	Real nu1_2=0.0;
	Real nu2_2=0.0;
	Real chi_2=0.0;

	scoring::rna::RNA_TorsionPotential const rna_torsion_potential;

	update_rotamer_ID(rotamer_ID, bins1, bins2, bins3, bins4);
	rotamer_ID.g2_std++; // lowest level torsion.
	rotamer_ID.subgroup_count++;


		if (rotamer_ID.d1 == 1) {
			delta1 = rna_torsion_potential.gaussian_parameter_set_delta_north()[1].center;
			chi_1 = rna_torsion_potential.gaussian_parameter_set_chi_north()[1].center+bin_size*(rotamer_ID.c1-1)-20;
			nu2_1 = rna_torsion_potential.gaussian_parameter_set_nu2_north()[1].center;
			nu1_1 = rna_torsion_potential.gaussian_parameter_set_nu1_north()[1].center;

		}	else {
			delta1 = rna_torsion_potential.gaussian_parameter_set_delta_south()[1].center;
			chi_1 = rna_torsion_potential.gaussian_parameter_set_chi_south()[1].center+bin_size*(rotamer_ID.c1-1)-20;
			nu2_1 = rna_torsion_potential.gaussian_parameter_set_nu2_south()[1].center;
			nu1_1 = rna_torsion_potential.gaussian_parameter_set_nu1_south()[1].center;
		}

		if (rotamer_ID.d1 == 1) gp = rna_torsion_potential.gaussian_parameter_set_epsilon_north()[rotamer_ID.e1];
		else         gp = rna_torsion_potential.gaussian_parameter_set_epsilon_south()[rotamer_ID.e1];
		epsilon1 = gp.center +	bin_size*(rotamer_ID.e1_std-1)-20;


		if(rotamer_ID.a2==1) alpha2= (240) + bin_size * rotamer_ID.a2_std;
		if(rotamer_ID.a2==2) alpha2=  0 + bin_size * rotamer_ID.a2_std;
		if(rotamer_ID.a2==3) alpha2= 120 + bin_size * rotamer_ID.a2_std;

		if(rotamer_ID.z1==1) zeta1= 0 + bin_size * rotamer_ID.z1_std;
		if(rotamer_ID.z1==2) zeta1= 180 + bin_size * rotamer_ID.z1_std;

		beta2 = bin_size * rotamer_ID.b2_std;

		if(rotamer_ID.g2==1) gamma2= (240) + bin_size * rotamer_ID.g2_std;
		if(rotamer_ID.g2==2) gamma2= 0 + bin_size * rotamer_ID.g2_std;
		if(rotamer_ID.g2==3) gamma2= 120 + bin_size * rotamer_ID.g2_std;

		current_rotamer.clear(); //Would be better to make this a class, RNA_SuiteToSuiteRotamer



		current_rotamer.push_back( delta1 );    //1
		current_rotamer.push_back( chi_1 );     //2
		current_rotamer.push_back( nu2_1 );     //3
		current_rotamer.push_back( nu1_1 );     //4
		current_rotamer.push_back( epsilon1 );  //5
		current_rotamer.push_back( zeta1 );     //6
		current_rotamer.push_back( alpha2 );    //7
		current_rotamer.push_back( beta2 );     //8
		current_rotamer.push_back( gamma2 );    //9
		current_rotamer.push_back( delta2 );    //10
		current_rotamer.push_back( chi_2 );     //11
		current_rotamer.push_back( nu2_2 );     //12
		current_rotamer.push_back( nu1_2 );     //13
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void
get_full_rotamers_include_syn_chi( utility::vector1< utility::vector1 <utility::vector1 <Real > > >& backbone_rotamers_groups, PuckerState const & pucker1, PuckerState const & pucker2) {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;



	scoring::rna::RNA_TorsionPotential const rna_torsion_potential;


	backbone_rotamers_groups.clear();

	Real delta1, nu2_1, nu1_1, epsilon1, zeta1;
	Real alpha2, beta2, gamma2, delta2, nu2_2, nu1_2, chi_1, chi_2;

	scoring::rna::Gaussian_parameter gp( 0.0, 0.0, 0.0);

	Real bin_size= option[torsion_bin];
	Size bins1=360/bin_size; //This is total bins
	Size bins2=180/bin_size; //This is total bins divided by 2;
	Size bins3=120/bin_size; //This is total bins divided by 3;
	Size bins4=40/bin_size+1; //This is the bin for chi and episilon, these two angle vary from -20+mean to 20+mean


	std::cout << "bin_size= " << bin_size << std::endl;
	std::cout << "bins1= " << bins1 << std::endl;
	std::cout << "bins2= " << bins2 << std::endl;
	std::cout << "bins3= " << bins3 << std::endl;
	std::cout << "bins4= " << bins4 << std::endl;



	for (int d1 = 1; d1 <= 2; d1++ ) {
	for (Size c1 = 1; c1 <= 2; c1++ ) { //New to sample both syn and anti chi angle
	for (Size c1_std = 1; c1_std <= bins4; c1_std++ ) {

	for (Size d2 = 1; d2 <= 2; d2++ ) {
	for (Size c2 = 1; c2 <= 2; c2++ ) { //New to sample both syn and anti chi angle
	for (Size c2_std = 1; c2_std <= bins4; c2_std++ ) {


	if(pucker1 != 0){ //if false, then value of chi1, nu2_1, nu1_1 and delta_1 are not used...just select c1_std=1 and c1=1 value. However value of pucker (d1) does matter
			if (pucker1 != d1) continue;
			if (c1_std != 1) continue;
			if (c1!=1) continue;
	}



	if(pucker2 != 0){ //if false, then value of chi2, nu2_2, nu1_2 are not used...just select c2_std=1 and c2=1 value. However value of pucker (d2) does matter
			if (pucker2 != d2) continue;
			if (c2_std != 1) continue;
			if (c2!=1) continue;
	}


		if (d1 == 1) { //3' endo
			//For 3' endo delta angle is not statistically different for anti and syn conformation
			delta1 = rna_torsion_potential.gaussian_parameter_set_delta_north()[1].center;

			if(c1==1) { //anti
				chi_1 = rna_torsion_potential.gaussian_parameter_set_chi_north()[1].center+bin_size*(c1_std-1)-20;
				nu2_1 = rna_torsion_potential.gaussian_parameter_set_nu2_north()[1].center;
				nu1_1 = rna_torsion_potential.gaussian_parameter_set_nu1_north()[1].center;
			} else { //syn
//			Syn_3prime  count 21
//			mean nu_1= 107.792  mean_square nu_1= 11655.6  standard_deviation_nu_1 6.04656
//			mean nu_2= 31.8163  mean_square nu_2= 1026.34  standard_deviation_nu_2 3.75109
//			mean chi= -52.4573  mean_square chi= 3105.06  standard_deviation_chi 18.796
//			mean delta= 84.702  mean_square delta= 7213.8  standard_deviation_delta 6.27459
				chi_1 = -52.4573 + bin_size*(c1_std-1)-20;
				nu2_1 = 31.8163;
				nu1_1 = 107.792;
			}

		}	else { //2' endo

			//For 2' endo delta, nu2_1 and nu1_1 angle are not statistically different for anti and syn conformation
			delta1 = rna_torsion_potential.gaussian_parameter_set_delta_south()[1].center;
			nu2_1 = rna_torsion_potential.gaussian_parameter_set_nu2_south()[1].center;
			nu1_1 = rna_torsion_potential.gaussian_parameter_set_nu1_south()[1].center;

			if(c1==1){ //anti
				chi_1 = rna_torsion_potential.gaussian_parameter_set_chi_south()[1].center+bin_size*(c1_std-1)-20;
			}else{
//			Syn_2prime  count 27
//			mean nu_1= 157.725  mean_square nu_1= 24903.2  standard_deviation_nu_1 5.09893
//			mean nu_2= -37.4744  mean_square nu_2= 1408.15  standard_deviation_nu_2 1.95549
//			mean chi= -49.8518  mean_square chi= 2594.73  standard_deviation_chi 10.4654
//			mean delta= 149.451  mean_square delta= 22393.5  standard_deviation_delta 7.60976
				chi_1 = -49.8518+ bin_size*(c1_std-1)-20;
			}
		}



			if (d2 == 1) { //3' endo
				//For 3' endo delta angle is not statistically different for anti and syn conformation
				delta2 = rna_torsion_potential.gaussian_parameter_set_delta_north()[1].center;

				if(c2==1) { //anti
					chi_2 = rna_torsion_potential.gaussian_parameter_set_chi_north()[1].center+bin_size*(c2_std-1)-20;
					nu2_2 = rna_torsion_potential.gaussian_parameter_set_nu2_north()[1].center;
					nu1_2 = rna_torsion_potential.gaussian_parameter_set_nu1_north()[1].center;
				} else { //syn
//				Syn_3prime  count 21
//				mean nu_1= 107.792  mean_square nu_1= 11655.6  standard_deviation_nu_1 6.04656
//				mean nu_2= 31.8163  mean_square nu_2= 1026.34  standard_deviation_nu_2 3.75109
//				mean chi= -52.4573  mean_square chi= 3105.06  standard_deviation_chi 18.796
//				mean delta= 84.702  mean_square delta= 7213.8  standard_deviation_delta 6.27459
					chi_2 = -52.4573 + bin_size*(c1_std-1)-20;
					nu2_2 = 31.8163;
					nu1_2 = 107.792;
				}
			}	else {

					//For 2' endo delta, nu2_1 and nu1_1 angle are not statistically different for anti and syn conformation
					delta2 = rna_torsion_potential.gaussian_parameter_set_delta_south()[1].center;
					nu2_2 = rna_torsion_potential.gaussian_parameter_set_nu2_south()[1].center;
					nu1_2 = rna_torsion_potential.gaussian_parameter_set_nu1_south()[1].center;

					if(c2==1){ //anti
						chi_2 = rna_torsion_potential.gaussian_parameter_set_chi_south()[1].center+bin_size*(c2_std-1)-20;
					}else{ //syn
//					Syn_2prime  count 27
//					mean nu_1= 157.725  mean_square nu_1= 24903.2  standard_deviation_nu_1 5.09893
//					mean nu_2= -37.4744  mean_square nu_2= 1408.15  standard_deviation_nu_2 1.95549
//					mean chi= -49.8518  mean_square chi= 2594.73  standard_deviation_chi 10.4654
//					mean delta= 149.451  mean_square delta= 22393.5  standard_deviation_delta 7.60976
						chi_2 = -49.8518+ bin_size*(c2_std-1)-20;
					}
			}

			Size spacing=12;

			std::cout << " delta1= " << std::setw(spacing) << delta1;
			std::cout << " chi_1= " << std::setw(spacing) << chi_1;
			std::cout << " nu2_1= " << std::setw(spacing) << nu2_1;
			std::cout << " nu1_1= " << std::setw(spacing) << nu1_1;
			std::cout << " delta2= " << std::setw(spacing) << delta2;
			std::cout << " chi_2= " << std::setw(spacing) << chi_2;
			std::cout << " nu2_2= " << std::setw(spacing) << nu2_2;
			std::cout << " nu1_2= " << std::setw(spacing) << nu1_2;
			std::cout << std::endl;

				//epsilon depends on delta of same residue...
				for (Size e1 = 1; e1 <= 1; e1++ ) {
					for (Size a2 = 1; a2 <= 3; a2++ ) {
					  //zeta depends on alpha of next residue...
						for (Size z1 = 1; z1 <= 2; z1++ ) {
							for (Size b2 = 1; b2 <= 1; b2++ ) {
								for (Size g2 = 1; g2 <= 3; g2++ ) {

									utility::vector1< utility::vector1 <Real > > backbone_rotamers_current;



//More Rotamers
/////////////////////////////////////////////////



			for ( int e1_std = 1; e1_std <= bins4; e1_std++ ) {


				if (d1 == 1) gp = rna_torsion_potential.gaussian_parameter_set_epsilon_north()[e1];
				else         gp = rna_torsion_potential.gaussian_parameter_set_epsilon_south()[e1];

				epsilon1 = gp.center +	bin_size*(e1_std-1)-20;

				for ( Size a2_std = 1; a2_std <= bins3; a2_std++ ) {


						if(a2==1) alpha2= (240) + bin_size * a2_std;
						if(a2==2) alpha2=  0 + bin_size * a2_std;
						if(a2==3) alpha2= 120 + bin_size * a2_std;

//						if(alpha2>=180) alpha2=alpha2-360;

						for ( Size z1_std = 1; z1_std <= bins2; z1_std++ ) {

							if(z1==1) zeta1= 0 + bin_size * z1_std;
							if(z1==2) zeta1= 180 + bin_size * z1_std;


								for ( Size b2_std = 1; b2_std <= bins1; b2_std++ ) {

										beta2 = bin_size * b2_std;

										for ( Size g2_std = 1; g2_std <= bins3; g2_std++ ) {

											if(g2==1) gamma2= (240) + bin_size * g2_std;
											if(g2==2) gamma2= 0 + bin_size * g2_std;
											if(g2==3) gamma2= 120 + bin_size * g2_std;

													utility::vector1 < Real >  backbone_rotamer; //Would be better to make this a class, RNA_SuiteToSuiteRotamer


													backbone_rotamer.push_back( delta1 );    //1
													backbone_rotamer.push_back( chi_1 );     //2
													backbone_rotamer.push_back( nu2_1 );     //3
													backbone_rotamer.push_back( nu1_1 );     //4
													backbone_rotamer.push_back( epsilon1 );  //5
													backbone_rotamer.push_back( zeta1 );     //6
													backbone_rotamer.push_back( alpha2 );    //7
													backbone_rotamer.push_back( beta2 );     //8
													backbone_rotamer.push_back( gamma2 );    //9
													backbone_rotamer.push_back( delta2 );    //10
													backbone_rotamer.push_back( chi_2 );     //11
													backbone_rotamer.push_back( nu2_2 );     //12
													backbone_rotamer.push_back( nu1_2 );     //13

/*
		Size spacing=8;

		std::cout <<  std::setw(18) << "Torsions= ";


		std::cout << std::setw(spacing)<< backbone_rotamer[1] << " ";
		std::cout << std::setw(spacing)<< backbone_rotamer[2] << " ";
		std::cout << std::setw(spacing)<< backbone_rotamer[3] << " ";
		std::cout << std::setw(spacing)<< backbone_rotamer[4] << " ";
		std::cout << std::setw(spacing)<< backbone_rotamer[5] << " ";
		std::cout << std::setw(spacing)<< backbone_rotamer[6] << " ";
		std::cout << std::setw(spacing)<< backbone_rotamer[7] << " ";
		std::cout << std::setw(spacing)<< backbone_rotamer[8] << " ";
		std::cout << std::setw(spacing)<< backbone_rotamer[9] << " ";
		std::cout << std::setw(spacing)<< backbone_rotamer[10] << " ";
		std::cout << std::setw(spacing)<< backbone_rotamer[11] << " ";
		std::cout << std::setw(spacing)<< backbone_rotamer[12] << " ";
		std::cout << std::setw(spacing)<< backbone_rotamer[13] << " ";
		std::cout << std::setw(spacing)<< backbone_rotamer[14] << " ";
		std::cout << std::endl;
*/


													backbone_rotamers_current.push_back( backbone_rotamer );

											  } // gamma2_samples
										  } // beta2_samples
									  } // zeta1_samples
								  } //  alpha2_samples
							  } // epsilon1_samples

							  	backbone_rotamers_groups.push_back(backbone_rotamers_current);

								} // gamma2
							} // beta2
						} // zeta1
					} // alpha2
				} // epsilon1
				} // chi2_std
				}// chi2
			} // delta2
		}	// chi1_std
		}// chi1
	} // delta1

	std::cout << "Exit get_full_rotamers_include_syn_chi function" << std::endl;
}





//Torsions:
//Bin size is always 20 degree
//alpha: Sample 360 degree in bin_size interval
//beta: Sample 360 degree in bin_size interval
//gamma: Sample 360 degree in bin_size interval
//delta: Sample two torsion North (3'-endo) or South (2'-endo)
//eplison: Sample center, center+bin_size, center-bin_size. Center is determine by sugar pucker (actaully the way code is written, as bin_size decrease with sample more point to cover same interval)
//zeta: Sample 360 degree in bin_size interval
//chi: Sample only anti conformation. Center of anti, center+bin_size, center-bin_size
//nu2 and nu1, sample center. Center is determine by sugar pucker
//chi_O2star: currently no sampled

void
get_full_rotamers_exclude_syn_chi( utility::vector1< utility::vector1 <utility::vector1 <Real > > >& backbone_rotamers_groups, PuckerState const & pucker1, PuckerState const & pucker2) {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;



	scoring::rna::RNA_TorsionPotential const rna_torsion_potential;


	backbone_rotamers_groups.clear();

	Real delta1, nu2_1, nu1_1, epsilon1, zeta1;
	Real alpha2, beta2, gamma2, delta2, nu2_2, nu1_2, chi_1, chi_2;

	scoring::rna::Gaussian_parameter gp( 0.0, 0.0, 0.0);

	Real bin_size= option[torsion_bin];
	Size bins1=360/bin_size; //This is total bins
	Size bins2=180/bin_size; //This is total bins divided by 2;
	Size bins3=120/bin_size; //This is total bins divided by 3;
	Size bins4=40/bin_size+1; //This is the bin for chi and episilon, these two angle vary from -20+mean to 20+mean


	std::cout << "bin_size= " << bin_size << std::endl;
	std::cout << "bins1= " << bins1 << std::endl;
	std::cout << "bins2= " << bins2 << std::endl;
	std::cout << "bins3= " << bins3 << std::endl;
	std::cout << "bins4= " << bins4 << std::endl;



	for (int d1 = 1; d1 <= 2; d1++ ) {
	for (Size c1 = 1; c1 <= bins4; c1++ ) {

		if ( pucker1 != 0 && (pucker1 != d1 || ((int)(bin_size*(c1-1)-20)) ) )  continue;


		if (d1 == 1) {
			delta1 = rna_torsion_potential.gaussian_parameter_set_delta_north()[1].center;
			chi_1 = rna_torsion_potential.gaussian_parameter_set_chi_north()[1].center+bin_size*(c1-1)-20;
			nu2_1 = rna_torsion_potential.gaussian_parameter_set_nu2_north()[1].center;
			nu1_1 = rna_torsion_potential.gaussian_parameter_set_nu1_north()[1].center;

		}	else {
			delta1 = rna_torsion_potential.gaussian_parameter_set_delta_south()[1].center;
			chi_1 = rna_torsion_potential.gaussian_parameter_set_chi_south()[1].center+bin_size*(c1-1)-20;
			nu2_1 = rna_torsion_potential.gaussian_parameter_set_nu2_south()[1].center;
			nu1_1 = rna_torsion_potential.gaussian_parameter_set_nu1_south()[1].center;
		}

		for (Size d2 = 1; d2 <= 2; d2++ ) {
		for (Size c2 = 1; c2 <= bins4; c2++ ) {

//			std::cout << " d2= " << d2;
//			std::cout << " c2= " << c2;
//			std::cout << " (c2-1)-20)= " << bin_size*(c2-1)-20;
//			std::cout << " Size((c2-1)-20)))= " << ((int)(bin_size*(c2-1)-20)) << std::endl;

			if ( pucker2 != 0 &&  (pucker2 != d2 || ((int)(bin_size*(c2-1)-20)) ) ) continue;


			if (d2 == 1) {
				delta2 = rna_torsion_potential.gaussian_parameter_set_delta_north()[1].center;
			  chi_2 = rna_torsion_potential.gaussian_parameter_set_chi_north()[1].center+bin_size*(c2-1)-20;
				nu2_2 = rna_torsion_potential.gaussian_parameter_set_nu2_north()[1].center;
				nu1_2 = rna_torsion_potential.gaussian_parameter_set_nu1_north()[1].center;
			}	else {
				delta2 = rna_torsion_potential.gaussian_parameter_set_delta_south()[1].center;
				chi_2 = rna_torsion_potential.gaussian_parameter_set_chi_south()[1].center+bin_size*(c2-1)-20;
				nu2_2 = rna_torsion_potential.gaussian_parameter_set_nu2_south()[1].center;
				nu1_2 = rna_torsion_potential.gaussian_parameter_set_nu1_south()[1].center;
			}

			Size spacing=6;

			std::cout << " delta1= " << std::setw(spacing) << delta1;
			std::cout << " chi_1= " << std::setw(spacing) << chi_1;
			std::cout << " nu2_1= " << std::setw(spacing) << nu2_1;
			std::cout << " nu1_1= " << std::setw(spacing) << nu1_1;
			std::cout << " delta2= " << std::setw(spacing) << delta2;
			std::cout << " chi_2= " << std::setw(spacing) << chi_2;
			std::cout << " nu2_2= " << std::setw(spacing) << nu2_2;
			std::cout << " nu1_2= " << std::setw(spacing) << nu1_2;
			std::cout << std::endl;

				//epsilon depends on delta of same residue...
				for (Size e1 = 1; e1 <= 1; e1++ ) {
					for (Size a2 = 1; a2 <= 3; a2++ ) {
					  //zeta depends on alpha of next residue...
						for (Size z1 = 1; z1 <= 2; z1++ ) {
							for (Size b2 = 1; b2 <= 1; b2++ ) {
								for (Size g2 = 1; g2 <= 3; g2++ ) {

									utility::vector1< utility::vector1 <Real > > backbone_rotamers_current;



//More Rotamers
/////////////////////////////////////////////////



			for ( int e1_std = 1; e1_std <= bins4; e1_std++ ) {


				if (d1 == 1) gp = rna_torsion_potential.gaussian_parameter_set_epsilon_north()[e1];
				else         gp = rna_torsion_potential.gaussian_parameter_set_epsilon_south()[e1];

				epsilon1 = gp.center +	bin_size*(e1_std-1)-20;

				for ( Size a2_std = 1; a2_std <= bins3; a2_std++ ) {


						if(a2==1) alpha2= (240) + bin_size * a2_std;
						if(a2==2) alpha2=  0 + bin_size * a2_std;
						if(a2==3) alpha2= 120 + bin_size * a2_std;

//						if(alpha2>=180) alpha2=alpha2-360;

						for ( Size z1_std = 1; z1_std <= bins2; z1_std++ ) {

							if(z1==1) zeta1= 0 + bin_size * z1_std;
							if(z1==2) zeta1= 180 + bin_size * z1_std;


								for ( Size b2_std = 1; b2_std <= bins1; b2_std++ ) {

										beta2 = bin_size * b2_std;

										for ( Size g2_std = 1; g2_std <= bins3; g2_std++ ) {

											if(g2==1) gamma2= (240) + bin_size * g2_std;
											if(g2==2) gamma2= 0 + bin_size * g2_std;
											if(g2==3) gamma2= 120 + bin_size * g2_std;

													utility::vector1 < Real >  backbone_rotamer; //Would be better to make this a class, RNA_SuiteToSuiteRotamer


													backbone_rotamer.push_back( delta1 );    //1
													backbone_rotamer.push_back( chi_1 );     //2
													backbone_rotamer.push_back( nu2_1 );     //3
													backbone_rotamer.push_back( nu1_1 );     //4
													backbone_rotamer.push_back( epsilon1 );  //5
													backbone_rotamer.push_back( zeta1 );     //6
													backbone_rotamer.push_back( alpha2 );    //7
													backbone_rotamer.push_back( beta2 );     //8
													backbone_rotamer.push_back( gamma2 );    //9
													backbone_rotamer.push_back( delta2 );    //10
													backbone_rotamer.push_back( chi_2 );     //11
													backbone_rotamer.push_back( nu2_2 );     //12
													backbone_rotamer.push_back( nu1_2 );     //13

/*
		Size spacing=8;

		std::cout <<  std::setw(18) << "Torsions= ";


		std::cout << std::setw(spacing)<< backbone_rotamer[1] << " ";
		std::cout << std::setw(spacing)<< backbone_rotamer[2] << " ";
		std::cout << std::setw(spacing)<< backbone_rotamer[3] << " ";
		std::cout << std::setw(spacing)<< backbone_rotamer[4] << " ";
		std::cout << std::setw(spacing)<< backbone_rotamer[5] << " ";
		std::cout << std::setw(spacing)<< backbone_rotamer[6] << " ";
		std::cout << std::setw(spacing)<< backbone_rotamer[7] << " ";
		std::cout << std::setw(spacing)<< backbone_rotamer[8] << " ";
		std::cout << std::setw(spacing)<< backbone_rotamer[9] << " ";
		std::cout << std::setw(spacing)<< backbone_rotamer[10] << " ";
		std::cout << std::setw(spacing)<< backbone_rotamer[11] << " ";
		std::cout << std::setw(spacing)<< backbone_rotamer[12] << " ";
		std::cout << std::setw(spacing)<< backbone_rotamer[13] << " ";
		std::cout << std::setw(spacing)<< backbone_rotamer[14] << " ";
		std::cout << std::endl;
*/


													backbone_rotamers_current.push_back( backbone_rotamer );

											  } // gamma2_samples
										  } // beta2_samples
									  } // zeta1_samples
								  } //  alpha2_samples
							  } // epsilon1_samples

							  	backbone_rotamers_groups.push_back(backbone_rotamers_current);

								} // gamma2
							} // beta2
						} // zeta1
					} // alpha2
				} // epsilon1
				}// chi2
			} // delta2
		}// chi1
	} // delta1

	std::cout << "Exit get_full_rotamer_exclude_syn_chi function" << std::endl;
}



void
get_rotamers( utility::vector1< utility::vector1 <utility::vector1 <Real > > >& backbone_rotamers_groups, PuckerState const & pucker1, PuckerState const & pucker2) {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	std::string rotamer_set_option = option[rotamer_set];

	if(rotamer_set_option=="exclude_syn_chi"){
		get_full_rotamers_exclude_syn_chi( backbone_rotamers_groups, pucker1, pucker2);
	} else {
		get_full_rotamers_include_syn_chi( backbone_rotamers_groups, pucker1, pucker2);
	}
}

//Actually create the backbone_rotamers_groups...computationally expensive don't call often
Size
get_total_rotamer_group(){
		std::cout << "get total rotamer group by creating an actual backbone_rotamers_group object !!! " << std::endl;

		utility::vector1< utility::vector1 <utility::vector1 <Real > > > backbone_rotamers_groups;
		get_rotamers( backbone_rotamers_groups, ALL, NORTH); //Assume that any (pucker1, pucker2) will create backbone_rotamers_groups with same total_rotamer_group value
		Size total_rotamer_group=backbone_rotamers_groups.size();
		std::cout << "total_rotamer_group= " << total_rotamer_group << std::endl;
		return total_rotamer_group;

}


///////////////////////////////////////////////////////////////////////////////

void
get_backbone_rotamers( utility::vector1< utility::vector1 <utility::vector1 <Real > > >&
											 backbone_rotamers_groups,
											 PuckerState const & pucker1,
											 PuckerState const & pucker2 ) {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	scoring::rna::RNA_TorsionPotential const rna_torsion_potential;

	Real bin_size=20;
	Size bins1=360/bin_size; //This is total bins
	Size bins2=180/bin_size; //This is total bins divided by 2;
	Size bins3=120/bin_size; //This is total bins divided by 3;

	std::cout << "bin_size= " << bin_size << std::endl;
	std::cout << "bins1= " << bins1 << std::endl;
	std::cout << "bins2= " << bins2 << std::endl;
	std::cout << "bins3= " << bins3 << std::endl;


	utility::vector1< Real > torsion_samples;
	torsion_samples.push_back( -0.5 );
	torsion_samples.push_back( 0.5 );


	backbone_rotamers_groups.clear();

	Real delta1, nu2_1, nu1_1, epsilon1, zeta1;
	Real alpha2, beta2, gamma2, delta2, nu2_2, nu1_2, chi_1, chi_2;

	scoring::rna::Gaussian_parameter gp( 0.0, 0.0, 0.0);




	for (int d1 = 1; d1 <= 2; d1++ ) {

		for (int c1 = -1; c1 <= 1; c1++ ) {

		if ( pucker1 > 0 &&  pucker1 != d1 ) continue;


		if (d1 == 1) {
			delta1 = rna_torsion_potential.gaussian_parameter_set_delta_north()[1].center;
			chi_1 = rna_torsion_potential.gaussian_parameter_set_chi_north()[1].center+bin_size*c1;
			nu2_1 = rna_torsion_potential.gaussian_parameter_set_nu2_north()[1].center;
			nu1_1 = rna_torsion_potential.gaussian_parameter_set_nu1_north()[1].center;

		}	else {
			delta1 = rna_torsion_potential.gaussian_parameter_set_delta_south()[1].center;
			chi_1 = rna_torsion_potential.gaussian_parameter_set_chi_south()[1].center+bin_size*c1;
			nu2_1 = rna_torsion_potential.gaussian_parameter_set_nu2_south()[1].center;
			nu1_1 = rna_torsion_potential.gaussian_parameter_set_nu1_south()[1].center;
		}

		for (int d2 = 1; d2 <= 2; d2++ ) {

			if ( pucker2 > 0 &&  pucker2 != d2 ) continue;

			if (d2 == 1) {
				delta2 = rna_torsion_potential.gaussian_parameter_set_delta_north()[1].center;
			  chi_2 = rna_torsion_potential.gaussian_parameter_set_chi_north()[1].center;
				nu2_2 = rna_torsion_potential.gaussian_parameter_set_nu2_north()[1].center;
				nu1_2 = rna_torsion_potential.gaussian_parameter_set_nu1_north()[1].center;
			}	else {
				delta2 = rna_torsion_potential.gaussian_parameter_set_delta_south()[1].center;
				chi_2 = rna_torsion_potential.gaussian_parameter_set_chi_south()[1].center;
				nu2_2 = rna_torsion_potential.gaussian_parameter_set_nu2_south()[1].center;
				nu1_2 = rna_torsion_potential.gaussian_parameter_set_nu1_south()[1].center;
			}


				//epsilon depends on delta of same residue...
				for (Size e1 = 1; e1 <= 1; e1++ ) {
					for (Size a2 = 1; a2 <= 3; a2++ ) {
					  //zeta depends on alpha of next residue...
						for (Size z1 = 1; z1 <= 2; z1++ ) {
							for (Size b2 = 1; b2 <= 1; b2++ ) {
								for (Size g2 = 1; g2 <= 3; g2++ ) {

									utility::vector1< utility::vector1 <Real > > backbone_rotamers_current;



//More Rotamers
/////////////////////////////////////////////////

			for ( int e1_std = -1; e1_std <= 1; e1_std++ ) {

				if (d1 == 1) gp = rna_torsion_potential.gaussian_parameter_set_epsilon_north()[e1];
				else         gp = rna_torsion_potential.gaussian_parameter_set_epsilon_south()[e1];

				epsilon1 = gp.center +	e1_std * bin_size;

				for ( Size a2_std = 1; a2_std <= torsion_samples.size(); a2_std++ ) {
						gp = rna_torsion_potential.gaussian_parameter_set_alpha()[a2];
						alpha2 = gp.center +	torsion_samples[ a2_std ] * gp.width;

						for ( Size z1_std = 1; z1_std <= torsion_samples.size(); z1_std++ ) {
								if (a2 == 1)    gp = rna_torsion_potential.gaussian_parameter_set_zeta_alpha_sc_minus()[z1];
								else if (a2==2) gp = rna_torsion_potential.gaussian_parameter_set_zeta_alpha_sc_plus()[z1];
								else            gp = rna_torsion_potential.gaussian_parameter_set_zeta_alpha_ap()[z1];
								zeta1 = gp.center +	torsion_samples[ z1_std ] * gp.width;

								for ( Size b2_std = 1; b2_std <= torsion_samples.size(); b2_std++ ) {

										gp = rna_torsion_potential.gaussian_parameter_set_beta()[b2];
										beta2 = gp.center +	torsion_samples[ b2_std ] * gp.width;

										for ( Size g2_std = 1; g2_std <= torsion_samples.size(); g2_std++ ) {
												gp = rna_torsion_potential.gaussian_parameter_set_gamma()[g2];
												gamma2 = gp.center +	torsion_samples[ g2_std ] * gp.width;




													utility::vector1 < Real >  backbone_rotamer; //Would be better to make this a class, RNA_SuiteToSuiteRotamer





													backbone_rotamer.push_back( delta1 );    //1
													backbone_rotamer.push_back( chi_1 );     //2
													backbone_rotamer.push_back( nu2_1 );     //3
													backbone_rotamer.push_back( nu1_1 );     //4
													backbone_rotamer.push_back( epsilon1 );  //5
													backbone_rotamer.push_back( zeta1 );     //6
													backbone_rotamer.push_back( alpha2 );    //7
													backbone_rotamer.push_back( beta2 );     //8
													backbone_rotamer.push_back( gamma2 );    //9
													backbone_rotamer.push_back( delta2 );    //10
													backbone_rotamer.push_back( chi_2 );     //11
													backbone_rotamer.push_back( nu2_2 );     //12
													backbone_rotamer.push_back( nu1_2 );     //13


/*

		Size spacing=8;

		std::cout <<  std::setw(18) << "Torsions= ";


		std::cout << std::setw(spacing)<< backbone_rotamer[7] << " ";
		std::cout << std::setw(spacing)<< backbone_rotamer[8] << " ";
		std::cout << std::setw(spacing)<< backbone_rotamer[9] << " ";
		std::cout << std::setw(spacing)<< backbone_rotamer[1] << " ";
		std::cout << std::setw(spacing)<< backbone_rotamer[5] << " ";
		std::cout << std::setw(spacing)<< backbone_rotamer[6] << " ";
		std::cout << std::setw(spacing)<< backbone_rotamer[2] << " ";
		std::cout << std::setw(spacing)<< backbone_rotamer[3] << " ";
		std::cout << std::setw(spacing)<< backbone_rotamer[4] << " ";
		std::cout << std::endl;
*/

													backbone_rotamers_current.push_back( backbone_rotamer );

											  } // gamma2_samples
										  } // beta2_samples
									  } // zeta1_samples
								  } //  alpha2_samples
							  } // epsilon1_samples

							  	backbone_rotamers_groups.push_back(backbone_rotamers_current);

								} // gamma2
							} // beta2
						} // zeta1
					} // alpha2
				} // epsilon1
			} // delta2
		}// chi1
	} // delta1

}



///////////////////////////////////////////////////////////////////////////////

void
apply_rotamer( pose::Pose & pose, Size const res_reb, Size current_rebuild_num ,utility::vector1< Real >& backbone_rotamer)
{
	using namespace core::id;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	bool const prepend_res = Is_prepend(current_rebuild_num);

	Size suite_upper_res_num, suite_lower_res_num;

	bool pucker1, pucker2;
		if (prepend_res) {
			pucker1=true;
			pucker2=false;
			suite_upper_res_num=res_reb+1;
			suite_lower_res_num=res_reb;
		} else {
			pucker1=false;
			pucker2=true;
			suite_upper_res_num=res_reb;
			suite_lower_res_num=res_reb-1;
		}


	if ( pucker1 ) {
		pose.set_torsion( TorsionID( suite_lower_res_num  , id::BB,  4 ), backbone_rotamer[1] );//delta i
		pose.set_torsion( TorsionID( suite_lower_res_num  , id::CHI, 1 ), backbone_rotamer[2] );//chi i
		pose.set_torsion( TorsionID( suite_lower_res_num  , id::CHI, 2 ), backbone_rotamer[3] );//nu2 i
		pose.set_torsion( TorsionID( suite_lower_res_num  , id::CHI, 3 ), backbone_rotamer[4] );//nu1 i
//		pose.set_torsion( TorsionID( suite_lower_res_num  , id::CHI, 4 ), 0.0 ); //Hydrogen of O2star.... currently not sampled, arbitrarily set to zero. June 13, 2009 (before this torsion angle was random)
	}

	pose.set_torsion( TorsionID( suite_lower_res_num  , id::BB,  5 ), backbone_rotamer[5] ); //epsilon i
	pose.set_torsion( TorsionID( suite_lower_res_num  , id::BB,  6 ), backbone_rotamer[6] ); //zetta i
	pose.set_torsion( TorsionID( suite_upper_res_num, id::BB,  1 ), backbone_rotamer[7] ); //alpha i+1
	pose.set_torsion( TorsionID( suite_upper_res_num, id::BB,  2 ), backbone_rotamer[8] ); //betta i+1
	pose.set_torsion( TorsionID( suite_upper_res_num, id::BB,  3 ), backbone_rotamer[9] ); //gamma i+1

	if ( pucker2 ) {
		pose.set_torsion( TorsionID( suite_upper_res_num, id::BB,  4 ), backbone_rotamer[10] ); //delta i+1
		pose.set_torsion( TorsionID( suite_upper_res_num, id::CHI, 1 ), backbone_rotamer[11] );//chi i+1
		pose.set_torsion( TorsionID( suite_upper_res_num, id::CHI, 2 ), backbone_rotamer[12] ); //nu2 i+1
		pose.set_torsion( TorsionID( suite_upper_res_num, id::CHI, 3 ), backbone_rotamer[13] ); //nu1 i+1
//		pose.set_torsion( TorsionID( suite_upper_res_num, id::CHI, 4 ), 0.0 ); //Hydrogen of O2star.... currently not sampled, arbitrarily set to zero. June 13, 2009 (before this torsion angle was random)
	}



/*
		Size spacing=8;

		std::cout <<  std::setw(18) << "Apply Torsions= ";


		std::cout << std::setw(spacing)<< backbone_rotamer[7] << " ";
		std::cout << std::setw(spacing)<< backbone_rotamer[8] << " ";
		std::cout << std::setw(spacing)<< backbone_rotamer[9] << " ";
		std::cout << std::setw(spacing)<< backbone_rotamer[(prepend_res ? 1: 10)] << " ";
		std::cout << std::setw(spacing)<< backbone_rotamer[5] << " ";
		std::cout << std::setw(spacing)<< backbone_rotamer[6] << " ";
		std::cout << std::setw(spacing)<< backbone_rotamer[(prepend_res ? 2: 11)] << " ";
		std::cout << std::setw(spacing)<< backbone_rotamer[(prepend_res ? 3: 12)] << " ";
		std::cout << std::setw(spacing)<< backbone_rotamer[(prepend_res ? 4: 13)] << " ";
		std::cout << std::endl;



	std::cout << "Check point apply rotamer 1 " << std::endl;
*/
}



///////////////////////////////////////////////////////////////////////////////
/*

void
dinucleotide_test()
{
	// Read in pdb.
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace core::kinematics;
	using namespace core::optimization;
	using namespace core::io::silent;
	using namespace core::id;
	using namespace protocols::rna;

	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "rna" );

	pose::Pose pose;
	std::string sequence = "cc";

	core::pose::make_pose_from_sequence(
																					pose,
																					sequence,
																					*rsd_set );
	dump_pdb( pose, "start.pdb" );

	scoring::rna::RNA_TorsionPotential const rna_torsion_potential;
	for (Size i = 1; i <=2; i++ ) { // starting values for torsions.

		pose.set_torsion( TorsionID( i, id::BB, 1 ), rna_torsion_potential.gaussian_parameter_set_alpha()[1].center );
		pose.set_torsion( TorsionID( i, id::BB, 2 ), rna_torsion_potential.gaussian_parameter_set_beta()[1].center );
		pose.set_torsion( TorsionID( i, id::BB, 3 ), rna_torsion_potential.gaussian_parameter_set_gamma()[1].center );
		pose.set_torsion( TorsionID( i, id::BB, 4 ), rna_torsion_potential.gaussian_parameter_set_delta_north()[1].center );
		pose.set_torsion( TorsionID( i, id::BB, 5 ), rna_torsion_potential.gaussian_parameter_set_epsilon_north()[1].center );
		pose.set_torsion( TorsionID( i, id::BB, 6 ), rna_torsion_potential.gaussian_parameter_set_zeta_alpha_sc_minus()[1].center );

		pose.set_torsion( TorsionID( i, id::CHI, 1 ), rna_torsion_potential.gaussian_parameter_set_chi_north()[1].center );
		pose.set_torsion( TorsionID( i, id::CHI, 2 ), rna_torsion_potential.gaussian_parameter_set_nu2_north()[1].center );
		pose.set_torsion( TorsionID( i, id::CHI, 3 ), rna_torsion_potential.gaussian_parameter_set_nu1_north()[1].center );
		pose.set_torsion( TorsionID( i, id::CHI, 4 ), 0.0 );


	}

	pose::Pose current_pose( pose );

	protocols::viewer::add_conformation_viewer( current_pose.conformation(), "current", 400, 400 );

	ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( RNA_HIRES_WTS );
	scorefxn->set_weight( dihedral_constraint, 1.0 );

	SilentFileData silent_file_data;
	std::string const silent_file = option[ out::file::silent  ]();

	// Now cycle through one of the torsions. Score and output.
	Size count( 1 );

	utility::vector1< utility::vector1 <Real > > backbone_rotamers;
// Since I change the get_backbone_rotamers function
//	get_backbone_rotamers( backbone_rotamers, ALL, ALL );

	for ( Size n = 1;  n <= backbone_rotamers.size(); n++ )		{
		apply_rotamer( pose, 1, backbone_rotamers[n] );

		current_pose = pose;

		ConstraintSetOP cst_set( new ConstraintSet() );
		create_dihedral_constraint( cst_set, current_pose, TorsionID( 1, id::BB, 4 ) );
		create_dihedral_constraint( cst_set, current_pose, TorsionID( 1, id::CHI, 1 ) );
		create_dihedral_constraint( cst_set, current_pose, TorsionID( 1, id::CHI, 2 ) );
		create_dihedral_constraint( cst_set, current_pose, TorsionID( 1, id::CHI, 3 ) );
		create_dihedral_constraint( cst_set, current_pose, TorsionID( 1, id::BB, 5 ) );
		create_dihedral_constraint( cst_set, current_pose, TorsionID( 1, id::BB, 6 ) );
		create_dihedral_constraint( cst_set, current_pose, TorsionID( 2, id::BB, 1 ) );
		create_dihedral_constraint( cst_set, current_pose, TorsionID( 2, id::BB, 2 ) );
		create_dihedral_constraint( cst_set, current_pose, TorsionID( 2, id::BB, 3 ) );
		create_dihedral_constraint( cst_set, current_pose, TorsionID( 2, id::BB, 4 ) );
		create_dihedral_constraint( cst_set, current_pose, TorsionID( 2, id::CHI, 1 ) );
		create_dihedral_constraint( cst_set, current_pose, TorsionID( 2, id::CHI, 2 ) );
		create_dihedral_constraint( cst_set, current_pose, TorsionID( 2, id::CHI, 3 ) );
		current_pose.constraint_set( cst_set );

		AtomTreeMinimizer minimizer;
		float const dummy_tol( 0.0000025);
		bool const use_nblist( true );
		MinimizerOptions options( "dfpmin", dummy_tol, use_nblist, false, false );
		options.nblist_auto_update( true );

		kinematics::MoveMap mm;
		mm.set_bb(  true );
		mm.set_chi( true );
		mm.set_jump( true );

		minimizer.run( current_pose, mm, *(scorefxn), options );
		//			(*scorefxn)(current_pose);

		std::string const tag( "S_"+lead_zero_string_of(count++, 3) );
		RNA_SilentStruct s( current_pose, tag );
		silent_file_data.write_silent_struct(s, silent_file, false);
		dump_pdb( current_pose, tag+".pdb" );
		dump_pdb( pose, "nomin_"+tag+".pdb" );
	}



}
*/

////////////////////////////////////////////////////
void
vary_bond_length( pose::Pose & pose,
									core::id::TorsionID & tor_id,
									core::kinematics::MoveMap & mm,
									core::scoring::constraints::ConstraintSetOP & cst_set )
{
	using namespace core::id;
	using namespace core::scoring;
	using namespace core::scoring::constraints;

	AtomID id1,id2,id3,id4, my_ID;
	bool const failure = pose.conformation().get_torsion_angle_atom_ids( tor_id, id1, id2, id3, id4 );
	if (failure) return;

 	core::kinematics::tree::Atom const * atom2 ( & pose.atom_tree().atom( id2 ) );
	core::kinematics::tree::Atom const * atom3 ( & pose.atom_tree().atom( id3 ) );

	DOF_ID dof_id;
	if ( atom2->parent() == atom3 ) {
		my_ID = id2;
	} else if ( atom3->parent() == atom2 ) {
		my_ID = id3;
	} else  {
		utility_exit_with_message( "Problem with atoms " );
	}

	dof_id = DOF_ID( my_ID, D );
	//	std::cout << "Attempt to vary bond length for resno " << my_ID.rsd() << " atom: " << pose.residue( my_ID.rsd() ).atom_name( my_ID.atomno() ) << std::endl;

	mm.set( dof_id, true );

	mm.set( DOF_ID( my_ID, PHI) , true );
	mm.set( DOF_ID( my_ID, THETA ), true );

	cst_set->add_dof_constraint( dof_id, new HarmonicFunc( (atom2->xyz() - atom3->xyz() ).length() , 0.01 ) );

}


/////////////////////////////////////////////////////////////////////////////////
void vary_geometry_RNA( pose::Pose & pose, kinematics::MoveMap & mm )
{

	using namespace core::id;
	using namespace core::scoring;
	using namespace core::scoring::rna;
	using namespace core::scoring::constraints;

	ConstraintSetOP cst_set( 	pose.constraint_set()->clone() ) ;

	//Change this to also include D DOF's for sidechains.
	for (Size i = 1; i <= pose.total_residue(); i++ ) {
		for (Size n = 1; n <= NUM_RNA_MAINCHAIN_TORSIONS; n++ ) {
			TorsionID tor_id( i, id::BB, n );
			vary_bond_length( pose, tor_id, mm, cst_set );
		}
		for (Size n = 1; n <= NUM_RNA_CHI_TORSIONS; n++ ) {
			TorsionID tor_id( i, id::CHI, n );
			vary_bond_length( pose, tor_id, mm, cst_set );
		}
	}

	pose.constraint_set( cst_set );

}

////////////////////////////////Start of Function coded by Parin///////////////////////////////////////////////

///////////////////////////////////////////Records///////////////////////////////////////////////////////////

struct job_file_data {
	Size rotamer_group_min;
	Size rotamer_group_max;
	Real base_score;
	std::string tag;
	Size job_number;
//	utility::vector1< Size > rotamers_group;
//	utility::vector1< Size > rotamers;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////


struct pose_data_struct_dummy {
	Real rmsd_wrt_min;
	Real rmsd;
	Real score;
	pose::Pose pose;
	std::string rotamers_string_input;
	utility::vector1< Size > rotamer_group_count;
	utility::vector1< Size > rotamer_count;
};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct pose_data_struct {
	Real rmsd_cor;
	Real rmsd_min;
	Real rmsd;
	Real loop_rmsd_cor;
	Real loop_rmsd_min;
	Real loop_rmsd;
	Real diff_torsions;
	Real score;
	Real new_score;
	Real norm_new_score;
	Real dist_to_close_chain;
	pose::PoseOP pose_OP;
//	utility::vector1< Size > rotamer_group_count;
//	utility::vector1< Size > rotamer_count;
	std::string tag;
	Size pose_num; //order of pose according to score;

	Real loop_rmsd_exclude_bulge;
	Real loop_rmsd_cor_exclude_bulge;
};

struct output_data_struct {

	Real rmsd_wrt_correct_minimized;
	Real rmsd_wrt_minimized;
	Real rmsd;
	Real loop_rmsd_wrt_correct_minimized;
	Real loop_rmsd_wrt_minimized;
	Real loop_rmsd;
	Real diff_torsions;
	Real current_score;
//	std::string tag;
//	Real O3_C5_distance;
//	utility::vector1< Size > rotamer_group_count;
//	utility::vector1< Size > rotamer_count;

};


struct count_struct{
	Size output_pose_count;
	Size good_rep_rotamer_count;
	Size good_atr_rotamer_count;
	Size good_angle_count;
	Size good_distance_count;
	Size C5_O3_distance_count;
	Size Near_Native_Rotamer_Count;
	Size Near_Native_Pose_Count;
	Size both_count;
	Size tot_rotamer_count;
	Size native_contact_count;
	Size fine_rmsd_count;
	Size rmsd_count;
	Size native_minimize_pass_count;
	Size native_minimize_fail_count;
	Size pass_before_last_res_reb_chain_break_U;
	Size total_before_last_res_reb_chain_break_U;
	Size pass_before_last_res_reb_chain_break_M;
	Size total_before_last_res_reb_chain_break_M;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

enum node_state{ NOT_INIT, CLOSE, OPEN };

struct temp_pose_holder_struct{
	pose::Pose pose;
	node_state state;
};
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct pose_data_struct2 {
	Real rmsd;
	Real score;
	Size group_rotamer;
	Size subgroup_rotamer;
	pose::Pose pose;
	std::vector< Real> diff_torsions;
	std::string tag;
	Real O3_C5_distance;
};

///////////////////////////////////////////////////////////////////////////////

pose::Pose
attach_next_nucleotide(pose::Pose base_pose, pose::Pose& original_pose, Size& current_rebuild_num)
{

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace core::kinematics;
	using namespace core::optimization;
	using namespace core::io::silent;
	using namespace core::id;
	using namespace protocols::rna;


	Size const full_pose_reb_res=get_full_pose_reb_res(current_rebuild_num);
	Size const reb_res=Convert_to_partial_pose_seq_num(full_pose_reb_res, current_rebuild_num);


  bool const prepend_res = Is_prepend(current_rebuild_num);
 	ResidueTypeSetCAP rsd_set; //This line is repeated
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "rna" ); //This line is repeated

	pose::Pose pose = base_pose;

	Output_fold_tree_info(pose, "rebuilt_pose (before appending/prepending residue)");

	Size const nres( pose.total_residue() ); //Does not count the residue that will be built in this function call.

  std::cout << "nres= " << nres << std::endl;

	chemical::AA res_aa = aa_from_name(RunTimeParameters::rebuild_residue_list[current_rebuild_num].name);
	std::cout << "res_aa string: " << RunTimeParameters::rebuild_residue_list[current_rebuild_num].name;
	std::cout << "  res_aa: " << res_aa << std::endl;

	ResidueOP new_rsd = conformation::ResidueFactory::create_residue( *(rsd_set->aa_map( res_aa )[1]) ) ;


//  core::pose::remove_lower_terminus_type_from_pose_residue(pose, reb_res );	 //need this if prepend from end, before building residue
//	core::pose::add_lower_terminus_type_to_pose_residue(pose, reb_res+1); //need this if prepend from end, after building residue
//	core::pose::remove_upper_terminus_type_from_pose_residue(pose, reb_res-1 );  //need this if prepend from end, before building residue
//	core::pose::add_upper_terminus_type_to_pose_residue( pose, reb_res); //need this if prepend from end, after building residue
//  Need to update the pose's seq_number (using Convert_to_partial_pose_seq_num function) for the function above

	if(prepend_res) {

		dump_pdb( pose, "before_removing_Virtual_posphate_on_previous_rebuild_res.pdb");
		core::pose::remove_variant_type_from_pose_residue( pose, "VIRTUAL_PHOSPHATE", reb_res ); //This is bad coding practice ...I thought I fixed this!!!
		dump_pdb( pose, "after_removing_Virtual_posphate_on_previous_rebuild_res.pdb");

		pose.prepend_polymer_residue_before_seqpos( *new_rsd, Convert_to_partial_pose_seq_num(full_pose_reb_res+1, current_rebuild_num-1), true);

		dump_pdb( pose, "before_adding_Virtual_posphate_on_current_rebuild_res.pdb");
		core::pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_PHOSPHATE", reb_res );
		dump_pdb( pose, "after_adding_Virtual_posphate_on_current_rebuild_res.pdb");


	}else {

    pose.append_polymer_residue_after_seqpos( *new_rsd, Convert_to_partial_pose_seq_num(full_pose_reb_res-1, current_rebuild_num-1), true);

//    pose.energies().clear(); (if use pose.conformation() functions, need to clear energy...not sure of the reason
//    but the version of the function in pose is equivalent to energies().clear() and then calling the version of the
//    the version of the function in pose.conformation
//  	pose.conformation().safely_append_polymer_residue_after_seqpos(*new_rsd, reb_res-1, true);
//  	pose.conformation().append_polymer_residue_after_seqpos( *new_rsd, reb_res-1, true)
	}


	Output_fold_tree_info(pose, "rebuilt_pose (after appending/prepending residue)");

	return pose;
}


// Debug virtual phosphate
/*
		std::cout << "After removing phosphate reb res";
		std::cout << " Charge of P= " << pose.residue(reb_res).atomic_charge( 1 );
		std::cout << " Charge of O1P= " << pose.residue(reb_res).atomic_charge( 2 );
		std::cout << " Charge of O2P= " << pose.residue(reb_res).atomic_charge( 3 );
		std::cout << " Charge of O5'= " << pose.residue(reb_res).atomic_charge( 4 ) << std::endl;

		std::cout << "After removing phosphate res+1";
		std::cout << " Charge of P = " << pose.residue(reb_res+1).atomic_charge( 1 );
		std::cout << " Charge of O1P= " << pose.residue(reb_res+1).atomic_charge( 2 );
		std::cout << " Charge of O2P= " << pose.residue(reb_res+1).atomic_charge( 3 );
		std::cout << " Charge of O5'= " << pose.residue(reb_res+1).atomic_charge( 4 ) << std::endl;
*/


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::string
Previous_Res_PuckerState(pose::Pose& base_pose, Size const & current_rebuild_num){
		using namespace core::scoring;
		using namespace core::scoring::rna;
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		std::cout << "In Previous_residue_PuckerState" << std::endl;
		std::cout << "  current_rebuild_num= " << current_rebuild_num << std::endl;

		bool const prepend_res = 	Is_prepend(current_rebuild_num);

		if(prepend_res) {
				std::cout << "Previous_Res_PuckerState should not be called in prepend mode ERROR" << std::endl;
				exit (1);
		}

		Size full_pose_reb_res=get_full_pose_reb_res(current_rebuild_num);
		Size previous_res=Convert_to_partial_pose_seq_num(full_pose_reb_res-1, current_rebuild_num-1);		 //Really want (full_pose_reb_res-1 of basepose)

		std::cout << "   full_pose_previous_res_seq_num= " << (full_pose_reb_res-1);
		std::cout << "   base_pose_previous_res_seq_num= " << previous_res;



		scoring::rna::RNA_TorsionPotential const rna_torsion_potential;
		Real const DELTA_CUTOFF( rna_torsion_potential.delta_cutoff() );
		std::cout << "  DELTA_CUTOFF angle=" << DELTA_CUTOFF;

		conformation::Residue const & rsd(base_pose.residue(previous_res));
		Real const & delta( rsd.mainchain_torsion( DELTA ) );

		std::cout << "  delta angle=" << delta << std::endl;

		if (delta <= DELTA_CUTOFF) {
			return "NORTH";
		} else {
			return "SOUTH";
		}
}


////////////////////////////////////////////////////////////////////////////////////////////////
void
o2star_minimize(pose::Pose& pose, core::scoring::ScoreFunctionOP & packer_scorefxn){



	//TR << "Repacking 2'-OH ... " << std::endl;

	pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
	task->initialize_from_command_line();

	for (Size i = 1; i <= pose.total_residue(); i++) {
		if ( !pose.residue(i).is_RNA() ) continue;
		task->nonconst_residue_task(i).and_extrachi_cutoff( 0 );
		//		task->nonconst_residue_task(i).or_ex4( true );
		task->nonconst_residue_task(i).or_include_current( true );
		// How about bump check?
	}

//	TR << "Orienting 2' hydroxyls..." << std::endl;

//		pack::pack_rotamers( pose, *packer_scorefxn, task);

		pack::rotamer_trials( pose, *packer_scorefxn, task);


}
/////////////////////////////////////////////////////////////////////////////////////////////////////
//We will delete the unrebuild residues in the order in which the residues will be rebuild
//For the first residue to be deleted of each chain_break, we could either delete by "unprepend" or by "unappend"
//Previously, we made the choice that for 1 residue gap, we will choose prepend. This means that for the first residue to be deleted of each chain_break we will choose to delete by "unprepend"

pose::Pose
create_partially_rebuild_pose(pose::Pose const & original_pose, Size& current_rebuild_num){

		utility::vector1 <Residue_info_struct> missing_residue_list;

		pose::Pose partially_rebuild_pose=original_pose;

		for(Size i=RunTimeParameters::rebuild_residue_list.size(); i>current_rebuild_num; i--){

				Size delete_res_full_pose_seq_num=RunTimeParameters::rebuild_residue_list[i].seq_num;

				Size delete_res_rebuild_pose_seq_num = Convert_to_partial_pose_seq_num(delete_res_full_pose_seq_num, i);
//				std::cout << "before delete_polymer_residue, i= " << i << std::endl;
//				std::cout << " delete_res_full_pose_seq_num " << delete_res_full_pose_seq_num << std::endl;
//				std::cout << " delete_res_rebuild_pose_seq_num " << delete_res_rebuild_pose_seq_num << std::endl;

//				test_pose.conformation().delete_residue_slow(i);  //Will need to use this if prepend residue to the first residue in pose (a lower_terminus lower)

				partially_rebuild_pose.delete_polymer_residue(delete_res_rebuild_pose_seq_num);

				//Would have had to include the special case for the first residue to be deleted of each chain_break if we had allowed first residue deleted of each chain_break to be "un-append" as well.
				if(RunTimeParameters::append_prepend_choice_list[i]==true){ //if the delete residue was deleted by "unprepend", then delete phosphate group of the loop's new three_prime_end residue.
					core::pose::add_variant_type_to_pose_residue( partially_rebuild_pose, "VIRTUAL_PHOSPHATE", Convert_to_partial_pose_seq_num(delete_res_full_pose_seq_num+1, i-1));
				}
				//In the future, might need virtual atom variant type for append case as well

		}

		return partially_rebuild_pose;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector< pose::Pose>
create_partially_rebuild_pose_list(pose::Pose const & original_pose, std::string output_pose_name){

	//0-> total unrebuild
	//num_residue_reb -> total rebuild

	 using namespace core::io::silent;
   using namespace core::chemical;
	 using namespace basic::options;
	 using namespace basic::options::OptionKeys;
   using namespace core::conformation;

	 std::cout << "enter create_partially_rebuild_pose_list function" << std::endl;

	 Size num_residue_reb= RunTimeParameters::rebuild_residue_list.size();
	 std::vector< pose::Pose> rebuild_pose_list;


	 for(Size current_rebuild_num=0; current_rebuild_num <= num_residue_reb; current_rebuild_num++){
//			std::cout << "current_rebuild_num= " << current_rebuild_num << std::endl;
			rebuild_pose_list.push_back(create_partially_rebuild_pose(original_pose, current_rebuild_num));
	 }

	 Size const node_number = option[node_number_input ] ;

   std::string name;
	 for(Size n=0; n<=num_residue_reb; n++) {
			name.clear();
			name.append(output_pose_name);
			name.append(lead_zero_string_of(n, 1));
			name.append("_");
			name.append(lead_zero_string_of(node_number, 3));
			dump_pdb( rebuild_pose_list[n], name+".pdb");
	 }

	 return rebuild_pose_list;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
void
UnFreeze_sugar_torsions(core::kinematics::MoveMap& mm, pose::Pose& pose){

	 using namespace core::id;
	 using namespace basic::options;
	 using namespace basic::options::OptionKeys;

	 std::cout << "UnFreeze pose sugar torsions" << std::endl;

	 Size const nres( pose.total_residue() );

	 for(Size i=1; i<=nres; i++){

			mm.set( TorsionID( i  , id::BB,  4 ), true ); //delta_i
			mm.set( TorsionID( i  , id::CHI, 2 ), true ); //nu2_i
			mm.set( TorsionID( i  , id::CHI, 3 ), true );	//nu1_i

	 }
}

*/

void
Freeze_sugar_torsions(core::kinematics::MoveMap& mm, pose::Pose& pose){

	 using namespace core::id;
	 using namespace basic::options;
	 using namespace basic::options::OptionKeys;

	 std::cout << "Freeze pose sugar torsions" << std::endl;

	 Size const nres( pose.total_residue() );

	 for(Size i=1; i<=nres; i++){

			mm.set( TorsionID( i  , id::BB,  4 ), false ); //delta_i
			mm.set( TorsionID( i  , id::CHI, 2 ), false ); //nu2_i
			mm.set( TorsionID( i  , id::CHI, 3 ), false );	//nu1_i

	 }
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



bool
Contains_residue(utility::vector1 <Residue_info_struct>& residue_list, Size residue_seq_num){

	for(Size i=1; i<=residue_list.size(); i++){

			if(residue_list[i].seq_num==residue_seq_num){
				return true;
			}
	}

	return false;
}



//////////////////////////////////////////////////////////////////////////////////////////////////////
void
Freeze_unrebuild_residues(core::kinematics::MoveMap& mm, Size current_rebuild_num, pose::Pose& current_pose){

		using namespace basic::options;
	  using namespace basic::options::OptionKeys;

		//using namespace RunTimeParameters;

		std::cout << "In Freeze_unrebuild_residues function " << std::endl;

		mm.set_bb(false);
	  mm.set_chi(false);



		utility::vector1 <Residue_info_struct> unfreeze_residue_list;
		utility::vector1 <Residue_info_struct> already_rebuild_residue_list;

		for(Size i=1; i<=current_rebuild_num; i++){
				already_rebuild_residue_list.push_back(RunTimeParameters::rebuild_residue_list[i]);
		}

		unfreeze_residue_list=already_rebuild_residue_list;



		//Unfreeze residues at the edge of the chain break as well. (only after at least one residue is rebuild from that edge)

		//June 13,2009
		//Reminder: In old mode, didn't unfreeze the residue before the first residue rebuild in the full_rebuild_mode ( did however implement this in the build_one_res_mode
		//The reason why you want to unfreeze the residue before the first residue rebuild is:
		//For prepend mode: The rebuild torsion alpha, beta and gamma is a part of the i+1 residue
		//For append mode: The rebuild torsion gamma and zeta (control the phospate group) is a part of the i-1 residue.

		utility::vector1 <Residue_info_struct> current_pose_residue_list;
		get_partial_pose_residue_list(current_pose_residue_list, current_rebuild_num);


		for(Size i=1; i<=current_pose_residue_list.size(); i++){

			//Check that not already in the list.
			if(Contains_residue(already_rebuild_residue_list, current_pose_residue_list[i].seq_num)==false){

				//Check if residue is at the edge of the loop.
				if(Contains_residue(already_rebuild_residue_list, current_pose_residue_list[i].seq_num+1) || Contains_residue(already_rebuild_residue_list, current_pose_residue_list[i].seq_num-1)){
						unfreeze_residue_list.push_back(current_pose_residue_list[i]);
				}

			}
		}
		std::cout << "current_rebuild_num= " << current_rebuild_num << std::endl;


		std::cout <<  std::setw(40)  << "rebuild_residue_list: "; Output_residue_list(RunTimeParameters::rebuild_residue_list);
		std::cout <<  std::setw(40)  << "already_rebuild_residue_list: "; Output_residue_list(already_rebuild_residue_list);
		std::cout <<  std::setw(40)  << "current_pose_residue_list: "; Output_residue_list(current_pose_residue_list);
		std::cout <<  std::setw(40)  << "unfreeze_residue_list: "; Output_residue_list(unfreeze_residue_list);

		//Apply mm.
		for(Size i=1; i<=unfreeze_residue_list.size(); i++){

				Size current_pose_seq_num = Convert_to_partial_pose_seq_num(unfreeze_residue_list[i].seq_num, current_rebuild_num);
				mm.set_bb(current_pose_seq_num, true);
				mm.set_chi(current_pose_seq_num, true );
		}


}

///////////////////////////////////////////////////////////////////////////////


std::string
Print_node_state(node_state& state){
	std::string string_state;
	switch(state) {
		case NOT_INIT:
			string_state="Not intialized";
			break;
		case OPEN:
			string_state="Open";
			break;
		case CLOSE:
			string_state="Close";
			break;
		default:
			string_state="Error";
			break;
	}
	return string_state;
}

///////////////////////////////////////////////////////////////////////////////
void
update_branch_position(utility::vector1< Size >& rotamer_group_count, utility::vector1< Size >& rotamer_count, utility::vector1< vector1 < node_state> >& rotamer_branch_state, Size rotamer_group_size, Size total_rotamer_group, Size num_residue_reb) {


  	for(Size n=num_residue_reb; n>=1; n--) {

		   if(rotamer_count[n]>rotamer_group_size) {
		   std::cout << "all rotamer in a rotamer group enumerated, move to the next group :";

		 	//After enumerating all rotamer in a rotamer group, move to the next group
		     rotamer_count[n]=1;
		     rotamer_group_count[n]=rotamer_group_count[n]+1;

		     rotamer_branch_state[n].assign(rotamer_group_size, NOT_INIT);
		   	 //The branch state of the nth residue is cleared before moving to the next branch on the nth residue level
		   	}
		   if(rotamer_group_count[n]>total_rotamer_group && n!=1) {

		  //After enumerating every groups, move to adjacent branch
       std::cout << "Every group enumerated, move to adjacent branch :";
			   rotamer_group_count[n]=1;
		   	 rotamer_count[n-1]=rotamer_count[n-1]+1;
		  	}
		}
}


////////////////////////////////////////////////////////////////////////////////////////////
Size
count_neighbor(utility::vector1< pose::Pose>& group_pose, vector1< vector1 <Real> >& rmsd_matrix, utility::vector1< node_state >& rotamer_group_state, Real neighbor_rsmd_cutoff, Size current_row) {
	Size neighbor_num=0;
//	std::cout <<" The current row is " << current_row << std::endl;
//	std::cout <<" neighbors_rsmd_cutoff " << neighbor_rsmd_cutoff << std::endl;
//	std::cout << "rsmd matrix of row " << current_row << std::endl;

	for(Size col =1; col<=group_pose.size(); col++) {
//		if (rotamer_group_state[current_row]==NOT_INIT && rotamer_group_state[col]==NOT_INIT){
			if(rmsd_matrix[current_row][col]<neighbor_rsmd_cutoff){
				neighbor_num=neighbor_num + 1;
				std::cout << "neighbor found!"<< std::endl;
			}
//		}
	}
	return neighbor_num;
}

void
create_filename_(std::string& filename, std::string name, Size job_number, Size current_rebuild_num, Size digits){

				std::string foldername="residueNum_";
				if(current_rebuild_num<10) {
					foldername.append(lead_zero_string_of(current_rebuild_num, 1));
				} else {
					foldername.append(lead_zero_string_of(current_rebuild_num, 2));
				}

				foldername.append("_");
				foldername.append(Get_one_letter_name(RunTimeParameters::rebuild_residue_list[current_rebuild_num].name));
		   	Size full_pose_res_reb = RunTimeParameters::rebuild_residue_list[current_rebuild_num].seq_num;

				if(full_pose_res_reb<10) {
					foldername.append(lead_zero_string_of(full_pose_res_reb, 1));
				} else {
					foldername.append(lead_zero_string_of(full_pose_res_reb, 2));
				}


				filename=foldername;
				filename.append("/");
				filename.append(name);
				filename.append("_");
				filename.append(lead_zero_string_of(job_number, digits)); //Allow for possibility that there will be more than 1000 total_jobs/num_branch
				filename.append(".txt");
}

void
create_filename(std::string& filename, std::string name, Size job_number, Size current_rebuild_num){

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if(option[old_file_format]==true){
		create_filename_(filename, name, job_number, current_rebuild_num, 3);
	} else {
		create_filename_(filename, name, job_number, current_rebuild_num, 4);
	}

}






void
create_communication_filename(std::string& filename, std::string name, Size const communication_num, Size node_number, Size current_rebuild_num){
		create_filename(filename, name, node_number, current_rebuild_num); //Overload create_filename
		filename.erase (filename.length()-4, 4); //Delete .txt
		filename.append("_"+lead_zero_string_of(communication_num, 3));
		filename.append(".txt");
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Tokenize(std::string& str, std::vector<std::string>& tokens)
{
	  using namespace std;
	  std::string delimiters= " \t\n\f\v";

    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
Size
get_column_pos(std::vector<std::string>& name_tokens, std::string col_name) {

		for( Size i = 0; i < name_tokens.size(); i++) {
			std::cout << "Column Name Token " << i << ": " << name_tokens[i] << "col_name :" << col_name <<std::endl;
			if (name_tokens[i]==col_name) {
			return i;
			std::cout << "Found match at column " << i << std::endl;
			}
		}
	return 99;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
void
create_column_pos_map(std::map< std::string, Size >& column_pos_map, std::ifstream& infile){
	using namespace std;

	std::vector<std::string> name_tokens;
	std::string name_string;
	getline(infile, name_string);
  Tokenize(name_string, name_tokens);

	for( Size i = 0; i < name_tokens.size(); i++) {
			std::cout << "Column Name Token " << i << ": " << name_tokens[i] << std::endl;
	}

//	column_pos_map.insert ( std::pair< string, Size>("rmsd", get_column_pos(name_tokens, "rmsd") ));
//	std::cout << "rmsd is located at: " << column_pos_map["rmsd"]<< std::endl;
//	column_pos_map.insert ( std::pair< string, Size>("rmsd_wrt_min", get_column_pos(name_tokens, "rmsd_wrt_min") ));
//	std::cout << "rmsd_wrt_min is located at: " << column_pos_map["rmsd_wrt_min"]<< std::endl;
	column_pos_map.insert ( std::pair< string, Size>("score", 4 );
	std::cout << "score is located at: " << column_pos_map["score"]<< std::endl;
	column_pos_map.insert ( std::pair< string, Size>("description", 5 );
	std::cout <<"description is located at: " << column_pos_map["description"]<< std::endl;
}
*/
//////////////////////////////////////////////////////////////////////////////////


void
cut_bulge_test(){

			using namespace basic::options;
			using namespace basic::options::OptionKeys;

			std::string tag = option[ in ::file::s ][1]; //S_000001.pdb

			std::ifstream infile;
		 	infile.open(tag.c_str());
			std::string data_string;



			std::string filename="Decoy_"+tag.substr(2,10);
			std::ofstream outfile;
	  	outfile.open(filename.c_str());


			while(getline(infile, data_string)) {

					std::vector<std::string> data_tokens;
					Tokenize(data_string, data_tokens);

					Size residueNum=0;
					if(data_tokens.size()>=5){
					residueNum= convert_string_to_int(data_tokens[5]);
					}

					if(residueNum==4 || residueNum==5 || residueNum==6){

					} else {
							outfile << data_string << "\n";
	//						std::cout << data_string << std::endl;
					}

//					for(Size i=0; i<data_tokens.size(); i++){
//						outfile << data_tokens[i] << "	";
//						std::cout << data_tokens[i] << "	";

//					}
//				outfile << "\n";


			}
		infile.clear();
  	infile.close();
		outfile.close();
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
convert_star_to_dash_test(){

			using namespace basic::options;
			using namespace basic::options::OptionKeys;

			std::string tag = option[ in ::file::s ][1]; //S_000001.pdb

			std::ifstream infile;
		 	infile.open(tag.c_str());
			std::string data_string;



			std::string filename="Q_"+tag.substr(2,12);
			std::ofstream outfile;
	  	outfile.open(filename.c_str());


			while(getline(infile, data_string)) {

					std::vector<std::string> data_tokens;
					Tokenize(data_string, data_tokens);

					for(int i=0; i<data_string.size(); i++){
						if(data_string[i]=='*'){
							data_string[i]='\'';
						}
					}

						outfile << data_string << "\n";


//					for(Size i=0; i<data_tokens.size(); i++){
//						outfile << data_tokens[i] << "	";
//						std::cout << data_tokens[i] << "	";

//					}
//				outfile << "\n";


			}
		infile.clear();
  	infile.close();
		outfile.close();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
split_pdb_file_test(){

			using namespace basic::options;
			using namespace basic::options::OptionKeys;

			std::string tag = option[ in ::file::s ][1]; // O1_000001.pdb or O2_000001.pdb

			std::ifstream infile;
		 	infile.open(tag.c_str());
			std::string data_string;

			std::string mode;
			if(tag[1]=='1'){
				mode="1";
				std::cout << "mode: " << mode << std::endl;
			} else if(tag[1]=='2'){
				mode="2";
				std::cout << "mode: " << mode << std::endl;
			} else {
				std::cout << "Error, no mode selected" << std::endl;
//				exit (1);
			}

			std::string filename;

			if(mode=="1"){
				filename="q_A_"+tag.substr(3,11);
			} else {
				filename="q_Z_"+tag.substr(3,11);
			}

			std::ofstream outfile;
			outfile.open(filename.c_str());


			while(getline(infile, data_string)) {

					std::vector<std::string> data_tokens;
					Tokenize(data_string, data_tokens);
					if(data_tokens[0]=="ENDMDL") break;

					outfile << data_string << "\n";
	//				std::cout << data_string << std::endl;

			}

			std::string filename2;
			if(mode=="1"){
				filename2="q_B_"+tag.substr(3,11);
			} else {
				filename2="q_C_"+tag.substr(3,11);
			}

			std::ofstream outfile2;
	  	outfile2.open(filename2.c_str());
			getline(infile, data_string);//Get thie MODEL1 line

			while(getline(infile, data_string)) {

					std::vector<std::string> data_tokens;
					Tokenize(data_string, data_tokens);

					outfile2 << data_string << "\n";
//				std::cout << data_string << std::endl;

			}


		infile.clear();
  	infile.close();
		outfile2.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
final_output_test(){

			using namespace basic::options;
			using namespace basic::options::OptionKeys;

			std::string pose_name = option[ in ::file::s ][1];

			std::string tag="blah.out";
			std::ifstream infile;
		 	infile.open(tag.c_str());
			std::string data_string;

			std::cout << "pose_name= " << pose_name << std::endl;

			std::string filename="final_data.txt";
			std::ofstream outfile;
	  	outfile.open(filename.c_str(), std::ios::app);


			while(getline(infile, data_string)) {


					std::vector<std::string> data_tokens;
					Tokenize(data_string, data_tokens);

					if(data_tokens.size()>=19){
						std::cout << "data_tokens[19]: " << data_tokens[19] << std::endl;

						if(data_tokens[19]==pose_name){
								infile.clear();
  							infile.close();
								break;
						}
//				std::cout << data_string << std::endl;
				}
			}

			///////////////////////////////////
			std::string filename2="temp.txt";
			std::ifstream infile2;
		 	infile2.open(filename2.c_str());
			std::string data_string2;
			getline(infile2, data_string2);
			infile2.clear();
  		infile2.close();
			///////////////////////////////////
			std::vector<std::string> data_tokens2;
			Tokenize(data_string2, data_tokens2);
			if(data_tokens2.size()!=0){
				for(int i=0; i< data_tokens2.size();i++){

				Real angle=convert_string_to_real(data_tokens2[i]);
				if(angle<-180) angle=angle+180;
				if(angle>180) angle=angle-180;

				outfile << angle << "		" ;
				std::cout << "data_tokens2[i]: " << angle;

				}
				std::cout << std::endl;
				outfile << data_string << "\n";
			}

//			std::cout <<"data_string2 " << data_string2 << std::endl;

//			outfile << data_string2;


		outfile.close();

}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
import_pose_data(utility::vector1<pose_data_struct>& all_poses_data, std::ifstream& infile, std::map <std::string, Size>& column_pos_map, char pose_type) {

/*
		  using namespace core::chemical;
			using namespace core::conformation;
			using namespace core::scoring;
			using namespace core::scoring::constraints;
			using namespace core::kinematics;
			using namespace core::optimization;
			using namespace core::io::silent;
			using namespace core::id;
			using namespace protocols::rna;

	  ResidueTypeSetCAP rsd_set; //This line is repeated
		rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "rna" ); //This line is repeated
*/
		std::string data_string;

		while(getline(infile, data_string)) {
//		for(Size blah=1; blah<=4; blah++){
//			getline(infile, data_string);

			pose_data_struct pose_data;
			std::vector<std::string> data_tokens;


			Tokenize(data_string, data_tokens);

			if(data_tokens[column_pos_map["description"]][0]==pose_type){
				std::cout << " " << data_tokens[column_pos_map["description"]] << std::endl;
				std::cout << " rsmd= " <<data_tokens[column_pos_map["rmsd"]];
				std::cout << " min= " <<data_tokens[column_pos_map["rmsd_min"]];
				std::cout << " cor= " <<data_tokens[column_pos_map["rmsd_cor"]];
				std::cout << " loop_rsmd= " <<data_tokens[column_pos_map["loop_rmsd"]];
				std::cout << " loop_min= " <<data_tokens[column_pos_map["loop_rmsd_min"]];
				std::cout << " loop_cor= " <<data_tokens[column_pos_map["loop_rmsd_cor"]];
				std::cout << " score= " <<data_tokens[column_pos_map["score"]];
				std::cout << " diff_tors= " <<data_tokens[column_pos_map["diff_torsions"]];
				std::cout << " dist_to_close_chain= " <<data_tokens[column_pos_map["dist_to_close_chain"]];
				std::cout<< std::endl;



//////////Convert the rsmd_wrt_min and score string into reals and store them in pose_data structure/////////
			  pose_data.tag = data_tokens[column_pos_map["description"]];
			  pose_data.rmsd=convert_string_to_real(data_tokens[column_pos_map["rmsd"]]);
			  pose_data.rmsd_min=convert_string_to_real(data_tokens[column_pos_map["rmsd_min"]]);
			  pose_data.rmsd_cor=convert_string_to_real(data_tokens[column_pos_map["rmsd_cor"]]);
			  pose_data.loop_rmsd=convert_string_to_real(data_tokens[column_pos_map["loop_rmsd"]]);
			  pose_data.loop_rmsd_min=convert_string_to_real(data_tokens[column_pos_map["loop_rmsd_min"]]);
			  pose_data.loop_rmsd_cor=convert_string_to_real(data_tokens[column_pos_map["loop_rmsd_cor"]]);
				pose_data.score=convert_string_to_real(data_tokens[column_pos_map["score"]]);
				pose_data.diff_torsions=convert_string_to_real(data_tokens[column_pos_map["diff_torsions"]]);
				pose_data.dist_to_close_chain=convert_string_to_real(data_tokens[column_pos_map["dist_to_close_chain"]]);

/*
 				pose::Pose pose;
		  	std::cout << "Import the following pose :" << pose_data.tag+".pdb" << std::endl;
				core::import_pose::pose_from_pdb( pose, *rsd_set, pose_data.tag+".pdb");
				pose_data.pose=pose;
*/

				all_poses_data.push_back( pose_data);
				data_string.clear();

			}
		} //End of getline() whileloop.

		std::cout << "# of pose_data imported into all_poses_data so far= " << all_poses_data.size() << std::endl;
}




////////The suspend time is a rough approximation/////////////////////////////////////////////////////////////////
void
suspend(Size seconds) {

		for(Size i=0; i<=25000; i++) {
			for(Size j=0; j<=seconds; j++) {
				pose_data_struct_dummy  dummy_pose_data;
				dummy_pose_data.score = 99999999.99;
			}
		}
}

//////////////////////////////////////////////////////////////////////////////////////////////
struct	job_info_struct{
		Size communication_num;
		Size node_num;
		Size job_number;
		bool is_done;
};

bool
Is_node_waiting_for_new_job(utility::vector1 <Size>& communication_num_list, Size const node_num, Size const current_rebuild_num){

		std::string status_filename;
    create_communication_filename(status_filename, "status", communication_num_list[node_num], node_num, current_rebuild_num);

		std::ifstream status_infile;
 	  status_infile.open(status_filename.c_str());

  	if (status_infile.fail()) {
	  	std::cout << "Could no open file named \"" << status_filename << "\".";
			std::cout << " node_num= " << std::setw(6) << std::left << node_num;
			std::cout <<  " communication_num= " << communication_num_list[node_num] << std::endl;
			suspend(1);
	  	status_infile.clear();
	  	status_infile.close(); //do I need this line?
	  return false;
		} else {
	  	status_infile.clear();
	  	status_infile.close(); //do I need this line?
		return true;
		}
}

bool
Is_job_finish(job_info_struct const & job_info, Size current_rebuild_num){
  using namespace std;

    std::string status_filename;

    create_communication_filename(status_filename, "status", job_info.communication_num, job_info.node_num, current_rebuild_num);

		std::ifstream status_infile;
 	  status_infile.open(status_filename.c_str());

  	if (status_infile.fail()) {
	  	std::cout << "Could no open file named \"" << status_filename << "\". Job is not finished yet" << std::endl;
	  	status_infile.clear();
	  	status_infile.close(); //do I need this line?
	  return false;
		} else {
	  	status_infile.clear();
	  	status_infile.close(); //do I need this line?
		return true;
		}

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
bool
Is_job_ready(std::string& filename, std::ifstream& infile){
  using namespace std;

 	  infile.open(filename.c_str());

  	if (infile.fail()) {
  	std::cout << "Could no open file named \"" << filename << "\". central_evaluation still running" << std::endl;
  	infile.clear();
  	infile.close(); //do I need this line?
	  return false;
		} else {
		return true; //This relies on the fact that job_file.txt have not been created yet.....consider changing.
		}
}
*/
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
sort_citeria(pose_data_struct pose_data_1, pose_data_struct pose_data_2){


	//Sort by score.
	return (pose_data_1.score < pose_data_2.score);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
sort_pose_data(utility::vector1< pose_data_struct >&	all_poses_data) {

		//Need to check if this work with vector1, if not switch to std::vector
		sort(all_poses_data.begin(), all_poses_data.end(), sort_citeria);

		for(Size i=1; i<= all_poses_data.size(); i++){
			all_poses_data[i].pose_num=i;
		}

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
Print_heavy_atoms(Size suite_num, pose::Pose const & pose1, pose::Pose const & pose2){

		using namespace conformation;
		using namespace basic::options;
		using namespace id;

		Size num_atoms;

		num_atoms=std::max(pose1.residue(suite_num).nheavyatoms(), pose2.residue(suite_num).nheavyatoms());

		std::cout << "num_atoms: " << num_atoms << std::endl;

		for(Size n=1; n<= num_atoms;  n++){

			std::cout << " atom num = " <<  n;
			std::cout << "  atom_name of the pose1 " <<  pose1.residue(suite_num).atom_name(n);
			std::cout << "  atom_name of the pose2 " <<  pose2.residue(suite_num).atom_name(n) << std::endl;

		}
}///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
suite_square_deviation(pose::Pose const & pose1, pose::Pose const & pose2, bool const & prepend_res, Size const & suite_num, Size& atom_count, Real& sum_sd, bool verbose){

  		chemical::AA res_aa =  pose1.residue(suite_num).aa();
  		chemical::AA res_aa2 =  pose2.residue(suite_num).aa();

  		Size first_sidechain_atom1=pose1.residue(suite_num).first_sidechain_atom();
  		Size first_sidechain_atom2=pose2.residue(suite_num).first_sidechain_atom();


			Size num_side_chain_atom;

	 		Size num_side_chain_atom_old_1=pose1.residue(suite_num).nheavyatoms()-(first_sidechain_atom1-1);
	 		Size num_side_chain_atom_old_2=pose2.residue(suite_num).nheavyatoms()-(first_sidechain_atom2-1);

			if(name_from_aa(res_aa)=="RAD") {
				if(verbose) std::cout << "name_from_aa: RAD" << std::endl;
				num_side_chain_atom=11;
			} else if(name_from_aa(res_aa)=="RCY") {
				if(verbose) std::cout << "name_from_aa: RCY" << std::endl;
				num_side_chain_atom=9;
			} else if(name_from_aa(res_aa)=="RGU") {
				if(verbose) std::cout << "name_from_aa: RGU" << std::endl;
				num_side_chain_atom=12;
			} else if(name_from_aa(res_aa)=="URA") {
				if(verbose) std::cout << "name_from_aa: URA" << std::endl;
				num_side_chain_atom=9;
			} else {
				std::cout << "Error, cannot identify residue type" << std::endl;
				num_side_chain_atom=0;
				exit (1);
		  }

//			Size first_sidechain_hydrogen1=pose1.residue(suite_num).first_sidechain_hydrogen();
//  		Size first_sidechain_hydrogen2=pose2.residue(suite_num).first_sidechain_hydrogen();


 			if(verbose){
  			std::cout << " residue type1= " <<  res_aa << " residue type2= " <<  res_aa2;
  			std::cout << " 1st_side_atom1= " <<  first_sidechain_atom1;
  			std::cout << " 1st_side_atom2= " <<  first_sidechain_atom2;
  			std::cout << " nheavyatoms1= " <<  pose1.residue(suite_num).nheavyatoms();
  			std::cout << " nheavyatoms2= " <<  pose2.residue(suite_num).nheavyatoms() << std::endl;
				std::cout << " num_side_chain_atom = " <<  num_side_chain_atom << std::endl;
  			std::cout << " num_side_chain_atom_old_1 = " <<  num_side_chain_atom_old_1 << std::endl;
  			std::cout << " num_side_chain_atom_old_2 = " <<  num_side_chain_atom_old_2 << std::endl;
				Print_heavy_atoms(suite_num,pose1, pose2);
			}

			for(Size n=1; n<= 11; n++){//RNA contain 11 heavy backbone atoms.
  			Size res_num;
  			if(prepend_res){
					if(n < 5) {
						res_num=suite_num+1;
					}else {
						res_num=suite_num;
					}
				} else {
					res_num=suite_num;
				}

				atom_count++;

  			Distance dist_squared = (pose1.residue(res_num).xyz( n ) - pose2.residue(res_num).xyz( n ) ).length_squared();


				sum_sd=sum_sd+dist_squared;

				if(verbose){
					std::cout << "atom_name of the pose1= " << pose1.residue(res_num).atom_name(n);
  				std::cout << " atom_name of the pose2= " << pose2.residue(res_num).atom_name(n);
					std::cout << " Backbone atom= " <<  n << " dist_squared= " << dist_squared << std::endl;
				}

  		}

			Size res_num;
			if(prepend_res){
				res_num=suite_num+1;
			}else{
				res_num=suite_num;
			}

			if(verbose){
				Distance dist_squared = (pose1.residue(res_num).xyz( 2 ) - pose2.residue(res_num).xyz( 3 ) ).length_squared();
				std::cout << "atom_name of the pose1= " << pose1.residue(res_num).atom_name(2);
  			std::cout << " atom_name of the pose2= " << pose2.residue(res_num).atom_name(3);
				std::cout << " Switch Phosphate1= " << " dist_squared= " << dist_squared << std::endl;
				dist_squared = (pose1.residue(res_num).xyz( 3 ) - pose2.residue(res_num).xyz( 2 ) ).length_squared();
				std::cout << "atom_name of the pose1= " << pose1.residue(res_num).atom_name(3);
  			std::cout << " atom_name of the pose2= " << pose2.residue(res_num).atom_name(2);
				std::cout << " Switch Phosphate2= " << " dist_squared= " << dist_squared << std::endl;
			}


  		//Need to use num_side_chain_atom from pose1 since a silly bug in Rosetta miscalculate num_heavy_atom by considering the virtaul O2star hydrogen to be heavy_atom when it is set to virtual in the current_pose_screen
  		for(Size n=0; n<= num_side_chain_atom-1; n++){ //Sidechain atoms include O2star
  			atom_count++;

  			Distance const dist_squared = (pose1.residue(suite_num).xyz( n+first_sidechain_atom1) - pose2.residue(suite_num).xyz( n+first_sidechain_atom2)).length_squared();

				sum_sd=sum_sd+dist_squared;

				if(verbose){
			  	std::cout << "atom_name of the atom1= " << pose1.residue(suite_num).atom_name(n+first_sidechain_atom1);
 		 			std::cout << " atom_name of the atom2= " << pose2.residue(suite_num).atom_name(n+first_sidechain_atom2);
					std::cout << "Side chain atom= " <<  n << " dist_squared= " << dist_squared << std::endl;
				}
 		 	}

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



//Real
//loop_rmsd(pose::Pose const)

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Sometime the last residue rebuild is a bulge residue...use this function when don't want to include bulge residue in calculation
Real
Loop_rmsd_exclude_last_residue_rebuild(pose::Pose const & pose1, pose::Pose const & pose2, Size& current_rebuild_num, bool verbose){
			using namespace scoring::constraints;
		using namespace conformation;
		using namespace basic::options;
		using namespace id;


		if(verbose){
			utility::vector1 <Residue_info_struct> already_rebuild_residue_list;
			for(Size i=1; i<=current_rebuild_num; i++){
					already_rebuild_residue_list.push_back(RunTimeParameters::rebuild_residue_list[i]);
			}
			utility::vector1 <Residue_info_struct> loop_rmsd_residue_list;
			for(Size i=1; i<=current_rebuild_num-1; i++){
					loop_rmsd_residue_list.push_back(RunTimeParameters::rebuild_residue_list[i]);
			}
			std::cout << "In loop rmsd_exclude_last_residue_rebuild function";
			std::cout << "  Number of residues build so far=" << current_rebuild_num << std::endl;
			std::cout <<  std::setw(40)  << "rebuild_residue_list: "; Output_residue_list(RunTimeParameters::rebuild_residue_list);
			std::cout <<  std::setw(40)  << "already_rebuild_residue_list: "; Output_residue_list(already_rebuild_residue_list);
			std::cout <<  std::setw(40)  << "loop_rmsd_residue_list: "; Output_residue_list(loop_rmsd_residue_list);
		}

  	Size atom_count=0;
  	Real sum_sd=0;

  	for(Size i=1; i<=current_rebuild_num-1; i++){ //current_rebuild_num-1 since exclude last residue rebuild

			Size suite_num = Convert_to_partial_pose_seq_num(RunTimeParameters::rebuild_residue_list[i].seq_num, current_rebuild_num);
			bool prepend_res = 	Is_prepend(i);

			if(verbose){
					std::cout << "Full_pose_seq_num= " << RunTimeParameters::rebuild_residue_list[i].seq_num << std::endl;
					std::cout << "Rebuild_pose_seq_num= " << suite_num << std::endl;
			}

			//add atom in the suites to atom_count
			//add sd of each atom to sum_sd
			suite_square_deviation(pose1, pose2, prepend_res, suite_num, atom_count, sum_sd, verbose);

		}

  	sum_sd=sum_sd/(atom_count);
  	Real rmsd=sqrt(sum_sd);

		if(verbose){
	 		std::cout << "sum_sd= " << sum_sd << std::endl;
 			std::cout << " atom_count= " << atom_count << std::endl;
 			std::cout << " rmsd= " << rmsd << std::endl;
 		}

  	return (std::max(0.01, rmsd));
//		std::cout << "Current distance: " << dist << std::endl << std::endl;
}




//Need to update this function!!!
Real
Loop_rmsd(pose::Pose const & pose1, pose::Pose const & pose2, Size& current_rebuild_num, bool verbose) {
		using namespace scoring::constraints;
		using namespace conformation;
		using namespace basic::options;
		using namespace id;


		if(verbose){
			utility::vector1 <Residue_info_struct> already_rebuild_residue_list;
			for(Size i=1; i<=current_rebuild_num; i++){
					already_rebuild_residue_list.push_back(RunTimeParameters::rebuild_residue_list[i]);
			}
			std::cout << "In loop rmsd function";
			std::cout << "  Number of residues build so far=" << current_rebuild_num << std::endl;
			std::cout <<  std::setw(40)  << "rebuild_residue_list: "; Output_residue_list(RunTimeParameters::rebuild_residue_list);
			std::cout <<  std::setw(40)  << "already_rebuild_residue_list: "; Output_residue_list(already_rebuild_residue_list);
		}

  	Size atom_count=0;
  	Real sum_sd=0;

  	for(Size i=1; i<=current_rebuild_num; i++){

			Size suite_num = Convert_to_partial_pose_seq_num(RunTimeParameters::rebuild_residue_list[i].seq_num, current_rebuild_num);
			bool prepend_res = 	Is_prepend(i);

			if(verbose){
					std::cout << "Full_pose_seq_num= " << RunTimeParameters::rebuild_residue_list[i].seq_num << std::endl;
					std::cout << "Rebuild_pose_seq_num= " << suite_num << std::endl;
			}

			//add atom in the suites to atom_count
			//add sd of each atom to sum_sd
			suite_square_deviation(pose1, pose2, prepend_res, suite_num, atom_count, sum_sd, verbose);

		}

  	sum_sd=sum_sd/(atom_count);
  	Real rmsd=sqrt(sum_sd);

		if(verbose){
	 		std::cout << "sum_sd= " << sum_sd << std::endl;
 			std::cout << " atom_count= " << atom_count << std::endl;
 			std::cout << " rmsd= " << rmsd << std::endl;
 		}

  	return (std::max(0.01, rmsd));
//		std::cout << "Current distance: " << dist << std::endl << std::endl;
}

//Size
//Get_rebuild_num_from_full_pose_seq_num(){
//}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//This now works for any rebuild residue.
Real
suite_rmsd_any_seq_num(pose::Pose const &  pose1,pose::Pose const &  pose2, Size const & seq_num, bool const prepend_res){

		using namespace scoring::constraints;
		using namespace conformation;
		using namespace basic::options;
		using namespace id;

		Size atom_count=0;
  	Real sum_sd=0;

		suite_square_deviation(pose1, pose2, prepend_res, seq_num, atom_count, sum_sd, false);

		sum_sd=sum_sd/(atom_count);
  	Real rmsd=sqrt(sum_sd);

		return (std::max(0.01, rmsd));

}

//This return of the suite_rmsd of the current rebuild residue
//In the future...change name of suite_rmsd to suite_rmsd of rebuild residue and make it call suite_rmsd_any_seq_num
Real
suite_rmsd(pose::Pose const &  pose1,pose::Pose const &  pose2, Size const & current_rebuild_num) {
		using namespace scoring::constraints;
		using namespace conformation;
		using namespace basic::options;
		using namespace id;

//		std::cout << "Enter Suite rmsd function" << std::endl;

		Size atom_count=0;
  	Real sum_sd=0;

		Size const full_pose_reb_res=get_full_pose_reb_res(current_rebuild_num);
		Size const suite_num = Convert_to_partial_pose_seq_num(full_pose_reb_res, current_rebuild_num);
		bool const prepend_res = 	Is_prepend(current_rebuild_num); //Once make generalize so that can return any residue suite_rmsd, this option will be determine by input full_pose_seq_num
		//Will have to decide on the definition of suite_rmsd for the unrebuild residues...

		suite_square_deviation(pose1, pose2, prepend_res, suite_num, atom_count, sum_sd, false);

		sum_sd=sum_sd/(atom_count);
  	Real rmsd=sqrt(sum_sd);

		return (std::max(0.01, rmsd));

}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
Calculate_rmsd(output_data_struct& output_data, pose::Pose const & current_pose, pose::Pose const & partial_original_pose, pose::Pose const & partial_minimized_original_pose, pose::Pose const & partial_minimized_then_delete_pose, Size current_rebuild_num, bool verbose){

			using namespace basic::options;
			using namespace basic::options::OptionKeys;

			bool const prepend_res = 	Is_prepend(current_rebuild_num);

			Size reb_res= get_rebuild_pose_res_num(current_rebuild_num);

		  output_data.rmsd = suite_rmsd(partial_original_pose, current_pose, current_rebuild_num);
      output_data.rmsd_wrt_minimized = suite_rmsd(partial_minimized_original_pose, current_pose, current_rebuild_num);
			output_data.rmsd_wrt_correct_minimized = suite_rmsd(partial_minimized_then_delete_pose, current_pose, current_rebuild_num);


			if(verbose==true){

				std::cout << "loop_rmsd_wrt_to_original_pose" << std::endl;
				output_data.loop_rmsd= Loop_rmsd(partial_original_pose, current_pose, current_rebuild_num, false);

				std::cout << "loop_rmsd_wrt_to_minimized" << std::endl;
				output_data.loop_rmsd_wrt_minimized=Loop_rmsd(partial_minimized_original_pose, current_pose, current_rebuild_num, false);

				std::cout << "loop_rmsd_wrt_correct_minimized" << std::endl;
				output_data.loop_rmsd_wrt_correct_minimized= Loop_rmsd(partial_minimized_then_delete_pose, current_pose, current_rebuild_num, true);

			} else {

				output_data.loop_rmsd= Loop_rmsd(partial_original_pose, current_pose, current_rebuild_num, false);

				output_data.loop_rmsd_wrt_minimized=Loop_rmsd(partial_minimized_original_pose, current_pose, current_rebuild_num, false);

				output_data.loop_rmsd_wrt_correct_minimized= Loop_rmsd(partial_minimized_then_delete_pose, current_pose, current_rebuild_num, false);

			}


}



//////////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1<node_state>
cluster_pose_with_score_priority(utility::vector1<pose::Pose>& cluster_poses, utility::vector1<node_state>& poses_state, Size current_rebuild_num){
	using namespace core::scoring;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	Real rmsd_cutoff= option[cluster_rmsd];
//	bool const all_atom = option[ all_atom_cluster] ;


	for(Size i=1; i<=cluster_poses.size(); i++){
		if(poses_state[i]==NOT_INIT){
			poses_state[i]=OPEN;
			for(Size j=i+1; j<=cluster_poses.size(); j++){



//				if(all_atom==false){
//					Size full_pose_reb_res=get_full_pose_reb_res(current_rebuild_num);
//					Size reb_res=Convert_to_partial_pose_seq_num(full_pose_reb_res, current_rebuild_num);
					if(suite_rmsd(cluster_poses[i], cluster_poses[j], current_rebuild_num) > rmsd_cutoff)	continue;
//				}

				std::cout << "pose" << j << " is a neighbor of pose " << i << std::endl;
				poses_state[j]=CLOSE;
			}
		}
	}

	return poses_state;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
cluster_last_pose_old(utility::vector1<temp_pose_holder_struct>& cluster_poses, Size& i, Size current_rebuild_num){

	using namespace core::scoring;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	Real rmsd_cutoff= option[cluster_rmsd];

	Size full_pose_reb_res=get_full_pose_reb_res(current_rebuild_num);
	Size reb_res=Convert_to_partial_pose_seq_num(full_pose_reb_res, current_rebuild_num);
	std::cout << "current_rebuild_num: " << current_rebuild_num;
	std::cout << " full_pose_reb_res= "<< full_pose_reb_res;
	std::cout << " reb_res= "<< reb_res << std::endl;
	std::cout << "In cluster_last_pose function, clustering pose: " << i << std::endl;

//	for(Size j=1; j<i; j++){
//		std::cout << "cluster_poses[" << j << "].state: " << Print_node_state(cluster_poses[j].state) << std::endl;
//	}

	for(Size j=1; j<i; j++){
		if(cluster_poses[j].state==OPEN){
			Real rmsd=suite_rmsd(cluster_poses[i].pose, cluster_poses[j].pose, current_rebuild_num);
//			std::cout << "The rmsd between pose "<<  i << " and pose "<< j << " is " << rmsd << std::endl;
			if(rmsd < rmsd_cutoff)	{
				std::cout << "pose " << i << " is a neighbor of pose " << j << " rmsd= " << rmsd << std::endl;
				cluster_poses[i].state=CLOSE;
				return false;
			}
    }
	}

	cluster_poses[i].state=OPEN;
	return true;

}

bool
Is_part_of_cluster(pose::Pose&  current_pose, pose::Pose&  cluster_center_pose, Size current_rebuild_num, Size current_pose_num, Size cluster_center_pose_num){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	//Check that every residue rebuilded of current_pose is individually within rmsd_cutoff of cluster_center_pose residues.
	Real rmsd_cutoff= option[cluster_rmsd];

	utility::vector1< Real > rmsd_list(current_rebuild_num, 9999.99);

	for(Size i=1; i<= current_rebuild_num; i++){
 		Size full_pose_seq_num=get_full_pose_reb_res(i);
		Size seq_num=Convert_to_partial_pose_seq_num(full_pose_seq_num , current_rebuild_num);

		bool prepend_res=Is_prepend(i);
//		rmsd_list[i]=suite_rmsd_any_seq_num(current_pose, cluster_center_pose, i, full_pose_seq_num);
		rmsd_list[i]=suite_rmsd_any_seq_num(current_pose, cluster_center_pose, seq_num, prepend_res);

		if(rmsd_list[i]>rmsd_cutoff) return false;
	}
	//			Real rmsd=suite_rmsd(current_pose, (*clustered_poses_data[j].pose_OP), current_rebuild_num);

	std::cout << "pose " << current_pose_num << " is a neighbor of pose " << cluster_center_pose_num << std::endl;

	for(Size i=1; i<= current_rebuild_num; i++){
 		Size full_pose_seq_num=get_full_pose_reb_res(i);
		Size seq_num=Convert_to_partial_pose_seq_num(full_pose_seq_num , current_rebuild_num);

		std::cout << "full_pose_seq_num= " << full_pose_seq_num << " seq_num= " << seq_num;
		std::cout << " rmsd_list[" << i << "]= " << rmsd_list[i] << std::endl;
	}

	return true;

}

bool
cluster_last_pose(pose::Pose&  current_pose,  pose_data_struct & current_pose_data, utility::vector1< pose_data_struct>& clustered_poses_data, Size current_rebuild_num){

		using namespace core::scoring;
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

//		Real rmsd_cutoff= option[cluster_rmsd];

		Size full_pose_reb_res=get_full_pose_reb_res(current_rebuild_num);
		Size reb_res=Convert_to_partial_pose_seq_num(full_pose_reb_res, current_rebuild_num);
		std::cout << "current_rebuild_num: " << current_rebuild_num;
		std::cout << " full_pose_reb_res= "<< full_pose_reb_res;
		std::cout << " reb_res= "<< reb_res << std::endl;
		std::cout << "In cluster_last_pose function, clustering pose: " << current_pose_data.pose_num << std::endl;

	for(Size j=1; j<=clustered_poses_data.size(); j++){

			if(Is_part_of_cluster(current_pose, (*clustered_poses_data[j].pose_OP), current_rebuild_num,  current_pose_data.pose_num, clustered_poses_data[j].pose_num)==true){
				return false; //False since already part of a cluster
			}
	}

	clustered_poses_data.push_back(current_pose_data);
	Size new_element_pos=clustered_poses_data.size();
	clustered_poses_data[new_element_pos].pose_OP = new pose::Pose;
	//Do I have to clean this up...or is it done automatically? Does garbage collection automatically delete the pose object once there is no Pose_OP pointing to it?
	(*clustered_poses_data[new_element_pos].pose_OP)=	current_pose;

	std::cout << "pose " << current_pose_data.pose_num << " was added to clustered_poses_data vector" << std::endl;

	return true;

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
central_cluster_pose(utility::vector1< pose_data_struct >& 	all_poses_data, Size num_branch, Size current_rebuild_num)	{

  using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace core::kinematics;
	using namespace core::optimization;
	using namespace core::io::silent;
	using namespace core::id;
	using namespace protocols::rna;

	std::cout << "Enter central_cluster_pose function"  << std::endl;

  ResidueTypeSetCAP rsd_set; //This line is repeated
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "rna" ); //This line is repeated



	Size num_pose_clustered=0;
	utility::vector1< pose_data_struct> clustered_poses_data;

		Size i=0;
	  while(true){
	  	//Break if enough pose is found for branching or have exhuasted all poses.
	  	i++;
 			pose::Pose current_pose;
			pose_data_struct current_pose_data=all_poses_data[i];
  		std::string pose_name="";
  		pose_name.append(current_pose_data.tag);
  		pose_name.append(".pdb");
			std::cout << "Importing pose: " << pose_name << std::endl;
			core::import_pose::pose_from_pdb(current_pose, *rsd_set, pose_name);

////////This setup the virtual phosphate...shouldn't change result since suite rmsd does not include the 3' phosphate group in the cluster////////////////////////////
////////Haven't test this yet...will do so in the future///////////////////////////////////////////////////////////////////////////////////////////////////////////
//			if(i==1){
//				Setup_rebuild_pose_virtual_phosphate(current_pose, current_rebuild_num, true);
//			} else {
//				Setup_rebuild_pose_virtual_phosphate(current_pose, current_rebuild_num, true);
//			}
//////////////////////////////////////////////////////////////////

			if(cluster_last_pose(current_pose, current_pose_data, clustered_poses_data, current_rebuild_num)==true) num_pose_clustered++;
				std::cout << "Currently: " << num_pose_clustered  << " pose clusters center " << "from " << i << " poses" << std::endl;
				std::cout << "Check pose_num, i= " << i << "  pose_num= " << current_pose_data.pose_num << std::endl;
			if(num_pose_clustered==num_branch || all_poses_data.size()==i) break;
		}


		std::cout << "Final: " << num_pose_clustered  << " pose clusters center " << "from " << i << " poses" << std::endl;
		std::cout << "Check: cluster_poses.size()= " <<  clustered_poses_data.size() << std::endl;


	//This is the new data_structure that will be returned

		for(Size n=1; n<=clustered_poses_data.size(); n++) {
				clustered_poses_data[n].tag[0]='C';
				dump_pdb( (*clustered_poses_data[n].pose_OP), clustered_poses_data[n].tag+".pdb" );
		}

//Manual clean up of memory
//	 for(Size i=1; i<=clustered_poses.data; i++){
//		int *p_var = NULL;     // new pointer declared
//		p_var = new int;       // memory dynamically allocated
//	 	delete p_var;          // freed up memory

//	clustered_poses_data[new_element_pos].pose_OP = new pose::Pose;
//		delete clustered_poses_data[i].pose_OP
//	 	clustered_poses_data[i].pose_OP = NULL;          // pointer changed to NULL
//	 }

   all_poses_data.clear();

   all_poses_data=clustered_poses_data;

	std::cout << "Exit central_cluster_pose  function"  << std::endl;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
bool
unfinish_jobs_left(utility::vector1< job_info_struct>& job_status_list){

	for(Size i=1; i<=job_status_list.size(); i++) {
		if(job_status_list[i].is_done==false) return true;
	}
	return false;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
Is_Chain_Closable(pose::Pose current_pose){

		//Pseudoknot parameter
		Real O3_C5_distance = (current_pose.residue(18).xyz(5) - current_pose.residue(17).xyz(9) ).length();
		std::cout << "O3_C5_distance= " << O3_C5_distance << std::endl;
		if(O3_C5_distance>11.0138) return false;
		else return true;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
outfile_column_name(std::ofstream& outfile, Size& current_rebuild_num){
		using namespace basic::options;
	  using namespace basic::options::OptionKeys;

  	Size spacing=8;

//			outfile << std::setw(18) << std::left << "description";

		outfile << std::setw(spacing+2) << std::left << "rmsd";
	  outfile << std::setw(spacing+2) << std::left << "rmsd_m";
	  outfile << std::setw(spacing+2) << std::left << "rmsd_c";
	 	outfile << std::setw(spacing+2) << std::left << "score";
	 	outfile << std::setw(spacing+2) << std::left << "lrmsd";
	  outfile << std::setw(spacing+2) << std::left << "lrmsd_m";
	  outfile << std::setw(spacing+2) << std::left << "lrmsd_c";
	 	outfile << std::setw(spacing+2) << std::left << "cc_dist"; //Distance to close chain
		outfile << std::setw(spacing+2) << std::left << "tot";


//	 		outfile << std::setw(spacing+2) << std::left << "O5";
//	 		outfile << std::setw(spacing+2) << std::left << "P";



	 	if(Is_prepend(current_rebuild_num)==false){

 			outfile << std::setw(spacing) << std::left << "ai";
 			outfile << std::setw(spacing) << std::left << "bi";
 			outfile << std::setw(spacing) << std::left << "gi";
 			outfile << std::setw(spacing) << std::left << "deli";
 			outfile << std::setw(spacing) << std::left << "ei-1";
 			outfile << std::setw(spacing) << std::left << "zi-1";
 			outfile << std::setw(spacing) << std::left << "chii";
 			outfile << std::setw(spacing) << std::left << "nu1i";
 			outfile << std::setw(spacing) << std::left << "nu2i";
 			outfile << std::setw(spacing) << std::left << "2OHi";

		} else {

 			outfile << std::setw(spacing) << std::left << "ai+1";
 			outfile << std::setw(spacing) << std::left << "bi+1";
 			outfile << std::setw(spacing) << std::left << "gi+1";
 			outfile << std::setw(spacing) << std::left << "deli";
 			outfile << std::setw(spacing) << std::left << "ei";
 			outfile << std::setw(spacing) << std::left << "zi";
 			outfile << std::setw(spacing) << std::left << "chii";
 			outfile << std::setw(spacing) << std::left << "nu1i";
 			outfile << std::setw(spacing) << std::left << "nu2i";
 			outfile << std::setw(spacing) << std::left << "2OHi";
 		}

		outfile << "\n";

}

void
create_column_pos_map(std::map< std::string, Size >& column_pos_map){
	using namespace std;

	column_pos_map.insert ( std::pair< string, Size>("description", 0 ));
	std::cout <<"description is located at: " << column_pos_map["description"]<< std::endl;
	column_pos_map.insert ( std::pair< string, Size>("rmsd", 1 ));
	std::cout << "rmsd is located at: " << column_pos_map["rmsd"]<< std::endl;
	column_pos_map.insert ( std::pair< string, Size>("rmsd_min", 2 ));
	std::cout << "rmsd_min is located at: " << column_pos_map["rmsd_min"]<< std::endl;
	column_pos_map.insert ( std::pair< string, Size>("rmsd_cor", 3 ));
	std::cout << "rmsd_correct is located at: " << column_pos_map["rmsd_cor"]<< std::endl;
	column_pos_map.insert ( std::pair< string, Size>("score", 4 ));
	std::cout << "score is located at: " << column_pos_map["score"]<< std::endl;
	column_pos_map.insert ( std::pair< string, Size>("loop_rmsd", 5 ));
	std::cout << "loop_rmsd is located at: " << column_pos_map["loop_rmsd"]<< std::endl;
	column_pos_map.insert ( std::pair< string, Size>("loop_rmsd_min", 6 ));
	std::cout << "loop_rmsd_min is located at: " << column_pos_map["loop_rmsd_min"]<< std::endl;
	column_pos_map.insert ( std::pair< string, Size>("loop_rmsd_cor", 7 ));
	std::cout << "loop_rmsd_correct is located at: " << column_pos_map["loop_rmsd_cor"]<< std::endl;
	column_pos_map.insert ( std::pair< string, Size>("dist_to_close_chain", 8 ));
	std::cout << "distance left to close chain: " << column_pos_map["dist_to_close_chain"]<< std::endl;
	column_pos_map.insert ( std::pair< string, Size>("diff_torsions", 9 ));
	std::cout << "diff_torsions is located at: " << column_pos_map["diff_torsions"]<< std::endl;

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
central_evaluation(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace core::kinematics;
	using namespace core::optimization;
	using namespace core::io::silent;
	using namespace core::id;
	using namespace protocols::rna;
	using namespace std;

///////////Create RunTimeParameters objects which is will referenced throughout the code/////////////////////////////////////////////////////////////////////////
	Initialize_RunTimeParameters();

////////////////////////////////////////Options///////////////////////////////////////////////////////////////////

	bool const apply_central_clustering = option[central_clustering];
  Size const total_nodes = option[total_nodes_input] ;
  Size const num_branch = option[num_branch_kept] ;
	Size const first_rebuild_num = option[hack_first_rebuild_num];
	bool const verbose_option = option[Verbose];
	Size num_residue_reb = RunTimeParameters::rebuild_residue_list.size();
	Size const total_rotamer_group=get_total_rotamer_group();

////////////////////////////////////////Residue rebuild for loop/////////////////////////////////////////////////////


	//////////create a map that associate column names with their position in data.txt/////////////

	std::map< std::string, Size > column_pos_map;
	create_column_pos_map(column_pos_map);

	/////////////////////////////////////////////////////////////////////////////////////////////////////
	//Jobs for rebuild_num=1
	utility::vector1< job_info_struct> job_status_list;
	for(Size i=1; i<= total_nodes; i++){
		job_info_struct job_info;
		job_info.node_num=i;
		job_info.job_number=i;
		job_info.communication_num=1;
		job_info.is_done=false;
		job_status_list.push_back(job_info);
	}

	for(Size current_rebuild_num=first_rebuild_num ; current_rebuild_num <= num_residue_reb; current_rebuild_num++) {


		bool const prepend_res = 	Is_prepend(current_rebuild_num);

		Size full_pose_reb_res=get_full_pose_reb_res(current_rebuild_num);
		Size reb_res= get_rebuild_pose_res_num(current_rebuild_num); //Reb_res_number in terms of the partially rebuilt pose.

		std::cout << "current_rebuild_num= " << current_rebuild_num << std::endl;
		std::cout << "full_pose_reb_res= " << full_pose_reb_res << std::endl;
		std::cout << "rebuild_pose_reb_res= " << reb_res << std::endl;

		utility::vector1< pose_data_struct > all_poses_data;

 		Size counter=0;

		while(unfinish_jobs_left(job_status_list)){

			counter++;
			if(counter==job_status_list.size()+1) counter=1; //iterate through all the members of job_status_list.

			if(Is_job_finish(job_status_list[counter], current_rebuild_num) && job_status_list[counter].is_done==false) {

  			std::cout << "Job: " << job_status_list[counter].job_number << "( counter= " << counter;
				std::cout << ", node_num= " << job_status_list[counter].node_num << ", communication_num= " << job_status_list[counter].communication_num << ") is finish" << std::endl;

  			job_status_list[counter].is_done=true;

			//////////////////////////////Open data.txt//////////////////////////////////////////
				std::string filename;
				create_filename(filename, "data", job_status_list[counter].job_number, current_rebuild_num);

				std::ifstream infile;
		 		infile.open(filename.c_str());
				if (infile.fail()) {
  				std::cout << "Error! \"" << filename << "\" could not be opened!" << std::endl;
  				exit (1);
				}
				std::cout << "Open \"" << filename << "\" successful!" << std::endl;

			/////////Import the data from data.txt into all_poses_data, each line in data.txt represent 1 pose_data struct////

			import_pose_data(all_poses_data,  infile, column_pos_map, 'M');

			infile.close();
			} else {
				suspend(1); //Suspend code for about 1 seconds
			} //End of if(Is_job_finish) statement

		}// End of while loop to check if all jobs are finish.

		sort_pose_data(all_poses_data);



//		Rescore_pose(all_poses_data, current_rebuild_num);



		std::string sort_filename;
		create_filename(sort_filename, "sort", 0 , current_rebuild_num);

		std::ofstream sort_outfile;
	  sort_outfile.open(sort_filename.c_str());

	  Size spacing=8;
	  Size tag_spacing;
		Size char_per_line=207;
	  if(all_poses_data[1].tag.length()<char_per_line){
  		tag_spacing=char_per_line;
  	} else {
  		tag_spacing=2*char_per_line;
  	}

    outfile_column_name(sort_outfile, current_rebuild_num);
	  for(Size i=1; i <= all_poses_data.size(); i++){
	  	sort_outfile << std::setw(tag_spacing) << std::left << all_poses_data[i].tag;
			sort_outfile << std::setw(spacing+2) << std::fixed << std::setprecision(3)  << std::left << all_poses_data[i].rmsd;
			sort_outfile << std::setw(spacing+2) << std::fixed << std::setprecision(3)  << std::left << all_poses_data[i].rmsd_min;
			sort_outfile << std::setw(spacing+2) << std::fixed << std::setprecision(3)  << std::left << all_poses_data[i].rmsd_cor;
 	 		sort_outfile << std::setw(spacing+2) << std::fixed << std::setprecision(2)  << std::left << all_poses_data[i].score;
 	 		sort_outfile << std::setw(spacing+2) << std::fixed << std::setprecision(3)  << std::left << all_poses_data[i].loop_rmsd;
			sort_outfile << std::setw(spacing+2) << std::fixed << std::setprecision(3)  << std::left << all_poses_data[i].loop_rmsd_min;
			sort_outfile << std::setw(spacing+2) << std::fixed << std::setprecision(3)  << std::left << all_poses_data[i].loop_rmsd_cor;
 	 		sort_outfile << std::setw(spacing+2) << std::fixed << std::setprecision(2)  << std::left << all_poses_data[i].dist_to_close_chain;
			sort_outfile << std::setw(spacing+2) << std::fixed << std::setprecision(1)  << std::left << all_poses_data[i].diff_torsions;
 	 		sort_outfile << "\n";
		}


		sort_outfile.close();

//		if(apply_central_clustering && (all_poses_data.size()>=4*num_branch)) { //Check if there is enough pose to perform clustering.

			central_cluster_pose(all_poses_data, num_branch, current_rebuild_num);
			std::cout << "Check: all_poses_data.size()= " <<  all_poses_data.size() << std::endl;

			std::string cluster_filename;
			create_filename(cluster_filename, "cluster", 0 , current_rebuild_num);

			std::ofstream cluster_outfile;
	  	cluster_outfile.open(cluster_filename.c_str());

			outfile_column_name(cluster_outfile, current_rebuild_num);
	  	for(Size i=1; i <= all_poses_data.size(); i++){

	  		cluster_outfile << std::setw(tag_spacing) << std::left << all_poses_data[i].tag;
				cluster_outfile << std::setw(spacing+2) << std::fixed << std::setprecision(3)  << std::left << all_poses_data[i].rmsd;
				cluster_outfile << std::setw(spacing+2) << std::fixed << std::setprecision(3)  << std::left << all_poses_data[i].rmsd_min;
				cluster_outfile << std::setw(spacing+2) << std::fixed << std::setprecision(3)  << std::left << all_poses_data[i].rmsd_cor;
 	 			cluster_outfile << std::setw(spacing+2) << std::fixed << std::setprecision(2)  << std::left << all_poses_data[i].score;
 	 			cluster_outfile << std::setw(spacing+2) << std::fixed << std::setprecision(3)  << std::left << all_poses_data[i].loop_rmsd;
				cluster_outfile << std::setw(spacing+2) << std::fixed << std::setprecision(3)  << std::left << all_poses_data[i].loop_rmsd_min;
				cluster_outfile << std::setw(spacing+2) << std::fixed << std::setprecision(3)  << std::left << all_poses_data[i].loop_rmsd_cor;
 	 			cluster_outfile << std::setw(spacing+2) << std::fixed << std::setprecision(2)  << std::left << all_poses_data[i].dist_to_close_chain;
				cluster_outfile << std::setw(spacing+2) << std::fixed << std::setprecision(1)  << std::left << all_poses_data[i].diff_torsions;
 	 			cluster_outfile << "\n";
			}

			cluster_outfile.close();
//		} else {
//			output_branching_pose(all_poses_data, num_branch, current_rebuild_num);
//		}





/////////////Print out sort.txt to screen///////////////////////////////////////////////////////////////

  for(Size i=1; i <= all_poses_data.size(); i++){
		std::cout << std::setw(3) << "pose " << i << "  ";
		std::cout << std::setw(spacing+2) << std::fixed << std::setprecision(3)  << std::left << all_poses_data[i].rmsd;
		std::cout << std::setw(spacing+2) << std::fixed << std::setprecision(3)  << std::left << all_poses_data[i].rmsd_min;
		std::cout << std::setw(spacing+2) << std::fixed << std::setprecision(3)  << std::left << all_poses_data[i].rmsd_cor;
 	 	std::cout << std::setw(spacing+2) << std::fixed << std::setprecision(2)  << std::left << all_poses_data[i].score;
 	 	std::cout << std::setw(spacing+2) << std::fixed << std::setprecision(3)  << std::left << all_poses_data[i].loop_rmsd;
		std::cout << std::setw(spacing+2) << std::fixed << std::setprecision(3)  << std::left << all_poses_data[i].loop_rmsd_min;
		std::cout << std::setw(spacing+2) << std::fixed << std::setprecision(3)  << std::left << all_poses_data[i].loop_rmsd_cor;
 	 	std::cout << std::setw(spacing+2) << std::fixed << std::setprecision(2)  << std::left << all_poses_data[i].dist_to_close_chain;
		std::cout << std::setw(spacing+2) << std::fixed << std::setprecision(1)  << std::left << all_poses_data[i].diff_torsions;
 	 	std::cout << std::endl;
	}

/////////////////////////////////Create job files for next_rebuild_num///////////////////////////////////////////////////////////////

		std::cout << "Jobs submission"  << std::endl;

		if(current_rebuild_num==num_residue_reb) continue; //No job_file for last_residue_rebuild!!

//		utility::vector1< job_info_struct> job_status_list;
		job_status_list.clear(); //clean jobs of previous residue rebuild

		//Need to submit num_branch jobs
		utility::vector1 <Size> communication_num_list(total_nodes, 0);

		Size job_number=1;
		Size node_num=1;
		while(job_number <= std::min( num_branch, all_poses_data.size() ) ){

			if(node_num==(total_nodes+1)) node_num=1;	//iterate over all the nodes

				//Check if node is ready to receive new job

				if(Is_node_waiting_for_new_job(communication_num_list, node_num, current_rebuild_num+1)){
					communication_num_list[node_num]++;
					std::cout << "Submitting new job to node_num= " << node_num << " communication_num= " << communication_num_list[node_num] << std::endl;

					job_info_struct job_info;
					job_info.node_num=node_num;
					job_info.job_number=job_number;
					job_info.communication_num=communication_num_list[node_num];
					job_info.is_done=false;
					job_status_list.push_back(job_info);

					std::string job_filename;
					create_communication_filename(job_filename, "job_file", job_info.communication_num, job_info.node_num, current_rebuild_num+1);
					std::ofstream job_outfile;
					job_outfile.open(job_filename.c_str());
					job_outfile << std::setw(tag_spacing) << std::left << all_poses_data[job_number].tag;
				  job_outfile << std::setw(20) << std::left << lead_zero_string_of(1, 3 ); //Obsolete, in the past use to split rotamer groups to many nodes
  				job_outfile << std::setw(20) << std::left << lead_zero_string_of(total_rotamer_group, 3 );  //Obsolete, in the past use to split rotamer groups to many nodes
					job_outfile << std::setw(20) << std::left << job_info.job_number;
					job_outfile << "\n";

					job_outfile.close();
					job_number++;
				}
			node_num++;
		}

		//Communicate to node that all jobs have been submitted-> move on to next rebuild num
		for(Size node_num=1; node_num<=total_nodes; node_num++){
				std::string job_filename;
				create_communication_filename(job_filename, "job_file", communication_num_list[node_num]+1, node_num, current_rebuild_num+1);
				std::ofstream job_outfile;
				job_outfile.open(job_filename.c_str());
				job_outfile << std::setw(tag_spacing) << std::left << "done";
				job_outfile.close();
		}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
	}//End of num_residue_rebuilt for loop

} //End of Central_evaluation function.

//////////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
get_rotamers_from_string(std::string& rotamers_string, utility::vector1< Size >& rotamers, Size current_rebuild_num, Size digits ){


		for(Size n=1; n < current_rebuild_num; n++){
//			std::cout << "n= " << n << std::endl;
//			std::cout << "rotamers_string=" << rotamers_string;
			std::string rotamer_string=rotamers_string.substr(digits*(n-1), digits);
//			std::cout << "rotamer_string= " << rotamer_string;

			rotamers[n]= convert_string_to_int(rotamer_string);

			std::cout << "rotamer= " << rotamers[n] << std::endl;
		}
	return rotamers;
}
///////////////////////////////////////////////////////////////////////////////////////////
void
Create_status_file(Size const  node_number, Size const communication_num, Size const current_rebuild_num){

  	std::string status_filename;
		create_communication_filename(status_filename, "status", communication_num, node_number, current_rebuild_num);
		std::ofstream status_outfile;
 		status_outfile.open(status_filename.c_str());
 		status_outfile << "finish" ;
 	 	status_outfile.close();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool
read_in_job_file(job_file_data& job_data, Size const node_number, Size const communication_num, Size const num_residue_reb, Size const current_rebuild_num) {

		using namespace std;

		Create_status_file(node_number, communication_num, current_rebuild_num); //Signal center computer that this node is ready to accept a new job.


		std::string job_filename;
		create_communication_filename(job_filename, "job_file", communication_num+1, node_number, current_rebuild_num);
		std::ifstream infile;
		while(true){

	 	  infile.open(job_filename.c_str());
  		if (infile.fail()) {
  			std::cout << "Could no open file named \"" << job_filename << "\". central_evaluation still running" << std::endl;
				suspend(10); //Suspend code for about 2.5 seconds before retry.
				infile.clear();
  			infile.close();
			} else {
				suspend(10); //Give time for central_computer to output job information after it has created the job_file.
				break;
			}
		}

		std::string job_string;
		getline(infile, job_string);
		infile.close(); //close infile after obtaining job_information

		std::vector<std::string> tokens;
  	Tokenize(job_string, tokens);

		for( Size i = 0; i < tokens.size(); i++) {
				std::cout << "Token " << i << ": " << tokens[i] << std::endl;
		}

		job_data.tag=tokens[0];

		if(job_data.tag=="done"){//No more jobs for this node ..skip loop to next rebuild_num
				std::cout << "no more jobs for this node ..skip to next residue rebuild" << std::endl;
//				Create_status_file(node_number, communication_num, current_rebuild_num); //not tested yet!
				return false;
		}

		job_data.rotamer_group_min= convert_string_to_int(tokens[1]);
		job_data.rotamer_group_max= convert_string_to_int(tokens[2]);
		job_data.job_number= convert_string_to_int(tokens[3]);


		std::cout << "rotamer_group_min " << job_data.rotamer_group_min << std::endl;
		std::cout << "rotamer_group_max " << job_data.rotamer_group_max << std::endl;
		std::cout << "job_number " << job_data.job_number << std::endl;

		return true;

}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//This function have been updated to handle multiple chain_breaks...will return the distance to close chain from the rebuild residue.
Real
Calculate_dist_to_close_chain(pose::Pose const & pose, Size  current_rebuild_num){
		using namespace basic::options;
	  using namespace basic::options::OptionKeys;
		using namespace core::chemical;
		using namespace core::conformation;
		using namespace core::kinematics;

		Size full_pose_reb_res=get_full_pose_reb_res(current_rebuild_num);

		utility::vector1 < utility::vector1 <Residue_info_struct> > rebuild_residue_group_list;
		create_rebuild_residue_groups(rebuild_residue_group_list); //Since calculate_dist_to_close_chain is called very often..creating this vector might slow down code.

		//These include only cutpoint create for loop rebuilding
		utility::vector1 <std::pair<Size, Size> > jump_points_list;
		get_jump_points(jump_points_list, rebuild_residue_group_list);

		Size chain_break_num=1;
		for(; chain_break_num<=jump_points_list.size(); chain_break_num++){
			if((full_pose_reb_res>jump_points_list[chain_break_num].first) && (full_pose_reb_res<jump_points_list[chain_break_num].second)){
					break;//if condition is satisfied, then we know that reb_res is part of this particular chainbreak/loop
			}
		}

		Size five_prime_jump_point=jump_points_list[chain_break_num].first;
		Size three_prime_jump_point=jump_points_list[chain_break_num].second;

		//The residues to be use to calculate the chain_break distance;
		Size five_prime_residue;
		Size three_prime_residue;

		if(Is_prepend(current_rebuild_num)){

			five_prime_residue=five_prime_jump_point;
			three_prime_residue=full_pose_reb_res;

			//Find five_prime_residue
					//Could create rebuild_residue_groups for the partially rebuild pose...but that will be computational expensive...
					//On the other hand would be more efficient in that only have to search within the group
			for(Size i=1; i<=current_rebuild_num; i++){

				if((RunTimeParameters::rebuild_residue_list[i].seq_num > five_prime_jump_point) && (RunTimeParameters::rebuild_residue_list[i].seq_num < three_prime_residue)){

					if(RunTimeParameters::rebuild_residue_list[i].seq_num > five_prime_residue) five_prime_residue=RunTimeParameters::rebuild_residue_list[i].seq_num;
				}
			}

		}else {


			five_prime_residue=full_pose_reb_res;
			three_prime_residue=three_prime_jump_point;

			//Find three_prime_residue
			for(Size i=1; i<=current_rebuild_num; i++){

				if((RunTimeParameters::rebuild_residue_list[i].seq_num > five_prime_residue) && (RunTimeParameters::rebuild_residue_list[i].seq_num < three_prime_jump_point)){

					if(RunTimeParameters::rebuild_residue_list[i].seq_num < three_prime_residue) three_prime_residue=RunTimeParameters::rebuild_residue_list[i].seq_num;
				}
			}
		}

	Size rebuild_pose_five_prime_residue= Convert_to_partial_pose_seq_num(five_prime_residue, current_rebuild_num);
	Size rebuild_pose_three_prime_residue= Convert_to_partial_pose_seq_num(three_prime_residue, current_rebuild_num);

	//5 is the C5' atom. 9 is the O3' atom.
	Real dist_to_close_chain=(pose.residue(rebuild_pose_three_prime_residue).xyz("C5*") - pose.residue(rebuild_pose_five_prime_residue).xyz("O3*") ).length();


/*
	utility::vector1 <Residue_info_struct> missing_residue_list;
	get_missing_residue_list(missing_residue_list, current_rebuild_num);
	std::cout << "In Calculate_dist_to_close_chain function" << std::endl;
	std::cout << "  Is_prepend: ";
	if(Is_prepend(current_rebuild_num)==true) std::cout << "true" << std::endl;
	if(Is_prepend(current_rebuild_num)==false) std::cout << "false" << std::endl;
	std::cout << "  full_pose_reb_res= " << full_pose_reb_res << std::endl;
	std::cout << "  current_rebuild_num= " << current_rebuild_num;
	std::cout << "  chain_break_num= " << chain_break_num << std::endl;
	std::cout << "  five_prime_jump_point= " << five_prime_jump_point;
	std::cout << "  three_prime_jump_point= " << three_prime_jump_point << std::endl;
	std::cout << "  full_pose_five_prime_residue: " << five_prime_residue << " full_pose_three_prime_residue: " << three_prime_residue <<std::endl;
	std::cout << "  rebuild_pose_five_prime_residue: " << rebuild_pose_five_prime_residue << " rebuild_pose_three_prime_residue: " << rebuild_pose_three_prime_residue <<std::endl;
	std::cout <<  std::setw(40)  << "  rebuild_residue_list: "; Output_residue_list(RunTimeParameters::rebuild_residue_list);
	std::cout <<  std::setw(40)  << "  missing_residue_list: "; Output_residue_list(missing_residue_list);
	std::cout << "  dist_to_close_chain: " << dist_to_close_chain << std::endl;
*/

	return dist_to_close_chain;

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Called when the before_last_residue is rebuilt. Check that distance is reasonable for chain closure.
//This function have been updated to handle multiple chain_breaks, June 14, 2009
bool
Check_chain_closable(pose::Pose const & pose, count_struct& count_data, std::string phase, Size& current_rebuild_num){

		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		Size* pass_count_P;
		Size*  total_count_P;
		Real cutoff_distance;

		if(phase=="Sampling"){
				pass_count_P=&count_data.pass_before_last_res_reb_chain_break_U;
				total_count_P=&count_data.total_before_last_res_reb_chain_break_U;
				cutoff_distance=11.0138;
				//Should make this value slightly larger to account for possible change in position due to minimization?
		} else if(phase=="Output"){
				pass_count_P=& count_data.pass_before_last_res_reb_chain_break_M;
				total_count_P=& count_data.total_before_last_res_reb_chain_break_M;
				cutoff_distance=11.0138;
				//11.0138 Angstrom is the distance of the fully strech chain from C5' to a O3' atom two residue before it
		} else {
			std::cout << "Error in Check chain closable function" << std::endl;
			exit (1);
		}


		//5 is the C5' atom. 9 is the O3' atom.
		Real dist_to_close_chain = Calculate_dist_to_close_chain(pose, current_rebuild_num);
		(*total_count_P)++;
		if(dist_to_close_chain<cutoff_distance){
				(*pass_count_P)++;
//				if(option[Verbose]){
//					std::cout << "before_last_res_reb_chain_closure_check: " << phase;
//					std::cout << " cutoff_distance= " << cutoff_distance;
//					std::cout << " dist_to_close_chain= " << dist_to_close_chain;
//					std::cout << " pass_count= " << (*pass_count_P);
//					std::cout << " total_count= " << (*total_count_P) << std::endl;
//				}
				return true;
		}
		return false;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
Output_data_general(core::io::silent::SilentFileData& silent_file_data, std::string const & silent_file, output_data_struct& output_data, pose::Pose const & current_pose, std::string const & tag, Size const & current_rebuild_num, std::ofstream& outfile, std::vector< Real> & diff_torsions , std::vector<Real> const *temp_vector = 0)
{

	using namespace core::io::silent;
  using namespace core::scoring;
  using namespace basic::options;
	using namespace basic::options::OptionKeys;


//	std::string rotamer_set_option = option[rotamer_set];
	Size const write_score_only = option[print_score_only] ;


  Size full_pose_reb_res=get_full_pose_reb_res(current_rebuild_num);
	Size reb_res=Convert_to_partial_pose_seq_num(full_pose_reb_res, current_rebuild_num);
	Real dist_to_close_chain=Calculate_dist_to_close_chain(current_pose, current_rebuild_num);


  ////////////////////////////write output/////////////////////////////////
  	Size spacing=8;
  	Size tag_spacing;
		Size char_per_line=207;
  	if(tag.length()<(char_per_line)){
  		tag_spacing=char_per_line;
  	} else {
  		tag_spacing=2*char_per_line;
  	}

  	if(output_data.current_score>-1) {
  		std::cout << "pose: " << tag << " has very bad score" << std::endl;
  		return; //Prevent code from outputting very bad score pose.
  	}


	  outfile << std::setw(tag_spacing) << std::left << tag;            // output description
    outfile << std::setw(spacing+2) << std::fixed << std::setprecision(3) << std::left << output_data.rmsd;
	  outfile << std::setw(spacing+2) << std::fixed << std::setprecision(3) << std::left << output_data.rmsd_wrt_minimized;
		outfile << std::setw(spacing+2) << std::fixed << std::setprecision(3) << std::left << output_data.rmsd_wrt_correct_minimized;
  	outfile << std::setw(spacing+2) << std::fixed << std::setprecision(2) << std::left << output_data.current_score;
 		outfile << std::setw(spacing+2) << std::fixed << std::setprecision(3) << std::left  << output_data.loop_rmsd;
	  outfile << std::setw(spacing+2) << std::fixed << std::setprecision(3) << std::left  << output_data.loop_rmsd_wrt_minimized;
		outfile << std::setw(spacing+2) << std::fixed << std::setprecision(3) << std::left  << output_data.loop_rmsd_wrt_correct_minimized;
		outfile << std::setw(spacing+2) << std::fixed << std::setprecision(2) << std::left << dist_to_close_chain;
		outfile << std::setw(spacing+2) << std::fixed << std::setprecision(1) << std::left  << std::left << diff_torsions[0];//Total difference

			if(Is_prepend(current_rebuild_num)){
				conformation::Residue const & Five_prime_res=current_pose.residue(reb_res);
				conformation::Residue const & Three_prime_res=current_pose.residue(reb_res+1);

				outfile << std::setw(spacing) << std::fixed << std::setprecision(1) << std::left << Three_prime_res.mainchain_torsion(1); //Alpha
				outfile << std::setw(spacing) << std::fixed << std::setprecision(1) << std::left << Three_prime_res.mainchain_torsion(2); //Beta
				outfile << std::setw(spacing) << std::fixed << std::setprecision(1) << std::left << Three_prime_res.mainchain_torsion(3); //Gamma
				outfile << std::setw(spacing) << std::fixed << std::setprecision(1) << std::left  << Five_prime_res.mainchain_torsion(4); //delta
				outfile << std::setw(spacing) << std::fixed << std::setprecision(1) << std::left  << Five_prime_res.mainchain_torsion(5); //epsilon
				outfile << std::setw(spacing) << std::fixed << std::setprecision(1) << std::left  << Five_prime_res.mainchain_torsion(6); //zeta
				outfile << std::setw(spacing) << std::fixed << std::setprecision(1) << std::left  << Five_prime_res.chi(1); //Chi
				outfile << std::setw(spacing) << std::fixed << std::setprecision(1) << std::left  << Five_prime_res.chi(2); //nu2
				outfile << std::setw(spacing) << std::fixed << std::setprecision(1) << std::left  << Five_prime_res.chi(3); //nu1
				outfile << std::setw(spacing) << std::fixed << std::setprecision(1) << std::left  << Five_prime_res.chi(4); //Chi_OH2

			}else {
				conformation::Residue const & Five_prime_res=current_pose.residue(reb_res-1);
				conformation::Residue const & Three_prime_res=current_pose.residue(reb_res);

				outfile << std::setw(spacing) << std::fixed << std::setprecision(1) << std::left  << Three_prime_res.mainchain_torsion(1); //Alpha
				outfile << std::setw(spacing) << std::fixed << std::setprecision(1) << std::left  << Three_prime_res.mainchain_torsion(2); //Beta
				outfile << std::setw(spacing) << std::fixed << std::setprecision(1) << std::left  << Three_prime_res.mainchain_torsion(3); //Gamma
				outfile << std::setw(spacing) << std::fixed << std::setprecision(1) << std::left  << Three_prime_res.mainchain_torsion(4); //Delta
				outfile << std::setw(spacing) << std::fixed << std::setprecision(1) << std::left << Five_prime_res.mainchain_torsion(5); //epsilon
				outfile << std::setw(spacing) << std::fixed << std::setprecision(1) << std::left  << Five_prime_res.mainchain_torsion(6); //zeta
				outfile << std::setw(spacing) << std::fixed << std::setprecision(1) << std::left  << Three_prime_res.chi(1); //Chi
				outfile << std::setw(spacing) << std::fixed << std::setprecision(1) << std::left  << Three_prime_res.chi(2); //nu2 //I might switch order of nu2 and nu1 before this.
				outfile << std::setw(spacing) << std::fixed << std::setprecision(1) << std::left  << Three_prime_res.chi(3); //nu1
				outfile << std::setw(spacing) << std::fixed << std::setprecision(1) << std::left  << Three_prime_res.chi(4); //Chi_OH2

			}

// 		if(temp_vector!=0){
//			for(Size i=0; i<(*temp_vector).size(); i++){
// 				outfile << std::setw(spacing+5) << std::setprecision(3) << std::left << (*temp_vector)[i];
//			}
// 		}
		outfile << "\n";

		//////////////////Create name of output pdb file////////////////////////////////

//		Real alpha=current_pose.residue(19).mainchain_torsion(1); //Alpha torsion  of rebuilt+1 residue*/
//		std::cout << "alpha 19=" << alpha << std::endl;

  	RNA_SilentStruct s( current_pose, tag );

		std::string name;
		name.append( "rmsd" );
		s.add_energy( name, output_data.rmsd );

		name.clear();
		name.append( "rms_min" );
		s.add_energy( name, output_data.rmsd_wrt_minimized);

		name.clear();
		name.append( "rmsd_cor" );
		s.add_energy( name, output_data.rmsd_wrt_correct_minimized);

		name.clear();
		name.append( "score2" );
		s.add_energy( name, output_data.current_score);

		name.clear();
		name.append( "dist_close_chain" );
		s.add_energy( name, dist_to_close_chain );

		if(temp_vector!=0){
			for(Size i=0; i<(*temp_vector).size(); i++){
				name.clear();
				name.append(lead_zero_string_of(i, 2));
				s.add_energy( name, (*temp_vector)[i]);
			}
		}

		silent_file_data.write_silent_struct(s, silent_file, write_score_only);


}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Updated June 14, however still need to check if this function is working properly
void
Calculate_diff_torsions(std::vector<Real>& diff_torsions, pose::Pose& current_pose, pose::Pose& original_pose, Size const & current_rebuild_num){

		using namespace core::scoring;
		using namespace core::scoring::rna;
		using namespace basic::options;
	  using namespace basic::options::OptionKeys;

//		std::cout << "In Calculate_diff_torsions function" << std::endl;


		bool const prepend_res = 	Is_prepend(current_rebuild_num);
		Size const full_pose_reb_res=get_full_pose_reb_res(current_rebuild_num);

//		std::cout << "current_rebuild_num=  " << current_rebuild_num << std::endl;
//		std::cout << "full_pose_reb_res=  " << full_pose_reb_res << std::endl;

		Size reb_res=Convert_to_partial_pose_seq_num(full_pose_reb_res, current_rebuild_num);

		utility::vector1 <Real> original_torsions(10 /*# torsion angles*/, 9999);
		utility::vector1 <Real> current_torsions(10 /*# torsion angles*/, 9999);

 		Size const nres( current_pose.total_residue() );
		Size const original_res= original_pose.total_residue();
		Size res_diff= original_res-nres;


		diff_torsions[0]=0;


		if(prepend_res==false){



			conformation::Residue const &  current_res(current_pose.residue(reb_res));
			conformation::Residue const &  original_res(original_pose.residue(full_pose_reb_res));

			conformation::Residue const &  previous_res(current_pose.residue(reb_res-1));
			conformation::Residue const &  previous_original_res(original_pose.residue(full_pose_reb_res-1));

			for(Size n=1; n<=4; n++){ //Alpha, Betta, Gramma, Delta of current res
				original_torsions[n] = original_res.mainchain_torsion(n);
				current_torsions[n]= current_res.mainchain_torsion(n);

				Real min_temp= std::min(std::abs(current_torsions[n]-original_torsions[n]+360),std::abs(current_torsions[n]-	original_torsions[n]-360));
				diff_torsions[n]= std::max(0.01,std::min(std::abs(current_torsions[n]-original_torsions[n]), min_temp));
				diff_torsions[0]=diff_torsions[0]+diff_torsions[n];
			}

			for(Size i=1; i<=2; i++){ //epsilon and zetta of previous res.
				Size n=i+4;

				original_torsions[n] = previous_original_res.mainchain_torsion(n);
				current_torsions[n]= previous_res.mainchain_torsion(n);
				Real min_temp= std::min(std::abs(current_torsions[n]-original_torsions[n]+360),std::abs(current_torsions[n]-original_torsions[n]-360));
				diff_torsions[n]= std::max(0.01,std::min(std::abs(current_torsions[n]-original_torsions[n]), min_temp));
				diff_torsions[0]=diff_torsions[0]+diff_torsions[n];
			}


			for(Size i=1; i<=4; i++){ // chi, nu1, nu2 and chi-OH2 of current res
				Size n=i+6;

				original_torsions[n]=original_res.chi(i);//careful i not n in RHS
				current_torsions[n]=current_res.chi(i);//careful i not n in RHS
				Real min_temp= std::min(std::abs(current_torsions[n]-original_torsions[n]+360),std::abs(current_torsions[n]-original_torsions[n]-360));
				diff_torsions[n]= std::max(0.01,std::min(std::abs(current_torsions[n]-original_torsions[n]), min_temp));

				if(i!=4){//don't include X OH.
					diff_torsions[0]=diff_torsions[0]+diff_torsions[n];
				}
			}

		} else {


			//Need to account for the fact that numbering change due to deletion.
			conformation::Residue const &  res_i(current_pose.residue(reb_res));
			conformation::Residue const &  res_iplus1(current_pose.residue(reb_res+1));

			conformation::Residue const &  original_res_i(original_pose.residue(full_pose_reb_res));
			conformation::Residue const &  original_res_iplus1(original_pose.residue(full_pose_reb_res+1));

 			for(Size n=1; n<=3; n++){ //alpha, beta, gamma res i+1 (the res prepended just before this one)
				original_torsions[n] = original_res_iplus1.mainchain_torsion(n);
				current_torsions[n]= res_iplus1.mainchain_torsion(n);

				Real min_temp= std::min(std::abs(current_torsions[n]-original_torsions[n]+360),std::abs(current_torsions[n]-	original_torsions[n]-360));
				diff_torsions[n]= std::max(0.01,std::min(std::abs(current_torsions[n]-original_torsions[n]), min_temp));
				diff_torsions[0]=diff_torsions[0]+diff_torsions[n];
			}


			for(Size n=4; n<=6; n++){ //delta, epsilon and zetta of res i.

				original_torsions[n] = original_res_i.mainchain_torsion(n);
				current_torsions[n]= res_i.mainchain_torsion(n);

				Real min_temp= std::min(std::abs(current_torsions[n]-original_torsions[n]+360),std::abs(current_torsions[n]-	original_torsions[n]-360));
				diff_torsions[n]= std::max(0.01,std::min(std::abs(current_torsions[n]-original_torsions[n]), min_temp));
				diff_torsions[0]=diff_torsions[0]+diff_torsions[n];
			}


			for(Size i=1; i<=4; i++){ // chi, nu1, nu2 and chi-OH2 of res i
				Size n=i+6;
				original_torsions[n]=original_res_i.chi(i);//careful i not n in RHS
				current_torsions[n]=res_i.chi(i);//careful i not n in RHS
				Real min_temp= std::min(std::abs(current_torsions[n]-original_torsions[n]+360),std::abs(current_torsions[n]-original_torsions[n]-360));
				diff_torsions[n]= std::max(0.01,std::min(std::abs(current_torsions[n]-original_torsions[n]), min_temp));

				if(i!=4){//don't include X OH.
					diff_torsions[0]=diff_torsions[0]+diff_torsions[n];
				}
			}


		}

//Print out the torsions of the rotamer
/*
		Size spacing=8;

		std::cout << std::setw(18)<<"original_torsion= ";
		for(Size n=1; n<=10; n++){
		std::cout << std::setw(spacing)<< original_torsions[n] << " ";
		}
		std::cout<< std::endl;

		std::cout << std::setw(18)<<"current_torsion= ";
		for(Size n=1; n<=10; n++){
		std::cout << std::setw(spacing)<< current_torsions[n] << " ";
		}
		std::cout<< std::endl;


		std::cout << std::setw(18)<<"diff_torsions= ";
		for(Size n=1; n<=10; n++){
		std::cout << std::setw(spacing)<< diff_torsions[n] << " ";
		}
		std::cout<< std::endl;

		std::cout << std::setw(18) << "sum_diff= "<< std::setw(spacing) << diff_torsions[0] << std::endl;
*/

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
Calculate_C3_O3_Position(pose::Pose& pose){

		using namespace scoring::constraints;
		using namespace conformation;
		using namespace basic::options;
		using namespace id;
//		using namespace numeric;

//   xyzVector const & a, xyzVector const & b, xyzVector & r )
	Size res_num=24;

	Vector C4_C5;
	Vector O3_C3;
	Vector diff_vector;

	subtract( pose.residue(res_num+1).xyz( 6 ), pose.residue(res_num+1).xyz( 5 ), C4_C5);
	subtract( pose.residue(res_num).xyz( 9 ), pose.residue(res_num).xyz( 8 ), O3_C3);
	subtract( O3_C3, C4_C5, diff_vector);


	Distance const C4_C5_length=	(C4_C5).length();
	Distance const O3_C3_length=	(O3_C3).length();
	Distance const diff_vector_length=	(diff_vector).length();


	std::cout << "03_C3_x= " << O3_C3[0];
	std::cout << " y= " << O3_C3[1];
	std::cout << " z= " << O3_C3[2];
	std::cout << " length= " << O3_C3_length << std::endl;

	std::cout << "C4_C5_x= " << C4_C5[0];
	std::cout << " y= " << C4_C5[1];
	std::cout << " z= " << C4_C5[2];
	std::cout << " length= " << C4_C5_length << std::endl;

	std::cout << "diff_vector_x= " << diff_vector[0];
	std::cout << " y= " << diff_vector[1];
	std::cout << " z= " << diff_vector[2];
	std::cout << " length= " << diff_vector_length << std::endl;

	diff_vector=diff_vector.normalize();
	//diff_vector_length=	;

	std::cout << "norm_diff_vector_x= " << diff_vector[0];
	std::cout << " y= " << diff_vector[1];
	std::cout << " z= " << diff_vector[2];
	std::cout << " length= " << (diff_vector).length() << std::endl;

	Real theta = std::acos(diff_vector[2])*180;
	Real phi = std::atan(diff_vector[1]/diff_vector[0])*180;

	std::cout << "theta= " << theta;
	std::cout << "  phi= " << phi << std::endl;






//	subtract( xyzVector const & a, xyzVector const & b, xyzVector & r )

}

///////////////////////////////////////////////////////////////////////////////////////////////////////
void GetSuiteBin(pose::Pose& pose , utility::vector1<utility::vector1<int> >& suiteBin, Size suite_num){
		using namespace conformation;
		using namespace basic::options;


//	conformation::Residue const & current_res=template_pose.residue(res_num);
//		conformation::Residue const & next_res=template_pose.residue(res_num+1);
		Real bin_size=0.5; //0.5 Amstrong
//		std::cout << "# atoms= " << pose.residue(suite_num).nheavyatoms() << std::endl;
//		Vector nullVector;
//		utility::vector1 <Vector> suiteCoordinates(pose.residue(suite_num).nheavyatoms() /*# suite atom coordinate*/, nullVector);

		std::cout << "atom_index of N1= " << pose.residue(suite_num).atom_index( "N1" ) << std::endl;

		Size n=0;
		for(Size m=1; m<= pose.residue(suite_num).nheavyatoms(); m++){



  		Size res_num;
			if(m <5) {
				res_num=suite_num+1;
			} else {
				res_num=suite_num;
			}

			if((m!=13) && (m!=11) && (m!=6)) continue;
  		n++;

			for(int i=0; i<=2 ; i++){
				suiteBin[n][i+1]=(pose.residue(res_num).xyz( m )[i])/bin_size;
				std::cout << "m= " << m  << " n= " << n << " i=" << i << " suiteBin[n][i]= " << suiteBin[n][i+1];
				std::cout << " pose.residue(res_num).xyz( m )[i] " << pose.residue(res_num).xyz( m )[i] << std::endl;
			}

//		suiteCoordinates[n]=pose.residue(res_num).xyz( n );
//  	std::cout << "atom number= " << n;
//  	std::cout << "atom_type_index of the atom= " << pose.residue(res_num).atom_type_index(n) << std::endl;
  	}
/*
		for(Size n=1; n<= pose.residue(suite_num).nheavyatoms(); n++){
			for(int i=0; i<=2 ; i++){
			suiteBin[n][i+1]=(suiteCoordinates[n][i])/bin_size;
			std::cout << "n= " << n << " i=" << i << " suiteBin[n][i]= " << suiteBin[n][i+1];
			std::cout << " suiteCoordinates[n][i]" << suiteCoordinates[n][i] << std::endl;
			}
		}
*/
}




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* This function serve the purpose of outputting the data of the three original poses (called just pose in function */
void
outfile_original_pose(core::io::silent::SilentFileData& silent_file_data, std::string const& silent_file, pose::Pose original_pose, pose::Pose minimized_original_pose, pose::Pose minimized_then_delete_pose, std::string tag , pose::Pose template_pose ,pose::Pose pose, std::ofstream& outfile, Size current_rebuild_num, core::scoring::ScoreFunctionOP & scorefxn){


			using namespace core::scoring;

			Size reb_res=get_rebuild_pose_res_num(current_rebuild_num);

			std::vector< Real> diff_torsions(11 /* # torsion angles+1*/, 9999.99);

			Calculate_diff_torsions(diff_torsions, pose, template_pose, current_rebuild_num);

			output_data_struct output_data={0,0,0,0,0,0,0,0};
			output_data.current_score=(*scorefxn)(pose); //Need to evaluate pose for energy score to show up in long.txt

/*
			protocols::rna::ensure_phosphate_nomenclature_matches_mini_parin( original_pose);
			protocols::rna::ensure_phosphate_nomenclature_matches_mini_parin( minimized_original_pose);
			protocols::rna::ensure_phosphate_nomenclature_matches_mini_parin( minimized_then_delete_pose);
			protocols::rna::ensure_phosphate_nomenclature_matches_mini_parin( pose);
*/
			Calculate_rmsd(output_data, pose, original_pose, minimized_original_pose, minimized_then_delete_pose, current_rebuild_num, false);


			Output_data_general(silent_file_data, silent_file, output_data, pose, tag, current_rebuild_num, outfile, diff_torsions);



}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool
sort_citeria2(pose_data_struct2 pose_data_1, pose_data_struct2 pose_data_2){

//Sort by score.
return (pose_data_1.score < pose_data_2.score);
}


void
sort_pose_data2(utility::vector1< pose_data_struct2 >&	pose_data_list) {

		//Need to check if this work with vector1, if not switch to std::vector
		sort(pose_data_list.begin(), pose_data_list.end(), sort_citeria2);


}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Best way to make this robust is to tokenize the tag.
std::string
create_tag(std::string prestring, Size group_rotamer, Size subgroup_rotamer, Size fine_rotamer_count, Size current_rebuild_num, std::string old_tag){


			std::string tag=old_tag;
			if(tag=="") tag.append(prestring);
			else tag[0]=prestring[0];

			tag.append("_");
		 	tag.append(lead_zero_string_of(group_rotamer, 3));
		 	tag.append("_");
		 	tag.append(lead_zero_string_of(subgroup_rotamer, 5));
			tag.append("_");
		 	tag.append(lead_zero_string_of(fine_rotamer_count, 4));
		 	tag.append("_");
		 	tag.append(lead_zero_string_of(current_rebuild_num, 2));
//		 	std::cout << "tag= " << tag << std::endl;
		 	return tag;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//name of function no longer correspond to what the  of the function does.
//Right now this function recreate the pose_data_list....which is not necessary...very inefficient...

void
cluster_pose_by_suite_rmsd(utility::vector1< pose_data_struct2 >& pose_data_list, Size& num_pose_kept, Size& current_rebuild_num){

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	utility::vector1<node_state> poses_state(pose_data_list.size(), OPEN);


	Real rmsd_cutoff=option[cluster_rmsd];
	Size num_clustered_pose=0;

	for(Size i=1; i<=pose_data_list.size(); i++){


		if(poses_state[i]==OPEN){
		num_clustered_pose++;
			for(Size j=i+1; j<=pose_data_list.size(); j++){
//				Didn't intend for function to work in this circumstance....but happen to work...
//				std::vector< Real> diff_torsions(11 /* # torsion angles+1*/, 9999.99);
//				Calculate_diff_torsions(diff_torsions, pose_data_list[i].pose, pose_data_list[j].pose, current_rebuild_num);


				Real rmsd=suite_rmsd(pose_data_list[i].pose, pose_data_list[j].pose, current_rebuild_num);

				if(rmsd < rmsd_cutoff) {
						std::cout << "rmsd= " << rmsd;
					  std::cout << "  pose" << pose_data_list[j].tag << " is a neighbor of pose " << pose_data_list[i].tag << std::endl;
					  poses_state[j]=CLOSE;
				}

				if(num_clustered_pose==num_pose_kept) break;
			}
		}
	}


	utility::vector1< pose_data_struct2> New_pose_data_list;
	Size pose_added=0;

	for(Size i=1; i<=pose_data_list.size(); i++) {
		if(poses_state[i]==OPEN){

			New_pose_data_list.push_back(pose_data_list[i]);
			pose_added++;

			if(pose_added==num_pose_kept) break;
	}
		}

	pose_data_list=New_pose_data_list;

}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////



void
Update_pose_data_list(Size& num_pose_kept, std::string& tag, utility::vector1< pose_data_struct2 >& pose_data_list, pose::Pose& current_pose, Real& current_score){

						//The order of evaluation of the two expression in the if statement is important!
						if(pose_data_list.size() < num_pose_kept || current_score < pose_data_list[num_pose_kept].score) {

							std::cout << "tag= " << tag;

							if(pose_data_list.size() >= num_pose_kept){
								std::cout << " cutoff score= " << pose_data_list[num_pose_kept].score;
							}
							std::cout<< " score= " << current_score;

								pose_data_struct2 current_pose_data;
								current_pose_data.pose=current_pose;
								current_pose_data.score = current_score;
//								current_pose_data.group_rotamer=group_rotamer;
//								current_pose_data.subgroup_rotamer=subgroup_rotamer;
								current_pose_data.tag=tag;
								pose_data_list.push_back(current_pose_data);
							std::cout << " pose_data_list.size= " << pose_data_list.size() << std::endl;
						}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
Cluster_poses_test(Size& num_pose_kept, utility::vector1< pose_data_struct2 >& pose_data_list, Size& current_rebuild_num){
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		sort_pose_data2(pose_data_list);
		if(option[apply_clustering]) cluster_pose_by_suite_rmsd(pose_data_list, num_pose_kept, current_rebuild_num);

		if(pose_data_list.size()>num_pose_kept){
				pose_data_list.erase(pose_data_list.begin()+num_pose_kept, pose_data_list.end());
		}

		std::cout<< "after erasing.. pose_data_list= " << pose_data_list.size() << std::endl;

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
//This function is not up-to-date
void
Fine_sampling(utility::vector1< pose_data_struct2 >& pose_data_list, Size& current_rebuild_num , Size& num_residue_reb, core::scoring::ScoreFunctionOP & scorefxn, core::scoring::ScoreFunctionOP & atr_rep_screening_scorefxn, pose::Pose& template_pose, Size& reb_res, std::string& silent_file, core::io::silent::SilentFileData& silent_file_data, std::ofstream& outfile, std::vector< pose::Pose>& original_poses, std::vector< pose::Pose>& minimized_original_poses, std::vector< pose::Pose>& minimized_then_delete_poses, Real& base_atr_score, Real& base_rep_score, Size num_pose_kept, Size multiplier, Real rep_cutoff){

	using namespace core::id;
	using namespace core::scoring;
	using namespace core::scoring::rna;



	utility::vector1< pose_data_struct2 > old_pose_data_list = pose_data_list;

	for(Size i=1; i<= pose_data_list.size(); i++){
		pose_data_list[i].tag.append("_");
		pose_data_list[i].tag.append(lead_zero_string_of(0, 3));
	}

	std::vector< Real> diff_torsions(11 , 9999.99);
	scoring::EMapVector	energy_map;

	std::cout << "residue rebuilt number= " << reb_res << std::endl;
	Size tot_rotamer_count=0;
 	Size good_rep_rotamer_count=0;
 	Size good_atr_rotamer_count=0;
 	Size both_count=0;

	for(Size i=1; i<=old_pose_data_list.size(); i++){

		pose::Pose current_pose=old_pose_data_list[i].pose;
		std::string	old_tag=old_pose_data_list[i].tag;
		old_tag[0]='F';

		pose::Pose current_pose_screen=current_pose;

 		Real current_score=pose_data_list[i].score;
 		Size group_rotamer=pose_data_list[i].group_rotamer;
 		Size subgroup_rotamer=pose_data_list[i].subgroup_rotamer;

		conformation::Residue const & current_res=current_pose.residue(reb_res);
		conformation::Residue const & next_res=current_pose.residue(reb_res+1);

		utility::vector1< Real > backbone_rotamer(13 , 9999.99);
		utility::vector1< Real > new_backbone_rotamer(13 , 9999.99);

		for(Size n=1; n<=1; n++){ //delta1
				Size i=n+3;
				backbone_rotamer[n]=current_res.mainchain_torsion(i);
		}

		for(Size n=2; n<=4; n++){ // chi_1, nu2_1, nu1_1 (recall that chi_1 2OH is determined seperately)
				Size i=n-1;
				backbone_rotamer[n]=current_res.chi(i);
		}

		for(Size n=5; n<=6; n++){ //epsilon1, zeta1
				Size i=n;
				backbone_rotamer[n]=current_res.mainchain_torsion(i);
		}


		for(Size n=7; n<=10; n++){ //alpha2, beta2, gamma2, delta2
				Size i=n-6;
				backbone_rotamer[n]=next_res.mainchain_torsion(i);
		}

		//Actually don't need this for prepend.
		for(Size n=11; n<=13; n++){ // chi_2, nu2_2, nu1_2
				Size i=n-10;
				backbone_rotamer[n]=next_res.chi(i);
		}
		int count=0;
		int min=-1;
		int max=1;
		Real delta_torsion=5;

		for(int a1 = min; a1 <= max; a1++ ) {
		for(int b1 = min; b1 <= max; b1++ ) {
		for(int g1 = min; g1 <= max; g1++ ) {
		for(int e1 = min; e1 <= max; e1++ ) {
		for(int z1 = min; z1 <= max; z1++ ) {
		for(int chi1 = min; chi1 <= max; chi1++ ) {
			if(a1==0 && b1==0 && g1==0 && e1==0 && z1==0 && chi1==0) continue;
			count++;

			//Work only for prepend right now.
			new_backbone_rotamer[1]= backbone_rotamer[1]; //delta i
			new_backbone_rotamer[2]= backbone_rotamer[2] + (chi1*delta_torsion);
			new_backbone_rotamer[3]= backbone_rotamer[3]; //nu2 i
			new_backbone_rotamer[4]= backbone_rotamer[4]; // nu1 i
			new_backbone_rotamer[5]= backbone_rotamer[5] +(e1*delta_torsion);
			new_backbone_rotamer[6]= backbone_rotamer[6] +(z1*delta_torsion);
			new_backbone_rotamer[7]= backbone_rotamer[7] +(a1*delta_torsion);
			new_backbone_rotamer[8]= backbone_rotamer[8] +(b1*delta_torsion);
			new_backbone_rotamer[9]= backbone_rotamer[9] +(g1*delta_torsion);

			apply_rotamer( current_pose_screen, reb_res, new_backbone_rotamer);


			(*atr_rep_screening_scorefxn)(current_pose_screen);
			energy_map=current_pose_screen.energies().total_energies();
			Real atr_score=atr_rep_screening_scorefxn->get_weight(fa_atr)*energy_map[scoring::fa_atr];
			Real rep_score=atr_rep_screening_scorefxn->get_weight(fa_rep)*energy_map[scoring::fa_rep];
			tot_rotamer_count++;

			Real delta_rep_score=rep_score-base_rep_score;
			Real delta_atr_score=atr_score-base_atr_score;

			if(delta_rep_score<rep_cutoff){
					good_rep_rotamer_count++;
			}

			if(delta_atr_score<-1){
					good_atr_rotamer_count++;
			}

			if((delta_rep_score<rep_cutoff) && (delta_atr_score<-1) && (delta_rep_score+delta_atr_score<0)){


					both_count++;
//					std::cout << " rep= " << rep_score-base_rep_score << " atr= " << atr_score-base_atr_score;
//					std::cout << " rep_n= " << good_rep_rotamer_count << " atr_n= " << good_atr_rotamer_count;
//					std::cout << " both= " << both_count << " tot= " << tot_rotamer_count << std::endl;

					apply_rotamer( current_pose, reb_res, new_backbone_rotamer);


					std::string tag=old_tag;
					tag.append("_");
					tag.append(lead_zero_string_of(count, 3));

					current_score=(*scorefxn)(current_pose);

					Calculate_diff_torsions(diff_torsions, current_pose, template_pose, current_rebuild_num);

    			Real rmsd = suite_rmsd(original_poses[current_rebuild_num], current_pose, reb_res);
   			 	Real rmsd_wrt_minimized = suite_rmsd(minimized_original_poses[current_rebuild_num], current_pose, reb_res);
					Real rmsd_wrt_correct_minimized = suite_rmsd(minimized_then_delete_poses[current_rebuild_num], current_pose, reb_res);
					if((rmsd<1.5) || (rmsd_wrt_minimized<1.5) || (rmsd_wrt_correct_minimized<1.5)){
//						Output_data_general(silent_file_data, silent_file, rmsd, rmsd_wrt_minimized, rmsd_wrt_correct_minimized, current_score, current_pose, tag, current_rebuild_num, outfile, diff_torsions);
					}

					Update_pose_data_list(num_pose_kept, tag, pose_data_list, current_pose, current_score);
					if((pose_data_list.size()==num_pose_kept*multiplier)) {
						Cluster_poses_test(num_pose_kept, pose_data_list, reb_res);
					}
			}

		}
		}
		}
		}
		}
		}
	}
	Cluster_poses_test(num_pose_kept, pose_data_list, reb_res);

}
*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
Remove_O2Star_hydrogen_variant_type(pose::Pose & pose){
	for(Size i=1; i<=pose.total_residue(); i++){
		core::pose::remove_variant_type_from_pose_residue( pose, "VIRTUAL_O2STAR_HYDROGEN", i);
	}
}

void
Delete_O2Star_hydrogen(pose::Pose & pose){

	for(Size i=1; i<=pose.total_residue(); i++){
		core::pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_O2STAR_HYDROGEN", i);
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
Setup_chain_break_scoring_for_minimize_original_pose(pose::Pose& current_pose){
 		std::cout << "Setup harmonic_chainbreak scoring " << std::endl;


	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace core::kinematics;
	using namespace core::optimization;
	using namespace core::io::silent;
	using namespace core::id;
	using namespace protocols::rna;

	using numeric::conversions::radians;



//    From RAD.param file
//	  ICOOR_INTERNAL  UPPER -175.907669   60.206192    1.607146   O3*   C3*   C4*   , Upper is P1
//    ICOOR_INTERNAL  LOWER  -64.027359   71.027062    1.593103   P     O5*   C5*   , Lower is O3'
//    Bug that bond distance is not the same. Rhiju suggest using 1.593103

//Original (Amber?) parameter 1.608, 119.8, 103.4

		Real const O3_P_distance( 1.593 ); //amber=1.608
		Real const O3_angle( 119.8 ); // 180-60.206192
		Real const  P_angle( 108.97 ); // Quite off from original (Amber?) ,180-71.027062

		Real const distance_stddev( 0.0659 ); // amber is 0.0659
		Real const angle_stddev_degrees_P( 8.54 ); // amber is 8.54 (P angle), 5.73 (O3 angle)
		Real const angle_stddev_degrees_O3( 5.73 );



		ConstraintSetOP cst_set( current_pose.constraint_set()->clone() );
		assert( cst_set ); //if ( !cst_set ) cst_set = new ConstraintSet();



		FuncOP const distance_func( new HarmonicFunc( O3_P_distance, distance_stddev ) );
		FuncOP const O3_angle_func( new HarmonicFunc( radians( O3_angle ), radians( angle_stddev_degrees_P ) ) );
		FuncOP const  P_angle_func( new HarmonicFunc( radians(  P_angle ), radians( angle_stddev_degrees_O3 ) ) );

//		bool const prepend_res = 	Is_prepend(current_rebuild_num);
//		Size full_pose_reb_res=get_full_pose_reb_res(current_rebuild_num);
//		Size reb_res=Convert_to_partial_pose_seq_num(full_pose_reb_res, current_rebuild_num);


		utility::vector1 <std::pair<Size, Size> > jump_points_list;
		utility::vector1 <Size> cut_points_list;

		utility::vector1 < utility::vector1 <Residue_info_struct> > rebuild_residue_group_list;
		create_rebuild_residue_groups(rebuild_residue_group_list); //A rebuild_residue_group is a list of Residue_info_struct sorted by seq_num. There is one group for each chain break.

		get_jump_points(jump_points_list, rebuild_residue_group_list); // One jump_point_pair per chain break, ordered by seq_num
		get_cut_points(cut_points_list, jump_points_list); // One cut_point per chain break, ordered by seq_num

		for(Size cut_point_num=1; cut_point_num<=cut_points_list.size(); cut_point_num++){

			Size five_prime_res=cut_points_list[cut_point_num];
			Size three_prime_res=cut_points_list[cut_point_num]+1;

			Residue const & rsd1( current_pose.residue(five_prime_res) );
			Residue const & rsd2( current_pose.residue(three_prime_res) );

			AtomID const C3_id( rsd1.atom_index( "C3*" ), five_prime_res);
			AtomID const O3_id( rsd1.atom_index( "O3*" ), five_prime_res);
			AtomID const  P_id( rsd2.atom_index( "P"   ), three_prime_res);
			AtomID const O5_id( rsd2.atom_index( "O5*" ), three_prime_res);

			// distance from O3* to P
			cst_set->add_constraint( new AtomPairConstraint( O3_id, P_id, distance_func ) );

			// angle at O3*
			cst_set->add_constraint( new AngleConstraint( C3_id, O3_id, P_id, O3_angle_func ) );

			// angle at P
			cst_set->add_constraint( new AngleConstraint( O3_id, P_id, O5_id,  P_angle_func ) );
		}

		current_pose.constraint_set( cst_set );

//		std::cout << "minimized_pose_cst_set" << std::endl;
//		cst_set->show(std::cout);

//		cst_set->remove_constraint( new AngleConstraint( O3_id, P_id, O5_id,  P_angle_func );
//		cst_set->remove_constraint( new AngleConstraint( C3_id, O3_id, P_id, O3_angle_func );
//		cst_set->remove_constraint( new AngleConstraint( C3_id, O3_id, P_id, O3_angle_func );

//		std::cout << "minimized_pose_cst_set" << std::endl;
//		cst_set->show(std::cout);


//		dump_pdb( current_pose, "after_setting_constraint.pdb" );

//		ConstraintSetOP cst_set( original_pose.constraint_set()->clone() );
//		assert( cst_set ); //if ( !cst_set ) cst_set = new ConstraintSet();
//		minimized_original_pose.constraint_set( cst_set );


}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//This assume that they is only 1 chain-break, will need to update the if statement once we allow multiple chain break.
void
Setup_chain_break_scoring(pose::Pose& current_pose, Size& num_residue_reb, Size& current_rebuild_num){
 		std::cout << "Setup harmonic_chainbreak scoring " << std::endl;


	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace core::kinematics;
	using namespace core::optimization;
	using namespace core::io::silent;
	using namespace core::id;
	using namespace protocols::rna;

	using numeric::conversions::radians;

//    From RAD.param file
//	  ICOOR_INTERNAL  UPPER -175.907669   60.206192    1.607146   O3*   C3*   C4*   , Upper is P1
//    ICOOR_INTERNAL  LOWER  -64.027359   71.027062    1.593103   P     O5*   C5*   , Lower is O3'
//    Bug that bond distance is not the same. Rhiju suggest using 1.593103

//Original (Amber?) parameter 1.608, 119.8, 103.4

		Real const O3_P_distance( 1.593 ); //amber=1.608
		Real const O3_angle( 119.8 ); // 180-60.206192
		Real const  P_angle( 108.97 ); // Quite off from original (Amber?) ,180-71.027062

		Real const distance_stddev( 0.0659 ); // amber is 0.0659
		Real const angle_stddev_degrees_P( 8.54 ); // amber is 8.54 (P angle), 5.73 (O3 angle)
		Real const angle_stddev_degrees_O3( 5.73 );



		ConstraintSetOP cst_set( current_pose.constraint_set()->clone() );
		assert( cst_set ); //if ( !cst_set ) cst_set = new ConstraintSet();



		FuncOP const distance_func( new HarmonicFunc( O3_P_distance, distance_stddev ) );
		FuncOP const O3_angle_func( new HarmonicFunc( radians( O3_angle ), radians( angle_stddev_degrees_P ) ) );
		FuncOP const  P_angle_func( new HarmonicFunc( radians(  P_angle ), radians( angle_stddev_degrees_O3 ) ) );

		bool const prepend_res = 	Is_prepend(current_rebuild_num);
		Size full_pose_reb_res=get_full_pose_reb_res(current_rebuild_num);
		Size reb_res=Convert_to_partial_pose_seq_num(full_pose_reb_res, current_rebuild_num);

		Size five_prime_res, three_prime_res;
		if(prepend_res){
			five_prime_res=reb_res-1;
			three_prime_res=reb_res;
		}else {
			five_prime_res=reb_res;
			three_prime_res=reb_res+1;
		}


		Residue const & rsd1( current_pose.residue(five_prime_res) );
		Residue const & rsd2( current_pose.residue(three_prime_res) );

		AtomID const C3_id( rsd1.atom_index( "C3*" ), five_prime_res);
		AtomID const O3_id( rsd1.atom_index( "O3*" ), five_prime_res);
		AtomID const  P_id( rsd2.atom_index( "P"   ), three_prime_res);
		AtomID const O5_id( rsd2.atom_index( "O5*" ), three_prime_res);



		// distance from O3* to P
		cst_set->add_constraint( new AtomPairConstraint( O3_id, P_id, distance_func ) );

		// angle at O3*
		cst_set->add_constraint( new AngleConstraint( C3_id, O3_id, P_id, O3_angle_func ) );

		// angle at P
		cst_set->add_constraint( new AngleConstraint( O3_id, P_id, O5_id,  P_angle_func ) );

		current_pose.constraint_set( cst_set );


//		dump_pdb( current_pose, "after_setting_constraint.pdb" );

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Satisfy_Chain_Closure(pose::Pose & pose){

//}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
Calculate_bond_angles(pose::Pose & pose, Size reb_res){

			using namespace scoring::constraints;
			using namespace conformation;
			using namespace basic::options;
			using namespace id;
			using numeric::conversions::degrees;

			Real const Ideal_O3_P_distance( 1.608 );
			Real const Ideal_O3_angle( 119.8 );
			Real const Ideal_P_angle( 103.4 );

			reb_res=19;

//			for(; reb_res<=24; reb_res++){

				Real O3_P_distance = (pose.residue(reb_res).xyz( 1 ) - pose.residue(reb_res-1).xyz(9) ).length();

				Vector C3_O3;
				Vector O3_P;
				Vector P_O5;

				subtract( pose.residue(reb_res-1).xyz( 8 ), pose.residue(reb_res-1).xyz( 9 ), C3_O3);
				subtract( pose.residue(reb_res-1).xyz( 9 ), pose.residue(reb_res).xyz( 1 ), O3_P);
				subtract( pose.residue(reb_res).xyz( 1 ), pose.residue(reb_res).xyz( 4 ), P_O5);

				Real O3_angle = degrees(angle_of( C3_O3, O3_P ));
				Real P_angle = degrees(angle_of( O3_P, P_O5 ));



				std::cout << "lower residue= " << (reb_res-1);
				std::cout << " P_O3_distance= " << std::setw(12) << std::setprecision(8) << O3_P_distance;
				std::cout << " C3_O3_P_angle= " << std::setw(12) << std::setprecision(8) << O3_angle;
				std::cout << " O3_P_O5_angle= " << std::setw(12) << std::setprecision(8) << P_angle << std::endl;
//			}

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
Import_base_pose(std::string base_pose_name, pose::Pose& base_pose, pose::Pose& template_pose){

		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace core::chemical;
		using namespace core::conformation;

		ResidueTypeSetCAP rsd_set;
		rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "rna" );

		base_pose_name.append(".pdb");

  	std::cout << "The following pose will be imported and use as the base pose :";
 	 	std::cout << base_pose_name << std::endl;
		core::import_pose::pose_from_pdb( base_pose, *rsd_set, base_pose_name);

    //////////////////////////Since base_pose is imported...need to reinitialize the fold tree////////

		std::cout << "Reinitialize fold_tree of imported base_pose with the fold_tree of the partially_rebuild_original_pose: " << std::endl;
		kinematics::FoldTree template_fold_tree=template_pose.fold_tree();
		base_pose.fold_tree( template_fold_tree );
		Output_fold_tree_info(base_pose, "Imported base_pose");
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Update June 14, still need to verify if the function works properly
void
Copy_torsion_angle(pose::Pose& pose, pose::Pose& template_pose, Size reb_res, Size current_rebuild_num){

 		using namespace core::chemical;
		using namespace core::conformation;
	  using namespace core::id;
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		bool const prepend_res = 	Is_prepend(current_rebuild_num);

//		std::cout << "In Copy_torsion_angle function " << std::endl;

			ResidueTypeSetCAP rsd_set; //This line is repeated
			rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "rna" ); //This line is repeated

			Size five_prime_res, three_prime_res;
			if(prepend_res){
				five_prime_res=reb_res-1;
				three_prime_res=reb_res;
			}else {
				five_prime_res=reb_res;
				three_prime_res=reb_res+1;
			}

			conformation::Residue const & lower_res=template_pose.residue(five_prime_res);
			conformation::Residue const & upper_res=template_pose.residue(three_prime_res);


			//Even through there is the chain_break, these torsions should be defined due to the existence of the upper and lower variant type atoms.

			for(Size n=1; n<=3; n++){ //alpha, beta, gamma of upper residue
				pose.set_torsion( TorsionID( three_prime_res, id::BB,  n ), upper_res.mainchain_torsion(n) );
/*
				AtomID id1,id2,id3,id4;
				pose.conformation().get_torsion_angle_atom_ids( TorsionID( three_prime_res, id::BB,  n ), id1, id2, id3, id4 );
				std::cout << "Torsion= " << n;
				std::cout << " atom1= " << upper_res.atom_name(id1.atomno());
				std::cout << " atom2= " << upper_res.atom_name(id2.atomno());
				std::cout << " atom3= " << upper_res.atom_name(id3.atomno());
				std::cout << " atom4= " << upper_res.atom_name(id4.atomno()) << std::endl;
*/

			}


			for(Size n=5; n<=6; n++){ //epsilon and zeta of lower residue
				pose.set_torsion( TorsionID( five_prime_res, id::BB,  n ), lower_res.mainchain_torsion(n) );
/*
				AtomID id1,id2,id3,id4;
				pose.conformation().get_torsion_angle_atom_ids( TorsionID( five_prime_res, id::BB,  n ), id1, id2, id3, id4 );
				std::cout << "Torsion= " << n;
				std::cout << " atom1= " << lower_res.atom_name(id1.atomno());
				std::cout << " atom2= " << lower_res.atom_name(id2.atomno());
				std::cout << " atom3= " << lower_res.atom_name(id3.atomno());
				std::cout << " atom4= " << lower_res.atom_name(id4.atomno()) << std::endl;
*/

			}


}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
Atr_rep_screening(pose::Pose& current_pose_screen, Size& reb_res, utility::vector1 <Real>& current_rotamer, count_struct& count_data, core::scoring::ScoreFunctionOP & atr_rep_screening_scorefxn, Real& base_rep_score, Real&  base_atr_score, std::string const & rotamer_set_option){

   			using namespace basic::options;
				using namespace basic::options::OptionKeys;
				using namespace core::scoring;

				Real rep_cutoff=option[fa_rep_cutoff];

				(*atr_rep_screening_scorefxn)(current_pose_screen);
				scoring::EMapVector energy_map=current_pose_screen.energies().total_energies();
				Real atr_score=atr_rep_screening_scorefxn->get_weight(fa_atr)*energy_map[scoring::fa_atr];
				Real rep_score=atr_rep_screening_scorefxn->get_weight(fa_rep)*energy_map[scoring::fa_rep];

				Real delta_rep_score=rep_score-base_rep_score;
				Real delta_atr_score=atr_score-base_atr_score;

				if(delta_rep_score<rep_cutoff) count_data.good_rep_rotamer_count++;

				if(delta_atr_score<(-1)) count_data.good_atr_rotamer_count++;

				//For both screening below took out condition: "(delta_rep_score<(rep_cutoff*option[Atr_rep_reweight_scaling]))" on May 17 since this condition didn't seem to add to screening.

				if( (delta_atr_score<(-1)) && ((delta_rep_score+delta_atr_score) < 0) ) {
						count_data.both_count++;
					  if(option[Verbose]){
							if(rotamer_set_option=="near_native_torsions") std::cout << "near_native_torsions_count= " << count_data.Near_Native_Rotamer_Count;
							std::cout << " rep= " << delta_rep_score << " atr= " << delta_atr_score;
							std::cout << " rep_n= " << count_data.good_rep_rotamer_count;
							std::cout << " atr_n= " << count_data.good_atr_rotamer_count;
							std::cout << " both= " << count_data.both_count << " tot= " << count_data.tot_rotamer_count << std::endl;
						}
						return true;
				} else {
						return false;
				}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//June 14, 2009 updated..still need to verify that function is working properly
bool
Chain_break_screening(pose::Pose& current_pose_screen, pose::Pose& current_pose_chain_break, Size& reb_res, utility::vector1 <Real>& current_rotamer, count_struct& count_data, core::scoring::ScoreFunctionOP & atr_rep_screening_scorefxn,  core::scoring::ScoreFunctionOP & constraint_scorefxn, protocols::rna::RNA_LoopCloser& rna_loop_closer, Real& base_rep_score, Size current_rebuild_num){

						using namespace basic::options;
						using namespace basic::options::OptionKeys;
						using namespace core::scoring;

						bool const prepend_res = 	Is_prepend(current_rebuild_num);

						apply_rotamer(current_pose_chain_break, reb_res, current_rebuild_num,  current_rotamer);

						Size five_prime_res, three_prime_res;
						if(prepend_res){
							five_prime_res=reb_res-1;
							three_prime_res=reb_res;
						}else {
							five_prime_res=reb_res;
							three_prime_res=reb_res+1;
						}

						//Distance between C5' of three_prime and O3' of three_prime residue at chain_break.
						Real const C5_O3_distance = (current_pose_chain_break.residue(three_prime_res).xyz("C5*") - current_pose_chain_break.residue(five_prime_res).xyz("O3*") ).length();

//						Calculate_bond_angles(current_pose_chain_break, reb_res);

							//4.4682=1.593+1.593+1.441 Angstrom which are the O3-P, P-O5, O5-C5 bond distances respectively
							//Where did I get the 4.627???,
							if(C5_O3_distance>4.627 || C5_O3_distance<2.0) return false; //basically cannot close chain if the C5_O3_distance is either too short or too long.
							count_data.C5_O3_distance_count++;

							///////////////Additional fa_rep screening//////////////////////
							(*atr_rep_screening_scorefxn)(current_pose_screen);
							scoring::EMapVector energy_map=current_pose_screen.energies().total_energies();
							Real rep_score=atr_rep_screening_scorefxn->get_weight(fa_rep)*energy_map[scoring::fa_rep];

							Real delta_rep_score=rep_score-base_rep_score;

							if(delta_rep_score<10) { //Very lenient rep_score screening here, no screen on atr_score.
								count_data.good_rep_rotamer_count++;
							} else {
								return false;
							}

//							(*constraint_scorefxn)(current_pose_chain_break);
//							std::cout << "Before CCD= " << std::endl;
//							constraint_scorefxn->show( std::cout, current_pose_chain_break);

							Real mean_dist_err=rna_loop_closer.apply( current_pose_chain_break, five_prime_res );
//							std::cout << "mean_dist_err= "<< mean_dist_err <<  std::endl;

							(*constraint_scorefxn)(current_pose_chain_break);
//							std::cout << "After CCD= " << std::endl;
//							constraint_scorefxn->show( std::cout, current_pose_chain_break);
							energy_map=current_pose_chain_break.energies().total_energies();
							Real angle_score = energy_map[scoring::angle_constraint];
							Real distance_score = energy_map[scoring::atom_pair_constraint];


//							temp_vector[0]=mean_dist_err;
//							temp_vector[1]=angle_score;
//							temp_vector[2]=distance_score;
//							temp_vector[3]=C5_O3_distance;

							if(angle_score<5) count_data.good_angle_count++;
							if(distance_score<5) count_data.good_distance_count++;
							if((angle_score<5) && (distance_score<5)){
								count_data.both_count++;
								if(option[Verbose]){
									std::cout << " rep= " << delta_rep_score << " rep_n= " << count_data.good_rep_rotamer_count;
									std::cout << "C5_O3= " << C5_O3_distance << " C5_O3_n= " << count_data.C5_O3_distance_count;
									std::cout << "  angle= " << angle_score << " dist= " << distance_score;
									std::cout << " angle_n= " << count_data.good_angle_count;
									std::cout << " dist_n= " << count_data.good_distance_count;
									std::cout << " both= " << count_data.both_count;
									std::cout << " tot= " << count_data.tot_rotamer_count << std::endl;
								}
								return true;
							} else {
								return false;
							}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Update June 14, 2009, still need to verify if function work properly
void
Pose_selection_by_full_score(utility::vector1< pose_data_struct2 >& pose_data_list, pose::Pose& current_pose, Size& reb_res, utility::vector1<Real>& current_rotamer, pose::Pose& template_pose, Size& current_rebuild_num, std::vector< pose::Pose>& original_poses, std::vector< pose::Pose>& minimized_original_poses, std::vector< pose::Pose>& minimized_then_delete_poses, Size& num_pose_kept, Size& multiplier, core::io::silent::SilentFileData& silent_file_data, std::string& silent_file, Size& group_rotamer, Size& subgroup_rotamer, std::string& old_tag, core::scoring::ScoreFunctionOP & scorefxn, Size & num_residue_reb, std::ofstream& outfile, fine_rotamer_ID_struct fine_rotamer_ID, count_struct& count_data){

			using namespace basic::options;
		  using namespace basic::options::OptionKeys;
			using namespace core::scoring;


			apply_rotamer(current_pose, reb_res, current_rebuild_num, current_rotamer);

			std::string tag = create_tag("U", group_rotamer, subgroup_rotamer, fine_rotamer_ID.count, current_rebuild_num , old_tag);

//		scorefxn->show( std::cout, current_pose );

			std::vector< Real> diff_torsions(11 /* # torsion angles+1*/, 9999.99);
			Calculate_diff_torsions(diff_torsions, current_pose, template_pose, current_rebuild_num);


//		dump_pdb( current_pose, tag+".pdb" );

			output_data_struct output_data={0,0,0,0,0,0,0,0};
			Calculate_rmsd(output_data, current_pose, original_poses[current_rebuild_num], minimized_original_poses[current_rebuild_num], minimized_then_delete_poses[current_rebuild_num], current_rebuild_num, false);

			output_data.current_score=(*scorefxn)(current_pose);

			if(output_data.current_score>-1) return; //Very bad score pose...don't bother to continue.

			if(option[Verbose]){
				std::cout << tag <<  std::endl;
	  		Output_data_general(silent_file_data, silent_file, output_data, current_pose, tag, current_rebuild_num, outfile, diff_torsions);
			}



			Update_pose_data_list(num_pose_kept, tag, pose_data_list, current_pose, output_data.current_score);

			if((pose_data_list.size()==num_pose_kept*multiplier)) {
					Cluster_poses_test(num_pose_kept, pose_data_list, current_rebuild_num);
			}

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool
Check_correct_contact(pose::Pose& pose, Size const res_shift){

	using namespace scoring::constraints;
	using namespace conformation;
	using namespace basic::options;
	using namespace id;

	//Valid for res 4 only
	//distance between 1H4

	Real distance=(pose.residue(21-res_shift).xyz(28)-pose.residue(15).xyz(12)).length();

//	std::cout << "distance= " << distance << std::endl;

	//3=2(ideal H-bond distance)+2(error from 20degree sampling, non-ideal distance)
	if(distance>4.0) return false;
	else return true;

}


struct atom{
	Size res_num;
	Size atom_num;
};


struct atom_pair{
	atom atom1;
	atom atom2;
};

//Add_atom_pair_constraint(std::vector<atom_pair>& atom_pair_constraint_list)
/*
void
Setup_atom_pair_constraint(std::vector<atom_pair>& atom_pair_constraint_list){

//	Add_atom_pair_constraint(atom_pair_constraint_list, ")

}
*/



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
Setup_rebuild_pose_virtual_phosphate(pose::Pose & pose, Size current_rebuild_num, bool verbose){

		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace core::chemical;
		using namespace core::conformation;


		if(verbose)	std::cout << "enter Setup_rebuild_pose_virtual_phosphate" << std::endl;

		Size num_residue_reb = RunTimeParameters::rebuild_residue_list.size();

		//Add virtual phosphate to every loop 3' end
		utility::vector1 < utility::vector1 <Residue_info_struct> > residue_group_list;
		utility::vector1 <Residue_info_struct> unrebuilded_residue_list;
		for(Size rebuild_num=current_rebuild_num+1; rebuild_num<=num_residue_reb; rebuild_num++){
 				unrebuilded_residue_list.push_back(RunTimeParameters::rebuild_residue_list[rebuild_num]);
		}

		if(verbose)	{
			for(Size j=1; j<=unrebuilded_residue_list.size(); j++){
				std::cout << Get_one_letter_name(unrebuilded_residue_list[j].name) << unrebuilded_residue_list[j].seq_num <<  std::endl;
			}
		}

		create_residue_group_list(residue_group_list, unrebuilded_residue_list);
		if(verbose)	std::cout << "residue_group_list.size()=" << residue_group_list.size() << std::endl;

		for(Size loop_num=1; loop_num<=residue_group_list.size(); loop_num++){

			Size full_pose_three_prime_end_res=residue_group_list[loop_num][residue_group_list[loop_num].size()].seq_num+1;
			if(verbose) std::cout << "full_pose_three_prime_end_res= " << full_pose_three_prime_end_res << std::endl;
			Size three_prime_end_res=Convert_to_partial_pose_seq_num(full_pose_three_prime_end_res, current_rebuild_num);
			if(verbose)	std::cout << "three_prime_end_res= " << three_prime_end_res << std::endl;
			dump_pdb(pose , "rescore_pose_before_Virtual_posphate.pdb");
			core::pose::add_variant_type_to_pose_residue(pose , "VIRTUAL_PHOSPHATE", three_prime_end_res);
			dump_pdb(pose , "rescore_pose_after_Virtual_posphate.pdb");

		}
//		std::cout << "exit Setup_rebuild_pose_virtual_phosphate" << std::endl;

}


void
output_data_rescore(core::io::silent::SilentFileData& silent_file_data, std::string const & silent_file, pose_data_struct& pose_data, std::ofstream& outfile, Size const current_rebuild_num){

	using namespace core::io::silent;
  using namespace core::scoring;
  using namespace basic::options;
	using namespace basic::options::OptionKeys;

	Size num_residue_reb = RunTimeParameters::rebuild_residue_list.size();

	Size spacing=8;
  Size tag_spacing;
	Size char_per_line=207;
  if(pose_data.tag.length()<char_per_line){
		tag_spacing=char_per_line;
 	} else {
 		tag_spacing=2*char_per_line;
 	}

	Size const write_score_only = option[print_score_only] ;

	outfile << std::setw(tag_spacing) << std::left << pose_data.tag;
	outfile << std::setw(spacing+2) << std::fixed << std::setprecision(3)  << std::left << pose_data.rmsd;
	outfile << std::setw(spacing+2) << std::fixed << std::setprecision(3)  << std::left << pose_data.rmsd_min;
	outfile << std::setw(spacing+2) << std::fixed << std::setprecision(3)  << std::left << pose_data.rmsd_cor;
 	outfile << std::setw(spacing+2) << std::fixed << std::setprecision(2)  << std::left << pose_data.score;
 	outfile << std::setw(spacing+2) << std::fixed << std::setprecision(3)  << std::left << pose_data.loop_rmsd;
	outfile << std::setw(spacing+2) << std::fixed << std::setprecision(3)  << std::left << pose_data.loop_rmsd_min;
	outfile << std::setw(spacing+2) << std::fixed << std::setprecision(3)  << std::left << pose_data.loop_rmsd_cor;
 	outfile << std::setw(spacing+2) << std::fixed << std::setprecision(2)  << std::left << pose_data.dist_to_close_chain;
	outfile << std::setw(spacing+2) << std::fixed << std::setprecision(1)  << std::left << pose_data.diff_torsions;
	outfile << std::setw(spacing+2) << std::fixed << std::setprecision(2)  << std::left << pose_data.new_score;
	outfile << std::setw(spacing+2) << std::fixed << std::setprecision(2)  << std::left << pose_data.norm_new_score;

	if(current_rebuild_num==num_residue_reb){
	 		outfile << std::setw(spacing+2) << std::fixed << std::setprecision(3)  << std::left << pose_data.loop_rmsd_exclude_bulge;
			outfile << std::setw(spacing+2) << std::fixed << std::setprecision(3)  << std::left << pose_data.loop_rmsd_cor_exclude_bulge;
	}
	outfile << "\n";

	RNA_SilentStruct s((*pose_data.pose_OP), pose_data.tag);


	std::string name;
	name.append( "rmsd" );
	s.add_energy( name, pose_data.rmsd );

	name.clear();
	name.append( "rmsd_cor" );
	s.add_energy( name, pose_data.rmsd_cor);

	name.clear();
	name.append( "loop_rmsd;" );
	s.add_energy( name, pose_data.loop_rmsd);

	name.clear();
	name.append( "loop_rmsd_cor" );
	s.add_energy( name, pose_data.loop_rmsd_cor);

	name.clear();
	name.append( "norm_new_score" );
	s.add_energy( name, pose_data.norm_new_score);

	name.clear();
	name.append( "dist_to_close_chain" );
	s.add_energy( name, pose_data.dist_to_close_chain );

	silent_file_data.write_silent_struct(s, silent_file, write_score_only);

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


Real
Rescore_pose_list(
std::vector< pose::Pose>& original_poses,
std::vector< pose::Pose>& minimized_original_poses,
std::vector< pose::Pose>& minimized_then_delete_poses,
core::scoring::ScoreFunctionOP const & scorefxn,
core::optimization::AtomTreeMinimizer& minimizer,
core::optimization::MinimizerOptions& options,
core::kinematics::MoveMap& mm,
std::string input_pose_type,
Real mean_ratio){

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace core::kinematics;
	using namespace core::optimization;
	using namespace core::io::silent;
	using namespace core::id;
	using namespace protocols::rna;
	using namespace std;

  ResidueTypeSetCAP rsd_set; //This line is repeated
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "rna" ); //This line is repeated
	std::string outfile_name_option = option[outfile_name];

	if(input_pose_type=="cluster"){
		outfile_name_option=outfile_name_option;
	} else if(input_pose_type=="N_data"){
		outfile_name_option="N_"+outfile_name_option;
	}else {
		std::cout << "Invalid input_pose choice " << std::endl;
		exit (1);
	}
	Size const write_score_only = option[print_score_only] ;
  SilentFileData silent_file_data;

	Size num_residue_reb = RunTimeParameters::rebuild_residue_list.size();

	//Analyze Multiple job folder at the same time?
	//For now assume that we are in the correct job folder


	for(Size current_rebuild_num=1; current_rebuild_num <= num_residue_reb; current_rebuild_num++) {
		Size full_pose_reb_res=get_full_pose_reb_res(current_rebuild_num);
		if(option[build_one_res_mode]!=0){
			if(option[build_one_res_mode]!=int(full_pose_reb_res)) continue;
		}


		std::string pose_filename;
		if(input_pose_type=="cluster"){
			create_filename(pose_filename, "cluster", 0 , current_rebuild_num);
			outfile_name_option=outfile_name_option;
		} else if(input_pose_type=="N_data"){
			create_filename(pose_filename, "N_data", 0 , current_rebuild_num);
			outfile_name_option="N_"+outfile_name_option;
		} else {
			std::cout << "Invalid input_pose choice " << std::endl;
			exit (1);
		}


		std::ifstream infile;
		infile.open(pose_filename.c_str());

		std::map< std::string, Size > column_pos_map;
		create_column_pos_map(column_pos_map);
		utility::vector1<pose_data_struct> all_poses_data;
		if(input_pose_type=="cluster"){
			import_pose_data(all_poses_data,  infile, column_pos_map, 'C');
		} else if(input_pose_type=="N_data"){
			import_pose_data(all_poses_data,  infile, column_pos_map, 'N');
		}	else {
			std::cout << "Invalid input_pose choice " << std::endl;
			exit (1);
		}
		infile.close();

		//import pose;
		for(Size i=1; i<=all_poses_data.size(); i++){

	  		std::string pose_name="";
	  		pose_name.append(all_poses_data[i].tag);
	  		pose_name.append(".pdb");
				std::cout << "Importing pose: " << pose_name << std::endl;
				all_poses_data[i].pose_OP = new pose::Pose;
				core::import_pose::pose_from_pdb((*all_poses_data[i].pose_OP), *rsd_set, pose_name);
				Setup_rebuild_pose_virtual_phosphate((*all_poses_data[i].pose_OP), current_rebuild_num, false);
				if(current_rebuild_num==num_residue_reb) Setup_chain_break_scoring((*all_poses_data[i].pose_OP), num_residue_reb, current_rebuild_num);
		}


		Real sum_ratio=0;
		Size count=0;

		for(Size i=1; i<=all_poses_data.size(); i++){
			all_poses_data[i].new_score=(*scorefxn)(*all_poses_data[i].pose_OP);

			if(option[exclude_bulge_residue_rmsd]==true && current_rebuild_num==num_residue_reb){
				all_poses_data[i].loop_rmsd_exclude_bulge=Loop_rmsd_exclude_last_residue_rebuild(original_poses[current_rebuild_num], (*all_poses_data[i].pose_OP), current_rebuild_num, true);
				all_poses_data[i].loop_rmsd_cor_exclude_bulge=Loop_rmsd_exclude_last_residue_rebuild(minimized_then_delete_poses[current_rebuild_num], (*all_poses_data[i].pose_OP), current_rebuild_num, true);
			}	else{
				all_poses_data[i].loop_rmsd_exclude_bulge=Loop_rmsd(original_poses[current_rebuild_num], (*all_poses_data[i].pose_OP), current_rebuild_num, false);
				all_poses_data[i].loop_rmsd_cor_exclude_bulge=Loop_rmsd(minimized_then_delete_poses[current_rebuild_num], (*all_poses_data[i].pose_OP), current_rebuild_num, false);
			}


			Real ratio=(all_poses_data[i].new_score/all_poses_data[i].score);
			if(ratio<3 && ratio>(1/3)){ //screen for outliers
					sum_ratio=sum_ratio+ratio;
					count++;
			}
		}


		if(input_pose_type!="N_data") {
			mean_ratio= sum_ratio/count;
		}

		std::cout << "mean_ratio= " << mean_ratio << std::endl;


		for(Size i=1; i<=all_poses_data.size(); i++){//Normalize the score
			all_poses_data[i].norm_new_score=(all_poses_data[i].new_score/mean_ratio);
		}


		std::string filename;
		std::string silent_file;
		std::ofstream outfile;


		create_filename(filename, outfile_name_option , 0 , current_rebuild_num);
		create_filename(silent_file, outfile_name_option+"_long" , 0 , current_rebuild_num);
	  outfile.open(filename.c_str());



	  for(Size i=1; i <= all_poses_data.size(); i++){
			output_data_rescore(silent_file_data, silent_file, all_poses_data[i], outfile, current_rebuild_num);
		}
		outfile.close();

		create_filename(filename, outfile_name_option+"_near_native" , 0 , current_rebuild_num);
		create_filename(silent_file, outfile_name_option+"_near_native"+"_long" , 0 , current_rebuild_num);
	  outfile.open(filename.c_str());

	  for(Size i=1; i <= all_poses_data.size(); i++){
			if(all_poses_data[i].loop_rmsd_cor_exclude_bulge< 1.5 || all_poses_data[i].loop_rmsd_cor< 1.5){
				output_data_rescore(silent_file_data, silent_file, all_poses_data[i], outfile, current_rebuild_num);
			}
		}
		outfile.close();

		create_filename(filename, outfile_name_option+"_non_native" , 0 , current_rebuild_num);
		create_filename(silent_file, outfile_name_option+"_non_native"+"_long" , 0 , current_rebuild_num);
	  outfile.open(filename.c_str());


	  for(Size i=1; i <= all_poses_data.size(); i++){
			if(all_poses_data[i].loop_rmsd_cor_exclude_bulge> 2.5 || all_poses_data[i].loop_rmsd_cor> 2.5){
				output_data_rescore(silent_file_data, silent_file, all_poses_data[i], outfile, current_rebuild_num);
			}
		}
		outfile.close();
	} //current_rebuild_num
	return mean_ratio;
}

void
Rescore_pose_list_evelope(
std::vector< pose::Pose>& original_poses,
std::vector< pose::Pose>& minimized_original_poses,
std::vector< pose::Pose>& minimized_then_delete_poses,
core::scoring::ScoreFunctionOP & scorefxn,
core::optimization::AtomTreeMinimizer& minimizer,
core::optimization::MinimizerOptions& options,
core::kinematics::MoveMap& mm){

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	std::string input_pose_type = option[input_pose];

	Real mean_ratio=Rescore_pose_list(original_poses, minimized_original_poses, minimized_then_delete_poses, scorefxn, minimizer, options, mm, "cluster", 999);

	if(input_pose_type=="N_data"){
		Rescore_pose_list(original_poses, minimized_original_poses, minimized_then_delete_poses, scorefxn, minimizer, options, mm, "N_data", mean_ratio);
	}

}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//build_next_nucleotide(original_poses, minimized_original_poses, minimized_then_delete_poses, original_pose, scorefxn, minimizer, options, mm,  silent_file_data, current_pose, current_pose_chain_break, base_pose, old_tag , rotamer_group_min, rotamer_group_max, current_rebuild_num, total_rotamer_group, job_number);

void
build_next_nucleotide(
std::vector< pose::Pose>& original_poses,
std::vector< pose::Pose>& minimized_original_poses,
std::vector< pose::Pose>& minimized_then_delete_poses,
pose::Pose& original_pose,
core::scoring::ScoreFunctionOP & scorefxn,
core::optimization::AtomTreeMinimizer& minimizer,
core::optimization::MinimizerOptions& options,
core::kinematics::MoveMap& mm,
core::io::silent::SilentFileData& silent_file_data,
pose::Pose& current_pose,
pose::Pose& current_pose_chain_break,
pose::Pose& base_pose,
std::string old_tag,
Size rotamer_group_min,
Size rotamer_group_max,
Size current_rebuild_num,
Size total_rotamer_group,
Size job_number){

		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace core::chemical;
		using namespace core::conformation;
		using namespace core::scoring;
		using namespace core::scoring::constraints;
		using namespace core::kinematics;
		using namespace core::optimization;
		using namespace core::io::silent;
		using namespace core::id;
		using namespace protocols::rna;

		ResidueTypeSetCAP rsd_set; //This line is repeated
		rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "rna" ); //This line is repeated

		bool const verbose_option = option[Verbose];
		Size num_pose_kept= option[pose_kept];
		Real rmsd_cutoff=option[suite_rmsd_cutoff];
		std::string const rotamer_set_option = option[rotamer_set];
		std::string const native_screen_option = option[native_screen];
		std::string const HO2star_sampling_option = option[HO2star_sampling];
//  	Size const node_number = option[node_number_input ];
		Size num_residue_reb = RunTimeParameters::rebuild_residue_list.size();
		std::string const starting_pose_option = option[starting_pose];
		bool const prepend_res = Is_prepend(current_rebuild_num);
		Size multiplier=2;
		Size full_pose_reb_res=get_full_pose_reb_res(current_rebuild_num);
		Size reb_res=Convert_to_partial_pose_seq_num(full_pose_reb_res, current_rebuild_num);

		std::cout << "Start rebuilding rebuild_num= " << current_rebuild_num ;
		std::cout << "  full_pose_reb_res= " << full_pose_reb_res ;
		std::cout << "  rebuild_pose_reb_res= " << reb_res << std::endl;


		std::cout << "Input HO2star_sampling_option : " << HO2star_sampling_option << std::endl;
		if((HO2star_sampling_option!="O2star_trial") && (HO2star_sampling_option!="Virtual_HOStar") && (HO2star_sampling_option!="Native")){
			std::cout << "Input HO2star_sampling_option : " << HO2star_sampling_option <<  " is invalid" << std::endl;
			exit (1);
		}

		std::vector< pose::Pose> *template_pose;
		if(starting_pose_option=="minimize") {
				template_pose=&minimized_then_delete_poses;
		} else {
				template_pose=&original_poses;
		}
		dump_pdb( (*template_pose)[current_rebuild_num], "template_" + lead_zero_string_of(current_rebuild_num, 2) +".pdb");

		utility::vector1 <Real> native_rotamer=get_native_rotamer((*template_pose)[num_residue_reb], current_rebuild_num);
		utility::vector1< pose_data_struct2 > pose_data_list;

		//Need to add virtual phosphate to every loop 3' end, don't actually need to be this if base pose is original_pose.
		Setup_rebuild_pose_virtual_phosphate(base_pose, current_rebuild_num-1, true);
		(*scorefxn)(base_pose);

		//If have allowed for more than one chain_break, then would have to reinitialize existing chemical::CUTPOINT_LOWER, chemical::CUTPOINT_UPPER, AtomPairConstraint and AngleConstraint

///////////////Determine whether to use rotamer_north or rotamer_south (for append case)////////////////////////
		std::cout << "before backbone_rotamers_group_creation" << std::endl;

		utility::vector1< utility::vector1 <utility::vector1 <Real > > > backbone_rotamers_groups;
		if( prepend_res) {
			//rotamer_groups_north and south are the same in this case
			get_rotamers( backbone_rotamers_groups, ALL, NORTH);
		} else {
			std::string previous_res_puckerstate=Previous_Res_PuckerState(base_pose, current_rebuild_num);
   		std::cout << "Previous res pucker state is " <<  previous_res_puckerstate << std::endl;
   		if (previous_res_puckerstate=="NORTH") {
   			get_rotamers( backbone_rotamers_groups, NORTH, ALL );
   		} else {
   			get_rotamers( backbone_rotamers_groups, SOUTH, ALL );
   		}
   	}
   	Size rotamer_group_size=backbone_rotamers_groups[1].size();

		if(total_rotamer_group==backbone_rotamers_groups.size()){
			std::cout << "Correct total_rotamer_group value" << std::endl;
		}else {
			std::cout << "Error: wrong total_rotamer_group value" << std::endl;
			exit (1);
		}

		std::cout << "The number of rotamer in each group=" << rotamer_group_size << std::endl;
 		std::cout << "The number of rotamers group=" << total_rotamer_group << std::endl;
		std::cout << "after backbone_rotamers_group_creation" << std::endl;


////////////////////////////Open output text file/////////////////////////////////////////

		std::string filename;
		create_filename(filename, "data", job_number, current_rebuild_num);
		std::ofstream outfile;
	  outfile.open(filename.c_str());
  	outfile_column_name(outfile, current_rebuild_num);

	  std::string silent_file;
	  create_filename(silent_file, "long", job_number, current_rebuild_num);

		std::cout << "Outfile original pose data" << std::endl;

		//Before only first node output original data, better this way since can check for consistency between nodes.

		 for(Size i=1; i <= num_residue_reb; i++){
				outfile_original_pose(silent_file_data, silent_file, original_poses[i], minimized_original_poses[i], minimized_then_delete_poses[i], "orig_pose"+lead_zero_string_of(i, 2), (*template_pose)[num_residue_reb], original_poses[i], outfile, i, scorefxn);
			}

			for(Size i=1; i <= num_residue_reb; i++){
				outfile_original_pose(silent_file_data, silent_file, original_poses[i], minimized_original_poses[i], minimized_then_delete_poses[i], "min_pose"+lead_zero_string_of(i, 2), (*template_pose)[num_residue_reb], minimized_original_poses[i], outfile, i, scorefxn);
			}

			for(Size i=1; i <= num_residue_reb; i++){
				outfile_original_pose(silent_file_data, silent_file, original_poses[i], minimized_original_poses[i], minimized_then_delete_poses[i], "corr_pose"+lead_zero_string_of(i, 2), (*template_pose)[num_residue_reb],  minimized_then_delete_poses[i], outfile, i, scorefxn);
			}


/////////////////////////////////The actual rebuilding code//////////////////////////////////////////////////////

		count_struct count_data={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
 		std::string tag;
 		Real current_score=0;
	  std::vector<Real> temp_vector(4 /*Number additional output variable*/, 9999.99);

//	  rotamer_ID_struct rotamer_ID={1,1,1,1,1,1,1,1,1,1,1,1,1,1};
//	  utility::vector1 <Real > current_rotamer(11 /* # torsion angles+1*/, 9999.99);;

//////////////////////////////////////Append/Prepend the rebuild residue////////////////////////////////////////////

 		current_pose=attach_next_nucleotide(base_pose, original_pose, current_rebuild_num);

 		Freeze_sugar_torsions(mm, current_pose); //Not sure what the purpose of this is, since there is no minimization during sampling
		protocols::rna::ensure_phosphate_nomenclature_matches_mini_parin( current_pose);

		if(HO2star_sampling_option=="Native"){
			//This mode is not well tested.
			std::cout << "get_O2_star_torsion_list of full template pose" << std::endl;
			utility::vector1 <Real> full_template_pose_O2star_torsion_list= get_O2star_torsion_list((*template_pose)[num_residue_reb], num_residue_reb);
			current_pose.set_torsion( TorsionID( reb_res, id::CHI, 4 ), full_template_pose_O2star_torsion_list[full_pose_reb_res] );
		}

/////////////////////////////////////Setup fa_atr ,fa_rep scoring///////////////////////////////////////////
		std::cout << "setup fa_atr, fa_rep scoring" << std::endl;
		scoring::EMapVector	energy_map;

		pose::Pose base_pose_screen=base_pose;

		ScoreFunctionOP atr_rep_screening_scorefxn = ScoreFunctionFactory::create_score_function( "fa_atr" ); //Should make this internal in the future

//		if(option[Virtual_H2star]){
			Delete_O2Star_hydrogen(base_pose_screen);
//		}

		dump_pdb( base_pose_screen, "base_pose_screen.pdb" );

		(*atr_rep_screening_scorefxn)(base_pose_screen);
		energy_map=base_pose_screen.energies().total_energies();
		Real base_atr_score=atr_rep_screening_scorefxn->get_weight(fa_atr)*energy_map[scoring::fa_atr];
		Real base_rep_score=atr_rep_screening_scorefxn->get_weight(fa_rep)*energy_map[scoring::fa_rep];
		std::cout << "base_rep= " << base_rep_score << " base_atr= " << base_atr_score << std::endl;

		pose::Pose current_pose_screen=current_pose;

// 		if(option[Virtual_H2star]){
 			Delete_O2Star_hydrogen(current_pose_screen);
// 		}
 		dump_pdb( current_pose, "current_pose_check.pdb" );
 		dump_pdb( current_pose_screen, "current_pose_screen.pdb" );


////////////////////Setup chain break scoring//////////////////////////////////////////////////////
 		ScoreFunctionOP constraint_scorefxn;
 		RNA_LoopCloser rna_loop_closer;

		//Does_rebuild_residue_close_chain(Size& current_rebuild_num){}
		//Assume one loop/chainbreak to be rebuild for now
 		if(current_rebuild_num==num_residue_reb){
 			if(prepend_res){
 				core::pose::remove_variant_type_from_pose_residue( current_pose, "VIRTUAL_PHOSPHATE", reb_res );
			 	protocols::rna::ensure_phosphate_nomenclature_matches_mini_parin( current_pose);
				//Non-Native torsion angles, does torsion value after CCD depend on input torsion choice?
 				current_pose.set_torsion( TorsionID( reb_res, id::BB,  2 ), 0.0  ); //beta i
 				current_pose.set_torsion( TorsionID( reb_res, id::BB,  3 ), 0.0 ); //gamma i
  		}

  		//This is for the CCD alogorithm.
  		Size five_prime_res, three_prime_res;
			if(prepend_res){
				five_prime_res=reb_res-1;
				three_prime_res=reb_res;
			}else {
				five_prime_res=reb_res;
				three_prime_res=reb_res+1;
			}
			core::pose::add_variant_type_to_pose_residue( current_pose , chemical::CUTPOINT_LOWER, five_prime_res );
			core::pose::add_variant_type_to_pose_residue( current_pose , chemical::CUTPOINT_UPPER, three_prime_res );

    	current_pose_chain_break=current_pose;

			//Setup up the chain break scoring for the current_pose_chain_break but not the
			//the current pose which will be done later.

 			Setup_chain_break_scoring(current_pose_chain_break, num_residue_reb, current_rebuild_num);
	 		constraint_scorefxn = ScoreFunctionFactory::create_score_function( "constraints" ); //Should make this internal in the future
 	  }
//////////////Sugar scoring screw up pose selection during rotamer sampling stage, so need to turn it off///////////////////

   	scorefxn->set_weight( rna_sugar_close, 0.0 );
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if(HO2star_sampling_option=="Virtual_HOStar"){
			current_pose.set_torsion( TorsionID( reb_res, id::CHI, 4 ), 0 );
			std::cout << "get_O2star_torsion_list(current_pose, current_rebuild_num), BEFORE Delete_O2Star_hydrogen" << std::endl;
			get_O2star_torsion_list(current_pose, current_rebuild_num);
			dump_pdb( current_pose, "Before_removing_O2Star_hydrogen_variant_type.pdb" );
 			Delete_O2Star_hydrogen(current_pose);
			std::cout << "get_O2star_torsion_list(current_pose, current_rebuild_num), AFTER Delete_O2Star_hydrogen" << std::endl;
			get_O2star_torsion_list(current_pose, current_rebuild_num);
			dump_pdb( current_pose, "After_adding_O2Star_hydrogen_variant_type.pdb" );
		}

////////////////////////////////////Rotamer sampling/////////////////////////////////////////////////////////////////

		std::time_t start_time=time(NULL);

		Size group_rotamer=rotamer_group_min;
		while( group_rotamer<= rotamer_group_max) {
				if(group_rotamer > total_rotamer_group) {
					std::cout << "group_rotamer > total_rotamer_group(108), break" << std::endl;
					break; //important just for first residue rebuilt
 				}

   			for(Size subgroup_rotamer=1; subgroup_rotamer<=rotamer_group_size; subgroup_rotamer++){
   					utility::vector1<Real> current_rotamer=backbone_rotamers_groups[group_rotamer][subgroup_rotamer];

   				//////Near Native Torsion screening////////////////////////////////////////////////////////////////////////
   					if(rotamer_set_option=="near_native_torsions"){
   						if(!Is_NearNative_Rotamer(current_rotamer, native_rotamer, prepend_res)) continue;
   						count_data.Near_Native_Rotamer_Count++;
   					}

					//////Near Native rmsd screening stage 1//////////////////////////////////////////////////////////////////
					if(native_screen_option=="rmsd"){
						apply_rotamer(current_pose_screen, reb_res, current_rebuild_num, current_rotamer);
						Real rmsd_compare=suite_rmsd((*template_pose)[current_rebuild_num], current_pose_screen, current_rebuild_num);

						if(rmsd_compare>(rmsd_cutoff+1)){
 								continue;
						}
						count_data.rmsd_count++;
					}

					////////////////////////////////////////////////////////////////////////////////////////////////////////

					fine_rotamer_ID_struct fine_rotamer_ID={0,0,0,0,0,0,false,0};

   					while(true){ //fine_sampling
/*
   						if(count_data.tot_rotamer_count==1000) {
   							  std::time_t end_time_temp=time(NULL);
   								std::cout<< "Sample rotamer= " << (long)(end_time_temp - start_time) << "seconds";
 							}
 							if(count_data.tot_rotamer_count==1001) break;
*/

   						//Refine current_rotamer
   						get_fine_rotamer(backbone_rotamers_groups[group_rotamer][subgroup_rotamer], current_rotamer, fine_rotamer_ID, current_rebuild_num);
   						if(fine_rotamer_ID.finish==true) break;
   						count_data.tot_rotamer_count++;

////////////////////////////////////////////////Screening//////////////////////////////////////////////////////////////

						apply_rotamer(current_pose_screen, reb_res, current_rebuild_num, current_rotamer);

            ////////////////Near Native rmsd screening stage 2(fine_sampling)//////////////////////////////////////////
						if(native_screen_option=="rmsd"){
							if(suite_rmsd((*template_pose)[current_rebuild_num], current_pose_screen, current_rebuild_num)>(rmsd_cutoff)){
 								continue;
							}
							count_data.fine_rmsd_count++;
							if(verbose_option){
								std::cout << "fine_rmsd_count = " << count_data.fine_rmsd_count;
								std::cout << " rmsd_count = " << count_data.rmsd_count;
								std::cout << " total count= " << count_data.tot_rotamer_count << std::endl;
							}
   					}

   					/////////////Check_chain_closable screening///////////////////////////////////////////////////////////////
						if(current_rebuild_num==(num_residue_reb-1)) { //Need to be updated if allow more than 1 chain_break.
							if(Check_chain_closable(current_pose_screen, count_data, "Sampling", current_rebuild_num)==false){
								continue;
							}
						}

   				 	/////////////Chain_break_screening or atr_rep_screening/////////////////////////////////////////////////////////////////////////
						if(current_rebuild_num==num_residue_reb){//Assume one loop/chainbreak to be rebuild for now

							if (Chain_break_screening(current_pose_screen, current_pose_chain_break, reb_res, current_rotamer, count_data, atr_rep_screening_scorefxn, constraint_scorefxn, rna_loop_closer, base_rep_score, current_rebuild_num)){

								Copy_torsion_angle(current_pose, current_pose_chain_break, reb_res, current_rebuild_num);
							} else continue;

						} else {
   						if (!Atr_rep_screening(current_pose_screen, reb_res, current_rotamer, count_data, atr_rep_screening_scorefxn, base_rep_score, base_atr_score, rotamer_set_option)){
   						   continue;
   						}
						}

						///o2star_trail, only rebuild residue
						if((HO2star_sampling_option=="O2star_trial")){
								current_pose.set_torsion( TorsionID( reb_res, id::CHI, 4 ), 0 );
								o2star_minimize(current_pose, scorefxn);
								//pack::rotamer_trials( current_pose, *(scorefxn), task);
						}
						////////////////Add pose to pose_data_list if pose have good score////////////////////////////////////////////

						Pose_selection_by_full_score(pose_data_list, current_pose, reb_res, current_rotamer, (*template_pose)[num_residue_reb], current_rebuild_num, original_poses, minimized_original_poses, minimized_then_delete_poses, num_pose_kept, multiplier,  silent_file_data, silent_file, group_rotamer, subgroup_rotamer, old_tag, scorefxn,  num_residue_reb, outfile, fine_rotamer_ID, count_data);

				}//fine_rotamer_while_loop
 				}//for loop

			group_rotamer++;
  	}//end of rebuilt while loop

  	Cluster_poses_test(num_pose_kept, pose_data_list, current_rebuild_num); //cluster pose one last time.

		std::time_t end_time1=time(NULL);
		std::time_t start_time_fine_sampling=time(NULL);

//This function is not up to date
// 		if((current_rebuild_num!=num_residue_reb) && option[Finer_sampling]){
//			Fine_sampling(pose_data_list, current_rebuild_num, num_residue_reb, scorefxn, atr_rep_screening_scorefxn, (*template_pose)[num_residue_reb], reb_res,  silent_file, silent_file_data, outfile, original_poses, minimized_original_poses, minimized_then_delete_poses, base_atr_score, base_rep_score, num_pose_kept, multiplier, rep_cutoff);
// 		}

 		std::time_t end_time_fine_sampling=time(NULL);

 		if(option[minimize_and_score_sugar]==true){
  		scorefxn->set_weight( rna_sugar_close, 0.7 );
  	} else {
    	scorefxn->set_weight( rna_sugar_close, 0.0 );
 		}


 		for(Size i=1; i<=pose_data_list.size(); i++){

 			tag=pose_data_list[i].tag;
 			current_pose=pose_data_list[i].pose;
 			current_score=pose_data_list[i].score;

			if(HO2star_sampling_option=="Virtual_HOStar"){
				if(i<=2 && job_number==1){
					std::cout << "get_O2star_torsion_list(current_pose, current_rebuild_num), BEFORE removing_O2Star_hydrogen variant_type/Before o2star minimize" << std::endl;
					get_O2star_torsion_list(current_pose, current_rebuild_num);
					dump_pdb( current_pose, "Before_Remove_O2Star_hydrogen_variant_type_"+tag+".pdb" );
				}

				Remove_O2Star_hydrogen_variant_type(current_pose);

				if(i<=2 && job_number==1){
					std::cout << "get_O2star_torsion_list(current_pose, current_rebuild_num), AFTER removing_O2Star_hydrogen variant_type/Before o2star minimize" << std::endl;
					get_O2star_torsion_list(current_pose, current_rebuild_num);
					dump_pdb( current_pose, "After_Remove_O2Star_hydrogen_variant_type_"+tag+".pdb" );
				}
			}

			o2star_minimize(current_pose, scorefxn);


 			//Discard really bad score pose that happend to get to this stage
 			if(current_score>0) continue;

 			//Setup up harmonic chain break scoring and CCD...not the best way...
			//This assume that they is only 1 chain-break, will need to update the if statement once we allow multiple chain break.
 			if(current_rebuild_num==num_residue_reb) Setup_chain_break_scoring(current_pose, num_residue_reb, current_rebuild_num);

			//Before Minimize////////////////////////////////////////////////////////////////////////////////
 			tag[0]='B';
    	std::cout << "tag= " << tag << std::endl;

			std::vector< Real> diff_torsions(11 /* # torsion angles+1*/, 9999.99);
 			Calculate_diff_torsions(diff_torsions, current_pose, (*template_pose)[num_residue_reb], current_rebuild_num);

			output_data_struct output_data={0,0,0,0,0,0,0,0};
			output_data.current_score=(*scorefxn)(current_pose);

			Calculate_rmsd(output_data, current_pose, original_poses[current_rebuild_num], minimized_original_poses[current_rebuild_num], minimized_then_delete_poses[current_rebuild_num], current_rebuild_num, false);

			Output_data_general(silent_file_data, silent_file, output_data, current_pose, tag, current_rebuild_num, outfile, diff_torsions);

			if(HO2star_sampling_option=="Virtual_HOStar"){
				if(i<=2 && job_number==1){
					std::cout << "get_O2_star_torsion_list of " << tag << std::endl;
			    get_O2star_torsion_list(current_pose, current_rebuild_num);
				}
			}
//			dump_pdb( current_pose, tag+".pdb" );


 			//////////////////////////////////////////////////////////////////////////////////////////////////


			if(option[ quick_test ]==false){

					if(option[full_rebuild_pose_minimize]==false){
						Freeze_unrebuild_residues(mm, current_rebuild_num , current_pose);
					}
					//Since the movemap was update, might need to refreeze the sugar_torsions
					if(option[minimize_and_score_sugar]==false){
						Freeze_sugar_torsions(mm, current_pose);
					} else {
//						UnFreeze_sugar_torsions(mm, current_pose); //This function doesn't work..forgot why it doesn't work though...
					}
					minimizer.run( current_pose, mm, *(scorefxn), options );
					o2star_minimize(current_pose, scorefxn);
			}

			tag[0]='M';

			if(HO2star_sampling_option=="Virtual_HOStar"){
				if(i<=2 && job_number==1){
					std::cout << "get_O2_star_torsion_list of " << tag << std::endl;
		    	get_O2star_torsion_list(current_pose, current_rebuild_num);
				}
			}
			////////////////////Final screening before output/////////////////////////////////////////////
			if(current_rebuild_num==num_residue_reb) {
						Size five_prime_res;
						if(prepend_res){
							five_prime_res=reb_res-1;
						}else {
							five_prime_res=reb_res;
						}
						rna_loop_closer.apply( current_pose, five_prime_res);
						//Need this for alpha torsion of 3 prime res to be correctly setup, important for O1P and O2P atom position
						Copy_torsion_angle(current_pose, current_pose, reb_res, current_rebuild_num);
					//Need to check if this work (May 5, 2009)
			}

			if(current_rebuild_num==(num_residue_reb-1)) {
				if(Check_chain_closable(current_pose, count_data, "Output", current_rebuild_num)==false){
					tag[0]='F'; //Mark as unclosable loop failure so that will not be use to rebuild next (last) residue
				}
			}

			//For native screen==rsmd mode
			//Basically, these will be discard and not continue rebuild, except last residue rebuild
			Real rmsd_compare=suite_rmsd((*template_pose)[current_rebuild_num], current_pose, current_rebuild_num);
			if((native_screen_option=="rmsd") && (rmsd_compare>(rmsd_cutoff+1)) && (current_rebuild_num!=num_residue_reb)){
				tag[0]='R'; //Mark as Rmsd failure
			}
			////////////////////////////////////////////////////////////////////////////////////////////////////////////


 			Calculate_diff_torsions(diff_torsions, current_pose, (*template_pose)[num_residue_reb], current_rebuild_num);

			output_data_struct output_data_M={0,0,0,0,0,0,0,0};
			output_data_M.current_score=(*scorefxn)(current_pose);

			//Verbose just for the first five pose...for debugging purposes
			if(i<=5){
				std::cout << "consistency check in minimize phase, i= " << i << std::endl;
				protocols::rna::ensure_phosphate_nomenclature_matches_mini_parin( original_poses[current_rebuild_num]);
				protocols::rna::ensure_phosphate_nomenclature_matches_mini_parin( minimized_original_poses[current_rebuild_num]);
				protocols::rna::ensure_phosphate_nomenclature_matches_mini_parin( minimized_then_delete_poses[current_rebuild_num]);
				protocols::rna::ensure_phosphate_nomenclature_matches_mini_parin( current_pose);

				Calculate_rmsd(output_data_M, current_pose, original_poses[current_rebuild_num], minimized_original_poses[current_rebuild_num], minimized_then_delete_poses[current_rebuild_num], current_rebuild_num, true);
			} else {
				Calculate_rmsd(output_data_M, current_pose, original_poses[current_rebuild_num], minimized_original_poses[current_rebuild_num], minimized_then_delete_poses[current_rebuild_num], current_rebuild_num, false);
			}

			std::cout << "tag= " << tag << std::endl;

			Output_data_general(silent_file_data, silent_file, output_data_M, current_pose, tag, current_rebuild_num, outfile, diff_torsions);
			dump_pdb( current_pose, tag+".pdb" );

			if(output_data_M.loop_rmsd_wrt_correct_minimized < 2.0 && (native_screen_option!="rmsd")){ //Select and output near native pose if NOT in rmsd native_screen_option mode.
				count_data.Near_Native_Pose_Count++;
				tag[0]='N';
				Output_data_general(silent_file_data, silent_file, output_data_M, current_pose, tag, current_rebuild_num, outfile, diff_torsions);
				dump_pdb( current_pose, tag+".pdb" );
			}
  	}

  	std::time_t end_time2=time(NULL);
  	std::cout << "pass_before_last_res_reb_chain_break sampling= ";
  	std::cout << count_data.pass_before_last_res_reb_chain_break_U;
  	std::cout << " out of " << count_data.total_before_last_res_reb_chain_break_U << std::endl;
  	std::cout << "pass_before_last_res_reb_chain_break minimize= ";
  	std::cout << count_data.total_before_last_res_reb_chain_break_M;
  	std::cout << " out of " << count_data.total_before_last_res_reb_chain_break_M << std::endl;
  	std::cout << "fine_rmsd_count = " << count_data.fine_rmsd_count;
		std::cout << " rmsd_count = " << count_data.rmsd_count;
  	std::cout << " Near_Native_Rotamer_Count= " << count_data.Near_Native_Rotamer_Count;
		std::cout << " Near_Native_Pose_Count= " << count_data.Near_Native_Pose_Count;
  	std::cout << " C5_O3_count= " << count_data.C5_O3_distance_count;
		std::cout << " angle_n= " << count_data.good_angle_count << " dist_n= " << count_data.good_distance_count;
		std::cout << " rep_n= " << count_data.good_rep_rotamer_count << " atr_n= " << count_data.good_atr_rotamer_count;
    std::cout << " both= " << count_data.both_count << " tot= " << count_data.tot_rotamer_count << std::endl;
  	std::cout<< "time taken" << std::endl;
  	std::cout<< "Sample rotamer= " << (long)(end_time1 - start_time) << "seconds";
  	std::cout<< "Fine sample old= " << (long)(end_time_fine_sampling - start_time_fine_sampling) << "seconds";
  	std::cout<< "  Minimize= " << (long)(end_time2 - end_time_fine_sampling) << "seconds." << std::endl;

  	outfile.close();

  	//This create the status.txt file which tells the central_computer that the job by this computer node is finish

//		Create_status_file(node_number, current_rebuild_num);

}//end of build_next_nucleotide


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
iterative_build_one_res_test(
std::vector< pose::Pose>& original_poses,
std::vector< pose::Pose>& minimized_original_poses,
std::vector< pose::Pose>& minimized_then_delete_poses,
pose::Pose& original_pose,
core::scoring::ScoreFunctionOP & scorefxn,
core::optimization::AtomTreeMinimizer& minimizer,
core::optimization::MinimizerOptions& options,
core::kinematics::MoveMap& mm,
core::io::silent::SilentFileData& silent_file_data) {

//current_poses, currents_scores, rotamer_branch_state all have dimension in the number of residue rebuild dimension...can take this out

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace core::kinematics;
	using namespace core::optimization;
	using namespace core::io::silent;
	using namespace core::id;
	using namespace protocols::rna;

////////////////////////////Options///////////////////////////////////////////////////////////////////////////
  Size const total_nodes = option[total_nodes_input ] ;
	Size const first_rebuild_num = option[ hack_first_rebuild_num];
	std::string const starting_pose_option = option[starting_pose];
  Size const node_number = option[node_number_input ] ;
	Size num_residue_reb = RunTimeParameters::rebuild_residue_list.size();
	std::string const rotamer_set_option = option[rotamer_set];


	//Need to define out here for conformation viewer to work
	pose::Pose current_pose = original_poses[num_residue_reb];
	pose::Pose current_pose_chain_break=original_poses[num_residue_reb];

	if(option[graphics]==true){
		protocols::viewer::add_conformation_viewer( current_pose.conformation(), "current", 400, 400 );
	}
///////////////// Declare type of residue as RNA type///////////////////////////////////////////////////////////

  ResidueTypeSetCAP rsd_set; //This line is repeated
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "rna" ); //This line is repeated
////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	Size const total_rotamer_group=get_total_rotamer_group();

//	std::vector<atom_pair> atom_pair_constraint_list;

	for(Size current_rebuild_num=first_rebuild_num; current_rebuild_num <= num_residue_reb; current_rebuild_num++) {

		Size full_pose_reb_res=get_full_pose_reb_res(current_rebuild_num);
		Size reb_res=Convert_to_partial_pose_seq_num(full_pose_reb_res, current_rebuild_num);

		pose::Pose base_pose;
		Size group_rotamer;
		Size rotamer_group_min;
		Size rotamer_group_max;
		Size job_number;
		std::string old_tag="";




/////////////////////////////////////////////////////////////////////////////////////
		std::cout << "current_rebuild_num= " << current_rebuild_num << std::endl;

		if(current_rebuild_num==1 || option[build_one_res_mode]!=0){ //first rebuild_num of full run or build_one_res_mode
				if((option[build_one_res_mode]!=0) && (option[build_one_res_mode]!=99)){ //!=0 means in build_one_res_mode...!=99 means that don't build every residue
					if(option[build_one_res_mode]!=int(full_pose_reb_res)) continue;
				}

				if(total_rotamer_group>=total_nodes){
					std::cout << "total_rotamer_group>=total_nodes case:";
					if(total_rotamer_group%total_nodes!=0) {
						std::cout << " total_rotamer_group is not divisible by total_nodes" << std::endl;
						exit (1);
					} else {
						std::cout << " total_rotamer_group is divisible by total_nodes" << std::endl;
						rotamer_group_min=(total_rotamer_group/total_nodes)*(node_number-1)+1;
  					rotamer_group_max=(total_rotamer_group/total_nodes)*(node_number);
  				}
				} else {
					std::cout << "total_rotamer_group<total_nodes case:";
						rotamer_group_min=node_number;
						rotamer_group_max=node_number;
				}
				//Note that rotamer_group_min and rotamer_group_max could exceed group_rotamer. This will be catch before rotamer sampling begins
				//The reason why I don't catch this here is so that the node does enter build_next_nucleotide and create a (blank) data.txt file which will be open by central_computer
				std::cout << "rotamer_group_min= " << rotamer_group_min;
				std::cout << " rotamer_group_max= " << rotamer_group_max << std::endl;

				group_rotamer=rotamer_group_min;

				if(starting_pose_option=="minimize") {
					base_pose=minimized_then_delete_poses[current_rebuild_num-1];
				} else {
					base_pose=original_poses[current_rebuild_num-1];
				}
				//Overide and use the selected as the base_pose.
				std::string base_pose_name = option[select_base_pose];
				if(base_pose_name!="default"){
					Import_base_pose(base_pose_name, base_pose, original_poses[current_rebuild_num-1]);
					old_tag=base_pose_name;
				}

				job_number=node_number;
				build_next_nucleotide(original_poses, minimized_original_poses, minimized_then_delete_poses, original_pose, scorefxn, minimizer, options, mm,  silent_file_data, current_pose, current_pose_chain_break, base_pose, old_tag , rotamer_group_min, rotamer_group_max, current_rebuild_num, total_rotamer_group, job_number);

				Create_status_file(job_number, 1, current_rebuild_num);

		} else { //Full (not first rebuild_num)


			Size communication_num=0;
			job_file_data job_data;
			while(read_in_job_file(job_data, node_number, communication_num, num_residue_reb, current_rebuild_num)){
				communication_num++;

				std::cout << "job_data.tag= " << job_data.tag << std::endl;
				rotamer_group_min=job_data.rotamer_group_min;
				rotamer_group_max=job_data.rotamer_group_max;
				job_number=job_data.job_number;

				old_tag=job_data.tag;
				Import_base_pose(old_tag, base_pose, original_poses[current_rebuild_num-1]);

				build_next_nucleotide(original_poses, minimized_original_poses, minimized_then_delete_poses, original_pose, scorefxn, minimizer, options, mm,  silent_file_data, current_pose, current_pose_chain_break, base_pose, old_tag , rotamer_group_min, rotamer_group_max, current_rebuild_num, total_rotamer_group, job_number);
			}
		}
	}//num_residue_rebuild for loop
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//If object are equal should return false?
bool comparefunction(utility::vector1<utility::vector1<int> > suiteBin1, utility::vector1<utility::vector1<int> >  suiteBin2){

		for(Size n=1; n<= suiteBin1.size(); n++){
			 for(Size i=1; i<= suiteBin1[1].size(); i++){

//				Size n=1;
//				Size i=1;
//					while(true){
//						if (i==suiteBin1[1].size()) {
//								i=1;
//								n++;
//								if(n>suiteBin1.size()) return false; //Equality case.
//						}

			 			if(suiteBin1[n][i]>suiteBin2[n][i]) return true;
						else if(suiteBin1[n][i]<suiteBin2[n][i]) return false;
//						else i++;
//					}
				}
	}
	return false; //Equality case.
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
void
Map_test(utility::vector1< utility::vector1 <utility::vector1 <Real > > >& backbone_rotamers_groups_north, pose::Pose base_pose, pose::Pose original_pose){

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace core::kinematics;
	using namespace core::optimization;
	using namespace core::io::silent;
	using namespace core::id;
	using namespace protocols::rna;

	ResidueTypeSetCAP rsd_set; //This line is repeated
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "rna" ); //This line is repeated
	std::string rotamer_set_option = option[rotamer_set];




		//Initialize rotamer_set
		if(rotamer_set_option=="cheat"){
			backbone_rotamers_groups_north.clear();
			get_cheat_rotamers( backbone_rotamers_groups_north, ALL, NORTH, original_pose, 1 );
			Size total_rotamer_group=backbone_rotamers_groups_north.size(); //Number of rotamer groups,
			Size rotamer_group_size=backbone_rotamers_groups_north[1].size();
			std::cout << "The number of rotamer in each group=" << rotamer_group_size << std::endl;
 			std::cout << "The number of rotamers group=" << total_rotamer_group << std::endl;
		}

/////////////////////////////Outfile////////////////////////////////////////

	std::string foldername;



		foldername="residueNum_1";

		std::string filename;

		filename=foldername;
		filename.append("/data");
		filename.append(".txt");

		std::ofstream outfile;
	  outfile.open(filename.c_str());

//  	outfile_column_name(outfile);


////////////////Setup basepose///////////////////////////////////
		pose::Pose pose=base_pose;
		dump_pdb( pose, "before.pdb" );

		for (Size n = 20; n >= 19; n-- ) {
		std::cout << n << std::endl;
		pose.conformation().delete_residue_slow(n);
		}

		for (Size n = 1; n < 17; n++ ) {
		std::cout << n << std::endl;
		pose.conformation().delete_residue_slow(1);
		}




		core::pose::remove_variant_type_from_pose_residue( pose, "VIRTUAL_PHOSPHATE", 2 );

		utility::vector1< std::pair< std::string, std::string> > atom_pairs;

		chemical::AA res_aa = aa_from_name( "RAD" );
		ResidueOP new_rsd = conformation::ResidueFactory::create_residue( *(rsd_set->aa_map( res_aa )[1]) ) ;
		atom_pairs.push_back( std::make_pair( " O3*", " O3*" ) );
		atom_pairs.push_back( std::make_pair( " C3*", " C3*" ) );
		atom_pairs.push_back( std::make_pair( " C4*", " C4*" ) );
		pose.replace_residue( 1, *new_rsd, atom_pairs );


		core::pose::remove_lower_terminus_type_from_pose_residue(pose, 1);	 //need this if prepend from end.
		pose.prepend_polymer_residue_before_seqpos( *new_rsd, 2, true);
		core::pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_PHOSPHATE", 2 );

		dump_pdb( pose, "after.pdb" );

/////////////////////////////////////////Setup Map///////////////////////////////////////////////////////

		Size suite_num=2;

		bool(*fn_pt)(utility::vector1<utility::vector1<int> >, utility::vector1<utility::vector1<int> >) = comparefunction;

		std::map<utility::vector1<utility::vector1<int> > , int, bool(*)(utility::vector1<utility::vector1<int> >, utility::vector1<utility::vector1<int> >) > suiteMap (fn_pt);

		std::map<utility::vector1<utility::vector1<int> > , int, bool(*)(utility::vector1<utility::vector1<int> >, utility::vector1<utility::vector1<int> >) > countMap (fn_pt);
		Size unique_bins=0;
		Size count=0;

		utility::vector1<int> atomBin(3, 0);
//		utility::vector1<utility::vector1<int> > suiteBin(pose.residue(suite_num).nheavyatoms()-7, atomBin);
		utility::vector1<utility::vector1<int> > suiteBin(3, atomBin);


/////////////////////////////////////////Map rotamers/////////////////////////////////////////////////////////
		Real current_score=0;

		for(Size group_rotamer=1; group_rotamer<= backbone_rotamers_groups_north.size(); group_rotamer++) {

   		for(Size subgroup_rotamer=1; subgroup_rotamer<=backbone_rotamers_groups_north[1].size(); subgroup_rotamer++){

				count++;

				apply_rotamer(pose, suite_num, backbone_rotamers_groups_north[group_rotamer][subgroup_rotamer]);
				std::string tag= create_tag2("M_", group_rotamer , subgroup_rotamer, 1);
				dump_pdb(pose, tag+".pdb");


				GetSuiteBin(pose, suiteBin, suite_num);

				if(countMap.find(suiteBin)==countMap.end()){
					 unique_bins++;
					 countMap[suiteBin]=unique_bins;

				} else {
				}

				/////////////////////////////////////////////////////////////////////////////////////////////////
				Size spacing=3;
				   outfile << std::setw(125) << std::left << tag;            // output description

   				for(Size n=1; n<= suiteBin.size(); n++){
						for(Size i=1; i<= suiteBin[1].size(); i++){
   		 				outfile << std::setw(spacing) << std::setprecision(3) << std::left << suiteBin[n][i];
						}
					}
					outfile << std::setw(spacing) << std::setprecision(3) << std::left << countMap[suiteBin];
					outfile << "\n";
				////////////////////////////////////////////////////////////////////////////////////////////////////

				if(suiteMap.find(suiteBin)==suiteMap.end()){
					 suiteMap[suiteBin]=1;
				} else {
					 suiteMap[suiteBin]++;
				}

//				std::cout << "countMap[suiteBin]= " << countMap[suiteBin];
//				std::cout << "  suiteMap[suiteBin]= " << suiteMap[suiteBin] << std::endl;

				std::cout << "unique_bins= " << unique_bins;
				std::cout << "  count= " << count << std::endl;
////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Must have accidentally delete num_pose_kept and pose_data_list
//	 			Update_pose_data_list(num_pose_kept, tag, pose_data_list, pose, current_score);
	 		}
		}
//		std::cout << "unique_bins= " << unique_bins;
//		std::cout << "  count= " << count << std::endl;

//  	Cluster_poses_test(num_pose_kept, pose_data_list, suite_num); //cluster pose one last time.





/////////////////Create 2 residue pose////////////////////
/*
	pose::Pose original_pose;
	std::string infile = option[ in ::file::s ][1];
	std::cout << infile << std::endl;
	core::import_pose::pose_from_pdb( original_pose, *rsd_set, infile );


	Size num_res= original_pose.total_residue();

	pose::Pose pose=original_pose;

	for (Size n = num_res; n > 1; n-- ) {
		std::cout << n << std::endl;
		pose.conformation().delete_residue_slow(n);

	}

		dump_pdb( pose, "before.pdb" );


		utility::vector1< std::pair< std::string, std::string> > atom_pairs;


//

		chemical::AA res_aa = aa_from_name( "RAD" );
		ResidueOP new_rsd = conformation::ResidueFactory::create_residue( *(rsd_set->aa_map( res_aa )[1]) ) ;
//		atom_pairs.push_back( std::make_pair( " O3*", " O3*" ) );
//		atom_pairs.push_back( std::make_pair( " C3*", " C3*" ) );
//		atom_pairs.push_back( std::make_pair( " C4*", " C4*" ) );
//		pose.replace_residue( 1, *new_rsd, atom_pairs );


		core::pose::remove_lower_terminus_type_from_pose_residue(pose, 1);	 //need this if prepend from end.
		pose.prepend_polymer_residue_before_seqpos( *new_rsd, 1, true);
		pose.conformation().delete_residue_slow(2);
//		core::pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_PHOSPHATE", 1 );
//		core::pose::remove_lower_terminus_type_from_pose_residue(pose, 1);	 //need this if prepend from end.
//		pose.prepend_polymer_residue_before_seqpos( *new_rsd, 1, true);
		pose.append_polymer_residue_after_seqpos( *new_rsd, 1, true);
//		core::pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_PHOSPHATE", 1 );


		dump_pdb( pose, "after.pdb" );

//		pose.append_residue_by_jump( *new_rsd, 1  );


		dump_pdb(pose, "two_res.pdb");
*/
/*
			Vector origin;
			origin[0]=0;
			origin[1]=0;
			origin[2]=0;


			std::cout << "x= " << origin[0] << std::endl;
			std::cout << "y= " << origin[1] << std::endl;
			std::cout << "z= " << origin[2] << std::endl;

			std::string const atom_name="C4*";

			pose.residue(1).set_xyz( atom_name , origin);
			conformation::Residue const & rsd( pose.residue(1) );

			pose.set_xyz( id::AtomID( rsd.atom_index( "C4*" ), 1 ) , origin);

			std::cout << "x= " << pose.residue(1).xyz(6)[0] << std::endl;
			std::cout << "y= " << pose.residue(1).xyz(6)[1] << std::endl;
			std::cout << "z= " << pose.residue(1).xyz(6)[2] << std::endl;


}

*/

/////////////////////////////////////////////////////////////////////////////////////////////
void
Output_data(pose::Pose& pose, std::string tag, Real& rmsd, core::io::silent::SilentFileData& silent_file_data){

	using namespace core::io::silent;
  using namespace core::scoring;
  using namespace basic::options;
	using namespace basic::options::OptionKeys;

	Size const write_score_only = option[print_score_only] ;

	std::string silent_file="long.txt";

  RNA_SilentStruct s(pose, tag );

	silent_file_data.write_silent_struct(s, silent_file, write_score_only);

}

void
Output_data_chi(pose::Pose& pose, std::string tag, Real& rmsd, core::io::silent::SilentFileData& silent_file_data, Real initial_chi_angle){

	using namespace core::io::silent;
  using namespace core::scoring;
  using namespace basic::options;
	using namespace basic::options::OptionKeys;

	Size const write_score_only = option[print_score_only] ;

	std::string silent_file="long.txt";

  RNA_SilentStruct s(pose, tag );

	conformation::Residue const & current_residue=pose.residue(1);
	std::string name="chi_angle";
	s.add_energy( name, current_residue.chi(1) );

	silent_file_data.write_silent_struct(s, silent_file, write_score_only);

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct rotamer_torsion_struct{
	std::string name;
	Real alpha;
	Real beta;
	Real gamma;
	Real delta;
	Real epsilon;
	Real zeta;
	Real chi;
	Real nu_1;
	Real nu_2;
};

struct rotamer_torsion_bin_struct{

	utility::vector1< Real > alpha;
	utility::vector1< Real > beta;
	utility::vector1< Real > gamma;
	utility::vector1< Real > delta;
	utility::vector1< Real > epsilon;
	utility::vector1< Real > zeta;
	utility::vector1< Real > chi;
	utility::vector1< Real > nu_1;
	utility::vector1< Real > nu_2;
};


Real
Get_torsion_from_bin_number(Size bin_number){
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		Real bin_size=option[torsion_bin];

		int num_bins=360/bin_size;
		int half_num_bins=180/bin_size;

		int int_bin_number=bin_number;

		Real torsion=((int_bin_number-half_num_bins)*bin_size)-(bin_size/2);



//	std::cout << "bin_number= " << bin_number << std::endl;

	return torsion;

}

Size
Get_bin_number(Real torsion_value){
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		Real bin_size=option[torsion_bin];

		int num_bins=360/bin_size;
		int half_num_bins=180/bin_size;

		Size bin_number;

		if(torsion_value < 0){
			bin_number=int(torsion_value/bin_size)+half_num_bins;
		}else {
			bin_number=int(torsion_value/bin_size)+half_num_bins+1;
		}
//	std::cout << "bin_number= " << bin_number << std::endl;

	return bin_number;

}

void
output_bin_data(rotamer_torsion_bin_struct& rotamer_torsion_bin, std::string name){

		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		Real bin_size=option[torsion_bin];
		int num_bins=360/bin_size;
		int half_num_bins=180/bin_size;

		std::string output_filename= name+"_torsion_bin_"+lead_zero_string_of(Size(bin_size), 2)+".txt";

		std::ofstream outfile;
	  outfile.open(output_filename.c_str());
		Size spacing=4;

//		outfile_column_name(outfile, current_rebuild_num);
	  for(Size i=1; i<=num_bins ; i++){

	  		outfile << std::setw(spacing+2) << std::setprecision(4) << std::left << i;
				outfile << std::setw(8) << std::setprecision(4) << std::left << Get_torsion_from_bin_number(i);
				outfile << std::setw(spacing+2) << std::setprecision(4) << std::left << rotamer_torsion_bin.alpha[i];
				outfile << std::setw(spacing+2) << std::setprecision(4) << std::left << rotamer_torsion_bin.beta[i];
 	 			outfile << std::setw(spacing+2) << std::setprecision(4) << std::left << rotamer_torsion_bin.gamma[i];
 	 			outfile << std::setw(spacing+2) << std::setprecision(4) << std::left << rotamer_torsion_bin.delta[i];
				outfile << std::setw(spacing+2) << std::setprecision(4) << std::left << rotamer_torsion_bin.epsilon[i];
				outfile << std::setw(spacing+2) << std::setprecision(4) << std::left << rotamer_torsion_bin.zeta[i];
 	 			outfile << std::setw(spacing+2) << std::setprecision(4) << std::left << rotamer_torsion_bin.chi[i];
				outfile << std::setw(spacing+2) << std::setprecision(4) << std::left << rotamer_torsion_bin.nu_2[i];
				outfile << std::setw(spacing+2) << std::setprecision(4) << std::left << rotamer_torsion_bin.nu_1[i];
 	 			outfile << "\n";

		}

		outfile.close();
}

void
output_torsion_data(utility::vector1< rotamer_torsion_struct >& rotamer_torsion_data, std::string name){
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		std::string output_filename= name+"_torsion_data.txt";

		std::ofstream outfile;
	  outfile.open(output_filename.c_str());
		Size spacing=8;

		outfile << std::setw(spacing+2) << std::setprecision(4) << std::left << "name";
		outfile << std::setw(spacing+2) << std::setprecision(4) << std::left << "alpha";
		outfile << std::setw(spacing+2) << std::setprecision(4) << std::left << "beta";
 	 	outfile << std::setw(spacing+2) << std::setprecision(4) << std::left << "gamma";
 	 	outfile << std::setw(spacing+2) << std::setprecision(4) << std::left << "delta";
		outfile << std::setw(spacing+2) << std::setprecision(4) << std::left << "epsilon";
		outfile << std::setw(spacing+2) << std::setprecision(4) << std::left << "zeta";
 	 	outfile << std::setw(spacing+2) << std::setprecision(4) << std::left << "chi";
		outfile << std::setw(spacing+2) << std::setprecision(4) << std::left << "nu_2";
		outfile << std::setw(spacing+2) << std::setprecision(4) << std::left << "nu_1";
 	 	outfile << "\n";


		for(Size i=1; i<=rotamer_torsion_data.size(); i++){

	  		outfile << std::setw(spacing+2) << std::fixed <<  std::setprecision(4) << std::left << rotamer_torsion_data[i].name;
				outfile << std::setw(spacing+2) << std::fixed <<  std::setprecision(4) << std::left << rotamer_torsion_data[i].alpha;
				outfile << std::setw(spacing+2) << std::fixed <<  std::setprecision(4) << std::left << rotamer_torsion_data[i].beta;
 	 			outfile << std::setw(spacing+2) << std::fixed <<  std::setprecision(4) << std::left << rotamer_torsion_data[i].gamma;
 	 			outfile << std::setw(spacing+2) << std::fixed <<  std::setprecision(4) << std::left << rotamer_torsion_data[i].delta;
				outfile << std::setw(spacing+2) << std::fixed <<  std::setprecision(4) << std::left << rotamer_torsion_data[i].epsilon;
				outfile << std::setw(spacing+2) << std::fixed <<  std::setprecision(4) << std::left << rotamer_torsion_data[i].zeta;
 	 			outfile << std::setw(spacing+2) << std::fixed <<  std::setprecision(4) << std::left << rotamer_torsion_data[i].chi;
				outfile << std::setw(spacing+2) << std::fixed <<  std::setprecision(4) << std::left << rotamer_torsion_data[i].nu_2;
				outfile << std::setw(spacing+2) << std::fixed <<  std::setprecision(4) << std::left << rotamer_torsion_data[i].nu_1;
 	 			outfile << "\n";

		}
	outfile.close();
}

void
Initialize_bin(rotamer_torsion_bin_struct& rotamer_torsion_bin){

		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		Real bin_size=option[torsion_bin];
		int num_bins=360/bin_size;
		int half_num_bins=180/bin_size;

		utility::vector1< Size > null_bin(num_bins , 0);

		rotamer_torsion_bin.alpha=null_bin;
		rotamer_torsion_bin.beta=null_bin;
		rotamer_torsion_bin.gamma=null_bin;
		rotamer_torsion_bin.delta=null_bin;
		rotamer_torsion_bin.epsilon=null_bin;
		rotamer_torsion_bin.zeta=null_bin;
		rotamer_torsion_bin.chi=null_bin;
		rotamer_torsion_bin.nu_2=null_bin;
		rotamer_torsion_bin.nu_1=null_bin;
}


void
Update_rotamer_bin(rotamer_torsion_struct& rotamer_torsion, rotamer_torsion_bin_struct& rotamer_torsion_bin){

		rotamer_torsion_bin.alpha[Get_bin_number(rotamer_torsion.alpha)]++;
		rotamer_torsion_bin.beta[Get_bin_number(rotamer_torsion.beta)]++;
		rotamer_torsion_bin.gamma[Get_bin_number(rotamer_torsion.gamma)]++;
		rotamer_torsion_bin.delta[Get_bin_number(rotamer_torsion.delta)]++;
		rotamer_torsion_bin.epsilon[Get_bin_number(rotamer_torsion.epsilon)]++;
		rotamer_torsion_bin.zeta[Get_bin_number(rotamer_torsion.zeta)]++;
		rotamer_torsion_bin.chi[Get_bin_number(rotamer_torsion.chi)]++;
		rotamer_torsion_bin.nu_2[Get_bin_number(rotamer_torsion.nu_2)]++;
		rotamer_torsion_bin.nu_1[Get_bin_number(rotamer_torsion.nu_1)]++;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
Torsion_analysis_test(){

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace core::kinematics;
	using namespace core::optimization;
	using namespace core::io::silent;
	using namespace core::id;
	using namespace protocols::rna;


			std::string filename="chemical/rna/1jj2.torsions";

			std::ifstream infile;
		 	infile.open(filename.c_str());
			if (infile.fail()) {
  				std::cout << "Error! \"" << filename << "\" could not be opened!" << std::endl;
  				infile.clear();
  				infile.close(); //do I need this line?
			}
			std::cout << "Open \"" << filename << "\" successful!" << std::endl;


		utility::vector1< rotamer_torsion_struct > rotamer_torsion_data;


		std::string data_string;
		while(getline(infile, data_string)) {

			std::vector<std::string> data_tokens;


			Tokenize(data_string, data_tokens);

//				std::cout << " name= " <<  data_tokens[0];
//				std::cout << " alpha= " << data_tokens[1];
//				std::cout << " beta= " << data_tokens[2];
//				std::cout << " gamma= " << data_tokens[3];
//				std::cout << " delta= " << data_tokens[4];
//				std::cout << " epsilon= " <<  data_tokens[5];
//				std::cout << " zeta= " << data_tokens[6];
//				std::cout << " chi= " << data_tokens[7];
//				std::cout << " nu_1= " << data_tokens[8];
//				std::cout << " nu_2= " << data_tokens[9];
//				std::cout<< std::endl;

				rotamer_torsion_struct rotamer_torsion;
				rotamer_torsion.name = data_tokens[0];
			  rotamer_torsion.alpha=convert_string_to_real(data_tokens[1]);
			 	rotamer_torsion.beta=convert_string_to_real(data_tokens[2]);
			  rotamer_torsion.gamma=convert_string_to_real(data_tokens[3]);
			  rotamer_torsion.delta=convert_string_to_real(data_tokens[4]);
			  rotamer_torsion.epsilon=convert_string_to_real(data_tokens[5]);
			  rotamer_torsion.zeta=convert_string_to_real(data_tokens[6]);
				rotamer_torsion.chi=convert_string_to_real(data_tokens[7]);
				rotamer_torsion.nu_2=convert_string_to_real(data_tokens[8]);
				rotamer_torsion.nu_1=convert_string_to_real(data_tokens[9]);
				rotamer_torsion_data.push_back(rotamer_torsion);
		}




		Real bin_size=option[torsion_bin];
		int num_bins=360/bin_size;
		int half_num_bins=180/bin_size;

		std::cout << "bin_size " << bin_size << std::endl;
		std::cout << "num_bins " << num_bins << std::endl;
		std::cout << "half_num_bins " << half_num_bins << std::endl;

		Real number =987.654321;
		Real number2=0.41312312;

		std::cout << std::fixed << std::setprecision(1) << number << std::endl;
		std::cout << std::fixed << std::setprecision(2) << number << std::endl;
		std::cout << std::fixed << std::setprecision(3) << number << std::endl;
		std::cout << std::fixed << std::setprecision(4) << number << std::endl;
		std::cout << std::fixed << std::setprecision(5) << number << std::endl;
		std::cout << std::fixed << std::setprecision(6) << number << std::endl;

		std::cout << std::fixed << std::setprecision(1) << number2 << std::endl;
		std::cout << std::fixed << std::setprecision(2) << number2 << std::endl;
		std::cout << std::fixed << std::setprecision(3) << number2 << std::endl;
		std::cout << std::fixed << std::setprecision(4) << number2 << std::endl;
		std::cout << std::fixed << std::setprecision(5) << number2 << std::endl;
		std::cout << std::fixed << std::setprecision(6) << number2 << std::endl;


		Size zero =0;

		std::cout << "Size zero = " << zero << std::endl;

		for(int i=1; i<=num_bins; i++){

//			std::cout << "(i-half_num_bins)=" << (i-half_num_bins) << std::endl;
//			std::cout << "((i-half_num_bins)*bin_size)=" << ((i-half_num_bins)*bin_size) << std::endl;

			Real torsion_value=((i-half_num_bins)*bin_size)-(bin_size/2);
			std::cout << "The bin of torsion: " << torsion_value << " is " << Get_bin_number(torsion_value) << std::endl;
		}

		Real sum_square_nu2=0;
		Real sum_nu2=0;
		Real count=0;
		Real sum_square_nu1=0;
		Real sum_nu1=0;
		Real sum_square_chi=0;
		Real sum_chi=0;
		Real sum_square_delta=0;
		Real sum_delta=0;


		rotamer_torsion_bin_struct rotamer_torsion_bin;
		Initialize_bin(rotamer_torsion_bin);

		for(Size i=1; i<=rotamer_torsion_data.size(); i++){
				Update_rotamer_bin(rotamer_torsion_data[i], rotamer_torsion_bin);
		}

		output_bin_data(rotamer_torsion_bin, "All");


		//Define syn chi angle cutoff as being in the interval [0;-100]
		//Define 3'endo as having delta angle outside the interval [0:115]

		rotamer_torsion_bin_struct syn_three_prime_rotamer_torsion_bin;
		Initialize_bin(syn_three_prime_rotamer_torsion_bin);
		utility::vector1< rotamer_torsion_struct > syn_three_prime_rotamer_torsion_data;

		sum_square_nu2=0;
		sum_nu2=0;
		count=0;
		sum_square_nu1=0;
		sum_nu1=0;
		sum_square_chi=0;
		sum_chi=0;
		sum_square_delta=0;
		sum_delta=0;


		for(Size i=1; i<=rotamer_torsion_data.size(); i++){

			if((rotamer_torsion_data[i].chi>-100) && (rotamer_torsion_data[i].chi<0) && (rotamer_torsion_data[i].delta<115) && (rotamer_torsion_data[i].delta>0.0)){
				count++;
				sum_nu2=sum_nu2+rotamer_torsion_data[i].nu_2;
				sum_nu1=sum_nu1+rotamer_torsion_data[i].nu_1;

				sum_square_nu2=sum_square_nu2+(rotamer_torsion_data[i].nu_2*rotamer_torsion_data[i].nu_2);
				sum_square_nu1=sum_square_nu1+(rotamer_torsion_data[i].nu_1*rotamer_torsion_data[i].nu_1);

				sum_chi=sum_chi+rotamer_torsion_data[i].chi;
				sum_square_chi=sum_square_chi+(rotamer_torsion_data[i].chi*rotamer_torsion_data[i].chi);

				sum_delta=sum_delta+rotamer_torsion_data[i].delta;
				sum_square_delta=sum_square_delta+(rotamer_torsion_data[i].delta*rotamer_torsion_data[i].delta);

				syn_three_prime_rotamer_torsion_data.push_back(rotamer_torsion_data[i]);
				Update_rotamer_bin(rotamer_torsion_data[i], syn_three_prime_rotamer_torsion_bin);

			}
		}

		std::cout << "Syn_3prime  " << "count " << count << std::endl;
		std::cout << "mean nu_1= " << (sum_nu1/count);
		std::cout << "  mean_square nu_1= " << (sum_square_nu1/count);
		std::cout << "  standard_deviation_nu_1 " << sqrt( (sum_square_nu1/count)-( (sum_nu1/count)*(sum_nu1/count) ) );
		std::cout << std::endl;
		std::cout << "mean nu_2= " << (sum_nu2/count);
		std::cout << "  mean_square nu_2= " << (sum_square_nu2/count);
		std::cout << "  standard_deviation_nu_2 " << sqrt( (sum_square_nu2/count)-( (sum_nu2/count)*(sum_nu2/count) ) );
		std::cout << std::endl;
		std::cout << "mean chi= " << (sum_chi/count);
		std::cout << "  mean_square chi= " << (sum_square_chi/count);
		std::cout << "  standard_deviation_chi " << sqrt( (sum_square_chi/count)-( (sum_chi/count)*(sum_chi/count) ) );
		std::cout << std::endl;
		std::cout << "mean delta= " << (sum_delta/count);
		std::cout << "  mean_square delta= " << (sum_square_delta/count);
		std::cout << "  standard_deviation_delta " << sqrt( (sum_square_delta/count)-( (sum_delta/count)*(sum_delta/count) ) );
		std::cout << std::endl;


		output_bin_data(syn_three_prime_rotamer_torsion_bin, "Syn_3prime");
		output_torsion_data(syn_three_prime_rotamer_torsion_data, "Syn_3prime");

		//Define syn chi angle cutoff as being in the interval [0;-100]
		//Define 2'endo as having delta angle in the interval [115:180] ...together there will be few points that will not be classified as both 2' and 3' endo

		rotamer_torsion_bin_struct syn_two_prime_rotamer_torsion_bin;
		Initialize_bin(syn_two_prime_rotamer_torsion_bin);
		utility::vector1< rotamer_torsion_struct > syn_two_prime_rotamer_torsion_data;

		sum_square_nu2=0;
		sum_nu2=0;
		count=0;
		sum_square_nu1=0;
		sum_nu1=0;
		sum_square_chi=0;
		sum_chi=0;
		sum_square_delta=0;
		sum_delta=0;

		for(Size i=1; i<=rotamer_torsion_data.size(); i++){

			if((rotamer_torsion_data[i].chi>-100) && (rotamer_torsion_data[i].chi<0) && (rotamer_torsion_data[i].delta>115) && (rotamer_torsion_data[i].delta<180)){
				count++;
				sum_nu2=sum_nu2+rotamer_torsion_data[i].nu_2;
				sum_nu1=sum_nu1+rotamer_torsion_data[i].nu_1;

				sum_square_nu2=sum_square_nu2+(rotamer_torsion_data[i].nu_2*rotamer_torsion_data[i].nu_2);
				sum_square_nu1=sum_square_nu1+(rotamer_torsion_data[i].nu_1*rotamer_torsion_data[i].nu_1);

				sum_chi=sum_chi+rotamer_torsion_data[i].chi;
				sum_square_chi=sum_square_chi+(rotamer_torsion_data[i].chi*rotamer_torsion_data[i].chi);

				sum_delta=sum_delta+rotamer_torsion_data[i].delta;
				sum_square_delta=sum_square_delta+(rotamer_torsion_data[i].delta*rotamer_torsion_data[i].delta);

				syn_two_prime_rotamer_torsion_data.push_back(rotamer_torsion_data[i]);
				Update_rotamer_bin(rotamer_torsion_data[i], syn_two_prime_rotamer_torsion_bin);

			}
		}

		std::cout << "Syn_2prime  " << "count " << count << std::endl;
		std::cout << "mean nu_1= " << (sum_nu1/count);
		std::cout << "  mean_square nu_1= " << (sum_square_nu1/count);
		std::cout << "  standard_deviation_nu_1 " << sqrt( (sum_square_nu1/count)-( (sum_nu1/count)*(sum_nu1/count) ) );
		std::cout << std::endl;
		std::cout << "mean nu_2= " << (sum_nu2/count);
		std::cout << "  mean_square nu_2= " << (sum_square_nu2/count);
		std::cout << "  standard_deviation_nu_2 " << sqrt( (sum_square_nu2/count)-( (sum_nu2/count)*(sum_nu2/count) ) );
		std::cout << std::endl;
		std::cout << "mean chi= " << (sum_chi/count);
		std::cout << "  mean_square chi= " << (sum_square_chi/count);
		std::cout << "  standard_deviation_chi " << sqrt( (sum_square_chi/count)-( (sum_chi/count)*(sum_chi/count) ) );
		std::cout << std::endl;
		std::cout << "mean delta= " << (sum_delta/count);
		std::cout << "  mean_square delta= " << (sum_square_delta/count);
		std::cout << "  standard_deviation_delta " << sqrt( (sum_square_delta/count)-( (sum_delta/count)*(sum_delta/count) ) );
		std::cout << std::endl;

		output_bin_data(syn_two_prime_rotamer_torsion_bin, "Syn_2prime");
		output_torsion_data(syn_two_prime_rotamer_torsion_data, "Syn_2prime");


		//three_endo test

		rotamer_torsion_bin_struct three_prime_rotamer_torsion_bin;
		Initialize_bin(three_prime_rotamer_torsion_bin);

		sum_square_nu2=0;
		sum_nu2=0;
		count=0;
		sum_square_nu1=0;
		sum_nu1=0;
		sum_square_chi=0;
		sum_chi=0;
		sum_square_delta=0;
		sum_delta=0;


		for(Size i=1; i<=rotamer_torsion_data.size(); i++){

			if((rotamer_torsion_data[i].delta<115.0) && (rotamer_torsion_data[i].delta>0.0)){
				count++;

				sum_nu2=sum_nu2+rotamer_torsion_data[i].nu_2;
				sum_nu1=sum_nu1+rotamer_torsion_data[i].nu_1;

				sum_square_nu2=sum_square_nu2+(rotamer_torsion_data[i].nu_2*rotamer_torsion_data[i].nu_2);
				sum_square_nu1=sum_square_nu1+(rotamer_torsion_data[i].nu_1*rotamer_torsion_data[i].nu_1);

				sum_chi=sum_chi+rotamer_torsion_data[i].chi;
				sum_square_chi=sum_square_chi+(rotamer_torsion_data[i].chi*rotamer_torsion_data[i].chi);

				sum_delta=sum_delta+rotamer_torsion_data[i].delta;
				sum_square_delta=sum_square_delta+(rotamer_torsion_data[i].delta*rotamer_torsion_data[i].delta);

				Update_rotamer_bin(rotamer_torsion_data[i], three_prime_rotamer_torsion_bin);
			}
		}

		std::cout << "3prime  " << "count " << count << std::endl;
		std::cout << "mean nu_1= " << (sum_nu1/count);
		std::cout << "  mean_square nu_1= " << (sum_square_nu1/count);
		std::cout << "  standard_deviation_nu_1 " << sqrt( (sum_square_nu1/count)-( (sum_nu1/count)*(sum_nu1/count) ) );
		std::cout << std::endl;
		std::cout << "mean nu_2= " << (sum_nu2/count);
		std::cout << "  mean_square nu_2= " << (sum_square_nu2/count);
		std::cout << "  standard_deviation_nu_2 " << sqrt( (sum_square_nu2/count)-( (sum_nu2/count)*(sum_nu2/count) ) );
		std::cout << std::endl;
		std::cout << "mean chi= " << (sum_chi/count);
		std::cout << "  mean_square chi= " << (sum_square_chi/count);
		std::cout << "  standard_deviation_chi " << sqrt( (sum_square_chi/count)-( (sum_chi/count)*(sum_chi/count) ) );
		std::cout << std::endl;
		std::cout << "mean delta= " << (sum_delta/count);
		std::cout << "  mean_square delta= " << (sum_square_delta/count);
		std::cout << "  standard_deviation_delta " << sqrt( (sum_square_delta/count)-( (sum_delta/count)*(sum_delta/count) ) );
		std::cout << std::endl;


		output_bin_data(three_prime_rotamer_torsion_bin, "3prime");

		//two_endo test

		rotamer_torsion_bin_struct two_prime_rotamer_torsion_bin;
		Initialize_bin(two_prime_rotamer_torsion_bin);

		sum_square_nu2=0;
		sum_nu2=0;
		count=0;
		sum_square_nu1=0;
		sum_nu1=0;
		sum_square_chi=0;
		sum_chi=0;
		sum_square_delta=0;
		sum_delta=0;

		for(Size i=1; i<=rotamer_torsion_data.size(); i++){

			if((rotamer_torsion_data[i].delta>115.0) && (rotamer_torsion_data[i].delta<180.0)){
				count++;

				sum_nu2=sum_nu2+rotamer_torsion_data[i].nu_2;
				sum_nu1=sum_nu1+rotamer_torsion_data[i].nu_1;

				sum_square_nu2=sum_square_nu2+(rotamer_torsion_data[i].nu_2*rotamer_torsion_data[i].nu_2);
				sum_square_nu1=sum_square_nu1+(rotamer_torsion_data[i].nu_1*rotamer_torsion_data[i].nu_1);

				sum_chi=sum_chi+rotamer_torsion_data[i].chi;
				sum_square_chi=sum_square_chi+(rotamer_torsion_data[i].chi*rotamer_torsion_data[i].chi);

				sum_delta=sum_delta+rotamer_torsion_data[i].delta;
				sum_square_delta=sum_square_delta+(rotamer_torsion_data[i].delta*rotamer_torsion_data[i].delta);

				Update_rotamer_bin(rotamer_torsion_data[i], two_prime_rotamer_torsion_bin);
			}
		}

		std::cout << "2prime  " << "count " << count << std::endl;
		std::cout << "mean nu_1= " << (sum_nu1/count);
		std::cout << "  mean_square nu_1= " << (sum_square_nu1/count);
		std::cout << "  standard_deviation_nu_1 " << sqrt( (sum_square_nu1/count)-( (sum_nu1/count)*(sum_nu1/count) ) );
		std::cout << std::endl;
		std::cout << "mean nu_2= " << (sum_nu2/count);
		std::cout << "  mean_square nu_2= " << (sum_square_nu2/count);
		std::cout << "  standard_deviation_nu_2 " << sqrt( (sum_square_nu2/count)-( (sum_nu2/count)*(sum_nu2/count) ) );
		std::cout << std::endl;
		std::cout << "mean chi= " << (sum_chi/count);
		std::cout << "  mean_square chi= " << (sum_square_chi/count);
		std::cout << "  standard_deviation_chi " << sqrt( (sum_square_chi/count)-( (sum_chi/count)*(sum_chi/count) ) );
		std::cout << std::endl;
		std::cout << "mean delta= " << (sum_delta/count);
		std::cout << "  mean_square delta= " << (sum_square_delta/count);
		std::cout << "  standard_deviation_delta " << sqrt( (sum_square_delta/count)-( (sum_delta/count)*(sum_delta/count) ) );
		std::cout << std::endl;

		output_bin_data(two_prime_rotamer_torsion_bin, "2prime");

}


///////////////////////////////////////////////////////////////////////////////////////////////////
void
Chi_angle_test(){

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace core::kinematics;
	using namespace core::optimization;
	using namespace core::io::silent;
	using namespace core::id;
	using namespace protocols::rna;
	using namespace core::scoring::rna;

  ///////////////// Declare type of residue as RNA type//////////////////////////////////////

  ResidueTypeSetCAP rsd_set; //This line is repeated
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "rna" ); //This line is repeated

	ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( RNA_HIRES_WTS );
	if ( option[ score::weights ].user() ) {
		scorefxn = ScoreFunctionFactory::create_score_function( option[ score::weights] );
	}


  //Might need this for decoy pose, at chain break.
  if(option[minimize_and_score_sugar]==true){
  	scorefxn->set_weight( rna_sugar_close, 0.7 );
  } else {
    scorefxn->set_weight( rna_sugar_close, 0.0 );
 	}



	AtomTreeMinimizer minimizer;
	float const dummy_tol( 0.0000025);
	bool const use_nblist( true );
	MinimizerOptions options( "dfpmin", dummy_tol, use_nblist, false, false );
	options.nblist_auto_update( true ); //update neighor atom position
	kinematics::MoveMap mm;   //Set whether the torsion angles of the residue is freeze or movable during minimization
	mm.set_bb( true );
	mm.set_chi( true );



	SilentFileData silent_file_data;  //This setup the silent files//
  ///////////////Read in pdb data to create original_pose ////////////////////////////////////

	pose::Pose pose;
	std::string tag = option[ in ::file::s ][1];

	std::cout << tag+".pdb" << std::endl;
	core::import_pose::pose_from_pdb( pose, *rsd_set, tag+".pdb");

  o2star_minimize(pose, scorefxn);

	if(option[minimize_and_score_sugar]==true){
//		UnFreeze_sugar_torsions(mm, pose);
		scorefxn->set_weight( rna_sugar_close, 0.7 );
  } else {
  	Freeze_sugar_torsions(mm, pose);
    scorefxn->set_weight( rna_sugar_close, 0.0 );
 	}

	Output_fold_tree_info(pose, "Original_pose");
 	dump_pdb( pose, "imported_pose.pdb");

	std::string pucker_choice = option[pucker];
	std::string base_choice=option[base];

	tag=pucker_choice+"_"+base_choice;

	if((base_choice!="RAD") && (base_choice!="RCY") && (base_choice!="URA")	&& (base_choice!="RGU")){
		std::cout << "Invalid sugar pucker choice" << std::endl;
		exit (1);
	}

	chemical::AA res_aa = aa_from_name( base_choice );
	std::cout << "  res_aa: " << res_aa << std::endl;


	ResidueOP new_rsd = conformation::ResidueFactory::create_residue( *(rsd_set->aa_map( res_aa )[1]) ) ;
	core::pose::remove_upper_terminus_type_from_pose_residue(pose, 1);
	pose.append_polymer_residue_after_seqpos( *new_rsd, 1, true);
	Size which_res=2;


	//Set to A form
	scoring::rna::RNA_TorsionPotential const rna_torsion_potential;
	pose.set_torsion( TorsionID( which_res, id::BB, 1 ), rna_torsion_potential.gaussian_parameter_set_alpha()[1].center );
	pose.set_torsion( TorsionID( which_res, id::BB, 2 ), rna_torsion_potential.gaussian_parameter_set_beta()[1].center );
	pose.set_torsion( TorsionID( which_res, id::BB, 3 ), rna_torsion_potential.gaussian_parameter_set_gamma()[1].center );
	pose.set_torsion( TorsionID( which_res-1, id::BB, 6 ), rna_torsion_potential.gaussian_parameter_set_zeta_alpha_sc_minus()[1].center );
	pose.set_torsion( TorsionID( which_res, id::CHI, 4 ), 0.0 );

 //nu_2=31.169
 //nu_1=110.123  Parameter from G of UUCG tetraloop

 //Statistical analysis of chemical/rna/1jj2.torsions
//  Syn_3prime  count 21
//	mean nu_1= 107.792  mean_square nu_1= 11655.6  standard_deviation_nu_1 6.04656
//  mean nu_2= 31.8163  mean_square nu_2= 1026.34  standard_deviation_nu_2 3.75109
//  3prime  count 2352
//  mean nu_1= 93.8932  mean_square nu_1= 8847.93  standard_deviation_nu_1 5.65716
//  mean nu_2= 38.1739  mean_square nu_2= 1468.51  standard_deviation_nu_2 3.35621
//  Syn_2prime  count 27
//  mean nu_1= 157.725  mean_square nu_1= 24903.2  standard_deviation_nu_1 5.09893
//  mean nu_2= -37.4744  mean_square nu_2= 1408.15  standard_deviation_nu_2 1.95549
//  2prime  count 394
//  mean nu_1= 156.316  mean_square nu_1= 24479  standard_deviation_nu_1 6.66222
//  mean nu_2= -36.5452  mean_square nu_2= 1353.85  standard_deviation_nu_2 4.27825



	if(pucker_choice=="Three_prime_endo"){
		pose.set_torsion( TorsionID( which_res-1, id::BB, 5 ), rna_torsion_potential.gaussian_parameter_set_epsilon_north()[1].center );
		pose.set_torsion( TorsionID( which_res, id::BB, 4 ), rna_torsion_potential.gaussian_parameter_set_delta_north()[1].center );
		pose.set_torsion( TorsionID( which_res, id::CHI, 1 ), rna_torsion_potential.gaussian_parameter_set_chi_north()[1].center );
//	pose.set_torsion( TorsionID( which_res, id::CHI, 2 ), rna_torsion_potential.gaussian_parameter_set_nu2_north()[1].center );
//	pose.set_torsion( TorsionID( which_res, id::CHI, 3 ), rna_torsion_potential.gaussian_parameter_set_nu1_north()[1].center );
		pose.set_torsion( TorsionID( which_res, id::CHI, 2 ), 31.8163 );
		pose.set_torsion( TorsionID( which_res, id::CHI, 3 ), 107.792 );
//		pose.set_torsion( TorsionID( which_res, id::CHI, 2 ), 31.54 );
//		pose.set_torsion( TorsionID( which_res, id::CHI, 3 ), 110.181 );


	} else if (pucker_choice=="Two_prime_endo"){

		apply_ideal_c2endo_sugar_coords( pose, 2 );
		pose.set_torsion( TorsionID( which_res, id::CHI, 4 ), 0.0 );
		pose.set_torsion( TorsionID( which_res-1, id::BB, 5 ), rna_torsion_potential.gaussian_parameter_set_epsilon_south()[1].center );
		pose.set_torsion( TorsionID( which_res, id::BB, 4 ), rna_torsion_potential.gaussian_parameter_set_delta_south()[1].center );
		pose.set_torsion( TorsionID( which_res, id::CHI, 1 ), rna_torsion_potential.gaussian_parameter_set_chi_south()[1].center );
//	pose.set_torsion( TorsionID( which_res, id::CHI, 2 ), rna_torsion_potential.gaussian_parameter_set_nu2_south()[1].center );
//	pose.set_torsion( TorsionID( which_res, id::CHI, 3 ), rna_torsion_potential.gaussian_parameter_set_nu1_south()[1].center );

		pose.set_torsion( TorsionID( which_res, id::CHI, 2 ), -37.4744 ); //Barely change compare to anti chi conformation case
		pose.set_torsion( TorsionID( which_res, id::CHI, 3 ), 157.725 );  //Barely change compare to anti chi conformation case
//c   -75.539327 -143.310093  -53.762169  155.328455 -101.175338   93.010944  -56.725585   -40.435535   161.391854    67.500550 0   319 (from 1JJ2)

//	pose.set_torsion( TorsionID( which_res, id::CHI, 2 ), -40.435535 );
//	pose.set_torsion( TorsionID( which_res, id::CHI, 3 ), 161.391854 );

//c   -51.413435  164.895453   51.203012  145.026688  -88.965989  136.661091  -44.211250   -37.908931   166.468108  -292.499450 0  2335 (from 1JJ2)

//	pose.set_torsion( TorsionID( which_res, id::CHI, 2 ), -37.908931  );
//	pose.set_torsion( TorsionID( which_res, id::CHI, 3 ), 166.468108  );


	} else {
		std::cout << "Invalid sugar pucker choice" << std::endl;
		exit (1);
	}


	for(Size i=0; i<=360; i++){

		Real chi_angle=i*1;
		std::cout << "chi_angle= " << chi_angle << std::endl;

		pose.set_torsion( TorsionID( 2, id::CHI, 1 ), chi_angle);//chi i+1

		pose::Pose new_pose=pose;

		core::pose::remove_lower_terminus_type_from_pose_residue(new_pose, 1);
		new_pose.conformation().delete_residue_slow(1);
		core::pose::add_lower_terminus_type_to_pose_residue(new_pose, 1);
		core::pose::add_variant_type_to_pose_residue( new_pose, "VIRTUAL_PHOSPHATE", 1 );

		(*scorefxn)(new_pose);
		scorefxn->show( std::cout, new_pose);

		Real rmsd=all_atom_rmsd( new_pose, new_pose);

		//Need this due to some scoring bug....discuss this will Rhiju tmr.
		if(i!=0){
			Output_data_chi(new_pose, tag+"_"+lead_zero_string_of(chi_angle, 3), rmsd, silent_file_data, chi_angle);
			dump_pdb( new_pose, tag+"_"+lead_zero_string_of(chi_angle, 3)+".pdb");
		}

	}


}


//Bug Occur only if
//1. o2star trail is turned on
//2. if choose pucker choice Two_prime_endo

void
RNA_torsion_bug(){

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace core::kinematics;
	using namespace core::optimization;
	using namespace core::io::silent;
	using namespace core::id;
	using namespace protocols::rna;
	using namespace core::scoring::rna;

  ///////////////// Declare type of residue as RNA type//////////////////////////////////////

  ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "rna" );

	ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( RNA_HIRES_WTS );
	SilentFileData silent_file_data;  //This setup the silent files//

	pose::Pose pose;
	std::string tag = option[ in ::file::s ][1];
	std::cout << tag+".pdb" << std::endl;
	core::import_pose::pose_from_pdb( pose, *rsd_set, tag+".pdb");

  scorefxn->set_weight( rna_sugar_close, 0.0 );


	///////////o2star trail/////////////Bug dissapear if this is turned off.//////////////////////////////////////////////
	pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
	task->initialize_from_command_line();

	for (Size i = 1; i <= pose.total_residue(); i++) {
		if ( !pose.residue(i).is_RNA() ) continue;
		task->nonconst_residue_task(i).and_extrachi_cutoff( 0 );
		task->nonconst_residue_task(i).or_include_current( true );
	}
	pack::rotamer_trials( pose, *scorefxn, task);
/////////////////////////////////////////////////////////////////////////////////////

	std::string pucker_choice = option[pucker];

	tag=pucker_choice;

	chemical::AA res_aa = aa_from_name("RAD"); //Bug occur regardless of base choice.
	std::cout << "Appending residue" << "  res_aa: " << res_aa << std::endl;
	ResidueOP new_rsd = conformation::ResidueFactory::create_residue( *(rsd_set->aa_map( res_aa )[1]) ) ;
	core::pose::remove_upper_terminus_type_from_pose_residue(pose, 1);
	pose.append_polymer_residue_after_seqpos( *new_rsd, 1, true);
	core::pose::add_upper_terminus_type_to_pose_residue( pose, 2);

	Size which_res=2;

	std::cout << "Applying torsion angle "<< std::endl;

	//Set to A form
	scoring::rna::RNA_TorsionPotential const rna_torsion_potential;
	pose.set_torsion( TorsionID( which_res, id::BB, 1 ), rna_torsion_potential.gaussian_parameter_set_alpha()[1].center );
	pose.set_torsion( TorsionID( which_res, id::BB, 2 ), rna_torsion_potential.gaussian_parameter_set_beta()[1].center );
	pose.set_torsion( TorsionID( which_res, id::BB, 3 ), rna_torsion_potential.gaussian_parameter_set_gamma()[1].center );
	pose.set_torsion( TorsionID( which_res-1, id::BB, 6 ), rna_torsion_potential.gaussian_parameter_set_zeta_alpha_sc_minus()[1].center );
	pose.set_torsion( TorsionID( which_res-1, id::BB, 5 ), rna_torsion_potential.gaussian_parameter_set_epsilon_north()[1].center ); //Previous residue is 3' endo
	pose.set_torsion( TorsionID( which_res, id::CHI, 4 ), 0.0 );


	if(pucker_choice=="Three_prime_endo"){
		pose.set_torsion( TorsionID( which_res, id::BB, 4 ), rna_torsion_potential.gaussian_parameter_set_delta_north()[1].center );
		pose.set_torsion( TorsionID( which_res, id::CHI, 1 ), rna_torsion_potential.gaussian_parameter_set_chi_north()[1].center );
		pose.set_torsion( TorsionID( which_res, id::CHI, 2 ), rna_torsion_potential.gaussian_parameter_set_nu2_north()[1].center );
		pose.set_torsion( TorsionID( which_res, id::CHI, 3 ), rna_torsion_potential.gaussian_parameter_set_nu1_north()[1].center );

	} else if (pucker_choice=="Two_prime_endo"){ //Bug if choose Two_prime_endo with o2star_minimize(pose, scorefxn);

		pose.set_torsion( TorsionID( which_res, id::BB, 4 ), rna_torsion_potential.gaussian_parameter_set_delta_south()[1].center );
		pose.set_torsion( TorsionID( which_res, id::CHI, 1 ), rna_torsion_potential.gaussian_parameter_set_chi_south()[1].center );
		pose.set_torsion( TorsionID( which_res, id::CHI, 2 ), rna_torsion_potential.gaussian_parameter_set_nu2_south()[1].center );
		pose.set_torsion( TorsionID( which_res, id::CHI, 3 ), rna_torsion_potential.gaussian_parameter_set_nu1_south()[1].center );
	} else {
		std::cout << "Invalid sugar pucker choice" << std::endl;
		exit (1);
	}


/*
	for(Size i=0; i<=5; i++){

		Real chi_angle=i*1;
		std::cout << "chi_angle= " << chi_angle << std::endl;

		pose.set_torsion( TorsionID( 2, id::CHI, 1 ), chi_angle);//chi i+1

		std::cout << "check point 1 " << std::endl;

//		pose::Pose new_pose=pose;

//		core::pose::remove_lower_terminus_type_from_pose_residue(new_pose, 1);
//		new_pose.conformation().delete_residue_slow(1);
//		core::pose::add_lower_terminus_type_to_pose_residue(new_pose, 1);
//		core::pose::add_variant_type_to_pose_residue( new_pose, "VIRTUAL_PHOSPHATE", 1 );

		(*scorefxn)(pose);
		std::cout << "check point 2 " << std::endl;

		scorefxn->show( std::cout, pose);

		Real rmsd=0;

		Output_data_chi(pose, tag+"_"+lead_zero_string_of(chi_angle, 3), rmsd, silent_file_data, chi_angle);
		dump_pdb( pose, tag+"_"+lead_zero_string_of(chi_angle, 3)+".pdb");

	}
*/



	for(Size i=0; i<=0; i++){

//		std::cout << " FIRST TRY (*scorefxn)(pose) "<< std::endl;
//		(*scorefxn)(pose);
//		scorefxn->show( std::cout, pose);

		Real chi_angle=i*1;
		std::cout << "chi_angle= " << chi_angle << std::endl;

		pose.set_torsion( TorsionID( 2, id::CHI, 1 ), chi_angle);//chi i+1

		pose::Pose new_pose=pose;

		std::cout << "removing first residue" << std::endl;
		core::pose::remove_lower_terminus_type_from_pose_residue(new_pose, 1);
		new_pose.conformation().delete_residue_slow(1);
		core::pose::add_lower_terminus_type_to_pose_residue(new_pose, 1);
		core::pose::add_variant_type_to_pose_residue( new_pose, "VIRTUAL_PHOSPHATE", 1 );

		std::cout << "(*scorefxn)(new_pose) "<< std::endl;
		(*scorefxn)(new_pose);
		std::cout << "scorefxn->show( std::cout, new_pose); "<< std::endl;
		scorefxn->show( std::cout, new_pose);

		Real rmsd=0;

		std::cout << "Ouput_data_chi "<< std::endl;
		Output_data_chi(new_pose, tag+"_"+lead_zero_string_of(chi_angle, 3), rmsd, silent_file_data, chi_angle);
		dump_pdb( new_pose, tag+"_"+lead_zero_string_of(chi_angle, 3)+".pdb");

	}


}

///////////////////////////////////////////////////////////////////////////////////////////////////

void
Full_minimize(){

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace core::kinematics;
	using namespace core::optimization;
	using namespace core::io::silent;
	using namespace core::id;
	using namespace protocols::rna;


  ///////////////// Declare type of residue as RNA type//////////////////////////////////////

  ResidueTypeSetCAP rsd_set; //This line is repeated
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "rna" ); //This line is repeated

	ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( RNA_HIRES_WTS );
	if ( option[ score::weights ].user() ) {
		scorefxn = ScoreFunctionFactory::create_score_function( option[ score::weights] );
	}


  //Might need this for decoy pose, at chain break.
  if(option[minimize_and_score_sugar]==true){
  	scorefxn->set_weight( rna_sugar_close, 0.7 );
  } else {
    scorefxn->set_weight( rna_sugar_close, 0.0 );
 	}



	AtomTreeMinimizer minimizer;
	float const dummy_tol( 0.0000025);
	bool const use_nblist( true );
	MinimizerOptions options( "dfpmin", dummy_tol, use_nblist, false, false );
	options.nblist_auto_update( true ); //update neighor atom position
	kinematics::MoveMap mm;   //Set whether the torsion angles of the residue is freeze or movable during minimization
	mm.set_bb( true );
	mm.set_chi( true );



	SilentFileData silent_file_data;  //This setup the silent files//
  ///////////////Read in pdb data to create original_pose ////////////////////////////////////

	pose::Pose pose;
	std::string tag = option[ in ::file::s ][1];

	std::cout << tag+".pdb" << std::endl;
	core::import_pose::pose_from_pdb( pose, *rsd_set, tag+".pdb");

  o2star_minimize(pose, scorefxn);

	if(option[minimize_and_score_sugar]==true){
//		UnFreeze_sugar_torsions(mm, pose);
		scorefxn->set_weight( rna_sugar_close, 0.7 );
  } else {
  	Freeze_sugar_torsions(mm, pose);
    scorefxn->set_weight( rna_sugar_close, 0.0 );
 	}


	Output_fold_tree_info(pose, "Original_pose");
 	Real rmsd=all_atom_rmsd( pose, pose);
	(*scorefxn)(pose);
  scorefxn->show( std::cout, pose);
	Output_data(pose, tag, rmsd, silent_file_data);
	dump_pdb( pose, "Input_"+tag);


/*
	conformation::Residue const & current_residue=pose.residue(2);
	std::cout << "input chi= " << current_residue.chi(1) << std::endl; //chi_O2star torsion


	for(Size i=1; i<=360; i++){

		Real chi_angle=i*1;
		std::cout << "chi_angle= " << chi_angle << std::endl;

		pose.set_torsion( TorsionID( 2, id::CHI, 1 ), chi_angle);//chi i+1

		pose::Pose new_pose=pose;

//		std::cout << "Check point 1 "<< std::endl;
		core::pose::remove_lower_terminus_type_from_pose_residue(new_pose, 1);
//		std::cout << "Check point 2 "<< std::endl;
		new_pose.conformation().delete_residue_slow(1);
//		std::cout << "Check point 3 "<< std::endl;
		core::pose::add_lower_terminus_type_to_pose_residue(new_pose, 1);
//		std::cout << "Check point 4 "<< std::endl;
		core::pose::add_variant_type_to_pose_residue( new_pose, "VIRTUAL_PHOSPHATE", 1 );

		(*scorefxn)(new_pose);
  	scorefxn->show( std::cout, new_pose);

		Real rmsd=all_atom_rmsd( new_pose, new_pose);



		Output_data_chi(new_pose, tag+"_"+lead_zero_string_of(chi_angle, 3), rmsd, silent_file_data, chi_angle);
		dump_pdb( new_pose, tag+"_"+lead_zero_string_of(chi_angle, 3));

	}
*/
//	dump_pdb( pose, tag+"_with_hydrogen.pdb" );

/*
  pose::Pose minimize_pose=pose;

  //This makes sure, the O2star hydrogen of the minimized pose are optimizely oriented.

	minimizer.run( minimize_pose, mm, *(scorefxn), options );
	o2star_minimize(minimize_pose, scorefxn);
	minimizer.run( minimize_pose, mm, *(scorefxn), options );
	o2star_minimize(minimize_pose, scorefxn);
	(*scorefxn)(minimize_pose);
	scorefxn->show( std::cout, minimize_pose );
	dump_pdb( minimize_pose, "min_"+tag+".pdb" );

	rmsd=all_atom_rmsd( minimize_pose, pose );


  Output_data(minimize_pose, "min_"+tag, rmsd, silent_file_data);
*/

}

///////////////////////////////////////////////////////////////////////////////////////////////

void
implement_rotamers_test()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	Size rotamer_group_max=72;

	rotamer_ID_struct rotamer_ID={1,1,1,1,1,1,1,1,1,1,1,1,1,1};
	utility::vector1 <Real > current_rotamer;

	std::time_t start_time=time(NULL);

	while( rotamer_ID.group_count<= rotamer_group_max) {
		get_next_rotamer(rotamer_ID, current_rotamer);

		if(option[Verbose]) {
			 std::cout << "group_count= " << rotamer_ID.group_count;
			 std::cout << "sub_group_count= " << rotamer_ID.subgroup_count << std::endl;
		}
	}

	std::time_t end_time=time(NULL);

	std::cout<< "Time taken= " << (long)(end_time - start_time) << "seconds." << std::endl;


}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
main_function_test()
{
	// Initialization
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace core::kinematics;
	using namespace core::optimization;
	using namespace core::io::silent;
	using namespace core::id;
	using namespace protocols::rna;

	std::cout.setf (std::ios_base::left , std::ios_base::adjustfield); //customize std::cout output format
	///////////////// Declare type of residue as RNA type//////////////////////////////////////

  ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "rna" );

	///////////Create RunTimeParameters objects which is will referenced throughout the code/////////////////////////////////////////////////////////////////////////
	Initialize_RunTimeParameters();

	///////////////// Options//////////////////////////////////////////////////////////////////
  std::string rebuild_algorithm_input = option[ rebuild_algorithm];
  std::string silent_file;
	std::string score_function_option = option[score_function];
  Size num_residue_reb = RunTimeParameters::rebuild_residue_list.size();

	std::cout << "score_function_option= " << score_function_option << std::endl;
  std::cout << "Number of residue rebuilt= " << num_residue_reb << std::endl;

  ////////////////Set the parameter for the energy score function and the minimizer//////////////////////

  ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( RNA_HIRES_WTS );
	if ( score_function_option != "default") {
		scorefxn = ScoreFunctionFactory::create_score_function( score_function_option );
	}


	//rna_fa_atr_base 0.23
	//rna_fa_rep_base 0.12
	//rna_base_stack?? This the angle dependence one?

//rna_stack_axis

	std::string bonus_stacking_score_option = option[bonus_stacking_score];
	if(bonus_stacking_score_option!="false"){

		if(bonus_stacking_score_option=="stack1"){
			scorefxn->set_weight( rna_fa_atr_base, 0.23);
			scorefxn->set_weight( rna_fa_rep_base, 0.12);
		}

		if(bonus_stacking_score_option=="stack2"){
			scorefxn->set_weight( rna_fa_atr_base, 0.46);
			scorefxn->set_weight( rna_fa_rep_base, 0.24);
		}

		if(bonus_stacking_score_option=="fa_stack"){
			scorefxn->set_weight( fa_stack, 0.50);
		}
	}

	scorefxn->set_weight( angle_constraint, 1.0 );
	scorefxn->set_weight( atom_pair_constraint, 1.0 );

//rna_hires_base_stacking.wts

//	Real fa_atr_weight=scorefxn->get_weight(fa_atr);
//	Real fa_rep_weight=scorefxn->get_weight(fa_rep);
//	std::cout << "fa_atr_weight= " << fa_atr_weight << std::endl;
//	std::cout << "fa_rep_weight= " << fa_rep_weight << std::endl;
//	std::cout << "option[Atr_rep_reweight_scaling]= " << option[Atr_rep_reweight_scaling] << std::endl;

//	scorefxn->set_weight(fa_atr, fa_atr_weight*option[Atr_rep_reweight_scaling]);
//	scorefxn->set_weight(fa_rep, fa_rep_weight*option[Atr_rep_reweight_scaling]);
//	std::cout << "new fa_atr_weight= " << scorefxn->get_weight(fa_atr) << std::endl;
//	std::cout << "new fa_rep_weight= " << scorefxn->get_weight(fa_rep) << std::endl;




	if(option[minimize_and_score_sugar]==true){
  	scorefxn->set_weight( rna_sugar_close, 0.7 );
  } else {
    scorefxn->set_weight( rna_sugar_close, 0.0 );
 	}

	AtomTreeMinimizer minimizer;
	float const dummy_tol( 0.0000025);
	bool const use_nblist( true );
	MinimizerOptions options( "dfpmin", dummy_tol, use_nblist, false, false );
	options.nblist_auto_update( true ); //update neighor atom position

  //////////////////////////////////////////////////////////////////////////////////////////////////////

  kinematics::MoveMap mm;
	mm.set_bb( true );
	mm.set_chi( true );
  SilentFileData silent_file_data;  //This setup the silent files//

  ///////////////Read in pdb data to create original_pose ////////////////////////////////////////////////////////////

	pose::Pose original_pose;

	std::string infile = option[ in ::file::s ][1];
	std::cout << infile << std::endl;
	core::import_pose::pose_from_pdb( original_pose, *rsd_set, infile );
	protocols::rna::ensure_phosphate_nomenclature_matches_mini_parin( original_pose);
  o2star_minimize(original_pose, scorefxn);	//Minimize o2star hydrogen to native position.

//////////////Set jump_point and cut_point in fold_tree/////////////////////////////////////////////


	Setup_jump_and_cut_point(original_pose);//ACTIVATE THIS
	std::cout << "finish Setup_jump_and_cut_point" << std::endl;


	///////////////////Miminize pose//////////////////////////////////////////////////////////////////////////////////////
	pose::Pose minimized_original_pose=original_pose;

	//Setup
	Setup_chain_break_scoring_for_minimize_original_pose(minimized_original_pose); //ACTIVATE THIS

	protocols::rna::ensure_phosphate_nomenclature_matches_mini_parin( minimized_original_pose);
	if(option[full_base_pose_minimize]==false){
			Freeze_unrebuild_residues(mm, num_residue_reb/*Number of res left unfreezed */, minimized_original_pose);
	}
	//Since the movemap was update, might need to refreeze the sugar_torsions
	if(option[minimize_and_score_sugar]==false){
			Freeze_sugar_torsions(mm, minimized_original_pose);
	} else {
//			UnFreeze_sugar_torsions(mm, minimized_original_pose); Problem since this unfreeze every residue sugar torsion
	}

	if(option[test_mode]==false){
		minimizer.run( minimized_original_pose, mm, *(scorefxn), options);
		o2star_minimize(minimized_original_pose, scorefxn);//This makes sure, the O2star hydrogen of the minimized pose are optimizely oriented.
		minimizer.run( minimized_original_pose, mm, *(scorefxn), options);
		o2star_minimize(minimized_original_pose, scorefxn);
	}


//For full minimize mode////
//	Setup_jump_and_cut_point(original_pose);
//	Setup_jump_and_cut_point(minimized_original_pose);
//	std::cout << "finish Setup_jump_and_cut_point" << std::endl;

	std::cout << "get_O2_star_torsion_list(original_pose)" << std::endl;
	utility::vector1 <Real> original_pose_O2star_torsion_list= get_O2star_torsion_list(original_pose, num_residue_reb);
	std::cout << "get_O2_star_torsion_list(minimized_original_pose)" << std::endl;
	utility::vector1 <Real> minimized_original_pose_O2star_torsion_list= get_O2star_torsion_list(minimized_original_pose, num_residue_reb);

//Initialize original  pose and  minimize original pose/////////////////////////////

	std::vector< pose::Pose> original_poses = create_partially_rebuild_pose_list(original_pose,  "original");
	std::vector< pose::Pose> minimized_then_delete_poses = create_partially_rebuild_pose_list(minimized_original_pose, "fully_minimized_then_delete");

//	for(Size i=1; i<= minimized_then_delete_poses.size(); i++){
//		std::cout << "minimized_then_delelets_poses[" << i << "] constraint_set" << std::endl;
//		ConstraintSetOP cst_set( minimized_then_delete_poses[i].constraint_set()->clone() );
//		assert( cst_set ); //if ( !cst_set ) cst_set = new ConstraintSet();
//		cst_set->show(std::cout);
//	}

	std::cout << "In main_function_test fold_tree of partially_rebuild_original_pose: " << std::endl;
	for(Size current_rebuild_num=0; current_rebuild_num<original_poses.size(); current_rebuild_num++){
			std::cout << "  "; Output_fold_tree_info(original_poses[current_rebuild_num], 	"partial_pose, current_rebuild_num_" + lead_zero_string_of(current_rebuild_num, 2));
	}

	//Obsolete
	std::vector< pose::Pose> minimized_original_poses = original_poses;

////////////////Call one of the rebuilt alogirithm///////////////////////////////////

	if(rebuild_algorithm_input=="Map"){
//		Map_test(backbone_rotamers_groups_north, original_poses[0], original_poses[num_residue_reb]);

	} else if(rebuild_algorithm_input=="iterative_build_one_res"){
	    iterative_build_one_res_test(original_poses, minimized_original_poses, minimized_then_delete_poses, original_pose, scorefxn, minimizer, options, mm,  silent_file_data);
	} else if(rebuild_algorithm_input=="Rescore_pose_list"){
	    Rescore_pose_list_evelope(original_poses, minimized_original_poses, minimized_then_delete_poses, scorefxn, minimizer, options, mm);
	} else {

		std::cout << "Error no rebuilt algorithm selected=" << std::endl;
	}

}
///////////////////////////////////////////////////////////////
void*
my_main( void* )
{

	using namespace basic::options;

  if ( option[ main_function ] ){  //added this option
	  main_function_test();
	} else if ( option[   Full_minimize_function ] ){  //added this option
		 Full_minimize();
	} else if ( option[   Chi_angle_function] ){  //added this option
		 Chi_angle_test();
	} else if ( option[   RNA_torsion_bug_function] ){  //added this option
		 RNA_torsion_bug();
//	} else if ( option[   Minimize_chi_angle_function] ){  //added this option
//		 Minimize_chi_angle_test();
//	} else if ( option[  Rescore_pose_list_function] ){  //added this option
//		 		Rescore_pose_list();
	} else if ( option[  Torsion_analysis_function] ){  //added this option
		 		Torsion_analysis_test();
	} else if ( option[   print_hbonds ] ){  //added this option
		 print_hbonds_test();
	} else if ( option[   cut_bulge ] ){  //added this option
		 cut_bulge_test();
	} else if ( option[   split_pdb_file ] ){  //added this option
		 split_pdb_file_test();
	} else if ( option[   convert_star_to_dash] ){  //added this option
		 convert_star_to_dash_test();
	} else if ( option[   final_output] ){  //added this option
		 final_output_test();
	} else if ( option[implement_rotamers] ){
		implement_rotamers_test();
	} else if ( option[ central_computer ] ){  //added this option
	  central_evaluation();
	} else if ( option[ color_by_geom_sol ] ){
	  color_by_geom_sol_RNA_test();
	}
	exit( 0 );
}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	using namespace basic::options;

	//Uh, options?

	NEW_OPT( RNA_torsion_bug_function, "RNA_torsion_bug_function", false);
	NEW_OPT( input_pose, "input_pose", "cluster");
	NEW_OPT( old_file_format, "old_file_format", false);
	NEW_OPT( exclude_bulge_residue_rmsd, "exclude_bulge_residue_rmsd", false);
	NEW_OPT( pucker, "For chi_angle_test function", "");
	NEW_OPT( base, "	For chi_angle_test function", "");
	NEW_OPT( 	vary_geometry , "		vary_geometry when minimize", false);
	NEW_OPT( 	outfile_name , "		outfile_name", "test");
	NEW_OPT( native_sequence_mode, "	Native_sequence_mode...disallow rebuild of non-native sequence", false);
	NEW_OPT( final_output, "	final_output()-RNA journal club", false);
	NEW_OPT( convert_star_to_dash, "	convert_star_to_dash_test()-RNA journal club", false);
	NEW_OPT( split_pdb_file, "	split_pdb_file()-RNA journal club", false);
	NEW_OPT( cut_bulge, "	cut_bulge-RNA journal club", false);
	NEW_OPT( print_hbonds, "print_hbonds_function", false);
	NEW_OPT( more_rotamers, "more rotamers",false);
	NEW_OPT( quick_test, "more rotamers",false);
	NEW_OPT( suite_rmsd_cutoff, "rmsd_cutoff", 2.0 );
	NEW_OPT( cluster_rmsd, "cluster_rmsd", 0.5 );
	NEW_OPT( fa_rep_cutoff, "fa_rep_cutoff", 4.0 );
	NEW_OPT( bonus_stacking_score, "bonus_stacking_score", "false" );
	NEW_OPT( main_function , "Rebuilt nucleotides",false);
	NEW_OPT( Full_minimize_function , "full minimize", false);
	NEW_OPT( 	Chi_angle_function , "Chi_angle_function", false);
//	NEW_OPT( Rescore_pose_list_function, " Rescore_pose_list_function" ,false);
	NEW_OPT( Torsion_analysis_function , "	Torsion_analysis_test", false);
	NEW_OPT( implement_rotamers , "quick function to determine function speed", false);
	NEW_OPT( central_computer  , "Indicate that this node is the central computer",false);
	NEW_OPT( read_job_file , "Read in job file which specify which rotamers to rebuild",false);
	NEW_OPT( HO2star_sampling, "option include O2star_trail, Virtual_HO2star and Native", "Virtual_HO2star");
	NEW_OPT( print_score_only , "Print_torision_angles as part of output if false", true);
	NEW_OPT( starting_pose, "Choose either original or minimize", "minimize");
	NEW_OPT( twoside_bound_original_poses , "rebuilt 1 res with residues both before and after it already built", false);
	NEW_OPT( build_pose , "Build a pose with a given rotamers combination input",false);
	NEW_OPT( full_base_pose_minimize, "Full minimize of the base pose", false);
	NEW_OPT( full_rebuild_pose_minimize, "Full minimize during the rebuild step", false);
	NEW_OPT( test_mode, "Make code faster during testing phrase.", false);
	NEW_OPT( Finer_sampling, "Implement Finer_sampling after initial sampling", false);
	NEW_OPT( all_atom_cluster, "Use all_rebuilt atom rmsd as cluster citeria instead of suite rmsd", false);
	NEW_OPT( minimize_and_score_sugar, "If false, sugar is fixed and its geometry is not score", false);
	NEW_OPT( build_one_res_mode, "This mode rebuilt 1 residue with the less of the structure given", 0);
	NEW_OPT( Verbose, "Verbose", false);
	NEW_OPT( build_loop, "build loop residue instead of end", true);
	NEW_OPT( score_function, "Score function choice, default if rna_hires", "default");
	NEW_OPT( rebuild_sequence, "Specify the rebuild nucleotides in the order you want to rebuild the", "");
	//Example: G5-C6-A7-A8
	NEW_OPT( rotamer_set, "Specify which rotamer set to use", "full");
	NEW_OPT( rebuild_algorithm, "Specify which rebuild_alogorith to use", "");
	NEW_OPT( rebuilt_pose_name, "rebuilt pose name", "01_01");
  NEW_OPT( apply_bound_function, "Apply bound function", " ");
  NEW_OPT( select_base_pose, "select a non-default pose as the base pose for the rebuilt process", "default");
  NEW_OPT( native_screen, "This is a cheat mode to screen only near native pose", "default");
  NEW_OPT( torsion_bin, "Bin size in normal mode , usually 20", 20.0);
  NEW_OPT( fine_torsion_bin, "Finer_sampling bin size", 20.0);
  NEW_OPT( torsion_cutoff, "use with the cheat rotamer option, torsion diff cutoff", 80.0);
	NEW_OPT( bound_cutoff_value, "value used to determine bound condition", 0);
	NEW_OPT( apply_clustering, "Apply group clustering", true);
	NEW_OPT( central_clustering, "Apply clustering in central_evaluation function", true);
	NEW_OPT( pose_kept, "number of pose kept for propagation in branch and bound", 0);
	NEW_OPT( num_branch_kept, "Number of branch to be expanded after rebuilding each residue", 1);
	NEW_OPT( hack_first_rebuild_num, "First_rebuild_num for for-loop in central evaluation and iterative_build_one_res_test..should be 1..aside from when debugging ", 1);
	NEW_OPT( num_residue_rebuilt, "number of residue to rebuilt", 1);
	NEW_OPT( total_nodes_input, "total number of computer nodes", 1);
	NEW_OPT( node_number_input, "the number of the computer node", 1);
	NEW_OPT( score_choice_option, "the number of the computer node", -1);
	NEW_OPT( color_by_geom_sol, "Color by geometric solvation",false );
	NEW_OPT(graphics, "Turn on molecular viewer", false);


	////////////////////////////////////////////////////////////////////////////
	// setup
	////////////////////////////////////////////////////////////////////////////
	devel::init(argc, argv);


	////////////////////////////////////////////////////////////////////////////
	// end of setup
	////////////////////////////////////////////////////////////////////////////

	protocols::viewer::viewer_main( my_main );

}


		////O2star trail for only rebuild nucleotide mode
//			pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( current_pose ));
//			task->initialize_from_command_line();
//			task->nonconst_residue_task(reb_res).and_extrachi_cutoff( 0 );
//			task->nonconst_residue_task(reb_res).or_include_current( true );
		///////////////////////////////////////////////////////////////////////////////////////////////////
/*

// 		protocols::rna::ensure_phosphate_nomenclature_matches_mini_parin( current_pose, reb_res );
//		protocols::rna::ensure_phosphate_nomenclature_matches_mini_parin( current_pose, reb_res+1 );

*/


