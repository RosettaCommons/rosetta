// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWise_ResidueInfo.cc
/// @brief Function used by the residue_info data structure
/// @detailed
/// @author Parin Sripakdeevong


//////////////////////////////////
#include <protocols/swa/rna/StepWiseRNA_Util.hh>
#include <protocols/swa/rna/StepWiseRNA_ResidueInfo.hh>
//////////////////////////////////
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <string>

#include <iostream>
// AUTO-REMOVED #include <fstream>
#include <sstream>
#include <utility/vector1.hh>
#include <map>


using namespace core;
using namespace ObjexxFCL;

namespace protocols {
namespace swa {
namespace rna {

	void
	Output_residue_struct(Residue_info const & residue){

		using namespace ObjexxFCL;
		using namespace ObjexxFCL::fmt;

		std::cout << Get_one_letter_name(residue.name);
		std::cout << lead_zero_string_of(residue.seq_num, 2);
		std::cout << A(1," ");
	}


	std::string
	Get_one_letter_name(std::string const & three_letter_name){
		if(three_letter_name=="RAD") return "A";
		if(three_letter_name=="RCY") return "C";
		if(three_letter_name=="URA") return "U";
		if(three_letter_name=="RGU") return "G";
		std::cout << "In get_one_letter_name_function, an invalid three_letter_name was passed into the function: " << three_letter_name << std::endl;
		exit (1);
	}

	std::string
	Get_three_letter_name(std::string const & one_letter_name){
		if(one_letter_name=="A") return "RAD";
		if(one_letter_name=="C") return "RCY";
		if(one_letter_name=="U") return "URA";
		if(one_letter_name=="G") return "RGU";
		std::cout << "In get_three_letter_name_function, an invalid one_letter_name was passed into the function: " << one_letter_name <<std::endl;
		exit (1);
	}




	Size
	get_max_seq_num_from_res_map(std::map< core::Size, core::Size > const & my_map){

		Size max_seq_num=0;
	  for (std::map< Size, Size>::const_iterator it=my_map.begin(); it!=my_map.end(); it++ ){
	    std::cout << it->first << " => " << it->second << std::endl;
	    if(it->first>=max_seq_num) max_seq_num=it->first;
		}
		return max_seq_num;
	}

	void
	output_res_map(std::map< core::Size, core::Size > const & my_map, Size const max_seq_num){

		using namespace ObjexxFCL;
		using namespace ObjexxFCL::fmt;
		using namespace std;
		Size spacing=4;
		std::cout << std::setw(30) << "full_pose_seq_num:";
		for(Size seq_num=1; seq_num<=max_seq_num; seq_num++){
			std::cout << std::setw(spacing) << seq_num;
		}
		std::cout << std::endl;

		std::cout << std::setw(30) << "partial_pose_seq_num:";
		for(Size seq_num=1; seq_num<=max_seq_num; seq_num++){
			if(my_map.find(seq_num)!=my_map.end()){
				std::cout << std::setw(spacing) << my_map.find(seq_num)->second;
			}else{
				std::cout << A(spacing,"-");
			}

		}
		std::cout << std::endl;

	}

	void
	output_is_prepend_map(std::map< core::Size, bool > const & my_map, Size const max_seq_num){

		using namespace ObjexxFCL;
		using namespace ObjexxFCL::fmt;
/*
		Size max_seq_num=0;
	  for (std::map<Size,bool>::const_iterator it=my_map.begin() ; it != my_map.end(); it++ ){
	    std::cout << it->first << " => " << it->second << std::endl;
	    if(it->first >= max_seq_num) max_seq_num=it->first;
		}
*/
		Size spacing=4;
		std::cout << std::setw(30) << "Is_residue_prepend:";
		for(Size seq_num=1; seq_num<=max_seq_num; seq_num++){
			char prepend_char;
			if(my_map.find(seq_num)!=my_map.end()){
				prepend_char = (my_map.find(seq_num)->second) ? 'P' : 'A';
			}else{
				prepend_char = '-';
			}
			std::cout << std::setw(spacing) << prepend_char;
		}
		std::cout << std::endl;

	}


	void
	Output_residue_list(utility::vector1<Residue_info> residue_list){
		using namespace ObjexxFCL;
		using namespace ObjexxFCL::fmt;


		sort_residue_list(residue_list); //maybe sure the list is sorted

		Size seq_num=1;
		for(Size n=1; n<=residue_list.size(); n++){
			Residue_info residue=residue_list[n];

			while(seq_num<residue.seq_num){
				std::cout << A(4," ");
				seq_num++;
			}

			Output_residue_struct(residue);
			seq_num++;
		}

		std::cout << std::endl;

	}

	utility::vector1< Residue_info >
	Get_residue_list_from_fasta(std::string const full_fasta_sequence){

		utility::vector1< Residue_info > full_residue_list;

		for(Size n=0; n<=full_fasta_sequence.size()-1; n++){
			std::string one_letter_name=string_of(full_fasta_sequence[n]);
			Residue_info residue;
			residue.seq_num=n+1;
			if(one_letter_name=="A" || one_letter_name=="a") residue.name="RAD";
			if(one_letter_name=="C" || one_letter_name=="c") residue.name="RCY";
			if(one_letter_name=="U" || one_letter_name=="u") residue.name="URA";
			if(one_letter_name=="G" || one_letter_name=="g") residue.name="RGU";

			full_residue_list.push_back(residue);
		}

		return full_residue_list;
	}


	Residue_info
	Get_residue_from_seq_num(Size const & seq_num, utility::vector1 <Residue_info> const & residue_list){

		for(Size i=1; i<=residue_list.size(); i++){
			if(seq_num==residue_list[i].seq_num){
				return residue_list[i];
			}
		}
		std::cout << "Error, in Get_residue_from_seq_num function. The seq_num " << seq_num << " does not exist in the residue_list" << 	std::endl;
		exit (1);
	}


	bool
	Contain_residue_at_seq_num(Size seq_num, utility::vector1 <Residue_info> const & residue_list){

		for(Size j=1; j<=residue_list.size(); j++){
			if(seq_num==residue_list[j].seq_num) {
				return true;
			}
		}
		return false;
	}

	utility::vector1 < utility::vector1 <Residue_info> >
	Create_strand_list(utility::vector1 <Residue_info> const & residue_list){

		utility::vector1 < utility::vector1 <Residue_info> >residue_group_list;

		utility::vector1 <Residue_info> Sorted_residue_list = residue_list;
		//Sort by seq_number, lowest sequence number at the top of the vector list.
		sort_residue_list(Sorted_residue_list);

		Size j=1;
		while(j<=Sorted_residue_list.size()){

			Size first_element=j;

			//Test if Sorted_residue_list[j] contain an adjacent residue at the three_prime_end
			while((j<Sorted_residue_list.size()) && ((Sorted_residue_list[j].seq_num+1)==Sorted_residue_list[j+1].seq_num)){
				j++;
			}
			Size last_element=j;

			utility::vector1 <Residue_info> residue_group;

			for(Size element=first_element; element<=last_element; element++){
				residue_group.push_back(Sorted_residue_list[element]);
			}
			residue_group_list.push_back(residue_group);

			j++;
		}

		return residue_group_list;
	}


	utility::vector1 <Residue_info>
	Set_Difference(utility::vector1 <Residue_info> const & residue_list_1, utility::vector1 <Residue_info> const & residue_list_2){

		utility::vector1 <Residue_info>	set_difference_residue_list;

		for(Size i=1; i<= residue_list_1.size(); i++){
			if(Contain_residue_at_seq_num(residue_list_1[i].seq_num, residue_list_2)) continue;
			set_difference_residue_list.push_back(residue_list_1[i]);
		}

		return set_difference_residue_list;

	}

	utility::vector1 <Residue_info>
	Set_Union(utility::vector1 <Residue_info> const & residue_list_1, utility::vector1 <Residue_info> const & residue_list_2){

		utility::vector1 <Residue_info> union_residue_list= residue_list_1;

		for(Size i=1; i<=residue_list_2.size(); i++){
			if(Contain_residue_at_seq_num(residue_list_2[i].seq_num, union_residue_list)==false){
				union_residue_list.push_back(residue_list_2[i]);
			}
		}

		return union_residue_list;

	}

	bool
	residue_list_sort_citeria(Residue_info residue_info_1, Residue_info residue_info_2){
		//Sort by seq_number, lowest sequence number at the top of the vector list.
		return (residue_info_1.seq_num < residue_info_2.seq_num);
	}

	void
	sort_residue_list(utility::vector1<Residue_info>& residue_list) {
		//Need to check if this work with vector1, if not switch to std::vector
		sort(residue_list.begin(), residue_list.end(), residue_list_sort_citeria);
	}

}
}
}
