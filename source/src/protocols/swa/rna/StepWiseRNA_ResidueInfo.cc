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
#include <core/pose/Pose.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/chemical/AtomType.hh>
#include <string>
//#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>
//#include <core/conformation/ResidueMatcher.hh>
#include <core/conformation/Conformation.hh>

#include <basic/Tracer.hh>

#include <iostream>
#include <fstream>
#include <sstream>
#include <utility/vector1.hh>
#include <map>


static basic::Tracer TR( "protocols.swa.rna.StepWiseRNA_ResidueInfo" );

using namespace core;
using namespace ObjexxFCL;

namespace protocols {
namespace swa {
namespace rna {

	void
	Print_torsion_info( core::pose::Pose const & pose, core::Size const seq_num, core::Size const rna_torsion_number, std::string const type ){

		using namespace core::id;
		using namespace core::chemical;
		using namespace core::conformation;
//		using namespace protocols::rna;

		id::TorsionID torsion_id;
		Real torsion_angle;
		if ( type == "back_bone" ){
			torsion_id = id::TorsionID( seq_num, id::BB, rna_torsion_number );
			torsion_angle = pose.residue( seq_num ).mainchain_torsion( rna_torsion_number );
		} else if ( type == "side_chain" ){
			torsion_id = id::TorsionID( seq_num, id::CHI, rna_torsion_number );
			torsion_angle = pose.residue( seq_num ).chi( rna_torsion_number );
		} else{
			TR << "In Print_torsion_info, invalid torsion_type: " << type << std::endl;
			exit( 1 );
		}

		TR << "Torsion angle = " << torsion_angle << std::endl;

		id::AtomID id1, id2, id3, id4;
		bool fail = pose.conformation().get_torsion_angle_atom_ids( torsion_id, id1, id2, id3, id4 );
		if ( fail ){
			TR << "In Print_torsion_info, would not get torsion angle atom ids: " <<  std::endl;
			exit( 1 );
		}

		conformation::Residue const & rsd_1 = pose.residue( id1.rsd() );
		conformation::Residue const & rsd_2 = pose.residue( id2.rsd() );
		conformation::Residue const & rsd_3 = pose.residue( id3.rsd() );
		conformation::Residue const & rsd_4 = pose.residue( id4.rsd() );

		TR << "atom " << id1.atomno()  << " " << 	"name = " << rsd_1.type().atom_name( id1.atomno() ) << " type = " << rsd_1.atom_type( id1.atomno() ).name()  << " " << rsd_1.atom_type_index( id1.atomno() ) 			<< " " << rsd_1.atomic_charge( id1.atomno() ) << std::endl;
		TR << "atom " << id2.atomno()  << " " << 	"name = " << rsd_2.type().atom_name( id2.atomno() ) << " type = " << rsd_2.atom_type( id2.atomno() ).name()  << " " << rsd_2.atom_type_index( id2.atomno() ) 			<< " " << rsd_2.atomic_charge( id2.atomno() ) << std::endl;
		TR << "atom " << id3.atomno()  << " " << 	"name = " << rsd_3.type().atom_name( id3.atomno() ) << " type = " << rsd_3.atom_type( id3.atomno() ).name()  << " " << rsd_3.atom_type_index( id3.atomno() ) 			<< " " << rsd_3.atomic_charge( id3.atomno() ) << std::endl;
		TR << "atom " << id4.atomno()  << " " << 	"name = " << rsd_4.type().atom_name( id4.atomno() ) << " type = " << rsd_4.atom_type( id4.atomno() ).name()  << " " << rsd_4.atom_type_index( id4.atomno() ) 			<< " " << rsd_4.atomic_charge( id4.atomno() ) << std::endl;

	}


	utility::vector1 < Residue_info >
	Convert_rebuild_residue_string_to_list( std::string const& rebuild_residue_string ){

		utility::vector1 < Residue_info > rebuild_residue_list;

//	if(rebuild_residue_string=="none"){
//		return rebuild_residue_list;
//	}

		 utility::vector1< std::string > nucleotides_token = Tokenize( rebuild_residue_string, "-" );

			for ( Size i = 1; i <= nucleotides_token.size(); i++ ){
//			TR<< nucleotides_token[i] << " ";

				Residue_info res_info;
				res_info.name = Get_three_letter_name( nucleotides_token[i].substr( 0, 1 ) );
				//std::string seq_num_string=nucleotides_token[i].substr(1,2);
				std::string seq_num_string = nucleotides_token[i].substr( 1, nucleotides_token[i].length() - 1 );
				res_info.seq_num = string_to_int( seq_num_string );
				rebuild_residue_list.push_back( res_info );

			}
//		TR << std::endl;

		return rebuild_residue_list;
	}


	void
	Output_residue_struct( Residue_info const & residue ){

		using namespace ObjexxFCL;
		using namespace ObjexxFCL::fmt;

		TR << Get_one_letter_name( residue.name );
		TR << lead_zero_string_of( residue.seq_num, 2 );
		TR << A( 1, " " );
	}


	std::string
	Get_one_letter_name( std::string const & three_letter_name ){
		if ( three_letter_name == "RAD" ) return "A";
		if ( three_letter_name == "RCY" ) return "C";
		if ( three_letter_name == "URA" ) return "U";
		if ( three_letter_name == "RGU" ) return "G";
		TR << "In get_one_letter_name_function, an invalid three_letter_name was passed into the function: " << three_letter_name << std::endl;
		exit ( 1 );
	}

	std::string
	Get_three_letter_name( std::string const & one_letter_name ){
		if ( one_letter_name == "A" ) return "RAD";
		if ( one_letter_name == "C" ) return "RCY";
		if ( one_letter_name == "U" ) return "URA";
		if ( one_letter_name == "G" ) return "RGU";
		TR << "In get_three_letter_name_function, an invalid one_letter_name was passed into the function: " << one_letter_name << std::endl;
		exit ( 1 );
	}




	Size
	get_max_seq_num_from_res_map( std::map< core::Size, core::Size > const & my_map ){

		Size max_seq_num = 0;
	  for ( std::map< Size, Size > ::const_iterator it = my_map.begin(); it != my_map.end(); it++ ){
	    TR << it->first << " =  > " << it->second << std::endl;
	    if ( it->first >= max_seq_num ) max_seq_num = it->first;
		}
		return max_seq_num;
	}

	void
	output_res_map( std::map< core::Size, core::Size > const & my_map, Size const max_seq_num ){

		using namespace ObjexxFCL;
		using namespace ObjexxFCL::fmt;
		using namespace std;
		Size spacing = 4;
		TR << std::setw( 30 ) << "full_pose_seq_num:";
		for ( Size seq_num = 1; seq_num <= max_seq_num; seq_num++ ){
			TR << std::setw( spacing ) << seq_num;
		}
		TR << std::endl;

		TR << std::setw( 30 ) << "partial_pose_seq_num:";
		for ( Size seq_num = 1; seq_num <= max_seq_num; seq_num++ ){
			if ( my_map.find( seq_num ) != my_map.end() ){
				TR << std::setw( spacing ) << my_map.find( seq_num )->second;
			} else{
				TR << A( spacing, "-" );
			}

		}
		TR << std::endl;

	}


	void
	Output_residue_list( utility::vector1< Residue_info > residue_list ){
		using namespace ObjexxFCL;
		using namespace ObjexxFCL::fmt;


		sort_residue_list( residue_list ); //maybe sure the list is sorted

		Size seq_num = 1;
		for ( Size n = 1; n <= residue_list.size(); n++ ){
			Residue_info residue = residue_list[n];

			while ( seq_num < residue.seq_num ){
				TR << A( 4, " " );
				seq_num++;
			}

			Output_residue_struct( residue );
			seq_num++;
		}

		TR << std::endl;

	}

	utility::vector1< Residue_info >
	Get_residue_list_from_fasta( std::string const full_fasta_sequence ){

		utility::vector1< Residue_info > full_residue_list;

		for ( Size n = 0; n <= full_fasta_sequence.size() - 1; n++ ){
			std::string one_letter_name = string_of( full_fasta_sequence[n] );
			Residue_info residue;
			residue.seq_num = n + 1;
			if ( one_letter_name == "A" || one_letter_name == "a" ) residue.name = "RAD";
			if ( one_letter_name == "C" || one_letter_name == "c" ) residue.name = "RCY";
			if ( one_letter_name == "U" || one_letter_name == "u" ) residue.name = "URA";
			if ( one_letter_name == "G" || one_letter_name == "g" ) residue.name = "RGU";

			full_residue_list.push_back( residue );
		}

		return full_residue_list;
	}


	Residue_info
	Get_residue_from_seq_num( Size const & seq_num, utility::vector1 < Residue_info > const & residue_list ){

		for ( Size i = 1; i <= residue_list.size(); i++ ){
			if ( seq_num == residue_list[i].seq_num ){
				return residue_list[i];
			}
		}
		TR << "Error, in Get_residue_from_seq_num function. The seq_num " << seq_num << " does not exist in the residue_list" << 	std::endl;
		exit ( 1 );
	}


	bool
	Contain_residue_at_seq_num( Size seq_num, utility::vector1 < Residue_info > const & residue_list ){

		for ( Size j = 1; j <= residue_list.size(); j++ ){
			if ( seq_num == residue_list[j].seq_num ) {
				return true;
			}
		}
		return false;
	}

	utility::vector1 < utility::vector1 < Residue_info > >
	Create_strand_list( utility::vector1 < Residue_info > const & residue_list ){

		utility::vector1 < utility::vector1 < Residue_info > > residue_group_list;

		utility::vector1 < Residue_info > Sorted_residue_list = residue_list;
		//Sort by seq_number, lowest sequence number at the top of the vector list.
		sort_residue_list( Sorted_residue_list );

		Size j = 1;
		while ( j <= Sorted_residue_list.size() ){

			Size first_element = j;

			//Test if Sorted_residue_list[j] contain an adjacent residue at the three_prime_end
			while ( ( j < Sorted_residue_list.size() ) && ( ( Sorted_residue_list[j].seq_num + 1 ) == Sorted_residue_list[j + 1].seq_num ) ){
				j++;
			}
			Size last_element = j;

			utility::vector1 < Residue_info > residue_group;

			for ( Size element = first_element; element <= last_element; element++ ){
				residue_group.push_back( Sorted_residue_list[element] );
			}
			residue_group_list.push_back( residue_group );

			j++;
		}

		return residue_group_list;
	}


	utility::vector1 < Residue_info >
	Set_Difference( utility::vector1 < Residue_info > const & residue_list_1, utility::vector1 < Residue_info > const & residue_list_2 ){

		utility::vector1 < Residue_info >	set_difference_residue_list;

		for ( Size i = 1; i <= residue_list_1.size(); i++ ){
			if ( Contain_residue_at_seq_num( residue_list_1[i].seq_num, residue_list_2 ) ) continue;
			set_difference_residue_list.push_back( residue_list_1[i] );
		}

		return set_difference_residue_list;

	}

	utility::vector1 < Residue_info >
	Set_Union( utility::vector1 < Residue_info > const & residue_list_1, utility::vector1 < Residue_info > const & residue_list_2 ){

		utility::vector1 < Residue_info > union_residue_list = residue_list_1;

		for ( Size i = 1; i <= residue_list_2.size(); i++ ){
			if ( Contain_residue_at_seq_num( residue_list_2[i].seq_num, union_residue_list ) == false ){
				union_residue_list.push_back( residue_list_2[i] );
			}
		}

		return union_residue_list;

	}

	bool
	residue_list_sort_citeria( Residue_info residue_info_1, Residue_info residue_info_2 ){
		//Sort by seq_number, lowest sequence number at the top of the vector list.
		return ( residue_info_1.seq_num < residue_info_2.seq_num );
	}

	void
	sort_residue_list( utility::vector1< Residue_info > & residue_list ) {
		//Need to check if this work with vector1, if not switch to std::vector
		sort( residue_list.begin(), residue_list.end(), residue_list_sort_citeria );
	}

}
}
}
