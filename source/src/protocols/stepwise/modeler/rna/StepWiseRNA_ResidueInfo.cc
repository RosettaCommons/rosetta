// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file StepWise_ResidueInfo.cc
/// @brief Function used by the residue_info data structure
/// @details
/// @author Parin Sripakdeevong


//////////////////////////////////
#include <protocols/stepwise/modeler/rna/StepWiseRNA_ResidueInfo.hh>
//////////////////////////////////
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/Pose.hh>
#include <string>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>

#include <basic/Tracer.hh>

#include <fstream>
#include <utility/vector1.hh>
#include <map>

#include <core/id/TorsionID.hh> // AUTO IWYU For TorsionID


static basic::Tracer TR( "protocols.stepwise.modeler.rna.StepWiseRNA_ResidueInfo" );

using namespace core;
using namespace ObjexxFCL;

namespace protocols {
namespace stepwise {
namespace modeler {
namespace rna {

void
print_torsion_info( core::pose::Pose const & pose, core::Size const seq_num, core::Size const rna_torsion_number, std::string const & type ) {

	using namespace core::id;
	using namespace core::chemical;
	using namespace core::conformation;

	id::TorsionID torsion_id;
	Real torsion_angle;
	if ( type == "back_bone" ) {
		torsion_id = id::TorsionID( seq_num, id::BB, rna_torsion_number );
		torsion_angle = pose.residue( seq_num ).mainchain_torsion( rna_torsion_number );
	} else if ( type == "side_chain" ) {
		torsion_id = id::TorsionID( seq_num, id::CHI, rna_torsion_number );
		torsion_angle = pose.residue( seq_num ).chi( rna_torsion_number );
	} else {
		TR << "In print_torsion_info, invalid torsion_type: " << type << std::endl;
		exit( 1 );
	}

	TR << "Torsion angle = " << torsion_angle << std::endl;

	id::AtomID id1, id2, id3, id4;
	bool fail = pose.conformation().get_torsion_angle_atom_ids( torsion_id, id1, id2, id3, id4 );
	if ( fail ) {
		TR << "In print_torsion_info, would not get torsion angle atom ids: " <<  std::endl;
		exit( 1 );
	}

	conformation::Residue const & rsd_1 = pose.residue( id1.rsd() );
	conformation::Residue const & rsd_2 = pose.residue( id2.rsd() );
	conformation::Residue const & rsd_3 = pose.residue( id3.rsd() );
	conformation::Residue const & rsd_4 = pose.residue( id4.rsd() );

	TR << "atom " << id1.atomno()  << " " <<  "name = " << rsd_1.type().atom_name( id1.atomno() ) << " type = " << rsd_1.atom_type( id1.atomno() ).name()  << " " << rsd_1.atom_type_index( id1.atomno() )    << " " << rsd_1.atomic_charge( id1.atomno() ) << std::endl;
	TR << "atom " << id2.atomno()  << " " <<  "name = " << rsd_2.type().atom_name( id2.atomno() ) << " type = " << rsd_2.atom_type( id2.atomno() ).name()  << " " << rsd_2.atom_type_index( id2.atomno() )    << " " << rsd_2.atomic_charge( id2.atomno() ) << std::endl;
	TR << "atom " << id3.atomno()  << " " <<  "name = " << rsd_3.type().atom_name( id3.atomno() ) << " type = " << rsd_3.atom_type( id3.atomno() ).name()  << " " << rsd_3.atom_type_index( id3.atomno() )    << " " << rsd_3.atomic_charge( id3.atomno() ) << std::endl;
	TR << "atom " << id4.atomno()  << " " <<  "name = " << rsd_4.type().atom_name( id4.atomno() ) << " type = " << rsd_4.atom_type( id4.atomno() ).name()  << " " << rsd_4.atom_type_index( id4.atomno() )    << " " << rsd_4.atomic_charge( id4.atomno() ) << std::endl;
}

void
output_res_map( std::map< core::Size, core::Size > const & my_map, core::Size const max_seq_num ){

	using namespace ObjexxFCL;
	using namespace ObjexxFCL::format;
	using namespace std;
	core::Size spacing = 4;
	TR << std::setw( 30 ) << "full_pose_seq_num:";
	for ( core::Size seq_num = 1; seq_num <= max_seq_num; seq_num++ ) {
		TR << std::setw( spacing ) << seq_num;
	}
	TR << std::endl;

	TR << std::setw( 30 ) << "partial_pose_seq_num:";
	for ( core::Size seq_num = 1; seq_num <= max_seq_num; seq_num++ ) {
		if ( my_map.find( seq_num ) != my_map.end() ) {
			TR << std::setw( spacing ) << my_map.find( seq_num )->second;
		} else {
			TR << A( spacing, "-" );
		}

	}
	TR << std::endl;
}

utility::vector1< Residue_info >
get_copy_dofs_from_fasta( std::string const & full_fasta_sequence ){

	utility::vector1< Residue_info > full_copy_dofs;

	for ( core::Size n = 0; n <= full_fasta_sequence.size() - 1; n++ ) {
		std::string one_letter_name = string_of( full_fasta_sequence[n] );
		Residue_info residue;
		residue.seq_num = n + 1;
		if ( one_letter_name == "A" || one_letter_name == "a" ) residue.name = "RAD";
		if ( one_letter_name == "C" || one_letter_name == "c" ) residue.name = "RCY";
		if ( one_letter_name == "U" || one_letter_name == "u" ) residue.name = "URA";
		if ( one_letter_name == "G" || one_letter_name == "g" ) residue.name = "RGU";

		full_copy_dofs.push_back( residue );
	}

	return full_copy_dofs;
}


bool
contain_residue_at_seq_num( core::Size seq_num, utility::vector1 < Residue_info > const & copy_dofs ){

	for ( auto const & copy_dof : copy_dofs ) {
		if ( seq_num == copy_dof.seq_num ) {
			return true;
		}
	}
	return false;
}

utility::vector1< utility::vector1< Residue_info > >
create_strand_list( utility::vector1< Residue_info > const & copy_dofs ){

	utility::vector1< utility::vector1< Residue_info > > residue_group_list;

	utility::vector1< Residue_info > Sorted_copy_dofs = copy_dofs;
	//Sort by seq_number, lowest sequence number at the top of the vector list.
	sort_copy_dofs( Sorted_copy_dofs );

	core::Size j = 1;
	while ( j <= Sorted_copy_dofs.size() ) {

		core::Size first_element = j;

		//Test if Sorted_copy_dofs[j] contain an adjacent residue at the three_prime_end
		while ( ( j < Sorted_copy_dofs.size() ) && ( ( Sorted_copy_dofs[j].seq_num + 1 ) == Sorted_copy_dofs[j + 1].seq_num ) ) {
			j++;
		}
		core::Size last_element = j;

		utility::vector1 < Residue_info > residue_group;

		for ( core::Size element = first_element; element <= last_element; element++ ) {
			residue_group.push_back( Sorted_copy_dofs[element] );
		}
		residue_group_list.push_back( residue_group );

		j++;
	}

	return residue_group_list;
}


utility::vector1 < Residue_info >
set_difference( utility::vector1 < Residue_info > const & copy_dofs_1, utility::vector1 < Residue_info > const & copy_dofs_2 ){

	utility::vector1 < Residue_info > set_difference_copy_dofs;

	for ( auto const & copy_dof : copy_dofs_1 ) {
		if ( contain_residue_at_seq_num( copy_dof.seq_num, copy_dofs_2 ) ) continue;
		set_difference_copy_dofs.push_back( copy_dof );
	}

	return set_difference_copy_dofs;
}

utility::vector1 < Residue_info >
set_union( utility::vector1 < Residue_info > const & copy_dofs_1, utility::vector1 < Residue_info > const & copy_dofs_2 ){

	utility::vector1 < Residue_info > union_copy_dofs = copy_dofs_1;

	for ( auto const & copy_dof : copy_dofs_2 ) {
		if ( contain_residue_at_seq_num( copy_dof.seq_num, union_copy_dofs ) == false ) {
			union_copy_dofs.push_back( copy_dof );
		}
	}

	return union_copy_dofs;
}

bool
copy_dofs_sort_criterion( Residue_info residue_info_1, Residue_info residue_info_2 ){
	//Sort by seq_number, lowest sequence number at the top of the vector list.
	return ( residue_info_1.seq_num < residue_info_2.seq_num );
}

void
sort_copy_dofs( utility::vector1< Residue_info > & copy_dofs ) {
	//Need to check if this work with vector1, if not switch to std::vector
	sort( copy_dofs.begin(), copy_dofs.end(), copy_dofs_sort_criterion );
}

} //rna
} //modeler
} //stepwise
} //protocols
