// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/sewing/extra_functions.cc
/// @brief  Utility functions for SEWING tests
/// @author Sharon Guffy (guffy@email.unc.edu)


#ifndef INCLUDED_protocols_sewing_extra_functions_HH
#define INCLUDED_protocols_sewing_extra_functions_HH

// Project Headers
#include <protocols/sewing/data_storage/SmartAssembly.hh>
#include <protocols/sewing/data_storage/SmartSegment.hh>
#include <protocols/sewing/data_storage/LigandSegment.hh>
#include <protocols/sewing/data_storage/LigandResidue.hh>
#include <protocols/sewing/hashing/ModelFileReader.hh>
#include <protocols/sewing/hashing/AlignmentFileGeneratorMover.hh>
#include <protocols/sewing/hashing/hasher_data.hh>


// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Protocol Headers
#include <basic/Tracer.hh>


namespace sewing_testing {
static basic::Tracer TR("sewing_testing");
///@brief Utility function to set up for testing deletes and switches
protocols::sewing::data_storage::SmartAssemblyOP create_simple_assembly( protocols::sewing::data_storage::SmartAssemblyOP assembly ){
	//Segments 4-6 will always be considered vital
	assembly->set_starting_segment( assembly->get_segment_vector()->at( 4 ), "all" );
	//Length of segment 4: 9
	//Length of segment 6: 13
	//Length of segment 13: 14
	//Length of segment 21: 7
	//Start with node 4-6
	//Add node 13-15 c-term
	bool n_term = false;
	core::Size first_res_id = 6;
	core::Size second_seg_id = 13; //All residue numbers will now be decreased by 1
	core::Size second_res_id = 7;
	assembly->add_segment( n_term, second_seg_id, first_res_id , second_res_id );
	//Add node 19-21 n-term
	n_term = true;
	//first_res_id is in segment 4
	//second_res_id is in segment 21
	first_res_id = 5;
	second_res_id = 3;
	second_seg_id = 21;
	//All resnums in segment 4 will be reduced by 2
	assembly->add_segment( n_term, second_seg_id, first_res_id , second_res_id );
	return assembly;
}
///@brief Utility function to set up for testing deletes and switches
protocols::sewing::data_storage::SmartAssemblyOP create_simple_non_vital_assembly( protocols::sewing::data_storage::SmartAssemblyOP assembly ){
	//Segments 4-6 will always be considered vital
	assembly->set_starting_segment( assembly->get_segment_vector()->at( 4 ), "all" );
	//Length of segment 4: 9
	//Length of segment 6: 13
	//Length of segment 13: 14
	//Length of segment 21: 7
	//Start with node 4-6
	//Add node 13-15 c-term
	bool n_term = false;
	core::Size first_res_id = 6;
	core::Size second_seg_id = 13; //All residue numbers will now be decreased by 1
	core::Size second_res_id = 7;
	assembly->add_segment( n_term, second_seg_id, first_res_id , second_res_id );
	//Add node 19-21 n-term
	n_term = true;
	//first_res_id is in segment 4
	//second_res_id is in segment 21
	first_res_id = 5;
	second_res_id = 3;
	second_seg_id = 21;
	//All resnums in segment 4 will be reduced by 2
	assembly->add_segment( n_term, second_seg_id, first_res_id , second_res_id );
	protocols::sewing::data_storage::SmartSegmentOP current_segment = assembly->get_n_terminal_segment();
	while ( current_segment ) {
		current_segment->set_is_vital( false );
		current_segment = current_segment->get_c_terminal_neighbor();
	}
	return assembly;
}

///@brief Utility function to make a double chimaera with the first add on the n terminus
protocols::sewing::data_storage::SmartAssemblyOP make_double_chimaera_nterm_first( protocols::sewing::data_storage::SmartAssemblyOP assembly ){
	assembly->set_starting_segment( assembly->get_segment_vector()->at( 31 ), "all" );
	bool nterm = true;
	core::Size segID_2 = 21;
	core::Size resID_1 = 2;
	core::Size resID_2 = 4; //Any residue numbers in seg1 will be increased by 2
	//Add segment 21 to the N terminus
	assembly->add_segment( nterm, segID_2, resID_1, resID_2 );
	//Add segment 1 to the C terminus
	//Segment 1 is 6 residues long
	//Chimaera has 9 residues after the basis residue, but residue numbering has changed (increased by 2)-> was 3 to 11, now 5 to 13
	nterm = false;
	segID_2 = 1;
	resID_1 = 8;
	resID_2 = 3;
	assembly->add_segment( nterm, segID_2, resID_1, resID_2 );
	return assembly;
}
///@brief Utility function to make a double chimaera with the first add on the c terminus
protocols::sewing::data_storage::SmartAssemblyOP  make_double_chimaera_cterm_first( protocols::sewing::data_storage::SmartAssemblyOP assembly ){
	assembly->set_starting_segment( assembly->get_segment_vector()->at( 31 ), "all" );
	bool nterm = false;
	core::Size segID_2 = 1;
	core::Size resID_1 = 6;
	core::Size resID_2 = 3; //Any residue numbers in seg2 will be increased by 3, seg1 unaffected
	//Add segment 1 to the C terminus
	assembly->add_segment( nterm, segID_2, resID_1, resID_2 );
	//Add segment 21 to the N terminus
	//Segment 21 is 7 residues long
	//Chimaera has 5 residues before the basis residue
	nterm = true;
	segID_2 = 21;
	resID_1 = 2;
	resID_2 = 4;
	assembly->add_segment( nterm, segID_2, resID_1, resID_2 );
	return assembly;
}

protocols::sewing::data_storage::SmartAssemblyOP
initial_assembly_with_auto_detected_ligand_contacts( protocols::sewing::data_storage::SmartAssemblyOP assembly){
	using namespace protocols::sewing;
	//First import the pose
	core::pose::Pose pose;
	core::import_pose::pose_from_file( pose, "protocols/sewing/inputs/zinc_site_10_12.pdb" );
	//Then convert into a SmartSegment
	std::map< core::Size, data_storage::SmartSegmentOP > pdbsegs; //This will be empty; we fill it in the function.


	utility::vector1< data_storage::LigandDescription > ligands;
	utility::vector1< data_storage::LigandDescription > expanded_ligands;
	std::map< core::Size, data_storage::LigandResidueCOP > partner_ligands;
	data_storage::LigandDescription ligand;
	ligand.ligand_id = 1;
	ligand.ligand_resnum_string = "26";
	ligand.auto_detect_contacts=true;
	ligands.push_back( ligand );
	std::string pose_segment_starts_string = "";
	std::string pose_segment_ends_string = "";
	std::string pose_segment_dssp = "";

	std::string required_resnums; //Empty b/c we included them in contacts
	core::select::residue_selector::ResidueSelectorCOP required_selector; //see above
	bool strict_dssp_changes = false;
	pdbsegs.clear();
	hashing::AlignmentFileGeneratorMover::add_pose_segments_to_segment_vector(
		pose,
		nullptr,
		assembly->get_segment_vector(),
		pdbsegs,
		pose_segment_starts_string,
		pose_segment_ends_string,
		pose_segment_dssp,
		ligands,
		partner_ligands,
		expanded_ligands,
		required_resnums,
		required_selector,
		strict_dssp_changes
	);
	assembly->pdb_segments( pdbsegs );
	//Now we set this as the starting segment
	data_storage::SmartSegmentOP starting_segment = assembly->local_segments().at( (assembly->pdb_segments().begin() )->second->get_segment_id() );
	assembly->set_starting_segment( starting_segment, "all" );
	return assembly;
}

protocols::sewing::data_storage::SmartAssemblyOP
initial_assembly_with_manually_detected_ligand_contacts( protocols::sewing::data_storage::SmartAssemblyOP assembly ){
	using namespace protocols::sewing;
	//First import the pose
	core::pose::Pose pose;
	core::import_pose::pose_from_file( pose, "protocols/sewing/inputs/zinc_site_10_12.pdb" );
	//Then convert into a SmartSegment
	std::map< core::Size, data_storage::SmartSegmentOP > pdbsegs; //This will be empty; we fill it in the function.



	utility::vector1< data_storage::LigandDescription > ligands;
	utility::vector1< data_storage::LigandDescription > expanded_ligands;
	std::map< core::Size, data_storage::LigandResidueCOP > partner_ligands;
	data_storage::LigandDescription ligand;
	ligand.ligand_id = 1;
	ligand.ligand_resnum_string = "26";
	ligand.auto_detect_contacts=false;

	data_storage::ContactDescription contact1;
	contact1.contact_resnum_string = "10";
	contact1.contact_atom_name="NE2";
	contact1.ligand_atom_name="ZN";
	data_storage::ContactDescription contact2;
	contact2.contact_resnum_string = "14";
	contact2.contact_atom_name="NE2";
	contact2.ligand_atom_name="ZN";
	ligand.ligand_contacts.push_back( contact1 );
	ligand.ligand_contacts.push_back( contact2 );
	ligands.push_back( ligand );
	std::string pose_segment_starts_string = "";
	std::string pose_segment_ends_string = "";
	std::string pose_segment_dssp = "";

	std::string required_resnums; //Empty b/c we included them in contacts
	core::select::residue_selector::ResidueSelectorCOP required_selector; //see above
	bool strict_dssp_changes = false;
	pdbsegs.clear();
	hashing::AlignmentFileGeneratorMover::add_pose_segments_to_segment_vector(
		pose,
		nullptr,
		assembly->get_segment_vector(),
		pdbsegs,
		pose_segment_starts_string,
		pose_segment_ends_string,
		pose_segment_dssp,
		ligands,
		partner_ligands,
		expanded_ligands,
		required_resnums,
		required_selector,
		strict_dssp_changes
	);
	assembly->pdb_segments( pdbsegs );
	//Now we set this as the starting segment
	data_storage::SmartSegmentOP starting_segment = assembly->local_segments().at( (assembly->pdb_segments().begin() )->second->get_segment_id() );
	assembly->set_starting_segment( starting_segment, "all" );
	return assembly;
}

}//end namespace sewing_testing

#endif
