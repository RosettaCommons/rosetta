// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/sewing/data_storage/SmartAssemblyTests.cxxtest.hh
/// @brief  Test sewing's SmartAssembly class
/// @author Minnie Langlois (minnie@email.unc.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

// Project Headers
#include <test/protocols/sewing/extra_functions.hh>
#include <protocols/sewing/data_storage/SmartAssembly.hh>
#include <protocols/sewing/data_storage/SmartSegment.hh>
#include <protocols/sewing/data_storage/LigandSegment.hh>
#include <protocols/sewing/data_storage/LigandResidue.hh>
#include <protocols/sewing/hashing/ModelFileReader.hh>
#include <protocols/sewing/hashing/AlignmentFileGeneratorMover.hh>
#include <protocols/sewing/hashing/hasher_data.hh>
#include <protocols/init_util.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Protocol Headers
#include <basic/Tracer.hh>
#include <utility/string_util.hh>
using namespace protocols::sewing;

static basic::Tracer TR("SmartAssemblyTests");


class SmartAssemblyTests : public CxxTest::TestSuite {
	//Define Variables
private:
	data_storage::SmartAssemblyOP assembly_;
	hashing::SegmentVectorOP segment_vector_;
	core::pose::Pose pose_;
public:

	void setUp(){
		protocols_init_with_additional_options( "-in:auto_setup_metals" );
		segment_vector_ = hashing::ModelFileReader::read_model_file( "protocols/sewing/inputs/test.segments" );
		assembly_ = data_storage::SmartAssemblyOP( new data_storage::SmartAssembly( segment_vector_ ) );
		core::pose::Pose pose_;
		core::import_pose::pose_from_file( pose_, "protocols/sewing/inputs/single_helix.pdb" );
		//SEGMENT 31 IS A SINGLE HELIX
		//ALL OTHERS ARE TRIPLETS (segid%3 = 1 nterm helix, 2 loop, 0 cterm helix)
	}

	void tearDown(){
		assembly_ = data_storage::SmartAssemblyOP( new data_storage::SmartAssembly( segment_vector_ ) );
	}

	void test_set_starting_segment(){
		TS_TRACE( "Begin set_starting_segment" );
		TS_ASSERT( 31 <= segment_vector_->size() );
		assembly_->set_starting_segment( segment_vector_->at( 28 ), "all" );
		TS_ASSERT( assembly_->get_n_terminal_segment() );
		TS_ASSERT_EQUALS( assembly_->get_n_terminal_segment()->get_segment_id(), 28 );
		TS_ASSERT( assembly_->get_c_terminal_segment() );
		TS_ASSERT_EQUALS( assembly_->get_c_terminal_segment()->get_segment_id(), 30 );
	}

	void test_can_make_chimera1(){
		TS_TRACE( "Begin can_make_chimera1" );
		assembly_->set_starting_segment( segment_vector_->at( 28 ), "all" );

		//first basis (c-term)
		//TS_TRACE( "first basis" );
		core::Size first_res_id = 15;
		bool n_term = false;

		//second basis (n-term)
		//TS_TRACE( "second basis" );
		core::Size second_seg_id = 7;
		TS_ASSERT( 7 <= segment_vector_->size() );
		core::Size second_res_id = 2;
		try{
			TS_TRACE( "Trying to Chimerize" );
			bool added = assembly_->add_segment( n_term, second_seg_id, first_res_id , second_res_id ) ;
			//   bool added = assembly_->add_segment( n_term );
			TS_ASSERT( added ); // assert seperately so as not to mess up the assembly trace
			TS_TRACE( "Checking if continuous" );
			TS_ASSERT( assembly_->is_continuous() );
			TS_TRACE( "... Done. Now exporting pose:" );
			core::pose::Pose assem_pose = assembly_->to_pose( "fa_standard" );
			TS_TRACE( "Making sure length is correct:" );
			TS_ASSERT( assem_pose.total_residue() == assembly_->get_length() );
			assem_pose.dump_pdb( "30_7_15_2_test.pdb" );
			TS_TRACE( "...Done.");
		}catch(...){
			TS_FAIL( "Problem with chimerize" );
		}
	}

	void test_can_make_chimera2(){
		TS_TRACE( "Begin can_make_chimera2" );
		assembly_->set_starting_segment( segment_vector_->at( 28 ), "all" );

		//first basis (c-term)
		TS_TRACE( "first basis" );
		core::Size first_res_id = 12;
		bool n_term = false;

		//second basis (n-term)
		TS_TRACE( "second basis" );
		core::Size second_seg_id = 19;
		TS_ASSERT( 19 <= segment_vector_->size() );
		core::Size second_res_id = 7;
		try{
			TS_TRACE( "Trying to Chimerize" );
			bool added = assembly_->add_segment( n_term, second_seg_id, first_res_id , second_res_id ) ;
			//   bool added = assembly_->add_segment( n_term );
			TS_ASSERT( added ); // assert seperately so as not to mess up the assembly trace
			TS_TRACE( "Checking if continuous" );
			TS_ASSERT( assembly_->is_continuous() );
			TS_TRACE( "... Done. Now exporting pose:" );
			core::pose::Pose assem_pose = assembly_->to_pose( "fa_standard" );
			TS_TRACE( "Making sure length is correct:" );
			TS_ASSERT( assem_pose.total_residue() == assembly_->get_length() );
			assem_pose.dump_pdb( "30_19_12_7_test.pdb" );
			TS_TRACE( "...Done.");
		}catch(...){
			TS_FAIL( "Problem with chimerize" );
		}
	}

	void test_can_make_chimera3(){
		TS_TRACE( "Begin can_make_chimera3" );
		assembly_->set_starting_segment( segment_vector_->at( 13 ), "all" );

		//first basis (c-term)
		TS_TRACE( "first basis" );
		core::Size first_res_id = 8;
		bool n_term = false;

		//second basis (n-term)
		TS_TRACE( "second basis" );
		core::Size second_seg_id = 7;
		TS_ASSERT( 7 <= segment_vector_->size() );
		core::Size second_res_id = 5;
		try{
			TS_TRACE( "Trying to Chimerize" );
			bool added = assembly_->add_segment( n_term, second_seg_id, first_res_id , second_res_id ) ;
			//   bool added = assembly_->add_segment( n_term );
			TS_ASSERT( added ); // assert seperately so as not to mess up the assembly trace
			TS_TRACE( "Checking if continuous" );
			TS_ASSERT( assembly_->is_continuous() );
			TS_TRACE( "... Done. Now exporting pose:" );
			core::pose::Pose assem_pose = assembly_->to_pose( "fa_standard" );
			TS_TRACE( "Making sure length is correct:" );
			TS_ASSERT( assem_pose.total_residue() == assembly_->get_length() );
			assem_pose.dump_pdb( "15_7_8_5_test.pdb" );
			TS_TRACE( "...Done.");
		}catch(...){
			TS_FAIL( "Problem with chimerize" );
		}
	}

	void test_simple_assembly_size(){
		TS_TRACE( "Begin simple_assembly_size" );
		assembly_ = sewing_testing::create_simple_assembly( assembly_ );
		TS_TRACE( assembly_->get_forward_assembly() );
		TS_ASSERT( assembly_->get_size() == 7 );
	}

	void test_simple_assembly_length(){
		TS_TRACE( "Begin simple_assembly_length" );
		assembly_ = sewing_testing::create_simple_assembly( assembly_ );
		TS_TRACE( assembly_->get_forward_assembly() );
		core::Size assembly_length = 0;
		data_storage::SmartSegmentOP current_seg = assembly_->get_n_terminal_segment();
		while ( current_seg != nullptr ) {
			assembly_length = assembly_length +  current_seg->get_length();
			current_seg = current_seg->get_c_terminal_neighbor();
		}
		TS_TRACE( "Reported assembly length: " + utility::to_string( assembly_->get_length() ) );
		TS_TRACE( "Calculated assembly length: " + utility::to_string( assembly_length ) );
		TS_ASSERT( assembly_->get_length() == assembly_length );
	}

	void test_simple_assembly_parent_connectivity(){
		TS_TRACE( "Begin simple_assembly_parent_connectivity" );
		assembly_ = sewing_testing::create_simple_assembly( assembly_ );
		TS_TRACE( assembly_->get_forward_assembly() );
		//Chimerized 4 and 21, 6 and 13
		TS_ASSERT( data_storage::SmartSegment::get_c_most_segment( assembly_->get_n_terminal_segment(), false )->is_chimaeric() );
		TS_ASSERT( data_storage::SmartSegment::get_c_most_segment( assembly_->get_n_terminal_segment(), true )->get_segment_id() == 15 );
		TS_ASSERT( data_storage::SmartSegment::get_n_most_segment( assembly_->get_c_terminal_segment(), false )->is_chimaeric() );
		TS_ASSERT( data_storage::SmartSegment::get_n_most_segment( assembly_->get_c_terminal_segment(), true )->get_segment_id() == 19 );
		//Now for the middle part
		TS_ASSERT( data_storage::SmartSegment::get_c_most_segment( assembly_->local_segments().at( 5 ), false )->is_chimaeric( ) );
		TS_ASSERT( data_storage::SmartSegment::get_c_most_segment( assembly_->local_segments().at( 5 ), true )->get_segment_id() == 15 );
		TS_ASSERT( data_storage::SmartSegment::get_n_most_segment( assembly_->local_segments().at( 5 ), false )->is_chimaeric() );
		TS_ASSERT( data_storage::SmartSegment::get_n_most_segment( assembly_->local_segments().at( 5 ), true )->get_segment_id() == 19 );
		TS_ASSERT( assembly_->local_segments().at( 4 )->get_c_terminal_neighbor()->get_segment_id() == 5 );
		TS_ASSERT( assembly_->local_segments().at( 6 )->get_n_terminal_neighbor()->get_segment_id() == 5 );
		TS_ASSERT( assembly_->local_segments().at( 13 )->get_c_terminal_neighbor()->get_segment_id() == 14 );
		TS_ASSERT( assembly_->local_segments().at( 21 )->get_n_terminal_neighbor()->get_segment_id() == 20 );

		//Check chimaera connectivity
		TS_ASSERT( assembly_->get_forward_assembly() == assembly_->get_reverse_assembly() );
	}

	void test_simple_delete_nterm(){
		TS_TRACE( "Begin simple_delete_nterm" );
		assembly_ = sewing_testing::create_simple_assembly( assembly_ );
		bool nterm = true;
		assembly_->delete_segment( nterm );
		TS_TRACE( assembly_->get_forward_assembly() );
		//Chimaera is 6 an 13
		TS_ASSERT( assembly_->get_last_change_was_n_terminal() );
		TS_ASSERT( assembly_->get_size() == 5 );
		TS_ASSERT( assembly_->get_n_terminal_segment()->get_segment_id() == 4 );
		TS_ASSERT( data_storage::SmartSegment::get_c_most_segment( assembly_->get_n_terminal_segment(), false )->is_chimaeric() );
		TS_ASSERT( data_storage::SmartSegment::get_c_most_segment( assembly_->get_n_terminal_segment(), true )->get_segment_id() == 15 );
		TS_ASSERT( data_storage::SmartSegment::get_n_most_segment( assembly_->get_c_terminal_segment(), false )->is_chimaeric() );
		TS_ASSERT( data_storage::SmartSegment::get_n_most_segment( assembly_->get_c_terminal_segment(), true )->get_segment_id() == 4 );
		TS_ASSERT( assembly_->local_segments().at( 6 )->get_n_terminal_neighbor()->get_segment_id() == 5 );
		TS_ASSERT( assembly_->local_segments().at( 13 )->get_c_terminal_neighbor()->get_segment_id() == 14 );
		TS_ASSERT( assembly_->get_forward_assembly() == assembly_->get_reverse_assembly() );

		//test connectivity of deleted part (19-21)
		TS_ASSERT( data_storage::SmartSegment::get_c_most_segment( assembly_->local_segments().at( 19), true )->get_segment_id() == 21 );
		TS_ASSERT( data_storage::SmartSegment::get_c_most_segment( assembly_->local_segments().at( 19), false )->get_segment_id() == 21 );
		TS_ASSERT( data_storage::SmartSegment::get_n_most_segment( assembly_->local_segments().at( 21), true )->get_segment_id() == 19 );
		TS_ASSERT( data_storage::SmartSegment::get_n_most_segment( assembly_->local_segments().at( 21), false )->get_segment_id() == 19 );
	}

	void test_simple_delete_cterm(){
		TS_TRACE( "Begin simple_delete_cterm" );
		assembly_ = sewing_testing::create_simple_assembly( assembly_ );
		bool nterm = false;
		assembly_->delete_segment( nterm );
		TS_TRACE( assembly_->get_forward_assembly() );
		//Chimaera is 4 and 21
		TS_ASSERT( !assembly_->get_last_change_was_n_terminal() );
		TS_ASSERT( assembly_->get_size() == 5 );
		TS_ASSERT( assembly_->get_c_terminal_segment()->get_segment_id() == 6 );
		TS_ASSERT( assembly_->local_segments().at( 4 )->get_c_terminal_neighbor()->get_segment_id() == 5 );
		TS_ASSERT( assembly_->local_segments().at( 21 )->get_n_terminal_neighbor()->get_segment_id() == 20 );
		TS_ASSERT( data_storage::SmartSegment::get_c_most_segment( assembly_->get_n_terminal_segment(), false )->is_chimaeric() );
		TS_ASSERT( data_storage::SmartSegment::get_c_most_segment( assembly_->get_n_terminal_segment(), true )->get_segment_id() == 6 );
		TS_ASSERT( data_storage::SmartSegment::get_n_most_segment( assembly_->get_c_terminal_segment(), false )->is_chimaeric() );
		TS_ASSERT( data_storage::SmartSegment::get_n_most_segment( assembly_->get_c_terminal_segment(), true )->get_segment_id() == 19 );
		TS_ASSERT( assembly_->get_forward_assembly() == assembly_->get_reverse_assembly() );


		//test connectivity of deleted part (13-15)
		TS_ASSERT( data_storage::SmartSegment::get_c_most_segment( assembly_->local_segments().at( 13), true )->get_segment_id() == 15 );
		TS_ASSERT( data_storage::SmartSegment::get_c_most_segment( assembly_->local_segments().at( 13), false )->get_segment_id() == 15 );
		TS_ASSERT( data_storage::SmartSegment::get_n_most_segment( assembly_->local_segments().at( 15), true )->get_segment_id() == 13 );
		TS_ASSERT( data_storage::SmartSegment::get_n_most_segment( assembly_->local_segments().at( 15), false )->get_segment_id() == 13 );
	}


	void test_do_not_delete_vital_segments_nterm(){
		TS_TRACE( "Begin do_not_delete_vital_segments_nterm" );
		assembly_ = sewing_testing::create_simple_assembly( assembly_ );
		bool nterm = true;
		bool deleted_first = assembly_->delete_segment( nterm );
		TS_ASSERT( deleted_first );
		bool deleted_second = assembly_->delete_segment( nterm );
		TS_ASSERT( !deleted_second );
	}
	void test_do_not_delete_vital_segments_cterm(){
		TS_TRACE( "Begin do_not_delete_vital_segments_cterm" );
		assembly_ = sewing_testing::create_simple_assembly( assembly_ );
		bool nterm = false;
		bool deleted_first = assembly_->delete_segment( nterm );
		TS_ASSERT( deleted_first );
		bool deleted_second = assembly_->delete_segment( nterm );
		TS_ASSERT( !deleted_second );
	}


	void test_do_not_hash_the_unhashable(){
		TS_TRACE( "Begin do_not_hash_the_unhashable" );
		//Make sure it won't hash a helix against a loop
		//Function should return 0
		core::Size bpnum;
		bpnum = ( assembly_->iterate_over_basis_pairs( assembly_->get_segment_vector()->at( 4 ),  assembly_->get_segment_vector()->at( 11 ), true ) ).second;
		TS_ASSERT( bpnum == 0 );
		bpnum = ( assembly_->iterate_over_basis_pairs( assembly_->get_segment_vector()->at( 6 ),  assembly_->get_segment_vector()->at( 11 ), false ) ).second;
		TS_ASSERT( bpnum == 0 );
		bpnum = ( assembly_->iterate_over_basis_pairs( assembly_->get_segment_vector()->at( 1 ),  assembly_->get_segment_vector()->at( 11 ), true ) ).second;
		TS_ASSERT( bpnum == 0 );
		bpnum = ( assembly_->iterate_over_basis_pairs( assembly_->get_segment_vector()->at( 30 ),  assembly_->get_segment_vector()->at( 11 ), false ) ).second;
		TS_ASSERT( bpnum == 0 );
		bpnum = ( assembly_->iterate_over_basis_pairs( assembly_->get_segment_vector()->at( 11 ),  assembly_->get_segment_vector()->at( 4 ), false ) ).second;
		TS_ASSERT( bpnum == 0 );
		bpnum = ( assembly_->iterate_over_basis_pairs( assembly_->get_segment_vector()->at( 11 ),  assembly_->get_segment_vector()->at( 6 ), true ) ).second;
		TS_ASSERT( bpnum == 0 );
		bpnum = ( assembly_->iterate_over_basis_pairs( assembly_->get_segment_vector()->at( 11 ),  assembly_->get_segment_vector()->at( 1 ), false ) ).second;
		TS_ASSERT( bpnum == 0 );
		bpnum = ( assembly_->iterate_over_basis_pairs( assembly_->get_segment_vector()->at( 11 ),  assembly_->get_segment_vector()->at( 30 ), true ) ).second;
		TS_ASSERT( bpnum == 0 );
	}




	void test_iterate_over_basis_pairs_basic(){
		TS_TRACE( "Begin iterate_over_basis_pairs_basic" );
		//Make sure it returns the right number of basis pairs
		//Make sure all the basis pairs it returns are in range (TODO)
		//Length of segment 4: 9
		//Length of segment 21: 7
		//Length of segment 1: 6
		//Length of segment 30: 24
		core::Size bpnum;
		//4 and 21 (and 21 and 4)
		bpnum = ( assembly_->iterate_over_basis_pairs( assembly_->get_segment_vector()->at( 4 ), assembly_->get_segment_vector()->at( 21 ), true ) ).second;
		TS_ASSERT( bpnum == (9 * 7) );
		bpnum = ( assembly_->iterate_over_basis_pairs( assembly_->get_segment_vector()->at( 21 ), assembly_->get_segment_vector()->at( 4 ), false ) ).second;
		TS_ASSERT( bpnum == (9 * 7) );
		//4 and 30 (and 30 and 4)
		bpnum = ( assembly_->iterate_over_basis_pairs( assembly_->get_segment_vector()->at( 4 ), assembly_->get_segment_vector()->at( 30 ), true ) ).second;
		TS_ASSERT( bpnum == (9 * 24) );
		bpnum = ( assembly_->iterate_over_basis_pairs( assembly_->get_segment_vector()->at( 30 ), assembly_->get_segment_vector()->at( 4 ), false ) ).second;
		TS_ASSERT( bpnum == (9 * 24) );
		//1 and 21 (and 21 and 1)
		bpnum = ( assembly_->iterate_over_basis_pairs( assembly_->get_segment_vector()->at( 1 ), assembly_->get_segment_vector()->at( 21 ), true ) ).second;
		TS_ASSERT( bpnum == (6 * 7) );
		bpnum = ( assembly_->iterate_over_basis_pairs( assembly_->get_segment_vector()->at( 21 ), assembly_->get_segment_vector()->at( 1 ), false ) ).second;
		TS_ASSERT( bpnum == (6 * 7) );
		//1 and 30 (and 30 and 1)
		bpnum = ( assembly_->iterate_over_basis_pairs( assembly_->get_segment_vector()->at( 1 ), assembly_->get_segment_vector()->at( 30 ), true ) ).second;
		TS_ASSERT( bpnum == (6 * 24) );
		bpnum = ( assembly_->iterate_over_basis_pairs( assembly_->get_segment_vector()->at( 30 ), assembly_->get_segment_vector()->at( 1 ), false ) ).second;
		TS_ASSERT( bpnum == (6 * 24) );
	}

	void test_simple_switch_nterm(){
		assembly_ = sewing_testing::create_simple_assembly( assembly_ );
		TS_TRACE( "Begin simple_switch_nterm" );
		TS_TRACE( assembly_->get_forward_assembly() );
		bool nterm = true;
		//Unfortunately switch currently does not have a non-random version
		//This will rely on iterate_over_basis_pairs
		assembly_->switch_segment( nterm );
		TS_TRACE( assembly_->get_forward_assembly() );
		TS_TRACE( assembly_->get_reverse_assembly() );
		//This should not fail with our starting assembly
		TS_ASSERT( assembly_->get_last_change_was_n_terminal() );
		//Make sure the length was handled properly
		core::Size assembly_length = 0;
		data_storage::SmartSegmentOP current_seg = assembly_->get_n_terminal_segment();
		while ( current_seg != nullptr ) {
			//assembly_length += current_seg->get_length();
			assembly_length = assembly_length +  current_seg->get_length();
			current_seg = current_seg->get_c_terminal_neighbor();
		}
		TS_TRACE( "Reported assembly length: " + utility::to_string( assembly_->get_length() ) );
		TS_TRACE( "Calculated assembly length: " + utility::to_string( assembly_length ) );
		TS_ASSERT( assembly_->get_length() == assembly_length );
		//Make sure the assembly reads the same forward and backward
		TS_ASSERT( assembly_->get_forward_assembly() == assembly_->get_reverse_assembly() );

	}

	void test_simple_switch_cterm(){
		TS_TRACE( "Begin simple_switch_cterm" );
		assembly_ = sewing_testing::create_simple_assembly( assembly_ );
		TS_TRACE( assembly_->get_forward_assembly() );
		bool nterm = false;
		//Unfortunately switch currently does not have a non-random version
		assembly_->switch_segment( nterm );
		TS_TRACE( assembly_->get_forward_assembly() );
		//This should not fail with our starting assembly
		TS_ASSERT( !assembly_->get_last_change_was_n_terminal() );
		//Make sure the length was handled properly
		core::Size assembly_length = 0;
		data_storage::SmartSegmentOP current_seg = assembly_->get_n_terminal_segment();
		while ( current_seg != nullptr ) {
			assembly_length = assembly_length +  current_seg->get_length();
			//assembly_length += current_seg->get_length();
			current_seg = current_seg->get_c_terminal_neighbor();
		}
		TS_TRACE( "Reported assembly length: " + utility::to_string( assembly_->get_length() ) );
		TS_TRACE( "Calculated assembly length: " + utility::to_string( assembly_length ) );
		TS_ASSERT( assembly_->get_length() == assembly_length );
		//Make sure the assembly reads the same forward and backward
		TS_ASSERT( assembly_->get_forward_assembly() == assembly_->get_reverse_assembly() );

	}

	void test_iterate_over_basis_pairs_chimaera_cterm_first(){
		TS_TRACE( "Begin iterate_over_basis_pairs_chimaera_cterm_first" );
		//Length of segment 4: 9
		//Length of segment 6: 13
		//Length of segment 1: 6
		//Length of segment 21: 7
		//Length of segment 30: 24
		//Length of segment 13: 14
		//Length of segment 31: 11
		//First make a basic chimaera between our single segment and a triplet to the c terminus
		assembly_->set_starting_segment( assembly_->get_segment_vector()->at( 31 ), "all" );
		bool nterm = false;
		core::Size segID_2 = 1;
		core::Size resID_1 = 6;
		core::Size resID_2 = 3; //Any residue numbers in seg2 will be increased by 3
		//Add segment 1 to the C terminus
		assembly_->add_segment( nterm, segID_2, resID_1, resID_2 );
		TS_TRACE( assembly_->get_forward_assembly() );
		TS_ASSERT( assembly_->get_n_terminal_segment()->is_chimaeric() );
		//Make sure parent connectivity is still maintained
		//Chimaera is 31 and 1
		TS_ASSERT( assembly_->local_segments().at( 1 )->get_c_terminal_neighbor()->get_segment_id() == 2 );

		TS_ASSERT( data_storage::SmartSegment::get_c_most_segment( assembly_->get_n_terminal_segment(), false )->get_segment_id() == assembly_->get_n_terminal_segment()->get_segment_id() );
		TS_ASSERT( data_storage::SmartSegment::get_c_most_segment( assembly_->get_n_terminal_segment(), true )->get_segment_id() == 3 );
		TS_ASSERT( data_storage::SmartSegment::get_n_most_segment( assembly_->get_c_terminal_segment(), false )->get_segment_id() == assembly_->get_n_terminal_segment()->get_segment_id() );
		TS_ASSERT( data_storage::SmartSegment::get_n_most_segment( assembly_->get_c_terminal_segment(), true )->get_segment_id() == assembly_->get_n_terminal_segment()->get_segment_id() );
		TS_ASSERT( assembly_->get_forward_assembly() == assembly_->get_reverse_assembly() );
		//Then call function on n_terminal_segment and a new triplet
		core::Size bpnum;
		bpnum = ( assembly_->iterate_over_basis_pairs( assembly_->get_n_terminal_segment(), assembly_->get_segment_vector()->at( 21 ), true ) ).second;
		//Segment 21 is 7 residues long
		//Chimaera has 5 residues before the basis residue
		TS_ASSERT( bpnum == (7 * 5 ) );
	}
	void test_iterate_over_basis_pairs_chimaera_nterm_first(){
		TS_TRACE( "Begin iterate_over_basis_pairs_chimaera_nterm_first" );
		//First make a basic chimaera between our single segment and a triplet to the n terminus
		//Then call function on c_terminal_segment and a new triplet
		//Length of segment 4: 9
		//Length of segment 6: 13
		//Length of segment 1: 6
		//Length of segment 21: 7
		//Length of segment 30: 24
		//Length of segment 13: 14
		//Length of segment 31: 11
		//First make a basic chimaera between our single segment and a triplet to the c terminus
		assembly_->set_starting_segment( assembly_->get_segment_vector()->at( 31 ), "all" );
		bool nterm = true;
		core::Size segID_2 = 21;
		core::Size resID_1 = 6;
		core::Size resID_2 = 4; //Any residue numbers in seg1 will be decreased by 2

		//Add segment 21 to the N terminus
		assembly_->add_segment( nterm, segID_2, resID_1, resID_2 );
		TS_TRACE( assembly_->get_forward_assembly() );
		TS_ASSERT( assembly_->get_c_terminal_segment()->is_chimaeric() );
		TS_ASSERT_EQUALS( assembly_->get_c_terminal_segment()->get_basis_pair().first.segment_id(), 31 );
		TS_ASSERT_EQUALS( assembly_->get_c_terminal_segment()->get_basis_pair().second.segment_id(), 21 );
		//Make sure parent connectivity is still maintained
		//Chimaera is 21 and 31
		TS_ASSERT( assembly_->local_segments().at( 21 )->get_n_terminal_neighbor()->get_segment_id() == 20 );
		TS_ASSERT( data_storage::SmartSegment::get_c_most_segment( assembly_->get_n_terminal_segment(), false )->get_segment_id() == assembly_->get_c_terminal_segment()->get_segment_id() );
		TS_ASSERT( data_storage::SmartSegment::get_c_most_segment( assembly_->get_n_terminal_segment(), true )->get_segment_id() == assembly_->get_c_terminal_segment()->get_segment_id() );
		TS_ASSERT( data_storage::SmartSegment::get_n_most_segment( assembly_->get_c_terminal_segment(), false )->get_segment_id() == assembly_->get_c_terminal_segment()->get_segment_id() );
		TS_ASSERT( data_storage::SmartSegment::get_n_most_segment( assembly_->get_c_terminal_segment(), true )->get_segment_id() == 19 );
		TS_ASSERT( assembly_->get_forward_assembly() == assembly_->get_reverse_assembly() );
		TS_TRACE( assembly_->get_comprehensive_forward_assembly() );
		//Then call function on n_terminal_segment and a new triplet
		core::Size bpnum;
		bpnum = ( assembly_->iterate_over_basis_pairs( assembly_->get_c_terminal_segment(), assembly_->get_segment_vector()->at( 1 ), false ) ).second;
		//Segment 1 is 6 residues long
		//Chimaera has 5 residues after the basis residue
		TS_TRACE( "Predicted number of available bp: 30" );
		TS_TRACE( "Calculated number: " + utility::to_string( bpnum ) );

		TS_ASSERT( bpnum == ( 6 * 5 ) );
	}

	void test_double_chimerize_single_helix_two_triplets_nterm_first(){
		TS_TRACE( "Begin double_chimerize_single_helix_two_triplets_nterm_first" );
		assembly_ = sewing_testing::make_double_chimaera_nterm_first( assembly_ );
		TS_TRACE( assembly_->get_forward_assembly() );
		//Test size
		TS_ASSERT( assembly_->get_size() == 5 );
		//Test length
		core::Size assembly_length = 0;
		data_storage::SmartSegmentOP current_seg = assembly_->get_n_terminal_segment();
		while ( current_seg != nullptr ) {
			assembly_length = assembly_length +  current_seg->get_length();
			//assembly_length += current_seg->get_length();
			current_seg = current_seg->get_c_terminal_neighbor();
		}
		TS_TRACE( "Reported assembly length: " + utility::to_string( assembly_->get_length() ) );
		TS_TRACE( "Calculated assembly length: " + utility::to_string( assembly_length ) );
		TS_ASSERT( assembly_->get_length() == assembly_length );
		//Test all connectivity (TODO)
		//19-20-(21|31||1)-2-3
		TS_ASSERT( data_storage::SmartSegment::get_c_most_segment( assembly_->get_n_terminal_segment(), false )->is_chimaeric() );
		TS_ASSERT( data_storage::SmartSegment::get_c_most_segment( assembly_->get_n_terminal_segment(), true )->get_segment_id() == 3 );
		TS_ASSERT( data_storage::SmartSegment::get_n_most_segment( assembly_->get_c_terminal_segment(), false )->is_chimaeric() );
		TS_ASSERT( data_storage::SmartSegment::get_n_most_segment( assembly_->get_c_terminal_segment(), true )->get_segment_id() == 19 );
		//Now for the middle part
		data_storage::SmartSegmentOP chimaera = assembly_->local_segments().at( 2 )->get_n_terminal_neighbor(); //Double chimaera, actually
		TS_ASSERT( chimaera->is_chimaeric() );
		//Make sure parents' connections to neighbors are maintained
		TS_ASSERT( chimaera->get_c_terminal_parent()->get_segment_id() == 1 );
		TS_ASSERT( chimaera->get_c_terminal_parent()->get_c_terminal_neighbor()->get_segment_id() == 2 );

		TS_ASSERT( chimaera->get_n_terminal_parent()->is_chimaeric() );

		TS_ASSERT( chimaera->get_n_terminal_parent()->get_n_terminal_parent()->get_segment_id() == 21 );
		TS_ASSERT( chimaera->get_n_terminal_parent()->get_c_terminal_parent()->get_segment_id() == 31 );
		TS_ASSERT( chimaera->get_n_terminal_parent()->get_c_terminal_parent()->get_c_terminal_neighbor() == nullptr );
		TS_ASSERT( chimaera->get_n_terminal_parent()->get_n_terminal_parent()->get_n_terminal_neighbor()->get_segment_id() == 20 );
		//Check chimaera connectivity
		TS_ASSERT( assembly_->get_forward_assembly() == assembly_->get_reverse_assembly() );
	}

	void test_double_chimerize_single_helix_two_triplets_cterm_first(){
		TS_TRACE( "Begin double_chimerize_single_helix_two_triplets_cterm_first" );
		assembly_ = sewing_testing::make_double_chimaera_cterm_first( assembly_ );
		TS_TRACE( assembly_->get_forward_assembly() );
		//Test size
		TS_ASSERT( assembly_->get_size() == 5 );
		//Test length
		core::Size assembly_length = 0;
		data_storage::SmartSegmentOP current_seg = assembly_->get_n_terminal_segment();
		while ( current_seg != nullptr ) {
			assembly_length = assembly_length +  current_seg->get_length();
			//assembly_length += current_seg->get_length();
			current_seg = current_seg->get_c_terminal_neighbor();
		}
		TS_TRACE( "Reported assembly length: " + utility::to_string( assembly_->get_length() ) );
		TS_TRACE( "Calculated assembly length: " + utility::to_string( assembly_length ) );
		TS_ASSERT( assembly_->get_length() == assembly_length );
		//Test all connectivity (TODO)
		//19-20-(21||31|1)-2-3
		TS_ASSERT( data_storage::SmartSegment::get_c_most_segment( assembly_->get_n_terminal_segment(), false )->is_chimaeric() );
		TS_ASSERT( data_storage::SmartSegment::get_c_most_segment( assembly_->get_n_terminal_segment(), true )->get_segment_id() == 3 );
		TS_ASSERT( data_storage::SmartSegment::get_n_most_segment( assembly_->get_c_terminal_segment(), false )->is_chimaeric() );
		TS_ASSERT( data_storage::SmartSegment::get_n_most_segment( assembly_->get_c_terminal_segment(), true )->get_segment_id() == 19 );
		//Now for the middle part
		data_storage::SmartSegmentOP chimaera = assembly_->local_segments().at( 2 )->get_n_terminal_neighbor(); //Double chimaera, actually
		TS_ASSERT( chimaera->is_chimaeric() );
		//Make sure parents' connections to neighbors are maintained
		TS_ASSERT( chimaera->get_n_terminal_parent()->get_n_terminal_neighbor()->get_segment_id() == 20 );
		TS_ASSERT( chimaera->get_c_terminal_parent()->is_chimaeric() );
		TS_ASSERT( chimaera->get_c_terminal_parent()->get_n_terminal_parent()->get_segment_id() == 31 );
		TS_ASSERT( chimaera->get_c_terminal_parent()->get_c_terminal_parent()->get_segment_id() == 1 );
		TS_ASSERT( chimaera->get_c_terminal_parent()->get_c_terminal_parent()->get_c_terminal_neighbor()->get_segment_id() == 2 );
		//Check chimaera connectivity
		TS_ASSERT( assembly_->get_forward_assembly() == assembly_->get_reverse_assembly() );

	}

	void test_basis_pair_calculations_recursive_unchimerize_n_term(){
		TS_TRACE( "Begin test_basis_pair_calculations_recursive_unchimerize_n_term" );
		assembly_ = sewing_testing::make_double_chimaera_nterm_first( assembly_ );
		bool nterm = true;
		assembly_->delete_segment( nterm );
		TS_TRACE( assembly_->get_forward_assembly() );
		//First check that the chimaera in the assembly has the correct basis pair
		//Correct answer will be (6, 3)
		core::Size seg1 = assembly_->get_n_terminal_segment()->get_basis_pair().first.segment_id();
		core::Size seg2 = assembly_->get_n_terminal_segment()->get_basis_pair().second.segment_id();
		core::Size res1 = assembly_->get_n_terminal_segment()->get_basis_pair().first.resnum();
		core::Size res2 = assembly_->get_n_terminal_segment()->get_basis_pair().second.resnum();
		TS_ASSERT_EQUALS( seg1, 31 );
		TS_ASSERT_EQUALS( seg2, 1 );
		TS_ASSERT_EQUALS( res1, 6 );
		TS_ASSERT_EQUALS( res2, 3 );


		//Then check that last_chimaera_deleted has the basis pair that we would need it to use to add it back
		//Correct answer will be (2, 4)
		core::Size deleted_seg1 = assembly_->get_last_chimaera_deleted().first.segment_id();
		core::Size deleted_seg2 = assembly_->get_last_chimaera_deleted().second.segment_id();
		core::Size deleted_res1 = assembly_->get_last_chimaera_deleted().first.resnum();
		core::Size deleted_res2 = assembly_->get_last_chimaera_deleted().second.resnum();
		TS_ASSERT_EQUALS( deleted_seg1, 31 );
		TS_ASSERT_EQUALS( deleted_seg2, 21 );
		TS_ASSERT_EQUALS( deleted_res1, 2 );
		TS_ASSERT_EQUALS( deleted_res2, 4 );
	}







	void test_basis_pair_calculations_recursive_unchimerize_c_term(){
		TS_TRACE( "Begin test_basis_pair_calculations_recursive_unchimerize_c_term" );
		assembly_ = sewing_testing::make_double_chimaera_cterm_first( assembly_ );
		TS_TRACE( assembly_->get_forward_assembly() );
		TS_TRACE( "Testing connectivity of parent chimaera" );
		data_storage::SmartSegmentOP chimaera = assembly_->local_segments().at( 2 )->get_n_terminal_neighbor(); //Double chimaera, actually
		TS_ASSERT( chimaera->get_n_terminal_parent()->get_segment_id() == 21 );
		TS_ASSERT( chimaera->get_n_terminal_parent()->get_n_terminal_neighbor()->get_segment_id() == 20 );
		TS_ASSERT( chimaera->get_c_terminal_parent()->is_chimaeric() );
		TS_ASSERT( chimaera->get_c_terminal_parent()->get_n_terminal_parent()->get_segment_id() == 31 );
		TS_ASSERT( chimaera->get_c_terminal_parent()->get_c_terminal_parent()->get_segment_id() == 1 );
		TS_ASSERT( chimaera->get_c_terminal_parent()->get_c_terminal_parent()->get_c_terminal_neighbor()->get_segment_id() == 2 );
		TS_TRACE( "Testing connectivity of parent chimaera" );
		TS_ASSERT( chimaera->get_c_terminal_parent()->get_c_terminal_neighbor() );
		TS_ASSERT( chimaera->get_c_terminal_parent()->get_c_terminal_neighbor()->get_segment_id() == 2 );
		bool nterm = false;
		assembly_->delete_segment( nterm );
		TS_TRACE( assembly_->get_forward_assembly() );
		//First check that the chimaera in the assembly has the correct basis pair
		//Correct answer will be (2, 4)
		core::Size seg1 = assembly_->get_c_terminal_segment()->get_basis_pair().first.segment_id();
		core::Size seg2 = assembly_->get_c_terminal_segment()->get_basis_pair().second.segment_id();
		core::Size res1 = assembly_->get_c_terminal_segment()->get_basis_pair().first.resnum();
		core::Size res2 = assembly_->get_c_terminal_segment()->get_basis_pair().second.resnum();
		TS_ASSERT_EQUALS( seg1, 31 );
		TS_ASSERT_EQUALS( seg2, 21 );
		TS_ASSERT_EQUALS( res1, 2 );
		TS_ASSERT_EQUALS( res2, 4 );
		//Then check that last_chimaera_deleted has the basis pair that we would need it to use to add it back
		//Correct answer will be (8, 3)
		core::Size deleted_seg1 = assembly_->get_last_chimaera_deleted().first.segment_id();
		core::Size deleted_seg2 = assembly_->get_last_chimaera_deleted().second.segment_id();
		core::Size deleted_res1 = assembly_->get_last_chimaera_deleted().first.resnum();
		core::Size deleted_res2 = assembly_->get_last_chimaera_deleted().second.resnum();
		TS_ASSERT_EQUALS( deleted_seg1, assembly_->get_c_terminal_segment()->get_segment_id() );
		TS_ASSERT_EQUALS( deleted_seg2, 1 );
		TS_ASSERT_EQUALS( deleted_res1, 8 );
		TS_ASSERT_EQUALS( deleted_res2, 3 );

	}


	void test_double_unchimerize_single_helix_two_triplets_nterm(){
		TS_TRACE( "Begin double_unchimerize_single_helix_two_triplets_nterm" );
		assembly_ = sewing_testing::make_double_chimaera_nterm_first( assembly_ );
		bool nterm = true;
		assembly_->delete_segment( nterm );
		TS_TRACE( assembly_->get_forward_assembly() );
		//Now we're going to test all the parameters of the assembly
		TS_ASSERT( assembly_->get_size() == 3 );
		//Test length
		core::Size assembly_length = 0;
		data_storage::SmartSegmentOP current_seg = assembly_->get_n_terminal_segment();
		while ( current_seg != nullptr ) {
			//assembly_length += current_seg->get_length();
			assembly_length = assembly_length +  current_seg->get_length();
			current_seg = current_seg->get_c_terminal_neighbor();
		}
		TS_TRACE( "Reported assembly length: " + utility::to_string( assembly_->get_length() ) );
		TS_TRACE( "Calculated assembly length: " + utility::to_string( assembly_length ) );
		TS_ASSERT( assembly_->get_length() == assembly_length );
		//Make sure the termini are right
		TS_ASSERT( assembly_->get_c_terminal_segment()->get_segment_id() == 3 );
		TS_ASSERT( assembly_->get_n_terminal_segment()->is_chimaeric() );
		TS_ASSERT( assembly_->get_n_terminal_segment()->get_c_terminal_parent()->get_segment_id() == 1 );
		TS_ASSERT( assembly_->get_n_terminal_segment()->get_n_terminal_parent()->get_segment_id() == 31 );
		//Test all connectivity
		//(31|1)-2-3
		TS_ASSERT( data_storage::SmartSegment::get_c_most_segment( assembly_->get_n_terminal_segment(), false )->is_chimaeric() );
		TS_ASSERT( data_storage::SmartSegment::get_c_most_segment( assembly_->get_n_terminal_segment(), true )->get_segment_id() == 3 );
		TS_ASSERT( data_storage::SmartSegment::get_n_most_segment( assembly_->get_c_terminal_segment(), false )->is_chimaeric() );
		TS_ASSERT( data_storage::SmartSegment::get_n_most_segment( assembly_->get_c_terminal_segment(), true )->is_chimaeric() );
		data_storage::SmartSegmentOP chimaera = assembly_->get_n_terminal_segment();
		TS_ASSERT( chimaera->is_chimaeric() );
		//Make sure parents' connections to neighbors are maintained
		TS_ASSERT( chimaera->get_n_terminal_parent()->get_segment_id() == 31 );
		TS_ASSERT( chimaera->get_n_terminal_parent()->get_c_terminal_neighbor() == nullptr );
		TS_ASSERT( chimaera->get_n_terminal_parent()->get_n_terminal_neighbor() == nullptr );
		TS_ASSERT( chimaera->get_c_terminal_parent()->get_c_terminal_neighbor()->get_segment_id() == 2 );

		//Make sure the deleted part's connectivity was restored (19-21
		TS_ASSERT( data_storage::SmartSegment::get_c_most_segment( assembly_->local_segments().at( 19), true )->get_segment_id() == 21 );
		TS_ASSERT( data_storage::SmartSegment::get_c_most_segment( assembly_->local_segments().at( 19), false )->get_segment_id() == 21 );
		TS_ASSERT( data_storage::SmartSegment::get_n_most_segment( assembly_->local_segments().at( 21), true )->get_segment_id() == 19 );
		TS_ASSERT( data_storage::SmartSegment::get_n_most_segment( assembly_->local_segments().at( 21), false )->get_segment_id() == 19 );

		//Check chimaera connectivity
		TS_ASSERT( assembly_->get_forward_assembly() == assembly_->get_reverse_assembly() );

		//Now for a real test
		//See if it's the same as if we had never added in the first place
		data_storage::SmartAssemblyOP test_assembly( new data_storage::SmartAssembly( segment_vector_ ) );
		test_assembly->set_starting_segment( test_assembly->get_segment_vector()->at( 31 ), "all" );
		test_assembly->add_segment( false, 1, 6, 3 );
		TS_ASSERT( assembly_->get_length() == test_assembly->get_length() );

	}


	//PROBLEM TEST??
	void test_double_unchimerize_single_helix_two_triplets_cterm(){
		TS_TRACE( "Begin double_unchimerize_single_helix_two_triplets_cterm" );
		assembly_ = sewing_testing::make_double_chimaera_cterm_first( assembly_ );
		bool nterm = false;
		assembly_->delete_segment( nterm );
		TS_TRACE( assembly_->get_forward_assembly() );
		//Now we're going to test all the parameters of the assembly
		TS_ASSERT( assembly_->get_size() == 3 );
		//Test length
		core::Size assembly_length = 0;
		data_storage::SmartSegmentOP current_seg = assembly_->get_n_terminal_segment();
		while ( current_seg != nullptr ) {
			assembly_length = assembly_length +  current_seg->get_length();
			//assembly_length += current_seg->get_length();
			current_seg = current_seg->get_c_terminal_neighbor();
		}
		TS_TRACE( "Reported assembly length: " + utility::to_string( assembly_->get_length() ) );
		TS_TRACE( "Calculated assembly length: " + utility::to_string( assembly_length ) );
		TS_ASSERT( assembly_->get_length() == assembly_length );
		//Make sure the termini are right
		TS_ASSERT( assembly_->get_n_terminal_segment()->get_segment_id() == 19 );
		TS_ASSERT( assembly_->get_c_terminal_segment()->is_chimaeric() );
		TS_ASSERT( assembly_->get_c_terminal_segment()->get_c_terminal_parent()->get_segment_id() == 31 );
		TS_ASSERT( assembly_->get_c_terminal_segment()->get_n_terminal_parent()->get_segment_id() == 21 );
		//Test all connectivity (TODO)
		//19-20-(21|31)
		TS_ASSERT( data_storage::SmartSegment::get_c_most_segment( assembly_->get_n_terminal_segment(), false )->is_chimaeric() );
		TS_ASSERT( data_storage::SmartSegment::get_c_most_segment( assembly_->get_n_terminal_segment(), true )->is_chimaeric() );
		TS_ASSERT( data_storage::SmartSegment::get_n_most_segment( assembly_->get_c_terminal_segment(), false )->is_chimaeric() );
		TS_ASSERT( data_storage::SmartSegment::get_n_most_segment( assembly_->get_c_terminal_segment(), true )->get_segment_id() == 19 );
		//Now for the middle part
		data_storage::SmartSegmentOP chimaera = assembly_->get_c_terminal_segment();
		TS_ASSERT( chimaera->is_chimaeric() );
		//Make sure parents' connections to neighbors are maintained
		TS_ASSERT( chimaera->get_n_terminal_parent()->get_n_terminal_neighbor()->get_segment_id() == 20 );
		TS_ASSERT( chimaera->get_c_terminal_parent()->get_segment_id() == 31 );
		TS_ASSERT( chimaera->get_c_terminal_parent()->get_c_terminal_neighbor() == nullptr );
		TS_ASSERT( chimaera->get_c_terminal_parent()->get_n_terminal_neighbor() == nullptr );
		//Make sure the deleted part's connectivity was restored (1-3)
		TS_ASSERT( data_storage::SmartSegment::get_c_most_segment( assembly_->local_segments().at( 1), true )->get_segment_id() == 3 );
		TS_ASSERT( data_storage::SmartSegment::get_c_most_segment( assembly_->local_segments().at( 1), false )->get_segment_id() == 3 );
		TS_ASSERT( data_storage::SmartSegment::get_n_most_segment( assembly_->local_segments().at( 3), true )->get_segment_id() == 1 );
		TS_ASSERT( data_storage::SmartSegment::get_n_most_segment( assembly_->local_segments().at( 3), false )->get_segment_id() == 1 );

		//Check chimaera connectivity
		TS_ASSERT( assembly_->get_forward_assembly() == assembly_->get_reverse_assembly() );

		//Now for a real test
		//See if it's the same as if we had never added in the first place
		data_storage::SmartAssemblyOP test_assembly( new data_storage::SmartAssembly( segment_vector_ ) );
		test_assembly->set_starting_segment( test_assembly->get_segment_vector()->at( 31 ), "all" );
		test_assembly->add_segment( true, 21, 2, 4 );
		TS_ASSERT( assembly_->get_length() == test_assembly->get_length() );
	}
	//In theory, if the delete works properly, the switch should work since it's just an add after the delete

	void test_reverting_double_unchimerize_nterm(){
		TS_TRACE( "Begin double_unchimerize_single_helix_two_triplets_nterm" );
		assembly_ = sewing_testing::make_double_chimaera_nterm_first( assembly_ );
		bool nterm = true;
		assembly_->delete_segment( nterm );
		//Make sure the segments are recognized as not in the assembly
		TS_ASSERT( assembly_->local_segments().at( 19 )->is_in_Assembly() == false );
		TS_ASSERT( assembly_->local_segments().at( 20 )->is_in_Assembly() == false );
		TS_ASSERT( assembly_->local_segments().at( 21 )->is_in_Assembly() == false );

		TS_TRACE( assembly_->get_forward_assembly() );
		assembly_->revert();

		//Now run all the tests on it that we ran when we first made it
		//Note the parent hierarchy/basis residues will be switched to as if we had added cterm first
		TS_ASSERT( assembly_->get_size() == 5 );
		//Test length
		core::Size assembly_length = 0;
		data_storage::SmartSegmentOP current_seg = assembly_->get_n_terminal_segment();
		while ( current_seg != nullptr ) {
			assembly_length = assembly_length +  current_seg->get_length();
			//assembly_length += current_seg->get_length();
			current_seg = current_seg->get_c_terminal_neighbor();
		}
		TS_TRACE( "Reported assembly length: " + utility::to_string( assembly_->get_length() ) );
		TS_TRACE( "Calculated assembly length: " + utility::to_string( assembly_length ) );
		TS_ASSERT( assembly_->get_length() == assembly_length );
		//Test all connectivity (TODO)
		//19-20-(21||31|1)-2-3
		TS_ASSERT( data_storage::SmartSegment::get_c_most_segment( assembly_->get_n_terminal_segment(), false )->is_chimaeric() );
		TS_ASSERT( data_storage::SmartSegment::get_c_most_segment( assembly_->get_n_terminal_segment(), true )->get_segment_id() == 3 );
		TS_ASSERT( data_storage::SmartSegment::get_n_most_segment( assembly_->get_c_terminal_segment(), false )->is_chimaeric() );
		TS_ASSERT( data_storage::SmartSegment::get_n_most_segment( assembly_->get_c_terminal_segment(), true )->get_segment_id() == 19 );
		//Now for the middle part
		data_storage::SmartSegmentOP chimaera = assembly_->local_segments().at( 2 )->get_n_terminal_neighbor(); //Double chimaera, actually
		TS_ASSERT( chimaera->is_chimaeric() );
		//Make sure parents' connections to neighbors are maintained
		TS_ASSERT( chimaera->get_n_terminal_parent()->get_n_terminal_neighbor()->get_segment_id() == 20 );
		TS_ASSERT( chimaera->get_c_terminal_parent()->is_chimaeric() );
		TS_ASSERT( chimaera->get_c_terminal_parent()->get_n_terminal_parent()->get_segment_id() == 31 );
		TS_ASSERT( chimaera->get_c_terminal_parent()->get_c_terminal_parent()->get_segment_id() == 1 );
		TS_ASSERT( chimaera->get_c_terminal_parent()->get_c_terminal_parent()->get_c_terminal_neighbor()->get_segment_id() == 2 );
		//Check chimaera connectivity
		TS_ASSERT( assembly_->get_forward_assembly() == assembly_->get_reverse_assembly() );

	}

	void test_reverting_double_unchimerize_cterm(){
		TS_TRACE( "Begin double_unchimerize_single_helix_two_triplets_cterm" );
		assembly_ = sewing_testing::make_double_chimaera_cterm_first( assembly_ );
		bool nterm = false;
		assembly_->delete_segment( nterm );
		TS_TRACE( assembly_->get_forward_assembly() );
		TS_ASSERT( assembly_->local_segments().at( 1 )->is_in_Assembly() == false );
		TS_ASSERT( assembly_->local_segments().at( 2 )->is_in_Assembly() == false );
		TS_ASSERT( assembly_->local_segments().at( 3 )->is_in_Assembly() == false );
		assembly_->revert();
		TS_TRACE( assembly_->get_forward_assembly() );
		TS_ASSERT( assembly_->get_size() == 5 );
		//Test length
		core::Size assembly_length = 0;
		data_storage::SmartSegmentOP current_seg = assembly_->get_n_terminal_segment();
		while ( current_seg != nullptr ) {
			assembly_length = assembly_length +  current_seg->get_length();
			//assembly_length += current_seg->get_length();
			current_seg = current_seg->get_c_terminal_neighbor();
		}
		TS_TRACE( "Reported assembly length: " + utility::to_string( assembly_->get_length() ) );
		TS_TRACE( "Calculated assembly length: " + utility::to_string( assembly_length ) );
		TS_ASSERT( assembly_->get_length() == assembly_length );
		//Test all connectivity (TODO)
		//19-20-(21|31||1)-2-3
		TS_ASSERT( data_storage::SmartSegment::get_c_most_segment( assembly_->get_n_terminal_segment(), false )->is_chimaeric() );
		TS_ASSERT( data_storage::SmartSegment::get_c_most_segment( assembly_->get_n_terminal_segment(), true )->get_segment_id() == 3 );
		TS_ASSERT( data_storage::SmartSegment::get_n_most_segment( assembly_->get_c_terminal_segment(), false )->is_chimaeric() );
		TS_ASSERT( data_storage::SmartSegment::get_n_most_segment( assembly_->get_c_terminal_segment(), true )->get_segment_id() == 19 );
		//Now for the middle part
		data_storage::SmartSegmentOP chimaera = assembly_->local_segments().at( 20 )->get_c_terminal_neighbor(); //Double chimaera, actually
		TS_ASSERT( chimaera->is_chimaeric() );
		//Make sure parents' connections to neighbors are maintained
		TS_ASSERT( chimaera->get_c_terminal_parent()->get_segment_id() == 1 );
		TS_ASSERT( chimaera->get_c_terminal_parent()->get_c_terminal_neighbor()->get_segment_id() == 2 );

		TS_ASSERT( chimaera->get_n_terminal_parent()->is_chimaeric() );

		TS_ASSERT( chimaera->get_n_terminal_parent()->get_n_terminal_parent()->get_segment_id() == 21 );
		TS_ASSERT( chimaera->get_n_terminal_parent()->get_c_terminal_parent()->get_segment_id() == 31 );
		TS_ASSERT( chimaera->get_n_terminal_parent()->get_c_terminal_parent()->get_c_terminal_neighbor() == nullptr );
		TS_ASSERT( chimaera->get_n_terminal_parent()->get_n_terminal_parent()->get_n_terminal_neighbor()->get_segment_id() == 20 );
		//Check chimaera connectivity
		TS_ASSERT( assembly_->get_forward_assembly() == assembly_->get_reverse_assembly() );
	}



	/////////////////////////////////////////////
	///////////pdbsegs testing//////////////////
	///////////////////////////////////////////


	void dont_test_starting_segment_is_vital(){
		std::map< core::Size, data_storage::SmartSegmentOP > pdbsegs;
		//utility::vector1< std::string > resnums;
		std::string pose_segment_starts_string = "";
		std::string pose_segment_ends_string = "";
		std::string pose_segment_dssp = "";
		utility::vector1< data_storage::LigandDescription > ligands;
		utility::vector1< data_storage::LigandDescription > expanded_ligands;
		std::map< core::Size, data_storage::LigandResidueCOP > partner_ligands;
		std::string required_resnums;
		core::select::residue_selector::ResidueSelectorCOP required_selector;
		bool strict_dssp_changes = true;
		pdbsegs.clear();
		hashing::AlignmentFileGeneratorMover::add_pose_segments_to_segment_vector(
			pose_,
			nullptr,
			assembly_->get_segment_vector(),
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
	}



	void dont_test_do_not_delete_starting_segment_simple_nterm(){
		std::map< core::Size, data_storage::SmartSegmentOP > pdbsegs;
		//utility::vector1< std::string > resnums;
		std::string pose_segment_starts_string = "";
		std::string pose_segment_ends_string = "";
		std::string pose_segment_dssp = "";
		utility::vector1< data_storage::LigandDescription > ligands;
		utility::vector1< data_storage::LigandDescription > expanded_ligands;
		std::map< core::Size, data_storage::LigandResidueCOP > partner_ligands;
		std::string required_resnums;
		core::select::residue_selector::ResidueSelectorCOP required_selector;
		bool strict_dssp_changes = true;
		pdbsegs.clear();
		hashing::AlignmentFileGeneratorMover::add_pose_segments_to_segment_vector(
			pose_,
			nullptr,
			assembly_->get_segment_vector(),
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

	}

	void dont_test_do_not_delete_starting_segment_simple_cterm(){
		std::map< core::Size, data_storage::SmartSegmentOP > pdbsegs;
		//utility::vector1< std::string > resnums;
		std::string pose_segment_starts_string = "";
		std::string pose_segment_ends_string = "";
		std::string pose_segment_dssp = "";
		utility::vector1< data_storage::LigandDescription > ligands;
		utility::vector1< data_storage::LigandDescription > expanded_ligands;
		std::map< core::Size, data_storage::LigandResidueCOP > partner_ligands;
		std::string required_resnums;
		core::select::residue_selector::ResidueSelectorCOP required_selector;
		bool strict_dssp_changes = true;
		pdbsegs.clear();
		hashing::AlignmentFileGeneratorMover::add_pose_segments_to_segment_vector(
			pose_,
			nullptr,
			assembly_->get_segment_vector(),
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


	}


	void dont_test_can_chimerize_chimaera(){
		std::map< core::Size, data_storage::SmartSegmentOP > pdbsegs;
		//utility::vector1< std::string > resnums;
		std::string pose_segment_starts_string = "";
		std::string pose_segment_ends_string = "";
		std::string pose_segment_dssp = "";
		utility::vector1< data_storage::LigandDescription > ligands;
		utility::vector1< data_storage::LigandDescription > expanded_ligands;
		std::map< core::Size, data_storage::LigandResidueCOP > partner_ligands;
		std::string required_resnums;
		core::select::residue_selector::ResidueSelectorCOP required_selector;
		bool strict_dssp_changes = true;
		pdbsegs.clear();
		core::Size num_added = hashing::AlignmentFileGeneratorMover::add_pose_segments_to_segment_vector(
			pose_,
			nullptr,
			assembly_->get_segment_vector(),
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

		TS_ASSERT_EQUALS( num_added, 1 );
		assembly_->pdb_segments( pdbsegs );
		data_storage::SmartSegmentOP starting_segment = assembly_->pdb_segments().at( 1 );
		assembly_->set_starting_segment( starting_segment, "all" );
		TS_ASSERT_EQUALS( assembly_->get_n_terminal_segment()->get_segment_id(), starting_segment->get_segment_id() );
		TS_ASSERT_EQUALS( assembly_->get_c_terminal_segment()->get_segment_id(), starting_segment->get_segment_id() );
	}



	void dont_test_do_not_delete_starting_segment_double_chimaera_nterm(){
		std::map< core::Size, data_storage::SmartSegmentOP > pdbsegs;
		//utility::vector1< std::string > resnums;
		std::string pose_segment_starts_string = "";
		std::string pose_segment_ends_string = "";
		std::string pose_segment_dssp = "";
		utility::vector1< data_storage::LigandDescription > ligands;
		utility::vector1< data_storage::LigandDescription > expanded_ligands;
		std::map< core::Size, data_storage::LigandResidueCOP > partner_ligands;
		std::string required_resnums;
		core::select::residue_selector::ResidueSelectorCOP required_selector;
		bool strict_dssp_changes = true;
		pdbsegs.clear();
		hashing::AlignmentFileGeneratorMover::add_pose_segments_to_segment_vector(
			pose_,
			nullptr,
			assembly_->get_segment_vector(),
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
	}


	void dont_test_do_not_delete_starting_segment_double_chimaera_cterm(){
		std::map< core::Size, data_storage::SmartSegmentOP > pdbsegs;
		//utility::vector1< std::string > resnums;
		std::string pose_segment_starts_string = "";
		std::string pose_segment_ends_string = "";
		std::string pose_segment_dssp = "";
		utility::vector1< data_storage::LigandDescription > ligands;
		utility::vector1< data_storage::LigandDescription > expanded_ligands;
		std::map< core::Size, data_storage::LigandResidueCOP > partner_ligands;
		std::string required_resnums;
		core::select::residue_selector::ResidueSelectorCOP required_selector;
		bool strict_dssp_changes = true;
		pdbsegs.clear();
		hashing::AlignmentFileGeneratorMover::add_pose_segments_to_segment_vector(
			pose_,
			nullptr,
			assembly_->get_segment_vector(),
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

	}

	////////////////////////////////////////////////////////////////////////////
	////////////////// ASSEMBLY RECOVERY TESTS /////////////////////////////////
	////////////////////////////////////////////////////////////////////////////


	void test_assemblies_are_continuous(){
		TS_TRACE( "Begin assemblies_are_continuous" );
		assembly_ = sewing_testing::make_double_chimaera_nterm_first( assembly_ );
		TS_ASSERT( assembly_->is_continuous() );
		assembly_ = data_storage::SmartAssemblyOP( new data_storage::SmartAssembly( segment_vector_ ) );
		assembly_ = sewing_testing::make_double_chimaera_cterm_first( assembly_ );
		TS_ASSERT( assembly_->is_continuous() );
		assembly_ = data_storage::SmartAssemblyOP( new data_storage::SmartAssembly( segment_vector_ ) );
		assembly_ = sewing_testing::create_simple_assembly( assembly_ );
		TS_ASSERT( assembly_->is_continuous() );
	}

	//////////BEGIN LIGAND-RELATED TESTS/////////////////
	void test_auto_detect_contacts(){
		assembly_ = sewing_testing::initial_assembly_with_auto_detected_ligand_contacts( assembly_ );
		//Check that the owner segment of the ligand is correct--segID should be 32
		TS_TRACE( "Ligand owner segment is " + utility::to_string( assembly_->get_local_ligands().at( 1 )->get_owner_segment()->get_segment_id() ) );
		TS_ASSERT( assembly_->get_local_ligands().at( 1 )->get_owner_segment()->get_segment_id() == 33 );
		TS_ASSERT( assembly_->get_local_ligands().at( 1 )->get_owner_segment()->get_const_owned_ligand_residues().size() > 0);
		TS_ASSERT( assembly_->get_local_ligands().at( 1 )->get_owner_segment()->get_const_owned_ligand_residues().at( 1 ) ==  assembly_->get_local_ligands().at( 1 ) );
		//Check that the contacts of the ligand are correct
		TS_ASSERT( assembly_->get_local_ligands().count( 1 ) != 0 );
		utility::vector1< data_storage::LigandContactOP > const & contacts = assembly_->get_local_ligands().at( 1 )->get_current_contacts();
		//All will have ligand_atom 1, segment_id 32, residue_atom 10
		//residue_number will be 10 for 1 and 12 for the other
		for ( data_storage::LigandContactOP contact: contacts ) {
			TS_TRACE( "Ligand contact is " + utility::to_string( contact->ligand_atom ) + ": " + utility::to_string( contact->segment_id ) + "-" + utility::to_string( contact->residue_number ) + "-" + utility::to_string( contact->residue_atom ) );
			TS_ASSERT( contact->ligand_atom == 1 );
			TS_ASSERT( contact->residue_atom == 10 );
			TS_ASSERT( contact->segment_id == 33 );
			TS_ASSERT( contact->residue_number == 10 || contact->residue_number == 14 );
		}
	}

	void test_maintain_contact_cterm_chimerization(){
		assembly_ = sewing_testing::initial_assembly_with_manually_detected_ligand_contacts( assembly_ );
		//Add to c terminus
		bool nterm = false;
		//Add segment 21 at residue 4
		//The basis residue on segment 32 (the C-terminal helix of the assembly) cannot be <= 14
		//Segment 32 is 25 residues long
		//Now do the add for real at basis residue 15 (this would be segment 13, which is 14 residues long)
		TS_ASSERT( assembly_->add_segment( nterm, 13, 15, 7 ) );
		core::Size chimaera_seg_id = assembly_->get_last_chimaera()->get_segment_id();
		TS_TRACE( "Chimaera is segment " + utility::to_string( chimaera_seg_id ) );
		//Contacts should be exactly the same except for the segID, which should be that of the new chimaera
		TS_ASSERT( assembly_->get_local_ligands().count( 1 ) != 0 );
		utility::vector1< data_storage::LigandContactOP > const & contacts = assembly_->get_local_ligands().at( 1 )->get_current_contacts();
		for ( data_storage::LigandContactOP contact: contacts ) {
			TS_TRACE( "Ligand contact is " + utility::to_string( contact->ligand_atom ) + ": " + utility::to_string( contact->segment_id ) + "-" + utility::to_string( contact->residue_number ) + "-" + utility::to_string( contact->residue_atom ) );
			TS_ASSERT( contact->ligand_atom == 1 );
			TS_ASSERT( contact->residue_atom == 10 );
			TS_ASSERT( contact->segment_id == chimaera_seg_id );
			TS_ASSERT( contact->residue_number == 10 || contact->residue_number == 14 );
		}
		//The new chimaera should know that it owns the ligand
		TS_ASSERT( std::dynamic_pointer_cast< data_storage::LigandSegment const >( assembly_->get_last_chimaera() ) != nullptr );
		TS_ASSERT( std::dynamic_pointer_cast< data_storage::LigandSegment const >( assembly_->get_last_chimaera() )->get_const_owned_ligand_residues().size() > 0 );
		//The ligand should know that its owner is the new chimaera
		TS_TRACE( "Ligand owner segment is " + utility::to_string( assembly_->get_local_ligands().at( 1 )->get_owner_segment()->get_segment_id() ) );
		TS_ASSERT( assembly_->get_local_ligands().at( 1 )->get_owner_segment()->get_segment_id() == chimaera_seg_id );
		//The old segment should still think that it owns the ligand
		TS_ASSERT( std::dynamic_pointer_cast< data_storage::LigandSegment const >( assembly_->local_segments().at( 33 ) )->get_const_owned_ligand_residues().size() > 0 );
	}
	void test_maintain_contact_nterm_chimerization(){
		assembly_ = sewing_testing::initial_assembly_with_manually_detected_ligand_contacts( assembly_ );
		//Add to n terminus
		bool nterm = true;
		//Add segment 21 at residue 4
		//The basis residue on segment 32 (the assembly) cannot be >= 10
		//Segment 32 is 25 residues long
		//Segment 21 is 7 residues long
		//Now do the add for real at basis residue 15
		TS_ASSERT(  assembly_->add_segment( nterm, 21, 5, 4 ) );
		//Residue numbering will change due to n terminal chimerization
		//Since seg1's basis residue was 1 and seg2's was 4, all resnums decrease by 1
		core::Size chimaera_seg_id = assembly_->get_last_chimaera()->get_segment_id();
		TS_TRACE( "Chimaera is segment " + utility::to_string( chimaera_seg_id ) );
		//Contacts should have the new segID, the same atoms, but their residue numbers will have changed--NEW = OLD + N - C
		//10 + 4 - 5 = 9; 14 + 4 - 5 = 13
		TS_ASSERT( assembly_->get_local_ligands().count( 1 ) != 0 );
		utility::vector1< data_storage::LigandContactOP > const & contacts = assembly_->get_local_ligands().at( 1 )->get_current_contacts();
		for ( data_storage::LigandContactOP contact: contacts ) {
			TS_TRACE( "Ligand contact is " + utility::to_string( contact->ligand_atom ) + ": " + utility::to_string( contact->segment_id ) + "-" + utility::to_string( contact->residue_number ) + "-" + utility::to_string( contact->residue_atom ) );
			TS_ASSERT( contact->ligand_atom == 1 );
			TS_ASSERT( contact->residue_atom == 10 );
			TS_ASSERT( contact->segment_id == chimaera_seg_id );
			TS_ASSERT( contact->residue_number == 9 || contact->residue_number == 13 );
		}
		//The new chimaera should know that it owns the ligand
		TS_ASSERT( std::dynamic_pointer_cast< data_storage::LigandSegment const >( assembly_->get_last_chimaera() ) != nullptr );
		TS_ASSERT( std::dynamic_pointer_cast< data_storage::LigandSegment const >( assembly_->get_last_chimaera() )->get_const_owned_ligand_residues().size() > 0 );
		//The ligand should know that its owner is the new chimaera
		TS_TRACE( "Ligand owner segment is " + utility::to_string( assembly_->get_local_ligands().at( 1 )->get_owner_segment()->get_segment_id() ) );
		TS_ASSERT( assembly_->get_local_ligands().at( 1 )->get_owner_segment()->get_segment_id() == chimaera_seg_id );
		//The old segment should still think that it owns the ligand
		TS_ASSERT( std::dynamic_pointer_cast< data_storage::LigandSegment const >( assembly_->local_segments().at( 33 ) )->get_const_owned_ligand_residues().size() > 0 );
	}
	void test_maintain_contact_double_chimerization_cterm_first(){
		assembly_ = sewing_testing::initial_assembly_with_manually_detected_ligand_contacts( assembly_ );
		//Add to c terminus
		bool nterm = false;
		//Add segment 13 at residue 7
		//Segment 13 is 14 residues long
		TS_ASSERT( assembly_->add_segment( nterm, 13, 15, 7 ) );
		//No resnum changes yet
		//Save chimaera
		data_storage::SmartSegmentOP initial_chimaera = assembly_->get_last_chimaera();
		//Now add to n terminus
		nterm = true;
		//Add segment 21 at residue 4
		//Segment 32 is 25 residues long, 21 is 7 residues long
		TS_ASSERT(  assembly_->add_segment( nterm, 21, 5, 4 ) );
		//All residue numbers are now decreased by one (contacts are now 9 and 13)
		data_storage::SmartSegmentOP double_chimaera = assembly_->get_last_chimaera();
		//All of these segments should think that they own the ligand
		TS_ASSERT( std::dynamic_pointer_cast< data_storage::LigandSegment const >( assembly_->local_segments().at( 33 ) )->get_const_owned_ligand_residues().size() > 0 );
		TS_ASSERT( std::dynamic_pointer_cast< data_storage::LigandSegment const >( initial_chimaera )->get_const_owned_ligand_residues().size() > 0 );
		TS_ASSERT( std::dynamic_pointer_cast< data_storage::LigandSegment const >( double_chimaera )->get_const_owned_ligand_residues().size() > 0 );
		//BUT the ligand should know that it's really the double chimaera who owns it
		TS_TRACE( "Ligand owner segment is " + utility::to_string( assembly_->get_local_ligands().at( 1 )->get_owner_segment()->get_segment_id() ) );
		TS_TRACE( "Double chimaera is segment " + utility::to_string( double_chimaera->get_segment_id() ) );
		TS_TRACE( "Initial chimaera is segment " + utility::to_string( initial_chimaera->get_segment_id() ) );
		TS_ASSERT( assembly_->get_local_ligands().at( 1 )->get_owner_segment()->get_segment_id() == double_chimaera->get_segment_id() );
		//The ligand contacts should reflect that the contacts are in the double chimaera and have the correct residue/atom numbers
		TS_ASSERT( assembly_->get_local_ligands().count( 1 ) != 0 );
		utility::vector1< data_storage::LigandContactOP > const & contacts = assembly_->get_local_ligands().at( 1 )->get_current_contacts();
		for ( data_storage::LigandContactOP contact: contacts ) {
			TS_TRACE( "Ligand contact is " + utility::to_string( contact->ligand_atom ) + ": " + utility::to_string( contact->segment_id ) + "-" + utility::to_string( contact->residue_number ) + "-" + utility::to_string( contact->residue_atom ) );
			TS_ASSERT( contact->ligand_atom == 1 );
			TS_ASSERT( contact->residue_atom == 10 );
			TS_ASSERT( contact->segment_id == double_chimaera->get_segment_id() );
			TS_ASSERT( contact->residue_number == 9 || contact->residue_number == 13 );
		}
	}
	void test_maintain_contact_double_chimerization_nterm_first(){
		assembly_ = sewing_testing::initial_assembly_with_manually_detected_ligand_contacts( assembly_ );
		//Add to n terminus
		bool nterm = true;
		//Add segment 21 at residue 4
		//Segment 32 is 25 residues long, 21 is 7 residues long
		TS_ASSERT(  assembly_->add_segment( nterm, 21, 5, 4 ) );
		//Save chimaera
		//All residue numbers are now decreased by one
		data_storage::SmartSegmentOP initial_chimaera = assembly_->get_last_chimaera();
		//Now add to c terminus
		nterm = false;
		//Add segment 13 at residue 7
		//Segment 13 is 14 residues long
		TS_ASSERT( assembly_->add_segment( nterm, 13, 15, 7 ) );
		data_storage::SmartSegmentOP double_chimaera = assembly_->get_last_chimaera();
		//All the same tests from cterm first apply
		//All of these segments should think that they own the ligand
		TS_ASSERT( std::dynamic_pointer_cast< data_storage::LigandSegment const >( assembly_->local_segments().at( 33 ) )->get_const_owned_ligand_residues().size() > 0 );
		TS_ASSERT( std::dynamic_pointer_cast< data_storage::LigandSegment const >( initial_chimaera )->get_const_owned_ligand_residues().size() > 0 );
		TS_ASSERT( std::dynamic_pointer_cast< data_storage::LigandSegment const >( double_chimaera )->get_const_owned_ligand_residues().size() > 0 );
		//BUT the ligand should know that it's really the double chimaera who owns it
		//TS_TRACE( "Ligand owner segment is " + utility::to_string( assembly_->get_local_ligands().at( 1 )->get_owner_segment()->get_segment_id() ) );
		TS_TRACE( "Ligand owner segment is " + utility::to_string( assembly_->get_local_ligands().at( 1 )->get_owner_segment()->get_segment_id() ) );
		TS_TRACE( "Double chimaera is segment " + utility::to_string( double_chimaera->get_segment_id() ) );
		TS_TRACE( "Initial chimaera is segment " + utility::to_string( initial_chimaera->get_segment_id() ) );

		TS_ASSERT( assembly_->get_local_ligands().at( 1 )->get_owner_segment()->get_segment_id() == double_chimaera->get_segment_id() );
		//The ligand contacts should reflect that the contacts are in the double chimaera and have the correct residue/atom numbers
		TS_ASSERT( assembly_->get_local_ligands().count( 1 ) != 0 );
		utility::vector1< data_storage::LigandContactOP > const & contacts = assembly_->get_local_ligands().at( 1 )->get_current_contacts();
		for ( data_storage::LigandContactOP contact: contacts ) {
			TS_TRACE( "Ligand contact is " + utility::to_string( contact->ligand_atom ) + ": " + utility::to_string( contact->segment_id ) + "-" + utility::to_string( contact->residue_number ) + "-" + utility::to_string( contact->residue_atom ) );
			TS_ASSERT( contact->ligand_atom == 1 );
			TS_ASSERT( contact->residue_atom == 10 );
			TS_ASSERT( contact->segment_id == double_chimaera->get_segment_id() );
			TS_ASSERT( contact->residue_number == 9 || contact->residue_number == 13 );
		}
	}
	void test_maintain_contact_unchimerize_cterm(){
		//This test is just to ensure that the assembly is the same before and after an add+delete
		assembly_ = sewing_testing::initial_assembly_with_manually_detected_ligand_contacts( assembly_ );
		//Add to c terminus
		bool nterm = false;
		TS_ASSERT( assembly_->add_segment( nterm, 13, 15, 7 ) );
		data_storage::SmartSegmentOP chimaera = assembly_->get_last_chimaera();
		assembly_->delete_segment( nterm );
		//The assembly should be as it was before any add/delete
		TS_TRACE( "Ligand owner segment is " + utility::to_string( assembly_->get_local_ligands().at( 1 )->get_owner_segment()->get_segment_id() ) );
		//TS_TRACE( "Double chimaera is segment " + utility::to_string( double_chimaera->get_segment_id() ) );
		//TS_TRACE( "Initial chimaera is segment " + utility::to_string( initial_chimaera->get_segment_id() ) );

		TS_ASSERT( assembly_->get_local_ligands().at( 1 )->get_owner_segment()->get_segment_id() == 33 );
		TS_ASSERT( assembly_->get_local_ligands().at( 1 )->get_owner_segment()->get_const_owned_ligand_residues().size() > 0);
		TS_ASSERT( assembly_->get_local_ligands().at( 1 )->get_owner_segment()->get_const_owned_ligand_residues().at( 1 ) ==  assembly_->get_local_ligands().at( 1 ) );
		//Check that the contacts of the ligand are correct
		TS_ASSERT( assembly_->get_local_ligands().count( 1 ) != 0 );
		utility::vector1< data_storage::LigandContactOP > const & contacts = assembly_->get_local_ligands().at( 1 )->get_current_contacts();
		//All will have ligand_atom 1, segment_id 32, residue_atom 10
		//residue_number will be 10 for 1 and 12 for the other
		for ( data_storage::LigandContactOP contact: contacts ) {
			TS_TRACE( "Ligand contact is " + utility::to_string( contact->ligand_atom ) + ": " + utility::to_string( contact->segment_id ) + "-" + utility::to_string( contact->residue_number ) + "-" + utility::to_string( contact->residue_atom ) );
			TS_ASSERT( contact->ligand_atom == 1 );
			TS_ASSERT( contact->residue_atom == 10 );
			TS_ASSERT( contact->segment_id == 33 );
			TS_ASSERT( contact->residue_number == 10 || contact->residue_number == 14 );
		}
		//BUT the chimaera should remember owning the ligand
		TS_ASSERT( std::dynamic_pointer_cast< data_storage::LigandSegment const >( chimaera )->get_const_owned_ligand_residues().size() > 0 );
	}
	void test_maintain_contact_unchimerize_nterm(){
		//This test is just to ensure that the assembly is the same before and after an add+delete
		assembly_ = sewing_testing::initial_assembly_with_manually_detected_ligand_contacts( assembly_ );
		//Add to c terminus
		bool nterm = true;
		TS_ASSERT(  assembly_->add_segment( nterm, 21, 5, 4 ) );
		data_storage::SmartSegmentOP chimaera = assembly_->get_last_chimaera();
		assembly_->delete_segment( nterm );
		//The assembly should be as it was before any add/delete
		TS_TRACE( "Ligand owner segment is " + utility::to_string( assembly_->get_local_ligands().at( 1 )->get_owner_segment()->get_segment_id() ) );
		//TS_TRACE( "Double chimaera is segment " + utility::to_string( double_chimaera->get_segment_id() ) );
		//TS_TRACE( "Initial chimaera is segment " + utility::to_string( initial_chimaera->get_segment_id() ) );

		TS_ASSERT( assembly_->get_local_ligands().at( 1 )->get_owner_segment()->get_segment_id() == 33 );
		TS_ASSERT( assembly_->get_local_ligands().at( 1 )->get_owner_segment()->get_const_owned_ligand_residues().size() > 0);
		TS_ASSERT( assembly_->get_local_ligands().at( 1 )->get_owner_segment()->get_const_owned_ligand_residues().at( 1 ) ==  assembly_->get_local_ligands().at( 1 ) );
		//Check that the contacts of the ligand are correct
		TS_ASSERT( assembly_->get_local_ligands().count( 1 ) != 0 );
		utility::vector1< data_storage::LigandContactOP > const & contacts = assembly_->get_local_ligands().at( 1 )->get_current_contacts();
		//All will have ligand_atom 1, segment_id 32, residue_atom 10
		//residue_number will be 10 for 1 and 12 for the other
		for ( data_storage::LigandContactOP contact: contacts ) {
			TS_TRACE( "Ligand contact is " + utility::to_string( contact->ligand_atom ) + ": " + utility::to_string( contact->segment_id ) + "-" + utility::to_string( contact->residue_number ) + "-" + utility::to_string( contact->residue_atom ) );
			TS_ASSERT( contact->ligand_atom == 1 );
			TS_ASSERT( contact->residue_atom == 10 );
			TS_ASSERT( contact->segment_id == 33 );
			TS_ASSERT( contact->residue_number == 10 || contact->residue_number == 14 );
		}
		//BUT the chimaera should remember owning the ligand
		TS_ASSERT( std::dynamic_pointer_cast< data_storage::LigandSegment const >( chimaera )->get_const_owned_ligand_residues().size() > 0 );
	}
	void test_maintain_contact_recurse_revert_cterm(){
		//Make the double chimaera as before
		assembly_ = sewing_testing::initial_assembly_with_manually_detected_ligand_contacts( assembly_ );
		bool nterm = false;
		TS_ASSERT( assembly_->add_segment( nterm, 13, 15, 7 ) );
		data_storage::SmartSegmentOP initial_chimaera = assembly_->get_last_chimaera();
		nterm = true;
		TS_ASSERT(  assembly_->add_segment( nterm, 21, 5, 4 ) );
		data_storage::SmartSegmentOP double_chimaera = assembly_->get_last_chimaera();
		//Now we're going to delete the original c-terminal addition
		assembly_->delete_segment( !nterm );
		//Since we deleted from the c term, all resnums are still decreased by 1 relative to the original (contacts are 9 and 13)
		data_storage::SmartSegmentOP final_chimaera = assembly_->get_last_chimaera();
		//The assembly should behave test_wise just like if we had only ever added to the n terminus
		TS_ASSERT( assembly_->get_local_ligands().count( 1 ) != 0 );
		utility::vector1< data_storage::LigandContactOP > const & contacts = assembly_->get_local_ligands().at( 1 )->get_current_contacts();


		TS_TRACE( "Ligand owner segment is " + utility::to_string( assembly_->get_local_ligands().at( 1 )->get_owner_segment()->get_segment_id() ) );
		TS_TRACE( "Double chimaera is segment " + utility::to_string( double_chimaera->get_segment_id() ) );
		TS_TRACE( "Initial chimaera is segment " + utility::to_string( initial_chimaera->get_segment_id() ) );
		TS_TRACE( "Final chimaera is segment " + utility::to_string( final_chimaera->get_segment_id() ) );

		for ( data_storage::LigandContactOP contact: contacts ) {
			TS_TRACE( "Ligand contact is " + utility::to_string( contact->ligand_atom ) + ": " + utility::to_string( contact->segment_id ) + "-" + utility::to_string( contact->residue_number ) + "-" + utility::to_string( contact->residue_atom ) );
			TS_ASSERT( contact->ligand_atom == 1 );
			TS_ASSERT( contact->residue_atom == 10 );
			TS_ASSERT( contact->segment_id == final_chimaera->get_segment_id() );
			TS_ASSERT( contact->residue_number == 9 || contact->residue_number == 13 );
		}
		//The new chimaera should know that it owns the ligand
		TS_ASSERT( std::dynamic_pointer_cast< data_storage::LigandSegment const >( final_chimaera )->get_const_owned_ligand_residues().size() > 0 );
		//The ligand should know that its owner is the new chimaera
		TS_ASSERT( assembly_->get_local_ligands().at( 1 )->get_owner_segment()->get_segment_id() == final_chimaera->get_segment_id() );
		//The old segment should still think that it owns the ligand
		TS_ASSERT( std::dynamic_pointer_cast< data_storage::LigandSegment const >( assembly_->local_segments().at( 33 ) )->get_const_owned_ligand_residues().size() > 0 );
		//All of the chimaerae will remember owning the ligand
		TS_ASSERT( std::dynamic_pointer_cast< data_storage::LigandSegment const >( initial_chimaera )->get_const_owned_ligand_residues().size() > 0 );
		TS_ASSERT( std::dynamic_pointer_cast< data_storage::LigandSegment const >( double_chimaera )->get_const_owned_ligand_residues().size() > 0 );
	}
	void test_maintain_contact_recurse_revert_nterm(){
		//Make the double chimaera as before
		assembly_ = sewing_testing::initial_assembly_with_manually_detected_ligand_contacts( assembly_ );
		bool nterm = true;
		TS_ASSERT(  assembly_->add_segment( nterm, 21, 5, 4 ) );
		data_storage::SmartSegmentOP initial_chimaera = assembly_->get_last_chimaera();
		nterm = false;
		TS_ASSERT( assembly_->add_segment( nterm, 13, 15, 7 ) );
		data_storage::SmartSegmentOP double_chimaera = assembly_->get_last_chimaera();
		//Now we're going to delete the original c-terminal addition
		assembly_->delete_segment( !nterm );
		data_storage::SmartSegmentOP final_chimaera = assembly_->get_last_chimaera();
		//The assembly should behave test_wise just like if we had only ever added to the c terminus
		TS_ASSERT( assembly_->get_local_ligands().count( 1 ) != 0 );
		utility::vector1< data_storage::LigandContactOP > const & contacts = assembly_->get_local_ligands().at( 1 )->get_current_contacts();
		TS_TRACE( "Ligand owner segment is " + utility::to_string( assembly_->get_local_ligands().at( 1 )->get_owner_segment()->get_segment_id() ) );
		TS_TRACE( "Double chimaera is segment " + utility::to_string( double_chimaera->get_segment_id() ) );
		TS_TRACE( "Initial chimaera is segment " + utility::to_string( initial_chimaera->get_segment_id() ) );
		TS_TRACE( "Final chimaera is segment " + utility::to_string( final_chimaera->get_segment_id() ) );
		for ( data_storage::LigandContactOP contact: contacts ) {
			TS_TRACE( "Ligand contact is " + utility::to_string( contact->ligand_atom ) + ": " + utility::to_string( contact->segment_id ) + "-" + utility::to_string( contact->residue_number ) + "-" + utility::to_string( contact->residue_atom ) );
			TS_ASSERT( contact->ligand_atom == 1 );
			TS_ASSERT( contact->residue_atom == 10 );
			TS_ASSERT( contact->segment_id == final_chimaera->get_segment_id() );
			TS_ASSERT( contact->residue_number == 10 || contact->residue_number == 14 );
		}
		//The new chimaera should know that it owns the ligand
		TS_ASSERT( std::dynamic_pointer_cast< data_storage::LigandSegment const >( final_chimaera )->get_const_owned_ligand_residues().size() > 0 );
		//The ligand should know that its owner is the new chimaera
		TS_ASSERT( assembly_->get_local_ligands().at( 1 )->get_owner_segment()->get_segment_id() == final_chimaera->get_segment_id() );
		//The old segment should still think that it owns the ligand
		TS_ASSERT( std::dynamic_pointer_cast< data_storage::LigandSegment const >( assembly_->local_segments().at( 33 ) )->get_const_owned_ligand_residues().size() > 0 );
		//All of the chimaerae will remember owning the ligand
		TS_ASSERT( std::dynamic_pointer_cast< data_storage::LigandSegment const >( initial_chimaera )->get_const_owned_ligand_residues().size() > 0 );
		TS_ASSERT( std::dynamic_pointer_cast< data_storage::LigandSegment const >( double_chimaera )->get_const_owned_ligand_residues().size() > 0 );
	}

};



