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
#include <protocols/sewing/data_storage/SmartSegment.hh>
#include <protocols/sewing/data_storage/Basis.hh>
#include <protocols/sewing/data_storage/SmartSewingResidue.hh>
#include <protocols/sewing/hashing/hasher_data.hh>
#include <protocols/sewing/hashing/Hasher.hh>
#include <protocols/sewing/hashing/ModelFileReader.hh>
#include <protocols/sewing/hashing/AlignmentGenerator.hh>
#include <protocols/init_util.hh>
#include <core/conformation/Atom.hh>
#include <numeric/xyzVector.hh>

// Protocol Headers
#include <basic/Tracer.hh>
using namespace protocols::sewing;

static basic::Tracer TR("HasherTests");


class NewSEWINGHasherTests : public CxxTest::TestSuite {
	//Define Variables
private:
	hashing::SegmentVectorOP segment_vector_;
	hashing::HasherSettings hasher_settings_;
	hashing::HasherOP hasher_;
public:

	void setUp(){
		core_init_with_additional_options( "-in:auto_setup_metals" );
		//Create segment vector
		segment_vector_ = hashing::ModelFileReader::read_model_file( "protocols/sewing/inputs/test.segments" );
		//Create hasher settings--no clashes, 4 atom( 1 residue) overlap, 3 boxes per dimension, hash between termini
		hasher_settings_ = hashing::HasherSettings( 0, 4, 3, true );
		hasher_ = hashing::HasherOP( new hashing::Hasher( hasher_settings_, segment_vector_ ) );
	}

	void tearDown(){
	}

	void test_alignment_generator(){
		//Create alignment generator
		//Only counts alignments that meet minimum hash score--since we're just testing total possible, disable that feature (at least 3 atoms will *always* overlap due to basis)
		hasher_settings_.min_hash_score = 0;
		hashing::AlignmentGenerator ag( hasher_settings_, segment_vector_ );
		//Only function to test is get_all_alignments--vector of all possible basis pairs for an edge given 2 segments involved
		//4 and 21 = 9 res  & 7 res
		utility::vector1< std::pair< core::Size, core::Size > > alignments_4_21 = ag.get_all_alignments( segment_vector_->at( 4 ), segment_vector_->at( 21 ) );
		utility::vector1< std::pair< core::Size, core::Size > > alignments_21_4 = ag.get_all_alignments( segment_vector_->at( 21 ), segment_vector_->at( 4 ) );
		TS_ASSERT_EQUALS( alignments_4_21.size(), 45); // segment_vector_->at(4)->get_length()*segment_vector_->at(21)->get_length());
		TS_ASSERT_EQUALS( alignments_21_4.size(),  49);// segment_vector_->at(4)->get_length()*segment_vector_->at(21)->get_length());
	}

	//Easy to test--true for opposite terminal segments, false for same termini or any non-terminal segments
	void test_can_hash(){

		TS_ASSERT( hasher_->can_hash( segment_vector_->at( 4 ), segment_vector_->at( 30) ) );
		TS_ASSERT( hasher_->can_hash( segment_vector_->at( 4 ), segment_vector_->at( 21) ) );
		TS_ASSERT( hasher_->can_hash( segment_vector_->at( 1 ), segment_vector_->at( 30) ) );
		TS_ASSERT( hasher_->can_hash( segment_vector_->at( 1 ), segment_vector_->at( 21) ) );
		TS_ASSERT( !hasher_->can_hash( segment_vector_->at( 4 ), segment_vector_->at( 1) ) );
		TS_ASSERT( !hasher_->can_hash( segment_vector_->at( 1 ), segment_vector_->at( 4) ) );
		TS_ASSERT( !hasher_->can_hash( segment_vector_->at( 21 ), segment_vector_->at( 30) ) );
		TS_ASSERT( !hasher_->can_hash( segment_vector_->at( 30 ), segment_vector_->at( 21) ) );

		TS_ASSERT( !hasher_->can_hash( segment_vector_->at( 4 ), segment_vector_->at( 11 ) ) );
		TS_ASSERT( !hasher_->can_hash( segment_vector_->at( 6 ), segment_vector_->at( 11 ) ) );

		TS_ASSERT( hasher_->can_hash( segment_vector_->at( 31 ), segment_vector_->at( 21) ) );
		TS_ASSERT( hasher_->can_hash( segment_vector_->at( 31 ), segment_vector_->at( 1 ) ) );
		TS_ASSERT( hasher_->can_hash( segment_vector_->at( 1 ), segment_vector_->at( 31) ) );
		TS_ASSERT( hasher_->can_hash( segment_vector_->at( 21 ), segment_vector_->at( 31 ) ) );
	}
	//iterate_over_basis_pairs (returns pair of iterators) was covered in the AlignmentGenerator test


	//transform_segments
	void test_transform_segments(){
		//Third atom of basis residue will be at the origin after transform
		TS_ASSERT_THROWS_NOTHING( hasher_->transform_segments( data_storage::Basis( 31, 5 ), true ) );
		TS_ASSERT_DELTA( hasher_->get_seg1()->get_const_residue_vector().at( 5 )->get_const_atom_vector().at( 3 ).xyz()[ 0 ], 0.0, 1e-5 );
		TS_ASSERT_DELTA( hasher_->get_seg1()->get_const_residue_vector().at( 5 )->get_const_atom_vector().at( 3 ).xyz()[ 1 ], 0.0, 1e-5 );
		TS_ASSERT_DELTA( hasher_->get_seg1()->get_const_residue_vector().at( 5 )->get_const_atom_vector().at( 3 ).xyz()[ 2 ], 0.0, 1e-5 );

		TS_ASSERT_THROWS_NOTHING( hasher_->transform_segments( data_storage::Basis( 31, 5 ), false ) );
		TS_ASSERT_DELTA( hasher_->get_seg2()->get_const_residue_vector().at( 5 )->get_const_atom_vector().at( 3 ).xyz()[ 0 ], 0.0, 1e-5 );
		TS_ASSERT_DELTA( hasher_->get_seg2()->get_const_residue_vector().at( 5 )->get_const_atom_vector().at( 3 ).xyz()[ 1 ], 0.0, 1e-5 );
		TS_ASSERT_DELTA( hasher_->get_seg2()->get_const_residue_vector().at( 5 )->get_const_atom_vector().at( 3 ).xyz()[ 2 ], 0.0, 1e-5 );
	}

	//hash_segments
	//Look at contents (size) of hashmap before and after function call, make sure everything is being inserted that is supposed to
	void test_hash_segments(){
		//First hash segment 31, length 11--basis residue not hashed
		TS_ASSERT_THROWS_NOTHING( hasher_->hash_segments( segment_vector_->at( 31 )->clone(), data_storage::Basis( 31, 5 ) ) );
		TS_ASSERT_EQUALS( hasher_->get_hash_map_size(),4*segment_vector_->at(31)->get_length() - 8 );
		//Clear hashmap
		hasher_ = hashing::HasherOP( new hashing::Hasher( hasher_settings_, segment_vector_ ) );
		//Then hash segments 4-6: seg4 length 9, seg6 length 13, seg5 length 13
		TS_ASSERT_THROWS_NOTHING( hasher_->hash_segments( segment_vector_->at( 4 )->clone(), data_storage::Basis( 4, 5 ) ) );
		TS_ASSERT_THROWS_NOTHING( hasher_->hash_segments( segment_vector_->at( 5 )->clone(), data_storage::Basis( 5, 5 ) ) );
		TS_ASSERT_THROWS_NOTHING( hasher_->hash_segments( segment_vector_->at( 6 )->clone(), data_storage::Basis( 6, 5 ) ) );
		TS_ASSERT_EQUALS( hasher_->get_hash_map_size(),( 4*segment_vector_->at(4)->get_length() + 4*segment_vector_->at(5)->get_length() + 4*segment_vector_->at(6)->get_length() - 24) );
	}

	//generate_key (pretty easy to test)
	//Manually set atom's xyz coords and then check the key it generates
	void test_generate_key(){
		//Copy an atom from a segment
		core::conformation::Atom test_atom( segment_vector_->at( 1 )->get_const_residue_vector().at( 1 )->get_const_atom_vector().at( 1 ) );
		//set its xyz
		numeric::xyzVector< core::Real> temp_coords( 0.0, 0.0, 0.0 );
		test_atom.xyz( temp_coords );

		//Test that the function works at all
		hashing::HashKey result = hasher_->generate_key( test_atom );
		TS_ASSERT_EQUALS( result[ 1 ],0 );
		TS_ASSERT_EQUALS( result[ 2 ],0 );
		TS_ASSERT_EQUALS( result[ 3 ],0 );

	}
	//find_aligned_atoms( Atom const & start_atom)
	void test_find_aligned_atoms(){
		//Set up hash map
		TS_ASSERT_THROWS_NOTHING( hasher_->transform_segments( data_storage::Basis( 4, 5 ), true ) );
		TS_ASSERT_THROWS_NOTHING( hasher_->hash_segments( hasher_->get_seg1()->clone(), data_storage::Basis( 4, 5 ) ) );

		//Set up test atom
		core::conformation::Atom test_atom( segment_vector_->at( 1 )->get_const_residue_vector().at( 1 )->get_const_atom_vector().at( 1 ) );
		//set its xyz to that of a non-basis residue in 31
		numeric::xyzVector< core::Real> ref_coords( hasher_->get_seg1()->get_const_residue_vector().at( 3 )->get_const_atom_vector().at( 1 ).xyz() );
		numeric::xyzVector< core::Real> temp_coords( ref_coords );

		for ( int i = -1; i <= 1; ++i ) {
			temp_coords[ 0 ] = ref_coords[ 0 ] + ( i * 0.25 );
			for ( int j = -1; j <= 1; ++j ) {
				temp_coords[ 1 ] = ref_coords[ 1 ] + ( j * 0.25 );
				for ( int k = -1; k <= 1; ++k ) {
					temp_coords[ 2 ] = ref_coords[ 2 ] + ( k * 0.25 );
					test_atom.xyz( temp_coords );
					TS_ASSERT_EQUALS( hasher_->find_aligned_atoms( test_atom ).size(), 1 );
				}
			}
		}
		//Check that it fails with wrong coordinates
		temp_coords[ 2 ] = ref_coords[ 2 ] + 0.5;
		test_atom.xyz( temp_coords );
		TS_ASSERT_EQUALS( hasher_->find_aligned_atoms( test_atom ).size(), 0 );
	}

};
