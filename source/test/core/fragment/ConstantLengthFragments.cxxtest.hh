// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/RotamerTrials.cxxtest.hh
/// @brief  test suite for rotamer_trials
/// @author Phil Bradley
/// @author Sergey Lyskov

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeFinder.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/FrameIteratorWorker_.hh>
#include <core/fragment/FragID_Iterator.hh>

#include <core/fragment/Frame.hh>
#include <core/fragment/FragCache.hh>
#include <core/fragment/util.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/types.hh>

#include <basic/Tracer.hh>
#include <ObjexxFCL/string.functions.hh>
#include <numeric/numeric.functions.hh>

//Auto Headers
#include <core/id/AtomID_Mask.hh>
#include <utility/fix_boinc_read.hh>
#include <utility/vector1.hh>


using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer tr("core.fragment.ConstantLengthFragments.cxxtest");
//MY_TRACERS("core.fragment.ConstantLengthFragments.cxxtest")

using namespace core;
using namespace fragment;
using namespace ObjexxFCL;

// // hacky test code ---  should live somewhere else
// void steal_constant_length_frag_set_from_pose ( pose::Pose const& pose, ConstantLengthFragSet& fragset ) {
//  Size nbb ( 3 ); // three backbone torsions for Protein
//  Size len = fragset.max_frag_length();
//  for ( Size pos = 1; pos <= pose.size() - len + 1; ++pos ) {
//   FragDataOP frag_raw = new FragData;
//   for ( Size i = 1; i<= len; i++ ) {
//    frag_raw->add_residue( new BBTorsionSRFD( nbb, pose.secstruct(pos), oneletter_code_from_aa(pose.residue( pos ).aa() ) ) );
//   };
//   FrameOP frame = new Frame( pos, len );
//   frag_raw->steal( pose, *frame );
//   frame->add_fragment ( frag_raw );
//   fragset.add( frame );
//  };
//  }


class FragmentConstantLengthTest : public CxxTest::TestSuite
{
	chemical::ResidueTypeSetCOP residue_set_;
	pose::Pose pose_random_, pose_;
public:
	FragmentConstantLengthTest() {};

	// Shared initialization goes here.
	void setUp() {
		core_init();
		residue_set_ = chemical::ChemicalManager::get_instance()->residue_type_set( chemical::CENTROID );

		//core::import_pose::pose_from_file( pose_, "core/pack/test_in.pdb" , core::import_pose::PDB_file);
		pose_ = create_test_in_pdb_pose();
	}
	//helpers
	void generate_random_pose();
	void sub_insertion( Size size, fragment::FragSet const& );

	// tests
	void test_fragment_stealing_insertion();
	void test_frag_cache();
	void test_frag_iterator();
	void test_insertmap();
	// Shared finalization goes here.
	void tearDown() {
	}

private:

};

void FragmentConstantLengthTest::generate_random_pose () {
	using namespace chemical;
	using namespace conformation;
	std::string sequence = pose_.sequence();
	//create extended pose from sequence
	for ( Size pos = 1; pos <= sequence.size(); pos++ ) {
		chemical::AA aa = aa_from_oneletter_code( sequence[ pos-1 ] );
		// Really, this could be changed to get_representative_type
		// But let's keep the unit test as-is
		ResidueTypeCOPs res_list = ResidueTypeFinder( *residue_set_).name3( name_from_aa ( aa ) ).get_all_possible_residue_types();
		ResidueOP new_rsd( ResidueFactory::create_residue( * ( res_list[ 1 ] ) ) );
		pose_random_.append_residue_by_bond( *new_rsd );
	}
	io::pdb::dump_pdb( pose_random_, "random_chain_pose" );
}

///////////////////////////////////////////////////////////////////////////////
// ------------------------------------------ //
/// @brief test for rotamer_trials
void FragmentConstantLengthTest::test_fragment_stealing_insertion()
{
	//using namespace core;
	using namespace conformation;
	using namespace chemical;
	using namespace scoring;
	using namespace pose;
	using namespace fragment;

	ConstantLengthFragSet fragset9mer( 9 );
	// test code is not linked against lib_protocls // nor should it for a test that lives in core!
	steal_constant_length_frag_set_from_pose ( pose_, fragset9mer );

	ConstantLengthFragSet fragset3mer( 3 );
	steal_constant_length_frag_set_from_pose ( pose_, fragset3mer );
	tr.Info << "stealing-test: run 3 mer insertions" << std::endl;
	sub_insertion( 3, fragset3mer );
	tr.Info << "stealing-test: run 9mer insertions" << std::endl;
	sub_insertion( 9, fragset9mer );
	tr.Info << "finished stealing test" << std::endl;
}

void FragmentConstantLengthTest::test_frag_cache() {
	using namespace pose;
	using namespace fragment;
	kinematics::MoveMap movemap;
	movemap.set_bb( true );
	Size len (9);
	ConstantLengthFragSet fragset( len );

	// test code is not linked against lib_protocls // nor should it for a test that lives in core!
	steal_constant_length_frag_set_from_pose ( pose_, fragset );

	{ // test FragCache
		FragCache< Size > silly_cache("ULTIMATE_SILLINESS");
		for ( Size pos = 1; pos<= pose_.size() ; pos++ ) {
			FrameList frames;
			if ( fragset.region( movemap, pos, pos, len, len, frames ) ) {
				Frame const& frame ( * ( frames[ 1 ] ) );
				silly_cache.store(frame, 1, pos);
			} else {
				TS_ASSERT( pos>pose_.size()-len+1 );
			};
		};
		FragCache< Real > empty_cache("NO_VALUES_HERE");
		FragCache< Size > another_silly_cache("ULTIMATE_SILLINESS");
		for ( Size pos = 1; pos <=  pose_.size(); pos++ ) {
			FrameList frames;
			if ( fragset.region( movemap, pos, pos, len, len, frames ) ) {
				Frame const& frame ( * ( frames[ 1 ] ) );
				Size val;
				TS_ASSERT( silly_cache.retrieve(frame, 1, val) ); //there should be a value ( return true )
				TS_ASSERT_EQUALS( val , pos );
				Real data;
				TS_ASSERT( !empty_cache.retrieve(frame, 1, data));
				Size val2;
				TS_ASSERT( another_silly_cache.retrieve(frame, 1, val2) ); //there should be a value ( return true )
				TS_ASSERT_EQUALS( val2 , pos );

				// for faster access to all values of the frame:
				FragCache< Size >::TCacheUnit& cache = silly_cache( frame );
				Size val3;
				TS_ASSERT( cache.retrieve( frame.frag_id( 1 ), val3 ) );
				TS_ASSERT_EQUALS( val3, pos );
			};
		}

#if 0
		// if compiled the following SHOULD throw a bad_cast exception
		// It did today! Oct 25 2007
		FragCache< Size > wrong_type_cache("NO_VALUES_HERE");
		FrameList frames;
		Size val4;
		fragset.region( movemap, 1, 1, 3, 3, frames);
		wrong_type_cache.retrieve( * (frames[1] ), 1, val4 );
#endif
	} // block
	////////////////////////////////////////////////////
	/// now test the FragStore
	{
		FragStore< Size > silly_cache("ULTIMATE_SILLINESS_STORE");
		for ( Size pos = 1; pos<= pose_.size() ; pos++ ) {
			FrameList frames;
			if ( fragset.region( movemap, pos, pos, len, len, frames ) ) {
				Frame const& frame ( * ( frames[ 1 ] ) );
				silly_cache.store(frame, 1, pos);
			} else {
				TS_ASSERT( pos>pose_.size()-len+1 );
			};
		};
		FragStore< Real > empty_cache("NO_VALUES_HERE_STORE");
		FragStore< Size > another_silly_cache("ULTIMATE_SILLINESS_STORE");
		for ( Size pos = 1; pos <=  pose_.size(); pos++ ) {
			FrameList frames;
			if ( fragset.region( movemap, pos, pos+len-1, len, len, frames ) ) {
				Frame const& frame ( * ( frames[ 1 ] ) );
				Size val;
				TS_ASSERT( silly_cache.retrieve(frame, 1, val) ); //there should be a value ( return true )
				TS_ASSERT_EQUALS( val , pos );
				//Real data;
				// TS_ASSERT( empty_cache.retrieve(frame, 1, data)); // this will break since no value has been stored
				Size val2;
				TS_ASSERT( another_silly_cache.retrieve(frame, 1, val2) ); //there should be a value ( return true )
				TS_ASSERT_EQUALS( val2 , pos );

				// for faster access to all values of the frame:
				FragStore< Size >::TCacheUnit& cache = silly_cache( frame );
				Size val3;
				TS_ASSERT( cache.retrieve( frame.frag_id( 1 ), val3 ) );
				TS_ASSERT_EQUALS( val3, pos );
			};
		}
	} // block for testing FragStore
}

void FragmentConstantLengthTest::test_frag_iterator() {
	using namespace pose;
	using namespace fragment;
	kinematics::MoveMap movemap;
	movemap.set_bb( true );
	Size len (9);
	ConstantLengthFragSet fragset( len );

	// test code is not linked against lib_protocls // nor should it for a test that lives in core!
	steal_constant_length_frag_set_from_pose ( pose_, fragset );

	{ // test FragCache
		FragCache< Size > silly_cache("ULTIMATE_SILLINESS");
		for ( Size pos = 1; pos<= pose_.size() ; pos++ ) {
			FrameList frames;
			if ( fragset.region( movemap, pos, pos, len, len, frames ) ) {
				Frame const& frame ( * ( frames[ 1 ] ) );
				silly_cache.store(frame, 1, pos);
			} else {
				TS_ASSERT( pos>pose_.size()-len+1 );
			};
		};
		FragSet& bfragset( fragset );
		ConstFrameIterator it = bfragset.begin();
		ConstFrameIterator eit= bfragset.end();

		for ( Size pos = 1; pos <=  pose_.size(); pos++ ) {
			FrameList frames;
			if ( fragset.region( movemap, pos, pos, len, len, frames ) ) {
				Frame const& frame ( * ( frames[ 1 ] ) );
				Size val;
				TS_ASSERT( silly_cache.retrieve(frame, 1, val) ); //there should be a value ( return true )
				TS_ASSERT_EQUALS( val , pos );

				Size val2;
				TS_ASSERT( it != eit );
				TS_ASSERT( silly_cache.retrieve( **it , 1, val2) ); //there should be a value ( return true )
				TS_ASSERT_EQUALS( val2 , val );
				++it;
			} // block for testing FragStore
		}
		// since we stopped one residue early there should be one left...
		TS_ASSERT( it == eit );
		{
			FragSet& bfragset( fragset );
			ConstFrameIterator it = bfragset.begin();
			ConstFrameIterator eit= bfragset.end();
			for ( Size pos=1 ; it!=eit; ++it ) {
				Size val;
				silly_cache.retrieve( **it, 1, val);
				TS_ASSERT_EQUALS( val, pos );
				++pos;
			};
		}
		{
			FragSet& bfragset( fragset );
			FragID_Iterator it = bfragset.begin();
			FragID_Iterator eit= bfragset.end();
			for ( Size pos=1 ; it!=eit; ++it ) {
				Size val;
				silly_cache.retrieve( *it, val);
				TS_ASSERT_EQUALS( val, pos );
				++pos;
			};
		}
		{

			FragCache< Size >::ScoredList scored_frags;
			silly_cache.scored_frag_ids( scored_frags, fragset.begin(), fragset.end() );
			FragID_List list;
			for ( int i = 1; i<= (int) scored_frags.size(); i++ ) {
				tr.Info << "score: "<< scored_frags[ i ].second << " frag_ID: " << scored_frags[ i ].first.first << std::endl;
				if ( numeric::mod( i,2 ) ) list.push_back(scored_frags[i].first );
			}
			ConstantLengthFragSet new_frag_set( len );
			new_frag_set.insert_fragID_list( list );

			FragCache< Size >::ScoredList scored_frags2;
			silly_cache.scored_frag_ids( scored_frags2, new_frag_set.begin(), new_frag_set.end() );
			for ( Size i = 1; i<=scored_frags2.size(); i++ ) {
				tr.Info << "score2: "<< scored_frags2[ i ].second << " frag_ID: " << scored_frags2[ i ].first.first << std::endl;
			}

		}
		{
			FrameList new_frame_list;
			std::copy(fragset.nonconst_begin(), fragset.nonconst_end(), std::back_inserter(new_frame_list) );
			FragID_Iterator it = new_frame_list.begin();
			FragID_Iterator eit = new_frame_list.end();
			for ( Size pos=1 ; it!=eit; ++it ) {
				Size val;
				silly_cache.retrieve( *it, val);
				TS_ASSERT_EQUALS( val, pos );
				++pos;
			}
		}
	}
}

void FragmentConstantLengthTest::test_insertmap() {
	using namespace pose;
	using namespace fragment;

	kinematics::MoveMap movemap;
	movemap.set_bb( false );

	Size len (3);
	ConstantLengthFragSet fragset( len );
	// test code is not linked against lib_protocls // nor should it for a test that lives in core!
	steal_constant_length_frag_set_from_pose ( pose_, fragset );

	for ( int i = 20; i<=30; i++ ) {
		movemap.set_bb( i, true );
	}
	InsertMap insert_map; InsertSize insert_size;
	fragset.generate_insert_map( movemap, insert_map, insert_size );
	int pos = 17;
	for ( InsertMap::const_iterator it=insert_map.begin(), eit=insert_map.end(); it!=eit; ++it ) {
		tr.Info << "active residue " << *it << std::endl;
		TS_ASSERT_EQUALS( *it, ++pos );
		TS_ASSERT_EQUALS( insert_size[ *it ], std::min( std::min( 3, 30-pos+1 ) , pos-17 ) )
			}
			TS_ASSERT_EQUALS( pos, 30 ); // that should be the last number

		//let's check that it doesn't find anything if we have less than 3 residues
		{
			movemap.set_bb( false );
			movemap.set_bb( 20, true ); //2mer
			movemap.set_bb( 21, true );

			movemap.set_bb( 30, true ); //1mer

			InsertMap insert_map; InsertSize insert_size;
			fragset.generate_insert_map( movemap, insert_map, insert_size );
			TS_ASSERT_EQUALS( insert_map.size(), 7 );
			TS_ASSERT_EQUALS( insert_size[ 20 ], 2 );
			TS_ASSERT_EQUALS( insert_size[ 21 ], 1 );
			TS_ASSERT_EQUALS( insert_size[ 30 ], 1 );
			TS_ASSERT_EQUALS( insert_size[ 14 ], 0 );
			TS_ASSERT_EQUALS( insert_size[ 19 ], 2 );
		}

		//let's check that it doesn't find anything if we have less than 3 residues
		{
			movemap.set_bb( false );
			movemap.set_bb( 20, true ); //3mer
			movemap.set_bb( 21, true );
			movemap.set_bb( 22, true );

			movemap.set_bb( 30, true ); //3mer
			movemap.set_bb( 31, true ); //
			movemap.set_bb( 32, true ); //

			InsertMap insert_map; InsertSize insert_size;
			fragset.generate_insert_map( movemap, insert_map, insert_size );

			tr.Info << "==============================================" << std::endl;
			for ( InsertMap::const_iterator it=insert_map.begin(), eit=insert_map.end(); it!=eit; ++it ) {
				tr.Info << "active residue " << *it << std::endl;
			}

			TS_ASSERT_EQUALS( insert_map.size(), 10 );
			if ( insert_map.size() == 6 ) {
				TS_ASSERT_EQUALS( insert_map[ 1 ], 18 );
				TS_ASSERT_EQUALS( insert_map[ 6 ], 28 );
			}
		}

	}


	void FragmentConstantLengthTest::sub_insertion( Size size, fragment::FragSet const& fragsetNmer ) {
		using namespace conformation;
	using namespace chemical;
	using namespace scoring;
	using namespace pose;
	using namespace fragment;
	kinematics::MoveMap movemap;
	movemap.set_bb( true );
	// fold via Nmers
	Pose pose_Nmers ( pose_ );
	for ( Size pos = 1; pos<= pose_Nmers.size() - size; pos++ ) {
		for ( Size i = 1; i <= 3; i++ ) {
			pose_Nmers.set_torsion( id::TorsionID( pos, id::BB, i ), numeric::random::uniform() );
		}
	};
	io::pdb::dump_pdb( pose_Nmers, "randomized.pdb");

	for ( Size pos = 1; pos<= pose_Nmers.size() - size; pos+= ( size / 2 ) ) {
		FrameList frames;
		fragsetNmer.region_all(  pos, pos, size, size, frames );
		TS_ASSERT_EQUALS( frames.size(), 1 );
		Frame const& frame ( * ( frames[ 1 ] ) );
		TS_ASSERT_EQUALS( frame.nr_frags(), 1 );
		TS_ASSERT_EQUALS( frame.length(), size );
		TS_ASSERT_EQUALS( frame.start(), pos );
		TS_ASSERT_EQUALS( frame.end(), pos + size - 1 );
		frame.apply( movemap, 1, pose_Nmers ); // apply the first fragment of this frame
	}
	io::pdb::dump_pdb( pose_Nmers, "pose_after_"+right_string_of(size,1)+"mer_insertion");
	for ( Size pos = 1; pos<= pose_Nmers.size(); pos ++ ) {
		for ( Size i = 1; i <= 3; i++ ) {
			TS_ASSERT_EQUALS( pose_Nmers.torsion( id::TorsionID( pos, id::BB, i) ), pose_.torsion( id::TorsionID( pos, id::BB, i) ) );
		}
	}
}
