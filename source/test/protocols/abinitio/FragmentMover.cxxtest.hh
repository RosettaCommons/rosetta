// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/RotamerTrials.cxxtest.hh
/// @brief  test suite for rotamer_trials
/// @author Oliver Lange

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>

#include <protocols/simple_moves/FragmentMover.hh>

#include <core/scoring/ScoreFunction.fwd.hh>
// AUTO-REMOVED #include <core/scoring/ScoreType.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunctionFactory.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/FragID_Iterator.hh>
#include <core/fragment/FrameIteratorWorker_.hh>

#include <core/fragment/util.hh>

// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>



#include <core/pose/Pose.hh>
#include <core/types.hh>

#include <basic/Tracer.hh>

//Auto Headers
#include <platform/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/signals/XYZEvent.fwd.hh>
#include <core/fragment/BaseCacheUnit.hh>
#include <core/fragment/ConstantLengthFragSet.fwd.hh>
#include <core/fragment/FragData.fwd.hh>
#include <core/fragment/FragID.fwd.hh>
#include <core/fragment/FragID.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/Frame.fwd.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameIterator.fwd.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/FrameIteratorWorker_.fwd.hh>
#include <core/fragment/FrameList.fwd.hh>
#include <core/fragment/FrameList.hh>
#include <core/fragment/SingleResidueFragData.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/id/AtomID_Mask.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/DOF_ID_Map.fwd.hh>
#include <core/id/DOF_ID_Mask.fwd.hh>
#include <core/id/JumpID.fwd.hh>
#include <core/id/JumpID.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/TorsionID.hh>
#include <core/id/types.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/pdb/file_data.fwd.hh>
#include <core/kinematics/AtomTree.fwd.hh>
#include <core/kinematics/Edge.fwd.hh>
#include <core/kinematics/Edge.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/kinematics/types.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/datacache/ObserverCache.fwd.hh>
#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
#include <core/pose/signals/ConformationEvent.fwd.hh>
#include <core/pose/signals/DestructionEvent.fwd.hh>
#include <core/pose/signals/EnergyEvent.fwd.hh>
#include <core/pose/signals/GeneralEvent.fwd.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <protocols/simple_moves/FragmentMover.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverStatus.hh>
#include <utility/py/PyAssert.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
#include <utility/fix_boinc_read.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/io/izstream.fwd.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/signals/BufferedSignalHub.fwd.hh>
#include <utility/signals/BufferedSignalHub.hh>
#include <utility/signals/Link.fwd.hh>
#include <utility/signals/Link.hh>
#include <utility/signals/LinkUnit.fwd.hh>
#include <utility/signals/LinkUnit.hh>
#include <utility/signals/SignalHub.fwd.hh>
#include <utility/signals/SignalHub.hh>
#include <utility/tag/Tag.fwd.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzVector.fwd.hh>
#include <ObjexxFCL/Dimension.fwd.hh>
#include <ObjexxFCL/Dimension.hh>
#include <ObjexxFCL/DimensionExpression.hh>
#include <ObjexxFCL/DynamicIndexRange.fwd.hh>
#include <ObjexxFCL/DynamicIndexRange.hh>
#include <ObjexxFCL/FArray.fwd.hh>
#include <ObjexxFCL/FArray.hh>
#include <ObjexxFCL/FArray1.fwd.hh>
#include <ObjexxFCL/FArray1.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2.fwd.hh>
#include <ObjexxFCL/FArray2.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArrayInitializer.fwd.hh>
#include <ObjexxFCL/FArrayInitializer.hh>
#include <ObjexxFCL/FArraySection.fwd.hh>
#include <ObjexxFCL/FArraySection.hh>
#include <ObjexxFCL/FArrayTraits.fwd.hh>
#include <ObjexxFCL/FArrayTraits.hh>
#include <ObjexxFCL/IndexRange.fwd.hh>
#include <ObjexxFCL/IndexRange.hh>
#include <ObjexxFCL/InitializerSentinel.hh>
#include <ObjexxFCL/Observer.fwd.hh>
#include <ObjexxFCL/Observer.hh>
#include <ObjexxFCL/ObserverMulti.hh>
#include <ObjexxFCL/ObserverSingle.hh>
#include <ObjexxFCL/ProxySentinel.hh>
#include <ObjexxFCL/SetWrapper.fwd.hh>
#include <ObjexxFCL/Star.fwd.hh>
#include <ObjexxFCL/Star.hh>
#include <ObjexxFCL/proxy_const_assert.hh>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <execinfo.h>
#include <iosfwd>
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <ostream>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <basic/MetricValue.fwd.hh>
#include <basic/Tracer.fwd.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/functional/hash.hpp>


using basic::T;
using basic::Error;
using basic::Warning;

//static basic::Tracer TR("core.fragment.ConstantLengthFragments.cxxtest");

using namespace core;
using namespace fragment;


static basic::Tracer tr("protocols.simple_moves.FragmentMover.cxxtest", basic::t_info );

class FragmentMoverTest : public CxxTest::TestSuite
{
  pose::Pose pose_;
public:
  FragmentMoverTest() {};

  // Shared initialization goes here.
  void setUp() {
    core_init();
    core::import_pose::pose_from_pdb( pose_, "protocols/abinitio/2GB3.pdb" );

    fragset3mer_ = ConstantLengthFragSetOP( new ConstantLengthFragSet );
    fragset3mer_->read_fragment_file( "protocols/abinitio/mfr_aa2GB3_03_05.200_v1_3" );

    TS_ASSERT_EQUALS( fragset3mer_->max_frag_length(), 3 );
  }
  void test_proposed_ss();
	void test_ss_check();
  void test_mover();
  void test_window_start();
	void test_steal();
  // Shared finalization goes here.
  void tearDown() {
  }

private:
  ConstantLengthFragSetOP fragset3mer_;
};


void FragmentMoverTest::test_mover() {
  using namespace pose;
  using namespace fragment;
  using namespace protocols;
  using namespace simple_moves;
  using namespace scoring;
  kinematics::MoveMapOP movemap( new kinematics::MoveMap ); //dummy ( functionality not used yet )
  //Size len (3);
  movemap->set_bb( true );
  //	for (Size ii=1; ii<=6; ii++ ) {
  //	movemap->set_bb( ii, false );
  //	movemap->set_bb( pose_.total_residue()-ii, false );
  //	TS_ASSERT( !movemap->get_bb( ii ) );
  //}
	// Pose pose = pose_;
  //FragmentMover mover( fragset3mer_, movemap, new GunnCost( 7.0 ) );
  //	moves::TrialMover wobble_min_trial

  //for ( Size i=1; i<= 100; i++ ) {
	// tr.Debug << " ----------------------------- " << i << " -------------------------------- " << std::endl;
    //wobbles.apply( pose );
  //  pose=pose_;
  //};
}

void FragmentMoverTest::test_steal() {
  using namespace pose;
  using namespace fragment;
  using namespace protocols;
  using namespace simple_moves;
  using namespace scoring;

	for ( ConstFrameIterator it = fragset3mer_->begin(), eit = fragset3mer_->end(); it!=eit; ++it ){
		TS_ASSERT_EQUALS( (*it)->nr_frags() , 200 );
	}
	steal_constant_length_frag_set_from_pose( pose_, *fragset3mer_ );
	for ( ConstFrameIterator it = fragset3mer_->begin(), eit = fragset3mer_->end(); it!=eit; ++it ){
		TS_ASSERT_EQUALS( (*it)->nr_frags() , 201 );
	}

}

//checks whether the proposed secondary structure is correct
// does not check whether the ss-check is applied correctly in FragmentMover
void FragmentMoverTest::test_proposed_ss() {
  using namespace pose;
  using namespace fragment;
  using namespace protocols;
  using namespace simple_moves;
  using namespace scoring;



  std::string proposed_ss;
  proposed_ss.reserve( pose_.total_residue() );
  proposed_ss = pose_.secstruct(); // full ss-string from pose
  std::string old_ss = proposed_ss;
  kinematics::MoveMapOP movemap( new kinematics::MoveMap ); //dummy ( functionality not used yet )
  //Size len (3);
  movemap->set_bb( true );
  ClassicFragmentMover mover ( fragset3mer_, movemap );
  TS_ASSERT( mover.valid_ss( proposed_ss ) );
  std::string full_result("");
	full_result.reserve( 130*200 );
  int ct = 1;
  for ( FragID_Iterator it = fragset3mer_->begin(), eit = fragset3mer_->end();
	it!=eit; ++it , ++ ct) {
		for ( int i = 1; i<=10 && it!=eit;  ++i ) ++it;
		if ( it == eit ) break;
    proposed_ss = old_ss;
    it->apply_ss( *movemap, proposed_ss );
    tr.Debug << "proposed_ss " << proposed_ss << std::endl;
    if ( mover.valid_ss( proposed_ss ) ) {
      full_result.push_back ( 'y' );
    } else {
      full_result.push_back( 'n' );
    }

    //each 100 fragment also try if it is the same result if we apply it to pose directly
    if ( ct % 1000  == 0 ) {
			tr.Info << "run full ss-check on fragment insertion" << std::endl;
      pose::Pose scratch_pose = pose_;
      it->apply( *movemap, scratch_pose );
      TS_ASSERT( scratch_pose.secstruct() == proposed_ss );
    }
  }

  tr.Info << "valid_ss history " << std::endl;
	//std::ofstream out( "valid_ss.dat" );
	//out << full_result << std::endl;
	std::string check_result;
	std::ifstream in( "protocols/abinitio/valid_ss.dat" );
	in >> check_result;
	TS_ASSERT( check_result == full_result );
}


//checks whether the proposed secondary structure is correct
// does not check whether the ss-check is applied correctly in FragmentMover
void FragmentMoverTest::test_ss_check() {
  using namespace pose;
  using namespace fragment;
  using namespace protocols;
  using namespace simple_moves;
  using namespace scoring;

	std::string const proposed_ss( pose_.secstruct() ); // full ss-string from pose
  kinematics::MoveMapOP movemap( new kinematics::MoveMap );
  movemap->set_bb( true );
  ClassicFragmentMover mover ( fragset3mer_, movemap );

	//check that we start with valid_ss
  TS_ASSERT( mover.valid_ss( proposed_ss ) );

	pose::Pose scratch_pose = pose_;
	tr.Info << "do a couple of fragment moves and check whether wrong ss is accepted" << std::endl;
	for ( Size ct = 1; ct <= 1000; ct++ ) {
		mover.apply( scratch_pose );
		//check that mover does not accept invalid ss
		TS_ASSERT( mover.valid_ss( scratch_pose.secstruct() ) );
	}
}


void FragmentMoverTest::test_window_start() {
  using namespace pose;
  using namespace fragment;
  using namespace protocols;
  using namespace simple_moves;
  using namespace scoring;

  // test if all insertable positions are chosen
  kinematics::MoveMapOP movemap( new kinematics::MoveMap ); //dummy ( functionality not used yet )
  movemap->set_bb( true );
	std::map< Size, Size > inserted_pos; // make it a set such that it will be sorted and no positions are doubled
  core::fragment::InsertMap insert_map;
	core::fragment::InsertSize insert_size;

	// disallow some residues to save sampling time ( and to make it more tricky )
	// total_residue/2.0 must be sampled, since it is used for the normalization further down
	Size center_pos = static_cast< int >( pose_.total_residue()/2.0 );
	for ( Size i= center_pos + 10; i<=pose_.total_residue(); i++ ) {
		movemap->set_bb( i, false );
	}
	for ( int i = 6; i<= static_cast< int > ( pose_.total_residue()/3.0 ); i++ ) {
		movemap->set_bb( i, false );
	}

	// there is already a check for generate_insert_map in core/fragment ... should be correct here
  fragset3mer_->generate_insert_map( *movemap, insert_map, insert_size );

	ClassicFragmentMover mover ( fragset3mer_, movemap );

	// this is a stochastic test, need to check the distribution --> enough sampling required.
	// thus we repeat with more sampling if threre is failure
	// in this way we don't spend too much time for daily testing but make sure that we don't fail
	// due to random sampling errors
	//
	Size pass = 1;
	Size success = false;
	Size nmoves = 20000; // starts with 50k -- relatively fast but some chance to succeed
	Real const end_bias ( 30.0 );

	for ( pass = 1; pass <= 3 && !success; pass ++ ) {
		inserted_pos.clear();
		nmoves *= 5;
		tr.Info << " test start distribution with " << nmoves << " trials... " << std::endl;
		for ( Size ct =1 ; ct <= nmoves; ct++ ) {
			Size pos;
			if ( mover.choose_window_start( pose_, fragset3mer_->max_frag_length(), pos ) ) {
				inserted_pos[ pos ] += 1;
			}
		}

		// all insertable positions are sampled ?
		TS_ASSERT_EQUALS( insert_map.size(), inserted_pos.size() );

		// find maximum
		Size max = inserted_pos[ center_pos ];
		//for ( Size ct = 1; ct <= insert_map.size(); ct++ ) {
		//	Size m = inserted_pos[ insert_map[ ct ] ];
		//	if ( m > max ) max = m;
		//}

		// check for bias
		success = true;
		for ( Size ct = 1; ct <= insert_map.size() && success; ct++ ) {
			Size begin = insert_map[ ct ];
			Size min_fixed_residues;
			Size const fixed_residues =
				pose_.fold_tree().count_fixed_residues( begin, insert_size[ begin ], min_fixed_residues );

			Real prob = std::exp( 1.0* ( (Real) min_fixed_residues - fixed_residues ) / end_bias );
			Real freq = 1.0/max * inserted_pos[ begin ];
			success = success && ( std::abs( prob - freq ) < 0.04 );
			if ( !success ) tr.Info << "inaccurate for position " << begin << " with " << prob << " vs " << freq << std::endl;
			// this checks for a correct distribution of positions: if slightly inaccurate -- maybe just more sampling required ? --> nmoves
		}
	} // passes
	//	TS_ASSERT( success );
	if ( !success ) {
		Size max = inserted_pos[ center_pos ];
		// make notification to test-system
			for ( Size ct = 1; ct <= insert_map.size() && success; ct++ ) {
				Size begin = insert_map[ ct ];
				Size min_fixed_residues;
				Size const fixed_residues =
					pose_.fold_tree().count_fixed_residues( begin, insert_size[ begin ], min_fixed_residues );

				Real prob = std::exp( 1.0* ( (Real) min_fixed_residues - fixed_residues ) / end_bias );
				Real freq = 1.0/max * inserted_pos[ begin ];
				TS_ASSERT_DELTA( prob, freq, 0.04 );
			}
	}
}
