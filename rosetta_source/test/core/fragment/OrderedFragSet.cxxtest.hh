// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/RotamerTrials.cxxtest.hh
/// @brief  test suite for rotamer_trials
/// @author Phil Bradley
/// @author Sergey Lyskov

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/OrderedFragSet.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/FrameIteratorWorker_.hh>
#include <core/fragment/FragID_Iterator.hh>

#include <core/fragment/Frame.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/FragCache.hh>
#include <core/fragment/BBTorsionSRFD.hh>
#include <core/fragment/util.hh>

#include <core/io/pdb/pose_io.hh>

#include <core/pack/rotamer_trials.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>



#include <core/pose/Pose.hh>

#include <core/scoring/ScoreFunction.hh>

#include <core/types.hh>

#include <basic/Tracer.hh>
#include <ObjexxFCL/string.functions.hh>
#include <numeric/random/random.hh>
#include <numeric/numeric.functions.hh>

//Auto Headers
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/VariantType.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/FragID_Iterator.fwd.hh>
#include <core/fragment/JumpingFrame.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID_Mask.hh>
#include <core/id/DOF_ID_Map.hh>
#include <core/id/DOF_ID_Mask.hh>
#include <core/id/NamedStubID.hh>
#include <core/io/pdb/file_data.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/signals/ConformationEvent.hh>
#include <core/pose/signals/DestructionEvent.hh>
#include <core/pose/signals/EnergyEvent.hh>
#include <core/scoring/types.hh>
#include <core/scoring/constraints/Constraints.fwd.hh>
#include <utility/stream_util.hh>
#include <utility/fix_boinc_read.hh>
#include <utility/keys/Key2Tuple.hh>
#include <ObjexxFCL/FArray.fwd.hh>
#include <set>


using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer tr("core.fragment.OrderedFragSet.cxxtest");
//MY_TRACERS("core.fragment.ConstantLengthFragments.cxxtest")

using namespace core;
using namespace fragment;
using namespace ObjexxFCL;

// hacky test code ---  should live somewhere else
// as defined in ConstantLengthFragments.cxxtest.hh
//extern void steal_constant_length_frag_set_from_pose ( pose::Pose const& pose, ConstantLengthFragSet& fragset );

class OrderedFragSetTest : public CxxTest::TestSuite
{
  chemical::ResidueTypeSetCAP residue_set_;
  pose::Pose pose_random_, pose_;
public:
  OrderedFragSetTest() {};

  // Shared initialization goes here.
  void setUp() {
		core_init();
		residue_set_ = chemical::ChemicalManager::get_instance()->residue_type_set( chemical::CENTROID );

		//core::import_pose::pose_from_pdb( pose_, "core/pack/test_in.pdb" );
		pose_ = create_test_in_pdb_pose();

  }
  //helpers
  void generate_random_pose();

  // tests
  void test_frag_iterator();

  // Shared finalization goes here.
  void tearDown() {
  }

private:

};

void OrderedFragSetTest::generate_random_pose() {
  using namespace chemical;
  using namespace conformation;
  std::string sequence = pose_.sequence();
  //create extended pose from sequence
  for ( Size pos = 1; pos <= sequence.size(); pos++ ) {
    chemical::AA aa = aa_from_oneletter_code( sequence[ pos-1 ] );
    ResidueTypeCAPs res_list = residue_set_->name3_map ( name_from_aa ( aa ) );
    ResidueOP new_rsd( ResidueFactory::create_residue( * ( res_list[ 1 ] ) ) );
    pose_random_.append_residue_by_bond( *new_rsd );
  }
  io::pdb::dump_pdb( pose_random_, "random_chain_pose" );
}


void OrderedFragSetTest::test_frag_iterator() {
  using namespace pose;
  using namespace fragment;
  kinematics::MoveMap movemap;
  movemap.set_bb( true );

  Size len9(9);
  ConstantLengthFragSet fragset9( len9 );
  steal_constant_length_frag_set_from_pose ( pose_, fragset9 );

  ConstantLengthFragSet fragset3( 3 );
  steal_constant_length_frag_set_from_pose ( pose_, fragset3 );

  OrderedFragSet fragset;
  // consolidate both fragsets in new OrderedFragSet
  for ( FrameIterator it = fragset3.begin(), eit = fragset3.end();
	it!=eit; ++it ) {
    fragset.add( *it );
  }

  for ( FrameIterator it = fragset9.begin(), eit = fragset9.end();
	it!=eit; ++it ) {
    fragset.add( *it );
  }

  // now there should be two frames with the same starting positions until we run out of
  // 9mers...  let's check this.
  Size ct = 1;
  for ( FrameIterator it = fragset.begin(), eit = fragset.end();
	it!=eit; ++it ) {
    //    tr.Info << " ct: " << ct << " " << (*it)->start() << "\n";
    TS_ASSERT( (*it)->start() == ct );
    if ( ct <= pose_.total_residue() - 8 ) ++it;
    TS_ASSERT( (*it)->start() == ct );
    ct++;
  }

  // a call to region should return us two Frames per position
  FrameList frames;
  fragset.region( movemap, 5, 5, 0, 0, frames);
  TS_ASSERT_EQUALS( frames.size() , 2 );

  // a call to region should return us two Frames per position
  frames.clear();
  fragset.region( movemap, 5, 6, 0, 0, frames);

  TS_ASSERT_EQUALS( frames.size(), 4 );

	//let's check the clone function
	FragSetOP cloned = fragset.clone();
	FrameIterator cit = cloned->begin();
	FrameIterator ecit = cloned->end();
  for ( FrameIterator it = fragset.begin(), eit = fragset.end();
				it!=eit; ++it,++cit ) {
    //    tr.Info << " ct: " << ct << " " << (*it)->start() << "\n";
    TS_ASSERT( cit != ecit );
		if (!( cit!=ecit )) break; //avoid run-time errors
		TS_ASSERT( (*it)->start() == (*cit)->start() );
  }
}

