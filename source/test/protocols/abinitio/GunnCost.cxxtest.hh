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
/// @author Oliver Lange

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>
// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
// AUTO-REMOVED #include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/conformation/ResidueFactory.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
// AUTO-REMOVED #include <core/fragment/FrameIterator.hh>
#include <core/fragment/FragID_Iterator.hh>

#include <core/fragment/Frame.hh>
// AUTO-REMOVED #include <core/fragment/FragData.hh>
#include <core/fragment/FragCache.hh>
// AUTO-REMOVED #include <core/fragment/BBTorsionSRFD.hh>
// AUTO-REMOVED #include <core/fragment/util.hh>

// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>

// AUTO-REMOVED #include <core/pack/rotamer_trials.hh>
// AUTO-REMOVED #include <core/pack/task/PackerTask.hh>
// AUTO-REMOVED #include <core/pack/task/TaskFactory.hh>



#include <core/pose/Pose.hh>

// AUTO-REMOVED #include <core/scoring/ScoreFunction.hh>

#include <core/types.hh>

#include <basic/Tracer.hh>
// AUTO-REMOVED #include <ObjexxFCL/string.functions.hh>

#include <protocols/simple_moves/GunnCost.hh>

//Auto Headers
#include <platform/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/AtomTypeSet.fwd.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/chemical/ElementSet.fwd.hh>
#include <core/chemical/MMAtomTypeSet.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/signals/XYZEvent.fwd.hh>
#include <core/fragment/BaseCacheUnit.hh>
#include <core/fragment/ConstantLengthFragSet.fwd.hh>
#include <core/fragment/FragCache.fwd.hh>
#include <core/fragment/FragData.fwd.hh>
#include <core/fragment/FragID.fwd.hh>
#include <core/fragment/FragID.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/fragment/Frame.fwd.hh>
#include <core/fragment/FrameIterator.fwd.hh>
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
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <protocols/simple_moves/FragmentMover.fwd.hh>
#include <protocols/simple_moves/FragmentMover.hh>
#include <protocols/simple_moves/SmoothFragmentMover.fwd.hh>
#include <protocols/simple_moves/SmoothFragmentMover.hh>
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
#include <utility/excn/EXCN_Base.hh>
#include <utility/excn/Exceptions.hh>
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
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <iosfwd>
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <ostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <basic/MetricValue.fwd.hh>
#include <basic/Tracer.fwd.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <boost/bind.hpp>
#include <boost/function.hpp>


using basic::T;
using basic::Error;
using basic::Warning;

//static basic::Tracer TR("core.fragment.ConstantLengthFragments.cxxtest");
std::string const TR("core.fragment.GunnCost.cxxtest");

using namespace core;
using namespace fragment;




class GunnCostTest : public CxxTest::TestSuite
{
	chemical::ResidueTypeSetCAP residue_set_;
	pose::Pose pose_random_, pose_;
public:
	GunnCostTest() {};

	// Shared initialization goes here.
	void setUp() {
		core_init();
		residue_set_ = chemical::ChemicalManager::get_instance()->residue_type_set( chemical::CENTROID );

		core::import_pose::pose_from_pdb( pose_, "protocols/abinitio/2GB3.pdb" );

		fragset3mer_ = ConstantLengthFragSetOP( new ConstantLengthFragSet( 3 ) );
		fragset3mer_->read_fragment_file( "protocols/abinitio/mfr_aa2GB3_03_05.200_v1_3" );
		fragset9mer_ = ConstantLengthFragSetOP( new ConstantLengthFragSet( 9 ) );
		fragset9mer_->read_fragment_file( "protocols/abinitio/mfr_aa2GB3_09_05.200_v1_3" );

		//generate_random_pose();
	}

	void test_cost();

	// Shared finalization goes here.
	void tearDown() {
	}

private:
	ConstantLengthFragSetOP fragset3mer_;
	ConstantLengthFragSetOP fragset9mer_;
};


void GunnCostTest::test_cost() {
	using namespace pose;
	using namespace fragment;
	using namespace protocols;
	using namespace simple_moves;
	kinematics::MoveMap movemap; //dummy ( functionality not used yet )
	//Size len (3);
	ConstantLengthFragSet& fragset=*fragset3mer_;
	FrameList frames;
	FrameOP aframe;
	GunnTuple data;
	GunnCost gc;
	movemap.set_bb( true );
	fragset.region(movemap,25,27,3,3,frames);
	aframe = frames[1];
	for (int i=1; i<=200; i++) {
		gc.compute_gunn( *aframe, i, data);
		T(TR, basic::t_error) << "25" << " " << i << " " << data.q1 << " " << data.q2 << " " << data.q3<< " " << data.q4<< " " << data.q5<< " " << data.q6 << "\n";
	}

	frames.clear();
	fragset.region(movemap,1,3,3,3,frames);
	aframe = frames[1];
	gc.compute_gunn( *aframe, 1, data);
	TS_ASSERT_DELTA( data.q1, -0.7429, 0.001 );
	TS_ASSERT_DELTA( data.q2, -0.6831, 0.001 );
	TS_ASSERT_DELTA( data.q3,  0.9926, 0.001 );
	TS_ASSERT_DELTA( data.q4,  0.7593, 0.001 );
	TS_ASSERT_DELTA( data.q5,  0.5410, 0.001 );
	TS_ASSERT_DELTA( data.q6,  5.9096, 0.001 );

	frames.clear();
	fragset.region(movemap,25,27,3,3,frames);
	aframe = frames[1];
	gc.compute_gunn( *aframe, 1, data);
	TS_ASSERT_DELTA( data.q1,  0.2617, 0.001 );
	TS_ASSERT_DELTA( data.q2, -0.8120, 0.001 );
	TS_ASSERT_DELTA( data.q3,  1.6389, 0.001 );
	TS_ASSERT_DELTA( data.q4,  0.4782, 0.001 );
	TS_ASSERT_DELTA( data.q5,  2.2299, 0.001 );
	TS_ASSERT_DELTA( data.q6,  5.4646, 0.001 );

	frames.clear();
	fragset.region(movemap,42,44,3,3,frames);
	aframe = frames[1];
	GunnTuple fr_data;
	float score_val[ 6 ] = { 10.964, 10.8948, 12.0666, 11.5839, 11.1111, 11.3222 };
	for (int i=1; i<=200; i++) {
		gc.compute_gunn( *aframe, i, fr_data);
		T(TR, basic::t_error) << "42" << " " << i << " " << fr_data.q1 << " " << fr_data.q2 << " " << fr_data.q3<< " " << fr_data.q4<< " " << fr_data.q5<< " " << fr_data.q6 << "score: " << gc.score_tuple(data,fr_data) << "\n";
		if (i<=6) {
			TS_ASSERT_DELTA( gc.score_tuple(data,fr_data), score_val[i-1], 0.01);
		}
	}
	// results that should come out:
	// tuple for 1 1 -0.742914 -0.683115 0.992633 0.759256 0.54098 5.90961
	// tuple for 1 2 -0.687352 -0.811139 0.672894 0.471779 0.154857 6.62757
	// tuple for 25 1 0.261736 -0.812008 1.63893 0.478188 2.22985 5.46585
	// these values are very close to those computed with rosetta++
	// the differences of 1e-2 can be explained with small (verified) differences  in bond-lengths and angles between my 2GB3.pdb (idealized) and the
	// standard geometry used in rosetta++.
}
