// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/moves/MoverContainer.cxxtest.hh
/// @brief  tests for container Movers classes.
/// @author Sergey Lyskov

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>

#include <core/pose/Pose.hh>
#include <core/types.hh>

// Unit headers
#include <protocols/moves/Mover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/moves/MoverContainer.hh>

#include <basic/Tracer.hh>

//Auto Headers
#include <platform/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/signals/XYZEvent.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/kinematics/AtomTree.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
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
#include <protocols/filters/Filter.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/MoverContainer.fwd.hh>
#include <protocols/moves/MoverStatus.hh>
#include <protocols/moves/RepeatMover.fwd.hh>
#include <utility/down_cast.hh>
#include <utility/vector0.fwd.hh>
#include <utility/vector0.hh>
#include <utility/vector0_bool.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
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
#include <iosfwd>
#include <limits>
#include <list>
#include <map>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>
#include <basic/MetricValue.fwd.hh>
#include <basic/Tracer.fwd.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <boost/bind.hpp>
#include <boost/function.hpp>


/// We want to isolate static instances of DummyMover so we put it inside privet namespace.
namespace MoveContainerCxxTest {
	#include <test/protocols/moves/DummyMover.hh>
}
using namespace MoveContainerCxxTest;

using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("protocols.moves.MoverContainer.cxxtest");

// using declarations
using namespace core;
using namespace protocols::moves;

///////////////////////////////////////////////////////////////////////////
/// @name ContainerMoversTest
/// @brief: test for: SequenceMover
/// @detailed
/// @author Sergey Lyskov Fri Nov 02 2007
///////////////////////////////////////////////////////////////////////////
class ContainerMoversTest : public CxxTest::TestSuite {

public:
	ContainerMoversTest() {}

	void setUp() {
		extern int command_line_argc; extern char ** command_line_argv;
		if( command_line_argc > 1 ) core::init::init(command_line_argc, command_line_argv);
		else {
			std::string commandline = "core.test -mute all";
 				initialize_from_commandline_w_db( commandline );
		}

		core::init::init_random_generators(1000, "mt19937");
	}

	void tearDown() {
	}

	void test_RepeatMover() {
		pose::Pose pose;

		for(Size i=0; i<256; i++) {
			DummyMover::reset();
			DummyMover * dm = new DummyMover;
			RepeatMover RM(dm, i);

			RM.apply(pose);

			TS_ASSERT_EQUALS(dm->call_count(), int(i));
		}
	}


	void test_SequenceMover() {
		//TR << "test_SequenceMover...\n";
		const Size N = 100; ///< number of movers to test in sequence.

		pose::Pose pose;

		DummyMover::reset();
		SequenceMover SM;
		for(Size i=0; i<N; i++) SM.add_mover(new DummyMover(i));

		SM.apply(pose);

		TS_ASSERT_EQUALS(DummyMover::call_records().size(), N);
		if( DummyMover::call_records().size() == N ) {
			for(Size i=0; i<N; i++) {
				TS_ASSERT_EQUALS(DummyMover::call_records()[i], int(i));
			}
		}
		//TR << "test_SequenceMover... Ok.\n";
	}


	void test_CycleMover() {
		const Size N = 114; ///< number of movers.
		const Size NSteps = 1514; ///< number of steps to test.

		pose::Pose pose;

		CycleMover CM;
		for(Size i=0; i<N; i++) CM.add_mover(new DummyMover(i));

		for(Size s=0; s<NSteps; s++) {
			DummyMover::reset();
			CM.apply(pose);

			TS_ASSERT_EQUALS(DummyMover::call_records().size(), 1u);
			if( DummyMover::call_records().size() == 1 ) {
				int k = s % N;
				TS_ASSERT_EQUALS(DummyMover::call_records()[0], k);
			}
		}
	}

};

