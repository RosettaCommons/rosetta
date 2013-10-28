// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   rosetta/benchmark/ShearMover.bench.cc
///
/// @brief  Varios moves benchmark
/// @author Sergey Lyskov


#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <protocols/simple_moves/BackboneMover.hh>

#include <core/kinematics/MoveMap.hh>

#include <apps/benchmark/performance/performance_benchmark.hh>

#include <utility/vector1.hh>

//Auto Headers
//#include <platform/types.hh>
//#include <core/types.hh>
//#include <core/chemical/AA.hh>
//#include <core/chemical/ResidueType.fwd.hh>
//#include <core/chemical/ResidueTypeSet.fwd.hh>
//#include <core/conformation/Conformation.fwd.hh>
//#include <core/conformation/Residue.fwd.hh>
//#include <core/conformation/signals/XYZEvent.fwd.hh>
//#include <core/id/AtomID.fwd.hh>
//#include <core/id/AtomID.hh>
//#include <core/id/AtomID_Map.fwd.hh>
//#include <core/id/AtomID_Map.hh>
//#include <core/id/AtomID_Mask.fwd.hh>
//#include <core/id/AtomID_Mask.hh>
//#include <core/id/DOF_ID.fwd.hh>
//#include <core/id/DOF_ID.hh>
//#include <core/id/DOF_ID_Map.fwd.hh>
//#include <core/id/DOF_ID_Map.hh>
//#include <core/id/DOF_ID_Mask.fwd.hh>
//#include <core/id/DOF_ID_Mask.hh>
//#include <core/id/JumpID.fwd.hh>
//#include <core/id/JumpID.hh>
//#include <core/id/NamedAtomID.fwd.hh>
//#include <core/id/NamedAtomID.hh>
//#include <core/id/NamedStubID.fwd.hh>
//#include <core/id/NamedStubID.hh>
//#include <core/id/TorsionID.fwd.hh>
//#include <core/id/TorsionID.hh>
//#include <core/id/types.hh>
//#include <core/import_pose/file_data.fwd.hh>
//#include <core/io/pdb/file_data.hh>
//#include <core/kinematics/AtomTree.fwd.hh>
//#include <core/kinematics/FoldTree.fwd.hh>
//#include <core/kinematics/Jump.fwd.hh>
//#include <core/kinematics/MoveMap.fwd.hh>
//#include <core/kinematics/Stub.fwd.hh>
//#include <core/kinematics/types.hh>
//#include <core/pose/PDBInfo.fwd.hh>
//#include <core/pose/Pose.fwd.hh>
//#include <core/pose/datacache/ObserverCache.fwd.hh>
//#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
//#include <core/pose/signals/ConformationEvent.fwd.hh>
//#include <core/pose/signals/ConformationEvent.hh>
//#include <core/pose/signals/DestructionEvent.fwd.hh>
//#include <core/pose/signals/DestructionEvent.hh>
//#include <core/pose/signals/EnergyEvent.fwd.hh>
//#include <core/pose/signals/EnergyEvent.hh>
//#include <core/pose/signals/GeneralEvent.fwd.hh>
//#include <core/pose/signals/GeneralEvent.hh>
//#include <core/scoring/Energies.fwd.hh>
//#include <core/scoring/ScoreFunction.fwd.hh>
//#include <core/scoring/ScoreFunctionInfo.fwd.hh>
//#include <core/scoring/constraints/Constraint.fwd.hh>
//#include <core/scoring/constraints/ConstraintSet.fwd.hh>
//#include <core/scoring/constraints/Constraints.fwd.hh>
//#include <basic/MetricValue.fwd.hh>
//// AUTO-REMOVED #include <basic/OStream.fwd.hh>
//#include <basic/Tracer.fwd.hh>
//#include <basic/datacache/BasicDataCache.fwd.hh>
//#include <protocols/filters/Filter.fwd.hh>
//#include <protocols/simple_moves/BackboneMover.fwd.hh>
//#include <basic/datacache/DataMap.fwd.hh>
//#include <protocols/moves/MonteCarlo.fwd.hh>
//#include <protocols/moves/Mover.fwd.hh>
//#include <protocols/moves/Mover.hh>
//#include <protocols/moves/MoverStatistics.hh>
//#include <protocols/moves/MoverStatus.hh>
//#include <utility/down_cast.hh>
//#include <utility/exit.hh>
//#include <utility/vector1.fwd.hh>
//#include <utility/vector1.hh>
//#include <utility/vector1_bool.hh>
//#include <utility/vectorL.fwd.hh>
//#include <utility/vectorL.hh>
//#include <utility/vectorL_Selector.hh>
//#include <utility/vectorL_bool.hh>
//#include <utility/tag/Tag.fwd.hh>
//#include <utility/pointer/ReferenceCount.fwd.hh>
//#include <utility/pointer/ReferenceCount.hh>
//#include <utility/pointer/access_ptr.fwd.hh>
//#include <utility/pointer/access_ptr.hh>
//#include <utility/pointer/owning_ptr.functions.hh>
//#include <utility/pointer/owning_ptr.fwd.hh>
//#include <utility/pointer/owning_ptr.hh>
//#include <utility/signals/BufferedSignalHub.fwd.hh>
//#include <utility/signals/BufferedSignalHub.hh>
//#include <utility/signals/Link.fwd.hh>
//#include <utility/signals/Link.hh>
//#include <utility/signals/LinkUnit.fwd.hh>
//#include <utility/signals/LinkUnit.hh>
//#include <utility/signals/SignalHub.fwd.hh>
//#include <utility/signals/SignalHub.hh>
//#include <numeric/numeric.functions.hh>
//#include <numeric/trig.functions.hh>
//#include <numeric/xyz.functions.fwd.hh>
//#include <numeric/xyzMatrix.fwd.hh>
//#include <numeric/xyzVector.fwd.hh>
//#include <numeric/xyzVector.hh>
//#include <ObjexxFCL/format.hh>
//#include <algorithm>
//#include <cassert>
//#include <cmath>
//#include <complex>
//#include <cstddef>
//#include <cstdlib>
//#include <iomanip>
//#include <iosfwd>
//#include <iostream>
//#include <limits>
//#include <list>
//#include <map>
//#include <ostream>
//#include <sstream>
//#include <string>
//#include <utility>
//#include <vector>
//#include <boost/bind.hpp>
//#include <boost/function.hpp>
//#include <boost/shared_ptr.hpp>


using namespace core;

class ShearMoverBenchmark : public PerformanceBenchmark
{
public:
	pose::PoseOP pose;
	kinematics::MoveMapOP movemap;
	protocols::simple_moves::ShearMover shear_mover;

	//ShearMoverBenchmark(std::string name) : PerformanceBenchmark(name), shear_mover(300., 100) {};
	ShearMoverBenchmark(std::string name) : PerformanceBenchmark(name) {};

	virtual void setUp() {
		pose = new pose::Pose();
		core::import_pose::pose_from_pdb(*pose, "test_in.pdb");

		movemap = new kinematics::MoveMap();
		movemap->set_chi( true );
		movemap->set_bb( true );

		//shear_mover.temperature(300.);
		shear_mover.nmoves(100);
		shear_mover.movemap(movemap);

		shear_mover.angle_max( 'H', 2.0 );
		shear_mover.angle_max( 'E', 2.0 );
		shear_mover.angle_max( 'L', 3.0 );

		// run once to trigger Ramachandran score calculation.
		shear_mover.apply(*pose);
	};

	virtual void run(core::Real scaleFactor) {
		//protocols::simple_moves::ShearMover shear_mover(movemap, 300., 100);
		for(int i=0; i<2000*scaleFactor; i++) {
			shear_mover.apply(*pose);
		}
	};

	virtual void tearDown() {
		pose = 0;
		movemap = 0;
	}
};
