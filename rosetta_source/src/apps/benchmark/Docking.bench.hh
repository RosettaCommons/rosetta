// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/apps/benchmark/Docking.bench.hh
///
/// @brief  Run a performance benchmark of docking
/// @author Monica Berrondo

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <protocols/docking/DockingProtocol.hh>

//#include <basic/Tracer.hh>
//#include <numeric/random/random.hh>


#include "benchmark.hh"
#include "init_util.hh"

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
//#include <core/id/AtomID_Mask.fwd.hh>
//#include <core/id/DOF_ID.fwd.hh>
//#include <core/id/DOF_ID.hh>
//#include <core/id/NamedAtomID.fwd.hh>
//#include <core/id/NamedAtomID.hh>
//#include <core/id/NamedStubID.fwd.hh>
//#include <core/id/NamedStubID.hh>
//#include <core/id/TorsionID.fwd.hh>
//#include <core/id/types.hh>
//#include <core/import_pose/file_data.fwd.hh>
//#include <core/io/pdb/file_data.hh>
//#include <core/kinematics/AtomTree.fwd.hh>
//#include <core/kinematics/FoldTree.fwd.hh>
//#include <core/kinematics/Jump.fwd.hh>
//#include <core/kinematics/MoveMap.fwd.hh>
//#include <core/kinematics/Stub.fwd.hh>
//#include <basic/options/option.hh>
//#include <basic/options/keys/OptionKeys.hh>
//#include <basic/options/keys/docking.OptionKeys.gen.hh>
//#include <core/pack/task/TaskFactory.fwd.hh>
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
//#include <protocols/docking/DockingHighRes.hh>
//#include <protocols/docking/DockingLowRes.hh>
//#include <protocols/docking/DockingProtocol.fwd.hh>
//#include <protocols/filters/Filter.fwd.hh>
//#include <protocols/jobdist/Jobs.fwd.hh>
//#include <protocols/loops/Loops.fwd.hh>
//#include <protocols/moves/DataMap.fwd.hh>
//#include <protocols/moves/MonteCarlo.fwd.hh>
//#include <protocols/moves/MonteCarlo.hh>
//#include <protocols/moves/MonteCarloStatus.hh>
//#include <protocols/moves/Mover.fwd.hh>
//#include <protocols/moves/Mover.hh>
//#include <protocols/moves/MoverContainer.fwd.hh>
//#include <protocols/moves/MoverStatistics.hh>
//#include <protocols/moves/MoverStatus.hh>
//#include <protocols/moves/RigidBodyMover.fwd.hh>
//#include <protocols/moves/mc_convergence_checks/ConvergenceCheck.fwd.hh>
//#include <protocols/moves/mc_convergence_checks/ConvergenceCheck.hh>
//#include <utility/Bound.fwd.hh>
//#include <utility/Bound.hh>
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
//#include <utility/excn/EXCN_Base.hh>
//#include <utility/excn/Exceptions.hh>
//#include <utility/file/FileName.fwd.hh>
//#include <utility/file/FileName.hh>
//#include <utility/file/PathName.fwd.hh>
//#include <utility/file/PathName.hh>
//#include <utility/keys/AutoKey.fwd.hh>
//#include <utility/keys/AutoKey.hh>
//#include <utility/keys/Key.fwd.hh>
//#include <utility/keys/Key.hh>
//#include <utility/keys/KeyLess.fwd.hh>
//#include <utility/keys/KeyLookup.fwd.hh>
//#include <utility/keys/KeyLookup.hh>
//#include <utility/keys/NoClient.fwd.hh>
//#include <utility/keys/NoClient.hh>
//#include <utility/keys/SmallKeyVector.fwd.hh>
//#include <utility/keys/SmallKeyVector.hh>
//#include <utility/keys/UserKey.fwd.hh>
//#include <utility/keys/VariantKey.fwd.hh>
//#include <utility/keys/VariantKey.hh>
//#include <utility/options/AnyOption.fwd.hh>
//#include <utility/options/AnyOption.hh>
//#include <utility/options/AnyVectorOption.fwd.hh>
//#include <utility/options/AnyVectorOption.hh>
//#include <utility/options/BooleanOption.fwd.hh>
//#include <utility/options/BooleanOption.hh>
//#include <utility/options/BooleanVectorOption.fwd.hh>
//#include <utility/options/BooleanVectorOption.hh>
//#include <utility/options/FileOption.fwd.hh>
//#include <utility/options/FileOption.hh>
//#include <utility/options/FileVectorOption.fwd.hh>
//#include <utility/options/FileVectorOption.hh>
//#include <utility/options/IntegerOption.fwd.hh>
//#include <utility/options/IntegerOption.hh>
//#include <utility/options/IntegerVectorOption.fwd.hh>
//#include <utility/options/IntegerVectorOption.hh>
//#include <utility/options/Option.fwd.hh>
//#include <utility/options/Option.hh>
//#include <utility/options/OptionCollection.fwd.hh>
//#include <utility/options/OptionCollection.hh>
//#include <utility/options/PathOption.fwd.hh>
//#include <utility/options/PathOption.hh>
//#include <utility/options/PathVectorOption.fwd.hh>
//#include <utility/options/PathVectorOption.hh>
//#include <utility/options/RealOption.fwd.hh>
//#include <utility/options/RealOption.hh>
//#include <utility/options/RealVectorOption.fwd.hh>
//#include <utility/options/RealVectorOption.hh>
//#include <utility/options/ScalarOption.fwd.hh>
//#include <utility/options/ScalarOption.hh>
//#include <utility/options/ScalarOption_T_.fwd.hh>
//#include <utility/options/ScalarOption_T_.hh>
//#include <utility/options/StringOption.fwd.hh>
//#include <utility/options/StringOption.hh>
//#include <utility/options/StringVectorOption.fwd.hh>
//#include <utility/options/StringVectorOption.hh>
//#include <utility/options/VariantOption.fwd.hh>
//#include <utility/options/VariantOption.hh>
//#include <utility/options/VectorOption.fwd.hh>
//#include <utility/options/VectorOption.hh>
//#include <utility/options/VectorOption_T_.fwd.hh>
//#include <utility/options/VectorOption_T_.hh>
//#include <utility/options/mpi_stderr.hh>
//#include <utility/options/keys/AnyOptionKey.fwd.hh>
//#include <utility/options/keys/AnyOptionKey.hh>
//#include <utility/options/keys/AnyVectorOptionKey.fwd.hh>
//#include <utility/options/keys/AnyVectorOptionKey.hh>
//#include <utility/options/keys/BooleanOptionKey.fwd.hh>
//#include <utility/options/keys/BooleanOptionKey.hh>
//#include <utility/options/keys/BooleanVectorOptionKey.fwd.hh>
//#include <utility/options/keys/BooleanVectorOptionKey.hh>
//#include <utility/options/keys/FileOptionKey.fwd.hh>
//#include <utility/options/keys/FileOptionKey.hh>
//#include <utility/options/keys/FileVectorOptionKey.fwd.hh>
//#include <utility/options/keys/FileVectorOptionKey.hh>
//#include <utility/options/keys/IntegerOptionKey.fwd.hh>
//#include <utility/options/keys/IntegerOptionKey.hh>
//#include <utility/options/keys/IntegerVectorOptionKey.fwd.hh>
//#include <utility/options/keys/IntegerVectorOptionKey.hh>
//#include <utility/options/keys/OptionKey.fwd.hh>
//#include <utility/options/keys/OptionKey.hh>
//#include <utility/options/keys/OptionKeys.hh>
//#include <utility/options/keys/PathOptionKey.fwd.hh>
//#include <utility/options/keys/PathOptionKey.hh>
//#include <utility/options/keys/PathVectorOptionKey.fwd.hh>
//#include <utility/options/keys/PathVectorOptionKey.hh>
//#include <utility/options/keys/RealOptionKey.fwd.hh>
//#include <utility/options/keys/RealOptionKey.hh>
//#include <utility/options/keys/RealVectorOptionKey.fwd.hh>
//#include <utility/options/keys/RealVectorOptionKey.hh>
//#include <utility/options/keys/ScalarOptionKey.fwd.hh>
//#include <utility/options/keys/ScalarOptionKey.hh>
//#include <utility/options/keys/StringOptionKey.fwd.hh>
//#include <utility/options/keys/StringOptionKey.hh>
//#include <utility/options/keys/StringVectorOptionKey.fwd.hh>
//#include <utility/options/keys/StringVectorOptionKey.hh>
//#include <utility/options/keys/VectorOptionKey.fwd.hh>
//#include <utility/options/keys/VectorOptionKey.hh>
//#include <utility/options/keys/all.hh>
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
//#include <ObjexxFCL/TypeTraits.hh>
//#include <ObjexxFCL/char.functions.hh>
//#include <ObjexxFCL/string.functions.hh>
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
//#include <set>
//#include <sstream>
//#include <string>
//#include <utility>
//#include <vector>
//#include <boost/bind.hpp>
//#include <boost/function.hpp>
//#include <boost/shared_ptr.hpp>


enum DockType {High, Low};

template  <DockType dock, int TScale>
class DockingBenchmark : public Benchmark
{
public:
	pose::PoseOP start_pose;
	protocols::docking::DockingProtocolOP docking;

	DockingBenchmark(std::string name) : Benchmark(name) {};

	virtual void setUp() {
		if( dock == Low ) core_init_with_additional_options( "-low_res_protocol_only -dock_pert 3 8 -run:constant_seed" );
		if( dock == High ) core_init_with_additional_options( "-docking_local_refine -run:constant_seed" );
		start_pose = new pose::Pose();
		core::import_pose::pose_from_pdb(*start_pose, "dock_in.pdb");
		docking = new protocols::docking::DockingProtocol();
		docking->set_native_pose(start_pose);
		docking->set_input_pose(start_pose);
	};


	virtual void run(int scaleFactor) {
		//for(int i=0; i<10; i++) {
		//	std::cout << "i="<< i << " R=" << numeric::random::uniform() << std::endl;
		//}

		for(int i=0; i<TScale*scaleFactor; i++) {
			core::pose::Pose pose;
			pose = *start_pose;
			docking->apply( pose );
		}
	};

	virtual void tearDown() {};
};

typedef DockingBenchmark<Low,  10> DockingBenchmark_low;
typedef DockingBenchmark<High, 3> DockingBenchmark_high;
