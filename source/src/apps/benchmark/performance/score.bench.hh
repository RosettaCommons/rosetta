// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   rosetta/benchmark/score.bench.cc
///
/// @brief  Scoring benchmark
/// @author Sergey Lyskov

#include <apps/benchmark/performance/performance_benchmark.hh>

#include <core/pose/Pose.hh>
//#include <core/import_pose/import_pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>

#include <utility/vector1.hh>

//Auto Headers
//#include <platform/types.hh>
//#include <core/types.hh>
//#include <core/chemical/AA.hh>
//#include <core/chemical/ResidueType.fwd.hh>
//#include <core/chemical/ResidueTypeSet.fwd.hh>
//#include <core/conformation/Atom.fwd.hh>
//#include <core/conformation/Conformation.fwd.hh>
//#include <core/conformation/Residue.fwd.hh>
//#include <core/conformation/signals/XYZEvent.fwd.hh>
//#include <core/graph/ArrayPool.hh>
//#include <core/graph/Graph.fwd.hh>
//#include <core/graph/Graph.hh>
//#include <core/conformation/PointGraph.fwd.hh>
//#include <core/conformation/PointGraphData.fwd.hh>
//#include <core/graph/UpperEdgeGraph.fwd.hh>
//#include <core/graph/unordered_object_pool.fwd.hpp>
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
//#include <core/kinematics/DomainMap.fwd.hh>
//#include <core/kinematics/FoldTree.fwd.hh>
//#include <core/kinematics/Jump.fwd.hh>
//#include <core/kinematics/Stub.fwd.hh>
////#include <core/optimization/MinimizerMap.fwd.hh>

//#include <core/pack/rotamer_set/RotamerSet.fwd.hh>
//#include <core/pack/task/PackerTask.fwd.hh>
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
//#include <core/scoring/ContextGraph.fwd.hh>
//#include <core/scoring/ContextGraphTypes.hh>
//#include <core/scoring/Energies.fwd.hh>
//#include <core/scoring/EnergiesCacheableDataType.hh>
//#include <core/scoring/EnergyGraph.fwd.hh>
//#include <core/scoring/EnergyGraph.hh>
//#include <core/scoring/EnergyMap.fwd.hh>
//#include <core/scoring/EnergyMap.hh>
//#include <core/scoring/LREnergyContainer.fwd.hh>
//#include <core/scoring/NeighborList.fwd.hh>
//#include <core/scoring/ScoreFunction.fwd.hh>
//#include <core/scoring/ScoreFunctionInfo.fwd.hh>
//#include <core/scoring/ScoreType.hh>
//#include <core/scoring/TenANeighborGraph.fwd.hh>
//#include <core/scoring/TwelveANeighborGraph.fwd.hh>
//#include <core/scoring/types.hh>
//#include <core/scoring/constraints/Constraint.fwd.hh>
//#include <core/scoring/constraints/ConstraintSet.fwd.hh>
//#include <core/scoring/constraints/Constraints.fwd.hh>
//#include <core/scoring/methods/ContextDependentLRTwoBodyEnergy.fwd.hh>
//#include <core/scoring/methods/ContextDependentOneBodyEnergy.fwd.hh>
//#include <core/scoring/methods/ContextDependentTwoBodyEnergy.fwd.hh>
//#include <core/scoring/methods/ContextIndependentLRTwoBodyEnergy.fwd.hh>
//#include <core/scoring/methods/ContextIndependentOneBodyEnergy.fwd.hh>
//#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.fwd.hh>
//#include <core/scoring/methods/EnergyMethod.fwd.hh>
//#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
//#include <core/scoring/methods/LongRangeTwoBodyEnergy.fwd.hh>
//#include <core/scoring/methods/Methods.hh>
//#include <core/scoring/methods/TwoBodyEnergy.fwd.hh>
//#include <core/scoring/methods/WholeStructureEnergy.fwd.hh>
//#include <basic/MetricValue.fwd.hh>
//// AUTO-REMOVED #include <basic/OStream.fwd.hh>
//#include <basic/Tracer.fwd.hh>
//#include <basic/datacache/BasicDataCache.fwd.hh>
//#include <basic/datacache/BasicDataCache.hh>
//#include <basic/datacache/CacheableData.fwd.hh>
//#include <basic/datacache/CacheableData.hh>
//#include <basic/datacache/DataCache.fwd.hh>
//#include <basic/datacache/DataCache.hh>
//#include <utility/down_cast.hh>
//#include <utility/exit.hh>
//#include <utility/string_util.hh>
//#include <utility/vector1.fwd.hh>
//#include <utility/vector1.hh>
//#include <utility/vector1_bool.hh>
//#include <utility/vectorL.fwd.hh>
//#include <utility/vectorL.hh>
//#include <utility/vectorL_Selector.hh>
//#include <utility/vectorL_bool.hh>
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
//#include <ObjexxFCL/CArray.fwd.hh>
//#include <ObjexxFCL/CArrayP.fwd.hh>
//#include <ObjexxFCL/Dimension.fwd.hh>
//#include <ObjexxFCL/Dimension.hh>
//#include <ObjexxFCL/DimensionExpression.hh>
//#include <ObjexxFCL/DynamicIndexRange.fwd.hh>
//#include <ObjexxFCL/DynamicIndexRange.hh>
//#include <ObjexxFCL/FArray.all.fwd.hh>
//#include <ObjexxFCL/FArray1D.fwd.hh>
//#include <ObjexxFCL/FArray1D.hh>
//#include <ObjexxFCL/FArray1.fwd.hh>
//#include <ObjexxFCL/FArray1.hh>
//#include <ObjexxFCL/FArray1A.fwd.hh>
//#include <ObjexxFCL/FArray1P.fwd.hh>
//#include <ObjexxFCL/FArray1.all.fwd.hh>
//#include <ObjexxFCL/FArray2D.fwd.hh>
//#include <ObjexxFCL/FArray2D.hh>
//#include <ObjexxFCL/FArray2.fwd.hh>
//#include <ObjexxFCL/FArray2.hh>
//#include <ObjexxFCL/FArray2A.fwd.hh>
//#include <ObjexxFCL/FArray2P.fwd.hh>
//#include <ObjexxFCL/FArray2.all.fwd.hh>
//#include <ObjexxFCL/FArray3D.fwd.hh>
//#include <ObjexxFCL/FArray3.fwd.hh>
//#include <ObjexxFCL/FArray3A.fwd.hh>
//#include <ObjexxFCL/FArray3P.fwd.hh>
//#include <ObjexxFCL/FArray3.all.fwd.hh>
//#include <ObjexxFCL/FArray4D.fwd.hh>
//#include <ObjexxFCL/FArray4.fwd.hh>
//#include <ObjexxFCL/FArray4A.fwd.hh>
//#include <ObjexxFCL/FArray4P.fwd.hh>
//#include <ObjexxFCL/FArray4.all.fwd.hh>
//#include <ObjexxFCL/FArray5D.fwd.hh>
//#include <ObjexxFCL/FArray5.fwd.hh>
//#include <ObjexxFCL/FArray5A.fwd.hh>
//#include <ObjexxFCL/FArray5P.fwd.hh>
//#include <ObjexxFCL/FArray5.all.fwd.hh>
//#include <ObjexxFCL/FArray.all.fwd.hh>
//#include <ObjexxFCL/FArray.hh>
//#include <ObjexxFCL/FArrayInitializer.fwd.hh>
//#include <ObjexxFCL/FArrayInitializer.hh>
//#include <ObjexxFCL/FArraySection.fwd.hh>
//#include <ObjexxFCL/FArraySection.hh>
//#include <ObjexxFCL/FArrayTraits.fwd.hh>
//#include <ObjexxFCL/FArrayTraits.hh>
//#include <ObjexxFCL/IndexRange.fwd.hh>
//#include <ObjexxFCL/IndexRange.hh>
//#include <ObjexxFCL/KeyFArray1D.fwd.hh>
//#include <ObjexxFCL/KeyFArray2D.fwd.hh>
//#include <ObjexxFCL/KeyFArray3D.fwd.hh>
//#include <ObjexxFCL/Observer.fwd.hh>
//#include <ObjexxFCL/Observer.hh>
//#include <ObjexxFCL/ObserverMulti.hh>
//#include <ObjexxFCL/ObserverSingle.hh>
//#include <ObjexxFCL/SetWrapper.fwd.hh>
//#include <ObjexxFCL/Star.fwd.hh>
//#include <ObjexxFCL/Star.hh>
//#include <algorithm>
//#include <cassert>
//#include <cmath>
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
//#include <vector>
//#include <boost/bind.hpp>
//#include <boost/config.hpp>
//#include <boost/function.hpp>
//#include <boost/pool/detail/mutex.hpp>
//#include <boost/pool/poolfwd.hpp>


using namespace core;

class ScoreBenchmark : public PerformanceBenchmark
{
public:
	pose::Pose pose;

	ScoreBenchmark(std::string name) : PerformanceBenchmark(name) {};

	virtual void setUp() {
		core::import_pose::pose_from_pdb(pose, "test_in.pdb");
	}

	virtual void run(core::Real scaleFactor) {
		scoring::ScoreFunction scorefxn;
		for(int i=0; i<45000*scaleFactor; i++) {
			scorefxn(pose);
			pose.energies().clear();
		}
	};

	virtual void tearDown() {};
};
