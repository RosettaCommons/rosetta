// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   rosetta/benchmark/Minimizer.bench.cc
///
/// @brief  Varios moves benchmark
/// @author Sergey Lyskov


#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include "benchmark.hh"

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
//#include <core/kinematics/DomainMap.fwd.hh>
//#include <core/kinematics/FoldTree.fwd.hh>
//#include <core/kinematics/Jump.fwd.hh>
//#include <core/kinematics/MoveMap.fwd.hh>
//#include <core/kinematics/Stub.fwd.hh>
//#include <core/kinematics/types.hh>
//#include <core/optimization/AtomTreeMinimizer.fwd.hh>
////#include <core/optimization/MinimizerMap.fwd.hh>
//#include <core/optimization/MinimizerOptions.fwd.hh>
//#include <core/optimization/types.hh>

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
//#include <core/scoring/Energies.fwd.hh>
//#include <core/scoring/EnergyMap.fwd.hh>
//#include <core/scoring/EnergyMap.hh>
//#include <core/scoring/ScoreFunction.fwd.hh>
//#include <core/scoring/ScoreFunctionInfo.fwd.hh>
//#include <core/scoring/ScoreType.hh>
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
//#include <core/scoring/methods/TwoBodyEnergy.fwd.hh>
//#include <core/scoring/methods/WholeStructureEnergy.fwd.hh>
//#include <core/scoring/symmetry/SymmetricScoreFunction.fwd.hh>
//#include <basic/MetricValue.fwd.hh>
//// AUTO-REMOVED #include <basic/OStream.fwd.hh>
//#include <basic/Tracer.fwd.hh>
//#include <basic/datacache/BasicDataCache.fwd.hh>
//#include <utility/down_cast.hh>
//#include <utility/exit.hh>
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
//#include <ObjexxFCL/FArray.all.fwd.hh>
//#include <ObjexxFCL/FArray1D.fwd.hh>
//#include <ObjexxFCL/FArray1.fwd.hh>
//#include <ObjexxFCL/FArray1A.fwd.hh>
//#include <ObjexxFCL/FArray1P.fwd.hh>
//#include <ObjexxFCL/FArray1.all.fwd.hh>
//#include <ObjexxFCL/FArray2D.fwd.hh>
//#include <ObjexxFCL/FArray2.fwd.hh>
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
//#include <ObjexxFCL/KeyFArray1D.fwd.hh>
//#include <ObjexxFCL/KeyFArray2D.fwd.hh>
//#include <ObjexxFCL/KeyFArray3D.fwd.hh>
//#include <algorithm>
//#include <cassert>
//#include <cmath>
//#include <cstddef>
//#include <cstdlib>
//#include <iomanip>
//#include <iosfwd>
//#include <iostream>
//#include <limits>
//#include <map>
//#include <ostream>
//#include <sstream>
//#include <string>
//#include <utility>
//#include <vector>
//#include <boost/bind.hpp>
//#include <boost/function.hpp>


using namespace core;

enum  ScoreFnType {SFT_dfpmin, SFT_dfpmin_armijo, SFT_dfpmin_armijo_nonmonotone};

template  <ScoreFnType sft, int TScale>
class MinimizerBenchmark : public Benchmark
{
public:
	pose::PoseOP start_pose;
	kinematics::MoveMap mm;
	core::scoring::ScoreFunctionOP scorefxn;
	core::optimization::AtomTreeMinimizer minimizer;

	MinimizerBenchmark(std::string name) : Benchmark(name) {};

	virtual void setUp() {
		start_pose = new pose::Pose();
		core::import_pose::pose_from_pdb(*start_pose, "test_in.pdb");

		scorefxn = core::scoring::ScoreFunctionFactory::create_score_function( "standard" );

				//kinematics::MoveMap mm;
				for ( int i=30; i<= 35; ++i ) {
						mm.set_bb ( i, true );
						mm.set_chi( i, true );
				}

		(*scorefxn)( *start_pose ); // to triger dunbrack loading/calcualtion
	};

	virtual void run(core::Real scaleFactor) {
		for(int i=0; i<TScale*scaleFactor; i++) {
			std::string stype = "unknow";
			if( sft == SFT_dfpmin ) stype = "dfpmin";
			if( sft == SFT_dfpmin_armijo ) stype = "dfpmin_armijo";
			if( sft == SFT_dfpmin_armijo_nonmonotone ) stype = "dfpmin_armijo_nonmonotone";
			core::optimization::MinimizerOptions options( stype/*"dfpmin"*/, 1e-1, true, true );
			core::pose::Pose pose;
			pose = *start_pose;
			minimizer.run( pose, mm, *scorefxn, options );
		}
	};

	virtual void tearDown() {};
};

typedef MinimizerBenchmark<SFT_dfpmin, 38> MinimizerBenchmark_dfpmin;
typedef MinimizerBenchmark<SFT_dfpmin_armijo, 4> MinimizerBenchmark_dfpmin_armijo;
typedef MinimizerBenchmark<SFT_dfpmin_armijo_nonmonotone, 4> MinimizerBenchmark_dfpmin_armijo_nonmonotone;


//class MinimizerBenchmark_dfpmin : public MinimizerBenchmark

/*
				{ // armijo
						MinimizerOptions options( "dfpmin_armijo", 1e-1, true, true );
						Pose pose;
						pose = start_pose;
						minimizer.run( pose, mm, *scorefxn, options );
						core::Real score = (*scorefxn)( pose );
						TR << "dfpmin_armijo/standard: " << score << "\n";
						TS_ASSERT_DELTA(score, 385.767, err_tol);
				}

				{ // non-monotone
						MinimizerOptions options( "dfpmin_armijo_nonmonotone", 1e-1, true, true );
						Pose pose;
						pose = start_pose;
						minimizer.run( pose, mm, *scorefxn, options );
						core::Real score = (*scorefxn)( pose );
						TR << "dfpmin_armijo_nonmonotone/standard: " << score << "\n";
						TS_ASSERT_DELTA(score, 385.945, err_tol);
				} */

