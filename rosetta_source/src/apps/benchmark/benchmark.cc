// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   rosetta/benchmark/benchmark.cc
///
/// @brief
/// @author Sergey Lyskov

#include <core/chemical/ChemicalManager.hh>
#include <core/scoring/Energies.hh>
#include <basic/Tracer.hh>
#include <devel/init.hh>
#include <core/types.hh>

#include <basic/options/option.hh>

#include <numeric/random/random.hh>

#include <apps/benchmark/benchmark.hh>

#include <time.h>
#include <fstream>

#if  !defined(WINDOWS) && !defined(WIN32)
	#include <sys/time.h>
	#include <sys/resource.h>
#endif


const char results_filename[] = "_performance_";

#include <apps/benchmark/score.bench.hh>
ScoreBenchmark Score_("core.scoring.Score");

#include <apps/benchmark/SmallMover.bench.hh>
SmallMoverBenchmark SmallMover_("protocols.moves.SmallMover");

#include <apps/benchmark/ShearMover.bench.hh>
ShearMoverBenchmark ShearMover_("protocols.moves.ShearMover");

#include <apps/benchmark/Minimizer.bench.hh>
//MinimizerBenchmark Minimizer_("protocols.optimization.Minimizer");
MinimizerBenchmark_dfpmin Minimizer_dfpmin_("protocols.optimization.Minimizer_dfpmin");
MinimizerBenchmark_dfpmin_armijo MinimizerBenchmark_dfpmin_armijo_("protocols.optimization.Minimizer_dfpmin_armijo");
MinimizerBenchmark_dfpmin_armijo_nonmonotone MinimizerBenchmark_dfpmin_armijo_nonmonotone_("protocols.optimization.Minimizer_dfpmin_armijo_nonmonotone");

#include <apps/benchmark/Docking.bench.hh>
DockingBenchmark_low DockingLow("protocols.docking.DockingLowRes");
DockingBenchmark_high DockingHigh("protocols.docking.DockingHighRes");

#include <apps/benchmark/Design.bench.hh>
//DesignBenchmark design("protocols.moves.PackRotamersMover");

#include <apps/benchmark/LigandDock.bench.hh>
//LigandDockBenchmark ligand_dock("protocols.ligand_docking.LigandDockProtocol");

#include <apps/benchmark/LigandDockScript.bench.hh>
//LigandDockScriptBenchmark ligand_dock_script("protocols.ligand_docking.LigandDockScript");

// option key includes

#include <basic/options/keys/james.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
//#include <platform/types.hh>
//#include <core/chemical/AA.hh>
//#include <core/chemical/AtomTypeSet.fwd.hh>
//#include <core/chemical/ChemicalManager.fwd.hh>
//#include <core/chemical/MMAtomTypeSet.fwd.hh>
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
//#include <core/grid/CartGrid.fwd.hh>
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
//#include <core/io/atom_tree_diffs/atom_tree_diff.hh>
//#include <core/import_pose/file_data.fwd.hh>
//#include <core/io/pdb/file_data.hh>
//#include <core/io/pdb/pose_io.hh>
//#include <core/kinematics/AtomTree.fwd.hh>
//#include <core/kinematics/DomainMap.fwd.hh>
//#include <core/kinematics/FoldTree.fwd.hh>
//#include <core/kinematics/Jump.fwd.hh>
//#include <core/kinematics/MoveMap.fwd.hh>
//#include <core/kinematics/MoveMap.hh>
//#include <core/kinematics/Stub.fwd.hh>
//#include <core/kinematics/types.hh>
//#include <core/optimization/AtomTreeMinimizer.fwd.hh>
//#include <core/optimization/AtomTreeMinimizer.hh>
////#include <core/optimization/MinimizerMap.fwd.hh>
//#include <core/optimization/MinimizerOptions.fwd.hh>
//#include <core/optimization/MinimizerOptions.hh>
//#include <core/optimization/types.hh>
//#include <basic/options/keys/OptionKeys.hh>
//#include <basic/options/keys/docking.OptionKeys.gen.hh>

//#include <core/pack/interaction_graph/InteractionGraphBase.fwd.hh>
//#include <core/pack/rotamer_set/RotamerSet.fwd.hh>
//#include <core/pack/rotamer_set/RotamerSets.fwd.hh>
//#include <core/pack/rotamer_set/UnboundRotamersOperation.fwd.hh>
//#include <core/pack/task/PackerTask.fwd.hh>
//#include <core/pack/task/TaskFactory.fwd.hh>
//#include <core/pose/PDBInfo.fwd.hh>
//#include <core/pose/Pose.fwd.hh>
//#include <core/pose/Pose.hh>
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
//#include <core/scoring/ScoreFunction.hh>
//#include <core/scoring/ScoreFunctionFactory.hh>
//#include <core/scoring/ScoreFunctionInfo.fwd.hh>
//#include <core/scoring/ScoreType.hh>
//#include <core/scoring/TenANeighborGraph.fwd.hh>
//#include <core/scoring/TwelveANeighborGraph.fwd.hh>
//#include <core/scoring/types.hh>
//#include <core/scoring/constraints/Constraint.fwd.hh>
//#include <core/scoring/constraints/ConstraintSet.fwd.hh>
//#include <core/scoring/constraints/Constraints.fwd.hh>
//#include <core/scoring/constraints/Constraints.hh>
//#include <core/scoring/constraints/Func.fwd.hh>
//#include <core/scoring/constraints/XYZ_Func.fwd.hh>
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
//#include <core/scoring/symmetry/SymmetricScoreFunction.fwd.hh>
//#include <core/id/SequenceMapping.fwd.hh>
//#include <basic/MetricValue.fwd.hh>
//// AUTO-REMOVED #include <basic/OStream.fwd.hh>
//#include <utility/stream_util.hh>
//#include <basic/Tracer.fwd.hh>
//#include <basic/datacache/BasicDataCache.fwd.hh>
//#include <basic/datacache/BasicDataCache.hh>
//#include <basic/datacache/CacheableData.fwd.hh>
//#include <basic/datacache/CacheableData.hh>
//#include <basic/datacache/DataCache.fwd.hh>
//#include <basic/datacache/DataCache.hh>
//#include <protocols/docking/DockingHighRes.hh>
//#include <protocols/docking/DockingLowRes.hh>
//#include <protocols/docking/DockingProtocol.fwd.hh>
//#include <protocols/docking/DockingProtocol.hh>
//#include <protocols/enzdes/EnzConstraintIO.fwd.hh>
//#include <protocols/enzdes/EnzConstraintIO.hh>
//#include <protocols/enzdes/EnzCstTemplateRes.hh>
//#include <protocols/enzdes/MatchConstraintFileInfo.fwd.hh>
//#include <protocols/enzdes/MatchConstraintFileInfo.hh>
//#include <protocols/filters/Filter.fwd.hh>
//#include <protocols/ligand_docking/LigandBaseProtocol.hh>
//#include <protocols/ligand_docking/LigandDockProtocol.fwd.hh>
//#include <protocols/ligand_docking/LigandDockProtocol.hh>
//#include <protocols/ligand_docking/ResidueTorsionRestraints.fwd.hh>
//#include <protocols/loops/Loops.fwd.hh>
//#include <protocols/match/downstream/ExternalGeomSampler.fwd.hh>
//#include <protocols/moves/BackboneMover.fwd.hh>
//#include <protocols/moves/BackboneMover.hh>
//#include <protocols/moves/DataMap.fwd.hh>
//#include <protocols/moves/MonteCarlo.fwd.hh>
//#include <protocols/moves/MonteCarlo.hh>
//#include <protocols/moves/MonteCarloStatus.hh>
//#include <protocols/moves/Mover.fwd.hh>
//#include <protocols/moves/Mover.hh>
//#include <protocols/moves/MoverContainer.fwd.hh>
//#include <protocols/moves/MoverStatistics.hh>
//#include <protocols/moves/MoverStatus.hh>
//#include <protocols/moves/PackRotamersMover.fwd.hh>
//#include <protocols/moves/PackRotamersMover.hh>
//#include <protocols/moves/RigidBodyMover.fwd.hh>
//#include <protocols/moves/mc_convergence_checks/ConvergenceCheck.fwd.hh>
//#include <protocols/moves/mc_convergence_checks/ConvergenceCheck.hh>
//#include <utility/Bound.fwd.hh>
//#include <utility/Bound.hh>
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
//#include <utility/tag/Tag.fwd.hh>
//#include <utility/excn/EXCN_Base.hh>
//#include <utility/excn/Exceptions.hh>
//#include <utility/file/FileName.fwd.hh>
//#include <utility/file/FileName.hh>
//#include <utility/file/PathName.fwd.hh>
//#include <utility/file/PathName.hh>
//#include <utility/io/izstream.fwd.hh>
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
//#include <numeric/random/random.fwd.hh>
//#include <numeric/random/uniform.hh>
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
//#include <boost/config.hpp>
//#include <boost/function.hpp>
//#include <boost/pool/detail/mutex.hpp>
//#include <boost/pool/poolfwd.hpp>
//#include <boost/shared_ptr.hpp>




using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("benchmark");

using namespace core;

std::vector<Benchmark *> &Benchmark::allBenchmarks()
{
	static std::vector<Benchmark*> * allBenchmarks = new std::vector<Benchmark*>;
	return *allBenchmarks;
}

double Benchmark::execute(int scaleFactor)
{
	/// Reseting RG system before each performance run.
	numeric::random::RandomGenerator::initializeRandomGenerators(
		 1000, numeric::random::_RND_TestRun_, "mt19937");

	TR << "Setting up "<< name() << "..." << std::endl;
	//for(int i=0; i<3; i++) {
	//	std::cout << "X i="<< i << " R=" << numeric::random::uniform() << std::endl;
	//}

	setUp();

	double t;

	#if  !defined(WINDOWS) && !defined(WIN32)
		TR << "Running(U) " << name() << "..." << std::endl;
		struct rusage R0, R1;

		getrusage(RUSAGE_SELF, &R0);
		run(scaleFactor);

		getrusage(RUSAGE_SELF, &R1);

		t = R1.ru_utime.tv_sec + R1.ru_utime.tv_usec*1e-6 - R0.ru_utime.tv_sec - R0.ru_utime.tv_usec*1e-6;
		TR << "Running(U) " << name() << "... Done. Time:" << t << std::endl;
	#else
		TR << "Running(W) " << name() << "..." << std::endl;
		t = clock();
		run(scaleFactor);
		t = clock() - t;
		t = t / CLOCKS_PER_SEC;
		TR << "Running(W) " << name() << "... Done. Time:" << t << std::endl;
	#endif


	TR << "Tear down "<< name() << "..." << std::endl;
	tearDown();
	TR << "Tear down "<< name() << "... Done." << std::endl << std::endl;

	result_ = t;
	return t;
}

void Benchmark::executeAllBenchmarks(int scaleFactor)
{
	TR << std::endl << "Executing all benchmarks..." << std::endl << std::endl;

	std::vector<Benchmark *> & all( allBenchmarks() );

	for(Size i=0; i<all.size(); i++) {
		Benchmark * B = all[i];
		B->execute(scaleFactor);
	}
	TR << std::endl << "Executing all benchmarks... Done." << std::endl;
}

///
/// Generting report file in python dict format: i.e: { 'Bench1': 1.5, 'Bench2': 1.6 }
///
std::string Benchmark::getReport()
{
	std::vector<Benchmark *> & all( allBenchmarks() );

	char buf[1024];

	std::string res = "{\n";
	for(Size i=0; i<all.size(); i++) {
		Benchmark * B = all[i];
		sprintf(buf, "%f", B->result_);
		res += "    '" + B->name_ + "':" + std::string(buf) + ",\n";
	}
	res += "}\n";
	return res;
}

int real_command_line_argc; char ** real_command_line_argv;
int command_line_argc; char ** command_line_argv;


int main( int argc, char *argv[])
{
	command_line_argc=argc; command_line_argv=argv;
	real_command_line_argc=argc; real_command_line_argv=argv;

	using namespace core;
	using namespace basic::options::OptionKeys;

	basic::options::option.add_relevant(run::benchmark_scale);
	basic::options::option.add_relevant(in::path::database);

	devel::init(argc, argv);

	//TR << "DB:"  << basic::options::option[ in::path::database ]() << "\n";

	chemical::ResidueTypeSetCAP residue_set
		( chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD ) );

	//TR << "Specified()=" << basic::options::option[ run::benchmark_scale ].user() << "\n";
	//TR << "Legal" << basic::options::option[ run::benchmark_scale ].legal() << "\n";
	//TR << "Active:" << basic::options::option[ run::benchmark_scale ].user() << "\n";
	//TR << "native:"  << basic::options::option[ james::native ]() << "\n";
	//TR << "DB:"  << basic::options::option[ in::path::database ]() << "\n";
	int scale = basic::options::option[ run::benchmark_scale ]();

	TR << "Mini Benchmark started! Scale factor: " << scale << " -------------" << std::endl;

	//TR << "CLOCKS_PER_SEC:" << CLOCKS_PER_SEC << "\n";

	Benchmark::executeAllBenchmarks(scale);
	std::string report = Benchmark::getReport();
	TR << "Results:" << std::endl << report;  TR.flush();

	/// Now, saving report to a file
	std::ofstream file(results_filename, std::ios::out | std::ios::binary);
	if(!file) {
		Error() << "Benchmark:: Unable to open file:" << results_filename << " for writing!!!" << std::endl;
		return 1;
	}
	file << report;

	file.close();

	TR << "Mini Benchmark ended.   --------------------------------" << std::endl;
	return 0;
}
