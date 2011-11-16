// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   protocols/jd2/MPIFileBufJobDistributor.cc
/// @brief  implementation of MPIFileBufJobDistributor
/// @author Oliver Lange olange@u.washington.edu


// MPI headers
#ifdef USEMPI
#include <mpi.h> //keep this first
#endif

// Unit headers
#include <protocols/jd2/archive/ArchiveManager.hh>
#include <protocols/jd2/archive/MPIArchiveJobDistributor.hh>
#include <protocols/jd2/archive/ArchiveBase.hh>
#include <protocols/jd2/BatchJobInputter.hh> //for BOGUS_BATCH_ID

// Package headers
#include <protocols/jd2/MpiFileBuffer.hh>

// to test setup-file
#include <protocols/topology_broker/TopologyBroker.hh>
#include <protocols/topology_broker/util.hh>

//for factory
#include <protocols/abinitio/IterativeAbrelax.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <basic/MemTracer.hh>
#include <basic/options/option.hh>
#include <utility/exit.hh>

#include <core/io/silent/SilentFileData.hh>

#include <ObjexxFCL/string.functions.hh>
#include <utility/file/file_sys_util.hh>

// Option headers
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/jd2.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/broker.OptionKeys.gen.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>

//#include <basic/options/keys/jd2.OptionKeys.gen.hh>
#include <basic/options/option.cc.gen.hh>
#include <basic/options/option_macros.hh>
#include <basic/prof.hh>

// C++ headers
#include <string>
// AUTO-REMOVED #include <ctime>
#include <sstream>
#include <iterator>

//Debug headers
// AUTO-REMOVED #include <protocols/abinitio/AbrelaxMover.hh>
#include <fstream> //testing
#include <utility/io/izstream.hh>

// Auto-header: duplicate removed #include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

#if (defined WIN32) //&& (!defined WIN_PYROSETTA)
// AUTO-REMOVED #include <windows.h>
#endif

//Auto Headers
#include <platform/types.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
// AUTO-REMOVED #include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
// AUTO-REMOVED #include <core/conformation/Conformation.fwd.hh>
// AUTO-REMOVED #include <core/conformation/Residue.fwd.hh>
// AUTO-REMOVED #include <core/conformation/RotamerSetBase.fwd.hh>
// AUTO-REMOVED #include <core/conformation/signals/ConnectionEvent.fwd.hh>
// AUTO-REMOVED #include <core/conformation/signals/LengthEvent.fwd.hh>
// AUTO-REMOVED #include <core/conformation/signals/XYZEvent.fwd.hh>
// AUTO-REMOVED #include <core/fragment/BaseCacheUnit.hh>
// AUTO-REMOVED #include <core/fragment/FragData.fwd.hh>
// AUTO-REMOVED #include <core/fragment/FragData.hh>
// AUTO-REMOVED #include <core/fragment/FragID.fwd.hh>
// AUTO-REMOVED #include <core/fragment/FragID_Iterator.fwd.hh>
#include <core/fragment/FragSet.fwd.hh>
// AUTO-REMOVED #include <core/fragment/FragSet.hh>
// AUTO-REMOVED #include <core/fragment/Frame.fwd.hh>
// AUTO-REMOVED #include <core/fragment/Frame.hh>
// AUTO-REMOVED #include <core/fragment/FrameIterator.fwd.hh>
// AUTO-REMOVED #include <core/fragment/FrameIterator.hh>
// AUTO-REMOVED #include <core/fragment/FrameIteratorWorker_.fwd.hh>
// AUTO-REMOVED #include <core/fragment/FrameIteratorWorker_.hh>
// AUTO-REMOVED #include <core/fragment/FrameList.fwd.hh>
// AUTO-REMOVED #include <core/fragment/FrameList.hh>
// AUTO-REMOVED #include <core/fragment/SingleResidueFragData.fwd.hh>
// AUTO-REMOVED #include <core/fragment/SingleResidueFragData.hh>
// AUTO-REMOVED #include <core/id/AtomID.fwd.hh>
// AUTO-REMOVED #include <core/id/AtomID.hh>
// AUTO-REMOVED #include <core/id/AtomID_Map.fwd.hh>
// AUTO-REMOVED #include <core/id/AtomID_Mask.fwd.hh>
// AUTO-REMOVED #include <core/id/DOF_ID.fwd.hh>
// AUTO-REMOVED #include <core/id/DOF_ID.hh>
// AUTO-REMOVED #include <core/id/DOF_ID_Map.fwd.hh>
// AUTO-REMOVED #include <core/id/DOF_ID_Mask.fwd.hh>
// AUTO-REMOVED #include <core/id/JumpID.fwd.hh>
// AUTO-REMOVED #include <core/id/JumpID.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedAtomID.hh>
// AUTO-REMOVED #include <core/id/NamedStubID.fwd.hh>
// AUTO-REMOVED #include <core/id/SequenceMapping.fwd.hh>
// AUTO-REMOVED #include <core/id/SequenceMapping.hh>
// AUTO-REMOVED #include <core/id/TorsionID.fwd.hh>
// AUTO-REMOVED #include <core/id/TorsionID.hh>
#include <core/id/types.hh>
#include <core/io/silent/SharedSilentData.hh>
#include <core/io/silent/SilentEnergy.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/io/silent/SilentStruct.fwd.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/silent.fwd.hh>
// AUTO-REMOVED #include <core/kinematics/AtomTree.fwd.hh>
// AUTO-REMOVED #include <core/kinematics/DomainMap.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
// AUTO-REMOVED #include <core/kinematics/Jump.fwd.hh>
// AUTO-REMOVED #include <core/kinematics/MinimizerMapBase.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
// AUTO-REMOVED #include <core/kinematics/MoveMap.hh>
#include <core/kinematics/ShortestPathInFoldTree.fwd.hh>
// AUTO-REMOVED #include <core/kinematics/Stub.fwd.hh>
// AUTO-REMOVED #include <core/kinematics/types.hh>
// AUTO-REMOVED #include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/Pose.fwd.hh>
// AUTO-REMOVED #include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/pose/datacache/ObserverCache.fwd.hh>
// AUTO-REMOVED #include <core/pose/metrics/PoseMetricContainer.fwd.hh>
// AUTO-REMOVED #include <core/pose/signals/ConformationEvent.fwd.hh>
// AUTO-REMOVED #include <core/pose/signals/DestructionEvent.fwd.hh>
// AUTO-REMOVED #include <core/pose/signals/EnergyEvent.fwd.hh>
// AUTO-REMOVED #include <core/pose/signals/GeneralEvent.fwd.hh>
// AUTO-REMOVED #include <core/scoring/Energies.fwd.hh>
// AUTO-REMOVED #include <core/scoring/EnergyGraph.fwd.hh>
// AUTO-REMOVED #include <core/scoring/EnergyMap.fwd.hh>
// AUTO-REMOVED #include <core/scoring/EnergyMap.hh>
// AUTO-REMOVED #include <core/scoring/LREnergyContainer.fwd.hh>
// AUTO-REMOVED #include <core/scoring/MinimizationData.fwd.hh>
// AUTO-REMOVED #include <core/scoring/MinimizationGraph.fwd.hh>
#include <core/scoring/ResidualDipolarCoupling.fwd.hh>
// AUTO-REMOVED #include <core/scoring/ResidualDipolarCoupling.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunction.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunctionInfo.fwd.hh>
// AUTO-REMOVED #include <core/scoring/ScoreType.hh>
// AUTO-REMOVED #include <core/scoring/constraints/AmbiguousNMRConstraint.fwd.hh>
// AUTO-REMOVED #include <core/scoring/constraints/AmbiguousNMRDistanceConstraint.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
// AUTO-REMOVED #include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
// AUTO-REMOVED #include <core/scoring/constraints/ConstraintSet.hh>
// AUTO-REMOVED #include <core/scoring/constraints/Constraints.fwd.hh>
// AUTO-REMOVED #include <core/scoring/constraints/Constraints.hh>
// AUTO-REMOVED #include <core/scoring/constraints/DOF_Constraint.fwd.hh>
// AUTO-REMOVED #include <core/scoring/constraints/DOF_Constraint.hh>
// AUTO-REMOVED #include <core/scoring/constraints/Func.fwd.hh>
// AUTO-REMOVED #include <core/scoring/constraints/Func.hh>
// AUTO-REMOVED #include <core/scoring/constraints/FuncFactory.fwd.hh>
// AUTO-REMOVED #include <core/scoring/constraints/HarmonicFunc.fwd.hh>
// AUTO-REMOVED #include <core/scoring/constraints/HarmonicFunc.hh>
// AUTO-REMOVED #include <core/scoring/constraints/XYZ_Func.fwd.hh>
// AUTO-REMOVED #include <core/scoring/methods/ContextDependentLRTwoBodyEnergy.fwd.hh>
// AUTO-REMOVED #include <core/scoring/methods/ContextDependentOneBodyEnergy.fwd.hh>
// AUTO-REMOVED #include <core/scoring/methods/ContextDependentTwoBodyEnergy.fwd.hh>
// AUTO-REMOVED #include <core/scoring/methods/ContextIndependentLRTwoBodyEnergy.fwd.hh>
// AUTO-REMOVED #include <core/scoring/methods/ContextIndependentOneBodyEnergy.fwd.hh>
// AUTO-REMOVED #include <core/scoring/methods/ContextIndependentTwoBodyEnergy.fwd.hh>
// AUTO-REMOVED #include <core/scoring/methods/EnergyMethod.fwd.hh>
// AUTO-REMOVED #include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
// AUTO-REMOVED #include <core/scoring/methods/LongRangeTwoBodyEnergy.fwd.hh>
// AUTO-REMOVED #include <core/scoring/methods/TwoBodyEnergy.fwd.hh>
// AUTO-REMOVED #include <core/scoring/methods/WholeStructureEnergy.fwd.hh>
// AUTO-REMOVED #include <protocols/abinitio/AbrelaxMover.fwd.hh>
#include <protocols/abinitio/FragmentSampler.fwd.hh>
// AUTO-REMOVED #include <protocols/abinitio/FragmentSampler.hh>
#include <protocols/abinitio/IterativeBase.hh>
#include <protocols/abinitio/IterativeCentroid.hh>
#include <protocols/abinitio/IterativeFullatom.hh>
#include <protocols/abinitio/PairingStatistics.fwd.hh>
#include <protocols/checkpoint/CheckPointer.fwd.hh>
#include <protocols/checkpoint/CheckPointer.hh>
#include <protocols/evaluation/PoseEvaluator.fwd.hh>
// AUTO-REMOVED #include <protocols/evaluation/PoseEvaluator.hh>
// AUTO-REMOVED #include <protocols/filters/Filter.fwd.hh>
#include <protocols/jd2/BatchJobInputter.fwd.hh>
#include <protocols/jd2/InnerJob.fwd.hh>
#include <protocols/jd2/Job.fwd.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.fwd.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobDistributorFactory.fwd.hh>
#include <protocols/jd2/JobInputter.fwd.hh>
#include <protocols/jd2/JobInputter.hh>
#include <protocols/jd2/JobOutputter.fwd.hh>
#include <protocols/jd2/MPIFileBufJobDistributor.fwd.hh>
#include <protocols/jd2/MPIFileBufJobDistributor.hh>
#include <protocols/jd2/MpiFileBuffer.fwd.hh>
#include <protocols/jd2/Parser.fwd.hh>
#include <protocols/jd2/SingleFileBuffer.fwd.hh>
#include <protocols/jd2/archive/ArchiveBase.fwd.hh>
#include <protocols/jd2/archive/ArchiveManager.fwd.hh>
#include <protocols/jd2/archive/EvaluatedArchive.hh>
#include <protocols/jd2/archive/MPIArchiveJobDistributor.fwd.hh>
#include <protocols/jobdist/Jobs.fwd.hh>
#include <protocols/jobdist/Jobs.hh>
#include <protocols/loops/Loop.fwd.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.fwd.hh>
#include <protocols/loops/Loops.hh>
// AUTO-REMOVED #include <protocols/loops/SlidingWindowLoopClosure.fwd.hh>
// AUTO-REMOVED #include <protocols/moves/DataMap.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
// AUTO-REMOVED #include <protocols/moves/MonteCarlo.hh>
// AUTO-REMOVED #include <protocols/moves/MonteCarloStatus.hh>
#include <protocols/moves/Mover.fwd.hh>
// AUTO-REMOVED #include <protocols/moves/Mover.hh>
// AUTO-REMOVED #include <protocols/moves/MoverContainer.fwd.hh>
// AUTO-REMOVED #include <protocols/moves/MoverStatus.hh>
// AUTO-REMOVED #include <protocols/moves/mc_convergence_checks/ConvergenceCheck.fwd.hh>
// AUTO-REMOVED #include <protocols/moves/mc_convergence_checks/ConvergenceCheck.hh>
#include <protocols/noesy_assign/CrossPeak.fwd.hh>
#include <protocols/noesy_assign/CrossPeak.hh>
#include <protocols/noesy_assign/CrossPeakInfo.fwd.hh>
#include <protocols/noesy_assign/CrossPeakInfo.hh>
// AUTO-REMOVED #include <protocols/noesy_assign/CrossPeakList.fwd.hh>
#include <protocols/noesy_assign/NoesyModule.fwd.hh>
// AUTO-REMOVED #include <protocols/noesy_assign/NoesyModule.hh>
#include <protocols/noesy_assign/PeakAssignment.fwd.hh>
// AUTO-REMOVED #include <protocols/noesy_assign/PeakAssignment.hh>
// AUTO-REMOVED #include <protocols/noesy_assign/PeakAssignmentParameters.hh>
// AUTO-REMOVED #include <protocols/noesy_assign/PeakAssignmentResidueMap.fwd.hh>
// AUTO-REMOVED #include <protocols/noesy_assign/PeakAssignmentResidueMap.hh>
#include <protocols/noesy_assign/PeakCalibrator.fwd.hh>
// AUTO-REMOVED #include <protocols/noesy_assign/PeakFileFormat.fwd.hh>
#include <protocols/noesy_assign/Resonance.hh>
#include <protocols/noesy_assign/ResonanceList.fwd.hh>
#include <protocols/noesy_assign/ResonanceList.hh>
// AUTO-REMOVED #include <protocols/relax/RelaxProtocolBase.fwd.hh>
#include <protocols/topology_broker/ClaimerMessage.fwd.hh>
#include <protocols/topology_broker/DofClaim.fwd.hh>
// AUTO-REMOVED #include <protocols/topology_broker/DofClaim.hh>
#include <protocols/topology_broker/TopologyBroker.fwd.hh>
#include <protocols/topology_broker/TopologyClaimer.fwd.hh>
// AUTO-REMOVED #include <protocols/topology_broker/TopologyClaimer.hh>
// AUTO-REMOVED #include <protocols/topology_broker/weights/AbinitioMoverWeight.hh>
// AUTO-REMOVED #include <protocols/topology_broker/weights/ConstAbinitioMoverWeight.hh>
#include <utility/Bound.fwd.hh>
#include <utility/Bound.hh>
// AUTO-REMOVED #include <utility/PyAssert.hh>
#include <utility/down_cast.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/excn/EXCN_Base.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/file/FileName.fwd.hh>
#include <utility/file/FileName.hh>
#include <utility/file/PathName.fwd.hh>
#include <utility/file/PathName.hh>
#include <utility/file/gzip_util.hh>
#include <utility/io/irstream.fwd.hh>
#include <utility/io/irstream.hh>
#include <utility/io/izstream.fwd.hh>
#include <utility/io/mpistream.hh>
#include <utility/io/orstream.fwd.hh>
#include <utility/io/orstream.hh>
#include <utility/io/ozstream.fwd.hh>
#include <utility/io/zipstream.hpp>
#include <utility/io/zipstream.ipp>
#include <utility/keys/AutoKey.fwd.hh>
#include <utility/keys/AutoKey.hh>
#include <utility/keys/Key.fwd.hh>
#include <utility/keys/Key.hh>
#include <utility/keys/KeyLess.fwd.hh>
#include <utility/keys/KeyLookup.fwd.hh>
#include <utility/keys/KeyLookup.hh>
#include <utility/keys/NoClient.fwd.hh>
#include <utility/keys/NoClient.hh>
#include <utility/keys/SmallKeyVector.fwd.hh>
#include <utility/keys/SmallKeyVector.hh>
#include <utility/keys/UserKey.fwd.hh>
#include <utility/keys/VariantKey.fwd.hh>
#include <utility/keys/VariantKey.hh>
#include <utility/options/AnyOption.fwd.hh>
#include <utility/options/AnyOption.hh>
#include <utility/options/AnyVectorOption.fwd.hh>
#include <utility/options/AnyVectorOption.hh>
#include <utility/options/BooleanOption.fwd.hh>
#include <utility/options/BooleanOption.hh>
#include <utility/options/BooleanVectorOption.fwd.hh>
#include <utility/options/BooleanVectorOption.hh>
#include <utility/options/FileOption.fwd.hh>
#include <utility/options/FileOption.hh>
#include <utility/options/FileVectorOption.fwd.hh>
#include <utility/options/FileVectorOption.hh>
#include <utility/options/IntegerOption.fwd.hh>
#include <utility/options/IntegerOption.hh>
#include <utility/options/IntegerVectorOption.fwd.hh>
#include <utility/options/IntegerVectorOption.hh>
#include <utility/options/Option.fwd.hh>
#include <utility/options/Option.hh>
#include <utility/options/OptionCollection.fwd.hh>
#include <utility/options/OptionCollection.hh>
#include <utility/options/PathOption.fwd.hh>
#include <utility/options/PathOption.hh>
#include <utility/options/PathVectorOption.fwd.hh>
#include <utility/options/PathVectorOption.hh>
#include <utility/options/RealOption.fwd.hh>
#include <utility/options/RealOption.hh>
#include <utility/options/RealVectorOption.fwd.hh>
#include <utility/options/RealVectorOption.hh>
#include <utility/options/ScalarOption.fwd.hh>
#include <utility/options/ScalarOption.hh>
#include <utility/options/ScalarOption_T_.fwd.hh>
#include <utility/options/ScalarOption_T_.hh>
#include <utility/options/StringOption.fwd.hh>
#include <utility/options/StringOption.hh>
#include <utility/options/StringVectorOption.fwd.hh>
#include <utility/options/StringVectorOption.hh>
#include <utility/options/VariantOption.fwd.hh>
#include <utility/options/VariantOption.hh>
#include <utility/options/VectorOption.fwd.hh>
#include <utility/options/VectorOption.hh>
#include <utility/options/VectorOption_T_.fwd.hh>
#include <utility/options/VectorOption_T_.hh>
#include <utility/options/mpi_stderr.hh>
#include <utility/options/keys/AnyOptionKey.fwd.hh>
#include <utility/options/keys/AnyOptionKey.hh>
#include <utility/options/keys/AnyVectorOptionKey.fwd.hh>
#include <utility/options/keys/AnyVectorOptionKey.hh>
#include <utility/options/keys/BooleanOptionKey.fwd.hh>
#include <utility/options/keys/BooleanOptionKey.hh>
#include <utility/options/keys/BooleanVectorOptionKey.fwd.hh>
#include <utility/options/keys/BooleanVectorOptionKey.hh>
#include <utility/options/keys/FileOptionKey.fwd.hh>
#include <utility/options/keys/FileOptionKey.hh>
#include <utility/options/keys/FileVectorOptionKey.fwd.hh>
#include <utility/options/keys/FileVectorOptionKey.hh>
#include <utility/options/keys/IntegerOptionKey.fwd.hh>
#include <utility/options/keys/IntegerOptionKey.hh>
#include <utility/options/keys/IntegerVectorOptionKey.fwd.hh>
#include <utility/options/keys/IntegerVectorOptionKey.hh>
#include <utility/options/keys/OptionKey.fwd.hh>
#include <utility/options/keys/OptionKey.hh>
#include <utility/options/keys/OptionKeys.hh>
#include <utility/options/keys/PathOptionKey.fwd.hh>
#include <utility/options/keys/PathOptionKey.hh>
#include <utility/options/keys/PathVectorOptionKey.fwd.hh>
#include <utility/options/keys/PathVectorOptionKey.hh>
#include <utility/options/keys/RealOptionKey.fwd.hh>
#include <utility/options/keys/RealOptionKey.hh>
#include <utility/options/keys/RealVectorOptionKey.fwd.hh>
#include <utility/options/keys/RealVectorOptionKey.hh>
#include <utility/options/keys/ScalarOptionKey.fwd.hh>
#include <utility/options/keys/ScalarOptionKey.hh>
#include <utility/options/keys/StringOptionKey.fwd.hh>
#include <utility/options/keys/StringOptionKey.hh>
#include <utility/options/keys/StringVectorOptionKey.fwd.hh>
#include <utility/options/keys/StringVectorOptionKey.hh>
#include <utility/options/keys/VectorOptionKey.fwd.hh>
#include <utility/options/keys/VectorOptionKey.hh>
#include <utility/options/keys/all.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
// AUTO-REMOVED #include <utility/signals/BufferedSignalHub.fwd.hh>
// AUTO-REMOVED #include <utility/signals/BufferedSignalHub.hh>
// AUTO-REMOVED #include <utility/signals/Link.fwd.hh>
// AUTO-REMOVED #include <utility/signals/Link.hh>
// AUTO-REMOVED #include <utility/signals/LinkUnit.fwd.hh>
// AUTO-REMOVED #include <utility/signals/LinkUnit.hh>
// AUTO-REMOVED #include <utility/signals/SignalHub.fwd.hh>
// AUTO-REMOVED #include <utility/signals/SignalHub.hh>
// AUTO-REMOVED #include <utility/tag/Tag.fwd.hh>
// AUTO-REMOVED #include <numeric/numeric.functions.hh>
// AUTO-REMOVED #include <numeric/sphericalVector.fwd.hh>
// AUTO-REMOVED #include <numeric/trig.functions.hh>
// AUTO-REMOVED #include <numeric/xyz.functions.fwd.hh>
// AUTO-REMOVED #include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzVector.fwd.hh>
// AUTO-REMOVED #include <numeric/xyzVector.hh>
// AUTO-REMOVED #include <numeric/random/random.fwd.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <ObjexxFCL/TypeTraits.hh>
#include <ObjexxFCL/char.functions.hh>
#include <algorithm>
#include <bitset>
#include <cassert>
// AUTO-REMOVED #include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <deque>
#include <functional>
#include <iomanip>
#include <ios>
#include <iosfwd>
#include <iostream>
#include <istream>
#include <limits>
#include <list>
#include <map>
#include <ostream>
#include <set>
#include <time.h>
#include <utility>
#include <vector>
// AUTO-REMOVED #include <basic/MetricValue.fwd.hh>
#include <basic/Tracer.fwd.hh>
// AUTO-REMOVED #include <basic/datacache/BasicDataCache.fwd.hh>
// AUTO-REMOVED #include <basic/datacache/CacheableData.fwd.hh>
// AUTO-REMOVED #include <basic/datacache/CacheableData.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/keys/AnchoredDesign.OptionKeys.gen.hh>
#include <basic/options/keys/DenovoProteinDesign.OptionKeys.gen.hh>
#include <basic/options/keys/DomainAssembly.OptionKeys.gen.hh>
#include <basic/options/keys/HelixAssembly.OptionKeys.gen.hh>
#include <basic/options/keys/LoopModel.OptionKeys.gen.hh>
#include <basic/options/keys/MM.OptionKeys.gen.hh>
#include <basic/options/keys/MonteCarlo.OptionKeys.gen.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/PCS.OptionKeys.gen.hh>
#include <basic/options/keys/RBSegmentRelax.OptionKeys.gen.hh>
#include <basic/options/keys/SSrbrelax.OptionKeys.gen.hh>
#include <basic/options/keys/abrelax.OptionKeys.gen.hh>
#include <basic/options/keys/antibody.OptionKeys.gen.hh>
#include <basic/options/keys/archive.OptionKeys.gen.hh>
#include <basic/options/keys/assembly.OptionKeys.gen.hh>
#include <basic/options/keys/backrub.OptionKeys.gen.hh>
#include <basic/options/keys/batch_relax.OptionKeys.gen.hh>
#include <basic/options/keys/bbg.OptionKeys.gen.hh>
#include <basic/options/keys/blivens.OptionKeys.gen.hh>
#include <basic/options/keys/boinc.OptionKeys.gen.hh>
#include <basic/options/keys/canonical_sampling.OptionKeys.gen.hh>
#include <basic/options/keys/casp.OptionKeys.gen.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/chrisk.OptionKeys.gen.hh>
#include <basic/options/keys/chunk.OptionKeys.gen.hh>
#include <basic/options/keys/cluster.OptionKeys.gen.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <basic/options/keys/crossmatch.OptionKeys.gen.hh>
#include <basic/options/keys/ddg.OptionKeys.gen.hh>
#include <basic/options/keys/dna.OptionKeys.gen.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <basic/options/keys/dwkulp.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/keys/enzdes.OptionKeys.gen.hh>
#include <basic/options/keys/evaluation.OptionKeys.gen.hh>
#include <basic/options/keys/evolution.OptionKeys.gen.hh>
#include <basic/options/keys/fast_loops.OptionKeys.gen.hh>
#include <basic/options/keys/filters.OptionKeys.gen.hh>
#include <basic/options/keys/fingerprint.OptionKeys.gen.hh>
#include <basic/options/keys/fldsgn.OptionKeys.gen.hh>
#include <basic/options/keys/flexPepDocking.OptionKeys.gen.hh>
#include <basic/options/keys/flexpack.OptionKeys.gen.hh>
#include <basic/options/keys/flxbb.OptionKeys.gen.hh>
#include <basic/options/keys/fold_and_dock.OptionKeys.gen.hh>
#include <basic/options/keys/fold_cst.OptionKeys.gen.hh>
#include <basic/options/keys/fold_from_loops.OptionKeys.gen.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>
#include <basic/options/keys/gpu.OptionKeys.gen.hh>
#include <basic/options/keys/holes.OptionKeys.gen.hh>
#include <basic/options/keys/hotspot.OptionKeys.gen.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/options/keys/james.OptionKeys.gen.hh>
#include <basic/options/keys/jumps.OptionKeys.gen.hh>
#include <basic/options/keys/krassk.OptionKeys.gen.hh>
#include <basic/options/keys/lh.OptionKeys.gen.hh>
#include <basic/options/keys/loopfcst.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/match.OptionKeys.gen.hh>
#include <basic/options/keys/matdes.OptionKeys.gen.hh>
#include <basic/options/keys/mc.OptionKeys.gen.hh>
#include <basic/options/keys/membrane.OptionKeys.gen.hh>
#include <basic/options/keys/motifs.OptionKeys.gen.hh>
#include <basic/options/keys/ms.OptionKeys.gen.hh>
#include <basic/options/keys/murphp.OptionKeys.gen.hh>
#include <basic/options/keys/mysql.OptionKeys.gen.hh>
#include <basic/options/keys/nonlocal.OptionKeys.gen.hh>
#include <basic/options/keys/optE.OptionKeys.gen.hh>
#include <basic/options/keys/optimization.OptionKeys.gen.hh>
#include <basic/options/keys/orbitals.OptionKeys.gen.hh>
#include <basic/options/keys/pH.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/packstat.OptionKeys.gen.hh>
#include <basic/options/keys/parser.OptionKeys.gen.hh>
#include <basic/options/keys/patterson.OptionKeys.gen.hh>
#include <basic/options/keys/pepspec.OptionKeys.gen.hh>
#include <basic/options/keys/phil.OptionKeys.gen.hh>
#include <basic/options/keys/pocket_grid.OptionKeys.gen.hh>
#include <basic/options/keys/pose_metrics.OptionKeys.gen.hh>
#include <basic/options/keys/qsar.OptionKeys.gen.hh>
#include <basic/options/keys/rdc.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/remodel.OptionKeys.gen.hh>
#include <basic/options/keys/resample.OptionKeys.gen.hh>
#include <basic/options/keys/rescore.OptionKeys.gen.hh>
#include <basic/options/keys/residues.OptionKeys.gen.hh>
#include <basic/options/keys/rigid.OptionKeys.gen.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>
#include <basic/options/keys/robert.OptionKeys.gen.hh>
#include <basic/options/keys/rot_anl.OptionKeys.gen.hh>
#include <basic/options/keys/rotamerdump.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/smhybrid.OptionKeys.gen.hh>
#include <basic/options/keys/swa.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/keys/templates.OptionKeys.gen.hh>
#include <basic/options/keys/threadsc.OptionKeys.gen.hh>
#include <basic/options/keys/ufv.OptionKeys.gen.hh>
#include <basic/options/keys/willmatch.OptionKeys.gen.hh>
#include <basic/options/keys/wum.OptionKeys.gen.hh>
#include <basic/options/option.cc.include.gen.hh>
// Auto-header: duplicate removed #include <basic/prof.hh>
#include <boost/bind.hpp>
// AUTO-REMOVED #include <boost/function.hpp>
#include <boost/functional/hash.hpp>
#include <zlib/zlib.h>
#include <zlib/zutil.h>

//Auto Headers


static basic::Tracer tr("protocols.jd2.Archive");
using basic::mem_tr;

OPT_1GRP_KEY( File, iterative, input_pool )
OPT_1GRP_KEY( String, iterative, input_pool_struct_type )

bool protocols::jd2::archive::ArchiveManager::options_registered_( false );

using namespace basic::options;
using namespace basic::options::OptionKeys;
//Mike: when you want to remove these Macros... leave them at least here as comment - since they provide documentation
void protocols::jd2::archive::ArchiveManager::register_options() {
	if ( !options_registered_ ) {
		options_registered_ = true;
		NEW_OPT( iterative::input_pool, "read these structures into pool", "" );
		NEW_OPT( iterative::input_pool_struct_type, "specifies the input-silent-struct type", "protein" );
	}
}

namespace protocols {
namespace jd2 {
namespace archive {

#ifdef WIN32
	void sleep(int seconds){
		//#if (defined WIN32) && (!defined WIN_PYROSETTA)
			Sleep( seconds * 1000 );
		//#endif
	}
#endif

using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core;




std::string
Batch::batch() const {
	return "batch_" + ObjexxFCL::lead_zero_string_of( id(), 6 );
}

std::string Batch::dir() const {
	return batch() + "/";
}

std::string Batch::silent_out() const {
	return batch() + "/decoys.out";
}

std::string Batch::silent_in() const {
	//	if ( has_silent_in() )
	return batch() + "/decoys.in";
	//	else
	//		return "";
}

std::string Batch::score_file() const {
	return batch() + "/score.fsc";
}

std::string Batch::flag_file() const {
	return batch() + "/flags";
}

std::string Batch::broker_file() const {
	return batch() + "/setup.tpb";
}

std::string Batch::extra_broker_files() const {
	return "";
}

void Batch::show( std::ostream& out, bool single_line ) const {
	std::string eol( single_line ? " " : "\n" );
	out << "ID " << id() << eol
			<< "INPUT " << ( has_silent_in() ? "yes" : "no" ) << eol
			<< "NSTRUCT " << nstruct() << eol
			<< "RETURNED " << decoys_returned() << eol
			<< "FINISHED " << ( has_finished() ? "yes" : "no" ) << eol
			<< "CANCELLED " << ( is_cancelled() ? "yes" : "no" ) << eol
			<< "ALLOW_READING_CANCELLED_DECOYS " << ( allow_reading_cancelled_decoys() ? "yes" : "no" ) << eol;
}

std::ostream& operator<< (std::ostream& out, Batch const& batch ) {
	batch.show( out, true );
	return out;
}

void Batch::write_info_file() const {
	utility::io::ozstream out( dir() + "BATCH_INFO" );
	tr.Debug << "write batch info " << dir() << "BATCH_INFO" << std::endl;
	show( out, false /*not single_line*/ );
}

void Batch::read_info_file() {
	core::Size this_id = id(); //to detec errors
	std::string this_batch = batch(); //for error report
	utility::io::izstream in( dir() + "BATCH_INFO" );
	if ( !in.good() ) throw( EXCN_Archive( "cannot find " + dir() + "BATCH_INFO" ) );
	in >> *this;
	if ( this_id != id() ) {
		throw( EXCN_Archive("Inconsistency detected when reading BATCH_INFO for "+ this_batch+" ID in BATCH_INFO is " + batch() ) );
	}
}

//instead of goto statements:   I think goto would be clearer... but there are coding guidlines to adhere...
void report_tag_error( Batch& batch, std::string const& expected_tag, std::string const& tag ) {
	throw( EXCN_Archive( "Error reading batch information for batch: "+batch.batch()+" expected_tag: "+expected_tag+ " found " + tag) );
}

void report_value_error( Batch& batch, std::string const& tag ) {
	throw( EXCN_Archive( "Error reading batch information for batch: "+batch.batch()+" wrong value for tag: "+tag ) );
}

std::istream& operator >> (std::istream& in, Batch &batch ) {
	std::string tag;
	std::string expected_tag;

	in >> tag;
	expected_tag = "ID";
	if ( tag == expected_tag ) {
		in >> batch.batch_id_;
		if ( !in.good() ) report_value_error( batch, tag );
	} else report_tag_error( batch, expected_tag, tag );;

	in >> tag;
	expected_tag = "INPUT";
	if ( tag == expected_tag ) {
		std::string yesno;
		in >> yesno;
		if ( !in.good() ) report_value_error( batch, tag );
		if ( yesno == "yes" ) batch.has_silent_in_ = true;
		else if ( yesno == "no" ) batch.has_silent_in_ = false;
		else report_value_error( batch, tag );
	} else report_tag_error( batch, expected_tag, tag );;

	in >> tag;
	expected_tag = "NSTRUCT";
	if ( tag == expected_tag ) {
		in >> batch.nstruct_;
		if ( !in.good() ) report_value_error( batch, tag );
	} else report_tag_error( batch, expected_tag, tag );;

	in >> tag;
	expected_tag = "RETURNED";
	if ( tag == expected_tag ) {
		in >> batch.decoys_returned_to_archive_;
		if ( !in.good() ) report_value_error( batch, tag );
	} else report_tag_error( batch, expected_tag, tag );;

	in >> tag;
	expected_tag = "FINISHED";
	if ( tag == expected_tag ) {
		std::string yesno;
		in >> yesno;
		if ( !in.good() ) report_value_error( batch, tag );
		if ( yesno == "yes" ) batch.has_finished_ = true;
		else if ( yesno == "no" ) batch.has_finished_ = false;
	} else report_tag_error( batch, expected_tag, tag );;

	in >> tag;
	expected_tag = "CANCELLED";
	if ( tag == expected_tag ) {
		std::string yesno;
		in >> yesno;
		if ( !in.good() ) report_value_error( batch, tag );
		if ( yesno == "yes" ) batch.is_cancelled_ = true;
		else if ( yesno == "no" ) batch.is_cancelled_ = false;
	} else report_tag_error( batch, expected_tag, tag );;

	return in;
}

//#ifndef WIN32

///@details constructor.  Notice it calls the parent class!  It also builds some internal variables for determining
///which processor it is in MPI land.
ArchiveManager::ArchiveManager( core::Size archive_rank, core::Size jd_master_rank, core::Size file_buf_rank ) :
  archive_rank_( archive_rank ),
	jd_master_rank_( jd_master_rank ),
	file_buf_rank_( file_buf_rank ),
	save_archive_time_interval_( 5 )
{
	runtime_assert( options_registered_ );
}

core::Size ArchiveManager::unfinished_batches() const {
	Size unfinished_batches( 0 );
	for ( BatchList::const_iterator it = batches().begin(); it != batches().end(); ++it ) {
		if ( !it->has_finished() && !it->is_cancelled() && it->valid() ) ++unfinished_batches;
	}
	return unfinished_batches;
}

void
ArchiveManager::go()
{
	tr.Debug << "starting ArchiveManager ..." << archive_rank_ << " " << jd_master_rank_ << " " << file_buf_rank_ << std::endl;
	theArchive_ = new abinitio::IterativeAbrelax( this );
	mem_tr << "initialized IterativeAbrelax" << std::endl;
	try {
		if ( !restore_archive() ) {
			if ( option[ OptionKeys::iterative::input_pool ].user() ) {
				std::string const& decoys( option[ OptionKeys::iterative::input_pool ]() );
				tr.Info << "reading decoys from " <<  decoys << " into archive " << std::endl;
				core::io::silent::SilentFileData sfd( decoys, false, false,  option[ OptionKeys::iterative::input_pool_struct_type ]() );
				sfd.read_file( decoys );
				theArchive_->init_from_decoy_set( sfd );
			}
		}
		save_archive();
		read_existing_batches();
	} catch ( utility::excn::EXCN_Base& excn ) {
		send_stop_to_jobdistributor();
		throw;
	}
	//	if ( batches_.size() == 0 ) theArchive_->generate_batch();
	sleep( 5 ); //give JobDistributor time to start up...
#ifdef USEMPI
	MPI_Status status;
	MPI_Request request;
#endif
	bool stop( false );
	bool print_status( true );
	while ( !stop || unfinished_batches() ) {
		if ( print_status ) {
			tr.Debug << "probing for message in ArchiveManager" << std::endl;
			tr.Debug << "\nSTATUS: " << (stop ? "STOP send: " : "" ) << "  ------ unfinished_batches: " << unfinished_batches() << std::endl;
			tr.Debug << "POOL_STATUS: " << std::endl;
			theArchive_->save_status( tr.Debug );
			tr.Debug << "END_STATUS\n\n"<< std::endl;
			print_status = false;
		}
		//is there a message ?
		int flag( -1 );
#ifdef USEMPI
		//no idea why... but 4 request seems to be the magical number... to receive the correct answer ... WEIRD
		MPI_Iprobe( jd_master_rank_, MPI_ARCHIVE_TAG, MPI_COMM_WORLD, &flag, &status );
		MPI_Iprobe( jd_master_rank_, MPI_ARCHIVE_TAG, MPI_COMM_WORLD, &flag, &status );
		MPI_Iprobe( jd_master_rank_, MPI_ARCHIVE_TAG, MPI_COMM_WORLD, &flag, &status );
		MPI_Iprobe( jd_master_rank_, MPI_ARCHIVE_TAG, MPI_COMM_WORLD, &flag, &status );
		MPI_Iprobe( jd_master_rank_, MPI_ARCHIVE_TAG, MPI_COMM_WORLD, &flag, &status );
#endif
		// 	if ( !flag ) { //nothing ...
		// 			//tell JobDistributor, that we are ready to receive message
		// 			int buf[ 4 ];
		// 			buf[ 0 ] = NOTIFICATION_QUERY;
		// 			MPI_Send( &buf, 1, MPI_INT, jd_master_rank_, MPI_JOB_DIST_TAG, MPI_COMM_WORLD );
		// 			sleep( 1 );
		// 			//check if there is something this time...
		// 			MPI_Iprobe( jd_master_rank_, MPI_ARCHIVE_TAG, MPI_COMM_WORLD, &flag, &status );
		// 		}
		try {
			//if there is a message -- go get it.
			int buf[ 6 ]={0,0,0,0,0,0};
			if ( flag ) {
#ifdef USEMPI
				int merrno = MPI_Recv( &buf, 6, MPI_INT, jd_master_rank_, MPI_ARCHIVE_TAG, MPI_COMM_WORLD, &status );
				if ( merrno != MPI_SUCCESS ) tr.Error << "MPI_Recv error " << std::endl;
#endif
			}	else { //nothing received
				idle();
				continue;
			}
			print_status = true;
			// here if we got a message
			Size const msg_tag( buf[ 0 ]);
			tr.Debug << "received message in ArchiveManager " << msg_tag << std::endl;

			switch( msg_tag ) {
			case JOB_COMPLETION: {
				Size const batch_id( buf[ 1 ] );
				bool const final( buf[ 2 ] == 1 );
				Size const bad( buf[ 3 ] );
				Size const good( buf[ 4 ] );
				Size const total( buf[ 5 ] ); //total nr of jobs
				tr.Debug << "ArchiveManager received JOB_COMPLETION " << batch_id << " " << bad << " " << good << " " << total << std::endl;
				jobs_completed_[ batch_id ] = CompletionMessage( batch_id, final, bad, good, total );
				break;
			}
			case QUEUE_EMPTY:	{
				Size const batch_id( buf[ 1 ] );

				//we ignore QUEUE_EMPTY if we know that a new batch has been submitted after issuing of this signal (i.e., the batch-number
				// coming with the message would be smaller than the currently highest batch number... however, there might be invalid batches...
				// find last valid and unfinished batch in list:
				Size max_working_batch_id( batches_.size() );
				if ( batches_.size() ) {
					while ( max_working_batch_id > 0
						&& ( !batches_[ max_working_batch_id ].valid() || batches_[ max_working_batch_id ].has_finished() ) )
						--max_working_batch_id;
					if ( batch_id <= max_working_batch_id ) {
						tr.Info << "ArchiveManager ignored outdated QUEUE_EMPTY with batch_id " << batch_id << " -- already submitted " << batches_.size() << std::endl;
						break;
					}
				}
				//any job-completions we should work thru before generating a new batch?
				PROF_START( basic::ARCHIVE_CRITICAL_JOBSCOMPLETE );
				while ( jobs_completed_.size() ) {
					jobs_completed(); //get thru these before making job decisions
					theArchive_->idle();
				}
				PROF_STOP( basic::ARCHIVE_CRITICAL_JOBSCOMPLETE );

				PROF_START( basic::ARCHIVE_GEN_BATCH );
				//this is a valid QUEUE_EMPTY request: do something about it
				tr.Info << "ArchiveManager received QUEUE_EMPTY" << std::endl;
				tr.Debug << "JD batch_id: " << batch_id << " max_working_batch_id: " << max_working_batch_id << std::endl;

				if ( !theArchive_->finished() ) {
					//if !finished Archive should always generate a batch...
					//but let's make sure by monitoring, since it would be bad if we hang in the communication...
					Size ct( batches_.size() );//monitor number of batches
					if ( !stop ) theArchive_->generate_batch();
					if ( ct == batches_.size() ) { //if generate_batch didn't create anything --- we still owe Jobdistributor a signal
						send_stop_to_jobdistributor(); //send stop
						stop = true;
					}
				} else {
					tr.Debug << "archive is finished ... spinning down" << std::endl;
					send_stop_to_jobdistributor();
					stop = true;
				}
				PROF_STOP( basic::ARCHIVE_GEN_BATCH );
				basic::prof_show();
				break;
			}
			default:
				utility_exit_with_message( "unknown msg in ArchiveManager " + ObjexxFCL::string_of( msg_tag ) );
			}
		} catch ( utility::excn::EXCN_Base &excn ) {
			tr.Error << "[ERROR] " << excn.msg() << std::endl;
			tr.Error << "spinning down" << std::endl;
			save_archive();
			//this usually doesn't work the jobs always run to completion ... let's hard exit for now.
			utility_exit_with_message( "error detected in ArchiveManager -- spinning down" );
			send_stop_to_jobdistributor();
			stop = true;
		}
	} //while loop
	save_archive();
	tr.Info << "ArchiveManager finished !!!" << std::endl;
}

void
ArchiveManager::idle() {

	{ //save archive
		static time_t last_save( time(NULL) );
		time_t now( time( NULL ) );
		Size const elapsedtime( now - last_save );
		if ( elapsedtime > save_archive_time_interval_ ) {
			save_archive();
			last_save = now;
		}
	}

	//	tr.Debug << "idle..." << std::endl;
	if ( jobs_completed_.size() ) {
		PROF_START( basic::ARCHIVE_JOBSCOMPLETE );
		jobs_completed();
		PROF_STOP( basic::ARCHIVE_JOBSCOMPLETE );
		return;
	};

	//	if ( !theArchive_->finished() && theArchive_->ready_for_batch() ) {
		//		theArchive_->generate_batch();
	//	} else {
		time_t before( time(NULL) );
		theArchive_->idle();
		time_t after( time( NULL ) );
		if ( after-before > 1 ) tr.Debug << "spend " << after-before << " seconds in archives idle method... " << std::endl;
		//sleep some more if idle didn't use much time
		if ( after-before < 5 ) sleep( (5 - ( after - before )) );
		//	}
}

void
ArchiveManager::jobs_completed() {// core::Size batch_id, bool final, core::Size bad ) {
	runtime_assert( jobs_completed_.begin() != jobs_completed_.end() );
	CompletionMessage msg = jobs_completed_.begin()->second;
	jobs_completed_.erase( jobs_completed_.begin() );
	Size batch_id( msg.batch_id );
	bool final( msg.final );
	Size bad( msg.bad );
	Size good_decoys( msg.good );
	Batch const& batch( batches_[ batch_id ] );

	// here if in integration-test mode, jump out if not final
	if ( option[ run::constant_seed ] && !final ) return;

	tr.Debug << "jobs_completed for " << batch.batch() << "..." << "already "
					 << batch.decoys_returned() << " decoys known" << std::endl;
	runtime_assert( batch.id() == batch_id );
	WriteOut_MpiFileBuffer file_buf( file_buf_rank_ );
	if ( bad ) {
		///there were some bad-jobs --- we might be at the end of this run... hard to tell
	}
	PROF_START( basic::ARCHIVE_BLOCK_FILE );
	if ( !final ) {
		tr.Debug << "not final ... block file" << std::endl;
		//careful if file isn't in FileBuf anymore... create it just for blocking ? don't block... ?
		file_buf.block_file( ".//"+batch.silent_out() ); //destructor will release file automatically
	} else {
		tr.Debug << "final ... close file " << std::endl;
		file_buf.close_file( ".//"+batch.silent_out() ); //that is not very nice, but file-buf isn't very smart with filenames...
		file_buf.close_file( ".//"+batch.score_file() );
	}
	if ( batch.is_cancelled() && !batch.allow_reading_cancelled_decoys() ) {
		tr.Debug << "returned decoys of cancelled batch.. ignore..." << std::endl;
		return;
	}
	PROF_STOP( basic::ARCHIVE_BLOCK_FILE );
	//sleep( 5 );
	PROF_START( basic::ARCHIVE_READ_DECOYS );
	utility::vector1< std::string > tags_in_file;


	if ( good_decoys ) {
		tr.Debug << "read file " << batch.silent_out() << std::endl;
		utility::io::izstream testin( batch.silent_out() );
		tr.Debug << "stream is " << ( testin.good() ? "good " : "bad" ) << std::endl;
		if ( !testin.good() ) { //this happens sometimes... usually it needs a little bit of waiting and then it works -- NFS lag ?
			//let's look at this later...
			jobs_completed_[ batch_id ] = msg;
			sleep( 5 );
			return;
		}

		using namespace core::io::silent;
		SilentFileData sfd;

		//this keeps order as in file... important since we skip already known tags by just keeping their number
		sfd.read_tags_fast( batch.silent_out(), tags_in_file );

		if ( !final ) {
			tr.Debug << "...and release file" << std::endl;
			file_buf.release_file( ".//"+batch.silent_out() );
		}

		tr.Debug << "found " << tags_in_file.size() << " decoys in " << batch.silent_out() << std::endl;

		utility::vector1< std::string >::iterator iter = tags_in_file.begin();
		for ( Size ct = 1;
					iter != tags_in_file.end() && ct <= batch.decoys_returned(); ++iter, ++ct ) { }; //just skipping...
		utility::vector1< std::string > tags_to_read;

		std::copy( iter, tags_in_file.end(), std::back_inserter( tags_to_read ) );
		if ( tags_to_read.size() ) {
			try {
				sfd.read_file( batch.silent_out(), tags_to_read );
			} catch ( utility::excn::EXCN_Base& excn ) { //or should we be more specific ?
				if ( final ) throw; //rethrow if it is the final version of the file...
				tr.Error << "[ignored ERROR] " << excn.msg() << std::endl;
				tr.Error << "this is not the final version of " << batch.silent_out() << "\n... maybe some data is still held in a cache of the filesystem..."
								 << " let's see if it works better the next time we have to read" << std::endl;
				//or sleep( 5 ) and retry as above ?
				return;
			}

			PROF_STOP( basic::ARCHIVE_READ_DECOYS );
			tr.Debug << "add " << tags_to_read.size() << " structures to archive " << std::endl;

			PROF_START( basic::ARCHIVE_EVAL_DECOYS );
			//read structures and add to archive
			theArchive_->read_structures( sfd, batch );
			PROF_STOP( basic::ARCHIVE_EVAL_DECOYS );
		} else {
			tr.Info << "no more decoys to read from file " << batch.silent_out() << std::endl;
			PROF_STOP( basic::ARCHIVE_READ_DECOYS );
		}


		PROF_START( basic::SAVE_ARCHIVE );
		if ( jobs_completed_.size() == 0 ) save_archive();
		PROF_STOP( basic::SAVE_ARCHIVE );
	} else { // no good decoys found
		tr.Debug << " no good decoys to read " << std::endl;
		throw EXCN_Archive( "all decoys returned with FAIL_BAD_INPUT" );
	}


	{ //now update our batch information and save to disck
		Batch& batch( batches_[ batch_id ] );
		batch.set_decoys_returned( tags_in_file.size() );
		if ( final ) {
			batch.mark_as_finished();
		}
		batch.write_info_file();
	}
}

void
ArchiveManager::queue_batch( Batch const& batch ) {
	tr.Debug << "queue new batch into MPIArchiveJobDistributor " << batch.flag_file() << std::endl;
	Size const size( 3 );
	int buf[ size ];
	buf[ 0 ] = ADD_BATCH;
	buf[ 1 ] = batch.id();
	buf[ 2 ] = batch.nstruct();
#ifdef USEMPI
	MPI_Send( &buf, size, MPI_INT, jd_master_rank_, MPI_JOB_DIST_TAG, MPI_COMM_WORLD );
	//need to have MPI_JOB_DIST_TAG... since it goes into main msg-loop of JobDist

	//send size of string
	std::string strbuf( batch.flag_file() );
	buf[ 0 ] = strbuf.size();
	buf[ 1 ] = batch.id();
	MPI_Send( buf, 2, MPI_INT, jd_master_rank_, MPI_JOB_DIST_TAG, MPI_COMM_WORLD );
	//send string
	MPI_Send(const_cast<char*> ( strbuf.data() ), strbuf.size(), MPI_CHAR, jd_master_rank_, MPI_JOB_DIST_TAG, MPI_COMM_WORLD );
#else
	protocols::jd2::JobDistributor::get_instance()->add_batch( batch.flag_file() );
#endif

}

void ArchiveManager::cancel_batches_previous_to( core::Size batch_id, bool allow_reading_of_decoys ) {
	for ( BatchList::iterator it = batches_.begin(); it!=batches_.end(); ++it) {
		if ( it->id() == batch_id ) break;
		cancel_batch( *it, allow_reading_of_decoys );
	}
}

void
ArchiveManager::cancel_batch( Batch& batch, bool allow_reading_of_decoys ) {
	if ( option[ OptionKeys::run::constant_seed ]() ) {
		tr.Warning << "asked to cancel batch, but ignore in constant_seed mode to enable integration test" << std::endl;
		return;
	}
	tr.Debug << "cancel batch  " << batch.flag_file() << std::endl;
	Size const size( 3 );
	int buf[ size ];
	buf[ 0 ] = CANCEL_BATCH;
	buf[ 1 ] = batch.id();
	buf[ 2 ] = batch.nstruct();
#ifdef USEMPI
	MPI_Send( &buf, size, MPI_INT, jd_master_rank_, MPI_JOB_DIST_TAG, MPI_COMM_WORLD );
	//need to have MPI_JOB_DIST_TAG... since it goes into main msg-loop of JobDist

	//send size of string
	std::string strbuf( BatchJobInputter::BOGUS_BATCH_ID );
	buf[ 0 ] = strbuf.size();
	buf[ 1 ] = batch.id();
	MPI_Send( buf, 2, MPI_INT, jd_master_rank_, MPI_JOB_DIST_TAG, MPI_COMM_WORLD );
	//send string
	MPI_Send(const_cast<char*> ( strbuf.data() ), strbuf.size(), MPI_CHAR, jd_master_rank_, MPI_JOB_DIST_TAG, MPI_COMM_WORLD );
#else
	protocols::jd2::JobDistributor::get_instance()->add_batch( BatchJobInputter::BOGUS_BATCH_ID );
#endif
	batch.mark_as_cancelled( allow_reading_of_decoys );
	batch.write_info_file();
}

void ArchiveManager::send_stop_to_jobdistributor() {

	// we can also do this the quick way:
	//   a) there is still some problem and sometimes jobs-just don't finish... haven't figured out why.
	//   b) we don't really care, if we are here there was either an error,  or the Archive is converged.
	//       at the current low-acceptance rate the remaining jobs are unlikely to yield anything useful...
	//   c) quick exits saves costly time on the cluster.
	if ( option[ OptionKeys::jd2::mpi_nowait_for_remaining_jobs ]() ) {
		save_archive();
		utility_exit_with_message("quick exit from job-distributor due to flag jd2::mpi_nowait_for_remaining_jobs --- this is not an error " );

	//we do this by sending empty batch.
	tr.Debug << "send STOP signal to JobDistributor " << std::endl;
	}
	Batch stop_batch( 0 );
	queue_batch( stop_batch );
}

void
ArchiveManager::read_existing_batches() {
	using utility::file::file_exists;
	using namespace basic::options::OptionKeys;

	//possible options:
	//reads all structures from batches as if they were newly coming in
	bool b_reread_all_structures( option[ OptionKeys::archive::reread_all_structures ]() /* default false */ );

	//don't know how to probe directory... take counter...
	core::Size id( 1 );
	Batch aBatch( id );
	batches_.clear();
	while ( file_exists( aBatch.flag_file() ) ) {
		Batch& new_batch( start_new_batch() );
		runtime_assert( new_batch.id() == id );
		tr.Info << "found existing batch " << new_batch.batch() << std::endl;
		try {
			finalize_batch( new_batch, true /*reread */ );
			tr.Debug << new_batch << std::endl;
		} catch ( EXCN_Archive& excn ) {
			//last started batch must have problems... ignore it
			tr.Warning << "[ WARNING ] "+new_batch.batch()+" is errorneous: " + excn.msg() << std::endl;
			tr.Warning << "[ WARNING ] ignoring this batch..." << std::endl;
			//fill batch list if exception state has left us without it
			if ( batches_.size() < id ) batches_.push_back( Batch( id ) );
			batches_.back().mark_as_invalid();
		}
		if ( b_reread_all_structures ) {
			if ( batches_[ id ].decoys_returned() ) {
				jobs_completed_[ id ] =
					CompletionMessage( id, batches_[id ].has_finished(), 0, batches_[ id ].decoys_returned(), batches_[ id ].nstruct() );
			}
			batches_[ id ].set_decoys_returned( 0 );
		}
		aBatch = Batch( ++id );
	}
}

Batch&
ArchiveManager::start_new_batch() {
	core::io::silent::SilentStructOPs empty;
	return start_new_batch( empty );
}

Batch&
ArchiveManager::start_new_batch( core::io::silent::SilentStructOPs const& start_decoys ) {
	using utility::file::file_exists;

	core::Size batch_id( batches_.size() + 1 );
	tr.Debug << "start new batch " << batch_id << std::endl;
	batches_.push_back( Batch( batch_id ) );
	Batch &new_batch( batches_.back() );

	new_batch.set_id( batch_id );
	//make directory:
	utility::file::create_directory( new_batch.dir() );
	if ( start_decoys.size() ) {
		new_batch.set_has_silent_in();
		core::io::silent::SilentFileData sfd;
		for ( core::io::silent::SilentStructOPs::const_iterator
						it = start_decoys.begin(); it != start_decoys.end(); ++it ) {
			sfd.add_structure( **it );
		}
		sfd.write_all( new_batch.silent_in() );
	}
	new_batch.user_options().add_built_in_options();
	add_all_rosetta_options( new_batch.user_options() );

//copy the system broker setup --- OBSOLET since Sept 20th 2010. now broker:setup is FileVector option.
// 	if ( !file_exists( new_batch.broker_file() ) && option[ OptionKeys::broker::setup ].user() ) {
// 		utility::io::ozstream batch_broker( new_batch.broker_file() );
// 		utility::io::izstream system_broker( option[ OptionKeys::broker::setup ]() );
// 		std::string line;
// 		while ( getline( system_broker, line ) ) batch_broker << line << std::endl;
// 	}
	new_batch.nstruct() = basic::options::option[ basic::options::OptionKeys::out::nstruct ];
	return batches_.back();
}

void report_batch_inconsistency( Batch& new_batch, std::string const &tag ) {
	throw( EXCN_Archive( "inconsistency detected when re-reading "+new_batch.batch()+" for " + tag) );
}

void
ArchiveManager::finalize_batch( Batch& new_batch, bool reread ) {
	using utility::file::file_exists;
	using namespace basic::options::OptionKeys;
	tr.Debug << "finalize_batch " << new_batch << std::endl;
	if ( !reread) new_batch.set_decoys_returned( 0 );

	if ( !utility::file::file_exists( new_batch.broker_file() ) ) {
		utility::io::ozstream broker( new_batch.broker_file() );
		broker << "# NO CLAIMERS PRESENT" << std::endl;
		broker.close();
	}

	if( file_exists( new_batch.flag_file() ) ) {
		tr.Debug << "checking aBatch.flag_file()... " << std::endl;
		utility::options::OptionCollection batch_opts;
		batch_opts.add_built_in_options();
		add_all_rosetta_options( batch_opts );
		try {
			tr.Debug << "load options from file" << std::endl;
			batch_opts.load_options_from_file_exception( new_batch.flag_file() );
		} catch ( utility::excn::EXCN_Msg_Exception &excn ) {
			tr.Error << "[ERROR] problems with flags in " << new_batch.flag_file() << " aborting... " << std::endl;
			// excn.show( tr.Error );
			batches_.pop_back();
			throw ( EXCN_Archive( new_batch.flag_file() + " contains errors: " + excn.msg() ) );
		}
		if ( !reread )  {
			//access all archive controlled options... so they are not in the "user_flags" anymore
			if ( batch_opts[ in::file::silent ].user() )
				tr.Warning << "option -in:file:silent will be overwritten by ArchiveMaster"
									 << " -- control directly via class Batch" << std::endl;
			if ( batch_opts[ out::nstruct ].user() )
				tr.Warning << "option -nstruct will be overwritten by ArchiveMaster "
									 << "-- control directly via class Batch" << std::endl;
			if ( batch_opts[ run::intermediate_structures ].user() )
				tr.Warning << "option -run::intermediate_structures will be overwritten by ArchiveMaster "
									 << "-- control directly via class Batch" << std::endl;
			if ( batch_opts[ out::file::silent ].user() )
				tr.Warning << "option -out:file:silent will be overwritten by ArchiveMaster "
									 << "-- control directly via class Batch" << std::endl;
			if ( batch_opts[ broker::setup ].user() )
				tr.Warning << "option -broker:setup will be overwritten by ArchiveMaster "
									 << "-- control directly via class Batch" << std::endl;
			if ( batch_opts[ out::file::scorefile ].user() )
				tr.Warning << "option -out:file:scorefile will be overwritten by ArchiveMaster "
									 << "-- control directly via class Batch" << std::endl;
		}

		bool has_silent( batch_opts[ in::file::silent ].user() );
		core::Size nstruct( batch_opts[ out::nstruct ]() );
		bool intermeds( batch_opts[ run::intermediate_structures ]() );
		std::string silent_out( batch_opts[ out::file::silent ]() );
		utility::vector1< std::string > broker( batch_opts[ broker::setup ]() );
		std::ostringstream broker_files;
		std::copy( broker.begin(), broker.end(), std::ostream_iterator<std::string>( broker_files, " "));
		std::string score_file( batch_opts[ out::file::scorefile ]() );

		// now the other options are "inaccessed options" and can be dumped to a stream
		std::stringstream user_flags;
		batch_opts.show_inaccessed_user_options( user_flags );
		tr.Debug << "user_options: \n" << user_flags.str() << std::endl;

		// and can be added to the batch-options
		new_batch.user_options().load_options_from_stream( user_flags, "USER_FLAGS" );
		if ( reread ) {
			new_batch.read_info_file();
			new_batch.set_intermediate_structs( intermeds ); //this is not read from BATCH_INFO

			// for all other values we just double-check consistency
			if ( new_batch.nstruct() != nstruct ) report_batch_inconsistency( new_batch, "NSTRUCT" );
			if ( new_batch.has_silent_in() !=  has_silent ) report_batch_inconsistency( new_batch, "INPUT" );
			if ( silent_out != new_batch.silent_out() ) report_batch_inconsistency( new_batch, "OUTPUT" );
			if ( broker_files.str() != new_batch.all_broker_files() ) report_batch_inconsistency( new_batch, "BROKER_FILE" );
			//TODO: determine how many decoys have been returned to archive...
		}
	}

	//now write the final flag-file
	utility::io::ozstream flag_out( new_batch.flag_file() );
	new_batch.user_options().show_user( flag_out );
	flag_out << "\n\n#Archive controlled flags" << std::endl;
	flag_out << "-out:file:silent " << new_batch.silent_out() << std::endl;
	if ( new_batch.has_silent_in() ) flag_out << "-in:file:silent " << new_batch.silent_in() << std::endl;

	flag_out << "-out:nstruct " << new_batch.nstruct() << std::endl;
	flag_out << "-out:file:scorefile " << new_batch.score_file() << std::endl;
	flag_out << "-broker:setup " << new_batch.all_broker_files() << std::endl;

	if ( new_batch.intermediate_structs() ) flag_out << "-run:intermediate_structures" << std::endl;
	if ( !reread ) {
		try {
			test_broker_settings( new_batch );
		} catch ( utility::excn::EXCN_Msg_Exception &excn ) {
			tr.Error << "[ERROR] problems with broker setup in " << new_batch.all_broker_files() << " aborting... " << std::endl;
			// excn.show( tr.Error );
			batches_.pop_back();
			throw ( EXCN_Archive( new_batch.all_broker_files() + " contains errors: " + excn.msg() ) );
		}
	}

	if ( !reread ) {
		new_batch.write_info_file();
	}

	if ( !new_batch.has_finished() && !new_batch.is_cancelled() && theArchive_->still_interested( new_batch ) ) {
		tr.Debug << "queue " << new_batch.batch() << " " << new_batch.flag_file() << std::endl;
		queue_batch( new_batch );
	} else {
		new_batch.mark_as_finished();
	}

	tr.Debug << "\n" << std::endl;

}

void
ArchiveManager::test_broker_settings( Batch const& batch ) {
	tr.Debug << "test broker settings...." << std::endl;
	OptionCollection vanilla_options( option );
  option.load_options_from_file( batch.flag_file() );
	try {
		topology_broker::TopologyBrokerOP topology_broker = new topology_broker::TopologyBroker();
		topology_broker::add_cmdline_claims( *topology_broker );
		tr.Debug << "setting of broker::setup  ";
		utility::vector1< std::string > files( option[ OptionKeys::broker::setup ]() );
		std::copy( files.begin(), files.end(), std::ostream_iterator<std::string>( tr.Debug, " "));
	} catch ( utility::excn::EXCN_Exception &excn ) {  // clean up options and rethrow
		option = vanilla_options;
		throw; //rethrow exception
	}
	option = vanilla_options;
}

void
ArchiveManager::save_archive() {
	theArchive_->save_to_file();
}


bool
ArchiveManager::restore_archive() {
	return theArchive_->restore_from_file();
}

//#endif //ndef WIN32

}//archive
}//jd2
}//protoco



