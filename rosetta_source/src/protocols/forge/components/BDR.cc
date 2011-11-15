// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/forge/components/BDR.cc
/// @brief
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// unit headers
#include <protocols/forge/components/BDR.hh>

// package headers
#include <protocols/forge/build/BuildInstruction.hh>
#include <protocols/forge/build/BuildManager.hh>
#include <protocols/forge/build/SegmentInsert.hh>
#include <protocols/forge/components/VarLengthBuild.hh>
#include <protocols/forge/methods/chainbreak_eval.hh>
#include <protocols/forge/methods/pose_mod.hh>
#include <protocols/forge/methods/util.hh>

// project headers
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <core/pose/Pose.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <basic/MetricValue.hh>
#include <basic/Tracer.hh>
#include <protocols/jumping/Dssp.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/LoopMover_CCD.hh>
#include <protocols/moves/PackRotamersMover.hh>
#include <protocols/toolbox/pose_manipulation.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/NeighborhoodByDistanceCalculator.hh>
#include <protocols/toolbox/task_operations/RestrictToNeighborhoodOperation.hh>
// AUTO-REMOVED #include <protocols/forge/remodel/RemodelConstraintGenerator.hh>

// C++ headers
#include <utility>

#include <platform/types.hh>
#include <core/types.hh>
#include <core/chemical/Adduct.fwd.hh>
#include <core/chemical/Adduct.hh>
#include <core/chemical/AtomICoor.fwd.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/AtomType.fwd.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.fwd.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ElementSet.fwd.hh>
#include <core/chemical/MMAtomType.fwd.hh>
#include <core/chemical/MMAtomTypeSet.fwd.hh>
#include <core/chemical/ResConnID.fwd.hh>
#include <core/chemical/ResConnID.hh>
#include <core/chemical/ResidueConnection.fwd.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/ResidueSelector.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/VariantType.fwd.hh>
#include <core/chemical/types.hh>
#include <core/chemical/orbitals/ICoorOrbitalData.hh>
#include <core/chemical/orbitals/OrbitalType.fwd.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.fwd.hh>
#include <core/chemical/sdf/MolData.fwd.hh>
#include <core/chemical/sdf/MolData.hh>
#include <core/conformation/Atom.fwd.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/PointGraph.fwd.hh>
#include <core/conformation/PointGraphData.fwd.hh>
#include <core/conformation/PseudoBond.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/RotamerSetBase.fwd.hh>
#include <core/conformation/orbitals/OrbitalXYZCoords.hh>
#include <core/conformation/signals/ConnectionEvent.fwd.hh>
#include <core/conformation/signals/ConnectionEvent.hh>
#include <core/conformation/signals/GeneralEvent.fwd.hh>
#include <core/conformation/signals/GeneralEvent.hh>
#include <core/conformation/signals/IdentityEvent.fwd.hh>
#include <core/conformation/signals/IdentityEvent.hh>
#include <core/conformation/signals/LengthEvent.fwd.hh>
#include <core/conformation/signals/LengthEvent.hh>
#include <core/conformation/signals/XYZEvent.fwd.hh>
#include <core/conformation/signals/XYZEvent.hh>
#include <core/fragment/BaseCacheUnit.hh>
#include <core/fragment/ConstantLengthFragSet.fwd.hh>
#include <core/fragment/FragData.fwd.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/FragID.fwd.hh>
#include <core/fragment/FragID_Iterator.fwd.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/Frame.fwd.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameIterator.fwd.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/FrameIteratorWorker_.fwd.hh>
#include <core/fragment/FrameIteratorWorker_.hh>
#include <core/fragment/FrameList.fwd.hh>
#include <core/fragment/FrameList.hh>
#include <core/fragment/OrderedFragSet.fwd.hh>
#include <core/fragment/SingleResidueFragData.fwd.hh>
#include <core/fragment/SingleResidueFragData.hh>
#include <core/graph/ArrayPool.hh>
#include <core/graph/DisjointSets.fwd.hh>
#include <core/graph/Graph.fwd.hh>
#include <core/graph/Graph.hh>
#include <core/graph/UpperEdgeGraph.fwd.hh>
#include <core/graph/unordered_object_pool.fwd.hpp>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/DOF_ID_Map.fwd.hh>
#include <core/id/DOF_ID_Map.hh>
#include <core/id/DOF_ID_Mask.fwd.hh>
#include <core/id/JumpID.fwd.hh>
#include <core/id/JumpID.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <core/id/SequenceMapping.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/TorsionID.hh>
#include <core/id/types.hh>
#include <core/kinematics/AtomPointer.fwd.hh>
#include <core/kinematics/AtomTree.fwd.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/AtomWithDOFChange.fwd.hh>
#include <core/kinematics/DomainMap.fwd.hh>
#include <core/kinematics/Edge.fwd.hh>
#include <core/kinematics/Edge.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MinimizerMapBase.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/RT.fwd.hh>
#include <core/kinematics/RT.hh>
#include <core/kinematics/ResidueCoordinateChangeList.fwd.hh>
#include <core/kinematics/ShortestPathInFoldTree.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/types.hh>
#include <core/kinematics/tree/Atom.fwd.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/pack/dunbrack/RotamerLibrary.fwd.hh>
#include <core/pack/interaction_graph/InteractionGraphBase.fwd.hh>
#include <core/pack/interaction_graph/InteractionGraphBase.hh>
#include <core/pack/interaction_graph/OnTheFlyInteractionGraph.fwd.hh>
#include <core/pack/interaction_graph/PrecomputedPairEnergiesInteractionGraph.fwd.hh>
#include <core/pack/rotamer_set/FixbbRotamerSets.fwd.hh>
#include <core/pack/rotamer_set/FixbbRotamerSets.hh>
#include <core/pack/rotamer_set/RotamerCouplings.fwd.hh>
#include <core/pack/rotamer_set/RotamerSet.fwd.hh>
#include <core/pack/rotamer_set/RotamerSetOperation.fwd.hh>
#include <core/pack/rotamer_set/RotamerSets.fwd.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSetsBase.fwd.hh>
#include <core/pack/rotamer_set/RotamerSetsBase.hh>
#include <core/pack/task/IGEdgeReweightContainer.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/ResfileReader.fwd.hh>
#include <core/pack/task/RotamerSampleOptions.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pack/task/operation/NoRepackDisulfides.fwd.hh>
#include <core/pack/task/operation/ResFilterFactory.fwd.hh>
#include <core/pack/task/operation/ResLvlTaskOperationFactory.fwd.hh>
#include <core/pack/task/operation/TaskOperation.fwd.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.fwd.hh>
#include <core/pose/MiniPose.fwd.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/util.hh>
#include <core/pose/datacache/ObserverCache.fwd.hh>
#include <core/pose/metrics/PoseMetricCalculatorBase.fwd.hh>
#include <core/pose/metrics/PoseMetricCalculatorBase.hh>
#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
#include <core/pose/signals/ConformationEvent.fwd.hh>
#include <core/pose/signals/DestructionEvent.fwd.hh>
#include <core/pose/signals/EnergyEvent.fwd.hh>
#include <core/pose/signals/GeneralEvent.fwd.hh>
#include <core/scoring/ContextGraph.fwd.hh>
#include <core/scoring/ContextGraphTypes.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/EnergiesCacheableDataType.hh>
#include <core/scoring/EnergyGraph.fwd.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/LREnergyContainer.fwd.hh>
#include <core/scoring/MinimizationGraph.fwd.hh>
#include <core/scoring/NeighborList.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionInfo.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/TenANeighborGraph.fwd.hh>
#include <core/scoring/TwelveANeighborGraph.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/Func.fwd.hh>
#include <core/scoring/constraints/Func.hh>
#include <core/scoring/constraints/FuncFactory.fwd.hh>
#include <core/scoring/constraints/HarmonicFunc.fwd.hh>
#include <core/scoring/constraints/HarmonicFunc.hh>
#include <core/scoring/constraints/XYZ_Func.fwd.hh>
#include <core/scoring/methods/ContextDependentLRTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextDependentOneBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextDependentTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentLRTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/EnergyMethod.fwd.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/methods/LongRangeTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/Methods.hh>
#include <core/scoring/methods/TwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/WholeStructureEnergy.fwd.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.fwd.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <protocols/checkpoint/CheckPointer.fwd.hh>
#include <protocols/checkpoint/CheckPointer.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/forge/build/BuildInstruction.fwd.hh>
#include <protocols/forge/build/BuildManager.fwd.hh>
#include <protocols/forge/build/Interval.fwd.hh>
#include <protocols/forge/build/Interval.hh>
#include <protocols/forge/build/SegmentInsert.fwd.hh>
#include <protocols/forge/components/BDR.fwd.hh>
#include <protocols/forge/components/VarLengthBuild.fwd.hh>
#include <protocols/forge/remodel/RemodelConstraintGenerator.fwd.hh>
#include <protocols/forge/remodel/RemodelData.hh>
#include <protocols/jobdist/Jobs.fwd.hh>
#include <protocols/jobdist/Jobs.hh>
#include <protocols/jumping/PairingsList.fwd.hh>
#include <protocols/jumping/StrandPairing.fwd.hh>
#include <protocols/jumping/StrandPairing.hh>
#include <protocols/loops/IndependentLoopMover.hh>
#include <protocols/loops/Loop.fwd.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/LoopMover.fwd.hh>
#include <protocols/loops/LoopMover.hh>
#include <protocols/loops/Loops.fwd.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverStatus.hh>
#include <protocols/moves/PackRotamersMover.fwd.hh>
#include <protocols/toolbox/pose_metric_calculators/NeighborhoodByDistanceCalculator.fwd.hh>
#include <protocols/toolbox/task_operations/RestrictOperationsBase.fwd.hh>
#include <protocols/toolbox/task_operations/RestrictOperationsBase.hh>
#include <protocols/toolbox/task_operations/RestrictToNeighborhoodOperation.fwd.hh>
#include <utility/Bound.fwd.hh>
#include <utility/Bound.hh>
#include <utility/PyAssert.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
#include <utility/stream_util.hh>
#include <utility/string_util.hh>
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
#include <utility/file/FileName.fwd.hh>
#include <utility/file/FileName.hh>
#include <utility/file/PathName.fwd.hh>
#include <utility/file/PathName.hh>
#include <utility/keys/AutoKey.fwd.hh>
#include <utility/keys/AutoKey.hh>
#include <utility/keys/Key.fwd.hh>
#include <utility/keys/Key.hh>
#include <utility/keys/Key2Tuple.fwd.hh>
#include <utility/keys/Key2Tuple.hh>
#include <utility/keys/Key3Tuple.fwd.hh>
#include <utility/keys/Key3Tuple.hh>
#include <utility/keys/Key4Tuple.fwd.hh>
#include <utility/keys/Key4Tuple.hh>
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
#include <utility/signals/BufferedSignalHub.fwd.hh>
#include <utility/signals/BufferedSignalHub.hh>
#include <utility/signals/Link.fwd.hh>
#include <utility/signals/Link.hh>
#include <utility/signals/LinkUnit.fwd.hh>
#include <utility/signals/LinkUnit.hh>
#include <utility/signals/PausableSignalHub.fwd.hh>
#include <utility/signals/PausableSignalHub.hh>
#include <utility/signals/SignalHub.fwd.hh>
#include <utility/signals/SignalHub.hh>
#include <utility/tag/Tag.fwd.hh>
#include <numeric/NumericTraits.hh>
#include <numeric/constants.hh>
#include <numeric/conversions.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/sphericalVector.fwd.hh>
#include <numeric/sphericalVector.hh>
#include <numeric/trig.functions.hh>
#include <numeric/types.hh>
#include <numeric/xyz.functions.fwd.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.fwd.hh>
#include <numeric/xyzVector.hh>
#include <numeric/internal/ColPointers.hh>
#include <numeric/internal/ColVectors.hh>
#include <numeric/internal/ColsPointer.hh>
#include <numeric/internal/RowPointers.hh>
#include <numeric/internal/RowVectors.hh>
#include <numeric/internal/RowsPointer.hh>
#include <numeric/random/random.fwd.hh>
#include <ObjexxFCL/Dimension.fwd.hh>
#include <ObjexxFCL/Dimension.hh>
#include <ObjexxFCL/DimensionExpression.hh>
#include <ObjexxFCL/DynamicIndexRange.fwd.hh>
#include <ObjexxFCL/DynamicIndexRange.hh>
#include <ObjexxFCL/FArray.all.fwd.hh>
#include <ObjexxFCL/FArray.fwd.hh>
#include <ObjexxFCL/FArray.hh>
#include <ObjexxFCL/FArray1.all.fwd.hh>
#include <ObjexxFCL/FArray1.fwd.hh>
#include <ObjexxFCL/FArray1.hh>
#include <ObjexxFCL/FArray1A.fwd.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray1P.fwd.hh>
#include <ObjexxFCL/FArray2.all.fwd.hh>
#include <ObjexxFCL/FArray2.fwd.hh>
#include <ObjexxFCL/FArray2.hh>
#include <ObjexxFCL/FArray2A.fwd.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray2P.fwd.hh>
#include <ObjexxFCL/FArray3.all.fwd.hh>
#include <ObjexxFCL/FArray3.fwd.hh>
#include <ObjexxFCL/FArray3A.fwd.hh>
#include <ObjexxFCL/FArray3D.fwd.hh>
#include <ObjexxFCL/FArray3P.fwd.hh>
#include <ObjexxFCL/FArray4.all.fwd.hh>
#include <ObjexxFCL/FArray4.fwd.hh>
#include <ObjexxFCL/FArray4A.fwd.hh>
#include <ObjexxFCL/FArray4D.fwd.hh>
#include <ObjexxFCL/FArray4P.fwd.hh>
#include <ObjexxFCL/FArray5.all.fwd.hh>
#include <ObjexxFCL/FArray5.fwd.hh>
#include <ObjexxFCL/FArray5A.fwd.hh>
#include <ObjexxFCL/FArray5D.fwd.hh>
#include <ObjexxFCL/FArray5P.fwd.hh>
#include <ObjexxFCL/FArray6.all.fwd.hh>
#include <ObjexxFCL/FArray6.fwd.hh>
#include <ObjexxFCL/FArray6A.fwd.hh>
#include <ObjexxFCL/FArray6D.fwd.hh>
#include <ObjexxFCL/FArray6P.fwd.hh>
#include <ObjexxFCL/FArrayInitializer.fwd.hh>
#include <ObjexxFCL/FArrayInitializer.hh>
#include <ObjexxFCL/FArraySection.fwd.hh>
#include <ObjexxFCL/FArraySection.hh>
#include <ObjexxFCL/FArrayTraits.fwd.hh>
#include <ObjexxFCL/FArrayTraits.hh>
#include <ObjexxFCL/IndexRange.fwd.hh>
#include <ObjexxFCL/IndexRange.hh>
#include <ObjexxFCL/InitializerSentinel.hh>
#include <ObjexxFCL/KeyFArray1D.fwd.hh>
#include <ObjexxFCL/KeyFArray2D.fwd.hh>
#include <ObjexxFCL/KeyFArray3D.fwd.hh>
#include <ObjexxFCL/KeyFArray4D.fwd.hh>
#include <ObjexxFCL/KeyFArray5D.fwd.hh>
#include <ObjexxFCL/KeyFArray6D.fwd.hh>
#include <ObjexxFCL/Observer.fwd.hh>
#include <ObjexxFCL/Observer.hh>
#include <ObjexxFCL/ObserverMulti.hh>
#include <ObjexxFCL/ObserverSingle.hh>
#include <ObjexxFCL/ProxySentinel.hh>
#include <ObjexxFCL/SetWrapper.fwd.hh>
#include <ObjexxFCL/Star.fwd.hh>
#include <ObjexxFCL/Star.hh>
#include <ObjexxFCL/TypeTraits.hh>
#include <ObjexxFCL/char.functions.hh>
#include <ObjexxFCL/proxy_const_assert.hh>
#include <ObjexxFCL/string.functions.hh>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <iomanip>
#include <iosfwd>
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <ostream>
#include <set>
#include <sstream>
#include <string>
#include <typeinfo>
#include <vector>
#include <basic/MetricValue.fwd.hh>
#include <basic/Tracer.fwd.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableData.fwd.hh>
#include <basic/datacache/CacheableData.hh>
#include <basic/datacache/DataCache.fwd.hh>
#include <basic/datacache/DataCache.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/pose_metrics.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <boost/algorithm/string/erase.hpp>
#include <boost/bind.hpp>
#include <boost/config.hpp>
#include <boost/function.hpp>
#include <boost/functional/hash.hpp>
#include <boost/pool/detail/mutex.hpp>
#include <boost/pool/poolfwd.hpp>



namespace protocols {
namespace forge {
namespace components {


static basic::Tracer TR( "protocols.forge.components.BDR" );


/// @brief default constructor
BDR::BDR() :
	Super( "BDR" ),
	use_fullmer_( false ),
	use_sequence_bias_( false ),
	max_linear_chainbreak_( 0.07 ),
	centroid_loop_mover_str_( "RemodelLoopMover" ),
	redesign_loop_neighborhood_( true ),
	dr_cycles_( 3 ),
	centroid_sfx_( core::scoring::ScoreFunctionFactory::create_score_function( "remodel_cen" ) ),
	fullatom_sfx_( core::scoring::ScoreFunctionFactory::create_score_function( core::scoring::STANDARD_WTS, core::scoring::SCORE12_PATCH ) )
{}


/// @brief copy constructor
BDR::BDR( BDR const & rval ) :
	//utility::pointer::ReferenceCount(),
	Super( rval ),
	manager_( rval.manager_ ),
	design_info_( rval.design_info_ ),
	use_fullmer_( rval.use_fullmer_ ),
	use_sequence_bias_( rval.use_sequence_bias_ ),
	max_linear_chainbreak_( rval.max_linear_chainbreak_ ),
	centroid_loop_mover_str_( rval.centroid_loop_mover_str_ ),
	redesign_loop_neighborhood_( rval.redesign_loop_neighborhood_ ),
	resfile_( rval.resfile_ ),
	dr_cycles_( rval.dr_cycles_ ),
	centroid_sfx_( rval.centroid_sfx_->clone() ),
	fullatom_sfx_( rval.fullatom_sfx_->clone() )
{
	if ( rval.vlb_.get() ) {
		vlb_ = new VarLengthBuild( *rval.vlb_ );
	}
}


/// @brief default destructor
BDR::~BDR() {}


/// @brief clone this object
BDR::MoverOP BDR::clone() {
	return new BDR( *this );
}


/// @brief create this type of object
BDR::MoverOP BDR::fresh_instance() {
	return new BDR();
}


/// @brief the centroid level score function, default "remodel_cen"
BDR::ScoreFunction const & BDR::centroid_scorefunction() const {
	return *centroid_sfx_;
}


/// @brief the full-atom level score function, default score12
BDR::ScoreFunction const & BDR::fullatom_scorefunction() const {
	return *fullatom_sfx_;
}


/// @brief add instruction to the manager of this BDR (no copy)
/// @param[in] bi BuildInstruction
/// @param[in] aa_during_design_refine The allowed amino acid sequence
///  during design.  Only applicable to BuildInstructions like
///  SegmentRebuild and SegmentInsert.  Make sure the length of this
///  string matches up properly.  Default empty string.
void BDR::add_instruction(
	BuildInstructionOP bi,
	String const & aa_during_design_refine
)
{
	manager_.add( bi );
	if ( !aa_during_design_refine.empty() ) {
		design_info_.push_back( std::make_pair( bi->original_interval(), aa_during_design_refine ) );
	}

	// additional instruction means we'll need a new re-init the VLB, so
	// go ahead and drop the existing one
	vlb_ = 0;
}


/// @brief create directed dependency between two instructions
void BDR::create_directed_dependency(
	BuildInstructionOP u,
	BuildInstructionOP v
)
{
	manager_.create_directed_dependency( u, v );
}


/// @brief set the centroid level score function
void BDR::centroid_scorefunction( ScoreFunction const & sfx ) {
	centroid_sfx_ = sfx.clone();
}


/// @brief set the centroid level score function
void BDR::centroid_scorefunction( ScoreFunctionOP sfx ) {
	centroid_sfx_ = sfx->clone();
}


/// @brief set the full-atom level score function
void BDR::fullatom_scorefunction( ScoreFunction const & sfx ) {
	fullatom_sfx_ = sfx.clone();
}


/// @brief set the full-atom level score function
void BDR::fullatom_scorefunction( ScoreFunctionOP sfx ) {
	fullatom_sfx_ = sfx->clone();
}


/// @brief apply defined moves to given Pose
void BDR::apply( Pose & pose ) {
	using core::pose::metrics::CalculatorFactory;
	using basic::MetricValue;
	using protocols::jumping::Dssp;
	using protocols::moves::MS_SUCCESS;
	using protocols::moves::FAIL_DO_NOT_RETRY;
	using protocols::moves::FAIL_RETRY;
	using protocols::toolbox::pose_metric_calculators::BuriedUnsatisfiedPolarsCalculator;
	using protocols::toolbox::pose_metric_calculators::NeighborhoodByDistanceCalculator;

	// assign secondary structure
	Dssp dssp( pose );
	dssp.insert_ss_into_pose( pose );

	// do centroid build
	if ( !centroid_build( pose ) ) { // build failed
		set_last_move_status( FAIL_RETRY );
		return;
	}

	// setup calculators
	CalculatorFactory::Instance().remove_calculator( neighborhood_calc_name() );
	CalculatorFactory::Instance().register_calculator(
		neighborhood_calc_name(),
		new NeighborhoodByDistanceCalculator( manager_.union_of_intervals_containing_undefined_positions() )
	);

	// do design-refine iteration
	if ( dr_cycles_ > 0 ) {

		if ( !design_refine( pose ) ) { // design-refine failed
			set_last_move_status( FAIL_RETRY );
			return;
		}

	}

	// if we've gotten to this point, then the structure has been
	// built properly
	set_last_move_status( MS_SUCCESS );

	// setup the PoseMetricCalculators and add them to the evaluators in the
	// JobOutputter
	CalculatorFactory::Instance().remove_calculator( loops_buns_polar_calc_name() );
	CalculatorFactory::Instance().remove_calculator( neighborhood_buns_polar_calc_name() );

	CalculatorFactory::Instance().register_calculator(
		loops_buns_polar_calc_name(),
		new BuriedUnsatisfiedPolarsCalculator(
			"default",
			"default",
			manager_.union_of_intervals_containing_undefined_positions()
		)
	);

	MetricValue< std::set< Size > > loops_neighborhood;
	pose.metric( neighborhood_calc_name(), "neighbors", loops_neighborhood );
	CalculatorFactory::Instance().register_calculator(
		neighborhood_buns_polar_calc_name(),
		new BuriedUnsatisfiedPolarsCalculator(
			"default",
			"default",
			loops_neighborhood.value()
		)
	);
}


std::string
BDR::get_name() const {
	return "BDR";
}

/// @brief run the centroid level build stage
/// @return true if loop closed, false otherwise
bool BDR::centroid_build(
	Pose & pose
) {
	using core::scoring::STANDARD_WTS;
	using core::scoring::SCORE12_PATCH;
	using core::scoring::ScoreFunctionOP;
	using core::scoring::ScoreFunctionFactory;
	using protocols::moves::MS_SUCCESS;

	using core::util::switch_to_residue_type_set;
	using protocols::forge::methods::restore_residues;
	using protocols::toolbox::pose_manipulation::construct_poly_ala_pose;

	// safety, clear the energies object
	pose.energies().clear();

	// make backup Pose for transferring sidechains
	Pose archive_pose = pose;
	Pose modified_archive_pose = archive_pose;
	manager_.modify( modified_archive_pose );

	// ensure modified_archive_pose is completely full-atom, otherwise mismatch
	// will occur when restoring sidechains at the end of the procedure
	bool mod_ap_is_full_atom = true;
	for ( Size i = 1, ie = modified_archive_pose.n_residue(); mod_ap_is_full_atom && i != ie; ++i ) {
		mod_ap_is_full_atom &= ( modified_archive_pose.residue( i ).residue_type_set().name() == core::chemical::FA_STANDARD );
	}

	if ( !mod_ap_is_full_atom ) {
		core::util::switch_to_residue_type_set( modified_archive_pose, core::chemical::FA_STANDARD );
	}

	// flip to poly-ala-gly-pro-disulf pose
	utility::vector1< Size > protein_residues;
	for ( Size i = 1, ie = pose.n_residue(); i <= ie; ++i ) {
		if ( pose.residue( i ).is_protein() ) {
			protein_residues.push_back( i );
		}
	}

	construct_poly_ala_pose( pose, protein_residues, true, true, true );

	// run VLB to build the new section, if no segments have been added/deleted
	// we use the same VLB so that fragment caching works properly
	if ( !vlb_.get() ) {
		vlb_ = new VarLengthBuild( manager_ );
	}

	vlb_->scorefunction( centroid_sfx_ );
	vlb_->vall_memory_usage( VLB_VallMemoryUsage::CLEAR_IF_CACHING_FRAGMENTS );
	vlb_->use_fullmer( use_fullmer_ );
	vlb_->max_linear_chainbreak( max_linear_chainbreak_ );
	vlb_->loop_mover_str( centroid_loop_mover_str_ );

	if ( use_sequence_bias_ ) {
		vlb_->original_sequence( archive_pose.sequence() );
	}

	vlb_->apply( pose );

	if ( vlb_->get_last_move_status() == MS_SUCCESS ) {

		// record the used manager w/ all mapping info
		manager_ = vlb_->manager();

		// safety, clear all the energies before restoring full-atom residues and
		// scoring
		pose.energies().clear();

		// Swap back original sidechains.  At the moment this is a two step process
		// in case any sidechains from SegmentInsert and the like that aren't in the
		// original archive pose need to be transferred.
		restore_residues( modified_archive_pose, pose );
		restore_residues( manager_.original2modified(), archive_pose, pose );

		// go ahead and score w/ full-atom here; we do this in case there are no
		// design-refine cycles -- it's useful to have e.g. rama in the output
		(*fullatom_sfx_)( pose );

		return true; // loop closed

	} else {

		pose = archive_pose;

	}

	return false; // false if loop not closed
}


/// @brief run the design-refine stage
/// @return currently always true
bool BDR::design_refine(
	Pose & pose
)
{
	using core::kinematics::FoldTree;
	using core::pack::task::operation::RestrictResidueToRepacking;
	using core::pack::task::operation::RestrictResidueToRepackingOP;
	using core::pack::task::operation::RestrictToRepacking;
	using core::scoring::STANDARD_WTS;
	using core::scoring::SCORE12_PATCH;
	using core::scoring::ScoreFunctionOP;
	using core::scoring::ScoreFunctionFactory;
	using protocols::forge::build::SegmentInsert;
	using protocols::loops::Loops;
	using protocols::loops::LoopMover_Refine_CCD;
	using protocols::moves::PackRotamersMover;
	using protocols::toolbox::task_operations::RestrictToNeighborhoodOperation;

	using core::pose::annotated_to_oneletter_sequence;
	using protocols::forge::methods::intervals_to_loops;
	using protocols::forge::methods::linear_chainbreak;
	using protocols::loops::remove_cutpoint_variants;

	typedef protocols::forge::build::BuildManager::Positions Positions;

	// collect new regions/positions
	std::set< Interval > loop_intervals = manager_.intervals_containing_undefined_positions();
	Original2Modified original2modified_interval_endpoints = manager_.original2modified_interval_endpoints();

	// collect loops
	Loops loops = intervals_to_loops( loop_intervals.begin(), loop_intervals.end() );

	// refine Mover used doesn't setup a fold tree, so do it here
	FoldTree loop_ft = protocols::forge::methods::fold_tree_from_loops( pose, loops );

	// save original fold tree
	FoldTree original_ft = pose.fold_tree();

	// define the score function
	ScoreFunctionOP sfx = fullatom_sfx_->clone();

	// setup the design TaskFactory
	TaskFactoryOP design_tf = generic_taskfactory();
	design_tf->push_back( new RestrictToNeighborhoodOperation( neighborhood_calc_name() ) );

	if ( !redesign_loop_neighborhood_ ) {
		// set repack only for non-loop positions
		RestrictResidueToRepackingOP repack_op = new RestrictResidueToRepacking();

		Positions new_positions = manager_.new_positions();
		for ( Size i = 1, ie = pose.n_residue(); i != ie; ++i ) {
			if ( new_positions.find( i ) == new_positions.end() ) {
				repack_op->include_residue( i );
			}
		}

		design_tf->push_back( repack_op );
	}

	// add explicit residue types to design task factory if requested
	for ( DesignInfo::const_iterator di = design_info_.begin(), die = design_info_.end(); di != die; ++di ) {
		if ( !di->second.empty() ) {
			String const aa = core::pose::annotated_to_oneletter_sequence( di->second );

			assert( original2modified_interval_endpoints.find( di->first.left ) != original2modified_interval_endpoints.end() );

			if ( aa.find( SegmentInsert::insertion_char() ) != String::npos ) { // SegmentInsert style
				process_insert_design_string( di->first, aa, original2modified_interval_endpoints, design_tf );
			} else { // SegmentRebuild style
				process_continuous_design_string( di->first, aa, original2modified_interval_endpoints, design_tf );
			}
		}
	}

	// setup the refine TaskFactory
	TaskFactoryOP refine_tf = generic_taskfactory();
	refine_tf->push_back( new RestrictToNeighborhoodOperation( neighborhood_calc_name() ) );
	refine_tf->push_back( new RestrictToRepacking() );

	// safety, clear the energies object
	pose.energies().clear();

	// run design-refine cycle
	for ( Size i = 0; i < dr_cycles_; ++i ) {

		// design the new section
		PackRotamersMover design( sfx );
		design.task_factory( design_tf );
		design.apply( pose );

		// set loop topology
		pose.fold_tree( loop_ft );

		// refine the new section
		LoopMover_Refine_CCD refine( loops, sfx );
		refine.false_movemap( manager_.movemap() );
		refine.set_task_factory( refine_tf );
		refine.apply( pose );

		// remove cutpoint variants -- shouldn't this happen at the end
		// of the refine Mover?
		remove_cutpoint_variants( pose );

		// set original topology
		pose.fold_tree( original_ft );
	}

	// must score one last time since we've removed variants and set
	// new topology, otherwise component energies not correct for
	// e.g. structure output
	(*sfx)( pose );

	// evaluate all chainbreaks using linear chainbreak
	bool cbreaks_pass = true;
	for ( Loops::const_iterator l = loops.begin(), le = loops.end(); l != le && cbreaks_pass; ++l ) {
		if ( l->cut() > 0 ) {
			Real const c = linear_chainbreak( pose, l->cut() );
			TR << "design_refine: final chainbreak = " << c << std::endl;
			cbreaks_pass = c <= max_linear_chainbreak_;
		}
	}

	return cbreaks_pass;
}


/// @brief return a TaskFactory useable as a starting point for either
///  design or refinement
BDR::TaskFactoryOP BDR::generic_taskfactory() {
	using core::pack::task::operation::IncludeCurrent;
	using core::pack::task::operation::InitializeFromCommandline;
	using core::pack::task::operation::ReadResfile;
	using core::pack::task::operation::ReadResfileOP;
	using core::pack::task::TaskFactory;
	using core::pack::task::operation::NoRepackDisulfides;

	TaskFactoryOP tf = new TaskFactory();

	tf->push_back( new InitializeFromCommandline() ); // also inits -ex options
	tf->push_back( new IncludeCurrent() ); // enforce keeping of input sidechains
	tf->push_back( new NoRepackDisulfides() );

	// load resfile op only if requested
	if ( !resfile_.empty() ) {
		ReadResfileOP rrf = new ReadResfile();
		rrf->filename( resfile_ );
		tf->push_back( rrf );
	}

	return tf;
}


/// @brief process a continuous design string, adding appropriate operations
///  to the TaskFactory
void BDR::process_continuous_design_string(
	Interval const & original_interval,
	String const & design_str,
	Original2Modified const & original2modified_interval_endpoints,
	TaskFactoryOP design_tf
)
{
	using core::pack::task::operation::RestrictAbsentCanonicalAAS;

	using core::chemical::aa_from_oneletter_code;

	Size const offset = original2modified_interval_endpoints.find( original_interval.left )->second;
	for ( Size i = 0, ie = design_str.length(); i < ie; ++i ) {
		utility::vector1< bool > allowed_aa_types( 20, false );

		switch ( design_str.at( i ) ) {
			case 's': // surface case, no CFWY
				allowed_aa_types = allowed_surface_aa();
				break;
			case '.': // protocol default design
				continue;
			default: // regular case, single aa type
				allowed_aa_types[ aa_from_oneletter_code( design_str.at( i ) ) ] = true;
				break;
		}

		design_tf->push_back( new RestrictAbsentCanonicalAAS( i + offset, allowed_aa_types ) );
	}
}


/// @brief process a design string containing an insert, adding appropriate
///  operations to the TaskFactory
void BDR::process_insert_design_string(
	Interval const & original_interval,
	String const & design_str,
	Original2Modified const & original2modified_interval_endpoints,
	TaskFactoryOP design_tf
)
{
	using core::pack::task::operation::RestrictAbsentCanonicalAAS;
	using core::pack::task::operation::RestrictResidueToRepacking;
	using core::pack::task::operation::RestrictResidueToRepackingOP;
	using protocols::forge::build::Interval;
	using protocols::forge::build::SegmentInsert;

	using core::chemical::aa_from_oneletter_code;

	char const insert_char = SegmentInsert::insertion_char();

	// Figure out the number of residues in each section.
	Interval const interval(
		original2modified_interval_endpoints.find( original_interval.left )->second,
		original2modified_interval_endpoints.find( original_interval.right )->second
	);

	Size const insert_char_idx = design_str.find( insert_char );
	Size const left_nres = insert_char_idx;
	Size const right_nres = design_str.size() - left_nres - 1;
	Size const insert_nres = interval.length() - left_nres - right_nres;

	// Make setup easy by building a new design string to expand the
	// insertion character into a series of the insertion character
	// the size of the insert.
	String aa = design_str;
	aa.replace( insert_char_idx, 1, insert_nres, insert_char );

	// setup TaskOperations
	RestrictResidueToRepackingOP repack_op = new RestrictResidueToRepacking();

	Size const left_offset = interval.left;
	for ( Size i = 0, ie = aa.size(); i < ie; ++i ) {
		utility::vector1< bool > allowed_aa_types( 20, false );

		if ( aa.at( i ) == insert_char ) { // repack only
			repack_op->include_residue( i + left_offset );
			continue;
		} else if ( aa.at( i ) == 's' ) { // surface case, no CFWY
			allowed_aa_types = allowed_surface_aa();
		} else if ( aa.at( i ) == '.' ) { // protocol default design
			continue;
		} else { // regular case, single aa type
			allowed_aa_types[ aa_from_oneletter_code( aa.at( i ) ) ] = true;
		}

		design_tf->push_back( new RestrictAbsentCanonicalAAS( i + left_offset, allowed_aa_types ) );
	}

	design_tf->push_back( repack_op );
}


/// @brief return a boolean vector specifying allowed a.a. when designing
///  on the surface
utility::vector1< bool > const & BDR::allowed_surface_aa() {
	using core::chemical::aa_from_oneletter_code;

	static String surface_aa = "ADEGHIKLMNPQRSTV";
	static utility::vector1< bool > v( 20, false );

	for ( Size i = 0, ie = surface_aa.length(); i < ie; ++i ) {
		v[ aa_from_oneletter_code( surface_aa.at( i ) ) ] = true;
	}

	return v;
}


} // namespace components
} // namespace forge
} // namespace protocols
