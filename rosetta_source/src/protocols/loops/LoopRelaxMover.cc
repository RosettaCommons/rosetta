// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief demo program for implementing loop relax + FA relax
/// @author Mike Tyka
/// @author Srivatsan Raman
/// @author James Thompson
/// @author Daniel J. Mandell

// include these first for building on Visual Studio

#include <protocols/loops/LoopRelaxMover.hh>

#ifdef BOINC_GRAPHICS
#include <protocols/boinc/boinc.hh>
#endif // BOINC_GRAPHICS

// Project Headers

#include <core/chemical/ResidueTypeSet.hh>

#include <core/conformation/ResidueFactory.hh>
#include <core/fragment/FragSet.hh>
#include <core/id/AtomID.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/io/raw_data/DisulfideFile.hh>
#include <core/kinematics/Jump.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/HarmonicFunc.hh>
#include <core/scoring/constraints/BoundConstraint.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/electron_density/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <protocols/checkpoint/CheckPointer.hh>
#include <protocols/comparative_modeling/util.hh>
#include <protocols/evaluation/PoseEvaluator.hh>
#include <protocols/evaluation/RmsdEvaluator.hh>
#include <protocols/evaluation/util.hh>
#include <protocols/idealize/IdealizeMover.hh>
#include <protocols/loops/LoopMover_CCD.hh>
#include <protocols/loops/LoopMoverFactory.hh>
#include <protocols/loops/LoopMover.fwd.hh>
#include <protocols/loops/LoopMover.hh>
#include <protocols/loops/LoopMover_KIC.hh>
#include <protocols/loops/LoopMover_QuickCCD.hh>
#include <protocols/loops/looprelax_protocols.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/util.hh>
#include <protocols/moves/MoverStatus.hh>
#include <protocols/moves/PackRotamersMover.hh>
#include <protocols/viewer/viewers.hh>
#include <utility/exit.hh>

// symmetry
#include <core/pose/symmetry/util.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/util.hh>

#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>
#include <protocols/moves/symmetry/SymPackRotamersMover.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>


// C++ headers
#include <fstream>
#include <iostream>
#include <string>

// option key includes

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
// Auto-header: duplicate removed #include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>

#include <platform/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/Adduct.fwd.hh>
#include <core/chemical/Adduct.hh>
#include <core/chemical/AtomICoor.fwd.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/AtomType.fwd.hh>
#include <core/chemical/AtomTypeSet.fwd.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
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
#include <core/conformation/symmetry/SymDof.fwd.hh>
#include <core/conformation/symmetry/SymmData.fwd.hh>
#include <core/conformation/symmetry/SymmetricConformation.fwd.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>
#include <core/fragment/BaseCacheUnit.hh>
#include <core/fragment/FragData.fwd.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/FragID.fwd.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/fragment/Frame.fwd.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameIterator.fwd.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/FrameIteratorWorker_.fwd.hh>
#include <core/fragment/FrameIteratorWorker_.hh>
#include <core/fragment/FrameList.fwd.hh>
#include <core/fragment/FrameList.hh>
#include <core/fragment/SecondaryStructure.fwd.hh>
#include <core/fragment/SecondaryStructure.hh>
#include <core/fragment/SingleResidueFragData.fwd.hh>
#include <core/fragment/SingleResidueFragData.hh>
#include <core/graph/ArrayPool.hh>
#include <core/graph/Graph.fwd.hh>
#include <core/graph/Graph.hh>
#include <core/graph/UpperEdgeGraph.fwd.hh>
#include <core/graph/unordered_object_pool.fwd.hpp>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID_Mask.fwd.hh>
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
#include <core/io/silent/SilentEnergy.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/io/silent/SilentStruct.fwd.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentStructCreator.fwd.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/kinematics/AtomPointer.fwd.hh>
#include <core/kinematics/AtomTree.fwd.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/AtomWithDOFChange.fwd.hh>
#include <core/kinematics/DomainMap.fwd.hh>
#include <core/kinematics/Edge.fwd.hh>
#include <core/kinematics/Edge.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.fwd.hh>
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
#include <core/pack/interaction_graph/InteractionGraphBase.fwd.hh>
#include <core/pack/rotamer_set/RotamerCouplings.fwd.hh>
#include <core/pack/rotamer_set/RotamerSet.fwd.hh>
#include <core/pack/rotamer_set/RotamerSetOperation.fwd.hh>
#include <core/pack/rotamer_set/RotamerSets.fwd.hh>
#include <core/pack/rotamer_set/symmetry/SymmetricRotamerSets.fwd.hh>
#include <core/pack/task/IGEdgeReweightContainer.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/RotamerSampleOptions.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pack/task/operation/TaskOperation.fwd.hh>
#include <core/pose/MiniPose.fwd.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/datacache/ObserverCache.fwd.hh>
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
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionInfo.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/TenANeighborGraph.fwd.hh>
#include <core/scoring/TwelveANeighborGraph.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/Func.fwd.hh>
#include <core/scoring/constraints/Func.hh>
#include <core/scoring/constraints/FuncFactory.fwd.hh>
#include <core/scoring/constraints/HarmonicFunc.fwd.hh>
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
#include <core/sequence/SequenceAlignment.fwd.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <protocols/evaluation/PoseEvaluator.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/jobdist/Jobs.fwd.hh>
#include <protocols/loops/Loop.fwd.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.fwd.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverCreator.hh>
#include <protocols/moves/PackRotamersMover.fwd.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/relax/MiniRelax.hh>
#include <protocols/relax/util.hh>
#include <protocols/viewer/GraphicsState.hh>
#include <protocols/viewer/triangle.hh>
#include <utility/Bound.fwd.hh>
#include <utility/Bound.hh>
#include <utility/PyAssert.hh>
#include <utility/down_cast.hh>
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
#include <utility/factory/WidgetRegistrator.hh>
#include <utility/file/FileName.fwd.hh>
#include <utility/file/FileName.hh>
#include <utility/file/PathName.fwd.hh>
#include <utility/file/PathName.hh>
#include <utility/io/ozstream.fwd.hh>
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
#include <numeric/HomogeneousTransform.hh>
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
#include <ObjexxFCL/FArray.fwd.hh>
#include <ObjexxFCL/FArray.hh>
#include <ObjexxFCL/FArray1.fwd.hh>
#include <ObjexxFCL/FArray1.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2.fwd.hh>
#include <ObjexxFCL/FArray2.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3.fwd.hh>
#include <ObjexxFCL/FArray3.hh>
#include <ObjexxFCL/FArray3D.fwd.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/FArrayInitializer.fwd.hh>
#include <ObjexxFCL/FArrayInitializer.hh>
#include <ObjexxFCL/FArraySection.fwd.hh>
#include <ObjexxFCL/FArraySection.hh>
#include <ObjexxFCL/FArrayTraits.fwd.hh>
#include <ObjexxFCL/FArrayTraits.hh>
#include <ObjexxFCL/IndexRange.fwd.hh>
#include <ObjexxFCL/IndexRange.hh>
#include <ObjexxFCL/InitializerSentinel.hh>
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
#include <complex>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iosfwd>
#include <limits>
#include <list>
#include <map>
#include <ostream>
#include <set>
#include <sstream>
#include <typeinfo>
#include <utility>
#include <vector>
#include <basic/MetricValue.fwd.hh>
#include <basic/Tracer.fwd.hh>
#include <basic/basic.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableData.fwd.hh>
#include <basic/datacache/CacheableData.hh>
#include <basic/datacache/DataCache.fwd.hh>
#include <basic/datacache/DataCache.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/option.hh>
#include <boost/algorithm/string/erase.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/functional/hash.hpp>
#include <boost/pool/poolfwd.hpp>
#include <boost/unordered/unordered_map.hpp>

//Auto using namespaces
namespace std { } using namespace std; // AUTO USING NS
namespace ObjexxFCL { } using namespace ObjexxFCL; // AUTO USING NS
namespace ObjexxFCL { namespace fmt { } } using namespace ObjexxFCL::fmt; // AUTO USING NS
//Auto using namespaces end
// Auto-header: duplicate removed #include <basic/options/option.hh>


#if defined(WIN32) || defined(__CYGWIN__)
	#include <ctime>
#endif


namespace protocols {
namespace loops {

//constructors
LoopRelaxMover::LoopRelaxMover() {
	set_defaults_();
}

// BE WARNED: THIS CONSTRUCTOR DOES NOT CALL SET_DEFAULTS().
// AS A RESULT, THE SCORE FUNCTIONS (AMONG OTHER THINGS) WILL
// NOT BE INITIALIZED
LoopRelaxMover::LoopRelaxMover(
	std::string const & remodel,
	std::string const & intermedrelax,
	std::string const & refine,
	std::string const & relax,
	loops::Loops loops
) :
	cmd_line_csts_( true ),
	copy_sidechains_( true ),
	n_rebuild_tries_( 1 ),
	rebuild_filter_( 999 ),
	remodel_( remodel ),
	intermedrelax_( intermedrelax ),
	refine_( refine ),
	relax_( relax ),
	loops_( loops )
{}

//destructor
LoopRelaxMover::~LoopRelaxMover() {}

void LoopRelaxMover::frag_libs(
	utility::vector1< core::fragment::FragSetOP > new_libs
) {
		frag_libs_ = new_libs;
}

utility::vector1< core::fragment::FragSetOP >
LoopRelaxMover::frag_libs() const {
	 return frag_libs_;
}

void LoopRelaxMover::scorefxns(
	core::scoring::ScoreFunctionOP centroid_scorefxn,
	core::scoring::ScoreFunctionOP fullatom_scorefxn
) {
	//cen_scorefxn_ = cen_scorefxn;
	//fa_scorefxn_ = fa_scorefxn;
	cen_scorefxn( centroid_scorefxn );
	fa_scorefxn ( fullatom_scorefxn );
}

void LoopRelaxMover::fa_scorefxn(
	core::scoring::ScoreFunctionOP fa_scorefxn
) {
	fa_scorefxn_ = fa_scorefxn;
}

void LoopRelaxMover::cen_scorefxn(
	core::scoring::ScoreFunctionOP cen_scorefxn
) {
	cen_scorefxn_ = cen_scorefxn;
}

void LoopRelaxMover::apply( core::pose::Pose & pose ) {
	int corelength = 0;
	std::string const curr_job_tag( get_current_tag() );

	// store the initial (possibly fullatom) pose
	//   used for computing RMS later
	//   we may steal sidechains from this pose as well
	core::pose::Pose start_pose = pose;

	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	// this typedef belongs in its own .fwd.hh file.
	typedef utility::pointer::owning_ptr< IndependentLoopMover > IndependentLoopMoverOP;

	basic::Tracer TR("protocols.looprelax");

	TR << "==== Loop protocol: ================================================="
		<< std::endl;
	TR << " remodel        " << remodel()        << std::endl;
	TR << " intermedrelax  " << intermedrelax()  << std::endl;
	TR << " refine         " << refine()         << std::endl;
	TR << " relax          " << relax()          << std::endl;

	// load native pose (if provided)
	core::pose::Pose native_pose;
	if ( option[ in::file::native ].user() ) {
		core::import_pose::pose_from_pdb( native_pose, option[ in::file::native ]() );
		core::pose::set_ss_from_phipsi( native_pose );
	} else	{
		native_pose = start_pose;
	}

	if ( start_pose.total_residue() != native_pose.total_residue() ) {
		// strip VRTs from the end, then compare lengths
		int nnonvrt_start = start_pose.total_residue();
		while ( nnonvrt_start>0 && start_pose.residue( nnonvrt_start ).aa() == core::chemical::aa_vrt ) nnonvrt_start--;

		int nnonvrt_native = native_pose.total_residue();
		while (  nnonvrt_native>0 && native_pose.residue( nnonvrt_native ).aa() == core::chemical::aa_vrt ) nnonvrt_native--;
		if ( nnonvrt_native != nnonvrt_start )
			utility_exit_with_message(
				"Start pose and native pose don't match in length"
			);
	}

	evaluation::MetaPoseEvaluatorOP evaluator = new evaluation::MetaPoseEvaluator;
	evaluation::read_common_evaluator_options(*evaluator);
	evaluator->add_evaluation(
		new evaluation::SelectRmsdEvaluator( native_pose, "_native" )
	);

#ifdef BOINC_GRAPHICS
	// set native for graphics
	boinc::Boinc::set_graphics_native_pose( native_pose );
#endif

	// pick loops if necessary
	protocols::loops::Loops loops = get_loops();
	// try to load loops from command line
	if ( loops.size() == 0 )
		protocols::loops::Loops loops = protocols::loops::get_loops_from_file();
	if ( loops.size() == 0 ) {
		TR.Debug << "picking loops by chainbreak score." << std::endl;
		loops = protocols::comparative_modeling::pick_loops_chainbreak(
			start_pose, option[ cm::min_loop_size ]()
		);

		if ( loops.size() == 0 ) {
			TR.Debug << "no loops found." << std::endl;
			remodel( "no" );
		}
	} // loops.size() == 0

	loops.verify_against( start_pose );
	TR.Debug << loops << std::endl;

	if ( option[ OptionKeys::loops::extended ]() ) loops.set_extended( true );
	bool const debug( option[ OptionKeys::loops::debug ]() );

	// superimpose native over core ?
	core::pose::Pose native_pose_super = native_pose;
	id::AtomID_Map< id::AtomID > atom_map;
	if(  option[ OptionKeys::loops::superimpose_native ]()  ){
		core::pose::initialize_atomid_map( atom_map, native_pose_super, core::id::BOGUS_ATOM_ID );
		for ( core::Size ir=1; ir <= native_pose.total_residue(); ++ir ) {
			if( !loops.is_loop_residue( ir ) ){
				id::AtomID const id1( native_pose_super.residue(ir).atom_index("CA"), ir );
				id::AtomID const id2( pose.residue(ir).atom_index("CA"), ir );
				atom_map.set(id1, id2);
			}
		}
		/*core::Real rms =*/ core::scoring::superimpose_pose( native_pose_super, pose, atom_map );
		if ( debug ) pose.dump_pdb(curr_job_tag + "_after_super_pose.pdb");
		if ( debug ) native_pose_super.dump_pdb(curr_job_tag + "_after_super_native.pdb");
	}

	bool const fullatom_input( pose.is_fullatom() );
	bool const fullatom_output(
		option[ out::file::fullatom ]() || refine() != "no" || relax() != "no"
	);

	if (fullatom_input) {
		core::util::switch_to_residue_type_set( pose, core::chemical::CENTROID );
	}


	checkpoint::CheckPointer checkpoints_("Loopbuild");
	// need to make sure that the bailout structure is set!
	checkpoints_.checkpoint( pose, curr_job_tag, "initial", true );

#ifdef BOINC_GRAPHICS
	// attach boinc graphics pose observer
	boinc::Boinc::attach_graphics_current_pose_observer( pose );
#endif

#ifdef GL_GRAPHICS
	protocols::viewer::add_conformation_viewer(
		pose.conformation(), "loops_pose"
	);
#endif
	////////////////////////////

	// loop rebuilding
	loops.auto_choose_cutpoints( pose );

	// if remove_extended_loops is specified, treat extended loops as missing density, by randomly placing the atoms
	// this is not great behavior, BUT the code is already tolerant of missing density treated in this fashion
	bool remove_extended_loops = option[ OptionKeys::loops::remove_extended_loops ]();
	if (remove_extended_loops) {
		for ( loops::Loops::const_iterator it = loops.begin(), it_end = loops.end();
					it != it_end; ++it
		) {
			if ( it->is_extended() && it->skip_rate() == 0.0 ) {
				TR << "Removing loop: " << *it << std::endl;
				int lstart = it->start(); if (lstart != 1) lstart++;

				for (core::Size r = lstart; r<= it->stop(); ++r) {
					for (core::Size k = 1; k<= pose.residue(r).natoms(); ++k) {
						numeric::xyzVector< core::Real > rnd_atm(
							900.000 + numeric::random::uniform()*100.000,
							900.000 + numeric::random::uniform()*100.000,
							900.000 + numeric::random::uniform()*100.000
						);
						pose.set_xyz( core::id::AtomID( k,r ) , rnd_atm );
					}
				}
			}
		}
	} // remove_extended_loops

	///////////////////////////////////////////////////////////////////////////////////////
	////
	////  	Remove Missing density
	////
	//// if this is an initial loop building exercise, build conservatively first, i.e. dont extend loop
	//// regions etc. Then repeat with a more aggressive approach to make sure the loops are closed etc.
	//// This makes sure all the missing density is gone before doing proper loop building.
	////
	if ( option[ OptionKeys::loops::build_initial ].user() ) {
		TR << "====================================================================================" << std::endl;
		TR << "===" << std::endl;
		TR << "===   Initial Building     " << std::endl;
		TR << "===" << std::endl;

		if ( !checkpoints_.recover_checkpoint( pose, curr_job_tag, "initial_build", false, true) ) {
			if ( debug ) pose.dump_pdb(curr_job_tag + "_before_initial_build.pdb");

			runtime_assert( frag_libs().size() > 0 );
			loops::LoopMover_Perturb_QuickCCD quick_ccd( loops );
			for ( Size i = 1; i <= frag_libs().size(); ++i ) {
				quick_ccd.add_fragments( frag_libs()[i] );
			}

			quick_ccd.get_checkpoints().set_type("InitialBuild");
			quick_ccd.set_current_tag( curr_job_tag );
			quick_ccd.set_native_pose( new core::pose::Pose ( native_pose ) );
			quick_ccd.set_scorefxn( cen_scorefxn_ );
			quick_ccd.set_build_attempts_( 1 );
			quick_ccd.set_grow_attempts_( 0 );
			quick_ccd.set_accept_aborted_loops_( true );
			quick_ccd.set_strict_loops( true );
			quick_ccd.set_random_order_( false );
			quick_ccd.set_build_all_loops_( true );
			quick_ccd.set_loop_combine_rate_( 0.0 );
			quick_ccd.apply( pose );
			loops::remove_cutpoint_variants( pose );

			//fpd if we care at all about the global coordinate frame
			//fpd (and thus have a root VRT) then don't recenter the pose
			if ( pose.residue_type( pose.fold_tree().root() ).aa() != core::chemical::aa_vrt )
				pose.center();
			(*cen_scorefxn_)(pose);
			if ( debug ) pose.dump_pdb(curr_job_tag + "_after_initial_build.pdb");

			checkpoints_.checkpoint( pose, curr_job_tag, "initial_build", true);
		}
		checkpoints_.debug( curr_job_tag, "initial_build", (*cen_scorefxn_)( pose ) );
	} // build_initial

	// Make sure loops can be grown in any protocol, not just QuickCCD
	if ( basic::options::option[ basic::options::OptionKeys::loops::random_grow_loops_by ].user() ) {
		loops.grow_all_loops( pose ,  basic::options::option[ basic::options::OptionKeys::loops::random_grow_loops_by ]() );
	}

	///////////////////////////////////////////////////////////////////////////////////////
	////
	//// 	add constraints if specified by user.
	////
	////
	if ( cmd_line_csts() ) {
		core::scoring::constraints::add_constraints_from_cmdline(
			pose, *cen_scorefxn_
		);
	} else {
		core::scoring::constraints::add_constraints_from_cmdline_to_scorefxn(
			*cen_scorefxn_
		);
	}

	// read in disulfides if specified by user
	if ( option[ in::fix_disulf ].user() ) {
		using std::pair;
		using utility::vector1;
		vector1< pair<Size, Size> > disulfides;

		core::io::raw_data::DisulfideFile ds_file( option[ in::fix_disulf ]() );
		ds_file.disulfides( disulfides, pose );
		pose.conformation().fix_disulfides( disulfides );
	}



	///////////////////////////////////////////////////////////////////////////////////////
	////
	//// 	Loop remodelling (centroid loop modelling)
	////
	////


	long starttime = time(NULL);
	if ( remodel() != "no" ) {
		TR << "====================================================================================" << std::endl;
		TR << "===" << std::endl;
		TR << "===   Remodel    " << std::endl;
		TR << "===" << std::endl;

		if ( !checkpoints_.recover_checkpoint( pose, curr_job_tag, "remodel", false, true) ) {

			if ( debug ) pose.dump_pdb(curr_job_tag + "_before_rebuild.pdb");

			using core::Size;
			using core::Real;
			using core::pose::Pose;
			for ( Size ii = 1; ii <= n_rebuild_tries(); ++ii ) {
				core::Real current_sc( rebuild_filter() + 1 );
				TR.Debug << "Remodeling attempt " << ii << "." << std::endl;
				if ( remodel() == "old_loop_relax") {
					LoopRebuild loop_rebuild( cen_scorefxn_, loops );
					loop_rebuild.apply( pose );
				} else {

/* // DJM: does this cause a crash if the only loop is terminal
          if ( remodel() == "perturb_kic" ) {
            // remove the terminal loops - perturb_kic doesnt seem to support these, sadly.
            protocols::loops::Loops newloops;
            for( core::Size i=1; i <= loops.size() ; i ++ ){
              if( loops[i].start() <= 1 ) continue;
              if( loops[i].stop() >= pose.total_residue() ) continue;
              newloops.add_loop( loops[i] );
            }
            loops = newloops;
          }
 */
          // DJM: need to cast this as IndependentLoopMover to set strict loops to true.
					IndependentLoopMoverOP remodel_mover( static_cast< loops::IndependentLoopMover * >
														 ( loops::get_loop_mover( remodel(), loops ).get() ) );
					core::kinematics::FoldTree f_orig=pose.fold_tree();
					if ( !remodel_mover ) {
						utility_exit_with_message( "Error: no remodel mover defined!" );
					}
					if ( ! ( remodel() == "perturb_kic" ) ) {
						runtime_assert( frag_libs().size() > 0 );
					}

					for ( Size i = 1; i <= frag_libs().size(); ++i ) {
						remodel_mover->add_fragments( frag_libs()[i] );
					}

					remodel_mover->set_scorefxn( cen_scorefxn_ );  // set cst score

					if ( remodel() == "perturb_kic" ) {
						core::kinematics::FoldTree f_new;
						protocols::loops::fold_tree_from_loops( pose, loops,  f_new, true );
						pose.fold_tree( f_new );
					}
					remodel_mover->get_checkpoints().set_type("Remodel");
					remodel_mover->set_current_tag( curr_job_tag );
					remodel_mover->set_native_pose( new Pose( native_pose ) );
					remodel_mover->apply( pose );

					if ( remodel() == "perturb_kic" ) { // DJM: skip this struct if initial closure fails
						if ( remodel_mover->get_last_move_status() != protocols::moves::MS_SUCCESS) {
							set_last_move_status(protocols::moves::FAIL_RETRY);
							TR << "Structure " << " failed initial kinematic closure. Skipping..." << std::endl;
							//	bool fail = true; // make this
							pose.fold_tree( f_orig );
							return;
						}
						pose.fold_tree( f_orig );
					} // if ( remodel() == "perturb_kic" )
				}
				current_sc = (*cen_scorefxn_)( pose );
				if ( current_sc <= rebuild_filter() ) break;  //fpd changed >= to <=
				                                              // ... don't we want to stop when our score is _less_ than the cutoff
			} // for ( in n_rebuild_tries )

			if ( debug ) pose.dump_pdb(curr_job_tag + "_after_rebuild.pdb");
			checkpoints_.checkpoint( pose, curr_job_tag, "remodel", true);
		} // recover checkpoint
		checkpoints_.debug( curr_job_tag, "remodel", (*cen_scorefxn_)( pose ) );
	} // if ( remodel != no )
	long endtime = time(NULL);

	TR << "Buildtime: " << endtime - starttime << std::endl;


	TR << pose.fold_tree() << std::endl;

	///////////////////////////////////////////////////////////////////////////////////////
	////
	////  Halfway stats
	////
	////
	setPoseExtraScores( pose, "cen_irms",  core::scoring::CA_rmsd( start_pose, pose ) );
	if ( option[ in::file::native ].user() ) {
		setPoseExtraScores( pose, "cen_rms",   core::scoring::native_CA_rmsd(native_pose, pose ) );
		setPoseExtraScores( pose, "cen_looprms",  loops::loop_rmsd(native_pose_super, pose, loops ) );
		setPoseExtraScores( pose, "cen_loopcarms",  loops::loop_rmsd(native_pose_super, pose, loops, true ) );
	}


	///////////////////////////////////////////////////////////////////////////////////////
	////
	////  Full atom part
	////
	////
	if ( fullatom_output ) {
		// make sure there aren't any cut point variants flying about

		loops::remove_cutpoint_variants( pose, true );

		// if no centroid modelling was done at all, then grab the original
		// fullatom pose.
		if ( remodel() == "no" &&
				!option[ OptionKeys::loops::build_initial ].user() &&
		     fullatom_input
		) {
			pose = start_pose;
		}

		TR << "===================================================================================="
			<< std::endl;
		TR << "===" << std::endl;
		TR << "===   Fullatom " << std::endl;
		TR << "===" << std::endl;

		if ( debug ) pose.dump_pdb(curr_job_tag + "_before_fullatom.pdb");
		TR << "Annotated sequence before fa switch: " << pose.annotated_sequence(true) << std::endl;

		//puts full-atom sidechains on loop regions
		core::util::switch_to_residue_type_set( pose, core::chemical::FA_STANDARD );
		pose.conformation().detect_bonds(); //apl fix this !

		utility::vector1< bool > needToRepack( pose.total_residue() , !fullatom_input );
		bool needToRepackAtAll = !fullatom_input;

		if ( debug ) pose.dump_pdb(curr_job_tag + "_after_fullatom.pdb");
		// copy sidechain torsions from input pose
		if ( fullatom_input && copy_sidechains() )  {

			for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
				// if remodelling was done, repack the loops - otherwise leave it.
				if ( remodel() != "no" ) {
					for ( int j = 1; j <= (int)loops.size(); ++j ) {
						if (    i >= core::Size( loops[j].start() ) - 3
						     && i <= core::Size( loops[j].stop() ) + 3 ) {
							// allow 3-residue leeway on either side for 'random_loops'
							// this kind of sucks.
							TR.Debug << "Repacking because in loop: " << i << std::endl;
							needToRepack[i] = true;
							break;
						}
					}
				}

				// if there is missing density in the sidechain, then we need to
				// repack check this my making sure that no SC atom is more than
				// 20A (?) away from CA
				if ( start_pose.residue_type(i).is_protein() && start_pose.residue_type(i).has("CA") ) {
					numeric::xyzVector< core::Real> ca_pos = start_pose.residue(i).atom("CA").xyz();
					for (int j=1; j<=(int)start_pose.residue(i).natoms(); ++j) {
						if ( (ca_pos - start_pose.residue(i).atom(j).xyz()).length() > 20 ) {
							TR.Debug << "Missing dens: " << i << std::endl;
							needToRepack[i] = true;
							break;
						}
					}
				} // missing density?

				//copy sidechains only for non-loop regions
				if ( !needToRepack[i] ) {
					if ( pose.residue_type(i).is_protein() )
						pose.replace_residue( i, start_pose.residue(i), true );
					TR.Debug << "Copying sidechain from template: " << i << std::endl;
				} else {
					needToRepackAtAll = true;
				}
			} // for ( i in pose.total_residue() )
		} // fa_input

		// create score function and add constraints for fullatom part
		if ( cmd_line_csts() ) {
			core::scoring::constraints::add_fa_constraints_from_cmdline(
				pose, *fa_scorefxn_
			);
		} else {
			core::scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn(
				*fa_scorefxn_
			);
		}



		// Add coordinate constraints to non-loop regions if desired
		core::Real constrain_rigid_segments_weight = option[ OptionKeys::loops::constrain_rigid_segments ]();
		if( constrain_rigid_segments_weight > 0.0 ){
			protocols::loops::Loops coordconstraint_segments;
			//core::pose::Pose coordconstrainted_pose = pose;
			core::pose::Pose constraint_target_pose = pose;


			coordconstraint_segments = loops.invert( pose.total_residue() );  // Invert the loops selection - i.e. the rigid areas are now defined
			std::cout << "Restraining the following segments: " << std::endl << coordconstraint_segments << std::endl;

			// ResidueTypeSet
			using namespace core;
			using namespace conformation;
			using namespace core::scoring;
			using namespace core::scoring::constraints;
			using namespace id;

			if ( pose.residue( pose.total_residue() ).aa() != core::chemical::aa_vrt ) {
				pose.append_residue_by_jump(
					*ResidueFactory::create_residue( pose.residue(1).residue_type_set().name_map( "VRT" ) ), pose.total_residue()/2 );
			}

			//fpd  nmonomerres is #residues in a single subunit (excluding virtuals)
			core::Size rootres = pose.fold_tree().root();
			core::Size nmonomerres = pose.total_residue()-1;
			core::conformation::symmetry::SymmetryInfoCOP symm_info;
			if ( core::pose::symmetry::is_symmetric(pose) ) {
				core::conformation::symmetry::SymmetricConformation & SymmConf (
					dynamic_cast<core::conformation::symmetry::SymmetricConformation &> ( pose.conformation()) );
				symm_info = SymmConf.Symmetry_Info();
				nmonomerres = symm_info->num_independent_residues();
			}

			if (!option[ OptionKeys::relax::coord_cst_width ].user() ) {
				Real const coord_sdev( option[ OptionKeys::relax::coord_cst_stdev ] );
					// default is 0.5 (from idealize) -- maybe too small
				for ( Size i = 1; i<=nmonomerres; ++i ) {
					if ( !pose.residue(i).is_polymer() ) continue;
					if ( coordconstraint_segments.is_loop_residue( i ) ) {
						Residue const & nat_i_rsd( constraint_target_pose.residue(i) );
						for ( Size ii = 1; ii<=nat_i_rsd.last_backbone_atom(); ++ii ) {
							pose.add_constraint( new CoordinateConstraint( AtomID(ii,i), AtomID(1,rootres), nat_i_rsd.xyz( ii ),
							                     new HarmonicFunc( 0.0, coord_sdev ) ) );
						}

						// now cst symmetry mates
						// if (symm_info) {
						// 	for ( core::conformation::symmetry::SymmetryInfo::Clones::const_iterator pos=symm_info->bb_clones( i ).begin(),
						// 		  epos=symm_info->bb_clones( i ).end(); pos != epos; ++pos ) {
						// 		for ( Size ii = 1; ii<= nat_i_rsd.last_backbone_atom(); ++ii ) {
						// 			pose.add_constraint( new CoordinateConstraint( AtomID(ii,*pos), AtomID(1,nres), nat_i_rsd.xyz( ii ),
						// 								 new HarmonicFunc( 0.0, coord_sdev ) ) );
						// 		}
						// 	}
						// }
					}
				}
			} else {
				Real const cst_width( option[ OptionKeys::relax::coord_cst_width ]() );
				Real const coord_sdev( option[ OptionKeys::relax::coord_cst_stdev ]() );
				for ( Size i = 1; i<=nmonomerres; ++i ) {
					if ( !pose.residue(i).is_polymer() ) continue;
					if( coordconstraint_segments.is_loop_residue( i ) ) {
						Residue const & nat_i_rsd( constraint_target_pose.residue(i) );
						for ( Size ii = 1; ii<= nat_i_rsd.last_backbone_atom(); ++ii ) {
							pose.add_constraint( new CoordinateConstraint( AtomID(ii,i), AtomID(1,rootres), nat_i_rsd.xyz( ii ),
							                     new BoundFunc( 0, cst_width, coord_sdev, "xyz" )) );
						}
						// now cst symmetry mates
						// if (symm_info) {
						// 	for ( core::conformation::symmetry::SymmetryInfo::Clones::const_iterator pos=symm_info->bb_clones( i ).begin(),
						// 		  epos=symm_info->bb_clones( i ).end(); pos != epos; ++pos ) {
						// 		for ( Size ii = 1; ii<= nat_i_rsd.last_backbone_atom(); ++ii ) {
						// 			pose.add_constraint( new CoordinateConstraint( AtomID(ii,*pos), AtomID(1,nres), nat_i_rsd.xyz( ii ),
						// 								 new BoundFunc( 0, cst_width, coord_sdev, "xyz" )) );
						// 		}
						// 	}
						// }
					}
				}
			}



			fa_scorefxn_->set_weight( coordinate_constraint, constrain_rigid_segments_weight );
		}




		// ----------------------------------------------------------





		// do the same (again) for fit-to-density
		if ( option[ edensity::mapfile ].user() ) {
			core::scoring::electron_density::add_dens_scores_from_cmdline_to_scorefxn( *fa_scorefxn_ );
		}

		if ( debug ) pose.dump_pdb(curr_job_tag + "_before_repack.pdb");

		if ( needToRepackAtAll ) { // kic refine does its own initial repacking
			TR << "Repacking required" << std::endl;
			TR << "Detecting disulfides" << std::endl;
			TR << "Annotated sequence before repack: " << pose.annotated_sequence(true) << std::endl;

			// repack loop + missing-density residues
			core::pack::task::PackerTaskOP taskstd
				= core::pack::task::TaskFactory::create_packer_task( pose );
			taskstd->restrict_to_repacking();
			taskstd->or_include_current(true);
			core::pose::symmetry::make_residue_mask_symmetric( pose, needToRepack );
			             // does nothing if pose is not symm
			taskstd->restrict_to_residues(needToRepack);

			fa_scorefxn_->show_line( TR, pose );
			core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap );  //fpd symmetrize this
			mm->set_bb( false );
			mm->set_chi( true );

			//fpd symmetrize this
			if ( core::pose::symmetry::is_symmetric( pose ) ) {
				moves::symmetry::SymPackRotamersMover pack1( fa_scorefxn_, taskstd );
				pack1.apply( pose );

				core::optimization::symmetry::SymAtomTreeMinimizer mzr;
				core::optimization::MinimizerOptions options("dfpmin_armijo_nonmonotone", 1e-5, true, false);
				core::pose::symmetry::make_symmetric_movemap( pose, *mm );

				mzr.run( pose, *mm, *fa_scorefxn_, options );
			} else {
				moves::PackRotamersMover pack1( fa_scorefxn_, taskstd );
				pack1.apply( pose );

				// quick SC minimization
				core::optimization::AtomTreeMinimizer mzr;
				core::optimization::MinimizerOptions options("dfpmin_armijo_nonmonotone", 1e-5, true, false);
				mzr.run( pose, *mm, *fa_scorefxn_, options );
			}

			fa_scorefxn_->show_line( TR, pose );
		} else {
			//fpd
			(*fa_scorefxn_)(pose);
			TR << "No repacking required" << std::endl;
			(*fa_scorefxn_)(pose);
		}

		if ( debug ) pose.dump_pdb(curr_job_tag + "_after_repack.pdb");
	} // fullatom_output

	//////////////////////////////////////////////////////////////////////////////////////
	////
	////  intermediate relax the structure
	////
	////
	if ( intermedrelax() != "no" ) {
		TR << "====================================================================================" << std::endl;
		TR << "===" << std::endl;
		TR << "===   Intermediate Relax  " << std::endl;
		TR << "===" << std::endl;

		core::kinematics::FoldTree f_new, f_orig=pose.fold_tree();
		if ( option[ OptionKeys::loops::relax_with_foldtree ].user() ){
			loops::fold_tree_from_loops( pose, loops, f_new );
			pose.fold_tree( f_new );
			loops::add_cutpoint_variants( pose );
			fa_scorefxn_->set_weight( core::scoring::chainbreak,        option[ OptionKeys::relax::chainbreak_weight ]() );
			fa_scorefxn_->set_weight( core::scoring::linear_chainbreak, option[ OptionKeys::relax::linear_chainbreak_weight ]() );
			fa_scorefxn_->set_weight( core::scoring::overlap_chainbreak, option[ OptionKeys::relax::overlap_chainbreak_weight ]() );
		}


		TR << pose.fold_tree() << std::endl;
		setPoseExtraScores( pose, "brlx_irms",  core::scoring::CA_rmsd( start_pose, pose ) );
		if ( option[ in::file::native ].user() ) {
			if(  option[ OptionKeys::loops::superimpose_native ]()  ){
				core::scoring::superimpose_pose( native_pose_super, pose, atom_map );
			}
			setPoseExtraScores( pose, "brlx_rms",   core::scoring::native_CA_rmsd(native_pose, pose ) );
			setPoseExtraScores( pose, "brlx_corerms", native_loop_core_CA_rmsd(native_pose, pose, loops, corelength ) );
			setPoseExtraScores( pose, "brlx_looprms",  loops::loop_rmsd(native_pose_super, pose, loops ) );
			setPoseExtraScores( pose, "brlx_loopcarms",  loops::loop_rmsd(native_pose_super, pose, loops,true ) );
		}

		if ( !checkpoints_.recover_checkpoint( pose, curr_job_tag, "relax", true , true) ) {
			if ( debug ) pose.dump_pdb(curr_job_tag + "_before_relax.pdb");
			if ( intermedrelax() == "relax" ) {
				relax::RelaxProtocolBaseOP relax_prot = relax::generate_relax_from_cmd();
				relax_prot->set_current_tag( curr_job_tag );
				relax_prot->apply( pose );
			} else if (( intermedrelax() == "fastrelax" ) ||
			           ( intermedrelax() == "seqrelax" )){
				relax::FastRelax seqrelax( fa_scorefxn_, option[ OptionKeys::relax::sequence_file ]() );
				seqrelax.set_current_tag( curr_job_tag );
				seqrelax.apply( pose );
			}
			checkpoints_.checkpoint( pose, curr_job_tag, "relax", true);
		}
		checkpoints_.debug( curr_job_tag, "relax", (*fa_scorefxn_)( pose ) );
	} // intermediate relax the structure

	///////////////////////////////////////////////////////////////////////////////////////
	////
	////  Loop refine (fullatom type loop modelling )
	////
	////

	if ( refine() != "no" ) {
		TR << "====================================================================================" << std::endl;
		TR << "===" << std::endl;
		TR << "===   Refine " << std::endl;
		TR << "===" << std::endl;

		long starttime = time(NULL);

		setPoseExtraScores( pose, "bref_irms",  core::scoring::CA_rmsd( start_pose, pose ) );
		if ( option[ in::file::native ].user() ) {
			if(  option[ OptionKeys::loops::superimpose_native ]()  ){
				core::scoring::superimpose_pose( native_pose_super, pose, atom_map );
			}
			setPoseExtraScores( pose, "bref_rms",   core::scoring::native_CA_rmsd(native_pose, pose ) );
			setPoseExtraScores( pose, "bref_corerms", native_loop_core_CA_rmsd(native_pose, pose, loops, corelength ) );
			setPoseExtraScores( pose, "bref_looprms",  loops::loop_rmsd(native_pose_super, pose, loops ) );
		}

		core::kinematics::FoldTree f_new, f_orig=pose.fold_tree();
		if ( refine() == "refine_kic" ) {
			loops::fold_tree_from_loops( pose, loops, f_new, true /* include terminal cutpoints */);
		}
		else {
			loops::fold_tree_from_loops( pose, loops, f_new);
		}
		pose.fold_tree( f_new );
		TR << "fold_tree_before_refine " << pose.fold_tree() << std::endl;

		if ( !checkpoints_.recover_checkpoint( pose, curr_job_tag, "refine", true , true) ) {

			if ( debug ) pose.dump_pdb(curr_job_tag + "_before_refine.pdb");
			if ( refine() == "refine_ccd" ) {
				loops::LoopMover_Refine_CCD refine_ccd( loops, fa_scorefxn_ );
				refine_ccd.set_native_pose( new core::pose::Pose ( native_pose ) );
				refine_ccd.apply( pose );
			} else
			if ( refine() == "refine_kic" ) {
				//loops.remove_terminal_loops( pose );
				loops::LoopMover_Refine_KIC refine_kic( loops, fa_scorefxn_ );
				refine_kic.set_native_pose( new core::pose::Pose ( native_pose ) );
				refine_kic.apply( pose );
			}
			if ( debug ) pose.dump_pdb(curr_job_tag + "_after_refine.pdb");
			checkpoints_.checkpoint( pose, curr_job_tag, "refine", true);
			// need to get the chainbreak score before the cutpoint variants are removed
			setPoseExtraScores(
				pose, "final_chainbreak",
				pose.energies().total_energies()[ core::scoring::chainbreak ]
			);
		}

		checkpoints_.debug( curr_job_tag, "refine", (*fa_scorefxn_)( pose ) );

		loops::remove_cutpoint_variants( pose );
		// restore simple fold tree
		pose.fold_tree( f_orig );

		endtime = time(NULL);

		TR << "Refinetime: " << endtime - starttime << std::endl;

	} // if ( refine != "no" )

	//////////////////////////////////////////////////////////////////////////////////////
	////
	////  Maybe idealize the structure before relax ?
	////
	////
	if ( option[ OptionKeys::loops::idealize_after_loop_close ].user() ) {

		if ( debug ){
			pose.dump_pdb(curr_job_tag + "_before_idealize.pdb");
		}
		protocols::idealize::IdealizeMover idealizer;
		idealizer.fast( false );
		idealizer.apply( pose );
		if ( debug ){
			pose.dump_pdb(curr_job_tag + "_after_idealize.pdb");
		}
	}


	//////////////////////////////////////////////////////////////////////////////////////
	////
	////  Relax the structure
	////
	////

	if ( relax() != "no" ) {
		TR << "====================================================================================" << std::endl;
		TR << "===" << std::endl;
		TR << "===   Relax  " << std::endl;
		TR << "===" << std::endl;

		core::kinematics::FoldTree f_new, f_orig=pose.fold_tree();
		if ( option[ OptionKeys::loops::relax_with_foldtree ].user() ){
			loops::fold_tree_from_loops( pose, loops, f_new );
			pose.fold_tree( f_new );
			loops::add_cutpoint_variants( pose );
			fa_scorefxn_->set_weight( core::scoring::chainbreak,        option[ OptionKeys::relax::chainbreak_weight ]() );
			fa_scorefxn_->set_weight( core::scoring::linear_chainbreak, option[ OptionKeys::relax::linear_chainbreak_weight ]() );
			fa_scorefxn_->set_weight( core::scoring::overlap_chainbreak, option[ OptionKeys::relax::overlap_chainbreak_weight ]() );
		}


		TR << pose.fold_tree() << std::endl;
		setPoseExtraScores( pose, "brlx_irms",  core::scoring::CA_rmsd( start_pose, pose ) );
		if ( option[ in::file::native ].user() ) {
			if(  option[ OptionKeys::loops::superimpose_native ]()  ){
				core::scoring::superimpose_pose( native_pose_super, pose, atom_map );
			}
			setPoseExtraScores( pose, "brlx_rms",   core::scoring::native_CA_rmsd(native_pose, pose ) );
			setPoseExtraScores( pose, "brlx_corerms", native_loop_core_CA_rmsd(native_pose, pose, loops, corelength ) );
			setPoseExtraScores( pose, "brlx_looprms",  loops::loop_rmsd(native_pose_super, pose, loops ) );
		}

		if ( !checkpoints_.recover_checkpoint( pose, curr_job_tag, "relax", true , true) ) {
			if ( debug ) pose.dump_pdb(curr_job_tag + "_before_relax.pdb");
			if ( relax() == "relax" ) {
				relax::RelaxProtocolBaseOP relax_prot = relax::generate_relax_from_cmd();
				relax_prot->set_current_tag( curr_job_tag );
				relax_prot->apply( pose );
			} else if (( relax() == "fastrelax" ) ||
			           ( relax() == "seqrelax" )) {
				relax::FastRelax seqrelax( fa_scorefxn_, option[ OptionKeys::relax::sequence_file ]() );
				seqrelax.set_current_tag( curr_job_tag );
				seqrelax.apply( pose );
			} else if ( relax() == "minirelax" ) {
				protocols::relax::MiniRelax mini_relax( fa_scorefxn_ );
				mini_relax.set_current_tag( curr_job_tag );
				mini_relax.apply( pose );
			}
			checkpoints_.checkpoint( pose, curr_job_tag, "relax", true);
		}
		checkpoints_.debug( curr_job_tag, "relax", (*fa_scorefxn_)( pose ) );

		// restore simple fold tree
		loops::remove_cutpoint_variants( pose );
		pose.fold_tree( f_orig );

		fa_scorefxn_->set_weight( core::scoring::chainbreak, 0.0 );
		fa_scorefxn_->set_weight( core::scoring::linear_chainbreak, 0.0 );
		fa_scorefxn_->set_weight( core::scoring::overlap_chainbreak, 0.0 );

		if ( debug ) pose.dump_pdb(curr_job_tag + "_after_relax.pdb");

		if ( option[ OptionKeys::loops::final_clean_fastrelax ]() ) {
			fa_scorefxn_->set_weight( core::scoring::atom_pair_constraint, 0.0 );
			fa_scorefxn_->set_weight( core::scoring::constant_constraint, 0.0 );
			fa_scorefxn_->set_weight( core::scoring::coordinate_constraint, 0.0 );
			fa_scorefxn_->set_weight( core::scoring::angle_constraint, 0.0 );
			fa_scorefxn_->set_weight( core::scoring::dihedral_constraint, 0.0 );
			fa_scorefxn_->set_weight( core::scoring::big_bin_constraint, 0.0 );
			if ( !checkpoints_.recover_checkpoint( pose, curr_job_tag, "ffrelax", true , true ) ) {
				relax::FastRelax fast_relax( fa_scorefxn_ );
				fast_relax.set_current_tag( curr_job_tag );
				fast_relax.apply( pose );
				checkpoints_.checkpoint( pose, curr_job_tag, "ffrelax", true);
			}
			checkpoints_.debug( curr_job_tag, "ffrelax", (*fa_scorefxn_)( pose ) );
		}
	} // relax the structure

	//////////////////////////////////////////////////////////////////////////////////////
	////
	////   Statistics
	////

	TR << "====================================================================================" << std::endl;
	TR << "===" << std::endl;
	TR << "===  Getting Statistics " << std::endl;
	TR << "===" << std::endl;
	TR << "===" << std::endl;

	setPoseExtraScores(
		pose, "irms",  core::scoring::CA_rmsd( start_pose, pose )
	);

	if ( option[ in::file::native ].user() ) {
		if(  option[ OptionKeys::loops::superimpose_native ]()  ){
			core::scoring::superimpose_pose( native_pose_super, pose, atom_map );
		}
		setPoseExtraScores(	pose, "rms",     core::scoring::native_CA_rmsd( native_pose, pose ));
		setPoseExtraScores(	pose, "looprms", loops::loop_rmsd( native_pose_super, pose, loops, false /*CA_only*/, true /*bb_only*/ ));
		setPoseExtraScores(	pose, "loop_heavy_rms", loops::loop_rmsd( native_pose_super, pose, loops, false /*CA_only*/, false /*bb_only*/ ));
		setPoseExtraScores(	pose, "loopcarms", loops::loop_rmsd( native_pose_super, pose, loops, true ));
		setPoseExtraScores(	pose, "corerms", native_loop_core_CA_rmsd( native_pose, pose, loops, corelength )	);
		setPoseExtraScores( pose, "corelen", corelength );
		//if( pose.is_fullatom() ) addScoresForLoopParts( pose, loops, (*fa_scorefxn_), native_pose, 10 );
		//else                     addScoresForLoopParts( pose, loops, (*cen_scorefxn_), native_pose, 10 );
	}

	core::Real final_score;

	if ( debug ) pose.dump_pdb(curr_job_tag + "_before_final_rescore.pdb");

	if (fullatom_output) final_score = (*fa_scorefxn_)(pose);  // may include constraint score
	else                 final_score = (*cen_scorefxn_)(pose); // may include constraint score

	core::pose::setPoseExtraScores(
		pose, std::string("final_looprelax_score"), final_score
	);

} // LoopRelaxMover

std::string
LoopRelaxMover::get_name() const {
	return "LoopRelaxMover";
}

void LoopRelaxMover::set_defaults_() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	cmd_line_csts_   = true;
	copy_sidechains_ = true;
	remodel_         = option[ OptionKeys::loops::remodel ]();
	intermedrelax_   = option[ OptionKeys::loops::intermedrelax ]();
	refine_          = option[ OptionKeys::loops::refine ]() ;
	relax_           = option[ OptionKeys::loops::relax ]();
	n_rebuild_tries_ = 1;
	rebuild_filter_  = 999.0;

	// use score4L by default (will be symm if needed)
	cen_scorefxn_ = get_cen_scorefxn();
	cen_scorefxn_->set_weight( core::scoring::chainbreak, 1.0*10.0/3.0 );

	// get cmd line scorefxn by default (will be symm if needed)
	fa_scorefxn_ = get_fa_scorefxn();
}

} // namespace loops
} // namespace protocols
