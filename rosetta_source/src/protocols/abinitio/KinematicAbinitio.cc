// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file FoldConstraints.cc
/// @brief ab-initio fragment assembly protocol for proteins under the influence of contraints (e.g., NOE)
/// @detailed
/// @author Oliver Lange
///


// Unit Headers
#include <protocols/abinitio/KinematicAbinitio.hh>

// Package Headers
#include <protocols/jumping/JumpSetup.hh>
#include <protocols/jumping/JumpSample.hh> //for diagnosis output
#include <protocols/abinitio/ResolutionSwitcher.hh>
#include <protocols/basic_moves/FragmentMover.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/ScoreFunction.hh>

#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>

#include <core/types.hh>

#include <core/fragment/OrderedFragSet.hh>
// AUTO-REMOVED #include <core/fragment/FragSetCollection.hh>
#include <core/fragment/FragmentIO.hh>

#include <core/kinematics/MoveMap.hh>

#include <protocols/moves/BackboneMover.hh>
#include <protocols/basic_moves/WobbleMover.hh>
#include <protocols/moves/MoverContainer.hh>
// AUTO-REMOVED #include <protocols/moves/RepeatMover.hh>

// Auto-header: duplicate removed #include <protocols/jumping/JumpSetup.hh>
#include <protocols/jumping/PairingLibrary.hh>
#include <protocols/jumping/util.hh>

#include <protocols/loops/SlidingWindowLoopClosure.hh>
#include <protocols/loops/Exceptions.hh>

// AUTO-REMOVED #include <protocols/evaluation/RmsdEvaluator.hh>
// AUTO-REMOVED #include <protocols/evaluation/ScoreEvaluator.hh>
#include <protocols/evaluation/JumpEvaluator.hh>
// AUTO-REMOVED #include <protocols/evaluation/EvaluatedTrialMover.hh>
#include <protocols/evaluation/ConstraintEvaluator.hh>

// Auto-header: duplicate removed #include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/SkipViolFunc.hh>

// Auto-header: duplicate removed #include <core/kinematics/MoveMap.hh>
#include <core/fragment/FragSet.hh>
// Auto-header: duplicate removed #include <core/fragment/OrderedFragSet.hh>


// ObjexxFCL Headers
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>

// Utility headers
#include <numeric/random/random.hh>
// AUTO-REMOVED #include <utility/io/izstream.hh> //for cheats
#include <utility/io/ozstream.hh> //for dump_frags

#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/jumps.OptionKeys.gen.hh>
#include <basic/options/keys/resample.OptionKeys.gen.hh>
#include <basic/options/keys/fold_cst.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/keys/in.OptionKeys.gen.hh>

#include <protocols/jobdist/Jobs.hh>
#include <protocols/idealize/IdealizeMover.hh>

//// C++ headers
//#include <cstdlib>
//#include <string>
#include <fstream>

#include <platform/types.hh>
#include <core/chemical/AA.hh>
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
#include <core/conformation/symmetry/SymDof.hh>
#include <core/conformation/symmetry/SymSlideInfo.fwd.hh>
#include <core/conformation/symmetry/SymSlideInfo.hh>
#include <core/conformation/symmetry/SymmData.fwd.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/fragment/BaseCacheUnit.hh>
#include <core/fragment/ConstantLengthFragSet.fwd.hh>
#include <core/fragment/FragData.fwd.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/FragID.fwd.hh>
#include <core/fragment/FragID.hh>
#include <core/fragment/FragID_Iterator.fwd.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/fragment/Frame.fwd.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameIterator.fwd.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/FrameIteratorWorker_.fwd.hh>
#include <core/fragment/FrameIteratorWorker_.hh>
#include <core/fragment/FrameList.fwd.hh>
#include <core/fragment/FrameList.hh>
#include <core/fragment/JumpingFrame.fwd.hh>
#include <core/fragment/OrderedFragSet.fwd.hh>
#include <core/fragment/SecondaryStructure.fwd.hh>
#include <core/fragment/SecondaryStructure.hh>
#include <core/fragment/SingleResidueFragData.fwd.hh>
#include <core/fragment/SingleResidueFragData.hh>
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
#include <core/id/DOF_ID_Range.fwd.hh>
#include <core/id/DOF_ID_Range.hh>
#include <core/id/JumpID.fwd.hh>
#include <core/id/JumpID.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/NamedStubID.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <core/id/SequenceMapping.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/TorsionID.hh>
#include <core/id/TorsionID_Range.fwd.hh>
#include <core/id/TorsionID_Range.hh>
#include <core/id/types.hh>
#include <core/io/raw_data/DisulfideFile.hh>
#include <core/io/silent/EnergyNames.hh>
#include <core/io/silent/ProteinSilentStruct.tmpl.hh>
#include <core/io/silent/SharedSilentData.hh>
#include <core/io/silent/SilentEnergy.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/io/silent/SilentStruct.fwd.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/silent.fwd.hh>
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
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MinimizerMapBase.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/kinematics/RT.fwd.hh>
#include <core/kinematics/RT.hh>
#include <core/kinematics/ResidueCoordinateChangeList.fwd.hh>
#include <core/kinematics/ShortestPathInFoldTree.fwd.hh>
#include <core/kinematics/ShortestPathInFoldTree.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/types.hh>
#include <core/kinematics/tree/Atom.fwd.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/pose/MiniPose.fwd.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/util.hh>
#include <core/pose/datacache/ObserverCache.fwd.hh>
#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
#include <core/pose/signals/ConformationEvent.fwd.hh>
#include <core/pose/signals/DestructionEvent.fwd.hh>
#include <core/pose/signals/EnergyEvent.fwd.hh>
#include <core/pose/signals/GeneralEvent.fwd.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/EnergyGraph.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/LREnergyContainer.fwd.hh>
#include <core/scoring/MinimizationData.fwd.hh>
#include <core/scoring/MinimizationGraph.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionInfo.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/Constraints.fwd.hh>
#include <core/scoring/constraints/Constraints.hh>
#include <core/scoring/constraints/DOF_Constraint.fwd.hh>
#include <core/scoring/constraints/DOF_Constraint.hh>
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
#include <core/scoring/methods/TwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/WholeStructureEnergy.fwd.hh>
#include <protocols/Protocol.fwd.hh>
#include <protocols/Protocol.hh>
#include <protocols/abinitio/ClassicAbinitio.hh>
#include <protocols/abinitio/FoldConstraints.fwd.hh>
#include <protocols/abinitio/FoldConstraints.hh>
#include <protocols/abinitio/KinematicAbinitio.fwd.hh>
#include <protocols/abinitio/KinematicControl.fwd.hh>
#include <protocols/abinitio/KinematicControl.hh>
#include <protocols/abinitio/MaxSeqSepConstraintSet.fwd.hh>
#include <protocols/abinitio/MaxSeqSepConstraintSet.hh>
#include <protocols/abinitio/ResolutionSwitcher.fwd.hh>
#include <protocols/basic_moves/FragmentMover.fwd.hh>
#include <protocols/basic_moves/SmoothFragmentMover.fwd.hh>
#include <protocols/basic_moves/SmoothFragmentMover.hh>
#include <protocols/checkpoint/CheckPointer.fwd.hh>
#include <protocols/checkpoint/CheckPointer.hh>
#include <protocols/evaluation/ConstraintEvaluator.fwd.hh>
#include <protocols/evaluation/PoseEvaluator.fwd.hh>
#include <protocols/evaluation/PoseEvaluator.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/idealize/IdealizeMover.fwd.hh>
#include <protocols/jobdist/Jobs.fwd.hh>
#include <protocols/jumping/JumpSetup.fwd.hh>
#include <protocols/jumping/PairingsList.fwd.hh>
#include <protocols/jumping/PairingsList.hh>
#include <protocols/loops/Loop.fwd.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.fwd.hh>
#include <protocols/loops/SlidingWindowLoopClosure.fwd.hh>
#include <protocols/moves/BackboneMover.fwd.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/moves/MetropolisHastingsMover.fwd.hh>
#include <protocols/moves/MinMover.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/MonteCarloStatus.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.fwd.hh>
#include <protocols/moves/MoverStatistics.hh>
#include <protocols/moves/MoverStatus.hh>
#include <protocols/moves/ThermodynamicMover.fwd.hh>
#include <protocols/moves/ThermodynamicMover.hh>
#include <protocols/moves/TrialMover.fwd.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/mc_convergence_checks/ConvergenceCheck.fwd.hh>
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
#include <utility/excn/EXCN_Base.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/file/FileName.fwd.hh>
#include <utility/file/FileName.hh>
#include <utility/file/PathName.fwd.hh>
#include <utility/file/PathName.hh>
#include <utility/file/gzip_util.hh>
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
#include <numeric/model_quality/rms.hh>
#include <numeric/random/random.fwd.hh>
#include <numeric/random/uniform.fwd.hh>
#include <numeric/random/uniform.hh>
#include <ObjexxFCL/Dimension.fwd.hh>
#include <ObjexxFCL/Dimension.hh>
#include <ObjexxFCL/DimensionExpression.hh>
#include <ObjexxFCL/DynamicIndexRange.fwd.hh>
#include <ObjexxFCL/DynamicIndexRange.hh>
#include <ObjexxFCL/FArray.fwd.hh>
#include <ObjexxFCL/FArray.hh>
#include <ObjexxFCL/FArray1.fwd.hh>
#include <ObjexxFCL/FArray1.hh>
#include <ObjexxFCL/FArray1A.fwd.hh>
#include <ObjexxFCL/FArray1A.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray1P.fwd.hh>
#include <ObjexxFCL/FArray1P.hh>
#include <ObjexxFCL/FArray2.fwd.hh>
#include <ObjexxFCL/FArray2.hh>
#include <ObjexxFCL/FArray2A.fwd.hh>
#include <ObjexxFCL/FArray2A.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <ObjexxFCL/FArray2P.fwd.hh>
#include <ObjexxFCL/FArray2P.hh>
#include <ObjexxFCL/FArrayInitializer.fwd.hh>
#include <ObjexxFCL/FArrayInitializer.hh>
#include <ObjexxFCL/FArraySection.fwd.hh>
#include <ObjexxFCL/FArraySection.hh>
#include <ObjexxFCL/FArrayTraits.fwd.hh>
#include <ObjexxFCL/FArrayTraits.hh>
#include <ObjexxFCL/Fstring.fwd.hh>
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
#include <ObjexxFCL/StaticIndexRange.fwd.hh>
#include <ObjexxFCL/StaticIndexRange.hh>
#include <ObjexxFCL/TypeTraits.hh>
#include <ObjexxFCL/byte.fwd.hh>
#include <ObjexxFCL/char.functions.hh>
#include <ObjexxFCL/proxy_const_assert.hh>
#include <ObjexxFCL/ubyte.fwd.hh>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
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
#include <sstream>
#include <string>
#include <time.h>
#include <typeinfo>
#include <utility>
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
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/prof.hh>
#include <boost/algorithm/string/erase.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/functional/hash.hpp>
#include <boost/lexical_cast.hpp>
#include <zlib/zlib.h>
#include <zlib/zutil.h>

//Auto using namespaces
namespace std { } using namespace std; // AUTO USING NS
namespace ObjexxFCL { } using namespace ObjexxFCL; // AUTO USING NS
namespace ObjexxFCL { namespace fmt { } } using namespace ObjexxFCL::fmt; // AUTO USING NS
//Auto using namespaces end


// option key includes



static numeric::random::RandomGenerator RG(198243);  // <- Magic number, do not change it!
static basic::Tracer tr("protocols.abinitio",basic::t_info);

using namespace core;
using scoring::constraints::ConstraintSet;
using scoring::constraints::ConstraintSetOP;
using kinematics::ShortestPathInFoldTree;

//@detail call this static routine before core::init to register the options
void protocols::abinitio::KinematicAbinitio::register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	Parent::register_options();

	option.add_relevant( resample::silent );
	option.add_relevant( resample::tag );
	option.add_relevant( resample::stage2 );
	option.add_relevant( resample::min_max_start_seq_sep );
 	option.add_relevant( jumps::bb_moves );
	option.add_relevant( jumps::no_wobble );
	option.add_relevant( jumps::no_shear );
	option.add_relevant( jumps::no_sample_ss_jumps );
	option.add_relevant( jumps::invrate_jump_move );
	option.add_relevant( jumps::chainbreak_weight_stage1 );
	option.add_relevant( jumps::chainbreak_weight_stage2 );
	option.add_relevant( jumps::chainbreak_weight_stage3 );
	option.add_relevant( jumps::chainbreak_weight_stage4 );
	option.add_relevant( jumps::increase_chainbreak );
	option.add_relevant( jumps::ramp_chainbreaks );
	option.add_relevant( jumps::overlap_chainbreak );
	// constraints/chainbreaks are phased in regarding their seq-sep. it grows from 15 (stage2 ) to MAX at end of stage4
	// seqsep_accelerate > 0 will increase MAX such that constraints/chainbreak are enforced earlier. seqsep_accelerate = 0 means
	// that MAX is exactly the longest seq-sep found in a given fold-tree.
	option.add_relevant( jumps::sep_switch_accelerate );
	option.add_relevant( jumps::dump_frags );
	option.add_relevant( OptionKeys::loops::idealize_after_loop_close );
	option.add_relevant( fold_cst::constraint_skip_rate );

}

namespace protocols {
namespace abinitio {

using namespace basic::options;
using namespace basic::options::OptionKeys;

KinematicAbinitio::~KinematicAbinitio() {}

//@detail c'stor
KinematicAbinitio::KinematicAbinitio(
	basic_moves::FragmentMoverOP brute_move_small,
	basic_moves::FragmentMoverOP brute_move_large,
	basic_moves::FragmentMoverOP smooth_move_small,
	int dummy /* otherwise the two constructors are ambigous */
) :
	FoldConstraints ( brute_move_small, brute_move_large, smooth_move_small, dummy ),
	bRampChainbreaks_( true )
{
  BaseClass::type( "KinematicAbinitio" );
	full_constraint_set_=NULL;
}


KinematicAbinitio::KinematicAbinitio(
   core::fragment::FragSetCOP fragset3mer,
   core::fragment::FragSetCOP fragset9mer,
   core::kinematics::MoveMapCOP movemap
) : FoldConstraints ( fragset3mer, fragset9mer, movemap ),
    bRampChainbreaks_( true )
{
  BaseClass::type( "KinematicAbinitio" );
	full_constraint_set_=NULL;
}

void
KinematicAbinitio::set_default_scores() {
  Parent::set_default_scores();
  // set low chainbreak energy for stage 1 + 2
  Real chainbreak_score_1 = option[ jumps::chainbreak_weight_stage1 ]();
  Real chainbreak_score_2 = option[ jumps::chainbreak_weight_stage2 ]();
  Real chainbreak_score_3 = option[ jumps::chainbreak_weight_stage3 ]();
  Real chainbreak_score_4 = option[ jumps::chainbreak_weight_stage4 ]();
  if ( !bRampChainbreaks_ ) {
    set_score_weight( scoring::linear_chainbreak, chainbreak_score_1, STAGE_1 );
    set_score_weight( scoring::linear_chainbreak, chainbreak_score_2, STAGE_2 );
    set_score_weight( scoring::linear_chainbreak, chainbreak_score_3, STAGE_3a );
    set_score_weight( scoring::linear_chainbreak, chainbreak_score_3, STAGE_3b );
    set_score_weight( scoring::linear_chainbreak, chainbreak_score_4, STAGE_4 );
  }
}

//@brief read out cmd-line options
void KinematicAbinitio::set_default_options() {
  Parent::set_default_options();
	bRampChainbreaks_ = option[ jumps::ramp_chainbreaks ]; //default is true
  bOverlapChainbreaks_ = option[ jumps::overlap_chainbreak ];
}

void KinematicAbinitio::replace_scorefxn( core::pose::Pose& pose, StageID stage, core::Real intra_stage_progress ) {
	Parent::replace_scorefxn( pose, stage, intra_stage_progress );
	scoring::ScoreFunction scorefxn( current_scorefxn() );
	kinematics().add_score_weights( scorefxn, 1.0 * stage / 2.0 + intra_stage_progress * 0.2 );
	current_scorefxn( scorefxn );
}

bool
KinematicAbinitio::prepare_stage1( core::pose::Pose &pose ) {
	//we set score term in here: we minimize in this one... need to set all score terms before that
	bool success = Parent::prepare_stage1( pose );
	if ( bRampChainbreaks_ ) {
    Real const setting( 0.25 / 3 * option[ jumps::increase_chainbreak ] );
    set_score_weight( scoring::linear_chainbreak, 0.0, ALL_STAGES );
    set_score_weight( scoring::linear_chainbreak, setting, STAGE_2);
  }
	return success;
}

bool
KinematicAbinitio::prepare_stage2( core::pose::Pose &pose ) {
  return Parent::prepare_stage2( pose );
}

bool
KinematicAbinitio::prepare_stage3( core::pose::Pose &pose ) {
  return Parent::prepare_stage3( pose );
}

bool
KinematicAbinitio::prepare_loop_in_stage3(
		core::pose::Pose& pose,
		Size iteration, /* loop_iteration*/
		Size total /* total_iterations */
) {
	Real progress( 1.0* iteration/total );
  if ( bRampChainbreaks_ ) {
    Real const fact(  progress * 1.0/3  * option[ jumps::increase_chainbreak ]);
    //score_stage3a ( score2 )
    set_score_weight( scoring::linear_chainbreak, 2.5 * fact, STAGE_3a );
    //score_stage3b ( score5 )
    set_score_weight( scoring::linear_chainbreak, 0.5 * fact, STAGE_3b );
  }
  return Parent::prepare_loop_in_stage3( pose, iteration, total ); //minimization happens here!
}


bool
KinematicAbinitio::prepare_loop_in_stage4(
  core::pose::Pose& pose,
  Size iteration, /* loop_iteration*/
  Size total /* total_iterations */
) {
	Real progress( 1.0* iteration/total );
  if ( bRampChainbreaks_ ) {
		Real const setting( ( 1.5*progress+2.5 ) * ( 1.0/3) * option[ jumps::increase_chainbreak ]);
    set_score_weight( scoring::linear_chainbreak, setting, STAGE_4 );
    set_current_weight( scoring::linear_chainbreak, setting );
    if ( bOverlapChainbreaks_ ) {
      set_score_weight( scoring::overlap_chainbreak, progress, STAGE_4 );
      set_current_weight( scoring::overlap_chainbreak, progress );
    }
  }
  return Parent::prepare_loop_in_stage4( pose, iteration, total );
}


void
KinematicAbinitio::set_max_seq_sep( core::pose::Pose& pose, Size max_dist ) {
  //if not already set, --> fix the chainbreaks dependent on distance
  kinematics().add_chainbreak_variants( pose, max_dist, constraints().shortest_path() );
  mc().reset( pose ); //necessary to avoid that new chainbreaks are immediatly removed by copying back the old pose
  Parent::set_max_seq_sep( pose, max_dist );
}

///@details the native_pose_ is used to determine orientation and pleating of jumps
// void
// KinematicAbinitio::set_native_pose( core::pose::Pose const& native_pose ) {
//   native_pose_ = new core::pose::Pose( native_pose );
// }

void
KinematicAbinitio::dump_jump_log( core::pose::Pose& pose, std::string const& filename ) {
  // diagnostics for jumps
  if ( get_native_pose() ) {
		pose::Pose native_pose( *get_native_pose() );
    native_pose.fold_tree( kinematics().sampling_fold_tree() );
    core::io::silent::ProteinSilentStruct pss;
    pss.fill_struct( pose, "jump_log" );
    evaluate_pose( pose, "jump_log", pss );

    utility::io::ozstream out( filename , std::ios_base::out | std::ios_base::app );
		using namespace ObjexxFCL::fmt;
    out << get_current_tag() << " "
				<< RJ(10, pss.get_energy ( "rms" ) ) << " "
			//				<< RJ(10, pss.get_energy ( "rms_native" ) ) << " "
				<< RJ(10, pss.get_energy ( "score" ) ) << " ";
    out << RJ(5, kinematics().sampling_fold_tree().num_jump() ) << " ";
		for ( Size i = 1; i<=kinematics().sampling_fold_tree().num_jump(); i++ )  {
			out << RJ(10, evaluation::JumpEvaluator( native_pose, i ).apply( pose ) ) << " ";
		}
		out << jumping::JumpSample( pose.fold_tree() ) << " " << std::endl;
  }
}

void
KinematicAbinitio::apply( core::pose::Pose& pose ) {
	if ( option[ resample::silent ].user() ) {
		std::string const filename( option[ resample::silent ]() );
		std::string tag("");
		if ( option[ resample::tag ].user()  ) {
			tag = option[ resample::tag ]();
		} else {
			jobdist::BasicJobCOP job = get_current_job();
			if ( job ) {
				tag = job->input_tag();
			}
		}

		if ( tag == "" ) {
			utility_exit_with_message("no resample tag found -- supply either via -resample:tag or via JobDistributor input_tag");
		}

		//retrieve pose
		io::silent::SilentFileData sfd;
		if ( !sfd.read_file( filename ) ) {
			utility_exit_with_message( "resampling for " + tag + " failed, problem reading silent file " + filename );
		}

		if ( !sfd.has_tag( tag ) ) {
			utility_exit_with_message( "resampling for " + tag + " failed, not found in silent file " + filename );
		}
		ConstraintSetOP cst_set = new ConstraintSet( *pose.constraint_set() ); //keep constraints from pose, (make new to remove constant)
		sfd.get_structure( tag ).fill_pose( pose ); //fill pose and fold-tree from file
		pose.constraint_set( cst_set ); // put constraint set back into place

		//retrieve fold_tree from pose
		KinematicControlOP recovered_control = new KinematicControl( kinematics() );
		if ( option[ resample::jumps ] ) {
			using namespace jumping;
			using namespace fragment;
			JumpSample target_jumps( pose.fold_tree() );
			target_jumps.steal_orientation_and_pleating( pose );
			tr.Warning << "enable JUMP MOVES in resamplin --- these are taken from SS-library"
								 <<" no matter where original structures came from" << std::endl;
			FragSetOP jump_frags = new OrderedFragSet;
			//generate fragments from all homologes
			core::fragment::FrameList jump_frames;
			target_jumps.generate_jump_frames( jump_frames, kinematics().movemap() );
			PairingsList library_pairings;
			for ( FrameList::iterator jump_frame = jump_frames.begin();
						jump_frame != jump_frames.end(); ++jump_frame ) {
				Pairing target_pairing( target_jumps.get_pairing( (*jump_frame)->start(), (*jump_frame)->stop() ) );
				library_pairings.push_back( target_pairing );
			}
			//fill remaining frames from ss-library
			jumping::StandardPairingLibrary::get_instance()->
				generate_jump_frags(
														library_pairings,
														kinematics().movemap(),
														true /*with Torsion*/,
														*jump_frags
				);
			using namespace protocols::basic_moves;
			basic_moves::ClassicFragmentMoverOP jump_mover = new ClassicFragmentMover( jump_frags, kinematics().movemap_ptr() );
			jump_mover->type( "JumpMoves" );
			jump_mover->set_check_ss( false ); // this doesn't make sense with jump fragments
			jump_mover->enable_end_bias_check( false ); //no sense for discontinuous fragments
			recovered_control->set_jump_mover( jump_mover );
		} else {
			tr.Warning << "disable JUMP MOVES in resampling... probably not much of difference" << std::endl;
			recovered_control->set_jump_mover( NULL ); // this mover will not have the correct jumps ... but jump-moves seem ineffective in stage3 an
		}
		recovered_control->set_sampling_fold_tree( pose.fold_tree() );
		set_kinematics( recovered_control );
		if ( !option[ resample::stage1 ] ) { //default true
			set_skip_stage1( true );
		}
		if ( !option[ resample::stage2 ] ) {
			set_skip_stage2( true );
		}
		if ( option[ resample::min_max_start_seq_sep ].user() ) {
			Real const min_sep( option[ resample::min_max_start_seq_sep ]()[ 1 ] );
			Real const max_sep( option[ resample::min_max_start_seq_sep ]()[ 2 ] );
			Real r = RG.uniform();
			Real val = r/(max_sep-min_sep)+min_sep;
			set_seq_sep_stage1( val );
			if ( val > 0.5 ) set_seq_sep_stage3( val ); //hard-coded number replace later
		}
	}

	//	evaluation::MetaPoseEvaluatorOP ori_evaluator = evaluator();
	tr.Debug << "set movemap from KinematicControl " << std::endl;
	set_movemap( kinematics().movemap_ptr() );
	//set fold-tree, apply a fragment from each jump_mover frame once
	kinematics().prepare_pose_for_sampling( pose );
	tr.Debug << "KinematicControl is implemented: ready to rumble " << std::endl;

	// set max_seq_separation fudge -factor such that distance contraints are not enforced far too late
	// let's have it depend on the number of jumps. -- this is totally arbitrary
	// maybe should check if it is actually better to keep the ramping speed fixed ( independen of max_possible_seq_separation ).
	max_seq_sep_fudge( 1.0 + kinematics().sampling_fold_tree().num_jump() * option[ jumps::sep_switch_accelerate ] );

	// diagnostics for jumps
	current_scorefxn()( pose );
	dump_jump_log( pose, "jumps_pre.log");
	if ( option[ jumps::dump_frags ] ) {
		fragment::FragmentIO().write_data(  "fragset_jumps.dump", *(kinematics().jump_mover()->fragments()) );
	}

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace scoring::constraints;
	ConstraintSetOP orig_constraints( NULL );

	if ( option[ fold_cst::keep_skipped_csts ] ) {
		// if this is active we let the pose leave this class with the modulated constraint set... thus we reset it each time again
		if ( !full_constraint_set_ ) {
			full_constraint_set_ = pose.constraint_set()->clone();
		} else {
			pose.constraint_set( full_constraint_set_ );
		}
	}

	ConstraintCOPs skipped_list;
	if ( pose.constraint_set()->has_residue_pair_constraints() ) {
		evaluator()->add_evaluation( new evaluation::ConstraintEvaluator( "total_cst", *pose.constraint_set() ) );
		if ( option[ fold_cst::constraint_skip_rate ].user() ) {
			orig_constraints = pose.constraint_set()->clone();
			Real const skip_rate( option[ fold_cst::constraint_skip_rate ]() );
			tr.Info << "Skip some constraints: " << skip_rate << std::endl;
			ConstraintCOPs cst_list = orig_constraints->get_all_constraints();
			ConstraintSetOP filtered_cst = new ConstraintSet; //empty
			for ( ConstraintCOPs::const_iterator it = cst_list.begin(), eit = cst_list.end();
						it != eit; ++it ) {
				Real local_skip( 1.0 );
				if ( option[ fold_cst::violation_skip_basis ].user() ) {
					Real const basis( option[ fold_cst::violation_skip_basis ] );
					Real const base_line( option[ fold_cst::violation_skip_ignore ] );
					try {
						SkipViolFunc const& cfunc = dynamic_cast< SkipViolFunc const& >( (*it)->get_func() );
						local_skip*=1.0*(cfunc.viols()-base_line)/basis;
						tr.Trace << "ponder constraint "; (*it)->show_def( tr.Trace, pose );
						tr.Trace << "skip prob: " << skip_rate*local_skip << " computed from, viols :" << cfunc.viols() << " - " << base_line << " base: "
										 << basis << " times overall skip_rate " << skip_rate << std::endl;
						if ( local_skip > 1.0 ) local_skip = 1.0;
					} catch ( std::bad_cast ){};
				}
				Real r = RG.uniform();
				if ( r > skip_rate*local_skip ) { //keep constraint
					tr.Trace << "keep constraint";
					(*it)->show_def( tr.Trace, pose );
					tr.Trace << std::endl;
					filtered_cst->add_constraint( (*it)->clone() ); //clone things so show_violation doesn't change the SkipViolFunc numbers
				} else skipped_list.push_back( *it );
			}
			pose.constraint_set( filtered_cst );
			if ( tr.Debug.visible() ) {
				filtered_cst->show_definition( std::cout, pose );
			}
		}
	}

  // run protocol
	Parent::apply( pose );
  bool success = ( get_last_move_status() == moves::MS_SUCCESS );

  // diagnostics for jumps
  dump_jump_log( pose, "jumps.log");

	if ( success && closure_protocol_ ) {
		closure_protocol_->scorefxn( new scoring::ScoreFunction( current_scorefxn() ) );
		closure_protocol_->movemap( movemap() );
		closure_protocol_->fragments( brute_move_small()->fragments() );
		closure_protocol_->set_current_tag( get_current_tag() );
		try {
			jumping::close_chainbreaks( closure_protocol_, pose, get_checkpoints(), get_current_tag(), kinematics().final_fold_tree() );
		} catch ( loops::EXCN_Loop_not_closed& excn ) {
			set_current_tag( "C_"+get_current_tag().substr(std::min(2,(int)get_current_tag().size())) );
			set_last_move_status( moves::FAIL_RETRY );
		}
		if( option[ OptionKeys::abinitio::debug ].user() ){
			output_debug_structure( pose, "loop_closed" );
		}
		//////////////////////////////////////////////////////////////////////////////////////
		////  Maybe idealize the structure before relax ?
		if ( option[ OptionKeys::loops::idealize_after_loop_close ]()){
				protocols::idealize::IdealizeMover idealizer;
				idealizer.fast( false );
				idealizer.apply( pose );
		}
	}

	{ //dump violated constraints into special file
		if ( orig_constraints && option[ basic::options::OptionKeys::out::file::silent ].user() ) {
			ConstraintSetCOP filtered_cst = pose.constraint_set();
			ConstraintCOPs cst_list = filtered_cst->get_all_constraints();
			std::string viol_file="_viols";
			if ( get_current_job() && get_current_job()->output_file_name() != "" ) {
				viol_file = get_current_job()->output_file_name()+viol_file;
			}
			utility::io::ozstream viol_stream;
			viol_stream.open_append( viol_file );
			for ( ConstraintCOPs::const_iterator it = cst_list.begin(), eit = cst_list.end();
						it != eit; ++it ) {
				viol_stream << get_current_tag() << " ";
				if ( (*it)->show_violations( tr.Debug, pose, 1, 1.1 /*threshold*/ ) ) {
					viol_stream << "AKTIV VIOL ";
				} else {
					viol_stream << "AKTIV PASS ";
				}
				(*it)->show_def( viol_stream, pose );
			} //for

			//no the constraint that were not active
			for ( ConstraintCOPs::const_iterator it = skipped_list.begin(), eit = skipped_list.end();
						it != eit; ++it ) {
				viol_stream << get_current_tag() << " ";
				if ( (*it)->show_violations( tr.Debug, pose, 1, 1.1 /*threshold*/ ) ) {
					viol_stream << "PASSIV VIOL ";
				} else {
					viol_stream << "PASSIV PASS ";
				}
				(*it)->show_def( viol_stream, pose );
			} //for

		}
	}//scope

	if ( orig_constraints && !option[ fold_cst::keep_skipped_csts ] ) {
		pose.constraint_set( orig_constraints );
	}

	if ( pose.constraint_set()->has_residue_pair_constraints() ) {
		evaluator()->pop_back(); //remove the added evaluation called total_cst
	}

}

std::string
KinematicAbinitio::get_name() const {
	return "KinematicAbinitio";
}

moves::MoverOP KinematicAbinitio::create_jump_moves( moves::MoverOP std_mover ) {
  basic_moves::FragmentMoverOP jump_mover = kinematics().jump_mover();
	if ( !jump_mover ) return std_mover;

  moves::RandomMoverOP combi_move ( new moves::RandomMover );
  Size nfrag_moves( option[ jumps::invrate_jump_move ] );
  for ( Size i = 1; i<=nfrag_moves; i++ ) { // a simple way of changing the relative frequencies
    combi_move->add_mover( std_mover );
  }
  combi_move->add_mover( jump_mover );

  return combi_move;
}


moves::MoverOP KinematicAbinitio::create_bb_moves(
  pose::Pose &,
  moves::MoverOP std_mover,
  bool bLargeWobble,
  Real crank_up_angle )
{
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  int const nmoves ( 5 ); //this is standard in relax
  Real temp = mc().temperature();

  // setup the move objects
	core::kinematics::MoveMapOP mm_temp( new core::kinematics::MoveMap( *movemap() ) );
  moves::SmallMoverOP small_mover( new moves::SmallMover( mm_temp, temp, nmoves ) );
  small_mover->angle_max( 'H', 2.0*crank_up_angle );
  small_mover->angle_max( 'E', 2.0*crank_up_angle );
  small_mover->angle_max( 'L', 3.0*crank_up_angle );

  // setup the move objects
  moves::ShearMoverOP shear_mover( new moves::ShearMover( mm_temp, temp, nmoves ) );
  shear_mover->angle_max( 'H', 2.0*crank_up_angle*2.0 );
  shear_mover->angle_max( 'E', 2.0*crank_up_angle*2.0 );
  shear_mover->angle_max( 'L', 3.0*crank_up_angle*2.0 );

  //setup cycle
  moves::RandomMoverOP cycle( new moves::RandomMover() );
  cycle->add_mover( std_mover ); //the current frag_trial
  cycle->add_mover( small_mover );
  if ( !option[ jumps::no_shear ]() ) cycle->add_mover( shear_mover );

  if ( !option[ jumps::no_wobble ]() ) {
    // and why not throwing in a wobble ?!
    if ( bLargeWobble ) {
      cycle->add_mover( new basic_moves::WobbleMover( brute_move_large()->fragments(), movemap() ) );
    } else {
      cycle->add_mover( new basic_moves::WobbleMover( brute_move_small()->fragments(), movemap() ) );
    }
  }
  return cycle;
}

moves::TrialMoverOP KinematicAbinitio::stage1_mover( pose::Pose &, moves::TrialMoverOP trials ) {
	if ( kinematics().jump_mover() ) {
		return new moves::TrialMover( create_jump_moves( trials->mover() ), mc_ptr() );
  } else return trials;
}


moves::TrialMoverOP KinematicAbinitio::stage2_mover( pose::Pose &, moves::TrialMoverOP trials ) {
	if ( kinematics().jump_mover() ) {
		return new moves::TrialMover( create_jump_moves( trials->mover() ), mc_ptr() );
  } else return trials;
}


moves::TrialMoverOP KinematicAbinitio::stage3_mover( pose::Pose &pose, int, int, moves::TrialMoverOP trials ) {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  moves::MoverOP moves = trials->mover();
  if ( kinematics().jump_mover() ) {
    //adds jump_moves to the already existant moves
    moves = create_jump_moves( moves );
  }

  // add also bb_moves to the existent moves
  if ( option[ jumps::bb_moves ] ) {
    // okay the bb_moves are asked for. We are ON AIR !
    bool bLargeWobble( true );
    Real crank_up_angle = 6.0; //factor to make the small_moves larger: after-all we are in centroid mode
    moves = create_bb_moves( pose, moves, bLargeWobble, crank_up_angle );
  }
	moves::TrialMoverOP retval = new moves::TrialMover( moves, mc_ptr() );
	retval->keep_stats_type( moves::accept_reject );
	return retval;
  //return new moves::TrialMover( moves, mc_ptr() );
}

moves::TrialMoverOP KinematicAbinitio::stage4_mover( pose::Pose &pose, int kk, moves::TrialMoverOP trials ) {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  moves::MoverOP moves = trials->mover();
  if (  kinematics().jump_mover() && kk <= 1 ) {
    moves = create_jump_moves( moves );
  }
  if ( option[ jumps::bb_moves ] ) {
    // okay the bb_moves are asked for. We are ON AIR !
    bool bLargeWobble( kk<=1 );
    Real crank_up_angle = 5.0; //factor to make the small_moves larger: after-all we are in centroid mode
    moves = create_bb_moves( pose, moves, bLargeWobble, crank_up_angle );
  }
  return new moves::TrialMover( moves, mc_ptr() );
}



void
JumpingFoldConstraintsWrapper::apply( core::pose::Pose& pose ) {
	ResolutionSwitcher res_switch(
			pose,
			!start_from_centroid(), //input pose /*is ignored now witch calls pose.is_fullatom() instead */
			true, // starts as centroid
			true //yeah we will apply it to a centroid structure
	);

  // initialize jumping
	jumping::JumpSample current_jumps;

	Size attempts( 10 );
	do {
		current_jumps = jump_def_->create_jump_sample();
	} while ( !current_jumps.is_valid() && attempts-- );

	if ( !current_jumps.is_valid() ) {
		utility_exit_with_message( "not able to build valid fold-tree in JumpingFoldConstraints::setup_foldtree" );
	}

	//allow jumps to move
	kinematics::MoveMapOP new_movemap = new	kinematics::MoveMap( *movemap() );
	new_movemap->set_jump( true );

	core::fragment::FragSetOP jump_frags;
	jump_frags = jump_def_->generate_jump_frags( current_jumps, *new_movemap );
	using namespace protocols::basic_moves;
	basic_moves::ClassicFragmentMoverOP jump_mover = new ClassicFragmentMover( jump_frags, new_movemap );
	jump_mover->type( "JumpMoves" );
	jump_mover->set_check_ss( false ); // this doesn't make sense with jump fragments
	jump_mover->enable_end_bias_check( false ); //no sense for discontinuous fragments

	abinitio::KinematicControlOP kc = new abinitio::KinematicControl;
	kc->set_sampling_fold_tree( current_jumps.fold_tree() );
	tr.Debug << "JumpingFoldConstraintsWrapper: sampling fold_tree " << current_jumps.fold_tree() << std::endl;
	kc->set_final_fold_tree( pose.fold_tree() );
	tr.Debug << "JumpingFoldConstraintsWrapper: final fold_tree " << pose.fold_tree() << std::endl;
	kc->set_jump_mover( jump_mover );
	kc->set_movemap( new_movemap );
  // run protocol
	if ( jump_mover() && option[ jumps::no_sample_ss_jumps ] ) {
		jump_mover()->apply_at_all_positions( pose ); //make sure each jump is initialized
		kc->set_jump_mover( NULL ); //but no sampling
	}

	set_kinematics( kc );

	Parent::apply( pose );
	if ( get_last_move_status() == moves::MS_SUCCESS && !return_centroid() ) {
		res_switch.apply( pose );
		set_last_move_status( res_switch.get_last_move_status() );
	}
}

std::string
JumpingFoldConstraintsWrapper::get_name() const {
	return "JumpingFoldConstraintsWrapper";
}

//@detail c'stor
JumpingFoldConstraintsWrapper::JumpingFoldConstraintsWrapper(
	basic_moves::FragmentMoverOP brute_move_small,
	basic_moves::FragmentMoverOP brute_move_large,
	basic_moves::FragmentMoverOP smooth_move_small,
	jumping::BaseJumpSetupOP jump_def,
	int dummy /* otherwise the two constructors are ambigous */
) : KinematicAbinitio ( brute_move_small, brute_move_large, smooth_move_small, dummy ),
		jump_def_ ( jump_def )
{
	BaseClass::type( "JumpingFoldConstraintsWrapper" );
}


JumpingFoldConstraintsWrapper::JumpingFoldConstraintsWrapper(
	core::fragment::FragSetCOP fragset3mer,
	core::fragment::FragSetCOP fragset9mer,
	core::kinematics::MoveMapCOP movemap,
	jumping::BaseJumpSetupOP jump_def
) : KinematicAbinitio ( fragset3mer, fragset9mer, movemap ),
		jump_def_( jump_def )
{
	BaseClass::type( "JumpingFoldConstraintsWrapper" );
}


} //abinitio
} //protocols
