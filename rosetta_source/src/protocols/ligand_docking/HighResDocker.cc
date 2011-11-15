// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/ResfileReader.cc
/// @brief  implementation of resfile reader and its command classes
/// @author Gordon Lemmon (glemmon@gmail.com), adapted from the ResfileReader code
/// by Steven Lewis (smlewi@unc.edu) and Andrew Leaver-Fay

// Unit Headers
#include <protocols/ligand_docking/HighResDocker.hh>
#include <protocols/ligand_docking/HighResDockerCreator.hh>
#include <protocols/ligand_docking/InterfaceBuilder.hh>
#include <protocols/ligand_docking/ligand_options/Interface.hh>
#include <protocols/ligand_docking/MinimizeLigand.hh>
#include <protocols/ligand_docking/MoveMapBuilder.hh>
#include <protocols/ligand_docking/TetherLigand.hh>
#include <core/pose/util.hh>

#include <protocols/ligand_docking/UnconstrainedTorsionsMover.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MinMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/PackRotamersMover.hh>
#include <protocols/moves/RotamerTrialsMover.hh>
#include <protocols/moves/RigidBodyMover.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Conformation.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pack/task/ResfileReader.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/rotamer_set/UnboundRotamersOperation.hh>

#include <core/scoring/ScoreFunction.hh>
// Package Headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/constraints/ResidueTypeConstraint.hh>
// Auto-header: duplicate removed #include <core/kinematics/MoveMap.hh>

// Project Headers
#include <core/chemical/ResidueType.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

// Utility Headers
#include <utility/string_util.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <core/kinematics/Edge.hh>
#include <core/pack/task/PackerTask.hh>

// Scripter Headers
#include <utility/tag/Tag.hh>
#include <protocols/moves/DataMap.hh>

// Boost Headers
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

//STL headers
#include <string>

#include <set>

#include <platform/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/Adduct.fwd.hh>
#include <core/chemical/Adduct.hh>
#include <core/chemical/AtomICoor.fwd.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/AtomType.fwd.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.fwd.hh>
#include <core/chemical/ElementSet.fwd.hh>
#include <core/chemical/MMAtomType.fwd.hh>
#include <core/chemical/MMAtomTypeSet.fwd.hh>
#include <core/chemical/ResConnID.fwd.hh>
#include <core/chemical/ResConnID.hh>
#include <core/chemical/ResidueConnection.fwd.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/ResidueSelector.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
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
#include <core/conformation/PseudoBond.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
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
#include <core/graph/Graph.fwd.hh>
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
#include <core/id/SequenceMapping.fwd.hh>
#include <core/id/SequenceMapping.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/TorsionID.hh>
#include <core/id/TorsionID_Range.fwd.hh>
#include <core/id/TorsionID_Range.hh>
#include <core/id/types.hh>
#include <core/kinematics/AtomPointer.fwd.hh>
#include <core/kinematics/AtomTree.fwd.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/AtomWithDOFChange.fwd.hh>
#include <core/kinematics/DomainMap.fwd.hh>
#include <core/kinematics/Edge.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MinimizerMapBase.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/kinematics/RT.fwd.hh>
#include <core/kinematics/RT.hh>
#include <core/kinematics/ResidueCoordinateChangeList.fwd.hh>
#include <core/kinematics/ShortestPathInFoldTree.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/types.hh>
#include <core/kinematics/tree/Atom.fwd.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/optimization/MinimizerOptions.fwd.hh>
#include <core/pack/dunbrack/ChiSet.fwd.hh>
#include <core/pack/interaction_graph/InteractionGraphBase.fwd.hh>
#include <core/pack/interaction_graph/InteractionGraphBase.hh>
#include <core/pack/interaction_graph/OnTheFlyInteractionGraph.fwd.hh>
#include <core/pack/interaction_graph/PrecomputedPairEnergiesInteractionGraph.fwd.hh>
#include <core/pack/rotamer_set/FixbbRotamerSets.fwd.hh>
#include <core/pack/rotamer_set/FixbbRotamerSets.hh>
#include <core/pack/rotamer_set/RotamerCouplings.fwd.hh>
#include <core/pack/rotamer_set/RotamerSet.fwd.hh>
#include <core/pack/rotamer_set/RotamerSetOperation.fwd.hh>
#include <core/pack/rotamer_set/RotamerSetOperation.hh>
#include <core/pack/rotamer_set/RotamerSets.fwd.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSetsBase.fwd.hh>
#include <core/pack/rotamer_set/RotamerSetsBase.hh>
#include <core/pack/rotamer_set/UnboundRotamersOperation.fwd.hh>
#include <core/pack/task/IGEdgeReweightContainer.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/ResfileReader.fwd.hh>
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
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/EnergyGraph.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/LREnergyContainer.fwd.hh>
#include <core/scoring/MinimizationGraph.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionInfo.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/Func.fwd.hh>
#include <core/scoring/constraints/Func.hh>
#include <core/scoring/constraints/FuncFactory.fwd.hh>
#include <core/scoring/constraints/HarmonicFunc.fwd.hh>
#include <core/scoring/constraints/HarmonicFunc.hh>
#include <core/scoring/constraints/ResidueTypeConstraint.fwd.hh>
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
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/jobdist/Jobs.fwd.hh>
#include <protocols/jobdist/Jobs.hh>
#include <protocols/ligand_docking/HighResDocker.fwd.hh>
#include <protocols/ligand_docking/InterfaceBuilder.fwd.hh>
#include <protocols/ligand_docking/LigandArea.fwd.hh>
#include <protocols/ligand_docking/LigandArea.hh>
#include <protocols/ligand_docking/MinimizeBackbone.fwd.hh>
#include <protocols/ligand_docking/MinimizeBackbone.hh>
#include <protocols/ligand_docking/MinimizeLigand.fwd.hh>
#include <protocols/ligand_docking/MoveMapBuilder.fwd.hh>
#include <protocols/ligand_docking/ResidueTorsionRestraints.fwd.hh>
#include <protocols/ligand_docking/ResidueTorsionRestraints.hh>
#include <protocols/ligand_docking/TetherLigand.fwd.hh>
#include <protocols/ligand_docking/UnconstrainedTorsionsMover.fwd.hh>
#include <protocols/ligand_docking/ligand_options/Interface.fwd.hh>
#include <protocols/loops/Loop.fwd.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/moves/MetropolisHastingsMover.fwd.hh>
#include <protocols/moves/MinMover.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/moves/MonteCarloStatus.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/MoverCreator.hh>
#include <protocols/moves/MoverStatus.hh>
#include <protocols/moves/PackRotamersMover.fwd.hh>
#include <protocols/moves/RigidBodyMover.fwd.hh>
#include <protocols/moves/RotamerTrialsMover.fwd.hh>
#include <protocols/moves/ThermodynamicMover.fwd.hh>
#include <protocols/moves/ThermodynamicMover.hh>
#include <protocols/moves/mc_convergence_checks/ConvergenceCheck.fwd.hh>
#include <protocols/moves/mc_convergence_checks/ConvergenceCheck.hh>
#include <utility/Bound.fwd.hh>
#include <utility/Bound.hh>
#include <utility/PyAssert.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
#include <utility/stream_util.hh>
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
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iosfwd>
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <ostream>
#include <sstream>
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
#include <basic/options/keys/enzdes.OptionKeys.gen.hh>
#include <boost/algorithm/string/erase.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/functional/hash.hpp>

namespace protocols {
namespace ligand_docking {

static basic::Tracer high_res_docker_tracer("protocols.ligand_docking.ligand_options.Protocol", basic::t_debug);

std::string
HighResDockerCreator::keyname() const
{
	return HighResDockerCreator::mover_name();
}

protocols::moves::MoverOP
HighResDockerCreator::create_mover() const {
	return new HighResDocker;
}

std::string
HighResDockerCreator::mover_name()
{
	return "HighResDocker";
}

///@brief
HighResDocker::HighResDocker():
		Mover("HighResDocker"),
		num_cycles_(0),
		repack_every_Nth_(0),
		score_fxn_(NULL),
		movemap_builder_(NULL),
		resfile_("")
{
	resfile_.clear();
	// Now use cycles and repack_every_Nth to replicate these options...
	//meiler2006: 50, 8;
	//abbreviated: 5, 4;
	//abbrev2: 6, 3;
}

HighResDocker::HighResDocker(HighResDocker const & that):
	    //utility::pointer::ReferenceCount(),
		protocols::moves::Mover( that ),
		num_cycles_(that.num_cycles_),
		repack_every_Nth_(that.repack_every_Nth_),
		//chains_(that.chains_),
		score_fxn_(that.score_fxn_),
		movemap_builder_(that.movemap_builder_)
{}

HighResDocker::~HighResDocker() {}

protocols::moves::MoverOP HighResDocker::clone() const {
	return new HighResDocker( *this );
}

protocols::moves::MoverOP HighResDocker::fresh_instance() const {
	return new HighResDocker;
}

std::string HighResDocker::get_name() const{
	return "HighResDocker";
}

///@brief parse XML (specifically in the context of the parser/scripting scheme)
void
HighResDocker::parse_my_tag(
		utility::tag::TagPtr const tag,
		protocols::moves::DataMap & datamap,
		protocols::filters::Filters_map const & /*filters*/,
		protocols::moves::Movers_map const & /*movers*/,
		core::pose::Pose const & /*pose*/
)
{
	if ( tag->getName() != "HighResDocker" ){
		utility_exit_with_message("This should be impossible");
	}

	// cycles and repack_every_Nth
	if ( ! tag->hasOption("cycles") ) utility_exit_with_message("'HighResDocker' mover requires cycles tag");
	if ( ! tag->hasOption("repack_every_Nth") ) utility_exit_with_message("'HighResDocker' mover requires repack_every_Nth tag");
	num_cycles_= tag->getOption<core::Size>("cycles");
	repack_every_Nth_= tag->getOption<core::Size>("repack_every_Nth");

	/// Score Function ///
	if ( ! tag->hasOption("scorefxn") ) utility_exit_with_message("'HighResDocker' requires 'scorefxn' tag");
	std::string scorefxn_name= tag->getOption<std::string>("scorefxn");
	score_fxn_= datamap.get< core::scoring::ScoreFunction * >( "scorefxns", scorefxn_name);

	/// MoveMapBuilder///
	if ( ! tag->hasOption("movemap_builder") ) utility_exit_with_message("'HighResDocker' requires 'movemap_builder' tag");
	std::string movemap_builder_name= tag->getOption<std::string>("movemap_builder");
	movemap_builder_= datamap.get< MoveMapBuilder * >( "movemap_builders", movemap_builder_name);

	/// Resfile ///
	if ( tag->hasOption("resfile") ){
		resfile_= tag->getOption<std::string>("resfile");
	}
}

MinimizeLigandOPs
HighResDocker::setup_ligands_to_minimize(core::pose::Pose pose){
	MinimizeLigandOPs minimize_ligands;

	LigandAreas const ligand_areas =
			movemap_builder_->get_sc_interface_builder()->get_ligand_areas();
	LigandAreas::const_iterator ligand_area_itr= ligand_areas.begin();
	LigandAreas::const_iterator const ligand_area_end= ligand_areas.end();
	//TODO Use BOOST_FOREACH
	for(; ligand_area_itr != ligand_area_end; ++ligand_area_itr){
		char const & chain= ligand_area_itr->first;
		LigandAreaOP const ligand_area = ligand_area_itr->second;
		core::Real const & degrees = ligand_area->minimize_ligand_;
		if(degrees > 0){
			MinimizeLigandOP minimize_ligand = new MinimizeLigand(chain, degrees);
			minimize_ligand->apply(pose);
			minimize_ligands.push_back(minimize_ligand);
		}
	}
	return minimize_ligands;
}

TetherLigandOPs
HighResDocker::tether_ligands(core::pose::Pose & pose){
	TetherLigandOPs ligand_tethers;

	LigandAreas const ligand_areas =
			movemap_builder_->get_sc_interface_builder()->get_ligand_areas();
	LigandAreas::const_iterator ligand_area_itr= ligand_areas.begin();
	LigandAreas::const_iterator const ligand_area_end= ligand_areas.end();

	for(; ligand_area_itr != ligand_area_end; ++ligand_area_itr){
		char const & chain= ligand_area_itr->first;
		LigandAreaOP const ligand_area = ligand_area_itr->second;
		core::Real const & tether_size = ligand_area->tether_ligand_;
		if(tether_size > 0){
			TetherLigandOP tether_ligand= new TetherLigand(chain, tether_size);
			tether_ligand->apply(pose);
			ligand_tethers.push_back(tether_ligand);
		}
	}
	return ligand_tethers;
}

void
HighResDocker::remove_ligand_tethers(core::pose::Pose pose, TetherLigandOPs ligand_tethers){
	TetherLigandOPs::const_iterator begin= ligand_tethers.begin();
	TetherLigandOPs::const_iterator const end= ligand_tethers.end();

	foreach(TetherLigandOP ligand_tether, ligand_tethers){
		ligand_tether->release(pose);
	}
}

void
HighResDocker::apply(core::pose::Pose & pose) {
	assert(num_cycles_ > 0);

	MinimizeLigandOPs minimized_ligands = setup_ligands_to_minimize(pose);

	TetherLigandOPs ligand_tethers= tether_ligands(pose);

	assert(movemap_builder_ && score_fxn_ ); // make sure the pointers point
	core::kinematics::MoveMapOP movemap = movemap_builder_->build(pose);

	protocols::moves::MonteCarloOP monteCarlo = new protocols::moves::MonteCarlo(pose, *score_fxn_, 2.0);/* temperature, from RosettaLigand paper */
	score_fxn_->score( pose ); // without this neither of the movers below were working
	// I believe that this may have been related to adding constraints incorrectly at other places in my code.
	// Rigid body exploration
	utility::vector1<protocols::moves::MoverOP> rigid_body_movers= create_rigid_body_movers(pose);

	for( core::Size cycle = 1; cycle <= num_cycles_; ++cycle ) {
		core::pack::task::PackerTaskOP packer_task = make_packer_task(pose);// has to be in the loop to be updated after each design cycle (w/resfiles)

		protocols::moves::MoverOP pack_mover;

		if(cycle % repack_every_Nth_ == 1){
			high_res_docker_tracer.Debug << "making PackRotamersMover" << std::endl;
			pack_mover= (protocols::moves::Mover *) new protocols::moves::PackRotamersMover(score_fxn_, packer_task);
		}
		else{
			high_res_docker_tracer.Debug << "making RotamerTrialsMover" << std::endl;
			pack_mover= (protocols::moves::Mover *) new protocols::moves::RotamerTrialsMover(score_fxn_, *packer_task);
		}

		// Wrap it in something to disable the torsion constraints before packing!
		pack_mover = new protocols::ligand_docking::UnconstrainedTorsionsMover( pack_mover, minimized_ligands );

		protocols::moves::MinMoverOP min_mover = new protocols::moves::MinMover( movemap, score_fxn_, "dfpmin_armijo_nonmonotone_atol", 1.0, true /*use_nblist*/ );
		min_mover->min_options()->nblist_auto_update(true); // does this cost us lots of time in practice?

		core::Real const score1 = (*score_fxn_)( pose );
		apply_rigid_body_moves(pose, rigid_body_movers);
		pack_mover->apply(pose);

		core::Real const score2 = (*score_fxn_)( pose );
		if(score2 - score1 < 15.0) {
			min_mover->apply(pose);
		}

		monteCarlo->boltzmann( pose );

	}

	remove_ligand_tethers(pose, ligand_tethers);
	// keep the best structure we found, not the current one
	monteCarlo->show_scores();
	monteCarlo->recover_low(pose);
}

core::pack::task::PackerTaskOP
HighResDocker::make_packer_task_from_vector(
		core::pose::Pose const & pose,
		ligand_options::Interface const allow_repack
) const{
	static bool pose_already_packed= false;
	core::pack::task::PackerTaskOP pack_task = core::pack::task::TaskFactory::create_packer_task(pose);
	pack_task->initialize_from_command_line(); // -ex1 -ex2  etc.

	core::pack::rotamer_set::UnboundRotamersOperationOP unboundrot_ = new core::pack::rotamer_set::UnboundRotamersOperation();
	unboundrot_->initialize_from_command_line();
	pack_task->append_rotamerset_operation( unboundrot_ );

	for(core::Size i = 1; i <= pose.total_residue(); ++i) {
		/// If several params files have the same name, allow switching among them
		/// This was previously only enabled with mutate_same_name3.  Now default.
		if( ! pose.residue(i).is_ligand()) continue;
		high_res_docker_tracer.Debug<<  "enabling packing for ligand residue "<< i << std::endl;
		enable_ligand_rotamer_packing(pose, i, pack_task);
	}

	if(resfile_.empty()){
		bool const use_resfile= basic::options::option[ basic::options::OptionKeys::packing::resfile ].user() ;
		if (use_resfile){
			high_res_docker_tracer<< "using OPTIONS resfile"<< std::endl;
			core::pack::task::parse_resfile(pose, *pack_task);
		}
		else{
			high_res_docker_tracer<< "restricting to repack"<< std::endl;
			for(core::Size i = 1; i <= pose.total_residue(); ++i) {
				if( ! pose.residue(i).is_ligand() )
				{
					pack_task->nonconst_residue_task( i ).restrict_to_repacking();
				}
			}
		}
	}
	else{
		high_res_docker_tracer<< "using XML resfile"<< std::endl;
		core::pack::task::parse_resfile(pose, *pack_task);
	}


	for(core::Size i = 1; i <= pose.total_residue(); ++i) {
		if ( allow_repack[i].type == ligand_options::InterfaceInfo::non_interface  ){
			pack_task->nonconst_residue_task( i ).prevent_repacking();
		}
	}

	// We always want the option (after the initial unbiased pack)
	// of sticking with our current nicely minimized conformation.
	if( pose_already_packed ){
		pack_task->or_include_current(true);
	}else {
		pose_already_packed=true;
	}

	return pack_task;
}

core::pack::task::PackerTaskOP
HighResDocker::make_packer_task(
		core::pose::Pose const & pose,
		bool all_residues
) const{
	if ( all_residues ){
		ligand_options::Interface interface(pose.n_residue(), ligand_options::InterfaceInfo(ligand_options::InterfaceInfo::is_interface)); // 0 is false, #
		return make_packer_task_from_vector(pose, interface);
	}else{ // the packer task interface should match the movemap interface
		InterfaceBuilderOP sc_interface_builder= movemap_builder_->get_sc_interface_builder();
		ligand_options::Interface side_chain_interface= sc_interface_builder->build(pose);
		return make_packer_task_from_vector(pose, side_chain_interface);
	}
}

void
HighResDocker::enable_ligand_rotamer_packing(
		core::pose::Pose const & pose,
		core::Size const ligand_residue_id,
		core::pack::task::PackerTaskOP & pack_task
) const{
	core::conformation::Residue const & this_residue= pose.residue(ligand_residue_id);
	core::chemical::ResidueTypeSet const & rsd_type_set = this_residue.residue_type_set();
	core::chemical::ResidueTypeCAPs allowed_types = rsd_type_set.name3_map( this_residue.name3() ); // a vector1

	assert(allowed_types.size() > 0);
	/// TODO consider removing this so resfiles can specify ligand mutations to allow
	if( allowed_types.size() == 1){
		pack_task->nonconst_residue_task( ligand_residue_id ).restrict_to_repacking();
		return;
	}
	// else
	for( core::Size j = 1; j <= allowed_types.size(); ++j ) {
		if( allowed_types[j]->name() == this_residue.name() ) continue; // already in the task's list
		///TODO figure out why this is nonconst.  Perhaps it could be const
		pack_task->nonconst_residue_task( ligand_residue_id ).allow_noncanonical_aa( allowed_types[j]->name() );
	}
}

utility::vector1<protocols::moves::MoverOP>
HighResDocker::create_rigid_body_movers(core::pose::Pose const & pose) const{
	utility::vector1<protocols::moves::MoverOP> rigid_body_movers;

	LigandAreas const ligand_areas =
			movemap_builder_->get_sc_interface_builder()->get_ligand_areas();

	foreach(LigandAreas::value_type ligand_area_pair, ligand_areas){
		char const & chain= ligand_area_pair.first;
		core::Size jump_id= core::pose::get_jump_id_from_chain(chain, pose);

		LigandAreaOP const ligand_area = ligand_area_pair.second;
		core::Real const & angstroms= ligand_area->high_res_angstroms_;
		core::Real const & degrees= ligand_area->high_res_degrees_;
		protocols::moves::MoverOP rigid_body_mover= new protocols::moves::RigidBodyPerturbMover( jump_id, degrees, angstroms);
		rigid_body_movers.push_back(rigid_body_mover);
	}
	return rigid_body_movers;
}

void HighResDocker::apply_rigid_body_moves(
		core::pose::Pose & pose,
		utility::vector1<protocols::moves::MoverOP> & rigid_body_movers
){
	utility::vector1<protocols::moves::MoverOP>::iterator rigid_body_mover= rigid_body_movers.begin();
	foreach(protocols::moves::MoverOP rigid_body_mover, rigid_body_movers){
		rigid_body_mover->apply(pose);
	}
}

/// Non-member functions

// Favor Native is part of the APPLY_TO_POSE section

} //namespace ligand_docking
} //namespace protocols
