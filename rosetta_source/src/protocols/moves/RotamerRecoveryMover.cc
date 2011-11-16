// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/moves/RotamerRecoveryMover.cc
/// @brief A wrapper that measures how similar the rotamers are between before and after running the child mover
/// @author Matthew O'Meara (mattjomeara@gmail.com)
/// Adapted from:
/// protocols::optimize_weights::IterativeOptEDriver::measure_rotamer_recovery()
/// and apps::pilot::doug::rotamer_prediction_benchmark()

#include <protocols/moves/RotamerRecoveryMover.hh>
#include <string>


// Setup Mover
#include <protocols/moves/RotamerRecoveryMoverCreator.hh>
namespace protocols{
namespace moves{

std::string
RotamerRecoveryMoverCreator::keyname() const
{
	return RotamerRecoveryMoverCreator::mover_name();
}

MoverOP
RotamerRecoveryMoverCreator::create_mover() const {
	return new RotamerRecoveryMover;
}

std::string
RotamerRecoveryMoverCreator::mover_name()
{
	return "RotamerRecoveryMover";
}

}
}


// Unit Headers
#include <protocols/moves/Mover.hh>
#include <protocols/moves/DataMap.hh>
#include <protocols/rotamer_recovery/RotamerRecovery.hh>

// Project Headers
#include <basic/Tracer.hh>
#include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/graph/Graph.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
// AUTO-REMOVED #include <core/pack/task/operation/TaskOperationFactory.hh>
// AUTO-REMOVED #include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
// AUTO-REMOVED #include <core/pose/util.hh>
// AUTO-REMOVED #include <core/pack/dunbrack/DunbrackRotamer.hh>
// AUTO-REMOVED #include <core/pack/dunbrack/RotamerLibrary.hh>
// AUTO-REMOVED #include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
// AUTO-REMOVED #include <core/scoring/TenANeighborGraph.hh>
#include <core/types.hh>
// AUTO-REMOVED #include <protocols/moves/PackRotamersMover.hh>
#include <protocols/rotamer_recovery/RotamerRecoveryFactory.hh>
#include <protocols/rotamer_recovery/RRProtocol.hh>
#include <protocols/rotamer_recovery/RRProtocolMover.hh>
#include <protocols/rotamer_recovery/RRComparer.hh>
#include <protocols/rotamer_recovery/RRReporter.hh>

// Option System Headers
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/exit.hh>
#include <utility/file/FileName.hh>
#include <utility/tag/Tag.hh>

// Numeric Headers
// AUTO-REMOVED #include <numeric/angle.functions.hh>

// C++ Headers
#include <algorithm>
#include <iostream>
#include <fstream>

#include <platform/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/Adduct.fwd.hh>
#include <core/chemical/Adduct.hh>
#include <core/chemical/AtomICoor.fwd.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/AtomType.fwd.hh>
// AUTO-REMOVED #include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.fwd.hh>
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
// AUTO-REMOVED #include <core/chemical/types.hh>
#include <core/chemical/orbitals/ICoorOrbitalData.hh>
#include <core/chemical/orbitals/OrbitalType.fwd.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.fwd.hh>
#include <core/chemical/sdf/MolData.fwd.hh>
#include <core/chemical/sdf/MolData.hh>
#include <core/conformation/Atom.fwd.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/Conformation.fwd.hh>
// AUTO-REMOVED #include <core/conformation/Conformation.hh>
// AUTO-REMOVED #include <core/conformation/PointGraph.fwd.hh>
// AUTO-REMOVED #include <core/conformation/PointGraphData.fwd.hh>
#include <core/conformation/PseudoBond.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/RotamerSetBase.fwd.hh>
#include <core/conformation/orbitals/OrbitalXYZCoords.hh>
#include <core/conformation/signals/ConnectionEvent.fwd.hh>
// AUTO-REMOVED #include <core/conformation/signals/ConnectionEvent.hh>
// AUTO-REMOVED #include <core/conformation/signals/GeneralEvent.fwd.hh>
// AUTO-REMOVED #include <core/conformation/signals/GeneralEvent.hh>
#include <core/conformation/signals/IdentityEvent.fwd.hh>
// AUTO-REMOVED #include <core/conformation/signals/IdentityEvent.hh>
#include <core/conformation/signals/LengthEvent.fwd.hh>
// AUTO-REMOVED #include <core/conformation/signals/LengthEvent.hh>
// AUTO-REMOVED #include <core/conformation/signals/XYZEvent.fwd.hh>
// AUTO-REMOVED #include <core/conformation/signals/XYZEvent.hh>
#include <core/conformation/symmetry/SymmetricConformation.fwd.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>
// AUTO-REMOVED #include <core/graph/ArrayPool.hh>
// AUTO-REMOVED #include <core/graph/Graph.fwd.hh>
// AUTO-REMOVED #include <core/graph/UpperEdgeGraph.fwd.hh>
// AUTO-REMOVED #include <core/graph/unordered_object_pool.fwd.hpp>
#include <core/id/AtomID.fwd.hh>
// AUTO-REMOVED #include <core/id/AtomID.hh>
// AUTO-REMOVED #include <core/id/AtomID_Map.fwd.hh>
// AUTO-REMOVED #include <core/id/AtomID_Map.hh>
// AUTO-REMOVED #include <core/id/AtomID_Mask.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
// AUTO-REMOVED #include <core/id/DOF_ID.hh>
// AUTO-REMOVED #include <core/id/DOF_ID_Map.fwd.hh>
// AUTO-REMOVED #include <core/id/DOF_ID_Map.hh>
// AUTO-REMOVED #include <core/id/DOF_ID_Mask.fwd.hh>
// AUTO-REMOVED #include <core/id/NamedAtomID.fwd.hh>
// AUTO-REMOVED #include <core/id/NamedAtomID.hh>
// AUTO-REMOVED #include <core/id/NamedStubID.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
// AUTO-REMOVED #include <core/id/types.hh>
// AUTO-REMOVED #include <core/kinematics/AtomPointer.fwd.hh>
// AUTO-REMOVED #include <core/kinematics/AtomTree.fwd.hh>
// AUTO-REMOVED #include <core/kinematics/AtomTree.hh>
// AUTO-REMOVED #include <core/kinematics/AtomWithDOFChange.fwd.hh>
#include <core/kinematics/DomainMap.fwd.hh>
// AUTO-REMOVED #include <core/kinematics/Edge.fwd.hh>
// AUTO-REMOVED #include <core/kinematics/Edge.hh>
// AUTO-REMOVED #include <core/kinematics/FoldTree.fwd.hh>
// AUTO-REMOVED #include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MinimizerMapBase.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/kinematics/RT.fwd.hh>
#include <core/kinematics/RT.hh>
// AUTO-REMOVED #include <core/kinematics/ResidueCoordinateChangeList.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/kinematics/Stub.hh>
// AUTO-REMOVED #include <core/kinematics/tree/Atom.fwd.hh>
// AUTO-REMOVED #include <core/kinematics/tree/Atom.hh>
// AUTO-REMOVED #include <core/pack/dunbrack/ChiSet.fwd.hh>
// AUTO-REMOVED #include <core/pack/dunbrack/DunbrackRotamer.fwd.hh>
// AUTO-REMOVED #include <core/pack/dunbrack/RotamerLibrary.fwd.hh>
// AUTO-REMOVED #include <core/pack/dunbrack/RotamerLibraryScratchSpace.fwd.hh>
// AUTO-REMOVED #include <core/pack/dunbrack/SemiRotamericSingleResidueDunbrackLibrary.fwd.hh>
#include <core/pack/interaction_graph/InteractionGraphBase.fwd.hh>
// AUTO-REMOVED #include <core/pack/interaction_graph/InteractionGraphBase.hh>
// AUTO-REMOVED #include <core/pack/interaction_graph/OnTheFlyInteractionGraph.fwd.hh>
// AUTO-REMOVED #include <core/pack/interaction_graph/PrecomputedPairEnergiesInteractionGraph.fwd.hh>
#include <core/pack/rotamer_set/FixbbRotamerSets.fwd.hh>
// AUTO-REMOVED #include <core/pack/rotamer_set/FixbbRotamerSets.hh>
#include <core/pack/rotamer_set/RotamerCouplings.fwd.hh>
// AUTO-REMOVED #include <core/pack/rotamer_set/RotamerSet.fwd.hh>
#include <core/pack/rotamer_set/RotamerSetOperation.fwd.hh>
#include <core/pack/rotamer_set/RotamerSets.fwd.hh>
// AUTO-REMOVED #include <core/pack/rotamer_set/RotamerSets.hh>
// AUTO-REMOVED #include <core/pack/rotamer_set/RotamerSetsBase.fwd.hh>
// AUTO-REMOVED #include <core/pack/rotamer_set/RotamerSetsBase.hh>
#include <core/pack/rotamer_set/symmetry/SymmetricRotamerSets.fwd.hh>
#include <core/pack/task/IGEdgeReweightContainer.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/ResfileReader.fwd.hh>
#include <core/pack/task/RotamerSampleOptions.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
// AUTO-REMOVED #include <core/pack/task/operation/ResFilterCreator.fwd.hh>
// AUTO-REMOVED #include <core/pack/task/operation/ResFilterFactory.fwd.hh>
// AUTO-REMOVED #include <core/pack/task/operation/ResLvlTaskOperationCreator.fwd.hh>
// AUTO-REMOVED #include <core/pack/task/operation/ResLvlTaskOperationFactory.fwd.hh>
#include <core/pack/task/operation/TaskOperation.fwd.hh>
#include <core/pack/task/operation/TaskOperation.hh>
// AUTO-REMOVED #include <core/pack/task/operation/TaskOperationCreator.fwd.hh>
// AUTO-REMOVED #include <core/pack/task/operation/TaskOperationFactory.fwd.hh>
#include <core/pack/task/operation/TaskOperations.fwd.hh>
// AUTO-REMOVED #include <core/pose/MiniPose.fwd.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/PDBPoseMap.fwd.hh>
#include <core/pose/PDBPoseMap.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/Remarks.fwd.hh>
#include <core/pose/Remarks.hh>
// AUTO-REMOVED #include <core/pose/datacache/ObserverCache.fwd.hh>
// AUTO-REMOVED #include <core/pose/metrics/PoseMetricContainer.fwd.hh>
// AUTO-REMOVED #include <core/pose/signals/ConformationEvent.fwd.hh>
// AUTO-REMOVED #include <core/pose/signals/DestructionEvent.fwd.hh>
// AUTO-REMOVED #include <core/pose/signals/EnergyEvent.fwd.hh>
// AUTO-REMOVED #include <core/pose/signals/GeneralEvent.fwd.hh>
// AUTO-REMOVED #include <core/scoring/ContextGraph.fwd.hh>
// AUTO-REMOVED #include <core/scoring/ContextGraph.hh>
// AUTO-REMOVED #include <core/scoring/ContextGraphTypes.hh>
// AUTO-REMOVED #include <core/scoring/Energies.fwd.hh>
// AUTO-REMOVED #include <core/scoring/EnergiesCacheableDataType.hh>
#include <core/scoring/EnergyGraph.fwd.hh>
// AUTO-REMOVED #include <core/scoring/EnergyGraph.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/LREnergyContainer.fwd.hh>
#include <core/scoring/MinimizationGraph.fwd.hh>
// AUTO-REMOVED #include <core/scoring/NeighborList.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionInfo.fwd.hh>
#include <core/scoring/ScoreType.hh>
// AUTO-REMOVED #include <core/scoring/TenANeighborGraph.fwd.hh>
// AUTO-REMOVED #include <core/scoring/TwelveANeighborGraph.fwd.hh>
// AUTO-REMOVED #include <core/scoring/constraints/Constraint.fwd.hh>
// AUTO-REMOVED #include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/methods/ContextDependentLRTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextDependentOneBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextDependentTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentLRTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/EnergyMethod.fwd.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/methods/LongRangeTwoBodyEnergy.fwd.hh>
// AUTO-REMOVED #include <core/scoring/methods/Methods.hh>
#include <core/scoring/methods/TwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/WholeStructureEnergy.fwd.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/jobdist/Jobs.fwd.hh>
// AUTO-REMOVED #include <protocols/jobdist/Jobs.hh>
#include <protocols/moves/DataMap.fwd.hh>
// Auto-header: duplicate removed #include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverCreator.hh>
#include <protocols/moves/MoverStatus.hh>
// AUTO-REMOVED #include <protocols/moves/PackRotamersMover.fwd.hh>
#include <protocols/moves/RotamerRecoveryMover.fwd.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/rotamer_recovery/RRComparer.fwd.hh>
#include <protocols/rotamer_recovery/RRReporter.fwd.hh>
#include <protocols/rotamer_recovery/RotamerRecovery.fwd.hh>
#include <utility/Bound.fwd.hh>
#include <utility/Bound.hh>
#include <utility/PyAssert.hh>
#include <utility/down_cast.hh>
// AUTO-REMOVED #include <utility/fixedsizearray1.fwd.hh>
// AUTO-REMOVED #include <utility/fixedsizearray1.hh>
// AUTO-REMOVED #include <utility/stream_util.hh>
// AUTO-REMOVED #include <utility/string_util.hh>
#include <utility/vector0.fwd.hh>
#include <utility/vector0.hh>
#include <utility/vector0_bool.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/file/FileName.fwd.hh>
#include <utility/file/PathName.fwd.hh>
#include <utility/file/PathName.hh>
// AUTO-REMOVED #include <utility/io/izstream.fwd.hh>
// AUTO-REMOVED #include <utility/io/ozstream.fwd.hh>
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
// AUTO-REMOVED #include <utility/signals/BufferedSignalHub.fwd.hh>
// AUTO-REMOVED #include <utility/signals/BufferedSignalHub.hh>
// AUTO-REMOVED #include <utility/signals/Link.fwd.hh>
// AUTO-REMOVED #include <utility/signals/Link.hh>
// AUTO-REMOVED #include <utility/signals/LinkUnit.fwd.hh>
// AUTO-REMOVED #include <utility/signals/LinkUnit.hh>
// AUTO-REMOVED #include <utility/signals/PausableSignalHub.fwd.hh>
// AUTO-REMOVED #include <utility/signals/PausableSignalHub.hh>
// AUTO-REMOVED #include <utility/signals/SignalHub.fwd.hh>
// AUTO-REMOVED #include <utility/signals/SignalHub.hh>
#include <utility/tag/Tag.fwd.hh>
// AUTO-REMOVED #include <numeric/NumericTraits.hh>
// AUTO-REMOVED #include <numeric/constants.hh>
// AUTO-REMOVED #include <numeric/conversions.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/sphericalVector.fwd.hh>
// AUTO-REMOVED #include <numeric/sphericalVector.hh>
#include <numeric/trig.functions.hh>
// AUTO-REMOVED #include <numeric/types.hh>
#include <numeric/xyz.functions.fwd.hh>
// AUTO-REMOVED #include <numeric/xyz.functions.hh>
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
// AUTO-REMOVED #include <numeric/random/random.fwd.hh>
// AUTO-REMOVED #include <ObjexxFCL/Dimension.fwd.hh>
// AUTO-REMOVED #include <ObjexxFCL/Dimension.hh>
// AUTO-REMOVED #include <ObjexxFCL/DimensionExpression.hh>
// AUTO-REMOVED #include <ObjexxFCL/DynamicIndexRange.fwd.hh>
// AUTO-REMOVED #include <ObjexxFCL/DynamicIndexRange.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArray.fwd.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArray.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArray1.fwd.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArray1.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArray1D.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArray2.fwd.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArray2.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArray2D.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArrayInitializer.fwd.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArrayInitializer.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArraySection.fwd.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArraySection.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArrayTraits.fwd.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArrayTraits.hh>
// AUTO-REMOVED #include <ObjexxFCL/IndexRange.fwd.hh>
// AUTO-REMOVED #include <ObjexxFCL/IndexRange.hh>
// AUTO-REMOVED #include <ObjexxFCL/InitializerSentinel.hh>
// AUTO-REMOVED #include <ObjexxFCL/Observer.fwd.hh>
// AUTO-REMOVED #include <ObjexxFCL/Observer.hh>
// AUTO-REMOVED #include <ObjexxFCL/ObserverMulti.hh>
// AUTO-REMOVED #include <ObjexxFCL/ObserverSingle.hh>
// AUTO-REMOVED #include <ObjexxFCL/ProxySentinel.hh>
// AUTO-REMOVED #include <ObjexxFCL/SetWrapper.fwd.hh>
// AUTO-REMOVED #include <ObjexxFCL/Star.fwd.hh>
// AUTO-REMOVED #include <ObjexxFCL/Star.hh>
#include <ObjexxFCL/TypeTraits.hh>
#include <ObjexxFCL/char.functions.hh>
// AUTO-REMOVED #include <ObjexxFCL/proxy_const_assert.hh>
#include <ObjexxFCL/string.functions.hh>
#include <cassert>
#include <cmath>
#include <cstddef>
// AUTO-REMOVED #include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iosfwd>
// AUTO-REMOVED #include <iterator>
#include <limits>
#include <list>
#include <map>
#include <ostream>
#include <set>
#include <sstream>
// AUTO-REMOVED #include <typeinfo>
#include <utility>
#include <vector>
// AUTO-REMOVED #include <basic/MetricValue.fwd.hh>
#include <basic/Tracer.fwd.hh>
// AUTO-REMOVED #include <basic/datacache/BasicDataCache.fwd.hh>
// AUTO-REMOVED #include <basic/datacache/BasicDataCache.hh>
// AUTO-REMOVED #include <basic/datacache/CacheableData.fwd.hh>
// AUTO-REMOVED #include <basic/datacache/CacheableData.hh>
// AUTO-REMOVED #include <basic/datacache/DataCache.fwd.hh>
// AUTO-REMOVED #include <basic/datacache/DataCache.hh>
#include <basic/options/keys/OptionKeys.hh>
// AUTO-REMOVED #include <boost/algorithm/string/erase.hpp>
// AUTO-REMOVED #include <boost/bind.hpp>
// AUTO-REMOVED #include <boost/config.hpp>
// AUTO-REMOVED #include <boost/function.hpp>
// AUTO-REMOVED #include <boost/functional/hash.hpp>
// AUTO-REMOVED #include <boost/pool/detail/mutex.hpp>
// AUTO-REMOVED #include <boost/pool/poolfwd.hpp>

//Auto Headers


//using std::ios::app;
using std::endl;
using std::max;
using std::ofstream;
using std::ostream;
using std::string;
using core::Real;
using core::Size;
using core::pack::pack_rotamers;
using core::pose::PoseOP;
using core::pack::task::PackerTaskOP;
using core::pack::task::TaskFactory;
using core::pack::task::TaskFactoryOP;
using core::pack::task::operation::InitializeFromCommandline;
using core::pack::task::operation::RestrictToRepacking;
using core::scoring::getScoreFunction;
using core::scoring::ScoreFunction;
using core::scoring::ScoreFunctionOP;
using basic::Tracer;
using utility::vector1;
using protocols::moves::MoverOP;
using protocols::rotamer_recovery::RRProtocolOP;
using protocols::rotamer_recovery::RRProtocolMover;
using protocols::rotamer_recovery::RRComparerOP;
using protocols::rotamer_recovery::RRReporterOP;
using protocols::rotamer_recovery::RotamerRecovery;
using protocols::rotamer_recovery::RotamerRecoveryFactory;
using protocols::rosetta_scripts::parse_mover;

namespace protocols {
namespace moves {

static Tracer TR("protocol.moves.RotamerRecoveryMover");

RotamerRecoveryMover::RotamerRecoveryMover() :
	rotamer_recovery_( NULL ),
	scfxn_(NULL),
	task_factory_(new TaskFactory),
	output_fname_("rotamer_recovery.out")
{
	task_factory_->push_back( new InitializeFromCommandline );
	task_factory_->push_back( new RestrictToRepacking );
}

RotamerRecoveryMover::RotamerRecoveryMover(
	string const & protocol,
	string const & comparer,
	string const & reporter,
	string const & output_fname,
	ScoreFunctionOP scfxn,
	TaskFactoryOP task_factory) :
	rotamer_recovery_(
		RotamerRecoveryFactory::get_instance()->get_rotamer_recovery(
			protocol, comparer, reporter )),
	scfxn_( scfxn ),
	task_factory_( task_factory ),
	output_fname_( output_fname )
	{ }

RotamerRecoveryMover::~RotamerRecoveryMover(){}

RotamerRecoveryMover::RotamerRecoveryMover( RotamerRecoveryMover const & src):
	//utility::pointer::ReferenceCount(),
	Mover( src ),
	rotamer_recovery_( src.rotamer_recovery_ ),
	scfxn_( src.scfxn_ ),
	task_factory_( src.task_factory_ )
{}


void
RotamerRecoveryMover::register_options() const {
	using basic::options::option;
	using namespace basic::options::OptionKeys;

	// Use RotamerRecovery to test new score functions eg. all the corrections
	option.add_relevant( corrections::correct );

	// Use full atom binary silent files for best io-performance
	option.add_relevant( in::file::fullatom );
	option.add_relevant( in::file::silent_struct_type );
	option.add_relevant( in::file::silent );

	// If using an outputter that writes to a database improve
	// io-performace by not writing out structures
	option.add_relevant( out::nooutput );

	rotamer_recovery_->register_options();

}

bool
RotamerRecoveryMover::reinitialize_for_each_job() const {
	return false;
}

bool
RotamerRecoveryMover::reinitialize_for_new_input() const {
	return false;
}

void
RotamerRecoveryMover::apply( Pose & pose
) {
	assert( rotamer_recovery_ );
	ScoreFunctionOP scfxn( score_function());
	scfxn->setup_for_scoring(pose);
	PackerTaskOP packer_task( task_factory_->create_task_and_apply_taskoperations( pose ));
	rotamer_recovery_->run(pose, *scfxn, *packer_task);
}

string
RotamerRecoveryMover::get_name() const {
	return "RotamerRecoveryMover";
}

MoverOP
RotamerRecoveryMover::fresh_instance() const {
	return new RotamerRecoveryMover;
}


MoverOP
RotamerRecoveryMover::clone() const {
	return new RotamerRecoveryMover( *this );
}

void
RotamerRecoveryMover::parse_my_tag(
	utility::tag::TagPtr const tag,
	DataMap & datamap,
	Filters_map const & /*filters*/,
	Movers_map const & movers,
	Pose const & /*pose*/ )
{
	string const & scorefxn_key( tag->getOption<std::string>("scorefxn", "score12" ));
	score_function( datamap.get< ScoreFunction * >("scorefxns", scorefxn_key) );

	if( rotamer_recovery_ ){
		TR << "WARNING: Attempting to redefine rotamer_recovery_ object from Parser Script" << endl;
		utility_exit();
	}

	if(tag->hasOption("protocol") && (tag->hasOption("mover") || tag->hasOption("mover_name"))){
		utility_exit_with_message("Please either the 'protocol' field or the 'mover' field but not both.");
	}

	RotamerRecoveryFactory * factory(RotamerRecoveryFactory::get_instance());

	RRProtocolOP protocol;
	if(tag->hasOption("mover") || tag->hasOption("mover_name")){
		MoverOP mover = parse_mover(tag->hasOption("mover") ?
			tag->getOption<string>("mover") : tag->getOption<string>("mover_name"), movers);
		protocol = new RRProtocolMover(mover);
	} else {
		protocol = factory->get_rotamer_recovery_protocol(tag->getOption<string>("protocol", "RRProtocolMinPack"));
	}
	RRComparerOP comparer(
		factory->get_rotamer_recovery_comparer(
			tag->getOption<string>("comparer", "RRComparerAutomorphicRMSD")));

	RRReporterOP reporter(
		factory->get_rotamer_recovery_reporter(
			tag->getOption<string>("reporter", "RRReporterSimple")));

	reporter->output_fname(tag->getOption<string>("output_fname"));

	rotamer_recovery_ = new RotamerRecovery(protocol, comparer, reporter);
}

ScoreFunctionOP
RotamerRecoveryMover::score_function(){
	if ( !scfxn_ )
		scfxn_ = getScoreFunction();

	return scfxn_;
}

void
RotamerRecoveryMover::score_function(
	ScoreFunctionOP scorefunction
) {
	scfxn_ = scorefunction;
}

void
RotamerRecoveryMover::show() {
	rotamer_recovery_->show();
}

void
RotamerRecoveryMover::show(
	ostream & out
) {
	rotamer_recovery_->show( out );
}


void
RotamerRecoveryMover::write_to_file(){
	write_to_file( output_fname_ );
}

void
RotamerRecoveryMover::write_to_file(
	string const & output_fname
) {
	ofstream fout;
	fout.open( output_fname.c_str(), std::ios::out );
	if( !fout.is_open() ){
		TR << "Unable to open output file '" << output_fname << "'." << endl;
		utility_exit();
	}

	rotamer_recovery_->show( fout );

	fout.close();

}



} // namespace moves
} // namespace protocols
