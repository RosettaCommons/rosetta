// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/etable/EtableEnergy.cc
/// @brief  Atom pair energy functions
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author Kevin P. Hinshaw (KevinHinshaw@gmail.com)
/// @author Andrew Leaver-Fay

#ifndef INCLUDED_core_scoring_etable_BaseEtableEnergy_tmpl_hh
#define INCLUDED_core_scoring_etable_BaseEtableEnergy_tmpl_hh

/// #define APL_TEMP_DEBUG


// Unit headers

// Package headers
#include <core/scoring/DerivVectorPair.hh>
#include <core/scoring/EnergiesCacheableDataType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/NeighborList.hh>
#include <core/scoring/NeighborList.tmpl.hh>
#include <core/scoring/ResidueNeighborList.hh>
#include <core/scoring/MinimizationData.hh>
#include <core/scoring/methods/Methods.hh>

#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>

//#include <core/scoring/etable/count_pair/CountPair1BC3.hh> // remove this
//#include <core/scoring/etable/count_pair/CountPair1BC4.hh> // remove this
//#include <core/scoring/etable/count_pair/CountPairIntraResC3.hh>
//#include <core/scoring/etable/count_pair/CountPairIntraResC4.hh>
#include <core/scoring/etable/count_pair/CountPairAll.hh>
#include <core/scoring/etable/count_pair/CountPairNone.hh>

#include <core/scoring/trie/CPDataCorrespondence.hh>
#include <core/scoring/trie/RotamerDescriptor.hh>
#include <core/scoring/trie/RotamerTrie.hh>
#include <core/scoring/trie/TrieCollection.hh>
#include <core/scoring/trie/TrieCountPairBase.hh>
#include <core/scoring/trie/trie.functions.hh>

#include <core/scoring/etable/etrie/CountPairData_1_1.hh>
#include <core/scoring/etable/etrie/CountPairData_1_2.hh>
#include <core/scoring/etable/etrie/CountPairData_1_3.hh>
#include <core/scoring/etable/etrie/CountPairDataGeneric.hh>
#include <core/scoring/etable/etrie/EtableAtom.hh>
#include <core/scoring/etable/etrie/TrieCountPair1BC3.hh>
#include <core/scoring/etable/etrie/TrieCountPair1BC4.hh>
#include <core/scoring/etable/etrie/TrieCountPairAll.hh>
#include <core/scoring/etable/etrie/TrieCountPairNone.hh>
#include <core/scoring/etable/etrie/TrieCountPairGeneric.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/kinematics/MinimizerMapBase.hh>
#include <core/conformation/RotamerSetBase.hh>

#include <core/chemical/AtomType.hh>
#include <basic/options/option.hh>

// option key includes

#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>

//Auto Headers
#include <platform/types.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/Adduct.fwd.hh>
#include <core/chemical/Adduct.hh>
#include <core/chemical/AtomICoor.fwd.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/AtomType.fwd.hh>
#include <core/chemical/AtomTypeSet.fwd.hh>
#include <core/chemical/AtomTypeSet.hh>
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
#include <core/chemical/types.hh>
#include <core/chemical/orbitals/ICoorOrbitalData.hh>
#include <core/chemical/orbitals/OrbitalType.fwd.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.fwd.hh>
#include <core/conformation/AbstractRotamerTrie.fwd.hh>
#include <core/conformation/AbstractRotamerTrie.hh>
#include <core/conformation/Atom.fwd.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/PointGraph.fwd.hh>
#include <core/conformation/PointGraph.hh>
#include <core/conformation/PointGraphData.fwd.hh>
#include <core/conformation/PointGraphData.hh>
#include <core/conformation/PseudoBond.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/RotamerSetBase.fwd.hh>
#include <core/conformation/find_neighbors.fwd.hh>
#include <core/conformation/find_neighbors.hh>
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
#include <utility/graph/ArrayPool.hh>
#include <utility/graph/Graph.fwd.hh>
#include <utility/graph/Graph.hh>
#include <utility/graph/UpperEdgeGraph.fwd.hh>
#include <utility/graph/UpperEdgeGraph.hh>
#include <utility/graph/unordered_object_pool.fwd.hpp>
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
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>
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
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MinimizerMapBase.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/kinematics/RT.fwd.hh>
#include <core/kinematics/RT.hh>
#include <core/kinematics/ResidueCoordinateChangeList.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/tree/Atom.fwd.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/pose/MiniPose.fwd.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/util.hh>
#include <core/pose/datacache/ObserverCache.fwd.hh>
#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
#include <core/pose/signals/ConformationEvent.fwd.hh>
#include <core/pose/signals/DestructionEvent.fwd.hh>
#include <core/pose/signals/EnergyEvent.fwd.hh>
#include <core/pose/signals/GeneralEvent.fwd.hh>
#include <core/scoring/ContextGraph.fwd.hh>
#include <core/scoring/ContextGraphTypes.hh>
#include <core/scoring/DerivVectorPair.fwd.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/EnergyGraph.fwd.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/LREnergyContainer.fwd.hh>
#include <core/scoring/MinimizationData.fwd.hh>
#include <core/scoring/MinimizationGraph.fwd.hh>
#include <core/scoring/NeighborList.fwd.hh>
#include <core/scoring/ResidueNeighborList.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionInfo.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/SecondaryStructureWeights.hh>
#include <core/scoring/TenANeighborGraph.fwd.hh>
#include <core/scoring/TwelveANeighborGraph.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/etable/BaseEtableEnergy.hh>
#include <core/scoring/etable/Etable.fwd.hh>
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/etable/EtableEnergy.fwd.hh>
#include <core/scoring/etable/EtableOptions.hh>
#include <core/scoring/etable/count_pair/CountPairAll.fwd.hh>
#include <core/scoring/etable/count_pair/CountPairCrossover3.hh>
#include <core/scoring/etable/count_pair/CountPairCrossover4.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.fwd.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.fwd.hh>
#include <core/scoring/etable/count_pair/CountPairNone.fwd.hh>
#include <core/scoring/etable/count_pair/types.hh>
#include <core/scoring/etable/etrie/CountPairDataGeneric.fwd.hh>
#include <core/scoring/etable/etrie/CountPairData_1_1.fwd.hh>
#include <core/scoring/etable/etrie/CountPairData_1_2.fwd.hh>
#include <core/scoring/etable/etrie/CountPairData_1_3.fwd.hh>
#include <core/scoring/etable/etrie/EtableAtom.fwd.hh>
#include <core/scoring/etable/etrie/EtableTrie.fwd.hh>
#include <core/scoring/etable/etrie/TrieCountPair1BC3.fwd.hh>
#include <core/scoring/etable/etrie/TrieCountPair1BC4.fwd.hh>
#include <core/scoring/etable/etrie/TrieCountPairAll.fwd.hh>
#include <core/scoring/etable/etrie/TrieCountPairGeneric.fwd.hh>
#include <core/scoring/etable/etrie/TrieCountPairNone.fwd.hh>
#include <core/scoring/elec/ElecAtom.fwd.hh>
#include <core/scoring/elec/FA_ElecEnergy.fwd.hh>
#include <core/scoring/hbonds/HBondEnergy.fwd.hh>
#include <core/scoring/hbonds/HBondOptions.fwd.hh>
#include <core/scoring/hbonds/hbtrie/HBAtom.fwd.hh>
#include <core/scoring/hbonds/hbtrie/HBCPData.fwd.hh>
#include <core/scoring/methods/ContextDependentLRTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextDependentOneBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextDependentTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentLRTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>
#include <core/scoring/methods/EnergyMethod.fwd.hh>
#include <core/scoring/methods/EnergyMethod.hh>
#include <core/scoring/methods/EnergyMethodCreator.fwd.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/methods/LongRangeTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/MMLJEnergyInter.fwd.hh>
#include <core/scoring/methods/ShortRangeTwoBodyEnergy.hh>
#include <core/scoring/methods/TwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/TwoBodyEnergy.hh>
#include <core/scoring/methods/WholeStructureEnergy.fwd.hh>
#include <core/scoring/mm/MMBondAngleResidueTypeParamSet.fwd.hh>
#include <core/scoring/mm/mmtrie/MMEnergyTableAtom.fwd.hh>
#include <core/scoring/trie/CPDataCorrespondence.fwd.hh>
#include <core/scoring/trie/RotamerTrie.fwd.hh>
#include <core/scoring/trie/RotamerTrieBase.fwd.hh>
#include <core/scoring/trie/RotamerTrieBase.hh>
#include <core/scoring/trie/TrieCollection.fwd.hh>
#include <core/scoring/trie/TrieCountPairBase.fwd.hh>
#include <core/scoring/trie/trie.functions.fwd.hh>
#include <utility/Bound.fwd.hh>
#include <utility/Bound.hh>
#include <utility/py/PyAssert.hh>
#include <utility/assert.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
#include <utility/in_place_list.fwd.hh>
#include <utility/in_place_list.hh>
#include <utility/stream_util.hh>
#include <utility/string_util.hh>
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
#include <numeric/xyzTriple.fwd.hh>
#include <numeric/xyzTriple.hh>
#include <numeric/xyzVector.fwd.hh>
#include <numeric/xyzVector.hh>
#include <numeric/internal/ColPointers.hh>
#include <numeric/internal/ColVectors.hh>
#include <numeric/internal/ColsPointer.hh>
#include <numeric/internal/RowPointers.hh>
#include <numeric/internal/RowVectors.hh>
#include <numeric/internal/RowsPointer.hh>
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
#include <ObjexxFCL/FArray3.hh>
#include <ObjexxFCL/FArray3A.fwd.hh>
#include <ObjexxFCL/FArray3D.fwd.hh>
#include <ObjexxFCL/FArray3D.hh>
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
#include <utility/assert.hh>
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
#include <set>
#include <sstream>
#include <string>
#include <time.h>
#include <typeinfo>
#include <utility>
#include <vector>
#include <basic/MetricValue.fwd.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableData.fwd.hh>
#include <basic/datacache/CacheableData.hh>
#include <basic/datacache/DataCache.fwd.hh>
#include <basic/datacache/DataCache.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/prof.hh>
#include <boost/algorithm/string/erase.hpp>
#include <boost/bind.hpp>
#include <boost/config.hpp>
#include <boost/function.hpp>
#include <boost/functional/hash.hpp>
#include <boost/pool/detail/mutex.hpp>
#include <boost/pool/poolfwd.hpp>
#include <boost/unordered_map.hpp>
#include <core/pose/util.tmpl.hh>

namespace core {
namespace scoring {
namespace etable {

#ifdef APL_TEMP_DEBUG
inline
Size & mingraph_n_atpairE_evals()
{
	static Size mingraph_n_atpairE_evals_;
	return mingraph_n_atpairE_evals_;
}
#endif

template< class Evaluator >
class ResResEnergyInvoker : public count_pair::Invoker
{
public:
	ResResEnergyInvoker(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		Evaluator const & evaluator,
		EnergyMap & emap
	) :
		rsd1_( rsd1 ),
		rsd2_( rsd2 ),
		evaluator_( evaluator ),
		emap_( emap )
	{}

	virtual ~ResResEnergyInvoker() {}

protected:

	conformation::Residue const & rsd1() { return rsd1_; }
	conformation::Residue const & rsd2() { return rsd2_; }
	Evaluator const & evaluator() { return evaluator_; }
	EnergyMap & emap() { return emap_; }

private:
	conformation::Residue const & rsd1_;
	conformation::Residue const & rsd2_;
	Evaluator const & evaluator_;
	EnergyMap & emap_;
};

template< class Evaluator >
class WholeWholeEnergyInvoker : public ResResEnergyInvoker< Evaluator >
{
public:
	typedef ResResEnergyInvoker< Evaluator > parent;
public:
	WholeWholeEnergyInvoker(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		Evaluator const & evaluator,
		EnergyMap & emap
	)
	: parent( rsd1, rsd2, evaluator, emap )
	{}

	virtual ~WholeWholeEnergyInvoker() {}

	virtual void invoke( count_pair::CountPairFunction const & cp )
	{
		parent::evaluator().residue_atom_pair_energy( parent::rsd1(), parent::rsd2(), cp, parent::emap() );
	}

};

template< class Evaluator >
class SC_BB_EnergyInvoker : public ResResEnergyInvoker< Evaluator >
{
public:
	typedef ResResEnergyInvoker< Evaluator > parent;
public:
	SC_BB_EnergyInvoker(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		Evaluator const & evaluator,
		EnergyMap & emap
	)
	: parent( rsd1, rsd2, evaluator, emap )
	{}

	virtual ~SC_BB_EnergyInvoker() {}

	virtual void invoke( count_pair::CountPairFunction const & cp )
	{
		parent::evaluator().residue_atom_pair_energy_sidechain_backbone( parent::rsd1(), parent::rsd2(), cp, parent::emap() );
	}

};

template< class Evaluator >
class SC_Whole_EnergyInvoker : public ResResEnergyInvoker< Evaluator >
{
public:
	typedef ResResEnergyInvoker< Evaluator > parent;
public:
	SC_Whole_EnergyInvoker(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		Evaluator const & evaluator,
		EnergyMap & emap
	)
	: parent( rsd1, rsd2, evaluator, emap )
	{}

	virtual ~SC_Whole_EnergyInvoker() {}

	virtual void invoke( count_pair::CountPairFunction const & cp )
	{
		parent::evaluator().residue_atom_pair_energy_sidechain_whole( parent::rsd1(), parent::rsd2(), cp, parent::emap() );
	}

};

template< class Evaluator >
class BB_BB_EnergyInvoker : public ResResEnergyInvoker< Evaluator >
{
public:
	typedef ResResEnergyInvoker< Evaluator > parent;
public:
	BB_BB_EnergyInvoker(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		Evaluator const & evaluator,
		EnergyMap & emap
	)
	: parent( rsd1, rsd2, evaluator, emap )
	{}

	virtual ~BB_BB_EnergyInvoker() {}

	virtual void invoke( count_pair::CountPairFunction const & cp )
	{
		parent::evaluator().residue_atom_pair_energy_backbone_backbone( parent::rsd1(), parent::rsd2(), cp, parent::emap() );
	}

};

template< class Evaluator >
class SC_SC_EnergyInvoker : public ResResEnergyInvoker< Evaluator >
{
public:
	typedef ResResEnergyInvoker< Evaluator > parent;
public:
	SC_SC_EnergyInvoker(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		Evaluator const & evaluator,
		EnergyMap & emap
	)
	: parent( rsd1, rsd2, evaluator, emap )
	{}

	virtual ~SC_SC_EnergyInvoker() {}

	virtual void invoke( count_pair::CountPairFunction const & cp )
	{
		parent::evaluator().residue_atom_pair_energy_sidechain_sidechain( parent::rsd1(), parent::rsd2(), cp, parent::emap() );
	}

};


using namespace etrie;
using namespace trie;
using namespace basic::options;

/// construction with an etable
template < class Derived >
BaseEtableEnergy< Derived >::BaseEtableEnergy(
	methods::EnergyMethodCreatorOP creator,
	Etable const & etable_in,
	methods::EnergyMethodOptions const & options,
	bool do_classic_intrares
):
parent( creator ),
etable_( etable_in ),
safe_max_dis2( etable_in.get_safe_max_dis2() ),
hydrogen_interaction_cutoff2_( option[ OptionKeys::score::fa_Hatr ] ?
std::pow( std::sqrt( etable_in.hydrogen_interaction_cutoff2()) + std::sqrt( safe_max_dis2 ), 2)
: etable_in.hydrogen_interaction_cutoff2() ),
exclude_DNA_DNA_( options.exclude_DNA_DNA() ),
do_classic_intrares_( do_classic_intrares ),
put_intra_into_total_( options.put_intra_into_total() ),
exclude_intra_res_protein_( options.exclude_intra_res_protein() )
{}

/// @details an explicit copy constructor is required so that the etable_evaluator_ instance,
/// which is held in an owning pointer, is not shared between multiple instances of this (or rather
/// the derived) class.
template < class Derived >
BaseEtableEnergy< Derived >::BaseEtableEnergy( BaseEtableEnergy< Derived > const & src ) :
parent( src ),
etable_( src.etable_ ),
safe_max_dis2( src.safe_max_dis2 ),
hydrogen_interaction_cutoff2_( src.hydrogen_interaction_cutoff2_ ),
exclude_DNA_DNA_( src.exclude_DNA_DNA_ ),
do_classic_intrares_( src.do_classic_intrares_ ),
put_intra_into_total_( src.put_intra_into_total_ ),
exclude_intra_res_protein_( src.exclude_intra_res_protein_ )
{
}


///////////////////////////////////////////////////////////////////////////////
template < class Derived >
void
BaseEtableEnergy< Derived >::finalize_total_energy(
	pose::Pose & pose,
	ScoreFunction const &,
	EnergyMap & totals
) const
{
#ifdef APL_TEMP_DEBUG
	if ( mingraph_n_atpairE_evals() != 0 ) {
		std::cout << "mingraph_n_atpairE_evals: " << mingraph_n_atpairE_evals() << std::endl;
	}
	mingraph_n_atpairE_evals() = 0;
#endif

	/// Only use this finalize_total_energies routine if we're using the long-distance minimization
	/// auto-update trick
	if ( pose.energies().use_nblist() && pose.energies().use_nblist_auto_update() ) {
#ifdef APL_TEMP_DEBUG
		Size n_d2( 0 );
#endif

		EnergyMap tbenergy_map;
		// add in contributions from the nblist atom-pairs
		NeighborList const & nblist
			( pose.energies().nblist( nblist_type() ) );

		nblist.check_domain_map( pose.energies().domain_map() );

		/// Trick to avoid calls to Conformation::residue()
		utility::vector1< conformation::Residue const * > resvect;
		resvect.reserve( pose.size() );
		for ( Size ii = 1; ii <= pose.size(); ++ii ) {
			resvect.push_back( & pose.residue( ii ) );
		}
		Real dsq(0.0);
		for ( Size i=1, i_end = pose.size(); i<= i_end; ++i ) {
			conformation::Residue const & ires( *resvect[i] );
			for ( Size ii=1, ii_end=ires.natoms(); ii<= ii_end; ++ii ) {
				/// 1. Iterate across intra-residue atom neighbors if there are any;
				/// 2. Iterate across inter-residue atom neighbors.
				//prepare_for_residue_pair( 1,1, pose ); // set intra-res
				AtomNeighbors const & intranbrs( nblist.intrares_upper_atom_neighbors(i,ii) );
				conformation::Atom const & iatom( ires.atom(ii) );
				for ( AtomNeighbors::const_iterator nbr=intranbrs.begin(),
						nbr_end=intranbrs.end(); nbr!= nbr_end; ++nbr ) {
					Size const jj( nbr->atomno() );
					Real const cp_weight( nbr->weight_func() * nbr->weight() );  //fpd

					conformation::Atom const & jatom( ires.atom(jj) );
					static_cast< Derived const & > (*this).intrares_evaluator().atom_pair_energy( iatom, jatom, cp_weight, tbenergy_map, dsq );
					if (!ires.is_protein()) {
						static_cast< Derived const & > (*this).nonprot_intrares_evaluator().atom_pair_energy( iatom, jatom, cp_weight, tbenergy_map, dsq );
					}
				}

				///prepare_for_residue_pair( 1,2, pose ); // set inter-res
				AtomNeighbors const & nbrs( nblist.upper_atom_neighbors(i,ii) );
				for ( AtomNeighbors::const_iterator nbr=nbrs.begin(),
						nbr_end=nbrs.end(); nbr!= nbr_end; ++nbr ) {
#ifdef APL_TEMP_DEBUG
					++n_d2;
#endif
					Size const  j( nbr->rsd() );
					Size const jj( nbr->atomno() );

					Real const cp_weight( nbr->weight_func()*nbr->weight() );
					conformation::Atom const & jatom( resvect[j]->atom(jj) );
					static_cast< Derived const & > (*this).interres_evaluator().atom_pair_energy( iatom, jatom, cp_weight, tbenergy_map, dsq );
				}
			}

		}
#ifdef APL_TEMP_DEBUG
		std::cout << "n_at_pairE_evals: " << n_d2 << std::endl;
#endif
		//std::cout << "totals[ fa_atr ] : " << totals[ fa_atr ] << std::endl;
		//std::cout << "totals[ fa_rep ] : " << totals[ fa_rep ] << std::endl;
		//std::cout << "totals[ fa_sol ] : " << totals[ fa_sol ] << std::endl;
		//std::cout << "totals[ fa_intra_atr ] : " << totals[ fa_intra_atr ] << std::endl;
		//std::cout << "totals[ fa_intra_rep ] : " << totals[ fa_intra_rep ] << std::endl;
		//std::cout << "totals[ fa_intra_sol ] : " << totals[ fa_intra_sol ] << std::endl;
		totals += tbenergy_map;
	}

}

template < class Derived >
bool
BaseEtableEnergy< Derived >::minimize_in_whole_structure_context( pose::Pose const & pose ) const
{
	return pose.energies().use_nblist_auto_update();
}

///////////////////////////////////////////////////////////////////////////////
template < class Derived >
void
BaseEtableEnergy< Derived >::setup_for_minimizing(
	pose::Pose & pose,
	ScoreFunction const & sfxn,
	kinematics::MinimizerMapBase const & min_map
) const
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( pose.energies().use_nblist() && pose.energies().use_nblist_auto_update() ) {
		// stash our nblist inside the pose's energies object
		Energies & energies( pose.energies() );

		// setup the atom-atom nblist
		NeighborListOP nblist;

		if ( pose.energies().use_nblist_auto_update() ) {
			Real const tolerated_narrow_nblist_motion = option[ run::nblist_autoupdate_narrow ];
			Real XX, HH, XH;

			if ( option[ score::fa_Hatr ]() ) {
				HH = XH = XX = etable_.max_dis() + 2 * tolerated_narrow_nblist_motion;
			} else {
				XX = etable_.max_dis() + 2 * tolerated_narrow_nblist_motion;
				XH = etable_.max_non_hydrogen_lj_radius() + etable_.max_hydrogen_lj_radius()
					+ 2 * tolerated_narrow_nblist_motion;
				HH = etable_.max_hydrogen_lj_radius() + etable_.max_hydrogen_lj_radius()
					+ 2 * tolerated_narrow_nblist_motion;
			}

			nblist = NeighborListOP( new NeighborList(
				min_map.domain_map(),
				XX*XX,
				XH*XH,
				HH*HH) );
			nblist->set_auto_update( tolerated_narrow_nblist_motion );
		} else {
			/// Use the default parameters
			nblist = NeighborListOP( new NeighborList(
				min_map.domain_map(),
				etable_.nblist_dis2_cutoff_XX(),
				etable_.nblist_dis2_cutoff_XH(),
				etable_.nblist_dis2_cutoff_HH()) );
		}
		// this partially becomes the EtableEnergy classes's responsibility
		nblist->setup( pose, sfxn, static_cast<Derived const&> (*this));

		energies.set_nblist( nblist_type(), nblist );
	}
}

// check compatibility with atomtypeset
template < class Derived >
void
BaseEtableEnergy< Derived >::setup_for_scoring(
	pose::Pose &pose,
	ScoreFunction const &scfxn
) const
{
	debug_assert( dynamic_cast< Derived const* > (this) );
	Derived const * ptr = static_cast< Derived const * > (this);
	ptr->setup_for_scoring_(pose,scfxn);
	if ( pose.energies().use_nblist() ) {
		NeighborList const & nblist( pose.energies().nblist( nblist_type() ) );
		//std::cout << "nblist autoupdate" << std::endl;
		nblist.prepare_for_scoring( pose, scfxn, *ptr );
	}
}

/// @details Make sure that the neighborlist is up-to-date bfore evaluating derivatives
template < class Derived >
void
BaseEtableEnergy< Derived >::setup_for_derivatives(
	pose::Pose &pose,
	ScoreFunction const &scfxn
) const
{
	//std::cout << "BaseEtable Setup for derivatives" << std::endl;
	debug_assert( dynamic_cast< Derived const* > (this) );
	Derived const * ptr = static_cast< Derived const* > (this);
	ptr->setup_for_scoring_(pose,scfxn);
	if ( pose.energies().use_nblist() ) {
		NeighborList const & nblist( pose.energies().nblist( nblist_type() ) );
		//std::cout << "nblist autoupdate2" << std::endl;
		nblist.prepare_for_scoring( pose, scfxn, *ptr );
	}
}


// The EtableEnergy method stores a vector of rotamer trie objects in the Energies
// object for use in rapid rotamer/background energy calculations.  Overrides default
// do-nothing behavior.
template < class Derived >
void
BaseEtableEnergy< Derived >::setup_for_packing(
	pose::Pose & pose,
	utility::vector1< bool > const &,
	utility::vector1< bool > const &
) const
{

	TrieCollectionOP tries( new TrieCollection );
	tries->total_residue( pose.size() );
	for ( Size ii = 1; ii <= pose.size(); ++ii ) {
		// Do not compute energy for virtual residues.
		if ( pose.residue(ii).aa() == core::chemical::aa_vrt ) continue;
		//iwd  Temporary hack: also skip ligand residues
		//if ( !pose.residue(ii).is_polymer() )
		// continue;

		EtableRotamerTrieOP one_rotamer_trie = create_rotamer_trie( pose.residue(ii), pose );
		tries->trie( ii, one_rotamer_trie );
	}
	pose.energies().data().set( EnergiesCacheableDataType::ETABLE_TRIE_COLLECTION, tries );

}

// @brief Creates a rotamer trie for the input set of rotamers and stores the trie
// in the rotamer set.
template < class Derived >
void
BaseEtableEnergy< Derived >::prepare_rotamers_for_packing(
	pose::Pose const & pose,
	conformation::RotamerSetBase & set
) const
{

	EtableRotamerTrieOP rottrie = create_rotamer_trie( set, pose );
	set.store_trie( methods::etable_method, rottrie );
}


// @brief Updates the cached rotamer trie for a residue if it has changed during the course of
// a repacking
template < class Derived >
void
BaseEtableEnergy< Derived >::update_residue_for_packing(
	pose::Pose & pose,
	Size resid
) const
{

	EtableRotamerTrieOP one_rotamer_trie = create_rotamer_trie( pose.residue( resid ), pose );

	// grab non-const & of the cached tries and replace resid's trie with a new one.
	TrieCollection & trie_collection
		( static_cast< TrieCollection & > (pose.energies().data().get( EnergiesCacheableDataType::ETABLE_TRIE_COLLECTION )));
	trie_collection.trie( resid, one_rotamer_trie );

}


// @brief convenience function
template < class Derived >
count_pair::CountPairFunctionCOP
BaseEtableEnergy< Derived >::get_count_pair_function(
	Size res1,
	Size res2,
	pose::Pose const & pose,
	ScoreFunction const & sfxn
) const
{
	return get_count_pair_function( pose.residue( res1 ), pose.residue( res2 ), pose, sfxn );
}

// @brief all of the residue level count pair logic goes here
/// and gets duplicated in a handful of other places within this class... this NEEDS to be reworked.
template < class Derived >
count_pair::CountPairFunctionCOP
BaseEtableEnergy< Derived >::get_count_pair_function(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	pose::Pose const & pose, // does the pose hold the connectivity information or do the residues?
	ScoreFunction const & sfxn
) const
{
	using namespace count_pair;

	if ( !calculate_interres( res1, res2 ) ) {
		return count_pair::CountPairFunctionCOP( count_pair::CountPairFunctionOP( new CountPairNone ) );
	}

	count_pair::CPCrossoverBehavior crossover = determine_crossover_behavior( res1, res2, pose, sfxn );
	return CountPairFactory::create_count_pair_function( res1, res2, crossover );
}

// @brief get a count pair object for intraresidue pair energies
template < class Derived >
count_pair::CountPairFunctionOP
BaseEtableEnergy< Derived >::get_intrares_countpair(
	conformation::Residue const & res,
	pose::Pose const & pose,
	ScoreFunction const & sfxn
) const
{
	using namespace count_pair;

	if ( !calculate_intrares( res ) ) {
		return count_pair::CountPairFunctionOP( new CountPairNone );
	}

	CPCrossoverBehavior crossover = determine_crossover_behavior( res, res, pose, sfxn );
	return CountPairFactory::create_intrares_count_pair_function( res, crossover );
}


/// @brief figure out the trie count pair function to use
/// Need to refactor this so that the decision "what kind of count pair behavior should I use" can be decoupled
/// from class instantiation, and therefore shared between the creation of the trie count pair classes and the regular
/// count pair classes
template < class Derived >
TrieCountPairBaseOP
BaseEtableEnergy< Derived >::get_count_pair_function_trie(
	conformation::RotamerSetBase const & set1,
	conformation::RotamerSetBase const & set2,
	pose::Pose const & pose,
	ScoreFunction const & sfxn
) const
{
	using namespace methods;
	conformation::Residue const & res1( pose.residue( set1.resid() ) );
	conformation::Residue const & res2( pose.residue( set2.resid() ) );
	trie::RotamerTrieBaseCOP trie1( utility::pointer::static_pointer_cast< trie::RotamerTrieBase const > ( set1.get_trie( etable_method ) ));
	trie::RotamerTrieBaseCOP trie2( utility::pointer::static_pointer_cast< trie::RotamerTrieBase const > ( set2.get_trie( etable_method ) ));

	return get_count_pair_function_trie( res1, res2, trie1, trie2, pose, sfxn );
}

template < class Derived >
TrieCountPairBaseOP
BaseEtableEnergy< Derived >::get_count_pair_function_trie(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	trie::RotamerTrieBaseCOP trie1,
	trie::RotamerTrieBaseCOP trie2,
	pose::Pose const & pose,
	ScoreFunction const & sfxn
) const
{
	using namespace count_pair;

	TrieCountPairBaseOP tcpfxn;
	if ( !calculate_interres( res1, res2 ) ) {
		tcpfxn = TrieCountPairBaseOP( new TrieCountPairNone() );
		return tcpfxn;
	}

	CPResidueConnectionType connection = count_pair::CountPairFactory::determine_residue_connection( res1, res2 );
	Size conn1 = trie1->get_count_pair_data_for_residue( res2.seqpos() );
	Size conn2 = trie2->get_count_pair_data_for_residue( res1.seqpos() );

	if ( connection == CP_ONE_BOND ) {
		CPCrossoverBehavior crossover = determine_crossover_behavior( res1, res2, pose, sfxn );
		switch ( crossover ) {
		case CP_CROSSOVER_3 :
			tcpfxn = TrieCountPairBaseOP( new TrieCountPair1BC3( conn1, conn2 ) );
			break;
		case CP_CROSSOVER_4 :
			tcpfxn = TrieCountPairBaseOP( new TrieCountPair1BC4( conn1, conn2 ) );
			break;
		default :
			utility_exit();
			break;
		}
	} else if ( connection == CP_MULTIPLE_BONDS_OR_PSEUDOBONDS ) {
		CPCrossoverBehavior crossover = determine_crossover_behavior( res1, res2, pose, sfxn );

		TrieCountPairGenericOP cpgen( new TrieCountPairGeneric( res1, res2, conn1, conn2 ) );
		if ( crossover == CP_CROSSOVER_3 ) {
			cpgen->crossover( 3 );
		} else if ( crossover == CP_CROSSOVER_4 ) {
			cpgen->crossover( 4 );
		} else {
			utility_exit();
		}
		cpgen->hard_crossover( false );
		tcpfxn = cpgen;
	} else {
		tcpfxn = TrieCountPairBaseOP( new TrieCountPairAll );
	}
	return tcpfxn;

}


/*
template < class Derived >
count_pair::CPResidueConnectionType
BaseEtableEnergy< Derived >::determine_residue_connection(
conformation::Residue const & res1,
conformation::Residue const & res2,
pose::Pose const &
) const
{
using namespace count_pair;
/// this code is incompatible with designing both disulfides and non-disulfies
/// at the same residue...
if ( res1.is_pseudobonded( res2.seqpos ) {
return CP_MULTIPLE_BONDS_OR_PSEUDOBONDS;
} else if ( res1.is_bonded(res2) ) {
if ( res1.connections_to_residue( res2 ).size() == 1 ) {
return CP_ONE_BOND;
} else {
return CP_MULTIPLE_BONDS_OR_PSEUDOBONDS;
}
} else {
return CP_NO_BONDS;
}
}
*/

// Note: The core::scoring::lkball::determine_crossover_behavior() function seems to copy code from this function.  There is a note there
// indicating that if this function is changed, that function should be changed correspondingly.  Why the author didn't add that note
// here is beyond me -- or, better yet, have both functions call a common utility function in order to avoid code duplication and drift.
// --VKM, 8 March 2016.
template < class Derived >
count_pair::CPCrossoverBehavior
BaseEtableEnergy< Derived >::determine_crossover_behavior(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	pose::Pose const &,
	ScoreFunction const & sfxn
) const
{
	using namespace count_pair;

	// Ask "are these two residues polymers and do they share a polymeric bond?"
	if (
			//res1.polymeric_sequence_distance(res2) == 1 //VKM, 20 Feb 2016: Doesn't handle cyclic geometry properly.  The following code does, however:
			res1.type().is_polymer() && res2.type().is_polymer() && res1.is_polymer_bonded(res2)
			) {
		if ( ( !sfxn.has_zero_weight( mm_twist ) && sfxn.has_zero_weight( rama )) ||
				( ( !res1.is_protein() || !res2.is_protein()) &&
				( !res1.is_RNA() || !res2.is_RNA())
				) ) {
			return CP_CROSSOVER_3;
		} else {
			return CP_CROSSOVER_4; // peptide bond w/ or w/o rama, but definately w/o mm_twist
		}
	} else if ( res1.seqpos() == res2.seqpos() ) {
		// logic for controlling intra-residue count pair behavior goes here; for now, default to crossover 3
		if ( do_classic_intrares_ ) {
			return CP_CROSSOVER_3; // used by EtableIntraClassicEnergy
		} else {
			return CP_CROSSOVER_4;
		}
	} else {
		return CP_CROSSOVER_3; // e.g. disulfides where seqsep != 1
	}
}


template < class Derived >
void
BaseEtableEnergy< Derived >::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap & emap
) const
{
	debug_assert( rsd1.seqpos() != rsd2.seqpos() );
	//std::cerr << __FILE__<< ' ' << __LINE__ << std::endl;

	if ( ! pose.energies().use_nblist() ) {
		prepare_for_residue_pair(rsd1.seqpos(),rsd2.seqpos(),pose);
		//count_pair::CountPairFunctionCOP cpfxn = get_count_pair_function( rsd1, rsd2, pose, sfxn );
		//cpfxn->residue_atom_pair_energy( rsd1, rsd2, static_cast<Derived const&> (*this), emap );
		if ( !calculate_interres( rsd1, rsd2 ) ) return;
		count_pair::CPCrossoverBehavior crossover = determine_crossover_behavior( rsd1, rsd2, pose, sfxn );
		WholeWholeEnergyInvoker< typename Derived::Evaluator > invoker( rsd1, rsd2, static_cast< Derived const & > (*this).interres_evaluator(), emap );
		count_pair::CountPairFactory::create_count_pair_function_and_invoke( rsd1, rsd2, crossover, invoker );
	}
}

template < class Derived >
bool
BaseEtableEnergy< Derived >::use_extended_residue_pair_energy_interface() const
{
	return true;
}

template < class Derived >
void
BaseEtableEnergy< Derived >::residue_pair_energy_ext(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResPairMinimizationData const & min_data,
	pose::Pose const & pose,
	ScoreFunction const &,
	EnergyMap & emap
) const
{

	if ( pose.energies().use_nblist_auto_update() ) return;

	///prepare_for_residue_pair( 1,2, pose ); // set inter-res
	debug_assert( utility::pointer::dynamic_pointer_cast< ResiduePairNeighborList const > (min_data.get_data( min_pair_data_type() ) ));
	ResiduePairNeighborList const & nblist( static_cast< ResiduePairNeighborList const & > ( min_data.get_data_ref( min_pair_data_type() ) ) );
	Real dsq;
	utility::vector1< SmallAtNb > const & neighbs( nblist.atom_neighbors() );
	for ( Size ii = 1, iiend = neighbs.size(); ii <= iiend; ++ii ) {
		conformation::Atom const & atom1( rsd1.atom( neighbs[ ii ].atomno1() ) );
		conformation::Atom const & atom2( rsd2.atom( neighbs[ ii ].atomno2() ) );
		static_cast< Derived const & > (*this).interres_evaluator().atom_pair_energy( atom1, atom2, neighbs[ ii ].weight(), emap, dsq );
#ifdef APL_TEMP_DEBUG
		++mingraph_n_atpairE_evals();
#endif
	}
}

template < class Derived >
void
BaseEtableEnergy< Derived >::setup_for_minimizing_for_residue_pair(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	kinematics::MinimizerMapBase const &, // is this necessary?
	ResSingleMinimizationData const &,//, r1_min_dat,
	ResSingleMinimizationData const &,// r2_min_dat,
	ResPairMinimizationData & pair_data
) const
{

	if ( pose.energies().use_nblist_auto_update() ) return;

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	count_pair::CountPairFunctionCOP count_pair = get_count_pair_function( rsd1, rsd2, pose, sfxn );

	// update the existing nblist if it's already present in the min_data object
	ResiduePairNeighborListOP nblist( utility::pointer::static_pointer_cast< core::scoring::ResiduePairNeighborList > ( pair_data.get_data( min_pair_data_type() ) ));
	if ( ! nblist ) nblist = ResiduePairNeighborListOP( new ResiduePairNeighborList );

	/// STOLEN CODE!
	Real const tolerated_narrow_nblist_motion = 0.75; //option[ run::nblist_autoupdate_narrow ];
	Real const XX = etable_.max_dis() + 2 * tolerated_narrow_nblist_motion;
	Real const XH = etable_.max_non_hydrogen_lj_radius() + etable_.max_hydrogen_lj_radius()
		+ 2 * tolerated_narrow_nblist_motion;
	Real const HH = etable_.max_hydrogen_lj_radius() + etable_.max_hydrogen_lj_radius()
		+ 2 * tolerated_narrow_nblist_motion;

	nblist->initialize_from_residues(
		XX*XX, XH*XH, HH*HH,
		rsd1, rsd2, count_pair );

	pair_data.set_data( min_pair_data_type(), nblist );

}

template < class Derived >
bool
BaseEtableEnergy< Derived >::requires_a_setup_for_scoring_for_residue_opportunity( pose::Pose const &  ) const
{
	// if ( pose.energies().use_nblist_auto_update() ) return false;

	//return true;
	return false; // TEMP -- disable autoupdate
}


template < class Derived >
void
BaseEtableEnergy< Derived >::setup_for_scoring_for_residue(
	conformation::Residue const &,// rsd,
	pose::Pose const &,
	ScoreFunction const &,
	ResSingleMinimizationData & // min_data
) const
{}
/*debug_assert( dynamic_cast< ResidueNblistData * > (min_data.get_data( min_single_data_type() )() ));
ResidueNblistData & nbdata( static_cast< ResidueNblistData & > ( min_data.get_data_ref( min_single_data_type() ) ));
nbdata.update( rsd );
}*/


template < class Derived >
bool
BaseEtableEnergy< Derived >::requires_a_setup_for_derivatives_for_residue_opportunity( pose::Pose const & ) const
{
	///if ( pose.energies().use_nblist_auto_update() ) return false;

	//return true;
	return false; // disable autoupdate (TEMP!)
}


template < class Derived >
void
BaseEtableEnergy< Derived >::setup_for_derivatives_for_residue(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	ResSingleMinimizationData & min_data
) const
{
	setup_for_scoring_for_residue( rsd, pose, sfxn, min_data );
}


template < class Derived >
bool
BaseEtableEnergy< Derived >::requires_a_setup_for_scoring_for_residue_pair_opportunity( pose::Pose const & ) const
{
	//return true;
	return false; /// TEMP -- disable autoupdate
}


template < class Derived >
void
BaseEtableEnergy< Derived >::setup_for_scoring_for_residue_pair(
	conformation::Residue const &,// rsd1,
	conformation::Residue const &,// rsd2,
	ResSingleMinimizationData const &,// minsingle_data1,
	ResSingleMinimizationData const &,// minsingle_data2,
	pose::Pose const &,
	ScoreFunction const &,
	ResPairMinimizationData &// data_cache
) const
{
	/*debug_assert( dynamic_cast< ResidueNblistData const * > (minsingle_data1.get_data( min_single_data_type() )() ));
	debug_assert( dynamic_cast< ResidueNblistData const * > (minsingle_data2.get_data( min_single_data_type() )() ));
	debug_assert( dynamic_cast< ResiduePairNeighborList * > (data_cache.get_data( min_pair_data_type() )() ));

	ResidueNblistData const & r1nbdat( static_cast< ResidueNblistData const & > ( minsingle_data1.get_data_ref( min_single_data_type() ) ));
	ResidueNblistData const & r2nbdat( static_cast< ResidueNblistData const & > ( minsingle_data2.get_data_ref( min_single_data_type() ) ));
	ResiduePairNeighborList & nblist( static_cast< ResiduePairNeighborList & > (  data_cache.get_data_ref( min_pair_data_type() ) ) );

	nblist.update(
	etable_.max_heavy_heavy_cutoff(),
	etable_.max_heavy_hydrogen_cutoff(),
	etable_.max_hydrogen_hydrogen_cutoff(),
	rsd1, rsd2, r1nbdat, r2nbdat );*/

}


template < class Derived >
bool
BaseEtableEnergy< Derived >::requires_a_setup_for_derivatives_for_residue_pair_opportunity( pose::Pose const & ) const
{
	return false;
}


template < class Derived >
void
BaseEtableEnergy< Derived >::setup_for_derivatives_for_residue_pair(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResSingleMinimizationData const & minsingle_data1,
	ResSingleMinimizationData const & minsingle_data2,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	ResPairMinimizationData & data_cache
) const
{
	setup_for_scoring_for_residue_pair( rsd1, rsd2, minsingle_data1, minsingle_data2, pose, sfxn, data_cache );
}

/*
template < class Derived >
void
BaseEtableEnergy< Derived >::eval_atom_derivative_for_residue_pair(
Size const atom_index,
conformation::Residue const & rsd1,
conformation::Residue const & rsd2,
ResSingleMinimizationData const &,
ResSingleMinimizationData const &,
ResPairMinimizationData const & min_data,
pose::Pose const & pose, // provides context
kinematics::DomainMap const &,
ScoreFunction const &,
EnergyMap const & weights,
Vector & F1,
Vector & F2
) const
{
debug_assert( dynamic_cast< ResiduePairNeighborList const * > (min_data.get_data( min_pair_data_type() )() ));
ResiduePairNeighborList const & nblist( static_cast< ResiduePairNeighborList const & > (min_data.get_data_ref( min_pair_data_type() )) );
utility::vector1< AtomNeighbor > const & atneighbs( rsd1.seqpos() < rsd2.seqpos() ?
nblist.r1_at_neighbors()[ atom_index ] :
nblist.r2_at_neighbors()[ atom_index ] );
if ( atneighbs.empty() ) return; // early exit

conformation::Atom const & atom1( rsd1.atom( atom_index ));
prepare_for_residue_pair( 1, 2, pose ); // set inter-res
Vector f1,f2;
//ResiduePairNeighborList::AtomNeighbors const & atneighbs( rsd1.seqpos() < rsd2.seqpos() ?
// nblist.r1_narrow_neighbors( atom_index ) :
// nblist.r2_narrow_neighbors( atom_index ) );
//for ( Size ii = atneighbs.head(), iiend = atneighbs.end(); ii != iiend; ii = atneighbs[ ii ].next() ) {
for ( Size ii = 1, iiend = atneighbs.size(); ii <= iiend; ++ii ) {
//Real const cp_weight( atneighbs[ ii ].data().weight() );  // do not use nbr->weight_func() here
Real const cp_weight( atneighbs[ ii ].weight() );  // do not use nbr->weight_func() here
conformation::Atom const & atom2( rsd2.atom( atneighbs[ ii ].atomno() ) );
Real const dE_dR_over_r( eval_dE_dR_over_r( atom1, atom2, weights, f1, f2 ) );
//std::cout << "  atom deriv: " << rsd1.seqpos() << " " << atom_index << " with " << rsd2.seqpos() << " " << nbr->atomno() << ". w= " << nbr->weight() << " dE_dR_over_r: " << dE_dR_over_r <<  std::endl;
if ( dE_dR_over_r != 0.0 ) {
F1 += dE_dR_over_r * cp_weight * f1;
F2 += dE_dR_over_r * cp_weight * f2;
}
}
}*/

template < class Derived >
void
BaseEtableEnergy< Derived >::eval_residue_pair_derivatives(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResSingleMinimizationData const &,
	ResSingleMinimizationData const &,
	ResPairMinimizationData const & min_data,
	pose::Pose const & pose, // provides context
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & r1_at_derivs,
	utility::vector1< DerivVectorPair > & r2_at_derivs
) const
{

	if ( pose.energies().use_nblist_auto_update() ) return;

	/*debug_assert( r1_at_derivs.size() >= rsd1.natoms() );
	debug_assert( r2_at_derivs.size() >= rsd2.natoms() );

	debug_assert( dynamic_cast< ResiduePairNeighborList const * > (min_data.get_data( min_pair_data_type() )() ));
	ResiduePairNeighborList const & nblist( static_cast< ResiduePairNeighborList const & > (min_data.get_data_ref( min_pair_data_type() )) );

	prepare_for_residue_pair( 1, 2, pose ); // set inter-res
	Vector f1,f2;
	utility::vector1< utility::vector1< AtomNeighbor > > const & r1_nbs = nblist.r1_at_neighbors();
	for ( Size ii = nblist.r1_nonempty().head(),  iiend = 0; ii != iiend; ii = nblist.r1_nonempty()[ ii ].next() ) {
	utility::vector1< AtomNeighbor > const & ii_nbs( r1_nbs[ ii ] );
	conformation::Atom const & iiatom( rsd1.atom( ii ) );
	for ( Size jj = 1, jj_end = ii_nbs.size(); jj <= jj_end; ++jj ) {
	conformation::Atom const & jjatom( rsd2.atom( ii_nbs[ jj ].atomno() ) );
	Real const dE_dR_over_r( eval_dE_dR_over_r( iiatom, jjatom, weights, f1, f2 ) );
	//std::cout << "  atom deriv: " << rsd1.seqpos() << " " << atom_index << " with " << rsd2.seqpos() << " " << nbr->atomno() << ". w= " << nbr->weight() << " dE_dR_over_r: " << dE_dR_over_r <<  std::endl;
	if ( dE_dR_over_r != 0.0 ) {
	f1 *= dE_dR_over_r * ii_nbs[ jj ].weight();
	f2 *= dE_dR_over_r * ii_nbs[ jj ].weight();
	r1_at_derivs[ ii ].f1() += f1;
	r1_at_derivs[ ii ].f2() += f2;
	r2_at_derivs[ ii_nbs[ jj ].atomno() ].f1() += -1*f1;
	r2_at_derivs[ ii_nbs[ jj ].atomno() ].f2() += -1*f2;
	}
	}
	}*/

	debug_assert( r1_at_derivs.size() >= rsd1.natoms() );
	debug_assert( r2_at_derivs.size() >= rsd2.natoms() );

	debug_assert( utility::pointer::dynamic_pointer_cast< ResiduePairNeighborList const > (min_data.get_data( min_pair_data_type() ) ));
	ResiduePairNeighborList const & nblist( static_cast< ResiduePairNeighborList const & > (min_data.get_data_ref( min_pair_data_type() )) );

	//prepare_for_residue_pair( 1, 2, pose ); // set inter-res
	typename Derived::Evaluator evaluator( static_cast< Derived const & > (*this).interres_evaluator() );
	evaluator.set_weights( weights );

	Vector f1,f2;
	utility::vector1< SmallAtNb > const & neighbs( nblist.atom_neighbors() );
	for ( Size ii = 1, iiend = neighbs.size(); ii <= iiend; ++ii ) {
		conformation::Atom const & atom1( rsd1.atom( neighbs[ ii ].atomno1() ) );
		conformation::Atom const & atom2( rsd2.atom( neighbs[ ii ].atomno2() ) );
		Real const dE_dR_over_r( evaluator.eval_dE_dR_over_r( atom1, atom2, weights, f1, f2 ) );
		//std::cout << "  atom deriv: " << rsd1.seqpos() << " " << atom_index << " with " << rsd2.seqpos() << " " << nbr->atomno() << ". w= " << nbr->weight() << " dE_dR_over_r: " << dE_dR_over_r <<  std::endl;
		if ( dE_dR_over_r != 0.0 ) {
			f1 *= dE_dR_over_r * neighbs[ ii ].weight();
			f2 *= dE_dR_over_r * neighbs[ ii ].weight();
			r1_at_derivs[ neighbs[ ii ].atomno1() ].f1() += f1;
			r1_at_derivs[ neighbs[ ii ].atomno1() ].f2() += f2;
			r2_at_derivs[ neighbs[ ii ].atomno2() ].f1() += -1*f1;
			r2_at_derivs[ neighbs[ ii ].atomno2() ].f2() += -1*f2;
		}
	}
}


/// @details do not use this during minimization
template < class Derived >
void
BaseEtableEnergy< Derived >::backbone_backbone_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap & emap
) const
{
	debug_assert( ! pose.energies().use_nblist() );
	//prepare_for_residue_pair(rsd1.seqpos(),rsd2.seqpos(),pose);
	//count_pair::CountPairFunctionCOP cpfxn = get_count_pair_function( rsd1, rsd2, pose, sfxn );
	//cpfxn->residue_atom_pair_energy_backbone_backbone( rsd1, rsd2, static_cast< Derived const&> (*this), emap );
	if ( !calculate_interres( rsd1, rsd2 ) ) return;
	count_pair::CPCrossoverBehavior crossover = determine_crossover_behavior( rsd1, rsd2, pose, sfxn );
	BB_BB_EnergyInvoker< typename Derived::Evaluator > invoker( rsd1, rsd2, static_cast< Derived const & > (*this).interres_evaluator(), emap );
	count_pair::CountPairFactory::create_count_pair_function_and_invoke( rsd1, rsd2, crossover, invoker );
}


/// @details do not use this during minimization
template < class Derived >
void
BaseEtableEnergy< Derived >::backbone_sidechain_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap & emap
) const
{
	debug_assert( ! pose.energies().use_nblist() );
	//prepare_for_residue_pair(rsd2.seqpos(),rsd1.seqpos(),pose);
	//count_pair::CountPairFunctionCOP cpfxn = get_count_pair_function( rsd2, rsd1, pose, sfxn );
	//cpfxn->residue_atom_pair_energy_sidechain_backbone( rsd2, rsd1, static_cast<Derived const&> (*this), emap );
	if ( !calculate_interres( rsd1, rsd2 ) ) return;
	count_pair::CPCrossoverBehavior crossover = determine_crossover_behavior( rsd2, rsd1, pose, sfxn );
	SC_BB_EnergyInvoker< typename Derived::Evaluator > invoker( rsd2, rsd1, static_cast< Derived const & > (*this).interres_evaluator(), emap );
	count_pair::CountPairFactory::create_count_pair_function_and_invoke( rsd2, rsd1, crossover, invoker );
}

//@details do not use this during minimization
template < class Derived >
void
BaseEtableEnergy< Derived >::sidechain_sidechain_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap & emap
) const
{
	debug_assert( ! pose.energies().use_nblist() );
	//prepare_for_residue_pair(rsd1.seqpos(),rsd2.seqpos(),pose);
	//count_pair::CountPairFunctionCOP cpfxn = get_count_pair_function( rsd1, rsd2, pose, sfxn );
	//cpfxn->residue_atom_pair_energy_sidechain_sidechain( rsd1, rsd2, static_cast<Derived const&> (*this), emap );
	if ( !calculate_interres( rsd1, rsd2 ) ) return;
	count_pair::CPCrossoverBehavior crossover = determine_crossover_behavior( rsd1, rsd2, pose, sfxn );
	SC_SC_EnergyInvoker< typename Derived::Evaluator > invoker( rsd1, rsd2, static_cast< Derived const & > (*this).interres_evaluator(), emap );
	count_pair::CountPairFactory::create_count_pair_function_and_invoke( rsd1, rsd2, crossover, invoker );
}


template < class Derived >
void
BaseEtableEnergy< Derived >::evaluate_rotamer_pair_energies(
	conformation::RotamerSetBase const & set1,
	conformation::RotamerSetBase const & set2,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap const & weights,
	ObjexxFCL::FArray2D< core::PackerEnergy > & energy_table
) const
{
	debug_assert( set1.resid() != set2.resid() );

	using namespace methods;
	ObjexxFCL::FArray2D< core::PackerEnergy > temp_table1( energy_table );
	ObjexxFCL::FArray2D< core::PackerEnergy > temp_table2( energy_table );

	temp_table1 = 0; temp_table2 = 0;

	// save weight information so that its available during tvt execution
	//weights_ = weights; <--- getting rid of this non-bitwise-const data.
	typename Derived::Evaluator evaluator( static_cast< Derived const & > (*this).interres_evaluator() );
	evaluator.set_weights( weights );

	EtableRotamerTrieCOP trie1( utility::pointer::static_pointer_cast< trie::RotamerTrieBase const > ( set1.get_trie( etable_method ) ));
	EtableRotamerTrieCOP trie2( utility::pointer::static_pointer_cast< trie::RotamerTrieBase const > ( set2.get_trie( etable_method ) ));

	// figure out which trie countPairFunction needs to be used for this set
	TrieCountPairBaseOP cp = get_count_pair_function_trie( set1, set2, pose, sfxn );

	/// now execute the trie vs trie algorithm.
	/// this steps through three rounds of type resolution before finally arriving at the
	/// actual trie_vs_trie method.  The type resolution calls allow the trie-vs-trie algorithm
	/// to be templated with full type knowledge (and therefore be optimized by the compiler for
	/// each variation on the count pair data used and the count pair funtions invoked.
	evaluator.trie_vs_trie( *trie1, *trie2, *cp, temp_table1, temp_table2 );
	//trie1->trie_vs_trie( *trie2, *cp, static_cast<Derived const&>(*this), temp_table1, temp_table2 );

	/// add in the energies calculated by the tvt alg.
	energy_table += temp_table1;
	//std::cout << "FINISHED evaluate_rotamer_pair_energies" << std::endl;

	/* // There should be a way to turn this on without recompiling...
	// debug
	ObjexxFCL::FArray2D< core::PackerEnergy > temp_table3( energy_table );
	temp_table3 = 0;
	EnergyMap emap;
	for ( Size ii = 1, ii_end = set1.num_rotamers(); ii <= ii_end; ++ii ) {
	for ( Size jj = 1, jj_end = set2.num_rotamers(); jj <= jj_end; ++jj ) {
	emap.zero();
	residue_pair_energy( *set1.rotamer( ii ), *set2.rotamer( jj ), pose, sfxn, emap );
	temp_table3( jj, ii ) += weights.dot( emap );
	if ( std::abs( temp_table1( jj, ii ) - temp_table3( jj, ii )) > 0.001 ) {
	std::cout << "Residues " << set1.resid() << " & " << set2.resid() << " rotamers: " << ii << " & " << jj;
	std::cout << " tvt/reg discrepancy: tvt= " <<  temp_table1( jj, ii ) << " reg= " << temp_table3( jj, ii );
	std::cout << " delta: " << temp_table1( jj, ii ) - temp_table3( jj, ii ) << std::endl;
	}
	}
	}
	//std::cout << "Finished RPE calcs for residues " << set1.resid() << " & " << set2.resid() << std::endl;
	*/
}

template < class Derived >
void
BaseEtableEnergy< Derived >::evaluate_rotamer_background_energies(
	conformation::RotamerSetBase const & set,
	conformation::Residue const & residue,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap const & weights,
	utility::vector1< core::PackerEnergy > & energy_vector
) const
{

	using namespace methods;
	// allocate space for the trie-vs-trie algorithm
	utility::vector1< core::PackerEnergy > temp_vector1( set.num_rotamers(), 0.0 );
	utility::vector1< core::PackerEnergy > temp_vector2( set.num_rotamers(), 0.0 );

	// save weight information so that its available during tvt execution
	typename Derived::Evaluator evaluator( static_cast< Derived const & > (*this).interres_evaluator() );
	evaluator.set_weights( weights );

	EtableRotamerTrieCOP trie1( utility::pointer::static_pointer_cast< trie::RotamerTrieBase const > ( set.get_trie( etable_method ) ));
	EtableRotamerTrieCOP trie2 = ( static_cast< TrieCollection const & >
		( pose.energies().data().get( EnergiesCacheableDataType::ETABLE_TRIE_COLLECTION )) ).trie( residue.seqpos() );

	//fpd
	if ( trie2 == NULL ) return;

	//prepare_for_residue_pair( set.resid(), residue.seqpos(), pose );

	// figure out which trie countPairFunction needs to be used for this set
	TrieCountPairBaseOP cp = get_count_pair_function_trie( pose.residue( set.resid() ), residue, trie1, trie2, pose, sfxn );

	/// now execute the trie vs path algorithm.
	/// this steps through three rounds of type resolution before finally arriving at the
	/// actual trie_vs_path method.  The type resolution calls allow the trie-vs-trie algorithm
	/// to be templated with full type knowledge (and therefore be optimized by the compiler for
	/// each variation on the count pair data used and the count pair funtions invoked.
	/// old! trie1->trie_vs_path( *trie2, *cp, static_cast<Derived const&> (*this), temp_vector1, temp_vector2 );
	evaluator.trie_vs_path( *trie1, *trie2, *cp, temp_vector1, temp_vector2 );

	/// add in the energies calculated by the tvt alg.
	for ( Size ii = 1; ii <= set.num_rotamers(); ++ii ) {
		energy_vector[ ii ] += temp_vector1[ ii ];
	}
	//std::cout << "FINISHED evaluate_rotamer_background_energies" << std::endl;

	/*
	//debug
	utility::vector1< Energy > temp_vector3( energy_vector.size(), 0.0f );
	EnergyMap emap;
	for ( Size ii = 1, ii_end = set.num_rotamers(); ii <= ii_end; ++ii ) {
	emap.zero();
	residue_pair_energy( *set.rotamer( ii ), residue, pose, sfxn, emap );
	temp_vector3[ ii ] += weights.dot( emap );
	if ( std::abs( temp_vector1[ ii ] - temp_vector3[ ii ]) > 0.001 ) {
	std::cout << "Residues " << set.resid() << " & " << residue.seqpos() << " rotamers: " << ii << " & bg";
	std::cout << " tvt/reg discrepancy: tvt= " <<  temp_vector1[ ii ] << " reg= " << temp_vector3[ ii ];
	std::cout << " delta: " << temp_vector1[ ii ] - temp_vector3[ ii ] << std::endl;
	}
	}
	std::cout << "Finished Rotamer BG calcs for residues " << set.resid() << " & " << residue.seqpos() << std::endl;
	*/
}


template < class Derived >
bool
BaseEtableEnergy< Derived >::use_extended_intrares_energy_interface() const
{
	return true;
}

template < class Derived >
void
BaseEtableEnergy< Derived >::eval_intrares_energy_ext(
	conformation::Residue const & rsd,
	ResSingleMinimizationData const & min_data,
	pose::Pose const & pose,
	ScoreFunction const &,
	EnergyMap & emap
) const
{
	if ( pose.energies().use_nblist_auto_update() ) return;

	debug_assert( utility::pointer::dynamic_pointer_cast< ResidueNblistData const > (min_data.get_data( min_single_data_type() ) ));
	ResidueNblistData const & nblist( static_cast< ResidueNblistData const & > ( min_data.get_data_ref( min_single_data_type() ) ) );

	//prepare_for_residue_pair( 1,1, pose ); // set intra-res
	utility::vector1< SmallAtNb > const & neighbs( nblist.atom_neighbors() );

	Real dsq;
	for ( Size ii = 1, iiend = neighbs.size(); ii <= iiend; ++ii ) {
		conformation::Atom const & atom1( rsd.atom( neighbs[ ii ].atomno1() ) );
		conformation::Atom const & atom2( rsd.atom( neighbs[ ii ].atomno2() ) );
		static_cast< Derived const & > (*this).intrares_evaluator().atom_pair_energy( atom1, atom2, neighbs[ ii ].weight(), emap, dsq );
		if (!rsd.is_protein()) {
			static_cast< Derived const & > (*this).nonprot_intrares_evaluator().atom_pair_energy( atom1, atom2, neighbs[ ii ].weight(), emap, dsq );
		}

		//std::cout << "evaluated " << neighbs[ ii ].atomno1() << " " << neighbs[ ii ].atomno2() << " " <<
		// emap[ fa_atr ] << " " << emap[ fa_rep ] << " " << emap[ fa_sol ] << " " <<
		// emap[ fa_intra_atr ] << " " << emap[ fa_intra_rep ] << " " << emap[ fa_intra_sol ] << std::endl;
	}
}

template < class Derived >
void
BaseEtableEnergy< Derived >::setup_for_minimizing_for_residue(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	kinematics::MinimizerMapBase const &,
	ResSingleMinimizationData & min_data
) const
{
	if ( pose.energies().use_nblist_auto_update() ) return;

	if ( (static_cast< Derived const & > (*this)).defines_intrares_energy( sfxn.weights() ) ) {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		// update the existing nblist if it's already present in the min_data object
		ResidueNblistDataOP nbdata( utility::pointer::static_pointer_cast< core::scoring::ResidueNblistData > ( min_data.get_data( min_single_data_type() ) ));
		if ( ! nbdata ) nbdata = ResidueNblistDataOP( new ResidueNblistData );

		/// STOLEN CODE!
		Real const tolerated_narrow_nblist_motion = option[ run::nblist_autoupdate_narrow ];
		Real const XX = etable_.max_dis() + 2 * tolerated_narrow_nblist_motion;
		Real const XH = etable_.max_non_hydrogen_lj_radius() + etable_.max_hydrogen_lj_radius()
			+ 2 * tolerated_narrow_nblist_motion;
		Real const HH = etable_.max_hydrogen_lj_radius() + etable_.max_hydrogen_lj_radius()
			+ 2 * tolerated_narrow_nblist_motion;

		count_pair::CountPairFunctionCOP count_pair = get_intrares_countpair( rsd, pose, sfxn );
		nbdata->initialize(
			rsd, count_pair, XX*XX, XH*XH, HH*HH );
		min_data.set_data( min_single_data_type(), nbdata );

	}
}


template < class Derived >
void
BaseEtableEnergy< Derived >::eval_intrares_derivatives(
	conformation::Residue const & rsd,
	ResSingleMinimizationData const & min_data,
	pose::Pose const & pose,
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & atom_derivs
) const
{
	if ( pose.energies().use_nblist_auto_update() ) return;

	debug_assert( utility::pointer::dynamic_pointer_cast< ResidueNblistData const > (min_data.get_data( min_single_data_type() ) ));
	ResidueNblistData const & nblist( static_cast< ResidueNblistData const & > ( min_data.get_data_ref( min_single_data_type() )) );

	utility::vector1< SmallAtNb > const & neighbs( nblist.atom_neighbors() );

	//prepare_for_residue_pair( 1, 1, pose ); // set intra-res
	typename Derived::Evaluator evaluator( static_cast< Derived const & > (*this).intrares_evaluator() );
	typename Derived::Evaluator nonprot_evaluator( static_cast< Derived const & > (*this).nonprot_intrares_evaluator() );
	evaluator.set_weights( weights );

	Vector f1(0.0),f2(0.0);
	for ( Size ii = 1, iiend = neighbs.size(); ii <= iiend; ++ii ) {
		conformation::Atom const & atom1( rsd.atom( neighbs[ ii ].atomno1() ) );
		conformation::Atom const & atom2( rsd.atom( neighbs[ ii ].atomno2() ) );
		Real dE_dR_over_r( evaluator.eval_dE_dR_over_r( atom1, atom2, weights, f1, f2 ) );
		if (!rsd.is_protein()) {
			dE_dR_over_r += nonprot_evaluator.eval_dE_dR_over_r( atom1, atom2, weights, f1, f2 );
		}
		if ( dE_dR_over_r != 0.0 ) {
			f1 *= dE_dR_over_r * neighbs[ ii ].weight();
			f2 *= dE_dR_over_r * neighbs[ ii ].weight();
			atom_derivs[ neighbs[ ii ].atomno1() ].f1() += f1;
			atom_derivs[ neighbs[ ii ].atomno1() ].f2() += f2;
			atom_derivs[ neighbs[ ii ].atomno2() ].f1() += -1 * f1;
			atom_derivs[ neighbs[ ii ].atomno2() ].f2() += -1 * f2;
		}
	}
}

/// @brief create a rotamer trie for a particular set, deciding upon the kind of count pair data that
/// needs to be contained by the trie.
///
template < class Derived >
EtableRotamerTrieOP
BaseEtableEnergy< Derived >::create_rotamer_trie(
	conformation::RotamerSetBase const & rotset,
	pose::Pose const & // will be need to create tries for disulfides
) const
{
	using namespace etrie;
	using namespace trie;

	CPDataCorrespondence cpdata_map( create_cpdata_correspondence_for_rotamerset( rotset ) );
	if ( cpdata_map.has_pseudobonds() ||
			cpdata_map.max_connpoints_for_residue() > 1 ||
			cpdata_map.n_entries() > 3 ) {
		EtableAtom at; CountPairDataGeneric cpdat;
		return create_trie( rotset, at, cpdat, cpdata_map, atomic_interaction_cutoff() );
	} else if ( cpdata_map.n_entries() == 1 || cpdata_map.n_entries() == 0 /* HACK! */ ) {
		EtableAtom at; CountPairData_1_1 cpdat;
		return create_trie( rotset, at, cpdat, cpdata_map, atomic_interaction_cutoff() );
	} else if ( cpdata_map.n_entries() == 2 ) {
		EtableAtom at; CountPairData_1_2 cpdat;
		return create_trie( rotset, at, cpdat, cpdata_map, atomic_interaction_cutoff() );
	} else if ( cpdata_map.n_entries() == 3 ) {
		EtableAtom at; CountPairData_1_3 cpdat;
		return create_trie( rotset, at, cpdat, cpdata_map, atomic_interaction_cutoff() );
	} else {
		/// As of 10/21, all count pair data combinations should be covered. This code should not execute.
		std::cerr << "Unsupported number of residue connections in trie construction." << std::endl;
		utility_exit();
		return 0;
	}
}

/// @details Create a one-residue trie.
template < class Derived >
EtableRotamerTrieOP
BaseEtableEnergy< Derived >::create_rotamer_trie(
	conformation::Residue const & residue,
	pose::Pose const & // will be need to create tries for disulfides
) const
{
	using namespace etrie;
	using namespace trie;

	CPDataCorrespondence cpdata_map( create_cpdata_correspondence_for_rotamer( residue ) );
	if ( cpdata_map.has_pseudobonds() ||
			cpdata_map.max_connpoints_for_residue() > 1 ||
			cpdata_map.n_entries() > 3 ) {
		EtableAtom at; CountPairDataGeneric cpdat;
		return create_trie( residue, at, cpdat, cpdata_map, atomic_interaction_cutoff() );
	} else if ( cpdata_map.n_entries() == 1 || cpdata_map.n_entries() == 0 /* HACK! */ ) {
		EtableAtom at; CountPairData_1_1 cpdat;
		return create_trie( residue, at, cpdat, cpdata_map, atomic_interaction_cutoff() );
	} else if ( cpdata_map.n_entries() == 2 ) {
		EtableAtom at; CountPairData_1_2 cpdat;
		return create_trie( residue, at, cpdat, cpdata_map, atomic_interaction_cutoff() );
	} else if ( cpdata_map.n_entries() == 3 ) {
		EtableAtom at; CountPairData_1_3 cpdat;
		return create_trie( residue, at, cpdat, cpdata_map, atomic_interaction_cutoff() );
	} else {
		/// As of 10/21, all count pair data combinations should be covered. This code should not execute.
		std::cerr << "Unsupported number of residue connections in trie construction." << std::endl;
		utility_exit();
		return 0;
	}
}

template < class Derived >
void
BaseEtableEnergy< Derived >::indicate_required_context_graphs( utility::vector1< bool > & /*context_graphs_required*/ ) const
{
	//context_graphs_required[ ten_A_neighbor_graph ] = true; // when did this get here?
}

template < class Derived >
void
BaseEtableEnergy< Derived >::bump_energy_full(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap & emap
) const
{
	EnergyMap tbemap;
	//prepare_for_residue_pair(rsd1.seqpos(),rsd2.seqpos(),pose);
	// count_pair::CountPairFunctionCOP cpfxn = get_count_pair_function( rsd1, rsd2, pose, sfxn );
	// cpfxn->residue_atom_pair_energy_sidechain_whole( rsd1, rsd2, static_cast<Derived const&> (*this), tbemap );
	if ( !calculate_interres( rsd1, rsd2 ) ) return;
	count_pair::CPCrossoverBehavior crossover = determine_crossover_behavior( rsd1, rsd2, pose, sfxn );
	typename Derived::Evaluator const & evaluator( static_cast< Derived const & > (*this).interres_evaluator() );
	SC_Whole_EnergyInvoker< typename Derived::Evaluator > invoker( rsd1, rsd2, evaluator, tbemap );
	count_pair::CountPairFactory::create_count_pair_function_and_invoke( rsd1, rsd2, crossover, invoker );

	emap[ evaluator.st_atr() ] += tbemap[ evaluator.st_atr() ]; // consider moving this into the derived class
	emap[ evaluator.st_rep() ] += tbemap[ evaluator.st_rep() ];

}

template < class Derived >
void
BaseEtableEnergy< Derived >::bump_energy_backbone(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap & emap
) const
{
	EnergyMap tbemap;
	//prepare_for_residue_pair(rsd1.seqpos(),rsd2.seqpos(),pose);
	//count_pair::CountPairFunctionCOP cpfxn = get_count_pair_function( rsd1, rsd2, pose, sfxn );
	//cpfxn->residue_atom_pair_energy_sidechain_backbone( rsd1, rsd2, static_cast<Derived const&> (*this), tbemap );
	if ( !calculate_interres( rsd1, rsd2 ) ) return;
	count_pair::CPCrossoverBehavior crossover = determine_crossover_behavior( rsd1, rsd2, pose, sfxn );
	typename Derived::Evaluator const & evaluator( static_cast< Derived const & > (*this).interres_evaluator() );
	SC_BB_EnergyInvoker< typename Derived::Evaluator > invoker( rsd1, rsd2, evaluator, tbemap );
	count_pair::CountPairFactory::create_count_pair_function_and_invoke( rsd1, rsd2, crossover, invoker );

	emap[ evaluator.st_atr() ] += tbemap[ evaluator.st_atr() ]; // consider moving this into the derived class
	emap[ evaluator.st_rep() ] += tbemap[ evaluator.st_rep() ];

}


///////////////////////////////////////////////////////////////////////////////

template < class Derived >
void
BaseEtableEnergy< Derived >::eval_atom_derivative(
	id::AtomID const & id,
	pose::Pose const & pose,
	kinematics::DomainMap const &, // domain_map,
	ScoreFunction const & /*sfxn*/, // needed for non-nblist minimization
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{
	Size const idresid = id.rsd();
	conformation::Atom const & atom1( pose.residue( idresid ).atom( id.atomno() ) );

	if ( pose.energies().use_nblist() ) {
		scoring::AtomNeighbors const & nbrs
			( pose.energies().nblist( nblist_type() ).atom_neighbors( id ) );

		typename Derived::Evaluator intrares_evaluator( static_cast< Derived const & > (*this).intrares_evaluator() );
		typename Derived::Evaluator nonprot_intrares_evaluator( static_cast< Derived const & > (*this).nonprot_intrares_evaluator() );
		typename Derived::Evaluator interres_evaluator( static_cast< Derived const & > (*this).interres_evaluator() );
		intrares_evaluator.set_weights( weights );
		interres_evaluator.set_weights( weights );
		Vector f1(0.0),f2(0.0);
		for ( scoring::AtomNeighbors::const_iterator it2=nbrs.begin(),
				it2e=nbrs.end(); it2 != it2e; ++it2 ) {
			scoring::AtomNeighbor const & nbr( *it2 );

			// Intra residue weights separate from interresidue weights; set which weights to use before scoring.
			// Comparison between idresid and nb.rsd performed by EtableEnergy but not by
			// CoarseEtableEnergy -> both are passed as arguments
			// static_cast<Derived const&> (*this).decide_scoretypes( idresid, nbr.rsd() );
			//prepare_for_residue_pair( idresid, nbr.rsd(), pose );
			if ( idresid == (Size) nbr.rsd() ) {
				Real const cp_weight( nbr.weight() );  // do not use nbr->weight_func() here
				conformation::Atom const & atom2( pose.residue( nbr.rsd() ).atom( nbr.atomno() ) );
				Real dE_dR_over_r( intrares_evaluator.eval_dE_dR_over_r( atom1, atom2, weights, f1, f2 ) );
				if (!pose.residue( idresid ).is_protein()) {
					dE_dR_over_r += nonprot_intrares_evaluator.eval_dE_dR_over_r( atom1, atom2, weights, f1, f2 );
				}
				//std::cout << "  gold atom deriv: " << idresid << " " << id.atomno() << " with " << nbr.rsd() << " " << nbr.atomno() << ". w= " << nbr.weight() << " dE_dR_over_r: " << dE_dR_over_r << std::endl ;
				if ( dE_dR_over_r != 0.0 ) {
					F1 += dE_dR_over_r * cp_weight * f1;
					F2 += dE_dR_over_r * cp_weight * f2;
				}
			} else {
				Real const cp_weight( nbr.weight() );  // do not use nbr->weight_func() here
				conformation::Atom const & atom2( pose.residue( nbr.rsd() ).atom( nbr.atomno() ) );
				Real const dE_dR_over_r( interres_evaluator.eval_dE_dR_over_r( atom1, atom2, weights, f1, f2 ) );
				//std::cout << "  gold atom deriv: " << idresid << " " << id.atomno() << " with " << nbr.rsd() << " " << nbr.atomno() << ". w= " << nbr.weight() << " dE_dR_over_r: " << dE_dR_over_r << std::endl ;
				if ( dE_dR_over_r != 0.0 ) {
					F1 += dE_dR_over_r * cp_weight * f1;
					F2 += dE_dR_over_r * cp_weight * f2;
				}
			}
		}
	} else {
		utility_exit_with_message("non-nblist minimize!");
	}
}


/// @brief return the Etables atomic distance cutoff
template < class Derived >
Distance
BaseEtableEnergy< Derived >::atomic_interaction_cutoff() const
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	return etable_.max_dis() +
		( option[ score::fa_Hatr ] ? chemical::MAX_CHEMICAL_BOND_TO_HYDROGEN_LENGTH * 2 : 0.0 ); /// HACK -- add in hydrogen/heavy max dist * 2
}


} // namespace etable
} // namespace scoring
} // namespace core

#endif
