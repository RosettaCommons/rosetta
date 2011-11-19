// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   protocols/features/ReportToDB.cc
///
/// @brief  report all data to a database
/// @author Matthew O'Meara


#ifdef USEMPI
#include <mpi.h>
#endif

#include <protocols/features/ReportToDB.hh>
#include <string>

// Setup Mover
#include <protocols/features/ReportToDBCreator.hh>

#include <basic/database/sql_utils.hh>

#include <platform/types.hh>
#include <core/svn_version.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/Adduct.fwd.hh>
#include <core/chemical/Adduct.hh>
#include <core/chemical/AtomICoor.fwd.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/AtomType.fwd.hh>
#include <core/chemical/AtomType.hh>
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
#include <core/chemical/VariantType.fwd.hh>
#include <core/chemical/types.hh>
#include <core/chemical/orbitals/ICoorOrbitalData.hh>
#include <core/chemical/orbitals/OrbitalType.fwd.hh>
#include <core/chemical/orbitals/OrbitalTypeMapper.fwd.hh>
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
#include <core/graph/ArrayPool.hh>
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
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/types.hh>
#include <core/io/pdb/pose_io.hh>
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
#include <core/pack/rotamer_set/RotamerCouplings.fwd.hh>
#include <core/pack/rotamer_set/RotamerSetOperation.fwd.hh>
#include <core/pack/task/IGEdgeReweightContainer.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/RotamerSampleOptions.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperation.fwd.hh>
#include <core/pose/MiniPose.fwd.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/PDBPoseMap.fwd.hh>
#include <core/pose/PDBPoseMap.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/pose/Remarks.fwd.hh>
#include <core/pose/Remarks.hh>
#include <core/pose/util.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/datacache/ObserverCache.fwd.hh>
#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
#include <core/pose/signals/ConformationEvent.fwd.hh>
#include <core/pose/signals/DestructionEvent.fwd.hh>
#include <core/pose/signals/EnergyEvent.fwd.hh>
#include <core/pose/signals/GeneralEvent.fwd.hh>
#include <core/scoring/AtomVDW.fwd.hh>
#include <core/scoring/CenHBPotential.fwd.hh>
#include <core/scoring/ContextGraph.fwd.hh>
#include <core/scoring/ContextGraph.hh>
#include <core/scoring/ContextGraphTypes.hh>
#include <core/scoring/DerivVectorPair.fwd.hh>
#include <core/scoring/DerivVectorPair.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergiesCacheableDataType.hh>
#include <core/scoring/EnergyGraph.fwd.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/EnvPairPotential.fwd.hh>
#include <core/scoring/GenBornPotential.fwd.hh>
#include <core/scoring/LREnergyContainer.fwd.hh>
#include <core/scoring/MembranePotential.fwd.hh>
#include <core/scoring/Membrane_FAPotential.fwd.hh>
#include <core/scoring/MinimizationData.fwd.hh>
#include <core/scoring/MinimizationGraph.fwd.hh>
#include <core/scoring/NeighborList.fwd.hh>
#include <core/scoring/OmegaTether.fwd.hh>
#include <core/scoring/P_AA.fwd.hh>
#include <core/scoring/PairEPotential.fwd.hh>
#include <core/scoring/PoissonBoltzmannPotential.fwd.hh>
#include <core/scoring/Ramachandran.fwd.hh>
#include <core/scoring/Ramachandran2B.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunctionInfo.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <core/scoring/ScoringManager.fwd.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/SecondaryStructurePotential.fwd.hh>
#include <core/scoring/SmoothEnvPairPotential.fwd.hh>
#include <core/scoring/TenANeighborGraph.fwd.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/TwelveANeighborGraph.fwd.hh>
#include <core/scoring/TwelveANeighborGraph.hh>
#include <core/scoring/UnfoldedStatePotential.fwd.hh>
#include <core/scoring/WaterAdductHBondPotential.fwd.hh>
#include <core/scoring/sasa.hh>
#include <core/scoring/carbon_hbonds/CarbonHBondPotential.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/disulfides/CentroidDisulfidePotential.fwd.hh>
#include <core/scoring/disulfides/DisulfideMatchingPotential.fwd.hh>
#include <core/scoring/disulfides/FullatomDisulfidePotential.fwd.hh>
#include <core/scoring/dna/DNA_BasePotential.fwd.hh>
#include <core/scoring/dna/DirectReadoutPotential.fwd.hh>
#include <core/scoring/etable/BaseEtableEnergy.hh>
#include <core/scoring/etable/Etable.fwd.hh>
#include <core/scoring/etable/EtableEnergy.fwd.hh>
#include <core/scoring/etable/EtableEnergy.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.fwd.hh>
#include <core/scoring/etable/count_pair/types.hh>
#include <core/scoring/etable/etrie/EtableAtom.fwd.hh>
#include <core/scoring/etable/etrie/EtableTrie.fwd.hh>
#include <core/scoring/geometric_solvation/DatabaseOccSolEne.fwd.hh>
#include <core/scoring/geometric_solvation/ExactOccludedHbondSolEnergy.hh>
#include <core/scoring/hbonds/HBondDatabase.fwd.hh>
#include <core/scoring/hbonds/HBondOptions.fwd.hh>
#include <core/scoring/hbonds/types.hh>
#include <core/scoring/interface/DDPlookup.fwd.hh>
#include <core/scoring/memb_etable/MembEtable.fwd.hh>
#include <core/scoring/methods/ContextDependentLRTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextDependentOneBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextDependentOneBodyEnergy.hh>
#include <core/scoring/methods/ContextDependentTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentLRTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>
#include <core/scoring/methods/EnergyMethod.fwd.hh>
#include <core/scoring/methods/EnergyMethod.hh>
#include <core/scoring/methods/EnergyMethodCreator.fwd.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/methods/LongRangeTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/Methods.hh>
#include <core/scoring/methods/OneBodyEnergy.fwd.hh>
#include <core/scoring/methods/OneBodyEnergy.hh>
#include <core/scoring/methods/ShortRangeTwoBodyEnergy.hh>
#include <core/scoring/methods/TwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/TwoBodyEnergy.hh>
#include <core/scoring/methods/WholeStructureEnergy.fwd.hh>
#include <core/scoring/mm/MMBondAngleLibrary.fwd.hh>
#include <core/scoring/mm/MMBondLengthLibrary.fwd.hh>
#include <core/scoring/mm/MMLJEnergyTable.fwd.hh>
#include <core/scoring/mm/MMLJLibrary.fwd.hh>
#include <core/scoring/mm/MMTorsionLibrary.fwd.hh>
#include <core/scoring/nv/NVlookup.fwd.hh>
#include <core/scoring/nv/NVlookup.hh>
#include <core/scoring/nv/NVscore.fwd.hh>
#include <core/scoring/nv/NVscore.hh>
#include <core/scoring/orbitals/OrbitalsLookup.hh>
#include <core/scoring/rna/RNA_AtomVDW.fwd.hh>
#include <core/scoring/rna/RNA_LowResolutionPotential.fwd.hh>
#include <core/scoring/rna/RNA_TorsionPotential.fwd.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.fwd.hh>
#include <core/scoring/trie/RotamerTrieBase.fwd.hh>
#include <core/scoring/trie/TrieCountPairBase.fwd.hh>
#include <protocols/features/FeaturesReporter.fwd.hh>
#include <protocols/features/FeaturesReporter.hh>
#include <protocols/features/FeaturesReporterCreator.fwd.hh>
#include <protocols/features/FeaturesReporterFactory.fwd.hh>
#include <protocols/features/FeaturesReporterFactory.hh>
#include <protocols/features/ProteinRMSDFeatures.fwd.hh>
#include <protocols/features/ProteinRMSDFeatures.hh>
#include <protocols/features/ProtocolFeatures.fwd.hh>
#include <protocols/features/ProtocolFeatures.hh>
#include <protocols/features/StructureFeatures.fwd.hh>
#include <protocols/features/StructureFeatures.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/jd2/InnerJob.fwd.hh>
#include <protocols/jd2/Job.fwd.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.fwd.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobInputter.fwd.hh>
#include <protocols/jd2/JobOutputter.fwd.hh>
#include <protocols/jd2/Parser.fwd.hh>
#include <protocols/jobdist/Jobs.fwd.hh>
#include <protocols/jobdist/Jobs.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/dssp/StrandPairing.fwd.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/moves/DataMap.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverCreator.hh>
#include <protocols/moves/MoverStatus.hh>
#include <protocols/features/ReportToDB.fwd.hh>
#include <protocols/rosetta_scripts/util.hh>
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
#include <utility/factory/WidgetRegistrator.hh>
#include <utility/file/FileName.fwd.hh>
#include <utility/file/FileName.hh>
#include <utility/file/PathName.fwd.hh>
#include <utility/file/PathName.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/file/gzip_util.hh>
#include <utility/io/irstream.fwd.hh>
#include <utility/io/irstream.hh>
#include <utility/io/izstream.fwd.hh>
#include <utility/io/izstream.hh>
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
#include <utility/sql_database/DatabaseSessionManager.fwd.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/Tag.hh>
#include <numeric/MathMatrix.hh>
#include <numeric/MathVector.hh>
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
#include <numeric/interpolation/spline/Bicubic_spline.hh>
#include <numeric/interpolation/spline/Cubic_spline.fwd.hh>
#include <numeric/interpolation/spline/Cubic_spline.hh>
#include <numeric/interpolation/spline/Interpolator.hh>
#include <numeric/interpolation/spline/SplineGenerator.hh>
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
#include <ObjexxFCL/ubyte.fwd.hh>
#include <ObjexxFCL/ubyte.hh>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iosfwd>
#include <iostream>
#include <istream>
#include <limits>
#include <list>
#include <map>
#include <math.h>
#include <numeric>
#include <ostream>
#include <set>
#include <sstream>
#include <typeinfo>
#include <utility>
#include <vector>
#include <basic/MetricValue.fwd.hh>
#include <basic/Tracer.fwd.hh>
#include <basic/Tracer.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableData.fwd.hh>
#include <basic/datacache/CacheableData.hh>
#include <basic/datacache/CacheableString.fwd.hh>
#include <basic/datacache/CacheableString.hh>
#include <basic/datacache/DataCache.fwd.hh>
#include <basic/datacache/DataCache.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/features.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/parser.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <boost/algorithm/string/erase.hpp>
#include <boost/bind.hpp>
#include <boost/config.hpp>
#include <boost/foreach.hpp>
#include <boost/function.hpp>
#include <boost/functional/hash.hpp>
#include <boost/pool/detail/mutex.hpp>
#include <boost/pool/poolfwd.hpp>
#include <boost/scoped_ptr.hpp>
#include <cppdb/errors.h>
#include <cppdb/frontend.h>
#include <zlib/zlib.h>
#include <zlib/zutil.h>



namespace protocols{
namespace features{

std::string
ReportToDBCreator::keyname() const
{
  return ReportToDBCreator::mover_name();
}

moves::MoverOP
ReportToDBCreator::create_mover() const {
  return new ReportToDB;
}

std::string
ReportToDBCreator::mover_name()
{
  return "ReportToDB";
}

/// Macros are not properly caught and passed along by my #inclusion
/// cleanup script
#define foreach BOOST_FOREACH

using basic::T;
using basic::Tracer;
using basic::Error;
using basic::Warning;
using basic::datacache::CacheableString;
using core::Size;
using core::pack::task::PackerTaskCOP;
using core::pack::task::TaskFactory;
using core::pose::Pose;
using core::pose::PoseOP;
using core::scoring::getScoreFunction;
using core::scoring::ScoreFunctionFactory;
using core::scoring::ScoreFunction;
using core::scoring::ScoreFunctionOP;
using core::scoring::ScoreFunctionCOP;
using core::scoring::STANDARD_WTS;
using cppdb::cppdb_error;
using cppdb::statement;
using cppdb::result;
using protocols::features::FeaturesReporterOP;
using protocols::features::ProteinRMSDFeatures;
using protocols::features::ProtocolFeatures;
using protocols::features::StructureFeatures;
using protocols::features::FeaturesReporterFactory;
using protocols::jd2::JobDistributor;
using protocols::moves::MoverOP;
using protocols::moves::DataMap;
using protocols::moves::Movers_map;
using protocols::rosetta_scripts::parse_task_operations;
using std::string;
using std::endl;
using std::accumulate;
using std::stringstream;
using utility::file::FileName;
using utility::vector0;
using utility::vector1;
using utility::tag::TagPtr;
using utility::sql_database::DatabaseSessionManager;
using utility::sql_database::session;

static Tracer TR("protocols.features.ReportToDB");

Size ReportToDB::struct_id_ = 0;
Size ReportToDB::protocol_id_ = 0;
bool ReportToDB::autoincrement_struct_id_ = true;
bool ReportToDB::protocol_table_initialized_ = false;


ReportToDB::ReportToDB():
	Mover("ReportToDB"),
	database_fname_("FeatureStatistics.db3"),
	database_mode_("sqlite3"),
	sample_source_("Rosetta: Unknown Protocol"),
	scfxn_(getScoreFunction()),
	use_transactions_(true),
	cache_size_(2000),
	task_factory_(new TaskFactory()),
	features_reporter_factory_(FeaturesReporterFactory::get_instance()),
	features_reporters_(),
	initialized( false )
{
	initialize_reporters();
}

ReportToDB::ReportToDB(string const & name):
	Mover(name),
	database_fname_("FeatureStatistics.db3"),
	database_mode_("sqlite3"),
	sample_source_("Rosetta: Unknown Protocol"),
	scfxn_( ScoreFunctionFactory::create_score_function( STANDARD_WTS ) ),
	cache_size_(2000),
	use_transactions_(true),
	task_factory_(new TaskFactory()),
	features_reporter_factory_(FeaturesReporterFactory::get_instance()),
	features_reporters_(),
	initialized( false )
{
	initialize_reporters();
}

ReportToDB::ReportToDB(
	string const & name,
	string const & database_fname,
	string const & sample_source,
	ScoreFunctionOP scfxn,
	bool use_transactions,
	Size cache_size) :
	Mover(name),
	database_fname_(database_fname),
	database_mode_("sqlite3"),
	sample_source_(sample_source),
	scfxn_(scfxn),
	use_transactions_(use_transactions),
	cache_size_(cache_size),
	task_factory_(new TaskFactory()),
	features_reporter_factory_(FeaturesReporterFactory::get_instance()),
	features_reporters_(),
	initialized( false )
{
	initialize_reporters();
}

ReportToDB::ReportToDB( ReportToDB const & src):
	Mover(src),
	database_fname_(src.database_fname_),
	database_mode_(src.database_mode_),
	sample_source_(src.sample_source_),
	scfxn_(new ScoreFunction(* src.scfxn_)),
	use_transactions_(src.use_transactions_),
	cache_size_(src.cache_size_),
	task_factory_(src.task_factory_),
	features_reporter_factory_(FeaturesReporterFactory::get_instance()),
	protocol_features_(src.protocol_features_),
	structure_features_(src.structure_features_),
	features_reporters_(src.features_reporters_),
	initialized(src.initialized)
{}

ReportToDB::~ReportToDB(){}

void
ReportToDB::register_options() const {
	using basic::options::option;
	using namespace basic::options::OptionKeys;

	// This mover is equiped to work with the Rosetta Scripts interface
	option.add_relevant( parser::protocol );

	//TODO call relevant_options on FeaturesMover objects
}

MoverOP ReportToDB::fresh_instance() const { return new ReportToDB; }

MoverOP ReportToDB::clone() const
{
	return new ReportToDB( *this );
}

void
ReportToDB::parse_db_tag_item(
	TagPtr const tag){

	if( tag->hasOption("db") ){
		database_fname_ = tag->getOption<string>("db");
	} else {
		TR << "Field 'db' required for use of ReportToDB mover in Rosetta Scripts." << endl;
	}
}

void
ReportToDB::parse_sample_source_tag_item(
	TagPtr const tag){
	if( tag->hasOption("sample_source") ){
		sample_source_ = tag->getOption<string>("sample_source");
	} else {
		TR << "Field 'sample_source' required for use of ReportToDB in Rosetta Scripts." << endl;
		TR << "The sample_sourOBce should describe where the samples came from. To access the description run \"sqlite3 'select description from sample_source' fname.db3\"" << endl;
		TR << "For example: Top4400 natives from Richardson Lab. Reduce placed hydrogens with -correct flag." << endl;
	}
}

void
ReportToDB::parse_protocol_id_tag_item(
	TagPtr const tag){

	if(tag->hasOption("protocol_id")){
#ifdef USEMPI
		int mpi_rank(0);
		MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
		protocol_id_ = tag->getOption<Size>("protocol_id") + mpi_rank;
#else
		protocol_id_ = tag->getOption<Size>("protocol_id");
#endif
	}
#ifdef USEMPI
	else {
		int mpi_rank(0);
		MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
		protocol_id_ = mpi_rank;
	}
#endif

}

void
ReportToDB::parse_first_struct_id_tag_item(
	TagPtr const tag) {

	if(tag->hasOption("first_struct_id")){
		Size first_struct_id = tag->getOption<Size>("first_struct_id");
		// Initialize struct_id_ to first_struct_id only if this is the
		// first time ReportToDB has been executed. It is sufficient to
		// test autoincrement_struct_id_.
		if (autoincrement_struct_id_){
#ifdef USEMPI
			int mpi_rank(0);
			MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
			struct_id_ = first_struct_id + mpi_rank - 1; // minus 1 to setup invariant
#else
			struct_id_ = first_struct_id - 1; // minus 1 to setup invariant
#endif
			autoincrement_struct_id_ = false;
		}
	}
#ifdef USEMPI
	else if (autoincrement_struct_id_){
		int mpi_rank(0);
		MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
		struct_id_ = mpi_rank - 1;
		autoincrement_struct_id_ = false;
	}
#endif

}

void
ReportToDB::parse_db_mode_tag_item(
	TagPtr const tag) {
	if(tag->hasOption("db_mode")){
		database_mode_ = tag->getOption<std::string>("db_mode");
	} else {
		database_mode_ = "sqlite3";
	}
}

void
ReportToDB::parse_use_transactions_tag_item(
	TagPtr const tag) {
	if(tag->hasOption("use_transactions")){
		use_transactions_ = tag->getOption<bool>("use_transactions");
	}
}

void
ReportToDB::parse_cache_size_tag_item(
	TagPtr const tag) {
	if(tag->hasOption("cache_size")){
		cache_size_ = tag->getOption<bool>("cache_size");
	}
}


/// Allow ReportToDB to be called from RosettaScripts
/// See
void
ReportToDB::parse_my_tag(
	TagPtr const tag,
	DataMap & data,
	Filters_map const & filters,
	Movers_map const & movers,
	Pose const & pose )
{

	// Name of output features database:
	// EXAMPLE: db=features_<sample_source>.db3
	// REQUIRED
	parse_db_tag_item(tag);

	// Description of features database
	// EXAMPLE: sample_source="This is a description of the sample source."
	// RECOMMENDED
	parse_sample_source_tag_item(tag);

	// Manually control the id of associated with this protocol
	// EXAMPLE: protocol_id=6
	// OPTIONAL
	parse_protocol_id_tag_item(tag);

	// Manually control the first structure id
	// EXAMPLE: first_struct_id=100
	// OPTIONAL
	parse_first_struct_id_tag_item(tag);

	// The database backend to use ('sqlite3', 'mysql' etc)
	// EXAMPLE: db_mode=sqlite3
	parse_db_mode_tag_item(tag);

	// Use transactions to group database i/o to be more efficient. Turning them off OBcan help debugging.
	// EXAMPLE: use_transactions=true
	// DEFAULT: TRUE
	parse_use_transactions_tag_item(tag);

	// Specify the maximum number 1k pages to keep in memory before writing to disk
	// EXAMPLE: cache_size=1000000  // this uses ~ 1GB of memory
	// DEFAULT: 2000
	parse_cache_size_tag_item(tag);

	task_factory_ = parse_task_operations(tag, data);

	vector0< TagPtr >::const_iterator begin=tag->getTags().begin();
	vector0< TagPtr >::const_iterator end=tag->getTags().end();

	for(; begin != end; ++begin){
		TagPtr feature_tag= *begin;
		//	foreach(TagPtr const & feature_tag, tag->getTags()){

		if(feature_tag->getName() != "feature"){
			TR.Error << "Please include only tags with name 'feature' as subtags of ReportToDB" << endl;
			TR.Error << "Tag with name '" << feature_tag->getName() << "' is invalid" << endl;
			utility_exit();
		}

		features_reporters_.push_back(
			features_reporter_factory_->get_features_reporter(
				feature_tag, data, filters, movers, pose));
	}

}

void
ReportToDB::initialize_reporters()
{
	// the protocol and structure features are special
	protocol_features_ = new ProtocolFeatures();
	structure_features_ = new StructureFeatures();
}

void
ReportToDB::initialize_database(
	utility::sql_database::sessionOP db_session
){
	if (!initialized){
		protocol_features_->write_schema_to_db(db_session);
		structure_features_->write_schema_to_db(db_session);
		foreach( FeaturesReporterOP const & reporter, features_reporters_ ){
			reporter->write_schema_to_db(db_session);
		}
		initialized = true;
	}

}

void
ReportToDB::apply( Pose& pose ){

	utility::sql_database::sessionOP db_session(basic::database::get_db_session(database_fname_, database_mode_, false, true));

	// Make sure energy objects are initialized
	(*scfxn_)(pose);

	PackerTaskCOP task(task_factory_->create_task_and_apply_taskoperations(pose));
	vector1< bool > relevant_residues(task->repacking_residues());

	TR << "Reporting features for "
		<< accumulate(relevant_residues.begin(), relevant_residues.end(), 0)
		<< " of the " << pose.total_residue()
		<< " total residues in the pose." << endl;

	if(use_transactions_) db_session->begin();
	initialize_database(db_session);
	if(use_transactions_) db_session->commit();

	if(use_transactions_) db_session->begin();

	stringstream stmt_ss; stmt_ss << "PRAGMA cache_size = " << cache_size_ << ";";
	statement stmt = (*db_session) << stmt_ss.str(); stmt.exec();

	if(!protocol_table_initialized_){
		protocol_id_ = protocol_features_->report_features(
			protocol_id_, db_session);
		protocol_table_initialized_ = true;

		statement stmt = (*db_session) <<
			"CREATE TABLE IF NOT EXISTS sample_source ( description TEXT );";
		stmt.exec();
		stmt = (*db_session)
			<< "INSERT INTO sample_source VALUES (?);"
			<< sample_source_;
		stmt.exec();
	}

	if(autoincrement_struct_id_){
		struct_id_ = structure_features_->report_features(
		  pose, relevant_residues, protocol_id_, db_session);
	} else {
#ifdef USEMPI
		int mpi_size(0);
		MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
		struct_id_+=mpi_size;
#else
		struct_id_++;
#endif
		structure_features_->report_features(
			pose, relevant_residues, struct_id_, protocol_id_, db_session);
	}

	for(Size i=1; i <= features_reporters_.size(); ++i){
		features_reporters_[i]->report_features(
			pose, relevant_residues, struct_id_, db_session);
	}

	if(use_transactions_) db_session->commit();
}

} // namespace
} // namespace
