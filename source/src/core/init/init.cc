// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/init/init.cc
/// @brief  options system initialization routines
/// @author Sergey Lyskov


#ifdef USEMPI
#include <mpi.h> // Must go first
#include <basic/TracerToFile.hh>
#endif

// Unit headers
#include <core/init/init.hh>

// Project Headers
#include <core/types.hh>
#include <core/id/bogus.hh>
#include <basic/options/option.hh>
#include <utility/sys_util.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/io/izstream.hh>
#include <basic/Tracer.hh>
#include <basic/prof.hh>

#include <core/chemical/bcl/util.hh>
#include <core/chemical/rdkit/util.hh>

#include <basic/random/init_random_generator.hh>
#include <numeric/random/random.hh>

// Factories (for registration)
#include <core/chemical/rotamers/RotamerLibrarySpecificationFactory.hh>
#include <core/pack/palette/PackerPaletteFactory.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibraryFactory.hh>
#include <core/select/jump_selector/JumpSelectorFactory.hh>
#include <core/select/residue_selector/ResidueSelectorFactory.hh>
#include <basic/resource_manager/ResourceLoaderFactory.hh>

// Classes in core that must register with factories
#include <core/init/score_function_corrections.hh>
#include <core/energy_methods/AACompositionEnergyCreator.hh>
#include <core/energy_methods/AARepeatEnergyCreator.hh>
#include <core/pack/guidance_scoreterms/approximate_buried_unsat_penalty/ApproximateBuriedUnsatPenaltyCreator.hh>
#include <core/pack/guidance_scoreterms/buried_unsat_penalty/BuriedUnsatPenaltyCreator.hh>
#include <core/scoring/constraints/ConstraintsEnergyCreator.hh>
#include <core/energy_methods/CentroidDisulfideEnergyCreator.hh>
#include <core/energy_methods/DisulfideMatchingEnergyCreator.hh>
#include <core/energy_methods/FullatomDisulfideEnergyCreator.hh>
#include <core/scoring/etable/EtableEnergyCreator.hh>
#include <core/pack/guidance_scoreterms/hbnet_energy/HBNetEnergyCreator.hh>
#include <core/scoring/hbonds/HBondEnergyCreator.hh>
#include <core/energy_methods/NPDHBondEnergyCreator.hh>
#include <core/scoring/methods/EnergyMethodRegistrator.hh>
#include <core/energy_methods/ArgCationPiEnergyCreator.hh>
#include <core/energy_methods/AspartimidePenaltyEnergyCreator.hh>
#include <core/energy_methods/CenPairEnergyCreator.hh>
#include <core/energy_methods/ChainbreakEnergyCreator.hh>
#include <core/energy_methods/BranchEnergyCreator.hh>
#include <core/energy_methods/ContactOrderEnergyCreator.hh>
#include <core/energy_methods/EnvEnergyCreator.hh>
#include <core/energy_methods/EnvSmoothEnergyCreator.hh>
#include <core/energy_methods/IntermolEnergyCreator.hh>
#include <core/energy_methods/MissingEnergyCreator.hh>
#include <core/scoring/nmer/NMerRefEnergyCreator.hh>
#include <core/scoring/nmer/NMerPSSMEnergyCreator.hh>
#include <core/scoring/nmer/NMerSVMEnergyCreator.hh>
#include <core/energy_methods/OmegaTetherEnergyCreator.hh>
#include <core/energy_methods/OtherPoseEnergyCreator.hh>
#include <core/energy_methods/P_AA_EnergyCreator.hh>
#include <core/energy_methods/P_AA_ss_EnergyCreator.hh>
#include <core/energy_methods/P_AA_pp_EnergyCreator.hh>
#include <core/energy_methods/AbegoEnergyCreator.hh>
#include <core/energy_methods/PackStatEnergyCreator.hh>
#include <core/energy_methods/PairEnergyCreator.hh>
#include <core/energy_methods/PeptideBondEnergyCreator.hh>
#include <core/energy_methods/ProClosureEnergyCreator.hh>
#include <core/energy_methods/BurialEnergyCreator.hh>
#include <core/energy_methods/Burial_v2EnergyCreator.hh>
#include <core/energy_methods/HRF_MSLabelingEnergyCreator.hh>
#include <core/energy_methods/HRFDynamicsEnergyCreator.hh>
#include <core/energy_methods/CCS_IMMSEnergyCreator.hh>
#include <core/energy_methods/CovalentLabelingEnergyCreator.hh>
#include <core/energy_methods/CovalentLabelingFAEnergyCreator.hh>
#include <core/energy_methods/DEPC_MS_EnergyCreator.hh>
#include <core/energy_methods/RG_Energy_FastCreator.hh>
#include <core/energy_methods/RG_LocalEnergyCreator.hh>
#include <core/energy_methods/SA_EnergyCreator.hh>
#include <core/energy_methods/SSElementMotifContactEnergyCreator.hh>
#include <core/energy_methods/D2H_SA_EnergyCreator.hh>
#include <core/energy_methods/ProQ_EnergyCreator.hh>
#include <core/energy_methods/RMS_EnergyCreator.hh>
#include <core/energy_methods/RamaPreProEnergyCreator.hh>
#include <core/energy_methods/RamachandranEnergy2BCreator.hh>
#include <core/energy_methods/RamachandranEnergyCreator.hh>
#include <core/energy_methods/ReferenceEnergyCreator.hh>
#include <core/energy_methods/ReferenceEnergyNoncanonicalCreator.hh>
#include <core/energy_methods/SecondaryStructureEnergyCreator.hh>
#include <core/energy_methods/SugarBackboneEnergyCreator.hh>
#include <core/energy_methods/DFIRE_EnergyCreator.hh>
#include <core/pack/guidance_scoreterms/sap/SapConstraintEnergyCreator.hh>
#include <core/scoring/vdwaals/VDW_EnergyCreator.hh>
#include <core/energy_methods/GoapEnergyCreator.hh>
#include <core/energy_methods/RingClosureEnergyCreator.hh>
#include <core/energy_methods/AromaticBackboneRestraintEnergyCreator.hh>
#include <core/energy_methods/MHCEpitopeEnergyCreator.hh>
#include <core/energy_methods/DumpTrajectoryEnergyCreator.hh>
#include <core/energy_methods/NetChargeEnergyCreator.hh>

#include <core/pack/dunbrack/DunbrackEnergyCreator.hh>
#include <core/pack/dunbrack/cenrot/CenRotDunEnergyCreator.hh>
#include <core/pack/guidance_scoreterms/voids_penalty_energy/VoidsPenaltyEnergyCreator.hh>

// define this for compiling a slimmed down version of mini libraries lacking about 3/4s of the code
// this is required for compiling a less memory hungry version of mini for Bluegene etc..
//#define MINI_SLIM
#ifndef MINI_SLIM
#include <core/scoring/carbon_hbonds/CarbonHBondEnergyCreator.hh>
#include <core/scoring/custom_pair_distance/FullatomCustomPairDistanceEnergyCreator.hh>
#include <core/energy_methods/CustomAtomPairEnergyCreator.hh>
#include <core/energy_methods/FastDensEnergyCreator.hh>
#include <core/energy_methods/ElecDensCenEnergyCreator.hh>
#include <core/energy_methods/ElecDensAllAtomCenEnergyCreator.hh>
#include <core/energy_methods/ElecDensEnergyCreator.hh>
#include <core/energy_methods/XtalMLEnergyCreator.hh>
#include <core/energy_methods/ElecDensAtomwiseEnergyCreator.hh>
//XRW_B_T1
//#include <core/scoring/etable/CoarseEtableEnergyCreator.hh>
//XRW_E_T1
#include <core/energy_methods/ExactOccludedHbondSolEnergyCreator.hh>
#include <core/energy_methods/ContextDependentGeometricSolEnergyCreator.hh>
#include <core/energy_methods/ContextIndependentGeometricSolEnergyCreator.hh>
#include <core/energy_methods/OccludedHbondSolEnergyCreator.hh>
#include <core/energy_methods/OccludedHbondSolEnergy_onebodyCreator.hh>
#include <core/energy_methods/DEEREnergyCreator.hh>
#include <core/scoring/elec/FA_ElecEnergyCreator.hh>
#include <core/energy_methods/ImplicitMembraneElecEnergyCreator.hh>
#include <core/energy_methods/FA_GrpElecEnergyCreator.hh>
#include <core/energy_methods/FA_ElecEnergyAroAroCreator.hh>
#include <core/energy_methods/FA_ElecEnergyAroAllCreator.hh>
#include <core/energy_methods/RNA_FA_ElecEnergyCreator.hh>
#include <core/energy_methods/HackAroEnergyCreator.hh>
#include <core/energy_methods/DNAChiEnergyCreator.hh>
#include <core/energy_methods/DNATorsionEnergyCreator.hh>
#include <core/energy_methods/DNA_BaseEnergyCreator.hh>
#include <core/energy_methods/DNA_DihedralEnergyCreator.hh>
#include <core/energy_methods/DirectReadoutEnergyCreator.hh>
#include <core/energy_methods/DistanceChainbreakEnergyCreator.hh>
#include <core/energy_methods/Fa_MbenvEnergyCreator.hh>
#include <core/energy_methods/Fa_MbsolvEnergyCreator.hh>
#include <core/energy_methods/FreeDOF_EnergyCreator.hh>
#include <core/energy_methods/HybridVDW_EnergyCreator.hh>
#include <core/energy_methods/GenBornEnergyCreator.hh>
#include <core/energy_methods/GenericBondedEnergyCreator.hh>
#include <core/energy_methods/VdWTinkerEnergyCreator.hh>
#include <core/energy_methods/MultipoleElecEnergyCreator.hh>
#include <core/energy_methods/SASAEnergyCreator.hh>
#include <core/energy_methods/FACTSEnergyCreator.hh>
#include <core/energy_methods/LK_PolarNonPolarEnergyCreator.hh>
#include <core/energy_methods/LK_hackCreator.hh>
#include <core/scoring/lkball/LK_BallEnergyCreator.hh>
#include <core/energy_methods/LinearChainbreakEnergyCreator.hh>
#include <core/energy_methods/LinearBranchEnergyCreator.hh>
#include <core/scoring/lkball/LK_DomeEnergyCreator.hh>
#include <core/pack/guidance_scoreterms/lk_dome/LK_DomePackEnergyCreator.hh>
#include <core/scoring/methods/MMBondAngleEnergyCreator.hh>
#include <core/scoring/methods/MMBondLengthEnergyCreator.hh>
#include <core/energy_methods/CartesianBondedEnergyCreator.hh>
#include <core/scoring/methods/MMTorsionEnergyCreator.hh>
#include <core/scoring/methods/MMLJEnergyIntraCreator.hh>
#include <core/scoring/methods/MMLJEnergyInterCreator.hh>
#include <core/energy_methods/MembraneCbetaEnergyCreator.hh>
#include <core/energy_methods/MembraneCenPairEnergyCreator.hh>
#include <core/energy_methods/MembraneEnvEnergyCreator.hh>
#include <core/energy_methods/MembraneEnvPenaltiesCreator.hh>
#include <core/energy_methods/MembraneLipoCreator.hh>
#include <core/energy_methods/MembraneEnvSmoothEnergyCreator.hh>
#include <core/energy_methods/MPPairEnergyCreator.hh>
#include <core/energy_methods/MPEnvEnergyCreator.hh>
#include <core/energy_methods/MPCBetaEnergyCreator.hh>
#include <core/scoring/membrane/MPNonHelixPenaltyCreator.hh>
#include <core/scoring/membrane/MPTerminiPenaltyCreator.hh>
#include <core/scoring/membrane/MPTMProjPenaltyCreator.hh>
#include <core/energy_methods/FaMPEnvEnergyCreator.hh>
#include <core/energy_methods/FaMPSolvEnergyCreator.hh>
#include <core/energy_methods/FaMPEnvSmoothEnergyCreator.hh>
#include <core/energy_methods/FaMPAsymEzCBEnergyCreator.hh>
#include <core/energy_methods/FaMPAsymEzCGEnergyCreator.hh>
#include <core/energy_methods/MPResidueLipophilicityEnergyCreator.hh>
#include <core/energy_methods/MPHelicalityEnergyCreator.hh>
#include <core/energy_methods/MPSpanInsertionEnergyCreator.hh>
#include <core/energy_methods/MPSpanAngleEnergyCreator.hh>
#include <core/energy_methods/pHEnergyCreator.hh>
#include <core/energy_methods/PoissonBoltzmannEnergyCreator.hh>
#include <core/energy_methods/ChemicalShiftAnisotropyEnergyCreator.hh>
#include <core/energy_methods/DipolarCouplingEnergyCreator.hh>
#include <core/energy_methods/ResidualDipolarCouplingEnergyCreator.hh>
#include <core/energy_methods/ResidualDipolarCouplingEnergy_RohlCreator.hh>
#include <core/energy_methods/SmoothCenPairEnergyCreator.hh>
#include <core/energy_methods/SmoothEnvEnergyCreator.hh>
#include <core/energy_methods/CenPairMotifEnergyCreator.hh>
#include <core/energy_methods/CenPairMotifDegreeEnergyCreator.hh>
#include <core/energy_methods/CenRotPairEnergyCreator.hh>
#include <core/energy_methods/CenRotEnvEnergyCreator.hh>
#include <core/energy_methods/CenHBEnergyCreator.hh>
#include <core/energy_methods/MotifDockEnergyCreator.hh>
#include <core/energy_methods/SuckerEnergyCreator.hh>
#include <core/energy_methods/GaussianOverlapEnergyCreator.hh>
#include <core/energy_methods/YHHPlanarityEnergyCreator.hh>
#include <core/energy_methods/HydroxylTorsionEnergyCreator.hh>
#include <core/pack/interaction_graph/SurfaceEnergyCreator.hh>
#include <core/pack/interaction_graph/HPatchEnergyCreator.hh>
#include <core/energy_methods/MgEnergyCreator.hh>
#include <core/energy_methods/RNA_MgPointEnergyCreator.hh>
#include <core/energy_methods/SymmetricLigandEnergyCreator.hh>
#include <core/energy_methods/UnfoldedStateEnergyCreator.hh>
#include <core/energy_methods/SplitUnfoldedTwoBodyEnergyCreator.hh>
#include <core/energy_methods/WaterAdductHBondEnergyCreator.hh>
#include <core/energy_methods/WaterAdductIntraEnergyCreator.hh>
#include <core/energy_methods/WaterSpecificEnergyCreator.hh>
#include <core/scoring/nv/NVscoreCreator.hh>
#include <core/scoring/orbitals/OrbitalsScoreCreator.hh>
#include <core/scoring/interface_/DDPscoreCreator.hh>
#include <core/energy_methods/HolesEnergyCreator.hh>
#include <core/energy_methods/SurfVolEnergyCreator.hh>
#include <core/energy_methods/SurfEnergyCreator.hh>
#include <core/energy_methods/RG_Energy_RNACreator.hh>
#include <core/energy_methods/RNA_BulgeEnergyCreator.hh>
#include <core/energy_methods/RNA_CoarseDistEnergyCreator.hh>
#include <core/energy_methods/RNA_FullAtomStackingEnergyCreator.hh>
#include <core/energy_methods/RNA_JR_SuiteEnergyCreator.hh>
#include <core/energy_methods/RNA_LJ_BaseEnergyCreator.hh>
#include <core/energy_methods/RNA_PairwiseLowResolutionEnergyCreator.hh>
#include <core/energy_methods/RNA_PartitionEnergyCreator.hh>
#include <core/energy_methods/RNP_LowResEnergyCreator.hh>
#include <core/energy_methods/RNP_LowResPairDistEnergyCreator.hh>
#include <core/energy_methods/RNP_LowResStackEnergyCreator.hh>
#include <core/energy_methods/RNA_SugarCloseEnergyCreator.hh>
#include <core/energy_methods/RNA_StubCoordinateEnergyCreator.hh>
#include <core/energy_methods/RNA_SuiteEnergyCreator.hh>
#include <core/energy_methods/TNA_SuiteEnergyCreator.hh>
#include <core/energy_methods/RNA_TorsionEnergyCreator.hh>
#include <core/energy_methods/RNA_VDW_EnergyCreator.hh>
#include <core/energy_methods/RNA_FullAtomVDW_BasePhosphateCreator.hh>
#include <core/energy_methods/StackElecEnergyCreator.hh>
#include <core/energy_methods/RNA_ChemicalShiftEnergyCreator.hh>
#include <core/energy_methods/RNA_ChemicalMappingEnergyCreator.hh>
#include <core/energy_methods/RNA_DataBackboneEnergyCreator.hh>
#include <core/energy_methods/LoopCloseEnergyCreator.hh>
#include <core/scoring/sym_e/symECreator.hh>
#include <core/energy_methods/FastSAXSEnergyCreator.hh>
#include <core/energy_methods/SAXSEnergyCreator.hh>
#include <core/energy_methods/SAXSEnergyCreatorFA.hh>
#include <core/energy_methods/SAXSEnergyCreatorCEN.hh>
#include <core/energy_methods/FiberDiffractionEnergyCreator.hh>
#include <core/energy_methods/FiberDiffractionEnergyDensCreator.hh>
#ifdef USECUDA
#include <core/energy_methods/FiberDiffractionEnergyGpuCreator.hh>
#endif
#include <core/energy_methods/PointWaterEnergyCreator.hh>

// Rotamer Library registration
#include <core/chemical/rotamers/RotamerLibrarySpecificationRegistrator.hh>
#include <core/chemical/rotamers/BasicRotamerLibrarySpecificationCreator.hh>
#include <core/chemical/rotamers/CenrotRotamerLibrarySpecificationCreator.hh>
#include <core/chemical/rotamers/DunbrackRotamerLibrarySpecificationCreator.hh>
#include <core/chemical/rotamers/NCAARotamerLibrarySpecificationCreator.hh>
#include <core/chemical/rotamers/PDBRotamerLibrarySpecificationCreator.hh>

#include <core/pack/dunbrack/SingleResidueDunbrackLibraryCreator.hh>
#include <core/pack/dunbrack/cenrot/SingleResidueCenrotLibraryCreator.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibraryRegistrator.hh>
#include <core/pack/rotamers/SingleBasicRotamerLibraryCreator.hh>
#include <core/pack/rotamers/SingleLigandRotamerLibraryCreator.hh>
#include <core/pack/rotamers/SingleNCAARotamerLibraryCreator.hh>
#include <core/pack/rotamers/StoredRotamerLibraryCreator.hh>

// Constraint registration
#include <core/scoring/constraints/ConstraintFactory.hh>
#include <core/scoring/constraints/BasicConstraintCreators.hh>

#include <core/pack/dunbrack/DunbrackConstraintCreator.hh>
#include <core/scoring/constraints/SequenceProfileConstraintCreator.hh>

// SilentStruct registration
#include <core/io/silent/SilentStructFactory.hh>
#include <core/io/silent/BasicSilentStructCreators.hh>
#include <core/import_pose/PDBSilentStructCreator.hh>

// Sequence registration
#include <core/sequence/SequenceFactory.hh>
#include <core/sequence/BasicSequenceCreators.hh>

// for registering TaskOperations, ResLvlTaskOperations, and ResFilters
#include <core/pack/task/operation/TaskOperationRegistrator.hh>
#include <core/pack/task/operation/TaskOperationCreators.hh>
#include <core/pack/task/operation/DesignRestrictionsCreator.hh>
#include <core/pack/task/operation/ClashBasedRepackShellCreator.hh>
#include <core/pack/task/operation/EnableMultiCoolAnnealerCreator.hh>
#include <core/pack/task/operation/EnableSmartAnnealerCreator.hh>
#include <core/pack/task/operation/KeepSequenceSymmetryCreator.hh>
#include <core/pack/task/operation/OperateOnCertainResiduesCreator.hh>
#include <core/pack/task/operation/OperateOnResidueSubsetCreator.hh>
#include <core/pack/task/operation/NoRepackDisulfidesCreator.hh>
#include <core/pack/task/operation/OptCysHGCreator.hh>
#include <core/pack/task/operation/OptHCreator.hh>
#include <core/pack/task/operation/ResLvlTaskOperationRegistrator.hh>
#include <core/pack/task/operation/ResLvlTaskOperationCreators.hh>
#include <core/pack/task/operation/ResFilterRegistrator.hh>
#include <core/pack/task/operation/ResFilterCreators.hh>
#include <core/pack/task/operation/RestrictInteractionGraphThreadsOperationCreator.hh>
// (end for registering TaskOperations, ResLvlTaskOperations, and ResFilters)

// ResidueSelectors
#include <core/select/residue_selector/AsymmetricUnitSelectorCreator.hh>
#include <core/select/residue_selector/BFactorSelectorCreator.hh>
#include <core/select/residue_selector/GlycanLayerSelectorCreator.hh>
#include <core/select/residue_selector/ResidueSelectorCreators.hh>
#include <core/select/residue_selector/ResiduePropertySelectorCreator.hh>
#include <core/select/residue_selector/PrimarySequenceNeighborhoodSelectorCreator.hh>
#include <core/select/residue_selector/SymmetricalResidueSelectorCreator.hh>
#include <core/pack/task/residue_selector/ClashBasedShellSelectorCreator.hh>
#include <core/select/residue_selector/CloseContactResidueSelectorCreator.hh>
#include <core/select/residue_selector/LogicResidueSelectorCreator.hh>

#include <core/select/residue_selector/ResidueSelectorRegistrator.hh>

// for creating and registering PackerPalettes
#include <core/pack/palette/PackerPaletteRegistrator.hh>
#include <core/pack/palette/CustomBaseTypePackerPaletteCreator.hh>
#include <core/pack/palette/DefaultPackerPaletteCreator.hh>
#include <core/pack/palette/NoDesignPackerPaletteCreator.hh>

// for creating and registering JumpSelectors:
#include <core/select/jump_selector/AndJumpSelectorCreator.hh>
#include <core/select/jump_selector/ExclusivelySharedJumpSelectorCreator.hh>
#include <core/select/jump_selector/InterchainJumpSelectorCreator.hh>
#include <core/select/jump_selector/JumpForResidueCreator.hh>
#include <core/select/jump_selector/JumpIndexSelectorCreator.hh>
#include <core/select/jump_selector/JumpSelectorRegistrator.hh>
#include <core/select/jump_selector/NotJumpSelectorCreator.hh>
#include <core/select/jump_selector/OrJumpSelectorCreator.hh>

#endif

#ifndef __native_client__
#ifndef WIN_PYROSETTA
#endif

#if defined(MAC) || defined(__APPLE__)  ||  defined(__OSX__)
#include <mach-o/dyld.h> // for _NSGetExecutablePath
#endif
#endif

#if defined( MAC ) || defined( __APPLE__ ) || defined( __OSX__ )
#include <sys/resource.h> // for getrlimit/setrlimit
#define PROCESS_STACK_SIZE (32 * 1024 * 1024) // 32 MB
#else // For linux, for example, we still need a bigger stack for ribosomes.
#ifndef WIN32
#include <sys/resource.h> // for getrlimit/setrlimit
#define PROCESS_STACK_SIZE (32 * 1024 * 1024) // 32 MB
#endif
#endif

// Windows headers
#if (defined WIN32) && (!defined WIN_PYROSETTA)
#include <windows.h>
#include <wincrypt.h>
#endif

// STL headers
#include <sstream>
#include <string>
#include <cstring>

// option key includes

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/testing.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/OptionKeys.hh>

// ResourceManager includes
#include <basic/resource_manager/ResourceLoaderRegistrator.hh>
#include <core/conformation/symmetry/SymmDataLoaderCreator.hh>
#include <core/import_pose/PoseResourceLoaderCreator.hh>

#include <core/init/init.ResourceLocatorCreators.ihh> // DO NOT AUTO-REMOVE
#include <core/init/init.ResourceLocatorRegistrators.ihh> // DO NOT AUTO-REMOVE


//#include <core/energy_methods/ElectronDensityLoaderCreator.hh>
//#include <core/energy_methods/FiberDiffractionLoaderCreator.hh>
//#include <core/chemical/ResidueLoaderCreator.hh>

//option key includes for deprecated pdbs
#include <basic/options/keys/LoopModel.OptionKeys.gen.hh>
#include <basic/options/keys/antibody.OptionKeys.gen.hh>

#include <core/pack/task/operation/ResFilterFactory.hh>
#include <core/pack/task/operation/ResLvlTaskOperationFactory.hh>
#include <core/pack/task/operation/TaskOperationFactory.hh>
#include <core/scoring/ScoringManager.hh>

#include <basic/init.hh>

#include <utility/version.hh>
#include <utility/vector1.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/CSI_Sequence.hh>
#include <utility/crash_report.hh>
#include <utility/options/OptionCollection.hh>


#ifdef UNICODE
typedef std::wostringstream ostringstream_t;
#else
using ostringstream_t = std::ostringstream;
#endif

namespace core {
namespace init {

static basic::Tracer TR("core.init");

/// The following global varialbles force the linker to always include
/// the EnergyMethodCreator files to be included in staticly linked
/// executables.  These variables will be initialized before main()
/// begins during the "dynamic initialization" phase of loading.
/// During this time, the Registrotor classes will register their templated
/// EnergyMethods with the ScoringManager before main() begins.

using namespace core::scoring::methods;
using namespace core::scoring::rna;

static EnergyMethodRegistrator< energy_methods::AARepeatEnergyCreator > AARepeatEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::AACompositionEnergyCreator > AACompositionEnergyCreator_registrator;
static EnergyMethodRegistrator< pack::guidance_scoreterms::approximate_buried_unsat_penalty::ApproximateBuriedUnsatPenaltyCreator > ApproximateBuriedUnsatPenaltyCreator_registrator;
static EnergyMethodRegistrator< scoring::constraints::ConstraintsEnergyCreator > ConstraintsEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods ::CentroidDisulfideEnergyCreator > CentroidDisulfideEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::DisulfideMatchingEnergyCreator > DisulfideMatchingEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::FullatomDisulfideEnergyCreator > FullatomDisulfideEnergyCreator_registrator;
static EnergyMethodRegistrator< scoring::etable::EtableEnergyCreator > EtableEnergyCreator_registrator;
static EnergyMethodRegistrator< scoring::etable::EtableClassicIntraEnergyCreator > EtableClassicIntraEnergyCreator_registrator;
static EnergyMethodRegistrator< pack::guidance_scoreterms::hbnet_energy::HBNetEnergyCreator > HBNetEnergyCreator_registrator;
static EnergyMethodRegistrator< scoring::hbonds::HBondEnergyCreator > HBondEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::MHCEpitopeEnergyCreator > MHCEpitopeEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::NPDHBondEnergyCreator > NPDHBondEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::ArgCationPiEnergyCreator > ArgCationPiEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::AspartimidePenaltyEnergyCreator > AspartimidePenaltyEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::ChainbreakEnergyCreator > ChainbreakEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::BranchEnergyCreator > BranchEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::CenPairEnergyCreator > CenPairEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::CenPairMotifEnergyCreator > CenPairMotifEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::CenPairMotifDegreeEnergyCreator > CenPairMotifDegreeEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::ContactOrderEnergyCreator > ContactOrderEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::IntermolEnergyCreator > IntermolEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::MissingEnergyCreator > MissingEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::EnvEnergyCreator > EnvEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::EnvSmoothEnergyCreator > EnvSmoothEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::OmegaTetherEnergyCreator > OmegaTetherEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::OtherPoseEnergyCreator > OtherPoseEnergyCreator_registrator;
static EnergyMethodRegistrator< core::scoring::methods::NMerRefEnergyCreator > NMerRefEnergyCreator_registrator;
static EnergyMethodRegistrator< core::scoring::methods::NMerPSSMEnergyCreator > NMerPSSMEnergyCreator_registrator;
static EnergyMethodRegistrator< core::scoring::methods::NMerSVMEnergyCreator > NMerSVMEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::P_AA_EnergyCreator > P_AA_EnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::P_AA_ss_EnergyCreator > P_AA_ss_EnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::P_AA_pp_EnergyCreator > P_AA_pp_EnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::AbegoEnergyCreator > AbegoEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::PackStatEnergyCreator > PackStatEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::PairEnergyCreator > PairEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::PeptideBondEnergyCreator > PeptideBondEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::ProClosureEnergyCreator > ProClosureEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::BurialEnergyCreator > BurialCreator_registrator;
static EnergyMethodRegistrator< energy_methods::Burial_v2EnergyCreator > Burial_v2Creator_registrator;
static EnergyMethodRegistrator< energy_methods::CovalentLabelingEnergyCreator> CovalentLabelingEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::CovalentLabelingFAEnergyCreator> CovalentLabelingFAEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::DEPC_MS_EnergyCreator > DEPC_MS_EnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::HRFDynamicsEnergyCreator > HRFDynamicsEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::HRF_MSLabelingEnergyCreator > HRF_MSLabelingEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::CCS_IMMSEnergyCreator > CCS_IMMSEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::RG_Energy_FastCreator > RG_Energy_FastCreator_registrator;
static EnergyMethodRegistrator< energy_methods::RG_LocalEnergyCreator > RG_LocalEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::SA_EnergyCreator > SA_EnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::SSElementMotifContactEnergyCreator > SSElementMotifContactEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::D2H_SA_EnergyCreator > D2H_SA_EnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::ProQ_EnergyCreator > ProQ_EnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::RMS_EnergyCreator > RMS_EnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::RamaPreProEnergyCreator > RamaPreProEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::RamachandranEnergy2BCreator > RamachandranEnergy2BCreator_registrator;
static EnergyMethodRegistrator< energy_methods::RamachandranEnergyCreator > RamachandranEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::ReferenceEnergyCreator > ReferenceEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::ReferenceEnergyNoncanonicalCreator > ReferenceEnergyNoncanonicalCreator_registrator;
static EnergyMethodRegistrator< energy_methods::SecondaryStructureEnergyCreator > SecondaryStructureEnergyCreator_registrator;
static EnergyMethodRegistrator< pack::guidance_scoreterms::sap::SapConstraintEnergyCreator > SapConstraintEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::dfire::DFIRE_EnergyCreator > DFIRE_EnergyCreator_registrator;
static EnergyMethodRegistrator< scoring::vdwaals::VDW_EnergyCreator > VDW_EnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::GoapEnergyCreator > GoapEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::RingClosureEnergyCreator > RingClosureEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::AromaticBackboneRestraintEnergyCreator > AromaticBackboneRestraintEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::SugarBackboneEnergyCreator > SugarBackboneEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::DumpTrajectoryEnergyCreator > DumpTrajectoryEnergy_registrator;
static EnergyMethodRegistrator< energy_methods::NetChargeEnergyCreator > NetChargeEnergyCreator_registrator;
static EnergyMethodRegistrator< pack::guidance_scoreterms::buried_unsat_penalty::BuriedUnsatPenaltyCreator > BuriedUnsatPenaltyCreator_registrator;
static EnergyMethodRegistrator< pack::guidance_scoreterms::voids_penalty_energy::VoidsPenaltyEnergyCreator > VoidsPenaltyEnergyCreator_registrator;

static EnergyMethodRegistrator< pack::dunbrack::DunbrackEnergyCreator > DunbrackEnergyCreator_registrator;
static EnergyMethodRegistrator< pack::dunbrack::cenrot::CenRotDunEnergyCreator > CenRotDunEnergyCreator_registrator;

// define this for compiling a slimmed down version of mini libraries lacking about 3/4s of the code
// this is required for compiling a less memory hungry version of mini for Bluegene etc..
#ifndef MINI_SLIM
static EnergyMethodRegistrator< scoring::nv::NVscoreCreator > NVscoreCreator_registrator;
static EnergyMethodRegistrator< scoring::orbitals::OrbitalsScoreCreator > OrbitalsScoreCreator_registrator;
static EnergyMethodRegistrator< scoring::interface_::DDPscoreCreator > DDPscoreCreator_registrator;
static EnergyMethodRegistrator< scoring::carbon_hbonds::CarbonHBondEnergyCreator > CarbonHBondEnergyCreator_registrator;
static EnergyMethodRegistrator< scoring::custom_pair_distance::FullatomCustomPairDistanceEnergyCreator > FullatomCustomPairDistanceEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::CustomAtomPairEnergyCreator > CustomAtomPairEnergy_registrator;
static EnergyMethodRegistrator< energy_methods::FastDensEnergyCreator > FastDensEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::ElecDensCenEnergyCreator > ElecDensCenEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::ElecDensAllAtomCenEnergyCreator > ElecDensAllAtomCenEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::ElecDensEnergyCreator > ElecDensEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::XtalMLEnergyCreator > XtalMLEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::ElecDensAtomwiseEnergyCreator > ElecDensAtomwiseEnergyCreator_registrator;
//XRW_B_T1
//static EnergyMethodRegistrator< scoring::etable::CoarseEtableEnergyCreator > CoarseEtableEnergyCreator_registrator;
//XRW_E_T1
static EnergyMethodRegistrator< energy_methods::ExactOccludedHbondSolEnergyCreator > ExactOccludedHbondSolEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::ContextDependentGeometricSolEnergyCreator > ContextDependentGeometricSolEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::ContextIndependentGeometricSolEnergyCreator > ContextIndependentGeometricSolEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::OccludedHbondSolEnergyCreator > OccludedHbondSolEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::OccludedHbondSolEnergy_onebodyCreator > OccludedHbondSolEnergy_onebodyCreator_registrator;
static EnergyMethodRegistrator< scoring::elec::FA_ElecEnergyCreator > FA_ElecEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::ImplicitMembraneElecEnergyCreator > ImplicitMembraneElecEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::RNA_FA_ElecEnergyCreator > RNA_FA_ElecEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::FA_ElecEnergyAroAroCreator > FA_ElecEnergyAroAroCreator_registrator;
static EnergyMethodRegistrator< energy_methods::FA_ElecEnergyAroAllCreator > FA_ElecEnergyAroAllCreator_registrator;
static EnergyMethodRegistrator< energy_methods::FA_GrpElecEnergyCreator > FA_GrpElecEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::HackAroEnergyCreator > HackAroEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::DEEREnergyCreator > DEEREnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::DNAChiEnergyCreator > DNAChiEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::DNATorsionEnergyCreator > DNATorsionEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::LoopCloseEnergyCreator > LoopCloseEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::MgEnergyCreator > MgEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::RNA_MgPointEnergyCreator > RNA_MgPointEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::DNA_BaseEnergyCreator > DNA_BaseEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::DNA_DihedralEnergyCreator > DNA_DihedralEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::DirectReadoutEnergyCreator > DirectReadoutEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::DistanceChainbreakEnergyCreator > DistanceChainbreakEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::Fa_MbenvEnergyCreator > Fa_MbenvEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::Fa_MbsolvEnergyCreator > Fa_MbsolvEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::FreeDOF_EnergyCreator > FreeDOF_EnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::GenBornEnergyCreator > GenBornEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::GenericBondedEnergyCreator > GenericBondedEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::VdWTinkerEnergyCreator > VdWTinkerEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::MultipoleElecEnergyCreator > MultipoleElecEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::SASAEnergyCreator > SASAEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::FACTSEnergyCreator > FACTSEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::HybridVDW_EnergyCreator > HybridVDW_EnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::LK_PolarNonPolarEnergyCreator > LK_PolarNonPolarEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::LK_hackCreator > LK_hackCreator_registrator;
static EnergyMethodRegistrator< scoring::lkball::LK_BallEnergyCreator > LK_BallEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::LinearChainbreakEnergyCreator > LinearChainbreakEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::LinearBranchEnergyCreator > LinearBranchEnergyCreator_registrator;
static EnergyMethodRegistrator< scoring::lkball::LK_DomeEnergyCreator > LK_DomeEnergyCreator_registrator;
static EnergyMethodRegistrator< pack::guidance_scoreterms::lk_dome::LK_DomePackEnergyCreator > LK_DomePackEnergyCreator_registrator;
static EnergyMethodRegistrator< scoring::methods::MMBondAngleEnergyCreator > MMBondAngleEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::CartesianBondedEnergyCreator > CartesianBondedEnergyCreator_registrator;
static EnergyMethodRegistrator< scoring::methods::MMBondLengthEnergyCreator > MMBondLengthEnergyCreator_registrator;
static EnergyMethodRegistrator< scoring::methods::MMTorsionEnergyCreator > MMTorsionEnergyCreator_registrator;
static EnergyMethodRegistrator< scoring::methods::MMLJEnergyInterCreator > MMLJEnergyInterCreator_registrator;
static EnergyMethodRegistrator< scoring::methods::MMLJEnergyIntraCreator > MMLJEnergyIntraCreator_registrator;
static EnergyMethodRegistrator< energy_methods::MembraneCbetaEnergyCreator > MembraneCbetaEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::MembraneCenPairEnergyCreator > MembraneCenPairEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::MembraneEnvEnergyCreator > MembraneEnvEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::MembraneEnvPenaltiesCreator > MembraneEnvPenaltiesCreator_registrator;
static EnergyMethodRegistrator< energy_methods::MPPairEnergyCreator >
	MPPairEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::MPEnvEnergyCreator >
	MPEnvEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::MPCbetaEnergyCreator >
	MPCbetaEnergyCreator_registrator;
static EnergyMethodRegistrator< scoring::membrane::MPNonHelixPenaltyCreator >
	MPNonHelixPenaltyCreator_registrator;
static EnergyMethodRegistrator< scoring::membrane::MPTerminiPenaltyCreator >
	MPTerminiPenaltyCreator_registrator;
static EnergyMethodRegistrator< scoring::membrane::MPTMProjPenaltyCreator >
	MPTMProjPenaltyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::FaMPEnvEnergyCreator >
	FaMPEnvEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::FaMPSolvEnergyCreator >
	FaMPSolvEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::FaMPEnvSmoothEnergyCreator >
	FaMPEnvSMoothEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::FaMPAsymEzCBEnergyCreator > FaMPAsymEzCBEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::FaMPAsymEzCGEnergyCreator > FaMPAsymEzCGEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::MPResidueLipophilicityEnergyCreator > MPResidueLipophilicityEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::MPHelicalityEnergyCreator > MPHelicalityEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::MPSpanInsertionEnergyCreator > MPSpanInsertionEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::MPSpanAngleEnergyCreator > MPSpanAngleEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::SmoothCenPairEnergyCreator > SmoothCenPairEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::SmoothEnvEnergyCreator > SmoothEnvEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::CenRotPairEnergyCreator > CenRotPairEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::CenRotEnvEnergyCreator > CenRotEnvEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::CenHBEnergyCreator > CenHBEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::MembraneLipoCreator > MembraneLipoCreator_registrator;
static EnergyMethodRegistrator< energy_methods::MembraneEnvSmoothEnergyCreator > MembraneEnvSmoothEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::pHEnergyCreator > pHEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::ChemicalShiftAnisotropyEnergyCreator > ChemicalShiftAnisotropyEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::DipolarCouplingEnergyCreator > DipolarCouplingEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::ResidualDipolarCouplingEnergyCreator > ResidualDipolarCouplingEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::ResidualDipolarCouplingEnergy_RohlCreator > ResidualDipolarCouplingEnergy_RohlCreator_registrator;
static EnergyMethodRegistrator< energy_methods::MotifDockEnergyCreator > MotifDockEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::SuckerEnergyCreator > SuckerEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::GaussianOverlapEnergyCreator > GaussianOverlapEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::YHHPlanarityEnergyCreator > YYHPlanarityEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::HydroxylTorsionEnergyCreator > HydroxylTorsionEnergyCreator_registrator;
static EnergyMethodRegistrator< pack::interaction_graph::SurfaceEnergyCreator > SurfaceEnergyCreator_registrator;
static EnergyMethodRegistrator< pack::interaction_graph::HPatchEnergyCreator > HPatchEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::SymmetricLigandEnergyCreator > SymmetricLigandEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::UnfoldedStateEnergyCreator > UnfoldedStateEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::SplitUnfoldedTwoBodyEnergyCreator > SplitUnfoldedTwoBodyEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::WaterAdductHBondEnergyCreator > WaterAdductHBondEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::WaterAdductIntraEnergyCreator > WaterAdductIntraEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::WaterSpecificEnergyCreator > WaterSpecificEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::HolesEnergyCreator > HolesEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::SurfVolEnergyCreator > SurfVolEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::SurfEnergyCreator > SurfEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::RG_Energy_RNACreator > RG_Energy_RNACreator_registrator;
static EnergyMethodRegistrator< energy_methods::RNA_BulgeEnergyCreator > RNA_BulgeEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::RNA_ChemicalMappingEnergyCreator > RNA_ChemicalMappingEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::RNA_DataBackboneEnergyCreator > RNA_DataBackboneEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::RNA_FullAtomStackingEnergyCreator > RNA_FullAtomStackingEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::RNA_JR_SuiteEnergyCreator > RNA_JR_SuiteEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::RNA_LJ_BaseEnergyCreator > RNA_LJ_BaseEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::RNA_PairwiseLowResolutionEnergyCreator > RNA_PairwiseLowResolutionEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::RNA_PartitionEnergyCreator > RNA_PartitionEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::RNP_LowResEnergyCreator > RNP_LowResEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::RNP_LowResPairDistEnergyCreator > RNP_LowResPairDistEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::RNP_LowResStackEnergyCreator > RNP_LowResStackEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::RNA_StubCoordinateEnergyCreator > RNA_StubCoordinateEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::RNA_SugarCloseEnergyCreator > RNA_SugarCloseEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::RNA_SuiteEnergyCreator > RNA_SuiteEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::TNA_SuiteEnergyCreator > TNA_SuiteEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::RNA_TorsionEnergyCreator > RNA_TorsionEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::RNA_VDW_EnergyCreator > RNA_VDW_EnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::RNA_FullAtomVDW_BasePhosphateCreator > RNA_FullAtomVDW_BasePhosphateCreator_registrator;
static EnergyMethodRegistrator< energy_methods::StackElecEnergyCreator > StackElecEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::RNA_ChemicalShiftEnergyCreator > NA_ChemicalShiftEnergyCreator_registrator;
//static EnergyMethodRegistrator< energy_methods::FreeResidueBonusEnergyCreator > FreeResidueBonusEnergyCreator_registrator;
static EnergyMethodRegistrator< scoring::sym_e::symECreator > symECreator_registrator;
static EnergyMethodRegistrator< energy_methods::PoissonBoltzmannEnergyCreator > PoissonBoltzmannEnergyCreator_registrator;

static EnergyMethodRegistrator< core::energy_methods::RNA_CoarseDistEnergyCreator > RNA_CoarseDistEnergyCreator_registrator;

static EnergyMethodRegistrator< energy_methods::FastSAXSEnergyCreator > FastSAXSEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::SAXSEnergyCreatorCEN > SAXSEnergyCreatorCEN_registrator;
static EnergyMethodRegistrator< energy_methods::SAXSEnergyCreatorFA > SAXSEnergyCreatorFA_registrator;
static EnergyMethodRegistrator< energy_methods::SAXSEnergyCreator > SAXSEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::PointWaterEnergyCreator > PointWaterEnergyCreator_registrator;

static EnergyMethodRegistrator< energy_methods::FiberDiffractionEnergyCreator > FiberDiffractionEnergyCreator_registrator;
static EnergyMethodRegistrator< energy_methods::FiberDiffractionEnergyDensCreator > FiberDiffractionEnergyDensCreator_registrator;
#ifdef USECUDA
static EnergyMethodRegistrator< energy_methods::FiberDiffractionEnergyGpuCreator > FiberDiffractionEnergyGpuCreator_registrator;
#endif

/// RotamerLibrary Registrators

using chemical::rotamers::RotamerLibrarySpecificationRegistrator;
static RotamerLibrarySpecificationRegistrator< chemical::rotamers::BasicRotamerLibrarySpecificationCreator > BasicRotamerLibrarySpecificationCreator_registrator;
static RotamerLibrarySpecificationRegistrator< chemical::rotamers::CenrotRotamerLibrarySpecificationCreator > CenrotRotamerLibrarySpecificationCreator_registrator;
static RotamerLibrarySpecificationRegistrator< chemical::rotamers::DunbrackRotamerLibrarySpecificationCreator > DunbrackRotamerLibrarySpecificationCreator_registrator;
static RotamerLibrarySpecificationRegistrator< chemical::rotamers::NCAARotamerLibrarySpecificationCreator > NCAARotamerLibrarySpecificationCreator_registrator;
static RotamerLibrarySpecificationRegistrator< chemical::rotamers::PDBRotamerLibrarySpecificationCreator > PDBRotamerLibrarySpecificationCreator_registrator;

using pack::rotamers::SingleResidueRotamerLibraryRegistrator;
static SingleResidueRotamerLibraryRegistrator< pack::dunbrack::SingleResidueDunbrackLibraryCreator > SingleResidueDunbrackLibraryCreator_registrator;
static SingleResidueRotamerLibraryRegistrator< pack::rotamers::SingleBasicRotamerLibraryCreator > SingleBasicRotamerLibraryCreator_registrator;
static SingleResidueRotamerLibraryRegistrator< pack::dunbrack::cenrot::SingleResidueCenrotLibraryCreator > SingleCenrotRotamerLibraryCreator_registrator;
static SingleResidueRotamerLibraryRegistrator< pack::rotamers::SingleLigandRotamerLibraryCreator > SingleLigandRotamerLibraryCreator_registrator;
static SingleResidueRotamerLibraryRegistrator< pack::rotamers::SingleNCAARotamerLibraryCreator > SingleNCAARotamerLibraryCreator_registrator;
static SingleResidueRotamerLibraryRegistrator< pack::rotamers::StoredRotamerLibraryCreator > StoredRotamerLibraryCreator_registrator;

/// Constraint Registrators
using namespace scoring::constraints;
static ConstraintRegistrator< core::scoring::constraints::AmbiguousConstraintCreator > AmbiguousConstraintCreator_registrator;
static ConstraintRegistrator< core::scoring::constraints::AmbiguousNMRConstraintCreator > AmbiguousNMRConstraintCreator_registrator;
static ConstraintRegistrator< core::scoring::constraints::AmbiguousNMRDistanceConstraintCreator > AmbiguousNMRDistanceConstraintCreator_registrator;
static ConstraintRegistrator< core::scoring::constraints::AngleConstraintCreator > AngleConstraintCreator_registrator;
static ConstraintRegistrator< core::scoring::constraints::AtomPairConstraintCreator > AtomPairConstraintCreator_registrator;
static ConstraintRegistrator< core::scoring::constraints::BasePairConstraintCreator > BasePairConstraintCreator_registrator;
static ConstraintRegistrator< core::scoring::constraints::BigBinConstraintCreator > BigBinConstraintCreator_registrator;
static ConstraintRegistrator< core::scoring::constraints::CoordinateConstraintCreator > CoordinateConstraintCreator_registrator;
static ConstraintRegistrator< core::scoring::constraints::DihedralConstraintCreator > DihedralConstraintCreator_registrator;
static ConstraintRegistrator< core::scoring::constraints::DihedralPairConstraintCreator > DihedralPairConstraintCreator_registrator;
static ConstraintRegistrator< core::scoring::constraints::KofNConstraintCreator > KofNConstraintCreator_registrator;
static ConstraintRegistrator< core::scoring::constraints::LocalCoordinateConstraintCreator > LocalCoordinateConstraintCreator_registrator;
static ConstraintRegistrator< core::scoring::constraints::MultiConstraintCreator > MultiConstraintCreator_registrator;
static ConstraintRegistrator< core::scoring::constraints::NamedAtomPairConstraintCreator > NamedAtomPairConstraintCreator_registrator;
static ConstraintRegistrator< core::scoring::constraints::NamedAngleConstraintCreator > NamedAngleConstraintCreator_registrator;
static ConstraintRegistrator< core::scoring::constraints::FabConstraintCreator > FabConstraintCreator_registrator;

static ConstraintRegistrator< core::pack::dunbrack::DunbrackConstraintCreator > DunbrackConstraintCreator_registrator;
static ConstraintRegistrator< core::scoring::constraints::SiteConstraintCreator > SiteConstraintCreator_registrator;
static ConstraintRegistrator< core::scoring::constraints::SiteConstraintResiduesCreator > SiteConstraintResiduesCreator_registrator;
static ConstraintRegistrator< core::scoring::constraints::SequenceProfileConstraintCreator > SequenceProfileConstraintCreator_registrator;


// SilentStruct registrators
using namespace core::io::silent;
static SilentStructRegistrator< core::io::silent::ProteinSilentStructCreator > ProteinSilentStructCreator_registrator;
static SilentStructRegistrator< core::io::silent::ProteinSilentStruct_SinglePrecCreator > ProteinSilentStruct_SinglePrecCreator_registrator;
static SilentStructRegistrator< core::io::silent::RNA_SilentStructCreator > RNA_SilentStructCreator_registrator;
static SilentStructRegistrator< core::import_pose::PDBSilentStructCreator > PDBSilentStructCreator_registrator;
static SilentStructRegistrator< core::io::silent::BinarySilentStructCreator > BinarySilentStructCreator_registrator;
static SilentStructRegistrator< core::io::silent::ScoreFileSilentStructCreator > ScoreFileSilentStructCreator_registrator;
static SilentStructRegistrator< core::io::silent::ScoreJumpFileSilentStructCreator > ScoreJumpFileSilentStructCreator_registrator;
static SilentStructRegistrator< core::io::silent::RigidBodySilentStructCreator > RigidBodySilentStructCreator_registrator;

// Sequence registrators
using namespace core::sequence;
static SequenceRegistrator< core::sequence::SimpleSequenceCreator > SimpleSequenceCreator_registrator;
static SequenceRegistrator< core::sequence::SequenceProfileCreator > SequenceProfileCreator_registrator;
static SequenceRegistrator< core::sequence::SequenceCouplingCreator > SequenceCouplingCreator_registrator;
static SequenceRegistrator< core::sequence::CompositeSequenceCreator > CompositeSequenceCreator_registrator;

// perhaps these long registration lists live in their own separate file?
using namespace pack::task::operation;
// register TaskOperationCreators
static TaskOperationRegistrator< RestrictYSDesignCreator > RestrictYSDesignCreator_registrator;
static TaskOperationRegistrator< PreventRepackingCreator > PreventRepackingCreator_registrator;
static TaskOperationRegistrator< PreserveCBetaCreator > PreserveCBetaCreator_registrator;
static TaskOperationRegistrator< AppendRotamerSetCreator > AppendRotamerSetCreator_registrator;
static TaskOperationRegistrator< AppendRotamerCreator > AppendRotamerCreator_registrator;
static TaskOperationRegistrator< DesignRestrictionsCreator > DesignRestrictionsCreator_registrator;
static TaskOperationRegistrator< EnableMultiCoolAnnealerCreator > EnableMultiCoolAnnealerCreator_registrator;
static TaskOperationRegistrator< EnableSmartAnnealerCreator > EnableSmartAnnealerCreator_registrator;
static TaskOperationRegistrator< ExtraRotamersCreator > ExtraRotamersCreator_registrator;
static TaskOperationRegistrator< ExtraChiCutoffCreator > ExtraChiCutoffCreator_registrator;
static TaskOperationRegistrator< SetRotamerCouplingsCreator > SetRotamerCouplingsCreator_registrator;
static TaskOperationRegistrator< SetRotamerLinksCreator > SetRotamerLinksCreator_registrator;
static TaskOperationRegistrator< ReadResfileCreator > ReadResfileCreator_registrator;
static TaskOperationRegistrator< ReadResfileAndObeyLengthEventsCreator > ReadResfileAndObeyLengthEventsCreator_registrator;
static TaskOperationRegistrator< IncludeCurrentCreator > IncludeCurrentCreator_registrator;
static TaskOperationRegistrator< InitializeExtraRotsFromCommandlineCreator > InitializeExtraRotsFromCommandlineCreator_registrator;
static TaskOperationRegistrator< InitializeFromCommandlineCreator > InitializeFromCommandlineCreator_registrator;
static TaskOperationRegistrator< UseMultiCoolAnnealerCreator > UseMultiCoolAnnealerCreator_registrator;
static TaskOperationRegistrator< InitializeFromOptionCollectionCreator > InitializeFromOptionCollectionCreator_registrator;
static TaskOperationRegistrator< ExtraRotamersGenericCreator > ExtraRotamersGenericCreator_registrator;
static TaskOperationRegistrator< RotamerExplosionCreator > RotamerExplosionCreator_registrator;
static TaskOperationRegistrator< RestrictAbsentCanonicalAASCreator > RestrictAbsentCanonicalAASCreator_registrator;
static TaskOperationRegistrator< DisallowIfNonnativeCreator > DisallowIfNonnativeCreator_registrator;
static TaskOperationRegistrator< RestrictToSpecifiedBaseResidueTypesCreator > RestrictToSpecifiedBaseResidueTypesCreator_registrator;
static TaskOperationRegistrator< ProhibitSpecifiedBaseResidueTypesCreator > ProhibitSpecifiedBaseResidueTypesCreator_registrator;
static TaskOperationRegistrator< RestrictToResiduePropertiesCreator > RestrictToResiduePropertiesCreator_registrator;
static TaskOperationRegistrator< ProhibitResiduePropertiesCreator > ProhibitResiduePropertiesCreator_registrator;
static TaskOperationRegistrator< RestrictResidueToRepackingCreator > RestrictResidueToRepackingCreator_registrator;
static TaskOperationRegistrator< RestrictToRepackingCreator > RestrictToRepackingCreator_registrator;
static TaskOperationRegistrator< RestrictInteractionGraphThreadsOperationCreator > RestrictInteractionGraphThreadsOperationCreator_registrator;
static TaskOperationRegistrator< OperateOnCertainResiduesCreator > OperateOnCertainResiduesCreator_registrator;
static TaskOperationRegistrator< OperateOnResidueSubsetCreator > OperateOnResidueSubsetCreator_registrator;
static TaskOperationRegistrator< NoRepackDisulfidesCreator > NoRepackDisulfidesCreator_registrator;
static TaskOperationRegistrator< ClashBasedRepackShellCreator > ClashBasedRepackShellCreator_registrator;
static TaskOperationRegistrator< OptCysHGCreator > OptCysHGCreator_registrator;
static TaskOperationRegistrator< OptHCreator > OptHCreator_registrator;
static TaskOperationRegistrator< KeepSequenceSymmetryCreator > KeepSequenceSymmetryCreator_registrator;
// register ResLvlTaskOperationCreators
static ResLvlTaskOperationRegistrator< RestrictToRepackingRLTCreator > RestrictToRepackingRLTCreator_registrator;
static ResLvlTaskOperationRegistrator< RestrictAbsentCanonicalAASRLTCreator > RestrictAbsentCanonicalAASRLTCreator_registrator;
static ResLvlTaskOperationRegistrator< RestrictAbsentCanonicalAASExceptNativeRLTCreator > RestrictAbsentCanonicalAASExceptNativeRLTCreator_registrator;
static ResLvlTaskOperationRegistrator< DisallowIfNonnativeRLTCreator > DisallowIfNonnativeRLTCreator_registrator;
static ResLvlTaskOperationRegistrator< PreventRepackingRLTCreator > PreventRepackingRLTCreator_registrator;
static ResLvlTaskOperationRegistrator< AddBehaviorRLTCreator > AddBehaviorRLTCreator_registrator;
static ResLvlTaskOperationRegistrator< IncludeCurrentRLTCreator > IncludeCurrentRLTCreator_registrator;
static ResLvlTaskOperationRegistrator< PreserveCBetaRLTCreator > PreserveCBetaRLTCreator_registrator;
static ResLvlTaskOperationRegistrator< ExtraChiCutoffRLTCreator > ExtraChiCutoffRLTCreator_registrator;
static ResLvlTaskOperationRegistrator< ExtraRotamersGenericRLTCreator > ExtraRotamersGenericRLTCreator_registrator;
// register ResFilterCreators
static ResFilterRegistrator< ResidueHasPropertyCreator > ResidueHasPropertyCreator_registrator;
static ResFilterRegistrator< ResiduePDBInfoHasLabelCreator > ResiduePDBInfoHasLabelCreator_registrator;
static ResFilterRegistrator< ResidueLacksPropertyCreator > ResidueLacksPropertyCreator_registrator;
static ResFilterRegistrator< ResidueName3IsCreator > ResidueName3IsCreator_registrator;
static ResFilterRegistrator< ResidueName3IsntCreator > ResidueName3IsntCreator_registrator;
static ResFilterRegistrator< ResidueIndexIsCreator > ResidueIndexIsCreator_registrator;
static ResFilterRegistrator< ResidueIndexIsntCreator > ResidueIndexIsntCreator_registrator;
//static ResFilterRegistrator< ResiduePDBIndexIsCreator > ResiduePDBIndexIsCreator_registrator;
//static ResFilterRegistrator< ResiduePDBIndexIsntCreator > ResiduePDBIndexIsntCreator_registrator;
static ResFilterRegistrator< ChainIsCreator > ChainIsCreator_registrator;
static ResFilterRegistrator< ChainIsntCreator > ChainIsntCreator_registrator;
static ResFilterRegistrator< AnyResFilterCreator > AnyResFilterCreator_registrator;
static ResFilterRegistrator< AllResFilterCreator > AllResFilterCreator_registrator;
static ResFilterRegistrator< NoResFilterCreator > NoResFilterCreator_registrator;
static ResFilterRegistrator< ResidueTypeFilterCreator > ResidueTypeResFilterCreator_registrator;

// register ResidueSelectorCreators
using namespace core::select::residue_selector;
static ResidueSelectorRegistrator< AndResidueSelectorCreator > reg_AndResidueSelectorCreator;
static ResidueSelectorRegistrator< BFactorSelectorCreator > reg_BfactorSelectorCreator;
static ResidueSelectorRegistrator< BinSelectorCreator > reg_BinSelectorCreator;
static ResidueSelectorRegistrator< BondedResidueSelectorCreator > reg_BondedResidueSelectorCreator;
static ResidueSelectorRegistrator< ChainSelectorCreator > reg_ChainSelectorCreator;
static ResidueSelectorRegistrator< core::pack::task::residue_selector::ClashBasedShellSelectorCreator > reg_ClashBasedShellSelectorCreator;
static ResidueSelectorRegistrator< CloseContactResidueSelectorCreator > reg_CloseContactResidueSelectorCreator;
static ResidueSelectorRegistrator< DensityFitResidueSelectorCreator > reg_DensityFitResidueSelectorCreator;
static ResidueSelectorRegistrator< GlycanResidueSelectorCreator > reg_GlycanResidueSelectorCreator;
static ResidueSelectorRegistrator< RandomGlycanFoliageSelectorCreator > reg_RandomGlycanFoliageSelectorCreator;
static ResidueSelectorRegistrator< GlycanLayerSelectorCreator > reg_GlycanLayerSelectorCreator;
static ResidueSelectorRegistrator< GlycanSequonsSelectorCreator > reg_GlycanSequonsSelectorCreator;
static ResidueSelectorRegistrator< InterGroupInterfaceByVectorSelectorCreator > reg_InterGroupInterfaceByVectorSelectorCreator;
static ResidueSelectorRegistrator< JumpDownstreamSelectorCreator > reg_JumpDownstreamSelectorCreator;
static ResidueSelectorRegistrator< JumpUpstreamSelectorCreator > reg_JumpUpstreamSelectorCreator;
static ResidueSelectorRegistrator< AsymmetricUnitSelectorCreator > reg_AsymmetricUnitSelectorCreator;
static ResidueSelectorRegistrator< LogicResidueSelectorCreator > reg_LogicResidueSelectorCreator;
static ResidueSelectorRegistrator< NeighborhoodResidueSelectorCreator > reg_NeighborhoodResidueSelectorCreator;
static ResidueSelectorRegistrator< NotResidueSelectorCreator > reg_NotResidueSelectorCreator;
static ResidueSelectorRegistrator< NumNeighborsSelectorCreator > reg_NumNeighborsSelectorCreator;
static ResidueSelectorRegistrator< LayerSelectorCreator > reg_LayerSelectorCreator;
static ResidueSelectorRegistrator< OrResidueSelectorCreator > reg_OrResidueSelectorCreator;
static ResidueSelectorRegistrator< PhiSelectorCreator > reg_PhiSelectorCreator;
static ResidueSelectorRegistrator< PrimarySequenceNeighborhoodSelectorCreator > reg_PrimarySequenceNeighborhoodSelectorCreator;
static ResidueSelectorRegistrator< RandomResidueSelectorCreator > reg_RandomResidueSelectorCreator;
static ResidueSelectorRegistrator< ResidueIndexSelectorCreator > reg_ResidueIndexSelectorCreator;
static ResidueSelectorRegistrator< ResidueNameSelectorCreator > reg_ResidueNameSelectorCreator;
static ResidueSelectorRegistrator< ResidueInSequenceMotifSelectorCreator > reg_ResidueInSequenceMotifSelectorCreator;
static ResidueSelectorRegistrator< ResiduePropertySelectorCreator > reg_ResiduePropertySelectorCreator;
static ResidueSelectorRegistrator< ResidueSpanSelectorCreator > reg_ResidueSpanSelectorCreator;
static ResidueSelectorRegistrator< TrueResidueSelectorCreator > reg_TrueResidueSelectorCreator;
static ResidueSelectorRegistrator< FalseResidueSelectorCreator > reg_FalseResidueSelectorCreator;
static ResidueSelectorRegistrator< ResidueInMembraneSelectorCreator > reg_ResidueInMembraneSelectorCreator;
static ResidueSelectorRegistrator< ResiduePDBInfoHasLabelSelectorCreator > reg_ResiduePDBInfoHasLabelSelectorCreator;
static ResidueSelectorRegistrator< SecondaryStructureSelectorCreator > reg_SecondaryStructureSelectorCreator;
static ResidueSelectorRegistrator< SSElementSelectorCreator > reg_SSElementSelectorCreator;
static ResidueSelectorRegistrator< SymmetricalResidueSelectorCreator > reg_SymmetricalResidueSelectorCreator;
static ResidueSelectorRegistrator< ScoreTermValueBasedSelectorCreator > reg_ScoreTermValueBasedSelectorCreator;
static ResidueSelectorRegistrator< SliceResidueSelectorCreator > reg_SliceResidueSelectorCreator;
static ResidueSelectorRegistrator< SimpleMetricSelectorCreator > reg_SimpleMetricSelectorCreator;

// Jump Selectors
using namespace core::select::jump_selector;
static JumpSelectorRegistrator< AndJumpSelectorCreator > reg_AndJumpSelectorCreator;
static JumpSelectorRegistrator< ExclusivelySharedJumpSelectorCreator > reg_ExclusivelySharedJumpSelectorCreator;
static JumpSelectorRegistrator< InterchainJumpSelectorCreator > reg_InterchainJumpSelectorCreator;
static JumpSelectorRegistrator< JumpForResidueCreator > reg_JumpForResidueCreator;
static JumpSelectorRegistrator< JumpIndexSelectorCreator > reg_JumpIndexSelectorCreator;
static JumpSelectorRegistrator< NotJumpSelectorCreator > reg_NotJumpSelectorCreator;
static JumpSelectorRegistrator< OrJumpSelectorCreator > reg_OrJumpSelectorCreator;

// register PackerPaletteCreators
using namespace core::pack::palette;
static PackerPaletteRegistrator< CustomBaseTypePackerPaletteCreator > reg_CustomBaseTypePackerPaletteCreator;
static PackerPaletteRegistrator< DefaultPackerPaletteCreator > reg_DefaultPackerPaletteCreator;
static PackerPaletteRegistrator< NoDesignPackerPaletteCreator > reg_NoDesignPackerPaletteCreator;

using basic::resource_manager::ResourceLoaderRegistrator;
static ResourceLoaderRegistrator< core::conformation::symmetry::SymmDataLoaderCreator > SymmDataLoaderCreator_registrator;
static ResourceLoaderRegistrator< core::import_pose::PoseResourceLoaderCreator > PoseResourceLoaderCreator_registrator;
//static ResourceLoaderRegistrator< core::energy_methods::ElectronDensityLoaderCreator > ElectronDensityLoaderCreator_registrator;
//static ResourceLoaderRegistrator< core::energy_methods::FiberDiffractionLoaderCreator > FiberDiffractionLoaderCreator_registrator;
//static ResourceLoaderRegistrator< core::chemical::ResidueLoaderCreator > ResidueLoaderCreator_registrator;

#endif

#ifndef WIN32
#ifndef NDEBUG
#ifndef PYROSETTA
//c++11: 201103L
//c++14: 201402L
//c++17: 201703L
//c++20: not set yet, please update
#ifdef CXX14
static_assert( __cplusplus > 201103L, "CXX14 is defined but the compiler version is <= c++11" );
#else
static_assert( __cplusplus != 201402L, "The compiler is set up to run c++14 but CXX14 is not defined" );
#endif

#ifdef CXX14_OR_LATER
//Using __cplusplus > 201103L instead of __cplusplus >= 201402L because there may be some early c++14 support in the versions between 201103L and 201402L
static_assert( __cplusplus > 201103L, "CXX14_OR_LATER is defined but the compiler version is <= c++11" );
#else
static_assert( __cplusplus < 201402L, "The compiler is set up to run a version of c++ at least as recent as c++14 but CXX14_OR_LATER is not defined" );
#endif

#ifdef CXX17
static_assert( __cplusplus > 201402L, "CXX17 is defined but the compiler version is <= c++14" );
#else
static_assert( __cplusplus != 201703L, "The compiler is set up to run c++17 but CXX17 is not defined" );
#endif

#ifdef CXX17_OR_LATER
static_assert( __cplusplus > 201402L, "CXX17_OR_LATER is defined but the compiler version is <= c++14" );
#else
static_assert( __cplusplus < 201703L, "The compiler is set up to run a version of c++ at least as recent as c++17 but CXX17_OR_LATER is not defined" );
#endif

#endif//PYROSETTA
#endif//NDEBUG
#endif//WIN32

using namespace basic::options;
using namespace basic::options::OptionKeys;

#ifdef USEMPI
void
init_mpi(int argc, char * argv []) {
	int already_initialized( 0 );
	MPI_Initialized( & already_initialized );
	if ( already_initialized == 0 ) MPI_Init(&argc, &argv);

}
#else
void
init_mpi(int, char **) {}
#endif


void
init_options(int argc, char * argv []) {

	// initialize options
	initialize().load( argc, argv, false /* no "free" cmd line args (just discarded anyway) */ );

	pre_tracer_process();
}

/// @brief Setup crash reporter
void init_crash_reporter(int argc, char * argv [])
{
	if ( argc >= 1 ) {
		utility::set_application_name(argv[0]);
	}

	std::stringstream option_stream;
	option_stream << basic::options::option;
	utility::set_options_string( option_stream.str() );

	if ( option[ run::crash_to_console ].user() ) { // If option is explicitly set, use that.
		utility::set_show_crash_report_on_console( option[ run::crash_to_console ]() );
	} else if ( utility::Version::public_release() ) {
		utility::set_show_crash_report_on_console( false ); // Don't show tracer backtrace for release versions
	} else {
		utility::set_show_crash_report_on_console( true ); // Developers seem to want the backtrace in the tracer.
	}

	if ( ! option[ run::nosignal ]() ) {
		utility::install_crash_handler(); // Install the signal listener
	}
}

/// @brief After the tracers have been initialized, now go back and modify some of the
/// values in the options system based on (hard coded) inter-flag relationships.
/// Some of these relationships are set in the basic::options::process() function, some
/// of them are handled in this .cc file.
void
init_complex_options()
{
	process();

	// Set option system global
	basic::options::option.set_show_accessed_options_flag( option[ out::show_accessed_options ].value() );
	basic::options::option.set_show_unused_options_flag( option[ out::show_unused_options ].value() );

	// Immediate stop if requested by options (used for testing purposes).
	if ( option[ testing::HCF ]() ) {
		user_fixable_issue_exit("Option -HCF set. Do not pass go. Do not collect $200.");
	}
}

void
init_tracers(){

#ifdef USEMPI
	if( option[ out::mpi_tracer_to_file ].user() ){
		int mpi_rank( 0 );
		MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank );

		std::stringstream outfilename;
		outfilename << option[ out::mpi_tracer_to_file ]() << "_" << mpi_rank;
		basic::otstreamOP redirect_tracer( new basic::TracerToFile( outfilename.str() ));
		basic::TracerImpl::set_ios_hook( redirect_tracer, basic::Tracer::get_all_channels_string(), false );
		basic::TracerImpl::super_mute( true );
		utility::CSI_Sequence::suppress_CSI_codes(); // We're redirecting output to a file - don't use CSI terminal codes
	}
#endif

	// set Tracer options
	basic::TracerOptions TO;

	if ( option[ out::mute ].active() )   TO.muted = option[ out::mute ]();

	if ( option[ out::unmute ].active() ) TO.unmuted = option[ out::unmute ]();
	if ( option[ out::level  ].active() ) TO.level   = option[ out::level ]();
	if ( option[ out::levels ].active() ) TO.levels  = option[ out::levels ]();
	if ( option[ out::chname ].active() ) TO.print_channel_name = option[ out::chname ]();
	if ( option[ out::chtimestamp ].active() ) TO.timestamp = option[ out::chtimestamp ]();

	if ( option[ out::no_color ]() ) utility::CSI_Sequence::suppress_CSI_codes();

	basic::TracerImpl::set_tracer_options( TO );

	// Initialize tracers for external libraries
	core::chemical::rdkit::initialize_rdkit_tracers();
	// Ifdef USBCL contained within initialization function
	core::chemical::bcl::initialize_bcl_tracers();
}


void
init_source_revision(){
	if ( option[ run::version ]() ) {
		TR << "Rosetta version information" << std::endl;
		TR << "Package:                 " << utility::Version::package() << std::endl;
		TR << "Testing Server Revision: " << utility::Version::revision() << std::endl;
		TR << "PEP-440 Version:         " << utility::Version::version() << std::endl;
		TR << "Git Commit:              " << utility::Version::commit() << std::endl;
		TR << "URL:                     " << utility::Version::url() << std::endl;
		TR << "Date:                    " << utility::Version::date() << std::endl;
	} else {
		TR << "Rosetta version:";
		if ( utility::Version::package() != "devel" ) TR << " " << utility::Version::package();
		if ( utility::Version::revision() != "None" ) TR << " r" << utility::Version::revision();
		TR << " " << utility::Version::version() << " " << utility::Version::commit() << " " << utility::Version::url() << " " << utility::Version::date() << std::endl;
		TR << "Rosetta extras: " << utility::Version::extras() << std::endl;
	}
}

void
init_paths(){
	if ( option[ in::path::path ].user() ) {
		utility::io::izstream::set_alternative_search_paths(
			option[ in::path::path ]());
	}
}

/// @detail If you deprecate a long standing flag, add a line to this function to describe what the deprecated flag has been replaced with.
/// If the user specifies one of these flags, Rosetta will utility exit with a helpful message directing the user towards the new functionality
void check_deprecated_flags(){

	utility::vector1<std::string> error_messages;

	// Add deprecated flags and corresponding helpful error messages here.
	// This is the only thing you need to do to deprecate a flag.
	if ( option[LoopModel::input_pdb].user() ) {
		error_messages.push_back("-LoopModel:input_pdb is no longer used.  Please use -s to input pdb files.");
	}

	if ( option[antibody::numbering_scheme].user() ) {
		error_messages.push_back("-numbering_scheme option no longer used.  Please use -input_ab_scheme instead. ");
	}


	if ( option[basic::options::OptionKeys::out::file::no_scores_in_pdb].user() ) {
		error_messages.push_back("-no_scores_in_pdb is no longer used.  Please use -output_pose_energies_table false");
	}

	// deprecated flags from antibody/snugdock code
	if ( option[antibody::constrain_cter].user() ) {
		error_messages.push_back("-constrain_cter has been renamed to -h3_loop_csts_lr");
	}

	if ( option[antibody::all_atom_mode_kink_constraint].user() ) {
		error_messages.push_back("-all_atom_mode_kink_constraint has been renamed to -h3_loop_csts_hr");
	}

	if ( option[antibody::auto_generate_kink_constraint].user() ) {
		error_messages.push_back("-auto_generate_kink_constraint has been renamed to -auto_generate_h3_kink_constraint");
	}

	if ( error_messages.size() > 0 ) {
		utility::vector1<std::string>::const_iterator error_it;
		for ( error_it = error_messages.begin(); error_it != error_messages.end(); ++error_it ) {
			TR.Fatal << "You have specified one or more deprecated flags:" <<std::endl;
			TR.Fatal << *error_it <<std::endl;
		}
		std::exit(1);
	}
}

void
report_application_command(int argc, char * argv []){
	TR << "command:";
	for ( int i=0; i< argc; ++i ) {
		TR << ' ' <<  argv[i];
	}
	TR << std::endl;
}

void
random_delay(){
	if ( !option[ run::nodelay ]() ) {
		// no silly waiting in DEBUG or BOINC builds
		// Test inside if statement so that nodelay option gets touched even in debug mode.
#ifdef NDEBUG
#ifndef BOINC
		if( option[ run::delay ]() > 0 ) {
			int waittime = option[ run::delay ]();
			TR << "Delaying start of mini for " << waittime << " seconds due to -delay option" << std::endl;
			utility::sys_sleep( waittime );
		} else
		if( option[ run::random_delay ]() > 0 ) {
			int waittime = (int) ( (Real)option[ run::random_delay ]() * numeric::random::uniform() );
			TR << "Delaying of mini for " << waittime << " seconds (maximum = "
			   <<  option[ run::random_delay ]()
				 << " )" << std::endl
				 << "This prevents extreme IO levels when multiple jobs start simultaneously on" << std::endl
				 << "large computer clusters  and is default now. To prevent this add the option -nodelay" << std::endl
				 << "To change the random wait time use -run::random_delay <int> " << std::endl;
			utility::sys_sleep( waittime );
		}
#endif
#endif
	}
}

void
locate_rosetta_database(){

#ifndef __native_client__
	if ( !option[ in::path::database ].user() ) {
		std::string database_path;
		char * descr = getenv("ROSETTA3_DB");
		if ( descr ) {
			TR << "found database environment variable ROSETTA3_DB: "<< descr << std::endl;
			database_path = std::string( descr );
		} else {

			char path[1024];
			uint32_t path_size = sizeof( path );
			std::string path_string;


#if defined(MAC) || defined(__APPLE__)  ||  defined(__OSX__)
			uint32_t result = _NSGetExecutablePath(path, &path_size );

			// _NSGetExecutablePath returns 0 if the path was successfully copied.
			// Copies a null-terminated string
			if ( result == 0 ) {
				path_string = std::string( path );
			}
#endif

#if defined(linux) || defined(__linux__) || defined(__linux)
			path_size = readlink( "/proc/self/exe", path, path_size ); // This works on plain Linux, but not FreeBSD or Solaris.

			// readlink returns -1 on error, otherwise number of bytes written.
			// Does not append null to string
			if ( path_size > 0 ) {
				path_string = std::string( path, path_size );
			}
#endif

#ifdef WIN32
			TR << "There is some way to automatically figure out the path to rosetta_database with GetModuleFileName in Windows. This is already set up for linux and mac, its probably just one line to change in core/init/init.cc" << std::endl;
			// I think its this -- can someone comment out and compile on Windows? -- rhiju
			// GetModuleFileName( NULL, path, path_size );
#endif

			TR << "Resolved executable path: " << path_string << std::endl;

			// Attempt to resolve from source/../database if 'source' is in path
			// or ../database if not.

			// Logic must be updated for windows paths if windows executable resolution is added.
			if ( path_string.length() > 0 ) {
				Size found = std::string::npos;

				if ( (found = path_string.find("source/")) && (found != std::string::npos) ) {
					std::string rosetta_exe_dir = path_string.substr(0,found);
					database_path = rosetta_exe_dir + "database/";
					TR << "Looking for database based on location of executable: " << database_path << std::endl;
				} else if ( (found = path_string.rfind("/")) && (found != std::string::npos) ) {
					std::string rosetta_exe_dir = path_string.substr(0,found);
					database_path = rosetta_exe_dir + "/../database/";
					TR << "Looking for database based on location of executable: " << database_path << std::endl;
				}
			} else {
				TR << "Could not determine location of executable." << std::endl;
			}
		}

		if ( database_path.size() > 0 ) {
			option[ in::path::database ].value( database_path );
		} else {
			TR << "Could not find database. Either specify -database or set environment variable ROSETTA3_DB." << std::endl;
		}
	}
#endif
}

// Check if common config file exists and load it/them as flags files.
//  These files can be in the home directory as .rosetta/flags or the working directory.
void
check_load_fconfig(){

	if ( option ( basic::options::OptionKeys::in::no_fconfig).value() ) {
		return;
	}
	std::string const homedir( utility::file::get_home_dir() );

	utility::vector1< std::string > flags = option( basic::options::OptionKeys::in::fconfig ).value();
	std::string const db_path( homedir+"/.rosetta/flags/" );

	TR << "Checking for fconfig files in pwd and ./rosetta/flags " << std::endl;

	for ( auto & flag : flags ) {
		TR.Debug <<"Checking for -" << flag << "- fconfig " << std::endl;
		std::string final_file_path;

		//Check if flags file exists, get path. Either in working dir or database and load.
		std::string full_db_path = db_path + flag;
		std::string cid;

		//Load options from fconfig file into options collections if present
		// Local overrides database.
		if ( utility::file::file_exists( flag ) ) {
			final_file_path = flag;
		} else if (  homedir != ""  && utility::file::file_exists( full_db_path ) ) {
			final_file_path = full_db_path;
		} else {
			continue;
		}

		try{
			TR << "Reading fconfig..." << final_file_path << std::endl;
			utility::io::izstream config_file( final_file_path );
			option.load_options_from_stream( config_file, final_file_path, cid, true /*print*/);
			TR << std::endl << std::endl;
		} catch ( utility::excn::Exception& excn ) {
			throw( CREATE_EXCEPTION(utility::excn::Exception,  "ERROR: " + excn.msg() ) );
		}
	}
}


void
init_profiling(){
	basic::prof_reset(); //reads option run::profile -- starts clock TOTAL
}

void
init_resources() {
#ifdef PROCESS_STACK_SIZE
	// set stack size of process at runtime
	struct rlimit rl;
	int error = getrlimit( RLIMIT_STACK, &rl );
	if ( !error && rl.rlim_cur < PROCESS_STACK_SIZE ) {
		rl.rlim_cur = PROCESS_STACK_SIZE;
		setrlimit( RLIMIT_STACK, &rl );
	}
#endif
}

/// @brief Init basic core systems: options system, random system.
void init(int argc, char * argv [])
{
	basic::init();

	core::id::initialize_core_id_globals();

	//Initialize MPI
	init_mpi(argc, argv);

	//The options system manages command line options
	init_options(argc, argv);

	//Setup the crash reporter system (if we're using it)
	// We do this after the option initialization so we can output option information
	init_crash_reporter(argc, argv);

	//Tracers control output to std::cout and std::cerr
	init_tracers();

#ifndef PYROSETTA
	// We want to print this in pretty much all cases (though silence it if we're completely muted.
	if ( TR.Error.visible() ) {
		std::cout << "********  (C) Copyright Rosetta Commons Member Institutions.  ***************" << std::endl;
		std::cout << "* Use of Rosetta for commercial purposes may require purchase of a license. *" << std::endl;
		std::cout << "********  See LICENSE.md or email license@uw.edu for more details. **********" << std::endl;
	}
#endif

	//Read flag config file (common options/custom setup)
	check_load_fconfig();

	// Invoke basic::options::process() which holds a set of complex logic
	// for option system modifications; this function requires that the
	// tracers first be initialized.
	init_complex_options();

	//Initialize the latest and greatest score function parameters
	init_score_function_corrections( basic::options::option );

	//Choose to output source version control information?
	init_source_revision();

	//Setup basic search paths
	init_paths();

	//Check for deprecated flags specified by the user and output error messages if necessary
	check_deprecated_flags();

	//Describe the application execution command
	report_application_command(argc, argv);

	//Initalize random number generators
	int rand_seed = basic::random::init_random_number_generators();
	core::chemical::rdkit::initialize_rdkit_random(rand_seed);


	//Choose to randomly delay execution to desyncronize parallel execution
	random_delay();

	//Locate rosetta_database
	locate_rosetta_database();

	// Initialize additional BCL external library settings
#ifdef USEBCL
	core::chemical::bcl::require_bcl();
	core::chemical::bcl::locate_bcl();
	core::chemical::bcl::initialize_bcl_main(); // beneficial to occur after locate_rosetta_database()
	core::chemical::bcl::initialize_bcl_random(rand_seed);
#endif

	//Profiling measures execution performance
	init_profiling();

	//Set up system resources
	init_resources();

	// help out user...
	if  ( argc == 1 )  TR << std::endl << "USEFUL TIP: Type -help to get the options for this Rosetta executable." << std::endl << std::endl;

}


/// @brief wrapper for core system Init
void init( utility::vector1<std::string> const & args )
{
	// create arguments in argc/argv format
	int argc = args.size();
	auto **argv = new char*[ argc ];
	for ( int ii = 0; ii < argc; ++ii ) {
		argv[ ii ] = new char[ args[ii+1].size()+1 ];
		strncpy( argv[ii], args[ii+1].c_str(), args[ii+1].size() );
		argv[ ii ][ args[ii+1].size() ] = 0; // ensure null termination
	}

	// call init
	init( argc, argv );

	// delete freestore
	for ( int ii = 0; ii < argc; ++ii ) {
		delete[] argv[ ii ];
	}
	delete[] argv;
}

} // namespace init
} // namespace core
