// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/protocols/jd3/StandardJobQueen.cxxtest.hh
/// @brief  test suite for the StandardJobQueen
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// basic headers
#include <basic/options/option.hh>

// Utility headers
#include <utility/string_util.hh>

// C++ headers
#include <sstream>

// Mover creator headers
#include <protocols/aa_composition/AddCompositionConstraintMoverCreator.hh>
#include <protocols/aa_composition/ClearCompositionConstraintsMoverCreator.hh>
#include <protocols/abinitio/DomainAssemblyCreator.hh>
#include <protocols/analysis/InterfaceAnalyzerMoverCreator.hh>
#include <protocols/analysis/LoopAnalyzerMoverCreator.hh>
#include <protocols/antibody/AntibodyCDRGrafterCreator.hh>
#include <protocols/antibody/AntibodyNumberingConverterMoverCreator.hh>
#include <protocols/antibody/constraints/CDRDihedralConstraintMoverCreator.hh>
#include <protocols/antibody/constraints/ParatopeEpitopeSiteConstraintMoverCreator.hh>
#include <protocols/antibody/constraints/ParatopeSiteConstraintMoverCreator.hh>
#include <protocols/antibody/design/AntibodyDesignMoverCreator.hh>
#include <protocols/antibody/design/AntibodyDesignProtocolCreator.hh>
#include <protocols/backrub/BackrubMoverCreator.hh>
#include <protocols/backrub/BackrubProtocolCreator.hh>
#include <protocols/backrub/BackrubSidechainMoverCreator.hh>
#include <protocols/canonical_sampling/MetricRecorderCreator.hh>
#include <protocols/canonical_sampling/MetropolisHastingsMoverCreator.hh>
#include <protocols/canonical_sampling/ParallelTemperingCreator.hh>
#include <protocols/canonical_sampling/HamiltonianExchangeCreator.hh>
#include <protocols/canonical_sampling/PDBTrajectoryRecorderCreator.hh>
#include <protocols/canonical_sampling/SilentTrajectoryRecorderCreator.hh>
#include <protocols/canonical_sampling/SimulatedTemperingCreator.hh>
#include <protocols/canonical_sampling/TrialCounterObserverCreator.hh>
#include <protocols/constraint_generator/AddConstraintsCreator.hh>
#include <protocols/constraint_generator/RemoveConstraintsCreator.hh>
#include <protocols/carbohydrates/LinkageConformerMoverCreator.hh>
#include <protocols/carbohydrates/GlycanRelaxMoverCreator.hh>
#include <protocols/carbohydrates/SimpleGlycosylateMoverCreator.hh>
#include <protocols/cryst/cryst_movers_creator.hh>
#include <protocols/cryst/refinable_lattice_creator.hh>
#include <protocols/comparative_modeling/LoopRelaxMoverCreator.hh>
#include <protocols/hybridization/HybridizeProtocolCreator.hh>
#include <protocols/hybridization/HybridizeProtocolCreator.hh>
#include <protocols/hybridization/BackboneTorsionSamplerCreator.hh>
#include <protocols/hybridization/BackboneTorsionPerturbationCreator.hh>
#include <protocols/cyclic_peptide/CreateAngleConstraintCreator.hh>
#include <protocols/cyclic_peptide/CreateDistanceConstraintCreator.hh>
#include <protocols/cyclic_peptide/CreateTorsionConstraintCreator.hh>
#include <protocols/cyclic_peptide/DeclareBondCreator.hh>
#include <protocols/cyclic_peptide/PeptideStubMoverCreator.hh>
#include <protocols/cyclic_peptide/TryDisulfPermutationsCreator.hh>
#include <protocols/cyclic_peptide/FlipChiralityMoverCreator.hh>
#include <protocols/contact_map/ContactMapCreator.hh>
#include <protocols/denovo_design/movers/DisulfidizeMoverCreator.hh>
#include <protocols/denovo_design/movers/AddSegmentDataMoverCreator.hh>
#include <protocols/denovo_design/movers/AlignResiduesMoverCreator.hh>
#include <protocols/denovo_design/movers/BridgeChainsMoverCreator.hh>
#include <protocols/denovo_design/movers/BuildDeNovoBackboneMoverCreator.hh>
#include <protocols/denovo_design/movers/ExtendChainMoverCreator.hh>
#include <protocols/denovo_design/movers/FastDesignCreator.hh>
#include <protocols/denovo_design/movers/SetResidueAliasMoverCreator.hh>
#include <protocols/denovo_design/movers/MakeAsymmetricStructureDataMoverCreator.hh>
#include <protocols/design_opt/GreedyOptMutationMoverCreator.hh>
#include <protocols/dna/DesignProteinBackboneAroundDNACreator.hh>
#include <protocols/dna/DnaInterfaceMinMoverCreator.hh>
#include <protocols/dna/DnaInterfaceMultiStateDesignCreator.hh>
#include <protocols/dna/DnaInterfacePackerCreator.hh>
#include <protocols/dna/SeparateDnaFromNonDnaCreator.hh>
#include <protocols/docking/ConformerSwitchMoverCreator.hh>
#include <protocols/docking/DockingProtocolCreator.hh>
#include <protocols/docking/DockSetupMoverCreator.hh>
#include <protocols/docking/DockingInitialPerturbationCreator.hh>
#include <protocols/generalized_kinematic_closure/GeneralizedKICCreator.hh>
#include <protocols/helical_bundle/BackboneGridSamplerCreator.hh>
#include <protocols/helical_bundle/BundleGridSamplerCreator.hh>
#include <protocols/helical_bundle/FitSimpleHelixCreator.hh>
#include <protocols/helical_bundle/MakeBundleCreator.hh>
#include <protocols/helical_bundle/MakeBundleHelixCreator.hh>
#include <protocols/helical_bundle/PerturbBundleCreator.hh>
#include <protocols/helical_bundle/PerturbBundleHelixCreator.hh>
#include <protocols/ncbb/SecStructMinimizeMoverCreator.hh>
#include <protocols/ncbb/NcbbDockDesignProtocolCreator.hh>
#include <protocols/ncbb/oop/OopDockDesignProtocolCreator.hh>
#include <protocols/ncbb/oop/OopCreatorMoverCreator.hh>
#include <protocols/symmetric_docking/SymDockProtocolCreator.hh>
#include <protocols/symmetric_docking/SymFoldandDockCreators.hh>
#include <protocols/electron_density/SetupForDensityScoringMoverCreator.hh>
#include <protocols/electron_density/BfactorFittingMoverCreator.hh>
#include <protocols/electron_density/ScaleMapIntensitiesCreator.hh>
#include <protocols/electron_density/VoxelSpacingRefinementMoverCreator.hh>
#include <protocols/electron_density/ReportFSCCreator.hh>
#include <protocols/enzdes/AddOrRemoveMatchCstsCreator.hh>
#include <protocols/enzdes/BackboneSamplerCreator.hh>
#include <protocols/enzdes/EnzRepackMinimizeCreator.hh>
#include <protocols/enzdes/PackRotamersMoverPartGreedyCreator.hh>
#include <protocols/enzdes/EnzdesMoversCreator.hh>
#include <protocols/enzdes/EnzdesMoversCreator.hh>
#include <protocols/evolution/EvolutionaryDynamicsMoverCreator.hh>
#include <protocols/evolution/NucleotideMutationCreator.hh>
#include <protocols/rna/movers/ErraserMinimizerMoverCreator.hh>
#include <protocols/rna/movers/RNAIdealizeMoverCreator.hh>
#include <protocols/features/ReportToDBCreator.hh>
#include <protocols/features/TrajectoryReportToDBCreator.hh>
#include <protocols/fldsgn/BluePrintBDRCreator.hh>
#include <protocols/fldsgn/CircularPermutationCreator.hh>
#include <protocols/fldsgn/MatchResiduesMoverCreator.hh>
#include <protocols/fldsgn/SheetRemodelConstraintGeneratorCreator.hh>
#include <protocols/fldsgn/potentials/SetAACompositionPotentialCreator.hh>
#include <protocols/fldsgn/potentials/SetSecStructEnergiesCreator.hh>
#include <protocols/flexpep_docking/FlexPepDockingProtocolCreator.hh>
#include <protocols/flxbb/FlxbbDesignCreator.hh>
#include <protocols/flxbb/InterlockAromaCreator.hh>
#include <protocols/forge/constraints/InverseRotamersCstGeneratorCreator.hh>
#include <protocols/forge/constraints/InvrotTreeCstGeneratorCreator.hh>
#include <protocols/forge/constraints/NtoCConstraintGeneratorCreator.hh>
#include <protocols/forge/constraints/RemoveRemodelCstsCreator.hh>
#include <protocols/forge/remodel/ResidueVicinityCstGeneratorCreator.hh>
#include <protocols/forge/remodel/RemodelMoverCreator.hh>
#include <protocols/grafting/AnchoredGraftMoverCreator.hh>
#include <protocols/grafting/CCDEndsGraftMoverCreator.hh>
#include <protocols/grafting/simple_movers/DeleteRegionMoverCreator.hh>
#include <protocols/grafting/simple_movers/InsertPoseIntoPoseMoverCreator.hh>
#include <protocols/grafting/simple_movers/ReplaceRegionMoverCreator.hh>
#include <protocols/grafting/simple_movers/KeepRegionMoverCreator.hh>
#include <protocols/idealize/IdealizeMoverCreator.hh>
#include <protocols/hotspot_hashing/movers/PlaceSurfaceProbeCreator.hh>
#include <protocols/kinematic_closure/KicMoverCreator.hh>
#include <protocols/ligand_docking/AddHydrogensCreator.hh>
#include <protocols/ligand_docking/CompoundTranslateCreator.hh>
#include <protocols/ligand_docking/ComputeLigandRDFCreator.hh>
#include <protocols/ligand_docking/FinalMinimizerCreator.hh>
#include <protocols/ligand_docking/GrowLigandCreator.hh>
#include <protocols/ligand_docking/HighResDockerCreator.hh>
#include <protocols/ligand_docking/InterfaceScoreCalculatorCreator.hh>
#include <protocols/ligand_docking/LigandDesignCreator.hh>
#include <protocols/ligand_docking/MinimizeBackboneCreator.hh>
#include <protocols/ligand_docking/RandomConformersCreator.hh>
#include <protocols/ligand_docking/RotateCreator.hh>
#include <protocols/ligand_docking/RotatesCreator.hh>
#include <protocols/ligand_docking/SlideTogetherCreator.hh>
#include <protocols/ligand_docking/StartFromCreator.hh>
#include <protocols/ligand_docking/TransformCreator.hh>
#include <protocols/ligand_docking/TranslateCreator.hh>
#include <protocols/ligand_docking/WriteLigandMolFileCreator.hh>
#include <protocols/loop_build/LoopmodelWrapperCreator.hh>
#include <protocols/loop_build/LoopMover_SlidingWindowCreator.hh>
#include <protocols/loop_modeler/LoopModelerCreator.hh>
#include <protocols/loop_modeling/LoopProtocolCreator.hh>
#include <protocols/loop_modeling/LoopBuilderCreator.hh>
#include <protocols/loop_modeling/refiners/MinimizationRefinerCreator.hh>
#include <protocols/loop_modeling/refiners/RepackingRefinerCreator.hh>
#include <protocols/loop_modeling/refiners/RotamerTrialsRefinerCreator.hh>
#include <protocols/loop_modeling/samplers/LegacyKicSamplerCreator.hh>
#include <protocols/loop_modeling/utilities/PrepareForCentroidCreator.hh>
#include <protocols/loop_modeling/utilities/PrepareForFullatomCreator.hh>
#include <protocols/loophash/LoopHashMoverWrapperCreator.hh>
#include <protocols/loophash/LoopHashDiversifierCreator.hh>
#include <protocols/loops/FoldTreeFromLoopsWrapperCreator.hh>
#include <protocols/loops/loop_closure/ccd/CCDLoopClosureMoverCreator.hh>
#include <protocols/loops/loop_mover/perturb/LoopMover_CCDCreator.hh>
#include <protocols/loops/loop_mover/perturb/LoopMover_KICCreator.hh>
#include <protocols/loops/loop_mover/perturb/LoopMover_QuickCCDCreator.hh>
#include <protocols/loops/loop_mover/perturb/LoopMover_QuickCCD_MovesCreator.hh>
#include <protocols/loops/loop_mover/refine/LoopMover_BackrubCreator.hh>
#include <protocols/loops/loop_mover/refine/LoopMover_CCDCreator.hh>
#include <protocols/loops/loop_mover/refine/LoopMover_KICCreator.hh>
#include <protocols/loops/loop_mover/refine/LoopRefineInnerCycleContainerCreator.hh>
#include <protocols/loops/loop_mover/refine/RepackTrialCreator.hh>
#include <protocols/loops/loop_mover/refine/ShearMinCCDTrialCreator.hh>
#include <protocols/loops/loop_mover/refine/SmallMinCCDTrialCreator.hh>
#include <protocols/loops/loop_mover/LoopCMCreator.hh>
#include <protocols/match/MatcherMoverCreator.hh>
#include <protocols/matdes/ExtractSubposeMoverCreator.hh>
#include <protocols/matdes/MatDesGreedyOptMutationMoverCreator.hh>
#include <protocols/matdes/SymDofMoverCreator.hh>
#include <protocols/matdes/SchemePlaceMotifsMoverCreator.hh>
#include <protocols/md/CartesianMDCreator.hh>
#include <protocols/simple_moves/BBGaussianMoverCreator.hh>
#include <protocols/simple_moves/PeriodicBoxMoverCreator.hh>
#include <protocols/motifs/MotifDnaPackerCreator.hh>
#include <protocols/motif_grafting/movers/MotifGraftCreator.hh>
#include <protocols/moves/DsspMoverCreator.hh>
#include <protocols/moves/IfMoverCreator.hh>
#include <protocols/moves/MoverContainerCreator.hh>
#include <protocols/moves/IteratedConvergenceMoverCreator.hh>
#include <protocols/moves/PyMOLMoverCreator.hh>
#include <protocols/moves/RampingMoverCreator.hh>
#include <protocols/moves/FilterReportAsPoseExtraScoresMoverCreator.hh>
#include <protocols/monte_carlo/MonteCarloResetCreator.hh>
#include <protocols/simple_moves/AlignChainMoverCreator.hh>
#include <protocols/monte_carlo/ResetBaselineMoverCreator.hh>
#include <protocols/simple_moves/RingConformationMoverCreator.hh>
#include <protocols/simple_moves/SimpleThreadingMoverCreator.hh>
#include <protocols/nonlocal/SingleFragmentMoverCreator.hh>
#include <protocols/normalmode/NormalModeRelaxMoverCreator.hh>
#include <protocols/normalmode/NormalModeMinimizerCreator.hh>
#include <protocols/pb_potential/SetupPoissonBoltzmannPotentialCreator.hh>
#include <protocols/pose_length_moves/FixAllLoopsMoverCreator.hh>
#include <protocols/pose_length_moves/InsertResMoverCreator.hh>
#include <protocols/pose_length_moves/NearNativeLoopCloserCreator.hh>
#include <protocols/protein_interface_design/movers/AddChainBreakCreator.hh>
#include <protocols/protein_interface_design/movers/AddSidechainConstraintsToHotspotsCreator.hh>
#include <protocols/protein_interface_design/movers/BackrubDDMoverCreator.hh>
#include <protocols/protein_interface_design/movers/BestHotspotCstMoverCreator.hh>
#include <protocols/protein_interface_design/movers/BuildAlaPoseCreator.hh>
#include <protocols/protein_interface_design/movers/DesignMinimizeHbondsCreator.hh>
#include <protocols/protein_interface_design/movers/DisulfideMoverCreator.hh>
#include <protocols/protein_interface_design/movers/DockAndRetrieveSidechainsCreator.hh>
#include <protocols/simple_moves/DumpPdbCreator.hh>
#include <protocols/protein_interface_design/movers/FavorNativeResiduePreCycleCreator.hh>
#include <protocols/protein_interface_design/movers/FavorNonNativeResiduePreCycleCreator.hh>
#include <protocols/protein_interface_design/movers/HotspotDisjointedFoldTreeMoverCreator.hh>
#include <protocols/protein_interface_design/movers/HotspotHasherMoverCreator.hh>
#include <protocols/protein_interface_design/movers/InterfaceRecapitulationMoverCreator.hh>
#include <protocols/protein_interface_design/movers/LoopFinderCreator.hh>
#include <protocols/protein_interface_design/movers/LoopLengthChangeCreator.hh>
#include <protocols/protein_interface_design/movers/LoopMoverFromCommandLineCreator.hh>
#include <protocols/protein_interface_design/movers/LoopOverCreator.hh>
#include <protocols/protein_interface_design/movers/LoopRemodelCreator.hh>
#include <protocols/protein_interface_design/movers/MapHotspotCreator.hh>
#include <protocols/protein_interface_design/movers/PatchdockTransformCreator.hh>
#include <protocols/protein_interface_design/movers/PeptideStapleDesignMoverCreator.hh>
#include <protocols/protein_interface_design/movers/PlacementAuctionMoverCreator.hh>
#include <protocols/protein_interface_design/movers/PlacementMinimizationMoverCreator.hh>
#include <protocols/protein_interface_design/movers/PlaceOnLoopCreator.hh>
#include <protocols/protein_interface_design/movers/PlaceSimultaneouslyMoverCreator.hh>
#include <protocols/protein_interface_design/movers/PlaceStubMoverCreator.hh>
#include <protocols/protein_interface_design/movers/PrepackMoverCreator.hh>
#include <protocols/protein_interface_design/movers/ProteinInterfaceMultiStateDesignMoverCreator.hh>
#include <protocols/protein_interface_design/movers/RandomMutationCreator.hh>
#include <protocols/protein_interface_design/movers/RepackMinimizeCreator.hh>
#include <protocols/protein_interface_design/movers/SaveAndRetrieveSidechainsCreator.hh>
#include <protocols/protein_interface_design/movers/SetAtomTreeCreator.hh>
#include <protocols/protein_interface_design/movers/SetTemperatureFactorCreator.hh>
#include <protocols/protein_interface_design/movers/SetupHotspotConstraintsMoverCreator.hh>
#include <protocols/protein_interface_design/movers/SetupHotspotConstraintsLoopsMoverCreator.hh>
#include <protocols/protein_interface_design/movers/SpinMoverCreator.hh>
#include <protocols/protein_interface_design/movers/TaskAwareCstsCreator.hh>
#include <protocols/protein_interface_design/movers/SubroutineMoverCreator.hh>
#include <protocols/protein_interface_design/movers/TryRotamersCreator.hh>
#include <protocols/protein_interface_design/movers/ShoveResidueMoverCreator.hh>
#include <protocols/protein_interface_design/movers/VLBCreator.hh>
#include <protocols/protein_interface_design/movers/DockWithHotspotMoverCreator.hh>
#include <protocols/protein_interface_design/movers/TopologyBrokerMoverCreator.hh>
#include <protocols/qsar/RenderGridsToKinemageCreator.hh>
#include <protocols/rbsegment_relax/MakeStarTopologyCreator.hh>
#include <protocols/rbsegment_relax/OptimizeThreadingCreator.hh>
#include <protocols/rbsegment_relax/IdealizeHelicesCreator.hh>
#include <protocols/relax/AtomCoordinateCstMoverCreator.hh>
#include <protocols/relax/FastRelaxCreator.hh>
#include <protocols/relax/LocalRelaxCreator.hh>
#include <protocols/residue_selectors/StoreResidueSubsetMoverCreator.hh>
#include <protocols/rigid/RigidBodyMoverCreator.hh>
#include <protocols/rigid/RollMoverCreator.hh>
#include <protocols/rigid/UniformRigidBodyCMCreator.hh>
#include <protocols/rigid/UniformRigidBodyMoverCreator.hh>
#include <protocols/rosetta_scripts/ParsedProtocolCreator.hh>
#include <protocols/rosetta_scripts/SavePoseMoverCreator.hh>
#include <protocols/rosetta_scripts/MultiplePoseMoverCreator.hh>
#include <protocols/rosetta_scripts/MultipleOutputWrapperCreator.hh>
#include <protocols/rotamer_recovery/RotamerRecoveryMoverCreator.hh>
#include <protocols/seeded_abinitio/CAcstGeneratorCreator.hh>
#include <protocols/seeded_abinitio/CloseFoldCreator.hh>
#include <protocols/seeded_abinitio/CoordinateCstCreator.hh>
#include <protocols/seeded_abinitio/DefineMovableLoopsCreator.hh>
#include <protocols/seeded_abinitio/GrowPeptidesCreator.hh>
#include <protocols/seeded_abinitio/SeedFoldTreeCreator.hh>
#include <protocols/seeded_abinitio/SeedSetupMoverCreator.hh>
#include <protocols/seeded_abinitio/SwapSegmentCreator.hh>
#include <protocols/seeded_abinitio/SegmentHybridizerCreator.hh>
#include <protocols/constraint_movers/AddConstraintsToCurrentConformationMoverCreator.hh>
#include <protocols/simple_moves/AddChainMoverCreator.hh>
#include <protocols/simple_moves/AddJobPairDataCreator.hh>
#include <protocols/simple_moves/ChangeAndResetFoldTreeMoverCreator.hh>
#include <protocols/minimization_packing/BoltzmannRotamerMoverCreator.hh>
#include <protocols/constraint_movers/ClearConstraintsMoverCreator.hh>
#include <protocols/calc_taskop_movers/ConsensusDesignMoverCreator.hh>
#include <protocols/constraint_movers/ConstraintSetMoverCreator.hh>
#include <protocols/simple_moves/ContingentAcceptMoverCreator.hh>
#include <protocols/simple_moves/CoupledMoverCreator.hh>
#include <protocols/simple_moves/DeleteChainMoverCreator.hh>
#include <protocols/simple_moves/DeleteChainsMoverCreator.hh>
#include <protocols/simple_ddg/ddGCreator.hh>
#include <protocols/simple_moves/DisulfideInsertionMoverCreator.hh>
#include <protocols/pose_creation/ExtendedPoseMoverCreator.hh>
#include <protocols/simple_moves/FavorSequenceProfileCreator.hh>
#include <protocols/simple_moves/FavorSymmetricSequenceCreator.hh>
#include <protocols/simple_moves/FindConsensusSequenceCreator.hh>
#include <protocols/calc_taskop_movers/ForceDisulfidesMoverCreator.hh>
#include <protocols/monte_carlo/GenericMonteCarloMoverCreator.hh>
#include <protocols/pose_creation/LoadPDBMoverCreator.hh>
#include <protocols/simple_moves/LoadUnboundRotMoverCreator.hh>
#include <protocols/pose_creation/MakePolyXMoverCreator.hh>
#include <protocols/simple_moves/MembraneTopologyCreator.hh>
#include <protocols/pose_creation/MergePDBMoverCreator.hh>
#include <protocols/minimization_packing/MinMoverCreator.hh>
#include <protocols/minimization_packing/MinPackMoverCreator.hh>
#include <protocols/simple_moves/ModifyVariantTypeMoverCreator.hh>
#include <protocols/monte_carlo/MonteCarloRecoverCreator.hh>
#include <protocols/monte_carlo/MonteCarloTestCreator.hh>
#include <protocols/simple_moves/MSDMoverCreator.hh>
#include <protocols/simple_moves/MutateResidueCreator.hh>
#include <protocols/minimization_packing/PackRotamersMoverCreator.hh>
#include <protocols/simple_moves/PSSM2BfactorMoverCreator.hh>
#include <protocols/pose_creation/RepeatPropagationMoverCreator.hh>
#include <protocols/constraint_movers/ResidueTypeConstraintMoverCreator.hh>
#include <protocols/simple_moves/ReportEffectivePKACreator.hh>
#include <protocols/minimization_packing/RotamerTrialsMinMoverCreator.hh>
#include <protocols/minimization_packing/RotamerTrialsMoverCreator.hh>
#include <protocols/simple_moves/RandomTorsionMoverCreator.hh>
#include <protocols/simple_moves/RandomOmegaFlipMoverCreator.hh>
#include <protocols/minimization_packing/SaneMinMoverCreator.hh>
#include <protocols/simple_moves/ScoreMoverCreator.hh>
#include <protocols/simple_moves/SequenceProfileMoverCreator.hh>
#include <protocols/simple_moves/SetTorsionCreator.hh>
#include <protocols/simple_moves/BackboneMoverCreator.hh>
#include <protocols/simple_moves/ShortBackrubMoverCreator.hh>
#include <protocols/simple_moves/StorePoseSnapshotCreator.hh>
#include <protocols/simple_moves/BackboneMoverCreator.hh>
#include <protocols/simple_moves/StructProfileMoverCreator.hh>
#include <protocols/simple_moves/SuperimposeMoverCreator.hh>
#include <protocols/simple_moves/PDBReloadMoverCreator.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMoverCreator.hh>
#include <protocols/simple_moves/SwitchChainOrderMoverCreator.hh>
#include <protocols/minimization_packing/TaskAwareMinMoverCreator.hh>
#include <protocols/simple_moves/VirtualRootMoverCreator.hh>
#include <protocols/simple_moves/TumbleCreator.hh>
#include <protocols/simple_moves/bin_transitions/InitializeByBinsCreator.hh>
#include <protocols/simple_moves/bin_transitions/PerturbByBinsCreator.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMCMoverCreator.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMoverCreator.hh>
#include <protocols/simple_moves/sidechain_moves/SetChiMoverCreator.hh>
#include <protocols/simple_moves/sidechain_moves/JumpRotamerSidechainMoverCreator.hh>
#include <protocols/simple_moves/sidechain_moves/PerturbRotamerSidechainMoverCreator.hh>
#include <protocols/simple_moves/sidechain_moves/PerturbChiSidechainMoverCreator.hh>
#include <protocols/symmetry/SetupForSymmetryMoverCreator.hh>
#include <protocols/symmetry/DetectSymmetryMoverCreator.hh>
#include <protocols/symmetry/SetupNCSMoverCreator.hh>
#include <protocols/minimization_packing/symmetry/SymMinMoverCreator.hh>
#include <protocols/minimization_packing/symmetry/SymPackRotamersMoverCreator.hh>
#include <protocols/minimization_packing/symmetry/SymRotamerTrialsMoverCreator.hh>
#include <protocols/minimization_packing/symmetry/TaskAwareSymMinMoverCreator.hh>
#include <protocols/task_operations/StoreCombinedStoredTasksMoverCreator.hh>
#include <protocols/task_operations/StoreCompoundTaskMoverCreator.hh>
#include <protocols/task_operations/StoreTaskMoverCreator.hh>
#include <protocols/environment/EnvMoverCreator.hh>
#include <protocols/environment/CoMTrackerCMCreator.hh>
#include <protocols/environment/ScriptCMCreator.hh>
#include <protocols/abinitio/abscript/AbscriptMoverCreator.hh>
#include <protocols/abinitio/abscript/ConstraintPreparerCreator.hh>
#include <protocols/abinitio/abscript/FragmentCMCreator.hh>
#include <protocols/abinitio/abscript/FragmentJumpCMCreator.hh>
#include <protocols/abinitio/abscript/AbscriptLoopCloserCMCreator.hh>
#include <protocols/abinitio/abscript/StructPerturberCMCreator.hh>
#include <protocols/abinitio/abscript/RigidChunkCMCreator.hh>
#include <protocols/membrane/AddMembraneMoverCreator.hh>
#include <protocols/membrane/FlipMoverCreator.hh>
#include <protocols/membrane/MembranePositionFromTopologyMoverCreator.hh>
#include <protocols/membrane/SetMembranePositionMoverCreator.hh>
#include <protocols/membrane/TransformIntoMembraneMoverCreator.hh>
#include <protocols/membrane/symmetry/SymmetricAddMembraneMoverCreator.hh>
#include <protocols/membrane/benchmark/SampleTiltAnglesCreator.hh>
#include <protocols/membrane/benchmark/MakeCanonicalHelixCreator.hh>
#include <protocols/membrane/visualize/VisualizeMembraneMoverCreator.hh>
#include <protocols/membrane/visualize/VisualizeEmbeddingMoverCreator.hh>
#include <protocols/docking/membrane/MPDockingMoverCreator.hh>
#include <protocols/docking/membrane/MPDockingSetupMoverCreator.hh>
#include <protocols/legacy_sewing/sampling/LegacyAppendAssemblyMoverCreator.hh>
#include <protocols/legacy_sewing/sampling/LegacyGivenPathAssemblyMoverCreator.hh>
#include <protocols/legacy_sewing/sampling/LegacyGreedyAssemblyMoverCreator.hh>
#include <protocols/legacy_sewing/sampling/LegacyMonteCarloAssemblyMoverCreator.hh>
#include <protocols/legacy_sewing/sampling/LegacyAssemblyConstraintsMoverCreator.hh>
#include <protocols/legacy_sewing/sampling/LegacyAddStartnodeFragmentsCreator.hh>
#include <protocols/legacy_sewing/sampling/LegacyRepeatAssemblyMoverCreator.hh>
#include <protocols/legacy_sewing/sampling/LegacyEnumerateAssemblyMoverCreator.hh>
#include <protocols/symmetric_docking/membrane/MPSymDockMoverCreator.hh>
#include <protocols/membrane/AddMPLigandMoverCreator.hh>
#include <protocols/relax/membrane/MPFastRelaxMoverCreator.hh>
#include <protocols/vardist_solaccess/LoadVarSolDistSasaCalculatorMover.hh>


class BackwardsProtocolsMoverCreatorTests : public CxxTest::TestSuite
{
public:

	void test_protocols_aa_composition_AddCompositionConstraintMoverCreator_name()
	{ protocols::aa_composition::AddCompositionConstraintMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "AddCompositionConstraintMover" ); }

	void test_protocols_aa_composition_ClearCompositionConstraintsMoverCreator_name()
	{ protocols::aa_composition::ClearCompositionConstraintsMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ClearCompositionConstraintsMover" ); }

	void test_protocols_abinitio_DomainAssemblyCreator_name()
	{ protocols::abinitio::DomainAssemblyCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "DomainAssembly" ); }

	void test_protocols_analysis_InterfaceAnalyzerMoverCreator_name()
	{ protocols::analysis::InterfaceAnalyzerMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "InterfaceAnalyzerMover" ); }

	void test_protocols_analysis_LoopAnalyzerMoverCreator_name()
	{ protocols::analysis::LoopAnalyzerMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "LoopAnalyzerMover" ); }

	void test_protocols_antibody_AntibodyCDRGrafterCreator_name()
	{ protocols::antibody::AntibodyCDRGrafterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "AntibodyCDRGrafter" ); }

	void test_protocols_antibody_AntibodyNumberingConverterMoverCreator_name()
	{ protocols::antibody::AntibodyNumberingConverterMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "AntibodyNumberingConverterMover" ); }

	void test_protocols_antibody_constraints_CDRDihedralConstraintMoverCreator_name()
	{ protocols::antibody::constraints::CDRDihedralConstraintMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "CDRDihedralConstraintMover" ); }

	void test_protocols_antibody_constraints_ParatopeEpitopeSiteConstraintMoverCreator_name()
	{ protocols::antibody::constraints::ParatopeEpitopeSiteConstraintMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ParatopeEpitopeConstraintMover" ); }

	void test_protocols_antibody_constraints_ParatopeSiteConstraintMoverCreator_name()
	{ protocols::antibody::constraints::ParatopeSiteConstraintMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ParatopeSiteConstraintMover" ); }

	void test_protocols_antibody_design_AntibodyDesignMoverCreator_name()
	{ protocols::antibody::design::AntibodyDesignMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "AntibodyDesignMover" ); }

	void test_protocols_antibody_design_AntibodyDesignProtocolCreator_name()
	{ protocols::antibody::design::AntibodyDesignProtocolCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "AntibodyDesignProtocol" ); }

	void test_protocols_backrub_BackrubMoverCreator_name()
	{ protocols::backrub::BackrubMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "Backrub" ); }

	void test_protocols_backrub_BackrubProtocolCreator_name()
	{ protocols::backrub::BackrubProtocolCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "BackrubProtocol" ); }

	void test_protocols_backrub_BackrubSidechainMoverCreator_name()
	{ protocols::backrub::BackrubSidechainMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "BackrubSidechain" ); }

	void test_protocols_canonical_sampling_MetricRecorderCreator_name()
	{ protocols::canonical_sampling::MetricRecorderCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "MetricRecorder" ); }

	void test_protocols_canonical_sampling_MetropolisHastingsMoverCreator_name()
	{ protocols::canonical_sampling::MetropolisHastingsMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "MetropolisHastings" ); }

	void test_protocols_canonical_sampling_ParallelTemperingCreator_name()
	{ protocols::canonical_sampling::ParallelTemperingCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ParallelTempering" ); }

	void test_protocols_canonical_sampling_HamiltonianExchangeCreator_name()
	{ protocols::canonical_sampling::HamiltonianExchangeCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "HamiltonianExchange" ); }

	void test_protocols_canonical_sampling_PDBTrajectoryRecorderCreator_name()
	{ protocols::canonical_sampling::PDBTrajectoryRecorderCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "PDBTrajectoryRecorder" ); }

	void test_protocols_canonical_sampling_SilentTrajectoryRecorderCreator_name()
	{ protocols::canonical_sampling::SilentTrajectoryRecorderCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SilentTrajectoryRecorder" ); }

	void test_protocols_canonical_sampling_SimulatedTemperingCreator_name()
	{ protocols::canonical_sampling::SimulatedTemperingCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SimulatedTempering" ); }

	void test_protocols_canonical_sampling_TrialCounterObserverCreator_name()
	{ protocols::canonical_sampling::TrialCounterObserverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "TrialCounterObserver" ); }

	void test_protocols_constraint_generator_AddConstraintsCreator_name()
	{ protocols::constraint_generator::AddConstraintsCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "AddConstraints" ); }

	void test_protocols_constraint_generator_RemoveConstraintsCreator_name()
	{ protocols::constraint_generator::RemoveConstraintsCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "RemoveConstraints" ); }

	void test_protocols_carbohydrates_LinkageConformerMoverCreator_name()
	{ protocols::carbohydrates::LinkageConformerMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "LinkageConformerMover" ); }

	void test_protocols_carbohydrates_GlycanRelaxMoverCreator_name()
	{ protocols::carbohydrates::GlycanRelaxMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "GlycanRelaxMover" ); }

	void test_protocols_carbohydrates_SimpleGlycosylateMoverCreator_name()
	{ protocols::carbohydrates::SimpleGlycosylateMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SimpleGlycosylateMover" ); }

	void test_protocols_cryst_ReportGradientsMoverCreator_name()
	{ protocols::cryst::ReportGradientsMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ReportGradients" ); }

	void test_protocols_cryst_SetCrystWeightMoverCreator_name()
	{ protocols::cryst::SetCrystWeightMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SetCrystWeight" ); }

	void test_protocols_cryst_RecomputeDensityMapMoverCreator_name()
	{ protocols::cryst::RecomputeDensityMapMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "RecomputeDensityMap" ); }

	void test_protocols_cryst_LoadDensityMapMoverCreator_name()
	{ protocols::cryst::LoadDensityMapMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "LoadDensityMap" ); }

	void test_protocols_cryst_FitBfactorsMoverCreator_name()
	{ protocols::cryst::FitBfactorsMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "FitBfactors" ); }

	void test_protocols_cryst_UpdateSolventMoverCreator_name()
	{ protocols::cryst::UpdateSolventMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "UpdateSolvent" ); }

	void test_protocols_cryst_TagPoseWithRefinementStatsMoverCreator_name()
	{ protocols::cryst::TagPoseWithRefinementStatsMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "TagPoseWithRefinementStats" ); }

	void test_protocols_cryst_SetRefinementOptionsMoverCreator_name()
	{ protocols::cryst::SetRefinementOptionsMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SetRefinementOptions" ); }

	void test_protocols_cryst_UpdateCrystInfoCreator_name()
	{ protocols::cryst::UpdateCrystInfoCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "UpdateCrystInfo" ); }

	void test_protocols_cryst_DockLatticeMoverCreator_name()
	{ protocols::cryst::DockLatticeMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "DockLatticeMover" ); }

	void test_protocols_cryst_MakeLatticeMoverCreator_name()
	{ protocols::cryst::MakeLatticeMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "MakeLatticeMover" ); }

	void test_protocols_cryst_MakeLayerMoverCreator_name()
	{ protocols::cryst::MakeLayerMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "MakeLayerMover" ); }

	void test_protocols_comparative_modeling_LoopRelaxMoverCreator_name()
	{ protocols::comparative_modeling::LoopRelaxMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "LoopRelaxMover" ); }

	void test_protocols_hybridization_HybridizeProtocolCreator_name()
	{ protocols::hybridization::HybridizeProtocolCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "Hybridize" ); }

	void test_protocols_hybridization_CartesianSamplerCreator_name()
	{ protocols::hybridization::CartesianSamplerCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "CartesianSampler" ); }

	void test_protocols_hybridization_BackboneTorsionSamplerCreator_name()
	{ protocols::hybridization::BackboneTorsionSamplerCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "BackboneTorsionSampler" ); }

	void test_protocols_hybridization_BackboneTorsionPerturbationCreator_name()
	{ protocols::hybridization::BackboneTorsionPerturbationCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "BackboneTorsionPerturbation" ); }

	void test_protocols_cyclic_peptide_CreateAngleConstraintCreator_name()
	{ protocols::cyclic_peptide::CreateAngleConstraintCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "CreateAngleConstraint" ); }

	void test_protocols_cyclic_peptide_CreateDistanceConstraintCreator_name()
	{ protocols::cyclic_peptide::CreateDistanceConstraintCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "CreateDistanceConstraint" ); }

	void test_protocols_cyclic_peptide_CreateTorsionConstraintCreator_name()
	{ protocols::cyclic_peptide::CreateTorsionConstraintCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "CreateTorsionConstraint" ); }

	void test_protocols_cyclic_peptide_DeclareBondCreator_name()
	{ protocols::cyclic_peptide::DeclareBondCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "DeclareBond" ); }

	void test_protocols_cyclic_peptide_PeptideStubMoverCreator_name()
	{ protocols::cyclic_peptide::PeptideStubMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "PeptideStubMover" ); }

	void test_protocols_cyclic_peptide_TryDisulfPermutationsCreator_name()
	{ protocols::cyclic_peptide::TryDisulfPermutationsCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "TryDisulfPermutations" ); }

	void test_protocols_cyclic_peptide_FlipChiralityMoverCreator_name()
	{ protocols::cyclic_peptide::FlipChiralityMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "FlipChiralityMover" ); }

	void test_protocols_contact_map_ContactMapCreator_name()
	{ protocols::contact_map::ContactMapCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ContactMap" ); }

	void test_protocols_denovo_design_DisulfidizeMoverCreator_name()
	{ protocols::denovo_design::DisulfidizeMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "Disulfidize" ); }

	void test_protocols_denovo_design_movers_AddSegmentDataMoverCreator_name()
	{ protocols::denovo_design::movers::AddSegmentDataMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "AddSegmentDataMover" ); }

	void test_protocols_denovo_design_movers_AlignResiduesMoverCreator_name()
	{ protocols::denovo_design::movers::AlignResiduesMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "AlignResiduesMover" ); }

	//void dont_test_protocols_denovo_design_movers_BridgeChainsMoverCreator_name()
	//{ protocols::denovo_design::movers::BridgeChainsMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "BridgeChainsMover" ); }

	void test_protocols_denovo_design_movers_BridgeChainsCreator_name()
	{ protocols::denovo_design::movers::BridgeChainsCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "BridgeChains" ); }

	void test_protocols_denovo_design_movers_BuildDeNovoBackboneMoverCreator_name()
	{ protocols::denovo_design::movers::BuildDeNovoBackboneMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "BuildDeNovoBackboneMover" ); }

	void test_protocols_denovo_design_movers_ExtendChainMoverCreator_name()
	{ protocols::denovo_design::movers::ExtendChainMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ExtendChain" ); }

	void test_protocols_denovo_design_movers_FastDesignCreator_name()
	{ protocols::denovo_design::movers::FastDesignCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "FastDesign" ); }

	void test_protocols_denovo_design_movers_SetResidueAliasMoverCreator_name()
	{ protocols::denovo_design::movers::SetResidueAliasMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SetResidueAlias" ); }

	void test_protocols_denovo_design_movers_MakeAsymmetricStructureDataMoverCreator_name()
	{ protocols::denovo_design::movers::MakeAsymmetricStructureDataMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "MakeAsymmetricStructureDataMover" ); }

	void test_protocols_design_opt_GreedyOptMutationMoverCreator_name()
	{ protocols::design_opt::GreedyOptMutationMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "GreedyOptMutationMover" ); }

	void test_protocols_dna_DesignProteinBackboneAroundDNACreator_name()
	{ protocols::dna::DesignProteinBackboneAroundDNACreator cr; TS_ASSERT_EQUALS( cr.keyname(), "DesignProteinBackboneAroundDNA" ); }

	void test_protocols_dna_DnaInterfaceMinMoverCreator_name()
	{ protocols::dna::DnaInterfaceMinMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "DnaInterfaceMinMover" ); }

	void test_protocols_dna_DnaInterfaceMultiStateDesignCreator_name()
	{ protocols::dna::DnaInterfaceMultiStateDesignCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "DnaInterfaceMultiStateDesign" ); }

	void test_protocols_dna_DnaInterfacePackerCreator_name()
	{ protocols::dna::DnaInterfacePackerCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "DnaInterfacePacker" ); }

	void test_protocols_dna_SeparateDnaFromNonDnaCreator_name()
	{ protocols::dna::SeparateDnaFromNonDnaCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SeparateDnaFromNonDna" ); }

	void test_protocols_docking_ConformerSwitchMoverCreator_name()
	{ protocols::docking::ConformerSwitchMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ConformerSwitchMover" ); }

	void test_protocols_docking_DockingProtocolCreator_name()
	{ protocols::docking::DockingProtocolCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "DockingProtocol" ); }

	void test_protocols_docking_DockSetupMoverCreator_name()
	{ protocols::docking::DockSetupMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "DockSetupMover" ); }

	void test_protocols_docking_DockingInitialPerturbationCreator_name()
	{ protocols::docking::DockingInitialPerturbationCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "DockingInitialPerturbation" ); }

	void test_protocols_generalized_kinematic_closure_GeneralizedKICCreator_name()
	{ protocols::generalized_kinematic_closure::GeneralizedKICCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "GeneralizedKIC" ); }

	void test_protocols_helical_bundle_BackboneGridSamplerCreator_name()
	{ protocols::helical_bundle::BackboneGridSamplerCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "BackboneGridSampler" ); }

	void test_protocols_helical_bundle_BundleGridSamplerCreator_name()
	{ protocols::helical_bundle::BundleGridSamplerCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "BundleGridSampler" ); }

	void test_protocols_helical_bundle_FitSimpleHelixCreator_name()
	{ protocols::helical_bundle::FitSimpleHelixCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "FitSimpleHelix" ); }

	void test_protocols_helical_bundle_MakeBundleCreator_name()
	{ protocols::helical_bundle::MakeBundleCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "MakeBundle" ); }

	void test_protocols_helical_bundle_MakeBundleHelixCreator_name()
	{ protocols::helical_bundle::MakeBundleHelixCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "MakeBundleHelix" ); }

	void test_protocols_helical_bundle_PerturbBundleCreator_name()
	{ protocols::helical_bundle::PerturbBundleCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "PerturbBundle" ); }

	void test_protocols_helical_bundle_PerturbBundleHelixCreator_name()
	{ protocols::helical_bundle::PerturbBundleHelixCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "PerturbBundleHelix" ); }

	void test_protocols_ncbb_SecStructMinimizeMoverCreator_name()
	{ protocols::ncbb::SecStructMinimizeMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SecStructMinimizeMover" ); }

	void test_protocols_ncbb_NcbbDockDesignProtocolCreator_name()
	{ protocols::ncbb::NcbbDockDesignProtocolCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "NcbbDockDesign" ); }

	void test_protocols_ncbb_oop_OopDockDesignProtocolCreator_name()
	{ protocols::ncbb::oop::OopDockDesignProtocolCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "OopDockDesign" ); }

	void test_protocols_ncbb_oop_OopCreatorMoverCreator_name()
	{ protocols::ncbb::oop::OopCreatorMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "OopCreatorMover" ); }

	void test_protocols_symmetric_docking_SymDockProtocolCreator_name()
	{ protocols::symmetric_docking::SymDockProtocolCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SymDockProtocol" ); }

	void test_protocols_symmetric_docking_SymFoldandDockMoveRbJumpMoverCreator_name()
	{ protocols::symmetric_docking::SymFoldandDockMoveRbJumpMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SymFoldandDockMoveRbJumpMover" ); }

	void test_protocols_symmetric_docking_SymFoldandDockSlideTrialMoverCreator_name()
	{ protocols::symmetric_docking::SymFoldandDockSlideTrialMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SymFoldandDockSlideTrialMover" ); }

	void test_protocols_symmetric_docking_SymFoldandDockRbTrialMoverCreator_name()
	{ protocols::symmetric_docking::SymFoldandDockRbTrialMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SymFoldandDockRbTrialMover" ); }

	void test_protocols_electron_density_SetupForDensityScoringMoverCreator_name()
	{ protocols::electron_density::SetupForDensityScoringMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SetupForDensityScoring" ); }

	void test_protocols_electron_density_BfactorFittingMoverCreator_name()
	{ protocols::electron_density::BfactorFittingMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "BfactorFitting" ); }

	void test_protocols_electron_density_ScaleMapIntensitiesCreator_name()
	{ protocols::electron_density::ScaleMapIntensitiesCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ScaleMapIntensities" ); }

	void test_protocols_electron_density_VoxelSpacingRefinementMoverCreator_name()
	{ protocols::electron_density::VoxelSpacingRefinementMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "VoxelSpacingRefinement" ); }

	void test_protocols_electron_density_ReportFSCCreator_name()
	{ protocols::electron_density::ReportFSCCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ReportFSC" ); }

	void test_protocols_enzdes_AddOrRemoveMatchCstsCreator_name()
	{ protocols::enzdes::AddOrRemoveMatchCstsCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "AddOrRemoveMatchCsts" ); }

	void test_protocols_enzdes_BackboneSamplerCreator_name()
	{ protocols::enzdes::BackboneSamplerCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "BackboneSampler" ); }

	void test_protocols_enzdes_EnzRepackMinimizeCreator_name()
	{ protocols::enzdes::EnzRepackMinimizeCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "EnzRepackMinimize" ); }

	void test_protocols_enzdes_PackRotamersMoverPartGreedyCreator_name()
	{ protocols::enzdes::PackRotamersMoverPartGreedyCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "PackRotamersMoverPartGreedy" ); }

	void test_protocols_enzdes_PredesignPerturbMoverCreator_name()
	{ protocols::enzdes::PredesignPerturbMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "PredesignPerturbMover" ); }

	void test_protocols_enzdes_UpdateEnzdesHeaderMoverCreator_name()
	{ protocols::enzdes::UpdateEnzdesHeaderMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "UpdateEnzdesHeader" ); }

	void test_protocols_evolution_EvolutionaryDynamicsMoverCreator_name()
	{ protocols::evolution::EvolutionaryDynamicsMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "EvolutionaryDynamics" ); }

	void test_protocols_evolution_NucleotideMutationCreator_name()
	{ protocols::evolution::NucleotideMutationCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "NucleotideMutation" ); }

	void test_protocols_farna_ErraserMinimizerMoverCreator_name()
	{ protocols::rna::movers::ErraserMinimizerMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ErraserMinimizerMover" ); }

	void test_protocols_farna_RNAIdealizeMoverCreator_name()
	{ protocols::rna::movers::RNAIdealizeMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "RNAIdealizeMover" ); }

	void test_protocols_features_ReportToDBCreator_name()
	{ protocols::features::ReportToDBCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ReportToDB" ); }

	void test_protocols_features_TrajectoryReportToDBCreator_name()
	{ protocols::features::TrajectoryReportToDBCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "TrajectoryReportToDB" ); }

	void test_protocols_fldsgn_BluePrintBDRCreator_name()
	{ protocols::fldsgn::BluePrintBDRCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "BluePrintBDR" ); }

	void test_protocols_fldsgn_CircularPermutationCreator_name()
	{ protocols::fldsgn::CircularPermutationCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "CircularPermutation" ); }

	void test_protocols_fldsgn_MatchResiduesMoverCreator_name()
	{ protocols::fldsgn::MatchResiduesMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "MatchResiduesMover" ); }

	void test_protocols_fldsgn_SheetRemodelConstraintGeneratorCreator_name()
	{ protocols::fldsgn::SheetRemodelConstraintGeneratorCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SheetCstGenerator" ); }

	void test_protocols_fldsgn_potentials_SetAACompositionPotentialCreator_name()
	{ protocols::fldsgn::potentials::SetAACompositionPotentialCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SetAACompositionPotential" ); }

	void test_protocols_fldsgn_potentials_SetSecStructEnergiesCreator_name()
	{ protocols::fldsgn::potentials::SetSecStructEnergiesCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SetSecStructEnergies" ); }

	void test_protocols_flexpep_docking_FlexPepDockingProtocolCreator_name()
	{ protocols::flexpep_docking::FlexPepDockingProtocolCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "FlexPepDock" ); }

	void test_protocols_flxbb_FlxbbDesignCreator_name()
	{ protocols::flxbb::FlxbbDesignCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "FlxbbDesign" ); }

	void test_protocols_flxbb_InterlockAromaCreator_name()
	{ protocols::flxbb::InterlockAromaCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "InterlockAroma" ); }

	void test_protocols_forge_constraints_InverseRotamersCstGeneratorCreator_name()
	{ protocols::forge::constraints::InverseRotamersCstGeneratorCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "InverseRotamersCstGenerator" ); }

	void test_protocols_forge_constraints_InvrotTreeCstGeneratorCreator_name()
	{ protocols::forge::constraints::InvrotTreeCstGeneratorCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "InvrotTreeCstGenerator" ); }

	void test_protocols_forge_constraints_NtoCConstraintGeneratorCreator_name()
	{ protocols::forge::constraints::NtoCConstraintGeneratorCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "NtoCConstraintGenerator" ); }

	void test_protocols_forge_constraints_RemoveRemodelCstsCreator_name()
	{ protocols::forge::constraints::RemoveRemodelCstsCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "RemoveRemodelCsts" ); }

	void test_protocols_forge_remodel_ResidueVicinityCstGeneratorCreator_name()
	{ protocols::forge::remodel::ResidueVicinityCstGeneratorCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ResidueVicinityCstCreator" ); }

	void test_protocols_forge_remodel_RemodelMoverCreator_name()
	{ protocols::forge::remodel::RemodelMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "RemodelMover" ); }

	void test_protocols_grafting_AnchoredGraftMoverCreator_name()
	{ protocols::grafting::AnchoredGraftMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "AnchoredGraftMover" ); }

	void test_protocols_grafting_CCDEndsGraftMoverCreator_name()
	{ protocols::grafting::CCDEndsGraftMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "CCDEndsGraftMover" ); }

	void test_protocols_grafting_simple_movers_DeleteRegionMoverCreator_name()
	{ protocols::grafting::simple_movers::DeleteRegionMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "DeleteRegionMover" ); }

	void test_protocols_grafting_simple_movers_InsertPoseIntoPoseMoverCreator_name()
	{ protocols::grafting::simple_movers::InsertPoseIntoPoseMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "InsertPoseIntoPoseMover" ); }

	void test_protocols_grafting_simple_movers_ReplaceRegionMoverCreator_name()
	{ protocols::grafting::simple_movers::ReplaceRegionMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ReplaceRegionMover" ); }

	void test_protocols_grafting_simple_movers_KeepRegionMoverCreator_name()
	{ protocols::grafting::simple_movers::KeepRegionMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "KeepRegionMover" ); }

	void test_protocols_idealize_IdealizeMoverCreator_name()
	{ protocols::idealize::IdealizeMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "Idealize" ); }

	void test_protocols_hotspot_hashing_movers_PlaceSurfaceProbeCreator_name()
	{ protocols::hotspot_hashing::movers::PlaceSurfaceProbeCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "PlaceSurfaceProbe" ); }

	void test_protocols_kinematic_closure_KicMoverCreator_name()
	{ protocols::kinematic_closure::KicMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "KicMover" ); }

	void test_protocols_ligand_docking_AddHydrogensCreator_name()
	{ protocols::ligand_docking::AddHydrogensCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "AddHydrogens" ); }

	void test_protocols_ligand_docking_CompoundTranslateCreator_name()
	{ protocols::ligand_docking::CompoundTranslateCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "CompoundTranslate" ); }

	void test_protocols_ligand_docking_ComputeLigandRDFCreator_name()
	{ protocols::ligand_docking::ComputeLigandRDFCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ComputeLigandRDF" ); }

	void test_protocols_ligand_docking_FinalMinimizerCreator_name()
	{ protocols::ligand_docking::FinalMinimizerCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "FinalMinimizer" ); }

	void test_protocols_ligand_docking_GrowLigandCreator_name()
	{ protocols::ligand_docking::GrowLigandCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "GrowLigand" ); }

	void test_protocols_ligand_docking_HighResDockerCreator_name()
	{ protocols::ligand_docking::HighResDockerCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "HighResDocker" ); }

	void test_protocols_ligand_docking_InterfaceScoreCalculatorCreator_name()
	{ protocols::ligand_docking::InterfaceScoreCalculatorCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "InterfaceScoreCalculator" ); }

	void test_protocols_ligand_docking_LigandDesignCreator_name()
	{ protocols::ligand_docking::LigandDesignCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "LigandDesign" ); }

	void test_protocols_ligand_docking_MinimizeBackboneCreator_name()
	{ protocols::ligand_docking::MinimizeBackboneCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "MinimizeBackbone" ); }

	void test_protocols_ligand_docking_RandomConformersCreator_name()
	{ protocols::ligand_docking::RandomConformersCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "RandomConformers" ); }

	void test_protocols_ligand_docking_RotateCreator_name()
	{ protocols::ligand_docking::RotateCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "Rotate" ); }

	void test_protocols_ligand_docking_RotatesCreator_name()
	{ protocols::ligand_docking::RotatesCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "Rotates" ); }

	void test_protocols_ligand_docking_SlideTogetherCreator_name()
	{ protocols::ligand_docking::SlideTogetherCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SlideTogether" ); }

	void test_protocols_ligand_docking_StartFromCreator_name()
	{ protocols::ligand_docking::StartFromCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "StartFrom" ); }

	void test_protocols_ligand_docking_TransformCreator_name()
	{ protocols::ligand_docking::TransformCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "Transform" ); }

	void test_protocols_ligand_docking_TranslateCreator_name()
	{ protocols::ligand_docking::TranslateCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "Translate" ); }

	void test_protocols_ligand_docking_WriteLigandMolFileCreator_name()
	{ protocols::ligand_docking::WriteLigandMolFileCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "WriteLigandMolFile" ); }

	void test_protocols_loop_build_LoopmodelWrapperCreator_name()
	{ protocols::loop_build::LoopmodelWrapperCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "LoopmodelWrapper" ); }

	void test_protocols_loop_build_LoopMover_SlidingWindowCreator_name()
	{ protocols::loop_build::LoopMover_SlidingWindowCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "LoopMover_SlidingWindow" ); }

	void test_protocols_loop_modeler_LoopModelerCreator_name()
	{ protocols::loop_modeler::LoopModelerCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "LoopModeler" ); }

	void test_protocols_loop_modeling_LoopProtocolCreator_name()
	{ protocols::loop_modeling::LoopProtocolCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "LoopProtocol" ); }

	void test_protocols_loop_modeling_LoopBuilderCreator_name()
	{ protocols::loop_modeling::LoopBuilderCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "LoopBuilder" ); }

	void test_protocols_loop_modeling_refiners_MinimizationRefinerCreator_name()
	{ protocols::loop_modeling::refiners::MinimizationRefinerCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "MinimizationRefiner" ); }

	void test_protocols_loop_modeling_refiners_RepackingRefinerCreator_name()
	{ protocols::loop_modeling::refiners::RepackingRefinerCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "RepackingRefiner" ); }

	void test_protocols_loop_modeling_refiners_RotamerTrialsRefinerCreator_name()
	{ protocols::loop_modeling::refiners::RotamerTrialsRefinerCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "RotamerTrialsRefiner" ); }

	void test_protocols_loop_modeling_samplers_LegacyKicSamplerCreator_name()
	{ protocols::loop_modeling::samplers::LegacyKicSamplerCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "LegacyKicSampler" ); }

	void test_protocols_loop_modeling_utilities_PrepareForCentroidCreator_name()
	{ protocols::loop_modeling::utilities::PrepareForCentroidCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "PrepareForCentroid" ); }

	void test_protocols_loop_modeling_utilities_PrepareForFullatomCreator_name()
	{ protocols::loop_modeling::utilities::PrepareForFullatomCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "PrepareForFullatom" ); }

	void test_protocols_loophash_LoopHashMoverWrapperCreator_name()
	{ protocols::loophash::LoopHashMoverWrapperCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "LoopHash" ); }

	void test_protocols_loophash_LoopHashDiversifierCreator_name()
	{ protocols::loophash::LoopHashDiversifierCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "LoopHashDiversifier" ); }

	void test_protocols_loops_FoldTreeFromLoopsCreator_name()
	{ protocols::loops::FoldTreeFromLoopsCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "FoldTreeFromLoops" ); }

	void test_protocols_loops_loop_closure_ccd_CCDLoopClosureMoverCreator_name()
	{ protocols::loops::loop_closure::ccd::CCDLoopClosureMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "CCDLoopClosureMover" ); }

	void test_protocols_loops_loop_mover_perturb_LoopMover_Perturb_CCDCreator_name()
	{ protocols::loops::loop_mover::perturb::LoopMover_Perturb_CCDCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "LoopMover_Perturb_CCD" ); }

	void test_protocols_loops_loop_mover_perturb_LoopMover_Perturb_KICCreator_name()
	{ protocols::loops::loop_mover::perturb::LoopMover_Perturb_KICCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "LoopMover_Perturb_KIC" ); }

	void test_protocols_loops_loop_mover_perturb_LoopMover_Perturb_QuickCCDCreator_name()
	{ protocols::loops::loop_mover::perturb::LoopMover_Perturb_QuickCCDCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "LoopMover_Perturb_QuickCCD" ); }

	void test_protocols_loops_loop_mover_perturb_LoopMover_Perturb_QuickCCD_MovesCreator_name()
	{ protocols::loops::loop_mover::perturb::LoopMover_Perturb_QuickCCD_MovesCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "LoopMover_Perturb_QuickCCD_Moves" ); }

	void test_protocols_loops_loop_mover_refine_LoopMover_Refine_BackrubCreator_name()
	{ protocols::loops::loop_mover::refine::LoopMover_Refine_BackrubCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "LoopMover_Refine_Backrub" ); }

	void test_protocols_loops_loop_mover_refine_LoopMover_Refine_CCDCreator_name()
	{ protocols::loops::loop_mover::refine::LoopMover_Refine_CCDCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "LoopMover_Refine_CCD" ); }

	void test_protocols_loops_loop_mover_refine_LoopMover_Refine_KICCreator_name()
	{ protocols::loops::loop_mover::refine::LoopMover_Refine_KICCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "LoopMover_Refine_KIC" ); }

	void test_protocols_loops_loop_mover_refine_LoopRefineInnerCycleContainerCreator_name()
	{ protocols::loops::loop_mover::refine::LoopRefineInnerCycleContainerCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "LoopRefineInnerCycleContainer" ); }

	void test_protocols_loops_loop_mover_refine_RepackTrialCreator_name()
	{ protocols::loops::loop_mover::refine::RepackTrialCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "RepackTrial" ); }

	void test_protocols_loops_loop_mover_refine_ShearMinCCDTrialCreator_name()
	{ protocols::loops::loop_mover::refine::ShearMinCCDTrialCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ShearMinCCDTrial" ); }

	void test_protocols_loops_loop_mover_refine_SmallMinCCDTrialCreator_name()
	{ protocols::loops::loop_mover::refine::SmallMinCCDTrialCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SmallMinCCDTrial" ); }

	void test_protocols_loops_loop_mover_LoopCMCreator_name()
	{ protocols::loops::loop_mover::LoopCMCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "LoopCM" ); }

	void test_protocols_match_MatcherMoverCreator_name()
	{ protocols::match::MatcherMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "MatcherMover" ); }

	void test_protocols_matdes_ExtractSubposeMoverCreator_name()
	{ protocols::matdes::ExtractSubposeMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ExtractSubposeMover" ); }

	void test_protocols_matdes_MatDesGreedyOptMutationMoverCreator_name()
	{ protocols::matdes::MatDesGreedyOptMutationMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "MatDesGreedyOptMutationMover" ); }

	void test_protocols_matdes_SymDofMoverCreator_name()
	{ protocols::matdes::SymDofMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SymDofMover" ); }

	void test_protocols_matdes_SchemePlaceMotifsMoverCreator_name()
	{ protocols::matdes::SchemePlaceMotifsMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SchemePlaceMotifs" ); }

	void test_protocols_md_CartesianMDCreator_name()
	{ protocols::md::CartesianMDCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "CartesianMD" ); }

	void test_protocols_simple_moves_BBGaussianMoverCreator_name()
	{ protocols::simple_moves::BBGaussianMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "BBGaussian" ); }

	void test_protocols_simple_moves_PeriodicBoxMoverCreator_name()
	{ protocols::simple_moves::PeriodicBoxMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "PeriodicBoxMover" ); }

	void test_protocols_motifs_MotifDnaPackerCreator_name()
	{ protocols::motifs::MotifDnaPackerCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "MotifDnaPacker" ); }

	void test_protocols_motif_grafting_movers_MotifGraftCreator_name()
	{ protocols::motif_grafting::movers::MotifGraftCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "MotifGraft" ); }

	void test_protocols_moves_DsspMoverCreator_name()
	{ protocols::moves::DsspMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "Dssp" ); }

	void test_protocols_moves_IfMoverCreator_name()
	{ protocols::moves::IfMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "If" ); }

	void test_protocols_moves_RandomMoverCreator_name()
	{ protocols::moves::RandomMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "RandomMover" ); }

	void test_protocols_moves_IteratedConvergenceMoverCreator_name()
	{ protocols::moves::IteratedConvergenceMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "IteratedConvergence" ); }

	void test_protocols_moves_PyMOLMoverCreator_name()
	{ protocols::moves::PyMOLMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "PyMOLMover" ); }

	void test_protocols_moves_RampingMoverCreator_name()
	{ protocols::moves::RampingMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "RampingMover" ); }

	void test_protocols_moves_FilterReportAsPoseExtraScoresMoverCreator_name()
	{ protocols::moves::FilterReportAsPoseExtraScoresMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "FilterReportAsPoseExtraScoresMover" ); }

	void test_protocols_simple_moves_MonteCarloResetCreator_name()
	{ protocols::monte_carlo::MonteCarloResetCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "MonteCarloReset" ); }

	void test_protocols_simple_moves_AlignChainMoverCreator_name()
	{ protocols::simple_moves::AlignChainMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "AlignChain" ); }

	void test_protocols_simple_moves_ResetBaselineMoverCreator_name()
	{ protocols::monte_carlo::ResetBaselineMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ResetBaseline" ); }

	void test_protocols_simple_moves_RingConformationMoverCreator_name()
	{ protocols::simple_moves::RingConformationMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "RingConformationMover" ); }

	void test_protocols_simple_moves_SimpleThreadingMoverCreator_name()
	{ protocols::simple_moves::SimpleThreadingMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SimpleThreadingMover" ); }

	void test_protocols_nonlocal_SingleFragmentMoverCreator_name()
	{ protocols::nonlocal::SingleFragmentMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SingleFragmentMover" ); }

	void test_protocols_normalmode_NormalModeRelaxMoverCreator_name()
	{ protocols::normalmode::NormalModeRelaxMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "NormalModeRelax" ); }

	void test_protocols_normalmode_NormalModeMinimizerCreator_name()
	{ protocols::normalmode::NormalModeMinimizerCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "NormalModeMinimizer" ); }

	void test_protocols_pb_potential_SetupPoissonBoltzmannPotentialCreator_name()
	{ protocols::pb_potential::SetupPoissonBoltzmannPotentialCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SetupPoissonBoltzmannPotential" ); }

	void test_protocols_pose_length_moves_FixAllLoopsMoverCreator_name()
	{ protocols::pose_length_moves::FixAllLoopsMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "FixAllLoopsMover" ); }

	void test_protocols_pose_length_moves_InsertResMoverCreator_name()
	{ protocols::pose_length_moves::InsertResMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "InsertResMover" ); }

	void test_protocols_pose_length_moves_NearNativeLoopCloserCreator_name()
	{ protocols::pose_length_moves::NearNativeLoopCloserCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "NearNativeLoopCloser" ); }

	void test_protocols_protein_interface_design_movers_AddChainBreakCreator_name()
	{ protocols::protein_interface_design::movers::AddChainBreakCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "AddChainBreak" ); }

	void test_protocols_protein_interface_design_movers_AddSidechainConstraintsToHotspotsCreator_name()
	{ protocols::protein_interface_design::movers::AddSidechainConstraintsToHotspotsCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "AddSidechainConstraintsToHotspots" ); }

	void test_protocols_protein_interface_design_movers_BackrubDDMoverCreator_name()
	{ protocols::protein_interface_design::movers::BackrubDDMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "BackrubDD" ); }

	void test_protocols_protein_interface_design_movers_BestHotspotCstMoverCreator_name()
	{ protocols::protein_interface_design::movers::BestHotspotCstMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "BestHotspotCst" ); }

	void test_protocols_protein_interface_design_movers_BuildAlaPoseCreator_name()
	{ protocols::protein_interface_design::movers::BuildAlaPoseCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "build_Ala_pose" ); }

	void test_protocols_protein_interface_design_movers_DesignMinimizeHbondsCreator_name()
	{ protocols::protein_interface_design::movers::DesignMinimizeHbondsCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "DesignMinimizeHbonds" ); }

	void test_protocols_protein_interface_design_movers_DisulfideMoverCreator_name()
	{ protocols::protein_interface_design::movers::DisulfideMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "DisulfideMover" ); }

	void test_protocols_protein_interface_design_movers_DockAndRetrieveSidechainsCreator_name()
	{ protocols::protein_interface_design::movers::DockAndRetrieveSidechainsCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "Docking" ); }

	void test_protocols_simple_moves_DumpPdbCreator_name()
	{ protocols::simple_moves::DumpPdbCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "DumpPdb" ); }

	void test_protocols_protein_interface_design_movers_FavorNativeResiduePreCycleCreator_name()
	{ protocols::protein_interface_design::movers::FavorNativeResiduePreCycleCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "FavorNativeResidue" ); }

	void test_protocols_protein_interface_design_movers_FavorNonNativeResiduePreCycleCreator_name()
	{ protocols::protein_interface_design::movers::FavorNonNativeResiduePreCycleCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "FavorNonNativeResidue" ); }

	void test_protocols_protein_interface_design_movers_HotspotDisjointedFoldTreeMoverCreator_name()
	{ protocols::protein_interface_design::movers::HotspotDisjointedFoldTreeMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "HotspotDisjointedFoldTree" ); }

	void test_protocols_protein_interface_design_movers_HotspotHasherMoverCreator_name()
	{ protocols::protein_interface_design::movers::HotspotHasherMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "HotspotHasher" ); }

	void test_protocols_protein_interface_design_movers_InterfaceRecapitulationMoverCreator_name()
	{ protocols::protein_interface_design::movers::InterfaceRecapitulationMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "InterfaceRecapitulation" ); }

	void test_protocols_protein_interface_design_movers_LoopFinderCreator_name()
	{ protocols::protein_interface_design::movers::LoopFinderCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "LoopFinder" ); }

	void test_protocols_protein_interface_design_movers_LoopLengthChangeCreator_name()
	{ protocols::protein_interface_design::movers::LoopLengthChangeCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "LoopLengthChange" ); }

	void test_protocols_protein_interface_design_movers_LoopMoverFromCommandLineCreator_name()
	{ protocols::protein_interface_design::movers::LoopMoverFromCommandLineCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "LoopMoverFromCommandLine" ); }

	void test_protocols_protein_interface_design_movers_LoopOverCreator_name()
	{ protocols::protein_interface_design::movers::LoopOverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "LoopOver" ); }

	void test_protocols_protein_interface_design_movers_LoopRemodelCreator_name()
	{ protocols::protein_interface_design::movers::LoopRemodelCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "LoopRemodel" ); }

	void test_protocols_protein_interface_design_movers_MapHotspotCreator_name()
	{ protocols::protein_interface_design::movers::MapHotspotCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "MapHotspot" ); }

	void test_protocols_protein_interface_design_movers_PatchdockTransformCreator_name()
	{ protocols::protein_interface_design::movers::PatchdockTransformCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "PatchdockTransform" ); }

	void test_protocols_protein_interface_design_movers_PeptideStapleDesignMoverCreator_name()
	{ protocols::protein_interface_design::movers::PeptideStapleDesignMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "StapleMover" ); }

	void test_protocols_protein_interface_design_movers_PlacementAuctionMoverCreator_name()
	{ protocols::protein_interface_design::movers::PlacementAuctionMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "Auction" ); }

	void test_protocols_protein_interface_design_movers_PlacementMinimizationMoverCreator_name()
	{ protocols::protein_interface_design::movers::PlacementMinimizationMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "PlacementMinimization" ); }

	void test_protocols_protein_interface_design_movers_PlaceOnLoopCreator_name()
	{ protocols::protein_interface_design::movers::PlaceOnLoopCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "PlaceOnLoop" ); }

	void test_protocols_protein_interface_design_movers_PlaceSimultaneouslyMoverCreator_name()
	{ protocols::protein_interface_design::movers::PlaceSimultaneouslyMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "PlaceSimultaneously" ); }

	void test_protocols_protein_interface_design_movers_PlaceStubMoverCreator_name()
	{ protocols::protein_interface_design::movers::PlaceStubMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "PlaceStub" ); }

	void test_protocols_protein_interface_design_movers_PrepackMoverCreator_name()
	{ protocols::protein_interface_design::movers::PrepackMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "Prepack" ); }

	void test_protocols_protein_interface_design_movers_ProteinInterfaceMultiStateDesignMoverCreator_name()
	{ protocols::protein_interface_design::movers::ProteinInterfaceMultiStateDesignMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ProteinInterfaceMS" ); }

	void test_protocols_protein_interface_design_movers_RandomMutationCreator_name()
	{ protocols::protein_interface_design::movers::RandomMutationCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "RandomMutation" ); }

	void test_protocols_protein_interface_design_movers_RepackMinimizeCreator_name()
	{ protocols::protein_interface_design::movers::RepackMinimizeCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "RepackMinimize" ); }

	void test_protocols_protein_interface_design_movers_SaveAndRetrieveSidechainsCreator_name()
	{ protocols::protein_interface_design::movers::SaveAndRetrieveSidechainsCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SaveAndRetrieveSidechains" ); }

	void test_protocols_protein_interface_design_movers_SetAtomTreeCreator_name()
	{ protocols::protein_interface_design::movers::SetAtomTreeCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "AtomTree" ); }

	void test_protocols_protein_interface_design_movers_SetTemperatureFactorCreator_name()
	{ protocols::protein_interface_design::movers::SetTemperatureFactorCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SetTemperatureFactor" ); }

	void test_protocols_protein_interface_design_movers_SetupHotspotConstraintsMoverCreator_name()
	{ protocols::protein_interface_design::movers::SetupHotspotConstraintsMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SetupHotspotConstraints" ); }

	void test_protocols_protein_interface_design_movers_SetupHotspotConstraintsLoopsMoverCreator_name()
	{ protocols::protein_interface_design::movers::SetupHotspotConstraintsLoopsMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SetupHotspotConstraintsLoops" ); }

	void test_protocols_protein_interface_design_movers_SpinMoverCreator_name()
	{ protocols::protein_interface_design::movers::SpinMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SpinMover" ); }

	void test_protocols_protein_interface_design_movers_TaskAwareCstsCreator_name()
	{ protocols::protein_interface_design::movers::TaskAwareCstsCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "TaskAwareCsts" ); }

	void test_protocols_protein_interface_design_movers_SubroutineMoverCreator_name()
	{ protocols::protein_interface_design::movers::SubroutineMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "Subroutine" ); }

	void test_protocols_protein_interface_design_movers_TryRotamersCreator_name()
	{ protocols::protein_interface_design::movers::TryRotamersCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "TryRotamers" ); }

	void test_protocols_protein_interface_design_movers_ShoveResidueMoverCreator_name()
	{ protocols::protein_interface_design::movers::ShoveResidueMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ShoveResidueMover" ); }

	void test_protocols_protein_interface_design_movers_VLBCreator_name()
	{ protocols::protein_interface_design::movers::VLBCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "VLB" ); }

	void test_protocols_protein_interface_design_movers_DockWithHotspotMoverCreator_name()
	{ protocols::protein_interface_design::movers::DockWithHotspotMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "DockWithHotspotMover" ); }

	void test_protocols_protein_interface_design_movers_TopologyBrokerMoverCreator_name()
	{ protocols::protein_interface_design::movers::TopologyBrokerMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "TopologyBrokerMover" ); }

	void test_protocols_qsar_RenderGridsToKinemageCreator_name()
	{ protocols::qsar::RenderGridsToKinemageCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "RenderGridsToKinemage" ); }

	void test_protocols_rbsegment_relax_MakeStarTopologyMoverCreator_name()
	{ protocols::rbsegment_relax::MakeStarTopologyMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "MakeStarTopology" ); }

	void test_protocols_rbsegment_relax_OptimizeThreadingMoverCreator_name()
	{ protocols::rbsegment_relax::OptimizeThreadingMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "OptimizeThreading" ); }

	void test_protocols_rbsegment_relax_IdealizeHelicesMoverCreator_name()
	{ protocols::rbsegment_relax::IdealizeHelicesMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "IdealizeHelices" ); }

	void test_protocols_relax_AtomCoordinateCstMoverCreator_name()
	{ protocols::relax::AtomCoordinateCstMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "AtomCoordinateCstMover" ); }

	void test_protocols_relax_FastRelaxCreator_name()
	{ protocols::relax::FastRelaxCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "FastRelax" ); }

	void test_protocols_relax_LocalRelaxCreator_name()
	{ protocols::relax::LocalRelaxCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "LocalRelax" ); }

	void test_protocols_residue_selectors_StoreResidueSubsetMoverCreator_name()
	{ protocols::residue_selectors::StoreResidueSubsetMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "StoreResidueSubset" ); }

	void test_protocols_rigid_RigidBodyPerturbNoCenterMoverCreator_name()
	{ protocols::rigid::RigidBodyPerturbNoCenterMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "RigidBodyPerturbNoCenter" ); }

	void test_protocols_rigid_RigidBodyTiltMoverCreator_name()
	{ protocols::rigid::RigidBodyTiltMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "RigidBodyTiltMover" ); }

	void test_protocols_rigid_RigidBodyTransMoverCreator_name()
	{ protocols::rigid::RigidBodyTransMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "RigidBodyTransMover" ); }

	void test_protocols_rigid_RollMoverCreator_name()
	{ protocols::rigid::RollMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "RollMover" ); }

	void test_protocols_rigid_UniformRigidBodyCMCreator_name()
	{ protocols::rigid::UniformRigidBodyCMCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "UniformRigidBodyCM" ); }

	void test_protocols_rigid_UniformRigidBodyMoverCreator_name()
	{ protocols::rigid::UniformRigidBodyMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "UniformRigidBodyMover" ); }

	void test_protocols_rosetta_scripts_ParsedProtocolCreator_name()
	{ protocols::rosetta_scripts::ParsedProtocolCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ParsedProtocol" ); }

	void test_protocols_rosetta_scripts_SavePoseMoverCreator_name()
	{ protocols::rosetta_scripts::SavePoseMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SavePoseMover" ); }

	void test_protocols_rosetta_scripts_MultiplePoseMoverCreator_name()
	{ protocols::rosetta_scripts::MultiplePoseMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "MultiplePoseMover" ); }

	void test_protocols_rosetta_scripts_MultipleOutputWrapperCreator_name()
	{ protocols::rosetta_scripts::MultipleOutputWrapperCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "MultipleOutputWrapper" ); }

	void test_protocols_rotamer_recovery_RotamerRecoveryMoverCreator_name()
	{ protocols::rotamer_recovery::RotamerRecoveryMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "RotamerRecoveryMover" ); }

	void test_protocols_seeded_abinitio_CAcstGeneratorCreator_name()
	{ protocols::seeded_abinitio::CAcstGeneratorCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "CAcstGenerator" ); }

	void test_protocols_seeded_abinitio_CloseFoldCreator_name()
	{ protocols::seeded_abinitio::CloseFoldCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "CloseFold" ); }

	void test_protocols_seeded_abinitio_CoordinateCstCreator_name()
	{ protocols::seeded_abinitio::CoordinateCstCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "CoordinateCst" ); }

	void test_protocols_seeded_abinitio_DefineMovableLoopsCreator_name()
	{ protocols::seeded_abinitio::DefineMovableLoopsCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "DefineMovableLoops" ); }

	void test_protocols_seeded_abinitio_GrowPeptidesCreator_name()
	{ protocols::seeded_abinitio::GrowPeptidesCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "GrowPeptides" ); }

	void test_protocols_seeded_abinitio_SeedFoldTreeCreator_name()
	{ protocols::seeded_abinitio::SeedFoldTreeCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SeedFoldTree" ); }

	void test_protocols_seeded_abinitio_SeedSetupMoverCreator_name()
	{ protocols::seeded_abinitio::SeedSetupMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SeedSetupMover" ); }

	void test_protocols_seeded_abinitio_SwapSegmentCreator_name()
	{ protocols::seeded_abinitio::SwapSegmentCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SwapSegment" ); }

	void test_protocols_seeded_abinitio_SegmentHybridizerCreator_name()
	{ protocols::seeded_abinitio::SegmentHybridizerCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SegmentHybridizer" ); }

	void test_protocols_simple_moves_AddConstraintsToCurrentConformationMoverCreator_name()
	{ protocols::constraint_movers::AddConstraintsToCurrentConformationMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "AddConstraintsToCurrentConformationMover" ); }

	void test_protocols_simple_moves_AddChainMoverCreator_name()
	{ protocols::simple_moves::AddChainMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "AddChain" ); }

	void test_protocols_simple_moves_AddJobPairDataCreator_name()
	{ protocols::simple_moves::AddJobPairDataCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "AddJobPairData" ); }

	void test_protocols_simple_moves_ChangeAndResetFoldTreeMoverCreator_name()
	{ protocols::simple_moves::ChangeAndResetFoldTreeMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ChangeAndResetFoldTreeMover" ); }

	void test_protocols_simple_moves_BoltzmannRotamerMoverCreator_name()
	{ protocols::minimization_packing::BoltzmannRotamerMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "BoltzmannRotamerMover" ); }

	void test_protocols_simple_moves_ClearConstraintsMoverCreator_name()
	{ protocols::constraint_movers::ClearConstraintsMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ClearConstraintsMover" ); }

	void test_protocols_simple_moves_ConsensusDesignMoverCreator_name()
	{ protocols::calc_taskop_movers::ConsensusDesignMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ConsensusDesignMover" ); }

	void test_protocols_simple_moves_ConstraintSetMoverCreator_name()
	{ protocols::constraint_movers::ConstraintSetMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ConstraintSetMover" ); }

	void test_protocols_simple_moves_ContingentAcceptMoverCreator_name()
	{ protocols::simple_moves::ContingentAcceptMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ContingentAccept" ); }

	void test_protocols_simple_moves_CoupledMoverCreator_name()
	{ protocols::simple_moves::CoupledMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "CoupledMover" ); }

	void test_protocols_simple_moves_DeleteChainMoverCreator_name()
	{ protocols::simple_moves::DeleteChainMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "DeleteChain" ); }

	void test_protocols_simple_moves_DeleteChainsMoverCreator_name()
	{ protocols::simple_moves::DeleteChainsMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "DeleteChainsMover" ); }

	void test_protocols_simple_moves_ddGCreator_name()
	{ protocols::simple_ddg::ddGCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ddG" ); }

	void test_protocols_simple_moves_DisulfideInsertionMoverCreator_name()
	{ protocols::simple_moves::DisulfideInsertionMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "DisulfideInsertion" ); }

	void test_protocols_simple_moves_ExtendedPoseMoverCreator_name()
	{ protocols::pose_creation::ExtendedPoseMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ExtendedPoseMover" ); }

	void test_protocols_simple_moves_FavorSequenceProfileCreator_name()
	{ protocols::simple_moves::FavorSequenceProfileCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "FavorSequenceProfile" ); }

	void test_protocols_simple_moves_FavorSymmetricSequenceCreator_name()
	{ protocols::simple_moves::FavorSymmetricSequenceCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "FavorSymmetricSequence" ); }

	void test_protocols_simple_moves_FindConsensusSequenceCreator_name()
	{ protocols::simple_moves::FindConsensusSequenceCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "FindConsensusSequence" ); }

	void test_protocols_simple_moves_ForceDisulfidesMoverCreator_name()
	{ protocols::calc_taskop_movers::ForceDisulfidesMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ForceDisulfides" ); }

	void test_protocols_simple_moves_GenericMonteCarloMoverCreator_name()
	{ protocols::monte_carlo::GenericMonteCarloMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "GenericMonteCarlo" ); }

	void test_protocols_simple_moves_LoadPDBMoverCreator_name()
	{ protocols::simple_moves::LoadPDBMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "LoadPDB" ); }

	void test_protocols_simple_moves_LoadUnboundRotMoverCreator_name()
	{ protocols::simple_moves::LoadUnboundRotMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "load_unbound_rot" ); }

	void test_protocols_simple_moves_MakePolyXMoverCreator_name()
	{ protocols::pose_creation::MakePolyXMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "MakePolyX" ); }

	void test_protocols_simple_moves_MembraneTopologyCreator_name()
	{ protocols::simple_moves::MembraneTopologyCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "MembraneTopology" ); }

	void test_protocols_pose_creation_MergePDBMoverCreator_name()
	{ protocols::pose_creation::MergePDBMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "MergePDB" ); }

	void test_protocols_simple_moves_MinMoverCreator_name()
	{ protocols::minimization_packing::MinMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "MinMover" ); }

	void test_protocols_simple_moves_MinPackMoverCreator_name()
	{ protocols::minimization_packing::MinPackMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "MinPackMover" ); }

	void test_protocols_simple_moves_ModifyVariantTypeMoverCreator_name()
	{ protocols::simple_moves::ModifyVariantTypeMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ModifyVariantType" ); }

	void test_protocols_simple_moves_MonteCarloRecoverCreator_name()
	{ protocols::monte_carlo::MonteCarloRecoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "MonteCarloRecover" ); }

	void test_protocols_simple_moves_MonteCarloTestCreator_name()
	{ protocols::monte_carlo::MonteCarloTestCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "MonteCarloTest" ); }

	void test_protocols_simple_moves_MSDMoverCreator_name()
	{ protocols::simple_moves::MSDMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "MSDMover" ); }

	void test_protocols_simple_moves_MutateResidueCreator_name()
	{ protocols::simple_moves::MutateResidueCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "MutateResidue" ); }

	void test_protocols_simple_moves_PackRotamersMoverCreator_name()
	{ protocols::minimization_packing::PackRotamersMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "PackRotamersMover" ); }

	void test_protocols_simple_moves_PSSM2BfactorMoverCreator_name()
	{ protocols::simple_moves::PSSM2BfactorMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "PSSM2Bfactor" ); }

	void test_protocols_simple_moves_RepeatPropagationMoverCreator_name()
	{ protocols::simple_moves::RepeatPropagationMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "RepeatPropagationMover" ); }

	void test_protocols_simple_moves_ResidueTypeConstraintMoverCreator_name()
	{ protocols::constraint_movers::ResidueTypeConstraintMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ResidueTypeConstraintMover" ); }

	void test_protocols_simple_moves_ReportEffectivePKACreator_name()
	{ protocols::simple_moves::ReportEffectivePKACreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ReportEffectivePKA" ); }

	void test_protocols_simple_moves_RotamerTrialsMinMoverCreator_name()
	{ protocols::minimization_packing::RotamerTrialsMinMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "RotamerTrialsMinMover" ); }

	void test_protocols_simple_moves_RotamerTrialsMoverCreator_name()
	{ protocols::minimization_packing::RotamerTrialsMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "RotamerTrialsMover" ); }

	void test_protocols_simple_moves_RandomTorsionMoverCreator_name()
	{ protocols::simple_moves::RandomTorsionMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "RandomTorsionMover" ); }

	void test_protocols_simple_moves_RandomOmegaFlipMoverCreator_name()
	{ protocols::simple_moves::RandomOmegaFlipMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "RandomOmegaFlipMover" ); }

	void test_protocols_simple_moves_SaneMinMoverCreator_name()
	{ protocols::minimization_packing::SaneMinMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SaneMinMover" ); }

	void test_protocols_simple_moves_ScoreMoverCreator_name()
	{ protocols::simple_moves::ScoreMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ScoreMover" ); }

	void test_protocols_simple_moves_SequenceProfileMoverCreator_name()
	{ protocols::simple_moves::SequenceProfileMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "profile" ); }

	void test_protocols_simple_moves_SetTorsionCreator_name()
	{ protocols::simple_moves::SetTorsionCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SetTorsion" ); }

	void test_protocols_simple_moves_ShearMoverCreator_name()
	{ protocols::simple_moves::ShearMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "Shear" ); }

	void test_protocols_simple_moves_ShortBackrubMoverCreator_name()
	{ protocols::simple_moves::ShortBackrubMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ShortBackrubMover" ); }

	void test_protocols_simple_moves_StorePoseSnapshotCreator_name()
	{ protocols::simple_moves::StorePoseSnapshotCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "StorePoseSnapshot" ); }

	void test_protocols_simple_moves_SmallMoverCreator_name()
	{ protocols::simple_moves::SmallMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "Small" ); }

	void test_protocols_simple_moves_StructProfileMoverCreator_name()
	{ protocols::simple_moves::StructProfileMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "StructProfileMover" ); }

	void test_protocols_simple_moves_SuperimposeMoverCreator_name()
	{ protocols::simple_moves::SuperimposeMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "Superimpose" ); }

	void test_protocols_simple_moves_PDBReloadMoverCreator_name()
	{ protocols::simple_moves::PDBReloadMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "PDBReload" ); }

	void test_protocols_simple_moves_SwitchResidueTypeSetMoverCreator_name()
	{ protocols::simple_moves::SwitchResidueTypeSetMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SwitchResidueTypeSetMover" ); }

	void test_protocols_simple_moves_SwitchChainOrderMoverCreator_name()
	{ protocols::simple_moves::SwitchChainOrderMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SwitchChainOrder" ); }

	void test_protocols_simple_moves_TaskAwareMinMoverCreator_name()
	{ protocols::minimization_packing::TaskAwareMinMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "TaskAwareMinMover" ); }

	void test_protocols_simple_moves_VirtualRootMoverCreator_name()
	{ protocols::simple_moves::VirtualRootMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "VirtualRoot" ); }

	void test_protocols_simple_moves_TumbleCreator_name()
	{ protocols::simple_moves::TumbleCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "Tumble" ); }

	void test_protocols_simple_moves_bin_transitions_InitializeByBinsCreator_name()
	{ protocols::simple_moves::bin_transitions::InitializeByBinsCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "InitializeByBins" ); }

	void test_protocols_simple_moves_bin_transitions_PerturbByBinsCreator_name()
	{ protocols::simple_moves::bin_transitions::PerturbByBinsCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "PerturbByBins" ); }

	void test_protocols_simple_moves_sidechain_moves_SidechainMCMoverCreator_name()
	{ protocols::simple_moves::sidechain_moves::SidechainMCMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SidechainMC" ); }

	void test_protocols_simple_moves_sidechain_moves_SidechainMoverCreator_name()
	{ protocols::simple_moves::sidechain_moves::SidechainMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "Sidechain" ); }

	void test_protocols_simple_moves_sidechain_moves_SetChiMoverCreator_name()
	{ protocols::simple_moves::sidechain_moves::SetChiMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SetChiMover" ); }

	void test_protocols_simple_moves_sidechain_moves_JumpRotamerSidechainMoverCreator_name()
	{ protocols::simple_moves::sidechain_moves::JumpRotamerSidechainMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "JumpRotamerSidechain" ); }

	void test_protocols_simple_moves_sidechain_moves_PerturbRotamerSidechainMoverCreator_name()
	{ protocols::simple_moves::sidechain_moves::PerturbRotamerSidechainMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "PerturbRotamerSidechain" ); }

	void test_protocols_simple_moves_sidechain_moves_PerturbChiSidechainMoverCreator_name()
	{ protocols::simple_moves::sidechain_moves::PerturbChiSidechainMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "PerturbChiSidechain" ); }

	void test_protocols_simple_moves_symmetry_ExtractAsymmetricUnitMoverCreator_name()
	{ protocols::symmetry::ExtractAsymmetricUnitMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ExtractAsymmetricUnit" ); }

	void test_protocols_simple_moves_symmetry_ExtractAsymmetricPoseMoverCreator_name()
	{ protocols::symmetry::ExtractAsymmetricPoseMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ExtractAsymmetricPose" ); }

	void test_protocols_simple_moves_symmetry_SetupForSymmetryMoverCreator_name()
	{ protocols::symmetry::SetupForSymmetryMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SetupForSymmetry" ); }

	void test_protocols_simple_moves_symmetry_DetectSymmetryMoverCreator_name()
	{ protocols::symmetry::DetectSymmetryMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "DetectSymmetry" ); }

	void test_protocols_simple_moves_symmetry_SetupNCSMoverCreator_name()
	{ protocols::symmetry::SetupNCSMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SetupNCS" ); }

	void test_protocols_simple_moves_symmetry_SymMinMoverCreator_name()
	{ protocols::minimization_packing::symmetry::SymMinMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SymMinMover" ); }

	void test_protocols_simple_moves_symmetry_SymPackRotamersMoverCreator_name()
	{ protocols::minimization_packing::symmetry::SymPackRotamersMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SymPackRotamersMover" ); }

	void test_protocols_simple_moves_symmetry_SymRotamerTrialsMoverCreator_name()
	{ protocols::minimization_packing::symmetry::SymRotamerTrialsMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SymRotamerTrialsMover" ); }

	void test_protocols_simple_moves_symmetry_TaskAwareSymMinMoverCreator_name()
	{ protocols::minimization_packing::symmetry::TaskAwareSymMinMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "TaskAwareSymMinMover" ); }

	void test_protocols_toolbox_task_operations_StoreCombinedStoredTasksMoverCreator_name()
	{ protocols::task_operations::StoreCombinedStoredTasksMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "StoreCombinedStoredTasksMover" ); }

	void test_protocols_toolbox_task_operations_StoreCompoundTaskMoverCreator_name()
	{ protocols::task_operations::StoreCompoundTaskMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "StoreCompoundTaskMover" ); }

	void test_protocols_toolbox_task_operations_StoreTaskMoverCreator_name()
	{ protocols::task_operations::StoreTaskMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "StoreTaskMover" ); }

	void test_protocols_environment_EnvMoverCreator_name()
	{ protocols::environment::EnvMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "Environment" ); }

	void test_protocols_environment_CoMTrackerCMCreator_name()
	{ protocols::environment::CoMTrackerCMCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "CoMTrackerCM" ); }

	void test_protocols_environment_ScriptCMCreator_name()
	{ protocols::environment::ScriptCMCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ScriptCM" ); }

	void test_protocols_abinitio_abscript_AbscriptMoverCreator_name()
	{ protocols::abinitio::abscript::AbscriptMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "AbscriptMover" ); }

	void test_protocols_abinitio_abscript_ConstraintPreparerCreator_name()
	{ protocols::abinitio::abscript::ConstraintPreparerCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ConstraintPreparer" ); }

	void test_protocols_abinitio_abscript_FragmentCMCreator_name()
	{ protocols::abinitio::abscript::FragmentCMCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "FragmentCM" ); }

	void test_protocols_abinitio_abscript_FragmentJumpCMCreator_name()
	{ protocols::abinitio::abscript::FragmentJumpCMCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "FragmentJumpCM" ); }

	void test_protocols_abinitio_abscript_AbscriptLoopCloserCMCreator_name()
	{ protocols::abinitio::abscript::AbscriptLoopCloserCMCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "AbscriptLoopCloserCM" ); }

	void test_protocols_abinitio_abscript_StructPerturberCMCreator_name()
	{ protocols::abinitio::abscript::StructPerturberCMCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "StructPerturberCM" ); }

	void test_protocols_abinitio_abscript_RigidChunkCMCreator_name()
	{ protocols::abinitio::abscript::RigidChunkCMCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "RigidChunkCM" ); }

	void test_protocols_membrane_AddMembraneMoverCreator_name()
	{ protocols::membrane::AddMembraneMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "AddMembraneMover" ); }

	void test_protocols_membrane_FlipMoverCreator_name()
	{ protocols::membrane::FlipMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "FlipMover" ); }

	void test_protocols_membrane_MembranePositionFromTopologyMoverCreator_name()
	{ protocols::membrane::MembranePositionFromTopologyMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "MembranePositionFromTopologyMover" ); }

	void test_protocols_membrane_SetMembranePositionMoverCreator_name()
	{ protocols::membrane::SetMembranePositionMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SetMembranePositionMover" ); }

	void test_protocols_membrane_TransformIntoMembraneMoverCreator_name()
	{ protocols::membrane::TransformIntoMembraneMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "TransformIntoMembraneMover" ); }

	void test_protocols_membrane_symmetry_SymmetricAddMembraneMoverCreator_name()
	{ protocols::membrane::symmetry::SymmetricAddMembraneMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SymmetricAddMembraneMover" ); }

	void test_protocols_membrane_benchmark_SampleTiltAnglesCreator_name()
	{ protocols::membrane::benchmark::SampleTiltAnglesCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SampleTiltAngles" ); }

	void test_protocols_membrane_benchmark_MakeCanonicalHelixCreator_name()
	{ protocols::membrane::benchmark::MakeCanonicalHelixCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "MakeCanonicalHelix" ); }

	void test_protocols_membrane_visualize_VisualizeMembraneMoverCreator_name()
	{ protocols::membrane::visualize::VisualizeMembraneMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "VisualizeMembraneMover" ); }

	void test_protocols_membrane_visualize_VisualizeEmbeddingMoverCreator_name()
	{ protocols::membrane::visualize::VisualizeEmbeddingMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "VisualizeEmbeddingMover" ); }

	void test_protocols_docking_membrane_MPDockingMoverCreator_name()
	{ protocols::docking::membrane::MPDockingMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "MPDockingMover" ); }

	void test_protocols_docking_membrane_MPDockingSetupMoverCreator_name()
	{ protocols::docking::membrane::MPDockingSetupMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "MPDockingSetupMover" ); }

	void test_protocols_legacy_sewing_LegacyAppendAssemblyMoverCreator_name()
	{ protocols::legacy_sewing::LegacyAppendAssemblyMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "LegacyAppendAssemblyMover" ); }

	void test_protocols_legacy_sewing_LegacyGivenPathAssemblyMoverCreator_name()
	{ protocols::legacy_sewing::LegacyGivenPathAssemblyMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "LegacyGivenPathAssemblyMover" ); }

	void test_protocols_legacy_sewing_LegacyGreedyAssemblyMoverCreator_name()
	{ protocols::legacy_sewing::LegacyGreedyAssemblyMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "LegacyGreedyAssemblyMover" ); }

	void test_protocols_legacy_sewing_LegacyMonteCarloAssemblyMoverCreator_name()
	{ protocols::legacy_sewing::LegacyMonteCarloAssemblyMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "LegacyMonteCarloAssemblyMover" ); }

	void test_protocols_legacy_sewing_LegacyAssemblyConstraintsMoverCreator_name()
	{ protocols::legacy_sewing::LegacyAssemblyConstraintsMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "LegacyAssemblyConstraintsMover" ); }

	void test_protocols_legacy_sewing_LegacyAddStartnodeFragmentsCreator_name()
	{ protocols::legacy_sewing::LegacyAddStartnodeFragmentsCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "LegacyAddStartnodeFragments" ); }

	void test_protocols_legacy_sewing_LegacyRepeatAssemblyMoverCreator_name()
	{ protocols::legacy_sewing::LegacyRepeatAssemblyMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "LegacyRepeatAssemblyMover" ); }

	void test_protocols_legacy_sewing_LegacyEnumerateAssemblyMoverCreator_name()
	{ protocols::legacy_sewing::LegacyEnumerateAssemblyMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "LegacyEnumerateAssemblyMover" ); }

	void test_protocols_symmetric_docking_membrane_MPSymDockMoverCreator_name()
	{ protocols::symmetric_docking::membrane::MPSymDockMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "MPSymDockMover" ); }

	void test_protocols_membrane_AddMPLigandMoverCreator_name()
	{ protocols::membrane::AddMPLigandMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "AddMPLigandMover" ); }

	void test_protocols_relax_membrane_MPFastRelaxMoverCreator_name()
	{ protocols::relax::membrane::MPFastRelaxMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "MPFastRelaxMover" ); }

	void test_protocols_vardist_solaccess_LoadVarSolDistSasaCalculatorMoverCreator_name()
	{ protocols::vardist_solaccess::LoadVarSolDistSasaCalculatorMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "LoadVarSolDistSasaCalculatorMover" ); }

	// void test_protocols_aa_composition_AddCompositionConstraintMoverCreator()
	// { protocols::aa_composition::AddCompositionConstraintMoverCreator cr; std::cout << "protocols::aa_composition::AddCompositionConstraintMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_aa_composition_ClearCompositionConstraintsMoverCreator()
	// { protocols::aa_composition::ClearCompositionConstraintsMoverCreator cr; std::cout << "protocols::aa_composition::ClearCompositionConstraintsMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_abinitio_DomainAssemblyCreator()
	// { protocols::abinitio::DomainAssemblyCreator cr; std::cout << "protocols::abinitio::DomainAssemblyCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_analysis_InterfaceAnalyzerMoverCreator()
	// { protocols::analysis::InterfaceAnalyzerMoverCreator cr; std::cout << "protocols::analysis::InterfaceAnalyzerMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_analysis_LoopAnalyzerMoverCreator()
	// { protocols::analysis::LoopAnalyzerMoverCreator cr; std::cout << "protocols::analysis::LoopAnalyzerMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_antibody_AntibodyCDRGrafterCreator()
	// { protocols::antibody::AntibodyCDRGrafterCreator cr; std::cout << "protocols::antibody::AntibodyCDRGrafterCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_antibody_AntibodyNumberingConverterMoverCreator()
	// { protocols::antibody::AntibodyNumberingConverterMoverCreator cr; std::cout << "protocols::antibody::AntibodyNumberingConverterMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_antibody_constraints_CDRDihedralConstraintMoverCreator()
	// { protocols::antibody::constraints::CDRDihedralConstraintMoverCreator cr; std::cout << "protocols::antibody::constraints::CDRDihedralConstraintMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_antibody_constraints_ParatopeEpitopeSiteConstraintMoverCreator()
	// { protocols::antibody::constraints::ParatopeEpitopeSiteConstraintMoverCreator cr; std::cout << "protocols::antibody::constraints::ParatopeEpitopeSiteConstraintMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_antibody_constraints_ParatopeSiteConstraintMoverCreator()
	// { protocols::antibody::constraints::ParatopeSiteConstraintMoverCreator cr; std::cout << "protocols::antibody::constraints::ParatopeSiteConstraintMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_antibody_design_AntibodyDesignMoverCreator()
	// { protocols::antibody::design::AntibodyDesignMoverCreator cr; std::cout << "protocols::antibody::design::AntibodyDesignMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_antibody_design_AntibodyDesignProtocolCreator()
	// { protocols::antibody::design::AntibodyDesignProtocolCreator cr; std::cout << "protocols::antibody::design::AntibodyDesignProtocolCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_backrub_BackrubMoverCreator()
	// { protocols::backrub::BackrubMoverCreator cr; std::cout << "protocols::backrub::BackrubMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_backrub_BackrubProtocolCreator()
	// { protocols::backrub::BackrubProtocolCreator cr; std::cout << "protocols::backrub::BackrubProtocolCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_backrub_BackrubSidechainMoverCreator()
	// { protocols::backrub::BackrubSidechainMoverCreator cr; std::cout << "protocols::backrub::BackrubSidechainMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_canonical_sampling_MetricRecorderCreator()
	// { protocols::canonical_sampling::MetricRecorderCreator cr; std::cout << "protocols::canonical_sampling::MetricRecorderCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_canonical_sampling_MetropolisHastingsMoverCreator()
	// { protocols::canonical_sampling::MetropolisHastingsMoverCreator cr; std::cout << "protocols::canonical_sampling::MetropolisHastingsMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_canonical_sampling_ParallelTemperingCreator()
	// { protocols::canonical_sampling::ParallelTemperingCreator cr; std::cout << "protocols::canonical_sampling::ParallelTemperingCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_canonical_sampling_HamiltonianExchangeCreator()
	// { protocols::canonical_sampling::HamiltonianExchangeCreator cr; std::cout << "protocols::canonical_sampling::HamiltonianExchangeCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_canonical_sampling_PDBTrajectoryRecorderCreator()
	// { protocols::canonical_sampling::PDBTrajectoryRecorderCreator cr; std::cout << "protocols::canonical_sampling::PDBTrajectoryRecorderCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_canonical_sampling_SilentTrajectoryRecorderCreator()
	// { protocols::canonical_sampling::SilentTrajectoryRecorderCreator cr; std::cout << "protocols::canonical_sampling::SilentTrajectoryRecorderCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_canonical_sampling_SimulatedTemperingCreator()
	// { protocols::canonical_sampling::SimulatedTemperingCreator cr; std::cout << "protocols::canonical_sampling::SimulatedTemperingCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_canonical_sampling_TrialCounterObserverCreator()
	// { protocols::canonical_sampling::TrialCounterObserverCreator cr; std::cout << "protocols::canonical_sampling::TrialCounterObserverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_constraint_generator_AddConstraintsCreator()
	// { protocols::constraint_generator::AddConstraintsCreator cr; std::cout << "protocols::constraint_generator::AddConstraintsCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_constraint_generator_RemoveConstraintsCreator()
	// { protocols::constraint_generator::RemoveConstraintsCreator cr; std::cout << "protocols::constraint_generator::RemoveConstraintsCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_carbohydrates_LinkageConformerMoverCreator()
	// { protocols::carbohydrates::LinkageConformerMoverCreator cr; std::cout << "protocols::carbohydrates::LinkageConformerMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_carbohydrates_GlycanRelaxMoverCreator()
	// { protocols::carbohydrates::GlycanRelaxMoverCreator cr; std::cout << "protocols::carbohydrates::GlycanRelaxMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_carbohydrates_SimpleGlycosylateMoverCreator()
	// { protocols::carbohydrates::SimpleGlycosylateMoverCreator cr; std::cout << "protocols::carbohydrates::SimpleGlycosylateMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_cryst_ReportGradientsMoverCreator()
	// { protocols::cryst::ReportGradientsMoverCreator cr; std::cout << "protocols::cryst::ReportGradientsMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_cryst_SetCrystWeightMoverCreator()
	// { protocols::cryst::SetCrystWeightMoverCreator cr; std::cout << "protocols::cryst::SetCrystWeightMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_cryst_RecomputeDensityMapMoverCreator()
	// { protocols::cryst::RecomputeDensityMapMoverCreator cr; std::cout << "protocols::cryst::RecomputeDensityMapMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_cryst_LoadDensityMapMoverCreator()
	// { protocols::cryst::LoadDensityMapMoverCreator cr; std::cout << "protocols::cryst::LoadDensityMapMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_cryst_FitBfactorsMoverCreator()
	// { protocols::cryst::FitBfactorsMoverCreator cr; std::cout << "protocols::cryst::FitBfactorsMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_cryst_UpdateSolventMoverCreator()
	// { protocols::cryst::UpdateSolventMoverCreator cr; std::cout << "protocols::cryst::UpdateSolventMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_cryst_TagPoseWithRefinementStatsMoverCreator()
	// { protocols::cryst::TagPoseWithRefinementStatsMoverCreator cr; std::cout << "protocols::cryst::TagPoseWithRefinementStatsMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_cryst_SetRefinementOptionsMoverCreator()
	// { protocols::cryst::SetRefinementOptionsMoverCreator cr; std::cout << "protocols::cryst::SetRefinementOptionsMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_cryst_UpdateCrystInfoCreator()
	// { protocols::cryst::UpdateCrystInfoCreator cr; std::cout << "protocols::cryst::UpdateCrystInfoCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_cryst_DockLatticeMoverCreator()
	// { protocols::cryst::DockLatticeMoverCreator cr; std::cout << "protocols::cryst::DockLatticeMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_cryst_MakeLatticeMoverCreator()
	// { protocols::cryst::MakeLatticeMoverCreator cr; std::cout << "protocols::cryst::MakeLatticeMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_cryst_MakeLayerMoverCreator()
	// { protocols::cryst::MakeLayerMoverCreator cr; std::cout << "protocols::cryst::MakeLayerMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_comparative_modeling_LoopRelaxMoverCreator()
	// { protocols::comparative_modeling::LoopRelaxMoverCreator cr; std::cout << "protocols::comparative_modeling::LoopRelaxMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_hybridization_HybridizeProtocolCreator()
	// { protocols::hybridization::HybridizeProtocolCreator cr; std::cout << "protocols::hybridization::HybridizeProtocolCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_hybridization_CartesianSamplerCreator()
	// { protocols::hybridization::CartesianSamplerCreator cr; std::cout << "protocols::hybridization::CartesianSamplerCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_hybridization_BackboneTorsionSamplerCreator()
	// { protocols::hybridization::BackboneTorsionSamplerCreator cr; std::cout << "protocols::hybridization::BackboneTorsionSamplerCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_hybridization_BackboneTorsionPerturbationCreator()
	// { protocols::hybridization::BackboneTorsionPerturbationCreator cr; std::cout << "protocols::hybridization::BackboneTorsionPerturbationCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_cyclic_peptide_CreateAngleConstraintCreator()
	// { protocols::cyclic_peptide::CreateAngleConstraintCreator cr; std::cout << "protocols::cyclic_peptide::CreateAngleConstraintCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_cyclic_peptide_CreateDistanceConstraintCreator()
	// { protocols::cyclic_peptide::CreateDistanceConstraintCreator cr; std::cout << "protocols::cyclic_peptide::CreateDistanceConstraintCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_cyclic_peptide_CreateTorsionConstraintCreator()
	// { protocols::cyclic_peptide::CreateTorsionConstraintCreator cr; std::cout << "protocols::cyclic_peptide::CreateTorsionConstraintCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_cyclic_peptide_DeclareBondCreator()
	// { protocols::cyclic_peptide::DeclareBondCreator cr; std::cout << "protocols::cyclic_peptide::DeclareBondCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_cyclic_peptide_PeptideStubMoverCreator()
	// { protocols::cyclic_peptide::PeptideStubMoverCreator cr; std::cout << "protocols::cyclic_peptide::PeptideStubMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_cyclic_peptide_TryDisulfPermutationsCreator()
	// { protocols::cyclic_peptide::TryDisulfPermutationsCreator cr; std::cout << "protocols::cyclic_peptide::TryDisulfPermutationsCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_cyclic_peptide_FlipChiralityMoverCreator()
	// { protocols::cyclic_peptide::FlipChiralityMoverCreator cr; std::cout << "protocols::cyclic_peptide::FlipChiralityMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_contact_map_ContactMapCreator()
	// { protocols::contact_map::ContactMapCreator cr; std::cout << "protocols::contact_map::ContactMapCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_denovo_design_DisulfidizeMoverCreator()
	// { protocols::denovo_design::DisulfidizeMoverCreator cr; std::cout << "protocols::denovo_design::DisulfidizeMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_denovo_design_movers_AddSegmentDataMoverCreator()
	// { protocols::denovo_design::movers::AddSegmentDataMoverCreator cr; std::cout << "protocols::denovo_design::movers::AddSegmentDataMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_denovo_design_movers_AlignResiduesMoverCreator()
	// { protocols::denovo_design::movers::AlignResiduesMoverCreator cr; std::cout << "protocols::denovo_design::movers::AlignResiduesMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_denovo_design_movers_BridgeChainsMoverCreator()
	// { protocols::denovo_design::movers::BridgeChainsMoverCreator cr; std::cout << "protocols::denovo_design::movers::BridgeChainsMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_denovo_design_movers_BridgeChainsCreator()
	// { protocols::denovo_design::movers::BridgeChainsCreator cr; std::cout << "protocols::denovo_design::movers::BridgeChainsCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_denovo_design_movers_BuildDeNovoBackboneMoverCreator()
	// { protocols::denovo_design::movers::BuildDeNovoBackboneMoverCreator cr; std::cout << "protocols::denovo_design::movers::BuildDeNovoBackboneMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_denovo_design_movers_ExtendChainMoverCreator()
	// { protocols::denovo_design::movers::ExtendChainMoverCreator cr; std::cout << "protocols::denovo_design::movers::ExtendChainMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_denovo_design_movers_FastDesignCreator()
	// { protocols::denovo_design::movers::FastDesignCreator cr; std::cout << "protocols::denovo_design::movers::FastDesignCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_denovo_design_movers_SetResidueAliasMoverCreator()
	// { protocols::denovo_design::movers::SetResidueAliasMoverCreator cr; std::cout << "protocols::denovo_design::movers::SetResidueAliasMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_denovo_design_movers_MakeAsymmetricStructureDataMoverCreator()
	// { protocols::denovo_design::movers::MakeAsymmetricStructureDataMoverCreator cr; std::cout << "protocols::denovo_design::movers::MakeAsymmetricStructureDataMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_design_opt_GreedyOptMutationMoverCreator()
	// { protocols::design_opt::GreedyOptMutationMoverCreator cr; std::cout << "protocols::design_opt::GreedyOptMutationMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_dna_DesignProteinBackboneAroundDNACreator()
	// { protocols::dna::DesignProteinBackboneAroundDNACreator cr; std::cout << "protocols::dna::DesignProteinBackboneAroundDNACreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_dna_DnaInterfaceMinMoverCreator()
	// { protocols::dna::DnaInterfaceMinMoverCreator cr; std::cout << "protocols::dna::DnaInterfaceMinMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_dna_DnaInterfaceMultiStateDesignCreator()
	// { protocols::dna::DnaInterfaceMultiStateDesignCreator cr; std::cout << "protocols::dna::DnaInterfaceMultiStateDesignCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_dna_DnaInterfacePackerCreator()
	// { protocols::dna::DnaInterfacePackerCreator cr; std::cout << "protocols::dna::DnaInterfacePackerCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_dna_SeparateDnaFromNonDnaCreator()
	// { protocols::dna::SeparateDnaFromNonDnaCreator cr; std::cout << "protocols::dna::SeparateDnaFromNonDnaCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_docking_ConformerSwitchMoverCreator()
	// { protocols::docking::ConformerSwitchMoverCreator cr; std::cout << "protocols::docking::ConformerSwitchMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_docking_DockingProtocolCreator()
	// { protocols::docking::DockingProtocolCreator cr; std::cout << "protocols::docking::DockingProtocolCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_docking_DockSetupMoverCreator()
	// { protocols::docking::DockSetupMoverCreator cr; std::cout << "protocols::docking::DockSetupMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_docking_DockingInitialPerturbationCreator()
	// { protocols::docking::DockingInitialPerturbationCreator cr; std::cout << "protocols::docking::DockingInitialPerturbationCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_generalized_kinematic_closure_GeneralizedKICCreator()
	// { protocols::generalized_kinematic_closure::GeneralizedKICCreator cr; std::cout << "protocols::generalized_kinematic_closure::GeneralizedKICCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_helical_bundle_BackboneGridSamplerCreator()
	// { protocols::helical_bundle::BackboneGridSamplerCreator cr; std::cout << "protocols::helical_bundle::BackboneGridSamplerCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_helical_bundle_BundleGridSamplerCreator()
	// { protocols::helical_bundle::BundleGridSamplerCreator cr; std::cout << "protocols::helical_bundle::BundleGridSamplerCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_helical_bundle_FitSimpleHelixCreator()
	// { protocols::helical_bundle::FitSimpleHelixCreator cr; std::cout << "protocols::helical_bundle::FitSimpleHelixCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_helical_bundle_MakeBundleCreator()
	// { protocols::helical_bundle::MakeBundleCreator cr; std::cout << "protocols::helical_bundle::MakeBundleCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_helical_bundle_MakeBundleHelixCreator()
	// { protocols::helical_bundle::MakeBundleHelixCreator cr; std::cout << "protocols::helical_bundle::MakeBundleHelixCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_helical_bundle_PerturbBundleCreator()
	// { protocols::helical_bundle::PerturbBundleCreator cr; std::cout << "protocols::helical_bundle::PerturbBundleCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_helical_bundle_PerturbBundleHelixCreator()
	// { protocols::helical_bundle::PerturbBundleHelixCreator cr; std::cout << "protocols::helical_bundle::PerturbBundleHelixCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_ncbb_SecStructMinimizeMoverCreator()
	// { protocols::ncbb::SecStructMinimizeMoverCreator cr; std::cout << "protocols::ncbb::SecStructMinimizeMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_ncbb_NcbbDockDesignProtocolCreator()
	// { protocols::ncbb::NcbbDockDesignProtocolCreator cr; std::cout << "protocols::ncbb::NcbbDockDesignProtocolCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_ncbb_oop_OopDockDesignProtocolCreator()
	// { protocols::ncbb::oop::OopDockDesignProtocolCreator cr; std::cout << "protocols::ncbb::oop::OopDockDesignProtocolCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_ncbb_oop_OopCreatorMoverCreator()
	// { protocols::ncbb::oop::OopCreatorMoverCreator cr; std::cout << "protocols::ncbb::oop::OopCreatorMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_symmetric_docking_SymDockProtocolCreator()
	// { protocols::symmetric_docking::SymDockProtocolCreator cr; std::cout << "protocols::symmetric_docking::SymDockProtocolCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_symmetric_docking_SymFoldandDockMoveRbJumpMoverCreator()
	// { protocols::symmetric_docking::SymFoldandDockMoveRbJumpMoverCreator cr; std::cout << "protocols::symmetric_docking::SymFoldandDockMoveRbJumpMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_symmetric_docking_SymFoldandDockSlideTrialMoverCreator()
	// { protocols::symmetric_docking::SymFoldandDockSlideTrialMoverCreator cr; std::cout << "protocols::symmetric_docking::SymFoldandDockSlideTrialMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_symmetric_docking_SymFoldandDockRbTrialMoverCreator()
	// { protocols::symmetric_docking::SymFoldandDockRbTrialMoverCreator cr; std::cout << "protocols::symmetric_docking::SymFoldandDockRbTrialMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_electron_density_SetupForDensityScoringMoverCreator()
	// { protocols::electron_density::SetupForDensityScoringMoverCreator cr; std::cout << "protocols::electron_density::SetupForDensityScoringMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_electron_density_BfactorFittingMoverCreator()
	// { protocols::electron_density::BfactorFittingMoverCreator cr; std::cout << "protocols::electron_density::BfactorFittingMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_electron_density_ScaleMapIntensitiesCreator()
	// { protocols::electron_density::ScaleMapIntensitiesCreator cr; std::cout << "protocols::electron_density::ScaleMapIntensitiesCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_electron_density_VoxelSpacingRefinementMoverCreator()
	// { protocols::electron_density::VoxelSpacingRefinementMoverCreator cr; std::cout << "protocols::electron_density::VoxelSpacingRefinementMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_electron_density_ReportFSCCreator()
	// { protocols::electron_density::ReportFSCCreator cr; std::cout << "protocols::electron_density::ReportFSCCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_enzdes_AddOrRemoveMatchCstsCreator()
	// { protocols::enzdes::AddOrRemoveMatchCstsCreator cr; std::cout << "protocols::enzdes::AddOrRemoveMatchCstsCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_enzdes_BackboneSamplerCreator()
	// { protocols::enzdes::BackboneSamplerCreator cr; std::cout << "protocols::enzdes::BackboneSamplerCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_enzdes_EnzRepackMinimizeCreator()
	// { protocols::enzdes::EnzRepackMinimizeCreator cr; std::cout << "protocols::enzdes::EnzRepackMinimizeCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_enzdes_PackRotamersMoverPartGreedyCreator()
	// { protocols::enzdes::PackRotamersMoverPartGreedyCreator cr; std::cout << "protocols::enzdes::PackRotamersMoverPartGreedyCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_enzdes_PredesignPerturbMoverCreator()
	// { protocols::enzdes::PredesignPerturbMoverCreator cr; std::cout << "protocols::enzdes::PredesignPerturbMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_enzdes_UpdateEnzdesHeaderMoverCreator()
	// { protocols::enzdes::UpdateEnzdesHeaderMoverCreator cr; std::cout << "protocols::enzdes::UpdateEnzdesHeaderMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_farna_ErraserMinimizerMoverCreator()
	// { protocols::rna::movers::ErraserMinimizerMoverCreator cr; std::cout << "protocols::rna::movers::ErraserMinimizerMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_farna_RNAIdealizeMoverCreator()
	// { protocols::rna::movers::RNAIdealizeMoverCreator cr; std::cout << "protocols::rna::movers::RNAIdealizeMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_features_ReportToDBCreator()
	// { protocols::features::ReportToDBCreator cr; std::cout << "protocols::features::ReportToDBCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_features_TrajectoryReportToDBCreator()
	// { protocols::features::TrajectoryReportToDBCreator cr; std::cout << "protocols::features::TrajectoryReportToDBCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_fldsgn_BluePrintBDRCreator()
	// { protocols::fldsgn::BluePrintBDRCreator cr; std::cout << "protocols::fldsgn::BluePrintBDRCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_fldsgn_CircularPermutationCreator()
	// { protocols::fldsgn::CircularPermutationCreator cr; std::cout << "protocols::fldsgn::CircularPermutationCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_fldsgn_MatchResiduesMoverCreator()
	// { protocols::fldsgn::MatchResiduesMoverCreator cr; std::cout << "protocols::fldsgn::MatchResiduesMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_fldsgn_SheetRemodelConstraintGeneratorCreator()
	// { protocols::fldsgn::SheetRemodelConstraintGeneratorCreator cr; std::cout << "protocols::fldsgn::SheetRemodelConstraintGeneratorCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_fldsgn_potentials_SetAACompositionPotentialCreator()
	// { protocols::fldsgn::potentials::SetAACompositionPotentialCreator cr; std::cout << "protocols::fldsgn::potentials::SetAACompositionPotentialCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_fldsgn_potentials_SetSecStructEnergiesCreator()
	// { protocols::fldsgn::potentials::SetSecStructEnergiesCreator cr; std::cout << "protocols::fldsgn::potentials::SetSecStructEnergiesCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_flexpep_docking_FlexPepDockingProtocolCreator()
	// { protocols::flexpep_docking::FlexPepDockingProtocolCreator cr; std::cout << "protocols::flexpep_docking::FlexPepDockingProtocolCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_flxbb_FlxbbDesignCreator()
	// { protocols::flxbb::FlxbbDesignCreator cr; std::cout << "protocols::flxbb::FlxbbDesignCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_flxbb_InterlockAromaCreator()
	// { protocols::flxbb::InterlockAromaCreator cr; std::cout << "protocols::flxbb::InterlockAromaCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_forge_constraints_InverseRotamersCstGeneratorCreator()
	// { protocols::forge::constraints::InverseRotamersCstGeneratorCreator cr; std::cout << "protocols::forge::constraints::InverseRotamersCstGeneratorCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_forge_constraints_InvrotTreeCstGeneratorCreator()
	// { protocols::forge::constraints::InvrotTreeCstGeneratorCreator cr; std::cout << "protocols::forge::constraints::InvrotTreeCstGeneratorCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_forge_constraints_NtoCConstraintGeneratorCreator()
	// { protocols::forge::constraints::NtoCConstraintGeneratorCreator cr; std::cout << "protocols::forge::constraints::NtoCConstraintGeneratorCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_forge_constraints_RemoveRemodelCstsCreator()
	// { protocols::forge::constraints::RemoveRemodelCstsCreator cr; std::cout << "protocols::forge::constraints::RemoveRemodelCstsCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_forge_remodel_ResidueVicinityCstGeneratorCreator()
	// { protocols::forge::remodel::ResidueVicinityCstGeneratorCreator cr; std::cout << "protocols::forge::remodel::ResidueVicinityCstGeneratorCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_forge_remodel_RemodelMoverCreator()
	// { protocols::forge::remodel::RemodelMoverCreator cr; std::cout << "protocols::forge::remodel::RemodelMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_grafting_AnchoredGraftMoverCreator()
	// { protocols::grafting::AnchoredGraftMoverCreator cr; std::cout << "protocols::grafting::AnchoredGraftMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_grafting_CCDEndsGraftMoverCreator()
	// { protocols::grafting::CCDEndsGraftMoverCreator cr; std::cout << "protocols::grafting::CCDEndsGraftMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_grafting_simple_movers_DeleteRegionMoverCreator()
	// { protocols::grafting::simple_movers::DeleteRegionMoverCreator cr; std::cout << "protocols::grafting::simple_movers::DeleteRegionMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_grafting_simple_movers_InsertPoseIntoPoseMoverCreator()
	// { protocols::grafting::simple_movers::InsertPoseIntoPoseMoverCreator cr; std::cout << "protocols::grafting::simple_movers::InsertPoseIntoPoseMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_grafting_simple_movers_ReplaceRegionMoverCreator()
	// { protocols::grafting::simple_movers::ReplaceRegionMoverCreator cr; std::cout << "protocols::grafting::simple_movers::ReplaceRegionMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_grafting_simple_movers_KeepRegionMoverCreator()
	// { protocols::grafting::simple_movers::KeepRegionMoverCreator cr; std::cout << "protocols::grafting::simple_movers::KeepRegionMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_idealize_IdealizeMoverCreator()
	// { protocols::idealize::IdealizeMoverCreator cr; std::cout << "protocols::idealize::IdealizeMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_hotspot_hashing_movers_PlaceSurfaceProbeCreator()
	// { protocols::hotspot_hashing::movers::PlaceSurfaceProbeCreator cr; std::cout << "protocols::hotspot_hashing::movers::PlaceSurfaceProbeCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_kinematic_closure_KicMoverCreator()
	// { protocols::kinematic_closure::KicMoverCreator cr; std::cout << "protocols::kinematic_closure::KicMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_ligand_docking_AddHydrogensCreator()
	// { protocols::ligand_docking::AddHydrogensCreator cr; std::cout << "protocols::ligand_docking::AddHydrogensCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_ligand_docking_CompoundTranslateCreator()
	// { protocols::ligand_docking::CompoundTranslateCreator cr; std::cout << "protocols::ligand_docking::CompoundTranslateCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_ligand_docking_ComputeLigandRDFCreator()
	// { protocols::ligand_docking::ComputeLigandRDFCreator cr; std::cout << "protocols::ligand_docking::ComputeLigandRDFCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_ligand_docking_FinalMinimizerCreator()
	// { protocols::ligand_docking::FinalMinimizerCreator cr; std::cout << "protocols::ligand_docking::FinalMinimizerCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_ligand_docking_GrowLigandCreator()
	// { protocols::ligand_docking::GrowLigandCreator cr; std::cout << "protocols::ligand_docking::GrowLigandCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_ligand_docking_HighResDockerCreator()
	// { protocols::ligand_docking::HighResDockerCreator cr; std::cout << "protocols::ligand_docking::HighResDockerCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_ligand_docking_InterfaceScoreCalculatorCreator()
	// { protocols::ligand_docking::InterfaceScoreCalculatorCreator cr; std::cout << "protocols::ligand_docking::InterfaceScoreCalculatorCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_ligand_docking_LigandDesignCreator()
	// { protocols::ligand_docking::LigandDesignCreator cr; std::cout << "protocols::ligand_docking::LigandDesignCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_ligand_docking_MinimizeBackboneCreator()
	// { protocols::ligand_docking::MinimizeBackboneCreator cr; std::cout << "protocols::ligand_docking::MinimizeBackboneCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_ligand_docking_RandomConformersCreator()
	// { protocols::ligand_docking::RandomConformersCreator cr; std::cout << "protocols::ligand_docking::RandomConformersCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_ligand_docking_RotateCreator()
	// { protocols::ligand_docking::RotateCreator cr; std::cout << "protocols::ligand_docking::RotateCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_ligand_docking_RotatesCreator()
	// { protocols::ligand_docking::RotatesCreator cr; std::cout << "protocols::ligand_docking::RotatesCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_ligand_docking_SlideTogetherCreator()
	// { protocols::ligand_docking::SlideTogetherCreator cr; std::cout << "protocols::ligand_docking::SlideTogetherCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_ligand_docking_StartFromCreator()
	// { protocols::ligand_docking::StartFromCreator cr; std::cout << "protocols::ligand_docking::StartFromCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_ligand_docking_TransformCreator()
	// { protocols::ligand_docking::TransformCreator cr; std::cout << "protocols::ligand_docking::TransformCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_ligand_docking_TranslateCreator()
	// { protocols::ligand_docking::TranslateCreator cr; std::cout << "protocols::ligand_docking::TranslateCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_ligand_docking_WriteLigandMolFileCreator()
	// { protocols::ligand_docking::WriteLigandMolFileCreator cr; std::cout << "protocols::ligand_docking::WriteLigandMolFileCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_loop_build_LoopmodelWrapperCreator()
	// { protocols::loop_build::LoopmodelWrapperCreator cr; std::cout << "protocols::loop_build::LoopmodelWrapperCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_loop_build_LoopMover_SlidingWindowCreator()
	// { protocols::loop_build::LoopMover_SlidingWindowCreator cr; std::cout << "protocols::loop_build::LoopMover_SlidingWindowCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_loop_modeler_LoopModelerCreator()
	// { protocols::loop_modeler::LoopModelerCreator cr; std::cout << "protocols::loop_modeling::LoopModelerCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_loop_modeling_LoopProtocolCreator()
	// { protocols::loop_modeling::LoopProtocolCreator cr; std::cout << "protocols::loop_modeling::LoopProtocolCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_loop_modeling_LoopBuilderCreator()
	// { protocols::loop_modeling::LoopBuilderCreator cr; std::cout << "protocols::loop_modeling::LoopBuilderCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_loop_modeling_refiners_MinimizationRefinerCreator()
	// { protocols::loop_modeling::refiners::MinimizationRefinerCreator cr; std::cout << "protocols::loop_modeling::refiners::MinimizationRefinerCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_loop_modeling_refiners_RepackingRefinerCreator()
	// { protocols::loop_modeling::refiners::RepackingRefinerCreator cr; std::cout << "protocols::loop_modeling::refiners::RepackingRefinerCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_loop_modeling_refiners_RotamerTrialsRefinerCreator()
	// { protocols::loop_modeling::refiners::RotamerTrialsRefinerCreator cr; std::cout << "protocols::loop_modeling::refiners::RotamerTrialsRefinerCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_loop_modeling_samplers_LegacyKicSamplerCreator()
	// { protocols::loop_modeling::samplers::LegacyKicSamplerCreator cr; std::cout << "protocols::loop_modeling::samplers::LegacyKicSamplerCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_loop_modeling_utilities_PrepareForCentroidCreator()
	// { protocols::loop_modeling::utilities::PrepareForCentroidCreator cr; std::cout << "protocols::loop_modeling::utilities::PrepareForCentroidCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_loop_modeling_utilities_PrepareForFullatomCreator()
	// { protocols::loop_modeling::utilities::PrepareForFullatomCreator cr; std::cout << "protocols::loop_modeling::utilities::PrepareForFullatomCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_loophash_LoopHashMoverWrapperCreator()
	// { protocols::loophash::LoopHashMoverWrapperCreator cr; std::cout << "protocols::loophash::LoopHashMoverWrapperCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_loophash_LoopHashDiversifierCreator()
	// { protocols::loophash::LoopHashDiversifierCreator cr; std::cout << "protocols::loophash::LoopHashDiversifierCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_loops_FoldTreeFromLoopsCreator()
	// { protocols::loops::FoldTreeFromLoopsCreator cr; std::cout << "protocols::loops::FoldTreeFromLoopsCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_loops_loop_closure_ccd_CCDLoopClosureMoverCreator()
	// { protocols::loops::loop_closure::ccd::CCDLoopClosureMoverCreator cr; std::cout << "protocols::loops::loop_closure::ccd::CCDLoopClosureMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_loops_loop_mover_perturb_LoopMover_Perturb_CCDCreator()
	// { protocols::loops::loop_mover::perturb::LoopMover_Perturb_CCDCreator cr; std::cout << "protocols::loops::loop_mover::perturb::LoopMover_Perturb_CCDCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_loops_loop_mover_perturb_LoopMover_Perturb_KICCreator()
	// { protocols::loops::loop_mover::perturb::LoopMover_Perturb_KICCreator cr; std::cout << "protocols::loops::loop_mover::perturb::LoopMover_Perturb_KICCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_loops_loop_mover_perturb_LoopMover_Perturb_QuickCCDCreator()
	// { protocols::loops::loop_mover::perturb::LoopMover_Perturb_QuickCCDCreator cr; std::cout << "protocols::loops::loop_mover::perturb::LoopMover_Perturb_QuickCCDCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_loops_loop_mover_perturb_LoopMover_Perturb_QuickCCD_MovesCreator()
	// { protocols::loops::loop_mover::perturb::LoopMover_Perturb_QuickCCD_MovesCreator cr; std::cout << "protocols::loops::loop_mover::perturb::LoopMover_Perturb_QuickCCD_MovesCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_loops_loop_mover_refine_LoopMover_Refine_BackrubCreator()
	// { protocols::loops::loop_mover::refine::LoopMover_Refine_BackrubCreator cr; std::cout << "protocols::loops::loop_mover::refine::LoopMover_Refine_BackrubCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_loops_loop_mover_refine_LoopMover_Refine_CCDCreator()
	// { protocols::loops::loop_mover::refine::LoopMover_Refine_CCDCreator cr; std::cout << "protocols::loops::loop_mover::refine::LoopMover_Refine_CCDCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_loops_loop_mover_refine_LoopMover_Refine_KICCreator()
	// { protocols::loops::loop_mover::refine::LoopMover_Refine_KICCreator cr; std::cout << "protocols::loops::loop_mover::refine::LoopMover_Refine_KICCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_loops_loop_mover_refine_LoopRefineInnerCycleContainerCreator()
	// { protocols::loops::loop_mover::refine::LoopRefineInnerCycleContainerCreator cr; std::cout << "protocols::loops::loop_mover::refine::LoopRefineInnerCycleContainerCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_loops_loop_mover_refine_RepackTrialCreator()
	// { protocols::loops::loop_mover::refine::RepackTrialCreator cr; std::cout << "protocols::loops::loop_mover::refine::RepackTrialCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_loops_loop_mover_refine_ShearMinCCDTrialCreator()
	// { protocols::loops::loop_mover::refine::ShearMinCCDTrialCreator cr; std::cout << "protocols::loops::loop_mover::refine::ShearMinCCDTrialCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_loops_loop_mover_refine_SmallMinCCDTrialCreator()
	// { protocols::loops::loop_mover::refine::SmallMinCCDTrialCreator cr; std::cout << "protocols::loops::loop_mover::refine::SmallMinCCDTrialCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_loops_loop_mover_LoopCMCreator()
	// { protocols::loops::loop_mover::LoopCMCreator cr; std::cout << "protocols::loops::loop_mover::LoopCMCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_match_MatcherMoverCreator()
	// { protocols::match::MatcherMoverCreator cr; std::cout << "protocols::match::MatcherMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_matdes_ExtractSubposeMoverCreator()
	// { protocols::matdes::ExtractSubposeMoverCreator cr; std::cout << "protocols::matdes::ExtractSubposeMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_matdes_MatDesGreedyOptMutationMoverCreator()
	// { protocols::matdes::MatDesGreedyOptMutationMoverCreator cr; std::cout << "protocols::matdes::MatDesGreedyOptMutationMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_matdes_SymDofMoverCreator()
	// { protocols::matdes::SymDofMoverCreator cr; std::cout << "protocols::matdes::SymDofMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_matdes_SchemePlaceMotifsMoverCreator()
	// { protocols::matdes::SchemePlaceMotifsMoverCreator cr; std::cout << "protocols::matdes::SchemePlaceMotifsMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_md_CartesianMDCreator()
	// { protocols::md::CartesianMDCreator cr; std::cout << "protocols::md::CartesianMDCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_BBGaussianMoverCreator()
	// { protocols::simple_moves::BBGaussianMoverCreator cr; std::cout << "protocols::simple_moves::BBGaussianMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_PeriodicBoxMoverCreator()
	// { protocols::simple_moves::PeriodicBoxMoverCreator cr; std::cout << "protocols::simple_moves::PeriodicBoxMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_motifs_MotifDnaPackerCreator()
	// { protocols::motifs::MotifDnaPackerCreator cr; std::cout << "protocols::motifs::MotifDnaPackerCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_motif_grafting_movers_MotifGraftCreator()
	// { protocols::motif_grafting::movers::MotifGraftCreator cr; std::cout << "protocols::motif_grafting::movers::MotifGraftCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_moves_DsspMoverCreator()
	// { protocols::moves::DsspMoverCreator cr; std::cout << "protocols::moves::DsspMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_moves_IfMoverCreator()
	// { protocols::moves::IfMoverCreator cr; std::cout << "protocols::moves::IfMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_moves_RandomMoverCreator()
	// { protocols::moves::RandomMoverCreator cr; std::cout << "protocols::moves::RandomMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_moves_IteratedConvergenceMoverCreator()
	// { protocols::moves::IteratedConvergenceMoverCreator cr; std::cout << "protocols::moves::IteratedConvergenceMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_moves_PyMOLMoverCreator()
	// { protocols::moves::PyMOLMoverCreator cr; std::cout << "protocols::moves::PyMOLMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_moves_RampingMoverCreator()
	// { protocols::moves::RampingMoverCreator cr; std::cout << "protocols::moves::RampingMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_moves_FilterReportAsPoseExtraScoresMoverCreator()
	// { protocols::moves::FilterReportAsPoseExtraScoresMoverCreator cr; std::cout << "protocols::moves::FilterReportAsPoseExtraScoresMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_MonteCarloResetCreator()
	// { protocols::monte_carlo::MonteCarloResetCreator cr; std::cout << "protocols::monte_carlo::MonteCarloResetCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_AlignChainMoverCreator()
	// { protocols::simple_moves::AlignChainMoverCreator cr; std::cout << "protocols::simple_moves::AlignChainMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_ResetBaselineMoverCreator()
	// { protocols::monte_carlo::ResetBaselineMoverCreator cr; std::cout << "protocols::monte_carlo::ResetBaselineMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_RingConformationMoverCreator()
	// { protocols::simple_moves::RingConformationMoverCreator cr; std::cout << "protocols::simple_moves::RingConformationMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_SimpleThreadingMoverCreator()
	// { protocols::simple_moves::SimpleThreadingMoverCreator cr; std::cout << "protocols::simple_moves::SimpleThreadingMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_nonlocal_SingleFragmentMoverCreator()
	// { protocols::nonlocal::SingleFragmentMoverCreator cr; std::cout << "protocols::nonlocal::SingleFragmentMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_normalmode_NormalModeRelaxMoverCreator()
	// { protocols::normalmode::NormalModeRelaxMoverCreator cr; std::cout << "protocols::normalmode::NormalModeRelaxMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_normalmode_NormalModeMinimizerCreator()
	// { protocols::normalmode::NormalModeMinimizerCreator cr; std::cout << "protocols::normalmode::NormalModeMinimizerCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_pb_potential_SetupPoissonBoltzmannPotentialCreator()
	// { protocols::pb_potential::SetupPoissonBoltzmannPotentialCreator cr; std::cout << "protocols::pb_potential::SetupPoissonBoltzmannPotentialCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_pose_length_moves_FixAllLoopsMoverCreator()
	// { protocols::pose_length_moves::FixAllLoopsMoverCreator cr; std::cout << "protocols::pose_length_moves::FixAllLoopsMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_pose_length_moves_InsertResMoverCreator()
	// { protocols::pose_length_moves::InsertResMoverCreator cr; std::cout << "protocols::pose_length_moves::InsertResMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_pose_length_moves_NearNativeLoopCloserCreator()
	// { protocols::pose_length_moves::NearNativeLoopCloserCreator cr; std::cout << "protocols::pose_length_moves::NearNativeLoopCloserCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_protein_interface_design_movers_AddChainBreakCreator()
	// { protocols::protein_interface_design::movers::AddChainBreakCreator cr; std::cout << "protocols::protein_interface_design::movers::AddChainBreakCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_protein_interface_design_movers_AddSidechainConstraintsToHotspotsCreator()
	// { protocols::protein_interface_design::movers::AddSidechainConstraintsToHotspotsCreator cr; std::cout << "protocols::protein_interface_design::movers::AddSidechainConstraintsToHotspotsCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_protein_interface_design_movers_BackrubDDMoverCreator()
	// { protocols::protein_interface_design::movers::BackrubDDMoverCreator cr; std::cout << "protocols::protein_interface_design::movers::BackrubDDMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_protein_interface_design_movers_BestHotspotCstMoverCreator()
	// { protocols::protein_interface_design::movers::BestHotspotCstMoverCreator cr; std::cout << "protocols::protein_interface_design::movers::BestHotspotCstMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_protein_interface_design_movers_BuildAlaPoseCreator()
	// { protocols::protein_interface_design::movers::BuildAlaPoseCreator cr; std::cout << "protocols::protein_interface_design::movers::BuildAlaPoseCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_protein_interface_design_movers_DesignMinimizeHbondsCreator()
	// { protocols::protein_interface_design::movers::DesignMinimizeHbondsCreator cr; std::cout << "protocols::protein_interface_design::movers::DesignMinimizeHbondsCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_protein_interface_design_movers_DisulfideMoverCreator()
	// { protocols::protein_interface_design::movers::DisulfideMoverCreator cr; std::cout << "protocols::protein_interface_design::movers::DisulfideMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_protein_interface_design_movers_DockAndRetrieveSidechainsCreator()
	// { protocols::protein_interface_design::movers::DockAndRetrieveSidechainsCreator cr; std::cout << "protocols::protein_interface_design::movers::DockAndRetrieveSidechainsCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_DumpPdbCreator()
	// { protocols::simple_moves::DumpPdbCreator cr; std::cout << "protocols::simple_moves::DumpPdbCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_protein_interface_design_movers_FavorNativeResiduePreCycleCreator()
	// { protocols::protein_interface_design::movers::FavorNativeResiduePreCycleCreator cr; std::cout << "protocols::protein_interface_design::movers::FavorNativeResiduePreCycleCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_protein_interface_design_movers_FavorNonNativeResiduePreCycleCreator()
	// { protocols::protein_interface_design::movers::FavorNonNativeResiduePreCycleCreator cr; std::cout << "protocols::protein_interface_design::movers::FavorNonNativeResiduePreCycleCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_protein_interface_design_movers_HotspotDisjointedFoldTreeMoverCreator()
	// { protocols::protein_interface_design::movers::HotspotDisjointedFoldTreeMoverCreator cr; std::cout << "protocols::protein_interface_design::movers::HotspotDisjointedFoldTreeMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_protein_interface_design_movers_HotspotHasherMoverCreator()
	// { protocols::protein_interface_design::movers::HotspotHasherMoverCreator cr; std::cout << "protocols::protein_interface_design::movers::HotspotHasherMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_protein_interface_design_movers_InterfaceRecapitulationMoverCreator()
	// { protocols::protein_interface_design::movers::InterfaceRecapitulationMoverCreator cr; std::cout << "protocols::protein_interface_design::movers::InterfaceRecapitulationMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_protein_interface_design_movers_LoopFinderCreator()
	// { protocols::protein_interface_design::movers::LoopFinderCreator cr; std::cout << "protocols::protein_interface_design::movers::LoopFinderCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_protein_interface_design_movers_LoopLengthChangeCreator()
	// { protocols::protein_interface_design::movers::LoopLengthChangeCreator cr; std::cout << "protocols::protein_interface_design::movers::LoopLengthChangeCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_protein_interface_design_movers_LoopMoverFromCommandLineCreator()
	// { protocols::protein_interface_design::movers::LoopMoverFromCommandLineCreator cr; std::cout << "protocols::protein_interface_design::movers::LoopMoverFromCommandLineCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_protein_interface_design_movers_LoopOverCreator()
	// { protocols::protein_interface_design::movers::LoopOverCreator cr; std::cout << "protocols::protein_interface_design::movers::LoopOverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_protein_interface_design_movers_LoopRemodelCreator()
	// { protocols::protein_interface_design::movers::LoopRemodelCreator cr; std::cout << "protocols::protein_interface_design::movers::LoopRemodelCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_protein_interface_design_movers_MapHotspotCreator()
	// { protocols::protein_interface_design::movers::MapHotspotCreator cr; std::cout << "protocols::protein_interface_design::movers::MapHotspotCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_evolution_NucleotideMutationCreator()
	// { protocols::evolution::NucleotideMutationCreator cr; std::cout << "protocols::evolution::NucleotideMutationCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_protein_interface_design_movers_PatchdockTransformCreator()
	// { protocols::protein_interface_design::movers::PatchdockTransformCreator cr; std::cout << "protocols::protein_interface_design::movers::PatchdockTransformCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_protein_interface_design_movers_PeptideStapleDesignMoverCreator()
	// { protocols::protein_interface_design::movers::PeptideStapleDesignMoverCreator cr; std::cout << "protocols::protein_interface_design::movers::PeptideStapleDesignMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_protein_interface_design_movers_PlacementAuctionMoverCreator()
	// { protocols::protein_interface_design::movers::PlacementAuctionMoverCreator cr; std::cout << "protocols::protein_interface_design::movers::PlacementAuctionMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_protein_interface_design_movers_PlacementMinimizationMoverCreator()
	// { protocols::protein_interface_design::movers::PlacementMinimizationMoverCreator cr; std::cout << "protocols::protein_interface_design::movers::PlacementMinimizationMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_protein_interface_design_movers_PlaceOnLoopCreator()
	// { protocols::protein_interface_design::movers::PlaceOnLoopCreator cr; std::cout << "protocols::protein_interface_design::movers::PlaceOnLoopCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_protein_interface_design_movers_PlaceSimultaneouslyMoverCreator()
	// { protocols::protein_interface_design::movers::PlaceSimultaneouslyMoverCreator cr; std::cout << "protocols::protein_interface_design::movers::PlaceSimultaneouslyMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_protein_interface_design_movers_PlaceStubMoverCreator()
	// { protocols::protein_interface_design::movers::PlaceStubMoverCreator cr; std::cout << "protocols::protein_interface_design::movers::PlaceStubMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_protein_interface_design_movers_PrepackMoverCreator()
	// { protocols::protein_interface_design::movers::PrepackMoverCreator cr; std::cout << "protocols::protein_interface_design::movers::PrepackMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_protein_interface_design_movers_ProteinInterfaceMultiStateDesignMoverCreator()
	// { protocols::protein_interface_design::movers::ProteinInterfaceMultiStateDesignMoverCreator cr; std::cout << "protocols::protein_interface_design::movers::ProteinInterfaceMultiStateDesignMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_protein_interface_design_movers_RandomMutationCreator()
	// { protocols::protein_interface_design::movers::RandomMutationCreator cr; std::cout << "protocols::protein_interface_design::movers::RandomMutationCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_protein_interface_design_movers_RepackMinimizeCreator()
	// { protocols::protein_interface_design::movers::RepackMinimizeCreator cr; std::cout << "protocols::protein_interface_design::movers::RepackMinimizeCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_protein_interface_design_movers_SaveAndRetrieveSidechainsCreator()
	// { protocols::protein_interface_design::movers::SaveAndRetrieveSidechainsCreator cr; std::cout << "protocols::protein_interface_design::movers::SaveAndRetrieveSidechainsCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_protein_interface_design_movers_SetAtomTreeCreator()
	// { protocols::protein_interface_design::movers::SetAtomTreeCreator cr; std::cout << "protocols::protein_interface_design::movers::SetAtomTreeCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_protein_interface_design_movers_SetTemperatureFactorCreator()
	// { protocols::protein_interface_design::movers::SetTemperatureFactorCreator cr; std::cout << "protocols::protein_interface_design::movers::SetTemperatureFactorCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_protein_interface_design_movers_SetupHotspotConstraintsMoverCreator()
	// { protocols::protein_interface_design::movers::SetupHotspotConstraintsMoverCreator cr; std::cout << "protocols::protein_interface_design::movers::SetupHotspotConstraintsMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_protein_interface_design_movers_SetupHotspotConstraintsLoopsMoverCreator()
	// { protocols::protein_interface_design::movers::SetupHotspotConstraintsLoopsMoverCreator cr; std::cout << "protocols::protein_interface_design::movers::SetupHotspotConstraintsLoopsMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_protein_interface_design_movers_SpinMoverCreator()
	// { protocols::protein_interface_design::movers::SpinMoverCreator cr; std::cout << "protocols::protein_interface_design::movers::SpinMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_protein_interface_design_movers_TaskAwareCstsCreator()
	// { protocols::protein_interface_design::movers::TaskAwareCstsCreator cr; std::cout << "protocols::protein_interface_design::movers::TaskAwareCstsCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_protein_interface_design_movers_SubroutineMoverCreator()
	// { protocols::protein_interface_design::movers::SubroutineMoverCreator cr; std::cout << "protocols::protein_interface_design::movers::SubroutineMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_protein_interface_design_movers_TryRotamersCreator()
	// { protocols::protein_interface_design::movers::TryRotamersCreator cr; std::cout << "protocols::protein_interface_design::movers::TryRotamersCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_protein_interface_design_movers_ShoveResidueMoverCreator()
	// { protocols::protein_interface_design::movers::ShoveResidueMoverCreator cr; std::cout << "protocols::protein_interface_design::movers::ShoveResidueMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_protein_interface_design_movers_VLBCreator()
	// { protocols::protein_interface_design::movers::VLBCreator cr; std::cout << "protocols::protein_interface_design::movers::VLBCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_protein_interface_design_movers_DockWithHotspotMoverCreator()
	// { protocols::protein_interface_design::movers::DockWithHotspotMoverCreator cr; std::cout << "protocols::protein_interface_design::movers::DockWithHotspotMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_protein_interface_design_movers_TopologyBrokerMoverCreator()
	// { protocols::protein_interface_design::movers::TopologyBrokerMoverCreator cr; std::cout << "protocols::protein_interface_design::movers::TopologyBrokerMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_qsar_RenderGridsToKinemageCreator()
	// { protocols::qsar::RenderGridsToKinemageCreator cr; std::cout << "protocols::qsar::RenderGridsToKinemageCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_rbsegment_relax_MakeStarTopologyMoverCreator()
	// { protocols::rbsegment_relax::MakeStarTopologyMoverCreator cr; std::cout << "protocols::rbsegment_relax::MakeStarTopologyMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_rbsegment_relax_OptimizeThreadingMoverCreator()
	// { protocols::rbsegment_relax::OptimizeThreadingMoverCreator cr; std::cout << "protocols::rbsegment_relax::OptimizeThreadingMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_rbsegment_relax_IdealizeHelicesMoverCreator()
	// { protocols::rbsegment_relax::IdealizeHelicesMoverCreator cr; std::cout << "protocols::rbsegment_relax::IdealizeHelicesMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_relax_AtomCoordinateCstMoverCreator()
	// { protocols::relax::AtomCoordinateCstMoverCreator cr; std::cout << "protocols::relax::AtomCoordinateCstMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_relax_FastRelaxCreator()
	// { protocols::relax::FastRelaxCreator cr; std::cout << "protocols::relax::FastRelaxCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_relax_LocalRelaxCreator()
	// { protocols::relax::LocalRelaxCreator cr; std::cout << "protocols::relax::LocalRelaxCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_residue_selectors_StoreResidueSubsetMoverCreator()
	// { protocols::residue_selectors::StoreResidueSubsetMoverCreator cr; std::cout << "protocols::residue_selectors::StoreResidueSubsetMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_rigid_RigidBodyPerturbNoCenterMoverCreator()
	// { protocols::rigid::RigidBodyPerturbNoCenterMoverCreator cr; std::cout << "protocols::rigid::RigidBodyPerturbNoCenterMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_rigid_RigidBodyTiltMoverCreator()
	// { protocols::rigid::RigidBodyTiltMoverCreator cr; std::cout << "protocols::rigid::RigidBodyTiltMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_rigid_RigidBodyTransMoverCreator()
	// { protocols::rigid::RigidBodyTransMoverCreator cr; std::cout << "protocols::rigid::RigidBodyTransMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_rigid_RollMoverCreator()
	// { protocols::rigid::RollMoverCreator cr; std::cout << "protocols::rigid::RollMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_rigid_UniformRigidBodyCMCreator()
	// { protocols::rigid::UniformRigidBodyCMCreator cr; std::cout << "protocols::rigid::UniformRigidBodyCMCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_rigid_UniformRigidBodyMoverCreator()
	// { protocols::rigid::UniformRigidBodyMoverCreator cr; std::cout << "protocols::rigid::UniformRigidBodyMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_rosetta_scripts_ParsedProtocolCreator()
	// { protocols::rosetta_scripts::ParsedProtocolCreator cr; std::cout << "protocols::rosetta_scripts::ParsedProtocolCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_rosetta_scripts_SavePoseMoverCreator()
	// { protocols::rosetta_scripts::SavePoseMoverCreator cr; std::cout << "protocols::rosetta_scripts::SavePoseMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_rosetta_scripts_MultiplePoseMoverCreator()
	// { protocols::rosetta_scripts::MultiplePoseMoverCreator cr; std::cout << "protocols::rosetta_scripts::MultiplePoseMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_rosetta_scripts_MultipleOutputWrapperCreator()
	// { protocols::rosetta_scripts::MultipleOutputWrapperCreator cr; std::cout << "protocols::rosetta_scripts::MultipleOutputWrapperCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_rotamer_recovery_RotamerRecoveryMoverCreator()
	// { protocols::rotamer_recovery::RotamerRecoveryMoverCreator cr; std::cout << "protocols::rotamer_recovery::RotamerRecoveryMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_seeded_abinitio_CAcstGeneratorCreator()
	// { protocols::seeded_abinitio::CAcstGeneratorCreator cr; std::cout << "protocols::seeded_abinitio::CAcstGeneratorCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_seeded_abinitio_CloseFoldCreator()
	// { protocols::seeded_abinitio::CloseFoldCreator cr; std::cout << "protocols::seeded_abinitio::CloseFoldCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_seeded_abinitio_CoordinateCstCreator()
	// { protocols::seeded_abinitio::CoordinateCstCreator cr; std::cout << "protocols::seeded_abinitio::CoordinateCstCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_seeded_abinitio_DefineMovableLoopsCreator()
	// { protocols::seeded_abinitio::DefineMovableLoopsCreator cr; std::cout << "protocols::seeded_abinitio::DefineMovableLoopsCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_seeded_abinitio_GrowPeptidesCreator()
	// { protocols::seeded_abinitio::GrowPeptidesCreator cr; std::cout << "protocols::seeded_abinitio::GrowPeptidesCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_seeded_abinitio_SeedFoldTreeCreator()
	// { protocols::seeded_abinitio::SeedFoldTreeCreator cr; std::cout << "protocols::seeded_abinitio::SeedFoldTreeCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_seeded_abinitio_SeedSetupMoverCreator()
	// { protocols::seeded_abinitio::SeedSetupMoverCreator cr; std::cout << "protocols::seeded_abinitio::SeedSetupMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_seeded_abinitio_SwapSegmentCreator()
	// { protocols::seeded_abinitio::SwapSegmentCreator cr; std::cout << "protocols::seeded_abinitio::SwapSegmentCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_seeded_abinitio_SegmentHybridizerCreator()
	// { protocols::seeded_abinitio::SegmentHybridizerCreator cr; std::cout << "protocols::seeded_abinitio::SegmentHybridizerCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_AddConstraintsToCurrentConformationMoverCreator()
	// { protocols::constraint_movers::AddConstraintsToCurrentConformationMoverCreator cr; std::cout << "protocols::simple_moves::AddConstraintsToCurrentConformationMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_AddChainMoverCreator()
	// { protocols::simple_moves::AddChainMoverCreator cr; std::cout << "protocols::simple_moves::AddChainMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_AddJobPairDataCreator()
	// { protocols::simple_moves::AddJobPairDataCreator cr; std::cout << "protocols::simple_moves::AddJobPairDataCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_ChangeAndResetFoldTreeMoverCreator()
	// { protocols::simple_moves::ChangeAndResetFoldTreeMoverCreator cr; std::cout << "protocols::simple_moves::ChangeAndResetFoldTreeMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_BoltzmannRotamerMoverCreator()
	// { protocols::minimization_packing::BoltzmannRotamerMoverCreator cr; std::cout << "protocols::minimization_packing::BoltzmannRotamerMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_ClearConstraintsMoverCreator()
	// { protocols::constraint_movers::ClearConstraintsMoverCreator cr; std::cout << "protocols::simple_moves::ClearConstraintsMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_ConsensusDesignMoverCreator()
	// { protocols::calc_taskop_movers::ConsensusDesignMoverCreator cr; std::cout << "protocols::calc_taskop_movers::ConsensusDesignMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_ConstraintSetMoverCreator()
	// { protocols::constraint_movers::ConstraintSetMoverCreator cr; std::cout << "protocols::constraint_movers::ConstraintSetMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_ContingentAcceptMoverCreator()
	// { protocols::simple_moves::ContingentAcceptMoverCreator cr; std::cout << "protocols::simple_moves::ContingentAcceptMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_CoupledMoverCreator()
	// { protocols::simple_moves::CoupledMoverCreator cr; std::cout << "protocols::simple_moves::CoupledMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_DeleteChainMoverCreator()
	// { protocols::simple_moves::DeleteChainMoverCreator cr; std::cout << "protocols::simple_moves::DeleteChainMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_DeleteChainsMoverCreator()
	// { protocols::simple_moves::DeleteChainsMoverCreator cr; std::cout << "protocols::simple_moves::DeleteChainsMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_ddGCreator()
	// { protocols::simple_ddg::ddGCreator cr; std::cout << "protocols::simple_ddg::ddGCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_DisulfideInsertionMoverCreator()
	// { protocols::simple_moves::DisulfideInsertionMoverCreator cr; std::cout << "protocols::simple_moves::DisulfideInsertionMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_evolution_EvolutionaryDynamicsMoverCreator()
	// { protocols::evolution::EvolutionaryDynamicsMoverCreator cr; std::cout << "protocols::evolution::EvolutionaryDynamicsMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_ExtendedPoseMoverCreator()
	// { protocols::pose_creation::ExtendedPoseMoverCreator cr; std::cout << "protocols::pose_creation::ExtendedPoseMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_FavorSequenceProfileCreator()
	// { protocols::simple_moves::FavorSequenceProfileCreator cr; std::cout << "protocols::simple_moves::FavorSequenceProfileCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_FavorSymmetricSequenceCreator()
	// { protocols::simple_moves::FavorSymmetricSequenceCreator cr; std::cout << "protocols::simple_moves::FavorSymmetricSequenceCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_FindConsensusSequenceCreator()
	// { protocols::simple_moves::FindConsensusSequenceCreator cr; std::cout << "protocols::simple_moves::FindConsensusSequenceCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_ForceDisulfidesMoverCreator()
	// { protocols::calc_taskop_movers::ForceDisulfidesMoverCreator cr; std::cout << "protocols::calc_taskop_movers::ForceDisulfidesMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_GenericMonteCarloMoverCreator()
	// { protocols::monte_carlo::GenericMonteCarloMoverCreator cr; std::cout << "protocols::monte_carlo::GenericMonteCarloMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_LoadPDBMoverCreator()
	// { protocols::simple_moves::LoadPDBMoverCreator cr; std::cout << "protocols::simple_moves::LoadPDBMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_LoadUnboundRotMoverCreator()
	// { protocols::simple_moves::LoadUnboundRotMoverCreator cr; std::cout << "protocols::simple_moves::LoadUnboundRotMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_MakePolyXMoverCreator()
	// { protocols::pose_creation::MakePolyXMoverCreator cr; std::cout << "protocols::pose_creation::MakePolyXMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_MembraneTopologyCreator()
	// { protocols::simple_moves::MembraneTopologyCreator cr; std::cout << "protocols::simple_moves::MembraneTopologyCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_MergePDBMoverCreator()
	// { protocols::simple_moves::MergePDBMoverCreator cr; std::cout << "protocols::simple_moves::MergePDBMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_MinMoverCreator()
	// { protocols::minimization_packing::MinMoverCreator cr; std::cout << "protocols::minimization_packing::MinMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_MinPackMoverCreator()
	// { protocols::minimization_packing::MinPackMoverCreator cr; std::cout << "protocols::minimization_packing::MinPackMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_ModifyVariantTypeMoverCreator()
	// { protocols::simple_moves::ModifyVariantTypeMoverCreator cr; std::cout << "protocols::simple_moves::ModifyVariantTypeMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_MonteCarloRecoverCreator()
	// { protocols::monte_carlo::MonteCarloRecoverCreator cr; std::cout << "protocols::monte_carlo::MonteCarloRecoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_MonteCarloTestCreator()
	// { protocols::monte_carlo::MonteCarloTestCreator cr; std::cout << "protocols::monte_carlo::MonteCarloTestCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_MSDMoverCreator()
	// { protocols::simple_moves::MSDMoverCreator cr; std::cout << "protocols::simple_moves::MSDMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_MutateResidueCreator()
	// { protocols::simple_moves::MutateResidueCreator cr; std::cout << "protocols::simple_moves::MutateResidueCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_PackRotamersMoverCreator()
	// { protocols::minimization_packing::PackRotamersMoverCreator cr; std::cout << "protocols::minimization_packing::PackRotamersMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_PSSM2BfactorMoverCreator()
	// { protocols::simple_moves::PSSM2BfactorMoverCreator cr; std::cout << "protocols::simple_moves::PSSM2BfactorMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_RepeatPropagationMoverCreator()
	// { protocols::simple_moves::RepeatPropagationMoverCreator cr; std::cout << "protocols::simple_moves::RepeatPropagationMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_ResidueTypeConstraintMoverCreator()
	// { protocols::constraint_movers::ResidueTypeConstraintMoverCreator cr; std::cout << "protocols::simple_moves::ResidueTypeConstraintMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_ReportEffectivePKACreator()
	// { protocols::simple_moves::ReportEffectivePKACreator cr; std::cout << "protocols::simple_moves::ReportEffectivePKACreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_RotamerTrialsMinMoverCreator()
	// { protocols::minimization_packing::RotamerTrialsMinMoverCreator cr; std::cout << "protocols::minimization_packing::RotamerTrialsMinMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_RotamerTrialsMoverCreator()
	// { protocols::minimization_packing::RotamerTrialsMoverCreator cr; std::cout << "protocols::minimization_packing::RotamerTrialsMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_RandomTorsionMoverCreator()
	// { protocols::simple_moves::RandomTorsionMoverCreator cr; std::cout << "protocols::simple_moves::RandomTorsionMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_RandomOmegaFlipMoverCreator()
	// { protocols::simple_moves::RandomOmegaFlipMoverCreator cr; std::cout << "protocols::simple_moves::RandomOmegaFlipMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_SaneMinMoverCreator()
	// { protocols::minimization_packing::SaneMinMoverCreator cr; std::cout << "protocols::minimization_packing::SaneMinMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_ScoreMoverCreator()
	// { protocols::simple_moves::ScoreMoverCreator cr; std::cout << "protocols::simple_moves::ScoreMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_SequenceProfileMoverCreator()
	// { protocols::simple_moves::SequenceProfileMoverCreator cr; std::cout << "protocols::simple_moves::SequenceProfileMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_SetTorsionCreator()
	// { protocols::simple_moves::SetTorsionCreator cr; std::cout << "protocols::simple_moves::SetTorsionCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_ShearMoverCreator()
	// { protocols::simple_moves::ShearMoverCreator cr; std::cout << "protocols::simple_moves::ShearMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_ShortBackrubMoverCreator()
	// { protocols::simple_moves::ShortBackrubMoverCreator cr; std::cout << "protocols::simple_moves::ShortBackrubMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_StorePoseSnapshotCreator()
	// { protocols::simple_moves::StorePoseSnapshotCreator cr; std::cout << "protocols::simple_moves::StorePoseSnapshotCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_SmallMoverCreator()
	// { protocols::simple_moves::SmallMoverCreator cr; std::cout << "protocols::simple_moves::SmallMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_StructProfileMoverCreator()
	// { protocols::simple_moves::StructProfileMoverCreator cr; std::cout << "protocols::simple_moves::StructProfileMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_SuperimposeMoverCreator()
	// { protocols::simple_moves::SuperimposeMoverCreator cr; std::cout << "protocols::simple_moves::SuperimposeMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_PDBReloadMoverCreator()
	// { protocols::simple_moves::PDBReloadMoverCreator cr; std::cout << "protocols::simple_moves::PDBReloadMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_SwitchResidueTypeSetMoverCreator()
	// { protocols::simple_moves::SwitchResidueTypeSetMoverCreator cr; std::cout << "protocols::simple_moves::SwitchResidueTypeSetMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_SwitchChainOrderMoverCreator()
	// { protocols::simple_moves::SwitchChainOrderMoverCreator cr; std::cout << "protocols::simple_moves::SwitchChainOrderMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_TaskAwareMinMoverCreator()
	// { protocols::minimization_packing::TaskAwareMinMoverCreator cr; std::cout << "protocols::minimization_packing::TaskAwareMinMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_VirtualRootMoverCreator()
	// { protocols::simple_moves::VirtualRootMoverCreator cr; std::cout << "protocols::simple_moves::VirtualRootMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_TumbleCreator()
	// { protocols::simple_moves::TumbleCreator cr; std::cout << "protocols::simple_moves::TumbleCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_bin_transitions_InitializeByBinsCreator()
	// { protocols::simple_moves::bin_transitions::InitializeByBinsCreator cr; std::cout << "protocols::simple_moves::bin_transitions::InitializeByBinsCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_bin_transitions_PerturbByBinsCreator()
	// { protocols::simple_moves::bin_transitions::PerturbByBinsCreator cr; std::cout << "protocols::simple_moves::bin_transitions::PerturbByBinsCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_sidechain_moves_SidechainMCMoverCreator()
	// { protocols::simple_moves::sidechain_moves::SidechainMCMoverCreator cr; std::cout << "protocols::simple_moves::sidechain_moves::SidechainMCMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_sidechain_moves_SidechainMoverCreator()
	// { protocols::simple_moves::sidechain_moves::SidechainMoverCreator cr; std::cout << "protocols::simple_moves::sidechain_moves::SidechainMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_sidechain_moves_SetChiMoverCreator()
	// { protocols::simple_moves::sidechain_moves::SetChiMoverCreator cr; std::cout << "protocols::simple_moves::sidechain_moves::SetChiMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_sidechain_moves_JumpRotamerSidechainMoverCreator()
	// { protocols::simple_moves::sidechain_moves::JumpRotamerSidechainMoverCreator cr; std::cout << "protocols::simple_moves::sidechain_moves::JumpRotamerSidechainMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_sidechain_moves_PerturbRotamerSidechainMoverCreator()
	// { protocols::simple_moves::sidechain_moves::PerturbRotamerSidechainMoverCreator cr; std::cout << "protocols::simple_moves::sidechain_moves::PerturbRotamerSidechainMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_sidechain_moves_PerturbChiSidechainMoverCreator()
	// { protocols::simple_moves::sidechain_moves::PerturbChiSidechainMoverCreator cr; std::cout << "protocols::simple_moves::sidechain_moves::PerturbChiSidechainMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_symmetry_ExtractAsymmetricUnitMoverCreator()
	// { protocols::symmetry::ExtractAsymmetricUnitMoverCreator cr; std::cout << "protocols::symmetry::ExtractAsymmetricUnitMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_symmetry_ExtractAsymmetricPoseMoverCreator()
	// { protocols::symmetry::ExtractAsymmetricPoseMoverCreator cr; std::cout << "protocols::symmetry::ExtractAsymmetricPoseMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_symmetry_SetupForSymmetryMoverCreator()
	// { protocols::symmetry::SetupForSymmetryMoverCreator cr; std::cout << "protocols::symmetry::SetupForSymmetryMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_symmetry_DetectSymmetryMoverCreator()
	// { protocols::symmetry::DetectSymmetryMoverCreator cr; std::cout << "protocols::symmetry::DetectSymmetryMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_symmetry_SetupNCSMoverCreator()
	// { protocols::symmetry::SetupNCSMoverCreator cr; std::cout << "protocols::symmetry::SetupNCSMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_symmetry_SymMinMoverCreator()
	// { protocols::minimization_packing::symmetry::SymMinMoverCreator cr; std::cout << "protocols::minimization_packing::symmetry::SymMinMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_symmetry_SymPackRotamersMoverCreator()
	// { protocols::minimization_packing::symmetry::SymPackRotamersMoverCreator cr; std::cout << "protocols::minimization_packing::symmetry::SymPackRotamersMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_symmetry_SymRotamerTrialsMoverCreator()
	// { protocols::minimization_packing::symmetry::SymRotamerTrialsMoverCreator cr; std::cout << "protocols::minimization_packing::symmetry::SymRotamerTrialsMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_simple_moves_symmetry_TaskAwareSymMinMoverCreator()
	// { protocols::minimization_packing::symmetry::TaskAwareSymMinMoverCreator cr; std::cout << "protocols::minimization_packing::symmetry::TaskAwareSymMinMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_toolbox_task_operations_StoreCombinedStoredTasksMoverCreator()
	// { protocols::task_operations::StoreCombinedStoredTasksMoverCreator cr; std::cout << "protocols::task_operations::StoreCombinedStoredTasksMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_toolbox_task_operations_StoreCompoundTaskMoverCreator()
	// { protocols::task_operations::StoreCompoundTaskMoverCreator cr; std::cout << "protocols::task_operations::StoreCompoundTaskMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_toolbox_task_operations_StoreTaskMoverCreator()
	// { protocols::task_operations::StoreTaskMoverCreator cr; std::cout << "protocols::task_operations::StoreTaskMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_environment_EnvMoverCreator()
	// { protocols::environment::EnvMoverCreator cr; std::cout << "protocols::environment::EnvMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_environment_CoMTrackerCMCreator()
	// { protocols::environment::CoMTrackerCMCreator cr; std::cout << "protocols::environment::CoMTrackerCMCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_environment_ScriptCMCreator()
	// { protocols::environment::ScriptCMCreator cr; std::cout << "protocols::environment::ScriptCMCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_abinitio_abscript_AbscriptMoverCreator()
	// { protocols::abinitio::abscript::AbscriptMoverCreator cr; std::cout << "protocols::abinitio::abscript::AbscriptMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_abinitio_abscript_ConstraintPreparerCreator()
	// { protocols::abinitio::abscript::ConstraintPreparerCreator cr; std::cout << "protocols::abinitio::abscript::ConstraintPreparerCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_abinitio_abscript_FragmentCMCreator()
	// { protocols::abinitio::abscript::FragmentCMCreator cr; std::cout << "protocols::abinitio::abscript::FragmentCMCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_abinitio_abscript_FragmentJumpCMCreator()
	// { protocols::abinitio::abscript::FragmentJumpCMCreator cr; std::cout << "protocols::abinitio::abscript::FragmentJumpCMCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_abinitio_abscript_AbscriptLoopCloserCMCreator()
	// { protocols::abinitio::abscript::AbscriptLoopCloserCMCreator cr; std::cout << "protocols::abinitio::abscript::AbscriptLoopCloserCMCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_abinitio_abscript_StructPerturberCMCreator()
	// { protocols::abinitio::abscript::StructPerturberCMCreator cr; std::cout << "protocols::abinitio::abscript::StructPerturberCMCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_abinitio_abscript_RigidChunkCMCreator()
	// { protocols::abinitio::abscript::RigidChunkCMCreator cr; std::cout << "protocols::abinitio::abscript::RigidChunkCMCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_membrane_AddMembraneMoverCreator()
	// { protocols::membrane::AddMembraneMoverCreator cr; std::cout << "protocols::membrane::AddMembraneMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_membrane_FlipMoverCreator()
	// { protocols::membrane::FlipMoverCreator cr; std::cout << "protocols::membrane::FlipMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_membrane_MembranePositionFromTopologyMoverCreator()
	// { protocols::membrane::MembranePositionFromTopologyMoverCreator cr; std::cout << "protocols::membrane::MembranePositionFromTopologyMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_membrane_SetMembranePositionMoverCreator()
	// { protocols::membrane::SetMembranePositionMoverCreator cr; std::cout << "protocols::membrane::SetMembranePositionMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_membrane_TransformIntoMembraneMoverCreator()
	// { protocols::membrane::TransformIntoMembraneMoverCreator cr; std::cout << "protocols::membrane::TransformIntoMembraneMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_membrane_symmetry_SymmetricAddMembraneMoverCreator()
	// { protocols::membrane::symmetry::SymmetricAddMembraneMoverCreator cr; std::cout << "protocols::membrane::symmetry::SymmetricAddMembraneMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_membrane_benchmark_SampleTiltAnglesCreator()
	// { protocols::membrane::benchmark::SampleTiltAnglesCreator cr; std::cout << "protocols::membrane::benchmark::SampleTiltAnglesCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_membrane_benchmark_MakeCanonicalHelixCreator()
	// { protocols::membrane::benchmark::MakeCanonicalHelixCreator cr; std::cout << "protocols::membrane::benchmark::MakeCanonicalHelixCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_membrane_visualize_VisualizeMembraneMoverCreator()
	// { protocols::membrane::visualize::VisualizeMembraneMoverCreator cr; std::cout << "protocols::membrane::visualize::VisualizeMembraneMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_membrane_visualize_VisualizeEmbeddingMoverCreator()
	// { protocols::membrane::visualize::VisualizeEmbeddingMoverCreator cr; std::cout << "protocols::membrane::visualize::VisualizeEmbeddingMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_docking_membrane_MPDockingMoverCreator()
	// { protocols::docking::membrane::MPDockingMoverCreator cr; std::cout << "protocols::docking::membrane::MPDockingMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_docking_membrane_MPDockingSetupMoverCreator()
	// { protocols::docking::membrane::MPDockingSetupMoverCreator cr; std::cout << "protocols::docking::membrane::MPDockingSetupMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_legacy_sewing_AppendAssemblyMoverCreator()
	// { protocols::legacy_sewing::AppendAssemblyMoverCreator cr; std::cout << "protocols::legacy_sewing::AppendAssemblyMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_legacy_sewing_GivenPathAssemblyMoverCreator()
	// { protocols::legacy_sewing::GivenPathAssemblyMoverCreator cr; std::cout << "protocols::legacy_sewing::GivenPathAssemblyMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_legacy_sewing_GreedyAssemblyMoverCreator()
	// { protocols::legacy_sewing::GreedyAssemblyMoverCreator cr; std::cout << "protocols::legacy_sewing::GreedyAssemblyMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_legacy_sewing_MonteCarloAssemblyMoverCreator()
	// { protocols::legacy_sewing::MonteCarloAssemblyMoverCreator cr; std::cout << "protocols::legacy_sewing::MonteCarloAssemblyMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_legacy_sewing_AssemblyConstraintsMoverCreator()
	// { protocols::legacy_sewing::AssemblyConstraintsMoverCreator cr; std::cout << "protocols::legacy_sewing::AssemblyConstraintsMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_legacy_sewing_AddStartnodeFragmentsCreator()
	// { protocols::legacy_sewing::AddStartnodeFragmentsCreator cr; std::cout << "protocols::legacy_sewing::AddStartnodeFragmentsCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_legacy_sewing_RepeatAssemblyMoverCreator()
	// { protocols::legacy_sewing::RepeatAssemblyMoverCreator cr; std::cout << "protocols::legacy_sewing::RepeatAssemblyMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_legacy_sewing_EnumerateAssemblyMoverCreator()
	// { protocols::legacy_sewing::EnumerateAssemblyMoverCreator cr; std::cout << "protocols::legacy_sewing::EnumerateAssemblyMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_symmetric_docking_membrane_MPSymDockMoverCreator()
	// { protocols::symmetric_docking::membrane::MPSymDockMoverCreator cr; std::cout << "protocols::symmetric_docking::membrane::MPSymDockMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_membrane_AddMPLigandMoverCreator()
	// { protocols::membrane::AddMPLigandMoverCreator cr; std::cout << "protocols::membrane::AddMPLigandMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_protocols_relax_membrane_MPFastRelaxMoverCreator()
	// { protocols::relax::membrane::MPFastRelaxMoverCreator cr; std::cout << "protocols::relax::membrane::MPFastRelaxMoverCreator " << cr.keyname() << std::endl; }



};
