// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/init.cc
/// @brief  Declare WidgetRegistrators as static (global) variables in this .cc file
///         so that at load time, they will be initialized, and the Creator classes
///         they register will be handed to the appropriate WidgetFactory.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#include <protocols/init.hh>
#include <core/init.hh>

#include <core/pack/task/operation/TaskOperationRegistrator.hh>
#include <protocols/moves/MoverFactory.hh>

#include <protocols/dna/RestrictDesignToProteinDNAInterfaceCreator.hh>
#include <protocols/dna/WatsonCrickRotamerCouplingsCreator.hh>
#include <protocols/protein_interface_design/movers/FavorNativeResiduePreCycleCreator.hh>
#include <protocols/protein_interface_design/movers/FavorNonNativeResiduePreCycleCreator.hh>
#include <protocols/toolbox/task_operations/SeqprofConsensusOperationCreator.hh>
#include <protocols/toolbox/task_operations/JointSequenceOperationCreator.hh>
#include <protocols/toolbox/task_operations/LimitAromaChi2OperationCreator.hh>
#include <protocols/toolbox/task_operations/LayerDesignOperationCreator.hh>
#include <protocols/toolbox/task_operations/PreventChainFromRepackingOperationCreator.hh>
#include <protocols/toolbox/task_operations/ProteinInterfaceDesignOperationCreator.hh>
#include <protocols/toolbox/task_operations/ReadResfileFromDBCreator.hh>
#include <protocols/toolbox/task_operations/RestrictByCalculatorsOperationCreator.hh>
#include <protocols/toolbox/task_operations/RestrictChainToRepackingOperationCreator.hh>
#include <protocols/toolbox/task_operations/RestrictToInterfaceOperationCreator.hh>
#include <protocols/toolbox/task_operations/RestrictToNeighborhoodOperationCreator.hh>
#include <protocols/toolbox/task_operations/DesignAroundOperationCreator.hh>
#include <protocols/toolbox/task_operations/ThreadSequenceOperationCreator.hh>
#include <protocols/toolbox/task_operations/RestrictNonSurfaceToRepackingOperationCreator.hh>
#include <protocols/toolbox/task_operations/RestrictToInterfaceVectorOperationCreator.hh>
#include <protocols/toolbox/task_operations/RestrictToInterfaceCreator.hh>
#include <protocols/enzdes/DetectProteinLigandInterfaceOperationCreator.hh>
#include <protocols/enzdes/AddLigandMotifRotamersOperationCreator.hh>
#include <protocols/enzdes/ProteinLigandInterfaceUpweighterOperationCreator.hh>
#include <protocols/enzdes/AddRigidBodyLigandConfsCreator.hh>
#include <protocols/enzdes/SetCatalyticResPackBehaviorCreator.hh>
#include <protocols/flexpep_docking/FlexPepDockingProtocolCreator.hh>

#include <protocols/fldsgn/potentials/sspot/NatbiasSecondaryStructureEnergyCreator.hh>

/// Constraint creators
#include <core/scoring/constraints/ConstraintFactory.hh>
#include <protocols/constraints_additional/AdditionalConstraintCreators.hh>
#include <protocols/constraints_additional/BindingSiteConstraint.hh>

/// Mover creators
#include <protocols/moves/InterfaceAnalyzerMoverCreator.hh>
#include <protocols/enzdes/AddOrRemoveMatchCstsCreator.hh>
#include <protocols/enzdes/BackboneSamplerCreator.hh>
#include <protocols/enzdes/EnzRepackMinimizeCreator.hh>
#include <protocols/enzdes/EnzdesMoversCreator.hh>
#include <protocols/moves/ParsedProtocolCreator.hh>
#include <protocols/idealize/IdealizeMoverCreator.hh>
#include <protocols/nonlocal/ExtendedPoseMoverCreator.hh>
#include <protocols/nonlocal/SingleFragmentMoverCreator.hh>
#include <protocols/protein_interface_design/movers/PeptideStapleDesignMoverCreator.hh>
#include <protocols/protein_interface_design/movers/SetTemperatureFactorCreator.hh>
#include <protocols/protein_interface_design/movers/AddSidechainConstraintsToHotspotsCreator.hh>
#include <protocols/protein_interface_design/movers/PlaceSimultaneouslyMoverCreator.hh>
#include <protocols/protein_interface_design/movers/SetupHotspotConstraintsMoverCreator.hh>
#include <protocols/protein_interface_design/movers/HotspotDisjointedFoldTreeMoverCreator.hh>
#include <protocols/protein_interface_design/movers/PlaceStubMoverCreator.hh>
#include <protocols/protein_interface_design/movers/LoopOverCreator.hh>
#include <protocols/protein_interface_design/movers/DesignMinimizeHbondsCreator.hh>
#include <protocols/protein_interface_design/movers/PlaceOnLoopCreator.hh>
#include <protocols/protein_interface_design/movers/RandomMutationCreator.hh>
#include <protocols/protein_interface_design/movers/LoopLengthChangeCreator.hh>
#include <protocols/protein_interface_design/movers/ProteinInterfaceMultiStateDesignMoverCreator.hh>
#include <protocols/protein_interface_design/movers/LoopRemodelCreator.hh>
#include <protocols/protein_interface_design/movers/SaveAndRetrieveSidechainsCreator.hh>
#include <protocols/protein_interface_design/movers/TryRotamersCreator.hh>
#include <protocols/protein_interface_design/movers/SubroutineMoverCreator.hh>
#include <protocols/protein_interface_design/movers/LoopFinderCreator.hh>
#include <protocols/protein_interface_design/movers/HotspotHasherMoverCreator.hh>
#include <protocols/protein_interface_design/movers/DomainAssemblyCreator.hh>
#include <protocols/protein_interface_design/movers/RepackMinimizeCreator.hh>
#include <protocols/protein_interface_design/movers/InterfaceRecapitulationMoverCreator.hh>
#include <protocols/protein_interface_design/movers/DisulfideMoverCreator.hh>
#include <protocols/protein_interface_design/movers/DumpPdbCreator.hh>
#include <protocols/protein_interface_design/movers/SetAtomTreeCreator.hh>
#include <protocols/protein_interface_design/movers/MapHotspotCreator.hh>
#include <protocols/protein_interface_design/movers/BestHotspotCstMoverCreator.hh>
#include <protocols/protein_interface_design/movers/BuildAlaPoseCreator.hh>
#include <protocols/protein_interface_design/movers/DockAndRetrieveSidechainsCreator.hh>
#include <protocols/protein_interface_design/movers/PlacementAuctionMoverCreator.hh>
#include <protocols/protein_interface_design/movers/PlacementMinimizationMoverCreator.hh>
#include <protocols/protein_interface_design/movers/VLBCreator.hh>
#include <protocols/protein_interface_design/movers/BackrubDDMoverCreator.hh>
#include <protocols/protein_interface_design/movers/PrepackMoverCreator.hh>
#include <protocols/protein_interface_design/movers/SpinMoverCreator.hh>
#include <protocols/protein_interface_design/movers/LoopMoverFromCommandLineCreator.hh>
#include <protocols/flxbb/FlxbbDesignCreator.hh>
#include <protocols/flxbb/InterlockAromaCreator.hh>
#include <protocols/fldsgn/BluePrintBDRCreator.hh>
#include <protocols/fldsgn/CircularPermutationCreator.hh>
#include <protocols/fldsgn/potentials/SetSecStructEnergiesCreator.hh>
#include <protocols/fldsgn/potentials/SetAACompositionPotentialCreator.hh>
#include <protocols/forge/remodel/RemodelLoopMoverCreator.hh>
#include <protocols/moves/IfMoverCreator.hh>
#include <protocols/moves/IteratedConvergenceMoverCreator.hh>
#include <protocols/moves/MutateResidueCreator.hh>
#include <protocols/moves/DsspMoverCreator.hh>
#include <protocols/moves/MatcherMoverCreator.hh>
#include <protocols/moves/MakePolyXMoverCreator.hh>
#include <protocols/moves/MinMoverCreator.hh>
#include <protocols/moves/MinPackMoverCreator.hh>
#include <protocols/moves/ScoreMoverCreator.hh>
#include <protocols/moves/ConsensusDesignMoverCreator.hh>
#include <protocols/moves/ConstraintSetMoverCreator.hh>
#include <protocols/moves/PackRotamersMoverCreator.hh>
#include <protocols/moves/RotamerTrialsMoverCreator.hh>
#include <protocols/moves/ConformerSwitchMoverCreator.hh>
#include <protocols/moves/SwitchResidueTypeSetMoverCreator.hh>
#include <protocols/moves/GenericMonteCarloMoverCreator.hh>
#include <protocols/moves/MonteCarloTestCreator.hh>
#include <protocols/moves/MonteCarloRecoverCreator.hh>
#include <protocols/moves/SidechainMCMoverCreator.hh>
#include <protocols/moves/TaskAwareMinMoverCreator.hh>
#include <protocols/moves/ReportToDBCreator.hh>
#include <protocols/moves/RotamerRecoveryMoverCreator.hh>
#include <protocols/moves/MetropolisHastingsMoverCreator.hh>
#include <protocols/moves/BackboneMoverCreator.hh>
#include <protocols/moves/BackrubMoverCreator.hh>
#include <protocols/moves/SidechainMoverCreator.hh>
#include <protocols/moves/BackrubSidechainMoverCreator.hh>
#include <protocols/moves/TrajectoryRecorderCreator.hh>
#include <protocols/moves/MetricRecorderCreator.hh>
#include <protocols/moves/RotamerTrialsMinMoverCreator.hh>
#include <protocols/moves/FavorSequenceProfileCreator.hh>
#include <protocols/relax/FastRelaxCreator.hh>
#include <protocols/dna/SeparateDnaFromNonDnaCreator.hh>
#include <protocols/dna/DnaInterfaceMinMoverCreator.hh>
#include <protocols/dna/DnaInterfacePackerCreator.hh>
#include <protocols/dna/DnaInterfaceMultiStateDesignCreator.hh>
#include <protocols/dna/DesignProteinBackboneAroundDNACreator.hh>
#include <protocols/motifs/MotifDnaPackerCreator.hh>
#include <protocols/ligand_docking/LigandDesignCreator.hh>
#include <protocols/ligand_docking/GrowLigandCreator.hh>
#include <protocols/ligand_docking/AddHydrogensCreator.hh>
#include <protocols/ligand_docking/TranslateCreator.hh>
#include <protocols/ligand_docking/TransformCreator.hh>
#include <protocols/ligand_docking/CompoundTranslateCreator.hh>
#include <protocols/ligand_docking/RotateCreator.hh>
#include <protocols/ligand_docking/StartFromCreator.hh>
#include <protocols/ligand_docking/RandomConformersCreator.hh>
#include <protocols/ligand_docking/SlideTogetherCreator.hh>
#include <protocols/ligand_docking/MinimizeBackboneCreator.hh>
#include <protocols/ligand_docking/HighResDockerCreator.hh>
#include <protocols/ligand_docking/FinalMinimizerCreator.hh>
#include <protocols/ligand_docking/InterfaceScoreCalculatorCreator.hh>
#include <protocols/loops/LoopMover_CCDCreator.hh>
#include <protocols/rosetta_scripts/SavePoseMoverCreator.hh>
#include <protocols/loophash/LoopHashMoverWrapperCreator.hh>
#include <protocols/docking/DockingProtocolCreator.hh>
#include <protocols/moves/RollMoverCreator.hh>

// Filter creators
#include <protocols/filters/FilterFactory.hh>

#include <protocols/enzdes/EnzFilterCreators.hh>
#include <protocols/protein_interface_design/dock_design_filter_creators.hh>
#include <protocols/protein_interface_design/filters/InterfaceHolesFilterCreator.hh>
#include <protocols/protein_interface_design/filters/DesignableResiduesFilterCreator.hh>
#include <protocols/protein_interface_design/filters/RmsdFilterCreator.hh>
#include <protocols/protein_interface_design/filters/DeltaFilterCreator.hh>
#include <protocols/protein_interface_design/filters/RelativePoseFilterCreator.hh>
#include <protocols/protein_interface_design/filters/BindingStrainFilterCreator.hh>
#include <protocols/protein_interface_design/filters/AverageDegreeFilterCreator.hh>
#include <protocols/protein_interface_design/filters/BoltzmannFilterCreator.hh>
#include <protocols/protein_interface_design/filters/FilterScanCreator.hh>
#include <protocols/protein_interface_design/filters/RotamerBoltzmannWeightFilterCreator.hh>
#include <protocols/protein_interface_design/filters/DisulfideFilterCreator.hh>
#include <protocols/protein_interface_design/filters/StubScoreFilterCreator.hh>
#include <protocols/protein_interface_design/filters/SequenceRecoveryFilterCreator.hh>
#include <protocols/ligand_docking/CompleteConnectionsFilterCreator.hh>
#include <protocols/ligand_docking/HeavyAtomFilterCreator.hh>
#include <protocols/filters/HolesFilterCreator.hh>
#include <protocols/filters/PackStatFilterCreator.hh>
#include <protocols/filters/BasicFilterCreators.hh>
#include <protocols/filters/ConservedPosMutationFilterCreator.hh>
#include <protocols/filters/ContingentFilterCreator.hh>
#include <protocols/filters/RGFilterCreator.hh>
#include <protocols/filters/ScoreCutoffFilterCreator.hh>
#include <protocols/filters/AtomicContactFilterCreator.hh>
#include <protocols/filters/TimeFilterCreator.hh>
#include <protocols/filters/AtomicDistanceFilterCreator.hh>
#include <protocols/filters/ScFilterCreator.hh>
#include <protocols/fldsgn/filters/CoreDunbrackFilterCreator.hh>
#include <protocols/fldsgn/filters/FragQualFilterCreator.hh>
#include <protocols/fldsgn/filters/HelixKinkFilterCreator.hh>
#include <protocols/fldsgn/filters/HelixPairingFilterCreator.hh>
#include <protocols/fldsgn/filters/HSSTripletFilterCreator.hh>
#include <protocols/fldsgn/filters/InterlockingAromaFilterCreator.hh>
#include <protocols/fldsgn/filters/NcontactsFilterCreator.hh>
#include <protocols/fldsgn/filters/ParallelBetaPairingPreferenceFilterCreator.hh>
#include <protocols/fldsgn/filters/SecondaryStructureFilterCreator.hh>
#include <protocols/fldsgn/filters/SheetTopologyFilterCreator.hh>

//Scoring Grid creators
#include <protocols/qsar/scoring_grid/GridFactory.hh>
#include <protocols/qsar/scoring_grid/AtrGridCreator.hh>
#include <protocols/qsar/scoring_grid/RepGridCreator.hh>
#include <protocols/qsar/scoring_grid/HbdGridCreator.hh>
#include <protocols/qsar/scoring_grid/HbaGridCreator.hh>
#include <protocols/qsar/scoring_grid/polarizGridCreator.hh>
#include <protocols/qsar/scoring_grid/VdwGridCreator.hh>

// DataLoader creators
#include <protocols/jd2/parser/DataLoaderFactory.hh>
#include <protocols/jd2/parser/StandardLoaderCreators.hh>
#include <protocols/ligand_docking/LigandDockingLoaderCreators.hh>

#include <core/scoring/methods/EnergyMethodRegisterer.hh>
#include <protocols/scoring/methods/ResidualDipolarCouplingEnergyRigidSegmentsCreator.hh>
#include <protocols/scoring/methods/SpecialRotamerEnergyCreator.hh>
#include <protocols/scoring/methods/pcs/PseudocontactShiftEnergyCreator.hh>
#include <protocols/scoring/methods/pcs2/PcsEnergyCreator.hh>
#include <protocols/scoring/methods/EnsembleEnergyCreator.hh>
#include <protocols/scoring/methods/InterchainEnvEnergyCreator.hh>
#include <protocols/scoring/methods/InterchainPairEnergyCreator.hh>
#include <protocols/scoring/methods/saxs/PDDFEnergyCreator.hh>

namespace protocols {

static core::scoring::methods::EnergyMethodRegistrator< scoring::methods::ResidualDipolarCouplingEnergyRigidSegmentsCreator > ResidualDipolarCouplingEnergyRigidSegmentsCreator_registrator;
static core::scoring::methods::EnergyMethodRegistrator< scoring::methods::SpecialRotamerEnergyCreator > SpecialRotamerEnergyCreator_registrator;
static core::scoring::methods::EnergyMethodRegistrator< scoring::methods::pcs::PseudocontactShiftEnergyCreator > PseudoconstactShiftEnergyCreator_registrator;
static core::scoring::methods::EnergyMethodRegistrator< scoring::methods::pcs2::PcsEnergyCreator > PcsEnergyCreator_registrator;
static core::scoring::methods::EnergyMethodRegistrator< scoring::methods::EnsembleEnergyCreator > EnsembleEnergyCreator_registrator;
static core::scoring::methods::EnergyMethodRegistrator< scoring::methods::InterchainEnvEnergyCreator > InterchainEnvEnergyCreator_registrator;
static core::scoring::methods::EnergyMethodRegistrator< scoring::methods::InterchainPairEnergyCreator > InterchainPairEnergyCreator_registrator;
static core::scoring::methods::EnergyMethodRegistrator< protocols::scoring::methods::saxs::PDDFEnergyCreator > PDDFEnergyCreator_registrator;
static core::scoring::methods::EnergyMethodRegistrator< protocols::fldsgn::potentials::sspot::NatbiasSecondaryStructureEnergyCreator > NatbiasSecondaryStructureEnergyCreator_registrator;

/// Constraint Registrators
static core::scoring::constraints::ConstraintRegistrator< protocols::constraints_additional::SequenceProfileConstraintCreator > SequenceProfileConstraintCreator_registrator;
static core::scoring::constraints::ConstraintRegistrator< protocols::constraints_additional::BindingSiteConstraintCreator > BindingSiteConstraintCreator_registrator;
static core::scoring::constraints::ConstraintRegistrator< protocols::constraints_additional::PocketConstraintCreator > PocketConstraintCreator_registrator;

using namespace core::pack::task::operation;

static TaskOperationRegistrator< protocols::dna::RestrictDesignToProteinDNAInterfaceCreator > RestrictDesignToProteinDNAInterfaceCreator_registrator;
static TaskOperationRegistrator< protocols::dna::WatsonCrickRotamerCouplingsCreator > WatsonCrickRotamerCouplingsCreator_registrator;
static TaskOperationRegistrator< protocols::enzdes::AddLigandMotifRotamersOperationCreator > AddLigandMotifRotamersOperationCreator_registrator;
static TaskOperationRegistrator< protocols::toolbox::task_operations::SeqprofConsensusOperationCreator > SeqprofConsensusOperationCreator_registrator;
static TaskOperationRegistrator< protocols::toolbox::task_operations::JointSequenceOperationCreator > JointSequenceOperationCreator_registrator;
static TaskOperationRegistrator< protocols::toolbox::task_operations::RestrictConservedLowDdgOperationCreator > RestrictConservedLowDdgOperationCreator_registrator;
static TaskOperationRegistrator< protocols::toolbox::task_operations::LimitAromaChi2OperationCreator > LimitAromaChi2OperationCreator_registrator;
static TaskOperationRegistrator< protocols::toolbox::task_operations::LayerDesignOperationCreator > LayerDesignOperationCreator_registrator;
static TaskOperationRegistrator< protocols::toolbox::task_operations::PreventChainFromRepackingOperationCreator > PreventChainFromRepackingOperationCreator_registrator;
static TaskOperationRegistrator< protocols::toolbox::task_operations::ProteinInterfaceDesignOperationCreator > ProteinInterfaceDesignOperationCreator_registrator;
static TaskOperationRegistrator< protocols::toolbox::task_operations::ReadResfileFromDBCreator > ReadResfileFromDBCreator_registrator;
static TaskOperationRegistrator< protocols::toolbox::task_operations::RestrictByCalculatorsOperationCreator > RestrictByCalculatorsOperationCreator_registrator;
static TaskOperationRegistrator< protocols::toolbox::task_operations::RestrictChainToRepackingOperationCreator > RestrictChainToRepackingOperationCreator_registrator;
static TaskOperationRegistrator< protocols::toolbox::task_operations::RestrictToInterfaceOperationCreator > RestrictToInterfaceOperationCreator_registrator;
static TaskOperationRegistrator< protocols::toolbox::task_operations::RestrictToInterfaceVectorOperationCreator > RestrictToInterfaceVectorOperationCreator_registrator;
static TaskOperationRegistrator< protocols::toolbox::task_operations::RestrictToNeighborhoodOperationCreator > RestrictToNeighborhoodOperationCreator_registrator;
static TaskOperationRegistrator< protocols::toolbox::task_operations::DesignAroundOperationCreator > DesignAroundOperationCreator_registrator;
static TaskOperationRegistrator< protocols::toolbox::task_operations::ThreadSequenceOperationCreator > ThreadSequenceOperationCreator_registrator;
static TaskOperationRegistrator< protocols::toolbox::task_operations::RestrictNonSurfaceToRepackingOperationCreator > RestrictNonSurfaceToRepackingOperationCreator_registrator;

static TaskOperationRegistrator< protocols::enzdes::DetectProteinLigandInterfaceOperationCreator > DetectProteinLigandInterfaceOperationCreator_registrator;
static TaskOperationRegistrator< protocols::enzdes::ProteinLigandInterfaceUpweighterOperationCreator > ProteinLigandInterfaceUpweighterOperationCreator_registrator;
static TaskOperationRegistrator< protocols::enzdes::AddRigidBodyLigandConfsCreator > AddRigidBodyLigandConfsCreator_registrator;
static TaskOperationRegistrator< protocols::enzdes::SetCatalyticResPackBehaviorCreator > SetCatalyticResPackBehaviorCreator_registrator;
static TaskOperationRegistrator< protocols::toolbox::task_operations::RestrictToInterfaceCreator > RestrictToInterfaceCreator_registrator;

using namespace moves;

static MoverRegistrator< protocols::moves::InterfaceAnalyzerMoverCreator > reg_InterfaceAnalyzerMoverCreator;
static MoverRegistrator< enzdes::AddOrRemoveMatchCstsCreator > reg_AddOrRemoveMatchCstsCreator;
static MoverRegistrator< enzdes::BackboneSamplerCreator > reg_BackboneSamplerCreator;
static MoverRegistrator< enzdes::EnzRepackMinimizeCreator > reg_EnzRepackMinimizeCreator;
static MoverRegistrator< enzdes::PredesignPerturbMoverCreator > reg_PredesignPerturbMoverCreator;
static MoverRegistrator< protocols::moves::ParsedProtocolCreator > reg_ParsedProtocolCreator;
static MoverRegistrator< protocols::idealize::IdealizeMoverCreator > reg_IdealizeMoverCreator;
static MoverRegistrator< protein_interface_design::movers::PeptideStapleDesignMoverCreator > reg_PeptideStapleDesignMoverCreator;
static MoverRegistrator< protein_interface_design::movers::AddSidechainConstraintsToHotspotsCreator > reg_AddSidechainConstraintsToHotspotsCreator;
static MoverRegistrator< protein_interface_design::movers::PlaceSimultaneouslyMoverCreator > reg_PlaceSimultaneouslyMoverCreator;
static MoverRegistrator< protein_interface_design::movers::SetTemperatureFactorCreator > reg_SetTemperatureFactorCreator;
static MoverRegistrator< protein_interface_design::movers::SetupHotspotConstraintsMoverCreator > reg_SetupHotspotConstraintsMoverCreator;
static MoverRegistrator< protein_interface_design::movers::PlaceStubMoverCreator > reg_PlaceStubMoverCreator;
static MoverRegistrator< protein_interface_design::movers::LoopOverCreator > reg_LoopOverCreator;
static MoverRegistrator< protein_interface_design::movers::SubroutineMoverCreator > reg_SubroutineCreator;
static MoverRegistrator< protein_interface_design::movers::DesignMinimizeHbondsCreator > reg_DesignMinimizeHbondsCreator;
static MoverRegistrator< protein_interface_design::movers::PlaceOnLoopCreator > reg_PlaceOnLoopCreator;
static MoverRegistrator< protein_interface_design::movers::ProteinInterfaceMultiStateDesignMoverCreator > reg_ProteinInterfaceMultiStateDesignMoverCreator;
static MoverRegistrator< protein_interface_design::movers::LoopRemodelCreator > reg_LoopRemodelCreator;
static MoverRegistrator< protein_interface_design::movers::SaveAndRetrieveSidechainsCreator > reg_SaveAndRetrieveSidechainsCreator;
static MoverRegistrator< protein_interface_design::movers::TryRotamersCreator > reg_TryRotamersCreator;
static MoverRegistrator< protein_interface_design::movers::LoopFinderCreator > reg_LoopFinderCreator;
static MoverRegistrator< protein_interface_design::movers::HotspotHasherMoverCreator > reg_HotspotHasherMoverCreator;
static MoverRegistrator< protein_interface_design::movers::DomainAssemblyCreator > reg_DomainAssemblyCreator;
static MoverRegistrator< protein_interface_design::movers::RepackMinimizeCreator > reg_RepackMinimizeCreator;
static MoverRegistrator< protein_interface_design::movers::InterfaceRecapitulationMoverCreator > reg_InterfaceRecapitulationMoverCreator;
static MoverRegistrator< protein_interface_design::movers::DisulfideMoverCreator > reg_DisulfideMoverCreator;
static MoverRegistrator< protein_interface_design::movers::DumpPdbCreator > reg_DumpPdbCreator;
static MoverRegistrator< protein_interface_design::movers::SetAtomTreeCreator > reg_SetAtomTreeCreator;
static MoverRegistrator< protein_interface_design::movers::MapHotspotCreator > reg_MapHotspotCreator;
static MoverRegistrator< protein_interface_design::movers::BestHotspotCstMoverCreator > reg_BestHotspotCstMoverCreator;
static MoverRegistrator< protein_interface_design::movers::BuildAlaPoseCreator > reg_BuildAlaPoseCreator;
static MoverRegistrator< protein_interface_design::movers::FavorNativeResiduePreCycleCreator > reg_FavorNativeResiduePreCycleCreator;
static MoverRegistrator< protein_interface_design::movers::FavorNonNativeResiduePreCycleCreator > reg_FavorNonNativeResiduePreCycleCreator;
static MoverRegistrator< protein_interface_design::movers::DockAndRetrieveSidechainsCreator > reg_DockAndRetrieveSidechainsCreator;
static MoverRegistrator< protein_interface_design::movers::PlacementAuctionMoverCreator > reg_PlacementAuctionMoverCreator;
static MoverRegistrator< protein_interface_design::movers::PlacementMinimizationMoverCreator > reg_PlacementMinimizationMoverCreator;
static MoverRegistrator< protein_interface_design::movers::RandomMutationCreator > reg_RandomMutationCreator;
static MoverRegistrator< protein_interface_design::movers::LoopLengthChangeCreator > reg_LoopLengthChangeCreator;
static MoverRegistrator< protein_interface_design::movers::VLBCreator > reg_VLBCreator;
static MoverRegistrator< protein_interface_design::movers::BackrubDDMoverCreator > reg_BackrubDDMoverCreator;
static MoverRegistrator< protein_interface_design::movers::PrepackMoverCreator > reg_PrepackMoverCreator;
static MoverRegistrator< protein_interface_design::movers::SpinMoverCreator > reg_SpinMoverCreator;
static MoverRegistrator< protein_interface_design::movers::HotspotDisjointedFoldTreeMoverCreator > reg_HotsotDisjointedFoldTreeMoverCreator;
static MoverRegistrator< protein_interface_design::movers::LoopMoverFromCommandLineCreator > reg_LoopMoverFromCommandLine;
static MoverRegistrator< flxbb::FlxbbDesignCreator > reg_FlxbbDesignCreator;
static MoverRegistrator< flxbb::InterlockAromaCreator > reg_InterlockAromaCreator;
static MoverRegistrator< fldsgn::BluePrintBDRCreator > reg_BluePrintBDRCreator;
static MoverRegistrator< fldsgn::CircularPermutationCreator > reg_CircularPermutationCreator;
static MoverRegistrator< fldsgn::potentials::SetSecStructEnergiesCreator > reg_SetSecStructEnergiesCreator;
static MoverRegistrator< fldsgn::potentials::SetAACompositionPotentialCreator > reg_AACompositionPotentialCreator;
static MoverRegistrator< forge::remodel::RemodelLoopMoverCreator > reg_RemodelLoopMoverCreator;
static MoverRegistrator< moves::IfMoverCreator > reg_IfMoverCreator;
static MoverRegistrator< moves::IteratedConvergenceMoverCreator > reg_IteratedConvergenceMoverCreator;
static MoverRegistrator< moves::MutateResidueCreator > reg_MutateResidueCreator;
static MoverRegistrator< moves::DsspMoverCreator > reg_DsspMoverCreator;
static MoverRegistrator< moves::MakePolyXMoverCreator > reg_MakePolyXMoverCreator;
static MoverRegistrator< moves::MatcherMoverCreator > reg_MatcherMoverCreator;
static MoverRegistrator< moves::MinMoverCreator > reg_MinMoverCreator;
static MoverRegistrator< moves::MinPackMoverCreator > reg_MinPackMoverCreator;
static MoverRegistrator< moves::ScoreMoverCreator > reg_ScoreMoverCreator;
static MoverRegistrator< moves::ConsensusDesignMoverCreator > reg_ConsensusDesignMoverCreator;
static MoverRegistrator< moves::ConstraintSetMoverCreator > reg_ConstraintSetMoverCreator;
static MoverRegistrator< moves::PackRotamersMoverCreator > reg_PackRotamersMoverCreator;
static MoverRegistrator< moves::ConformerSwitchMoverCreator > reg_ConformerSwitchMoverCreator;
static MoverRegistrator< moves::SwitchResidueTypeSetMoverCreator > reg_SwitchResidueTypeSetMoverCreator;
static MoverRegistrator< moves::GenericMonteCarloMoverCreator > reg_GenericMonteCarloMoverCreator;
static MoverRegistrator< moves::MonteCarloTestCreator > reg_MonteCarloTestCreator;
static MoverRegistrator< moves::MonteCarloRecoverCreator > reg_GenericMonteCarloRecoverCreator;
static MoverRegistrator< moves::SidechainMCMoverCreator > reg_SidechainMCMoverCreator;
static MoverRegistrator< moves::TaskAwareMinMoverCreator > reg_TaskAwareMinMoverCreator;
static MoverRegistrator< moves::ReportToDBCreator > reg_ReportToDBCreator;
static MoverRegistrator< moves::RotamerRecoveryMoverCreator > reg_RotamerRecoveryMoverCreator;
static MoverRegistrator< moves::MetropolisHastingsMoverCreator > reg_MetropolisHastingsMoverCreator;
static MoverRegistrator< moves::SmallMoverCreator > reg_SmallMoverCreator;
static MoverRegistrator< moves::ShearMoverCreator > reg_ShearMoverCreator;
static MoverRegistrator< moves::BackrubMoverCreator > reg_BackrubMoverCreator;
static MoverRegistrator< moves::SidechainMoverCreator > reg_SidechainMoverCreator;
static MoverRegistrator< moves::BackrubSidechainMoverCreator > reg_BackrubSidechainMoverCreator;
static MoverRegistrator< moves::TrajectoryRecorderCreator > reg_TrajectoryRecorderCreator;
static MoverRegistrator< moves::MetricRecorderCreator > reg_MetricRecorderCreator;
static MoverRegistrator< moves::RotamerTrialsMoverCreator > reg_RotamerTrialsMoverCreator;
static MoverRegistrator< moves::RotamerTrialsMinMoverCreator > reg_RotamerTrialsMinMoverCreator;
static MoverRegistrator< moves::FavorSequenceProfileCreator > reg_FavorSequenceProfileCreator;
static MoverRegistrator< flexpep_docking::FlexPepDockingProtocolCreator > reg_FlexPepDockingProtocolCreator;
static MoverRegistrator< nonlocal::ExtendedPoseMoverCreator > reg_ExtendedPoseMoverCreator;
static MoverRegistrator< nonlocal::SingleFragmentMoverCreator > reg_SingleFragmentMoverCreator;
static MoverRegistrator< relax::FastRelaxCreator> reg_FastRelaxCreator;
static MoverRegistrator< dna::SeparateDnaFromNonDnaCreator > reg_SeparateDnaFromNonDnaCreator;
static MoverRegistrator< dna::DnaInterfaceMinMoverCreator > reg_DnaInterfaceMinMoverCreator;
static MoverRegistrator< dna::DnaInterfacePackerCreator > reg_DnaInterfacePackerCreator;
static MoverRegistrator< dna::DnaInterfaceMultiStateDesignCreator > reg_DnaInterfaceMultiStateDesignCreator;
static MoverRegistrator< dna::DesignProteinBackboneAroundDNACreator > reg_DesignProteinBackboneAroundDNACreator;
static MoverRegistrator< motifs::MotifDnaPackerCreator > reg_MotifDnaPackerCreator;
static MoverRegistrator< ligand_docking::LigandDesignCreator > reg_LigandDesignCreator;
static MoverRegistrator< ligand_docking::GrowLigandCreator > reg_GrowLigandCreator;
static MoverRegistrator< ligand_docking::AddHydrogensCreator > reg_AddHydrogensCreator;
static MoverRegistrator< ligand_docking::TranslateCreator > reg_TranslateCreator;
static MoverRegistrator< ligand_docking::TransformCreator > reg_TransformCreator;
static MoverRegistrator< ligand_docking::CompoundTranslateCreator > reg_CompoundTranslateCreator;
static MoverRegistrator< ligand_docking::RotateCreator > reg_RotateCreator;
static MoverRegistrator< ligand_docking::StartFromCreator > reg_StartFromCreator;
static MoverRegistrator< ligand_docking::RandomConformersCreator > reg_RandomConformersCreator;
static MoverRegistrator< ligand_docking::SlideTogetherCreator > reg_SlideTogetherCreator;
static MoverRegistrator< ligand_docking::MinimizeBackboneCreator > reg_MinimizeBackboneCreator;
static MoverRegistrator< ligand_docking::HighResDockerCreator > reg_HighResDockerCreator;
static MoverRegistrator< ligand_docking::FinalMinimizerCreator > reg_FinalMinimizerCreator;
static MoverRegistrator< ligand_docking::InterfaceScoreCalculatorCreator > reg_InterfaceScoreCalculatorCreator;
static MoverRegistrator< loops::LoopMover_Refine_CCDCreator > reg_LoopMover_Refine_CCDCreator;
static MoverRegistrator< rosetta_scripts::SavePoseMoverCreator > SavePoseMoverCreator;
static MoverRegistrator< loophash::LoopHashMoverWrapperCreator > reg_LoopHashMoverWrapperCreator;
static MoverRegistrator< docking::DockingProtocolCreator > DockingProtocolCreator;
static MoverRegistrator< moves::RollMoverCreator > reg_RollMoverCreator;

using namespace filters;

static FilterRegistrator< protocols::enzdes::DiffAtomSasaFilterCreator > reg_DiffAtomSasaFilterCreator;
static FilterRegistrator< protocols::enzdes::EnzScoreFilterCreator > reg_EnzScoreFilterCreator;
static FilterRegistrator< protocols::enzdes::LigBurialFilterCreator > reg_LigBurialFilterCreator;
static FilterRegistrator< protocols::enzdes::LigDSasaFilterCreator > reg_LigDSasaFilterCreator;
static FilterRegistrator< protocols::enzdes::LigInterfaceEnergyFilterCreator > reg_LigInterfaceEnergyFilterCreator;
static FilterRegistrator< protocols::enzdes::RepackWithoutLigandFilterCreator > reg_RepackWithoutLigandFilterCreator;
static FilterRegistrator< protocols::enzdes::EnzdesScorefileFilterCreator > reg_EnzdesScorefileFilterCreator;
static FilterRegistrator< protocols::filters::AtomicContactFilterCreator > reg_AtomicContactFilterCreator;
static FilterRegistrator< protocols::filters::TimeFilterCreator > reg_TimeFilterCreator;
static FilterRegistrator< protocols::filters::AtomicDistanceFilterCreator > reg_AtomicDistanceFilterCreator;
static FilterRegistrator< protocols::ligand_docking::CompleteConnectionsFilterCreator > reg_CompleteConnectionsFilterCreator;
static FilterRegistrator< protocols::ligand_docking::HeavyAtomFilterCreator > reg_HeavyAtomFilterCreator;
static FilterRegistrator< protocols::filters::CombinedFilterCreator > reg_CombinedFilterCreator;
static FilterRegistrator< protocols::filters::CompoundFilterCreator > reg_CompoundFilterCreator;
static FilterRegistrator< protocols::filters::ConservedPosMutationFilterCreator > reg_ConservedPosMutationFilterCreator;
static FilterRegistrator< protocols::filters::ContingentFilterCreator > reg_ContingentFilterCreator;
static FilterRegistrator< protocols::filters::FalseFilterCreator > reg_FalseFilterCreator;
static FilterRegistrator< protocols::filters::HolesFilterCreator > reg_HolesFilterCreator;
static FilterRegistrator< protocols::filters::MoveBeforeFilterCreator > reg_MoveBeforeFilterCreator;
static FilterRegistrator< protocols::filters::PackStatFilterCreator > reg_PackStatFilterCreator;
static FilterRegistrator< protocols::filters::StochasticFilterCreator > reg_StochasticFilterCreator;
static FilterRegistrator< protocols::filters::ScoreCutoffFilterCreator > reg_ScoreCutoffFilterCreator;
static FilterRegistrator< protocols::filters::ScFilterCreator > reg_ScFilterCreator;
static FilterRegistrator< protocols::fldsgn::filters::CoreDunbrackFilterCreator > reg_CoreDunbrackFilterCreator;
static FilterRegistrator< protocols::fldsgn::filters::FragQualFilterCreator > reg_FragQualFilterCreator;
static FilterRegistrator< protocols::fldsgn::filters::HelixPairingFilterCreator > reg_HelixPairingFilterCreator;
static FilterRegistrator< protocols::fldsgn::filters::HelixKinkFilterCreator > reg_HelixKinkFilterCreator;
static FilterRegistrator< protocols::fldsgn::filters::HSSTripletFilterCreator > reg_HSSTripletFilterCreator;
static FilterRegistrator< protocols::fldsgn::filters::InterlockingAromaFilterCreator > reg_InterlockingAromaFilterCreator;
static FilterRegistrator< protocols::fldsgn::filters::NcontactsFilterCreator > reg_NcontactsFilterCreator;
static FilterRegistrator< protocols::fldsgn::filters::SecondaryStructureFilterCreator > reg_SecondaryStructureFilterCreator;
static FilterRegistrator< protocols::fldsgn::filters::ParallelBetaPairingPreferenceFilterCreator > reg_ParallelBetaPairingPreferenceFilterCreator;
static FilterRegistrator< protocols::fldsgn::filters::SheetTopologyFilterCreator > reg_SheetTopologyFilterCreator;
static FilterRegistrator< protocols::protein_interface_design::AlaScanFilterCreator > reg_AlaScanFilterCreator;
static FilterRegistrator< protocols::protein_interface_design::BuriedUnsatHbondFilterCreator > reg_BuriedUnsatHbondFilterCreator;
static FilterRegistrator< protocols::protein_interface_design::DdgFilterCreator > reg_DdgFilterCreator;
static FilterRegistrator< protocols::protein_interface_design::EnergyPerResidueFilterCreator > reg_EnergyPerResidueFilterCreator;
static FilterRegistrator< protocols::protein_interface_design::filters::DesignableResiduesFilterCreator > reg_DesignableResiduesFilterCreator;
static FilterRegistrator< protocols::protein_interface_design::filters::DisulfideFilterCreator > reg_DisulfideFilterCreator;
static FilterRegistrator< protocols::protein_interface_design::filters::InterfaceHolesFilterCreator > reg_InterfaceHolesFilterCreator;
static FilterRegistrator< protocols::protein_interface_design::filters::RmsdFilterCreator > reg_RmsdFilterCreator;
static FilterRegistrator< protocols::protein_interface_design::filters::AverageDegreeFilterCreator > reg_AverageFilterCreator;
static FilterRegistrator< protocols::protein_interface_design::filters::BoltzmannFilterCreator > reg_BoltzmannCreator;
static FilterRegistrator< protocols::protein_interface_design::filters::DeltaFilterCreator > reg_DeltaFilterCreator;
static FilterRegistrator< protocols::protein_interface_design::filters::RelativePoseFilterCreator > reg_RelativePoseFilterCreator;
static FilterRegistrator< protocols::protein_interface_design::filters::BindingStrainFilterCreator > reg_BindingStrainFilterCreator;
static FilterRegistrator< protocols::protein_interface_design::filters::FilterScanFilterCreator > reg_FilterScanFilterCreator;
static FilterRegistrator< protocols::protein_interface_design::filters::RotamerBoltzmannWeightFilterCreator > reg_RotamerBoltzmannWeightFilterCreator;
static FilterRegistrator< protocols::protein_interface_design::filters::SequenceRecoveryFilterCreator > reg_SequenceRecoveryFilterCreator;
static FilterRegistrator< protocols::protein_interface_design::filters::StubScoreFilterCreator > reg_StubScoreFilterCreator;
static FilterRegistrator< protocols::protein_interface_design::HbondsToResidueFilterCreator > reg_HbondsToResidueFilterCreator;
static FilterRegistrator< protocols::protein_interface_design::InterfaceSasaFilterCreator > reg_InterfaceSasaFilterCreator;
static FilterRegistrator< protocols::protein_interface_design::NeighborTypeFilterCreator > reg_NeighborTypeFilterCreator;
static FilterRegistrator< protocols::protein_interface_design::ResidueBurialFilterCreator > reg_ResidueBurialFilterCreator;
static FilterRegistrator< protocols::protein_interface_design::ResidueDistanceFilterCreator > reg_ResidueDistanceFilterCreator;
static FilterRegistrator< protocols::protein_interface_design::ResiduesInInterfaceFilterCreator > reg_ResiduesInInterfaceFilterCreator;
static FilterRegistrator< protocols::protein_interface_design::ScoreTypeFilterCreator > reg_ScoreTypeFilterCreator;
static FilterRegistrator< protocols::protein_interface_design::TerminusDistanceFilterCreator > reg_TerminusDistanceFilterCreator;

using namespace qsar::scoring_grid;
static GridRegistrator<AtrGridCreator> reg_AtrGridCreator;
static GridRegistrator<RepGridCreator> reg_RepGridCreator;
static GridRegistrator<HbaGridCreator> reg_HbaGridCreator;
static GridRegistrator<HbdGridCreator> reg_HbdGridCreator;
static GridRegistrator<polarizGridCreator> reg_polarizGridCreator;
static GridRegistrator<VdwGridCreator> reg_VdwGridCreator;

using namespace jd2::parser;
static DataLoaderRegistrator< ScoreFunctionLoaderCreator > reg_ScoreFunctionLoaderCreator;
static DataLoaderRegistrator< TaskOperationLoaderCreator > reg_TaskOperationLoaderCreator;
static DataLoaderRegistrator< ScoringGridLoaderCreator > reg_ScoringGridLoaderCreator;
static DataLoaderRegistrator< FragSetLoaderCreator > reg_FragSetLoaderCreator;
static DataLoaderRegistrator< MonteCarloLoaderCreator > reg_MonteCarloLoaderCreator;
static DataLoaderRegistrator< ligand_docking::InterfaceBuilderLoaderCreator > reg_InterfaceBuilderLoaderCreator;
static DataLoaderRegistrator< ligand_docking::MoveMapBuilderLoaderCreator > reg_MoveMapBuilderLoaderCreator;
static DataLoaderRegistrator< ligand_docking::LigandAreaLoaderCreator > reg_LigandAreaLoaderCreator;

void init( int argc, char * argv [] )
{
	core::init( argc, argv );
}

void init( utility::vector1< std::string > const & args )
{
	core::init( args );
}

}

