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

// Filter creator headers
#include <protocols/analysis/PeptideDeriverFilterCreator.hh>
#include <protocols/cyclic_peptide/OversaturatedHbondAcceptorFilterCreator.hh>
#include <protocols/enzdes/EnzFilterCreators.hh>
#include <protocols/enzdes/EnzFilterCreators.hh>
#include <protocols/enzdes/EnzFilterCreators.hh>
#include <protocols/enzdes/EnzFilterCreators.hh>
#include <protocols/enzdes/EnzFilterCreators.hh>
#include <protocols/enzdes/EnzFilterCreators.hh>
#include <protocols/enzdes/EnzFilterCreators.hh>
#include <protocols/enzdes/EnzFilterCreators.hh>
#include <protocols/enzdes/RemoveLigandFilterCreator.hh>
#include <protocols/filters/CalculatorFilterCreator.hh>
#include <protocols/filters/BasicFilterCreators.hh>
#include <protocols/filters/BasicFilterCreators.hh>
#include <protocols/filters/ContingentFilterCreator.hh>
#include <protocols/filters/BasicFilterCreators.hh>
#include <protocols/filters/BasicFilterCreators.hh>
#include <protocols/filters/BasicFilterCreators.hh>
#include <protocols/filters/ReplicateFilterCreator.hh>
#include <protocols/filters/BasicFilterCreators.hh>
#include <protocols/filters/TimeFilterCreator.hh>
#include <protocols/simple_filters/RelativePoseFilterCreator.hh>
#include <protocols/simple_filters/StemFinderFilterCreator.hh>
#include <protocols/simple_filters/SSMotifFinderFilterCreator.hh>
#include <protocols/simple_filters/AngleToVectorFilterCreator.hh>
#include <protocols/fldsgn/filters/CoreDunbrackFilterCreator.hh>
#include <protocols/fldsgn/filters/FragQualFilterCreator.hh>
#include <protocols/fldsgn/filters/HelixKinkFilterCreator.hh>
#include <protocols/fldsgn/filters/HelixPairingFilterCreator.hh>
#include <protocols/fldsgn/filters/HSSTripletFilterCreator.hh>
#include <protocols/fldsgn/filters/InterlockingAromaFilterCreator.hh>
#include <protocols/fldsgn/filters/NcontactsFilterCreator.hh>
#include <protocols/fldsgn/filters/ParallelBetaPairingPreferenceFilterCreator.hh>
#include <protocols/fldsgn/filters/SecondaryStructureFilterCreator.hh>
#include <protocols/fldsgn/filters/SecondaryStructureCountFilterCreator.hh>
#include <protocols/fldsgn/filters/SecondaryStructureHasResidueFilterCreator.hh>
#include <protocols/fldsgn/filters/SheetTopologyFilterCreator.hh>
#include <protocols/denovo_design/filters/ExposedHydrophobicsFilterCreator.hh>
#include <protocols/denovo_design/filters/PreProlineFilterCreator.hh>
#include <protocols/denovo_design/filters/SSPredictionFilterCreator.hh>
#include <protocols/helical_bundle/BundleReporterFilterCreator.hh>
#include <protocols/indexed_structure_store/filters/FragmentLookupFilterCreator.hh>
#include <protocols/ligand_docking/AtomCountFilterCreator.hh>
#include <protocols/ligand_docking/CompleteConnectionsFilterCreator.hh>
#include <protocols/ligand_docking/ChainExistsFilterCreator.hh>
#include <protocols/ligand_docking/HeavyAtomFilterCreator.hh>
#include <protocols/ligand_docking/HBondAcceptorFilterCreator.hh>
#include <protocols/ligand_docking/HBondDonorFilterCreator.hh>
#include <protocols/ligand_docking/MolecularMassFilterCreator.hh>
#include <protocols/ligand_docking/MolarMassFilterCreator.hh>
#include <protocols/loops/filters/LoopAnalyzerFilterCreator.hh>
#include <protocols/matdes/ClashCheckFilterCreator.hh>
#include <protocols/matdes/GetRBDOFValuesCreator.hh>
#include <protocols/matdes/InterfacePackingFilterCreator.hh>
#include <protocols/matdes/OligomericAverageDegreeFilterCreator.hh>
#include <protocols/matdes/SymUnsatHbondFilterCreator.hh>
#include <protocols/protein_interface_design/filters/AtomicContactCountFilterCreator.hh>
#include <protocols/protein_interface_design/filters/AverageDegreeFilterCreator.hh>
#include <protocols/protein_interface_design/filters/BindingStrainFilterCreator.hh>
#include <protocols/protein_interface_design/filters/BoltzmannFilterCreator.hh>
#include <protocols/protein_interface_design/filters/DesignableResiduesFilterCreator.hh>
#include <protocols/protein_interface_design/filters/DisulfideFilterCreator.hh>
#include <protocols/protein_interface_design/filters/FilterScanCreator.hh>
#include <protocols/protein_interface_design/filters/HbondsToResidueFilterCreator.hh>
#include <protocols/protein_interface_design/filters/HbondsToAtomFilterCreator.hh>
#include <protocols/protein_interface_design/filters/InterfaceHolesFilterCreator.hh>
#include <protocols/protein_interface_design/filters/RelativeSegmentFilterCreator.hh>
#include <protocols/protein_interface_design/filters/RmsdFilterCreator.hh>
#include <protocols/protein_interface_design/filters/RmsdSimpleFilterCreator.hh>
#include <protocols/protein_interface_design/filters/ClashWithTargetFilterCreator.hh>
#include <protocols/protein_interface_design/filters/LRmsdFilterCreator.hh>
#include <protocols/protein_interface_design/filters/IRmsdFilterCreator.hh>
#include <protocols/protein_interface_design/filters/FNatFilterCreator.hh>
#include <protocols/protein_interface_design/filters/SequenceRecoveryFilterCreator.hh>
#include <protocols/protein_interface_design/filters/SpecificResiduesNearInterfaceFilterCreator.hh>
#include <protocols/protein_interface_design/filters/SSamountFilterCreator.hh>
#include <protocols/protein_interface_design/filters/StubScoreFilterCreator.hh>
#include <protocols/protein_interface_design/filters/StubScoreLoopsFilterCreator.hh>
#include <protocols/protein_interface_design/filters/TorsionFilterCreator.hh>
#include <protocols/simple_filters/AlaScanCreator.hh>
#include <protocols/simple_filters/AtomicContactFilterCreator.hh>
#include <protocols/simple_filters/AtomicDistanceFilterCreator.hh>
#include <protocols/simple_filters/AveragePathLengthFilterCreator.hh>
#include <protocols/simple_filters/BuriedUnsatHbondFilterCreator.hh>
#include <protocols/simple_filters/ConservedPosMutationFilterCreator.hh>
#include <protocols/simple_filters/ConstraintScoreFilterCreator.hh>
#include <protocols/simple_filters/DdgFilterCreator.hh>
#include <protocols/simple_filters/DeltaFilterCreator.hh>
#include <protocols/simple_filters/DisulfideEntropyFilterCreator.hh>
#include <protocols/simple_filters/EnergyPerResidueFilterCreator.hh>
#include <protocols/simple_filters/ExpiryFilterCreator.hh>
#include <protocols/simple_filters/NonSequentialNeighborsFilterCreator.hh>
#include <protocols/simple_filters/FileExistFilterCreator.hh>
#include <protocols/simple_filters/FileRemoveFilterCreator.hh>
#include <protocols/simple_filters/GeometryFilterCreator.hh>
#include <protocols/simple_filters/HolesFilterCreator.hh>
#include <protocols/simple_filters/InterfaceSasaFilterCreator.hh>
#include <protocols/simple_filters/InterfaceBindingEnergyDensityFilterCreator.hh>
#include <protocols/simple_filters/InterRepeatContactFilterCreator.hh>
#include <protocols/simple_filters/IntraRepeatContactFilterCreator.hh>
#include <protocols/simple_filters/LeastNativeLike9merFilterCreator.hh>
#include <protocols/simple_filters/MotifScoreFilterCreator.hh>
#include <protocols/simple_filters/MultipleSigmoidsFilterCreator.hh>
#include <protocols/simple_filters/MutationsFilterCreator.hh>
#include <protocols/simple_filters/NeighborTypeFilterCreator.hh>
#include <protocols/simple_filters/NetChargeFilterCreator.hh>
#include <protocols/simple_filters/NMerPSSMEnergyFilterCreator.hh>
#include <protocols/simple_filters/NMerSVMEnergyFilterCreator.hh>
#include <protocols/simple_filters/OperatorFilterCreator.hh>
#include <protocols/simple_filters/PackStatFilterCreator.hh>
#include <protocols/simple_filters/PoseCommentFilterCreator.hh>
#include <protocols/simple_filters/PoseInfoFilterCreator.hh>
#include <protocols/simple_filters/RangeFilterCreator.hh>
#include <protocols/simple_filters/ReportFilterCreator.hh>
#include <protocols/simple_filters/RepeatParameterFilterCreator.hh>
#include <protocols/simple_filters/ResidueCountFilterCreator.hh>
#include <protocols/simple_filters/ResidueDistanceFilterCreator.hh>
#include <protocols/simple_filters/ResidueDepthFilterCreator.hh>
#include <protocols/simple_filters/ResidueIEFilterCreator.hh>
#include <protocols/simple_filters/ResiduesInInterfaceFilterCreator.hh>
#include <protocols/simple_filters/ResidueSetChainEnergyFilterCreator.hh>
#include <protocols/simple_filters/RotamerBoltzmannWeightFilterCreator.hh>
#include <protocols/simple_filters/RotamerBoltzmannWeight2Creator.hh>
#include <protocols/simple_filters/SavePoseConstraintToFileFilterCreator.hh>
#include <protocols/simple_filters/SSElementMotifContactFilterCreator.hh>
#include <protocols/simple_filters/SaveResfileToDiskFilterCreator.hh>
#include <protocols/simple_filters/ScoreCutoffFilterCreator.hh>
#include <protocols/simple_filters/ScoreTypeFilterCreator.hh>
#include <protocols/simple_filters/ShapeComplementarityFilterCreator.hh>
#include <protocols/simple_filters/SigmoidFilterCreator.hh>
#include <protocols/simple_filters/SymmetricMotifFilterCreator.hh>
#include <protocols/simple_filters/SidechainRmsdFilterCreator.hh>
#include <protocols/simple_filters/DdGScanCreator.hh>
#include <protocols/simple_filters/TaskAwareSASAFilterCreator.hh>
#include <protocols/simple_filters/TaskAwareScoreTypeFilterCreator.hh>
#include <protocols/simple_filters/TerminusDistanceFilterCreator.hh>
#include <protocols/simple_filters/TotalSasaFilterCreator.hh>

class BackwardsProtocolsFilterCreatorTests : public CxxTest::TestSuite
{
public:

	void test_protocols_analysis_PeptideDeriverFilterCreator_name()
	{ protocols::analysis::PeptideDeriverFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "PeptideDeriver" ); }

	void test_protocols_cyclic_peptide_OversaturatedHbondAcceptorFilterCreator_name()
	{ protocols::cyclic_peptide::OversaturatedHbondAcceptorFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "OversaturatedHbondAcceptorFilter" ); }

	void test_protocols_enzdes_DiffAtomSasaFilterCreator_name()
	{ protocols::enzdes::DiffAtomSasaFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "DiffAtomBurial" ); }

	void test_protocols_enzdes_EnzdesScorefileFilterCreator_name()
	{ protocols::enzdes::EnzdesScorefileFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "EnzdesScorefileFilter" ); }

	void test_protocols_enzdes_EnzScoreFilterCreator_name()
	{ protocols::enzdes::EnzScoreFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "EnzScore" ); }

	void test_protocols_enzdes_ResidueConformerFilterCreator_name()
	{ protocols::enzdes::ResidueConformerFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ResidueConformerFilter" ); }

	void test_protocols_enzdes_LigBurialFilterCreator_name()
	{ protocols::enzdes::LigBurialFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "LigBurial" ); }

	void test_protocols_enzdes_LigDSasaFilterCreator_name()
	{ protocols::enzdes::LigDSasaFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "DSasa" ); }

	void test_protocols_enzdes_LigInterfaceEnergyFilterCreator_name()
	{ protocols::enzdes::LigInterfaceEnergyFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "LigInterfaceEnergy" ); }

	void test_protocols_enzdes_RepackWithoutLigandFilterCreator_name()
	{ protocols::enzdes::RepackWithoutLigandFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "RepackWithoutLigand" ); }

	void test_protocols_enzdes_RemoveLigandFilterCreator_name()
	{ protocols::enzdes::RemoveLigandFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "RemoveLigandFilter" ); }

	void test_protocols_filters_CalculatorFilterCreator_name()
	{ protocols::filters::CalculatorFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "CalculatorFilter" ); }

	void test_protocols_filters_CombinedFilterCreator_name()
	{ protocols::filters::CombinedFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "CombinedValue" ); }

	void test_protocols_filters_CompoundFilterCreator_name()
	{ protocols::filters::CompoundFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "CompoundStatement" ); }

	void test_protocols_filters_ContingentFilterCreator_name()
	{ protocols::filters::ContingentFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ContingentFilter" ); }

	void test_protocols_filters_FalseFilterCreator_name()
	{ protocols::filters::FalseFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "FalseFilter" ); }

	void test_protocols_filters_IfThenFilterCreator_name()
	{ protocols::filters::IfThenFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "IfThenFilter" ); }

	void test_protocols_filters_MoveBeforeFilterCreator_name()
	{ protocols::filters::MoveBeforeFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "MoveBeforeFilter" ); }

	void test_protocols_filters_ReplicateFilterCreator_name()
	{ protocols::filters::ReplicateFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ReplicateFilter" ); }

	void test_protocols_filters_StochasticFilterCreator_name()
	{ protocols::filters::StochasticFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "Stochastic" ); }

	void test_protocols_filters_TimeFilterCreator_name()
	{ protocols::filters::TimeFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "Time" ); }

	void test_protocols_simple_filters_RelativePoseFilterCreator_name()
	{ protocols::simple_filters::RelativePoseFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "RelativePose" ); }

	void test_protocols_simple_filters_StemFinderFilterCreator_name()
	{ protocols::simple_filters::StemFinderFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "StemFinder" ); }

	void test_protocols_simple_filters_SSMotifFinderFilterCreator_name()
	{ protocols::simple_filters::SSMotifFinderFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SSMotifFinder" ); }

	void test_protocols_simple_filters_AngleToVectorFilterCreator_name()
	{ protocols::simple_filters::AngleToVectorFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "AngleToVector" ); }

	void test_protocols_fldsgn_filters_CoreDunbrackFilterCreator_name()
	{ protocols::fldsgn::filters::CoreDunbrackFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "CoreDunbrack" ); }

	void test_protocols_fldsgn_filters_FragQualFilterCreator_name()
	{ protocols::fldsgn::filters::FragQualFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "FragQual" ); }

	void test_protocols_fldsgn_filters_HelixKinkFilterCreator_name()
	{ protocols::fldsgn::filters::HelixKinkFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "HelixKink" ); }

	void test_protocols_fldsgn_filters_HelixPairingFilterCreator_name()
	{ protocols::fldsgn::filters::HelixPairingFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "HelixPairing" ); }

	void test_protocols_fldsgn_filters_HSSTripletFilterCreator_name()
	{ protocols::fldsgn::filters::HSSTripletFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "HSSTriplet" ); }

	void test_protocols_fldsgn_filters_InterlockingAromaFilterCreator_name()
	{ protocols::fldsgn::filters::InterlockingAromaFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "InterlockingAroma" ); }

	void test_protocols_fldsgn_filters_NcontactsFilterCreator_name()
	{ protocols::fldsgn::filters::NcontactsFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "Ncontacts" ); }

	void test_protocols_fldsgn_filters_ParallelBetaPairingPreferenceFilterCreator_name()
	{ protocols::fldsgn::filters::ParallelBetaPairingPreferenceFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ParallelBetaPairingPreference" ); }

	void test_protocols_fldsgn_filters_SecondaryStructureFilterCreator_name()
	{ protocols::fldsgn::filters::SecondaryStructureFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SecondaryStructure" ); }

	void test_protocols_fldsgn_filters_SecondaryStructureCountFilterCreator_name()
	{ protocols::fldsgn::filters::SecondaryStructureCountFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SecondaryStructureCount" ); }

	void test_protocols_fldsgn_filters_SecondaryStructureHasResidueFilterCreator_name()
	{ protocols::fldsgn::filters::SecondaryStructureHasResidueFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SecondaryStructureHasResidue" ); }

	void test_protocols_fldsgn_filters_SheetTopologyFilterCreator_name()
	{ protocols::fldsgn::filters::SheetTopologyFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SheetTopology" ); }

	void test_protocols_denovo_design_filters_ExposedHydrophobicsFilterCreator_name()
	{ protocols::denovo_design::filters::ExposedHydrophobicsFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ExposedHydrophobics" ); }

	void test_protocols_denovo_design_filters_PreProlineFilterCreator_name()
	{ protocols::denovo_design::filters::PreProlineFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "PreProline" ); }

	void test_protocols_denovo_design_filters_SSPredictionFilterCreator_name()
	{ protocols::denovo_design::filters::SSPredictionFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SSPrediction" ); }

	void test_protocols_helical_bundle_BundleReporterFilterCreator_name()
	{ protocols::helical_bundle::BundleReporterFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "BundleReporter" ); }

	void test_protocols_indexed_structure_store_filters_FragmentLookupFilterCreator_name()
	{ protocols::indexed_structure_store::filters::FragmentLookupFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "FragmentLookupFilter" ); }

	void test_protocols_ligand_docking_AtomCountFilterCreator_name()
	{ protocols::ligand_docking::AtomCountFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "AtomCount" ); }

	void test_protocols_ligand_docking_CompleteConnectionsFilterCreator_name()
	{ protocols::ligand_docking::CompleteConnectionsFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "CompleteConnections" ); }

	void test_protocols_ligand_docking_ChainExistsFilterCreator_name()
	{ protocols::ligand_docking::ChainExistsFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ChainExists" ); }

	void test_protocols_ligand_docking_HeavyAtomFilterCreator_name()
	{ protocols::ligand_docking::HeavyAtomFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "HeavyAtom" ); }

	void test_protocols_ligand_docking_HBondAcceptorFilterCreator_name()
	{ protocols::ligand_docking::HBondAcceptorFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "HBondAcceptor" ); }

	void test_protocols_ligand_docking_HBondDonorFilterCreator_name()
	{ protocols::ligand_docking::HBondDonorFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "HBondDonor" ); }

	void test_protocols_ligand_docking_MolecularMassFilterCreator_name()
	{ protocols::ligand_docking::MolecularMassFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "MolecularMass" ); }

	void test_protocols_ligand_docking_MolarMassFilterCreator_name()
	{ protocols::ligand_docking::MolarMassFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "MolarMass" ); }

	void test_protocols_loops_filters_LoopAnalyzerFilterCreator_name()
	{ protocols::loops::filters::LoopAnalyzerFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "LoopAnalyzerFilter" ); }

	void test_protocols_matdes_ClashCheckFilterCreator_name()
	{ protocols::matdes::ClashCheckFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ClashCheck" ); }

	void test_protocols_matdes_GetRBDOFValuesCreator_name()
	{ protocols::matdes::GetRBDOFValuesCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "GetRBDOFValues" ); }

	void test_protocols_matdes_InterfacePackingFilterCreator_name()
	{ protocols::matdes::InterfacePackingFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "InterfacePacking" ); }

	void test_protocols_matdes_OligomericAverageDegreeFilterCreator_name()
	{ protocols::matdes::OligomericAverageDegreeFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "OligomericAverageDegree" ); }

	void test_protocols_matdes_SymUnsatHbondFilterCreator_name()
	{ protocols::matdes::SymUnsatHbondFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SymUnsatHbonds" ); }

	void test_protocols_protein_interface_design_filters_AtomicContactCountFilterCreator_name()
	{ protocols::protein_interface_design::filters::AtomicContactCountFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "AtomicContactCount" ); }

	void test_protocols_protein_interface_design_filters_AverageDegreeFilterCreator_name()
	{ protocols::protein_interface_design::filters::AverageDegreeFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "AverageDegree" ); }

	void test_protocols_protein_interface_design_filters_BindingStrainFilterCreator_name()
	{ protocols::protein_interface_design::filters::BindingStrainFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "BindingStrain" ); }

	void test_protocols_protein_interface_design_filters_BoltzmannFilterCreator_name()
	{ protocols::protein_interface_design::filters::BoltzmannFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "Boltzmann" ); }

	void test_protocols_protein_interface_design_filters_DesignableResiduesFilterCreator_name()
	{ protocols::protein_interface_design::filters::DesignableResiduesFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "DesignableResidues" ); }

	void test_protocols_protein_interface_design_filters_DisulfideFilterCreator_name()
	{ protocols::protein_interface_design::filters::DisulfideFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "DisulfideFilter" ); }

	void test_protocols_protein_interface_design_filters_FilterScanFilterCreator_name()
	{ protocols::protein_interface_design::filters::FilterScanFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "FilterScan" ); }

	void test_protocols_protein_interface_design_filters_HbondsToResidueFilterCreator_name()
	{ protocols::protein_interface_design::filters::HbondsToResidueFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "HbondsToResidue" ); }

	void test_protocols_protein_interface_design_filters_HbondsToAtomFilterCreator_name()
	{ protocols::protein_interface_design::filters::HbondsToAtomFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "HbondsToAtom" ); }

	void test_protocols_protein_interface_design_filters_InterfaceHolesFilterCreator_name()
	{ protocols::protein_interface_design::filters::InterfaceHolesFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "InterfaceHoles" ); }

	void test_protocols_protein_interface_design_filters_RelativeSegmentFilterCreator_name()
	{ protocols::protein_interface_design::filters::RelativeSegmentFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "RelativeSegment" ); }

	void test_protocols_protein_interface_design_filters_RmsdFilterCreator_name()
	{ protocols::protein_interface_design::filters::RmsdFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "Rmsd" ); }

	void test_protocols_protein_interface_design_filters_RmsdSimpleFilterCreator_name()
	{ protocols::protein_interface_design::filters::RmsdSimpleFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "RmsdSimple" ); }

	void test_protocols_protein_interface_design_filters_ClashWithTargetFilterCreator_name()
	{ protocols::protein_interface_design::filters::ClashWithTargetFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ClashWithTarget" ); }

	void test_protocols_protein_interface_design_filters_LRmsdFilterCreator_name()
	{ protocols::protein_interface_design::filters::LRmsdFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "LRmsd" ); }

	void test_protocols_protein_interface_design_filters_IRmsdFilterCreator_name()
	{ protocols::protein_interface_design::filters::IRmsdFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "IRmsd" ); }

	void test_protocols_protein_interface_design_filters_FNatFilterCreator_name()
	{ protocols::protein_interface_design::filters::FNatFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "FNat" ); }

	void test_protocols_protein_interface_design_filters_SequenceRecoveryFilterCreator_name()
	{ protocols::protein_interface_design::filters::SequenceRecoveryFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SequenceRecovery" ); }

	void test_protocols_protein_interface_design_filters_SpecificResiduesNearInterfaceFilterCreator_name()
	{ protocols::protein_interface_design::filters::SpecificResiduesNearInterfaceFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SpecificResiduesNearInterface" ); }

	void test_protocols_protein_interface_design_filters_SSamountFilterCreator_name()
	{ protocols::protein_interface_design::filters::SSamountFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SSamount" ); }

	void test_protocols_protein_interface_design_filters_StubScoreFilterCreator_name()
	{ protocols::protein_interface_design::filters::StubScoreFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "StubScore" ); }

	void test_protocols_protein_interface_design_filters_StubScoreLoopsFilterCreator_name()
	{ protocols::protein_interface_design::filters::StubScoreLoopsFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "StubScoreLoops" ); }

	void test_protocols_protein_interface_design_filters_TorsionCreator_name()
	{ protocols::protein_interface_design::filters::TorsionCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "Torsion" ); }

	void test_protocols_simple_filters_AlaScanFilterCreator_name()
	{ protocols::simple_filters::AlaScanFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "AlaScan" ); }

	void test_protocols_simple_filters_AtomicContactFilterCreator_name()
	{ protocols::simple_filters::AtomicContactFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "AtomicContact" ); }

	void test_protocols_simple_filters_AtomicDistanceFilterCreator_name()
	{ protocols::simple_filters::AtomicDistanceFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "AtomicDistance" ); }

	void test_protocols_simple_filters_AveragePathLengthFilterCreator_name()
	{ protocols::simple_filters::AveragePathLengthFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "AveragePathLength" ); }

	void test_protocols_simple_filters_BuriedUnsatHbondFilterCreator_name()
	{ protocols::simple_filters::BuriedUnsatHbondFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "BuriedUnsatHbonds" ); }

	void test_protocols_simple_filters_ConservedPosMutationFilterCreator_name()
	{ protocols::simple_filters::ConservedPosMutationFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ConservedPosMutationFilter" ); }

	void test_protocols_simple_filters_ConstraintScoreFilterCreator_name()
	{ protocols::simple_filters::ConstraintScoreFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ConstraintScore" ); }

	void test_protocols_simple_filters_DdgFilterCreator_name()
	{ protocols::simple_filters::DdgFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "Ddg" ); }

	void test_protocols_simple_filters_DeltaFilterCreator_name()
	{ protocols::simple_filters::DeltaFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "Delta" ); }

	void test_protocols_simple_filters_DisulfideEntropyFilterCreator_name()
	{ protocols::simple_filters::DisulfideEntropyFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "DisulfideEntropy" ); }

	void test_protocols_simple_filters_EnergyPerResidueFilterCreator_name()
	{ protocols::simple_filters::EnergyPerResidueFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "EnergyPerResidue" ); }

	void test_protocols_simple_filters_ExpiryFilterCreator_name()
	{ protocols::simple_filters::ExpiryFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "Expiry" ); }

	void test_protocols_simple_filters_NonSequentialNeighborsFilterCreator_name()
	{ protocols::simple_filters::NonSequentialNeighborsFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "NonSequentialNeighbors" ); }

	void test_protocols_simple_filters_FileExistFilterCreator_name()
	{ protocols::simple_filters::FileExistFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "FileExist" ); }

	void test_protocols_simple_filters_FileRemoveFilterCreator_name()
	{ protocols::simple_filters::FileRemoveFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "FileRemove" ); }

	void test_protocols_simple_filters_GeometryFilterCreator_name()
	{ protocols::simple_filters::GeometryFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "Geometry" ); }

	void test_protocols_simple_filters_HolesFilterCreator_name()
	{ protocols::simple_filters::HolesFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "Holes" ); }

	void test_protocols_simple_filters_InterfaceSasaFilterCreator_name()
	{ protocols::simple_filters::InterfaceSasaFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "Sasa" ); }

	void test_protocols_simple_filters_InterfaceBindingEnergyDensityFilterCreator_name()
	{ protocols::simple_filters::InterfaceBindingEnergyDensityFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "InterfaceBindingEnergyDensityFilter" ); }

	void test_protocols_simple_filters_InterRepeatContactFilterCreator_name()
	{ protocols::simple_filters::InterRepeatContactFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "InterRepeatContactsPerResidue" ); }

	void test_protocols_simple_filters_IntraRepeatContactFilterCreator_name()
	{ protocols::simple_filters::IntraRepeatContactFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "IntraRepeatContactsPerResidue" ); }

	void test_protocols_simple_filters_LeastNativeLike9merFilterCreator_name()
	{ protocols::simple_filters::LeastNativeLike9merFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "worst9mer" ); }

	void test_protocols_simple_filters_MotifScoreFilterCreator_name()
	{ protocols::simple_filters::MotifScoreFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "MotifScore" ); }

	void test_protocols_simple_filters_MultipleSigmoidsFilterCreator_name()
	{ protocols::simple_filters::MultipleSigmoidsFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "MultipleSigmoids" ); }

	void test_protocols_simple_filters_MutationsFilterCreator_name()
	{ protocols::simple_filters::MutationsFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "Mutations" ); }

	void test_protocols_simple_filters_NeighborTypeFilterCreator_name()
	{ protocols::simple_filters::NeighborTypeFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "NeighborType" ); }

	void test_protocols_simple_filters_NetChargeFilterCreator_name()
	{ protocols::simple_filters::NetChargeFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "NetCharge" ); }

	void test_protocols_simple_filters_NMerPSSMEnergyFilterCreator_name()
	{ protocols::simple_filters::NMerPSSMEnergyFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "NMerPSSMEnergy" ); }

	void test_protocols_simple_filters_NMerSVMEnergyFilterCreator_name()
	{ protocols::simple_filters::NMerSVMEnergyFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "NMerSVMEnergy" ); }

	void test_protocols_simple_filters_OperatorFilterCreator_name()
	{ protocols::simple_filters::OperatorFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "Operator" ); }

	void test_protocols_simple_filters_PackStatFilterCreator_name()
	{ protocols::simple_filters::PackStatFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "PackStat" ); }

	void test_protocols_simple_filters_PoseCommentFilterCreator_name()
	{ protocols::simple_filters::PoseCommentFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "PoseComment" ); }

	void test_protocols_simple_filters_PoseInfoFilterCreator_name()
	{ protocols::simple_filters::PoseInfoFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "PoseInfo" ); }

	void test_protocols_simple_filters_RangeFilterCreator_name()
	{ protocols::simple_filters::RangeFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "Range" ); }

	void test_protocols_simple_filters_ReportFilterCreator_name()
	{ protocols::simple_filters::ReportFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "Report" ); }

	void test_protocols_simple_filters_RepeatParameterFilterCreator_name()
	{ protocols::simple_filters::RepeatParameterFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "RepeatParameter" ); }

	void test_protocols_simple_filters_ResidueCountFilterCreator_name()
	{ protocols::simple_filters::ResidueCountFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ResidueCount" ); }

	void test_protocols_simple_filters_ResidueDistanceFilterCreator_name()
	{ protocols::simple_filters::ResidueDistanceFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ResidueDistance" ); }

	void test_protocols_simple_filters_ResidueDepthFilterCreator_name()
	{ protocols::simple_filters::ResidueDepthFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ResidueDepth" ); }

	void test_protocols_simple_filters_ResidueIEFilterCreator_name()
	{ protocols::simple_filters::ResidueIEFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ResidueIE" ); }

	void test_protocols_simple_filters_ResiduesInInterfaceFilterCreator_name()
	{ protocols::simple_filters::ResiduesInInterfaceFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ResInInterface" ); }

	void test_protocols_simple_filters_ResidueSetChainEnergyFilterCreator_name()
	{ protocols::simple_filters::ResidueSetChainEnergyFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ResidueSetChainEnergy" ); }

	void test_protocols_simple_filters_RotamerBoltzmannWeightFilterCreator_name()
	{ protocols::simple_filters::RotamerBoltzmannWeightFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "RotamerBoltzmannWeight" ); }

	void test_protocols_simple_filters_RotamerBoltzmannWeight2Creator_name()
	{ protocols::simple_filters::RotamerBoltzmannWeight2Creator cr; TS_ASSERT_EQUALS( cr.keyname(), "RotamerBoltzmannWeight2" ); }

	void test_protocols_simple_filters_SavePoseConstraintToFileFilterCreator_name()
	{ protocols::simple_filters::SavePoseConstraintToFileFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SavePoseConstraintToFile" ); }

	void test_protocols_simple_filters_SSElementMotifContactFilterCreator_name()
	{ protocols::simple_filters::SSElementMotifContactFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SSDegree" ); }

	void test_protocols_simple_filters_SaveResfileToDiskFilterCreator_name()
	{ protocols::simple_filters::SaveResfileToDiskFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SaveResfileToDisk" ); }

	void test_protocols_simple_filters_ScoreCutoffFilterCreator_name()
	{ protocols::simple_filters::ScoreCutoffFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ScoreCutoffFilter" ); }

	void test_protocols_simple_filters_ScoreTypeFilterCreator_name()
	{ protocols::simple_filters::ScoreTypeFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ScoreType" ); }

	void test_protocols_simple_filters_ShapeComplementarityFilterCreator_name()
	{ protocols::simple_filters::ShapeComplementarityFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ShapeComplementarity" ); }

	void test_protocols_simple_filters_SigmoidFilterCreator_name()
	{ protocols::simple_filters::SigmoidFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "Sigmoid" ); }

	void test_protocols_simple_filters_SymmetricMotifFilterCreator_name()
	{ protocols::simple_filters::SymmetricMotifFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SymmetricMotif" ); }

	void test_protocols_simple_filters_SidechainRmsdFilterCreator_name()
	{ protocols::simple_filters::SidechainRmsdFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SidechainRmsd" ); }

	void test_protocols_simple_filters_DdGScanCreator_name()
	{ protocols::simple_filters::DdGScanCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "DdGScan" ); }

	void test_protocols_simple_filters_TaskAwareSASAFilterCreator_name()
	{ protocols::simple_filters::TaskAwareSASAFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "TaskAwareSASA" ); }

	void test_protocols_simple_filters_TaskAwareScoreTypeFilterCreator_name()
	{ protocols::simple_filters::TaskAwareScoreTypeFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "TaskAwareScoreType" ); }

	void test_protocols_simple_filters_TerminusDistanceFilterCreator_name()
	{ protocols::simple_filters::TerminusDistanceFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "TerminusDistance" ); }

	void test_protocols_simple_filters_TotalSasaFilterCreator_name()
	{ protocols::simple_filters::TotalSasaFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "TotalSasa" ); }


	//void test_protocols_analysis_PeptideDeriverFilterCreator()
	//{ protocols::analysis::PeptideDeriverFilterCreator cr; std::cout << "protocols::analysis::PeptideDeriverFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_cyclic_peptide_OversaturatedHbondAcceptorFilterCreator()
	//{ protocols::cyclic_peptide::OversaturatedHbondAcceptorFilterCreator cr; std::cout << "protocols::cyclic_peptide::OversaturatedHbondAcceptorFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_enzdes_DiffAtomSasaFilterCreator()
	//{ protocols::enzdes::DiffAtomSasaFilterCreator cr; std::cout << "protocols::enzdes::DiffAtomSasaFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_enzdes_EnzdesScorefileFilterCreator()
	//{ protocols::enzdes::EnzdesScorefileFilterCreator cr; std::cout << "protocols::enzdes::EnzdesScorefileFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_enzdes_EnzScoreFilterCreator()
	//{ protocols::enzdes::EnzScoreFilterCreator cr; std::cout << "protocols::enzdes::EnzScoreFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_enzdes_ResidueConformerFilterCreator()
	//{ protocols::enzdes::ResidueConformerFilterCreator cr; std::cout << "protocols::enzdes::ResidueConformerFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_enzdes_LigBurialFilterCreator()
	//{ protocols::enzdes::LigBurialFilterCreator cr; std::cout << "protocols::enzdes::LigBurialFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_enzdes_LigDSasaFilterCreator()
	//{ protocols::enzdes::LigDSasaFilterCreator cr; std::cout << "protocols::enzdes::LigDSasaFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_enzdes_LigInterfaceEnergyFilterCreator()
	//{ protocols::enzdes::LigInterfaceEnergyFilterCreator cr; std::cout << "protocols::enzdes::LigInterfaceEnergyFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_enzdes_RepackWithoutLigandFilterCreator()
	//{ protocols::enzdes::RepackWithoutLigandFilterCreator cr; std::cout << "protocols::enzdes::RepackWithoutLigandFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_enzdes_RemoveLigandFilterCreator()
	//{ protocols::enzdes::RemoveLigandFilterCreator cr; std::cout << "protocols::enzdes::RemoveLigandFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_filters_CalculatorFilterCreator()
	//{ protocols::filters::CalculatorFilterCreator cr; std::cout << "protocols::filters::CalculatorFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_filters_CombinedFilterCreator()
	//{ protocols::filters::CombinedFilterCreator cr; std::cout << "protocols::filters::CombinedFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_filters_CompoundFilterCreator()
	//{ protocols::filters::CompoundFilterCreator cr; std::cout << "protocols::filters::CompoundFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_filters_ContingentFilterCreator()
	//{ protocols::filters::ContingentFilterCreator cr; std::cout << "protocols::filters::ContingentFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_filters_FalseFilterCreator()
	//{ protocols::filters::FalseFilterCreator cr; std::cout << "protocols::filters::FalseFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_filters_IfThenFilterCreator()
	//{ protocols::filters::IfThenFilterCreator cr; std::cout << "protocols::filters::IfThenFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_filters_MoveBeforeFilterCreator()
	//{ protocols::filters::MoveBeforeFilterCreator cr; std::cout << "protocols::filters::MoveBeforeFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_filters_ReplicateFilterCreator()
	//{ protocols::filters::ReplicateFilterCreator cr; std::cout << "protocols::filters::ReplicateFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_filters_StochasticFilterCreator()
	//{ protocols::filters::StochasticFilterCreator cr; std::cout << "protocols::filters::StochasticFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_filters_TimeFilterCreator()
	//{ protocols::filters::TimeFilterCreator cr; std::cout << "protocols::filters::TimeFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_RelativePoseFilterCreator()
	//{ protocols::simple_filters::RelativePoseFilterCreator cr; std::cout << "protocols::simple_filters::RelativePoseFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_StemFinderFilterCreator()
	//{ protocols::simple_filters::StemFinderFilterCreator cr; std::cout << "protocols::simple_filters::StemFinderFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_SSMotifFinderFilterCreator()
	//{ protocols::simple_filters::SSMotifFinderFilterCreator cr; std::cout << "protocols::simple_filters::SSMotifFinderFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_AngleToVectorFilterCreator()
	//{ protocols::simple_filters::AngleToVectorFilterCreator cr; std::cout << "protocols::simple_filters::AngleToVectorFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_fldsgn_filters_CoreDunbrackFilterCreator()
	//{ protocols::fldsgn::filters::CoreDunbrackFilterCreator cr; std::cout << "protocols::fldsgn::filters::CoreDunbrackFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_fldsgn_filters_FragQualFilterCreator()
	//{ protocols::fldsgn::filters::FragQualFilterCreator cr; std::cout << "protocols::fldsgn::filters::FragQualFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_fldsgn_filters_HelixKinkFilterCreator()
	//{ protocols::fldsgn::filters::HelixKinkFilterCreator cr; std::cout << "protocols::fldsgn::filters::HelixKinkFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_fldsgn_filters_HelixPairingFilterCreator()
	//{ protocols::fldsgn::filters::HelixPairingFilterCreator cr; std::cout << "protocols::fldsgn::filters::HelixPairingFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_fldsgn_filters_HSSTripletFilterCreator()
	//{ protocols::fldsgn::filters::HSSTripletFilterCreator cr; std::cout << "protocols::fldsgn::filters::HSSTripletFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_fldsgn_filters_InterlockingAromaFilterCreator()
	//{ protocols::fldsgn::filters::InterlockingAromaFilterCreator cr; std::cout << "protocols::fldsgn::filters::InterlockingAromaFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_fldsgn_filters_NcontactsFilterCreator()
	//{ protocols::fldsgn::filters::NcontactsFilterCreator cr; std::cout << "protocols::fldsgn::filters::NcontactsFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_fldsgn_filters_ParallelBetaPairingPreferenceFilterCreator()
	//{ protocols::fldsgn::filters::ParallelBetaPairingPreferenceFilterCreator cr; std::cout << "protocols::fldsgn::filters::ParallelBetaPairingPreferenceFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_fldsgn_filters_SecondaryStructureFilterCreator()
	//{ protocols::fldsgn::filters::SecondaryStructureFilterCreator cr; std::cout << "protocols::fldsgn::filters::SecondaryStructureFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_fldsgn_filters_SecondaryStructureCountFilterCreator()
	//{ protocols::fldsgn::filters::SecondaryStructureCountFilterCreator cr; std::cout << "protocols::fldsgn::filters::SecondaryStructureCountFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_fldsgn_filters_SecondaryStructureHasResidueFilterCreator()
	//{ protocols::fldsgn::filters::SecondaryStructureHasResidueFilterCreator cr; std::cout << "protocols::fldsgn::filters::SecondaryStructureHasResidueFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_fldsgn_filters_SheetTopologyFilterCreator()
	//{ protocols::fldsgn::filters::SheetTopologyFilterCreator cr; std::cout << "protocols::fldsgn::filters::SheetTopologyFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_denovo_design_filters_ExposedHydrophobicsFilterCreator()
	//{ protocols::denovo_design::filters::ExposedHydrophobicsFilterCreator cr; std::cout << "protocols::denovo_design::filters::ExposedHydrophobicsFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_denovo_design_filters_PreProlineFilterCreator()
	//{ protocols::denovo_design::filters::PreProlineFilterCreator cr; std::cout << "protocols::denovo_design::filters::PreProlineFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_denovo_design_filters_SSPredictionFilterCreator()
	//{ protocols::denovo_design::filters::SSPredictionFilterCreator cr; std::cout << "protocols::denovo_design::filters::SSPredictionFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_helical_bundle_BundleReporterFilterCreator()
	//{ protocols::helical_bundle::BundleReporterFilterCreator cr; std::cout << "protocols::helical_bundle::BundleReporterFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_indexed_structure_store_filters_FragmentLookupFilterCreator()
	//{ protocols::indexed_structure_store::filters::FragmentLookupFilterCreator cr; std::cout << "protocols::indexed_structure_store::filters::FragmentLookupFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_ligand_docking_AtomCountFilterCreator()
	//{ protocols::ligand_docking::AtomCountFilterCreator cr; std::cout << "protocols::ligand_docking::AtomCountFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_ligand_docking_CompleteConnectionsFilterCreator()
	//{ protocols::ligand_docking::CompleteConnectionsFilterCreator cr; std::cout << "protocols::ligand_docking::CompleteConnectionsFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_ligand_docking_ChainExistsFilterCreator()
	//{ protocols::ligand_docking::ChainExistsFilterCreator cr; std::cout << "protocols::ligand_docking::ChainExistsFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_ligand_docking_HeavyAtomFilterCreator()
	//{ protocols::ligand_docking::HeavyAtomFilterCreator cr; std::cout << "protocols::ligand_docking::HeavyAtomFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_ligand_docking_HBondAcceptorFilterCreator()
	//{ protocols::ligand_docking::HBondAcceptorFilterCreator cr; std::cout << "protocols::ligand_docking::HBondAcceptorFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_ligand_docking_HBondDonorFilterCreator()
	//{ protocols::ligand_docking::HBondDonorFilterCreator cr; std::cout << "protocols::ligand_docking::HBondDonorFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_ligand_docking_MolecularMassFilterCreator()
	//{ protocols::ligand_docking::MolecularMassFilterCreator cr; std::cout << "protocols::ligand_docking::MolecularMassFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_ligand_docking_MolarMassFilterCreator()
	//{ protocols::ligand_docking::MolarMassFilterCreator cr; std::cout << "protocols::ligand_docking::MolarMassFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_loops_filters_LoopAnalyzerFilterCreator()
	//{ protocols::loops::filters::LoopAnalyzerFilterCreator cr; std::cout << "protocols::loops::filters::LoopAnalyzerFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_matdes_ClashCheckFilterCreator()
	//{ protocols::matdes::ClashCheckFilterCreator cr; std::cout << "protocols::matdes::ClashCheckFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_matdes_GetRBDOFValuesCreator()
	//{ protocols::matdes::GetRBDOFValuesCreator cr; std::cout << "protocols::matdes::GetRBDOFValuesCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_matdes_InterfacePackingFilterCreator()
	//{ protocols::matdes::InterfacePackingFilterCreator cr; std::cout << "protocols::matdes::InterfacePackingFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_matdes_OligomericAverageDegreeFilterCreator()
	//{ protocols::matdes::OligomericAverageDegreeFilterCreator cr; std::cout << "protocols::matdes::OligomericAverageDegreeFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_matdes_SymUnsatHbondFilterCreator()
	//{ protocols::matdes::SymUnsatHbondFilterCreator cr; std::cout << "protocols::matdes::SymUnsatHbondFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_protein_interface_design_filters_AtomicContactCountFilterCreator()
	//{ protocols::protein_interface_design::filters::AtomicContactCountFilterCreator cr; std::cout << "protocols::protein_interface_design::filters::AtomicContactCountFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_protein_interface_design_filters_AverageDegreeFilterCreator()
	//{ protocols::protein_interface_design::filters::AverageDegreeFilterCreator cr; std::cout << "protocols::protein_interface_design::filters::AverageDegreeFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_protein_interface_design_filters_BindingStrainFilterCreator()
	//{ protocols::protein_interface_design::filters::BindingStrainFilterCreator cr; std::cout << "protocols::protein_interface_design::filters::BindingStrainFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_protein_interface_design_filters_BoltzmannFilterCreator()
	//{ protocols::protein_interface_design::filters::BoltzmannFilterCreator cr; std::cout << "protocols::protein_interface_design::filters::BoltzmannFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_protein_interface_design_filters_DesignableResiduesFilterCreator()
	//{ protocols::protein_interface_design::filters::DesignableResiduesFilterCreator cr; std::cout << "protocols::protein_interface_design::filters::DesignableResiduesFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_protein_interface_design_filters_DisulfideFilterCreator()
	//{ protocols::protein_interface_design::filters::DisulfideFilterCreator cr; std::cout << "protocols::protein_interface_design::filters::DisulfideFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_protein_interface_design_filters_FilterScanFilterCreator()
	//{ protocols::protein_interface_design::filters::FilterScanFilterCreator cr; std::cout << "protocols::protein_interface_design::filters::FilterScanFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_protein_interface_design_filters_HbondsToResidueFilterCreator()
	//{ protocols::protein_interface_design::filters::HbondsToResidueFilterCreator cr; std::cout << "protocols::protein_interface_design::filters::HbondsToResidueFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_protein_interface_design_filters_HbondsToAtomFilterCreator()
	//{ protocols::protein_interface_design::filters::HbondsToAtomFilterCreator cr; std::cout << "protocols::protein_interface_design::filters::HbondsToAtomFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_protein_interface_design_filters_InterfaceHolesFilterCreator()
	//{ protocols::protein_interface_design::filters::InterfaceHolesFilterCreator cr; std::cout << "protocols::protein_interface_design::filters::InterfaceHolesFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_protein_interface_design_filters_RelativeSegmentFilterCreator()
	//{ protocols::protein_interface_design::filters::RelativeSegmentFilterCreator cr; std::cout << "protocols::protein_interface_design::filters::RelativeSegmentFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_protein_interface_design_filters_RmsdFilterCreator()
	//{ protocols::protein_interface_design::filters::RmsdFilterCreator cr; std::cout << "protocols::protein_interface_design::filters::RmsdFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_protein_interface_design_filters_RmsdSimpleFilterCreator()
	//{ protocols::protein_interface_design::filters::RmsdSimpleFilterCreator cr; std::cout << "protocols::protein_interface_design::filters::RmsdSimpleFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_protein_interface_design_filters_ClashWithTargetFilterCreator()
	//{ protocols::protein_interface_design::filters::ClashWithTargetFilterCreator cr; std::cout << "protocols::protein_interface_design::filters::ClashWithTargetFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_protein_interface_design_filters_LRmsdFilterCreator()
	//{ protocols::protein_interface_design::filters::LRmsdFilterCreator cr; std::cout << "protocols::protein_interface_design::filters::LRmsdFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_protein_interface_design_filters_IRmsdFilterCreator()
	//{ protocols::protein_interface_design::filters::IRmsdFilterCreator cr; std::cout << "protocols::protein_interface_design::filters::IRmsdFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_protein_interface_design_filters_FNatFilterCreator()
	//{ protocols::protein_interface_design::filters::FNatFilterCreator cr; std::cout << "protocols::protein_interface_design::filters::FNatFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_protein_interface_design_filters_SequenceRecoveryFilterCreator()
	//{ protocols::protein_interface_design::filters::SequenceRecoveryFilterCreator cr; std::cout << "protocols::protein_interface_design::filters::SequenceRecoveryFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_protein_interface_design_filters_SpecificResiduesNearInterfaceFilterCreator()
	//{ protocols::protein_interface_design::filters::SpecificResiduesNearInterfaceFilterCreator cr; std::cout << "protocols::protein_interface_design::filters::SpecificResiduesNearInterfaceFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_protein_interface_design_filters_SSamountFilterCreator()
	//{ protocols::protein_interface_design::filters::SSamountFilterCreator cr; std::cout << "protocols::protein_interface_design::filters::SSamountFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_protein_interface_design_filters_StubScoreFilterCreator()
	//{ protocols::protein_interface_design::filters::StubScoreFilterCreator cr; std::cout << "protocols::protein_interface_design::filters::StubScoreFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_protein_interface_design_filters_StubScoreLoopsFilterCreator()
	//{ protocols::protein_interface_design::filters::StubScoreLoopsFilterCreator cr; std::cout << "protocols::protein_interface_design::filters::StubScoreLoopsFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_protein_interface_design_filters_TorsionCreator()
	//{ protocols::protein_interface_design::filters::TorsionCreator cr; std::cout << "protocols::protein_interface_design::filters::TorsionCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_AlaScanFilterCreator()
	//{ protocols::simple_filters::AlaScanFilterCreator cr; std::cout << "protocols::simple_filters::AlaScanFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_AtomicContactFilterCreator()
	//{ protocols::simple_filters::AtomicContactFilterCreator cr; std::cout << "protocols::simple_filters::AtomicContactFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_AtomicDistanceFilterCreator()
	//{ protocols::simple_filters::AtomicDistanceFilterCreator cr; std::cout << "protocols::simple_filters::AtomicDistanceFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_AveragePathLengthFilterCreator()
	//{ protocols::simple_filters::AveragePathLengthFilterCreator cr; std::cout << "protocols::simple_filters::AveragePathLengthFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_BuriedUnsatHbondFilterCreator()
	//{ protocols::simple_filters::BuriedUnsatHbondFilterCreator cr; std::cout << "protocols::simple_filters::BuriedUnsatHbondFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_ConservedPosMutationFilterCreator()
	//{ protocols::simple_filters::ConservedPosMutationFilterCreator cr; std::cout << "protocols::simple_filters::ConservedPosMutationFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_ConstraintScoreFilterCreator()
	//{ protocols::simple_filters::ConstraintScoreFilterCreator cr; std::cout << "protocols::simple_filters::ConstraintScoreFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_DdgFilterCreator()
	//{ protocols::simple_filters::DdgFilterCreator cr; std::cout << "protocols::simple_filters::DdgFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_DeltaFilterCreator()
	//{ protocols::simple_filters::DeltaFilterCreator cr; std::cout << "protocols::simple_filters::DeltaFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_DisulfideEntropyFilterCreator()
	//{ protocols::simple_filters::DisulfideEntropyFilterCreator cr; std::cout << "protocols::simple_filters::DisulfideEntropyFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_EnergyPerResidueFilterCreator()
	//{ protocols::simple_filters::EnergyPerResidueFilterCreator cr; std::cout << "protocols::simple_filters::EnergyPerResidueFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_ExpiryFilterCreator()
	//{ protocols::simple_filters::ExpiryFilterCreator cr; std::cout << "protocols::simple_filters::ExpiryFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_NonSequentialNeighborsFilterCreator()
	//{ protocols::simple_filters::NonSequentialNeighborsFilterCreator cr; std::cout << "protocols::simple_filters::NonSequentialNeighborsFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_FileExistFilterCreator()
	//{ protocols::simple_filters::FileExistFilterCreator cr; std::cout << "protocols::simple_filters::FileExistFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_FileRemoveFilterCreator()
	//{ protocols::simple_filters::FileRemoveFilterCreator cr; std::cout << "protocols::simple_filters::FileRemoveFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_GeometryFilterCreator()
	//{ protocols::simple_filters::GeometryFilterCreator cr; std::cout << "protocols::simple_filters::GeometryFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_HolesFilterCreator()
	//{ protocols::simple_filters::HolesFilterCreator cr; std::cout << "protocols::simple_filters::HolesFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_InterfaceSasaFilterCreator()
	//{ protocols::simple_filters::InterfaceSasaFilterCreator cr; std::cout << "protocols::simple_filters::InterfaceSasaFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_InterfaceBindingEnergyDensityFilterCreator()
	//{ protocols::simple_filters::InterfaceBindingEnergyDensityFilterCreator cr; std::cout << "protocols::simple_filters::InterfaceBindingEnergyDensityFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_InterRepeatContactFilterCreator()
	//{ protocols::simple_filters::InterRepeatContactFilterCreator cr; std::cout << "protocols::simple_filters::InterRepeatContactFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_IntraRepeatContactFilterCreator()
	//{ protocols::simple_filters::IntraRepeatContactFilterCreator cr; std::cout << "protocols::simple_filters::IntraRepeatContactFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_LeastNativeLike9merFilterCreator()
	//{ protocols::simple_filters::LeastNativeLike9merFilterCreator cr; std::cout << "protocols::simple_filters::LeastNativeLike9merFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_MotifScoreFilterCreator()
	//{ protocols::simple_filters::MotifScoreFilterCreator cr; std::cout << "protocols::simple_filters::MotifScoreFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_MultipleSigmoidsFilterCreator()
	//{ protocols::simple_filters::MultipleSigmoidsFilterCreator cr; std::cout << "protocols::simple_filters::MultipleSigmoidsFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_MutationsFilterCreator()
	//{ protocols::simple_filters::MutationsFilterCreator cr; std::cout << "protocols::simple_filters::MutationsFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_NeighborTypeFilterCreator()
	//{ protocols::simple_filters::NeighborTypeFilterCreator cr; std::cout << "protocols::simple_filters::NeighborTypeFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_NetChargeFilterCreator()
	//{ protocols::simple_filters::NetChargeFilterCreator cr; std::cout << "protocols::simple_filters::NetChargeFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_NMerPSSMEnergyFilterCreator()
	//{ protocols::simple_filters::NMerPSSMEnergyFilterCreator cr; std::cout << "protocols::simple_filters::NMerPSSMEnergyFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_NMerSVMEnergyFilterCreator()
	//{ protocols::simple_filters::NMerSVMEnergyFilterCreator cr; std::cout << "protocols::simple_filters::NMerSVMEnergyFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_OperatorFilterCreator()
	//{ protocols::simple_filters::OperatorFilterCreator cr; std::cout << "protocols::simple_filters::OperatorFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_PackStatFilterCreator()
	//{ protocols::simple_filters::PackStatFilterCreator cr; std::cout << "protocols::simple_filters::PackStatFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_PoseCommentFilterCreator()
	//{ protocols::simple_filters::PoseCommentFilterCreator cr; std::cout << "protocols::simple_filters::PoseCommentFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_PoseInfoFilterCreator()
	//{ protocols::simple_filters::PoseInfoFilterCreator cr; std::cout << "protocols::simple_filters::PoseInfoFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_RangeFilterCreator()
	//{ protocols::simple_filters::RangeFilterCreator cr; std::cout << "protocols::simple_filters::RangeFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_ReportFilterCreator()
	//{ protocols::simple_filters::ReportFilterCreator cr; std::cout << "protocols::simple_filters::ReportFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_RepeatParameterFilterCreator()
	//{ protocols::simple_filters::RepeatParameterFilterCreator cr; std::cout << "protocols::simple_filters::RepeatParameterFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_ResidueCountFilterCreator()
	//{ protocols::simple_filters::ResidueCountFilterCreator cr; std::cout << "protocols::simple_filters::ResidueCountFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_ResidueDistanceFilterCreator()
	//{ protocols::simple_filters::ResidueDistanceFilterCreator cr; std::cout << "protocols::simple_filters::ResidueDistanceFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_ResidueDepthFilterCreator()
	//{ protocols::simple_filters::ResidueDepthFilterCreator cr; std::cout << "protocols::simple_filters::ResidueDepthFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_ResidueIEFilterCreator()
	//{ protocols::simple_filters::ResidueIEFilterCreator cr; std::cout << "protocols::simple_filters::ResidueIEFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_ResiduesInInterfaceFilterCreator()
	//{ protocols::simple_filters::ResiduesInInterfaceFilterCreator cr; std::cout << "protocols::simple_filters::ResiduesInInterfaceFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_ResidueSetChainEnergyFilterCreator()
	//{ protocols::simple_filters::ResidueSetChainEnergyFilterCreator cr; std::cout << "protocols::simple_filters::ResidueSetChainEnergyFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_RotamerBoltzmannWeightFilterCreator()
	//{ protocols::simple_filters::RotamerBoltzmannWeightFilterCreator cr; std::cout << "protocols::simple_filters::RotamerBoltzmannWeightFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_RotamerBoltzmannWeight2Creator()
	//{ protocols::simple_filters::RotamerBoltzmannWeight2Creator cr; std::cout << "protocols::simple_filters::RotamerBoltzmannWeight2Creator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_SavePoseConstraintToFileFilterCreator()
	//{ protocols::simple_filters::SavePoseConstraintToFileFilterCreator cr; std::cout << "protocols::simple_filters::SavePoseConstraintToFileFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_SSElementMotifContactFilterCreator()
	//{ protocols::simple_filters::SSElementMotifContactFilterCreator cr; std::cout << "protocols::simple_filters::SSElementMotifContactFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_SaveResfileToDiskFilterCreator()
	//{ protocols::simple_filters::SaveResfileToDiskFilterCreator cr; std::cout << "protocols::simple_filters::SaveResfileToDiskFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_ScoreCutoffFilterCreator()
	//{ protocols::simple_filters::ScoreCutoffFilterCreator cr; std::cout << "protocols::simple_filters::ScoreCutoffFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_ScoreTypeFilterCreator()
	//{ protocols::simple_filters::ScoreTypeFilterCreator cr; std::cout << "protocols::simple_filters::ScoreTypeFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_ShapeComplementarityFilterCreator()
	//{ protocols::simple_filters::ShapeComplementarityFilterCreator cr; std::cout << "protocols::simple_filters::ShapeComplementarityFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_SigmoidFilterCreator()
	//{ protocols::simple_filters::SigmoidFilterCreator cr; std::cout << "protocols::simple_filters::SigmoidFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_SymmetricMotifFilterCreator()
	//{ protocols::simple_filters::SymmetricMotifFilterCreator cr; std::cout << "protocols::simple_filters::SymmetricMotifFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_SidechainRmsdFilterCreator()
	//{ protocols::simple_filters::SidechainRmsdFilterCreator cr; std::cout << "protocols::simple_filters::SidechainRmsdFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_DdGScanCreator()
	//{ protocols::simple_filters::DdGScanCreator cr; std::cout << "protocols::simple_filters::DdGScanCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_TaskAwareSASAFilterCreator()
	//{ protocols::simple_filters::TaskAwareSASAFilterCreator cr; std::cout << "protocols::simple_filters::TaskAwareSASAFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_TaskAwareScoreTypeFilterCreator()
	//{ protocols::simple_filters::TaskAwareScoreTypeFilterCreator cr; std::cout << "protocols::simple_filters::TaskAwareScoreTypeFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_TerminusDistanceFilterCreator()
	//{ protocols::simple_filters::TerminusDistanceFilterCreator cr; std::cout << "protocols::simple_filters::TerminusDistanceFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_protocols_simple_filters_TotalSasaFilterCreator()
	//{ protocols::simple_filters::TotalSasaFilterCreator cr; std::cout << "protocols::simple_filters::TotalSasaFilterCreator " << cr.keyname() << std::endl; }


};
