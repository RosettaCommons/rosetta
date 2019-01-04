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
#include <devel/domain_insertion/FusePosesNtoCMoverCreator.hh>
#include <devel/domain_insertion/InsertionSiteTestMoverCreator.hh>
#include <devel/domain_insertion/FusePosesNtoCMoverCreator.hh>
#include <devel/enzdes/EnzdesRemodelMoverCreator.hh>
//#include <devel/vardist_solaccess/LoadVarSolDistSasaCalculatorMover.hh>
#include <devel/znhash/SymmZnMoversAndTaskOpsCreators.hh>
#include <devel/znhash/SymmZnMoversAndTaskOpsCreators.hh>
#include <devel/denovo_design/ConnectJumpsCreator.hh>
#include <devel/denovo_design/DumpStatsSSCreator.hh>
#include <devel/denovo_design/RestrictRegionCreator.hh>
#include <protocols/monte_carlo/GenericSimulatedAnnealerCreator.hh>
#include <devel/loop_creation/LoopCreationMoverCreator.hh>
#include <devel/loop_creation/FragmentLoopInserterCreator.hh>
#include <devel/loop_creation/CCDLoopCloserCreator.hh>
#include <devel/loop_creation/LoophashLoopInserterCreator.hh>
#include <devel/loop_creation/IterativeLoophashLoopInserterCreator.hh>
#include <devel/matdes/SymmetrizerMoverCreator.hh>
#include <devel/matdes/StoreQuasiSymmetricTaskMoverCreator.hh>
#include <devel/matdes/GenericSymmetricSamplerCreator.hh>
#include <devel/replica_docking/AddEncounterConstraintMoverCreator.hh>
#include <devel/replica_docking/ModulatedMoverCreator.hh>
#include <devel/loophash_loopclosure/LoopHashLoopClosureMoverCreator.hh>


class BackwardsDevelMoverCreatorTests : public CxxTest::TestSuite
{
public:

	void test_devel_domain_insertion_FusePosesNtoCMoverCreator_name()
	{ devel::domain_insertion::FusePosesNtoCMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "FusePosesNtoCMover" ); }

	void test_devel_domain_insertion_InsertionSiteTestMoverCreator_name()
	{ devel::domain_insertion::InsertionSiteTestMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "InsertionSiteTestMover" ); }

	void test_devel_domain_insertion_SetupCoiledCoilFoldTreeMoverCreator_name()
	{ devel::domain_insertion::SetupCoiledCoilFoldTreeMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "SetupCoiledCoilFoldTreeMover" ); }

	void test_devel_enzdes_EnzdesRemodelMoverCreator_name()
	{ devel::enzdes::EnzdesRemodelMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "EnzdesRemodelMover" ); }

	//void test_devel_vardist_solaccess_LoadVarSolDistSasaCalculatorMoverCreator_name()
	//{ devel::vardist_solaccess::LoadVarSolDistSasaCalculatorMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "LoadVarSolDistSasaCalculatorMover" ); }

	void test_devel_znhash_InsertZincCoordinationRemarkLinesCreator_name()
	{ devel::znhash::InsertZincCoordinationRemarkLinesCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "InsertZincCoordinationRemarkLines" ); }

	void test_devel_znhash_LoadZnCoordNumHbondCalculatorMoverCreator_name()
	{ devel::znhash::LoadZnCoordNumHbondCalculatorMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "LoadZnCoordNumHbondCalculatorMover" ); }

	void test_devel_denovo_design_ConnectJumpsCreator_name()
	{ devel::denovo_design::ConnectJumpsCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ConnectJumps" ); }

	void test_devel_denovo_design_DumpStatsSSCreator_name()
	{ devel::denovo_design::DumpStatsSSCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "DumpStatsSS" ); }

	void test_devel_denovo_design_RestrictRegionCreator_name()
	{ devel::denovo_design::RestrictRegionCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "RestrictRegion" ); }

	void test_devel_denovo_design_GenericSimulatedAnnealerCreator_name()
	{ protocols::monte_carlo::GenericSimulatedAnnealerCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "GenericSimulatedAnnealer" ); }

	void test_devel_loop_creation_LoopCreationMoverCreator_name()
	{ devel::loop_creation::LoopCreationMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "LoopCreationMover" ); }

	void test_devel_loop_creation_FragmentLoopInserterCreator_name()
	{ devel::loop_creation::FragmentLoopInserterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "FragmentLoopInserter" ); }

	void test_devel_loop_creation_CCDLoopCloserCreator_name()
	{ devel::loop_creation::CCDLoopCloserCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "CCDLoopCloser" ); }

	void test_devel_loop_creation_LoophashLoopInserterCreator_name()
	{ devel::loop_creation::LoophashLoopInserterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "LoophashLoopInserter" ); }

	void test_devel_loop_creation_IterativeLoophashLoopInserterCreator_name()
	{ devel::loop_creation::IterativeLoophashLoopInserterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "IterativeLoophashLoopInserter" ); }

	void test_devel_matdes_SymmetrizerMoverCreator_name()
	{ devel::matdes::SymmetrizerMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "Symmetrizer" ); }

	void test_devel_matdes_StoreQuasiSymmetricTaskMoverCreator_name()
	{ devel::matdes::StoreQuasiSymmetricTaskMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "StoreQuasiSymmetricTaskMover" ); }

	void test_devel_matdes_GenericSymmetricSamplerCreator_name()
	{ devel::matdes::GenericSymmetricSamplerCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "GenericSymmetricSampler" ); }

	void test_devel_replica_docking_AddEncounterConstraintMoverCreator_name()
	{ devel::replica_docking::AddEncounterConstraintMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "AddEncounterConstraintMover" ); }

	void test_devel_replica_docking_ModulatedMoverCreator_name()
	{ devel::replica_docking::ModulatedMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "ModulatedMover" ); }

	void test_devel_loophash_loopclosure_LoopHashLoopClosureMoverCreator_name()
	{ devel::loophash_loopclosure::LoopHashLoopClosureMoverCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "LoopHashLoopClosureMover" ); }





	// void test_devel_domain_insertion_FusePosesNtoCMoverCreator()
	// { devel::domain_insertion::FusePosesNtoCMoverCreator cr; std::cout << "devel::domain_insertion::FusePosesNtoCMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_devel_domain_insertion_InsertionSiteTestMoverCreator()
	// { devel::domain_insertion::InsertionSiteTestMoverCreator cr; std::cout << "devel::domain_insertion::InsertionSiteTestMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_devel_domain_insertion_SetupCoiledCoilFoldTreeMoverCreator()
	// { devel::domain_insertion::SetupCoiledCoilFoldTreeMoverCreator cr; std::cout << "devel::domain_insertion::SetupCoiledCoilFoldTreeMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_devel_enzdes_EnzdesRemodelMoverCreator()
	// { devel::enzdes::EnzdesRemodelMoverCreator cr; std::cout << "devel::enzdes::EnzdesRemodelMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_devel_vardist_solaccess_LoadVarSolDistSasaCalculatorMoverCreator()
	// { devel::vardist_solaccess::LoadVarSolDistSasaCalculatorMoverCreator cr; std::cout << "devel::vardist_solaccess::LoadVarSolDistSasaCalculatorMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_devel_znhash_InsertZincCoordinationRemarkLinesCreator()
	// { devel::znhash::InsertZincCoordinationRemarkLinesCreator cr; std::cout << "devel::znhash::InsertZincCoordinationRemarkLinesCreator " << cr.keyname() << std::endl; }
	//
	// void test_devel_znhash_LoadZnCoordNumHbondCalculatorMoverCreator()
	// { devel::znhash::LoadZnCoordNumHbondCalculatorMoverCreator cr; std::cout << "devel::znhash::LoadZnCoordNumHbondCalculatorMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_devel_denovo_design_ConnectJumpsCreator()
	// { devel::denovo_design::ConnectJumpsCreator cr; std::cout << "devel::denovo_design::ConnectJumpsCreator " << cr.keyname() << std::endl; }
	//
	// void test_devel_denovo_design_DumpStatsSSCreator()
	// { devel::denovo_design::DumpStatsSSCreator cr; std::cout << "devel::denovo_design::DumpStatsSSCreator " << cr.keyname() << std::endl; }
	//
	// void test_devel_denovo_design_RestrictRegionCreator()
	// { devel::denovo_design::RestrictRegionCreator cr; std::cout << "devel::denovo_design::RestrictRegionCreator " << cr.keyname() << std::endl; }
	//
	// void test_devel_denovo_design_GenericSimulatedAnnealerCreator()
	// { devel::denovo_design::GenericSimulatedAnnealerCreator cr; std::cout << "devel::denovo_design::GenericSimulatedAnnealerCreator " << cr.keyname() << std::endl; }
	//
	// void test_devel_loop_creation_LoopCreationMoverCreator()
	// { devel::loop_creation::LoopCreationMoverCreator cr; std::cout << "devel::loop_creation::LoopCreationMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_devel_loop_creation_FragmentLoopInserterCreator()
	// { devel::loop_creation::FragmentLoopInserterCreator cr; std::cout << "devel::loop_creation::FragmentLoopInserterCreator " << cr.keyname() << std::endl; }
	//
	// void test_devel_loop_creation_CCDLoopCloserCreator()
	// { devel::loop_creation::CCDLoopCloserCreator cr; std::cout << "devel::loop_creation::CCDLoopCloserCreator " << cr.keyname() << std::endl; }
	//
	// void test_devel_loop_creation_LoophashLoopInserterCreator()
	// { devel::loop_creation::LoophashLoopInserterCreator cr; std::cout << "devel::loop_creation::LoophashLoopInserterCreator " << cr.keyname() << std::endl; }
	//
	// void test_devel_loop_creation_IterativeLoophashLoopInserterCreator()
	// { devel::loop_creation::IterativeLoophashLoopInserterCreator cr; std::cout << "devel::loop_creation::IterativeLoophashLoopInserterCreator " << cr.keyname() << std::endl; }
	//
	// void test_devel_matdes_SymmetrizerMoverCreator()
	// { devel::matdes::SymmetrizerMoverCreator cr; std::cout << "devel::matdes::SymmetrizerMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_devel_matdes_StoreQuasiSymmetricTaskMoverCreator()
	// { devel::matdes::StoreQuasiSymmetricTaskMoverCreator cr; std::cout << "devel::matdes::StoreQuasiSymmetricTaskMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_devel_matdes_GenericSymmetricSamplerCreator()
	// { devel::matdes::GenericSymmetricSamplerCreator cr; std::cout << "devel::matdes::GenericSymmetricSamplerCreator " << cr.keyname() << std::endl; }
	//
	// void test_devel_replica_docking_AddEncounterConstraintMoverCreator()
	// { devel::replica_docking::AddEncounterConstraintMoverCreator cr; std::cout << "devel::replica_docking::AddEncounterConstraintMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_devel_replica_docking_ModulatedMoverCreator()
	// { devel::replica_docking::ModulatedMoverCreator cr; std::cout << "devel::replica_docking::ModulatedMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_devel_loophash_loopclosure_LoopHashLoopClosureMoverCreator()
	// { devel::loophash_loopclosure::LoopHashLoopClosureMoverCreator cr; std::cout << "devel::loophash_loopclosure::LoopHashLoopClosureMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_devel_splice_SpliceCreator()
	// { devel::splice::SpliceCreator cr; std::cout << "devel::splice::SpliceCreator " << cr.keyname() << std::endl; }
	//
	// void test_devel_splice_RBOutMoverCreator()
	// { devel::splice::RBOutMoverCreator cr; std::cout << "devel::splice::RBOutMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_devel_splice_RBInMoverCreator()
	// { devel::splice::RBInMoverCreator cr; std::cout << "devel::splice::RBInMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_devel_splice_AlignEndsMoverCreator()
	// { devel::splice::AlignEndsMoverCreator cr; std::cout << "devel::splice::AlignEndsMoverCreator " << cr.keyname() << std::endl; }
	//
	// void test_devel_cutoutdomain_CutOutDomainCreator()
	// { devel::cutoutdomain::CutOutDomainCreator cr; std::cout << "devel::cutoutdomain::CutOutDomainCreator " << cr.keyname() << std::endl; }

};
