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
#include <devel/denovo_design/filters/CavityVolumeFilterCreator.hh>
#include <devel/denovo_design/filters/CoreResiduesPerElementFilterCreator.hh>
#include <devel/denovo_design/filters/FoldabilityFilterCreator.hh>
#include <devel/replica_docking/InteractionScoreFilterCreator.hh>
#include <devel/replica_docking/IrmsdFilterCreator.hh>
#include <devel/replica_docking/FnatFilterCreator.hh>
#include <devel/replica_docking/LrmsdFilterCreator.hh>
#include <devel/replica_docking/FnonnatFilterCreator.hh>
#include <devel/replica_docking/CaIrmsdFilterCreator.hh>
#include <devel/buns/BuriedUnsatHbondFilter2Creator.hh>

class BackwardsDevelFilterCreatorTests : public CxxTest::TestSuite
{
public:

	void test_devel_denovo_design_filters_CavityVolumeFilterCreator_name()
	{ devel::denovo_design::filters::CavityVolumeFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "CavityVolume" ); }

	void test_devel_denovo_design_filters_CoreResiduesPerElementFilterCreator_name()
	{ devel::denovo_design::filters::CoreResiduesPerElementFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "CoreResiduesPerElement" ); }

	void test_devel_denovo_design_filters_FoldabilityFilterCreator_name()
	{ devel::denovo_design::filters::FoldabilityFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "Foldability" ); }

	void test_devel_replica_docking_InteractionScoreFilterCreator_name()
	{ devel::replica_docking::InteractionScoreFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "I_sc" ); }

	void test_devel_replica_docking_IrmsdFilterCreator_name()
	{ devel::replica_docking::IrmsdFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "Irms" ); }

	void test_devel_replica_docking_FnatFilterCreator_name()
	{ devel::replica_docking::FnatFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "Fnat" ); }

	void test_devel_replica_docking_LrmsdFilterCreator_name()
	{ devel::replica_docking::LrmsdFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "Lrmsd" ); }

	void test_devel_replica_docking_FnonnatFilterCreator_name()
	{ devel::replica_docking::FnonnatFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "Fnonnat" ); }

	void test_devel_replica_docking_CaIrmsdFilterCreator_name()
	{ devel::replica_docking::CaIrmsdFilterCreator cr; TS_ASSERT_EQUALS( cr.keyname(), "Ca_Irms" ); }

	void test_devel_buns_BuriedUnsatHbondFilter2Creator_name()
	{ devel::buns::BuriedUnsatHbondFilter2Creator cr; TS_ASSERT_EQUALS( cr.keyname(), "BuriedUnsatHbonds2" ); }


	//void test_devel_denovo_design_filters_CavityVolumeFilterCreator()
	//{ devel::denovo_design::filters::CavityVolumeFilterCreator cr; std::cout << "devel::denovo_design::filters::CavityVolumeFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_devel_denovo_design_filters_CoreResiduesPerElementFilterCreator()
	//{ devel::denovo_design::filters::CoreResiduesPerElementFilterCreator cr; std::cout << "devel::denovo_design::filters::CoreResiduesPerElementFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_devel_denovo_design_filters_FoldabilityFilterCreator()
	//{ devel::denovo_design::filters::FoldabilityFilterCreator cr; std::cout << "devel::denovo_design::filters::FoldabilityFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_devel_replica_docking_InteractionScoreFilterCreator()
	//{ devel::replica_docking::InteractionScoreFilterCreator cr; std::cout << "devel::replica_docking::InteractionScoreFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_devel_replica_docking_IrmsdFilterCreator()
	//{ devel::replica_docking::IrmsdFilterCreator cr; std::cout << "devel::replica_docking::IrmsdFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_devel_replica_docking_FnatFilterCreator()
	//{ devel::replica_docking::FnatFilterCreator cr; std::cout << "devel::replica_docking::FnatFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_devel_replica_docking_LrmsdFilterCreator()
	//{ devel::replica_docking::LrmsdFilterCreator cr; std::cout << "devel::replica_docking::LrmsdFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_devel_replica_docking_FnonnatFilterCreator()
	//{ devel::replica_docking::FnonnatFilterCreator cr; std::cout << "devel::replica_docking::FnonnatFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_devel_replica_docking_CaIrmsdFilterCreator()
	//{ devel::replica_docking::CaIrmsdFilterCreator cr; std::cout << "devel::replica_docking::CaIrmsdFilterCreator " << cr.keyname() << std::endl; }
	//
	//void test_devel_buns_BuriedUnsatHbondFilter2Creator()
	//{ devel::buns::BuriedUnsatHbondFilter2Creator cr; std::cout << "devel::buns::BuriedUnsatHbondFilter2Creator " << cr.keyname() << std::endl; }


};
