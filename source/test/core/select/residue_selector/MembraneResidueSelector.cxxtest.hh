// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  core/select/residue_selector/MembraneResidueSeelctor.cxxtest.hh
/// @brief  Unit test for membrane residue selector
/// @author Rebecca Alford (rfalford12@gmail.com)

// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

// Project Headers
#include <core/select/residue_selector/MembraneResidueSelector.hh> 
#include <protocols/membrane/AddMembraneMover.hh> 
#include <core/conformation/Conformation.hh> 
#include <core/conformation/membrane/MembraneInfo.hh> 

// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Protocol Headers
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR("MembraneResidueSeelctor");

class MembraneResidueSeelctorTest : public CxxTest::TestSuite {

public:

	void setUp(){

		using namespace core::pose; 
		using namespace core::import_pose; 
		using namespace protocols::membrane; 

		core_init();

		// Setup DsbB (4 TM segments) and define transmembrane spans during add
		// membrane mover setup
		core::pose::PoseOP temp_pose = pose_from_file("core/select/residue_selector/DsbB_4tm_segs.pdb", core::import_pose::PDB_file);
		AddMembraneMoverOP add_memb( new AddMembraneMover( "core/select/residue_selector/DsbB_4tm_segs_B.span" ) );
		add_memb->apply( *temp_pose );
		pose_ = core::pose::PoseOP( new Pose( *temp_pose ) ); 

	}

	void tearDown(){

	}

	void test_membrane_residue_selection() {

		using namespace core::select::residue_selector; 
		TS_TRACE( "Testing selection of all transmembrane segments in the pose" );
		core::Size rsdnum( pose_->conformation().membrane_info()->membrane_rsd_num() );

		MembraneResidueSelectorOP memb_rsd_selector = MembraneResidueSelectorOP( new MembraneResidueSelector );
		ResidueSubset subset = memb_rsd_selector->apply( *pose_ ); 
		TS_ASSERT( subset[ rsdnum ] == true ); 
	}

private: 

	core::pose::PoseCOP pose_; 

};
