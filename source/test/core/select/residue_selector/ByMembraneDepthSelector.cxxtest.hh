// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  core/select/residue_selector/ByMembraneDepthSelector.cxxtest.hh
/// @brief  Select residues based upon their depth in the membrane, relative to the center
/// @author Rebecca Alford (rfalford12@gmail.com)

// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

// Project Headers
#include <core/select/residue_selector/ByMembraneDepthSelector.hh> 
#include <core/select/residue_selector/ResidueSelector.hh> 

#include <protocols/membrane/AddMembraneMover.hh>
#include <core/conformation/Conformation.hh> 
#include <core/conformation/membrane/MembraneInfo.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Utility Headers
#include <basic/Tracer.hh>

// C++ Headers
#include <string> 

static THREAD_LOCAL basic::Tracer TR( "ByMembraneDepthSelector" );

class ByMembraneDepthSelector : public CxxTest::TestSuite {

public:

	void setUp(){
	
		using namespace core::pose;
		using namespace core::import_pose; 
		using namespace protocols::membrane; 

		core_init();

		// Setting up the pose (using residues in 1AFO)
		core::pose::PoseOP temp_pose = pose_from_file("core/conformation/membrane/1AFO_AB.pdb", core::import_pose::PDB_file);
		AddMembraneMoverOP add_memb( new AddMembraneMover( "core/conformation/membrane/1AFO_AB.span" ) );
		add_memb->apply( *temp_pose );
		pose_ = core::pose::PoseOP( new Pose( *temp_pose ) );
		
	}

	void tearDown(){

	}
	
public: // test functions

	// Test that under default conditions, all residues in 1AFO are selected
	void test_complete_selection() {
	
		using namespace core::select::residue_selector;
		TS_TRACE( "Testing selection of all residues in the pose based on full range (unless the pose is inordinately large, which it is not" );
 
		ByMembraneDepthSelectorOP membrane_rsd_selector( new core::select::residue_selector::ByMembraneDepthSelector() );
		ResidueSubset subset = membrane_rsd_selector->apply( *pose_ );
		for ( core::Size ii = 1; ii <= subset.size(); ++ii ) {
			TS_ASSERT( subset[ii] );
		}
	}
	
	// Test that bounds are set correctly by getters and setters
	void test_correct_getter_setter_behavior() {
	
		using namespace core::select::residue_selector;
		TS_TRACE( "Testing that getters and setters result in reasonable values" );
		
		 
		ByMembraneDepthSelectorOP membrane_rsd_selector( new core::select::residue_selector::ByMembraneDepthSelector() );
	
		// Check a set of default values
		TS_ASSERT( membrane_rsd_selector->minimum_depth() == 0 );
		TS_ASSERT( membrane_rsd_selector->maximum_depth() == 100000 );
		
		// Check a set of mofied values
		membrane_rsd_selector->minimum_depth( 50 );
		membrane_rsd_selector->maximum_depth( 51 );
		TS_ASSERT( membrane_rsd_selector->minimum_depth() == 50 );
		TS_ASSERT( membrane_rsd_selector->maximum_depth() == 51 );

	}
	
	// Test that residues on both sides of the center plane are selected given a zero bounded selection range
	void test_zero_bounded_selection() {
	
		using namespace core::select::residue_selector;
		
		ByMembraneDepthSelectorOP membrane_rsd_selector( new core::select::residue_selector::ByMembraneDepthSelector() );
		membrane_rsd_selector->minimum_depth( 0 );
		membrane_rsd_selector->maximum_depth( 5 );
		ResidueSubset subset = membrane_rsd_selector->apply( *pose_ );
		
		// ones that should be true
		TS_ASSERT( subset[1] );
		TS_ASSERT( subset[20] );
		TS_ASSERT( subset[21] );
		TS_ASSERT( subset[22] );
		TS_ASSERT( subset[23] );
		TS_ASSERT( subset[24] );
		TS_ASSERT( subset[25] );
		TS_ASSERT( subset[26] );
		TS_ASSERT( subset[60] );
		TS_ASSERT( subset[61] );
		TS_ASSERT( subset[62] );
		TS_ASSERT( subset[63] );
		TS_ASSERT( subset[64] );
		TS_ASSERT( subset[65] );
		TS_ASSERT( subset[66] );
		TS_ASSERT( subset[67] );
		TS_ASSERT( subset[68] );
		
		// Random ones that should not be true
		TS_ASSERT( !subset[10] );
		TS_ASSERT( !subset[18] );
		TS_ASSERT( !subset[30] );
		TS_ASSERT( !subset[40] );
		TS_ASSERT( !subset[46] );
		TS_ASSERT( !subset[52] );
		TS_ASSERT( !subset[56] );
		TS_ASSERT( !subset[80] );
	
	}

private:

	core::pose::PoseCOP pose_;

};
