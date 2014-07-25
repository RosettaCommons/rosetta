// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 	 core/conformation/membrane/MembraneInfo.cxxtest.hh
///
/// @brief 	 Unit Test: MembraneInfo Object
/// @details The membrane info object is the central information store for
///			 data describing the conformation of a membrane protein. This includes:
///				(1) Transmembrane Spans
///				(2) Per-residue lipophilicity
///				(3) Position of the membrane residue & memrbane jump
///				(4) Optionally describing position of the membrane planes for visualization
///			 Last Modified: 7/8/14
///
/// @author  Rebecca Alford (rfalford12@gmail.com)

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Package Headers
#include <core/kinematics/FoldTree.hh>

#include <core/conformation/membrane/SpanningTopology.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/Conformation.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <core/types.hh>

// Utility Headers
#include <utility/vector1.hh>

using namespace core::kinematics;
using namespace core::conformation;
using namespace core::conformation::membrane;

class MembraneInfoTest : public CxxTest::TestSuite {
	
public: // test functions
    
    // Test Setup Functions ///////////////////////////
    
    /// @brief Setup Test
    void setUp() {
	
		using namespace core::import_pose;
        using namespace core::pose;
		
        // Initialize core & options system
        core_init();
        
        // Load in pose from pdb
        pose_ = new Pose();
        pose_from_pdb( *pose_, "protocols/membrane/1C3W_TR_A.pdb" );
        
		// Initialize Spans from spanfile
		std::string spanfile = "protocols/membrane/1C3W_A.span";
		SpanningTopologyOP topology = new SpanningTopology( spanfile, pose_->total_residue() );
		
		// Setup membrane info object
		pose_->conformation().setup_membrane(
			223, // position of the membrane residue
			topology, // spanning topology
			1 // membrane_jump
		);
	}
	
	/// @brief Tear Down Test
	void tearDown() {}
	
	// Test Methods /////////////////////////////////
	
	/// @brief Check conformation invariant passes when membraneInfo object is present
	void test_conformation_invariant() {
		
		TS_TRACE("Testing membrane conformation invariants");
        TS_ASSERT( pose_->conformation().is_membrane() );

	}
		
	/// @brief Check that MembraneInfo contains a reasonable setup
	void test_reasonable_setup() {
		
		TS_TRACE( "Testing a reasonable setup in MembraneInfo" );
		
		// Readability - grab object straight from conformation
		MembraneInfoOP membrane_info = pose_->conformation().membrane();
		
		// Thickness & steepness
		core::Real thickness = membrane_info->membrane_thickness();
		core::Real steepness = membrane_info->membrane_steepness();
		
		TS_ASSERT_EQUALS( thickness, 15 );
		TS_ASSERT_EQUALS( steepness, 10 );
		
		// Residue number and jump position
		core::Size resnum = membrane_info->membrane_rsd_num();
		core::Size jump = membrane_info->membrane_jump();
		
		TS_ASSERT_EQUALS( resnum, 223 );
		TS_ASSERT_EQUALS( jump, 1 );
		
		// Spanning topology (no lips info yet)
		core::Size total_spans = membrane_info->spanning_topology()->total_spans();
		TS_ASSERT_EQUALS( total_spans, 7 );
		
	}
		
private:

	core::pose::PoseOP pose_;
	
}; // MembraneInfo unit test
