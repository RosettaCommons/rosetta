// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 	 protocols/membrane/UniformPositionMover.cxxtest.hh
///
/// @brief 	 Unit Test: Uniform Position Mover
/// @details The add membrane mover sets up an anchored membrane fold tree, is responsible for
///			 adding the membrane residue and passing correct information to MembraneInfo for setup.
///			 Last Modified: 7/8/14
///
/// @author  Rebecca Alford (rfalford12@gmail.com)

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Package Headers
#include <protocols/simple_moves/UniformPositionMover.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/kinematics/FoldTree.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <core/types.hh>

// Utility Headers
#include <numeric/xyzVector.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/membrane_new.OptionKeys.gen.hh>

using namespace core;

class UniformPositionMoverTest : public CxxTest::TestSuite {
	
public: // test functions
    
    // Test Setup Functions ///////////////////////////
    
    /// @brief Setup Test
    void setUp() {
		
		using namespace core::pose;
		using namespace core::import_pose;
		using namespace basic::options;
        
        // Initialize core & options system
        core_init();
		
		// Set Anchored FoldTree Option
		option[ OptionKeys::membrane_new::anchored_foldtree ](true);
        
        // Load in pose from pdb
        pose_ = new Pose();
        pose_from_pdb( *pose_, "protocols/membrane/1C3W_TR_A.pdb" );
		
		// Add virtual atom to the root of the pose
		setup_virtual( *pose_ );
		
	}
	
	/// @brief Tear Down Test
	void tearDown() {}
	
	// Test Methods /////////////////////////////////
	
	/// @brief Test Uniform Translation (Check first coord)
	void test_uniform_translation() {
		
		TS_TRACE( "Testing uniform translation move" );
		
		using namespace protocols::simple_moves;
		
		// Move the CA of the first residue to the center
		Vector new_position(0, 0, 0);
		UniformPositionTranslationMoverOP translate = new UniformPositionTranslationMover( new_position, 1 );
		translate->apply( *pose_ );
		
		// Grab the first CA of the first pose residue
		Vector xyz1( pose_->residue( 1 ).atom( 2 ).xyz() );
		
		// Check first position stub center was shifted correctly
		TS_ASSERT( position_equal_within_delta( new_position, xyz1, 0.0001 ) );
		
		// Grab the second CA of the first pose residue
		Vector xyz2( pose_->residue( 2 ).atom( 2 ).xyz() );
		Vector new_position2(-0.071, 2.588, -2.753);
		TS_ASSERT( position_equal_within_delta( new_position2, xyz2, 0.0001 ) );
		
	}
	
	/// @brief Test Uniform Rotation (Check 2nd coord)
	void test_uniform_rotation() {
		
		TS_TRACE( "Testing uniform rotation move" );
		
		using namespace protocols::simple_moves;
		
		// Move the CA of the first residue to the center
		Vector axis(0, 0, 1);
		Real theta(2.289);
		UniformPositionRotationMoverOP rotate = new UniformPositionRotationMover( theta, axis, 1 );
		rotate->apply( *pose_ );
		
		// Grab the first CA of the first pose residue
		Vector new_position(29.5255, -15.4427, -3.8636 );
		Vector xyz1( pose_->residue( 1 ).atom( 2 ).xyz() );
		
		// Check first position stub center was shifted correctly
		TS_ASSERT( position_equal_within_delta( new_position, xyz1, 0.0001 ) );
		
		// Grab the second CA of the first pose residue
		Vector xyz2( pose_->residue( 2 ).atom( 2 ).xyz() );
		Vector new_position2(28.5557,  -18.1667, -6.2970);
		TS_ASSERT( position_equal_within_delta( new_position2, xyz2, 0.0001 ) );
		
	}
	
	
private:
	
	/// @brief Helper method - add virtual residue as the root of the pose
	void
	setup_virtual( core::pose::Pose & pose ) {
		
		using namespace core::conformation;
		using namespace core::chemical;
		using namespace core::kinematics;
		
		// Grab the current residue typeset and create a new residue
		ResidueTypeSetCAP const & residue_set(
											  ChemicalManager::get_instance()->residue_type_set( pose.is_fullatom() ? core::chemical::FA_STANDARD : core::chemical::CENTROID )
											  );
		
		// Create a new Residue from rsd typeset of type VRT
		ResidueTypeCOPs const & rsd_type_list( residue_set->name3_map("VRT") );
		ResidueType const & virt( *rsd_type_list[1] );
		ResidueOP rsd( ResidueFactory::create_residue( virt ) );
		
		// Append residue by jump, creating a new chain
		pose.append_residue_by_jump( *rsd, 1, "", "", true );
		
		// Make the anchoring residue the root of the fold tree
		core::Size const nres = pose.total_residue();
		FoldTree newF( pose.fold_tree() );
		newF.reorder( nres );
		pose.fold_tree( newF );
		
	}
	
	/// @brief Position equal within delta (helper method)
	bool position_equal_within_delta( Vector a, Vector b, Real delta ) {
		
		TS_ASSERT_DELTA( a.x(), b.x(), delta );
		TS_ASSERT_DELTA( a.y(), b.y(), delta );
		TS_ASSERT_DELTA( a.z(), b.z(), delta );
		
		return true;
	}
	
private:
	
	core::pose::PoseOP pose_;
	
}; // Uniform Position mover Unit Test
