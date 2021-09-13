// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/protocols/scoring/InterfaceInfo.cxxtest.hh
/// @brief  test suite for protocols/scoring/InterfaceInfo.cc
/// @author Brahm Yachnin (byachnin@visterrainc.com)


// Test headers
#include <cxxtest/TestSuite.h>

#include <test/protocols/init_util.hh>


#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>

// Unit headers

// Package headers


//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <protocols/scoring/InterfaceInfo.hh>
#include <core/conformation/Residue.hh>

#include <basic/Tracer.hh> // AUTO IWYU For Error, Warning


using basic::Error;
using basic::Warning;

// using declarations
using namespace core;
using namespace core::pose;
using namespace core::conformation;
using namespace protocols::scoring;

/// @name InterfaceInfo
/// @brief: test the InterfaceInfo object, which holds Interfaces of a Pose

class InterfaceInfoTests : public CxxTest::TestSuite
{
public:
	PoseOP pose;
	PoseOP onechain_pose;

	void setUp() {
		protocols_init();

		// Need to score poses to get a neighbour graph
		core::scoring::ScoreFunction sfxn;
		sfxn.add_weights_from_file("ref2015");

		pose = utility::pointer::make_shared< Pose >();
		onechain_pose = utility::pointer::make_shared< Pose >();
		core::import_pose::pose_from_file( *pose, "protocols/scoring/dock_in.pdb" );
		core::import_pose::pose_from_file( *onechain_pose, "core/pose/onechain.pdb" );
		sfxn(*pose);
		sfxn(*onechain_pose);
	}

	void tearDown() {
		pose.reset();
		onechain_pose.reset();
	}

	/// @brief Tests for a normal, two-chain InterfaceInfo object
	void test_InterfaceInfoTest() {
		// Make an InterfaceInfo object
		InterfaceInfoOP iface_info( utility::pointer::make_shared<InterfaceInfo>() );

		// Check the interface distance.  Should be 6.0 (the default)
		TS_ASSERT_DELTA(iface_info->distance(), 6.0, 1e-6);
		// Check how many jumps we have.  Should be 1 (the default)
		TS_ASSERT_EQUALS(iface_info->num_jump(), 1);
		// Clear the jumps
		iface_info->clear_jumps();
		TS_ASSERT_EQUALS(iface_info->num_jump(), 0);
		// Add back jump 1
		iface_info->add_jump(1);
		TS_ASSERT_EQUALS(iface_info->num_jump(), 1);

		// Calculate the interfaces
		iface_info->calculate(*pose);

		// Get the number of interface residues for jump 1
		core::Size iface_nres = iface_info->interface_nres(1);
		TS_ASSERT_EQUALS(iface_nres, 32);

		// Two residues that are in contact should be in the interface
		core::Size res1_posenum(194);  // Chain E, Asp 194
		core::Size res2_posenum(263);  // Chain I, Tyr 18
		Residue res1 = pose->residue(res1_posenum);
		Residue res2 = pose->residue(res2_posenum);
		TS_ASSERT(iface_info->is_interface(res1));
		TS_ASSERT(iface_info->is_interface(res2));

		// They should also be a pair
		TS_ASSERT(iface_info->is_pair(res1, res2));

		// Pick some residues further away
		core::Size res_nearby_posenum(195);  // Chain E, Ser 195, at interface
		core::Size res_border_posenum(96);  // Chain E, Ser 96, just outside 6 A interface
		core::Size res_far_posenum(204);  // Chain E, Asn 204, far from interface

		Residue res_nearby = pose->residue(res_nearby_posenum);
		Residue res_border = pose->residue(res_border_posenum);
		Residue res_far = pose->residue(res_far_posenum);

		// res_nearby should be in the interface, while res_border and res_far should not
		TS_ASSERT(iface_info->is_interface(res_nearby));
		TS_ASSERT(!iface_info->is_interface(res_border));
		TS_ASSERT(!iface_info->is_interface(res_far));

		// Extend the interface to 8 A, and recalculate
		iface_info->distance(8);
		TS_ASSERT_DELTA(iface_info->distance(), 8.0, 1e-6);
		iface_info->calculate(*pose);

		// Should be more residues in the interface now
		core::Size bigger_iface_nres = iface_info->interface_nres(1);
		TS_ASSERT_EQUALS(bigger_iface_nres, 51);
		TS_ASSERT_LESS_THAN(iface_nres, bigger_iface_nres);

		// Now res_border should be in the interface, while res_far still shouldn't
		TS_ASSERT(iface_info->is_interface(res1));
		TS_ASSERT(iface_info->is_interface(res2));
		TS_ASSERT(iface_info->is_interface(res_nearby));
		TS_ASSERT(iface_info->is_interface(res_border));
		TS_ASSERT(!iface_info->is_interface(res_far));
	}

	/// @brief Tests for a one-chain InterfaceInfo object (i.e. there are no real Interfaces)
	void test_OneChainInterfaceInfoTest() {
		// Make an InterfaceInfo object
		InterfaceInfoOP iface_info( utility::pointer::make_shared<InterfaceInfo>() );

		// Check the interface distance.  Should be 6.0 (the default)
		TS_ASSERT_DELTA(iface_info->distance(), 6.0, 1e-6);
		// Check how many jumps we have.  Should be 1 (the default)
		TS_ASSERT_EQUALS(iface_info->num_jump(), 1);
		// Clear the jumps
		iface_info->clear_jumps();
		TS_ASSERT_EQUALS(iface_info->num_jump(), 0);
		// Add back jump 1
		iface_info->add_jump(1);
		TS_ASSERT_EQUALS(iface_info->num_jump(), 1);

		// Calculate the interfaces
		iface_info->calculate(*onechain_pose);

		// Jump 1 (which doesn't exist) shouldn't have any residues in the interface
		core::Size iface_nres = iface_info->interface_nres(1);
		TS_ASSERT_EQUALS(iface_nres, 0);

		// If a protein has no interfaces, its residues shouldn't be interface residues
		core::Size res1_posenum(20);
		Residue res1 = onechain_pose->residue(res1_posenum);
		TS_ASSERT(!iface_info->is_interface(res1));
	}

};
