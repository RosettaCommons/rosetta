// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/chemical/VariantTypeTests.cxxtest.hh
/// @brief  Testing VariantType code.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers


// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ResidueType.hh>

// Protocols Headers
#include <protocols/simple_moves/ModifyVariantTypeMover.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("VariantTypeTests");


class VariantTypeTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp() {
		core_init();

	}

	void tearDown() {

	}

	void test_general_residue_virtualization() {
		core::pose::PoseOP pose( core::import_pose::pose_from_file( "core/pack/guidance_scoreterms/voids_penalty_energy/1ubq.pdb" ) );
		protocols::simple_moves::ModifyVariantTypeMover add_virt;
		add_virt.set_additional_type_to_add( "VIRTUAL_RESIDUE_VARIANT" );
		add_virt.apply( *pose );
		for ( core::Size ir(1), irmax(pose->total_residue()); ir<=irmax; ++ir ) {
			for ( core::Size ia(1), iamax(pose->residue(ir).natoms()); ia<=iamax; ++ia ) {
				TS_ASSERT( pose->residue(ir).atom_type(ia).is_virtual() );
			}
			TS_ASSERT( pose->residue_type(ir).is_virtual_residue() || pose->residue_type(ir).is_inverted_virtual_residue() );
		}
	}


};
