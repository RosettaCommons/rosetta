// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/scoring/RamaPrePro.cxxtest.hh
/// @brief  Unit tests for the RamaPrePro energy.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

// Project Headers
#include <core/scoring/RamaPrePro.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueProperties.hh>
#include <core/chemical/ResidueTypeSet.hh>

// Protocol Headers
#include <protocols/cyclic_peptide/FlipChiralityMover.hh>

// Basic Headers
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR("RamaPreProTests_symmgly");


class RamaPreProTests_symmgly : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init_with_additional_options("-symmetric_gly_tables true");
	}

	void tearDown(){
	}

	void do_test(
		std::string const &seq,
		bool const add_nmethyl,
		bool const flip_chirality
	) {
		core::pose::PoseOP pose( new core::pose::Pose );
		core::pose::make_pose_from_sequence(*pose, seq, "fa_standard", false);

		if ( add_nmethyl ) {
			core::chemical::ResidueTypeSetCOP rsd_set(pose->residue(2).residue_type_set());
			core::chemical::ResidueTypeCOP rsd_type( pose->residue_type(2).get_self_ptr() );
			core::chemical::ResidueTypeCOP new_rsd_type( rsd_set->get_residue_type_with_variant_added( *rsd_type,
				core::chemical::ResidueProperties::get_variant_from_string( "N_METHYLATION" ) ).get_self_ptr() );
			core::pose::replace_pose_residue_copying_existing_coordinates( *pose, 2, *new_rsd_type );
			//pose->dump_pdb("vtemp_ramaprepro.pdb"); //DELETE ME
		}

		if ( flip_chirality ) {
			protocols::cyclic_peptide::FlipChiralityMover flipper;
			flipper.apply( *pose );
		}

		core::scoring::RamaPrePro const & rama( core::scoring::ScoringManager::get_instance()->get_RamaPrePro() );
		TR << "\nPHI\tPSI\n";
		core::Size leftcount(0);
		for ( core::Size i=1; i<=1000; ++i ) {
			utility::vector1 < core::Real > phipsi;
			rama.random_mainchain_torsions( pose->residue_type(2).get_self_ptr(), pose->residue_type(3).get_self_ptr(), phipsi);
			TR << phipsi[1] << "\t" << phipsi[2] << "\n";
			if ( phipsi[1] <= 0 ) ++leftcount;
		}
		TR << "LEFT: " << leftcount << "\tRIGHT: " << 1000-leftcount << std::endl;
		TS_ASSERT( leftcount < 600 && leftcount > 400 );
	}

	/// @brief Test the drawing of random mainchain torsion values from the
	/// Ramachandran probability distribution of gly.
	void test_random_phipsi_gly() {
		do_test("AGAA", false, false);
	}

	/// @brief Test the drawing of random mainchain torsion values from the
	/// Ramachandran probability distribution for a pre-proline gly.
	void test_random_phipsi_gly_prepro() {
		do_test("AGPA", false, false);
	}

	/// @brief Test the drawing of random mainchain torsion values from the
	/// Ramachandran probability distribution of gly, with flipping.
	void test_random_phipsi_gly_invert() {
		do_test("AGAA", false, true);
	}

	/// @brief Test the drawing of random mainchain torsion values from the
	/// Ramachandran probability distribution for a pre-proline gly, with flipping.
	void test_random_phipsi_gly_prepro_invert() {
		do_test("AGPA", false, true);
	}

};



