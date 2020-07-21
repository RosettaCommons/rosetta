// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/backbone_moves/RandomizeBBByRamaPreProTests.cxxtest.hh
/// @brief  Unit tests for the RandomizeBBByRamaPrePro mover.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers
#include <protocols/backbone_moves/RandomizeBBByRamaPrePro.hh>

// Protocols Headers
#include <protocols/cyclic_peptide/PeptideStubMover.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>

// Numeric Headers
#include <numeric/angle.functions.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("RandomizeBBByRamaPreProTests");


class RandomizeBBByRamaPreProTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init();

	}

	void tearDown(){

	}

	/// @brief Apply the mover to a thousand alanines and ensure that negative phi gets sampled more than positive phi.
	void test_rama_distribution_l_ala() {
		TR << "Starting RandomizeBBByRamaPreProTests:test_rama_distribution_l_ala()" << std::endl;

		using namespace protocols::cyclic_peptide;
		using namespace core::select::residue_selector;
		using namespace protocols::backbone_moves;

		// Selector for alanines:
		ResidueIndexSelectorOP index_selector( utility::pointer::make_shared< ResidueIndexSelector >() );
		index_selector->set_index_range( 2, 21 );

		// Mover for building peptide:
		PeptideStubMover stubmover;
		stubmover.add_residue( PSM_append, "GLY", 1, true, "", 1, 0, nullptr, "" );
		for ( core::Size j(2); j<=21; ++j ) {
			stubmover.add_residue( PSM_append, "ALA", 1, false, "", 1, 0, nullptr, "" );
		}
		stubmover.add_residue( PSM_append, "GLY", 1, false, "", 1, 0, nullptr, "" );

		// Mover to test:
		RandomizeBBByRamaPrePro rand_bb;
		rand_bb.set_residue_selector( index_selector );

		// Counters for test:
		core::Size pos_phicount(0), neg_phicount(0);

		// Iterations of applying mover:
		for ( core::Size i(1); i<=50; ++i ) {
			core::pose::Pose pose;

			//Build a 20-residue ala pose, with glys on the end.
			{
				stubmover.apply( pose );
			}

			//Torsions to 180:
			for ( core::Size j(1); j<=21; ++j ) {
				pose.set_omega(j, 180);
			}
			pose.update_residue_neighbors();

			//Apply the mover:
			rand_bb.apply( pose );

			for ( core::Size j(2); j<=21; ++j ) {
				if ( numeric::principal_angle_degrees( pose.phi(j) ) < 0 ) {
					++neg_phicount;
				} else {
					++pos_phicount;
				}
			}
		}

		TR << "Positive phi count:\t" << pos_phicount << std::endl;
		TR << "Negative phi count:\t" << neg_phicount << std::endl;

		TS_ASSERT_EQUALS( pos_phicount + neg_phicount, 1000 ); //We should have done 1000 samples (50*20).
		TS_ASSERT_LESS_THAN( pos_phicount, neg_phicount ); //The negative phi count should be greater.
		TS_ASSERT_LESS_THAN( 20, pos_phicount ); //The positive phi count should be nonzero.
		TS_ASSERT_LESS_THAN( static_cast<core::Real>( pos_phicount ) / static_cast<core::Real>( pos_phicount + neg_phicount ), 0.15 ); //The positive phi count should be less than 15% of the total.

		TR << "Completed RandomizeBBByRamaPreProTests:test_rama_distribution_l_ala()" << std::endl;
	}

	/// @brief Apply the mover to a thousand D-alanines and ensure that positive phi gets sampled more than negative phi.
	void test_rama_distribution_d_ala() {
		TR << "Starting RandomizeBBByRamaPreProTests:test_rama_distribution_d_ala()" << std::endl;

		using namespace protocols::cyclic_peptide;
		using namespace core::select::residue_selector;
		using namespace protocols::backbone_moves;

		// Selector for alanines:
		ResidueIndexSelectorOP index_selector( utility::pointer::make_shared< ResidueIndexSelector >() );
		index_selector->set_index_range( 2, 21 );

		// Mover for building peptide:
		PeptideStubMover stubmover;
		stubmover.add_residue( PSM_append, "GLY", 1, true, "", 1, 0, nullptr, "" );
		for ( core::Size j(2); j<=21; ++j ) {
			stubmover.add_residue( PSM_append, "DALA", 1, false, "", 1, 0, nullptr, "" );
		}
		stubmover.add_residue( PSM_append, "GLY", 1, false, "", 1, 0, nullptr, "" );

		// Mover to test:
		RandomizeBBByRamaPrePro rand_bb;
		rand_bb.set_residue_selector( index_selector );

		// Counters for test:
		core::Size pos_phicount(0), neg_phicount(0);

		// Iterations of applying mover:
		for ( core::Size i(1); i<=50; ++i ) {
			core::pose::Pose pose;

			//Build a 20-residue ala pose, with glys on the end.
			{
				stubmover.apply( pose );
			}

			//Torsions to 180:
			for ( core::Size j(1); j<=21; ++j ) {
				pose.set_omega(j, 180);
			}
			pose.update_residue_neighbors();

			//Apply the mover:
			rand_bb.apply( pose );

			for ( core::Size j(2); j<=21; ++j ) {
				if ( numeric::principal_angle_degrees( pose.phi(j) ) < 0 ) {
					++neg_phicount;
				} else {
					++pos_phicount;
				}
			}
		}

		TR << "Positive phi count:\t" << pos_phicount << std::endl;
		TR << "Negative phi count:\t" << neg_phicount << std::endl;

		TS_ASSERT_EQUALS( pos_phicount + neg_phicount, 1000 ); //We should have done 1000 samples (50*20).
		TS_ASSERT_LESS_THAN( neg_phicount, pos_phicount ); //The positive phi count should be greater.
		TS_ASSERT_LESS_THAN( 20, neg_phicount ); //The negative phi count should be nonzero.
		TS_ASSERT_LESS_THAN( static_cast<core::Real>( neg_phicount ) / static_cast<core::Real>( pos_phicount + neg_phicount ), 0.15 ); //The negative phi count should be less than 15% of the total.

		TR << "Completed RandomizeBBByRamaPreProTests:test_rama_distribution_d_ala()" << std::endl;
	}

};
