// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/protocols/denovo_design/filters/PreProlineFilter.cxxtest.hh
/// @brief  test suite for protocols::denovo_design::filters::PreProlineFilter
/// @author Tom Linsky (tlinsky@uw.edu)


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Unit headers
#include <protocols/denovo_design/filters/PreProlineFilter.hh>

// Protocol headers
#include <protocols/denovo_design/util.hh>
#include <protocols/simple_moves/MutateResidue.hh>

// Core headers
#include <core/io/pdb/build_pose_as_is.hh>
#include <core/select/residue_selector/ChainSelector.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

// Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

// Boost headers
#include <boost/assign.hpp>

// C++ headers

static THREAD_LOCAL basic::Tracer TR( "protocols.denovo_design.PreProlineFilterTests.cxxtest" );

// --------------- Test Class --------------- //
class PreProlineFilterTests : public CxxTest::TestSuite {
public:

	// Shared initialization goes here.
	void setUp() {
		// load params for ligand
		protocols_init();

		// set preserve header always for "tomponent"
		basic::options::option[basic::options::OptionKeys::run::preserve_header].value(true);
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	void test_compute() {
		using namespace protocols::denovo_design;
		using namespace protocols::denovo_design::filters;

		core::pose::Pose input_pose;
		core::io::pdb::build_pose_from_pdb_as_is( input_pose, "protocols/optimize_weights/1UAD.wt_complex.pdb" );
		// obliterate pdb info so that residues are numbered sanely
		input_pose.pdb_info( core::pose::PDBInfoOP( new core::pose::PDBInfo( input_pose ) ) );

		PreProlineFilter prepro;
		prepro.set_use_statistical_potential( false );

		core::Real const bad_residues = prepro.compute( input_pose );
		// should be one violation -- a CIS-peptide between SER and PRO
		TS_ASSERT_DELTA( bad_residues, 1.0, 1e-3 );

		// now mutate a terminal residue to proline
		protocols::simple_moves::MutateResidue mut( 169, "PRO" );
		mut.apply( input_pose );
		core::Real const prochain_bad = prepro.compute( input_pose );

		// show be one violation still even though the residue before pro has O torsion
		TS_ASSERT_DELTA( prochain_bad, 1.0, 1e-3 );

		// now introduce a proline in a helix -- residue 91
		protocols::simple_moves::MutateResidue mut2( 91, "PRO" );
		mut2.apply( input_pose );
		core::Real const prohelix_bad = prepro.compute( input_pose );

		// should now be 2 violations
		TS_ASSERT_DELTA( prohelix_bad, 2.0, 1e-3 );

		// now select chain a
		core::select::residue_selector::ChainSelectorOP selchain(
			new core::select::residue_selector::ChainSelector() );
		utility::vector1< std::string > const chains = boost::assign::list_of ("A");
		selchain->set_chain_strings( chains );
		prepro.set_selector( selchain );
		core::Real const sel_bad = prepro.compute( input_pose );
		// should be one bad residue again
		TS_ASSERT_DELTA( sel_bad, 1.0, 1e-3 );
	}

	void test_spline() {
		using namespace protocols::denovo_design;
		using namespace protocols::denovo_design::filters;

		core::pose::Pose input_pose;
		core::io::pdb::build_pose_from_pdb_as_is( input_pose, "protocols/optimize_weights/1UAD.wt_complex.pdb" );
		// obliterate pdb info so that residues are numbered sanely
		input_pose.pdb_info( core::pose::PDBInfoOP( new core::pose::PDBInfo( input_pose ) ) );

		PreProlineFilter prepro;
		prepro.set_use_statistical_potential( true );
		core::Real const prepro_score = prepro.compute( input_pose );
		// score should be somewhere around 0
		TS_ASSERT_DELTA( prepro_score, -3.0, 1.0 );

		// mutate residue 91 to a pro, and we should get bad score
		protocols::simple_moves::MutateResidue mut2( 91, "PRO" );
		mut2.apply( input_pose );
		core::Real const prohelix_bad = prepro.compute( input_pose );
		TR << "before = " << prepro_score << " now score is " << prohelix_bad << std::endl;
		TS_ASSERT_DELTA( prohelix_bad, 1.0, 1.0 );
	}
};
