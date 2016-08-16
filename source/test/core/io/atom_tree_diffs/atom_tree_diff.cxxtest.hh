// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/io/atom_tree_diffs/atom_tree_diff.cxxtest.hh
/// @brief  test suite for atom_tree_diff file format
/// @author Ian Davis

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>
#include <test/UTracer.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
// Auto-header: duplicate removed #include <core/conformation/Residue.hh>
#include <core/import_pose/atom_tree_diffs/atom_tree_diff.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

//Auto Headers
#include <utility/vector1.hh>


static basic::Tracer TR("core.import_pose.atom_tree_diffs.atom_tree_diff.cxxtest");

using namespace core;

class AtomTreeDiffTests : public CxxTest::TestSuite
{
public:
	AtomTreeDiffTests() {};

	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	// Shared finalization goes here.
	void tearDown() {
	}


	void test_save_and_restore()
	{
		//test::UTracer UT("core/io/atom_tree_diffs/atom_tree_diff_test.u", &TR);

		/// Init random generators system to insure that random number sequence is the same for each run.
		core::init::init_random_generators(1000, "mt19937");


		pose::Pose start_pose( create_test_in_pdb_pose()), modified_pose, restored_pose;
		//core::import_pose::pose_from_file( start_pose, "core/io/test_in.pdb" , core::import_pose::PDB_file);
		//UTRACE << pose.fold_tree() << std::endl;

		// Randomize the pose a little
		using numeric::random::rg;
		modified_pose = start_pose; // make a copy
		for ( int i = 0; i < 20; ++i ) {
			Size rsd_no = rg().random_range( 1, modified_pose.total_residue() );
			modified_pose.set_phi( rsd_no, 10*rg().gaussian() + modified_pose.phi(rsd_no) );
			modified_pose.set_psi( rsd_no, 10*rg().gaussian() + modified_pose.psi(rsd_no) );
			Size num_chi = modified_pose.residue(rsd_no).nchi();
			for ( Size j = 1; j < num_chi; ++j ) {
				modified_pose.set_chi( j, rsd_no, 10*rg().gaussian() + modified_pose.chi(j, rsd_no) );
			}
		}

		// Verify that it's not like the input structure any more
		Real rms_to_orig = scoring::rmsd_no_super(start_pose, modified_pose, scoring::is_heavyatom);
		TR << "RMS to original: " << rms_to_orig << std::endl;
		TS_ASSERT( rms_to_orig > 1.0 );

		// Now mutate it a little for good measure (can't do RMS to orig after this)
		for ( int i = 0; i < 10; ++i ) {
			// Don't mutate the ends because they're variant residue types
			Size rsd_no = rg().random_range( 2, modified_pose.total_residue()-1 );
			using namespace core::conformation;
			ResidueOP newres = ResidueFactory::create_residue(
				modified_pose.residue(rsd_no).residue_type_set()->name_map("LYS"),
				modified_pose.residue(rsd_no), modified_pose.conformation());
			modified_pose.replace_residue(rsd_no, *newres, true /*orient backbone*/);
			// Change chi angles for mutated res away from their default values
			Size num_chi = modified_pose.residue(rsd_no).nchi();
			for ( Size j = 1; j < num_chi; ++j ) {
				modified_pose.set_chi( j, rsd_no, 90*rg().gaussian() + modified_pose.chi(j, rsd_no) );
			}
		}

		// Serialize the modified structure as a atom_tree_diff file
		std::ostringstream outss;
		std::map< std::string, core::Real > my_scores; // empty
		scoring::ScoreFunctionOP sfxn = scoring::get_score_function();
		core::import_pose::atom_tree_diffs::map_of_weighted_scores(modified_pose, *sfxn, my_scores);
		core::import_pose::atom_tree_diffs::dump_atom_tree_diff(outss, "tag", my_scores, start_pose, modified_pose);

		// Restore headers and structure from atom_tree_diff file, assert minimal error introduced.
		std::string tag_out;
		std::map< std::string, core::Real > my_scores_out; // empty
		std::istringstream inss( outss.str() );
		TS_ASSERT( core::import_pose::atom_tree_diffs::header_from_atom_tree_diff(inss, tag_out, my_scores_out) );
		TS_ASSERT( tag_out == "tag" );
		for ( std::map< std::string, core::Real >::iterator pair = my_scores_out.begin(), pair_end = my_scores_out.end(); pair != pair_end; ++pair ) {
			TS_ASSERT_DELTA( my_scores[ pair->first ], pair->second, (1e-4)*std::abs(pair->second) );
		}

		TS_ASSERT( core::import_pose::atom_tree_diffs::pose_from_atom_tree_diff(inss, start_pose, restored_pose) );
		Real rms_to_modified = scoring::rmsd_no_super(modified_pose, restored_pose, scoring::is_heavyatom);
		TR << "RMS error from save/restore: " << rms_to_modified << std::endl;
		TS_ASSERT( rms_to_modified < 1e-3 );
	}


};
