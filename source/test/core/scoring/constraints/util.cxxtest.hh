// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/core/scoring/constraints/util.cxxtest.hh
/// @brief  test suite for constraints utilities
/// @author James Thompson (math)
/// @author Steven Lewis smlewi@gmail.com (cst_file utilities)

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>
#include <test/UTracer.hh>
#include <util/pose_funcs.hh>

#include <core/types.hh>

#include <core/scoring/constraints/util.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>

#include <core/scoring/ScoreFunction.hh>

#include <core/pose/Pose.hh>

#include <core/conformation/Residue.hh>

#include <core/id/AtomID.hh>

#include <utility/io/ozstream.hh>

#include <utility/vector1.hh>


using namespace core;

class ConstraintUtilTests : public CxxTest::TestSuite {
public:
	ConstraintUtilTests() {}

	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	//IO UTILITIES TESTS
	//IO test recipe:
	//write a "constraints" string to a file
	//read the constraints in
	//use ConstraintSet show to match strings against what constraints ought to be

	void test_get_cst_file_option() {
		std::string const cst_file("cst_file");

		core_init_with_additional_options( "-constraints:cst_file " + cst_file );

		TS_ASSERT_EQUALS(core::scoring::constraints::get_cst_file_option(), cst_file);
		//std::cout << "cst file is " << cst_file << " " << core::scoring::constraints::get_cst_file_option() << std::endl;

		return;
	}

	void test_get_cst_fa_file_option() {
		std::string const cst_fa_file("cst_fa_file");

		core_init_with_additional_options( "-constraints:cst_fa_file " + cst_fa_file );

		TS_ASSERT_EQUALS(core::scoring::constraints::get_cst_fa_file_option(), cst_fa_file);
		//std::cout << "cst file is " << cst_fa_file << " " << core::scoring::constraints::get_cst_fa_file_option() << std::endl;

		return;
	}

	void test_add_constraints_from_cmdline_to_scorefxn(){

		//prepare constants and option system
		core::Real const weight_val(12.345);
		std::ostringstream convert;
		convert << weight_val;
		std::string const weight_str(convert.str());
		//std::cout << weight_str << std::endl;

		core_init_with_additional_options( "-constraints:cst_weight " + weight_str );

		//apply to scorefunction
		core::scoring::ScoreFunction scorefunction;

		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::atom_pair_constraint ), 0, 1e-12);
		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::angle_constraint     ), 0, 1e-12);
		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::dihedral_constraint  ), 0, 1e-12);
		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::coordinate_constraint), 0, 1e-12);

		TS_ASSERT(scorefunction.has_zero_weight(core::scoring::atom_pair_constraint ));
		TS_ASSERT(scorefunction.has_zero_weight(core::scoring::angle_constraint     ));
		TS_ASSERT(scorefunction.has_zero_weight(core::scoring::dihedral_constraint  ));
		TS_ASSERT(scorefunction.has_zero_weight(core::scoring::coordinate_constraint));

		core::scoring::constraints::add_constraints_from_cmdline_to_scorefxn(scorefunction);

		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::atom_pair_constraint ), weight_val, 1e-12);
		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::angle_constraint     ), weight_val, 1e-12);
		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::dihedral_constraint  ), weight_val, 1e-12);
		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::coordinate_constraint), weight_val, 1e-12);

		TS_ASSERT(scorefunction.has_nonzero_weight(core::scoring::atom_pair_constraint ));
		TS_ASSERT(scorefunction.has_nonzero_weight(core::scoring::angle_constraint     ));
		TS_ASSERT(scorefunction.has_nonzero_weight(core::scoring::dihedral_constraint  ));
		TS_ASSERT(scorefunction.has_nonzero_weight(core::scoring::coordinate_constraint));

		//Ideally we would confirm no OTHER scorefunction terms are changed...
		return;
	}

	void test_add_fa_constraints_from_cmdline_to_scorefxn(){

		//prepare constants and option system
		core::Real const weight_val(4.5678);
		std::ostringstream convert;
		convert << weight_val;
		std::string const weight_str(convert.str());
		//std::cout << weight_str << std::endl;

		core_init_with_additional_options( "-constraints:cst_fa_weight " + weight_str );

		//apply to scorefunction
		core::scoring::ScoreFunction scorefunction;

		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::atom_pair_constraint ), 0, 1e-12);
		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::angle_constraint     ), 0, 1e-12);
		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::dihedral_constraint  ), 0, 1e-12);
		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::coordinate_constraint), 0, 1e-12);

		TS_ASSERT(scorefunction.has_zero_weight(core::scoring::atom_pair_constraint ));
		TS_ASSERT(scorefunction.has_zero_weight(core::scoring::angle_constraint     ));
		TS_ASSERT(scorefunction.has_zero_weight(core::scoring::dihedral_constraint  ));
		TS_ASSERT(scorefunction.has_zero_weight(core::scoring::coordinate_constraint));

		core::scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn(scorefunction);

		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::atom_pair_constraint ), weight_val, 1e-12);
		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::angle_constraint     ), weight_val, 1e-12);
		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::dihedral_constraint  ), weight_val, 1e-12);
		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::coordinate_constraint), weight_val, 1e-12);

		TS_ASSERT(scorefunction.has_nonzero_weight(core::scoring::atom_pair_constraint ));
		TS_ASSERT(scorefunction.has_nonzero_weight(core::scoring::angle_constraint     ));
		TS_ASSERT(scorefunction.has_nonzero_weight(core::scoring::dihedral_constraint  ));
		TS_ASSERT(scorefunction.has_nonzero_weight(core::scoring::coordinate_constraint));

		//Ideally we would confirm no OTHER scorefunction terms are changed...
		return;
	}

	void test_add_constraints_from_cmdline_to_pose() {

		//create pose and constraints file
		std::string const cst_file("cst_file");
		core_init_with_additional_options( "-constraints:cst_file " + cst_file );

		//Pose is:
		//seq chain num
		//L A 1
		//D A 2
		//N B 3
		//L B 4
		core::pose::Pose pose(create_pdb_string_2res_1ten_2res_trp_cage_pose());

		//write constraints file (locally)
		//yes, this is a strange multiline string literal
		std::string const cst_file_string(
			"AtomPair  CA     1  CA     3 HARMONIC 10 0.2\n"
			"AtomPair  CA     2  CA     4 HARMONIC 10 0.2\n");
		//formatted to match show_definition output from ConstraintSet

		//std::cout << cst_file_string << std::endl;

		utility::io::ozstream out( cst_file );
		out << cst_file_string;
		out.close();

		//check that constraint set is empty before we start
		TS_ASSERT(pose.constraint_set()->is_empty());
		TS_ASSERT(!pose.constraint_set()->residue_pair_constraint_exists(1,3));
		TS_ASSERT(!pose.constraint_set()->residue_pair_constraint_exists(2,4));

		//add constraints...
		core::scoring::constraints::add_constraints_from_cmdline_to_pose(pose);

		//simple tests
		TS_ASSERT(!pose.constraint_set()->is_empty());
		TS_ASSERT(pose.constraint_set()->residue_pair_constraint_exists(1,3));
		TS_ASSERT(pose.constraint_set()->residue_pair_constraint_exists(2,4));

		// std::cout << "show" << std::endl;
		// pose.constraint_set()->show(std::cout);
		// std::cout << "show_definition" << std::endl;
		// pose.constraint_set()->show_definition(std::cout, pose);
		// std::cout << "show_numbers" << std::endl;
		// pose.constraint_set()->show_numbers(std::cout);
		// std::cout << "end" << std::endl;

		//complex show string style test
		std::ostringstream show_capture;
		pose.constraint_set()->show_definition(show_capture, pose);

		TS_ASSERT_EQUALS(show_capture.str(), cst_file_string);

		return;
	}

	void test_add_fa_constraints_from_cmdline_to_pose() {

		//create pose and constraints file
		std::string const cst_fa_file("cst_fa_file");
		core_init_with_additional_options( "-constraints:cst_fa_file " + cst_fa_file );

		//Pose is:
		//seq chain num
		//L A 1
		//D A 2
		//N B 3
		//L B 4
		core::pose::Pose pose(create_pdb_string_2res_1ten_2res_trp_cage_pose());

		//write constraints file (locally)
		//yes, this is a strange multiline string literal
		std::string const cst_fa_file_string(
			"AtomPair  CA     1  CA     3 HARMONIC 11 0.2\n"
			"AtomPair  CA     2  CA     4 HARMONIC 11 0.2\n");
		//formatted to match show_definition output from ConstraintSet

		//std::cout << cst_fa_file_string << std::endl;

		utility::io::ozstream out( cst_fa_file );
		out << cst_fa_file_string;
		out.close();

		//check that constraint set is empty before we start
		TS_ASSERT(pose.constraint_set()->is_empty());
		TS_ASSERT(!pose.constraint_set()->residue_pair_constraint_exists(1,3));
		TS_ASSERT(!pose.constraint_set()->residue_pair_constraint_exists(2,4));

		//add constraints...
		core::scoring::constraints::add_fa_constraints_from_cmdline_to_pose(pose);

		//simple tests
		TS_ASSERT(!pose.constraint_set()->is_empty());
		TS_ASSERT(pose.constraint_set()->residue_pair_constraint_exists(1,3));
		TS_ASSERT(pose.constraint_set()->residue_pair_constraint_exists(2,4));

		// std::cout << "show" << std::endl;
		// pose.constraint_set()->show(std::cout);
		// std::cout << "show_definition" << std::endl;
		// pose.constraint_set()->show_definition(std::cout, pose);
		// std::cout << "show_numbers" << std::endl;
		// pose.constraint_set()->show_numbers(std::cout);
		// std::cout << "end" << std::endl;

		//complex show string style test
		std::ostringstream show_capture;
		pose.constraint_set()->show_definition(show_capture, pose);

		TS_ASSERT_EQUALS(show_capture.str(), cst_fa_file_string);

		return;
	}

	//this just runs both test_add_constraints_from_cmdline_to_pose and _to_scorefunction
	//it's naked code duplication, which is appropriate for unit testing
	void test_add_constraints_from_cmdline(){

		//prepare constants and option system
		core::Real const weight_val(45.678);
		std::ostringstream convert;
		convert << weight_val;
		std::string const weight_str(convert.str());
		//std::cout << weight_str << std::endl;

		std::string const cst_file("cst_file");

		core_init_with_additional_options( "-constraints:cst_file " + cst_file + " -constraints:cst_weight " + weight_str );

		//create Pose, scorefunction, cst file
		//Pose is:
		//seq chain num
		//L A 1
		//D A 2
		//N B 3
		//L B 4
		core::pose::Pose pose(create_pdb_string_2res_1ten_2res_trp_cage_pose());

		//write constraints file (locally)
		//yes, this is a strange multiline string literal
		std::string const cst_file_string(
			"AtomPair  CA     1  CA     3 HARMONIC 14 0.2\n"
			"AtomPair  CA     2  CA     4 HARMONIC 14 0.2\n");
		//formatted to match show_definition output from ConstraintSet

		//std::cout << cst_file_string << std::endl;

		utility::io::ozstream out( cst_file );
		out << cst_file_string;
		out.close();

		core::scoring::ScoreFunction scorefunction;


		//pre-tests

		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::atom_pair_constraint ), 0, 1e-12);
		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::angle_constraint     ), 0, 1e-12);
		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::dihedral_constraint  ), 0, 1e-12);
		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::coordinate_constraint), 0, 1e-12);

		TS_ASSERT(scorefunction.has_zero_weight(core::scoring::atom_pair_constraint ));
		TS_ASSERT(scorefunction.has_zero_weight(core::scoring::angle_constraint     ));
		TS_ASSERT(scorefunction.has_zero_weight(core::scoring::dihedral_constraint  ));
		TS_ASSERT(scorefunction.has_zero_weight(core::scoring::coordinate_constraint));

		TS_ASSERT(pose.constraint_set()->is_empty());
		TS_ASSERT(!pose.constraint_set()->residue_pair_constraint_exists(1,3));
		TS_ASSERT(!pose.constraint_set()->residue_pair_constraint_exists(2,4));

		core::scoring::constraints::add_constraints_from_cmdline(pose, scorefunction);

		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::atom_pair_constraint ), weight_val, 1e-12);
		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::angle_constraint     ), weight_val, 1e-12);
		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::dihedral_constraint  ), weight_val, 1e-12);
		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::coordinate_constraint), weight_val, 1e-12);

		TS_ASSERT(scorefunction.has_nonzero_weight(core::scoring::atom_pair_constraint ));
		TS_ASSERT(scorefunction.has_nonzero_weight(core::scoring::angle_constraint     ));
		TS_ASSERT(scorefunction.has_nonzero_weight(core::scoring::dihedral_constraint  ));
		TS_ASSERT(scorefunction.has_nonzero_weight(core::scoring::coordinate_constraint));

		//Ideally we would confirm no OTHER scorefunction terms are changed...

		TS_ASSERT(!pose.constraint_set()->is_empty());
		TS_ASSERT(pose.constraint_set()->residue_pair_constraint_exists(1,3));
		TS_ASSERT(pose.constraint_set()->residue_pair_constraint_exists(2,4));

		// std::cout << "show" << std::endl;
		// pose.constraint_set()->show(std::cout);
		// std::cout << "show_definition" << std::endl;
		// pose.constraint_set()->show_definition(std::cout, pose);
		// std::cout << "show_numbers" << std::endl;
		// pose.constraint_set()->show_numbers(std::cout);
		// std::cout << "end" << std::endl;

		//complex show string style test
		std::ostringstream show_capture;
		pose.constraint_set()->show_definition(show_capture, pose);

		TS_ASSERT_EQUALS(show_capture.str(), cst_file_string);

		return;
	}

	//this just runs both test_add_fa_constraints_from_cmdline_to_pose and _to_scorefunction
	//it's naked code duplication, which is appropriate for unit testing
	void test_add_fa_constraints_from_cmdline(){

		//prepare constants and option system
		core::Real const weight_val(90.123);
		std::ostringstream convert;
		convert << weight_val;
		std::string const weight_str(convert.str());
		//std::cout << weight_str << std::endl;

		std::string const cst_fa_file("cst_fa_file");

		core_init_with_additional_options( "-constraints:cst_fa_file " + cst_fa_file + " -constraints:cst_fa_weight " + weight_str );

		//create Pose, scorefunction, cst file
		//Pose is:
		//seq chain num
		//L A 1
		//D A 2
		//N B 3
		//L B 4
		core::pose::Pose pose(create_pdb_string_2res_1ten_2res_trp_cage_pose());

		//write constraints file (locally)
		//yes, this is a strange multiline string literal
		std::string const cst_fa_file_string(
			"AtomPair  CA     1  CA     3 HARMONIC 10 0.25\n"
			"AtomPair  CA     2  CA     4 HARMONIC 10 0.25\n");
		//formatted to match show_definition output from ConstraintSet

		//std::cout << cst_fa_file_string << std::endl;

		utility::io::ozstream out( cst_fa_file );
		out << cst_fa_file_string;
		out.close();

		core::scoring::ScoreFunction scorefunction;


		//pre-tests

		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::atom_pair_constraint ), 0, 1e-12);
		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::angle_constraint     ), 0, 1e-12);
		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::dihedral_constraint  ), 0, 1e-12);
		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::coordinate_constraint), 0, 1e-12);

		TS_ASSERT(scorefunction.has_zero_weight(core::scoring::atom_pair_constraint ));
		TS_ASSERT(scorefunction.has_zero_weight(core::scoring::angle_constraint     ));
		TS_ASSERT(scorefunction.has_zero_weight(core::scoring::dihedral_constraint  ));
		TS_ASSERT(scorefunction.has_zero_weight(core::scoring::coordinate_constraint));

		TS_ASSERT(pose.constraint_set()->is_empty());
		TS_ASSERT(!pose.constraint_set()->residue_pair_constraint_exists(1,3));
		TS_ASSERT(!pose.constraint_set()->residue_pair_constraint_exists(2,4));

		core::scoring::constraints::add_fa_constraints_from_cmdline(pose, scorefunction);

		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::atom_pair_constraint ), weight_val, 1e-12);
		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::angle_constraint     ), weight_val, 1e-12);
		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::dihedral_constraint  ), weight_val, 1e-12);
		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::coordinate_constraint), weight_val, 1e-12);

		TS_ASSERT(scorefunction.has_nonzero_weight(core::scoring::atom_pair_constraint ));
		TS_ASSERT(scorefunction.has_nonzero_weight(core::scoring::angle_constraint     ));
		TS_ASSERT(scorefunction.has_nonzero_weight(core::scoring::dihedral_constraint  ));
		TS_ASSERT(scorefunction.has_nonzero_weight(core::scoring::coordinate_constraint));

		//Ideally we would confirm no OTHER scorefunction terms are changed...

		TS_ASSERT(!pose.constraint_set()->is_empty());
		TS_ASSERT(pose.constraint_set()->residue_pair_constraint_exists(1,3));
		TS_ASSERT(pose.constraint_set()->residue_pair_constraint_exists(2,4));

		// std::cout << "show" << std::endl;
		// pose.constraint_set()->show(std::cout);
		// std::cout << "show_definition" << std::endl;
		// pose.constraint_set()->show_definition(std::cout, pose);
		// std::cout << "show_numbers" << std::endl;
		// pose.constraint_set()->show_numbers(std::cout);
		// std::cout << "end" << std::endl;

		//complex show string style test
		std::ostringstream show_capture;
		pose.constraint_set()->show_definition(show_capture, pose);

		TS_ASSERT_EQUALS(show_capture.str(), cst_fa_file_string);

		return;
	}

	//to test merging, add one constraint to the pose directly, test that only it exists, add one more from commandline, and test that both coexist neatly.
	void test_merge_fa_constraints_from_cmdline_to_pose(){

		//create pose with one constraint
		//Pose is:
		//seq chain num
		//L A 1
		//D A 2
		//N B 3
		//L B 4
		core::pose::Pose pose(create_pdb_string_2res_1ten_2res_trp_cage_pose());

		//"AtomPair  CA     1  CA     3 HARMONIC 10 0.29\n"
		core::scoring::func::FuncOP fx( new core::scoring::func::HarmonicFunc( 10.0, 0.29 ) );
		core::scoring::constraints::ConstraintCOP cst13( new core::scoring::constraints::AtomPairConstraint(
				core::id::AtomID( pose.residue( 1 ).atom_index( "CA" ), 1 ),
				core::id::AtomID( pose.residue( 3 ).atom_index( "CA" ), 3 ),
				fx ) );

		pose.add_constraint(cst13);

		//create constraints file
		std::string const cst_fa_file("cst_fa_file");
		core_init_with_additional_options( "-constraints:cst_fa_file " + cst_fa_file );

		//write constraints file (locally)
		std::string const cst_fa_file_string(
			//"AtomPair  CA     1  CA     3 HARMONIC 10 0.29\n"
			"AtomPair  CA     2  CA     4 HARMONIC 10 0.27\n");
		//formatted to match show_definition output from ConstraintSet

		std::string const cst_notfile_string(
			"AtomPair  CA     1  CA     3 HARMONIC 10 0.29\n");
		//"AtomPair  CA     2  CA     4 HARMONIC 10 0.27\n");

		std::string const cst_bothfile_string(cst_notfile_string + cst_fa_file_string);

		//std::cout << "cst_fa_file_string \n" << cst_fa_file_string << std::endl;
		//std::cout << "cst_notfile_string \n" << cst_notfile_string << std::endl;
		//std::cout << "cst_bothfile_string \n" << cst_bothfile_string << std::endl;

		utility::io::ozstream out( cst_fa_file );
		out << cst_fa_file_string;
		out.close();

		//check that constraint set is correct before we start
		TS_ASSERT(!pose.constraint_set()->is_empty());
		TS_ASSERT(pose.constraint_set()->residue_pair_constraint_exists(1,3));
		TS_ASSERT(!pose.constraint_set()->residue_pair_constraint_exists(2,4));

		//complex show string style test
		std::ostringstream show_capture_pre;
		pose.constraint_set()->show_definition(show_capture_pre, pose);

		TS_ASSERT_EQUALS(show_capture_pre.str(), cst_notfile_string);

		//add constraints...
		core::scoring::constraints::merge_fa_constraints_from_cmdline_to_pose(pose);

		//simple tests
		TS_ASSERT(!pose.constraint_set()->is_empty());
		TS_ASSERT(pose.constraint_set()->residue_pair_constraint_exists(1,3));
		TS_ASSERT(pose.constraint_set()->residue_pair_constraint_exists(2,4));

		// std::cout << "show" << std::endl;
		// pose.constraint_set()->show(std::cout);
		// std::cout << "show_definition" << std::endl;
		// pose.constraint_set()->show_definition(std::cout, pose);
		// std::cout << "show_numbers" << std::endl;
		// pose.constraint_set()->show_numbers(std::cout);
		// std::cout << "end" << std::endl;

		//complex show string style test
		std::ostringstream show_capture_post;
		pose.constraint_set()->show_definition(show_capture_post, pose);

		TS_ASSERT_EQUALS(show_capture_post.str(), cst_bothfile_string);

		return;
	}

	//merge test: some weights pre-exist, only set the non-pre-set ones
	void test_merge_fa_constraints_from_cmdline_to_scorefxn(){

		//prepare constants and option system
		core::Real const weight_val(6.78);
		std::ostringstream convert;
		convert << weight_val;
		std::string const weight_str(convert.str());
		//std::cout << weight_str << std::endl;

		core_init_with_additional_options( "-constraints:cst_fa_weight " + weight_str );

		//set up scorefunction
		core::scoring::ScoreFunction scorefunction;
		core::Real const preset_weight(20.13);
		scorefunction.set_weight(core::scoring::coordinate_constraint, preset_weight);

		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::atom_pair_constraint ), 0, 1e-12);
		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::angle_constraint     ), 0, 1e-12);
		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::dihedral_constraint  ), 0, 1e-12);
		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::coordinate_constraint), preset_weight, 1e-12);

		TS_ASSERT(scorefunction.has_zero_weight(core::scoring::atom_pair_constraint ));
		TS_ASSERT(scorefunction.has_zero_weight(core::scoring::angle_constraint     ));
		TS_ASSERT(scorefunction.has_zero_weight(core::scoring::dihedral_constraint  ));
		TS_ASSERT(scorefunction.has_nonzero_weight(core::scoring::coordinate_constraint));

		core::scoring::constraints::merge_fa_constraints_from_cmdline_to_scorefxn(scorefunction);
		//this should be noisy too...but even if we can TS_ASSERT makes noise, -mute all will still fail it

		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::atom_pair_constraint ), weight_val, 1e-12);
		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::angle_constraint     ), weight_val, 1e-12);
		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::dihedral_constraint  ), weight_val, 1e-12);
		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::coordinate_constraint), preset_weight, 1e-12);

		TS_ASSERT(scorefunction.has_nonzero_weight(core::scoring::atom_pair_constraint ));
		TS_ASSERT(scorefunction.has_nonzero_weight(core::scoring::angle_constraint     ));
		TS_ASSERT(scorefunction.has_nonzero_weight(core::scoring::dihedral_constraint  ));
		TS_ASSERT(scorefunction.has_nonzero_weight(core::scoring::coordinate_constraint));

		//Ideally we would confirm no OTHER scorefunction terms are changed...
		return;
	}

	void test_merge_fa_constraints_from_cmdline(){

		//create pose with one constraint
		//Pose is:
		//seq chain num
		//L A 1
		//D A 2
		//N B 3
		//L B 4
		core::pose::Pose pose(create_pdb_string_2res_1ten_2res_trp_cage_pose());

		//"AtomPair  CA     1  CA     3 HARMONIC 14 0.2\n"
		core::scoring::func::FuncOP fx( new core::scoring::func::HarmonicFunc( 14.0, 0.2 ) );
		core::scoring::constraints::ConstraintCOP cst13( new core::scoring::constraints::AtomPairConstraint(
				core::id::AtomID( pose.residue( 1 ).atom_index( "CA" ), 1 ),
				core::id::AtomID( pose.residue( 3 ).atom_index( "CA" ), 3 ),
				fx ) );

		pose.add_constraint(cst13);

		//create constraints file
		std::string const cst_fa_file("cst_fa_file");

		//write constraints file (locally)
		std::string const cst_fa_file_string(
			//"AtomPair  CA     1  CA     3 HARMONIC 14 0.2\n"
			"AtomPair  CA     2  CA     4 HARMONIC 12 0.2\n");
		//formatted to match show_definition output from ConstraintSet

		std::string const cst_notfile_string(
			"AtomPair  CA     1  CA     3 HARMONIC 14 0.2\n");
		//"AtomPair  CA     2  CA     4 HARMONIC 12 0.2\n");

		std::string const cst_bothfile_string(cst_notfile_string + cst_fa_file_string);

		//std::cout << "cst_fa_file_string \n" << cst_fa_file_string << std::endl;
		//std::cout << "cst_notfile_string \n" << cst_notfile_string << std::endl;
		//std::cout << "cst_bothfile_string \n" << cst_bothfile_string << std::endl;

		utility::io::ozstream out( cst_fa_file );
		out << cst_fa_file_string;
		out.close();

		//set up scorefunction
		core::Real const weight_val(4567.8);
		std::ostringstream convert;
		convert << weight_val;
		std::string const weight_str(convert.str());
		//std::cout << weight_str << std::endl;

		core::scoring::ScoreFunction scorefunction;
		core::Real const preset_weight(20.15);
		scorefunction.set_weight(core::scoring::coordinate_constraint, preset_weight);

		//load this into options
		core_init_with_additional_options( "-constraints:cst_fa_file " + cst_fa_file + " -constraints:cst_fa_weight " + weight_str );

		//check that constraint set is correct before we start
		TS_ASSERT(!pose.constraint_set()->is_empty());
		TS_ASSERT(pose.constraint_set()->residue_pair_constraint_exists(1,3));
		TS_ASSERT(!pose.constraint_set()->residue_pair_constraint_exists(2,4));

		//complex show string style test
		std::ostringstream show_capture_pre;
		pose.constraint_set()->show_definition(show_capture_pre, pose);
		TS_ASSERT_EQUALS(show_capture_pre.str(), cst_notfile_string);

		//check that weights are correct
		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::atom_pair_constraint ), 0, 1e-12);
		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::angle_constraint     ), 0, 1e-12);
		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::dihedral_constraint  ), 0, 1e-12);
		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::coordinate_constraint), preset_weight, 1e-12);

		TS_ASSERT(scorefunction.has_zero_weight(core::scoring::atom_pair_constraint ));
		TS_ASSERT(scorefunction.has_zero_weight(core::scoring::angle_constraint     ));
		TS_ASSERT(scorefunction.has_zero_weight(core::scoring::dihedral_constraint  ));
		TS_ASSERT(scorefunction.has_nonzero_weight(core::scoring::coordinate_constraint));

		//add constraints...
		core::scoring::constraints::merge_fa_constraints_from_cmdline(pose, scorefunction);
		//this should be noisy too...but even if we can TS_ASSERT makes noise, -mute all will still fail it

		//simple tests
		TS_ASSERT(!pose.constraint_set()->is_empty());
		TS_ASSERT(pose.constraint_set()->residue_pair_constraint_exists(1,3));
		TS_ASSERT(pose.constraint_set()->residue_pair_constraint_exists(2,4));

		// std::cout << "show" << std::endl;
		// pose.constraint_set()->show(std::cout);
		// std::cout << "show_definition" << std::endl;
		// pose.constraint_set()->show_definition(std::cout, pose);
		// std::cout << "show_numbers" << std::endl;
		// pose.constraint_set()->show_numbers(std::cout);
		// std::cout << "end" << std::endl;

		//complex show string style test
		std::ostringstream show_capture_post;
		pose.constraint_set()->show_definition(show_capture_post, pose);

		TS_ASSERT_EQUALS(show_capture_post.str(), cst_bothfile_string);

		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::atom_pair_constraint ), weight_val, 1e-12);
		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::angle_constraint     ), weight_val, 1e-12);
		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::dihedral_constraint  ), weight_val, 1e-12);
		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::coordinate_constraint), preset_weight, 1e-12);

		TS_ASSERT(scorefunction.has_nonzero_weight(core::scoring::atom_pair_constraint ));
		TS_ASSERT(scorefunction.has_nonzero_weight(core::scoring::angle_constraint     ));
		TS_ASSERT(scorefunction.has_nonzero_weight(core::scoring::dihedral_constraint  ));
		TS_ASSERT(scorefunction.has_nonzero_weight(core::scoring::coordinate_constraint));

		//Ideally we would confirm no OTHER scorefunction terms are changed...
		return;
	}

	//to test merging, add one constraint to the pose directly, test that only it exists, add one more from commandline, and test that both coexist neatly.
	void test_merge_constraints_from_cmdline_to_pose(){

		//create pose with one constraint
		//Pose is:
		//seq chain num
		//L A 1
		//D A 2
		//N B 3
		//L B 4
		core::pose::Pose pose(create_pdb_string_2res_1ten_2res_trp_cage_pose());

		//"AtomPair  CA     1  CA     3 HARMONIC 10 0.29\n"
		core::scoring::func::FuncOP fx( new core::scoring::func::HarmonicFunc( 10.5, 0.29 ) );
		core::scoring::constraints::ConstraintCOP cst13( new core::scoring::constraints::AtomPairConstraint(
				core::id::AtomID( pose.residue( 1 ).atom_index( "CA" ), 1 ),
				core::id::AtomID( pose.residue( 3 ).atom_index( "CA" ), 3 ),
				fx ) );

		pose.add_constraint(cst13);

		//create constraints file
		std::string const cst_file("cst_file");
		core_init_with_additional_options( "-constraints:cst_file " + cst_file );

		//write constraints file (locally)
		std::string const cst_file_string(
			//"AtomPair  CA     1  CA     3 HARMONIC 10.5 0.29\n"
			"AtomPair  CA     2  CA     4 HARMONIC 17 0.27\n");
		//formatted to match show_definition output from ConstraintSet

		std::string const cst_notfile_string(
			"AtomPair  CA     1  CA     3 HARMONIC 10.5 0.29\n");
		//"AtomPair  CA     2  CA     4 HARMONIC 17 0.27\n");

		std::string const cst_bothfile_string(cst_notfile_string + cst_file_string);

		//std::cout << "cst_file_string \n" << cst_file_string << std::endl;
		//std::cout << "cst_notfile_string \n" << cst_notfile_string << std::endl;
		//std::cout << "cst_bothfile_string \n" << cst_bothfile_string << std::endl;

		utility::io::ozstream out( cst_file );
		out << cst_file_string;
		out.close();

		//check that constraint set is correct before we start
		TS_ASSERT(!pose.constraint_set()->is_empty());
		TS_ASSERT(pose.constraint_set()->residue_pair_constraint_exists(1,3));
		TS_ASSERT(!pose.constraint_set()->residue_pair_constraint_exists(2,4));

		//complex show string style test
		std::ostringstream show_capture_pre;
		pose.constraint_set()->show_definition(show_capture_pre, pose);

		TS_ASSERT_EQUALS(show_capture_pre.str(), cst_notfile_string);

		//add constraints...
		core::scoring::constraints::merge_constraints_from_cmdline_to_pose(pose);

		//simple tests
		TS_ASSERT(!pose.constraint_set()->is_empty());
		TS_ASSERT(pose.constraint_set()->residue_pair_constraint_exists(1,3));
		TS_ASSERT(pose.constraint_set()->residue_pair_constraint_exists(2,4));

		// std::cout << "show" << std::endl;
		// pose.constraint_set()->show(std::cout);
		// std::cout << "show_definition" << std::endl;
		// pose.constraint_set()->show_definition(std::cout, pose);
		// std::cout << "show_numbers" << std::endl;
		// pose.constraint_set()->show_numbers(std::cout);
		// std::cout << "end" << std::endl;

		//complex show string style test
		std::ostringstream show_capture_post;
		pose.constraint_set()->show_definition(show_capture_post, pose);

		TS_ASSERT_EQUALS(show_capture_post.str(), cst_bothfile_string);

		return;
	}

	//merge test: some weights pre-exist, only set the non-pre-set ones
	void test_merge_constraints_from_cmdline_to_scorefxn(){

		//prepare constants and option system
		core::Real const weight_val(67.34);
		std::ostringstream convert;
		convert << weight_val;
		std::string const weight_str(convert.str());
		//std::cout << weight_str << std::endl;

		core_init_with_additional_options( "-constraints:cst_weight " + weight_str );

		//set up scorefunction
		core::scoring::ScoreFunction scorefunction;
		core::Real const preset_weight(2.13);
		scorefunction.set_weight(core::scoring::coordinate_constraint, preset_weight);

		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::atom_pair_constraint ), 0, 1e-12);
		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::angle_constraint     ), 0, 1e-12);
		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::dihedral_constraint  ), 0, 1e-12);
		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::coordinate_constraint), preset_weight, 1e-12);

		TS_ASSERT(scorefunction.has_zero_weight(core::scoring::atom_pair_constraint ));
		TS_ASSERT(scorefunction.has_zero_weight(core::scoring::angle_constraint     ));
		TS_ASSERT(scorefunction.has_zero_weight(core::scoring::dihedral_constraint  ));
		TS_ASSERT(scorefunction.has_nonzero_weight(core::scoring::coordinate_constraint));

		core::scoring::constraints::merge_constraints_from_cmdline_to_scorefxn(scorefunction);
		//this should be noisy too...but even if we can TS_ASSERT makes noise, -mute all will still fail it

		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::atom_pair_constraint ), weight_val, 1e-12);
		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::angle_constraint     ), weight_val, 1e-12);
		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::dihedral_constraint  ), weight_val, 1e-12);
		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::coordinate_constraint), preset_weight, 1e-12);

		TS_ASSERT(scorefunction.has_nonzero_weight(core::scoring::atom_pair_constraint ));
		TS_ASSERT(scorefunction.has_nonzero_weight(core::scoring::angle_constraint     ));
		TS_ASSERT(scorefunction.has_nonzero_weight(core::scoring::dihedral_constraint  ));
		TS_ASSERT(scorefunction.has_nonzero_weight(core::scoring::coordinate_constraint));

		//Ideally we would confirm no OTHER scorefunction terms are changed...
		return;
	}

	void test_merge_constraints_from_cmdline(){

		//create pose with one constraint
		//Pose is:
		//seq chain num
		//L A 1
		//D A 2
		//N B 3
		//L B 4
		core::pose::Pose pose(create_pdb_string_2res_1ten_2res_trp_cage_pose());

		//"AtomPair  CA     1  CA     3 HARMONIC 14 0.1\n"
		core::scoring::func::FuncOP fx( new core::scoring::func::HarmonicFunc( 14.0, 0.1 ) );
		core::scoring::constraints::ConstraintCOP cst13( new core::scoring::constraints::AtomPairConstraint(
				core::id::AtomID( pose.residue( 1 ).atom_index( "CA" ), 1 ),
				core::id::AtomID( pose.residue( 3 ).atom_index( "CA" ), 3 ),
				fx ) );

		pose.add_constraint(cst13);

		//create constraints file
		std::string const cst_file("cst_file");

		//write constraints file (locally)
		std::string const cst_file_string(
			//"AtomPair  CA     1  CA     3 HARMONIC 14 0.1\n"
			"AtomPair  CA     2  CA     4 HARMONIC 12 1.2\n");
		//formatted to match show_definition output from ConstraintSet

		std::string const cst_notfile_string(
			"AtomPair  CA     1  CA     3 HARMONIC 14 0.1\n");
		//"AtomPair  CA     2  CA     4 HARMONIC 12 1.2\n");

		std::string const cst_bothfile_string(cst_notfile_string + cst_file_string);

		//std::cout << "cst_file_string \n" << cst_file_string << std::endl;
		//std::cout << "cst_notfile_string \n" << cst_notfile_string << std::endl;
		//std::cout << "cst_bothfile_string \n" << cst_bothfile_string << std::endl;

		utility::io::ozstream out( cst_file );
		out << cst_file_string;
		out.close();

		//set up scorefunction
		core::Real const weight_val(385.5);
		std::ostringstream convert;
		convert << weight_val;
		std::string const weight_str(convert.str());
		//std::cout << weight_str << std::endl;

		core::scoring::ScoreFunction scorefunction;
		core::Real const preset_weight(24.15);
		scorefunction.set_weight(core::scoring::coordinate_constraint, preset_weight);

		//load this into options
		core_init_with_additional_options( "-constraints:cst_file " + cst_file + " -constraints:cst_weight " + weight_str );

		//check that constraint set is correct before we start
		TS_ASSERT(!pose.constraint_set()->is_empty());
		TS_ASSERT(pose.constraint_set()->residue_pair_constraint_exists(1,3));
		TS_ASSERT(!pose.constraint_set()->residue_pair_constraint_exists(2,4));

		//complex show string style test
		std::ostringstream show_capture_pre;
		pose.constraint_set()->show_definition(show_capture_pre, pose);
		TS_ASSERT_EQUALS(show_capture_pre.str(), cst_notfile_string);

		//check that weights are correct
		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::atom_pair_constraint ), 0, 1e-12);
		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::angle_constraint     ), 0, 1e-12);
		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::dihedral_constraint  ), 0, 1e-12);
		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::coordinate_constraint), preset_weight, 1e-12);

		TS_ASSERT(scorefunction.has_zero_weight(core::scoring::atom_pair_constraint ));
		TS_ASSERT(scorefunction.has_zero_weight(core::scoring::angle_constraint     ));
		TS_ASSERT(scorefunction.has_zero_weight(core::scoring::dihedral_constraint  ));
		TS_ASSERT(scorefunction.has_nonzero_weight(core::scoring::coordinate_constraint));

		//add constraints...
		core::scoring::constraints::merge_constraints_from_cmdline(pose, scorefunction);
		//this should be noisy too...but even if we can TS_ASSERT makes noise, -mute all will still fail it

		//simple tests
		TS_ASSERT(!pose.constraint_set()->is_empty());
		TS_ASSERT(pose.constraint_set()->residue_pair_constraint_exists(1,3));
		TS_ASSERT(pose.constraint_set()->residue_pair_constraint_exists(2,4));

		// std::cout << "show" << std::endl;
		// pose.constraint_set()->show(std::cout);
		// std::cout << "show_definition" << std::endl;
		// pose.constraint_set()->show_definition(std::cout, pose);
		// std::cout << "show_numbers" << std::endl;
		// pose.constraint_set()->show_numbers(std::cout);
		// std::cout << "end" << std::endl;

		//complex show string style test
		std::ostringstream show_capture_post;
		pose.constraint_set()->show_definition(show_capture_post, pose);

		TS_ASSERT_EQUALS(show_capture_post.str(), cst_bothfile_string);

		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::atom_pair_constraint ), weight_val, 1e-12);
		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::angle_constraint     ), weight_val, 1e-12);
		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::dihedral_constraint  ), weight_val, 1e-12);
		TS_ASSERT_DELTA(scorefunction.get_weight(core::scoring::coordinate_constraint), preset_weight, 1e-12);

		TS_ASSERT(scorefunction.has_nonzero_weight(core::scoring::atom_pair_constraint ));
		TS_ASSERT(scorefunction.has_nonzero_weight(core::scoring::angle_constraint     ));
		TS_ASSERT(scorefunction.has_nonzero_weight(core::scoring::dihedral_constraint  ));
		TS_ASSERT(scorefunction.has_nonzero_weight(core::scoring::coordinate_constraint));

		//Ideally we would confirm no OTHER scorefunction terms are changed...
		return;
	}

	//MATH FUNCTIONS TESTS

	void test_gaussian_functions() {
		using core::Size;
	 	using core::Real;
		using utility::vector1;
		using namespace core::scoring::constraints;

		Real const TOLERATED_ERROR( 1e-2 );

		Real const mean( 0 );
		Real const sdev( 1 );
		Real const weight( 1 );

		TS_ASSERT_DELTA(
			dgaussian( 0, mean, sdev, weight ), 0.3969, TOLERATED_ERROR
		);
		TS_ASSERT_DELTA(
			logdgaussian( 0, mean, sdev, weight ), -0.9239, TOLERATED_ERROR
		);

	}

	void test_exponential_functions() {
		using core::Size;
	 	using core::Real;
		using utility::vector1;
		using namespace core::scoring::constraints;

		Real const TOLERATED_ERROR( 1e-2 );

		Real const anchor( 0 );
		Real const rate  ( 1 );
		Real const weight( 1 );

		TS_ASSERT_DELTA(
			dexponential( 1, anchor, rate, weight ), 0.3678, TOLERATED_ERROR
		);
	}

}; // ConstraintUtilTests
