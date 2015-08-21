// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/ValidateDunbrackBinaries.cxxtest.hh
/// @brief  test suites to validate Dunbrack binaries for unit tests
/// @author Rocco Moretti

/// @details These are unit tests intended to test if the Dunbrack binaries match the version
/// in the ASCII library. If these unit tests fail, your first course of action should be
/// to re-run them with a clean checkout of master. (DO NOT SKIP THIS STEP) If they still fail,
/// this means that the Dunbrack binaries in your database are corrupt. Delete the *.bin files under
/// $ROSETTA_DATABASE/rotamers (there will likely be three of them, - exact names are printed
/// to the (non-muted) tracer in the course of a failing test).
///
/// If the clean checkout tests fine, then the issue is with your modifications - you somehow changed
/// Dunbrack rotamer reading/writing. If this was intentional, you'll need to increment the
/// current_binary_format_version_id() functions in src/core/pack/dunbrack/RotamerLibrary.cc,
/// and possibly update the operator==() function of SingleResidueDunbrackLibrary and associated
/// classes/subclasses.
///
/// *DO NOT DELETE THE BINARY FILES UNLESS YOU'RE GETTING FAILURES WITH A *CLEAN* CHECKOUT.*
/// - Doing can mask potentially serious changes to the Rosetta codebase.

// Test headers
#include <core/pack/dunbrack/RotamerLibrary.hh>

#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <core/types.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR("core.pack.dunbrack.ValidateDunbrackBinaries.cxxtest");

void print_relevant_info() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	core::pack::dunbrack::RotamerLibrary* rotamer_library(  core::pack::dunbrack::RotamerLibrary::get_instance() );

	TR.Error << "---------------------- Settings: --------------------------------" << std::endl;
	TR.Error << "Database Directory(s): " << std::endl;
	for ( core::Size ii(1); ii <= option[ in::path::database ]().size(); ++ii ) {
		TR.Error << "\t\t" << option[ in::path::database ](ii).name() << std::endl;
	}
	TR.Error << "No binary Dunlib : " << (option[ in::file::no_binary_dunlib ] ? " true " : " false " ) << std::endl;
	TR.Error << "Dun10: " << (option[ corrections::score::dun10 ] ? " true " : " false " ) << std::endl;
	TR.Error << "-correct " << (option[ corrections::correct ] ? " true " : " false " ) << std::endl;
	if ( option[ corrections::score::dun10 ] ) {
		TR.Error << "Dunbrack 2010 directory: " << option[ corrections::score::dun10_dir ].value() << std::endl;
	} else {
		TR.Error << "Dunbrack 2002 file: " << option[ corrections::score::dun02_file ].value() << std::endl;
	}
	TR.Error << "Dunbrack library binary file: " << rotamer_library->get_binary_name() << std::endl;
	TR.Error << "-----------------------------------------------------------------" << std::endl;
}
class ValidateDun10BinariesTests : public CxxTest::TestSuite
{
public:

	// Shared initialization goes here.
	void setUp() {
		// Standard Dun10 binary loading.
		core_init_with_additional_options( "-out:levels core.pack.dunbrack:debug" );
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	void test_dun10_binaries() {

		core::pack::dunbrack::RotamerLibrary* rotamer_library(  core::pack::dunbrack::RotamerLibrary::get_instance() );

		if ( ! rotamer_library->validate_dunbrack_binary() ) {
			TR << "Failure validating the Dunbrack2010 binary" << std::endl;
			print_relevant_info();
			TS_FAIL("Dunbrack2010 ASCII/binary inconsistency");
		}
	}
};

class ValidateDun02BinariesTests : public CxxTest::TestSuite
{
public:

	// Shared initialization goes here.
	void setUp() {
		core_init_with_additional_options( "-out:levels core.pack.dunbrack:debug -restore_pre_talaris_2013_behavior -override_rsd_type_limit" );
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	void test_dun02_binaries() {

		core::pack::dunbrack::RotamerLibrary* rotamer_library(  core::pack::dunbrack::RotamerLibrary::get_instance() );

		if ( ! rotamer_library->validate_dunbrack_binary() ) {
			TR << "Failure validating the Dunbrack2002 binary" << std::endl;
			print_relevant_info();
			TS_FAIL("Dunbrack2002 ASCII/binary inconsistency");
		}
	}
};

class ValidateDun02CorrectBinariesTests : public CxxTest::TestSuite
{
public:

	// Shared initialization goes here.
	void setUp() {
		core_init_with_additional_options( "-out:levels core.pack.dunbrack:debug -restore_pre_talaris_2013_behavior -correct -override_rsd_type_limit" );
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	void test_dun02_correct_binaries() {

		core::pack::dunbrack::RotamerLibrary* rotamer_library(  core::pack::dunbrack::RotamerLibrary::get_instance() );

		if ( ! rotamer_library->validate_dunbrack_binary() ) {
			TR << "Failure validating the Dunbrack2002 -correct binary" << std::endl;
			print_relevant_info();
			TS_FAIL("Dunbrack2002 -correct ASCII/binary inconsistency");
		}
	}
};
