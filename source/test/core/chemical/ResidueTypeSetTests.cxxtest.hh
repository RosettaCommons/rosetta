// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/chemical/ResidueTypeSetTests.cxxtest.hh
/// @brief unit test for ResidueTypeSet, not really complete yet
/// @author Florian Richter, jan 11,


// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit Headers
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>

// Project Headers
#include <core/chemical/AA.hh>
#include <core/chemical/ChemicalManager.hh>

// Platform Headers
#include <basic/Tracer.hh>
#include <basic/database/open.hh>
#
// ObjexxFCL Headers
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

// C++ Headers
#include <string>
#include <ostream>

using std::endl;
using std::string;

static basic::Tracer TR("core.chemical.ResidueTypeSetTests.cxxtest");

class ResidueTypeSetTests : public CxxTest::TestSuite {

public:
	core::Real delta_percent;

	void setUp() {
		core_init();
		delta_percent=0.0001;
	}

	void tearDown() {}

	void test_residue_type_sets() {
		using namespace core::chemical;
		using namespace ObjexxFCL::format;

		int width = 15;
		TR << A(width,"ResidueTypeSet") << A(width,"NumResTypes") << endl;

		string rss;
		ResidueTypeSetOP rs;
		rss = FA_STANDARD;
		rs = ChemicalManager::get_instance()->nonconst_residue_type_set(rss );
		TR << A(width, rss) << I(width,rs->residue_types().size()) << endl;

		//rss = CENTROID;
		//rs = ChemicalManager::get_instance()->residue_type_set(rss );
		//TR << A(width, rss) << I(width,rs->residue_types().size()) << endl;

		ResidueType const & serine = rs->name_map( "SER" );
		TS_ASSERT_DELTA(serine.molecular_mass(), 87, delta_percent);

		ResidueTypeOP modser = serine.clone();
		modser->nbr_radius( 15.0);
		modser->name( "bigser" );

		//get some stuff from the residue type set
		core::Size n_res_types = rs->residue_types().size();
		core::Size n_ser_types = rs->name3_map( "SER" ).size();
		core::Size n_gln_types = rs->name3_map( "GLN" ).size();
		core::Size n_ser_aa = rs->aa_map( aa_ser ).size();
		//ResidueTypeCOP pointer10 = rs->residue_types()[10];

		//now change the residue type set
		rs->add_residue_type( modser );

		//now make sure everything is as should be
		TS_ASSERT( n_res_types + 1 == rs->residue_types().size());
		TS_ASSERT( n_ser_types + 1 == rs->name3_map( "SER" ).size()  );
		TS_ASSERT( n_gln_types == rs->name3_map( "GLN" ).size() );
		TS_ASSERT( n_ser_aa + 1 == rs->aa_map( aa_ser ).size() );
		TS_ASSERT( rs->has_name("bigser") );

		//ResidueTypeCOP * newpointer10 = &(rs->residue_types()[10]);
		//TR << "old pointer addr is " << pointer10 << ", new pointer addr is " << newpointer10 << std::endl;
		//TS_ASSERT( pointer10 == newpointer10 );
	}


	void test_extra_params() {
		using namespace core::chemical;
		using namespace ObjexxFCL::format;

		std::vector< std::string > extra_params_files( 1, "core/chemical/1pqc.params");
		std::vector< std::string > extra_patch_files( 1, "core/chemical/1pqc_test.patch");

		ResidueTypeSet rtset( FA_STANDARD, basic::database::full_name( "chemical/residue_type_sets/"+FA_STANDARD+"/" ), extra_params_files, extra_patch_files );

		TS_ASSERT( rtset.has_name3("QC1") );
		TS_ASSERT( rtset.has_name("QC1") );
		TS_ASSERT( rtset.has_name("QC1_p:1pqcTestPatch") );

		ResidueType const & plain( rtset.name_map("QC1") );
		ResidueType const & decorated( rtset.name_map("QC1_p:1pqcTestPatch") );

		// This is here mainly to make sure it doesn't crash.
		ResidueType const & dec_gen( rtset.get_residue_type_with_variant_added( plain, "TESTPATCH" ) );

		TS_ASSERT_EQUALS( plain.natoms(), 43 );
		TS_ASSERT_EQUALS( decorated.natoms(), 48 );
		TS_ASSERT_EQUALS( dec_gen.natoms(), 48 );

		TS_ASSERT( plain.has("F9") );
		TS_ASSERT( plain.has("S1") );
		TS_ASSERT( plain.has("H12") );
		TS_ASSERT( ! plain.has("OP1") );
		TS_ASSERT( ! plain.has("3HP2") );

		TS_ASSERT( decorated.has("F9") );
		TS_ASSERT( dec_gen.has("S1") );
		TS_ASSERT( ! decorated.has("H12") );
		TS_ASSERT( dec_gen.has("OP1") );
		TS_ASSERT( decorated.has("3HP2") );

	}

	void test_delete_residue() {
		using namespace core::chemical;

		std::vector< std::string > extra_params_files( 1, "core/chemical/1pqc.params");
		std::vector< std::string > extra_patch_files;
		ResidueTypeSet rtset( FA_STANDARD, basic::database::full_name( "chemical/residue_type_sets/"+FA_STANDARD+"/" ), extra_params_files, extra_patch_files );

		TS_ASSERT(rtset.has_name("QC1"));
		TS_ASSERT(rtset.has_name3("QC1"));

		rtset.remove_residue_type("QC1");

		TS_ASSERT(!rtset.has_name("QC1"));
		TS_ASSERT(!rtset.has_name3("QC1"));
	}


};
