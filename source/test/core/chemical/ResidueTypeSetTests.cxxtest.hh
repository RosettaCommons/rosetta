// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/chemical/ResidueTypeSetTests.cxxtest.hh
/// @brief unit test for ResidueTypeSet, not really complete yet
/// @author Florian Richter, jan 11,


// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit Headers
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeFinder.hh>

// Project Headers
#include <core/chemical/AA.hh>
#include <core/chemical/ChemicalManager.hh>

// Platform Headers
#include <basic/Tracer.hh>
#include <basic/database/open.hh>

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

		string rss;
		ResidueTypeSetOP rs;
		rss = FA_STANDARD;
		rs = ChemicalManager::get_instance()->nonconst_residue_type_set(rss ).get_self_ptr();
		TR << A(width,"ResidueTypeSet") << A(width,"NumBaseResTypes") << endl;
		TR << A(width, rss) << I(width,rs->base_residue_types().size()) << endl;
		TR << A(width,"ResidueTypeSet") << A(width,"NumCustomResTypes") << endl;
		TR << A(width, rss) << I(width,rs->custom_residue_types().size()) << endl;

		ResidueType const & serine = rs->name_map( "SER" );
		TS_ASSERT_DELTA(serine.mass(), 87.0900, delta_percent);

		ResidueTypeOP modser = serine.clone();
		modser->nbr_radius( 15.0);
		modser->name( "bigser" );

		//get some stuff from the residue type set
		core::Size n_base_res_types   = rs->base_residue_types().size();
		core::Size n_custom_res_types = rs->custom_residue_types().size();
		core::Size n_ser_types = ResidueTypeFinder( *rs ).name3( "SER" ).get_all_possible_residue_types().size();
		core::Size n_gln_types = ResidueTypeFinder( *rs ).name3( "GLN" ).get_all_possible_residue_types().size();
		core::Size n_ser_aa = ResidueTypeFinder( *rs ).aa( aa_ser ).get_all_possible_residue_types().size();

		//now change the residue type set
		rs->add_custom_residue_type( modser );

		//now make sure everything is as should be
		TS_ASSERT( n_base_res_types == rs->base_residue_types().size());
		TS_ASSERT( n_custom_res_types + 1 == rs->custom_residue_types().size());
		TS_ASSERT( n_ser_types + 1 == ResidueTypeFinder( *rs ).name3( "SER" ).get_all_possible_residue_types().size() );
		TS_ASSERT( n_gln_types == ResidueTypeFinder( *rs ).name3( "GLN" ).get_all_possible_residue_types().size() );
		TS_ASSERT( n_ser_aa + 1 == ResidueTypeFinder( *rs ).aa( aa_ser ).get_all_possible_residue_types().size() );
		TS_ASSERT( rs->has_name("bigser") );
	}


	void test_extra_params() {
		using namespace core::chemical;
		using namespace ObjexxFCL::format;

		std::vector< std::string > extra_params_files( 1, "core/chemical/1pqc.params");
		std::vector< std::string > extra_patch_files( 1, "core/chemical/1pqc_test.patch");

		ResidueTypeSetOP rtset( new ResidueTypeSet( FA_STANDARD, basic::database::full_name( "chemical/residue_type_sets/"+FA_STANDARD+"/" ) ) );
		rtset->init( extra_params_files, extra_patch_files );

		TS_ASSERT( rtset->has_name3("QC1") );
		TS_ASSERT( rtset->has_name("QC1") );
		TS_ASSERT( rtset->has_name("QC1:1pqcTestPatch") );

		ResidueType const & plain( rtset->name_map("QC1") );
		ResidueType const & decorated( rtset->name_map("QC1:1pqcTestPatch") );

		// This is here mainly to make sure it doesn't crash.
		ResidueType const & dec_gen( rtset->get_residue_type_with_variant_added( plain, SPECIAL_ROT ) );

		TR << "BASE\tFULL" << std::endl;
		TR << plain.base_name() << "\t" << plain.name() << std::endl;
		TR << decorated.base_name() << "\t" << decorated.name() << std::endl;
		TR << std::endl; TR.flush();

		TS_ASSERT_EQUALS( decorated.name(), "QC1:1pqcTestPatch" );
		TS_ASSERT_EQUALS( decorated.base_name(), "QC1" );
		TS_ASSERT_EQUALS( plain.name(), "QC1" );
		TS_ASSERT_EQUALS( plain.base_name(), "QC1" );

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

		ResidueTypeSetOP rtset( new ResidueTypeSet( FA_STANDARD, basic::database::full_name( "chemical/residue_type_sets/"+FA_STANDARD+"/" ) ) );
		rtset->init( extra_params_files, extra_patch_files );

		TS_ASSERT(rtset->has_name("QC1"));
		TS_ASSERT(rtset->has_name3("QC1"));

		// Approved use: we're literally testing if this dangerous function works as expected.
		rtset->remove_base_residue_type_DO_NOT_USE("QC1");

		TS_ASSERT(!rtset->has_name("QC1"));
		TS_ASSERT(!rtset->has_name3("QC1"));
	}

	/// @brief Tests the ability of the ResidueTypeSet to give me the mirror-image ResidueType to a given type.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	void test_get_mirror_type() {
		using namespace core::chemical;

		ResidueTypeSetOP restypeset( new ResidueTypeSet( FA_STANDARD, basic::database::full_name( "chemical/residue_type_sets/"+FA_STANDARD+"/" ) ) );
		std::vector< std::string > extra_params_files;
		std::vector< std::string > extra_patch_files;
		restypeset->init( extra_params_files, extra_patch_files );

		//Test types:
		ResidueTypeCOP lcys_nterm_cterm( restypeset->name_mapOP("CYS:NtermProteinFull:CtermProteinFull") );
		ResidueTypeCOP dile_nacetyl( restypeset->name_mapOP("DILE:AcetylatedNtermProteinFull") );
		ResidueTypeCOP plain_c53( restypeset->name_mapOP("C53") );

		//Return types:
		ResidueTypeCOP mirrored_lcys_nterm_cterm( restypeset->get_mirrored_type( lcys_nterm_cterm ) );
		ResidueTypeCOP mirrored_dile_nacetyl( restypeset->get_mirrored_type( dile_nacetyl ) );
		ResidueTypeCOP mirrored_plain_c53( restypeset->get_mirrored_type( plain_c53 ) );

		TS_ASSERT_EQUALS( mirrored_lcys_nterm_cterm->name3(), "DCS" );
		TS_ASSERT_EQUALS( mirrored_lcys_nterm_cterm->name(), "DCYS:CtermProteinFull:NtermProteinFull" );
		TS_ASSERT_EQUALS( mirrored_lcys_nterm_cterm->base_name(), "DCYS" );

		TS_ASSERT_EQUALS( mirrored_dile_nacetyl->name3(), "ILE" );
		TS_ASSERT_EQUALS( mirrored_dile_nacetyl->name(), "ILE:AcetylatedNtermProteinFull" );
		TS_ASSERT_EQUALS( mirrored_dile_nacetyl->base_name(), "ILE" );

		TS_ASSERT_EQUALS( mirrored_plain_c53->name3(), "C53" ); //Currently, the D-patches don't update the name3 for noncanonicals.  Fix this if/when they do.
		TS_ASSERT_EQUALS( mirrored_plain_c53->name(), "DC53" );
		TS_ASSERT_EQUALS( mirrored_plain_c53->base_name(), "DC53" );

	}

};
