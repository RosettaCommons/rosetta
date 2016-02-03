// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/io/pdb/HeaderInformationTests.cxxtest.hh
/// @brief  test suite for classes associated with core::io::pdb::Field
/// @author Matthew O'Meara (mattjomeara@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit headers
#include <core/io/HeaderInformation.hh>

// Program headers
#include <core/io/StructFileRep.hh>
#include <core/io/pdb/pdb_reader.hh>
#include <core/io/pdb/pdb_writer.hh>

// Basic headers
#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/string_util.hh>
#include <utility/vector1.hh>

// ObjexxFCL headers
#include <ObjexxFCL/format.hh>


// C++ headers
#include <sstream>
#include <vector>

static basic::Tracer TR("core.io.HeaderInformationTests.cxxtest");

using namespace core;

class HeaderInformationTests : public CxxTest::TestSuite
{

public:
	bool cache_option__preserve_header_;


	// Shared initialization goes here.
	void setUp() {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		core_init();

		if ( option.has(run::preserve_header) ) {
			cache_option__preserve_header_ = option[ run::preserve_header ]();
			option[ run::preserve_header ](true);
		}

	}

	void tearDown() {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		if ( option.has(run::preserve_header) ) {
			option[ run::preserve_header ].value(cache_option__preserve_header_);
		}
	}


	// ------------------------------------------ //
	/// @brief test how PDB IO for header info
	void test_header_information_basic() {
		using std::string;
		using std::endl;

		std::stringstream header_in;
		header_in
			<< "HEADER    PEPTIDASE                               13-JAN-98   1A2Z" << endl
			<< "TITLE     LIFE IS LIKE A HURRICANE HERE IN DUCKBURG. RACE CARS," << endl
			<< "TITLE    2 LASERS, AEROPLANES, IT'S A DUCK-BLUR! MIGHT SOLVE A MYSTERY" << endl
			<< "TITLE    3 OR REWRITE HISTORY! DUCKTALES! WOO-OO!" << endl
			<< "KEYWDS    X-RAY STRUCTURE, NESG, Q5ZSV0, STRUCTURAL GENOMICS, PSI," << endl
			<< "KEYWDS   2 PROTEIN STRUCTURE INITIATIVE, NORTHEAST STRUCTURAL GENOMICS" << endl
			<< "KEYWDS   3 CONSORTIUM, PROTEIN BINDING" << endl
			<< "COMPND    MOL_ID: 1;" << endl
			<< "COMPND   2 MOLECULE: PYRROLIDONE CARBOXYL PEPTIDASE;" << endl
			<< "COMPND   3 CHAIN: A, B, C, D;" << endl
			<< "COMPND   4 FRAGMENT: RESIDUES 36-337;" << endl
			<< "COMPND   5 SYNONYM: PYROGLUTAMYL AMINOPEPTIDASE;" << endl
			<< "COMPND   6 EC: 3.4.19.3;" << endl
			<< "COMPND   7 ENGINEERED: YES;" << endl
			<< "COMPND   8 MUTATION: YES;" << endl
			<< "COMPND   9 OTHER_DETAILS: CHIMERA CONSISTING OF RESIDUES 1-366" << endl
			<< "COMPND  10 (SEQUENCE DATABASE RESIDUES 27-392) OF MALTOSE-BINDING" << endl
			<< "COMPND  11 PERIPLASMIC PROTEIN, A 5 ALANINE LINKER, AND RESIDUES 477-" << endl
			<< "COMPND  12 526 (SEQUENCE DATABASE RESIDUES 77-126) OF HOMEOBOX OF" << endl
			<< "COMPND  13 MATING-TYPE PROTEIN A-1." << endl
			<< "EXPDTA    X-RAY DIFFRACTION; FIBER DIFFRACTION; NEUTRON DIFFRACTION;" << endl
			<< "EXPDTA   2 ELECTRON CRYSTALLOGRAPHY; ELECTRON MICROSCOPY; SOLID-STATE" << endl
			<< "EXPDTA   3 NMR; SOLUTION NMR; SOLUTION SCATTERING; THEORETICAL MODEL" << endl
			<< "TER" << endl;

		string header_str_in(header_in.str());

		using utility::string_split;
		using utility::vector1;
		using core::io::StructFileRep;
		using core::io::StructFileRepOP;
		using ObjexxFCL::rstripped_whitespace;

		StructFileRep file_data( core::io::pdb::create_sfr_from_pdb_file_contents( header_str_in ) );
		string header_str_out( core::io::pdb::create_pdb_contents_from_sfr( file_data ) );

		vector1<string> h_in(string_split(header_str_in, '\n'));
		vector1<string> h_out(string_split(header_str_out, '\n'));

		TS_ASSERT_EQUALS(h_in.size(),h_out.size());
		for ( Size i = 1, ie= h_in.size(); i != ie; ++i ) {
			TS_ASSERT_EQUALS(h_in[i], rstripped_whitespace(h_out[i]));
		}
	}
};

