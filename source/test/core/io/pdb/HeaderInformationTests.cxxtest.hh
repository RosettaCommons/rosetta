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
#include <core/io/pdb/HeaderInformation.hh>

// Program headers
#include <core/io/pdb/file_data.hh>
#include <core/io/pdb/pdb_dynamic_reader.hh>

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

static basic::Tracer TR("core.io.pdb.HeaderInformationTests.cxxtest");

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

		if(option.has(run::preserve_header)){
			cache_option__preserve_header_ = option[ run::preserve_header ]();
			option[ run::preserve_header ](true);
		}

	}

	void tearDown() {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		if(option.has(run::preserve_header)){
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
			<< "TITLE     SITTING ON A PARK BENCH EYEING LITTLE GIRLS WITH BAD INTENT." << endl
			<< "TITLE    2 SNOT RUNNING DOWN HIS NOSE GREASY FINGERS SMEARING SHABBY" << endl
			<< "TITLE    3 CLOTHES. DRYING IN THE COLD SUN WATCHING AS THE FRILLY" << endl
			<< "TITLE    4 PANTIES RUN. FEELING LIKE A DEAD DUCK SPITTING OUT PIECES" << endl
			<< "TITLE    5 OF HIS BROKEN LUCK." << endl
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
		using core::io::pdb::FileData;
		using core::io::pdb::PDB_DReader;
		using ObjexxFCL::rstripped_whitespace;

		PDB_DReader pdb_reader;
		FileData file_data(pdb_reader.createFileData(header_str_in));
		string header_str_out(pdb_reader.createPDBData(file_data));

		vector1<string> h_in(string_split(header_str_in, '\n'));
		vector1<string> h_out(string_split(header_str_out, '\n'));

		TS_ASSERT_EQUALS(h_in.size(),h_out.size());
		for(Size i = 1, ie= h_in.size(); i != ie; ++i){
			TS_ASSERT_EQUALS(h_in[i], rstripped_whitespace(h_out[i]));
		}
	}
};

