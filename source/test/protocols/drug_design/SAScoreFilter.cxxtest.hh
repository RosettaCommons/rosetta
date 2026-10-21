// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/drug_design/SAScoreFilter.cxxtest.hh
/// @brief  test for SAScoreFilter
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>

// Project Headers
#include <core/types.hh>

#include <protocols/drug_design/SAScoreFilter.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/residue_io.hh>

// Utility Headers
#include <utility/file/FileName.hh>
#include <basic/Tracer.hh>
#include <utility/io/izstream.hh>

#include <string>
#include <sstream>

#include <core/chemical/rdkit/RDKit.fwd.hh>

#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>

static basic::Tracer TR("protocols.drug_design.SAScoreFilter.cxxtest.hh");

// --------------- Test Class --------------- //

class SAScoreFilterTests : public CxxTest::TestSuite {

private:
public:

	void setUp() {
		core_init();
	}

	void tearDown() {
	}

	void test_larger() {
	}

	void test_sascores() {
		protocols::drug_design::SAScoreFilter sascorer;

		// These examples and values come from the RDKit source
		utility::io::izstream testfile("protocols/drug_design/zim.100.sascore.txt");


		std::string line;
		while ( getline(testfile, line) ) {
			std::string smiles, name;
			core::Real ref_score;
			std::stringstream linestream( line );
			linestream >> smiles >> name >> ref_score;
			if ( smiles == "smiles" ) { continue; }

			RDKit::ROMolOP rdmol( RDKit::SmilesToMol( smiles) );

			// I decided not to bother with testing the back-and-forth through Rosetta
			core::Real computed_score( sascorer.calculate_rdkit(*rdmol) );
			TR.Debug << "Testing " << smiles << std::endl;
			TSM_ASSERT_DELTA( smiles, computed_score, ref_score, 0.001  ); //listed in file to 3 decimal places
		}
	}


};
