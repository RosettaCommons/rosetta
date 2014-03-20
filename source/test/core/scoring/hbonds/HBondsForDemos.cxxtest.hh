// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/hbonds/HBonds.cxxtest.hh
/// @brief  Test evaluation of hbond potentials for each HBEvalType across parameter ranges.
/// @author Matthew O'Meara

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/UTracer.hh>
#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>

// Package headers
#include <core/scoring/hbonds/HBondSet.hh>

// Project headers
#include <core/types.hh>

//Auto Headers
#include <core/scoring/hbonds/HBondOptions.hh>


using namespace core;
  using namespace conformation;
  using namespace scoring;
    using namespace hbonds;

using pose::Pose;

static basic::Tracer TR("core.scoring.hbonds.HBondsForDemos.cxxtest");

class HBondsForDemos : public CxxTest::TestSuite {

public:
  void setUp() {
		core_init();
	}

	void tearDown(){}


	void test_easy_construction_printing_hbond_set(){
		using std::endl;


		// Construct an hbond set
		Pose pose = create_trpcage_ideal_pose();
		HBondSet hbond_set(pose);

		TR << "Print minimal hbond set information:" << endl;
		TR << hbond_set;
		TR << endl;

		TR << "Same info just using the .show() this time.";
		hbond_set.show(TR);
		TR << endl;

		TR << "Print interprable hbond set information:" << endl;
		hbond_set.show(pose, true, TR);
		TR << endl;

		TR << "Print all the hbonds involving a few specific residues: (3,4,5,6)" << endl;
		hbond_set.show(pose, 3, /*print header*/ true,  TR);
		hbond_set.show(pose, 4, /*print header*/ false, TR);
		hbond_set.show(pose, 5, /*print header*/ false, TR);
		hbond_set.show(pose, 6, /*print header*/ false, TR);

		/* Expected output last set of show statements

    #Dch Dn Dres Da  Ach An Ares Aa  length AHDang BAHang  BAtor weight energy
    #A    7  LEU  N  A    3  TYR  O    2.03  146.5  141.5 -122.2  1.000 -1.266
    #A    8  LYS  N  A    4  ILE  O    2.07  156.6  148.9 -115.6  1.000 -1.442
    #A    5  GLN  N  A    1  ASN  O    2.03  153.9  154.4 -143.8  1.000 -1.424
    #A    6  TRP  N  A    2  LEU  O    2.19  142.9  142.5 -149.9  1.000 -0.980
    #A    9  ASP  N  A    6  TRP  O    2.55  141.1   95.4 -106.5  1.000 -0.137
    #A   11  GLY  N  A    6  TRP  O    2.37  144.0  141.8  162.9  1.000 -0.393

		*/

	}
};
