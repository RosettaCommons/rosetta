// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/make_rot_lib/MakeRotLib.cxxtest.hh
/// @brief  test for MakeRotLib
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/UTracer.hh>

// Unit header
#include <protocols/make_rot_lib/MakeRotLib.hh>


#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>

#include <utility/io/izstream.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR("protocols.make_rot_lib.cxxtest.hh");

// returns true on error.
bool load_RotData( std::istream & in, protocols::MakeRotLib::RotVec & rotvec ) {
	using namespace protocols::MakeRotLib;
	char next;
	std::string line;
	while ( in.good() ) {
		while ( next=in.peek(), next == ' ' || next == '\n' || next == '\t' ) { in.get(); } // Discard leading whitespace
		if ( ! in.good() ) { break; }
		if ( in.peek() == '#' ) {
			getline( in,line ); // Discard the comment line
			continue;
		}
		RotData rd(0,0);
		if( rd.load( in ) ) { return true; }
		rotvec.push_back( rd );
		TR.Debug << "Loaded " << rotvec.size() << " RotData objects." << std::endl;
	}
	return false;
}

// returns true on error.
bool load_RotData( std::string const & filename, protocols::MakeRotLib::RotVec & rotvec ) {
	utility::io::izstream stream( filename );
	return load_RotData( stream, rotvec );
}

// --------------- Test Class --------------- //

class MakeRotLibTests : public CxxTest::TestSuite {

public:

	// --------------- Fixtures --------------- //

	void setUp() {
		core_init();
	}

	void tearDown() {
	}

	// ------------- Helper Functions ------------- //


	// --------------- Test Cases --------------- //

	// Note: These test results were made from existing behavior, and may differ from "ideal" behavior.

	void test_init_rotamers_centroids() {
		using namespace protocols::MakeRotLib;
		using namespace core;

		test::UTracer UT_rot("protocols/make_rot_lib/init_rotamers.txt");
		test::UTracer UT_cen("protocols/make_rot_lib/init_centroids.txt");

		RotVec rotamers, centroids;
		Size ncluster( 0 );
		std::string aa_name;

		init_rotamers_centroids( rotamers, centroids, ncluster, "protocols/make_rot_lib/NVL_rot_lib_options_-60_-40.in", aa_name, false, 180, 180 );

		TS_ASSERT_EQUALS(ncluster, 9);
		TS_ASSERT_EQUALS(aa_name, "NVL:MethylatedCtermProteinFull:AcetylatedNtermProteinFull");

		for(core::Size ii(1); ii <= rotamers.size(); ++ii ) {
			UT_rot << "# Rotamer " << ii << std::endl;
			rotamers[ii].show(UT_rot);
		}

		for(core::Size jj(1); jj <= centroids.size(); ++jj ) {
			UT_cen << "# Centroid " << jj << std::endl;
			centroids[jj].show(UT_cen);
		}

	}

	// This should get cut down -- right now I'm doing this to test issues with the integration test.
	// As might be expected from minimization, this shows sensitivity to platform, \
	// having different results on the test server Mac clang build.
	void DISABLED_min_rotamers() {
		using namespace protocols::MakeRotLib;
		using namespace core;
		using namespace core::scoring;

		test::UTracer UT("protocols/make_rot_lib/min_rotamers.txt");

		ScoreFunctionOP scrfxn( getScoreFunctionLegacy( PRE_TALARIS_2013_STANDARD_WTS ) );
		scrfxn->set_weight( fa_dun, 0.0 );
		scrfxn->set_weight( p_aa_pp, 0.0 );
		scrfxn->set_weight( rama, 0.0 );
		scrfxn->set_weight( fa_intra_rep, 0.0 );
		scrfxn->set_weight( fa_intra_atr, 0.0 );
		scrfxn->set_weight( fa_rep, 0.0 );
		scrfxn->set_weight( fa_atr, 0.0 );
		scrfxn->set_weight( mm_twist, 1.0 );
		scrfxn->set_weight( mm_lj_inter_rep, 1.0 );
		scrfxn->set_weight( mm_lj_inter_atr, 1.0 );
		scrfxn->set_weight( mm_lj_intra_rep, 1.0 );
		scrfxn->set_weight( mm_lj_intra_atr, 1.0 );

		RotVec rotamers;
		TS_ASSERT_EQUALS( load_RotData( "protocols/make_rot_lib/init_rotamers.txt", rotamers ), false);

		min_rotamers( rotamers, scrfxn, "NVL:MethylatedCtermProteinFull:AcetylatedNtermProteinFull" );

		for(core::Size ii(1); ii <= rotamers.size(); ++ii ) {
			UT << "# Min Rotamer " << ii << std::endl;
			rotamers[ii].show(UT);
		}
	}

	void test_calc_all_dist() {
		using namespace protocols::MakeRotLib;
		using namespace core;
		using namespace core::scoring;

		test::UTracer UT("protocols/make_rot_lib/all_dist_rotamers.txt");
		RotVec rotamers, centroids;
		TS_ASSERT_EQUALS( load_RotData( "protocols/make_rot_lib/min_rotamers.txt", rotamers ), false);
		TS_ASSERT_EQUALS( load_RotData( "protocols/make_rot_lib/init_centroids.txt", centroids ), false);

		calc_all_dist( rotamers, centroids );

		for(core::Size ii(1); ii <= rotamers.size(); ++ii ) {
			UT << "# All dist Rotamer " << ii << std::endl;
			rotamers[ii].show(UT);
    }
	}

	void test_calc_rotamer_clusters() {
		using namespace protocols::MakeRotLib;
		using namespace core;
		using namespace core::scoring;

		test::UTracer UT("protocols/make_rot_lib/cluster_rotamers.txt");
		RotVec rotamers;
		TS_ASSERT_EQUALS( load_RotData( "protocols/make_rot_lib/all_dist_rotamers.txt", rotamers ), false);

		calc_rotamer_clusters( rotamers );

		for(core::Size ii(1); ii <= rotamers.size(); ++ii ) {
			UT << "# Calc Rotamer Clusters Rotamer " << ii << std::endl;
			rotamers[ii].show(UT);
    }
	}


	void test_calc_centroids() {
		using namespace protocols::MakeRotLib;
		using namespace core;
		using namespace core::scoring;

		test::UTracer UT("protocols/make_rot_lib/calc_centroids.txt");
		RotVec rotamers, centroids;
		TS_ASSERT_EQUALS( load_RotData( "protocols/make_rot_lib/cluster_rotamers.txt", rotamers ), false);
		TS_ASSERT_EQUALS( load_RotData( "protocols/make_rot_lib/init_centroids.txt", centroids ), false);

		calc_centroids( rotamers, centroids );

		for(core::Size ii(1); ii <= centroids.size(); ++ii ) {
			UT << "# Calc centroids centroid " << ii << std::endl;
			centroids[ii].show(UT);
    }
	}

};//end class
