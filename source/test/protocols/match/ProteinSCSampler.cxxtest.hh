// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/protocols/match/ProteinSCSampler.cxxtest.hh
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>

#include <protocols/match/upstream/ProteinSCSampler.hh>
#include <protocols/match/upstream/OriginalScaffoldBuildPoint.hh>

#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Utility headers
#include <utility/exit.hh>

/// Project headers
#include <core/types.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>

// C++ headers
#include <string>
#include <iostream>

//Auto Headers
#include <core/pack/dunbrack/DunbrackRotamer.hh>
#include <utility/vector1.hh>


using namespace protocols::match;
using namespace protocols::match::upstream;


// --------------- Test Class --------------- //

class ProteinSCSamplerTests : public CxxTest::TestSuite {

	public:


	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.


	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	// Shared finalization goes here.
	void tearDown() {
	}


	// --------------- Test Cases --------------- //
	void test_sc_sampler_ctor() {
		using namespace core;

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		OriginalBackboneBuildPointOP res2bp( new OriginalBackboneBuildPoint( trpcage.residue( 2 ), 1 ) );

		DunbrackSCSampler dunsampler;
		DunbrackSCSampler::DunbrackRotamerSampleDataVector samps1(
			dunsampler.samples( *res2bp, trpcage.residue_type( 2 ) ));

		//std::cout << "Sample 1: " << trpcage.residue_type( 2 ).name() << " " << samps1.size() << std::endl;
		TS_ASSERT( samps1.size() == 9 );

		/*Real prob_cummulative( 0.0 );
		for ( Size ii = 1; ii <= samps1.size(); ++ii ) {
			std::cout << ii << ": prob= " << samps1[ ii ].probability();
			prob_cummulative +=  samps1[ ii ].probability();
			std::cout << " probcummulative= " << prob_cummulative;

			for ( Size jj = 1; jj <= samps1[ ii ].nchi(); ++jj ) {
				std::cout << " ( " << samps1[ ii ].rot_well()[ jj ] << ", " <<  samps1[ ii ].chi_mean()[ jj ] << ", "  << samps1[ ii ].chi_sd()[ jj ] << " )";
			}
			std::cout << std::endl;
		}*/


		DunbrackSCSampler::DunbrackRotamerSampleDataVector samps2(
			dunsampler.samples( *res2bp, trpcage.residue_type( 3 ) ));

		//std::cout << "Sample 2: " << trpcage.residue_type( 3 ).name() << " " << samps2.size() << std::endl;
		TS_ASSERT( samps2.size() == 18 );

		/*prob_cummulative = 0.0;
		for ( Size ii = 1; ii <= samps2.size(); ++ii ) {
			std::cout << ii << ": prob= " << samps2[ ii ].probability();
			prob_cummulative +=  samps2[ ii ].probability();
			std::cout << " probcummulative= " << prob_cummulative;
			for ( Size jj = 1; jj <= samps2[ ii ].nchi(); ++jj ) {
				std::cout << " ( " << samps2[ ii ].rot_well()[ jj ] << ", " <<  samps2[ ii ].chi_mean()[ jj ] << ", "  << samps2[ ii ].chi_sd()[ jj ] << " )";
			}
			std::cout << std::endl;
		}*/


		DunbrackSCSampler::DunbrackRotamerSampleDataVector samps3(
			dunsampler.samples( *res2bp, trpcage.residue_type( 4 ) ));

		//std::cout << "Sample 3: " << trpcage.residue_type( 4 ).name() << " " << samps3.size() << std::endl;
		TS_ASSERT( samps3.size() == 9 );

		/*prob_cummulative = 0.0;
		for ( Size ii = 1; ii <= samps3.size(); ++ii ) {
			std::cout << ii << ": prob= " << samps3[ ii ].probability();
			prob_cummulative +=  samps3[ ii ].probability();
			std::cout << " probcummulative= " << prob_cummulative;
			for ( Size jj = 1; jj <= samps3[ ii ].nchi(); ++jj ) {
				std::cout << " ( " << samps3[ ii ].rot_well()[ jj ] << ", " <<  samps3[ ii ].chi_mean()[ jj ] << ", "  << samps3[ ii ].chi_sd()[ jj ] << " )";
			}
			std::cout << std::endl;
		}*/

	}

	void test_sc_sampler_desymmeterize() {
		using namespace core;

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		OriginalBackboneBuildPointOP res2bp( new OriginalBackboneBuildPoint( trpcage.residue( 2 ), 1 ) );

		DunbrackSCSampler dunsampler_symm, dunsampler_desymm;
		dunsampler_desymm.set_desymmeterize( true );

		TS_ASSERT( dunsampler_symm.desymmeterize() == false );
		TS_ASSERT( dunsampler_desymm.desymmeterize() == true );

		/// residue 9 on trpcage is ASP
		DunbrackSCSampler::DunbrackRotamerSampleDataVector asp_samps_symm(
			dunsampler_symm.samples( *res2bp, trpcage.residue_type( 9 ) ));

		DunbrackSCSampler::DunbrackRotamerSampleDataVector asp_samps_desymm(
			dunsampler_desymm.samples( *res2bp, trpcage.residue_type( 9 ) ));

		TS_ASSERT( asp_samps_symm.size() * 2 == asp_samps_desymm.size() );

		for ( core::Size ii = 1; ii <= asp_samps_symm.size(); ++ii ) {
			core::Size ii1 = 2*ii - 1;
			core::Size ii2 = 2*ii;

			TS_ASSERT_DELTA( asp_samps_symm[ ii ].probability() * 0.5, asp_samps_desymm[ ii1 ].probability(), 1e-6 );
			TS_ASSERT_DELTA( asp_samps_symm[ ii ].probability() * 0.5, asp_samps_desymm[ ii2 ].probability(), 1e-6 );

			TS_ASSERT_DELTA( asp_samps_symm[ ii ].chi_mean()[2]      , asp_samps_desymm[ ii1 ].chi_mean()[2], 1e-6 );
			TS_ASSERT_DELTA( asp_samps_symm[ ii ].chi_mean()[2] + 180, asp_samps_desymm[ ii2 ].chi_mean()[2], 1e-6 );

			TS_ASSERT( asp_samps_symm[ ii ].nrchi_sample() );
			if ( asp_samps_symm[ ii ].nrchi_sample() ) {
				TS_ASSERT_DELTA( asp_samps_symm[ ii ].nrchi_lower_boundary()      , asp_samps_desymm[ ii1 ].nrchi_lower_boundary(), 1e-6 );
				TS_ASSERT_DELTA( asp_samps_symm[ ii ].nrchi_lower_boundary() + 180, asp_samps_desymm[ ii2 ].nrchi_lower_boundary(), 1e-6 );

				TS_ASSERT_DELTA( asp_samps_symm[ ii ].nrchi_upper_boundary()      , asp_samps_desymm[ ii1 ].nrchi_upper_boundary(), 1e-6 );
				TS_ASSERT_DELTA( asp_samps_symm[ ii ].nrchi_upper_boundary() + 180, asp_samps_desymm[ ii2 ].nrchi_upper_boundary(), 1e-6 );
			}
		}

	}


};
