// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief Test core::scoring::disulfides::FullatomDisulfidePotential
/// @author Spencer Bliven <blivens@u.washington.edu>
/// @created October 2008

// Test headers
#include <cxxtest/TestSuite.h>

// Core headers
#include <core/types.hh>
#include <core/scoring/disulfides/FullatomDisulfidePotential.hh>
#include <core/scoring/func/Func.hh>

//Utility
#include <numeric/constants.hh>
#include <test/UTracer.hh>

//Auto Headers
#include <ObjexxFCL/FArray1D.hh>


//C++

using namespace core;
using namespace std;
using namespace core::scoring::disulfides;
using namespace numeric::constants::d;
using test::UTracer;


class FullatomDisulfidePotentialTest : public CxxTest::TestSuite {
public:

	FullatomDisulfidePotentialTest() {}


	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.

	// Shared initialization goes here.
	void setUp() {
	}

	// Shared finalization goes here.
	void tearDown() {
	}


	// --------------- Test Cases --------------- //

	void test_ss_dst()
	{
		UTracer UT("core/scoring/dslf_ss_dst.u");

		Real dslf_ss_dst_min(1.5);
		Real dslf_ss_dst_max(3.0);
		Real dslf_ss_dst_step(.0025);

		scoring::constraints::FuncOP f( new SG_Dist_Func );

		UT << "X\tf(X)\tdf/dx\n";
		for(Real x(dslf_ss_dst_min); x<=dslf_ss_dst_max; x+=dslf_ss_dst_step) {
			Energy score( f->func(x) );
			Energy deriv( f->dfunc(x));
			UT << x << "\t"
				<< score << "\t"
				<< deriv << "\n";
		}

	}

	void test_cs_ang() {
		UTracer UT("core/scoring/dslf_cs_ang.u");

		Real dslf_cs_ang_min(0);
		Real dslf_cs_ang_max(180);
		Real dslf_cs_ang_step(.5);

		scoring::constraints::FuncOP f(new CB_Angle_Func);

		UT << "X\tf(X)\tdf/dx\n";
		for(Real x(dslf_cs_ang_min); x<=dslf_cs_ang_max; x+=dslf_cs_ang_step) {
			Energy score( f->func(x*degrees_to_radians) );
			Energy deriv( f->dfunc(x*degrees_to_radians));
			UT << x << "\t"
				<< score << "\t"
				<< deriv << std::endl;
		}

	}

	void test_ca_dih() {
		UTracer UT("core/scoring/dslf_ca_dih.u");

		Real dslf_dih_min(-180);
		Real dslf_dih_step(1.8);
		Real dslf_dih_max(180);


		scoring::constraints::FuncOP f(new CBSG_Dihedral_Func);

		UT << "X\tf(X)\tdf/dx\n";
		for(Real x(dslf_dih_min); x<=dslf_dih_max; x+=dslf_dih_step) {
			Energy score( f->func(x*degrees_to_radians) );
			Energy deriv( f->dfunc(x*degrees_to_radians));
			UT << x << "\t"
				<< score << "\t"
				<< deriv << std::endl;
		}
	}

	void test_ss_dih() {
		UTracer UT("core/scoring/dslf_ss_dih.u");

		Real dslf_dih_min(-180);
		Real dslf_dih_step(1.8);
		Real dslf_dih_max(180);

		scoring::constraints::FuncOP f(new SGSG_Dihedral_Func);
		UT << "X\tf(X)\tdf/dx\n";
		for(Real x(dslf_dih_min); x<=dslf_dih_max; x+=dslf_dih_step) {
			Energy score( f->func(x*degrees_to_radians) );
			Energy deriv( f->dfunc(x*degrees_to_radians));
			UT << x << "\t"
				<< score << "\t"
				<< deriv << std::endl;
		}

	}

};



