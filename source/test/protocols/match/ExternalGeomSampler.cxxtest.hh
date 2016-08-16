// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/protocols/match/ExternalGeomSampler.cxxtest.hh
/// @brief  test suite for protocols::match::ExternalGeomSampler
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>

#include <protocols/toolbox/match_enzdes_util/ExternalGeomSampler.hh>

// Utility headers
#include <utility/exit.hh>

/// Project headers
#include <core/types.hh>

// C++ headers
#include <string>
#include <iostream>

//Auto Headers
#include <utility/vector1.hh>


//using namespace protocols::match;
//using namespace protocols::match::downstream;
using namespace protocols::toolbox::match_enzdes_util;

// --------------- Test Class --------------- //

class ExternalGeomSamplerTests : public CxxTest::TestSuite {

public:
	typedef core::Real Real;
	typedef ExternalGeomSampler::HTReal HTReal;

public:

	utility::vector1< Real >  d_1;
	utility::vector1< Real >  ang_U2D1_1;
	utility::vector1< Real >  tor_U3D1_1;
	utility::vector1< Real >  ang_U1D2_1;
	utility::vector1< Real >  tor_U1D3_1;
	utility::vector1< Real >  tor_U2D2_1;

	Real dis_D1D2;
	Real dis_D2D3;
	Real ang_D1D2D3;

	// --------------- Fixtures --------------- //


	// Shared initialization goes here.
	void setUp() {

		d_1.clear();
		ang_U2D1_1.clear();
		tor_U3D1_1.clear();
		ang_U1D2_1.clear();
		tor_U1D3_1.clear();
		tor_U2D2_1.clear();


		d_1.push_back( 3.2 );
		ang_U2D1_1.push_back(  90.0 );
		tor_U3D1_1.push_back( 90.0  );
		ang_U1D2_1.push_back( 100.0 ); ang_U1D2_1.push_back( 110.0 ); ang_U1D2_1.push_back( 90.0 );
		tor_U1D3_1.push_back( 0.0 ); tor_U1D3_1.push_back( 30.0 ); tor_U1D3_1.push_back( 60.0 ); tor_U1D3_1.push_back( 90.0 ); tor_U1D3_1.push_back( 120.0 ); tor_U1D3_1.push_back( 150.0 );
		tor_U1D3_1.push_back(180.0); tor_U1D3_1.push_back(-30.0 ); tor_U1D3_1.push_back(-60.0 ); tor_U1D3_1.push_back(-90.0 ); tor_U1D3_1.push_back(-120.0 ); tor_U1D3_1.push_back(-150.0 );
		tor_U2D2_1.push_back( 90.0 );

		dis_D1D2 = 1.5;
		dis_D2D3 = 1.8;
		ang_D1D2D3 = 109.5;

	}

	// Shared finalization goes here.
	void tearDown() {
	}


	// --------------- Test Cases --------------- //
	void test_sampler_ctor() {

		ExternalGeomSampler sam;

		TS_ASSERT( sam.n_tor_U3D1_samples() == 0 );
		TS_ASSERT( sam.n_dis_U1D1_samples() == 0 );
		TS_ASSERT( sam.n_ang_U2D1_samples() == 0 );
		TS_ASSERT( sam.n_ang_U1D2_samples() == 0 );
		TS_ASSERT( sam.n_tor_U2D2_samples() == 0 );
		TS_ASSERT( sam.n_tor_U1D3_samples() == 0 );

	}

	void test_sampler_initialization() {

		ExternalGeomSampler sam;

		sam.set_dis_U1D1_samples( d_1 );
		sam.set_ang_U2D1_samples( ang_U2D1_1 );
		sam.set_tor_U3D1_samples( tor_U3D1_1 );
		sam.set_ang_U1D2_samples( ang_U1D2_1 );
		sam.set_tor_U1D3_samples( tor_U1D3_1 );
		sam.set_tor_U2D2_samples( tor_U2D2_1 );


		TS_ASSERT( sam.n_tor_U3D1_samples() == 1 );
		TS_ASSERT( sam.n_dis_U1D1_samples() == 1 );
		TS_ASSERT( sam.n_ang_U2D1_samples() == 1 );
		TS_ASSERT( sam.n_ang_U1D2_samples() == 3 );
		TS_ASSERT( sam.n_tor_U2D2_samples() == 1 );
		TS_ASSERT( sam.n_tor_U1D3_samples() == 12 );

	}

	void test_external_geom_sampler_transforms() {
		ExternalGeomSampler sam;
		sam.set_dis_U1D1_samples( d_1 );
		sam.set_ang_U2D1_samples( ang_U2D1_1 );
		sam.set_tor_U3D1_samples( tor_U3D1_1 );
		sam.set_ang_U1D2_samples( ang_U1D2_1 );
		sam.set_tor_U1D3_samples( tor_U1D3_1 );
		sam.set_tor_U2D2_samples( tor_U2D2_1 );
		sam.set_dis_D1D2(dis_D1D2);
		sam.set_dis_D2D3(dis_D2D3);
		sam.set_ang_D1D2D3(ang_D1D2D3);

		sam.precompute_transforms();


		core::Vector p1(0.0,0.0,0.0), p2(0.0,0.0,1.5), p3(0.0,0.9,2.4);

		HTReal start( p1, p2, p3 );

		HTReal zrotor_U3D1, xrotor_U3D1, xrotor_U1D3, zrotor_U2D2, zrotor_U1D3, ztransd;
		HTReal ang_D1D2D3xrot, ztrans45, ztrans56;
		zrotor_U3D1.set_zaxis_rotation_deg( 90 );
		xrotor_U3D1.set_xaxis_rotation_deg( -90 );
		xrotor_U1D3.set_xaxis_rotation_deg( -1 * (180.0 - 100) );
		zrotor_U2D2.set_zaxis_rotation_deg( 90 );
		ang_D1D2D3xrot.set_xaxis_rotation_deg( -1 * (180 - ang_D1D2D3 ));
		ztransd.walk_along_z( 3.2 );
		ztrans45.walk_along_z( dis_D1D2 );
		ztrans56.walk_along_z( dis_D2D3 );

		HTReal ht4 = start * zrotor_U3D1 * xrotor_U3D1 * ztransd;
		HTReal ht5 = ht4 * zrotor_U2D2 * xrotor_U1D3 * ztrans45;
		HTReal ht6 = ht5 * zrotor_U1D3 * ang_D1D2D3xrot * ztrans56;

		core::Vector p4( ht4.point()), p5( ht5.point()), p6( ht6.point());

		HTReal ex_ht4 = start * sam.transform( HT_tor_U3D1, 1 ) * sam.transform( HT_ang_U2D1, 1 );
		ex_ht4.walk_along_z( sam.dis_U1D1_samples()[ 1 ] );

		HTReal ex_ht5 = ex_ht4 * sam.transform( HT_tor_U2D2, 1 ) * sam.transform( HT_ang_U1D2, 1 );
		HTReal ex_ht6 = ex_ht5 * sam.transform( HT_tor_U1D3, 1 );


		core::Vector xp4( ex_ht4.point()), xp5( ex_ht5.point()), xp6( ex_ht6.point());

		/*
		std::cout << "ht4:" << std::endl;
		std::cout << "  " << ht4.xx() << " " << ht4.yx() << " " << ht4.zx() << " " << ht4.px() << std::endl;
		std::cout << "  " << ht4.xy() << " " << ht4.yy() << " " << ht4.zy() << " " << ht4.py() << std::endl;
		std::cout << "  " << ht4.xz() << " " << ht4.yz() << " " << ht4.zz() << " " << ht4.pz() << std::endl;

		std::cout << "ex_ht4:" << std::endl;
		std::cout << "  " << ex_ht4.xx() << " " << ex_ht4.yx() << " " << ex_ht4.zx() << " " << ex_ht4.px() << std::endl;
		std::cout << "  " << ex_ht4.xy() << " " << ex_ht4.yy() << " " << ex_ht4.zy() << " " << ex_ht4.py() << std::endl;
		std::cout << "  " << ex_ht4.xz() << " " << ex_ht4.yz() << " " << ex_ht4.zz() << " " << ex_ht4.pz() << std::endl;


		std::cout << "ht5:" << std::endl;
		std::cout << "  " << ht5.xx() << " " << ht5.yx() << " " << ht5.zx() << " " << ht5.px() << std::endl;
		std::cout << "  " << ht5.xy() << " " << ht5.yy() << " " << ht5.zy() << " " << ht5.py() << std::endl;
		std::cout << "  " << ht5.xz() << " " << ht5.yz() << " " << ht5.zz() << " " << ht5.pz() << std::endl;

		std::cout << "ex_ht5:" << std::endl;
		std::cout << "  " << ex_ht5.xx() << " " << ex_ht5.yx() << " " << ex_ht5.zx() << " " << ex_ht5.px() << std::endl;
		std::cout << "  " << ex_ht5.xy() << " " << ex_ht5.yy() << " " << ex_ht5.zy() << " " << ex_ht5.py() << std::endl;
		std::cout << "  " << ex_ht5.xz() << " " << ex_ht5.yz() << " " << ex_ht5.zz() << " " << ex_ht5.pz() << std::endl;


		std::cout << "ht6:" << std::endl;
		std::cout << "  " << ht6.xx() << " " << ht6.yx() << " " << ht6.zx() << " " << ht6.px() << std::endl;
		std::cout << "  " << ht6.xy() << " " << ht6.yy() << " " << ht6.zy() << " " << ht6.py() << std::endl;
		std::cout << "  " << ht6.xz() << " " << ht6.yz() << " " << ht6.zz() << " " << ht6.pz() << std::endl;

		std::cout << "ex_ht6:" << std::endl;
		std::cout << "  " << ex_ht6.xx() << " " << ex_ht6.yx() << " " << ex_ht6.zx() << " " << ex_ht6.px() << std::endl;
		std::cout << "  " << ex_ht6.xy() << " " << ex_ht6.yy() << " " << ex_ht6.zy() << " " << ex_ht6.py() << std::endl;
		std::cout << "  " << ex_ht6.xz() << " " << ex_ht6.yz() << " " << ex_ht6.zz() << " " << ex_ht6.pz() << std::endl;


		HTReal blah = zrotor_U2D2;// * xrotor_U1D3 * ztrans45;
		std::cout << "blah:" << std::endl;
		std::cout << "  " << blah.xx() << " " << blah.yx() << " " << blah.zx() << " " << blah.px() << std::endl;
		std::cout << "  " << blah.xy() << " " << blah.yy() << " " << blah.zy() << " " << blah.py() << std::endl;
		std::cout << "  " << blah.xz() << " " << blah.yz() << " " << blah.zz() << " " << blah.pz() << std::endl;

		HTReal blah2 = sam.transform( HT_tor_U2D2, 1 );// * sam.transform( HT_ang_U1D2, 1 );
		std::cout << "blah2:" << std::endl;
		std::cout << "  " << blah2.xx() << " " << blah2.yx() << " " << blah2.zx() << " " << blah2.px() << std::endl;
		std::cout << "  " << blah2.xy() << " " << blah2.yy() << " " << blah2.zy() << " " << blah2.py() << std::endl;
		std::cout << "  " << blah2.xz() << " " << blah2.yz() << " " << blah2.zz() << " " << blah2.pz() << std::endl;

		std::cout << "p4: " << p4.x() << " " << p4.y() << " " << p4.z() << std::endl;
		std::cout << "p5: " << p5.x() << " " << p5.y() << " " << p5.z() << std::endl;
		std::cout << "p6: " << p6.x() << " " << p6.y() << " " << p6.z() << std::endl;

		std::cout << "xp4: " << xp4.x() << " " << xp4.y() << " " << xp4.z() << std::endl;
		std::cout << "xp5: " << xp5.x() << " " << xp5.y() << " " << xp5.z() << std::endl;
		std::cout << "xp6: " << xp6.x() << " " << xp6.y() << " " << xp6.z() << std::endl;

		HTReal t1 = start * zrotor_U3D1;
		HTReal xt1 = start * sam.transform( HT_tor_U3D1, 1 );

		core::Vector pt1( t1.point() ), xpt1( xt1.point() );
		std::cout << "pt1: " << pt1.x() << " " << pt1.y() << " " << pt1.z() << std::endl;
		std::cout << "xpt1: " << xpt1.x() << " " << xpt1.y() << " " << xpt1.z() << std::endl;

		std::cout << "t1:" << std::endl;
		std::cout << "  " << t1.xx() << " " << t1.yx() << " " << t1.zx() << " " << t1.px() << std::endl;
		std::cout << "  " << t1.xy() << " " << t1.yy() << " " << t1.zy() << " " << t1.py() << std::endl;
		std::cout << "  " << t1.xz() << " " << t1.yz() << " " << t1.zz() << " " << t1.pz() << std::endl;

		std::cout << "xt1:" << std::endl;
		std::cout << "  " << xt1.xx() << " " << xt1.yx() << " " << xt1.zx() << " " << xt1.px() << std::endl;
		std::cout << "  " << xt1.xy() << " " << xt1.yy() << " " << xt1.zy() << " " << xt1.py() << std::endl;
		std::cout << "  " << xt1.xz() << " " << xt1.yz() << " " << xt1.zz() << " " << xt1.pz() << std::endl;

		HTReal sampA = sam.transform( HT_ang_U2D1, 1 );

		std::cout << "xrotor_U3D1:" << std::endl;
		std::cout << "  " << xrotor_U3D1.xx() << " " << xrotor_U3D1.yx() << " " << xrotor_U3D1.zx() << " " << xrotor_U3D1.px() << std::endl;
		std::cout << "  " << xrotor_U3D1.xy() << " " << xrotor_U3D1.yy() << " " << xrotor_U3D1.zy() << " " << xrotor_U3D1.py() << std::endl;
		std::cout << "  " << xrotor_U3D1.xz() << " " << xrotor_U3D1.yz() << " " << xrotor_U3D1.zz() << " " << xrotor_U3D1.pz() << std::endl;

		std::cout << "sampA:" << std::endl;
		std::cout << "  " << sampA.xx() << " " << sampA.yx() << " " << sampA.zx() << " " << sampA.px() << std::endl;
		std::cout << "  " << sampA.xy() << " " << sampA.yy() << " " << sampA.zy() << " " << sampA.py() << std::endl;
		std::cout << "  " << sampA.xz() << " " << sampA.yz() << " " << sampA.zz() << " " << sampA.pz() << std::endl;
		*/


		TS_ASSERT( p4.distance_squared( ex_ht4.point() ) < 1e-6 );
		TS_ASSERT( p5.distance_squared( ex_ht5.point() ) < 1e-6 );
		TS_ASSERT( p6.distance_squared( ex_ht6.point() ) < 1e-6 );
	}

};


