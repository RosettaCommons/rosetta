// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/interpolation/spline
/// @brief  test suite for numeric::interpolation::spline::Bicubic_spline
/// @author Steven Combs (steven.combs@vanderbilt.edu)
/// This tests the functions that are in the bicubic spline class except for
/// the e_periodic steps.



// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit headers
#include <numeric/types.hh>
//#include <numeric/interpolation/spline/Bicubic_spline.hh>
//#include <numeric/interpolation/spline/Cubic_spline.hh>
//#include <numeric/interpolation/spline/Cubic_spline.fwd.hh>
//#include <numeric/MathVector_operations.hh>

#include <basic/Tracer.hh>

#include <cassert>


//static basic::Tracer TR( "numeric.interpolation.spline.bicubic_interpolation_function_cxxtest_hh" );

// --------------- Test Class --------------- //

using namespace numeric;
class bicubic_interpolation_function_tests : public CxxTest::TestSuite {

public:
	//shared data


	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.

	// Shared initialization goes here.
	void setUp()
	{

	} //Match contents of Histogram_sample.hist


	// Shared finalization goes here.
	void tearDown() {

	}


/// @brief Interpolate in a grid with the values, and second derivatives given, and
/// simultaneously evaluate the derivative.  No option for working with periodic ranges.
/// Instead, make sure that interpolation doesn't need to span > 180 degrees.
void bicubic_interpolation(
	Real v00, Real d2dx200, Real d2dy200, Real d4dx2y200,
	Real v01, Real d2dx201, Real d2dy201, Real d4dx2y201,
	Real v10, Real d2dx210, Real d2dy210, Real d4dx2y210,
	Real v11, Real d2dx211, Real d2dy211, Real d4dx2y211,
	Real dxp, // in the range [0..1) representing the distance to the left bin boundary
	Real dyp, // in the range [0..1) representing the distance to the lower bin boundary
	Real binwx, // the size of the bin witdh for x
	Real binwy, // the size of the bin width for y
	Real & val,
	Real & dvaldx,
	Real & dvaldy
)
{
	assert( dxp >= 0 && dxp < 1.0 );
	assert( dyp >= 0 && dyp < 1.0 );
	Real dxm = 1-dxp;
	Real dym = 1-dyp;
	Real binwx_over6 = binwx/6;
	Real binwy_over6 = binwy/6;
	Real dx3p = ( dxp*dxp*dxp - dxp) * binwx*binwx_over6;
	Real dx3m = ( dxm*dxm*dxm - dxm) * binwx*binwx_over6;
	Real dy3p = ( dyp*dyp*dyp - dyp) * binwy*binwy_over6;
	Real dy3m = ( dym*dym*dym - dym) * binwy*binwy_over6;
	Real invbinwx = 1/binwx;
	Real invbinwy = 1/binwy;

	/*	val =
		dxm *   ( dym * values_( (   i - 1) % dimx, ( j - 1) % dimy) + dyp  * values_( (  i - 1) % dimx, j % dimy))
		+ dxp * ( dym * values_(    i      % dimx, ( j - 1) % dimy) + dyp  * values_(   i      % dimx, j % dimy))
		+dx3m * ( dym * dsecox_( (   i - 1) % dimx, ( j - 1) % dimy) + dyp  * dsecox_( (  i - 1) % dimx, j % dimy))
		+dx3p * ( dym * dsecox_(    i      % dimx, ( j - 1) % dimy) + dyp  * dsecox_(   i      % dimx, j % dimy))
		+ dxm * ( dy3m * dsecoy_( (  i - 1) % dimx, ( j - 1) % dimy) + dy3p * dsecoy_( (  i - 1) % dimx, j % dimy))
		+ dxp * ( dy3m * dsecoy_(   i      % dimx, ( j - 1) % dimy) + dy3p * dsecoy_(   i      % dimx, j % dimy))
		+dx3m * ( dy3m * dsecoxy_( ( i - 1) % dimx, ( j - 1) % dimy) + dy3p * dsecoxy_( ( i - 1) % dimx, j % dimy))
		+dx3p * ( dy3m * dsecoxy_(  i      % dimx, ( j - 1) % dimy) + dy3p * dsecoxy_(  i      % dimx, j % dimy))*/
	val =
		dxm *   (  dym *     v00   +  dyp *     v01 )
		+ dxp * (  dym *     v10   +  dyp *     v11 )
		+dx3m * (  dym * d2dx200   +  dyp * d2dx201 )
		+dx3p * (  dym * d2dx210   +  dyp * d2dx211 )
		+ dxm * ( dy3m * d2dy200   + dy3p * d2dy201 )
		+ dxp * ( dy3m * d2dy210   + dy3p * d2dy211 )
		+dx3m * ( dy3m * d4dx2y200 + dy3p * d4dx2y201 )
		+dx3p * ( dy3m * d4dx2y210 + dy3p * d4dx2y211 );

	/*dvaldx = -( dym * values_( ( i - 1) % dimx, ( j - 1) % dimy) + dyp * values_( ( i - 1) % dimx  , j % dimy)) / delta_[ 0]
		+( dym * values_( i % dimx      , ( j - 1) % dimy) + dyp * values_( i % dimx        , j % dimy)) / delta_[ 0]
		- ( 3 * dxm * dxm - 1) * delta_[ 0] / 6 *( dym*dsecox_( ( i - 1) % dimx, ( j - 1) % dimy) + dyp * dsecox_( ( i - 1) % dimx, j%dimy))
		+ ( 3 * dxp * dxp - 1) * delta_[ 0] / 6 *( dym*dsecox_( i % dimx      , ( j - 1) % dimy) + dyp * dsecox_( i % dimx      , j%dimy))
		-( dy3m * dsecoy_( ( i - 1) % dimx , ( j - 1) % dimy) + dy3p * dsecoy_( ( i-1) % dimx , j % dimy)) / delta_[ 0]
		+( dy3m * dsecoy_( i % dimx       , ( j - 1) % dimy) + dy3p * dsecoy_( i % dimx     , j % dimy)) / delta_[ 0]
		- ( 3 * dxm * dxm - 1) * delta_[ 0] / 6 *( dy3m*dsecoxy_( ( i - 1)%dimx, ( j - 1) % dimy) + dy3p * dsecoxy_( ( i - 1) % dimx, j % dimy))
		+ ( 3 * dxp * dxp - 1) * delta_[ 0] / 6 *( dy3m*dsecoxy_( i % dimx    , ( j - 1) % dimy) + dy3p * dsecoxy_( i % dimx      , j % dimy));*/

	dvaldx =
		-( dym * v00 + dyp * v01 ) * invbinwx
		+( dym * v10 + dyp * v11 ) * invbinwx
		- ( 3 * dxm * dxm - 1) * binwx_over6 *( dym  * d2dx200   + dyp  * d2dx201 )
		+ ( 3 * dxp * dxp - 1) * binwx_over6 *( dym  * d2dx210   + dyp  * d2dx211 )
		-( dy3m * d2dy200 + dy3p * d2dy201 ) * invbinwx
		+( dy3m * d2dy210 + dy3p * d2dy211 ) * invbinwx
		- ( 3 * dxm * dxm - 1) * binwx_over6 *( dy3m * d4dx2y200 + dy3p * d4dx2y201 )
		+ ( 3 * dxp * dxp - 1) * binwx_over6 *( dy3m * d4dx2y210 + dy3p * d4dx2y211 );

	/*dvaldy =
		dxm *( -values_( ( i-1)%dimx  , (j-1)%dimy)+values_( (i-1)%dimx  , j%dimy))/delta_[ 1]
		+ dxp *( -values_( i%dimx      , (j-1)%dimy)+values_(i%dimx      , j%dimy))/delta_[ 1]
		+dx3m *( -dsecox_( ( i-1)%dimx  , (j-1)%dimy)+dsecox_( (i-1)%dimx  , j%dimy))/delta_[ 1]
		+dx3p *( -dsecox_( i%dimx      , (j-1)%dimy)+dsecox_(i%dimx      , j%dimy))/delta_[ 1]
		+ dxm *( -( 3 * dym * dym - 1) * dsecoy_( ( i-1)%dimx , ( j - 1)% dimy) +( 3 * dyp * dyp - 1) * dsecoy_( ( i-1)%dimx , j % dimy)) * delta_[ 1]/ 6
		+ dxp *( -( 3 * dym * dym - 1) * dsecoy_( i%dimx     , ( j - 1)% dimy) +( 3 * dyp * dyp - 1) * dsecoy_( i%dimx     , j % dimy)) * delta_[ 1]/ 6
		+ dx3m *( -( 3 * dym * dym - 1) * dsecoxy_( ( i-1)%dimx, ( j - 1)% dimy) +( 3 * dyp * dyp - 1) * dsecoxy_( ( i-1)%dimx, j % dimy)) * delta_[ 1]/ 6
		+dx3p *( -( 3 * dym * dym - 1) * dsecoxy_( i%dimx    , ( j - 1)% dimy) +( 3 * dyp * dyp - 1) * dsecoxy_( i%dimx    , j % dimy)) * delta_[ 1]/ 6;*/

	dvaldy =
		dxm   *( -v00 + v01 ) * invbinwy
		+ dxp *( -v10 + v11 ) * invbinwy
		+dx3m *( -d2dx200 + d2dx201) * invbinwy
		+dx3p *( -d2dx210 + d2dx211) * invbinwy
		+ dxm *( -( 3 * dym * dym - 1) * d2dy200   + ( 3 * dyp * dyp - 1) * d2dy201 ) * binwy_over6
		+ dxp *( -( 3 * dym * dym - 1) * d2dy210   + ( 3 * dyp * dyp - 1) * d2dy211 ) * binwy_over6
		+dx3m *( -( 3 * dym * dym - 1) * d4dx2y200 + ( 3 * dyp * dyp - 1) * d4dx2y201 ) * binwy_over6
		+dx3p *( -( 3 * dym * dym - 1) * d4dx2y210 + ( 3 * dyp * dyp - 1) * d4dx2y211 ) * binwy_over6;

	//TR << "bicubic interpolation " << val << " " << dvaldx << " " << dvaldy << std::endl;
}


  void test_bicubic_interpolation_function() {
    Real v00( 11.86960029602051), d2dx200( -0.0112837441265583 ),  d2dy200( -0.04301836714148521 ), d4dx2y200  ( 0.005927249323576689 );
    Real v01( 9.21034049987793),  d2dx201( 0.1994424015283585),    d2dy201( 0.00607264693826437),   d4dx2y201(-0.00843812245875597);
    Real v10( 11.25056076049805), d2dx210( -0.002053411677479744), d2dy210(  0.06962345540523529),  d4dx2y210(  -0.004980601370334625);
    Real v11( 12.20607280731201), d2dx211( -0.1664439588785172),   d2dy211( -0.1316054910421371),   d4dx2y211(  0.006051852367818356 );
    Real dxp( 0.9034520028470581 );
    Real dyp( 0.9236237223457749 );
    Real binwx( 10.);
    Real binwy( 10.);

    Real val, dvaldx, dvaldy;

    bicubic_interpolation(
       v00,  d2dx200,  d2dy200,  d4dx2y200,
       v01,  d2dx201,  d2dy201,  d4dx2y201,
       v10,  d2dx210,  d2dy210,  d4dx2y210,
       v11,  d2dx211,  d2dy211,  d4dx2y211,
       dxp, dyp, binwx, binwy,
       val, dvaldx, dvaldy
    );

    Real delta=1e-6;

    Real valxm, dvaldxxm, dvaldyxm;
    bicubic_interpolation(
       v00,  d2dx200,  d2dy200,  d4dx2y200,
       v01,  d2dx201,  d2dy201,  d4dx2y201,
       v10,  d2dx210,  d2dy210,  d4dx2y210,
       v11,  d2dx211,  d2dy211,  d4dx2y211,
       dxp-delta/binwx, dyp, binwx, binwy,
       valxm, dvaldxxm, dvaldyxm
    );

    Real valxp, dvaldxxp, dvaldyxp;
    bicubic_interpolation(
       v00,  d2dx200,  d2dy200,  d4dx2y200,
       v01,  d2dx201,  d2dy201,  d4dx2y201,
       v10,  d2dx210,  d2dy210,  d4dx2y210,
       v11,  d2dx211,  d2dy211,  d4dx2y211,
       dxp+delta/binwx, dyp, binwx, binwy,
       valxp, dvaldxxp, dvaldyxp
    );

    //TR << "dvaldx analytic " << dvaldx << " and numeric " << (valxp-valxm)/(2*delta) << std::endl;
    TS_ASSERT_DELTA( dvaldx, (valxp-valxm)/(2*delta), 1e-6 );


    Real valym, dvaldxym, dvaldyym;
    bicubic_interpolation(
       v00,  d2dx200,  d2dy200,  d4dx2y200,
       v01,  d2dx201,  d2dy201,  d4dx2y201,
       v10,  d2dx210,  d2dy210,  d4dx2y210,
       v11,  d2dx211,  d2dy211,  d4dx2y211,
       dxp, dyp-delta/binwy, binwx, binwy,
       valym, dvaldxym, dvaldyym
    );

    Real valyp, dvaldxyp, dvaldyyp;
    bicubic_interpolation(
       v00,  d2dx200,  d2dy200,  d4dx2y200,
       v01,  d2dx201,  d2dy201,  d4dx2y201,
       v10,  d2dx210,  d2dy210,  d4dx2y210,
       v11,  d2dx211,  d2dy211,  d4dx2y211,
       dxp, dyp+delta/binwy, binwx, binwy,
       valyp, dvaldxyp, dvaldyyp
    );

    //TR << "dvaldy analytic " << dvaldy << " and numeric " << (valyp-valym)/(2*delta) << std::endl;
    TS_ASSERT_DELTA( dvaldy, (valyp-valym)/(2*delta), 1e-6 );


  }

};
