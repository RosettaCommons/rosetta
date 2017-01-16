
//////////////////////////////////////////////////////////////////////
///
/// @brief
/// read the header file!
///
/// @references
/// Numerical Recipes in c++ 2nd edition
/// Ralf Mueller
///
///
/// @author Steven Combs, Ralf Mueller, Jens Meiler
/// ported to Rosetta by Andrew Leaver-Fay
/// generalized to N dimensions by Andrew Watkins
///
/////////////////////////////////////////////////////////////////////////

// Unit headers
#include <numeric/interpolation/spline/TricubicSpline.hh>
#include <numeric/interpolation/spline/PolycubicSpline.tmpl.hh>

// Package headers
#include <numeric/types.hh>
#include <numeric/interpolation/spline/CubicSpline.hh>
#include <numeric/interpolation/spline/BicubicSpline.hh>
#include <numeric/MathVector.hh>
#include <numeric/MathMatrix.hh>
#include <numeric/MathNTensor.hh>
#include <utility/fixedsizearray1.hh>
#include <utility/vector1.hh>

namespace numeric {
namespace interpolation {
namespace spline {

/// @brief Dummy function, never to be called.  This is only here to ensure that the compiler creates
/// PolycubicSpline<3> through PolycubicSpline<9> classes.
void hokey_template_workaround() {
	using numeric::interpolation::spline::BorderFlag;

	numeric::interpolation::spline::PolycubicSpline< 3 > pcs3;
	pcs3.train( utility::fixedsizearray1< BorderFlag, 3 >(),
		utility::fixedsizearray1< double, 3 >(),
		utility::fixedsizearray1< double, 3 >(),
		numeric::MathNTensor< Real, 3 >(),
		utility::fixedsizearray1< bool, 3 >(),
		utility::fixedsizearray1< std::pair< Real, Real >, 3 >( std::pair< Real, Real >( 0,0 ) ) );
	pcs3.F( utility::fixedsizearray1< Real, 3 >() );
	pcs3.dFdxi( 1, utility::fixedsizearray1< Real, 3 >() );

	numeric::interpolation::spline::PolycubicSpline< 4 > pcs4;
	pcs4.train( utility::fixedsizearray1< BorderFlag, 4 >(),
		utility::fixedsizearray1< double, 4 >(),
		utility::fixedsizearray1< double, 4 >(),
		numeric::MathNTensor< Real, 4 >(),
		utility::fixedsizearray1< bool, 4 >(),
		utility::fixedsizearray1< std::pair< Real, Real >, 4 >( std::pair< Real, Real >( 0,0 ) ) );
	pcs4.F( utility::fixedsizearray1< Real, 4 >() );
	pcs4.dFdxi( 1, utility::fixedsizearray1< Real, 4 >() );

	numeric::interpolation::spline::PolycubicSpline< 5 > pcs5;
	pcs5.train( utility::fixedsizearray1< BorderFlag, 5 >(),
		utility::fixedsizearray1< double, 5 >(),
		utility::fixedsizearray1< double, 5 >(),
		numeric::MathNTensor< Real, 5 >(),
		utility::fixedsizearray1< bool, 5 >(),
		utility::fixedsizearray1< std::pair< Real, Real >, 5 >( std::pair< Real, Real >( 0,0 ) ) );
	pcs5.F( utility::fixedsizearray1< Real, 5 >() );
	pcs5.dFdxi( 1, utility::fixedsizearray1< Real, 5 >() );

	numeric::interpolation::spline::PolycubicSpline< 6 > pcs6;
	pcs6.train( utility::fixedsizearray1< BorderFlag, 6 >(),
		utility::fixedsizearray1< double, 6 >(),
		utility::fixedsizearray1< double, 6 >(),
		numeric::MathNTensor< Real, 6 >(),
		utility::fixedsizearray1< bool, 6 >(),
		utility::fixedsizearray1< std::pair< Real, Real >, 6 >( std::pair< Real, Real >( 0,0 ) ) );
	pcs6.F( utility::fixedsizearray1< Real, 6 >() );
	pcs6.dFdxi( 1, utility::fixedsizearray1< Real, 6 >() );

	numeric::interpolation::spline::PolycubicSpline< 7 > pcs7;
	pcs7.train( utility::fixedsizearray1< BorderFlag, 7 >(),
		utility::fixedsizearray1< double, 7 >(),
		utility::fixedsizearray1< double, 7 >(),
		numeric::MathNTensor< Real, 7 >(),
		utility::fixedsizearray1< bool, 7 >(),
		utility::fixedsizearray1< std::pair< Real, Real >, 7 >( std::pair< Real, Real >( 0,0 ) ) );
	pcs7.F( utility::fixedsizearray1< Real, 7 >() );
	pcs7.dFdxi( 1, utility::fixedsizearray1< Real, 7 >() );

	numeric::interpolation::spline::PolycubicSpline< 8 > pcs8;
	pcs8.train( utility::fixedsizearray1< BorderFlag, 8 >(),
		utility::fixedsizearray1< double, 8 >(),
		utility::fixedsizearray1< double, 8 >(),
		numeric::MathNTensor< Real, 8 >(),
		utility::fixedsizearray1< bool, 8 >(),
		utility::fixedsizearray1< std::pair< Real, Real >, 8 >( std::pair< Real, Real >( 0,0 ) ) );
	pcs8.F( utility::fixedsizearray1< Real, 8 >() );
	pcs8.dFdxi( 1, utility::fixedsizearray1< Real, 8 >() );

	numeric::interpolation::spline::PolycubicSpline< 9 > pcs9;
	pcs9.train( utility::fixedsizearray1< BorderFlag, 9 >(),
		utility::fixedsizearray1< double, 9 >(),
		utility::fixedsizearray1< double, 9 >(),
		numeric::MathNTensor< Real, 9 >(),
		utility::fixedsizearray1< bool, 9 >(),
		utility::fixedsizearray1< std::pair< Real, Real >, 9 >( std::pair< Real, Real >( 0,0 ) ) );
	pcs9.F( utility::fixedsizearray1< Real, 9 >() );
	pcs9.dFdxi( 1, utility::fixedsizearray1< Real, 9 >() );

}

/// @brief Given a PolycubicSplineBase and a set of coordinates, call PolycubicSpline<N>::F and return the value.
/// @details Convenience function to hide the switch/case logic.  Only works for PolycubicSplines of dimensionality 3 through 9.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
Real
get_PolycubicSpline_F(
	PolycubicSplineBaseCOP splinebase,
	utility::vector1< Real > const &coords
) {
	runtime_assert( splinebase );
	Size const dim( splinebase->dimensionality() );
	switch( dim ) {
	case 3 :
		{
		PolycubicSplineCOP<3> spline( utility::pointer::dynamic_pointer_cast<PolycubicSpline<3> const>( splinebase ) );
		runtime_assert(spline);
		utility::fixedsizearray1< Real, 3 > fixedcoords;
		for ( Size i=1; i<=dim; ++i ) fixedcoords[i] = coords[i];
		return spline->F(fixedcoords);
	}
		break;
	case 4 :
		{
		PolycubicSplineCOP<4> spline( utility::pointer::dynamic_pointer_cast<PolycubicSpline<4> const>( splinebase ) );
		runtime_assert(spline);
		utility::fixedsizearray1< Real, 4 > fixedcoords;
		for ( Size i=1; i<=dim; ++i ) fixedcoords[i] = coords[i];
		return spline->F(fixedcoords);
	}
		break;
	case 5 :
		{
		PolycubicSplineCOP<5> spline( utility::pointer::dynamic_pointer_cast<PolycubicSpline<5> const>( splinebase ) );
		runtime_assert(spline);
		utility::fixedsizearray1< Real, 5 > fixedcoords;
		for ( Size i=1; i<=dim; ++i ) fixedcoords[i] = coords[i];
		return spline->F(fixedcoords);
	}
		break;
	case 6 :
		{
		PolycubicSplineCOP<6> spline( utility::pointer::dynamic_pointer_cast<PolycubicSpline<6> const>( splinebase ) );
		runtime_assert(spline);
		utility::fixedsizearray1< Real, 6 > fixedcoords;
		for ( Size i=1; i<=dim; ++i ) fixedcoords[i] = coords[i];
		return spline->F(fixedcoords);
	}
		break;
	case 7 :
		{
		PolycubicSplineCOP<7> spline( utility::pointer::dynamic_pointer_cast<PolycubicSpline<7> const>( splinebase ) );
		runtime_assert(spline);
		utility::fixedsizearray1< Real, 7 > fixedcoords;
		for ( Size i=1; i<=dim; ++i ) fixedcoords[i] = coords[i];
		return spline->F(fixedcoords);
	}
		break;
	case 8 :
		{
		PolycubicSplineCOP<8> spline( utility::pointer::dynamic_pointer_cast<PolycubicSpline<8> const>( splinebase ) );
		runtime_assert(spline);
		utility::fixedsizearray1< Real, 8 > fixedcoords;
		for ( Size i=1; i<=dim; ++i ) fixedcoords[i] = coords[i];
		return spline->F(fixedcoords);
	}
		break;
	case 9 :
		{
		PolycubicSplineCOP<9> spline( utility::pointer::dynamic_pointer_cast<PolycubicSpline<9> const>( splinebase ) );
		runtime_assert(spline);
		utility::fixedsizearray1< Real, 9 > fixedcoords;
		for ( Size i=1; i<=dim; ++i ) fixedcoords[i] = coords[i];
		return spline->F(fixedcoords);
	}
		break;
	default :
		utility_exit_with_message( "Error in numeric::interpolation::spline::get_PolycubicSpline_F(): The provided spline dimensionality must be between 3 and 9, inclusive." );
		break;
	}

	return 0.0; //This line is never reached, and is only here to keep compilers happy.
}

/// @brief Given a PolycubicSplineBase and a set of coordinates, call PolycubicSpline<N>::dFdall and return the value.
/// @details Convenience function to hide the switch/case logic.  Only works for PolycubicSplines of dimensionality 3 through 9.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
void
get_PolycubicSpline_gradient(
	PolycubicSplineBaseCOP splinebase,
	utility::vector1< Real > const &coords,
	utility::vector1< Real > &gradient_out
) {
	runtime_assert( splinebase );
	Size const dim( splinebase->dimensionality() );
	switch( dim ) {
	case 3 :
		{
		PolycubicSplineCOP<3> spline( utility::pointer::dynamic_pointer_cast<PolycubicSpline<3> const>( splinebase ) );
		runtime_assert(spline);
		utility::fixedsizearray1< Real, 3 > fixedcoords;
		for ( Size i=1; i<=dim; ++i ) fixedcoords[i] = coords[i];
		utility::fixedsizearray1< Real, 3 > const result ( spline->dFdall(fixedcoords) );
		gradient_out.resize(dim);
		for ( Size i=1; i<=dim; ++i ) gradient_out[i] = result[i];
	}
		break;
	case 4 :
		{
		PolycubicSplineCOP<4> spline( utility::pointer::dynamic_pointer_cast<PolycubicSpline<4> const>( splinebase ) );
		runtime_assert(spline);
		utility::fixedsizearray1< Real, 4 > fixedcoords;
		for ( Size i=1; i<=dim; ++i ) fixedcoords[i] = coords[i];
		utility::fixedsizearray1< Real, 4 > const result ( spline->dFdall(fixedcoords) );
		gradient_out.resize(dim);
		for ( Size i=1; i<=dim; ++i ) gradient_out[i] = result[i];
	}
		break;
	case 5 :
		{
		PolycubicSplineCOP<5> spline( utility::pointer::dynamic_pointer_cast<PolycubicSpline<5> const>( splinebase ) );
		runtime_assert(spline);
		utility::fixedsizearray1< Real, 5 > fixedcoords;
		for ( Size i=1; i<=dim; ++i ) fixedcoords[i] = coords[i];
		utility::fixedsizearray1< Real, 5 > const result ( spline->dFdall(fixedcoords) );
		gradient_out.resize(dim);
		for ( Size i=1; i<=dim; ++i ) gradient_out[i] = result[i];
	}
		break;
	case 6 :
		{
		PolycubicSplineCOP<6> spline( utility::pointer::dynamic_pointer_cast<PolycubicSpline<6> const>( splinebase ) );
		runtime_assert(spline);
		utility::fixedsizearray1< Real, 6 > fixedcoords;
		for ( Size i=1; i<=dim; ++i ) fixedcoords[i] = coords[i];
		utility::fixedsizearray1< Real, 6 > const result ( spline->dFdall(fixedcoords) );
		gradient_out.resize(dim);
		for ( Size i=1; i<=dim; ++i ) gradient_out[i] = result[i];
	}
		break;
	case 7 :
		{
		PolycubicSplineCOP<7> spline( utility::pointer::dynamic_pointer_cast<PolycubicSpline<7> const>( splinebase ) );
		runtime_assert(spline);
		utility::fixedsizearray1< Real, 7 > fixedcoords;
		for ( Size i=1; i<=dim; ++i ) fixedcoords[i] = coords[i];
		utility::fixedsizearray1< Real, 7 > const result ( spline->dFdall(fixedcoords) );
		gradient_out.resize(dim);
		for ( Size i=1; i<=dim; ++i ) gradient_out[i] = result[i];
	}
		break;
	case 8 :
		{
		PolycubicSplineCOP<8> spline( utility::pointer::dynamic_pointer_cast<PolycubicSpline<8> const>( splinebase ) );
		runtime_assert(spline);
		utility::fixedsizearray1< Real, 8 > fixedcoords;
		for ( Size i=1; i<=dim; ++i ) fixedcoords[i] = coords[i];
		utility::fixedsizearray1< Real, 8 > const result ( spline->dFdall(fixedcoords) );
		gradient_out.resize(dim);
		for ( Size i=1; i<=dim; ++i ) gradient_out[i] = result[i];
	}
		break;
	case 9 :
		{
		PolycubicSplineCOP<9> spline( utility::pointer::dynamic_pointer_cast<PolycubicSpline<9> const>( splinebase ) );
		runtime_assert(spline);
		utility::fixedsizearray1< Real, 9 > fixedcoords;
		for ( Size i=1; i<=dim; ++i ) fixedcoords[i] = coords[i];
		utility::fixedsizearray1< Real, 9 > const result ( spline->dFdall(fixedcoords) );
		gradient_out.resize(dim);
		for ( Size i=1; i<=dim; ++i ) gradient_out[i] = result[i];
	}
		break;
	default :
		utility_exit_with_message( "Error in numeric::interpolation::spline::get_PolycubicSpline_gradient(): The provided spline dimensionality must be between 3 and 9, inclusive." );
		break;
	}
}



}//end namespace spline
}//end namespace interpolation
}//end namespace numeric

