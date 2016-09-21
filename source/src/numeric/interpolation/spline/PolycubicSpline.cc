
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
#include <numeric/interpolation/spline/Cubic_spline.hh>
#include <numeric/interpolation/spline/Bicubic_spline.hh>
#include <numeric/MathVector.hh>
#include <numeric/MathMatrix.hh>
#include <numeric/MathNTensor.hh>
#include <utility/fixedsizearray1.hh>

namespace numeric {
namespace interpolation {
namespace spline {

#ifdef WIN32
inline double pow(Size x, Size y)
{
	return pow( (double)x, y);
}
#endif

void hokey_template_workaround() {
	using numeric::interpolation::spline::BorderFlag;

	/*
	Use cubic and bicubic for these.
	numeric::interpolation::spline::PolycubicSpline< 1 > pcs1;
	pcs1.train( utility::fixedsizearray1< BorderFlag, 1 >(),
	utility::fixedsizearray1< double, 1 >(),
	utility::fixedsizearray1< double, 1 >(),
	numeric::MathNTensor< Real, 1 >(),
	utility::fixedsizearray1< bool, 1 >(),
	utility::fixedsizearray1< std::pair< Real, Real >, 1 >( std::pair< Real, Real >( 0,0 ) ) );

	numeric::interpolation::spline::PolycubicSpline< 2 > pcs2;
	pcs2.train( utility::fixedsizearray1< BorderFlag, 2 >(),
	utility::fixedsizearray1< double, 2 >(),
	utility::fixedsizearray1< double, 2 >(),
	numeric::MathNTensor< Real, 2 >(),
	utility::fixedsizearray1< bool, 2 >(),
	utility::fixedsizearray1< std::pair< Real, Real >, 2 >( std::pair< Real, Real >( 0,0 ) ) );
	*/

	//Not strictly necessary... except maybe?

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

	//numeric::interpolation::spline::PolycubicSpline< 7 > pcs7;
	//numeric::interpolation::spline::PolycubicSpline< 8 > pcs8;
	//numeric::interpolation::spline::PolycubicSpline< 9 > pcs9;
	//numeric::interpolation::spline::PolycubicSpline< 10 > pcs10;
	//numeric::interpolation::spline::PolycubicSpline< 11 > pcs11;
	//numeric::interpolation::spline::PolycubicSpline< 12 > pcs12;
}

}//end namespace spline
}//end namespace interpolation
}//end namespace numeric

