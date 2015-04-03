// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


//////////////////////////////////////////////////////////////////////
///
/// @brief
/// Polycubic spline for smoothly interpolating a function in n dimensions
///
/// @details
///
///
/// @references
/// Numerical Recipes in c++ 2nd edition
/// Ralf Mueller
///
///
/// @author Steven Combs, Ralf Mueller, Jens Meiler
/// ported to Rosetta by Andrew Leaver-Fay
/// generalized to n dimensions by Andrew Watkins
///
/////////////////////////////////////////////////////////////////////////

#ifndef INCLUDED_numeric_interpolation_spline_PolycubicSpline_hh
#define INCLUDED_numeric_interpolation_spline_PolycubicSpline_hh

#include <numeric/types.hh>
#include <numeric/interpolation/spline/Cubic_spline.fwd.hh>
#include <numeric/MathTensor.hh>
#include <numeric/MathNTensor.hh>
#include <numeric/MathMatrix.hh>
#include <numeric/MathVector.hh>

#include <utility>

namespace numeric {
namespace interpolation {
namespace spline {

class PolycubicSpline
{
public:

	//////////////////////////////////
	// construction and destruction //
	//////////////////////////////////

	/// construct generic PolycubicSpline
	PolycubicSpline(){n_xs_ = 3;}
    PolycubicSpline(Size n_xs) : n_xs_( n_xs ) {
        border_.resize( n_xs );
        start_.resize( n_xs );
        delta_.resize( n_xs );
        //n_derivs_.resize( static_cast< Size > ( std::pow( 2, n_xs ) ) );
        firstbe_.resize( n_xs );
        LinCont_.resize( n_xs );
    }

	/// copy constructor
	PolycubicSpline* Clone() const
	{
	  return new PolycubicSpline( *this);
	}

	/////////////////
	// data access //
	/////////////////

    /// get the second order derivatives of the spline
    /// for 3 dimensions, you would pass 1-8 and get
    /// values, z, y, yz, x, xz, xy, xyz
    MathNTensor< Real> const & get_deriv( Size n ) const
	{
        return n_derivs_[ n ];
	}
    
    
    utility::vector1< Real > get_all_derivs( utility::vector1< Size > indices )
	{
        utility::vector1< Real > ret;
        for ( Size i = 1; i <= n_derivs_.size(); ++i ) ret.push_back( n_derivs_[ i ]( indices ) );
        return ret;
	}

	////////////////
	// operations //
	////////////////

	/// @return value at (x1, x2, ... xn)
	Real F( utility::vector1< Real > xs ) const;

	/// @return partial derivative at (x1, x2, ... xn) for var i
	Real dFdxi( Size n, utility::vector1< Real > xs ) const;

	/// @return partial derivatives at (x1, x2, ... xn)
	utility::vector1< Real > dFdall( utility::vector1< Real > xs ) const;

	/// @return value and derivative at (x, y)
	//void FdF( utility::vector1< Real > xs, Real & val, utility::vector1< Real > & dvaldxs ) const;

	/// train PolycubicSpline
	void train
	(
        const utility::vector1< BorderFlag > BORDER,//[3],
		const utility::vector1< double > START,//[3],
		const utility::vector1< double > DELTA,//[3],
		const MathNTensor< Real > &RESULTS,
		const utility::vector1< bool > LINCONT,//[3],
		const utility::vector1< std::pair< Real, Real > > FIRSTBE//[3]
	);


private:
    Size n_xs_; ///< number of dimensions
    utility::vector1< BorderFlag > border_;   ///< controls the behavior at x/y_0 and x/y_dim-1

	utility::vector1< Real > start_;
	utility::vector1< Real > delta_;    ///< gives the arguments as a sequence of equidistant points

    utility::vector1< MathNTensor< Real> > n_derivs_; // has 000 = values_,
                                                     // 001 = z deriv, etc.
	/*MathTensor< Real> values_;     ///< f(x,y,z)
	MathTensor< Real> dsecox_;     ///< second order derivatives for x -- d**2/dx**2 f(x,y,z)
	MathTensor< Real> dsecoy_;     ///< second order derivatives for y
	MathTensor< Real> dsecoxy_;    ///< second order derivatives for x and y
	MathTensor< Real> dsecoz_;     ///< second order derivatives for z
	MathTensor< Real> dsecoxz_;    ///< second order derivatives for xz
	MathTensor< Real> dsecoyz_;    ///< second order derivatives for yz
	MathTensor< Real> dsecoxyz_;   ///< second order derivatives for x y and z*/

	utility::vector1< std::pair< Real, Real> > firstbe_; ///< first order derivative at x_0/dim-1, y_0/dim-1, z_0/dim-1 can be set for BorderFlag FIRSTDER

	utility::vector1< bool > LinCont_;    ///< if the argument x is outside the range decide if the spline should be continued linearly


};

}//end namespace spline
}//end namespace interpolation
}//end namespace numeric


#endif /* POLYCUBIC_SPLINE_HH_ */
