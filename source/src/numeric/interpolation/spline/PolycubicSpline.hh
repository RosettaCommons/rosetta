// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


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
#include <numeric/interpolation/spline/PolycubicSpline.fwd.hh>
#include <numeric/interpolation/spline/PolycubicSplineBase.hh>
#include <numeric/interpolation/spline/Cubic_spline.fwd.hh> //For BorderFlag enum definition.
#include <numeric/MathNTensor.hh>

#include <utility>

namespace numeric {
namespace interpolation {
namespace spline {

template< Size N >
class PolycubicSpline : public PolycubicSplineBase
{
public:

	typedef PolycubicSplineBase parent;

	//////////////////////////////////
	// construction and destruction //
	//////////////////////////////////

	/// @brief construct generic PolycubicSpline
	PolycubicSpline():
		parent(N),
		border_( utility::fixedsizearray1< BorderFlag, ( N ) >( e_Natural ) ),
		start_( utility::fixedsizearray1< Real, ( N ) >( 0.0 ) ),
		delta_( utility::fixedsizearray1< Real, ( N ) >( 0.0 ) ),
		n_derivs_( utility::fixedsizearray1< MathNTensor< Real, N >, ( 1 << N ) >( MathNTensor< Real, N >() ) ),
		firstbe_( utility::fixedsizearray1< std::pair< Real, Real >, N >( std::pair< Real, Real >(0.0,0.0) ) ),
		LinCont_( utility::fixedsizearray1< bool, N >( false ) )
	{
		/*
		utility::fixedsizearray1< BorderFlag, N > border_;   ///< controls the behavior at x/y_0 and x/y_dim-1

		utility::fixedsizearray1< Real, N > start_;
		utility::fixedsizearray1< Real, N > delta_;    ///< gives the arguments as a sequence of equidistant points

		utility::fixedsizearray1< MathNTensor< Real, N >, ( 1 << N ) > n_derivs_; // has 000 = values_,
		// 001 = z deriv, etc.
		utility::fixedsizearray1< std::pair< Real, Real>, N > firstbe_; ///< first order derivative at x_0/dim-1, y_0/dim-1, z_0/dim-1 can be set for BorderFlag FIRSTDER

		utility::fixedsizearray1< bool, N > LinCont_;    ///< if the argument x is outside the range decide if the spline should be continued */
	}

	/// @brief copy constructor
	PolycubicSpline( PolycubicSpline const &src ) :
		parent( src.dimensionality() ),
		border_(src.border_),
		start_(src.start_),
		delta_(src.delta_),
		n_derivs_(src.n_derivs_),
		firstbe_(src.firstbe_),
		LinCont_(src.LinCont_)
	{}

	/////////////////
	// data access //
	/////////////////

	/// get the second order derivatives of the spline
	/// for 3 dimensions, you would pass 1-8 and get
	/// values, z, y, yz, x, xz, xy, xyz
	MathNTensor< Real, N > const & get_deriv( Size n ) const
	{
		return n_derivs_[ n ];
	}

	Size dimensionality() const { return N; }


	utility::fixedsizearray1< Real, ( 1 << N ) > get_all_derivs( utility::fixedsizearray1< Size, N > const & indices )
	{
		utility::fixedsizearray1< Real, ( 1 << N ) > ret;
		for ( Size i = 1; i <= ( 1 << N ); ++i ) ret[ i ] = n_derivs_[ i ]( indices );
		return ret;
	}

	////////////////
	// operations //
	////////////////

	/// @return value at (x1, x2, ... xn)
	Real F( utility::fixedsizearray1< Real, N > const & xs ) const;

	/// @return partial derivative at (x1, x2, ... xn) for var i
	Real dFdxi( Size n, utility::fixedsizearray1< Real, N > const & xs ) const;

	/// @return partial derivatives at (x1, x2, ... xn)
	utility::fixedsizearray1< Real, N > dFdall( utility::fixedsizearray1< Real, N > const & xs ) const;

	/// @return value and derivative at (x, y)
	//void FdF( utility::vector1< Real > xs, Real & val, utility::vector1< Real > & dvaldxs ) const;

	/// train PolycubicSpline
	void train
	(
		utility::fixedsizearray1< BorderFlag, N > const & BORDER,//[3],
		utility::fixedsizearray1< double, N > const & START,//[3],
		utility::fixedsizearray1< double, N > const & DELTA,//[3],
		MathNTensor< Real, N > const & RESULTS,
		utility::fixedsizearray1< bool, N > const & LINCONT,//[3],
		utility::fixedsizearray1< std::pair< Real, Real >, N > const & FIRSTBE//[3]
	);


private:
	utility::fixedsizearray1< BorderFlag, N > border_;   ///< controls the behavior at x/y_0 and x/y_dim-1

	utility::fixedsizearray1< Real, N > start_;
	utility::fixedsizearray1< Real, N > delta_;    ///< gives the arguments as a sequence of equidistant points

	utility::fixedsizearray1< MathNTensor< Real, N >, ( 1 << N ) > n_derivs_; // has 000 = values_,
	// xyz
	// 001 = z deriv, etc.

	utility::fixedsizearray1< std::pair< Real, Real>, N > firstbe_; ///< first order derivative at x_0/dim-1, y_0/dim-1, z_0/dim-1 can be set for BorderFlag FIRSTDER

	utility::fixedsizearray1< bool, N > LinCont_;    ///< if the argument x is outside the range decide if the spline should be continued linearly


};

// Stub implementations for 1 and 2 to prevent problems!
template<>
class PolycubicSpline< 1 > : public PolycubicSplineBase
{
public:

	typedef PolycubicSplineBase parent;

	//////////////////////////////////
	// construction and destruction //
	//////////////////////////////////

	/// @brief construct generic PolycubicSpline
	PolycubicSpline():
		parent(1),
		border_( utility::fixedsizearray1< BorderFlag, ( 1 ) >( e_Natural ) ),
		start_( utility::fixedsizearray1< Real, ( 1 ) >( 0.0 ) ),
		delta_( utility::fixedsizearray1< Real, ( 1 ) >( 0.0 ) ),
		n_derivs_( utility::fixedsizearray1< MathNTensor< Real, 1 >, 2 >( MathNTensor< Real, 1 >() ) ),
		firstbe_( utility::fixedsizearray1< std::pair< Real, Real >, 1 >( std::pair< Real, Real >(0.0,0.0) ) ),
		LinCont_( utility::fixedsizearray1< bool, 1 >( false ) )
	{}

	/// @brief copy constructor
	PolycubicSpline( PolycubicSpline const &src ) :
		parent( 1 ),
		border_(src.border_),
		start_(src.start_),
		delta_(src.delta_),
		n_derivs_(src.n_derivs_),
		firstbe_(src.firstbe_),
		LinCont_(src.LinCont_)
	{}

	/////////////////
	// data access //
	/////////////////

	/// @brief get the second order derivatives of the spline
	/// for 3 dimensions, you would pass 1-8 and get
	/// values, z, y, yz, x, xz, xy, xyz
	MathNTensor< Real, 1 > const & get_deriv( Size n ) const
	{
		return n_derivs_[ n ];
	}

	/// @brief Get the dimensionality.
	///
	Size dimensionality() const { return 1; }


	utility::fixedsizearray1< Real, 2 > get_all_derivs( utility::fixedsizearray1< Size, 1 > & /*indices*/ )
	{
		utility::fixedsizearray1< Real, 2 > ret;
		return ret;
	}

	////////////////
	// operations //
	////////////////

	/// @return value at (x1, x2, ... xn)
	Real F( utility::fixedsizearray1< Real, 1 > const & /*xs*/ ) const { return 0; }

	/// @return partial derivative at (x1, x2, ... xn) for var i
	Real dFdxi( Size , utility::fixedsizearray1< Real, 1 > const & /*xs*/ ) const { return 0; }

	/// @return partial derivatives at (x1, x2, ... xn)
	utility::fixedsizearray1< Real, 1 > dFdall( utility::fixedsizearray1< Real, 1 > const & /*xs*/ ) const { return 0; }

	void train
	(
		utility::fixedsizearray1< BorderFlag, 1 > const & ,//[3],
		utility::fixedsizearray1< double, 1 > const & ,//[3],
		utility::fixedsizearray1< double, 1 > const & ,//[3],
		MathNTensor< Real, 1 > const & ,
		utility::fixedsizearray1< bool, 1 > const & ,//[3],
		utility::fixedsizearray1< std::pair< Real, Real >, 1 > const & //[3]
	) {
		utility_exit_with_message( "Error in numeric/interpolation/spline/PolycubicSpline<1>::train(): This function has not yet been implemented!  PolycubicSplines do not currently support the 1D or 2D cases; for these, use Spline or BicubicSpline classes." );
	}


private:
	utility::fixedsizearray1< BorderFlag, 1 > border_;   ///< controls the behavior at x/y_0 and x/y_dim-1
	utility::fixedsizearray1< Real, 1 > start_;
	utility::fixedsizearray1< Real, 1 > delta_;    ///< gives the arguments as a sequence of equidistant points
	utility::fixedsizearray1< MathNTensor< Real, 1 >, 2 > n_derivs_; // has 000 = values_,
	utility::fixedsizearray1< std::pair< Real, Real>, 1 > firstbe_;
	utility::fixedsizearray1< bool, 1 > LinCont_;
};

template<>
class PolycubicSpline< 2 >  : public PolycubicSplineBase
{
public:

	typedef PolycubicSplineBase parent;

	//////////////////////////////////
	// construction and destruction //
	//////////////////////////////////

	/// construct generic PolycubicSpline
	PolycubicSpline():
		parent(2),
		border_( utility::fixedsizearray1< BorderFlag, ( 2 ) >( e_Natural ) ),
		start_( utility::fixedsizearray1< Real, ( 2 ) >( 0.0 ) ),
		delta_( utility::fixedsizearray1< Real, ( 2 ) >( 0.0 ) ),
		n_derivs_( utility::fixedsizearray1< MathNTensor< Real, 2 >, 4 >( MathNTensor< Real, 2 >() ) ),
		firstbe_( utility::fixedsizearray1< std::pair< Real, Real >, 2 >( std::pair< Real, Real >(0.0,0.0) ) ),
		LinCont_( utility::fixedsizearray1< bool, 2 >( false ) )
	{}


	/// @brief Get the dimensionality.
	///
	Size dimensionality() const { return 2; }

	/// @brief copy constructor
	PolycubicSpline( PolycubicSpline const &src ) :
		parent( 2 ),
		border_( src.border_ ),
		n_derivs_(src.n_derivs_),
		firstbe_(src.firstbe_),
		LinCont_(src.LinCont_)
	{}

	/////////////////
	// data access //
	/////////////////

	/// get the second order derivatives of the spline
	/// for 3 dimensions, you would pass 1-8 and get
	/// values, z, y, yz, x, xz, xy, xyz
	MathNTensor< Real, 2 > const & get_deriv( Size n ) const
	{
		return n_derivs_[ n ];
	}


	utility::fixedsizearray1< Real, 4 > get_all_derivs( utility::fixedsizearray1< Size, 2 > &  )
	{
		utility::fixedsizearray1< Real, 4 > ret;
		return ret;
	}

	////////////////
	// operations //
	////////////////

	/// @return value at (x1, x2, ... xn)
	Real F( utility::fixedsizearray1< Real, 2 > const &  ) const { return 0; }

	/// @return partial derivative at (x1, x2, ... xn) for var i
	Real dFdxi( Size , utility::fixedsizearray1< Real, 2 > const &  ) const { return 0; }

	/// @return partial derivatives at (x1, x2, ... xn)
	utility::fixedsizearray1< Real, 2 > dFdall( utility::fixedsizearray1< Real, 2 > const &  ) const { return 0; }

	void train
	(
		utility::fixedsizearray1< BorderFlag, 2 > const & ,//[3],
		utility::fixedsizearray1< double, 2 > const & ,//[3],
		utility::fixedsizearray1< double, 2 > const & ,//[3],
		MathNTensor< Real, 2 > & ,
		utility::fixedsizearray1< bool, 2 > const & ,//[3],
		utility::fixedsizearray1< std::pair< Real, Real >, 2 > const & //[3]
	) {
		utility_exit_with_message( "Error in numeric/interpolation/spline/PolycubicSpline<2>::train(): This function has not yet been implemented!  PolycubicSplines do not currently support the 1D or 2D cases; for these, use Spline or BicubicSpline classes." );
	}


private:
	utility::fixedsizearray1< BorderFlag, 2 > border_;   ///< controls the behavior at x/y_0 and x/y_dim-1
	utility::fixedsizearray1< Real, 2 > start_;
	utility::fixedsizearray1< Real, 2 > delta_;    ///< gives the arguments as a sequence of equidistant points
	utility::fixedsizearray1< MathNTensor< Real, 2 >, 4 > n_derivs_; // has 000 = values_,
	utility::fixedsizearray1< std::pair< Real, Real>, 2 > firstbe_;
	utility::fixedsizearray1< bool, 2 > LinCont_;
};

void hokey_template_workaround();

}//end namespace spline
}//end namespace interpolation
}//end namespace numeric


#endif /* POLYCUBIC_SPLINE_HH_ */
