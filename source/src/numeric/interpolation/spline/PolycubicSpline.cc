// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
#include <numeric/interpolation/spline/PolycubicSpline.hh>

// Package headers
#include <numeric/types.hh>
#include <numeric/interpolation/spline/Cubic_spline.hh>
#include <numeric/interpolation/spline/Bicubic_spline.hh>
#include <numeric/MathVector.hh>
#include <numeric/MathMatrix.hh>
#include <numeric/MathNTensor.hh>

namespace numeric {
namespace interpolation {
namespace spline {

#ifdef WIN32
inline double pow(Size x, Size y)
{
	return pow( (double)x, y);
}
#endif

//! train PolycubicSpline
void PolycubicSpline::train
(
	const utility::vector1< BorderFlag > & BORDER,//[3],
	const utility::vector1< double > & START,//[3],
	const utility::vector1< double > & DELTA,//[3],
	const MathNTensor< Real > & RESULTS,
	const utility::vector1< bool > & LINCONT,//[3],
	const utility::vector1< std::pair< Real, Real > > & FIRSTBE//[3]
)
{
	// this does not work for the base bicubic case
	assert( n_xs_ > 2 );

	//check, if the points are given in positive direction
	for ( Size i = 1; i <= n_xs_; ++i ) {
		assert( DELTA[ i ] > 0 );
	}

	//determine values for all n dimensions
	utility::vector1< Size > dims;
	for ( Size i = 1; i <= n_xs_; ++i ) {
		dims.push_back( RESULTS.n_dimensions( i ) );
	}

	//assigning values
	for ( Size i = 1; i <= n_xs_; ++i ) {
		border_[ i ]  = BORDER[ i ];
		start_[ i ]  = START[ i ];
		delta_[ i ]  = DELTA[ i ];
		LinCont_[ i ] = LINCONT[ i ];
		firstbe_[ i ] = FIRSTBE[ i ];
	}

	Size number_of_derivs = 1 << n_xs_;

	for ( Size i = 1; i <= number_of_derivs; ++i ) {
		n_derivs_.push_back( RESULTS );
	}

	//train seven times for fxx, fyy, fzz, fxxyy, fxxzz, fyyzz, fxxyyzz
	//reduction to Spline(n-1)D by training only (n-1)D-layers at the same time

	if ( n_xs_ == 3 ) { // base case must be treated separately because of mathmatrices
		// loop over first index
		for ( Size first_index = 0; first_index < dims[ 1 ]; ++first_index ) {
			MathMatrix< Real> values( dims[ n_xs_ - 1 ], dims[ n_xs_ ] );
			for ( Size im1 = 1; im1 <= dims[ n_xs_ - 1 ]; ++im1 ) {
				for ( Size i = 1; i <= dims[ n_xs_ ]; ++i ) {
					utility::vector1< Size > position( n_xs_, 0 );
					position[ 1 ] = first_index + 1;
					position[ n_xs_-1 ] = im1;
					position[ n_xs_   ] = i;
					values( im1-1, i-1 )= n_derivs_[ 1 ]( position ); // because values is a mathmatrix we have to use numbers one lesser
				}
			}
			BicubicSpline bs;
			BorderFlag border[2]               = { BORDER[ n_xs_-1 ],  BORDER[ n_xs_ ]  };
			Real start[2]                      = { START[ n_xs_-1 ],   START[ n_xs_ ]   };
			Real delta[2]                      = { DELTA[ n_xs_-1 ],   DELTA[ n_xs_ ]   };
			bool lin_cont[2]                   = { LINCONT[ n_xs_-1 ], LINCONT[ n_xs_ ] };
			std::pair< Real, Real > firstbe[2] = { FIRSTBE[ n_xs_-1 ], FIRSTBE[ n_xs_ ] };

			bs.train( border, start, delta, values, lin_cont, firstbe);
			n_derivs_[3].replace_layer( first_index, bs.get_dsecox());
			n_derivs_[2].replace_layer( first_index, bs.get_dsecoy());
			n_derivs_[4].replace_layer( first_index, bs.get_dsecoxy());
		}

		Size const half_n_derivs = n_derivs_.size() / 2;

		// loop over all but 1 indices (for base case: two indices)
		utility::vector1< Size > indices2( n_xs_, 0 );
		Size p = 1;
		while ( indices2[ n_xs_ ] == 0 ) {

			utility::vector1< Size > position( n_xs_, 1 ); // to be "i, indices2"
			for ( Size k = 1; k <= n_xs_ - 1; k++ ) position[ k + 1 ] = indices2[ k ] + 1;

			// This holds the 4 MathVectors of that have no derivative in the most significant dimension (e.g. x in xyz)
			utility::vector1< MathVector< Real > > n_derivs( half_n_derivs, MathVector< Real > ( dims[ n_xs_ ] ) );

			for ( Size i = 0; i < dims[ 1 ]; ++i ) {
				position[ 1 ] = i + 1;
				// amw all the n_derivs that have no x component, i.e. the FIRST HALF of n_derivs_
				for ( Size ii = 1; ii <= half_n_derivs; ++ii ) {
					n_derivs[ ii ]( i ) = n_derivs_[ ii ](  position );
				}
			}

			utility::vector1< CubicSpline > cs( half_n_derivs );
			for ( Size ii = 1; ii <= half_n_derivs; ++ii ) {
				cs[ ii ].train(   BORDER[ 1], START[1], DELTA[ 1],  n_derivs[ ii ], FIRSTBE[ 1]);
			}

			for ( Size i = 0; i < dims[ n_xs_ ]; ++i ) {

				position[ 1 ] = i + 1;
				for ( Size ii = 1; ii <= half_n_derivs; ++ii ) {
					n_derivs_[ half_n_derivs + ii ]( position ) = cs[ ii ].get_dsecox()( i );
				}
			}

			indices2[ 1 ]++;
			while ( indices2[ p ] == dims[ p ] ) {
				indices2[ p ] = 0;
				indices2[ ++p ]++;
				if ( indices2[ p ] != dims[ p ] ) p = 1;
			}
		}
	} else { // reduction of dimension, but just by one.
		utility::vector1< Real > smalldims;
		for ( Size dimi = 1; dimi <= n_xs_-1; dimi++ ) {
			smalldims.push_back( dims[ dimi ] );
		}
		MathNTensor< Real> values( smalldims );

		utility::vector1< Size > indices( n_xs_, 0 );
		Size p = 1;
		for ( Size layer = 0; layer < dims[ 1 ]; ++layer ) {
			while ( indices[ n_xs_ ] == 0 ) {
				utility::vector1< Size > position( n_xs_, 0 );
				for ( Size k = 2; k <= n_xs_; k++ ) position[ k ] = indices[ k ] + 1;
				position[ 1 ] = layer+1;

				values( smalldims )= n_derivs_[ 1 ]( position ); // because values is a mathmatrix we have to use numbers one lesser
				PolycubicSpline bs;
				utility::vector1< BorderFlag > border( n_xs_ - 1 );
				utility::vector1< Real > start( n_xs_ - 1 );
				utility::vector1< Real > delta( n_xs_ - 1 );
				utility::vector1< bool > lin_cont( n_xs_ - 1 );
				utility::vector1< std::pair< Real, Real > > firstbe( n_xs_ - 1 );
				for ( Size initi = 2; initi <= n_xs_; ++initi ) {
					border[ initi ] = BORDER[ initi ];
					start[ initi ] = START[ initi ];
					delta[ initi ] = DELTA[ initi ];
					lin_cont[ initi ] = LINCONT[ initi ];
					firstbe[ initi ] = FIRSTBE[ initi ];
				}

				// amw std::cout << "Going to train bicubic spline on " << border[0] << " " << border[1] << " start " << start[0] << " " << start[1] << " delta "<<delta[0]<<" " << delta[1] << " lin_cont " << lin_cont[0] << " " << lin_cont[ 1] << " first_be " << firstbe[0].first << ", " << firstbe[0].second << " " << firstbe[1].first << ", " << firstbe[1].second << std::endl;
				bs.train( border, start, delta, values, lin_cont, firstbe);

				for ( Size di = 2; di <= n_derivs_.size() / 2; ++di ) {
					n_derivs_[ di ].replace_layer( layer, bs.get_deriv( di ) );
				}

				indices[ p ]++;
				while ( indices[ n_xs_ ] == dims[ n_xs_ ] ) {
					indices[ p ] = 0;
					indices[ ++p ]++;
					if ( indices[ n_xs_ ] != dims[ n_xs_ ] ) p = 1;
				}
			}

			utility::vector1< Size > indices2( n_xs_ - 1, 0 );
			p = 1;
			while ( indices2[ n_xs_ - 1 ] == 0 ) {

				utility::vector1< MathVector< Real > > n_derivs;
				for ( Size ii = 1; ii <= n_derivs_.size() / 2; ++ii ) {
					n_derivs.push_back( MathVector< Real > ( dims[ n_xs_ ] ) );
				}

				//( n_derivs_.size()/2, dims[ n_xs_ ] ), dsecoy( dims[ n_xs_ ] ), dsecoz( dims[ n_xs_ ] ), dsecoyz( dims[ n_xs_ ] );
				for ( Size i( 0); i < dims[ n_xs_ ]; ++i ) {
					utility::vector1< Size > position( n_xs_, 0 );
					for ( Size k = 2; k <= n_xs_; k++ ) {
						position[ k ] = indices2[ k ] + 1;
					}
					position[ 1 ] = i+1;
					// amw all the n_derivs that have no x component, i.e. the FIRST HALF of n_derivs_

					for ( Size ii = 1; ii <= n_derivs_.size() / 2; ++ii ) {
						//n_derivs_[ ii ](   indices[n_xs_]) = n_derivs_[ ii ](  position);
						n_derivs[ ii ]( indices2[ n_xs_ ] ) = n_derivs_[ ii ](  position);
					}
				}

				utility::vector1< CubicSpline > cs( n_derivs_.size() / 2 );
				for ( Size ii = 1; ii <= n_derivs_.size() / 2; ++ii ) {
					cs[ ii ].train( BORDER[ 1], START[1], DELTA[ 1],  n_derivs[ ii ], FIRSTBE[ 1]);
				}

				for ( Size i( 0); i < dims[ n_xs_ ]; ++i ) {
					utility::vector1< Size > position( n_xs_, 0 );
					for ( Size k = 2; k <= n_xs_; k++ ) {
						position[ k ] = indices2[ k ] + 1;
					}
					position[ 1 ] = i+1;

					// amw I need to figure out what these indices should become!
					for ( Size ii = 1; ii <= n_derivs_.size() / 2; ++ii ) {
						n_derivs_[n_derivs_.size()/2 + ii]( position ) = cs[ ii ].get_dsecox()(    i);
					}
				}

				indices2[ 1 ]++;
				while ( indices2[ p ] == dims[ p ] ) {
					indices2[ p ] = 0;
					indices2[ ++p ]++;
					if ( indices2[ p ] != dims[ p ] ) p = 1;
				}
			}
		}
	}
	return;
}


//! return value at certain (x, y, z)
double PolycubicSpline::F( const utility::vector1 < double > xs ) const
{
	utility::vector1< Size > dims;
	for ( Size i = 1; i <= n_xs_; ++i ) {
		dims.push_back( n_derivs_[ 1 ].n_dimensions( i ) );
	}

	//check if argument is in range for non-periodic splines
	for ( Size i = 0; i < n_xs_; ++i ) {
		if ( ( border_[ i ] != e_Periodic)
				&& ( ( xs[ i+1 ] < start_[ i ]) || ( start_[ i ] + ( dims[ i+1 ]-1) * delta_[ i ] < xs[ i+1 ])) ) {
			//BCL_Assert( LinCont_[ 0], "argument out of range for non-periodic spline!");
			//BCL_Message( util::Message::e_Debug, "argument out of range, using linear continuation");
			if ( xs[ i ] < start_[ i ] ) {
				utility::vector1< int > pos( n_xs_, 0 );
				for ( Size si = 1; si <= n_xs_; ++si ) {
					pos[ si ] = ( si == i+1 ) ? start_[ si ] : xs[ si ];
				}
				return F( pos ) + ( xs[ i ] - start_[ i ] ) * dFdxi( i, pos );
			}
			if ( xs[ i ] > start_[ i ] + ( dims[ i ] - 1) * delta_[ i ] ) {
				utility::vector1< int > pos( n_xs_, 0 );
				for ( Size si = 1; si <= n_xs_; ++si ) {
					pos[ si ] = ( si == i+1 ) ? start_[ i ] + ( dims[ si ] - 1 ) * delta_[ i ] : xs[ si ];
				}
				return F( pos ) + ( xs[ i ] - start_[ i ] - ( dims[ i ] - 1 ) * delta_[ i ] ) * dFdxi( i, pos );
			}
		}
	}

	//determine i with start_[ 0]+(i-1)*delta_[ 0] < x < start_[ 0]+i*delta_[ 0] for the correct supporting points
	utility::vector1< int > ijk;
	for ( Size i = 0; i < n_xs_; ++i ) {
		ijk.push_back( int( floor( ( xs[ i+1 ]-start_[ i ]) / delta_[ i ] ) ) );
		while ( start_[ i ] + ijk[ i+1 ] * delta_[ i] < xs[ i+1 ] ) ijk[ i+1 ]++;
	}

	// this formula was derived from the Spline2D formula
	// the idea is to 'combine every part of the 2D formula
	// with dzm, dzp, dz3m, dz3p'

	utility::vector1< double > delta_akts;
	utility::vector1< double > dps;
	utility::vector1< double > dms;
	utility::vector1< double > d3ps;
	utility::vector1< double > d3ms;

	for ( Size i = 0; i < n_xs_; ++i ) {
		delta_akts.push_back( xs[ i+1 ] - start_[ i ] - (ijk[ i+1 ] - 1) * delta_[ i ] );
		dps.push_back( delta_akts[ i+1 ] / delta_[ i ] );
		dms.push_back( 1 - dps[ i+1 ] );
		d3ps.push_back( ( dps[ i+1 ] * dps[ i+1 ] * dps[ i+1 ] - dps[ i+1 ] ) * sqr( delta_[ i ] ) / 6 );
		d3ms.push_back( ( dms[ i+1 ] * dms[ i+1 ] * dms[ i+1 ] - dms[ i+1 ] ) * sqr( delta_[ i ] ) / 6 );
		while ( ijk[ i+1 ] < 1 ) ijk[ i+1 ] += dims[ i+1 ];
	}

	double value = 0;
	Size n_points = static_cast< int >( pow( 2, n_xs_ ) );
	for ( Size posi = 1; posi <= n_points; ++posi ) {
		// if posi - 1 &
		double coeff = 1;
		utility::vector1< int > position( n_xs_, 0 );
		for ( Size iid = 1; iid <= n_points; ++iid ) {
			for ( Size xsi = 1; xsi <= n_xs_; ++xsi ) {
				if ( ( iid - 1 ) & static_cast< int >( pow( 2, xsi - 1 ) ) ) {
					coeff *= ( ( posi - 1 ) & static_cast< int >( pow( 2, xsi - 1 ) ) ) ? d3ps[ xsi ] : d3ms[ xsi ];
				} else {
					coeff *= ( ( posi - 1 ) & static_cast< int >( pow( 2, xsi - 1 ) ) ) ? dps[ xsi ] : dms[ xsi ];
				}
				position.push_back( ( ( posi - 1 ) & static_cast< int >( pow( 2, xsi ) ) ) ? ( ijk[ xsi ] - 1 % dims[ xsi ] ) : ( ijk[ xsi ] - 1 % dims[ xsi ] ) );
			}
			value += coeff * n_derivs_[ iid ]( position );

		}
	}

	return value;
}

//! return partial derivative at certain (x, y, z) for x
double PolycubicSpline::dFdxi( Size n, const utility::vector1 < double > xs ) const
{
	utility::vector1< Size > dims;
	for ( Size i = 1; i <= n_xs_; ++i ) {
		dims.push_back( n_derivs_[ 1 ].n_dimensions( i ) );
	}

	//check if argument is in range for non-periodic splines
	for ( Size i = 0; i < n_xs_; ++i ) {
		if ( ( border_[ i ] != e_Periodic)
				&& ( ( xs[ i+1 ] < start_[ i ]) || ( start_[ i ] + ( dims[ i+1 ]-1) * delta_[ i ] < xs[ i+1 ])) ) {
			//BCL_Assert( LinCont_[ 0], "argument out of range for non-periodic spline!");
			//BCL_Message( util::Message::e_Debug, "argument out of range, using linear continuation");
			if ( xs[ i ] < start_[ i ] ) {
				utility::vector1< int > pos( n_xs_, 0 );
				for ( Size si = 1; si <= n_xs_; ++si ) {
					pos[ si ] = ( si == i+1 ) ? start_[ si ] : xs[ si ];
				}
				return F( pos ) + ( xs[ i ] - start_[ i ] )*dFdxi( n, pos );
			}
			if ( xs[ i ] > start_[ i ] + ( dims[ i ] - 1) * delta_[ i ] ) {
				utility::vector1< int > pos( n_xs_, 0 );
				for ( Size si = 1; si <= n_xs_; ++si ) {
					pos[ si ] = ( si == i+1 ) ? start_[ i ] + ( dims[ si ] - 1 ) * delta_[ i ] : xs[ si ];
				}
				return F( pos ) + ( xs[ i ] - start_[ i ] - ( dims[ i ] - 1 ) * delta_[ i ] ) * dFdxi( n, pos );
			}
		}
	}

	//determine i with start_+(i-1)*delta_ < x < start_+i*delta_ for the correct supporting points
	//determine i with start_[ 0]+(i-1)*delta_[ 0] < x < start_[ 0]+i*delta_[ 0] for the correct supporting points
	utility::vector1< int > ijk;
	for ( Size i = 0; i < n_xs_; ++i ) {
		ijk.push_back( int( floor( ( xs[ i+1 ]-start_[ i ]) / delta_[ i ] ) ) );
		while ( start_[ i ] + ijk[i+1] * delta_[ i ] < xs[ i+1 ] ) ijk[ i+1 ]++;
	}

	utility::vector1< double > delta_akts;
	utility::vector1< double > dps;
	utility::vector1< double > dms;
	utility::vector1< double > d3ps;
	utility::vector1< double > d3ms;

	for ( Size i = 0; i < n_xs_; ++i ) {
		delta_akts.push_back( xs[ i+1 ] - start_[ i ] - (ijk[ i+1 ] - 1) * delta_[ i ] );
		dps.push_back( delta_akts[ i+1 ] / delta_[ i ] );
		dms.push_back( 1 - dps[ i+1 ] );
		d3ps.push_back( ( dps[ i+1 ] * dps[ i+1 ] * dps[ i+1 ] - dps[ i+1 ] ) * sqr( delta_[ i ] ) / 6 );
		d3ms.push_back( ( dms[ i+1 ] * dms[ i+1 ] * dms[ i+1 ] - dms[ i+1 ] ) * sqr( delta_[ i ] ) / 6 );
		while ( ijk[ i+1 ] < 1 ) ijk[ i+1 ] += dims[ i+1 ];
	}

	// We are finding the derivative d xn. We never should find ourselves using d3ps[ n ] or d3ms[ n ]

	double value = 0;
	Size n_points = static_cast< int >( pow( 2, n_xs_ ) );
	for ( Size posi = 1; posi <= n_points; ++posi ) {
		// if posi - 1 &
		double coeff = 1;
		utility::vector1< Size > position( n_xs_, 0 );
		for ( Size iid = 1; iid <= n_points; ++iid ) {
			for ( Size xsi = 1; xsi <= n_xs_; ++xsi ) {
				if ( xsi != n ) {
					if ( ( iid - 1 ) & static_cast< int >( pow( 2, xsi - 1 ) ) ) {
						coeff *= ( ( posi - 1 ) & static_cast< int >( pow( 2, xsi - 1 ) ) ) ? d3ps[ xsi ] : d3ms[ xsi ];
					} else {
						coeff *= ( ( posi - 1 ) & static_cast< int >( pow( 2, xsi - 1 ) ) ) ? dps[ xsi ] : dms[ xsi ];
					}
				}
				if ( ( posi - 1 ) & static_cast< int >( pow( 2, xsi - 1 ) ) && xsi == n ) {
					coeff *= -1;
				}
				position.push_back( ( ( posi - 1 ) & static_cast< int >( pow( 2, xsi ) ) ) ? ( ijk[ xsi ] - 1 % dims[ xsi ] ) : ( ijk[ xsi ] - 1 % dims[ xsi ] ) );
			}
			if ( ( iid - 1 ) & static_cast< int >( pow( 2, n - 1 ) ) ) {
				value += coeff * n_derivs_[ iid ]( position ) * (3 * dms[n]*dms[n] - 1) * delta_[ n-1] / 6;
			} else {
				value += coeff * n_derivs_[ iid ]( position ) / delta_[ n-1 ];
			}
		}
	}

	return value;
}

utility::vector1< Real > PolycubicSpline::dFdall( utility::vector1< Real > xs ) const
{
	utility::vector1< Real > values;
	for ( Size n = 1; n <= n_xs_; ++n ) {
		values.push_back( dFdxi( n, xs ) );
	}
	return values;
}

}//end namespace spline
}//end namespace interpolation
}//end namespace numeric

