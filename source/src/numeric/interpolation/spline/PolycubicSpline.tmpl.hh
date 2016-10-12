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

#ifndef INCLUDED_numeric_interpolation_spline_PolycubicSpline_tmpl_hh
#define INCLUDED_numeric_interpolation_spline_PolycubicSpline_tmpl_hh

// Unit headers
#include <numeric/interpolation/spline/TricubicSpline.hh>
#include <numeric/interpolation/spline/PolycubicSpline.hh>

// Package headers
#include <numeric/types.hh>
#include <numeric/interpolation/spline/CubicSpline.hh>
#include <numeric/interpolation/spline/BicubicSpline.hh>
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

//! train PolycubicSpline
template< Size N >
void PolycubicSpline< N >::train
(
	utility::fixedsizearray1< BorderFlag, N > const & BORDER,//[3],
	utility::fixedsizearray1< double, N > const & START,//[3],
	utility::fixedsizearray1< double, N > const & DELTA,//[3],
	MathNTensor< Real, N > const & RESULTS,
	utility::fixedsizearray1< bool, N > const & LINCONT,//[3],
	utility::fixedsizearray1< std::pair< Real, Real >, N > const & FIRSTBE//[3]
) {
	// this does not work for the base bicubic case
	assert( N > 2 );

	//check, if the points are given in positive direction
	for ( Size i = 1; i <= N; ++i ) {
		assert( DELTA[ i ] > 0 );
	}

	//determine values for all n dimensions
	utility::fixedsizearray1< Size, N > dims;
	for ( Size i = 1; i <= N; ++i ) {
		dims[ i ] = RESULTS.n_bins( i );
	}

	//assigning values
	for ( Size i = 1; i <= N; ++i ) {
		border_[ i ]  = BORDER[ i ];
		start_[ i ]  = START[ i ];
		delta_[ i ]  = DELTA[ i ];
		LinCont_[ i ] = LINCONT[ i ];
		firstbe_[ i ] = FIRSTBE[ i ];
	}

	Size number_of_derivs = 1 << N;

	for ( Size i = 1; i <= number_of_derivs; ++i ) {
		n_derivs_[ i ] = RESULTS;
	}

	//train seven times for fxx, fyy, fzz, fxxyy, fxxzz, fyyzz, fxxyyzz
	//reduction to Spline(n-1)D by training only (n-1)D-layers at the same time

	if ( N == 3 ) { // base case must be treated separately because of mathmatrices
		// loop over first index
		for ( Size layer = 0; layer < dims[ 1 ]; ++layer ) {
			MathMatrix< Real> values( dims[ 2 ], dims[ 3 ] );
			for ( Size row = 0; row < dims[ 2 ]; ++row ) {
				for ( Size col = 0; col < dims[ 3]; ++col ) {
					utility::vector1< Size > position( 3, 0 );
					position[ 1 ] = layer;// + 1;
					position[ 2 ] = row;
					position[ 3 ] = col;
					values( row, col )= n_derivs_[ 1 ]( position );
				}
			}
			BicubicSpline bs;
			BorderFlag border[2]               = { BORDER[ 2 ],  BORDER[ 3 ]  };
			Real start[2]                      = { START[ 2 ],   START[ 3 ]   };
			Real delta[2]                      = { DELTA[ 2 ],   DELTA[ 3 ]   };
			bool lin_cont[2]                   = { LINCONT[ 2 ], LINCONT[ 3 ] };
			std::pair< Real, Real > firstbe[2] = { FIRSTBE[ 2 ], FIRSTBE[ 3 ] };

			bs.train( border, start, delta, values, lin_cont, firstbe);
			n_derivs_[3].replace_layer( layer, bs.get_dsecox());
			n_derivs_[2].replace_layer( layer, bs.get_dsecoy());
			n_derivs_[4].replace_layer( layer, bs.get_dsecoxy());
		}

		for ( Size row = 0; row < dims[ 2 ]; ++row ) {
			for ( Size col = 0; col < dims[ 3 ]; ++col ) {
				MathVector< Real > values( dims[ 1 ] ), dsecoz( dims[ 1 ] ), dsecoy( dims[ 1 ] ), dsecoyz( dims[ 1 ] );
				for ( Size layer = 0; layer < dims[ 1 ]; ++layer ) {
					values(  layer ) = n_derivs_[ 1 ]( layer, row, col );
					dsecoy(  layer ) = n_derivs_[ 3 ]( layer, row, col );
					dsecoz(  layer ) = n_derivs_[ 2 ]( layer, row, col );
					dsecoyz( layer ) = n_derivs_[ 4 ]( layer, row, col );
				}
				CubicSpline cs, csz, csy, csyz;
				cs.train(   BORDER[ 1 ], START[ 1 ], DELTA[ 1 ], values , FIRSTBE[ 1 ] );
				csy.train(  BORDER[ 1 ], START[ 1 ], DELTA[ 1 ], dsecoy , FIRSTBE[ 1 ] );
				csz.train(  BORDER[ 1 ], START[ 1 ], DELTA[ 1 ], dsecoz , FIRSTBE[ 1 ] );
				csyz.train( BORDER[ 1 ], START[ 1 ], DELTA[ 1 ], dsecoyz, FIRSTBE[ 1 ] );
				for ( Size layer = 0; layer < dims[ 1 ]; ++layer ) {
					// We can't use a fixedsizearray1 here because the template
					// engine can't prove that N isn't 5 and nonetheless in this
					// branch.
					// Oh wait - we can just use N. I'm dumb.
					utility::fixedsizearray1< Size, N > position;
					position[ 1 ] = layer;// + 1;
					position[ 2 ] = row;
					position[ 3 ] = col;
					n_derivs_[ 5 ]( position ) = cs.get_dsecox()(   layer );
					n_derivs_[ 7 ]( position ) = csy.get_dsecox()(  layer );
					n_derivs_[ 6 ]( position ) = csz.get_dsecox()(  layer );
					n_derivs_[ 8 ]( position ) = csyz.get_dsecox()( layer );
				}
			}
		}
	} else { // reduction of dimension, but just by one.

		constexpr Size half_derivs = 1 << (N-1);

		for ( Size dim1 = 0; dim1 < dims[ 1 ]; ++dim1 ) {
			utility::fixedsizearray1< Size, N-1 > rest_of_dimensions;
			for ( Size ii = 1; ii <= N-1; ++ii ) { rest_of_dimensions[ ii ] = dims[ ii + 1 ]; }
			MathNTensor< Real, N-1 > values( rest_of_dimensions );

			// Loop over dim2 to dimN (rod1 to rodN-1)
			// assume "indices" <N> will hold bin indices in 1 to N-1
			utility::fixedsizearray1< Size, N > indices;
			utility::fixedsizearray1< Size, N > maxes;
			// initialize maxes to dims 2 to N in 1 to N-1
			for ( Size ii = 1; ii <= N-1; ++ii ) maxes[ ii ] = dims[ ii+1 ];

			Size p = 1;
			while ( indices[ N ] == 0 ) {

				utility::fixedsizearray1< Size, N > position;
				position[ 1 ] = dim1;// + 1;

				// position[1] is dim1, so position 2 maps to indices 1, indices N maps to nothing
				for ( Size ii = 1; ii <= N-1; ++ii ) { position[ ii + 1 ] = indices[ ii ]; }

				// All but first index goes to this MathMatrix
				utility::fixedsizearray1< Size, N-1 > rest_of_indices;
				for ( Size ii = 1; ii <= N-1; ++ii ) { rest_of_indices[ ii ] = indices[ ii ]; }
				values( rest_of_indices ) = n_derivs_[ 1 ]( position );


				// Incrementation
				++indices[ 1 ];
				while ( indices[ p ] == maxes[ p ] ) {
					indices[ p ] = 0;
					indices[ ++p ]++;
					if ( indices[ p ] != maxes[ p ] )  p = 1;
				}
			}

			// Create fixedsizearrays with all but the first dim for each trait
			PolycubicSpline< N-1 > ps;
			utility::fixedsizearray1< BorderFlag, N-1 > border;
			utility::fixedsizearray1< Real, N-1 > start;
			utility::fixedsizearray1< Real, N-1 > delta;
			utility::fixedsizearray1< bool, N-1 > lin_cont;
			utility::fixedsizearray1< std::pair< Real, Real >, N-1 > firstbe( std::pair< Real, Real >( 0, 0 ) );

			for ( Size ii = 2; ii <= N; ++ii ) {
				border[ ii - 1 ] = BORDER[ ii ];
				start[ ii - 1 ] = START[ ii ];
				delta[ ii - 1 ] = DELTA[ ii ];
				lin_cont[ ii - 1 ] = LINCONT[ ii ];
				firstbe[ ii - 1 ] = FIRSTBE[ ii ];
			}

			ps.train( border, start, delta, values, lin_cont, firstbe);

			// Replace layers in n_derivs_ from this spline, ceux qui correspondent
			for ( Size di = 2; di <= half_derivs; ++di ) {
				n_derivs_[ di ].replace_layer( dim1, ps.get_deriv( di ) );
			}
			/*
			n_derivs_[3].replace_layer( dim1, bs.get_dsecox());
			n_derivs_[2].replace_layer( dim1, bs.get_dsecoy());
			n_derivs_[4].replace_layer( dim1, bs.get_dsecoxy());
			*/
		}

		// all but first dim in indices2< N >, last index means nothing
		utility::fixedsizearray1< Size, N > indices2;
		utility::fixedsizearray1< Size, N > maxes2;
		// initialize maxes to dims 2 to N in 1 to N-1
		for ( Size ii = 1; ii <= N-1; ++ii ) maxes2[ ii ] = dims[ ii+1 ];

		Size p2 = 1;
		while ( indices2[ N ] == 0 ) {

			// I THINK that what this means is basically, once you've trained
			// half the dimension N-1 layers, you train N/2 layers of dimension
			// 1 to bridge from the "no x deriv" layers to the "x deriv" layers.
			// HELP HELP I AM HERE HERE HERE

			utility::fixedsizearray1< MathVector< Real >, half_derivs > first_half_derivs;
			for ( Size ii = 1; ii <= half_derivs; ++ii ) {
				first_half_derivs[ ii ] = MathVector< Real >( dims[ 1 ] );
			}

			// Initialize each layer for each one
			for ( Size layer = 0; layer < dims[ 1 ]; ++layer ) {
				for ( Size ii = 1; ii <= half_derivs; ++ii ) {
					utility::fixedsizearray1< Size, N > position;
					position[ 1 ] = layer;
					for ( Size jj = 1; jj <= N-1; ++jj ) {
						position[ jj+1 ] = indices2[ jj ];
					}
					first_half_derivs[ ii ]( layer ) = n_derivs_[ ii ]( position );
				}
			}

			// 1 << N-1 splines
			utility::vector1< CubicSpline > first_half_splines; //, half_derivs > first_half_splines;
			for ( Size ii = 1; ii <= half_derivs; ++ii ) {
				CubicSpline cs;
				cs.train( BORDER[ 1 ], START[ 1 ], DELTA[ 1 ], first_half_derivs[ ii ], FIRSTBE[ 1 ] );
				first_half_splines.push_back( cs );
				//first_half_splines[ ii ].train( BORDER[ 1 ], START[ 1 ], DELTA[ 1 ], first_half_derivs[ ii ], FIRSTBE[ 1 ] );
			}

			for ( Size layer = 0; layer < dims[ 1 ]; ++layer ) {
				utility::fixedsizearray1< Size, N > position;
				position[ 1 ] = layer;
				for ( Size jj = 1; jj <= N-1; ++jj ) {
					position[ jj+1 ] = indices2[ jj ];
				}

				for ( Size ii = 1; ii <= half_derivs; ++ii ) {
					n_derivs_[ ii + half_derivs ]( position ) = first_half_splines[ ii ].get_dsecox()( layer );
				}
			}


			// Incrementation
			++indices2[ 1 ];
			while ( indices2[ p2 ] == maxes2[ p2 ] ) {
				indices2[ p2 ] = 0;
				indices2[ ++p2 ]++;
				if ( indices2[ p2 ] != maxes2[ p2 ] )  p2 = 1;
			}
		}
	}
	return;
}


//! return value at certain (x, y, z)
template< Size N >
double PolycubicSpline< N >::F( utility::fixedsizearray1< double, N > const & xs ) const
{
	utility::fixedsizearray1< Size, N > dims;
	for ( Size ii = 1; ii <= N; ++ii ) {
		dims[ ii ] = n_derivs_[ 1 ].n_bins( ii );
	}

	//check if argument is in range for non-periodic splines
	// Note: we should unit test this out-of-range behavior and confirm it's reasonable;
	// if you're out of range for two axes you don't interpolate both?
	for ( Size ii = 1; ii <= N; ++ii ) {
		if ( ( border_[ ii ] != e_Periodic)
				&& ( ( xs[ ii ] < start_[ ii ]) || ( start_[ ii ] + ( dims[ ii ] - 1) * delta_[ ii ] < xs[ ii ])) ) {

			// Before start
			if ( xs[ ii ] < start_[ ii ] ) {
				utility::fixedsizearray1< Real, N > pos;
				for ( Size si = 1; si <= N; ++si ) {
					pos[ si ] = ( si == ii ) ? start_[ si ] : xs[ si ];
				}
				return F( pos ) + ( xs[ ii ] - start_[ ii ] ) * dFdxi( ii, pos );
			}

			// After end
			if ( xs[ ii ] > start_[ ii ] + ( dims[ ii ] - 1 ) * delta_[ ii ] ) {
				utility::fixedsizearray1< Real, N > pos;
				for ( Size si = 1; si <= N; ++si ) {
					pos[ si ] = ( si == ii ) ? start_[ ii ] + ( dims[ si ] - 1 ) * delta_[ ii ] : xs[ si ];
				}
				return F( pos ) + ( xs[ ii ] - start_[ ii ] - ( dims[ ii ] - 1 ) * delta_[ ii ] ) * dFdxi( ii, pos );
			}
		}
	}

	//determine i with start_[ 0]+(i-1)*delta_[ 0] < x < start_[ 0]+i*delta_[ 0] for the correct supporting points
	// pos: the one-indexed left-lower-corner coordinates (bins, not next_bins)
	utility::fixedsizearray1< int, N > bins;
	for ( Size ii = 1; ii <= N; ++ii ) {
		bins[ ii ] = int( floor( ( xs[ ii ]-start_[ ii ] ) / delta_[ ii ] ) );
		while ( start_[ ii ] + bins[ ii ] * delta_[ ii ] < xs[ ii ] ) bins[ ii ]++;
	}

	// this formula was derived from the Spline2D formula
	// the idea is to 'combine every part of the 2D formula
	// with dzm, dzp, dz3m, dz3p'

	utility::fixedsizearray1< double, N > delta_akts;
	utility::fixedsizearray1< double, N > dps;
	utility::fixedsizearray1< double, N > dms;
	utility::fixedsizearray1< double, N > d3ps;
	utility::fixedsizearray1< double, N > d3ms;

	for ( Size ii = 1; ii <= N; ++ii ) {
		delta_akts[ ii ] = xs[ ii ] - start_[ ii ] - ( bins[ ii ] - 1 ) * delta_[ ii ];
		dps[ ii ] = delta_akts[ ii ] / delta_[ ii ];
		dms[ ii ] = 1.0 - dps[ ii ];
		d3ps[ ii ] = ( dps[ ii ] * dps[ ii ] * dps[ ii ] - dps[ ii ] ) * sqr( delta_[ ii ] ) / 6;
		d3ms[ ii ] = ( dms[ ii ] * dms[ ii ] * dms[ ii ] - dms[ ii ] ) * sqr( delta_[ ii ] ) / 6;

		//std::cout << ii << ": delta_akt " << delta_akts[ ii ] << " dps " << dps[ ii ] << " dms " << dms[ ii ] << " d3ps " << d3ps[ ii ] << " d3ms " << d3ms[ ii ] << "\n";
		while ( bins[ ii ] < 1 ) bins[ ii ] += dims[ ii ];
	}

	Real value = 0;
	for ( Size posi = 1; posi <= ( 1 << N ); ++posi ) {
		// if posi - 1 &
		utility::fixedsizearray1< Size, N > position;
		for ( Size ii = 1; ii <= N; ++ii ) {
			position[ ii ] = ( ( posi - 1 ) & static_cast< int >( ( 1 << (N - ii) ) ) ) ? bins[ ii ] % dims[ ii ] : ( bins[ ii ] - 1 ) % dims[ ii ];
		}
		//std::cout << "One position relevant to ( ";
		//for ( Size ii = 1; ii <= N; ++ii ) std::cout << xs[ ii ] << " ";
		//std::cout << " ) is ( ";
		//for ( Size ii = 1; ii <= N; ++ii ) std::cout << position[ ii ] << " ";
		//std::cout << " )\n";

		for ( Size iid = 1; iid <= ( 1 << N ); ++iid ) {
			double coeff = n_derivs_[ iid ]( position );
			for ( Size ii = 1; ii <= N; ++ii ) {
				bool const deriv_being_taken_ii = ( ( iid - 1 ) & static_cast< int >( 1 << (N - ii) ) );
				bool const m_vs_p_ii = ( ( posi - 1 ) & static_cast< int >( 1 << (N - ii) ) );

				if ( deriv_being_taken_ii ) {
					coeff *= m_vs_p_ii ? d3ps[ ii ] : d3ms[ ii ];
				} else {
					coeff *= m_vs_p_ii ? dps[ ii ] : dms[ ii ];
				}
			}

			//if ( coeff != 0 ) {
			// std::cout << "iid " << iid << " coeff " << coeff << std::endl;
			//}
			value += coeff;
		}
	}

	return value;
}

//! return partial derivative at certain (x, y, z) for x
template< Size N >
double PolycubicSpline< N >::dFdxi( Size n, const utility::fixedsizearray1< double, N > & xs ) const
{
	utility::fixedsizearray1< Size, N > dims;
	for ( Size i = 1; i <= N; ++i ) {
		dims[ i ] = n_derivs_[ 1 ].n_bins( i );
	}

	//check if argument is in range for non-periodic splines
	for ( Size ii = 1; ii <= N; ++ii ) {
		if ( ( border_[ ii ] != e_Periodic)
				&& ( ( xs[ ii ] < start_[ ii ]) || ( start_[ ii ] + ( dims[ ii ]-1) * delta_[ ii ] < xs[ ii ])) ) {
			if ( xs[ ii ] < start_[ ii ] ) {
				utility::fixedsizearray1< Real, N > pos;
				for ( Size si = 1; si <= N; ++si ) {
					pos[ si ] = ( si == ii ) ? start_[ si ] : xs[ si ];
				}
				return dFdxi( n, pos );
			}
			if ( xs[ ii ] > start_[ ii ] + ( dims[ ii ] - 1 ) * delta_[ ii ] ) {
				utility::fixedsizearray1< Real, N > pos;
				for ( Size si = 1; si <= N; ++si ) {
					pos[ si ] = ( si == ii ) ? start_[ ii ] + ( dims[ si ] - 1 ) * delta_[ ii ] : xs[ si ];
				}
				return dFdxi( n, pos );
			}
		}
	}

	//determine i with start_+(i-1)*delta_ < x < start_+i*delta_ for the correct supporting points
	//determine i with start_[ 0]+(i-1)*delta_[ 0] < x < start_[ 0]+i*delta_[ 0] for the correct supporting points
	utility::fixedsizearray1< int, N > bins;
	for ( Size ii = 1; ii <= N; ++ii ) {
		bins[ ii ] = int( floor( ( xs[ ii ]-start_[ ii ]) / delta_[ ii ] ) );
		while ( start_[ ii ] + bins[ ii ] * delta_[ ii ] < xs[ ii ] ) bins[ ii ]++;
	}

	utility::fixedsizearray1< double, N > delta_akts;
	utility::fixedsizearray1< double, N > dps;
	utility::fixedsizearray1< double, N > dms;
	utility::fixedsizearray1< double, N > d3ps;
	utility::fixedsizearray1< double, N > d3ms;

	for ( Size ii = 1; ii <= N; ++ii ) {
		delta_akts[ ii ] = xs[ ii ] - start_[ ii ] - ( bins[ ii ] - 1) * delta_[ ii ];
		dps[ ii ] = delta_akts[ ii ] / delta_[ ii ];
		dms[ ii ] = 1 - dps[ ii ];
		d3ps[ ii ] = ( dps[ ii ] * dps[ ii ] * dps[ ii ] - dps[ ii ] ) * sqr( delta_[ ii ] ) / 6;
		d3ms[ ii ] = ( dms[ ii ] * dms[ ii ] * dms[ ii ] - dms[ ii ] ) * sqr( delta_[ ii ] ) / 6;
		//std::cout << ii << ": delta_akt " << delta_akts[ ii ] << " dps " << dps[ ii ] << " dms " << dms[ ii ] << " d3ps " << d3ps[ ii ] << " d3ms " << d3ms[ ii ] << "\n";

		while ( bins[ ii ] < 1 ) bins[ ii ] += dims[ ii ];
	}

	// We are finding the derivative d xn. We never should find ourselves using d3ps[ n ] or d3ms[ n ]

	double value = 0;

	for ( Size posi = 1; posi <= ( 1 << N ); ++posi ) {

		// Set up position from point index
		utility::fixedsizearray1< Size, N > position;
		for ( Size ii = 1; ii <= N; ++ii ) {
			position[ ii ] = ( ( posi - 1 ) & static_cast< int >( ( 1 << (N - ii) ) ) ) ? bins[ ii ] % dims[ ii ] : ( bins[ ii ] - 1 ) % dims[ ii ];
		}

		// Each derivative's term obeys particular rules for each position.
		// 1. Is the derivative taken wrt that variable?
		//    a. If so... Is the position bin - 1?
		//       i.  If so, multiply coefficient by d3ms[ ii ]
		//       ii. If not, multiply by d3ps[ ii ]
		//    b. If not.. Is the position bin - 1?
		//       i.  If so, multiply coefficient by dms[ ii ]
		//       ii. If not, multiply by dps[ ii ]
		// 2. Oh, but is this also n? n works differently.
		//    a. if not deriv: 1/delta[ii]
		//    b. if deriv, (3 * dxm*dxm - 1) * delta_[ 0] / 6, m or p based on above

		for ( Size iid = 1; iid <= ( 1 << N ); ++iid ) {

			// This is set to match TricubicSpline. Initially the opposite sign
			// convention made sense to me.
			double coeff = -1;
			if ( ( posi - 1 ) & static_cast< int >( 1 << (N - n) ) ) {
				coeff *= -1;
			}

			for ( Size ii = 1; ii <= N; ++ii ) {
				bool const deriv_being_taken_ii = ( ( iid - 1 ) & static_cast< int >( 1 << (N - ii) ) );
				bool const m_vs_p_ii = ( ( posi - 1 ) & static_cast< int >( 1 << (N - ii) ) );

				if ( ii == n ) {
					if ( deriv_being_taken_ii ) {
						Real const multiplicand = ( m_vs_p_ii ) ? dps[ii] : dms[ii];
						coeff *= n_derivs_[ iid ]( position ) * (3 * multiplicand * multiplicand - 1) * delta_[ n ] / 6;
					} else {
						coeff *= n_derivs_[ iid ]( position ) / delta_[ n ];
					}
				} else {

					if ( deriv_being_taken_ii ) {
						coeff *= m_vs_p_ii ? d3ps[ ii ] : d3ms[ ii ];
					} else {
						coeff *= m_vs_p_ii ? dps[ ii ] : dms[ ii ];
					}
				}
			}
			value += coeff;
		}
	}

	return value;
}

template< Size N >
utility::fixedsizearray1< Real, N >
PolycubicSpline< N >::dFdall( utility::fixedsizearray1< Real, N > const & xs ) const
{
	utility::fixedsizearray1< Real, N > values;
	for ( Size n = 1; n <= N; ++n ) {
		values[ n ] = dFdxi( n, xs );
	}
	return values;
}


}//end namespace spline
}//end namespace interpolation
}//end namespace numeric

#endif

