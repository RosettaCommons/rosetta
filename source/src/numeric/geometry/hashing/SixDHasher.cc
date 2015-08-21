// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   /numeric/geometry/hashing/SixDHasher.cc
/// @brief  Implementation of 6D hasher classes
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini


// Unit headers
#include <numeric/geometry/hashing/SixDHasher.hh>
//#include <basic/Tracer.hh>

#include <boost/foreach.hpp>


namespace numeric {
namespace geometry {
namespace hashing {

/// @details Auto-generated virtual destructor
SixDCoordinateBinner::~SixDCoordinateBinner() {}

//static basic::Tracer TR("numeric.geometry.hashing.SixDHasher");


SixDCoordinateBinner::SixDCoordinateBinner()
{}

SixDCoordinateBinner::SixDCoordinateBinner(
	BoundingBox const & bounding_box,
	Size3 const & euler_offsets,
	utility::fixedsizearray1< Real, 6 > bin_widths
) :
	bounding_box_( bounding_box ),
	bin_widths_( bin_widths )
{
	Vector span = bounding_box_.upper() - bounding_box_.lower();
	Vector new_upper = bounding_box_.upper();

	for ( Size ii = 1; ii <= 3; ++ii ) {
		dimsizes_[ ii ] = static_cast< Size > ( span( ii ) / bin_widths_[ ii ] );
		if ( dimsizes_[ ii ] * bin_widths_[ ii ] < span( ii ) ) {
			dimsizes_[ ii ] += 1;
			new_upper( ii ) = bounding_box_.lower()( ii ) + dimsizes_[ ii ] * bin_widths_[ ii ];
		}
	}
	bounding_box_.set_upper( new_upper );

	for ( Size ii = 4; ii <= 5; ++ii ) {
		dimsizes_[ ii ] = static_cast< Size > ( 360.0 / bin_widths_[ ii ] );
	}
	dimsizes_[ 6 ] = static_cast< Size > ( 180.0 / bin_widths_[ 6 ] );
	if ( dimsizes_[6] == 0 ) dimsizes_[6] = 1;

	/// Add an extra bin so that values near 180 ( wi binwidth/2) and near 0 ( wi binwidth/2 ) can be joined.
	if ( euler_offsets[ 3 ] != 0 ) {
		++dimsizes_[ 6 ];
	}

	dimprods_[ 6 ] = 1;
	for ( Size ii = 5; ii >= 1; --ii ) {
		dimprods_[ ii ] = dimprods_[ ii + 1 ] * dimsizes_[ ii + 1 ];
	}

	for ( Size ii = 1; ii <= 3; ++ii ) {
		if ( euler_offsets[ ii ] != 0 ) {
			euler_offsets_[ ii ] = bin_widths_[ ii + 3 ] / 2;
		} else {
			euler_offsets_[ ii ] = 0.0;
		}
	}
	for ( Size ii = 1; ii <= 6; ++ii ) { halfbin_widths_[ ii ] = bin_widths_[ ii ] / 2; }
}


/// @details When floating point comparison breaks down, it is possible to have a point
/// hash outside of the bounding volume:  359.999999 / 10.0 -> 36 instead of 35.
/// For this reason, it's important to mod the resulting bin index by the number of bins.
Bin6D
SixDCoordinateBinner::bin6( Real6 const & values ) const {

	numeric::xyzVector< numeric::Real > xyzcoord( values[ 1 ], values[ 2 ], values[ 3 ] );

	assert( bounding_box_.contains( xyzcoord ) );
	assert( values[ 4 ] >= 0.0 && values[ 4 ] < 360.0 );
	assert( values[ 5 ] >= 0.0 && values[ 5 ] < 360.0 );
	assert( values[ 6 ] >= 0 && values[ 6 ] <= 180.0 );

	numeric::xyzVector< numeric::Real > from_corner = xyzcoord - bounding_box_.lower();
	Bin6D bins;

	bins[ 1 ] = static_cast< Size > ( from_corner.x() / bin_widths_[ 1 ] );
	if ( bins[ 1 ] == dimsizes_[ 1 ] ) --bins[ 1 ];

	bins[ 2 ] = static_cast< Size > ( from_corner.y() / bin_widths_[ 2 ] );
	if ( bins[ 2 ] == dimsizes_[ 2 ] ) --bins[ 2 ];

	bins[ 3 ] = static_cast< Size > ( from_corner.z() / bin_widths_[ 3 ] );
	if ( bins[ 3 ] == dimsizes_[ 3 ] ) --bins[ 3 ];

	Real3 euler = wrap_euler_angles( values );

	bins[ 4 ] = static_cast< Size > ( euler[ 1 ] / bin_widths_[ 4 ] ) % dimsizes_[ 4 ];
	bins[ 5 ] = static_cast< Size > ( euler[ 2 ] / bin_widths_[ 5 ] ) % dimsizes_[ 5 ];
	bins[ 6 ] = static_cast< Size > ( euler[ 3 ] / bin_widths_[ 6 ] ) % dimsizes_[ 6 ];

	return bins;
}

std::vector < boost::uint64_t >
SixDCoordinateBinner::radial_bin_index( Size radius, Real6 const & center ) const {
	std::vector < boost::uint64_t > bin_indices;
	Bin6D center_bin = bin6( center );
	std::vector < SBin6D > offsets = offset_tree_.lookup(radius, center_bin, dimsizes_);
	for ( Size j = 0; j < offsets.size(); ++j ) {
		Bin6D candidate;
		// apparently fixedarray1 (bin6d) has no + operator to another fixedarray1
		for ( Size k = 1; k <= 6; ++k ) {
			//go through and add
			candidate[k] = center_bin[k] + offsets[j][k];
			//TR<< candidate[k] << "\t" << center_bin[k] << "\t" << offsets[j][k] << std::endl;
		}
		// TR << std::endl;
		bin_indices.push_back( bin_index( candidate) );
	}
	return bin_indices;
}


/// @details
Bin6D
SixDCoordinateBinner::halfbin6( Real6 const & values ) const
{
	Bin6D bin = bin6( values );

	numeric::xyzVector< numeric::Real > xyzcoord( values[ 1 ], values[ 2 ], values[ 3 ] );
	numeric::xyzVector< numeric::Real > from_corner = xyzcoord - bounding_box_.lower();

	numeric::xyzVector< numeric::Real > r3lower_corner = bounding_box_.lower() +
		Vector( bin[1]*bin_widths_[1], bin[2]*bin_widths_[2], bin[3]*bin_widths_[3] );
	numeric::xyzVector< numeric::Real > remnant = xyzcoord - r3lower_corner;

	Bin6D halfbin( 0 );
	for ( Size ii = 1; ii <= 3; ++ii ) {
		if ( remnant[ ii-1 ] >= halfbin_widths_[ ii ] ) {
			halfbin[ii] = 1;
		}
	}

	Real3 euler = wrap_euler_angles( values );
	for ( Size ii = 1; ii <= 3; ++ii ) {
		Real lower_boundary = bin[ii+3]*bin_widths_[ii+3];
		Real dist = euler[ii] - lower_boundary;
		//std::cout << "Euler " << ii << " " << values[ii+3] <<  " " << euler[ ii ] << " " << lower_boundary << " " << dist << " " << euler_offsets_[ii] << std::endl;
		if ( dist >= halfbin_widths_[ii+3] ) {
			halfbin[ii+3] = 1;
		}
	}
	return halfbin;
}

Real6
SixDCoordinateBinner::bin_center_point( Bin6D const & bin ) const
{
	Real6 center;
	for ( Size ii = 1; ii <= 3; ++ii ) center[ ii ] = bounding_box_.lower()( ii ) + bin[ ii ] * bin_widths_[ ii ] + halfbin_widths_[ ii ];
	for ( Size ii = 1; ii <= 3; ++ii ) center[ ii + 3 ] = euler_offsets_[ ii ] + bin[ ii + 3 ] * bin_widths_[ ii + 3 ] + halfbin_widths_[ ii + 3 ];
	return center;
}


Real3
SixDCoordinateBinner::wrap_euler_angles( Real6 const & values ) const
{
	Real3 euler;

	if ( values[ 6 ] > euler_offsets_[ 3 ] ) {
		euler[ 3 ] = values[ 6 ] + euler_offsets_[ 3 ];

	} else {
		euler[ 3 ] = values[ 6 ];
	}

	if ( ( values[ 6 ] < euler_offsets_[ 3 ] || values[ 6 ] >= 180.0 - euler_offsets_[ 3 ] ) &&
			( values[ 4 ] - euler_offsets_[ 1 ] > 180.0 ) ) {
		/// WRAP if phi > 180 to the region of negative theta rotation.
		/// The idea is that, when we're wrapping theta, half of the points have to stay fixed so they end up in the same
		/// bin: if all points wrapped, then none would land in the same bin.
		euler[ 1 ] = values[ 4 ] - euler_offsets_[ 1 ] - 180.0;
		euler[ 2 ] = values[ 5 ] - euler_offsets_[ 2 ] - 180.0;
		if ( euler[ 2 ] < 0 ) euler[ 2 ] += 360;
	} else {
		/// Leave phi/psi in their usual positions
		euler[ 1 ] = values[ 4 ] - euler_offsets_[ 1 ];
		if ( euler[ 1 ] < 0 ) euler[ 1 ] += 360;

		euler[ 2 ] = values[ 5 ] - euler_offsets_[ 2 ];
		if ( euler[ 2 ] < 0 ) euler[ 2 ] += 360;
	}
	return euler;
}


SixDOffsetTree::SixDOffsetTree(){
}

void SixDOffsetTree::init( Size max_radius ){
	SBin6D test ( 0 );
	// completely stupid way of doing
	// just testing on how fast it is
	SSize smax_radius = max_radius;
	for ( SSize i = -smax_radius; i <= smax_radius; i++ ) {
		test[1] = i;
		for ( SSize j = -smax_radius; j <= smax_radius; j++ ) {
			test[2] = j;
			if ( sum_radius(test, 2) > max_radius ) continue;
			for ( SSize k = -smax_radius; k <= smax_radius; k++ ) {
				test[3] = k;
				if ( sum_radius(test, 3) > max_radius ) continue;
				for ( SSize m = -smax_radius; m <= smax_radius; m++ ) {
					test[4] = m;
					if ( sum_radius(test, 4) > max_radius ) continue;
					for ( SSize n = -smax_radius; n <= smax_radius; n++ ) {
						test[5] = n;
						if ( sum_radius(test, 5) > max_radius ) continue;
						for ( SSize p = -smax_radius; p <= smax_radius; p++ ) {
							test[6] = p;
							if ( sum_radius(test) > max_radius ) continue;
							insert( test );
						}
					}
				}
			}
		}
	}
}

std::vector< SBin6D > SixDOffsetTree::lookup( Size radius, Bin6D const & center, Bin6D const & bounds ) const {
	typedef boost::unordered_map<SSize, Size>::value_type map_itr;
	std::vector < SBin6D > offset_list;
	if ( radius > data_.size() ) {
		std::cerr << "numeric.geometry.hashing.SixDHasher  Tree: Radius out of bounds" << std::endl;
	} else {
		Bin6D idx ( 0 );
		SBin6D depth ( 0 );
		BOOST_FOREACH ( map_itr i, data_[radius][0] ) {
			idx[1] = i.second;
			depth[1] = i.first;
			if ( depth[1] + center[1] >= bounds[1] ) continue;

			BOOST_FOREACH ( map_itr j, data_[radius][idx[1]] ) {
				idx[2] = j.second;
				depth[2] = j.first;
				if ( depth[2] + center[2] >= bounds[2] ) continue;

				BOOST_FOREACH ( map_itr k, data_[radius][idx[2]] ) {
					idx[3] = k.second;
					depth[3] = k.first;
					if ( depth[3] + center[3] >= bounds[3] ) continue;

					BOOST_FOREACH ( map_itr m, data_[radius][idx[3]] ) {
						idx[4] = m.second;
						depth[4] = m.first;
						if ( depth[4] + center[4] >= bounds[4] ) continue;

						BOOST_FOREACH ( map_itr n, data_[radius][idx[4]] ) {
							idx[5] = n.second;
							depth[5] = n.first;
							if ( depth[5] + center[5] >= bounds[5] ) continue;

							BOOST_FOREACH ( map_itr p, data_[radius][idx[5]] ) {
								idx[6] = p.second;
								depth[6] = p.first;
								if ( idx[6] == 1 && depth[6]+center[6] < bounds[6] ) offset_list.push_back(depth);
							}
						}
					}
				}
			}
		}
	}
	return offset_list;
}

numeric::Size SixDOffsetTree::sum_radius( SBin6D & input, numeric::Size range ){
	Size radius = 0;
	if ( range > 6 ) {
		std::cerr << "numeric.geometry.hashing.SixDHasher Bad range in sum_radius" << std::endl;
		return 1000;
	}
	for ( Size i = 1; i<= range; i++ ) {
		radius += std::abs(input[i]);
	}
	return radius;
}

bool SixDOffsetTree::insert( SBin6D & input, Size depth, Size caller ) {
	Size radius = sum_radius(input);
	while ( data_.size() <= radius ) {
		std::vector< boost::unordered_map < SSize, Size > > tmp;
		boost::unordered_map < SSize, Size > depth_1;
		tmp.push_back(depth_1);
		data_.push_back(tmp);
	}

	if ( depth == 1 ) {
		if ( data_[radius][0].find( input[depth] ) == data_[radius][0].end() ) {
			data_[radius][0][ input[depth] ] = 0;
		}
		insert( input, depth + 1, 0 );
		return true;
	}
	if ( data_[radius][caller][ input[depth-1] ] == 0 ) {
		// map doesn't exist
		boost::unordered_map < SSize,Size > tmp;
		data_[radius].push_back(tmp);
		data_[radius][caller][ input[depth-1] ] = data_[radius].size() - 1;
	}
	Size idx = data_[radius][caller][ input[depth-1]];
	if ( data_[radius][idx].find( input[depth] ) == data_[radius][idx].end() ) {
		data_[radius][idx][ input[depth] ] = 0;
	}
	if ( depth == 6 ) {
		data_[radius][idx][ input[depth] ] = 1;
		return true;
	} else {
		// and then go ONE LEVEL DEEPER! BRAWWWWMMMMMMM!
		return insert( input, depth + 1, idx );
	}
}

} //hashing
} //geometry
} //numeric

