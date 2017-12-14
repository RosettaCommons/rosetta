// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/toolbox/match_enzdes_util/ExternalGeomSampler.cc
/// @brief
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

#include <protocols/toolbox/match_enzdes_util/ExternalGeomSampler.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace toolbox {
namespace match_enzdes_util {

ExternalGeomSampler::~ExternalGeomSampler() = default;

ExternalGeomSampler::ExternalGeomSampler() :
	parent(),
	dis_D1D2_( 0.0 ),
	dis_D2D3_( 0.0 ),
	ang_D1D2D3_( 0.0 ),
	transforms_uptodate_( false )
{}


ExternalGeomSampler::ExternalGeomSampler( ExternalGeomSampler const & other )
:
	parent(),
	dis_U1D1_samples_( other.dis_U1D1_samples_ ),
	ang_U2D1_samples_( other.ang_U2D1_samples_ ),
	tor_U3D1_samples_( other.tor_U3D1_samples_ ),
	ang_U1D2_samples_( other.ang_U1D2_samples_ ),
	tor_U2D2_samples_( other.tor_U2D2_samples_ ),
	tor_U1D3_samples_( other.tor_U1D3_samples_ ),
	dis_D1D2_( other.dis_D1D2_ ),
	dis_D2D3_( other.dis_D2D3_ ),
	ang_D1D2D3_( other.ang_D1D2D3_ ),
	transforms_uptodate_( other.transforms_uptodate_ ),
	transforms_( other.transforms_ )
{}

ExternalGeomSampler &
ExternalGeomSampler::operator = ( ExternalGeomSampler const & rhs )
{
	if ( & rhs != this ) {
		dis_U1D1_samples_ = rhs.dis_U1D1_samples_;
		ang_U2D1_samples_ = rhs.ang_U2D1_samples_;
		tor_U3D1_samples_ = rhs.tor_U3D1_samples_;
		ang_U1D2_samples_ = rhs.ang_U1D2_samples_;
		tor_U2D2_samples_ = rhs.tor_U2D2_samples_;
		tor_U1D3_samples_ = rhs.tor_U1D3_samples_;
		dis_D1D2_ = rhs.dis_D1D2_;
		dis_D2D3_ = rhs.dis_D2D3_;
		ang_D1D2D3_ = rhs.ang_D1D2D3_;
		transforms_uptodate_ = rhs.transforms_uptodate_;
		transforms_ = rhs.transforms_;
	}

	return *this;
}


void ExternalGeomSampler::set_tor_U3D1_samples( utility::vector1< Real > const & tor_U3D1_samples )
{
	transforms_uptodate_ = false;
	tor_U3D1_samples_ = tor_U3D1_samples;
}

void ExternalGeomSampler::set_dis_U1D1_samples( utility::vector1< Real > const & dis_U1D1_samples )
{
	transforms_uptodate_ = false;
	dis_U1D1_samples_ = dis_U1D1_samples;
}

void ExternalGeomSampler::set_ang_U2D1_samples( utility::vector1< Real > const & ang_U2D1_samples )
{
	transforms_uptodate_ = false;
	ang_U2D1_samples_ = ang_U2D1_samples;
}

void ExternalGeomSampler::set_ang_U1D2_samples( utility::vector1< Real > const & ang_U1D2_samples )
{
	transforms_uptodate_ = false;
	ang_U1D2_samples_ = ang_U1D2_samples;
}

void ExternalGeomSampler::set_tor_U2D2_samples( utility::vector1< Real > const & tor_U2D2_samples )
{
	transforms_uptodate_ = false;
	tor_U2D2_samples_ = tor_U2D2_samples;
}

void ExternalGeomSampler::set_tor_U1D3_samples( utility::vector1< Real > const & tor_U1D3_samples )
{
	transforms_uptodate_ = false;
	tor_U1D3_samples_ = tor_U1D3_samples;
}

void ExternalGeomSampler::set_tor_U3D1_samples( std::list< Real > const & tor_U3D1_samples )
{
	transforms_uptodate_ = false;
	tor_U3D1_samples_.resize( tor_U3D1_samples.size() );
	std::copy( tor_U3D1_samples.begin(), tor_U3D1_samples.end(), tor_U3D1_samples_.begin() );
}

void ExternalGeomSampler::set_dis_U1D1_samples( std::list< Real > const & dis_U1D1_samples )
{
	transforms_uptodate_ = false;
	dis_U1D1_samples_.resize( dis_U1D1_samples.size() );
	std::copy( dis_U1D1_samples.begin(), dis_U1D1_samples.end(), dis_U1D1_samples_.begin() );
}

void ExternalGeomSampler::set_ang_U2D1_samples( std::list< Real > const & ang_U2D1_samples )
{
	transforms_uptodate_ = false;
	ang_U2D1_samples_.resize( ang_U2D1_samples.size() );
	std::copy( ang_U2D1_samples.begin(), ang_U2D1_samples.end(), ang_U2D1_samples_.begin() );
}

void ExternalGeomSampler::set_ang_U1D2_samples( std::list< Real > const & ang_U1D2_samples )
{
	transforms_uptodate_ = false;
	ang_U1D2_samples_.resize( ang_U1D2_samples.size() );
	std::copy( ang_U1D2_samples.begin(), ang_U1D2_samples.end(), ang_U1D2_samples_.begin() );
}

void ExternalGeomSampler::set_tor_U2D2_samples( std::list< Real > const & tor_U2D2_samples )
{
	transforms_uptodate_ = false;
	tor_U2D2_samples_.resize( tor_U2D2_samples.size() );
	std::copy( tor_U2D2_samples.begin(), tor_U2D2_samples.end(), tor_U2D2_samples_.begin() );
}

void ExternalGeomSampler::set_tor_U1D3_samples( std::list< Real > const & tor_U1D3_samples )
{
	transforms_uptodate_ = false;
	tor_U1D3_samples_.resize( tor_U1D3_samples.size() );
	std::copy( tor_U1D3_samples.begin(), tor_U1D3_samples.end(), tor_U1D3_samples_.begin() );
}

void ExternalGeomSampler::set_dis_D1D2( Real distance ) {
	dis_D1D2_ = distance;
	transforms_uptodate_ = false;
}

void ExternalGeomSampler::set_dis_D2D3( Real distance ) {
	dis_D2D3_ = distance;
	transforms_uptodate_ = false;
}

void ExternalGeomSampler::set_ang_D1D2D3( Real ang_in_degrees ) {
	ang_D1D2D3_ = ang_in_degrees;
	transforms_uptodate_ = false;
}

/// @brief Must be called after the samples are set, and the
/// internal geometry of the three downstream coordinates
/// (point 4, 5, and 6) are described.  Does nothing if the transforms
/// are up to date.
void
ExternalGeomSampler::precompute_transforms()
{
	transforms_.resize( n_external_transforms );
	transforms_[ HT_tor_U3D1 ].resize( tor_U3D1_samples_.size() );
	transforms_[ HT_ang_U2D1 ].resize( ang_U2D1_samples_.size() );
	transforms_[ HT_tor_U2D2 ].resize( tor_U2D2_samples_.size() );
	transforms_[ HT_ang_U1D2 ].resize( ang_U1D2_samples_.size() );
	transforms_[ HT_tor_U1D3 ].resize( tor_U1D3_samples_.size() );

	for ( Size ii = 1; ii <= tor_U3D1_samples_.size(); ++ii ) {
		transforms_[ HT_tor_U3D1 ][ ii ].set_zaxis_rotation_deg( tor_U3D1_samples_[ ii ] );
	}

	for ( Size ii = 1; ii <= ang_U2D1_samples_.size(); ++ii ) {
		// negative rotation about the x axis of a magnitude 180 - angle 2-3-4
		transforms_[ HT_ang_U2D1 ][ ii ].set_xaxis_rotation_deg( -1 * ( 180 - ang_U2D1_samples_[ ii ] ));
	}

	for ( Size ii = 1; ii <= tor_U2D2_samples_.size(); ++ii ) {
		transforms_[ HT_tor_U2D2 ][ ii ].set_zaxis_rotation_deg( tor_U2D2_samples_[ ii ] );
	}


	for ( Size ii = 1; ii <= ang_U1D2_samples_.size(); ++ii ) {
		// negative rotation about the x axis of a magnitude 180 - angle B
		transforms_[ HT_ang_U1D2 ][ ii ].set_xaxis_rotation_deg( -1 * ( 180 - ang_U1D2_samples_[ ii ] ) );
		// pre-multiply by the stride from atom 4 to atom 5.
		transforms_[ HT_ang_U1D2 ][ ii ].walk_along_z( dis_D1D2_ );
	}

	HTReal ht_ang_d1d2d3;
	ht_ang_d1d2d3.set_xaxis_rotation_deg( -1 * ( 180 - ang_D1D2D3_ ) );
	for ( Size ii = 1; ii <= tor_U1D3_samples_.size(); ++ii ) {
		/// pre-multiply by the bond angle and the step along z.
		HTReal ht_tor_U1D3_zrot;
		ht_tor_U1D3_zrot.set_zaxis_rotation_deg( tor_U1D3_samples_[ ii ] );
		transforms_[ HT_tor_U1D3 ][ ii ] = ht_tor_U1D3_zrot * ht_ang_d1d2d3;
		transforms_[ HT_tor_U1D3 ][ ii ].walk_along_z( dis_D2D3_ );
	}

	transforms_uptodate_ = true;
}


void
ExternalGeomSampler::flip_upstream_downstream_samples()
{
	transforms_uptodate_ = false;

	utility::vector1< Real > new_angU2D1_samples = ang_U1D2_samples_;
	ang_U1D2_samples_ = ang_U2D1_samples_;
	ang_U2D1_samples_ = new_angU2D1_samples;

	utility::vector1< Real > new_torU3D1_samples = tor_U1D3_samples_;
	tor_U1D3_samples_ = tor_U3D1_samples_;
	tor_U3D1_samples_ = new_torU3D1_samples;
}

}
}
}

