// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author


// Unit headers
#include <core/pack/dunbrack/DunbrackRotamer.hh>

// Package headers
#include <core/pack/dunbrack/ChiSet.hh>

// Project headers
#include <basic/interpolate.hh>

// Numeric headers
#include <numeric/constants.hh>
#include <numeric/angle.functions.hh>
#include <numeric/random/random.hh>

#include <cmath>
#include <iostream>

#include <basic/Tracer.hh>
#include <basic/basic.hh>

#include <core/chemical/ResidueType.hh>
#include <utility/vector1.hh>

using basic::T;
using basic::Error;
using basic::Warning;

static THREAD_LOCAL basic::Tracer TR( "core.pack.dunbrack" );

namespace core {
namespace pack {
namespace dunbrack {

Size positive_pow( Size mantissa, Size exponent )
{
	if ( exponent == 0 ) return 1;
	if ( exponent == 1 ) return mantissa;
	if ( exponent == 2 ) return mantissa * mantissa;
	if ( exponent == 3 ) return mantissa * mantissa * mantissa;

	int tmp = positive_pow( mantissa, exponent/2 );
	if ( exponent % 2 == 0 ) return tmp * tmp;
	else return mantissa * tmp * tmp;
}

/// @details Fun Fact: virtual destructor must still be defined even if it's abstract
RotamerBuildingData::~RotamerBuildingData() {}

/// @details alternative interpolate_rotamers that uses polylinear interpolation
template < Size N >
void interpolate_rotamers(
	utility::fixedsizearray1< DunbrackRotamer< FOUR, N >, ( 1 << N ) > const & rot,
	utility::fixedsizearray1< Real, N > bb_err, Real binrange,
	Size nchi_aa,
	DunbrackRotamer< FOUR, N, Real > & interpolated_rotamer
)
{
	utility::vector1< Real > tmp;

	for ( Size i = 1; i <= nchi_aa; ++i ) {
		// get the dunbrack chi angle means
		Real interpolated_value;
		utility::fixedsizearray1< Real, ( 1 << N ) > chi_mean;
		utility::fixedsizearray1< Real, ( 1 << N ) > chi_sd;
		for ( Size roti = 1; roti <= rot.size(); ++roti ) {
			chi_mean.push_back( static_cast< Real > ( rot[ roti ].chi_mean( i ) ) );
			chi_sd.push_back( static_cast< Real > ( rot[ roti ].chi_sd( i ) ) );
		}
		interpolate_polylinear_by_value( chi_mean, bb_err, binrange, true /*treat_as_angles*/, interpolated_value, tmp );

		interpolated_rotamer.chi_mean( i, interpolated_value );

		// get the dunbrack chi angle sdevs
		interpolate_polylinear_by_value( chi_sd, bb_err, binrange, false /*don't treat_as_angles */, interpolated_value, tmp );
		interpolated_rotamer.chi_sd( i, interpolated_value );
		interpolated_rotamer.rotwell( i, rot[ 1 ].rotwell( i ) );

		// ctsa - check validity of result
		if ( interpolated_rotamer.chi_sd(i) < 0.0 ) {
			utility_exit_with_message( "interpolated_rotamer.chi_sd < 0 in fill_chi_set" );
		}

	} // i=1, i<= nchi_aa_

	utility::fixedsizearray1< Real, N > rot_prob;
	for ( Size roti = 1; roti <= rot.size(); ++roti ) {
		rot_prob[ roti ] = static_cast< Real > ( rot[ roti ].rotamer_probability() );
	}

	Real interpolated_prob;
	interpolate_polylinear_by_value( rot_prob, bb_err, binrange, false /*dont' treat_as_angles*/ , interpolated_prob, tmp );
	interpolated_rotamer.rotamer_probability( interpolated_prob );
}


/* OL: I thought copying the whole chi_set_vector is unnecessary and made a new version of this function */
void
expand_proton_chi_oldversion(
	pack::task::ExtraRotSample ex_samp_level,
	chemical::ResidueTypeCOP concrete_residue,
	Size proton_chi,
	utility::vector1< ChiSetOP > & chi_set_vector
)
{
	using namespace pack::task;
	// count the number of extra hydroxyl samples -- n
	// copy the chi_set_vector n times into a temporary
	// set the hydroxyl_chi value for these n samples in the temporary
	// assign the temporary to the input chi_set_vector

	utility::vector1< Real > const & samples = concrete_residue->proton_chi_samples( proton_chi );
	// i.e., -60, 60, 180

	// i.e., 10, 20  to model -40 -50 -60 -70 -80 etc.
	utility::vector1< Real > const & extra_samples = concrete_residue->proton_chi_extra_samples( proton_chi );

	Size chi_id = concrete_residue->proton_chi_2_chi( proton_chi );

	bool const include_extra( ex_samp_level != NO_EXTRA_CHI_SAMPLES );

	Size nsamples = samples.size() * ( 1 + ( include_extra ? extra_samples.size() * 2 : 0 ) );
	utility::vector1< ChiSetOP > newchi_vect( nsamples * chi_set_vector.size() );

	// copy old chi_set_vector nsample times
	for ( Size ii = 1; ii <= nsamples; ++ii ) {
		Size offset = (ii-1) * chi_set_vector.size();
		for ( Size jj = 1; jj <= chi_set_vector.size(); ++jj ) {
			newchi_vect[ jj + offset ] = ChiSetOP( new pack::dunbrack::ChiSet( *(chi_set_vector[ jj ]) ) );
		}
	}

	// add extra chi samples
	Size count( 1 );
	for ( Size ii = 1; ii <= samples.size(); ++ii ) {
		Real ii_sample = samples[ ii ]; //chi-angle
		for ( Size jj = 1; jj <= chi_set_vector.size(); ++jj ) {
			newchi_vect[ count ]->chi[ chi_id ] = ii_sample;
			++count;
			if ( include_extra ) {
				for ( Size kk = 1; kk <= extra_samples.size(); ++kk ) {
					newchi_vect[ count ]->chi[ chi_id ] = ii_sample + extra_samples[ kk ];
					++count;
					newchi_vect[ count ]->chi[ chi_id ] = ii_sample - extra_samples[ kk ];
					++count;
				}
			}
		}
	}

	debug_assert( count - 1 == nsamples * chi_set_vector.size() );
	chi_set_vector = newchi_vect;
}

/// olli -- I think this is slightly simpler to read than the original version
/// apl -- needs to find a new residence
void
expand_proton_chi(
	pack::task::ExtraRotSample ex_samp_level,
	chemical::ResidueTypeCOP concrete_residue,
	Size proton_chi,
	utility::vector1< ChiSetOP > & chi_set_vector
)
{
	using namespace pack::task;
	// count the number of extra hydroxyl samples -- n
	// copy the chi_set_vector n times into a temporary
	// set the hydroxyl_chi value for these n samples in the temporary
	// assign the temporary to the input chi_set_vector

	utility::vector1< Real > const & samples = concrete_residue->proton_chi_samples( proton_chi );

	// i.e., -60, 60, 180
	debug_assert( samples.size() > 0 ); // or less harsh and just a return ?

	// i.e., 10, 20  to model -40 -50 -60 -70 -80 etc.
	utility::vector1< Real > const & extra_samples = concrete_residue->proton_chi_extra_samples( proton_chi );

	Size chi_id = concrete_residue->proton_chi_2_chi( proton_chi );

	bool const include_extra( ex_samp_level != NO_EXTRA_CHI_SAMPLES );

	Size nsamples = samples.size() * ( 1 + ( include_extra ? extra_samples.size() * 2 : 0 ) );
	chi_set_vector.reserve( nsamples * chi_set_vector.size() ); //preallocate the necessary memory

	// add extra chi samples
	ChiSetOP new_chi_vec; ChiSetOP base_chi_vec;
	Size nr_of_old_elem = chi_set_vector.size();
	for ( Size jj = 1; jj <= nr_of_old_elem; ++jj ) {
		for ( Size ii = 1; ii <= samples.size(); ++ii ) {
			Real ii_sample = samples[ ii ]; //chi-angle
			if ( ii == 1 ) {  //change first chi in place:
				base_chi_vec = new_chi_vec = chi_set_vector[ jj ];
			} else { // make copies for all others
				new_chi_vec = ChiSetOP( new pack::dunbrack::ChiSet( *base_chi_vec ) );
				chi_set_vector.push_back( new_chi_vec );
			}
			new_chi_vec->chi[ chi_id ] = ii_sample;

			if ( include_extra ) {
				for ( Size kk = 1; kk <= extra_samples.size(); ++kk ) {
					chi_set_vector.push_back( new_chi_vec = ChiSetOP( new pack::dunbrack::ChiSet( *base_chi_vec  ) ) );
					new_chi_vec->chi[ chi_id ] = ii_sample  + extra_samples[ kk ];
					chi_set_vector.push_back( new_chi_vec = ChiSetOP( new pack::dunbrack::ChiSet( *base_chi_vec  ) ) );
					new_chi_vec->chi[ chi_id ] = ii_sample  - extra_samples[ kk ];
				} // for extra_samples
			} // include extra
		} // for sample.size()
	} // for jj (chi_set_vector)

	debug_assert( chi_set_vector.size()  == nsamples * nr_of_old_elem );
}


DunbrackRotamerSampleData::DunbrackRotamerSampleData() :
	nrchi_sample_( false ),
	nchi_( 0 ),
	probability_( 0.0 ),
	nrchi_lower_boundary_( 0.0 ),
	nrchi_upper_boundary_( 0.0 ),
	nrchi_probability_( 0.0 )
{}

DunbrackRotamerSampleData::DunbrackRotamerSampleData( bool is_nrchi_sample ) :
	nrchi_sample_( is_nrchi_sample ),
	nchi_( 0 ),
	probability_( 0.0 ),
	nrchi_lower_boundary_( 0.0 ),
	nrchi_upper_boundary_( 0.0 ),
	nrchi_probability_( 0.0 )
{}

DunbrackRotamerSampleData::~DunbrackRotamerSampleData() {}

void DunbrackRotamerSampleData::set_nrchi_sample( bool setting )
{
	nrchi_sample_ = setting;
}

void DunbrackRotamerSampleData::set_nchi( Size nchi ) {
	debug_assert( nchi_ <= DUNBRACK_MAX_SCTOR );
	nchi_ = nchi;
}

void DunbrackRotamerSampleData::set_rotwell(  Size chi_index, Size rotwell )
{
	debug_assert( chi_index > 0 && chi_index <= nchi_ );
	rot_well_[ chi_index ] = rotwell;
}

void DunbrackRotamerSampleData::set_rotwell( utility::vector1< Size > const & rotwell )
{
	debug_assert( ( ! nrchi_sample_ && rotwell.size() == nchi_ ) ||
		( nrchi_sample_ && rotwell.size() == nchi_ - 1 ) );
	std::copy( rotwell.begin(), rotwell.end(), rot_well_.begin() );
}

void DunbrackRotamerSampleData::set_chi_mean( Size chi_index, Real mean )
{
	debug_assert( chi_index > 0 && chi_index <= nchi_ );
	chi_mean_[ chi_index ] = mean;
}

void DunbrackRotamerSampleData::set_chi_sd( Size chi_index, Real sd )
{
	debug_assert( chi_index > 0 && chi_index <= nchi_ );
	chi_sd_[ chi_index ] = sd;
}

void DunbrackRotamerSampleData::set_prob( Real probability )
{
	probability_ = probability;
}

void DunbrackRotamerSampleData::set_nrchi_lower_boundary( Real low )
{
	debug_assert( nrchi_sample_ );
	nrchi_lower_boundary_ = low;
}

void DunbrackRotamerSampleData::set_nrchi_upper_boundary( Real high )
{
	debug_assert( nrchi_sample_ );
	nrchi_upper_boundary_ = high;
}

void DunbrackRotamerSampleData::set_nrchi_probability( Real nrchi_prob )
{
	debug_assert( nrchi_sample_ );
	nrchi_probability_ = nrchi_prob;
}

void
DunbrackRotamerSampleData::assign_random_chi(
	utility::vector1< Real > & chi_angles,
	numeric::random::RandomGenerator & RG,
	core::Real temperature /* scale distributions like a temperature,  default 1.0 = is xray temperature */
) const
{
	debug_assert( chi_angles.size() >= nchi() );

	for ( core::Size ii = 1; ii <= nchi(); ++ii ) {
		//if ( chi_is_nonrotameric( ii ) ) {
		// nonrotameric chi angles not currently supported
		//runtime_assert( false );
		//} else {
		chi_angles[ ii ] = basic::periodic_range( chi_mean_[ii] + RG.gaussian() * chi_sd_[ii] * temperature, 360.0 );
		//}
	}

	// set any remaining chi uniformly (proton chi)
	for ( core::Size ii = nchi()+1; ii <= chi_angles.size(); ++ii ) {
		chi_angles[ ii ] = basic::periodic_range( RG.uniform()*360.0 - 180.0, 360.0 );
	}
}

Real
DunbrackRotamerSampleData::chi_probability(
	utility::vector1< Real > const & chi_angles,
	core::Real temperature /* scale distributions like a temperature, default 1.0 = is xray temperature */
) const
{
	debug_assert( chi_angles.size() >= nchi() );
	Real const norm_gauss( std::sqrt( numeric::constants::r::pi_2 ) );
	Real prob(1);

	for ( Size ii = 1; ii <= nchi(); ++ii ) {
		//if ( chi_is_nonrotameric( ii ) ) {
		// nonrotameric chi angles not currently supported
		//runtime_assert( false );
		//} else {
		// Gaussian function with area 1 for rotameric angles
		Real const angle_diff( chi_mean_[ii] - numeric::nearest_angle_degrees( chi_angles[ii], chi_mean_[ii] ) );
		Real const sd ( chi_sd_[ii] * temperature );
		Real const variance( sd*sd );
		prob *= std::exp( -(angle_diff*angle_diff)/(2*variance) ) / sd / norm_gauss;
		//}
	}

	for ( Size ii = nchi()+1; ii <= chi_angles.size(); ++ii ) {
		// uniform function with area 1 for all other angles
		prob *= 1.0 / 360.0;
	}

	return prob;
}

} // namespace dunbrack
} // namespace scoring
} // namespace core
