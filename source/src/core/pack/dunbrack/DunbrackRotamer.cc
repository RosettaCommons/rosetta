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

static thread_local basic::Tracer TR( "core.pack.dunbrack" );

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
	//std::cout.precision(16);
	//std::cout << "  v00  " << v00 << "  d2dx200  " << d2dx200 << "  d2dy200  " << d2dy200 << "  d4dx2y200  " << d4dx2y200 << std::endl;
	//std::cout << "  v01  " << v01 << "  d2dx201  " << d2dx201 << "  d2dy201  " << d2dy201 << "  d4dx2y201  " << d4dx2y201 << std::endl;
	//std::cout << "  v10  " << v10 << "  d2dx210  " << d2dx210 << "  d2dy210  " << d2dy210 << "  d4dx2y210  " << d4dx2y210 << std::endl;
	//std::cout << "  v11  " << v11 << "  d2dx211  " << d2dx211 << "  d2dy211  " << d2dy211 << "  d4dx2y211  " << d4dx2y211 << std::endl;
	//std::cout << "  dxp " <<  dxp << std::endl;
	//std::cout << "  dyp " <<  dyp << std::endl;
	//std::cout << "  binwx " <<  binwx << std::endl;
	//std::cout << "  binwy " <<  binwy << std::endl;
	//std::cout << std::endl;

debug_assert( dxp >= 0 && dxp < 1.0 );
debug_assert( dyp >= 0 && dyp < 1.0 );
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

	val =
		dxm *   (  dym *     v00   +  dyp *     v01 )
		+ dxp * (  dym *     v10   +  dyp *     v11 )
		+dx3m * (  dym * d2dx200   +  dyp * d2dx201 )
		+dx3p * (  dym * d2dx210   +  dyp * d2dx211 )
		+ dxm * ( dy3m * d2dy200   + dy3p * d2dy201 )
		+ dxp * ( dy3m * d2dy210   + dy3p * d2dy211 )
		+dx3m * ( dy3m * d4dx2y200 + dy3p * d4dx2y201 )
		+dx3p * ( dy3m * d4dx2y210 + dy3p * d4dx2y211 );

	dvaldx =
		-( dym * v00 + dyp * v01 ) * invbinwx
		+( dym * v10 + dyp * v11 ) * invbinwx
		- ( 3 * dxm * dxm - 1) * binwx_over6 *( dym  * d2dx200   + dyp  * d2dx201 )
		+ ( 3 * dxp * dxp - 1) * binwx_over6 *( dym  * d2dx210   + dyp  * d2dx211 )
		-( dy3m * d2dy200 + dy3p * d2dy201 ) * invbinwx
		+( dy3m * d2dy210 + dy3p * d2dy211 ) * invbinwx
		- ( 3 * dxm * dxm - 1) * binwx_over6 *( dy3m * d4dx2y200 + dy3p * d4dx2y201 )
		+ ( 3 * dxp * dxp - 1) * binwx_over6 *( dy3m * d4dx2y210 + dy3p * d4dx2y211 );

	dvaldy =
		dxm   *( -v00 + v01 ) * invbinwy
		+ dxp *( -v10 + v11 ) * invbinwy
		+dx3m *( -d2dx200 + d2dx201) * invbinwy
		+dx3p *( -d2dx210 + d2dx211) * invbinwy
		+ dxm *( -( 3 * dym * dym - 1) * d2dy200   + ( 3 * dyp * dyp - 1) * d2dy201 ) * binwy_over6
		+ dxp *( -( 3 * dym * dym - 1) * d2dy210   + ( 3 * dyp * dyp - 1) * d2dy211 ) * binwy_over6
		+dx3m *( -( 3 * dym * dym - 1) * d4dx2y200 + ( 3 * dyp * dyp - 1) * d4dx2y201 ) * binwy_over6
		+dx3p *( -( 3 * dym * dym - 1) * d4dx2y210 + ( 3 * dyp * dyp - 1) * d4dx2y211 ) * binwy_over6;

}

void
tricubic_interpolation(
	Real v000, Real dvdx000, Real dvdy000, Real dvdz000, Real dvdxy000, Real dvdxz000, Real dvdyz000, Real dvdxyz000,
	Real v001, Real dvdx001, Real dvdy001, Real dvdz001, Real dvdxy001, Real dvdxz001, Real dvdyz001, Real dvdxyz001,
	Real v010, Real dvdx010, Real dvdy010, Real dvdz010, Real dvdxy010, Real dvdxz010, Real dvdyz010, Real dvdxyz010,
	Real v011, Real dvdx011, Real dvdy011, Real dvdz011, Real dvdxy011, Real dvdxz011, Real dvdyz011, Real dvdxyz011,
	Real v100, Real dvdx100, Real dvdy100, Real dvdz100, Real dvdxy100, Real dvdxz100, Real dvdyz100, Real dvdxyz100,
	Real v101, Real dvdx101, Real dvdy101, Real dvdz101, Real dvdxy101, Real dvdxz101, Real dvdyz101, Real dvdxyz101,
	Real v110, Real dvdx110, Real dvdy110, Real dvdz110, Real dvdxy110, Real dvdxz110, Real dvdyz110, Real dvdxyz110,
	Real v111, Real dvdx111, Real dvdy111, Real dvdz111, Real dvdxy111, Real dvdxz111, Real dvdyz111, Real dvdxyz111,
	Real dxp, Real dyp, Real dzp,
	Real binwx, Real binwy, Real binwz,
	Real & val,
	Real & dvaldx,
	Real & dvaldy,
	Real & dvaldz
)
{

	Real const invbinwx(1/binwx), invbinwy(1/binwy), invbinwz(1/binwz);
	Real const binwx_over_6(binwx/6), binwy_over_6(binwy/6), binwz_over_6(binwz/6);
	Real const dxm( 1-dxp), dym(1-dyp), dzm(1-dzp);

	Real const dx3p( ( dxp * dxp * dxp - dxp) * binwx * binwx_over_6 );
	Real const dx3m( ( dxm * dxm * dxm - dxm) * binwx * binwx_over_6 );
	Real const dy3p( ( dyp * dyp * dyp - dyp) * binwy * binwy_over_6 );
	Real const dy3m( ( dym * dym * dym - dym) * binwy * binwy_over_6 );
	Real const dz3p( ( dzp * dzp * dzp - dzp) * binwz * binwz_over_6 );
	Real const dz3m( ( dzm * dzm * dzm - dzm) * binwz * binwz_over_6 );

	val=
		dzm
		*(dxm*(dym*v000+dyp*v010)
		+dxp*(dym*v100+dyp*v110)
		+dx3m*(dym*dvdx000+dyp*dvdx010)
		+dx3p*(dym*dvdx100+dyp*dvdx110)
		+dxm*(dy3m*dvdy000+dy3p*dvdy010)
		+dxp*(dy3m*dvdy100+dy3p*dvdy110)
		+dx3m*(dy3m*dvdxy000+dy3p*dvdxy010)
		+dx3p*(dy3m*dvdxy100+dy3p*dvdxy110))

		+dzp
		*(dxm*(dym*v001+dyp*v011)
		+dxp*(dym*v101+dyp*v111)
		+dx3m*(dym*dvdx001+dyp*dvdx011)
		+dx3p*(dym*dvdx101+dyp*dvdx111)
		+dxm*(dy3m*dvdy001+dy3p*dvdy011)
		+dxp*(dy3m*dvdy101+dy3p*dvdy111)
		+dx3m*(dy3m*dvdxy001+dy3p*dvdxy011)
		+dx3p*(dy3m*dvdxy101+dy3p*dvdxy111))

		+dz3m
		*(dxm*(dym*dvdz000+dyp*dvdz010)
		+dxp*(dym*dvdz100+dyp*dvdz110)
		+dx3m*(dym*dvdxz000+dyp*dvdxz010)
		+dx3p*(dym*dvdxz100+dyp*dvdxz110)
		+dxm*(dy3m*dvdyz000+dy3p*dvdyz010)
		+dxp*(dy3m*dvdyz100+dy3p*dvdyz110)
		+dx3m*(dy3m*dvdxyz000+dy3p*dvdxyz010)
		+dx3p*(dy3m*dvdxyz100+dy3p*dvdxyz110))

		+dz3p
		*(dxm*(dym*dvdz001+dyp*dvdz011)
		+dxp*(dym*dvdz101+dyp*dvdz111)
		+dx3m*(dym*dvdxz001+dyp*dvdxz011)
		+dx3p*(dym*dvdxz101+dyp*dvdxz111)
		+dxm*(dy3m*dvdyz001+dy3p*dvdyz011)
		+dxp*(dy3m*dvdyz101+dy3p*dvdyz111)
		+dx3m*(dy3m*dvdxyz001+dy3p*dvdxyz011)
		+dx3p*(dy3m*dvdxyz101+dy3p*dvdxyz111));

	dvaldx=
		dzm
		*(
		-(dym*v000+dyp*v010)*invbinwx
		+(dym*v100+dyp*v110)*invbinwx
		-(3*dxm*dxm-1)*binwx_over_6*(dym*dvdx000+dyp*dvdx010)
		+(3*dxp*dxp-1)*binwx_over_6*(dym*dvdx100+dyp*dvdx110)
		-(dy3m*dvdy000+dy3p*dvdy010)*invbinwx
		+(dy3m*dvdy100+dy3p*dvdy110)*invbinwx
		-(3*dxm*dxm-1)*binwx_over_6*(dy3m*dvdxy000+dy3p*dvdxy010)
		+(3*dxp*dxp-1)*binwx_over_6*(dy3m*dvdxy100+dy3p*dvdxy110)
		)

		+dzp
		*(
		-(dym*v001+dyp*v011)*invbinwx
		+(dym*v101+dyp*v111)*invbinwx
		-(3*dxm*dxm-1)*binwx_over_6*(dym*dvdx001+dyp*dvdx011)
		+(3*dxp*dxp-1)*binwx_over_6*(dym*dvdx101+dyp*dvdx111)
		-(dy3m*dvdy001+dy3p*dvdy011)*invbinwx
		+(dy3m*dvdy101+dy3p*dvdy111)*invbinwx
		-(3*dxm*dxm-1)*binwx_over_6*(dy3m*dvdxy001+dy3p*dvdxy011)
		+(3*dxp*dxp-1)*binwx_over_6*(dy3m*dvdxy101+dy3p*dvdxy111)
		)

		+dz3m
		*(
		-(dym*dvdz000+dyp*dvdz010)*invbinwx
		+(dym*dvdz100+dyp*dvdz110)*invbinwx
		-(3*dxm*dxm-1)*binwx_over_6*(dym*dvdxz000+dyp*dvdxz010)
		+(3*dxp*dxp-1)*binwx_over_6*(dym*dvdxz100+dyp*dvdxz110)
		-(dy3m*dvdyz000+dy3p*dvdyz010)*invbinwx
		+(dy3m*dvdyz100+dy3p*dvdyz110)*invbinwx
		-(3*dxm*dxm-1)*binwx_over_6*(dy3m*dvdxyz000+dy3p*dvdxyz010)
		+(3*dxp*dxp-1)*binwx_over_6*(dy3m*dvdxyz100+dy3p*dvdxyz110)
		)

		+dz3p
		*(
		-(dym*dvdz001+dyp*dvdz011)*invbinwx
		+(dym*dvdz101+dyp*dvdz111)*invbinwx
		-(3*dxm*dxm-1)*binwx_over_6*(dym*dvdxz001+dyp*dvdxz011)
		+(3*dxp*dxp-1)*binwx_over_6*(dym*dvdxz101+dyp*dvdxz111)
		-(dy3m*dvdyz001+dy3p*dvdyz011)*invbinwx
		+(dy3m*dvdyz101+dy3p*dvdyz111)*invbinwx
		-(3*dxm*dxm-1)*binwx_over_6*(dy3m*dvdxyz001+dy3p*dvdxyz011)
		+(3*dxp*dxp-1)*binwx_over_6*(dy3m*dvdxyz101+dy3p*dvdxyz111));

	dvaldy=
		dzm
		*(
		dxm*(-v000+v010)*invbinwy
		+dxp*(-v100+v110)*invbinwy
		+dx3m*(-dvdx000+dvdx010)*invbinwy
		+dx3p*(-dvdx100+dvdx110)*invbinwy
		+dxm*(-(3*dym*dym-1)*dvdy000+(3*dyp*dyp-1)*dvdy010)*binwy_over_6
		+dxp*(-(3*dym*dym-1)*dvdy100+(3*dyp*dyp-1)*dvdy110)*binwy_over_6
		+dx3m*(-(3*dym*dym-1)*dvdxy000+(3*dyp*dyp-1)*dvdxy010)*binwy_over_6
		+dx3p*(-(3*dym*dym-1)*dvdxy100+(3*dyp*dyp-1)*dvdxy110)*binwy_over_6
		)

		+dzp
		*(
		dxm*(-v001+v011)*invbinwy
		+dxp*(-v101+v111)*invbinwy
		+dx3m*(-dvdx001+dvdx011)*invbinwy
		+dx3p*(-dvdx101+dvdx111)*invbinwy
		+dxm*(-(3*dym*dym-1)*dvdy001+(3*dyp*dyp-1)*dvdy011)*binwy_over_6
		+dxp*(-(3*dym*dym-1)*dvdy101+(3*dyp*dyp-1)*dvdy111)*binwy_over_6
		+dx3m*(-(3*dym*dym-1)*dvdxy001+(3*dyp*dyp-1)*dvdxy011)*binwy_over_6
		+dx3p*(-(3*dym*dym-1)*dvdxy101+(3*dyp*dyp-1)*dvdxy111)*binwy_over_6
		)

		+dz3m
		*(
		dxm*(-dvdz000+dvdz010)*invbinwy
		+dxp*(-dvdz100+dvdz110)*invbinwy
		+dx3m*(-dvdxz000+dvdxz010)*invbinwy
		+dx3p*(-dvdxz100+dvdxz110)*invbinwy
		+dxm*(-(3*dym*dym-1)*dvdyz000+(3*dyp*dyp-1)*dvdyz010)*binwy_over_6
		+dxp*(-(3*dym*dym-1)*dvdyz100+(3*dyp*dyp-1)*dvdyz110)*binwy_over_6
		+dx3m*(-(3*dym*dym-1)*dvdxyz000+(3*dyp*dyp-1)*dvdxyz010)*binwy_over_6
		+dx3p*(-(3*dym*dym-1)*dvdxyz100+(3*dyp*dyp-1)*dvdxyz110)*binwy_over_6
		)

		+dz3p
		*(
		dxm*(-dvdz001+dvdz011)*invbinwy
		+dxp*(-dvdz101+dvdz111)*invbinwy
		+dx3m*(-dvdxz001+dvdxz011)*invbinwy
		+dx3p*(-dvdxz101+dvdxz111)*invbinwy
		+dxm*(-(3*dym*dym-1)*dvdyz001+(3*dyp*dyp-1)*dvdyz011)*binwy_over_6
		+dxp*(-(3*dym*dym-1)*dvdyz101+(3*dyp*dyp-1)*dvdyz111)*binwy_over_6
		+dx3m*(-(3*dym*dym-1)*dvdxyz001+(3*dyp*dyp-1)*dvdxyz011)*binwy_over_6
		+dx3p*(-(3*dym*dym-1)*dvdxyz101+(3*dyp*dyp-1)*dvdxyz111)*binwy_over_6);

	dvaldz=
		-(dxm*(dym*v000+dyp*v010)
		+dxp*(dym*v100+dyp*v110)
		+dx3m*(dym*dvdx000+dyp*dvdx010)
		+dx3p*(dym*dvdx100+dyp*dvdx110)
		+dxm*(dy3m*dvdy000+dy3p*dvdy010)
		+dxp*(dy3m*dvdy100+dy3p*dvdy110)
		+dx3m*(dy3m*dvdxy000+dy3p*dvdxy010)
		+dx3p*(dy3m*dvdxy100+dy3p*dvdxy110))
		*invbinwz

		+(dxm*(dym*v001+dyp*v011)
		+dxp*(dym*v101+dyp*v111)
		+dx3m*(dym*dvdx001+dyp*dvdx011)
		+dx3p*(dym*dvdx101+dyp*dvdx111)
		+dxm*(dy3m*dvdy001+dy3p*dvdy011)
		+dxp*(dy3m*dvdy101+dy3p*dvdy111)
		+dx3m*(dy3m*dvdxy001+dy3p*dvdxy011)
		+dx3p*(dy3m*dvdxy101+dy3p*dvdxy111))
		*invbinwz

		-(3*dzm*dzm-1)*binwz_over_6
		*(dxm*(dym*dvdz000+dyp*dvdz010)
		+dxp*(dym*dvdz100+dyp*dvdz110)
		+dx3m*(dym*dvdxz000+dyp*dvdxz010)
		+dx3p*(dym*dvdxz100+dyp*dvdxz110)
		+dxm*(dy3m*dvdyz000+dy3p*dvdyz010)
		+dxp*(dy3m*dvdyz100+dy3p*dvdyz110)
		+dx3m*(dy3m*dvdxyz000+dy3p*dvdxyz010)
		+dx3p*(dy3m*dvdxyz100+dy3p*dvdxyz110))

		+(3*dzp*dzp-1)*binwz_over_6
		*(dxm*(dym*dvdz001+dyp*dvdz011)
		+dxp*(dym*dvdz101+dyp*dvdz111)
		+dx3m*(dym*dvdxz001+dyp*dvdxz011)
		+dx3p*(dym*dvdxz101+dyp*dvdxz111)
		+dxm*(dy3m*dvdyz001+dy3p*dvdyz011)
		+dxp*(dy3m*dvdyz101+dy3p*dvdyz111)
		+dx3m*(dy3m*dvdxyz001+dy3p*dvdxyz011)
		+dx3p*(dy3m*dvdxyz101+dy3p*dvdxyz111));

}

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
			if (ii == 1) { 	//change first chi in place:
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
