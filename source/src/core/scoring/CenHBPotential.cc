// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/CenHBPotential.cc
/// @brief  Smooth, differentiable version of centroid hbond term
/// @author Frank DiMaio


#include <core/scoring/CenHBPotential.hh>

#include <basic/database/open.hh>
#include <utility/io/izstream.hh>

#include <utility/vector1.hh>
#include <boost/bind.hpp>


namespace core {
namespace scoring {

///////////////////////////////////////////////////////////////////////////////////////////////

Real CenHBPotential::func( Size seqsep, Real d, Real xd, Real xh ) const {
	using numeric::constants::f::pi;

	utility::vector1< Real > const & As       = (seqsep<=4) ? sr_As_ : lr_As_;
	utility::vector1< Vector > const & mus    = (seqsep<=4) ? sr_mus_ : lr_mus_;
	utility::vector1< Vector > const & sigmas = (seqsep<=4) ? sr_sigmas_ : lr_sigmas_;

	Real y = 0;
	for ( Size i=1; i<=As.size(); ++i ) {
		y += As[i] * exp( -(d-mus[i][0])*(d-mus[i][0]) / (2*sigmas[i][0]*sigmas[i][0]) )
			* ( exp( sigmas[i][1]*cos( pi/180*(xd-mus[i][1]) ) ) + exp( sigmas[i][1]*cos( pi/180*(xd+mus[i][1]) ) ) )
			* ( exp( sigmas[i][2]*cos( pi/180*(xh-mus[i][2]) ) ) + exp( sigmas[i][2]*cos( pi/180*(xh+mus[i][2]) ) ) );
	}
	return y;
}

Vector CenHBPotential::dfunc( Size seqsep, Real d, Real xd, Real xh ) const {
	using numeric::constants::f::pi;

	utility::vector1< Real > const & As       = (seqsep<=4) ? sr_As_ : lr_As_;
	utility::vector1< Vector > const & mus    = (seqsep<=4) ? sr_mus_ : lr_mus_;
	utility::vector1< Vector > const & sigmas = (seqsep<=4) ? sr_sigmas_ : lr_sigmas_;

	Vector dy = Vector(0,0,0);
	for ( Size i=1; i<=As.size(); ++i ) {
		Real s1 = exp( -(d-mus[i][0])*(d-mus[i][0]) / (2*sigmas[i][0]*sigmas[i][0]) );
		Real s2a = exp( sigmas[i][1]*cos( pi/180*(xd-mus[i][1]) ) );
		Real s2b = exp( sigmas[i][1]*cos( pi/180*(xd+mus[i][1]) ) );
		Real s3a = exp( sigmas[i][2]*cos( pi/180*(xh-mus[i][2]) ) );
		Real s3b = exp( sigmas[i][2]*cos( pi/180*(xh+mus[i][2]) ) );
		Real s2 = s2a+s2b;
		Real s3 = s3a+s3b;

		dy[0] += -As[i] * s1 * s2 * s3 * (d-mus[i][0]) / (sigmas[i][0]*sigmas[i][0]);
		dy[1] += -As[i] * s1 * s3 * pi/180 * sigmas[i][1] * ( sin( pi/180*(xd-mus[i][1])) * s2a + sin( pi/180*(xd+mus[i][1])) * s2b );
		dy[2] += -As[i] * s1 * s2 * pi/180 * sigmas[i][2] * ( sin( pi/180*(xh-mus[i][2])) * s3a + sin( pi/180*(xh+mus[i][2])) * s3b );
	}
	return dy;
}

///////////////////////////////////////////////////////////////////////////////////////////////

CenHBPotential::CenHBPotential() {
	// defaults
	cutoff_sr_ = 6.0;
	cutoff_lr_ = 6.0;

	// load the data
	std::string tag,line;
	utility::io::izstream stream;
	basic::database::open( stream, "scoring/score_functions/centroid_smooth/cen_hb_params.txt");

	while ( getline( stream, line ) ) {
		std::istringstream l(line);
		l >> tag;
		if ( tag == "HBOND_BB_LR:" ) {
			l >> tag;
			if ( tag == "CUTOFF" ) {
				l >> cutoff_lr_;
			} else if ( tag == "GAUSSIAN3D" ) {
				Size ngauss; l >> ngauss;
				Real A;
				numeric::xyzVector<Real> mu, sigma;
				for ( Size i=1; i<=ngauss; ++i ) {
					l >> A >> mu[0] >> mu[1] >> mu[2] >> sigma[0] >> sigma[1] >> sigma[2];
					add_lr_gaussian(A,mu,sigma);
				}
			}
		} else if ( tag == "HBOND_BB_SR:" ) {
			l >> tag;
			if ( tag == "CUTOFF" ) {
				l >> cutoff_sr_;
			} else if ( tag == "GAUSSIAN3D" ) {
				Size ngauss; l >> ngauss;
				Real A;
				numeric::xyzVector<Real> mu, sigma;
				for ( Size i=1; i<=ngauss; ++i ) {
					l >> A >> mu[0] >> mu[1] >> mu[2] >> sigma[0] >> sigma[1] >> sigma[2];
					add_sr_gaussian(A,mu,sigma);
				}
			}
		} else if ( tag != "#" ) {
			utility_exit_with_message("bad format for cen_smooth_params.txt");
		}

		if ( l.fail() ) utility_exit_with_message("bad format for cen_smooth_params.txt");
	}
}

// implementation of soft version
Real
CenHBPotential::func_soft( Vector a1, Vector a2, Vector b1, Vector b2, Vector dv ) const
{
	using numeric::constants::f::pi;

	Real f, fa, fb, fd, fdv1, fdv2, dota, dotb, dotdv1, dotdv2, d;
	Real const dv0( sin(pi/12.0) ), dvs( sin(pi/6.0) );

	d = dv.length();
	fd = fade( d, 6.0, 3.0, true ); // 1 at < 6.0, fade until 9.0 Ang

	// normalize
	dv.normalize();
	a1.normalize();
	a2.normalize();
	b1.normalize();
	b2.normalize();

	dota = a1.dot( a2 );
	dotb = b1.dot( b2 );
	dotdv1 = dv.dot( a1 );
	dotdv2 = dv.dot( b1 );

	fa  = fade(       dota,    0.75,    0.75, false  ); // 1 at < 30', fade until 90'
	fb  = fade(  dotb*dotb,    0.75,    0.50, false ); // 1 at > 150 or < -150', fade until 120/-120'
	fdv1 = fade(dotdv1*dotdv1, dv0*dv0, dvs*dvs-dv0*dv0, true ); // 1 at 75~105', fade until 60~120'
	fdv2 = fade(dotdv2*dotdv2, dv0*dv0, dvs*dvs-dv0*dv0, true ); // 1 at 75~105', fade until 60~120'
	f = -fd*fa*fb*fdv1*fdv2;
	return f;
}

void
CenHBPotential::dfunc_soft( Vector a1, Vector a2, Vector b1, Vector b2, Vector dv,
	utility::vector1< Vector > &df_dABNC_1,
	utility::vector1< Vector > &df_dABNC_2
) const
{
	using numeric::constants::f::pi;

	//utility::vector1< Vector > df_dABNC_1( 4, Vector(0.0) ); // in order of A1, B1, N1, C1
	//utility::vector1< Vector > df_dABNC_2( 4, Vector(0.0) ); // in order of A2, B2, N2, C2

	Real f, fa, fb, fd, fdv1, fdv2;
	Real dota, dotb, dotdv1, dotdv2, d, a1norm, a2norm, b1norm, b2norm;
	Real dfd_dd, dfa_ddota, dfb_ddotb, dfdv1_ddotdv1, dfdv2_ddotdv2;
	Real const dv0( sin(pi/12.0) ), dvs( sin(pi/6.0) );

	d = dv.length();
	a1norm = a1.length();
	a2norm = a2.length();
	b1norm = b1.length();
	b2norm = b2.length();

	// normalize
	dv /= d;
	a1 /= a1norm;
	a2 /= a2norm;
	b1 /= b1norm;
	b2 /= b2norm;

	dota = a1.dot( a2 );
	dotb = b1.dot( b2 );
	dotdv1 = dv.dot( a1 );
	dotdv2 = dv.dot( b1 );

	fd = fade( d, 6.0, 3.0, true ); // 1 at < 6.0, fade until 9.0 Ang
	fa  = fade(       dota,    0.75,    0.75, false  ); // 1 at < 30', fade until 90'
	fb  = fade(  dotb*dotb,    0.75,    0.50, false ); // 1 at > 150 or < -150', fade until 120/-120'
	fdv1 = fade(dotdv1*dotdv1, dv0*dv0, dvs*dvs-dv0*dv0, true ); // 1 at 75~105', fade until 60~120'
	fdv2 = fade(dotdv2*dotdv2, dv0*dv0, dvs*dvs-dv0*dv0, true ); // 1 at 75~105', fade until 60~ dfd_dd     = dfade( d, 6.0, 3.0, true ); // 1 at < 6.0, fade until 9.0 Ang

	f = -fd*fa*fb*fdv1*fdv2;

	if ( f > -1e-6 ) return;

	dfd_dd     = dfade( d, 6.0, 3.0, true );
	dfa_ddota  = dfade(       dota,    0.75,    0.75, false  );
	dfb_ddotb  = dfade(  dotb*dotb,    0.75,    0.50, false );
	dfdv1_ddotdv1 = dfade(dotdv1*dotdv1, dv0*dv0, dvs*dvs-dv0*dv0, true );
	dfdv2_ddotdv2 = dfade(dotdv2*dotdv2, dv0*dv0, dvs*dvs-dv0*dv0, true );

	dfb_ddotb     *= 2*dotb;
	dfdv1_ddotdv1 *= 2*dotdv1;
	dfdv2_ddotdv2 *= 2*dotdv2;

	Vector dd_dA1 = -dv; // normalized
	Vector dd_dA2 = -dd_dA1;

	Vector ddota_dA1 = (-a2 + dota*a1)/a1norm;
	Vector ddota_dA2 = (-a1 + dota*a2)/a2norm;
	//Vector ddota_dB1 = -ddota_dA1;
	//Vector ddota_dB2 = -ddota_dB2;

	Vector ddotb_dN1 = (-b2 + dotb*b1)/b1norm;
	Vector ddotb_dN2 = (-b1 + dotb*b2)/b2norm;
	//Vector ddotb_dC1 = -ddotb_dN1;
	//Vector ddotb_dC2 = -ddotb_dN2;

	//dv1 = dv_dot_a1
	Vector ddotdv1_dA2 = (a1 - dotdv1*a1)/d;
	Vector ddotdv1_dB1 = (dv - dotdv1*dv)/a1norm;
	//Vector ddotdv1_dA1 = -ddotdv1_dA2 - ddotdv1_dB1;

	// dv2 = dv_dot_b1
	Vector ddotdv2_dN1 = (-dv + dotdv2*b1)/b1norm;
	Vector ddotdv2_dA1 = (-b1 + dotdv2*dv)/d;
	//Vector ddotdv2_dC1 = -ddotdv2_dN1;
	//Vector ddotdv2_dA2 = -ddotdv2_dA1;

	// contribution from f(d): A1, A2
	if ( fd > 1e-6 ) {
		df_dABNC_1[1] += f/fd * dfd_dd*dd_dA1;
		df_dABNC_2[1] += f/fd * dfd_dd*dd_dA2;
	}

	// contribution from f(dota): A1, A2, B1, B2
	if ( fa > 1e-6 ) {
		df_dABNC_1[1] += f/fa * dfa_ddota*ddota_dA1;
		df_dABNC_2[1] += f/fa * dfa_ddota*ddota_dA2;
		df_dABNC_1[2] -= f/fa * dfa_ddota*ddota_dA1;
		df_dABNC_2[2] -= f/fa * dfa_ddota*ddota_dA2;
	}

	// contribution from f(dotb): N1, C1, N2, C2
	if ( fb > 1e-6 ) {
		df_dABNC_1[3] += f/fb * dfb_ddotb*ddotb_dN1;
		df_dABNC_2[3] += f/fb * dfb_ddotb*ddotb_dN2;
		df_dABNC_1[4] -= f/fb * dfb_ddotb*ddotb_dN1;
		df_dABNC_2[4] -= f/fb * dfb_ddotb*ddotb_dN2;
	}

	// contribution from f(dv1): A1, B1, A2
	if ( fdv1 > 1e-6 ) {
		df_dABNC_1[2] += f/fdv1 * dfdv1_ddotdv1*ddotdv1_dB1;
		df_dABNC_2[1] += f/fdv1 * dfdv1_ddotdv1*ddotdv1_dA2;
		df_dABNC_1[1] -= f/fdv1 * dfdv1_ddotdv1*(ddotdv1_dB1 + ddotdv1_dA2);
	}

	// contribution from f(dv2): A1, A2, N1, C1
	if ( fdv2 > 1e-6 ) {
		df_dABNC_1[1] += f/fdv2 * dfdv2_ddotdv2*ddotdv2_dA1;
		df_dABNC_2[1] -= f/fdv2 * dfdv2_ddotdv2*ddotdv2_dA1;
		df_dABNC_1[3] += f/fdv2 * dfdv2_ddotdv2*ddotdv2_dN1;
		df_dABNC_1[4] -= f/fdv2 * dfdv2_ddotdv2*ddotdv2_dN1;
	}
}

}
}
