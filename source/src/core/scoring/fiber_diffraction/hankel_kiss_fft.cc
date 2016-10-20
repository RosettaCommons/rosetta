// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/fiber_diffractiobn/hankel_kiss_fft.cc
/// @brief  Fast Hankel transform needed for low-resolution scoring of fiber diffraction data
/// @author Wojciech Potrzebowski and Ingemar Andre


#include <core/scoring/fiber_diffraction/hankel_kiss_fft.hh>

#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/FArray1D.hh>

#include <basic/Tracer.hh>
#include <utility/basic_sys_util.hh>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>

#ifdef WIN32
#define _USE_MATH_DEFINES
#include <math.h>
#endif

static basic::Tracer TR("core.scoring.fiber_diffraction.hankel_kiss_fft");

namespace core {
namespace scoring {
namespace fiber_diffraction {

core::Real alpha_func (
	core::Size n,
	core::Real k1,
	core::Real k2,
	core::Real alpha
){
	debug_assert( k2 != 0 );
	return alpha * exp( n * alpha ) - k1 / k2;
}

core::Real  alpha_deriv_func (
	core::Size n,
	core::Real alpha
){
	return ( 1 + n * alpha ) * exp( n * alpha );
}

void hankel_set_alpha
(
	Hankel *p_hankel
){
	debug_assert( p_hankel->n != 0);
	core::Real X_MIN(0.0);
#ifdef WIN32
	p_hankel->alpha = NRbisafe(p_hankel->n, p_hankel->k1, p_hankel->k2, X_MIN, 30 * log(10.) / p_hankel->n);
#else
	p_hankel->alpha = NRbisafe(p_hankel->n, p_hankel->k1, p_hankel->k2, X_MIN, 30 * log(10) / p_hankel->n);
#endif
}

void c_array_mult
(
	core::Size length,
	core::Real *dp_1,
	core::Real *dp_2
){
	for ( core::Size count = 0 ; count < length ; ++count ) {
		c_mult_ip( dp_1 , dp_2 );
		dp_1 += 2;
		dp_2 += 2;
	}
}

void hankel_free (
	Hankel *p_hankel
){
	free( p_hankel->f );
	free( p_hankel->f0 );
	free( p_hankel->j );
	free( p_hankel->rp );
	if ( p_hankel->snn != NULL ) {
		free( p_hankel->snn );
	}
	free( p_hankel->data_in_fft );
	free( p_hankel );
}


void c_mult_ip (
	core::Real *one,
	core::Real *two
){
	core::Real temp;
	temp = *one * *two - *( one + 1 ) * *( two + 1 );
	*( one + 1 ) = *one * *( two + 1 ) + *( one + 1 ) * *two;
	*one = temp;
}


void hankel_r_mult (
	Hankel *p_hankel
){
	core::Real *f;
	core::Real *rp;
	f = p_hankel->f;
	rp = p_hankel->rp;

	for ( core::Size count = 0 ; count < p_hankel->n ; count++ ) {
		*f++ *= p_hankel->rp0 * *rp;
		*f++ *= p_hankel->rp0 * *rp;
		rp++;
	}
}

void hankel_r_div
(
	Hankel *p_hankel
){
	core::Real *f;
	core::Real *rp;
	f = p_hankel->f;
	rp = p_hankel->rp;

	for ( core::Size count = 0 ; count < p_hankel->n ; count++ ) {
		*f++ /= p_hankel->rp0 * *rp;
		*f++ /= p_hankel->rp0 * *rp;
		rp++;
	}
}


void hankel_trans_no_lec (
	Hankel *p_hankel
){
	hankel_r_mult( p_hankel );
	dfour1_plan( p_hankel->f - 1 , 2 * p_hankel->n , 1 , p_hankel->data_in_fft);
	c_array_mult( 2 * p_hankel->n , p_hankel->f , p_hankel->j );
	dfour1_plan( p_hankel->f - 1 , 2 * p_hankel->n , 1 , p_hankel->data_in_fft);
	memset( p_hankel->f + 2 * p_hankel->n , 0 , 2 * p_hankel->n *sizeof( core::Real ) );
	p_hankel->rp0 = hankel_get_p0( p_hankel );
	hankel_r_div( p_hankel );
}


core::Real hankel_get_p0 (
	Hankel *p_hankel
){
	debug_assert( p_hankel->k1 != 0 && p_hankel->rp0 != 0 );
	return ( p_hankel->k2 * p_hankel->alpha /
		( p_hankel->k1 * p_hankel->k1 * p_hankel->rp0 ) );
}


void hankel_make_snn (
	Hankel *p_hankel
){
	core::Real *snn;
	core::Real factor;
	core::Real pre;
	if ( p_hankel->n == 0 ) {
		utility_exit_with_message( "Empty malloc in hankel_copy" );
	}
	snn = p_hankel->snn = ( core::Real * ) malloc( p_hankel->n *
		sizeof( core::Real ) );
	if ( p_hankel->snn == NULL ) {
		utility_exit_with_message( "snn malloc failed in hankel_copy" );
	}

	pre = pow( hankel_get_rc( p_hankel ) , p_hankel->l + 1 );

	factor = 2 * M_PI * hankel_get_p0( p_hankel ) * hankel_get_rc( p_hankel );

	for ( core::Size count = 0 ; count < p_hankel->n ; count++ ) {
		*snn++ = pre * jn( p_hankel->l + 1, factor * exp( p_hankel->alpha * count ) );
	}
}


void hankel_in_machine (
	Hankel *p_hankel
){
	hankel_make_j( p_hankel );
	hankel_make_rp( p_hankel );
	if ( p_hankel->l <= p_hankel->lec_order ) {
		hankel_make_snn( p_hankel );
	} else {
		p_hankel->snn = NULL;
	}
	memset( p_hankel->f + 2 * p_hankel->n , 0 ,
		2 * p_hankel->n * sizeof( core::Real ) );
}


core::Real hankel_get_rc (
	Hankel *p_hankel
){
	return p_hankel->rp0 * exp(-p_hankel->alpha / 2);
}


void hankel_make_j (
	Hankel *p_hankel
){
	core::Size count_max;
	core::Real *j;
	core::Real factor;
	core::Real arg;
	if ( p_hankel->n == 0 ) {
		utility_exit_with_message( "Empty malloc in hankel_copy" );
	}
	j = p_hankel->j = ( core::Real * ) malloc( 4 * p_hankel->n *
		sizeof( core::Real ) );
	if ( p_hankel->j == NULL ) {
		utility_exit_with_message( "malloc failed in make_j");
	}
	factor = 2 * M_PI * p_hankel->k2 * p_hankel->alpha /
		( p_hankel->k1 * p_hankel->k1 );
	count_max = 2 * p_hankel->n;
	for ( core::Size count = 0 ; count < count_max ; count++ ) {
		arg = factor * exp( p_hankel->alpha * count );
		*j++ = p_hankel->alpha * arg * jn( p_hankel->l, arg);
		*j++ = 0;
	}
	dfour1_plan( p_hankel->j - 1 , 2 * p_hankel->n , -1, p_hankel->data_in_fft );
	d_array_scale( 4 * p_hankel->n , 1.0 / ( 2 * p_hankel->n ), p_hankel->j );
}


void hankel_make_rp (
	Hankel *p_hankel
){
	core::Real *rp;

	rp = p_hankel->rp = ( core::Real * ) malloc( p_hankel->n * sizeof( core::Real ) );
	if ( p_hankel->rp == NULL ) {
		utility_exit_with_message( "fourth malloc failed in top_hat" );
	}
	for ( core::Size count = 0 ; count < p_hankel->n ;
			*rp++ = exp( p_hankel->alpha * count++ ) ) ; // fpd is this semicolon correct?
}


Hankel * hankel_make_input (
	core::Size length,
	core::Real k1,
	core::Real k2,
	core::Real b,
	int lec_order,
	core::Real *f0,
	ObjexxFCL::FArray3D < std::complex<float> > & hankel_in,
	core::Size & lindex,
	core::Size & nindex,
	int l
){
	core::Real *f;
	//core::Real r;
	Hankel *p_hankel;
	p_hankel = ( Hankel * ) malloc( sizeof( Hankel ) );
	if ( p_hankel == NULL ) {
		utility_exit_with_message( "first malloc failed in make_input" );
	}

	p_hankel->n = length;
	p_hankel->k1 = k1;
	p_hankel->k2 = k2;
	p_hankel->lec_order = lec_order;
	p_hankel->l = l;
	f = p_hankel->f = ( core::Real * ) malloc( 4 * p_hankel->n *  sizeof( core::Real ) );

	if ( p_hankel->f == NULL ) {
		utility_exit_with_message( "Second malloc failed in make_input" );
	}

	p_hankel->f0 = ( core::Real * ) malloc( 2 * p_hankel->lec_order *
		sizeof( core::Real ) );
	if ( p_hankel->f0 == NULL ) {
		utility_exit_with_message( "Third malloc failed in make_input" );
	}

	memcpy( p_hankel->f0 , f0 , 2 * p_hankel->lec_order * sizeof( core::Real ) );

	p_hankel->data_in_fft = ( std::complex<double> *) malloc(sizeof(std::complex<double> ) * 2 * p_hankel->n );
	if ( p_hankel->data_in_fft == NULL ) {
		utility_exit_with_message( "Fourth malloc failed in make_input" );
	}
	hankel_set_alpha( p_hankel );
	if ( b == 0 ) {
		p_hankel->rp0 = sqrt( p_hankel->k2 * p_hankel->alpha) / p_hankel->k1;
	} else {
		p_hankel->rp0 = b * exp( -p_hankel->alpha * p_hankel->n );
	}

	for ( core::Size count = 0 ; count < p_hankel->n ; count++ ) {
		//r = p_hankel->rp0 * exp( p_hankel->alpha * count );
		*f = hankel_in(count+1,lindex,nindex).real();
		*(f+1) = hankel_in(count+1,lindex,nindex).imag();
		f += 2;
	}

	hankel_in_machine( p_hankel);

	return p_hankel;
}

void d_array_scale (
	core::Size length,
	core::Real factor,
	core::Real *dp_in
){
	core::Size count;
	for ( count = 0 ; count < length ; count++ ) {
		*dp_in++ *= factor;
	}
}

void set_r_array (
	core::Size num_r_points,
	core::Real k1,
	core::Real k2,
	core::Real max_r,
	ObjexxFCL::FArray1D<float> & rc
) {
	core::Real X_MIN(0.0);
#ifdef WIN32
	core::Real alpha = NRbisafe(num_r_points, k1, k2, X_MIN , 30 * log(10.) / num_r_points);
#else
	core::Real alpha = NRbisafe(num_r_points, k1, k2, X_MIN , 30 * log(10) / num_r_points);
#endif
	core::Real rp0 = max_r * exp( -alpha * num_r_points );

	debug_assert( num_r_points <= rc.size() );
	for ( core::Size count = 0 ; count < num_r_points; count++ ) {
		core::Real r = rp0 * exp( alpha * count );
		rc(count+1) = r;
	}
}

void set_r_inv_array (
	core::Size num_r_points,
	core::Real k1,
	core::Real k2,
	core::Real max_r,
	ObjexxFCL::FArray1D<float> & Rinv
) {
	core::Real X_MIN(0.0);
#ifdef WIN32
	core::Real alpha = NRbisafe(num_r_points, k1, k2, X_MIN , 30 * log(10.) / num_r_points );
#else
	core::Real alpha = NRbisafe(num_r_points, k1, k2, X_MIN , 30 * log(10) / num_r_points );
#endif
	core::Real rp0 = max_r * exp( -alpha * num_r_points );
	core::Real Rp0_new ( k2 * alpha /
		( k1 * k1 * rp0 ) );

	debug_assert( num_r_points <= Rinv.size() );
	for ( core::Size count = 0 ; count < num_r_points; count++ ) {
		core::Real R = Rp0_new * exp( alpha * count );
		Rinv(count+1) = R;
	}
}

void dfour1_plan(double *data, core::Size nn, int isign, std::complex<double> *in)
{
	double *f;
	f=data+1;
	for ( core::Size i=0; i<nn; ++i ) {
		std::complex<double> cmpx( *f , *( f+1 ) );
		in[i] = cmpx;
		f+=2;
	}
	if ( isign == -1 ) {
		numeric::fourier::kiss_fft_state fft_params;
		fft_params.resize( nn, 0 );
		numeric::fourier::kiss_fft(&fft_params, in, in );
	}
	if ( isign == 1  ) {
		numeric::fourier::kiss_fft_state ifft_params;
		ifft_params.resize( nn, 1 );
		numeric::fourier::kiss_fft(&ifft_params, in, in );
	}
	f=data+1;
	for ( core::Size i=0; i< nn; ++i ) {
		*f = in[i].real();
		*( f + 1 ) =  in[i].imag();
		f+=2;
	}
}

core::Real NRbisafe(core::Size n, core::Real k1, core::Real k2, core::Real X1, core::Real X2 ) {
	core::Real ACC( 1e-9);
	core::Size max_iterations(100);
	core::Real dF,dX,dXold,F,Fh,Fl;
	core::Real temp,Xh,Xl,rts;

	Fl = alpha_func( n, k1, k2, X1 );
	Fh = alpha_func( n, k1, k2, X2 );

	if ( Fl*Fh >= 0.0 ) {
		utility_exit_with_message( "Cannot derive roots of alpha function" );
	}
	if ( Fl < 0.0 ) {
		Xl=X1;
		Xh=X2;
	} else {
		Xh=X1;
		Xl=X2;
	}
	rts=0.5*(X1+X2);
	dXold=fabs(X2-X1);
	dX=dXold;

	F = alpha_func( n, k1, k2, rts );
	dF = alpha_deriv_func( n, rts );

	for ( core::Size j=1; j<=max_iterations; j++ ) {
		if ( (((rts-Xh)*dF-F)*((rts-Xl)*dF-F) >= 0.0) || (fabs(2.0*F) > fabs(dXold*dF)) ) {
			dXold=dX;
			dX=0.5*(Xh-Xl);
			rts=Xl+dX;
			if ( Xl == rts ) return rts;
		} else {
			dXold=dX;
			dX=F/dF;
			temp=rts;
			rts -= dX;
			if ( temp == rts ) return rts;
		}
		if ( fabs(dX) < ACC ) return rts;
		F = alpha_func( n, k1, k2, rts );
		dF = alpha_deriv_func( n, rts );

		if ( F < 0.0 ) Xl=rts;
		else Xh=rts;
	}
	utility_exit_with_message("Newton-Raphson algorithm didn't converge in 100 iterations");
}

}
}
}
