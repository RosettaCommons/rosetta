// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/fiber_diffractiobn/hankel_kiss_fft.hh
/// @brief  Fast Hankel transform needed for low-resolution scoring of fiber diffraction data
/// @author Wojciech Potrzebowski and Ingemar Andre

#ifndef INCLUDED_core_scoring_fiber_diffraction_hankel_kiss_hh
#define INCLUDED_core_scoring_fiber_diffraction_hankel_kiss_hh

#include <ObjexxFCL/FArray3D.fwd.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <core/types.hh>
#include <complex>
#include <numeric/fourier/kiss_fft.hh>
#include <numeric/fourier/FFT.hh>
#include <utility/vector0.hh>

namespace core {
namespace scoring {
namespace fiber_diffraction {

//TODO:Can be revritten to class
typedef struct Hankel{
	core::Size n;
	core::Real k1;
	core::Real k2;
	core::Real rp0;
	core::Real alpha;
	core::Real *f;
	core::Real *j;
	core::Real *rp;
	//utility::vector0< core::Real > rp_vec;
	int lec_order;
	core::Real *f0;
	int l;
	core::Real *snn;
	std::complex<double> * data_in_fft;
} Hankel;

void hankel_free (
	Hankel *p_hankel
);

void hankel_trans_no_lec (
	Hankel *p_hankel
);

Hankel * hankel_make_input (
	core::Size length,
	core::Real k1,
	core::Real k2,
	core::Real b,
	int lec_order,
	core::Real *f0,
	ObjexxFCL::FArray3D < std::complex<float> >  & hankel_in,
	core::Size & lindex,
	core::Size & nindex,
	int l
);

core::Real hankel_get_rc (
	Hankel *p_hankel
);

void c_array_mult (
	core::Size length,
	core::Real *dp_1,
	core::Real *dp_2
);

void d_array_scale (
	core::Size length,
	core::Real factor,
	core::Real *dp_in
);

void set_r_array (
	core::Size num_r_points,
	core::Real k1,
	core::Real k2,
	core::Real max_r,
	ObjexxFCL::FArray1D<float> & rc
);

void set_r_inv_array (
	core::Size num_r_points,
	core::Real k1,
	core::Real k2,
	core::Real max_r,
	ObjexxFCL::FArray1D<float> & Rinv
);

core::Real alpha_func(
	core::Size n,
	core::Real k1,
	core::Real k2,
	core::Real alpha
);

core::Real alpha_deriv_func(
	core::Size n,
	core::Real alpha
);


void hankel_set_alpha(
	Hankel *p_hankel
);

void c_mult_ip(
	core::Real *one ,
	const core::Real *two
);

void hankel_r_mult(
	Hankel *p_hankel
);

void hankel_r_div(
	Hankel *p_hankel
);

core::Real hankel_get_p0(
	Hankel *p_hankel
);

void hankel_make_snn(
	Hankel *p_hankel
);

void hankel_in_machine(
	Hankel *p_hankel
);

core::Real hankel_get_rc(
	Hankel *p_hankel
);

void hankel_make_j(
	Hankel *p_hankel
);

void hankel_make_rp (
	Hankel *p_hankel
);

void dfour1_plan(
	double *data,
	core::Size nn,
	int isign,
	std::complex<double> *in
	//numeric::fourier::kiss_fft_cpx *in
);

core::Real NRbisafe(
	core::Size n,
	core::Real k1,
	core::Real k2,
	core::Real X1,
	core::Real X2 );
}
}
}
#endif
