// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


/// @file   rgg.hh
/// @brief  Header file for the EISPACK rgg routine.
/// @author Kale Kundert

#ifndef INCLUDED_numeric_linear_algebra_rgg_HH
#define INCLUDED_numeric_linear_algebra_rgg_HH

// Fortran Emulation Headers
#include <fem/fem.hpp>

namespace numeric {
namespace linear_algebra {

double epslon(
	double const & x);

void qzhes(
	int const & nm,
	int const & n,
	fem::major_types::arr_ref<double,2> a,
	fem::major_types::arr_ref<double,2> b,
	bool const & matz,
	fem::major_types::arr_ref<double,2> z);

void qzit(
	int const & nm,
	int const & n,
	fem::major_types::arr_ref<double,2> a,
	fem::major_types::arr_ref<double,2> b,
	double const & eps1,
	bool const & matz,
	fem::major_types::arr_ref<double,2> z,
	int & ierr);

void qzval(
	int const & nm,
	int const & n,
	fem::major_types::arr_ref<double,2> a,
	fem::major_types::arr_ref<double,2> b,
	fem::major_types::arr_ref<double> alfr,
	fem::major_types::arr_ref<double> alfi,
	fem::major_types::arr_ref<double> beta,
	bool const & matz,
	fem::major_types::arr_ref<double,2> z);

void qzvec(
	int const & nm,
	int const & n,
	fem::major_types::arr_cref<double,2> a,
	fem::major_types::arr_ref<double,2> b,
	fem::major_types::arr_cref<double> alfr,
	fem::major_types::arr_cref<double> alfi,
	fem::major_types::arr_cref<double> beta,
	fem::major_types::arr_ref<double,2> z);

void rgg(
	int const & nm,
	int const & n,
	fem::major_types::arr_ref<double,2> a,
	fem::major_types::arr_ref<double,2> b,
	fem::major_types::arr_ref<double> alfr,
	fem::major_types::arr_ref<double> alfi,
	fem::major_types::arr_ref<double> beta,
	int const & matz,
	fem::major_types::arr_ref<double,2> z,
	int & ierr);

} // namespace linear_algebra
} // namespace numeric

#endif
