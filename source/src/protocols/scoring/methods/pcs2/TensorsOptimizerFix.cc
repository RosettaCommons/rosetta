// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

//////////////////////////////////////////////
///
/// @file protocols/scoring/methods/pcs2/TensorsOptimizerFix.cc
///
/// @brief
///
/// @details
///
/// @param
///
/// @return
///
/// @remarks
///
/// @references
///
/// @authorv Christophe Schmitz
///
////////////////////////////////////////////////


// Unit headers
#include <protocols/scoring/methods/pcs2/TensorsOptimizerFix.hh>
#include <protocols/scoring/methods/pcs2/PcsDataCenter.hh>
// Package headers

// Project headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/exit.hh>

// Numeric headers
#include <numeric/constants.hh>

// Objexx headers

// C++ headers

#include <utility/vector1.hh>


namespace protocols {
namespace scoring {
namespace methods {
namespace pcs2 {

static thread_local basic::Tracer TR_TensorsOptimizerFix( "protocols.scoring.methods.pcs.TensorsOptimizerFix" );

TensorsOptimizerFix::TensorsOptimizerFix(PcsDataCenter const & pcs_d_c/*, core::Real xM, core::Real yM, core::Real zM*/):
	pcs_d_c_(pcs_d_c)//, xM_(xM), yM_(yM), zM_(zM)
{
	/*
	core::Size i;
	const utility::vector1<core::Real> & X_all(pcs_d_c_.get_X_all());
	const utility::vector1<core::Real> & Y_all(pcs_d_c_.get_Y_all());
	const utility::vector1<core::Real> & Z_all(pcs_d_c_.get_Z_all());

	utility::vector1<core::Real> A(5, 0);

	for(i = 1; i <= X_all.size(); i++){
	fill_A_line(A, xM_, yM_, zM_, X_all[i], Y_all[i], Z_all[i]);
	Xxx_coef_vect_.push_back(A[1]);
	Xxy_coef_vect_.push_back(A[2]);
	Xxz_coef_vect_.push_back(A[3]);
	Xyy_coef_vect_.push_back(A[4]);
	Xyz_coef_vect_.push_back(A[5]);
	}
	*/
}

TensorsOptimizerFix::~TensorsOptimizerFix(){
}


core::Real
TensorsOptimizerFix::func( core::optimization::Multivec const & vars ) const{
	core::Size i, j, idx;
	core::Size n_la, n_pcs;
	// core::Real x, y, z;
	core::Real PCS_calc, PCS_exp, score, score_la;;
	utility::vector1<core::Real> A(5, 0);

	score = 0;
	n_la = pcs_d_c_.get_n_lanthanides();

	if ( n_la != (vars.size()) / 5 ) {
		TR_TensorsOptimizerFix << "n_la: " << n_la << " (vars.size()) / 5):" << (vars.size()) / 5 << " vars.size(): " << vars.size() <<std::endl;
		utility_exit_with_message("The number of lanthanides is inconsistent with the derivatives size of vars");
	}


	const utility::vector1<PcsDataLanthanide> & pcs_d_l_vec(pcs_d_c_.get_pcs_data_per_lanthanides_all());
	utility::vector1< utility::vector1<core::Real> >  const & A_all( pcs_d_c_.get_A_all());
	/*
	const utility::vector1<core::Real> & X_all(pcs_d_c_.get_X_all());
	const utility::vector1<core::Real> & Y_all(pcs_d_c_.get_Y_all());
	const utility::vector1<core::Real> & Z_all(pcs_d_c_.get_Z_all());
	*/

	if ( n_la != pcs_d_l_vec.size() ) {
		utility_exit_with_message("The number of lanthanides is inconsistent");
	}

	for ( i = 1; i <= n_la; i++ ) {
		core::Real Xxx(vars[5*(i-1) + 1]);
		core::Real Xxy(vars[5*(i-1) + 2]);
		core::Real Xxz(vars[5*(i-1) + 3]);
		core::Real Xyy(vars[5*(i-1) + 4]);
		core::Real Xyz(vars[5*(i-1) + 5]);

		PcsDataLanthanide const & pcs_d_l(pcs_d_l_vec[i]);
		utility::vector1<core::Size> const & A_index( pcs_d_l.get_A_index());
		utility::vector1< core::Real > const & cstyle_b(pcs_d_l.get_cstyle_b());

		core::Real weight(pcs_d_l.get_weight());
		//if weight = 0???

		n_pcs = pcs_d_l.get_n_pcs();

		score_la = 0;
		for ( j = 1; j <= n_pcs; j++ ) {
			idx = A_index[j]; //index on the vector X_all_, Y_all_, Z_all_
			/*
			x = X_all[idx];
			y = Y_all[idx];
			z = Z_all[idx];
			fill_A_line(A, xM_, yM_, zM_, x, y, z);
			PCS_calc = A[1]*Xxx + A[2]*Xxy + A[3]*Xxz + A[4]*Xyy + A[5]*Xyz;
			*/
			PCS_calc = A_all[idx][1]*Xxx + A_all[idx][2]*Xxy + A_all[idx][3]*Xxz + A_all[idx][4]*Xyy + A_all[idx][5]*Xyz;
			/*
			PCS_calc =
			Xxx_coef_vect_[idx]*Xxx +
			Xxy_coef_vect_[idx]*Xxy +
			Xxz_coef_vect_[idx]*Xxz +
			Xyy_coef_vect_[idx]*Xyy +
			Xyz_coef_vect_[idx]*Xyz;
			*/

			PCS_exp = cstyle_b[j];
			score_la += (PCS_calc-PCS_exp) * (PCS_calc-PCS_exp);
		}

		score += sqrt(score_la) * weight / pcs_d_l.get_normalization_factor();
	}

	// std::cerr << score << std::endl;
	return(score);
}

void
TensorsOptimizerFix::dfunc_numeric(core::optimization::Multivec const & vars,
	core::optimization::Multivec & dE_dvars
) const{

	core::Size i;
	core::Real delta(0.0001);

	for ( i = 1; i <= vars.size(); i++ ) {
		core::optimization::Multivec vars_bis_pd (vars);
		core::optimization::Multivec vars_bis_md (vars);
		vars_bis_pd[i] += delta;
		vars_bis_md[i] -= delta;
		core::Real value_pd (func(vars_bis_pd) / delta / 2);
		core::Real value_md (func(vars_bis_md) / delta / 2);
		dE_dvars[i] = value_pd - value_md;
	}
}

void
TensorsOptimizerFix::dfunc(core::optimization::Multivec const & vars,
	core::optimization::Multivec & dE_dvars
) const{

	dfunc_exact(vars, dE_dvars);
	/*
	core::Size i;
	std::cerr << "df_exact ";
	for(i = 1; i <= vars.size(); i++){
	std::cerr << std::setw(8) << dE_dvars[i] << " " ;
	}
	std::cerr << std::endl;

	dfunc_numeric(vars, dE_dvars);
	std::cerr << "df_numer ";
	for(i = 1; i <= vars.size(); i++){
	std::cerr  << std::setw(8)<< dE_dvars[i] << " " ;
	}
	std::cerr << std::endl;
	*/
}


core::Real
TensorsOptimizerFix::operator ()( core::optimization::Multivec const & vars ) const{
	return(func(vars));
}

void
TensorsOptimizerFix::dfunc_exact(core::optimization::Multivec const & vars,
	core::optimization::Multivec & dE_dvars
) const{
	core::Size i, j, k, idx;
	core::Size n_la, n_pcs;
	// core::Real x, y, z, xS, yS, zS;
	core::Real PCS_calc, PCS_exp;
	// core::Real common2, common3;
	// core::Real r5, r2, r3;

	utility::vector1<core::Real> A(5, 0);

	n_la = pcs_d_c_.get_n_lanthanides();

	if ( n_la != (dE_dvars.size()) / 5 ) {
		utility_exit_with_message("The number of lanthanides is inconsistent with the derivativessize of dE_dvars");
	}

	//lets' init everyone to zero because we are going to use +=
	for ( i = 1; i <= dE_dvars.size(); i++ ) {
		dE_dvars[i] = 0;
	}

	utility::vector1< core::Real > vec_temp;
	vec_temp.resize(dE_dvars.size(), 0);

	//those are the lanthanides current coordinates.

	const utility::vector1<PcsDataLanthanide> & pcs_d_l_vec(pcs_d_c_.get_pcs_data_per_lanthanides_all());
	utility::vector1< utility::vector1<core::Real> >  const & A_all( pcs_d_c_.get_A_all());
	/*
	const utility::vector1<core::Real> & X_all(pcs_d_c_.get_X_all());
	const utility::vector1<core::Real> & Y_all(pcs_d_c_.get_Y_all());
	const utility::vector1<core::Real> & Z_all(pcs_d_c_.get_Z_all());
	*/

	if ( n_la != pcs_d_l_vec.size() ) {
		utility_exit_with_message("The number of lanthanides is inconsistent");
	}


	for ( i = 1; i <= n_la; i++ ) {
		core::Size p;
		for ( p = 1; p <= vec_temp.size(); p++ ) {
			vec_temp[p] = 0;
		}

		core::Real Xxx(vars[5*(i-1) + 1]);
		core::Real Xxy(vars[5*(i-1) + 2]);
		core::Real Xxz(vars[5*(i-1) + 3]);
		core::Real Xyy(vars[5*(i-1) + 4]);
		core::Real Xyz(vars[5*(i-1) + 5]);


		PcsDataLanthanide const & pcs_d_l(pcs_d_l_vec[i]);
		utility::vector1<core::Size> const & A_index( pcs_d_l.get_A_index());
		utility::vector1< core::Real > const & cstyle_b(pcs_d_l.get_cstyle_b());

		core::Real weight(pcs_d_l.get_weight());
		n_pcs = pcs_d_l.get_n_pcs();

		core::Real square_sum = 0;

		for ( j = 1; j <= n_pcs; j++ ) {
			idx = A_index[j]; //index on the vector X_all_, Y_all_, Z_all_
			/*
			xS = X_all[idx];
			yS = Y_all[idx];
			zS = Z_all[idx];
			fill_A_line(A, xM_, yM_, zM_, xS, yS, zS);
			PCS_calc = A[1]*Xxx + A[2]*Xxy + A[3]*Xxz + A[4]*Xyy + A[5]*Xyz;
			*/
			PCS_calc = A_all[idx][1]*Xxx + A_all[idx][2]*Xxy + A_all[idx][3]*Xxz + A_all[idx][4]*Xyy + A_all[idx][5]*Xyz;
			/*
			PCS_calc =
			Xxx_coef_vect_[idx]*Xxx +
			Xxy_coef_vect_[idx]*Xxy +
			Xxz_coef_vect_[idx]*Xxz +
			Xyy_coef_vect_[idx]*Xyy +
			Xyz_coef_vect_[idx]*Xyz;
			*/

			PCS_exp = cstyle_b[j];

			core::Real deviation(PCS_calc - PCS_exp);

			square_sum += deviation * deviation;

			//This is for the derivatve of Xxx, Xxy, Xxz, Xyy, Xyz

			for ( k = 1; k <= 5; k++ ) {
				// vec_temp[5*(i-1) + k] += deviation * A[k];
				vec_temp[5*(i-1) + k] += deviation * A_all[idx][k];
			}

			/*
			vec_temp[5*(i-1) + 1] += deviation * Xxx_coef_vect_[idx];
			vec_temp[5*(i-1) + 2] += deviation * Xxy_coef_vect_[idx];
			vec_temp[5*(i-1) + 3] += deviation * Xxz_coef_vect_[idx];
			vec_temp[5*(i-1) + 4] += deviation * Xyy_coef_vect_[idx];
			vec_temp[5*(i-1) + 5] += deviation * Xyz_coef_vect_[idx];
			*/
		}

		core::Real one_over_RMSD(1.0/sqrt(square_sum));
		core::Real common_temp( one_over_RMSD * weight / pcs_d_l.get_normalization_factor() );
		for ( k = 1; k <= 5; k++ ) {
			dE_dvars[5*(i-1) + k] += vec_temp[5*(i-1) + k] * common_temp;
		}
	}

}

}//namespace pcs2
}//namespace methods
}//namespace scoring
}//namespace protocols
