// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

//////////////////////////////////////////////
///
/// @file protocols/scoring/methods/pcs2/TensorsOptimizerSvd.cc
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
#include <protocols/scoring/methods/pcs2/TensorsOptimizerSvd.hh>
#include <protocols/scoring/methods/pcs2/PcsDataCenter.hh>
#include <protocols/scoring/methods/pcs2/PcsTensor.hh>

// Package headers

// Project headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/exit.hh>

// Numeric headers
#include <numeric/constants.hh>

// Objexx headers

// C++ headers

#include <protocols/scoring/methods/pcs2/TensorsOptimizer.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace scoring {
namespace methods {
namespace pcs2 {

using namespace core::optimization;

static basic::Tracer TR_TensorsOptimizerSvd( "protocols.scoring.methods.pcs.TensorsOptimizerSvd" );

/*
TensorsOptimizerSvd::TensorsOptimizerSvd():
pcs_d_c_(PcsDataCenter())
{

utility_exit_with_message("You shouldn't call the empty constructor for class TensorsOptimizerSvd");
}
*/

TensorsOptimizerSvd::TensorsOptimizerSvd(PcsDataCenter /*const*/ & pcs_d_c):
	pcs_d_c_(pcs_d_c)
{
}

TensorsOptimizerSvd::~TensorsOptimizerSvd()= default;


bool
TensorsOptimizerSvd::abort_min(core::optimization::Multivec const & vars ) const{

	// return false;

	core::Real limit(4000);
	core::Size i;
	core::Size n_la;

	n_la = pcs_d_c_.get_n_lanthanides();

	core::Real xM(vars[1]);
	core::Real yM(vars[1]);
	core::Real zM(vars[1]);
	if ( (xM*xM + yM*yM + zM*zM) > 90000 ) {
		return true;
	}

	for ( i = 1; i <= n_la; i++ ) {
		core::Real Xxx(vars[3 + 5*(i-1) + 1]);
		core::Real Xxy(vars[3 + 5*(i-1) + 2]);
		core::Real Xxz(vars[3 + 5*(i-1) + 3]);
		core::Real Xyy(vars[3 + 5*(i-1) + 4]);
		core::Real Xyz(vars[3 + 5*(i-1) + 5]);
		if ( (fabs(Xxx) + fabs(Xxy) + fabs(Xxz) + fabs(Xyy) + fabs(Xyz)) > limit ) {
			return true;
		}
	}
	return false;
}

core::Real
TensorsOptimizerSvd::func( core::optimization::Multivec const & vars ) const{
	core::Size i; //, j, idx;
	core::Size n_la;//, n_pcs;
	core::Real score; //PCS_calc, PCS_exp, score_la;;
	utility::vector1<core::Real> A(5, 0);


	score = 0;
	n_la = pcs_d_c_.get_n_lanthanides();

	//those are the lanthanides current coordinates.
	core::Real xM(vars[1]);
	core::Real yM(vars[2]);
	core::Real zM(vars[3]);

	PcsTensor pcs_t(0, 0, 0, 0, 0, "");

	pcs_d_c_.update_matrix_A_all_for_svd(xM, yM, zM);
	// std::cerr << "updating from func" << std::endl;
	// std::cerr  << std::setprecision(20) << xM << " " << yM << " " << zM  << std::endl;
	utility::vector1<PcsDataLanthanide> & pcs_d_l_vec(pcs_d_c_.get_pcs_data_per_lanthanides_all());

	if ( n_la != pcs_d_l_vec.size() ) {
		utility_exit_with_message("The number of lanthanides is inconsistent");
	}


	// utility::vector1< utility::vector1<core::Real> >  const & A_all(pcs_d_c_.get_A_all());


	for ( i = 1; i <= n_la; i++ ) {
		PcsDataLanthanide & pcs_d_l(pcs_d_l_vec[i]);

		//  core::Real weight(pcs_d_l.get_weight());
		//  score += pcs_d_l.calculate_cost_only_with_svd();

		score += pcs_d_l.calculate_tensor_and_cost_with_svd(pcs_t);
		//score += pcs_d_l.calculate_cost_only_with_svd();// * weight * pcs_d_l.get_normalization_factor_inversed();
		//weight and normalization_factor allready applied in cost_with_svd function!

	}

	return(score);
}

void
TensorsOptimizerSvd::dfunc_numeric(core::optimization::Multivec const & vars,
	core::optimization::Multivec & dE_dvars
) const{

	core::Size i;
	core::Real delta(0.0001);

	for ( i = 1; i <= vars.size(); i++ ) {
		core::optimization::Multivec vars_bis_pd (vars);
		core::optimization::Multivec vars_bis_md (vars);
		vars_bis_pd[i] += delta;
		vars_bis_md[i] -= delta;
		core::Real value_pd (func(vars_bis_pd) / delta / 2.0);
		core::Real value_md (func(vars_bis_md) / delta / 2.0);
		dE_dvars[i] = value_pd - value_md;
	}
}

void
TensorsOptimizerSvd::dfunc(core::optimization::Multivec const & vars,
	core::optimization::Multivec & dE_dvars
) const{

	dfunc_exact(vars, dE_dvars);
	//dfunc_numeric(vars, dE_dvars);
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
TensorsOptimizerSvd::operator ()( core::optimization::Multivec const & vars ) const{
	return(func(vars));
}

void
TensorsOptimizerSvd::dfunc_exact(core::optimization::Multivec const & vars,
	core::optimization::Multivec & dE_dvars
) const{
	core::Size i, j, idx;
	core::Size n_la, n_pcs;
	core::Real x, y, z, xS, yS, zS;
	core::Real PCS_calc, PCS_exp, common2, common3;
	core::Real r5, r2, r3;

	utility::vector1<core::Real> A(5, 0);

	n_la = pcs_d_c_.get_n_lanthanides();

	//lets' init everyone to zero because we are going to use +=
	for ( i = 1; i <= dE_dvars.size(); i++ ) {
		dE_dvars[i] = 0;
	}

	utility::vector1< core::Real > vec_temp;
	vec_temp.resize(dE_dvars.size(), 0);

	//those are the lanthanides current coordinates.
	core::Real xM(vars[1]);
	core::Real yM(vars[2]);
	core::Real zM(vars[3]);

	utility::vector1<PcsDataLanthanide> & pcs_d_l_vec(pcs_d_c_.get_pcs_data_per_lanthanides_all());

	utility::vector1< utility::vector1<core::Real> >  const & A_all(pcs_d_c_.get_A_all());

	const utility::vector1<core::Real> & X_all(pcs_d_c_.get_X_all());
	const utility::vector1<core::Real> & Y_all(pcs_d_c_.get_Y_all());
	const utility::vector1<core::Real> & Z_all(pcs_d_c_.get_Z_all());


	if ( n_la != pcs_d_l_vec.size() ) {
		utility_exit_with_message("The number of lanthanides is inconsistent");
	}


	PcsTensor pcs_t(0, 0, 0, 0, 0, "");
	for ( i = 1; i <= n_la; i++ ) {
		PcsDataLanthanide & pcs_d_l(pcs_d_l_vec[i]);
		utility::vector1<core::Size> const & A_index( pcs_d_l.get_A_index());
		utility::vector1< core::Real > const & cstyle_b(pcs_d_l.get_cstyle_b());
		core::Real weight(pcs_d_l.get_weight());

		pcs_d_l.retrieve_tensor_from_svd(pcs_t);

		core::Real Xxx(pcs_t.get_chi_xx());
		core::Real Xxy(pcs_t.get_chi_xy());
		core::Real Xxz(pcs_t.get_chi_xz());
		core::Real Xyy(pcs_t.get_chi_yy());
		core::Real Xyz(pcs_t.get_chi_yz());

		core::Size p;
		for ( p = 1; p <= vec_temp.size(); p++ ) {
			vec_temp[p] = 0;
		}


		n_pcs = pcs_d_l.get_n_pcs();

		core::Real square_sum = 0;

		for ( j = 1; j <= n_pcs; j++ ) {
			idx = A_index[j]; //index on the vector X_all_, Y_all_, Z_all_
			xS = X_all[idx];
			yS = Y_all[idx];
			zS = Z_all[idx];
			PCS_calc = A_all[idx][1]*Xxx + A_all[idx][2]*Xxy + A_all[idx][3]*Xxz + A_all[idx][4]*Xyy + A_all[idx][5]*Xyz;
			PCS_exp = cstyle_b[j];

			core::Real deviation(PCS_calc - PCS_exp);

			square_sum += deviation * deviation;

			//This is for the derivatve of xM, yM, zM; could be a little bit optimized?
			x=xS-xM;
			y=yS-yM;
			z=zS-zM;
			r2=x*x+y*y+z*z;
			r3=sqrt(r2)*r2;
			r5=r3*r2;

			common2 = FACT_20_PI_OVER_10000*PCS_calc*r3;
			common3 = deviation/r5;

			vec_temp[1] += (x*common2 - 2*( x*Xxx + y*Xxy + z*Xxz)) * common3;
			vec_temp[2] += (y*common2 - 2*( x*Xxy + y*Xyy + z*Xyz)) * common3;
			vec_temp[3] += (z*common2 - 2*(-z*Xxx + x*Xxz - z*Xyy + y*Xyz)) * common3;
		}
		//  std::cerr << "square_sum" <<square_sum << std::endl;
		core::Real one_over_RMSD(1.0/sqrt(square_sum));
		core::Real common_temp( one_over_RMSD * weight * pcs_d_l.get_normalization_factor_inversed() );

		dE_dvars[1] += common_temp * vec_temp[1];
		dE_dvars[2] += common_temp * vec_temp[2];
		dE_dvars[3] += common_temp * vec_temp[3];
	}

	dE_dvars[1] *= FACT_10000_OVER_4PI;
	dE_dvars[2] *= FACT_10000_OVER_4PI;
	dE_dvars[3] *= FACT_10000_OVER_4PI;

}

}//namespace pcs2
}//namespace methods
}//namespace scoring
}//namespace protocols
