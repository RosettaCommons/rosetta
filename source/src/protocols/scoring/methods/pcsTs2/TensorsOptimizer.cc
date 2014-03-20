// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

 //////////////////////////////////////////////
 /// @begin
 ///
 /// @file protocols/scoring/TensorsOptimizer.cc
 ///
 /// @brief
 ///
 /// @detailed
 ///
 /// @param
 ///
 /// @return
 ///
 /// @remarks
 ///
 /// @references C Schmitz et.al. J Mol Biol. Mar 9, 2012; 416(5): 668–677 ; Yagi H et.al Structure, 2013, 21(6):883-890
 ///
 /// @authorsv Christophe Schmitz , Kala Bharath Pilla
 ///
 /// @last_modified Mar 2014
 ////////////////////////////////////////////////


// Unit headers
#include <protocols/scoring/methods/pcsTs2/TensorsOptimizer.hh>
#include <protocols/scoring/methods/pcsTs2/PseudocontactShiftData.hh>
#include <protocols/scoring/methods/pcsTs2/PseudocontactShiftInput.hh> // REQUIRED FOR WINDOWS
// Package headers

// Project headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/exit.hh>

// Numeric headers
#include <numeric/constants.hh>

// Objexx headers

// C++ headers
#include <sstream>
#include <iostream>

namespace protocols{
namespace scoring{
namespace methods{
namespace pcsTs2{

basic::Tracer TR_tsr_opt_Ts2("protocols.scoring.methods.pcsTs2.TensorsOptimizer_Ts2");


TensorsOptimizer_Ts2::TensorsOptimizer_Ts2(PCS_data_Ts2 const & pcs_d):
	pcs_d_(pcs_d)
{
}

TensorsOptimizer_Ts2::~TensorsOptimizer_Ts2(){
}

core::Real
TensorsOptimizer_Ts2::operator ()( core::optimization::Multivec const & vars ) const{
	core::Size i, j, idx;
	core::Size n_la, n_pcs;
	core::Real x, y, z;
	core::Real PCS_calc, PCS_exp, score, score_la;;
	utility::vector1<core::Real> A(5, 0);

	score = 0;
	n_la = pcs_d_.get_n_lanthanides();

	if(n_la != (vars.size() - 3) / 5){
		TR_tsr_opt_Ts2 << "n_la: " << n_la << " (vars.size() - 3) / 5):" << (vars.size() - 3) / 5 << " vars.size(): " << vars.size() <<std::endl;
		utility_exit_with_message("The number of lanthanides is inconsistent with the derivatives size of vars");
	}

	//those are the lanthanides current coordinates.
	core::Real xM(vars[1]);
	core::Real yM(vars[2]);
	core::Real zM(vars[3]);

	const utility::vector1<PCS_data_per_lanthanides_Ts2> & pcs_d_p_l_a(pcs_d_.get_pcs_data_per_lanthanides_all());

  const utility::vector1<core::Real> & X_all(pcs_d_.get_X_all());
  const utility::vector1<core::Real> & Y_all(pcs_d_.get_Y_all());
  const utility::vector1<core::Real> & Z_all(pcs_d_.get_Z_all());

	if(n_la != pcs_d_p_l_a.size()){
		utility_exit_with_message("The number of lanthanides is inconsistent");
	}

	for (i = 1; i <= n_la; i++){
		core::Real Xxx(vars[3 + 5*(i-1) + 1]);
		core::Real Xxy(vars[3 + 5*(i-1) + 2]);
		core::Real Xxz(vars[3 + 5*(i-1) + 3]);
		core::Real Xyy(vars[3 + 5*(i-1) + 4]);
		core::Real Xyz(vars[3 + 5*(i-1) + 5]);

		PCS_data_per_lanthanides_Ts2 const & pcs_d_p_l(pcs_d_p_l_a[i]);
		utility::vector1<core::Size> const & A_index( pcs_d_p_l.get_A_index());
		ObjexxFCL::FArray1D< core::Real > const & fstyle_b(pcs_d_p_l.get_fstyle_b());

		core::Real weight(pcs_d_p_l.get_weight());
		n_pcs = pcs_d_p_l.get_n_pcs();

		//		std::cerr << weight << " applyed to minimization" << std::endl;
		score_la = 0;
		for (j = 1; j <= n_pcs; j++){
			idx = A_index[j]; //index on the vector X_all_, Y_all_, Z_all_
			x = X_all[idx];
			y = Y_all[idx];
			z = Z_all[idx];
			fill_A_line(A, xM, yM, zM, x, y, z);
			PCS_calc = A[1]*Xxx + A[2]*Xxy + A[3]*Xxz + A[4]*Xyy + A[5]*Xyz;
			PCS_exp = fstyle_b(j);
			score_la += (PCS_calc-PCS_exp) * (PCS_calc-PCS_exp);
		}
		score += score_la * weight / pcs_d_p_l.get_normalization_factor() / pcs_d_p_l.get_normalization_factor();
	}
	/*
	std::cerr << "current: ";
	for (i = 1; i<= vars.size(); i++){
		std::cerr << vars[i] << " ";
	}
	std::cerr << std::endl;
	*/

  return(sqrt(score));
}

void
TensorsOptimizer_Ts2::dfunc(core::optimization::Multivec const & vars,
												core::optimization::Multivec & dE_dvars
												) const{
	core::Size i, j, k, idx;
	core::Size n_la, n_pcs;
	core::Real x, y, z, xS, yS, zS;
	core::Real PCS_calc, PCS_exp, common, common2, common3;
	core::Real r5, r2, r3;

	utility::vector1<core::Real> A(5, 0);

	n_la = pcs_d_.get_n_lanthanides();

	if(n_la != (dE_dvars.size() - 3) / 5){
		utility_exit_with_message("The number of lanthanides is inconsistent with the derivativessize of dE_dvars");
	}

	//lets' init everyone to zero because we are going to use +=
	for (i = 1; i <= dE_dvars.size(); i++){
		dE_dvars[i] = 0;
	}

	//those are the lanthanides current coordinates.
	core::Real xM(vars[1]);
	core::Real yM(vars[2]);
	core::Real zM(vars[3]);

	const utility::vector1<PCS_data_per_lanthanides_Ts2> & pcs_d_p_l_a(pcs_d_.get_pcs_data_per_lanthanides_all());

  const utility::vector1<core::Real> & X_all(pcs_d_.get_X_all());
  const utility::vector1<core::Real> & Y_all(pcs_d_.get_Y_all());
  const utility::vector1<core::Real> & Z_all(pcs_d_.get_Z_all());

	if(n_la != pcs_d_p_l_a.size()){
		utility_exit_with_message("The number of lanthanides is inconsistent");
	}

	for (i = 1; i <= n_la; i++){
		core::Real Xxx(vars[3 + 5*(i-1) + 1]);
		core::Real Xxy(vars[3 + 5*(i-1) + 2]);
		core::Real Xxz(vars[3 + 5*(i-1) + 3]);
		core::Real Xyy(vars[3 + 5*(i-1) + 4]);
		core::Real Xyz(vars[3 + 5*(i-1) + 5]);

		PCS_data_per_lanthanides_Ts2 const & pcs_d_p_l(pcs_d_p_l_a[i]);
		utility::vector1<core::Size> const & A_index( pcs_d_p_l.get_A_index());
		ObjexxFCL::FArray1D< core::Real > const & fstyle_b(pcs_d_p_l.get_fstyle_b());

		core::Real weight(pcs_d_p_l.get_weight());
		n_pcs = pcs_d_p_l.get_n_pcs();
		for (j = 1; j <= n_pcs; j++){
			idx = A_index[j]; //index on the vector X_all_, Y_all_, Z_all_
			xS = X_all[idx];
			yS = Y_all[idx];
			zS = Z_all[idx];
			fill_A_line(A, xM, yM, zM, xS, yS, zS);
			PCS_calc = A[1]*Xxx + A[2]*Xxy + A[3]*Xxz + A[4]*Xyy + A[5]*Xyz;
			PCS_exp = fstyle_b(j);
			common = 2 * (PCS_calc - PCS_exp) * weight;

			//This is for the derivatve of Xxx, Xxy, Xxz, Xyy, Xyz
			for (k = 1; k <= 5; k++){
				dE_dvars[3 + 5*(i-1) + k] += common * A[k];
			}

			//This is for the derivatve of xM, yM, zM
			//could be a little bit optimized?
			x=xS-xM;
			y=yS-yM;
			z=zS-zM;
			r2=x*x+y*y+z*z;
			r3=sqrt(r2)*r2;
			r5=r3*r2;
			common2 = 5*r3*4*PCS_calc * numeric::constants::d::pi / 10000.0  ;
			common3 = common/(4*r5 * numeric::constants::d::pi)  *10000.0 ;
			dE_dvars[1] += (x*common2 - 2*(x*Xxx + y*Xxy + z*Xxz)) * common3;
			dE_dvars[2] += (y*common2 - 2*(x*Xxy + y*Xyy + z*Xyz)) * common3;
			dE_dvars[3] += (z*common2 - 2*(-z*Xxx + x*Xxz -z*Xyy + y*Xyz)) * common3;

		}
	}

}


}//namespace pcsTs2
}//namespace methods
}//namespace scoring
}//namespace protocols
