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
/// @file protocols/scoring/TensorsOptimizer.cc
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
#include <protocols/scoring/methods/pcs/TensorsOptimizer.hh>
#include <protocols/scoring/methods/pcs/PseudocontactShiftData.hh>
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

#include <utility/vector1.hh>


namespace protocols {
namespace scoring {
namespace methods {
namespace pcs {

static thread_local basic::Tracer TR_tsr_opt( "protocols.scoring.methods.pcs.TensorsOptimizer" );

TensorsOptimizer::TensorsOptimizer(PCS_data const & pcs_d):
	pcs_d_(pcs_d)
{
}

TensorsOptimizer::~TensorsOptimizer(){
}

core::Real
TensorsOptimizer::operator ()( core::optimization::Multivec const & vars ) const{
	core::Size i, j, idx;
	core::Size n_la, n_pcs;
	core::Real x, y, z;
	core::Real PCS_calc, PCS_exp, score, score_la;;
	utility::vector1<core::Real> A(5, 0);

	score = 0;
	n_la = pcs_d_.get_n_lanthanides();

	if ( n_la != (vars.size() - 3) / 5 ) {
		TR_tsr_opt << "n_la: " << n_la << " (vars.size() - 3) / 5):" << (vars.size() - 3) / 5 << " vars.size(): " << vars.size() <<std::endl;
		utility_exit_with_message("The number of lanthanides is inconsistent with the derivatives size of vars");
	}

	//those are the lanthanides current coordinates.
	core::Real xM(vars[1]);
	core::Real yM(vars[2]);
	core::Real zM(vars[3]);

	const utility::vector1<PCS_data_per_lanthanides> & pcs_d_p_l_a(pcs_d_.get_pcs_data_per_lanthanides_all());

	const utility::vector1<core::Real> & X_all(pcs_d_.get_X_all());
	const utility::vector1<core::Real> & Y_all(pcs_d_.get_Y_all());
	const utility::vector1<core::Real> & Z_all(pcs_d_.get_Z_all());

	if ( n_la != pcs_d_p_l_a.size() ) {
		utility_exit_with_message("The number of lanthanides is inconsistent");
	}

	for ( i = 1; i <= n_la; i++ ) {
		core::Real Xxx(vars[3 + 5*(i-1) + 1]);
		core::Real Xxy(vars[3 + 5*(i-1) + 2]);
		core::Real Xxz(vars[3 + 5*(i-1) + 3]);
		core::Real Xyy(vars[3 + 5*(i-1) + 4]);
		core::Real Xyz(vars[3 + 5*(i-1) + 5]);

		PCS_data_per_lanthanides const & pcs_d_p_l(pcs_d_p_l_a[i]);
		utility::vector1<core::Size> const & A_index( pcs_d_p_l.get_A_index());
		ObjexxFCL::FArray1D< core::Real > const & fstyle_b(pcs_d_p_l.get_fstyle_b());

		core::Real weight(pcs_d_p_l.get_weight());
		n_pcs = pcs_d_p_l.get_n_pcs();

		//  std::cerr << weight << " applyed to minimization" << std::endl;
		score_la = 0;
		for ( j = 1; j <= n_pcs; j++ ) {
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
TensorsOptimizer::dfunc(core::optimization::Multivec const & vars,
	core::optimization::Multivec & dE_dvars
) const{
	core::Size i, j, k, idx;
	core::Size n_la, n_pcs;
	core::Real x, y, z, xS, yS, zS;
	core::Real PCS_calc, PCS_exp, common, common2, common3;
	core::Real r5, r2, r3;

	utility::vector1<core::Real> A(5, 0);

	n_la = pcs_d_.get_n_lanthanides();

	if ( n_la != (dE_dvars.size() - 3) / 5 ) {
		utility_exit_with_message("The number of lanthanides is inconsistent with the derivativessize of dE_dvars");
	}

	//lets' init everyone to zero because we are going to use +=
	for ( i = 1; i <= dE_dvars.size(); i++ ) {
		dE_dvars[i] = 0;
	}

	//those are the lanthanides current coordinates.
	core::Real xM(vars[1]);
	core::Real yM(vars[2]);
	core::Real zM(vars[3]);

	const utility::vector1<PCS_data_per_lanthanides> & pcs_d_p_l_a(pcs_d_.get_pcs_data_per_lanthanides_all());

	const utility::vector1<core::Real> & X_all(pcs_d_.get_X_all());
	const utility::vector1<core::Real> & Y_all(pcs_d_.get_Y_all());
	const utility::vector1<core::Real> & Z_all(pcs_d_.get_Z_all());

	if ( n_la != pcs_d_p_l_a.size() ) {
		utility_exit_with_message("The number of lanthanides is inconsistent");
	}

	for ( i = 1; i <= n_la; i++ ) {
		core::Real Xxx(vars[3 + 5*(i-1) + 1]);
		core::Real Xxy(vars[3 + 5*(i-1) + 2]);
		core::Real Xxz(vars[3 + 5*(i-1) + 3]);
		core::Real Xyy(vars[3 + 5*(i-1) + 4]);
		core::Real Xyz(vars[3 + 5*(i-1) + 5]);

		PCS_data_per_lanthanides const & pcs_d_p_l(pcs_d_p_l_a[i]);
		utility::vector1<core::Size> const & A_index( pcs_d_p_l.get_A_index());
		ObjexxFCL::FArray1D< core::Real > const & fstyle_b(pcs_d_p_l.get_fstyle_b());

		core::Real weight(pcs_d_p_l.get_weight());
		n_pcs = pcs_d_p_l.get_n_pcs();
		for ( j = 1; j <= n_pcs; j++ ) {
			idx = A_index[j]; //index on the vector X_all_, Y_all_, Z_all_
			xS = X_all[idx];
			yS = Y_all[idx];
			zS = Z_all[idx];
			fill_A_line(A, xM, yM, zM, xS, yS, zS);
			PCS_calc = A[1]*Xxx + A[2]*Xxy + A[3]*Xxz + A[4]*Xyy + A[5]*Xyz;
			PCS_exp = fstyle_b(j);
			common = 2 * (PCS_calc - PCS_exp) * weight;

			//This is for the derivatve of Xxx, Xxy, Xxz, Xyy, Xyz
			for ( k = 1; k <= 5; k++ ) {
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


	/*
	for (i = 1; i<= dE_dvars.size(); i++){
	std::cerr << dE_dvars[i] << " ";
	}

	std::cerr <<" should be equalt to:" << std::endl;
	dfunc_test(vars);
	std::cerr <<" current: ";

	for (i = 1; i<= dE_dvars.size(); i++){
	std::cerr << vars[i] << " ";
	}
	std::cerr << std::endl;
	*/
}


/*
//This function is meant to calculate numeric derivatives to test.
void
TensorsOptimizer::dfunc_test(optimization::Multivec const & vars) const{

core::Size i, j, k, idx;
core::Size n_la, n_pcs;
core::Real x, y, z, xS, yS, zS;
core::Real PCS_calc, PCS_exp, common, common2, common3, PCS_calc_p_d, PCS_calc_m_d;
core::Real r5, r2, r3;
core::Real delta(0.00001);

utility::vector1<core::Real> A(5, 0);

utility::vector1<core::Real> A_p_d(5, 0);
utility::vector1<core::Real> A_m_d(5, 0);

n_la = pcs_d_.get_n_lanthanides();

optimization::Multivec dE_dvars(3+5*n_la);

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

const utility::vector1<PCS_data_per_lanthanides> & pcs_d_p_l_a(pcs_d_.get_pcs_data_per_lanthanides_all());

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
const PCS_data_per_lanthanides & pcs_d_p_l(pcs_d_p_l_a[i]);
const utility::vector1<core::Size> & A_index( pcs_d_p_l.get_A_index());
const FArray1D< core::Real > & fstyle_b(pcs_d_p_l.get_fstyle_b());
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
common = 2.0 * (PCS_calc - PCS_exp) * weight;


PCS_calc_p_d = A[1]*(Xxx+delta) + A[2]*Xxy + A[3]*Xxz + A[4]*Xyy + A[5]*Xyz;
PCS_calc_m_d = A[1]*(Xxx-delta) + A[2]*Xxy + A[3]*Xxz + A[4]*Xyy + A[5]*Xyz;
dE_dvars[3 + 5*(i-1) + 1] += common * (PCS_calc_p_d - PCS_calc_m_d) / (2.0*delta);

PCS_calc_p_d = A[1]*Xxx + A[2]*(Xxy+delta) + A[3]*Xxz + A[4]*Xyy + A[5]*Xyz;
PCS_calc_m_d = A[1]*Xxx + A[2]*(Xxy-delta) + A[3]*Xxz + A[4]*Xyy + A[5]*Xyz;
dE_dvars[3 + 5*(i-1) + 2] += common * (PCS_calc_p_d - PCS_calc_m_d) / (2.0*delta);

PCS_calc_p_d = A[1]*Xxx + A[2]*Xxy + A[3]*(Xxz+delta) + A[4]*Xyy + A[5]*Xyz;
PCS_calc_m_d = A[1]*Xxx + A[2]*Xxy + A[3]*(Xxz-delta) + A[4]*Xyy + A[5]*Xyz;
dE_dvars[3 + 5*(i-1) + 3] += common * (PCS_calc_p_d - PCS_calc_m_d) / (2.0*delta);

PCS_calc_p_d = A[1]*Xxx + A[2]*Xxy + A[3]*Xxz + A[4]*(Xyy+delta) + A[5]*Xyz;
PCS_calc_m_d = A[1]*Xxx + A[2]*Xxy + A[3]*Xxz + A[4]*(Xyy-delta) + A[5]*Xyz;
dE_dvars[3 + 5*(i-1) + 4] += common * (PCS_calc_p_d - PCS_calc_m_d) / (2.0*delta);

PCS_calc_p_d = A[1]*Xxx + A[2]*Xxy + A[3]*Xxz + A[4]*Xyy + A[5]*(Xyz+delta);
PCS_calc_m_d = A[1]*Xxx + A[2]*Xxy + A[3]*Xxz + A[4]*Xyy + A[5]*(Xyz-delta);
dE_dvars[3 + 5*(i-1) + 5] += common * (PCS_calc_p_d - PCS_calc_m_d) / (2.0*delta);


fill_A_line(A_p_d, xM + delta, yM, zM, xS, yS, zS);
PCS_calc_p_d = A_p_d[1]*Xxx + A_p_d[2]*Xxy + A_p_d[3]*Xxz + A_p_d[4]*Xyy + A_p_d[5]*Xyz;
fill_A_line(A_m_d, xM - delta, yM, zM, xS, yS, zS);
PCS_calc_m_d = A_m_d[1]*Xxx + A_m_d[2]*Xxy + A_m_d[3]*Xxz + A_m_d[4]*Xyy + A_m_d[5]*Xyz;
dE_dvars[1] += common * (PCS_calc_p_d - PCS_calc_m_d) / (2.0*delta);

fill_A_line(A_p_d, xM, yM + delta, zM, xS, yS, zS);
PCS_calc_p_d = A_p_d[1]*Xxx + A_p_d[2]*Xxy + A_p_d[3]*Xxz + A_p_d[4]*Xyy + A_p_d[5]*Xyz;
fill_A_line(A_m_d, xM, yM - delta, zM, xS, yS, zS);
PCS_calc_m_d = A_m_d[1]*Xxx + A_m_d[2]*Xxy + A_m_d[3]*Xxz + A_m_d[4]*Xyy + A_m_d[5]*Xyz;
dE_dvars[2] += common * (PCS_calc_p_d - PCS_calc_m_d) / (2.0*delta);


fill_A_line(A_p_d, xM, yM, zM + delta, xS, yS, zS);
PCS_calc_p_d = A_p_d[1]*Xxx + A_p_d[2]*Xxy + A_p_d[3]*Xxz + A_p_d[4]*Xyy + A_p_d[5]*Xyz;
fill_A_line(A_m_d, xM, yM, zM - delta, xS, yS, zS);
PCS_calc_m_d = A_m_d[1]*Xxx + A_m_d[2]*Xxy + A_m_d[3]*Xxz + A_m_d[4]*Xyy + A_m_d[5]*Xyz;
dE_dvars[3] += common * (PCS_calc_p_d - PCS_calc_m_d) / (2.0*delta);

}
}


for (i = 1; i<= dE_dvars.size(); i++){
std::cerr << dE_dvars[i] << " ";
}
std::cerr << std::endl;

}
*/


}//namespace pcs
}//namespace methods
}//namespace scoring
}//namespace protocols
