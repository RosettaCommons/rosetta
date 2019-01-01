// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/nmr/pcs/PCSTensorOptimizer.cc
/// @brief   Implementation of class PCSTensorOptimizer
/// @details last Modified: 07/05/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

// Unit headers
#include <protocols/nmr/pcs/PCSTensorOptimizer.hh>

// Package headers
#include <core/scoring/nmr/pcs/PCSSingle.hh>
#include <core/scoring/nmr/pcs/PCSSingleSet.hh>
#include <core/scoring/nmr/util.hh>

// Project headers
#include <core/optimization/Multifunc.hh>
#include <core/optimization/types.hh>
#include <core/types.hh>

// Basic headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/exit.hh>

// Numeric headers
#include <numeric/constants.hh>
#include <numeric/xyzVector.hh>

// Objexx headers
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray1D.hh>

// C++ headers
#include <cmath>

namespace protocols {
namespace nmr {
namespace pcs {

static basic::Tracer TR( "protocols.nmr.pcs.PCSTensorOptimizer" );

/// @brief constructor with a vector of PCSSingleSet pointers as argument
PCSTensorOptimizer::PCSTensorOptimizer( utility::vector1< PCSSingleSetOP > const & singleset_vec ) :
	singleset_vec_(singleset_vec)
{ }

/// @brief destructor
PCSTensorOptimizer::~PCSTensorOptimizer() {}

/// @brief error function used in optimization of the PCS tensor parameter
core::Real
PCSTensorOptimizer::operator()( Multivec const & tensor_params ) const {
	Size number_singlesets_to_optimize = singleset_vec_.size();

	// We are optimizing the whole vector of PCSSingleSets at once because the metals should share similar coordinates.
	// Thus we assign them all the same xyz-coordinates and expect that the number of input parameters is N*5 + 3 whereas N = number_singlesets_to_optimize
	if ( number_singlesets_to_optimize != (tensor_params.size()-3) / 5 ) {
		utility_exit_with_message("ERROR in PCSTensorOptimizer operator(). Number of tensor parameters passed is inconsistent with the number of PCSSingleSets to optimize. You need to provide N*5 + 3 parameters (N = number PCSSingleSets).");
	}

	Vector metal_coords(tensor_params[1], tensor_params[2], tensor_params[3]);
	Real total_score(0);

	for ( Size i = 1; i <= number_singlesets_to_optimize; ++i ) {
		Real chiT_xx(tensor_params[3 + 5*(i-1) + 1]);
		Real chiT_xy(tensor_params[3 + 5*(i-1) + 2]);
		Real chiT_xz(tensor_params[3 + 5*(i-1) + 3]);
		Real chiT_yy(tensor_params[3 + 5*(i-1) + 4]);
		Real chiT_yz(tensor_params[3 + 5*(i-1) + 5]);

		singleset_vec_[i]->update_matrix_A(metal_coords);
		ObjexxFCL::FArray2D<Real> const & matrix_A =  singleset_vec_[i]->get_matrix_A();
		ObjexxFCL::FArray1D<Real> const & pcs_values = singleset_vec_[i]->get_pcs_values();
		Size n_pcs(singleset_vec_[i]->get_number_pcs());
		ObjexxFCL::FArray1D<Real> const & single_pcs_weights = singleset_vec_[i]->get_pcs_single_weights();
		Real singleset_score(0);

		for ( Size j = 1; j <= n_pcs; ++j ) {
			Real pcs_calc = matrix_A(j,1)*chiT_xx + matrix_A(j,2)*chiT_xy + matrix_A(j,3)*chiT_xz + matrix_A(j,4)*chiT_yy + matrix_A(j,5)*chiT_yz;
			Real pcs_exp  = pcs_values(j);
			singleset_score += (pcs_calc - pcs_exp) * (pcs_calc - pcs_exp) * single_pcs_weights(j);
		}
		//total_score += std::sqrt(singleset_score) * singleset_vec_[i]->get_weight();
		total_score += singleset_score * singleset_vec_[i]->get_weight();
	}
	//return std::sqrt(total_score);
	return total_score;
}

/// @brief gradient function used in optimization of the PCS tensor parameter
void
PCSTensorOptimizer::dfunc(
	Multivec const & tensor_params,
	Multivec & dPCS_dparams
) const
{
	Size number_singlesets_to_optimize = singleset_vec_.size();

	// We are optimizing the whole vector of PCSSingleSets at once because the metals should share similar coordinates.
	// Thus we assign them all the same xyz-coordinates and expect that the number of input parameters is N*5 + 3 whereas N = number_singlesets_to_optimize
	if ( number_singlesets_to_optimize != (tensor_params.size()-3) / 5 ) {
		utility_exit_with_message("ERROR in PCSTensorOptimizer dfunc(). Number of tensor parameters passed is inconsistent with the number of PCSSingleSets to optimize. You need to provide N*5 + 3 parameters (N = number PCSSingleSets).");
	}

	// Set deviations to 0 at the beginning because we will sum up their individual contributions
	for ( Size n = 1; n <= dPCS_dparams.size(); ++n ) {
		dPCS_dparams[n] = 0;
	}

	Vector metal_coords(tensor_params[1], tensor_params[2], tensor_params[3]);

	// loop over different lanthanides
	for ( Size i = 1; i <= number_singlesets_to_optimize; ++i ) {
		Real chiT_xx(tensor_params[3 + 5*(i-1) + 1]);
		Real chiT_xy(tensor_params[3 + 5*(i-1) + 2]);
		Real chiT_xz(tensor_params[3 + 5*(i-1) + 3]);
		Real chiT_yy(tensor_params[3 + 5*(i-1) + 4]);
		Real chiT_yz(tensor_params[3 + 5*(i-1) + 5]);

		// Get some data relevant for PCS calculation for one lanthanide dataset (i.e. PCSSingleSet)
		singleset_vec_[i]->update_matrix_A(metal_coords);
		ObjexxFCL::FArray2D<Real> const & matrix_A =  singleset_vec_[i]->get_matrix_A();
		ObjexxFCL::FArray1D<Real> const & pcs_values = singleset_vec_[i]->get_pcs_values();
		utility::vector1< utility::vector1< utility::vector1< Vector > > > const & spin_coordinates = singleset_vec_[i]->get_spin_coordinates();
		Size n_pcs(singleset_vec_[i]->get_number_pcs());
		ObjexxFCL::FArray1D<Real> const & single_pcs_weights = singleset_vec_[i]->get_pcs_single_weights();
		Real scal(1.0 / singleset_vec_[i]->get_scaling_factor());

		// loop over number pcs for one given lanthanide
		for ( Size j = 1; j <= n_pcs; ++j ) {
			Real pcs_calc = matrix_A(j,1)*chiT_xx + matrix_A(j,2)*chiT_xy + matrix_A(j,3)*chiT_xz + matrix_A(j,4)*chiT_yy + matrix_A(j,5)*chiT_yz;
			Real pcs_exp  = pcs_values(j);
			Real diff = pcs_calc - pcs_exp;

			// PCS gradient with respect to chi tensor values
			for ( Size k = 1; k <= 5; ++k ) {
				dPCS_dparams[3 + 5*(i-1)+k] += 2.0 * diff * single_pcs_weights(j) * matrix_A(j,k);
			}

			// If symmetric_pcs_calc_ is set to true for given PCSSingleSet,
			// middle vector of spin_coordinates_ has dimension of the number
			// of symmetric subunits, otherwise the dimension is 1
			for ( Size su = 1; su <= spin_coordinates[j].size(); ++su ) {
				Real dPCS_dx_one_su(0);
				Real dPCS_dy_one_su(0);
				Real dPCS_dz_one_su(0);
				Size num_eq_spins(spin_coordinates[j][su].size());

				// loop over equivalent spins (e.g. CH3 protons, or spins on multiple subunits) that give rise to one observed PCS
				for ( Size l = 1; l <= num_eq_spins; ++l ) {
					Real x(spin_coordinates[j][su][l].x() - metal_coords.x());
					Real y(spin_coordinates[j][su][l].y() - metal_coords.y());
					Real z(spin_coordinates[j][su][l].z() - metal_coords.z());
					Real x2(x * x);
					Real y2(y * y);
					Real z2(z * z);

					// PCS prefactor
					Real r2( x2 + y2 + z2);
					Real r3( r2 * std::sqrt(r2));
					Real r5( r2 * r3 );
					Real value_1_4_PI_r5(10000.0 / (4.0 * numeric::constants::d::pi * r5));
					Real value_4_PI_r3((4.0 * numeric::constants::d::pi * r3)/10000.0);
					Real pcs_calc_one_su = scal * value_1_4_PI_r5 * ( chiT_xx*(x2-z2) + chiT_xy*(2.0*x*y) + chiT_xz*(2.0*x*z) + chiT_yy*(y2-z2) + chiT_yz*(2.0*y*z) );

					// PCS gradient with respect to metal coordinates
					dPCS_dx_one_su  += (5.0 * pcs_calc_one_su * value_4_PI_r3 * x - scal * 2.0 * (chiT_xx*x + chiT_xy*y + chiT_xz*z))
						* 2.0 * diff * single_pcs_weights(j) * value_1_4_PI_r5;
					dPCS_dy_one_su  += (5.0 * pcs_calc_one_su * value_4_PI_r3 * y - scal * 2.0 * (chiT_xy*x + chiT_yy*y + chiT_yz*z))
						* 2.0 * diff * single_pcs_weights(j) * value_1_4_PI_r5;
					dPCS_dz_one_su  += (5.0 * pcs_calc_one_su * value_4_PI_r3 * z - scal * 2.0 * (chiT_xz*x + chiT_yz*y + (-chiT_xx-chiT_yy)*z))
						* 2.0 * diff * single_pcs_weights(j) * value_1_4_PI_r5;
				} // number equivalent spins
				if ( singleset_vec_[i]->get_averaging_type() == core::scoring::nmr::MEAN ) {
					dPCS_dx_one_su /= num_eq_spins;
					dPCS_dy_one_su /= num_eq_spins;
					dPCS_dz_one_su /= num_eq_spins;
				}
				dPCS_dparams[1] += dPCS_dx_one_su;
				dPCS_dparams[2] += dPCS_dy_one_su;
				dPCS_dparams[3] += dPCS_dz_one_su;
			} // number symmetric subunits
		} // number pcs per singleset
	} // number lanthanides
}

} // namespace pcs
} // namespace nmr
} // namespace protocols

