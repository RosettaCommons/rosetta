// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/rdc/RDCMultiSet.cc
/// @brief   Implementation of class RDCMultiSet
/// @details last Modified: 07/29/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

// Unit headers
#include <core/scoring/nmr/rdc/RDCMultiSet.hh>

// Package headers
#include <core/scoring/nmr/rdc/RDCSingle.hh>
#include <core/scoring/nmr/rdc/RDCSingleSet.hh>
#include <core/scoring/nmr/rdc/RDCTensor.hh>
#include <core/scoring/nmr/util.hh>
#include <core/scoring/nmr/rdc/parameters.hh>
#include <core/io/nmr/util.hh>
#include <basic/svd/SVD_Solver.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/id/AtomID.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/pose/symmetry/util.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/nmr.OptionKeys.gen.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/fixedsizearray1.hh>

// Numeric headers
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/HomogeneousTransform.hh>
#include <numeric/constants.hh>
#include <numeric/nls/lmmin.hh>
#include <numeric/random/random.hh>

// ObjexxFCL headers
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray1D.hh>

// C++ headers
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include <limits>

// Boost headers
#include <boost/algorithm/string.hpp>

namespace core {
namespace scoring {
namespace nmr {
namespace rdc {

static basic::Tracer TR( "core.scoring.nmr.rdc.RDCMultiSet" );


/// @brief construct from data files
///        set default values for RDCMultiSet weight and computation_type
RDCMultiSet::RDCMultiSet(
	utility::vector1<std::string> const & filenames,
	std::string const & medium_label,
	pose::Pose const & pose
) :
	alignment_medium_(medium_label),
	weight_(1.0),
	number_experiments_(filenames.size()),
	svd_solver_(nullptr),
	tensor_( new RDCTensor ),
	computation_type_(RDCMultiSet::SVD),
	norm_type_(NORM_TYPE_NH),
	symmetric_rdc_calc_(false),
	correct_sign_(false),
	ave_type_(MEAN),
	fixed_tensor_(false)
{
	register_options();
	init_from_cml();
	init_from_rdc_filedata(filenames, pose);
}

/// @brief constructor with full argument list
RDCMultiSet::RDCMultiSet(
	utility::vector1<std::string> const & filenames,
	std::string const & medium_label,
	pose::Pose const & pose,
	Real const weight,
	std::string computation_type
) :
	alignment_medium_(medium_label),
	weight_(weight),
	number_experiments_(filenames.size()),
	svd_solver_(nullptr),
	tensor_( new RDCTensor ),
	norm_type_(NORM_TYPE_NH),
	symmetric_rdc_calc_(false),
	correct_sign_(false),
	ave_type_(MEAN),
	fixed_tensor_(false)
{
	register_options();
	init_from_cml();
	convert_string_to_computation_type(computation_type);
	init_from_rdc_filedata(filenames, pose);
}

/// @brief copy constructor
RDCMultiSet::RDCMultiSet(RDCMultiSet const & other) :
	alignment_medium_(other.alignment_medium_),
	weight_(other.weight_),
	total_number_rdc_(other.total_number_rdc_),
	number_experiments_(other.number_experiments_),
	matrix_A_(other.matrix_A_),
	rdc_values_(other.rdc_values_),
	rdc_single_weights_(other.rdc_single_weights_),
	svd_solver_( other.svd_solver_ ? new basic::svd::SVD_Solver(*(other.svd_solver_)) : nullptr ),
	tensor_( other.tensor_ ? new RDCTensor(*(other.tensor_)) : nullptr ),
	computation_type_(other.computation_type_),
	spin_coordinates_(other.spin_coordinates_),
	norm_type_(other.norm_type_),
	symmetric_rdc_calc_(other.symmetric_rdc_calc_),
	correct_sign_(other.correct_sign_),
	ave_type_(other.ave_type_),
	nls_repeats_(other.nls_repeats_),
	fixed_tensor_(other.fixed_tensor_)
{
	// Make deep copy of RDCSingleSet vector
	deep_copy_rdc_single_set_vec(other.rdc_singleset_vec_);
}

/// @brief assignment operator
RDCMultiSet&
RDCMultiSet::operator=(RDCMultiSet const & rhs) {
	if ( this != &rhs ) {
		alignment_medium_ = rhs.alignment_medium_;
		weight_ = rhs.weight_;
		total_number_rdc_ = rhs.total_number_rdc_;
		number_experiments_ = rhs.number_experiments_;
		matrix_A_ = rhs.matrix_A_;
		rdc_values_ = rhs.rdc_values_;
		rdc_single_weights_ = rhs.rdc_single_weights_;
		svd_solver_ = rhs.svd_solver_ ? basic::svd::SVD_SolverOP( new basic::svd::SVD_Solver(*(rhs.svd_solver_)) ) : nullptr;
		tensor_ = rhs.tensor_ ? RDCTensorOP( new RDCTensor(*(rhs.tensor_)) ) : nullptr;
		computation_type_ = rhs.computation_type_;
		spin_coordinates_ = rhs.spin_coordinates_;
		norm_type_ = rhs.norm_type_;
		symmetric_rdc_calc_ = rhs.symmetric_rdc_calc_;
		correct_sign_ = rhs.correct_sign_;
		ave_type_ = rhs.ave_type_;
		nls_repeats_ = rhs.nls_repeats_;
		fixed_tensor_ = rhs.fixed_tensor_;

		// Make deep copy of PCSSingleSet vector
		rdc_singleset_vec_.clear();
		deep_copy_rdc_single_set_vec(rhs.rdc_singleset_vec_);
	}
	return *this;
}

/// @brief destructor
RDCMultiSet::~RDCMultiSet() {}

/// @brief utility function used in constructor to initialize RDCMultiSet object from data files.
void
RDCMultiSet::init_from_rdc_filedata(
	utility::vector1<std::string> const & filenames,
	pose::Pose const & pose
)
{
	using pose::symmetry::is_symmetric;
	using pose::symmetry::symmetry_info;

	// We split behavior depending on the value of num_subunits
	// If symmetric_rdc_calc_ is true, we consider only the RDCs for the asymmetric subunit
	// and they are filled in the matrix_A_, rdc_values_ etc. times num_subunits
	// If symmetric_rdc_calc_ is false, RDCs for all subunits (i.e. for symmetric poses also for the
	// symmetric subunits, in case of monomeric protein only for the ASU) have to be provided and
	// we fill them all into matrix_A_, rdc_values_ etc. but only once

	Size num_subunits(1);
	if ( symmetric_rdc_calc_ && is_symmetric( pose ) ) {
		conformation::symmetry::SymmetryInfoCOP syminfo_ptr = symmetry_info( pose );
		num_subunits = syminfo_ptr->subunits();
	}

	number_experiments_ = filenames.size();
	rdc_singleset_vec_.reserve(filenames.size());
	total_number_rdc_ = 0;
	for ( auto const & f : filenames ) {
		RDCSingleSetOP singleset_ptr( new RDCSingleSet(f, pose) );
		rdc_singleset_vec_.push_back(singleset_ptr);
		total_number_rdc_ += singleset_ptr->get_number_rdc();
	}

	// Resize matrix, vectors and SVDSolver to the total number of RDCs (times the number of symmetric subunits if true)
	Size total_length = total_number_rdc_ * num_subunits;
	matrix_A_.dimension(total_length, 5);
	// initialize SVD solver later, only when we use SVD and update matrix A
	// in this way, we can still create a RDCMultiSet from a sparse RDC dataset
	// and calculate the RDC from a fixed tensor
	//svd_solver_ = basic::svd::SVD_SolverOP( new basic::svd::SVD_Solver(total_length, 5) );
	rdc_values_.dimension(total_length);
	rdc_single_weights_.dimension(total_length);
	spin_coordinates_.resize(total_length);

	// Fill the rdc value array
	// Calculate and fill the single weights array
	Size index_offset(0);
	for ( Size i = 1; i <= number_experiments_; ++i ) {
		Real max_rdc(0);          // largest RDC of current experiment i
		Size no_rdcs = rdc_singleset_vec_[i]->get_number_rdc(); // number of RDCs of current experiment i; increments the index_offset (times the number of symmetric subunits)
		utility::vector1<RDCSingle> single_rdcs = rdc_singleset_vec_[i]->get_single_rdc_vec();

		// Fill rdc value vector
		for ( Size su = 1; su <= num_subunits; ++su ) {
			for ( Size j = 1; j <= no_rdcs; ++j ) {

				// Normalize RDCs relative to NH or CH if needed and update values in RDCSingleSet.
				// Otherwise we assume RDC values that have already been normalized.
				Real single_rdc_val(single_rdcs[j].get_rdc_exp());
				if ( norm_type_ == NORM_TYPE_NH ) {
					single_rdc_val /= rdc_scaling_factor_toNH(rdc_singleset_vec_[i]->get_rdc_type());
				} else if ( norm_type_ == NORM_TYPE_CH ) {
					single_rdc_val /= rdc_scaling_factor_toCH(rdc_singleset_vec_[i]->get_rdc_type());
				}
				rdc_values_( index_offset + (su-1)*no_rdcs + j ) = single_rdc_val;
				rdc_singleset_vec_[i]->rdc_single_vec_[j].set_rdc_exp(single_rdc_val);
				if ( std::abs(single_rdc_val) > max_rdc ) { max_rdc = std::abs(single_rdc_val); }

				// Resize the inner vector of the spin_coordinates array to the number of equivalent spins
				spin_coordinates_[ index_offset + (su-1)*no_rdcs + j ].resize(single_rdcs[j].get_spinsAB().size());
			} // number RDCs

			// Now calculate the single rdc weights and fill single weights vector
			Real single_rdc_err(1.0);
			Real single_rdc_weight(1.0);
			for ( Size j = 1; j <= no_rdcs; ++j ) {

				single_rdc_err = single_rdcs[j].get_rdc_err();
				// Normalize RDC errors relative to NH or CH if needed and update values in RDCSingleSet.
				if ( norm_type_ == NORM_TYPE_NH ) {
					single_rdc_err /= rdc_scaling_factor_toNH(rdc_singleset_vec_[i]->get_rdc_type());
				} else if ( norm_type_ == NORM_TYPE_CH ) {
					single_rdc_err /= rdc_scaling_factor_toCH(rdc_singleset_vec_[i]->get_rdc_type());
				}
				rdc_singleset_vec_[i]->rdc_single_vec_[j].set_rdc_err(single_rdc_err);

				if ( rdc_singleset_vec_[i]->get_single_rdc_weighting_scheme() == CONST ) {
					single_rdc_weight = 1.0;
				} else {
					runtime_assert_msg(single_rdc_err > 1e-6, "ERROR in calculating single RDC weights. Experimental RDC error is unreasonably small (< 1e-6) and would produce a large weight (> 10e+12) for the chosen weighting scheme. Check RDC datafile.");
					if ( rdc_singleset_vec_[i]->get_single_rdc_weighting_scheme() == SIGMA ) {
						single_rdc_weight = 1.0/(single_rdc_err * single_rdc_err);
					} else if ( rdc_singleset_vec_[i]->get_single_rdc_weighting_scheme() == OBSIG ) {
						single_rdc_weight = std::abs(rdc_values_( index_offset + (su-1)*no_rdcs + j ))/(single_rdc_err * single_rdc_err * max_rdc);
					}
				}
				rdc_single_weights_( index_offset + (su-1)*no_rdcs + j ) = single_rdc_weight;
				rdc_singleset_vec_[i]->rdc_single_vec_[j].set_weight(single_rdc_weight);
			} // number RDCs

		} // number subunits

		// increment the index offset
		index_offset += no_rdcs * num_subunits;

	} // number of experiments (i.e. number of RDCSingleSets)
	//svd_solver_->set_vector_b(rdc_values_); // lazy creation of SVDSolver (see above)
}

/// @brief updates the spin coordinates every time the pose is changed
///        make sure that this function is called before update_matrix_A() is called
void
RDCMultiSet::update_spin_coordinates(pose::Pose const & pose) {
	using pose::symmetry::is_symmetric;
	using pose::symmetry::symmetry_info;

	Size index_offset(0);

	// If pose is symmetric and calculation with automatic deduction of symmetric RDCs is true,
	// then we pull out automatically the coordinates of the symmetric clones.
	Size num_subunits(1);
	conformation::symmetry::SymmetryInfoCOP syminfo_ptr;
	if ( symmetric_rdc_calc_ && is_symmetric( pose ) ) {
		syminfo_ptr = symmetry_info( pose );
		num_subunits = syminfo_ptr->subunits();
	}

	// The outer vector runs over No. experiments * No. RDCs (* No. subunits)
	for ( Size i = 1; i <= number_experiments_; ++i ) {
		Size no_rdcs = rdc_singleset_vec_[i]->get_number_rdc();

		for ( Size su = 1; su <= num_subunits; ++su ) {

			for ( Size j = 1; j <= no_rdcs; ++j ) {
				// Number of equivalent spins over which we perform averaging of the RDC
				Size num_eq_spins(rdc_singleset_vec_[i]->get_single_rdc_vec()[j].get_spinsAB().size());

				for ( Size k = 1; k <= num_eq_spins; ++k ) {

					Size rsdA  = rdc_singleset_vec_[i]->get_single_rdc_vec()[j].get_spinsAB()[k].first.rsd();
					Size rsdB  = rdc_singleset_vec_[i]->get_single_rdc_vec()[j].get_spinsAB()[k].second.rsd();
					Size atomA = rdc_singleset_vec_[i]->get_single_rdc_vec()[j].get_spinsAB()[k].first.atomno();
					Size atomB = rdc_singleset_vec_[i]->get_single_rdc_vec()[j].get_spinsAB()[k].second.atomno();
					// We prepare a vector in case we pull out the symmetric clones too.
					utility::vector1< Size > spinA_rsds_per_rdc_all_subunits(1, rsdA);
					utility::vector1< Size > spinB_rsds_per_rdc_all_subunits(1, rsdB);

					if ( symmetric_rdc_calc_ && is_symmetric( pose ) ) {
						utility::vector1< Size > symm_spinA_rsds = syminfo_ptr->bb_clones(rsdA);
						utility::vector1< Size > symm_spinB_rsds = syminfo_ptr->bb_clones(rsdB);
						runtime_assert_msg(symm_spinA_rsds.size() == symm_spinB_rsds.size(),
							"ERROR while updating RDC spin coordinates. Vectors of symmetric spins A and B have unequal length.");
						spinA_rsds_per_rdc_all_subunits.resize(num_subunits);
						spinB_rsds_per_rdc_all_subunits.resize(num_subunits);
						std::copy(symm_spinA_rsds.begin(), symm_spinA_rsds.end(), spinA_rsds_per_rdc_all_subunits.begin()+1);
						std::copy(symm_spinB_rsds.begin(), symm_spinB_rsds.end(), spinB_rsds_per_rdc_all_subunits.begin()+1);

					}

					spin_coordinates_[index_offset + (su-1)*no_rdcs + j][k] = std::make_pair(pose.residue(spinA_rsds_per_rdc_all_subunits[su]).atom(atomA).xyz(),
						pose.residue(spinB_rsds_per_rdc_all_subunits[su]).atom(atomB).xyz());

				} // number equivalent spins

			} // number rdcs per experiment

		} // number subunits

		// increment the index offset
		index_offset += no_rdcs * num_subunits;

	} // number experiments
}


/// @brief builds matrix_A from the spin coordinates and hands it over to the SVD solver
void
RDCMultiSet::update_matrix_A() {
	if ( !svd_solver_ ) {
		svd_solver_ = basic::svd::SVD_SolverOP( new basic::svd::SVD_Solver(rdc_values_.size(), 5) );
		svd_solver_->set_vector_b(rdc_values_);
	}
	utility::fixedsizearray1<Real,5> A_row;
	Size index_offset(0);
	Size num_subunits = spin_coordinates_.size() / total_number_rdc_;

	// loop over the different RDCSingleSets
	// and use their specific dipolar coupling constant
	// as prefactor in the RDC expression
	// concatenate all RDC equations in a single matrix A
	Real Dmax;
	for ( Size i = 1; i <= number_experiments_; ++i ) {
		Size no_rdcs = rdc_singleset_vec_[i]->get_number_rdc();
		Dmax = rdc_D_max(rdc_singleset_vec_[i]->get_rdc_type(), correct_sign_);

		// If RDCs are normalized to NH use Dmax(NH) otherwise Dmax(CAHA)
		// this seems a bit cumbersome, but we keep track of the correct sign of Dmax
		if ( norm_type_ == NORM_TYPE_NH || norm_type_ == NORM_TYPE_NONE ) {
			Dmax /= rdc_scaling_factor_toNH(rdc_singleset_vec_[i]->get_rdc_type());
		} else if ( norm_type_ == NORM_TYPE_CH ) {
			Dmax /= rdc_scaling_factor_toCH(rdc_singleset_vec_[i]->get_rdc_type());
		}
		for ( Size su = 1; su <= num_subunits; ++su ) {
			for ( Size j = 1; j <= no_rdcs; ++j ) {
				fill_matrix_A_row(A_row, spin_coordinates_[index_offset + (su-1)*no_rdcs + j], Dmax); // loop over spins in symmetric subunits
				for ( Size k = 1; k <= 5; ++k ) {              // calculate average rdc for degenerate spins (e.g. CH3)
					matrix_A_(index_offset + (su-1)*no_rdcs + j, k) = A_row[k];
				}
			} // number RDCs per experiment

		} // number subunits

		// increment the index offset
		index_offset += no_rdcs * num_subunits;

	} // number of experiments
	svd_solver_->set_matrix_A(matrix_A_);
}

/// @brief updates the array holding the single RDC weights that are used during score calculation
///        needs to be called e.g. when the weighting scheme of one of the RDCSingleSets was changed
void
RDCMultiSet::update_single_rdc_weighting() {
	Size index_offset(0);
	Size num_subunits = spin_coordinates_.size() / total_number_rdc_;

	for ( Size i = 1; i <= number_experiments_; ++i ) {
		Size no_rdcs = rdc_singleset_vec_[i]->get_number_rdc();
		for ( Size su = 1; su <= num_subunits; ++su ) {
			// Get max RDC
			Real max_rdc(0);
			for ( Size j = 1; j <= no_rdcs; ++ j ) {
				if ( std::abs(rdc_values_( index_offset + (su-1)*no_rdcs + j )) > max_rdc ) {
					max_rdc = std::abs(rdc_values_( index_offset + (su-1)*no_rdcs + j ));
				}
			}

			// Now calculate the single RDC weights and fill single weights vector
			Real single_rdc_err(1.0);
			Real single_rdc_weight(1.0);
			for ( Size j = 1; j <= no_rdcs; ++j ) {
				if ( rdc_singleset_vec_[i]->get_single_rdc_weighting_scheme() == CONST ) {
					single_rdc_weight = 1.0;
				} else {
					single_rdc_err = rdc_singleset_vec_[i]->rdc_single_vec_[j].get_rdc_err();
					runtime_assert_msg(single_rdc_err > 1e-6,
						"ERROR in calculating single RDC weights. Experimental RDC error is unreasonably small (< 1e-6) and would produce a large weight (> 10e+12) for the chosen weighting scheme. Check RDC datafile.");
					if ( rdc_singleset_vec_[i]->get_single_rdc_weighting_scheme() == SIGMA ) {
						single_rdc_weight = 1.0/(single_rdc_err * single_rdc_err);
					} else if ( rdc_singleset_vec_[i]->get_single_rdc_weighting_scheme() == OBSIG ) {
						single_rdc_weight = std::abs(rdc_values_( index_offset + (su-1)*no_rdcs + j ))/(single_rdc_err * single_rdc_err * max_rdc);
					}
				}
				rdc_single_weights_( index_offset + (su-1)*no_rdcs + j ) = single_rdc_weight;
				rdc_singleset_vec_[i]->rdc_single_vec_[j].set_weight(single_rdc_weight);

			} // number RDCs per experiment

		} // number subunits

		// increment the index offset and clear error vector
		index_offset += no_rdcs * num_subunits;
	} // number experiments
}

/// @brief utility function to fill matrix_A that is used for SVD
void
RDCMultiSet::fill_matrix_A_row(
	utility::fixedsizearray1<Real,5> & A_row,
	utility::vector1< RDCMultiSet::SpinPairCoordinates > const & spin_coord,
	Real const & Dmax
)
{
	// Set all elements in row to zero because we use operator+=
	std::fill(A_row.begin(), A_row.end(), 0.0);
	Size num_eq_spins(spin_coord.size());

	for ( Size i = 1; i <= num_eq_spins; ++i ) {
		Real x(spin_coord[i].second.x() - spin_coord[i].first.x());
		Real y(spin_coord[i].second.y() - spin_coord[i].first.y());
		Real z(spin_coord[i].second.z() - spin_coord[i].first.z());

		if ( norm_type_ == NORM_TYPE_NH || norm_type_ == NORM_TYPE_NONE ) {
			// Scale bond vector to length of NH bond (1.041 Ang.)
			Real d = std::sqrt(x * x + y * y + z * z);
			x *= (1.041 / d);
			y *= (1.041 / d);
			z *= (1.041 / d);
		} else if ( norm_type_ == NORM_TYPE_CH ) {
			// Scale bond vector to length of CAHA bond (1.107 Ang.)
			Real d = std::sqrt(x * x + y * y + z * z);
			x *= (1.107 / d);
			y *= (1.107 / d);
			z *= (1.107 / d);
		}
		Real x2(x * x);
		Real y2(y * y);
		Real z2(z * z);

		A_row[1] += 1.5 * Dmax * (x2 - z2);
		A_row[2] += 1.5 * Dmax * (2.0 * x * y);
		A_row[3] += 1.5 * Dmax * (2.0 * x * z);
		A_row[4] += 1.5 * Dmax * (y2 - z2);
		A_row[5] += 1.5 * Dmax * (2.0 * y * z);
	}
	if ( ave_type_ == MEAN ) {
		for ( Size j = 1; j <= 5; ++j ) {
			A_row[j] /= num_eq_spins;
		}
	}
}

/// @brief utility function that calculates one single RDC value given the input arguments
///        Da, R, the spin coordinates, the maximal dipolar coupling constant and a rotation matrix
Real
RDCMultiSet::basic_rdc_equation(
	utility::fixedsizearray1<Real,2> const & par,
	SpinPairCoordinates const & spin_coord,
	Real const & Dmax,
	Matrix const & rotM
)
{
	// vector between spins A and B
	Real x(spin_coord.second.x() - spin_coord.first.x());
	Real y(spin_coord.second.y() - spin_coord.first.y());
	Real z(spin_coord.second.z() - spin_coord.first.z());

	if ( norm_type_ == NORM_TYPE_NH || norm_type_ == NORM_TYPE_NONE ) {
		// Scale bond vector to length of NH bond (1.041 Ang.)
		Real d = std::sqrt(x * x + y * y + z * z);
		x *= (1.041 / d);
		y *= (1.041 / d);
		z *= (1.041 / d);
	} else if ( norm_type_ == NORM_TYPE_CH ) {
		// Scale bond vector to length of CAHA bond (1.107 Ang.)
		Real d = std::sqrt(x * x + y * y + z * z);
		x *= (1.107 / d);
		y *= (1.107 / d);
		z *= (1.107 / d);
	}

	// transformed vector after rotation
	Real x_t(rotM(1,1)*x + rotM(1,2)*y + rotM(1,3)*z);
	Real y_t(rotM(2,1)*x + rotM(2,2)*y + rotM(2,3)*z);
	Real z_t(rotM(3,1)*x + rotM(3,2)*y + rotM(3,3)*z);

	// Da = par[1], R = par[2]
	// If RDCs are normalized to NH we report the alignment magnitude relative to NH.
	// If RDCs are normalized to CH we report the alignment magnitude relative to CH.
	// Thus, we have to convert it back.
	Real Dmax_ = norm_type_ == NORM_TYPE_CH ? rdc_D_max(RDC_TYPE_CAHA, correct_sign_) : rdc_D_max(RDC_TYPE_NH, correct_sign_);
	Real A_zz = (4.0 * par[1]) / (3.0 * Dmax_);
	Real A_xx = par[1]/Dmax_ * (-2.0/3.0 + par[2]);
	Real A_yy = par[1]/Dmax_ * (-2.0/3.0 - par[2]);

	Real rdc(1.5 * Dmax * (A_xx * (x_t * x_t) + A_yy * (y_t * y_t) + A_zz * (z_t * z_t)));
	return rdc;
}

/// @brief utility functions to calculate the RDC from different sets of input arguments
///        arguments are a pointer to the parameters (Da, R), the coordinates of equivalent spins,
///        the maximal dipolar coupling constant and a rotation matrix that must be previously
///        constructed from the Euler angles. Fixed values of Da and R can be provided.
Real
RDCMultiSet::frdc(
	Real const *par,
	utility::vector1< SpinPairCoordinates > const & spin_coord,
	Real const & Dmax,
	Matrix const & rotM
)
{
	utility::fixedsizearray1<Real,2> params = { par[0], par[1] };

	Size num_eq_spins(spin_coord.size());
	Real rdc(0);
	for ( Size i = 1; i <= num_eq_spins; ++i ) {
		rdc += basic_rdc_equation(params, spin_coord[i], Dmax, rotM);
	}
	if ( ave_type_ == MEAN ) {
		rdc /= num_eq_spins;
	}
	return rdc;
}

Real
RDCMultiSet::frdc_Da(
	Real const *par,
	Real const & Da,
	utility::vector1< SpinPairCoordinates > const & spin_coord,
	Real const & Dmax,
	Matrix const & rotM
)
{
	utility::fixedsizearray1<Real,2> params = { Da, par[0] };

	Size num_eq_spins(spin_coord.size());
	Real rdc(0);
	for ( Size i = 1; i <= num_eq_spins; ++i ) {
		rdc += basic_rdc_equation(params, spin_coord[i], Dmax, rotM);
	}
	if ( ave_type_ == MEAN ) {
		rdc /= num_eq_spins;
	}
	return rdc;
}

Real
RDCMultiSet::frdc_R(
	Real const *par,
	Real const & R,
	utility::vector1< SpinPairCoordinates > const & spin_coord,
	Real const & Dmax,
	Matrix const & rotM
)
{
	utility::fixedsizearray1<Real,2> params = { par[0], R };

	Size num_eq_spins(spin_coord.size());
	Real rdc(0);
	for ( Size i = 1; i <= num_eq_spins; ++i ) {
		rdc += basic_rdc_equation(params, spin_coord[i], Dmax, rotM);
	}
	if ( ave_type_ == MEAN ) {
		rdc /= num_eq_spins;
	}
	return rdc;
}

Real
RDCMultiSet::frdc_Da_R(
	Real const & Da,
	Real const & R,
	utility::vector1< SpinPairCoordinates > const & spin_coord,
	Real const & Dmax,
	Matrix const & rotM
)
{
	utility::fixedsizearray1<Real,2> params = { Da, R };

	Size num_eq_spins(spin_coord.size());
	Real rdc(0);
	for ( Size i = 1; i <= num_eq_spins; ++i ) {
		rdc += basic_rdc_equation(params, spin_coord[i], Dmax, rotM);
	}
	if ( ave_type_ == MEAN ) {
		rdc /= num_eq_spins;
	}
	return rdc;
}

/// @brief rdc error function used in the lmmin function
///        * par is an array of fit parameters [alpha, beta, gamma, (Da, R)]
///        * data is a pointer to the the RDCMultiSet object i.e. to all data needed
///        for RDC calculation and NLS fitting
///        * fvc is an array holding the residuals of the fit calculation
void rdc_erf(
	Real const *par,
	int m_dat,
	void const *data,
	Real *fvec,
	int */*info*/
)
{
	RDCMultiSet * rdc_multiset_ptr = static_cast< RDCMultiSet * >(const_cast< void* >(data));
	Real * nonconst_par = const_cast< Real* >(par);

	// Constrain the fit parameter to reasonable ranges.
	// Test if alpha is in range [0,360) and constrain it to that range otherwise
	if ( !( 0.0 <= par[0] ) || !( par[0] < 360.0 ) ) {
		nonconst_par[0] = 180.0*std::tanh(par[0])+180.0;
	}

	// Test if beta is in range [0,360) and constrain it to that range otherwise
	if ( !( 0.0 <= par[1] ) || !( par[1] < 360.0 ) ) {
		nonconst_par[1] = 180.0*std::tanh(par[1])+180.0;
	}

	// Test if alpha is in range [0,360) and constrain it to that range otherwise
	if ( !( 0.0 <= par[2] ) || !( par[2] < 360.0 ) ) {
		nonconst_par[2] = 180.0*std::tanh(par[2])+180.0;
	}

	// Test if R is in range (0, 2/3] and constrain it to that range otherwise
	if ( rdc_multiset_ptr->computation_type_ == RDCMultiSet::NLS ) {
		if ( !( 0.0 < par[4] ) || !( par[4] <= 2.0/3.0 ) ) {
			nonconst_par[4] = (1.0/3.0)*std::tanh(par[4])+(1.0/3.0);
		}
	} else if ( rdc_multiset_ptr->computation_type_ == RDCMultiSet::NLSDA ) {
		if ( !( 0.0 < par[3] ) || !( par[3] <= 2.0/3.0 ) ) {
			nonconst_par[3] = (1.0/3.0)*std::tanh(par[3])+(1.0/3.0);
		}
	}

	Vector euler_angles(nonconst_par[0], nonconst_par[1], nonconst_par[2]); // euler angles
	Matrix rotM = rotation_matrix_from_euler_angles(euler_angles, rdc_multiset_ptr->tensor_->get_euler_convention());

	// upper index range for the individual experiments in the rdc value and single weights array
	// upper index is the partial sum of the number of rdcs per experiment (times the number of subunits)
	utility::vector1< Size > upper_index(rdc_multiset_ptr->number_experiments_);
	Size num_subunits = rdc_multiset_ptr->spin_coordinates_.size() / rdc_multiset_ptr->total_number_rdc_;
	for ( Size i = 1; i <= rdc_multiset_ptr->number_experiments_; ++i ) {
		upper_index[i] = rdc_multiset_ptr->rdc_singleset_vec_[i]->get_number_rdc() * num_subunits;
		if ( i > 1 ) { upper_index[i] +=  upper_index[i-1]; }
	}
	runtime_assert_msg(upper_index.back() == rdc_multiset_ptr->spin_coordinates_.size(), "ERROR in RDC error function. Upper index of RDC array is unequal the total number of RDCs.");
	Real Dmax;

	for ( Size i = 1, j = 1, i_end = static_cast< Size >(m_dat); i <= i_end; ++i ) { // loop over rdc_values_ array (i.e. all rdc values)
		// update maximal dipolar coupling constant for next experiment
		if ( i > upper_index[j] ) { ++j; }
		Dmax = rdc_D_max(rdc_multiset_ptr->rdc_singleset_vec_[j]->get_rdc_type(), rdc_multiset_ptr->correct_sign_);

		// if RDCs are normalized to N-H use Dmax(NH) otherwise Dmax(CH)
		if ( rdc_multiset_ptr->norm_type_ == NORM_TYPE_NH || rdc_multiset_ptr->norm_type_ == NORM_TYPE_NONE ) {
			Dmax /= rdc_scaling_factor_toNH(rdc_multiset_ptr->rdc_singleset_vec_[j]->get_rdc_type());
		} else if ( rdc_multiset_ptr->norm_type_ == NORM_TYPE_CH ) {
			Dmax /= rdc_scaling_factor_toCH(rdc_multiset_ptr->rdc_singleset_vec_[j]->get_rdc_type());
		}

		fvec[i-1] = rdc_multiset_ptr->rdc_values_(i);
		// calculate average rdc for a group of degenerate spins (e.g. methyl protons) and calculate the residual
		if ( rdc_multiset_ptr->computation_type_ == RDCMultiSet::NLS ) {
			fvec[i-1] -= rdc_multiset_ptr->frdc(&nonconst_par[3], rdc_multiset_ptr->spin_coordinates_[i], Dmax, rotM);
		} else if ( rdc_multiset_ptr->computation_type_ == RDCMultiSet::NLSDA ) {
			fvec[i-1] -= rdc_multiset_ptr->frdc_Da(&nonconst_par[3], rdc_multiset_ptr->tensor_->get_Da(), rdc_multiset_ptr->spin_coordinates_[i], Dmax, rotM);
		} else if ( rdc_multiset_ptr->computation_type_ == RDCMultiSet::NLSR ) {
			fvec[i-1] -= rdc_multiset_ptr->frdc_R(&nonconst_par[3], rdc_multiset_ptr->tensor_->get_R(), rdc_multiset_ptr->spin_coordinates_[i], Dmax, rotM);
		} else if ( rdc_multiset_ptr->computation_type_ == RDCMultiSet::NLSDAR ) {
			fvec[i-1] -= rdc_multiset_ptr->frdc_Da_R(rdc_multiset_ptr->tensor_->get_Da(), rdc_multiset_ptr->tensor_->get_R(),
				rdc_multiset_ptr->spin_coordinates_[i], Dmax, rotM);
		}
		// Weight by single RDC value weight
		//fvec[i-1] *= std::sqrt(rdc_multiset_ptr->rdc_single_weights_(i));
		fvec[i-1] *= rdc_multiset_ptr->rdc_single_weights_(i);
	}
}

/// @brief solves the RDC tensor using SVD and returns the weighted RDC score
///        according to the single RDC weighting scheme and the weights of the individual experiments
Real
RDCMultiSet::solve_tensor_and_compute_score_by_svd() {
	// update matrix A to be sure that it is set
	update_matrix_A();

	// decompose matrix A and get solution x for SLE Ax = b
	//runtime_assert_msg(svd_solver_, "ERROR while trying to decompose matrix A in RDCMultiSet. SVD_Solver has not been initialized yet.");
	svd_solver_->run_decomp_svd();
	svd_solver_->run_solve_svd();
	utility::vector1<Real> Saupe(svd_solver_->get_svd_solution());

	// tensor elements are:
	// Saupe[1] = S_xx
	// Saupe[2] = S_xy
	// Saupe[3] = S_xz
	// Saupe[4] = S_yy
	// Saupe[5] = S_yz

	tensor_->set_tensor_in_arbitrary_frame(Saupe);

	// Calculate RDC score taking into account the weights of the single RDC values and the individual experiments
	Real total_score(0);
	Size index_offset(0);
	Size num_subunits = spin_coordinates_.size() / total_number_rdc_;
	for ( Size i = 1; i <= number_experiments_; ++i ) {
		Size no_rdcs(rdc_singleset_vec_[i]->get_number_rdc());
		Real singleset_score(0);

		for ( Size su = 1; su <= num_subunits; ++su ) {
			for ( Size j = 1; j <= no_rdcs; ++j ) {
				Real calc_rdc(0);
				for ( Size k = 1; k <= 5; ++k ) {
					calc_rdc += matrix_A_( index_offset + (su-1)*no_rdcs + j, k ) * Saupe[k];
				}
				rdc_singleset_vec_[i]->rdc_single_vec_[j].set_rdc_calc(calc_rdc);
				Real diff = calc_rdc - rdc_values_( index_offset + (su-1)*no_rdcs + j );
				singleset_score += diff * diff * rdc_single_weights_( index_offset + (su-1)*no_rdcs + j );
			} // number RDCs per experiment
		} // number subunits

		singleset_score /= num_subunits;
		//total_score += std::sqrt(singleset_score) * rdc_singleset_vec_[i]->get_weight();
		total_score += singleset_score * rdc_singleset_vec_[i]->get_weight();
		index_offset += no_rdcs * num_subunits;
	} // number experiments

	// Show score and tensor after svd
	if ( TR.Trace.visible() ) {
		TR.Trace << "RDC score and tensor for alignment medium " << alignment_medium_ << " after SVD: " << total_score << std::endl;
		tensor_->show_tensor_stats(TR.Trace, false);
	}

	return total_score;
}

/// @brief solves the RDC tensor using NLS and returns the weighted RDC score
///        according to the single RDC weighting scheme and the weights of the individual experiments
Real
RDCMultiSet::solve_tensor_and_compute_score_by_nls() {

	// Set initial values for parameters
	Size number_params(3);
	utility::vector1<Real> params;
	utility::vector1<Real> params_best;

	for ( Size i = 1; i <= 3; ++i ) {
		params.push_back(360.0 * numeric::random::uniform());
	}

	if ( computation_type_ == RDCMultiSet::NLS ) {
		number_params = 5;
		params.push_back(-9999);
		params.push_back(-9999);
		params_best.resize(5);
	} else if ( (computation_type_ == RDCMultiSet::NLSDA) || (computation_type_ == RDCMultiSet::NLSR) ) {
		number_params = 4;
		params.push_back(-9999);
		params_best.resize(4);
	} else if ( computation_type_ == RDCMultiSet::NLSDAR ) {
		number_params = 3;
		params_best.resize(3);
	} else if ( computation_type_ == RDCMultiSet::SVD ) {
		utility_exit_with_message("ERROR in RDC tensor and score calculation. Cannot use NLS fitting because specified computation type is SVD.");
	}
	params_best=params; // Set best fit params to initial guess to avoid that they are not initialized in case lmmin fails

	// definition of auxiliary parameters for lmmin
	numeric::nls::lm_status_struct status;
	Real bestnorm = std::numeric_limits< Real >::infinity();

	// Perform nonlinear least squares fitting
	for ( Size i = 1; i <= nls_repeats_; ++i ) {
		TR.Trace << "RDC NLS fitting" << std::endl;
		// Random starting values
		TR.Trace << "Initialize angles to random values:";
		for ( Size j = 1; j <= 3; ++j ) {
			params[j] = 360.0 * numeric::random::uniform();
			TR.Trace << " " << params[j];
		}
		TR.Trace << std::endl;
		numeric::nls::lmmin( number_params, &params[1], rdc_values_.size1(), (const void*) this, rdc_erf, &status, numeric::nls::lm_printout_std);
		TR.Trace << "Iteration: " << i << " status.fnorm: " << status.fnorm << " bestnorm: "<< bestnorm << std::endl;
		//save to best fitting parameter
		if ( status.fnorm < bestnorm ) {
			TR.Trace << "status.fnorm: " << status.fnorm << " replaced bestnorm: " << bestnorm << std::endl;
			bestnorm=status.fnorm;
			params_best = params;
		}
	}

	// Set tensor with best parameters obtained from fit
	// Alignment magnitude is scaled relative to NH or CH.
	// Thus, we add the maximal dipolar coupling constant to the tensor parameters.
	if ( computation_type_ == RDCMultiSet::NLS ) {
		params_best.resize(6);
	} else if ( computation_type_ == RDCMultiSet::NLSDA ) {
		params_best.resize(6);
		params_best[5] = tensor_->get_Da();
		std::swap(params_best[4], params_best[5]);
	} else if ( computation_type_ == RDCMultiSet::NLSR ) {
		params_best.resize(6);
		params_best[5] = tensor_->get_R();
	} else if ( computation_type_ == RDCMultiSet::NLSDAR ) {
		params_best.resize(6);
		params_best[4] = tensor_->get_Da();
		params_best[5] = tensor_->get_R();
	}
	params_best[6] = (norm_type_ == NORM_TYPE_CH) ? rdc_D_max(RDC_TYPE_CAHA, correct_sign_) : rdc_D_max(RDC_TYPE_NH, correct_sign_);
	tensor_->set_tensor_in_pas(params_best);

	Vector euler_angles(params_best[1], params_best[2], params_best[3]); // euler angles
	Matrix rotM = rotation_matrix_from_euler_angles(euler_angles, tensor_->get_euler_convention());

	// Calculate RDC score taking into account the weights of the single RDC values and the individual experiments
	Real total_score(0);
	Size index_offset(0);
	Size num_subunits = spin_coordinates_.size() / total_number_rdc_;
	Real Dmax;
	for ( Size i = 1; i <= number_experiments_; ++i ) {
		Size no_rdcs(rdc_singleset_vec_[i]->get_number_rdc());
		Real singleset_score(0);
		Dmax = rdc_D_max(rdc_singleset_vec_[i]->get_rdc_type(), correct_sign_);
		// if RDCs are normalized to N-H use Dmax(NH) otherwise Dmax(CH)
		if ( norm_type_ == NORM_TYPE_NH || norm_type_ == NORM_TYPE_NONE ) {
			Dmax /= rdc_scaling_factor_toNH(rdc_singleset_vec_[i]->get_rdc_type());
		} else if ( norm_type_ == NORM_TYPE_CH ) {
			Dmax /= rdc_scaling_factor_toCH(rdc_singleset_vec_[i]->get_rdc_type());
		}

		for ( Size su = 1; su <= num_subunits; ++su ) {    // loop over spins in symmetric subunits
			for ( Size j = 1; j <= no_rdcs; ++j ) {     // calculate average rdc for a group of degenerate spins
				Real calc_rdc = frdc(&params_best[4], spin_coordinates_[ index_offset + (su-1)*no_rdcs + j ], Dmax, rotM);
				rdc_singleset_vec_[i]->rdc_single_vec_[j].set_rdc_calc(calc_rdc);
				Real diff = calc_rdc - rdc_values_(index_offset + (su-1)*no_rdcs + j);
				singleset_score += diff * diff * rdc_single_weights_(index_offset + (su-1)*no_rdcs + j);
			} // number RDCs per experiment
		} // number subunits

		singleset_score /= num_subunits;
		//total_score += std::sqrt(singleset_score) * rdc_singleset_vec_[i]->get_weight();
		total_score += singleset_score * rdc_singleset_vec_[i]->get_weight();
		index_offset += no_rdcs * num_subunits;
	} // number of experiments

	// Show score and tensor after fit
	if ( TR.Trace.visible() ) {
		TR.Trace << "RDC score and tensor for alignment medium " << alignment_medium_ << " after NLS: " << total_score << std::endl;
		tensor_->show_tensor_stats(TR.Trace, true);
	}

	return total_score;
}

/// @brief sets the xyz derivative of the RDC
///        RDCTensor must be determined before
///        first call solve_tensor_and_compute_score_by_svd() or
///        solve_tensor_and_compute_score_by_nls() before setting derivatives
void
RDCMultiSet::set_atom_derivatives(pose::Pose const & pose) {
	using pose::symmetry::is_symmetric;
	using pose::symmetry::symmetry_info;

	if ( tensor_->is_rdc_tensor_in_arbitrary_frame() ) {

		// Get the current alignment tensor
		Matrix saupe(Matrix::rows(tensor_->get_T_xx(), tensor_->get_T_xy(),  tensor_->get_T_xz(),
			tensor_->get_T_xy(), tensor_->get_T_yy(),  tensor_->get_T_yz(),
			tensor_->get_T_xz(), tensor_->get_T_yz(), -tensor_->get_T_xx()-tensor_->get_T_yy()));

		// Update matrix A
		update_matrix_A();

		// Calculate and set xyz derivatives
		Size index_offset(0);
		Size num_subunits(1);
		if ( symmetric_rdc_calc_ && is_symmetric( pose ) ) {
			conformation::symmetry::SymmetryInfoCOP syminfo_ptr = symmetry_info( pose );
			num_subunits = syminfo_ptr->subunits();
		}
		Real Dmax;
		for ( Size i = 1; i <= number_experiments_; ++i ) {
			Size no_rdcs(rdc_singleset_vec_[i]->get_number_rdc());
			// RDC prefactor
			Dmax = rdc_D_max(rdc_singleset_vec_[i]->get_rdc_type(), correct_sign_);
			// if RDCs are normalized to NH use Dmax(NH) otherwise Dmax(CH)
			if ( norm_type_ == NORM_TYPE_NH || norm_type_ == NORM_TYPE_NONE ) {
				Dmax /= rdc_scaling_factor_toNH(rdc_singleset_vec_[i]->get_rdc_type());
			} else if ( norm_type_ == NORM_TYPE_CH ) {
				Dmax /= rdc_scaling_factor_toCH(rdc_singleset_vec_[i]->get_rdc_type());
			}
			for ( Size su = 1; su <= num_subunits; ++su ) {
				// Calculate RDC per experiment
				for ( Size j = 1; j <= no_rdcs; ++j ) {
					Real calc_rdc = matrix_A_(index_offset + (su-1)*no_rdcs + j, 1) * saupe.xx()
						+ matrix_A_(index_offset + (su-1)*no_rdcs + j, 2) * saupe.xy()
						+ matrix_A_(index_offset + (su-1)*no_rdcs + j, 3) * saupe.xz()
						+ matrix_A_(index_offset + (su-1)*no_rdcs + j, 4) * saupe.yy()
						+ matrix_A_(index_offset + (su-1)*no_rdcs + j, 5) * saupe.yz();
					Real diff = calc_rdc - rdc_values_(index_offset + (su-1)*no_rdcs + j);

					// Here we iterate over the number of equivalent (degenerate) spins that produce the same rdc
					for ( Size k = 1, k_end = rdc_singleset_vec_[i]->get_single_rdc_vec()[j].get_spinsAB().size(); k <= k_end; ++k ) {

						Real x = spin_coordinates_[index_offset + (su-1)*no_rdcs + j][k].second.x()
							- spin_coordinates_[index_offset + (su-1)*no_rdcs + j][k].first.x();
						Real y = spin_coordinates_[index_offset + (su-1)*no_rdcs + j][k].second.y()
							- spin_coordinates_[index_offset + (su-1)*no_rdcs + j][k].first.y();
						Real z = spin_coordinates_[index_offset + (su-1)*no_rdcs + j][k].second.z()
							- spin_coordinates_[index_offset + (su-1)*no_rdcs + j][k].first.z();

						if ( norm_type_ == NORM_TYPE_NH || norm_type_ == NORM_TYPE_NONE ) {
							// Scale bond vector to length of NH bond (1.041 Ang.)
							Real d = std::sqrt(x * x + y * y + z * z);
							x *= (1.041 / d);
							y *= (1.041 / d);
							z *= (1.041 / d);
						} else if ( norm_type_ == NORM_TYPE_CH ) {
							// Scale bond vector to length of CAHA bond (1.107 Ang.)
							Real d = std::sqrt(x * x + y * y + z * z);
							x *= (1.107 / d);
							y *= (1.107 / d);
							z *= (1.107 / d);
						}
						Real x2(x * x);
						Real y2(y * y);
						Real z2(z * z);
						Real r2( x2 + y2 + z2);

						// atom derivatives
						Real dx = ( ( 3.0 * Dmax * (saupe.xx()*x + saupe.xy()*y + saupe.xz()*z) - calc_rdc * 3.0 * (x / r2))
							* diff * rdc_single_weights_(index_offset + (su-1)*no_rdcs + j) );
						Real dy = ( ( 3.0 * Dmax * (saupe.xy()*x + saupe.yy()*y + saupe.yz()*z) - calc_rdc * 3.0 * (y / r2))
							* diff * rdc_single_weights_(index_offset + (su-1)*no_rdcs + j) );
						Real dz = ( ( 3.0 * Dmax * (saupe.xz()*x + saupe.yz()*y + (-saupe.xx()-saupe.yy())*z) - calc_rdc * 3.0 * (z / r2))
							* diff * rdc_single_weights_(index_offset + (su-1)*no_rdcs + j) );

						rdc_singleset_vec_[i]->rdc_single_vec_[j].set_atom_derivatives(k, dx, dy, dz);

					} // equivalent spins

				} //number RDCs per experiment

			} // number subunits

			// Increment index offset
			index_offset += no_rdcs * num_subunits;
		}

	} else if ( tensor_->is_rdc_tensor_in_pas() ) {
		// Get current alignment tensor
		Vector euler_angles(tensor_->get_alpha(), tensor_->get_beta(), tensor_->get_gamma());
		Matrix rotM = rotation_matrix_from_euler_angles(euler_angles, tensor_->get_euler_convention());
		utility::fixedsizearray1<Real,2> T = { tensor_->get_Da(), tensor_->get_R() };

		// convert Da and R to tensor eigenvalues
		// Da is scaled to NH or CH
		Real Dmax_ = norm_type_ == NORM_TYPE_CH ? rdc_D_max(RDC_TYPE_CAHA, correct_sign_) : rdc_D_max(RDC_TYPE_NH, correct_sign_);
		Real A_zz = (4.0 * T[1]) / (3.0 * Dmax_);
		Real A_xx = T[1]/Dmax_ * (-2.0/3.0 + T[2]);
		Real A_yy = T[1]/Dmax_ * (-2.0/3.0 - T[2]);

		Size index_offset(0);
		Size num_subunits(1);
		if ( symmetric_rdc_calc_ && is_symmetric( pose ) ) {
			conformation::symmetry::SymmetryInfoCOP syminfo_ptr = symmetry_info( pose );
			num_subunits = syminfo_ptr->subunits();
		}
		Real Dmax;
		for ( Size i = 1; i <= number_experiments_; ++i ) {
			Size no_rdcs(rdc_singleset_vec_[i]->get_number_rdc());
			// RDC prefactor
			Dmax = rdc_D_max(rdc_singleset_vec_[i]->get_rdc_type(), correct_sign_);
			// if RDCs are normalized to NH use Dmax(NH) otherwise Dmax(CH)
			if ( norm_type_ == NORM_TYPE_NH || norm_type_ == NORM_TYPE_NONE ) {
				Dmax /= rdc_scaling_factor_toNH(rdc_singleset_vec_[i]->get_rdc_type());
			} else if ( norm_type_ == NORM_TYPE_CH ) {
				Dmax /= rdc_scaling_factor_toCH(rdc_singleset_vec_[i]->get_rdc_type());
			}
			for ( Size su = 1; su <= num_subunits; ++su ) {
				// Calculate RDC per experiment
				for ( Size j = 1; j <= no_rdcs; ++j ) {
					Real calc_rdc = frdc(&T[1], spin_coordinates_[index_offset + (su-1)*no_rdcs + j], Dmax, rotM);
					Real diff = calc_rdc - rdc_values_(index_offset + (su-1)*no_rdcs + j);

					// Here we iterate over the number of equivalent (degenerate) spins that produce the same rdc
					for ( Size k = 1, k_end = rdc_singleset_vec_[i]->get_single_rdc_vec()[j].get_spinsAB().size(); k <= k_end ; ++k ) {

						Real x = spin_coordinates_[index_offset + (su-1)*no_rdcs + j][k].second.x()
							- spin_coordinates_[index_offset + (su-1)*no_rdcs + j][k].first.x();
						Real y = spin_coordinates_[index_offset + (su-1)*no_rdcs + j][k].second.y()
							- spin_coordinates_[index_offset + (su-1)*no_rdcs + j][k].first.y();
						Real z = spin_coordinates_[index_offset + (su-1)*no_rdcs + j][k].second.z()
							- spin_coordinates_[index_offset + (su-1)*no_rdcs + j][k].first.z();

						if ( norm_type_ == NORM_TYPE_NH || norm_type_ == NORM_TYPE_NONE ) {
							// Scale bond vector to length of NH bond (1.041 Ang.)
							Real d = std::sqrt(x * x + y * y + z * z);
							x *= (1.041 / d);
							y *= (1.041 / d);
							z *= (1.041 / d);
						} else if ( norm_type_ == NORM_TYPE_CH ) {
							// Scale bond vector to length of CAHA bond (1.107 Ang.)
							Real d = std::sqrt(x * x + y * y + z * z);
							x *= (1.107 / d);
							y *= (1.107 / d);
							z *= (1.107 / d);
						}
						// transformed vector after rotation
						Real x_t(rotM(1,1)*x + rotM(1,2)*y + rotM(1,3)*z);
						Real y_t(rotM(2,1)*x + rotM(2,2)*y + rotM(2,3)*z);
						Real z_t(rotM(3,1)*x + rotM(3,2)*y + rotM(3,3)*z);
						Real r2( x_t*x_t + y_t*y_t + z_t*z_t);

						// atom derivatives
						Real dx = ( ( 3.0 * Dmax * (A_xx * rotM(1,1) * x_t + A_yy * rotM(2,1) * y_t + A_zz * rotM(3,1) * z_t)
							- calc_rdc * 3.0 * ((rotM(1,1) * x_t + rotM(2,1) * y_t + rotM(3,1) * z_t) / r2))
							* diff * rdc_single_weights_(index_offset + (su-1)*no_rdcs + j) );
						Real dy = ( ( 3.0 * Dmax * (A_xx * rotM(1,2) * x_t + A_yy * rotM(2,2) * y_t + A_zz * rotM(3,2) * z_t)
							- calc_rdc * 3.0 * ((rotM(1,2) * x_t + rotM(2,2) * y_t + rotM(3,2) * z_t) / r2))
							* diff * rdc_single_weights_(index_offset + (su-1)*no_rdcs + j) );
						Real dz = ( ( 3.0 * Dmax * (A_xx * rotM(1,3) * x_t + A_yy * rotM(2,3) * y_t + A_zz * rotM(3,3) * z_t)
							- calc_rdc * 3.0 * ((rotM(1,3) * x_t + rotM(2,3) * y_t + rotM(3,3) * z_t) / r2))
							* diff * rdc_single_weights_(index_offset + (su-1)*no_rdcs + j) );

						rdc_singleset_vec_[i]->rdc_single_vec_[j].set_atom_derivatives(k, dx, dy, dz);

					} // equivalent spins

				} // number RDCs per experiment

			}  // number subunits

			// Increment index offset
			index_offset += no_rdcs * num_subunits;
		}
	} else {
		utility_exit_with_message( "ERROR when setting RDC xyz derivative. RDCTensor is not set in arbitrary or principal axis frame. First call \"solve_tensor_and_compute_score_by_svd()\" or \"solve_tensor_and_compute_score_by_nls()\"." );
	}
}

/// @brief calculate RDC values from a given tensor, set values in the RDCSingle vector
///        and return the RDC score
Real
RDCMultiSet::compute_rdc_values_and_score_from_tensor(RDCTensor const & tensor) {

	if ( tensor.is_rdc_tensor_in_arbitrary_frame() ) {

		utility::fixedsizearray1<Real,5> Saupe = { tensor.get_T_xx(), tensor.get_T_xy(), tensor.get_T_xz(), tensor.get_T_yy(), tensor.get_T_yz() };

		utility::fixedsizearray1<Real,5> Arow;
		Real total_score(0);
		Size index_offset(0);
		Size num_subunits = spin_coordinates_.size() / total_number_rdc_;
		Real Dmax;

		for ( Size i = 1; i <= number_experiments_; ++i ) {
			Size no_rdcs(rdc_singleset_vec_[i]->get_number_rdc());
			Real singleset_score(0);
			Dmax = rdc_D_max(rdc_singleset_vec_[i]->get_rdc_type(), correct_sign_);
			if ( norm_type_ == NORM_TYPE_NH || norm_type_ == NORM_TYPE_NONE ) {
				Dmax /= rdc_scaling_factor_toNH(rdc_singleset_vec_[i]->get_rdc_type());
			} else if ( norm_type_ == NORM_TYPE_CH ) {
				Dmax /= rdc_scaling_factor_toCH(rdc_singleset_vec_[i]->get_rdc_type());
			}
			for ( Size su = 1; su <= num_subunits; ++su ) {

				for ( Size j = 1; j <= no_rdcs; ++j ) {
					Real calc_rdc(0);
					fill_matrix_A_row(Arow, spin_coordinates_[index_offset + (su-1)*no_rdcs + j], Dmax);

					for ( Size k = 1; k <= 5; ++k ) {
						calc_rdc += Arow[k] * Saupe[k];
					}
					rdc_singleset_vec_[i]->rdc_single_vec_[j].set_rdc_calc(calc_rdc);
					Real diff = calc_rdc - rdc_values_( index_offset + (su-1)*no_rdcs + j );
					singleset_score += diff * diff * rdc_single_weights_( index_offset + (su-1)*no_rdcs + j );
				}
			}
			singleset_score /= num_subunits;
			//total_score += std::sqrt(singleset_score) * rdc_singleset_vec_[i]->get_weight();
			total_score += singleset_score * rdc_singleset_vec_[i]->get_weight();
			index_offset += no_rdcs * num_subunits;
		}
		return total_score;

	} else if ( tensor.is_rdc_tensor_in_pas() ) {

		Vector euler_angles(tensor.get_alpha(), tensor.get_beta(), tensor.get_gamma()); // euler angles
		Matrix rotM = rotation_matrix_from_euler_angles(euler_angles, tensor.get_euler_convention());
		utility::fixedsizearray1<Real,2> Da_R = { tensor.get_Da(), tensor.get_R() };

		Real total_score(0);
		Size index_offset(0);
		Size num_subunits = spin_coordinates_.size() / total_number_rdc_;
		Real Dmax;
		for ( Size i = 1; i <= number_experiments_; ++i ) {
			Size no_rdcs(rdc_singleset_vec_[i]->get_number_rdc());
			Real singleset_score(0);
			Dmax = rdc_D_max(rdc_singleset_vec_[i]->get_rdc_type(), correct_sign_);
			if ( norm_type_ == NORM_TYPE_NH || norm_type_ == NORM_TYPE_NONE ) {
				Dmax /= rdc_scaling_factor_toNH(rdc_singleset_vec_[i]->get_rdc_type());
			} else if ( norm_type_ == NORM_TYPE_CH ) {
				Dmax /= rdc_scaling_factor_toCH(rdc_singleset_vec_[i]->get_rdc_type());
			}
			for ( Size su = 1; su <= num_subunits; ++su ) {
				for ( Size j = 1; j <= no_rdcs; ++j ) {
					Real calc_rdc = frdc(&Da_R[1], spin_coordinates_[ index_offset + (su-1)*no_rdcs + j ], Dmax, rotM);
					rdc_singleset_vec_[i]->rdc_single_vec_[j].set_rdc_calc(calc_rdc);
					Real diff = calc_rdc - rdc_values_(index_offset + (su-1)*no_rdcs + j);
					singleset_score += diff * diff * rdc_single_weights_(index_offset + (su-1)*no_rdcs + j);
				}
			}
			singleset_score /= num_subunits;
			//total_score += std::sqrt(singleset_score) * rdc_singleset_vec_[i]->get_weight();
			total_score += singleset_score * rdc_singleset_vec_[i]->get_weight();
			index_offset += no_rdcs * num_subunits;
		}
		return total_score;

	} else {
		utility_exit_with_message( "ERROR when trying to calculate RDC values from given tensor. RDCTensor is not set in arbitrary or principal axis frame. First call \"solve_tensor_and_compute_score_by_svd()\" or \"solve_tensor_and_compute_score_by_nls()\"." );
	}
}

/// @brief calculate RDC values and the RDC score from the dataset's current tensor
Real
RDCMultiSet::compute_rdc_values_and_score_from_tensor() {
	if ( !tensor_ ) {
		utility_exit_with_message("ERROR  while trying to calculate RDCMultiSet score from current tensor. No RDCTensor object set.");
	}
	return compute_rdc_values_and_score_from_tensor( *tensor_ );
}

/// @brief register options
void
RDCMultiSet::register_options() {
	using namespace basic::options;
	option.add_relevant(OptionKeys::nmr::rdc::nls_repeats);
	option.add_relevant(OptionKeys::nmr::rdc::normalization_type);
	option.add_relevant(OptionKeys::nmr::rdc::use_symmetry_calc);
	option.add_relevant(OptionKeys::nmr::rdc::correct_sign);
}

void
RDCMultiSet::init_from_cml() {
	using namespace basic::options;
	norm_type_          = convert_string_to_normalization_type( option[ basic::options::OptionKeys::nmr::rdc::normalization_type ]() );
	symmetric_rdc_calc_ = option[ basic::options::OptionKeys::nmr::rdc::use_symmetry_calc ]();
	correct_sign_       = option[ basic::options::OptionKeys::nmr::rdc::correct_sign ]();
	nls_repeats_        = option[ basic::options::OptionKeys::nmr::rdc::nls_repeats ]();
}

/// @brief utility function to convert string to class computation type enum
void
RDCMultiSet::convert_string_to_computation_type(std::string const & computation_type) {
	if ( computation_type == "SVD" || computation_type == "svd" ) {
		computation_type_ = RDCMultiSet::SVD;
	} else if ( computation_type == "NLS" || computation_type == "nls" ) {
		computation_type_ = RDCMultiSet::NLS;
	} else if ( computation_type == "NLSDA" || computation_type == "nlsda" ) {
		computation_type_ = RDCMultiSet::NLSDA;
	} else if ( computation_type == "NLSR" || computation_type == "nlsr" ) {
		computation_type_ = RDCMultiSet::NLSR;
	} else if ( computation_type == "NLSDAR" || computation_type == "nlsdar" ) {
		computation_type_ = RDCMultiSet::NLSDAR;
	} else {
		utility_exit_with_message( "ERROR: Provided string does not match any RDCMultiSet computation type. Possible options are \"SVD\", \"NLS\", \"NLSDA\", \"NLSR\" and \"NLSDAR\"." );
	}
}

void
RDCMultiSet::deep_copy_rdc_single_set_vec(utility::vector1< RDCSingleSetOP > const & other_vec) {
	rdc_singleset_vec_.resize(other_vec.size());
	for ( Size i = 1; i <= rdc_singleset_vec_.size(); ++i ) {
		rdc_singleset_vec_[i] = RDCSingleSetOP( new RDCSingleSet( *(other_vec[i]) ) );
	}
}

void
RDCMultiSet::set_computation_type(std::string const & type) {
	convert_string_to_computation_type(type);
}

void
RDCMultiSet::set_averaging_type(std::string const & type) {
	ave_type_ = convert_string_to_averaging_type(type);
}

void
RDCMultiSet::show(std::ostream & tracer) const {
	tracer << "   * * * RDCMultiSet Summary Report * * *   " << std::endl;
	tracer << "Alignment Medium: " << alignment_medium_ << std::endl;
	tracer << "No experiments:   " << number_experiments_ << std::endl;
	tracer << "No RDCs:          " << total_number_rdc_ << std::endl;
	tracer << "Experimental conditions: " << std::endl;
	for ( Size i = 1; i <= number_experiments_; ++i ) {
		tracer << " - " << rdc_singleset_vec_[i]->get_dataset_name() << " : [ " << rdc_singleset_vec_[i]->get_number_rdc() << " RDCs, "
			<< convert_rdc_type_to_string(rdc_singleset_vec_[i]->get_rdc_type()) << " ]" << std::endl;
	}
	tracer << "Alignment Tensor: " << std::endl;
	computation_type_ == SVD ? tensor_->show_tensor_stats(tracer, false) : tensor_->show_tensor_stats(tracer, true);
}

} // namespace rdc
} // namespace nmr
} // namespace scoring
} // namespace core
