// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/pcs/PCSSingleSet.cc
/// @brief   Implementation of class PCSSingleSet
/// @details last Modified: 06/22/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

// Unit headers
#include <core/scoring/nmr/pcs/PCSSingleSet.hh>

// Package headers
#include <core/scoring/nmr/pcs/PCSSingle.hh>
#include <core/scoring/nmr/pcs/PCSTensor.hh>
#include <core/scoring/nmr/util.hh>
#include <core/io/nmr/AtomSelection.hh>
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
#include <utility/string_util.hh>
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
#include <limits>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <numeric>

// Boost headers
#include <boost/algorithm/string.hpp>

namespace core {
namespace scoring {
namespace nmr {
namespace pcs {

static basic::Tracer TR( "core.scoring.nmr.pcs.PCSSingleSet" );

/// @brief constructor with arguments
///        set default values for PCSSingleSet weight, single_pcs_weighting_scheme and computation type
PCSSingleSet::PCSSingleSet(
	std::string const & filename,
	std::string const & metal_ion_label,
	pose::Pose const & pose
) :
	dataset_name_(utility::file_basename(filename)),
	metal_ion_label_(metal_ion_label),
	weight_(1.0),
	scaling_factor_(1.0),
	tensor_( new PCSTensor ),
	svd_solver_(nullptr),
	normalized_data_(false),
	symmetric_pcs_calc_(false),
	ave_type_(MEAN)
{
	register_options();
	init_from_cml();
	computation_type_ = PCSSingleSet::SVD;
	single_pcs_weighting_scheme_ = CONST;
	init_from_pcs_filedata(filename, pose);
}

/// @brief constructor with full argument list
PCSSingleSet::PCSSingleSet(
	std::string const & filename,
	std::string const & metal_ion_label,
	pose::Pose const & pose,
	Real const weight,
	std::string single_pcs_weigting,
	std::string computation_type
) :
	dataset_name_(utility::file_basename(filename)),
	metal_ion_label_(metal_ion_label),
	weight_(weight),
	scaling_factor_(1.0),
	tensor_( new PCSTensor ),
	svd_solver_(nullptr),
	normalized_data_(false),
	symmetric_pcs_calc_(false),
	ave_type_(MEAN)

{
	register_options();
	init_from_cml();
	convert_string_to_computation_type(computation_type);
	single_pcs_weighting_scheme_ = convert_string_to_weighting_scheme(single_pcs_weigting);
	init_from_pcs_filedata(filename, pose);
}

/// @brief copy constructor
PCSSingleSet::PCSSingleSet(PCSSingleSet const & other) :
	dataset_name_(other.dataset_name_),
	metal_ion_label_(other.metal_ion_label_),
	pcs_single_vec_(other.pcs_single_vec_),
	weight_(other.weight_),
	scaling_factor_(other.scaling_factor_),
	number_pcs_(other.number_pcs_),
	matrix_A_(other.matrix_A_),
	pcs_values_(other.pcs_values_),
	pcs_single_weights_(other.pcs_single_weights_),
	tensor_( other.tensor_ ? new PCSTensor(*(other.tensor_)) : nullptr ),
	svd_solver_( other.svd_solver_ ? new basic::svd::SVD_Solver( *(other.svd_solver_) ) : nullptr ),
	single_pcs_weighting_scheme_(other.single_pcs_weighting_scheme_),
	computation_type_(other.computation_type_),
	spin_coordinates_(other.spin_coordinates_),
	metal_coord_bounds_(other.metal_coord_bounds_),
	normalized_data_(other.normalized_data_),
	symmetric_pcs_calc_(other.symmetric_pcs_calc_),
	ave_type_(other.ave_type_),
	nls_repeats_(other.nls_repeats_)
{}

/// @brief assignment operator
PCSSingleSet&
PCSSingleSet::operator=(PCSSingleSet const & rhs) {
	if ( this != &rhs ) {
		dataset_name_ = rhs.dataset_name_;
		metal_ion_label_ = rhs.metal_ion_label_;
		pcs_single_vec_ = rhs.pcs_single_vec_;
		weight_ = rhs.weight_;
		scaling_factor_ = rhs.scaling_factor_;
		number_pcs_ = rhs.number_pcs_;
		matrix_A_ = rhs.matrix_A_;
		pcs_values_ = rhs.pcs_values_;
		pcs_single_weights_ = rhs.pcs_single_weights_;
		tensor_ = rhs.tensor_ ? PCSTensorOP( new PCSTensor(*(rhs.tensor_)) ) : nullptr;
		svd_solver_ = rhs.svd_solver_ ? basic::svd::SVD_SolverOP( new basic::svd::SVD_Solver( *(rhs.svd_solver_))) : nullptr;
		single_pcs_weighting_scheme_ = rhs.single_pcs_weighting_scheme_;
		computation_type_ = rhs.computation_type_;
		spin_coordinates_ = rhs.spin_coordinates_;
		metal_coord_bounds_ = rhs.metal_coord_bounds_;
		normalized_data_ = rhs.normalized_data_;
		symmetric_pcs_calc_ = rhs.symmetric_pcs_calc_;
		ave_type_ = rhs.ave_type_;
		nls_repeats_ = rhs.nls_repeats_;
	}
	return *this;
}

/// @brief destructor
PCSSingleSet::~PCSSingleSet() {}

/// @brief utility function used in constructor to initialize PCSSingelSet object from filedata and pose.
void
PCSSingleSet::init_from_pcs_filedata(
	std::string const & filename,
	pose::Pose const & pose
)
{
	using pose::symmetry::is_symmetric;
	using pose::symmetry::symmetry_info;

	// Get the pcs filedata and create a vector of PCSSingle objects
	// resize the spin coordinates vector to the required size
	utility::vector1< utility::vector1< core::io::nmr::AtomSelection > > spinsA_all;
	utility::vector1< Real > values;
	utility::vector1< Real > errors;
	core::io::nmr::read_pcs_datafile(filename, spinsA_all, values, errors);
	dataset_name_ = utility::file_basename(filename);
	number_pcs_ = values.size();
	spin_coordinates_.resize(number_pcs_);
	// Set up things for SVD
	matrix_A_.dimension(number_pcs_, 5);
	// initialize SVD solver later, only when we use SVD and update matrix A
	// in this way, we can still create a PCSSingleSet from a sparse PCS dataset
	// and calculate the PCS from a fixed tensor
	//svd_solver_ = basic::svd::SVD_SolverOP( new basic::svd::SVD_Solver(number_pcs_, 5) );
	pcs_values_.dimension(number_pcs_);
	pcs_single_weights_.dimension(number_pcs_);

	TR.Info << "Creating PCSSingleSet from file " << dataset_name_ << ". Number of PCSs: " << number_pcs_ << std::endl;

	if ( normalized_data_ ) {
		// Normalize PCS data if needed and set the scaling factor
		Real mean = std::accumulate(values.begin(), values.end(), 0.0) / values.size();
		Real sq_sum = 0.0;
		for ( const Real & v : values ) {
			sq_sum += (v - mean) * (v - mean);
		}
		scaling_factor_ = std::sqrt(sq_sum / values.size());
	} else {
		scaling_factor_ = 1.0;
	}
	// If we perform PCS calculation with automatic deduction of protein symmetry,
	// we resize the middle vector of the spin_coordinates_ array to the number of subunits
	Size num_subunits(1);
	conformation::symmetry::SymmetryInfoCOP syminfo_ptr = NULL;
	if ( symmetric_pcs_calc_ && is_symmetric( pose ) ) {
		syminfo_ptr = symmetry_info( pose );
		num_subunits = syminfo_ptr->subunits();
	}

	Real max_pcs(0);
	for ( Size i = 1; i <= number_pcs_; ++i ) {
		Real norm_pcs = values[i] / scaling_factor_;
		Real norm_err = errors[i] / scaling_factor_;
		PCSSingle single_pcs(spinsA_all[i], pose, norm_pcs, norm_err);
		// If we perform PCS calculation with automatic deduction of protein symmetry,
		// make sure that only PCSs for the ASU are provided
		if ( symmetric_pcs_calc_ && is_symmetric( pose ) ) {
			for ( Size j = 1, j_end = single_pcs.get_protein_spins().size(); j <= j_end; ++j ) {
				// We test that the spins belong to subunit with index 1
				if ( syminfo_ptr->subunit_index(single_pcs.get_protein_spins()[j].rsd())  != 1 ) {
					utility_exit_with_message("ERROR in creation of PCSSingleSet for dataset " + dataset_name_
						+ ". For PCS calculation with automatic symmetry deduction the input datafile must contain only residue selections for the asymmetric subunit.");
				}
			}
		}
		pcs_single_vec_.push_back(single_pcs);
		pcs_values_(i) = norm_pcs;
		if ( std::abs(norm_pcs) > max_pcs ) {
			max_pcs = std::abs(norm_pcs);
		}

		// Resize spin_coordinates array
		spin_coordinates_[i].resize(num_subunits);
		for ( Size su = 1; su <= num_subunits; ++su ) {
			spin_coordinates_[i][su].resize(single_pcs.get_protein_spins().size());
		}
	}
	//svd_solver_->set_vector_b(pcs_values_); // lazy evaluation of SVD solver (see above)

	// Calculate weights of single pcs values
	for ( Size i = 1; i <= number_pcs_; ++i ) {
		if ( single_pcs_weighting_scheme_ == CONST ) {
			pcs_single_weights_(i) = 1.0;
			pcs_single_vec_[i].set_weight(1.0);
		} else {
			Real pcs_weight(1.0);
			Real pcs_err = pcs_single_vec_[i].get_pcs_err();
			runtime_assert_msg(pcs_err > 1e-6, "ERROR in calculating single PCS weights. Experimental PCS error is unreasonably small (< 1e-6) and would produce a large weight (> 10e+12) for the chosen weighting scheme. Check PCS datafile.");
			if ( single_pcs_weighting_scheme_ == SIGMA ) {
				pcs_weight = (1.0/(pcs_err * pcs_err));
			} else if ( single_pcs_weighting_scheme_ == OBSIG ) {
				pcs_weight = (std::abs(pcs_values_(i))/(pcs_err * pcs_err * max_pcs));
			}
			pcs_single_weights_(i) = pcs_weight;
			pcs_single_vec_[i].set_weight(pcs_weight);
		}
	}
}

/// @brief register options
void
PCSSingleSet::register_options() {
	using namespace basic::options;
	option.add_relevant(OptionKeys::nmr::pcs::nls_repeats);
	option.add_relevant(OptionKeys::nmr::pcs::use_symmetry_calc);
	option.add_relevant(OptionKeys::nmr::pcs::normalize_data);
}

void
PCSSingleSet::init_from_cml() {
	using namespace basic::options;
	symmetric_pcs_calc_ = option[ basic::options::OptionKeys::nmr::pcs::use_symmetry_calc ]();
	nls_repeats_ = option[ basic::options::OptionKeys::nmr::pcs::nls_repeats ]();
	normalized_data_ = option[ basic::options::OptionKeys::nmr::pcs::normalize_data ]();
}

/// @brief utility function to fill matrix_A that is used for SVD in PCSSingleSet
///        We average the PCS here over the number of equivalent spins (e.g. methyl protons).
void
PCSSingleSet::fill_matrix_A_row(
	utility::fixedsizearray1<Real,5 > & A_row,
	utility::vector1< Vector > const & spin_coord,
	Vector const & metal_coord,
	Real const & scal
)
{
	// Set all elements in row to zero because we use operator+=
	std::fill(A_row.begin(), A_row.end(), 0.0);
	Size num_eq_spins = spin_coord.size();

	for ( Size i = 1; i <= num_eq_spins; ++i ) {
		Real x(spin_coord[i].x() - metal_coord.x());
		Real y(spin_coord[i].y() - metal_coord.y());
		Real z(spin_coord[i].z() - metal_coord.z());
		Real x2(x * x);
		Real y2(y * y);
		Real z2(z * z);

		// PCS prefactor
		Real r2( x2 + y2 + z2);
		Real r5( r2 * r2 * std::sqrt(r2) );
		Real value_1_4_PI_r5(10000.0 / (4.0 * numeric::constants::d::pi * r5));

		A_row[1] += scal * value_1_4_PI_r5 * (x2 - z2);
		A_row[2] += scal * value_1_4_PI_r5 * (2.0 * x * y);
		A_row[3] += scal * value_1_4_PI_r5 * (2.0 * x * z);
		A_row[4] += scal * value_1_4_PI_r5 * (y2 - z2);
		A_row[5] += scal * value_1_4_PI_r5 * (2.0 * y * z);
	}

	if ( ave_type_ == MEAN ) {
		for ( Size j = 1; j <= 5; ++j ) {
			A_row[j] /= num_eq_spins;
		}
	}
}

/// @brief utility function that calculates one single PCS value given the input arguments
///        xM, yM, zM, Xax, Xrh, the spin coordinates, a rotation matrix and a scaling factor
Real
PCSSingleSet::basic_pcs_equation(
	utility::fixedsizearray1<Real,5> const & par,
	Vector const & spin_coord,
	Matrix const & rotM,
	Real const & scal
)
{
	// vector between spin and metal center
	Real x(spin_coord.x() - par[1]);
	Real y(spin_coord.y() - par[2]);
	Real z(spin_coord.z() - par[3]);

	// transformed vector after rotation
	Real x_t(rotM(1,1)*x + rotM(1,2)*y + rotM(1,3)*z);
	Real y_t(rotM(2,1)*x + rotM(2,2)*y + rotM(2,3)*z);
	Real z_t(rotM(3,1)*x + rotM(3,2)*y + rotM(3,3)*z);

	Real r2(x_t*x_t + y_t*y_t + z_t*z_t);
	Real r5(r2 * r2 * std::sqrt(r2));
	Real value_1_12_PI_r5(10000.0 / (12.0 * numeric::constants::d::pi * r5));

	Real pcs(scal * value_1_12_PI_r5 * (par[4] * (3.0 * z_t*z_t - r2) + par[5] * 1.5 * (x_t*x_t - y_t*y_t)));
	return pcs;
}

/// @brief utility functions to calculate the PCS from different sets of input arguments
///        arguments are a pointer to the parameters (xM, yM, zM, Xax, Xrh), the coordinates of equivalent spins,
///        a rotation matrix that must be previously constructed from the Euler angles and a scaling factor.
///        Fixed values of Xax and Xrh can be provided.
///        We average the PCS here over the number of equivalent spins (e.g. methyl protons).
Real
PCSSingleSet::fpcs(
	Real const *par,
	utility::vector1< Vector > const & spin_coord,
	Matrix const & rotM,
	Real const & scal
)
{
	utility::fixedsizearray1<Real,5> params = { par[0], par[1], par[2], par[3], par[4] };

	Size num_eq_spins = spin_coord.size();
	Real pcs(0);
	for ( Size i = 1; i <= num_eq_spins; ++i ) {
		pcs += basic_pcs_equation(params, spin_coord[i], rotM, scal);
	}
	if ( ave_type_ == MEAN ) {
		pcs /= num_eq_spins;
	}
	return pcs;
}

Real
PCSSingleSet::fpcs_Xax(
	Real const *par,
	Real const & Xax,
	utility::vector1< Vector > const & spin_coord,
	Matrix const & rotM,
	Real const & scal
)
{
	utility::fixedsizearray1<Real,5> params = { par[0], par[1], par[2], Xax, par[3] };

	Size num_eq_spins = spin_coord.size();
	Real pcs(0);
	for ( Size i = 1; i <= num_eq_spins; ++i ) {
		pcs += basic_pcs_equation(params, spin_coord[i], rotM, scal);
	}
	if ( ave_type_ == MEAN ) {
		pcs /= num_eq_spins;
	}
	return pcs;
}

Real
PCSSingleSet::fpcs_Xrh(
	Real const *par,
	Real const & Xrh,
	utility::vector1< Vector > const & spin_coord,
	Matrix const & rotM,
	Real const & scal
)
{
	utility::fixedsizearray1<Real,5> params = { par[0], par[1], par[2], par[3], Xrh };

	Size num_eq_spins = spin_coord.size();
	Real pcs(0);
	for ( Size i = 1; i <= num_eq_spins; ++i ) {
		pcs += basic_pcs_equation(params, spin_coord[i], rotM, scal);
	}
	if ( ave_type_ == MEAN ) {
		pcs /= num_eq_spins;
	}
	return pcs;
}

Real
PCSSingleSet::fpcs_Xax_Xrh(
	Real const *par,
	Real const & Xax,
	Real const & Xrh,
	utility::vector1< Vector > const & spin_coord,
	Matrix const & rotM,
	Real const & scal
)
{
	utility::fixedsizearray1<Real,5> params = { par[0], par[1], par[2], Xax, Xrh };

	Size num_eq_spins = spin_coord.size();
	Real pcs(0);
	for ( Size i = 1; i <= num_eq_spins; ++i ) {
		pcs += basic_pcs_equation(params, spin_coord[i], rotM, scal);
	}
	if ( ave_type_ == MEAN ) {
		pcs /= num_eq_spins;
	}
	return pcs;
}

/// @brief pcs error function used in the lmmin function
///        * par is an array of fit parameters [alpha, beta, gamma, xM, yM, zM, (Xax, Xrh)]
///        * data is a pointer to the the PCSSingleSet object i.e. to all data needed
///        for PCS calculation and NLS fitting
///        * fvc is an array holding the residuals of the fit calculation
void
pcs_erf(
	Real const *par,
	int m_dat,
	void const * data,
	Real *fvec,
	int */*info*/
)
{
	PCSSingleSet * pcs_singleset_ptr = static_cast< PCSSingleSet * >(const_cast< void* >(data));
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

	// Test if xM is in range [x_min, x_max] and constrain it to that range otherwise
	if ( !( pcs_singleset_ptr->metal_coord_bounds_[1] <= par[3] ) || !( par[3] <= pcs_singleset_ptr->metal_coord_bounds_[4] ) ) {
		nonconst_par[3] = (std::abs(pcs_singleset_ptr->metal_coord_bounds_[4] - pcs_singleset_ptr->metal_coord_bounds_[1])*0.5)*std::tanh(par[3])
			+ (pcs_singleset_ptr->metal_coord_bounds_[4] + pcs_singleset_ptr->metal_coord_bounds_[1])/2.0;
	}

	// Test if yM is in range [y_min, y_max] and constrain it to that range otherwise
	if ( !( pcs_singleset_ptr->metal_coord_bounds_[2] <= par[4] ) || !( par[4] <= pcs_singleset_ptr->metal_coord_bounds_[5] ) ) {
		nonconst_par[4] = (std::abs(pcs_singleset_ptr->metal_coord_bounds_[5] - pcs_singleset_ptr->metal_coord_bounds_[2])*0.5)*std::tanh(par[4])
			+ (pcs_singleset_ptr->metal_coord_bounds_[5] + pcs_singleset_ptr->metal_coord_bounds_[2])/2.0;
	}

	// Test if zM is in range [z_min, z_max] and constrain it to that range otherwise
	if ( !( pcs_singleset_ptr->metal_coord_bounds_[3] <= par[5] ) || !( par[5] <= pcs_singleset_ptr->metal_coord_bounds_[6] ) ) {
		nonconst_par[5] = (std::abs(pcs_singleset_ptr->metal_coord_bounds_[6] - pcs_singleset_ptr->metal_coord_bounds_[3])*0.5)*std::tanh(par[5])
			+ (pcs_singleset_ptr->metal_coord_bounds_[6] + pcs_singleset_ptr->metal_coord_bounds_[3])/2.0;
	}

	Vector euler_angles(nonconst_par[0], nonconst_par[1], nonconst_par[2]); // euler angles
	Matrix rotM = rotation_matrix_from_euler_angles(euler_angles, pcs_singleset_ptr->tensor_->get_euler_convention());
	Real scal(1.0 / pcs_singleset_ptr->scaling_factor_);

	for ( int i = 1; i <= m_dat; ++i ) {                // loop over the number of pcs
		fvec[i-1] = pcs_singleset_ptr->pcs_values_(i);
		for ( uint j = 1, j_end = pcs_singleset_ptr->spin_coordinates_[i].size(); j <= j_end; ++j ) { // loop over spins in symmetric subunits (if symmetric_pcs_calc_ is true)
			// otherwise middle vector in spin_coordinates_ has dimension 1
			// calculate average pcs for a group of degenerate spins (e.g. methyl protons
			// or symmetric spins if we don't use the automatic deduction of symmetric spins
			// but provide them explicitly as one selection in the input file) and calculate the residual
			if ( pcs_singleset_ptr->computation_type_ == PCSSingleSet::NLS ) {
				fvec[i-1] -= pcs_singleset_ptr->fpcs(&nonconst_par[3], pcs_singleset_ptr->spin_coordinates_[i][j], rotM, scal);
			} else if ( pcs_singleset_ptr->computation_type_ == PCSSingleSet::NLSAX ) {
				fvec[i-1] -= pcs_singleset_ptr->fpcs_Xax(&nonconst_par[3], pcs_singleset_ptr->tensor_->get_ax(), pcs_singleset_ptr->spin_coordinates_[i][j], rotM, scal);
			} else if ( pcs_singleset_ptr->computation_type_ == PCSSingleSet::NLSRH ) {
				fvec[i-1] -= pcs_singleset_ptr->fpcs_Xrh(&nonconst_par[3], pcs_singleset_ptr->tensor_->get_rh(), pcs_singleset_ptr->spin_coordinates_[i][j], rotM, scal);
			} else if ( pcs_singleset_ptr->computation_type_ == PCSSingleSet::NLSAXRH ) {
				fvec[i-1] -= pcs_singleset_ptr->fpcs_Xax_Xrh(&nonconst_par[3], pcs_singleset_ptr->tensor_->get_ax(), pcs_singleset_ptr->tensor_->get_rh(),
					pcs_singleset_ptr->spin_coordinates_[i][j], rotM, scal);
			}
		}
		//fvec[i-1] *= std::sqrt(pcs_singleset_ptr->pcs_single_weights_(i));
		fvec[i-1] *= pcs_singleset_ptr->pcs_single_weights_(i);
	}
}

/// @brief updates matrix_A using the metal coordinates it gets from
///        the grid search of the PCSMultiSet object
///        hands matrix_A over to the SVD solver too
void
PCSSingleSet::update_matrix_A(Vector const & metal_coord) {
	if ( !svd_solver_ ) {
		svd_solver_ = basic::svd::SVD_SolverOP( new basic::svd::SVD_Solver(number_pcs_, 5) );
		svd_solver_->set_vector_b(pcs_values_);
	}
	utility::fixedsizearray1<Real,5> A_row;
	Real scal(1.0/scaling_factor_);
	matrix_A_ = 0; // set matrix A to zero since we are using operator +=
	for ( Size i = 1; i <= number_pcs_; ++i ) {              // loop over number of pcs
		for ( Size su = 1, su_end = spin_coordinates_[i].size(); su <= su_end; ++su ) {  // loop over spins in symmetric subunits (if symmetric_pcs_calc_ is true)
			// otherwise middle vector in spin_coordinates_ has dimension 1
			fill_matrix_A_row(A_row, spin_coordinates_[i][su], metal_coord, scal);    // calculate average pcs for a group of degenerate spins (e.g. methyl protons)
			for ( Size k = 1; k <= 5; ++k ) {            // and sum it up
				matrix_A_(i,k) += A_row[k];
			}
		}
	}
	//runtime_assert_msg(svd_solver_, "ERROR while trying to set matrix A in PCSSingleSet. SVD_Solver has not been initialized yet.");
	svd_solver_->set_matrix_A(matrix_A_);
}

/// @brief updates the spin coordinates every time the pose is changed
///        make sure that this function is called before update_matrix_A() is called
void
PCSSingleSet::update_spin_coordinates(pose::Pose const & pose) {
	using pose::symmetry::is_symmetric;
	using pose::symmetry::symmetry_info;

	// Here we split behavior. If summetric_pcs_calc_ is true and the pose is indeed symmetric
	// we deduce the symmetric spins from the pose and fill their coordinates in the middle vector of spin_coordinates_
	Size num_subunits(1);
	conformation::symmetry::SymmetryInfoCOP syminfo_ptr = NULL;
	if ( symmetric_pcs_calc_ && is_symmetric( pose ) ) {
		syminfo_ptr = symmetry_info( pose );
		num_subunits = syminfo_ptr->subunits();
	}
	for ( Size i = 1; i <= number_pcs_; ++i ) {
		Size num_eq_spins(pcs_single_vec_[i].get_protein_spins().size());

		for ( Size j = 1; j <= num_eq_spins; ++j ) {
			Size spin_rsd_asu(pcs_single_vec_[i].get_protein_spins()[j].rsd());
			utility::vector1< Size > rsds_for_pcs_all_subunits( 1, spin_rsd_asu );

			if ( symmetric_pcs_calc_ && is_symmetric( pose ) && syminfo_ptr ) {
				utility::vector1< Size > symm_spin_rsds = syminfo_ptr->bb_clones(spin_rsd_asu);
				rsds_for_pcs_all_subunits.resize(num_subunits);
				std::copy(symm_spin_rsds.begin(), symm_spin_rsds.end(), rsds_for_pcs_all_subunits.begin()+1);
			}

			for ( Size k = 1; k <= num_subunits; ++k ) {
				spin_coordinates_[i][k][j] = pose.residue(rsds_for_pcs_all_subunits[k]).atom(pcs_single_vec_[i].get_protein_spins()[j].atomno()).xyz();

			} // number subunits

		} //number equivalent spins

	} // number pcs
}

void
PCSSingleSet::set_computation_type(std::string const & type) {
	convert_string_to_computation_type(type);
}

void
PCSSingleSet::set_averaging_type(std::string const & type) {
	ave_type_ = convert_string_to_averaging_type(type);
}

void
PCSSingleSet::show(std::ostream & tracer) const {
	tracer << "   * * * PCSSingleSet Summary Report * * *   " << std::endl;
	tracer << " PCS dataset " << dataset_name_ << " contains " << number_pcs_ << " pcs values and has weight " << weight_ << "." << std::endl;
	tracer << " PCS values: " << std::endl;
	tracer << " * * * * * * PCS values * * * * * * " << std::endl;
	for ( Size i = 1; i <= number_pcs_; ++i ) {
		pcs_single_vec_[i].show(tracer);
	}
	tracer << " * * * * * * PCS tensor * * * * * * " << std::endl;
	computation_type_ == SVD ? tensor_->show_tensor_stats(tracer, false) : tensor_->show_tensor_stats(tracer, true);
}

/// @brief solves the PCS tensor using SVD and returns the weighted PCS score
///        according to the single PCS weighting scheme
///        matrix A and spin coordinates must have been updated before calling
///        this function for the first time
Real
PCSSingleSet::solve_tensor_and_compute_score_by_svd(Vector const & metal_coord) {
	// update matrix A to be sure that it is set
	update_matrix_A(metal_coord);

	// decompose matrix A and get solution x for SLE Ax = b
	//runtime_assert_msg(svd_solver_, "ERROR while trying to decompose matrix A in PCSSingleSet. SVD_Solver has not been initialized yet.");
	svd_solver_->run_decomp_svd();
	svd_solver_->run_solve_svd();
	utility::vector1<Real>chiT(svd_solver_->get_svd_solution());
	chiT.push_back(metal_coord.x());
	chiT.push_back(metal_coord.y());
	chiT.push_back(metal_coord.z());

	// tensor elements are:
	// chiT[1] = chi_xx
	// chiT[2] = chi_xy
	// chiT[3] = chi_xz
	// chiT[4] = chi_yy
	// chiT[5] = chi_yz
	// chiT[6] = xM
	// chiT[7] = yM
	// chiT[8] = zM

	tensor_->set_tensor_in_arbitrary_frame(chiT);

	// Calculate PCS score taking into account the single pcs weights
	Real score(0);
	for ( Size i = 1; i <= number_pcs_; ++i ) {
		Real calc_pcs(0);

		for ( Size j = 1; j <= 5; ++j ) {
			calc_pcs += matrix_A_(i,j) * chiT[j];
		}
		pcs_single_vec_[i].set_pcs_calc(calc_pcs);
		Real diff = calc_pcs - pcs_values_(i);
		score += diff * diff * pcs_single_weights_(i);
	}

	// Show score and tensor after svd
	if ( TR.Trace.visible() ) {
		TR.Trace << "PCS score and tensor for dataset " << dataset_name_ << " after SVD: " << /*std::sqrt(score)*/ score << std::endl;
		tensor_->show_tensor_stats(TR.Trace, false);
	}

	//return std::sqrt(score);
	return score;
}

/// @brief solves the PCS tensor using NLS and returns the weighted PCS score
///        according to the single PCS weighting scheme
Real
PCSSingleSet::solve_tensor_and_compute_score_by_nls(Vector const & metal_coord) {

	// Set initial values for parameters
	Size number_params(6);
	utility::vector1<Real> params;
	utility::vector1<Real> params_best;

	for ( Size i = 1; i <= 3; ++i ) {
		params.push_back(360.0 * numeric::random::uniform());
	}
	params.push_back(metal_coord.x());
	params.push_back(metal_coord.y());
	params.push_back(metal_coord.z());

	if ( computation_type_ == PCSSingleSet::NLS ) {
		number_params = 8;
		params.push_back(-9999);
		params.push_back(-9999);
		params_best.resize(8);
	} else if ( (computation_type_ == PCSSingleSet::NLSAX) || (computation_type_ == PCSSingleSet::NLSRH) ) {
		number_params = 7;
		params.push_back(-9999);
		params_best.resize(7);
	} else if ( computation_type_ == PCSSingleSet::NLSAXRH ) {
		number_params = 6;
		params_best.resize(6);
	} else if ( computation_type_ == PCSSingleSet::SVD ) {
		utility_exit_with_message("ERROR in PCS tensor and score calculation. Cannot use NLS fitting because specified computation type is SVD.");
	}
	params_best=params; // Set best fit params to initial guess to avoid that they are not initialized in case lmmin fails

	// definition of auxiliary parameters for lmmin
	numeric::nls::lm_status_struct status;
	Real bestnorm = std::numeric_limits< Real >::infinity();

	// Perform nonlinear least squares fitting
	for ( Size i = 1; i <= nls_repeats_; ++i ) {
		// Random starting values
		TR.Trace << "PCS NLS fitting" << std::endl;
		TR.Trace << "Initialize angles to random values:";
		for ( Size j = 1; j <= 3; ++j ) {
			params[j] = 360.0 * numeric::random::uniform();
			TR.Trace << " " << params[j];
		}
		TR.Trace << std::endl;
		numeric::nls::lmmin( number_params, &params[1], number_pcs_, (const void*) this, pcs_erf, &status, numeric::nls::lm_printout_std);
		TR.Trace << "Iteration: " << i << " status.fnorm: " << status.fnorm << " bestnorm: "<< bestnorm << std::endl;
		//save to best fitting parameter
		if ( status.fnorm < bestnorm ) {
			TR.Trace << "status.fnorm: " << status.fnorm << " replaced bestnorm: " << bestnorm << std::endl;
			bestnorm=status.fnorm;
			params_best = params;
		}
	}

	// Set tensor with best parameters obtained from fit
	if ( computation_type_ == PCSSingleSet::NLSAX ) {
		params_best.resize(8);
		params_best[8] = tensor_->get_ax();
		std::swap(params_best[7], params_best[8]);
	} else if ( computation_type_ == PCSSingleSet::NLSRH ) {
		params_best.resize(8);
		params_best[8] = tensor_->get_rh();
	} else if ( computation_type_ == PCSSingleSet::NLSAXRH ) {
		params_best.resize(8);
		params_best[7] = tensor_->get_ax();
		params_best[8] = tensor_->get_rh();
	}
	tensor_->set_tensor_in_pas(params_best);

	// Calculate PCS score taking into account the single pcs weights
	Real score(0);

	Vector euler_angles(params_best[1], params_best[2], params_best[3]); // euler angles
	Matrix rotM = rotation_matrix_from_euler_angles(euler_angles, tensor_->get_euler_convention());
	Real scal(1.0 / scaling_factor_);

	for ( Size i = 1; i <= number_pcs_; ++i ) {            // loop over the number of pcs
		Real calc_pcs(0);
		for ( Size j = 1, j_end = spin_coordinates_[i].size(); j <= j_end; ++j ) {   // loop over spins in symmetric subunits (if symmetric_pcs_calc_ is true)
			// otherwise middle vector in spin_coordinates_ has dimension 1
			calc_pcs += fpcs(&params_best[4], spin_coordinates_[i][j], rotM, scal);   // calculate average pcs for a group of degenerate spins (e.g. methyl protons) and sum it up
		}
		pcs_single_vec_[i].set_pcs_calc(calc_pcs);
		Real diff = calc_pcs - pcs_values_(i);
		score += diff * diff * pcs_single_weights_(i);
	}

	// Show tensor after fit
	if ( TR.Trace.visible() ) {
		TR.Trace << "PCS score and tensor for dataset " << dataset_name_ << " after NLS: " << /*std::sqrt(score)*/ score << std::endl;
		tensor_->show_tensor_stats(TR.Trace, true);
	}

	//return std::sqrt(score);
	return score;
}

/// @brief sets the xyz derivative of the PCS
///        PCSTensor must be determined previously
///        call solve_tensor_and_compute_score_by_svd() or
///        solve_tensor_and_compute_score_by_nls() first
void
PCSSingleSet::set_atom_derivatives(pose::Pose & pose) {
	using pose::symmetry::is_symmetric;
	using pose::symmetry::symmetry_info;

	if ( tensor_->is_pcs_tensor_in_arbitrary_frame() ) {

		// Get the current tensor for the asymmetric unit
		Matrix chi_tensor_asu(Matrix::rows(
			tensor_->get_T_xx(), tensor_->get_T_xy(), tensor_->get_T_xz(),
			tensor_->get_T_xy(), tensor_->get_T_yy(), tensor_->get_T_yz(),
			tensor_->get_T_xz(), tensor_->get_T_yz(), -tensor_->get_T_xx()-tensor_->get_T_yy()));

		Vector metal_coord_asu(tensor_->get_metal_center());
		Real scal(1.0 / scaling_factor_);

		// Update matrix A with the current metal position
		update_matrix_A(metal_coord_asu);

		// Is pose symmetric?
		conformation::symmetry::SymmetryInfoCOP syminfo_ptr = NULL;
		Size num_subunits(1);
		if ( is_symmetric( pose ) ) {
			syminfo_ptr = symmetry_info( pose );
			debug_assert(syminfo_ptr);
			num_subunits = syminfo_ptr->subunits();
		}

		// Calculate and set xyz derivatives
		for ( Size i = 1; i <= number_pcs_; ++i ) {
			Real calc_pcs = matrix_A_(i,1) * chi_tensor_asu.xx()
				+ matrix_A_(i,2) * chi_tensor_asu.xy()
				+ matrix_A_(i,3) * chi_tensor_asu.xz()
				+ matrix_A_(i,4) * chi_tensor_asu.yy()
				+ matrix_A_(i,5) * chi_tensor_asu.yz();
			Real diff = calc_pcs - pcs_values_(i);

			// Here we iterate over the number of equivalent (degenerate) spins that produce the same pcs
			for ( Size j = 1, j_end = pcs_single_vec_[i].get_protein_spins().size(); j <= j_end; ++j ) {
				Real dx(0);
				Real dy(0);
				Real dz(0);

				// If symmetric_pcs_calc_ is true i.e. we deduce the symmetric spins automatically
				// we iterate over the number of symmetric subunits and transform the tensor
				// into the frame of the other subunits to calculate the derivative with
				// respect to all tensors in the whole structure that give rise to one pcs value
				Size residue_from = pcs_single_vec_[i].get_protein_spins()[j].rsd();
				utility::vector1< Size > residues_to(1, residue_from);

				if ( is_symmetric( pose ) && syminfo_ptr ) {
					residues_to.resize(num_subunits);
					if ( !syminfo_ptr->is_asymmetric_seqpos(residue_from) ) {
						residue_from = syminfo_ptr->get_asymmetric_seqpos(residue_from);
						residues_to[1] = residue_from;
					}
					utility::vector1< Size > symm_rsds = syminfo_ptr->bb_clones(residue_from);
					std::copy(symm_rsds.begin(), symm_rsds.end(), residues_to.begin()+1);
				}
				for ( Size su = 1; su <= num_subunits; ++su ) {
					Matrix chi_tensor_symm = apply_tensor_transformation( pose, chi_tensor_asu, residue_from, residues_to[su]);
					Vector metal_coord_symm = apply_vector_transformation( pose, metal_coord_asu, residue_from, residues_to[su]);

					// We calculate only the xyz derivatives for the asymmetric subunit.
					// Therefore, index of middle vector is 1.
					Real x(spin_coordinates_[i][1][j].x() - metal_coord_symm.x());
					Real y(spin_coordinates_[i][1][j].y() - metal_coord_symm.y());
					Real z(spin_coordinates_[i][1][j].z() - metal_coord_symm.z());
					Real x2(x * x);
					Real y2(y * y);
					Real z2(z * z);

					// PCS prefactor
					Real r2( x2 + y2 + z2);
					Real r3( r2 * std::sqrt(r2));
					Real r5( r2 * r3 );
					Real value_1_4_PI_r5(10000.0 / (4.0 * numeric::constants::d::pi * r5));
					Real value_4_PI_r3((4.0 * numeric::constants::d::pi * r3)/10000.0);
					Real calc_pcs_this_su = scal * value_1_4_PI_r5 * ( chi_tensor_symm.xx()*(x2-z2) + chi_tensor_symm.xy()*(2.0*x*y)
						+ chi_tensor_symm.xz()*(2.0*x*z) + chi_tensor_symm.yy()*(y2-z2) + chi_tensor_symm.yz()*(2.0*y*z) );

					// atom derivatives
					dx += ( (5.0 * calc_pcs_this_su * value_4_PI_r3 * x - scal*2.0*(chi_tensor_symm.xx()*x + chi_tensor_symm.xy()*y + chi_tensor_symm.xz()*z))
						* diff * pcs_single_weights_(i) * value_1_4_PI_r5);
					dy += ( (5.0 * calc_pcs_this_su * value_4_PI_r3 * y - scal*2.0*(chi_tensor_symm.xy()*x + chi_tensor_symm.yy()*y + chi_tensor_symm.yz()*z))
						* diff * pcs_single_weights_(i) * value_1_4_PI_r5);
					dz += ( (5.0 * calc_pcs_this_su * value_4_PI_r3 * z - scal*2.0*(chi_tensor_symm.xz()*x + chi_tensor_symm.yz()*y + (-chi_tensor_symm.xx()-chi_tensor_symm.yy())*z))
						* diff * pcs_single_weights_(i) * value_1_4_PI_r5);

				} // number subunits
				pcs_single_vec_[i].set_atom_derivatives(j, dx, dy, dz );
			} // equivalent spins

		} // all spins with assigned pcs value per asu

	} else if ( tensor_->is_pcs_tensor_in_pas() ) {

		// Get current tensor for the asymmetric unit
		Vector euler_angles_asu(tensor_->get_alpha(), tensor_->get_beta(), tensor_->get_gamma());
		Matrix rotM_asu = rotation_matrix_from_euler_angles(euler_angles_asu, tensor_->get_euler_convention());
		Vector metal_coord_asu(tensor_->get_metal_center());
		utility::fixedsizearray1< Real, 5 > chiT = { metal_coord_asu.x(), metal_coord_asu.y(), metal_coord_asu.z(), tensor_->get_ax(), tensor_->get_rh() };
		Real scal(1.0 / scaling_factor_);

		// Is pose symmetric
		conformation::symmetry::SymmetryInfoCOP syminfo_ptr = NULL;
		Size num_subunits(1);
		if ( is_symmetric( pose ) ) {
			syminfo_ptr = symmetry_info( pose );
			debug_assert(syminfo_ptr);
			num_subunits = syminfo_ptr->subunits();
		}

		// Calculate and set xyz derivatives
		for ( Size i = 1; i <= number_pcs_; ++i ) {
			Real calc_pcs(0);
			for ( Size j = 1, j_end = spin_coordinates_[i].size(); j <= j_end; ++j ) {
				calc_pcs += fpcs(&chiT[1], spin_coordinates_[i][j], rotM_asu, scal);
			}
			Real diff = calc_pcs - pcs_values_(i);

			// Here we iterate over the number of equivalent (degenerate) spins that produce the same pcs
			for ( Size j = 1, j_end = pcs_single_vec_[i].get_protein_spins().size(); j <= j_end; ++j ) {
				Real dx(0);
				Real dy(0);
				Real dz(0);

				// If symmetric_pcs_calc_ is true i.e. we deduce the symmetric spins automatically
				// we iterate over the number of symmetric subunits and transform the tensor
				// into the frame of the other subunits to calculate the derivative with
				// respect to all tensors in the whole structure that give rise to one pcs value
				Size residue_from = pcs_single_vec_[i].get_protein_spins()[j].rsd();
				utility::vector1< Size > residues_to(1, residue_from);

				if ( is_symmetric( pose ) && syminfo_ptr ) {
					residues_to.resize(num_subunits);
					if ( !syminfo_ptr->is_asymmetric_seqpos(residue_from) ) {
						residue_from = syminfo_ptr->get_asymmetric_seqpos(residue_from);
						residues_to[1] = residue_from;
					}
					utility::vector1< Size > symm_rsds = syminfo_ptr->bb_clones(residue_from);
					std::copy(symm_rsds.begin(), symm_rsds.end(), residues_to.begin()+1);
				}
				for ( core ::Size su = 1; su <= num_subunits; ++su ) {
					// Euler rotation matrix and metal coordinates in the new frame
					Matrix rotM_symm = apply_tensor_transformation( pose, rotM_asu, residue_from, residues_to[su]);
					Vector metal_coord_symm = apply_vector_transformation( pose, metal_coord_asu, residue_from, residues_to[su]);

					// We calculate only the xyz derivatives for the asymmetric subunit.
					// Therefore, index of middle vector is 1.
					Real x(spin_coordinates_[i][1][j].x() - metal_coord_symm.x());
					Real y(spin_coordinates_[i][1][j].y() - metal_coord_symm.y());
					Real z(spin_coordinates_[i][1][j].z() - metal_coord_symm.z());

					Real R11(rotM_symm(1,1));
					Real R12(rotM_symm(1,2));
					Real R13(rotM_symm(1,3));
					Real R21(rotM_symm(2,1));
					Real R22(rotM_symm(2,2));
					Real R23(rotM_symm(2,3));
					Real R31(rotM_symm(3,1));
					Real R32(rotM_symm(3,2));
					Real R33(rotM_symm(3,3));

					// transformed vector after rotation
					Real x_t(R11*x + R12*y + R13*z);
					Real y_t(R21*x + R22*y + R23*z);
					Real z_t(R31*x + R32*y + R33*z);

					Real r2(x_t*x_t + y_t*y_t + z_t*z_t);
					Real r3(r2 * std::sqrt(r2));
					Real r5(r3 * r2);
					Real r7(r5 * r2);
					Real value_1_12_PI(10000.0 / 12.0 * numeric::constants::d::pi);

					// precompute some terms of the derivative
					Real dAdx(-5.0 * (x*(R11*R11+R21*R21+R31*R31) + y*(R11*R12+R21*R22+R31*R32) + z*(R11*R13+R21*R23+R31*R33)) / r7);
					Real dAdy(-5.0 * (x*(R11*R12+R21*R22+R31*R32) + y*(R12*R12+R22*R22+R32*R32) + z*(R12*R13+R22*R23+R32*R33)) / r7);
					Real dAdz(-5.0 * (x*(R11*R13+R21*R23+R31*R33) + y*(R12*R13+R22*R23+R32*R33) + z*(R13*R13+R23*R23+R33*R33)) / r7);

					Real dBdx(2.0 * x * (chiT[4]*(2*R31*R31-R11*R11-R21*R21) + 1.5*chiT[5]*(R11*R11-R21*R21))
						+ 2.0 * y * (chiT[4]*(2*R31*R32-R11*R12-R21*R22) + 1.5*chiT[5]*(R11*R12-R21*R22))
						+ 2.0 * z * (chiT[4]*(2*R31*R33-R11*R13-R21*R23) + 1.5*chiT[5]*(R11*R13-R21*R23)));
					Real dBdy(2.0 * x * (chiT[4]*(2*R31*R32-R11*R12-R21*R22) + 1.5*chiT[5]*(R11*R12-R21*R22))
						+ 2.0 * y * (chiT[4]*(2*R32*R32-R12*R12-R22*R22) + 1.5*chiT[5]*(R12*R12-R22*R22))
						+ 2.0 * z * (chiT[4]*(2*R32*R33-R12*R13-R22*R23) + 1.5*chiT[5]*(R12*R13-R22*R23)));
					Real dBdz(2.0 * x * (chiT[4]*(2*R31*R33-R11*R13-R21*R23) + 1.5*chiT[5]*(R11*R13-R21*R23))
						+ 2.0 * y * (chiT[4]*(2*R32*R33-R12*R13-R22*R23) + 1.5*chiT[5]*(R12*R13-R22*R23))
						+ 2.0 * z * (chiT[4]*(2*R33*R33-R13*R13-R23*R23) + 1.5*chiT[5]*(R13*R13-R23*R23)));

					// atom derivatives
					dx += -value_1_12_PI * scal * (dAdx * (chiT[4]*(3.0*z_t*z_t-r2) + 1.5*chiT[5]*(x_t*x_t-y_t*y_t)) + 1.0/r5 * dBdx) * diff * pcs_single_weights_(i);
					dy += -value_1_12_PI * scal * (dAdy * (chiT[4]*(3.0*z_t*z_t-r2) + 1.5*chiT[5]*(x_t*x_t-y_t*y_t)) + 1.0/r5 * dBdy) * diff * pcs_single_weights_(i);
					dz += -value_1_12_PI * scal * (dAdz * (chiT[4]*(3.0*z_t*z_t-r2) + 1.5*chiT[5]*(x_t*x_t-y_t*y_t)) + 1.0/r5 * dBdz) * diff * pcs_single_weights_(i);

				} // number subunits
				pcs_single_vec_[i].set_atom_derivatives(j, dx, dy, dz );
			} // equivalent spins

		} // all spins with assigned pcs value per asu

	} else {
		utility_exit_with_message( "ERROR when setting PCS xyz derivative. PCSTensor is not set in arbitrary or principal axis frame. First call \"solve_tensor_and_compute_score_by_svd()\" or \"solve_tensor_and_compute_score_by_nls()\"." );
	}
}

/// @brief calculate PCS values from a given tensor and set values in the PCSSingle vector
Real
PCSSingleSet::compute_pcs_values_and_score_from_tensor(PCSTensor const & tensor) {
	if ( tensor.is_pcs_tensor_in_arbitrary_frame() ) {

		utility::fixedsizearray1<Real,5> chiT = { tensor.get_T_xx(), tensor.get_T_xy(), tensor.get_T_xz(), tensor.get_T_yy(), tensor.get_T_yz()};

		Real scal(1.0 / scaling_factor_);
		utility::fixedsizearray1<Real,5> Arow;
		Real score(0);

		for ( Size i = 1; i <= number_pcs_; ++i ) {
			Real calc_pcs(0);
			for ( Size su = 1, su_end = spin_coordinates_[i].size(); su <= su_end; ++su ) {
				fill_matrix_A_row(Arow, spin_coordinates_[i][su], tensor.get_metal_center(), scal);
				for ( Size k = 1; k <= 5; ++k ) {
					calc_pcs += Arow[k] * chiT[k];
				}
			}
			pcs_single_vec_[i].set_pcs_calc(calc_pcs);
			Real diff = calc_pcs - pcs_values_(i);
			score += diff * diff * pcs_single_weights_(i);
		}
		//return std::sqrt(score);
		return score;

	} else if ( tensor.is_pcs_tensor_in_pas() ) {

		Vector euler_angles(tensor.get_alpha(), tensor.get_beta(), tensor.get_gamma()); // euler angles
		Matrix rotM = rotation_matrix_from_euler_angles(euler_angles, tensor.get_euler_convention());
		utility::fixedsizearray1<Real,5> chiT = { tensor.get_metal_center().x(), tensor.get_metal_center().y(), tensor.get_metal_center().z(),
			tensor.get_ax(), tensor.get_rh() };

		Real scal(1.0 / scaling_factor_);
		Real score(0);

		for ( Size i = 1; i <= number_pcs_; ++i ) {
			Real calc_pcs(0);
			for ( Size su = 1, su_end = spin_coordinates_[i].size(); su <= su_end; ++su ) {
				calc_pcs += fpcs(&chiT[1], spin_coordinates_[i][su], rotM, scal);
			}
			pcs_single_vec_[i].set_pcs_calc(calc_pcs);
			Real diff = calc_pcs - pcs_values_(i);
			score += diff * diff * pcs_single_weights_(i);
		}
		//return std::sqrt(score);
		return score;

	} else {
		utility_exit_with_message( "ERROR when trying to calculate PCS values from given tensor. PCSTensor is not set in arbitrary or principal axis frame. First call \"solve_tensor_and_compute_score_by_svd()\" or \"solve_tensor_and_compute_score_by_nls()\"." );
	}
}

/// @brief calculate PCS values and the PCS score from the dataset's current tensor
Real
PCSSingleSet::compute_pcs_values_and_score_from_tensor() {
	if ( !tensor_ ) {
		utility_exit_with_message("ERROR  while trying to calculate PCSSingleSet score from current tensor. No PCSTensor object set.");
	}
	return compute_pcs_values_and_score_from_tensor( *tensor_ );
}

/// @brief utility functions to convert string to class specific enums
void
PCSSingleSet::convert_string_to_computation_type(std::string const & computation_type) {
	std::string str = boost::to_upper_copy<std::string>(computation_type);
	if ( str == "SVD" ) {
		computation_type_ = PCSSingleSet::SVD;
	} else if ( str == "NLS" ) {
		computation_type_ = PCSSingleSet::NLS;
	} else if ( str == "NLSAX" ) {
		computation_type_ = PCSSingleSet::NLSAX;
	} else if ( str == "NLSRH" ) {
		computation_type_ = PCSSingleSet::NLSRH;
	} else if ( str == "NLSAXRH" ) {
		computation_type_ = PCSSingleSet::NLSAXRH;
	} else {
		utility_exit_with_message( "ERROR: Provided string does not match any PCSSingleSet computation type. Possible options are \"SVD\", \"NLS\", \"NLSAX\", \"NLSRH\" and \"NLSAXRH\"." );
	}
}

} // namespace pcs
} // namespace nmr
} // namespace scoring
} // namespace core
