// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/rdc/RDCMultiSet.hh
/// @brief   class that stores and handles RDC data collected for one alignment medium
///          can contain multiple RDCSingleSets, each obtained for a different type of dipolar couplings (i.e. a different NMR experiment)
/// @details last Modified: 07/27/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

#ifndef INCLUDED_core_scoring_nmr_rdc_RDCMultiSet_HH
#define INCLUDED_core_scoring_nmr_rdc_RDCMultiSet_HH

// Unit headers
#include <core/scoring/nmr/rdc/RDCMultiSet.fwd.hh>

// Package headers
#include <core/io/nmr/AtomSelection.fwd.hh>
#include <core/scoring/nmr/rdc/RDCSingleSet.fwd.hh>
#include <core/scoring/nmr/rdc/RDCTensor.fwd.hh>
#include <core/scoring/nmr/util.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Basic headers
#include <basic/Tracer.fwd.hh>
#include <basic/svd/SVD_Solver.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/fixedsizearray1.hh>

// Numeric headers
#include <numeric/xyzVector.fwd.hh>
#include <numeric/xyzMatrix.fwd.hh>

// ObjexxFCL headers
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray1D.hh>

// C++ headers
#include <iostream>
#include <string>

namespace core {
namespace scoring {
namespace nmr {
namespace rdc {

class RDCMultiSet {
public:

	typedef std::pair< Vector, Vector > SpinPairCoordinates;

	/// @brief type of tensor and score calculation
	///        SVD     = singular value decomposition (coupled with grid search)
	///        NLS     = non-linear least squares fitting
	///        NLSDA   = non-linear least squares fitting with known alignment magnitude
	///        NLSR    = non-linear least squares fitting with known rhombicity
	///        NLSDAR  = non-linear least squares fitting with known alignment magnitude and rhombicity
	enum COMPUTATION_TYPE {
		SVD    = 1,
		NLS    = 2,
		NLSDA  = 3,
		NLSR   = 4,
		NLSDAR = 5
	};

public: // Methods

	/// @brief construct from data files
	///        set default values for RDCMultiSet weight and computation_type
	RDCMultiSet(
		utility::vector1<std::string> const & filenames,
		std::string const & medium_label,
		pose::Pose const & pose
	);

	/// @brief constructor with full argument list
	RDCMultiSet(
		utility::vector1<std::string> const & filenames,
		std::string const & medium_label,
		pose::Pose const & pose,
		Real const weight,
		std::string computation_type = "SVD"
	);

	/// @brief copy constructor
	RDCMultiSet(RDCMultiSet const & other);

	/// @brief assignment operator
	RDCMultiSet &
	operator=(RDCMultiSet const & rhs);

	/// @brief destructor
	~RDCMultiSet();

	/// @brief updates the spin coordinates every time the pose is changed
	///        make sure that this function is called before update_matrix_A() is called
	void
	update_spin_coordinates(pose::Pose const & pose);

	/// @brief builds matrix_A from the spin coordinates and hands it over to the SVD solver
	void
	update_matrix_A();

	/// @brief updates the array holding the single RDC weights that are used during score calculation
	///        needs to be called e.g. when the weighting scheme of one of the RDCSingleSets was changed
	void
	update_single_rdc_weighting();

	/// @brief solves the RDC tensor using SVD and returns the weighted RDC score
	///        according to the single RDC weighting scheme and the weights of the individual experiments
	Real
	solve_tensor_and_compute_score_by_svd();

	/// @brief solves the RDC tensor using NLS and returns the weighted RDC score
	///        according to the single RDC weighting scheme and the weights of the individual experiments
	Real
	solve_tensor_and_compute_score_by_nls();

	/// @brief calculate RDC values from a given tensor, set values in the RDCSingle vector
	///        and return the RDC score
	Real
	compute_rdc_values_and_score_from_tensor(RDCTensor const & tensor);

	/// @brief calculate RDC values and the RDC score from the dataset's current tensor
	Real
	compute_rdc_values_and_score_from_tensor();

	/// @brief sets the xyz derivative of the RDC
	///        RDCTensor must be determined before
	///        first call solve_tensor_and_compute_score_by_svd() or
	///        solve_tensor_and_compute_score_by_nls() before setting derivatives
	void
	set_atom_derivatives(pose::Pose const & pose);

	// Getters
	Real get_weight() const { return weight_; }
	Size get_total_number_rdc() const { return total_number_rdc_; }
	Size get_number_experiments() const { return number_experiments_; }
	std::string get_alignment_medium_label() const { return alignment_medium_; }
	utility::vector1< RDCSingleSetOP > & get_rdc_singleset_vec() { return rdc_singleset_vec_; }
	utility::vector1< RDCSingleSetOP > const & get_rdc_singleset_vec() const { return rdc_singleset_vec_; }
	ObjexxFCL::FArray2D<Real> const & get_matrix_A() const { return matrix_A_; }
	ObjexxFCL::FArray1D<Real> const & get_rdc_values() const { return rdc_values_; }
	ObjexxFCL::FArray1D<Real> const & get_rdc_single_weights() const { return rdc_single_weights_; }
	RDCTensorCOP get_tensor_const() const { return tensor_; }
	RDCTensorOP get_tensor() { return tensor_; }
	COMPUTATION_TYPE get_computation_type() const { return computation_type_; }
	NMR_VALUE_AVERAGING_TYPE get_averaging_type() const { return ave_type_; }
	utility::vector1< utility::vector1< RDCMultiSet::SpinPairCoordinates > > const & get_spin_coordinates() const { return spin_coordinates_; }
	RDC_NORM_TYPE get_normalization_type() const { return norm_type_; }
	bool symmetric_rdc_calc() const { return symmetric_rdc_calc_; }
	bool correct_sign() const { return correct_sign_; }
	bool tensor_fixed() const { return fixed_tensor_; }

	// Setters
	void set_weight(Real weight) { weight_ = weight; }
	void set_tensor(RDCTensorOP const & tensor) { tensor_ = tensor; }
	void set_computation_type(std::string const & type);
	void set_averaging_type(std::string const & type);
	void fix_tensor() { fixed_tensor_ = true; }
	void unfix_tensor() { fixed_tensor_ = false; }

	void show(std::ostream & TR) const;

private: // Methods

	/// @brief default constructor
	RDCMultiSet();

	/// @brief utility function used in constructor to initialize RDCMultiSet object from data files.
	void init_from_rdc_filedata(
		utility::vector1<std::string> const & filenames,
		pose::Pose const & pose
	);

	/// @brief register options
	void register_options();
	void init_from_cml();

	/// @brief utility function to convert string to class computation type enum
	void convert_string_to_computation_type(std::string const & computation_type);
	void deep_copy_rdc_single_set_vec(utility::vector1< RDCSingleSetOP > const & other_vec);

	/// @brief utility function to fill matrix_A that is used for SVD
	void
	fill_matrix_A_row(
		utility::fixedsizearray1<Real,5> & A_row,
		utility::vector1< SpinPairCoordinates > const & spin_coord,
		Real const & Dmax
	);

	/// @brief utility function that calculates one single RDC value given the input arguments
	///        Da, R, the spin coordinates, the maximal dipolar coupling constant and a rotation matrix
	Real
	basic_rdc_equation(
		utility::fixedsizearray1<Real,2> const & par,
		SpinPairCoordinates const & spin_coord,
		Real const & Dmax,
		Matrix const & rotM
	);

	/// @brief utility functions to calculate the RDC from different sets of input arguments
	///        arguments are a pointer to the parameters (Da, R), the coordinates of equivalent spins,
	///        the maximal dipolar coupling constant and a rotation matrix that must be previously
	///        constructed from the Euler angles. Fixed values of Da and R can be provided.
	Real
	frdc(
		Real const *par,
		utility::vector1< SpinPairCoordinates > const & spin_coord,
		Real const & Dmax,
		Matrix const & rotM
	);

	Real
	frdc_Da(
		Real const *par,
		Real const & Da,
		utility::vector1< SpinPairCoordinates > const & spin_coord,
		Real const & Dmax,
		Matrix const & rotM
	);

	Real
	frdc_R(
		Real const *par,
		Real const & R,
		utility::vector1< SpinPairCoordinates > const & spin_coord,
		Real const & Dmax,
		Matrix const & rotM
	);

	Real
	frdc_Da_R(
		Real const & Da,
		Real const & R,
		utility::vector1< SpinPairCoordinates > const & spin_coord,
		Real const & Dmax,
		Matrix const & rotM
	);

	/// @brief rdc error function used in the lmmin function
	///        * par is an array of fit parameters [alpha, beta, gamma, (Da, R)]
	///        * data is a pointer to the the RDCMultiSet object i.e. to all data needed
	///        for RDC calculation and NLS fitting
	///        * fvc is an array holding the residuals of the fit calculation
	friend
	void
	rdc_erf(
		Real const *par,
		int m_dat,
		void const *data,
		Real *fvec,
		int */*info*/
	);

private: // Data

	std::string alignment_medium_;
	utility::vector1< RDCSingleSetOP > rdc_singleset_vec_;
	Real weight_;
	Size total_number_rdc_;
	Size number_experiments_;

	// Since we are searching one single alignment tensor for all experiments
	// these arrays hold the RDC expression for all RDCSingleSets together.
	// In addition, the spin coordinates of all experiments are stored in one single matrix
	// and all experiments should share the same computation type and tensor (as well as fixed values for Da and R).
	ObjexxFCL::FArray2D<Real> matrix_A_;
	ObjexxFCL::FArray1D<Real> rdc_values_;
	ObjexxFCL::FArray1D<Real> rdc_single_weights_;
	basic::svd::SVD_SolverOP svd_solver_;
	RDCTensorOP tensor_;
	COMPUTATION_TYPE computation_type_;

	// When we run several repeats of NLS, we don't want to retrieve the spin coordinates every time.
	// Therefore, we look them up only once and store them as long as the pose hasn't changed.
	// The outer vector runs over the number of unique RDCs times the number of symmetric subunits,
	// and the inner vector runs over ambiguous spins (e.g. methyl protons or CA-HA1, CA-HA2 pairs for glycin)
	// over which the RDC is averaged
	utility::vector1< utility::vector1< SpinPairCoordinates > > spin_coordinates_;

	// data the user may provide optionally
	RDC_NORM_TYPE norm_type_;
	bool symmetric_rdc_calc_;
	bool correct_sign_; // correct sign of the 15N gyromagnetic ratio
	NMR_VALUE_AVERAGING_TYPE ave_type_;
	Size nls_repeats_;
	// Perform RDC calculation from fixed input tensor values
	// and keep tensor values fixed during the score calculation
	bool fixed_tensor_;

};

} // namespace rdc
} // namespace nmr
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_scoring_nmr_rdc_RDCMultiSet_HH
