// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/pcs/PCSSingleSet.hh
/// @brief   class that stores and handles data of one single PCS dataset (i.e. for one lanthanide ion)
/// @details last Modified: 06/22/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

#ifndef INCLUDED_core_scoring_nmr_pcs_PCSSingleSet_HH
#define INCLUDED_core_scoring_nmr_pcs_PCSSingleSet_HH

// Unit headers
#include <core/scoring/nmr/pcs/PCSSingleSet.fwd.hh>

// Package headers
#include <core/scoring/nmr/pcs/PCSSingle.fwd.hh>
#include <core/scoring/nmr/pcs/PCSTensor.fwd.hh>
#include <core/scoring/nmr/util.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Basic headers
#include <basic/Tracer.fwd.hh>
#include <basic/svd/SVD_Solver.fwd.hh>

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
namespace pcs {

class PCSSingleSet {

public: // Enums

	/// @brief type of tensor and score calculation
	///        SVD     = singular value decomposition (coupled with grid search)
	///        NLS     = non-linear least squares fitting
	///        NLSAX   = non-linear least squares fitting with fixed axial tensor component
	///        NLSRH   = non-linear least squares fitting with fixed rhombic tensor component
	///        NLSAXRH = non-linear least squares fitting with fixed axial and rhombic tensor component
	enum COMPUTATION_TYPE {
		SVD     = 1,
		NLS     = 2,
		NLSAX   = 3,
		NLSRH   = 4,
		NLSAXRH = 5
	};

public: // Methods

	/// @brief construct from pcs datafile
	///        set default values for PCSSingleSet weight, single_pcs_weighting_scheme and computation type
	PCSSingleSet(
		std::string const & filename,
		std::string const & metal_ion_label,
		pose::Pose const & pose
	);

	/// @brief constructor with full argument list
	PCSSingleSet(
		std::string const & filename,
		std::string const & metal_ion_label,
		pose::Pose const & pose,
		Real const weight,
		std::string single_pcs_weigting = "CONST",
		std::string computation_type = "SVD"
	);

	/// @brief copy constructor
	PCSSingleSet(PCSSingleSet const & other);

	/// @brief assignment operator
	PCSSingleSet&
	operator=(PCSSingleSet const & rhs);

	/// @brief destructor
	~PCSSingleSet();

	/// @brief updates the spin coordinates every time the pose is changed
	///        make sure that this function is called before update_matrix_A() is called
	void
	update_spin_coordinates(pose::Pose const & pose);

	/// @brief updates matrix_A using the metal coordinates it gets from
	///        the grid search of the PCSMultiSet object
	///        hands matrix_A over to the SVD solver too
	void
	update_matrix_A(Vector const & metal_coord);

	/// @brief solves the PCS tensor using SVD and returns the weighted PCS score
	///        according to the single PCS weighting scheme
	Real
	solve_tensor_and_compute_score_by_svd(Vector const & metal_coord);

	/// @brief solves the PCS tensor using NLS and returns the weighted PCS score
	///        according to the single PCS weighting scheme
	Real
	solve_tensor_and_compute_score_by_nls(Vector const & metal_coord);

	/// @brief calculate PCS values from a given tensor, set values in the PCSSingle vector
	///        and return the PCS score
	Real
	compute_pcs_values_and_score_from_tensor(PCSTensor const & tensor);

	/// @brief calculate PCS values and the PCS score from the dataset's current tensor
	Real
	compute_pcs_values_and_score_from_tensor();

	/// @brief sets the xyz derivative of the PCS
	///        PCSTensor must be determined before
	///        first call solve_tensor_and_compute_score_by_svd() or
	///        solve_tensor_and_compute_score_by_nls() before setting derivatives
	void
	set_atom_derivatives(pose::Pose & pose);

	// Getters
	Real get_weight() const { return weight_; }
	Real get_scaling_factor() const { return scaling_factor_; }
	Size get_number_pcs() const { return number_pcs_; }
	std::string get_dataset_name() const { return dataset_name_; }
	std::string get_metal_ion_label() const { return metal_ion_label_; }
	utility::vector1<PCSSingle> const & get_single_pcs_vec() const { return pcs_single_vec_; }
	ObjexxFCL::FArray2D<Real> const & get_matrix_A() const { return matrix_A_; }
	ObjexxFCL::FArray1D<Real> const & get_pcs_values() const { return pcs_values_; }
	ObjexxFCL::FArray1D<Real> const & get_pcs_single_weights() const { return pcs_single_weights_; }
	PCSTensorCOP get_tensor_const() const { return tensor_; }
	PCSTensorOP get_tensor() { return tensor_; }
	COMPUTATION_TYPE get_computation_type() const { return computation_type_; }
	NMR_VALUE_AVERAGING_TYPE get_averaging_type() const { return ave_type_; }
	utility::vector1< utility::vector1< utility::vector1< Vector > > > const & get_spin_coordinates() const { return spin_coordinates_; }
	utility::fixedsizearray1<Real,6> const & get_metal_coord_bounds() const { return metal_coord_bounds_; }
	bool normalized_data() const { return normalized_data_; }
	bool symmetric_pcs_calc() const { return symmetric_pcs_calc_; }

	// Setters
	void set_weight(Real weight) { weight_ = weight; }
	void set_tensor(PCSTensorOP const & tensor) { tensor_ = tensor; }
	void set_metal_coord_bounds(utility::fixedsizearray1<Real,6> const & metal_coord_bounds) { metal_coord_bounds_ = metal_coord_bounds; }
	void set_computation_type(std::string const & type);
	void set_averaging_type(std::string const & type);

	void show(std::ostream & TR) const;

private: // Methods

	/// @brief default constructor
	PCSSingleSet();

	/// @brief utility function used in constructor to initialize PCSSingelSet
	///        object from pcs data file and pose.
	void
	init_from_pcs_filedata(
		std::string const & filename,
		pose::Pose const & pose
	);

	/// @brief register options
	void register_options();
	void init_from_cml();

	/// @brief utility function to fill matrix_A that is used for SVD in PCSSingleSet
	///        We average the PCS here over the number of equivalent spins (e.g. methyl protons).
	void
	fill_matrix_A_row(
		utility::fixedsizearray1<Real,5> & A_row,
		utility::vector1< Vector > const & spin_coords,
		Vector const & metal_coord,
		Real const & scal = 1.0
	);

	/// @brief utility functions to convert string to class specific enums
	void convert_string_to_computation_type(std::string const & computation_type);

	/// @brief utility function that calculates one single PCS value given the input arguments
	///        xM, yM, zM, Xax, Xrh, the spin coordinates, a rotation matrix and a scaling factor
	Real
	basic_pcs_equation(
		utility::fixedsizearray1<Real,5> const & par,
		Vector const & spin_coord,
		Matrix const & rotM,
		Real const & scal = 1.0
	);

	/// @brief utility functions to calculate the PCS from different sets of input arguments
	///        arguments are a pointer to the parameters (xM, yM, zM, Xax, Xrh), the coordinates of equivalent spins,
	///        a rotation matrix that must be previously constructed from the Euler angles and a scaling factor.
	///        Fixed values of Xax and Xrh can be provided.
	Real
	fpcs(
		Real const *par,
		utility::vector1< Vector > const & spin_coord,
		Matrix const & rotM,
		Real const & scal = 1.0
	);

	Real
	fpcs_Xax(
		Real const *par,
		Real const & Xax,
		utility::vector1< Vector > const & spin_coord,
		Matrix const & rotM,
		Real const & scal = 1.0
	);

	Real
	fpcs_Xrh(
		Real const *par,
		Real const & Xrh,
		utility::vector1< Vector > const & spin_coord,
		Matrix const & rotM,
		Real const & scal = 1.0
	);

	Real
	fpcs_Xax_Xrh(
		Real const *par,
		Real const & Xax,
		Real const & Xrh,
		utility::vector1< Vector > const & spin_coord,
		Matrix const & rotM,
		Real const & scal = 1.0
	);

	/// @brief pcs error function used in the lmmin function
	///        * par is an array of fit parameters [alpha, beta, gamma, xM, yM, zM, scal, (Xax, Xrh)]
	///        * data is a pointer to the the PCSSingleSet object i.e. to all data needed
	///        for PCS calculation and NLS fitting
	///        * fvc is an array holding the residuals of the fit calculation
	friend
	void
	pcs_erf(
		Real const *par,
		int m_dat,
		void const *data,
		Real *fvec,
		int */*info*/
	);

private: // Data

	std::string dataset_name_;
	std::string metal_ion_label_;
	utility::vector1<PCSSingle> pcs_single_vec_;
	Real weight_;
	Real scaling_factor_;
	Size number_pcs_;
	ObjexxFCL::FArray2D<Real> matrix_A_;
	ObjexxFCL::FArray1D<Real> pcs_values_;
	ObjexxFCL::FArray1D<Real> pcs_single_weights_;
	PCSTensorOP tensor_;
	basic::svd::SVD_SolverOP svd_solver_;
	SINGLE_NMR_VALUE_WEIGHTING single_pcs_weighting_scheme_;
	COMPUTATION_TYPE computation_type_;

	// The matrix A gets computed many times during the grid search. During that time
	// the pose stays unchanged. Thus, to prevent retrieving the spin coordinates every time
	// we test another metal ion position we store them here in a lookup table.
	// Also, for several repeats of NLS we don't want to retrieve the spin coordinates every time.
	// The outer vector runs over the number of pcs, the middle vector over symmetrical spins
	// (in case the pose is not symmetric this vector has only size 1) and the inner vector runs
	// over equivalent spins (e.g. methyl protons)
	utility::vector1< utility::vector1< utility::vector1< Vector > > > spin_coordinates_;

	// data the user may provide optionally
	// vector with metal coordinate bounds
	// must be set to a reasonable range when NLS fitting is used
	utility::fixedsizearray1<Real,6> metal_coord_bounds_;
	bool normalized_data_;
	bool symmetric_pcs_calc_;
	NMR_VALUE_AVERAGING_TYPE ave_type_;
	Size nls_repeats_;
};

} // namespace pcs
} // namespace nmr
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_scoring_nmr_pcs_PCSSingleSet_HH
