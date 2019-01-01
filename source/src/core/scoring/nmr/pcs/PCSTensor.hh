// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/pcs/PCSTensor.hh
/// @brief   derived class of NMRTensor specific for the treatment of PCSs
///          the tensor we are using here is the magnetic susceptibility anisotropy tensor (delta chi tensor)
/// @details last Modified: 06/15/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

#ifndef INCLUDED_core_scoring_nmr_pcs_PCSTensor_HH
#define INCLUDED_core_scoring_nmr_pcs_PCSTensor_HH

// Unit headers
#include <core/scoring/nmr/NMRTensor.hh>
#include <core/scoring/nmr/pcs/PCSTensor.fwd.hh>
#include <core/scoring/nmr/util.hh>

// Package headers

// Project headers
#include <core/types.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.fwd.hh>

// Basic headers
#include <basic/Tracer.fwd.hh>

// Utility headers
#include <utility/vector1.hh>

// C++ headers
#include <string>
#include <iostream>

namespace core {
namespace scoring {
namespace nmr {
namespace pcs {

class PCSTensor : public NMRTensor {
public: // Methods

	/// @brief default constructor
	PCSTensor();

	/// @brief constructor with arguments
	///        initializes tensor with upper diagonal elements and the coordinates of the metal center
	///        parameters in vector should be in order chi__xx, chi_xy, chi_xz, chi_yy, chi_yz, xM, yM, zM
	PCSTensor(utility::vector1<Real> const & tensor_params);

	/// @brief copy constructor
	PCSTensor(PCSTensor const & other);

	/// @brief assignment operator
	PCSTensor & operator=(PCSTensor const & rhs);

	/// @brief virtual destructor
	~PCSTensor() override;

	/// @brief every tensor should offer a function that returns its name
	std::string tensor_name() const override;

	/// @brief sets the delta chi tensor as it exists in the arbitrary frame
	///        resets the tensor upper diagonal elements and the coordinates of the metal center
	///        parameters in vector should be in order chi__xx, chi_xy, chi_xz, chi_yy, chi_yz, xM, yM, zM
	void set_tensor_in_arbitrary_frame(utility::vector1<Real> const & tensor_params) override;

	/// @brief sets the delta chi tensor as it exists in the principal axis system (PAS)
	///        resets the tensor principal values, the axial and rhombic component,
	///        the three Euler angles and the coordinates of the metal center
	///        parameters in vector should be in order alpha, beta, gamma, xM, yM, zM, ax, rh
	void set_tensor_in_pas(utility::vector1<Real> const & tensor_params) override;

	/// @brief brings tensor in principal axis frame by diagonalization
	///        make sure to call either set_tensor_in_arbitrary_frame() or set_tensor_in_pas() first
	void diagonalize_tensor() override;

	/// @brief brings delta chi tensor principal values in correct order
	///        and reconfigures tensor into unique tensor representation
	///        make sure to call set_tensor_in_pas() or diagonalize_tensor() first
	void reorder_tensor() override;

	/// @brief shows summary of tensor statistics
	/// @details If show_in_pas is true and tensor has been set or transformed in the
	///          principal axis system, show tensor parameters in PAS. If either of
	///          the two conditions is not fulfilled show tensor matrix in arbitrary frame.
	void show_tensor_stats(std::ostream & TR, bool show_in_pas=true) const override;

	/// @brief serialize a PCSTensor to a json_spirit object
	utility::json_spirit::Value serialize() const override;

	/// @brief deserialize a json_spirit object to a PCSTensor
	void deserialize(utility::json_spirit::mObject data) override;

	// Accessors
	bool is_pcs_tensor_in_arbitrary_frame() const { return pcs_tensor_in_arbitrary_frame_; }
	bool is_pcs_tensor_in_pas() const { return pcs_tensor_in_pas_; }
	bool is_pcs_tensor_diagonalized() const { return pcs_tensor_diagonalized_; }
	bool is_pcs_tensor_reconfigured() const { return pcs_tensor_reconfigured_; }
	Vector const & get_metal_center() const { return metal_center_; }
	void set_metal_center(Real xM, Real yM, Real zM) { metal_center_.x(xM); metal_center_.y(yM); metal_center_.z(zM); }


private: // Methods

	/// @brief converts tensor principal values into its axial and rhombic component
	void ax_rh_from_eigenvalues();

	/// @brief converts axial and rhombic component into the tensor principal values
	void ax_rh_to_eigenvalues();

	/// @brief returns the delta chi tensor in 3x3 matrix form
	Matrix build_tensor_matrix();

private: // Data

	// boolean variables to store the state of the tensor
	bool pcs_tensor_in_arbitrary_frame_;
	bool pcs_tensor_in_pas_;
	bool pcs_tensor_diagonalized_;
	bool pcs_tensor_reconfigured_;

	// metal center
	Vector metal_center_;

};

} // namespace nmr
} // namespace nmr
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_scoring_nmr_pcs_PCSTensor_HH
