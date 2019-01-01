// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/rdc/RDCTensor.hh
/// @brief   derived class of NMRTensor specific for the treatment of RDCs
///          the tensor we are using here is the alignment tensor A that
///          is related to the Saupe order matrix according to S = 3/2 A
/// @details last Modified: 06/12/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

#ifndef INCLUDED_core_scoring_nmr_rdc_RDCTensor_HH
#define INCLUDED_core_scoring_nmr_rdc_RDCTensor_HH

// Unit headers
#include <core/scoring/nmr/NMRTensor.hh>
#include <core/scoring/nmr/rdc/RDCTensor.fwd.hh>
#include <core/scoring/nmr/util.hh>

// Project headers
#include <core/types.hh>

// Basic headers
#include <basic/Tracer.fwd.hh>

// Utility headers
#include <utility/vector1.hh>

// Numeric headers
#include <numeric/xyzVector.fwd.hh>
#include <numeric/xyzMatrix.fwd.hh>

// C++ headers
#include <string>
#include <iostream>

namespace core {
namespace scoring {
namespace nmr {
namespace rdc {

class RDCTensor : public NMRTensor {

public: // Methods

	/// @brief default constructor
	RDCTensor();

	/// @brief constructor with arguments
	///        initializes alignment tensor with upper diagonal elements
	///        vector should contain tensor parameters in the order A_xx, A_xy, A_xz, A_yy, A_yz
	RDCTensor(utility::vector1<Real> const & tensor_params);

	/// @brief copy constructor
	RDCTensor(RDCTensor const & other);

	/// @brief assignment operator
	RDCTensor &
	operator=(RDCTensor const & rhs);

	/// @brief virtual destructor
	~RDCTensor() override;

	/// @brief every tensor should offer a function that returns its name
	std::string tensor_name() const override;

	/// @brief sets the alignment tensor as it exists in the arbitrary frame
	///        resets the tensor upper diagonal elements
	///        vector should contain tensor parameters in the order A_xx, A_xy, A_xz, A_yy, A_yz
	void set_tensor_in_arbitrary_frame(utility::vector1<Real> const & tensor_params) override;

	/// @brief sets the alignment tensor as it exists in the principal axis system (PAS)
	///        resets the alignment tensor principal values, the alignment order, the rhombicity,
	///        the three Euler angles and the dipolar coupling constant
	///        vector should contain tensor parameters in the order alpha, beta, gamma, Da, R, Dmax
	void set_tensor_in_pas(utility::vector1<Real> const & tensor_params) override;

	/// @brief brings alignment tensor in principal axis frame by diagonalization
	///        make sure to call either set_tensor_in_arbitrary_frame() and set_Dmax() or
	///        set_tensor_in_pas() first
	void diagonalize_tensor() override;

	/// @brief brings alignment tensor principal values in correct order
	///        and reconfigure tensor into unique tensor representation
	///        make sure to call set_tensor_in_pas() or diagonalize_tensor() first
	void reorder_tensor() override;

	/// @brief shows summary of tensor statistics
	/// @details If show_in_pas is true and tensor has been set or transformed in the
	///          principal axis system, show tensor parameters in PAS. If either of
	///          the two conditions is not fulfilled show tensor matrix in arbitrary frame.
	void show_tensor_stats(std::ostream & TR, bool show_in_pas=true) const override;

	/// @brief serialize an RDCTensor to a json_spirit object
	utility::json_spirit::Value serialize() const override;

	/// @brief deserialize a json_spirit object to an RDCTensor
	void deserialize(utility::json_spirit::mObject data) override;

	// Accessors
	Real get_Da() const { return Da_; }
	Real get_R() const { return R_; }
	Real get_Dmax() const { return Dmax_; }
	bool is_rdc_tensor_in_arbitrary_frame() const { return rdc_tensor_in_arbitrary_frame_; }
	bool is_rdc_tensor_in_pas() const { return rdc_tensor_in_pas_; }
	bool is_tensor_diagonalized() const { return rdc_tensor_diagonalized_; }
	bool is_tensor_reconfigured() const { return rdc_tensor_reconfigured_; }

	// Setters
	void set_Dmax(Real Dmax) { Dmax_ = Dmax; dmax_set_ = true; }

private: //Methods

	/// @brief converts alignment tensor principal values into axial and rhombic component, alignment order and rhombicity
	void Da_R_ax_rh_from_eigenvalues();

	/// @brief calculates alignment tensor principal values, its axial and rhombic component from the alignment order and rhombicity
	void Da_R_to_ax_rh_and_eigenvalues();

	/// @brief returns the alignment tensor in 3x3 matrix form
	Matrix build_tensor_matrix();

private: // Data

	// Alignment order and rhombicity
	Real Da_;
	Real R_;

	// Maximal dipolar coupling constant
	Real Dmax_;

	// boolean variables to store the state of the tensor
	bool rdc_tensor_in_arbitrary_frame_;
	bool rdc_tensor_in_pas_;
	bool dmax_set_;
	bool rdc_tensor_diagonalized_;
	bool rdc_tensor_reconfigured_;

};

} // namespace rdc
} // namespace nmr
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_scoring_nmr_rdc_RDCTensor_HH
