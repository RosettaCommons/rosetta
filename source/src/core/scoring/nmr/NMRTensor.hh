// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/NMRTensor.hh
/// @brief   base tensor class for the representation of different NMR interactions like PCS, RDC or CSA
/// @details last Modified: 06/08/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

#ifndef INCLUDED_core_scoring_nmr_NMRTensor_HH
#define INCLUDED_core_scoring_nmr_NMRTensor_HH

// Unit headers
#include <core/scoring/nmr/NMRTensor.fwd.hh>
#include <core/scoring/nmr/util.hh>

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/json_spirit/json_spirit_value.h>

// Numeric headers
#include <numeric/angle.functions.hh>

// C++ headers
#include <string>
#include <iostream>

namespace core {
namespace scoring {
namespace nmr {

class NMRTensor : public utility::pointer::ReferenceCount {
public: // Methods

	/// @brief default constructor
	NMRTensor();

	/// @brief copy constructor
	NMRTensor(NMRTensor const & other);

	/// @brief assignment operator
	NMRTensor & operator=(NMRTensor const & rhs);

	/// @brief virtual destructor
	virtual ~NMRTensor();

	/// @brief every tensor should offer a function that returns its name
	virtual std::string tensor_name() const = 0;

	/// @brief every NMR tensor should offer a function to set and update its upper diagonal elements
	///        in this way the NMR tensor is represented in the arbitrary (molecular) frame
	virtual void set_tensor_in_arbitrary_frame(utility::vector1<Real> const & tensor_params) = 0;

	/// @brief every NMR tensor should offer a function to set and update its axial and rhombic component and three Euler angles
	///        in this way the NMR tensor is represented in its principal axis system (PAS)
	virtual void set_tensor_in_pas(utility::vector1<Real> const & tensor_params) = 0;

	/// @brief bring NMR tensor in principal axis frame by diagonalization
	virtual void diagonalize_tensor();

	/// @brief Bring the NMR tensor principal values in correct order
	///        and reconfigure tensor into unique tensor representation
	virtual void reorder_tensor();

	/// @brief show summary of tensor statistics
	virtual void show_tensor_stats(std::ostream &, bool ) const;

	/// @brief serialize an NMRTensor to a json_spirit object
	virtual utility::json_spirit::Value serialize() const;

	/// @brief deserialize a json_spirit object to an NMRTensor
	virtual void deserialize(utility::json_spirit::mObject data);

	// Accessors
	Real get_T_xx() const { return T_xx_; }
	Real get_T_xy() const { return T_xy_; }
	Real get_T_xz() const { return T_xz_; }
	Real get_T_yy() const { return T_yy_; }
	Real get_T_yz() const { return T_yz_; }
	Real get_ax() const { return ax_; }
	Real get_rh() const { return rh_; }
	Real get_Eig_xx() const { return Eig_xx_; }
	Real get_Eig_yy() const { return Eig_yy_; }
	Real get_Eig_zz() const { return Eig_zz_; }
	Real get_alpha() const { return alpha_; }
	Real get_beta() const { return beta_; }
	Real get_gamma() const { return gamma_; }
	EULER_CONVENTION get_euler_convention() const { return convention_; }

	// Mutators
	void set_T_xx( Real T_xx ) { T_xx_ = T_xx; }
	void set_T_xy( Real T_xy ) { T_xy_ = T_xy; }
	void set_T_xz( Real T_xz ) { T_xz_ = T_xz; }
	void set_T_yy( Real T_yy ) { T_yy_ = T_yy; }
	void set_T_yz( Real T_yz ) { T_yz_ = T_yz; }
	void set_ax( Real ax ) { ax_ = ax; }
	void set_rh( Real rh ) { rh_ = rh; }
	void set_Eig_xx( Real Eig_xx ) { Eig_xx_ = Eig_xx; }
	void set_Eig_yy( Real Eig_yy ) { Eig_yy_ = Eig_yy; }
	void set_Eig_zz( Real Eig_zz ) { Eig_zz_ = Eig_zz; }
	void set_alpha( Real alpha ) { alpha_ = numeric::nonnegative_principal_angle_degrees(alpha); }
	void set_beta( Real beta ) { beta_ = numeric::nonnegative_principal_angle_degrees(beta); }
	void set_gamma( Real gamma) { gamma_ = numeric::nonnegative_principal_angle_degrees(gamma); }

private: // Data

	// Euler angles and convention
	Real alpha_;
	Real beta_;
	Real gamma_;
	EULER_CONVENTION convention_;

	// axial and rhombic component
	Real ax_;
	Real rh_;

	// 5 unique tensor values
	Real T_xx_;
	Real T_xy_;
	Real T_xz_;
	Real T_yy_;
	Real T_yz_;

	// 3 principal components of the tensor in the principal axis system (PAS)
	Real Eig_xx_;
	Real Eig_yy_;
	Real Eig_zz_;

};

} // namespace nmr
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_scoring_nmr_NMRTensor_HH
