// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/PCSTensor.cc
/// @brief   Implementation of class PCSTensor
/// @details last Modified: 06/15/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

// Unit headers
#include <core/scoring/nmr/pcs/PCSTensor.hh>
#include <core/scoring/nmr/util.hh>

// Package headers

// Project headers
#include <core/types.hh>

// Basic headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/tools/make_vector.hh>

// Numeric headers
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.json.hh>

// ObjexxFCL headers
#include <ObjexxFCL/format.hh>

// C++ headers
#include <string>
#include <iostream>

namespace core {
namespace scoring {
namespace nmr {
namespace pcs {

static basic::Tracer TR( "core.scoring.nmr.pcs.PCSTensor" );

/// @brief default constructor
PCSTensor::PCSTensor() :
	NMRTensor(),
	pcs_tensor_in_arbitrary_frame_(false),
	pcs_tensor_in_pas_(false),
	pcs_tensor_diagonalized_(false),
	pcs_tensor_reconfigured_(false),
	metal_center_(0.0)
{}

/// @brief constructor with arguments
///        initializes tensor with upper diagonal elements and the coordinates of the metal center
///        parameters in vector should be in order chi__xx, chi_xy, chi_xz, chi_yy, chi_yz, xM, yM, zM
PCSTensor::PCSTensor(utility::vector1<Real> const & tensor_params) :
	NMRTensor(),
	pcs_tensor_in_arbitrary_frame_(false),
	pcs_tensor_in_pas_(false),
	pcs_tensor_diagonalized_(false),
	pcs_tensor_reconfigured_(false)
{
	set_tensor_in_arbitrary_frame(tensor_params);
}

/// @brief copy constructor
PCSTensor::PCSTensor(PCSTensor const & other) :
	NMRTensor(other),
	pcs_tensor_in_arbitrary_frame_(other.pcs_tensor_in_arbitrary_frame_),
	pcs_tensor_in_pas_(other.pcs_tensor_in_pas_),
	pcs_tensor_diagonalized_(other.pcs_tensor_diagonalized_),
	pcs_tensor_reconfigured_(other.pcs_tensor_reconfigured_),
	metal_center_(other.metal_center_)
{}

/// @brief assignment operator
PCSTensor &
PCSTensor::operator=(PCSTensor const & rhs) {
	if ( this != &rhs ) {
		NMRTensor::operator =(rhs);
		pcs_tensor_in_arbitrary_frame_ = rhs.pcs_tensor_in_arbitrary_frame_;
		pcs_tensor_in_pas_ = rhs.pcs_tensor_in_pas_;
		pcs_tensor_diagonalized_ = rhs.pcs_tensor_diagonalized_;
		pcs_tensor_reconfigured_ = rhs.pcs_tensor_reconfigured_;
		metal_center_ = rhs.metal_center_;
	}
	return *this;
}

/// @brief virtual destructor
PCSTensor::~PCSTensor() {}

/// @brief every tensor should offer a function that returns its name
std::string
PCSTensor::tensor_name() const {
	return "PCSTensor";
}

/// @brief sets the delta chi tensor as it exists in the arbitrary frame
///        resets the tensor upper diagonal elements and the coordinates of the metal center
///        parameters in vector should be in order chi__xx, chi_xy, chi_xz, chi_yy, chi_yz, xM, yM, zM
void
PCSTensor::set_tensor_in_arbitrary_frame(utility::vector1<Real> const & tensor_params) {
	runtime_assert_msg(tensor_params.size() >= 8, "ERROR: Vector contains not enough parameter to instantiate PCSTensor object.");
	set_T_xx(tensor_params[1]);
	set_T_xy(tensor_params[2]);
	set_T_xz(tensor_params[3]);
	set_T_yy(tensor_params[4]);
	set_T_yz(tensor_params[5]);
	metal_center_ = Vector(tensor_params[6], tensor_params[7], tensor_params[8]);
	pcs_tensor_in_arbitrary_frame_ = true;
}

/// @brief sets the delta chi tensor as it exists in the principal axis system (PAS)
///        resets the tensor principal values, the axial and rhombic component,
///        the three Euler angles and the coordinates of the metal center
///        parameters in vector should be in order alpha, beta, gamma, xM, yM, zM, ax, rh
void
PCSTensor::set_tensor_in_pas(utility::vector1<Real> const & tensor_params) {
	runtime_assert_msg(tensor_params.size() >= 8, "ERROR: Vector contains not enough parameter to instantiate PCSTensor object.");
	set_ax(tensor_params[7]);
	set_rh(tensor_params[8]);
	set_alpha(tensor_params[1]);
	set_beta(tensor_params[2]);
	set_gamma(tensor_params[3]);
	metal_center_ = Vector(tensor_params[4], tensor_params[5], tensor_params[6]);
	ax_rh_to_eigenvalues();
	pcs_tensor_in_pas_ = true;
	pcs_tensor_diagonalized_ = true;
}

/// @brief brings tensor in principal axis frame by diagonalization
///        make sure to call either set_tensor_in_arbitrary_frame() or set_tensor_in_pas() first
void
PCSTensor::diagonalize_tensor() {
	if ( !pcs_tensor_in_arbitrary_frame_ && !pcs_tensor_in_pas_ ) {
		utility_exit_with_message("ERROR: PCSTensor is not set. Call \"set_tensor_in_arbitray_frame()\" or \"set_tensor_in_pas()\" before diagonalization.");
	}
	if ( pcs_tensor_in_arbitrary_frame_ ) {
		Matrix chi_ = build_tensor_matrix();
		Real tol(0.0);
		Matrix evecs;
		Vector evals = numeric::eigenvector_jacobi(chi_, tol, evecs);
		evecs.transpose();
		Vector angles = euler_angles_from_rotation_matrix(evecs, get_euler_convention());
		set_Eig_xx(evals(1));
		set_Eig_yy(evals(2));
		set_Eig_zz(evals(3));
		ax_rh_from_eigenvalues();
		pcs_tensor_diagonalized_ = true;
		set_alpha(angles(1));
		set_beta(angles(2));
		set_gamma(angles(3));
		pcs_tensor_in_pas_ = true;
	} else {
		return;
	}
}

/// @brief brings delta chi tensor principal values in correct order
///        and reconfigures tensor into unique tensor representation
///        make sure to call set_tensor_in_pas() or diagonalize_tensor() first
void
PCSTensor::reorder_tensor() {
	if ( !pcs_tensor_reconfigured_ ) {

		if ( !pcs_tensor_diagonalized_ ) {
			TR.Warning << "PCSTensor is not diagonalized. Call \"diagonalize_tensor()\" first" << std::endl;
			diagonalize_tensor();
		}

		Real Eig_xx = get_Eig_xx();
		Real Eig_yy = get_Eig_yy();
		Real Eig_zz = get_Eig_zz();
		Real alpha  = get_alpha();
		Real beta   = get_beta();
		Real gamma  = get_gamma();
		utility::vector1<Real> tensor_params;

		tensor_params.push_back(Eig_xx);
		tensor_params.push_back(Eig_yy);
		tensor_params.push_back(Eig_zz);
		tensor_params.push_back(alpha);
		tensor_params.push_back(beta);
		tensor_params.push_back(gamma);

		// Call to utility function that reconfigures tensor
		order_tensor_parameters(tensor_params, get_euler_convention());

		set_Eig_xx(tensor_params[1]);
		set_Eig_yy(tensor_params[2]);
		set_Eig_zz(tensor_params[3]);
		set_alpha(tensor_params[4]);
		set_beta(tensor_params[5]);
		set_gamma(tensor_params[6]);

		ax_rh_from_eigenvalues();

		pcs_tensor_reconfigured_ = true;
	} else {
		return;
	}
}

/// @brief shows summary of tensor statistics
/// @details If show_in_pas is true and tensor has been set or transformed in the
///          principal axis system, show tensor parameters in PAS. If either of
///          the two conditions is not fulfilled show tensor matrix in arbitrary frame.
void
PCSTensor::show_tensor_stats(std::ostream & tracer, bool show_in_pas) const {
	using namespace ObjexxFCL;
	using namespace ObjexxFCL::format;
	Size width(10);
	Size precision(3);
	tracer << "   * * * PCSTensor parameters * * *   "  << std::endl;
	if ( show_in_pas && ( pcs_tensor_in_pas_ || (pcs_tensor_diagonalized_ && pcs_tensor_reconfigured_ ) ) ) {
		tracer << "Xax    (10-32 m3) = " << F(width, precision, get_ax()) << std::endl;
		tracer << "Xrh    (10-32 m3) = " << F(width, precision, get_rh()) << std::endl;
		tracer << "alpha  (degrees)  = " << F(width, precision, get_alpha()) << std::endl;
		tracer << "beta   (degrees)  = " << F(width, precision, get_beta()) << std::endl;
		tracer << "gamma  (degrees)  = " << F(width, precision, get_gamma()) << std::endl;
		tracer << "xM     (Ang.)     = " << F(width, precision, get_metal_center().x()) << std::endl;
		tracer << "yM     (Ang.)     = " << F(width, precision, get_metal_center().y()) << std::endl;
		tracer << "zM     (Ang.)     = " << F(width, precision, get_metal_center().z()) << std::endl;
	} else {
		tracer << "[[ " << F(width, precision, get_T_xx()) << ", "
			<< F(width, precision, get_T_xy()) << ", "
			<< F(width, precision, get_T_xz()) << " ]," << std::endl;
		tracer << " [ " << F(width, precision, get_T_xy()) << ", "
			<< F(width, precision, get_T_yy()) << ", "
			<< F(width, precision, get_T_yz()) << " ]," << std::endl;
		tracer << " [ " << F(width, precision, get_T_xz()) << ", "
			<< F(width, precision, get_T_yz()) << ", "
			<< F(width, precision, -get_T_xx()-get_T_yy()) << " ]]" << std::endl;
		tracer << "xM     = " << F(width, precision, get_metal_center().x()) << std::endl;
		tracer << "yM     = " << F(width, precision, get_metal_center().y()) << std::endl;
		tracer << "zM     = " << F(width, precision, get_metal_center().z()) << std::endl;
	}
}

/// @brief serialize a PCSTensor to a json_spirit object
utility::json_spirit::Value
PCSTensor::serialize() const
{
	using utility::json_spirit::Value;
	using utility::json_spirit::Pair;

	Pair tensor_in_arbitrary("pcs_tensor_in_arbitrary_frame", Value(pcs_tensor_in_arbitrary_frame_));
	Pair tensor_in_pas("pcs_tensor_in_pas", Value(pcs_tensor_in_pas_));
	Pair tensor_diagonalized("pcs_tensor_diagonalized", Value(pcs_tensor_diagonalized_));
	Pair tensor_reconfigured("pcs_tensor_reconfigured", Value(pcs_tensor_reconfigured_));
	Pair metal_center("metal_center",numeric::serialize(metal_center_));
	Pair nmr_tensor_data("nmr_tensor_data",NMRTensor::serialize());

	return Value(utility::tools::make_vector(tensor_in_arbitrary,tensor_in_pas,tensor_diagonalized,tensor_reconfigured,metal_center,nmr_tensor_data));
}

/// @brief deserialize a json_spirit object to a PCSTensor
void
PCSTensor::deserialize(utility::json_spirit::mObject data)
{
	pcs_tensor_in_arbitrary_frame_ = data["pcs_tensor_in_arbitrary_frame"].get_bool();
	pcs_tensor_in_pas_             = data["pcs_tensor_in_pas"].get_bool();
	pcs_tensor_diagonalized_       = data["pcs_tensor_diagonalized"].get_bool();
	pcs_tensor_reconfigured_       = data["pcs_tensor_reconfigured"].get_bool();
	metal_center_                  = numeric::deserialize<core::Real>(data["metal_center"].get_array());
	NMRTensor::deserialize(data["nmr_tensor_data"].get_obj());
}

/// @brief converts delta chi tensor principal values into its axial and rhombic component
void
PCSTensor::ax_rh_from_eigenvalues() {
	Real Eig_zz = get_Eig_zz();
	Real Eig_yy = get_Eig_yy();
	Real Eig_xx = get_Eig_xx();
	Real ax = Eig_zz - (Eig_xx + Eig_yy)/2.0;
	Real rh = Eig_xx - Eig_yy;
	set_ax(ax);
	set_rh(rh);
}

/// @brief converts axial and rhombic component into the delta chi tensor principal values
void
PCSTensor::ax_rh_to_eigenvalues() {
	Real ax = get_ax();
	Real rh = get_rh();
	Real Eig_xx = -ax/3.0 + rh/2.0;
	Real Eig_yy = -ax/3.0 - rh/2.0;
	Real Eig_zz = 2.0/3.0 * ax;
	set_Eig_xx(Eig_xx);
	set_Eig_yy(Eig_yy);
	set_Eig_zz(Eig_zz);
}

/// @brief return the delta chi tensor in 3x3 matrix form
Matrix
PCSTensor::build_tensor_matrix() {
	Matrix chi_;
	chi_(1,1) = get_T_xx(); chi_(1,2) = get_T_xy(); chi_(1,3) = get_T_xz();
	chi_(2,1) = get_T_xy(); chi_(2,2) = get_T_yy(); chi_(2,3) = get_T_yz();
	chi_(3,1) = get_T_xz(); chi_(3,2) = get_T_yz(); chi_(3,3) = -get_T_xx() - get_T_yy();
	return chi_;
}

} // namespace pcs
} // namespace nmr
} // namespace scoring
} // namespace core

