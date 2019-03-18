// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/RDCTensor.cc
/// @brief   Implementation of class NMRTensor
/// @details last Modified: 06/08/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

// Unit headers
#include <core/scoring/nmr/rdc/RDCTensor.hh>
#include <core/scoring/nmr/util.hh>

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

// ObjexxFCL headers
#include <ObjexxFCL/format.hh>

// C++ headers
#include <string>
#include <iostream>

namespace core {
namespace scoring {
namespace nmr {
namespace rdc {

static basic::Tracer TR( "core.scoring.nmr.rdc.RDCTensor" );

/// @brief default constructor
RDCTensor::RDCTensor() :
	NMRTensor(),
	Da_(0),
	R_(0),
	Dmax_(0),
	rdc_tensor_in_arbitrary_frame_(false),
	rdc_tensor_in_pas_(false),
	dmax_set_(false),
	rdc_tensor_diagonalized_(false),
	rdc_tensor_reconfigured_(false)
{}

/// @brief constructor with arguments
///        initializes alignment tensor with upper diagonal elements
///        vector should contain tensor parameters in the order A_xx, A_xy, A_xz, A_yy, A_yz
RDCTensor::RDCTensor(utility::vector1<Real> const & tensor_params) :
	Da_(0),
	R_(0),
	Dmax_(0),
	rdc_tensor_in_arbitrary_frame_(false),
	rdc_tensor_in_pas_(false),
	dmax_set_(false),
	rdc_tensor_diagonalized_(false),
	rdc_tensor_reconfigured_(false)
{
	set_tensor_in_arbitrary_frame(tensor_params);
}

/// @brief copy constructor
RDCTensor::RDCTensor(RDCTensor const & other) :
	NMRTensor(other),
	Da_(other.Da_),
	R_(other.R_),
	Dmax_(other.Dmax_),
	rdc_tensor_in_arbitrary_frame_(other.rdc_tensor_in_arbitrary_frame_),
	rdc_tensor_in_pas_(other.rdc_tensor_in_pas_),
	dmax_set_(other.dmax_set_),
	rdc_tensor_diagonalized_(other.rdc_tensor_diagonalized_),
	rdc_tensor_reconfigured_(other.rdc_tensor_reconfigured_)
{}

/// @brief assignment operator
RDCTensor &
RDCTensor::operator=(RDCTensor const & rhs) {
	if ( this != &rhs ) {
		NMRTensor::operator =(rhs);
		Da_ = rhs.Da_;
		R_ = rhs.R_;
		Dmax_ = rhs.Dmax_;
		rdc_tensor_in_arbitrary_frame_ = rhs.rdc_tensor_in_arbitrary_frame_;
		rdc_tensor_in_pas_ = rhs.rdc_tensor_in_pas_;
		dmax_set_ = rhs.dmax_set_;
		rdc_tensor_diagonalized_ = rhs.rdc_tensor_diagonalized_;
		rdc_tensor_reconfigured_ = rhs.rdc_tensor_reconfigured_;
	}
	return *this;
}

/// @brief virtual destructor
RDCTensor::~RDCTensor() {}

/// @brief every tensor should offer a function that returns its name
std::string
RDCTensor::tensor_name() const {
	return "RDCTensor";
}

/// @brief sets the alignment tensor as it exists in the arbitrary frame
///        resets the tensor upper diagonal elements
///        vector should contain tensor parameters in the order A_xx, A_xy, A_xz, A_yy, A_yz
void
RDCTensor::set_tensor_in_arbitrary_frame(utility::vector1<Real> const & tensor_params) {
	runtime_assert_msg(tensor_params.size() >= 5, "ERROR: Vector contains not enough parameter to instantiate RDCTensor object.");
	set_T_xx(tensor_params[1]);
	set_T_xy(tensor_params[2]);
	set_T_xz(tensor_params[3]);
	set_T_yy(tensor_params[4]);
	set_T_yz(tensor_params[5]);
	rdc_tensor_in_arbitrary_frame_ = true;
}

/// @brief sets the alignment tensor as it exists in the principal axis system (PAS)
///        resets the alignment tensor principal values, the alignment order, the rhombicity,
///        the three Euler angles and the dipolar coupling constant
///        vector should contain tensor parameters in the order alpha, beta, gamma, Da, R, Dmax
void
RDCTensor::set_tensor_in_pas(utility::vector1<Real> const & tensor_params) {
	runtime_assert_msg(tensor_params.size() >= 6, "ERROR: Vector contains not enough parameter to instantiate RDCTensor object.");
	Da_ = tensor_params[4];
	runtime_assert_msg( ( (0.0 <= tensor_params[5]) && (tensor_params[5] <= 2.0/3.0) ), "ERROR: Rhombicity must be in range (0, 2/3].");
	R_ = tensor_params[5];
	set_Dmax(tensor_params[6]);
	set_alpha(tensor_params[1]);
	set_beta(tensor_params[2]);
	set_gamma(tensor_params[3]);
	Da_R_to_ax_rh_and_eigenvalues();
	rdc_tensor_in_pas_ = true;
	rdc_tensor_diagonalized_ = true;
}

/// @brief brings alignment tensor in principal axis frame by diagonalization
///        make sure to call either set_tensor_in_arbitrary_frame() and set_Dmax() or
///        set_tensor_in_pas() first
void
RDCTensor::diagonalize_tensor() {
	if ( !rdc_tensor_in_arbitrary_frame_ && !rdc_tensor_in_pas_ ) {
		utility_exit_with_message("ERROR: RDCTensor is not set. Call \"set_tensor_in_arbitray_frame()\" or \"set_tensor_in_pas()\" before diagonalization.");
	}
	if ( !dmax_set_ ) {
		utility_exit_with_message("ERROR: Dmax is not set. Call \"set_Dmax()\" before diagonalization.");
	}
	if ( rdc_tensor_in_arbitrary_frame_ ) {
		Matrix A_ = build_tensor_matrix();
		Real tol(0.0);
		Matrix evecs;
		Vector evals = numeric::eigenvector_jacobi(A_, tol, evecs);
		evecs.transpose();
		Vector angles = euler_angles_from_rotation_matrix(evecs, get_euler_convention());
		set_Eig_xx(evals(1));
		set_Eig_yy(evals(2));
		set_Eig_zz(evals(3));
		Da_R_ax_rh_from_eigenvalues();
		rdc_tensor_diagonalized_ = true;
		set_alpha(angles(1));
		set_beta(angles(2));
		set_gamma(angles(3));
		rdc_tensor_in_pas_ = true;
	} else {
		return;
	}
}

/// @brief brings alignment tensor principal values in correct order
///        and reconfigure tensor into unique tensor representation
///        make sure to call set_tensor_in_pas() or diagonalize_tensor() first
void
RDCTensor::reorder_tensor() {
	if ( !rdc_tensor_reconfigured_ ) {

		if ( !rdc_tensor_diagonalized_ ) {
			TR.Warning << "RDCTensor is not diagonalized. Call \"diagonalize_tensor()\" first" << std::endl;
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

		Da_R_ax_rh_from_eigenvalues();

		rdc_tensor_reconfigured_ = true;
	} else {
		return;
	}
}

/// @brief shows summary of tensor statistics
/// @details If show_in_pas is true and tensor has been set or transformed in the
///          principal axis system, show tensor parameters in PAS. If either of
///          the two conditions is not fulfilled show tensor matrix in arbitrary frame.
void
RDCTensor::show_tensor_stats(std::ostream & tracer, bool show_in_pas) const {
	using namespace ObjexxFCL;
	using namespace ObjexxFCL::format;
	Size width(10);
	Size precision1(6);
	Size precision2(3);
	tracer << "   * * * RDCTensor parameters * * *   "  << std::endl;
	if ( show_in_pas && ( rdc_tensor_in_pas_ || ( rdc_tensor_diagonalized_ && rdc_tensor_reconfigured_ ) ) ) {
		tracer << "Da    (Hz)   = " << F(width, precision2, get_Da()) << std::endl;
		tracer << "R            = " << F(width, precision2, get_R()) << std::endl;
		tracer << "Aa    (Hz)   = " << E(width, precision1, get_ax()) << std::endl;
		tracer << "Ar    (Hz)   = " << E(width, precision1, get_rh()) << std::endl;
		tracer << "alpha (deg.) = " << F(width, precision2, get_alpha()) << std::endl;
		tracer << "beta  (deg.) = " << F(width, precision2, get_beta()) << std::endl;
		tracer << "gamma (deg.) = " << F(width, precision2, get_gamma()) << std::endl;
	} else {
		tracer << "[[ " << E(width, precision1, get_T_xx()) << ", "
			<< E(width, precision1, get_T_xy()) << ", "
			<< E(width, precision1, get_T_xz()) << " ]," << std::endl;
		tracer << " [ " << E(width, precision1, get_T_xy()) << ", "
			<< E(width, precision1, get_T_yy()) << ", "
			<< E(width, precision1, get_T_yz()) << " ]," << std::endl;
		tracer << " [ " << E(width, precision1, get_T_xz()) << ", "
			<< E(width, precision1, get_T_yz()) << ", "
			<< E(width, precision1, -get_T_xx()-get_T_yy()) << " ]]" << std::endl;
	}
}

/// @brief serialize a PCSTensor to a json_spirit object
utility::json_spirit::Value
RDCTensor::serialize() const
{
	using utility::json_spirit::Value;
	using utility::json_spirit::Pair;

	Pair magnitude("magnitude", Value(Da_));
	Pair rhombicity("rhombicity", Value(R_));
	Pair Dmax("Dmax", Value(Dmax_));
	Pair tensor_in_arbitrary("rdc_tensor_in_arbitrary_frame", Value(rdc_tensor_in_arbitrary_frame_));
	Pair tensor_in_pas("rdc_tensor_in_pas", Value(rdc_tensor_in_pas_));
	Pair tensor_diagonalized("rdc_tensor_diagonalized", Value(rdc_tensor_diagonalized_));
	Pair tensor_reconfigured("rdc_tensor_reconfigured", Value(rdc_tensor_reconfigured_));
	Pair nmr_tensor_data("nmr_tensor_data",NMRTensor::serialize());
	return Value(utility::tools::make_vector(magnitude,rhombicity,Dmax,tensor_in_arbitrary,tensor_in_pas,tensor_diagonalized,tensor_reconfigured,nmr_tensor_data));
}

/// @brief deserialize a json_spirit object to a PCSTensor
void
RDCTensor::deserialize(utility::json_spirit::mObject data)
{
	Da_                            = data["magnitude"].get_real();
	R_                             = data["rhombicity"].get_real();
	Dmax_                          = data["Dmax"].get_real();
	rdc_tensor_in_arbitrary_frame_ = data["rdc_tensor_in_arbitrary_frame"].get_bool();
	rdc_tensor_in_pas_             = data["rdc_tensor_in_pas"].get_bool();
	rdc_tensor_diagonalized_       = data["rdc_tensor_diagonalized"].get_bool();
	rdc_tensor_reconfigured_       = data["rdc_tensor_reconfigured"].get_bool();
	NMRTensor::deserialize(data["nmr_tensor_data"].get_obj());
}

/// @brief converts alignment tensor principal values into axial and rhombic component, alignment order and rhombicity
void
RDCTensor::Da_R_ax_rh_from_eigenvalues() {
	if ( !dmax_set_ ) {
		utility_exit_with_message("ERROR: Dmax for RDCTensor is not set. Call \"set_Dmax()\" first.");
	}
	Real Eig_xx = get_Eig_xx();
	Real Eig_yy = get_Eig_yy();
	Real Eig_zz = get_Eig_zz();
	Real ax = Eig_zz - (Eig_xx + Eig_yy)/2.0;
	Real rh = (Eig_xx - Eig_yy);
	set_ax(ax);
	set_rh(rh);
	Da_ = 0.5 * Dmax_ * ax;
	R_ = rh/ax;
}

/// @brief calculates alignment tensor principal values, its axial and rhombic component from the alignment order and rhombicity
void
RDCTensor::Da_R_to_ax_rh_and_eigenvalues() {
	if ( !dmax_set_ ) {
		utility_exit_with_message("ERROR: Dmax for RDCTensor is not set. Call \"set_Dmax()\" first.");
	}
	Real Eig_zz = (4.0 * Da_) / (3.0 * Dmax_);
	Real Eig_xx = Da_/Dmax_ * (-2.0/3.0 + R_);
	Real Eig_yy = Da_/Dmax_ * (-2.0/3.0 - R_);
	set_Eig_xx(Eig_xx);
	set_Eig_yy(Eig_yy);
	set_Eig_zz(Eig_zz);
	Real ax = 2.0 * Da_ / Dmax_;
	Real rh = R_ * 2.0 * Da_ / Dmax_;
	set_ax(ax);
	set_rh(rh);
}

/// @brief returns the alignment tensor in 3x3 matrix form
Matrix
RDCTensor::build_tensor_matrix() {
	Matrix A_;
	A_(1,1) = get_T_xx(); A_(1,2) = get_T_xy(); A_(1,3) = get_T_xz();
	A_(2,1) = get_T_xy(); A_(2,2) = get_T_yy(); A_(2,3) = get_T_yz();
	A_(3,1) = get_T_xz(); A_(3,2) = get_T_yz(); A_(3,3) = -get_T_xx() - get_T_yy();
	return A_;
}

} // namespace rdc
} // namespace nmr
} // namespace scoring
} // namespace core

