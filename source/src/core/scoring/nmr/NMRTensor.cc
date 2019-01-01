// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/NMRTensor.cc
/// @brief   Implementation of class NMRTensor
/// @details last Modified: 06/08/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

// Unit headers
#include <core/scoring/nmr/NMRTensor.hh>

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/tools/make_vector.hh>

namespace core {
namespace scoring {
namespace nmr {

/// @brief default constructor
NMRTensor::NMRTensor() :
	utility::pointer::ReferenceCount(),
	alpha_(0.0),
	beta_(0.0),
	gamma_(0.0),
	convention_(ZYZ_CONVENTION),
	ax_(0.0),
	rh_(0.0),
	T_xx_(0.0),
	T_xy_(0.0),
	T_xz_(0.0),
	T_yy_(0.0),
	T_yz_(0.0),
	Eig_xx_(0.0),
	Eig_yy_(0.0),
	Eig_zz_(0.0)
{}

/// @brief copy constructor
NMRTensor::NMRTensor(NMRTensor const & other) :
	utility::pointer::ReferenceCount(other),
	alpha_(other.alpha_),
	beta_(other.beta_),
	gamma_(other.gamma_),
	convention_(other.convention_),
	ax_(other.ax_),
	rh_(other.rh_),
	T_xx_(other.T_xx_),
	T_xy_(other.T_xy_),
	T_xz_(other.T_xz_),
	T_yy_(other.T_yy_),
	T_yz_(other.T_yz_),
	Eig_xx_(other.Eig_xx_),
	Eig_yy_(other.Eig_yy_),
	Eig_zz_(other.Eig_zz_)
{}

/// @brief assignment operator
NMRTensor &
NMRTensor::operator=(NMRTensor const & rhs) {
	if ( this != &rhs ) {
		utility::pointer::ReferenceCount::operator=(rhs);
		alpha_ = rhs.alpha_;
		beta_ = rhs.beta_;
		gamma_ = rhs.gamma_;
		convention_ = rhs.convention_;
		ax_ = rhs.ax_;
		rh_ = rhs.rh_;
		T_xx_ = rhs.T_xx_;
		T_xy_ = rhs.T_xy_;
		T_xz_ = rhs.T_xz_;
		T_yy_ = rhs.T_yy_;
		T_yz_ = rhs.T_yz_;
		Eig_xx_ = rhs.Eig_xx_;
		Eig_yy_ = rhs.Eig_yy_;
		Eig_zz_ = rhs.Eig_zz_;
	}
	return *this;
}
/// @brief virtual destructor
NMRTensor::~NMRTensor() {}

/// @brief bring NMR tensor in principal axis frame by diagonalization
void
NMRTensor::diagonalize_tensor() {}

/// @brief Bring the NMR tensor principal values in correct order
///        and reconfigure tensor into unique tensor representation
void
NMRTensor::reorder_tensor() {}

/// @brief show summary of tensor statistics
void
NMRTensor::show_tensor_stats(std::ostream &, bool ) const {}

/// @brief serialize an NMRTensor to a json_spirit object
utility::json_spirit::Value
NMRTensor::serialize() const
{
	using utility::json_spirit::Value;
	using utility::json_spirit::Pair;

	Pair alpha("alpha", Value(alpha_));
	Pair beta("beta", Value(beta_));
	Pair gamma("gamma", Value(gamma_));
	Pair euler("euler", Value(static_cast<int>(convention_)));
	Pair axial("axial", Value(ax_));
	Pair rhombic("rhombic", Value(rh_));
	Pair T_xx("T_xx", Value(T_xx_));
	Pair T_xy("T_xy", Value(T_xy_));
	Pair T_xz("T_xz", Value(T_xz_));
	Pair T_yy("T_yy", Value(T_yy_));
	Pair T_yz("T_yz", Value(T_yz_));
	Pair Eig_xx("Eig_xx", Value(Eig_xx_));
	Pair Eig_yy("Eig_yy", Value(Eig_yy_));
	Pair Eig_zz("Eig_zz", Value(Eig_zz_));
	return Value(utility::tools::make_vector(alpha,beta,gamma,euler,axial,rhombic,T_xx,T_xy,T_xz,T_yy,T_yz,Eig_xx,Eig_yy,Eig_zz));
}

/// @brief deserialize a json_spirit object to an NMRTensor
void
NMRTensor::deserialize(utility::json_spirit::mObject data)
{
	alpha_      = data["alpha"].get_real();
	beta_       = data["beta"].get_real();
	gamma_      = data["gamma"].get_real();
	convention_ = static_cast<EULER_CONVENTION>(data["euler"].get_int());
	ax_         = data["axial"].get_real();
	rh_         = data["rhombic"].get_real();
	T_xx_       = data["T_xx"].get_real();
	T_xy_       = data["T_xy"].get_real();
	T_xz_       = data["T_xz"].get_real();
	T_yy_       = data["T_yy"].get_real();
	T_yz_       = data["T_yz"].get_real();
	Eig_xx_     = data["Eig_xx"].get_real();
	Eig_yy_     = data["Eig_yy"].get_real();
	Eig_zz_     = data["Eig_zz"].get_real();
}

} // namespace nmr
} // namespace scoring
} // namespace core

