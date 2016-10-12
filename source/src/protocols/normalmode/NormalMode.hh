// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    protocols/normalmode/NormalMode.hh
/// @brief   initialization for NormalMode
/// @details
/// @author  Hahnbeom Park

#ifndef INCLUDED_protocols_normalmode_NormalMode_hh
#define INCLUDED_protocols_normalmode_NormalMode_hh

// Unit headers
#include <numeric/xyzVector.hh>
#include <numeric/MathMatrix.hh>
//#include <numeric/xyzMatrix.hh>

// Package headers

// Project headers
#include <core/pose/Pose.hh>
#include <core/id/AtomID.hh>
#include <core/id/TorsionID.hh>

//// C++ headers
#include <string>

#include <core/kinematics/MoveMap.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace normalmode {

struct ContactStruct
{
	core::Real k;
	core::Vector r21;
	core::Vector c21;
	core::Real d2;
	core::Size tor1, tor2;
};

class NormalMode
{
public:

	NormalMode();

	NormalMode( std::string const mode,
		core::Real distcut );

	~NormalMode();

	// Option setups
	void set_harmonic_constants( core::Real const &k_uniform );

	void set_harmonic_constants( core::Real const &k_short,
		core::Real const &k_SS,
		core::Real const &k_long );

	void set_torsions_using( core::Size const seqpos ){
		torsions_using_assigned_ = true;
		torsions_using_.push_back( seqpos );
	}

	void clear_torsions_using(){
		torsions_using_assigned_ = false;
		torsions_using_.resize( 0 );
	}

	void torsion( bool const bool_in,
		bool const use_phi = true,
		bool const use_psi = true,
		bool const eckart_correction = true );

	void eckart_correction( bool const bool_in ){ eckart_correction_ = bool_in; }

	void solve( core::pose::Pose const & pose );

	// Accessors
	bool torsion() const { return torsion_; }
	core::Size natm() const { return xyz_.size(); }
	core::Size ntor() const { return e_.size(); }

	core::Size nmode() const {
		if ( torsion() ) {
			return ntor();
		} else {
			return 3*natm()-6;
		}
	}

	std::string mode() const { return mode_; }
	core::Real dist2_cut() const { return dist2_cut_; }

	// Index
	utility::vector1< core::id::AtomID > get_atomID( ) const { return atomID_; }
	utility::vector1< core::id::TorsionID > get_torID( ) const { return torID_; }

	// Cartesian eigenvector
	utility::vector1< utility::vector1< core::Vector > > get_eigvec_cart() const { return eigvec_cart_; }
	utility::vector1< core::Vector > get_eigvec_cart( core::Size const imode ) const { return eigvec_cart_[imode+6]; }

	// Torsion eigenvector
	utility::vector1< utility::vector1< core::Real > > get_eigvec_tor() const { return eigvec_tor_; }
	utility::vector1< core::Real > get_eigvec_tor( core::Size const imode ) const { return eigvec_tor_[imode]; }

	utility::vector1< core::Real > get_eigval() const { return eigval_; }
	core::Real get_eigval( core::Size const imode ) const {
		if ( torsion() ) {
			return eigval_[imode];
		} else {
			return eigval_[imode+6];
		}
	}

	utility::vector1< core::Real > get_importance() const { return importance_; }
	core::Real get_importance( core::Size const imode ) const { return importance_[imode]; }

	core::Real get_k( core::Size const i, core::Size const j ) const { return k_[i][j]; }

	core::Vector xyz( core::Size const i ) const { return xyz_[i]; }

	core::Size a_to_i( core::Size const a ) const { return a_to_i_[a]; }
	core::Size i_to_a( core::Size const i ) const { return i_to_a_[i]; }

private:

	void prepare_coord( core::pose::Pose const &pose );

	void set_harmonic_constant_map( core::pose::Pose const &pose );

	utility::vector1< utility::vector1< core::Real > > make_Hessian_ANM();
	utility::vector1< utility::vector1< core::Real > > make_Hessian_TNM();

	utility::vector1< utility::vector1< core::Real > >
	convert_to_torsion_crd( utility::vector1< utility::vector1< core::Real > > const &U );

	void
	update_inertia_tensor( core::Size const a,
		core::Real &Ma,
		core::Vector &Ra,
		numeric::MathMatrix< core::Real > &Icorr,
		numeric::MathMatrix< core::Real > &Ia,
		utility::vector1< core::Vector > const xyz_com );

	void
	calculate_Jacobi_correction( core::Size const a,
		core::Real const &Ma,
		core::Vector const &Ra,
		core::Real const &Msum,
		numeric::MathMatrix< core::Real > const &Ia,
		numeric::MathMatrix< core::Real > const &i_Isum,
		core::Vector &Aa,
		core::Vector &ta );

	numeric::MathMatrix< core::Real >
	calc_inertia_tensor( numeric::MathMatrix< core::Real > const Icorr,
		core::Vector const &Ra,
		core::Real const Msum );

	void
	eigsrt( utility::vector1< utility::vector1 < core::Real > > &eigvec,
		utility::vector1< core::Real > &eigval );

public:

private:
	bool torsion_;
	std::string mode_;
	core::Real dist2_cut_;

	bool use_uniform_k_;
	core::Real k_uniform_;
	core::Real k_short_;
	core::Real k_SS_;
	core::Real k_long_;

	// TNM
	bool use_phi_, use_psi_;
	bool eckart_correction_;
	utility::vector1< core::Size > torsions_using_;
	bool torsions_using_assigned_;

	utility::vector1< utility::vector1< core::Real > > k_; // harmonic constant

	utility::vector1< core::Size > a_to_i_;
	utility::vector1< core::Size > i_to_a_;
	utility::vector1< core::Vector > e_;
	utility::vector1< core::Vector > s_;
	utility::vector1< core::Vector > tau_;
	utility::vector1< core::Vector > xyz_;
	utility::vector1< core::Real > importance_;

	// Index
	utility::vector1< core::id::AtomID > atomID_;
	utility::vector1< core::id::TorsionID > torID_;

	utility::vector1< core::Real > eigval_;
	utility::vector1< utility::vector1< core::Vector > > eigvec_cart_;
	utility::vector1< utility::vector1< core::Real > > eigvec_tor_;
}; // end class NormalMode

} // moves
} // protocols

#endif
