// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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

using namespace core;

struct ContactStruct
{
	Real k;
	Vector r21;
	Vector c21;
	Real d2;
	Size tor1, tor2;
};

class NormalMode
{
public:

	NormalMode();

	NormalMode( std::string const mode,
		Real distcut );

	~NormalMode();

	// Option setups
	void set_harmonic_constants( Real const &k_uniform );

	void set_harmonic_constants( Real const &k_short,
		Real const &k_SS,
		Real const &k_long );

	void set_torsions_using( Size const seqpos ){
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
	Size natm() const { return xyz_.size(); }
	Size ntor() const { return e_.size(); }

	Size nmode() const {
		if ( torsion() ) {
			return ntor();
		} else {
			return 3*natm()-6;
		}
	}

	std::string mode() const { return mode_; }
	Real dist2_cut() const { return dist2_cut_; }

	// Index
	utility::vector1< id::AtomID > get_atomID( ) const { return atomID_; }
	utility::vector1< id::TorsionID > get_torID( ) const { return torID_; }

	// Cartesian eigenvector
	utility::vector1< utility::vector1< Vector > > get_eigvec_cart() const { return eigvec_cart_; }
	utility::vector1< Vector > get_eigvec_cart( Size const imode ) const { return eigvec_cart_[imode+6]; }

	// Torsion eigenvector
	utility::vector1< utility::vector1< Real > > get_eigvec_tor() const { return eigvec_tor_; }
	utility::vector1< Real > get_eigvec_tor( Size const imode ) const { return eigvec_tor_[imode]; }

	utility::vector1< Real > get_eigval() const { return eigval_; }
	Real get_eigval( Size const imode ) const {
		if ( torsion() ) {
			return eigval_[imode];
		} else {
			return eigval_[imode+6];
		}
	}

	utility::vector1< Real > get_importance() const { return importance_; }
	Real get_importance( Size const imode ) const { return importance_[imode]; }

	Real get_k( Size const i, Size const j ) const { return k_[i][j]; }

	Vector xyz( Size const i ) const { return xyz_[i]; }

	Size a_to_i( Size const a ) const { return a_to_i_[a]; }
	Size i_to_a( Size const i ) const { return i_to_a_[i]; }

private:

	void prepare_coord( pose::Pose const &pose );

	void set_harmonic_constant_map( pose::Pose const &pose );

	utility::vector1< utility::vector1< Real > > make_Hessian_ANM();
	utility::vector1< utility::vector1< Real > > make_Hessian_TNM();

	utility::vector1< utility::vector1< Real > >
	convert_to_torsion_crd( utility::vector1< utility::vector1< Real > > const &U );

	void
	update_inertia_tensor( Size const a,
		Real &Ma,
		Vector &Ra,
		numeric::MathMatrix< Real > &Icorr,
		numeric::MathMatrix< Real > &Ia,
		utility::vector1< Vector > const xyz_com );

	void
	calculate_Jacobi_correction( Size const a,
		Real const &Ma,
		Vector const &Ra,
		Real const &Msum,
		numeric::MathMatrix< Real > const &Ia,
		numeric::MathMatrix< Real > const &i_Isum,
		Vector &Aa,
		Vector &ta );

	numeric::MathMatrix< Real >
	calc_inertia_tensor( numeric::MathMatrix< Real > const Icorr,
		Vector const &Ra,
		Real const Msum );


	Real
	pythag( Real a, Real b );

	inline
	Real
	sign( Real a, Real b )
	{ return b >= 0.0 ? std::fabs(a) : -std::fabs(a); }

	void
	svdcmp( utility::vector1< utility::vector1< Real > > &a,
		Size const m, Size const n,
		utility::vector1< Real > &w,
		utility::vector1< utility::vector1< Real > > &v );

	void
	eigsrt( utility::vector1< utility::vector1 < Real > > &eigvec,
		utility::vector1< Real > &eigval );

public:

private:
	bool torsion_;
	std::string mode_;
	Real dist2_cut_;

	bool use_uniform_k_;
	Real k_uniform_;
	Real k_short_;
	Real k_SS_;
	Real k_long_;

	// TNM
	bool use_phi_, use_psi_;
	bool eckart_correction_;
	utility::vector1< Size > torsions_using_;
	bool torsions_using_assigned_;

	utility::vector1< utility::vector1< Real > > k_; // harmonic constant

	utility::vector1< Size > a_to_i_;
	utility::vector1< Size > i_to_a_;
	utility::vector1< Vector > e_;
	utility::vector1< Vector > s_;
	utility::vector1< Vector > tau_;
	utility::vector1< Vector > xyz_;
	utility::vector1< Real > importance_;

	// Index
	utility::vector1< id::AtomID > atomID_;
	utility::vector1< id::TorsionID > torID_;

	utility::vector1< Real > eigval_;
	utility::vector1< utility::vector1< Vector > > eigvec_cart_;
	utility::vector1< utility::vector1< Real > > eigvec_tor_;
}; // end class NormalMode

} // moves
} // protocols

#endif
