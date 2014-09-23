// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief

// Unit headers
#include <core/scoring/constraints/BigBinConstraint.hh>

// Package headers
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/constraints/BoundConstraint.hh>

// Project headers
#include <core/chemical/ResidueType.hh>
#include <core/pose/Pose.hh>

#include <core/scoring/func/Func.hh>
#include <core/scoring/func/Func.fwd.hh>

#include <numeric/trig.functions.hh>
#include <numeric/conversions.hh>

#include <utility/exit.hh>
#include <ObjexxFCL/string.functions.hh>

#include <utility/vector1.hh>


// C++ Headers

namespace core {
namespace scoring {
namespace constraints {

/// @brief one line definition "BigBin res_number bin_char sdev"
void BigBinConstraint::read_def(
	std::istream & in,
	core::pose::Pose const & pose,
	func::FuncFactory const & /*func_factory*/
) {
	my_csts_.clear();

  using core::id::AtomID;
	using numeric::conversions::radians;
	AtomID C0, N1, CA1, C1, N2, CA2;

	std::string dummy;
	in >> dummy >> res_ >> bin_ >> sdev_;

	//chu skip the rest of line since this is a single line defintion.
	while( in.good() && (in.get() != '\n') ) {}

	if ( res_ > 1 && res_ < pose.total_residue() - 1 ) {
		C0_  = AtomID(pose.residue_type(res_ - 1).atom_index("C" ), res_ - 1);
		N1_  = AtomID(pose.residue_type(res_    ).atom_index("N" ), res_    );
		CA1_ = AtomID(pose.residue_type(res_    ).atom_index("CA"), res_    );
		C1_  = AtomID(pose.residue_type(res_    ).atom_index("C" ), res_    );
		N2_  = AtomID(pose.residue_type(res_ + 1).atom_index("N" ), res_ + 1);
		CA2_ = AtomID(pose.residue_type(res_ + 1).atom_index("CA"), res_ + 1);
	}

	Real const two_pi( 6.283185 );
	Real omega_lower_ = radians( -175.0 );
	Real omega_upper_ = radians(  175.0 );

	using ObjexxFCL::string_of;

	if ( bin_ == 'O' ) {
		omega_lower_ = radians( -10.0 );
		omega_upper_ = radians(  10.0 );
	} else {
		Real psi_upper_( 0.0 ), psi_lower_( 0.0 );
		Real phi_upper_( 0.0 ), phi_lower_( 0.0 );
		if ( bin_ == 'G' ) {
			phi_lower_ = radians(    0.0 );
			phi_upper_ = radians(  180.0 );
			psi_lower_ = radians( -100.0 );
			psi_upper_ = radians(  100.0 );
		} else if ( bin_ == 'E' ) {
			phi_lower_ = radians(    0.0 );
			phi_upper_ = radians(  180.0 );
			psi_lower_ = radians(  100.0 );
			psi_upper_ = radians(  -90.0 );
		} else if ( bin_ == 'A' ) {
			phi_lower_ = radians( -180.0 );
			phi_upper_ = radians(    0.0 );
			psi_lower_ = radians(  -50.0 );
			psi_upper_ = radians(  -30.0 );
		} else if ( bin_ == 'B') {
			phi_lower_ = radians( -180.0 );
			phi_upper_ = radians(    0.0 );
			psi_lower_ = radians(  100.0 );
			psi_upper_ = radians(  175.0 );
		} else {
			utility_exit_with_message(
				"Error: don't recognize bin " + string_of(bin_) + "!"
			);
		}

		func::FuncOP phi_func( new PeriodicBoundFunc(
			phi_lower_, phi_upper_, sdev_, "phi_" + string_of(bin_), two_pi
		) );
		func::FuncOP psi_func( new PeriodicBoundFunc(
				psi_lower_, psi_upper_, sdev_, "psi_" + string_of(bin_), two_pi
		) );
		func::FuncOP omega_func( new PeriodicBoundFunc(
				omega_lower_, omega_upper_, sdev_, "omega_" + string_of(bin_), two_pi
		) );

		ConstraintOP phi_cst( new DihedralConstraint( C0_, N1_, CA1_, C1_, phi_func ) );
		ConstraintOP psi_cst( new DihedralConstraint( N1_, CA1_, C1_, N2_, psi_func ) );
		ConstraintOP omega_cst( new DihedralConstraint( CA1_, C1_, N2_, CA2_, omega_func ) );

		my_csts_.push_back( phi_cst );
		my_csts_.push_back( psi_cst );
		my_csts_.push_back( omega_cst );
	}
}

void BigBinConstraint::show( std::ostream & out ) const {
	out << res_ << ' ' << bin_ << ' ' << sdev_ << std::endl;
}

/////////////////////////////////////////////////////////////////////////////
void
BigBinConstraint::score(
	func::XYZ_Func const & xyz, EnergyMap const & weights, EnergyMap & emap
) const {
	using utility::vector1;
	typedef vector1< ConstraintOP >::const_iterator iter;
	for ( iter it = my_csts_.begin(), end = my_csts_.end(); it != end; ++it ) {
		(*it)->score( xyz, weights, emap );
	}
}

void
BigBinConstraint::fill_f1_f2(
	AtomID const & id,
	func::XYZ_Func const & xyz,
	Vector & F1,
	Vector & F2,
	EnergyMap const & weights
) const {
	// call fill_f1_f2 on my_csts_
	using utility::vector1;
	typedef vector1< ConstraintOP >::const_iterator iter;
	for ( iter it = my_csts_.begin(), end = my_csts_.end(); it != end; ++it ) {
		(*it)->fill_f1_f2( id, xyz, F1, F2, weights );
	}
}

} // constraints
} // scoring
} // core
