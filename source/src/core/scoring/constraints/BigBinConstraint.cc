// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>
#include <utility/vector1.srlz.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace constraints {

std::string BigBinConstraint::type() const {
	return "BigBin";
}

ConstraintOP
BigBinConstraint::clone() const {
	return ConstraintOP( new BigBinConstraint( *this ) );
}


bool BigBinConstraint::operator == ( Constraint const & other ) const
{
	if ( !       same_type_as_me( other ) ) return false;
	if ( ! other.same_type_as_me( *this ) ) return false;

	auto const & other_bbcst( static_cast< BigBinConstraint const & > (other) );
	if ( C0_   != other_bbcst.C0_   ) return false;
	if ( N1_   != other_bbcst.N1_   ) return false;
	if ( CA1_  != other_bbcst.CA1_  ) return false;
	if ( C1_   != other_bbcst.C1_   ) return false;
	if ( N2_   != other_bbcst.N2_   ) return false;
	if ( CA1_  != other_bbcst.CA1_  ) return false;
	if ( res_  != other_bbcst.res_  ) return false;
	if ( bin_  != other_bbcst.bin_  ) return false;
	if ( sdev_ != other_bbcst.sdev_ ) return false;

	if ( my_csts_.size() != other_bbcst.my_csts_.size() ) return false;
	for ( core::Size ii = 1; ii <= my_csts_.size(); ++ii ) {
		if ( *my_csts_[ii] != *other_bbcst.my_csts_[ii] ) return false;
	}
	return true;
}

bool BigBinConstraint::same_type_as_me( Constraint const & other ) const {
	return dynamic_cast< BigBinConstraint const * > ( &other );
}

///c-tor
BigBinConstraint::BigBinConstraint(
	AtomID C0,
	AtomID N1,
	AtomID CA1,
	AtomID C1,
	AtomID N2,
	AtomID CA2,
	char bin,
	ScoreType scotype /* = dihedral_constraint */
):
	Constraint( scotype ),
	C0_( C0 ),
	N1_( N1 ),
	CA1_( CA1 ),
	C1_( C1 ),
	N2_( N2 ),
	CA2_( CA2 ),
	bin_(bin)
{}

BigBinConstraint::BigBinConstraint() :
	Constraint( dihedral_constraint ),
	res_( 0 ),
	bin_( 'A' ),
	sdev_( 0.5 )
{}

BigBinConstraint::BigBinConstraint( Size const res, char const bin, core::Real const sdev ) :
	Constraint( dihedral_constraint ),
	res_( res ),
	bin_( bin ),
	sdev_( sdev )
{}

///
Size
BigBinConstraint::natoms() const {
	return 6;
}

id::AtomID const &
BigBinConstraint::atom( Size const n ) const
{
	switch( n ) {
	case 1 :
		return C0_;
	case 2 :
		return N1_;
	case 3 :
		return CA1_;
	case 4 :
		return C1_;
	case 5 :
		return N2_;
	case 6 :
		return CA2_;
	default :
		utility_exit_with_message( "BigBinConstraint::atom() bad argument" );
	}
	return C0_;
}


/// @brief one line definition "BigBin res_number bin_char sdev"
void BigBinConstraint::read_def(
	std::istream & in,
	core::pose::Pose const & pose,
	func::FuncFactory const & /*func_factory*/
) {
	my_csts_.clear();

	using core::id::AtomID;
	using numeric::conversions::radians;

	std::string dummy;
	in >> dummy >> res_ >> bin_ >> sdev_;

	//chu skip the rest of line since this is a single line defintion.
	while ( in.good() && (in.get() != '\n') ) {}

	if ( res_ > 1 && res_ < pose.size() - 1 ) {
		C0_  = AtomID(pose.residue_type(res_ - 1).atom_index("C" ), res_ - 1);
		N1_  = AtomID(pose.residue_type(res_    ).atom_index("N" ), res_    );
		CA1_ = AtomID(pose.residue_type(res_    ).atom_index("CA"), res_    );
		C1_  = AtomID(pose.residue_type(res_    ).atom_index("C" ), res_    );
		N2_  = AtomID(pose.residue_type(res_ + 1).atom_index("N" ), res_ + 1);
		CA2_ = AtomID(pose.residue_type(res_ + 1).atom_index("CA"), res_ + 1);
	}

	Real const two_pi( 6.283185 );
	Real omega_lower = radians( -175.0 );
	Real omega_upper = radians(  175.0 );

	using ObjexxFCL::string_of;

	if ( bin_ == 'O' ) {
		//The following assignments are never used, and were flagged by clang analysis.
		//Uncomment if ever they're to be used for something.
		/*omega_lower = radians( -10.0 );
		omega_upper = radians(  10.0 );*/
	} else {
		Real psi_upper( 0.0 ), psi_lower( 0.0 );
		Real phi_upper( 0.0 ), phi_lower( 0.0 );
		if ( bin_ == 'G' ) {
			phi_lower = radians(    0.0 );
			phi_upper = radians(  180.0 );
			psi_lower = radians( -100.0 );
			psi_upper = radians(  100.0 );
		} else if ( bin_ == 'E' ) {
			phi_lower = radians(    0.0 );
			phi_upper = radians(  180.0 );
			psi_lower = radians(  100.0 );
			psi_upper = radians(  -90.0 );
		} else if ( bin_ == 'A' ) {
			phi_lower = radians( -180.0 );
			phi_upper = radians(    0.0 );
			psi_lower = radians(  -50.0 );
			psi_upper = radians(  -30.0 );
		} else if ( bin_ == 'B' ) {
			phi_lower = radians( -180.0 );
			phi_upper = radians(    0.0 );
			psi_lower = radians(  100.0 );
			psi_upper = radians(  175.0 );
		} else {
			utility_exit_with_message(
				"Error: don't recognize bin " + string_of(bin_) + "!"
			);
		}

		func::FuncOP phi_func( new PeriodicBoundFunc(
			phi_lower, phi_upper, sdev_, "phi_" + string_of(bin_), two_pi
			) );
		func::FuncOP psi_func( new PeriodicBoundFunc(
			psi_lower, psi_upper, sdev_, "psi_" + string_of(bin_), two_pi
			) );
		func::FuncOP omega_func( new PeriodicBoundFunc(
			omega_lower, omega_upper, sdev_, "omega_" + string_of(bin_), two_pi
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
	for ( auto const & my_cst : my_csts_ ) {
		my_cst->score( xyz, weights, emap );
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
	for ( auto const & my_cst : my_csts_ ) {
		my_cst->fill_f1_f2( id, xyz, F1, F2, weights );
	}
}

BigBinConstraint::BigBinConstraint( BigBinConstraint const & src ) :
	Constraint( src ),
	C0_  ( src.C0_   ),
	N1_  ( src.N1_   ),
	CA1_ ( src.CA1_  ),
	C1_  ( src.C1_   ),
	N2_  ( src.N2_   ),
	CA2_ ( src.CA2_  ),
	res_ ( src.res_  ),
	bin_ ( src.bin_  ),
	sdev_( src.sdev_ ),
	my_csts_( src.my_csts_.size() )
{
	for ( Size ii = 1; ii <= src.my_csts_.size(); ++ii ) {
		my_csts_[ ii ] = src.my_csts_[ ii ]->clone();
	}
}

} // constraints
} // scoring
} // core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::constraints::BigBinConstraint::save( Archive & arc ) const {
	arc( cereal::base_class< Constraint >( this ) );
	arc( CEREAL_NVP( C0_ ) ); // AtomID
	arc( CEREAL_NVP( N1_ ) ); // AtomID
	arc( CEREAL_NVP( CA1_ ) ); // AtomID
	arc( CEREAL_NVP( C1_ ) ); // AtomID
	arc( CEREAL_NVP( N2_ ) ); // AtomID
	arc( CEREAL_NVP( CA2_ ) ); // AtomID
	arc( CEREAL_NVP( res_ ) ); // Size
	arc( CEREAL_NVP( bin_ ) ); // char
	arc( CEREAL_NVP( sdev_ ) ); // core::Real
	arc( CEREAL_NVP( my_csts_ ) ); // utility::vector1<ConstraintOP>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::constraints::BigBinConstraint::load( Archive & arc ) {
	arc( cereal::base_class< Constraint >( this ) );
	arc( C0_ ); // AtomID
	arc( N1_ ); // AtomID
	arc( CA1_ ); // AtomID
	arc( C1_ ); // AtomID
	arc( N2_ ); // AtomID
	arc( CA2_ ); // AtomID
	arc( res_ ); // Size
	arc( bin_ ); // char
	arc( sdev_ ); // core::Real
	arc( my_csts_ ); // utility::vector1<ConstraintOP>
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::constraints::BigBinConstraint );
CEREAL_REGISTER_TYPE( core::scoring::constraints::BigBinConstraint )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_constraints_BigBinConstraint )
#endif // SERIALIZATION
