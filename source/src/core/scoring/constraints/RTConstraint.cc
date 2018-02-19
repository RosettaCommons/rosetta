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
#include <core/scoring/constraints/RTConstraint.hh>

// Package Headers
#include <core/scoring/func/FuncFactory.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/func/XYZ_Func.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/EnergyMap.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/id/NamedAtomID.hh>
#include <core/id/NamedStubID.hh>

#include <core/kinematics/Stub.hh>
#include <core/kinematics/RT.hh>

#include <core/conformation/Conformation.hh>

// Utility Headers
#include <basic/Tracer.hh>

// Numeric Headers
#include <numeric/deriv/distance_deriv.hh>
#include <numeric/xyz.functions.hh>

// C++ Headers
#include <cstdlib>
#include <iostream>


#include <core/id/SequenceMapping.hh>
#include <utility>
#include <utility/vector1.hh>


static basic::Tracer tr( "core.scoring.constraints" );

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Numeric serialization headers
#include <numeric/xyz.serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace constraints {

RTConstraint::RTConstraint() :
	Constraint( atom_pair_constraint ),
	stub1_( id::StubID::BOGUS_STUB_ID() ),
	stub2_( id::StubID::BOGUS_STUB_ID() ),
	func_( /* NULL */ ) {}

RTConstraint::RTConstraint(
	id::StubID const & stub1,
	id::StubID const & stub2,
	kinematics::RT const & rt_target,
	func::FuncOP func,
	ScoreType const & score_type
):
	Constraint( score_type ),
	stub1_( stub1 ),
	stub2_( stub2 ),
	rt_target_( rt_target ),
	func_(std::move( func ))
{}

RTConstraint::RTConstraint(
	RTConstraint const & src
):
	Constraint( src ),
	stub1_( src.stub1_ ),
	stub2_( src.stub2_ ),
	rt_target_( src.rt_target_ ),
	func_( src.func_ ? src.func_->clone() : src.func_ )
{}

RTConstraint::~RTConstraint() = default;

ConstraintOP
RTConstraint::clone() const {
	return ConstraintOP( new RTConstraint( *this ) );
}

ConstraintOP
RTConstraint::clone( func::FuncOP new_func ) const {
	return ConstraintOP( new RTConstraint( stub1_, stub2_, rt_target_, new_func, this->score_type() ));
}

kinematics::RT
RTConstraint::calculate_rt( func::XYZ_Func const & xyz ) const
{
	kinematics::Stub stub1(
		xyz( stub1_.atom( 1 ) ),
		xyz( stub1_.atom( 2 ) ),
		xyz( stub1_.atom( 3 ) )
	);

	kinematics::Stub stub2(
		xyz( stub2_.atom( 1 ) ),
		xyz( stub2_.atom( 2 ) ),
		xyz( stub2_.atom( 3 ) )
	);

	kinematics::RT rt(
		stub1,
		stub2
	);

	return rt;
}

Real
RTConstraint::dist( func::XYZ_Func const & xyz ) const
{
	kinematics::RT rt_instant( calculate_rt( xyz ) );

	return std::sqrt( rt_target_.distance_squared( rt_instant ) );
}

void
RTConstraint::score( func::XYZ_Func const & xyz, EnergyMap const &, EnergyMap & emap ) const
{
	emap[ this->score_type() ] += func_->func( dist( xyz ) );
}

std::string
RTConstraint::type() const {
	return "RTConstraint";
}

Size
RTConstraint::natoms() const
{
	return 6;
}

///
id::AtomID const &
RTConstraint::atom( Size const i ) const
{
	if ( i > 0 && i < 4 ) {
		return stub1_.atom( i );
	} else if ( i > 3 && i < 7 ) {
		return stub2_.atom( i - 3 );
	} else {
		utility_exit_with_message( "RTConstraint::atom() -- index is out of bounds" );
	}
}

// atom deriv
void
RTConstraint::fill_f1_f2(
	AtomID const & atom,
	func::XYZ_Func const & xyz,
	Vector & F1,
	Vector & F2,
	EnergyMap const & weights
) const
{
	utility_exit_with_message( " derivative of RTConstraint not supported yet " );

	// The following code is an attempt at calculating a derivative
	// Tests are already written
	tr.Debug << "atom:" << atom << std::endl;

	kinematics::RT rt_instant( calculate_rt(xyz) );
	kinematics::RT rt_target( rt_target_ );

	Vector xyz_instant(0.0), xyz_target(0.0);

	if ( atom == stub1_.atom( 1 ) || atom == stub1_.atom( 2 ) || atom == stub1_.atom( 3 ) ) {
		rt_instant.inverse(); rt_target.inverse();

		kinematics::Stub stub_instant(rt_instant);
		kinematics::Stub stub_target(rt_target);

		Real length1( xyz( stub1_.atom( 1 ) ).distance( xyz( stub1_.atom( 2 ) ) ));
		Real length2( xyz( stub1_.atom( 2 ) ).distance( xyz( stub1_.atom( 3 ) ) ));
		Real angle(numeric::angle_degrees(xyz( stub1_.atom( 1 ) ),
			xyz( stub1_.atom( 2 ) ),
			xyz( stub1_.atom( 3 ) )));

		if ( atom == stub1_.atom( 1 ) ) {
			tr.Debug << "stub1 atom1" << std::endl;
			xyz_instant = stub_instant.build_fake_xyz( 1, length1, length2, angle );
			xyz_target = stub_target.build_fake_xyz( 1, length1, length2, angle );
		} else if ( atom == stub1_.atom( 2 ) ) {
			tr.Debug << "stub1 atom2" << std::endl;
			xyz_instant = stub_instant.build_fake_xyz( 2, length1, length2, angle );
			xyz_target = stub_target.build_fake_xyz( 2, length1, length2, angle );
		} else {
			tr.Debug << "stub1 atom3" << std::endl;
			xyz_instant = stub_instant.build_fake_xyz( 3, length1, length2, angle );
			xyz_target = stub_target.build_fake_xyz( 3, length1, length2, angle );
		}
	} else if ( atom == stub2_.atom( 1 ) || atom == stub2_.atom( 2 ) || atom == stub2_.atom( 3 ) ) {
		kinematics::Stub stub_instant(rt_instant);
		kinematics::Stub stub_target(rt_target_);

		Real length1( xyz( stub2_.atom( 1 ) ).distance( xyz( stub2_.atom( 2 ) ) ));
		Real length2( xyz( stub2_.atom( 2 ) ).distance( xyz( stub2_.atom( 3 ) ) ));
		Real angle(numeric::angle_degrees(xyz( stub2_.atom( 1 ) ),
			xyz( stub2_.atom( 2 ) ),
			xyz( stub2_.atom( 3 ) )));

		if ( atom == stub2_.atom( 1 ) ) {
			tr.Debug << "stub2 atom1" << std::endl;
			xyz_instant = stub_instant.build_fake_xyz( 1, length1, length2, angle );
			xyz_target = stub_target.build_fake_xyz( 1, length1, length2, angle );
		} else if ( atom == stub2_.atom( 2 ) ) {
			tr.Debug << "stub2 atom2" << std::endl;
			xyz_instant = stub_instant.build_fake_xyz( 2, length1, length2, angle );
			xyz_target = stub_target.build_fake_xyz( 2, length1, length2, angle );
		} else {
			tr.Debug << "stub2 atom3" << std::endl;
			xyz_instant = stub_instant.build_fake_xyz( 3, length1, length2, angle );
			xyz_target = stub_target.build_fake_xyz( 3, length1, length2, angle );
		}
	} else {
		tr.Debug << "no atom" << std::endl;
		return;
	}
	tr.Debug << "xyz_instant:" << xyz_instant.to_string() << std::endl;
	tr.Debug << "xyz_target:" << xyz_target.to_string() << std::endl;

	Real dist(0.0); Vector f1(0.0), f2(0.0);

	numeric::deriv::distance_f1_f2_deriv( xyz_instant, xyz_target, dist, f1, f2 );
	tr.Debug << "dist:" << dist << std::endl;
	tr.Debug << "f1:" << f1.to_string() << std::endl;
	tr.Debug << "f2:" << f2.to_string() << std::endl;

	Real wderiv( weights[ this->score_type() ] * func_->dfunc( dist ));
	tr.Debug << "wderiv:" << wderiv << std::endl;

	F1 += wderiv * f1;
	tr.Debug << "F1:" << F1.to_string() << std::endl;
	F2 += wderiv * f2;
	tr.Debug << "F2:" << F2.to_string() << std::endl;
}

bool RTConstraint::operator == ( Constraint const & rhs ) const
{
	if ( !     same_type_as_me(   rhs ) ) return false;
	if ( ! rhs.same_type_as_me( *this ) ) return false;

	auto const & rhs_rtc( static_cast< RTConstraint const & > ( rhs ));
	if ( stub1_ != rhs_rtc.stub1_ ) return false;
	if ( stub2_ != rhs_rtc.stub2_ ) return false;
	if ( rt_target_.get_rotation() != rhs_rtc.rt_target_.get_rotation() ) return false;
	if ( rt_target_.get_translation() != rhs_rtc.rt_target_.get_translation() ) return false;
	if ( score_type() != rhs.score_type() ) return false;

	return func_ == rhs_rtc.func_ || ( func_ && rhs_rtc.func_ && *func_ == *rhs_rtc.func_ );
}

bool RTConstraint::same_type_as_me( Constraint const & other ) const
{
	return dynamic_cast< RTConstraint const * > ( &other );
}

ConstraintOP
RTConstraint::remapped_clone( pose::Pose const& src, pose::Pose const& dest, id::SequenceMappingCOP smap ) const {
	id::NamedAtomID atom1( core::pose::atom_id_to_named_atom_id( atom(1), src ));
	id::NamedAtomID atom2( core::pose::atom_id_to_named_atom_id( atom(2), src ));
	id::NamedAtomID atom3( core::pose::atom_id_to_named_atom_id( atom(3), src ));
	id::NamedAtomID atom4( core::pose::atom_id_to_named_atom_id( atom(4), src ));
	id::NamedAtomID atom5( core::pose::atom_id_to_named_atom_id( atom(5), src ));
	id::NamedAtomID atom6( core::pose::atom_id_to_named_atom_id( atom(6), src ));

	if ( smap ) {
		atom1.rsd() = (*smap)[ atom(1).rsd() ];
		atom2.rsd() = (*smap)[ atom(2).rsd() ];
		atom3.rsd() = (*smap)[ atom(3).rsd() ];
		atom4.rsd() = (*smap)[ atom(4).rsd() ];
		atom5.rsd() = (*smap)[ atom(5).rsd() ];
		atom6.rsd() = (*smap)[ atom(6).rsd() ];
	}

	id::AtomID id1( core::pose::named_atom_id_to_atom_id( atom1, dest ));
	id::AtomID id2( core::pose::named_atom_id_to_atom_id( atom2, dest ));
	id::AtomID id3( core::pose::named_atom_id_to_atom_id( atom3, dest ));
	id::AtomID id4( core::pose::named_atom_id_to_atom_id( atom4, dest ));
	id::AtomID id5( core::pose::named_atom_id_to_atom_id( atom5, dest ));
	id::AtomID id6( core::pose::named_atom_id_to_atom_id( atom6, dest ));

	if ( id1.valid() && id2.valid() && id3.valid() && id4.valid() ) {
		return ConstraintOP( new RTConstraint( id::StubID( id1, id2, id3 ),
			id::StubID( id4, id5, id6 ),
			rt_target_,
			func_ ? func_->clone() : func_,
			this->score_type() ) );
	} else {
		return nullptr;
	}
}

void RTConstraint::show( std::ostream& out ) const
{
	out << "--- RTConstraint ---" << std::endl
		<< "stub1:" << stub1_ << std::endl
		<< "stub2:" << stub2_ << std::endl
		<< "rt_target" << rt_target_ << std::endl
		<< "func:" << func_ << std::endl
		<< "--------------------" << std::endl;
}

func::Func const & RTConstraint::get_func() const {
	return *func_;
}

}
}
}


#ifdef    SERIALIZATION
/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::constraints::RTConstraint::save( Archive & arc ) const {
	arc( cereal::base_class< Constraint >( this ) );
	arc( CEREAL_NVP( stub1_ ) );
	arc( CEREAL_NVP( stub2_ ) );
	arc( CEREAL_NVP( rt_target_ ) );
	arc( CEREAL_NVP( func_ ) );
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::constraints::RTConstraint::load( Archive & arc ) {
	arc( cereal::base_class< Constraint >( this ) );
	arc( stub1_ );
	arc( stub2_ );
	arc( rt_target_ );
	arc( func_ );
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::constraints::RTConstraint );
CEREAL_REGISTER_TYPE( core::scoring::constraints::RTConstraint )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_constraints_RTConstraint )
#endif // SERIALIZATION
