// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief Inserts a Fragment into a Pose, similar to old Rosetta++ main_frag_trial algorithm.
/// @author Oliver Lange
/// @author James Thompson

// Unit Headers
#include <core/scoring/constraints/BoundConstraint.hh>

// Package Headers
#include <core/scoring/func/Func.hh>

// Project Headers
#include <basic/Tracer.hh>
#include <utility/exit.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

// Utility headers

static THREAD_LOCAL basic::Tracer tr( "core.constraints.BoundFunc", basic::t_info );

// C++ headers

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace constraints {

BoundFunc::BoundFunc( Real const lb, Real const ub, Real sd, std::string type ):
	lb_( lb ), ub_( ub ), sd_ ( sd ), rswitch_( 0.5 ), type_( type )
{}

BoundFunc::BoundFunc( Real const lb, Real const ub, Real sd, Real rswitch, std::string type ) :
	lb_( lb ), ub_( ub ), sd_ ( sd ), rswitch_( rswitch ), type_( type )
{}

func::FuncOP
BoundFunc::clone() const { return func::FuncOP( new BoundFunc( *this ) ); }

bool BoundFunc::operator==( Func const & rhs ) const
{
	if ( !     same_type_as_me(   rhs ) ) return false;
	if ( ! rhs.same_type_as_me( *this ) ) return false;

	BoundFunc const & rhs_bound( static_cast< BoundFunc const & > (rhs) );
	if ( lb_ != rhs_bound.lb_ ) return false;
	if ( ub_ != rhs_bound.ub_ ) return false;
	if ( sd_ != rhs_bound.sd_ ) return false;
	if ( rswitch_ != rhs_bound.rswitch_ ) return false;
	return type_ == rhs_bound.type_;
}

bool BoundFunc::same_type_as_me( Func const & other ) const {
	return dynamic_cast< BoundFunc const * > (&other);
}

Real
BoundFunc::func( Real const x ) const
{

	//Real const rswitch_offset (( rswitch_ * rswitch_ ) - rswitch_ );
	Real delta;
	if ( x > ub_ ) {
		delta = x - ub_;
	} else if ( lb_ <= x ) {
		delta = 0;
	} else if ( x < lb_ ) {
		delta = lb_ - x;
	} else {
		delta = 0;
	}
	//  tr.Trace << "evaluate x in [ lb_ ub_ ]: delta " << x << " " << lb_ << " " << ub_ << " " << delta << std::endl;

	delta/=sd_;

	if ( x > ub_ + rswitch_*sd_ ) {
		return 2 * rswitch_ * delta - rswitch_ * rswitch_;
	} else {
		return delta * delta;
	}
}

Real
BoundFunc::dfunc( Real const x ) const {

	Real delta( 0 );
	if ( x > ub_ ) {
		delta = x - ub_;
	} else if ( x < lb_ ) {
		delta = lb_ - x;
	}
	if ( x > ub_ + rswitch_*sd_ ) {
		return 2.0*rswitch_/sd_;
	}
	return 2.0 * (delta) / ( sd_ * sd_ ) * (( x < lb_ ) ? (-1.0) : 1.0);
}


Size
BoundFunc::show_violations( std::ostream& out, Real x, Size verbose_level, Real threshold) const {
	if ( verbose_level > 100 ) {
		out << " " << type_ << " " ;
	}
	if ( verbose_level > 75  ) {
		out << x << " [ " << lb_ << " , " << ub_ << "] ";
		if ( x < lb_ && ( this->func(x) > threshold ) ) out << "TOO SHORT";
		else if ( x > ub_ && ( this->func(x) > threshold ) ) out << "VIOLATED";
		else out << "PASSED";
		if ( verbose_level > 80 ) {
			out << " " << type_ << "\n";
		} else out << "\n";
	} else if ( verbose_level > 100 ) {
		out << x << " in[ " << lb_ << " , " << ub_ << " ] ";
		if ( x < lb_ ) out << (x - lb_) / sd_ << "\n";
		else if ( x > ub_ ) out << (x - ub_) / sd_ << "\n";
		else out << "0.0\n";
	} else if ( verbose_level > 70 ) {
		if ( x < lb_  && ( this->func(x) > threshold ) ) out << "-";
		else if ( x > ub_ && ( this->func(x) > threshold ) ) out << "+";
		else out << ".";
	}


	if ( this->func(x) <= threshold ) return 0;
	return 1;
}

void
BoundFunc::show_definition( std::ostream &out ) const {
	using namespace ObjexxFCL::format;
	std::streamsize const input_precision(out.precision()); // bug #0000005; SML
	out << "BOUNDED " << std::setprecision( 4 ) << RJ(7, lb_) << " " << RJ(7, ub_) << " " << RJ(3,sd_) << " ";
	if ( rswitch_ != 0.5 ) out << RJ(5,rswitch_ ) << " ";
	out << type_;
	out << std::setprecision(input_precision) << "\n";
}

void
BoundFunc::read_data( std::istream& in ) {
	using namespace ObjexxFCL;
	basic::Tracer trInfo( "core.io.constraints", basic::t_info );
	std::string tag;
	in >> lb_ >> ub_ >> sd_ >> tag;
	if ( !in.good() ) {
		in.setstate( std::ios_base::failbit );
	}
	if ( is_float( tag ) ) {
		rswitch_ = float_of( tag );
		in >> type_;
	} else {
		//std::string line;
		//  getline( in, line );
		type_ = tag;//+line;
	}
}


PeriodicBoundFunc::PeriodicBoundFunc(
	Real const lb,
	Real const ub,
	Real sd,
	std::string type,
	Real const periodicity_in
) :
	BoundFunc(
	basic::periodic_range(lb, periodicity_in),
	basic::periodic_range(ub,periodicity_in),
	sd, type
	),
	periodicity_( periodicity_in )
{}

func::FuncOP
PeriodicBoundFunc::clone() const { return func::FuncOP( new PeriodicBoundFunc( *this ) ); }

bool PeriodicBoundFunc::operator == ( Func const & rhs ) const {
	if ( !     same_type_as_me(   rhs ) ) return false;
	if ( ! rhs.same_type_as_me( *this ) ) return false;

	PeriodicBoundFunc const & rhs_pbf( static_cast< PeriodicBoundFunc const & > ( rhs ) );
	if ( periodicity_ != rhs_pbf.periodicity_ ) return false;

	return parent::operator == ( rhs );
}

bool PeriodicBoundFunc::same_type_as_me( Func const & rhs ) const {
	return dynamic_cast< PeriodicBoundFunc const * > (&rhs);
}

void
PeriodicBoundFunc::read_data( std::istream& in )
{
	in >> periodicity_ ;
	parent::read_data( in );
}

Real PeriodicBoundFunc::func(Real const x ) const
{
	return parent::func( basic::periodic_range(x , periodicity_ ) );
}

Real
PeriodicBoundFunc::dfunc( Real const x ) const
{
	return parent::dfunc( basic::periodic_range(x , periodicity_ ) );
}

void
PeriodicBoundFunc::show_definition( std::ostream &out ) const
{
	using namespace ObjexxFCL::format;
	out << "PERIODICITYBOUNDED " << RJ(7, periodicity_) << " ";
	parent::show_definition( out );
}

Size
PeriodicBoundFunc::show_violations( std::ostream& out, Real x, Size verbose_level, Real threshold ) const
{
	return parent::show_violations( out, basic::periodic_range( x, periodicity_ ), verbose_level, threshold );
}

OffsetPeriodicBoundFunc::OffsetPeriodicBoundFunc(
	Real const lb,
	Real const ub,
	Real sd,
	std::string type,
	Real const periodicity_in,
	Real const offset_in
) :
	BoundFunc( lb, ub, sd, type ),
	periodicity_( periodicity_in ),
	offset_( offset_in )
{}

func::FuncOP OffsetPeriodicBoundFunc::clone() const { return func::FuncOP( new OffsetPeriodicBoundFunc( *this ) ); }

bool OffsetPeriodicBoundFunc::operator == ( Func const & rhs ) const {
	if ( !     same_type_as_me(   rhs ) ) return false;
	if ( ! rhs.same_type_as_me( *this ) ) return false;

	OffsetPeriodicBoundFunc const & rhs_opbf( static_cast< OffsetPeriodicBoundFunc const & > ( rhs ) );
	if ( periodicity_ != rhs_opbf.periodicity_ ) return false;
	if ( offset_ != rhs_opbf.offset_ ) return false;

	return parent::operator == ( rhs );
}

bool OffsetPeriodicBoundFunc::same_type_as_me( Func const & rhs ) const
{
	return dynamic_cast< OffsetPeriodicBoundFunc const * > (&rhs);
}


void
OffsetPeriodicBoundFunc::read_data( std::istream& in )
{
	in >> offset_ >> periodicity_ ;
	parent::read_data( in );
}


Real OffsetPeriodicBoundFunc::func(Real const x ) const
{
	return parent::func( basic::periodic_range(x - offset_, periodicity_ ) );
}

Real OffsetPeriodicBoundFunc::dfunc( Real const x ) const
{
	return parent::dfunc( basic::periodic_range(x - offset_, periodicity_ ) );
}

void
OffsetPeriodicBoundFunc::show_definition( std::ostream &out ) const
{
	using namespace ObjexxFCL::format;
	out << "OFFSETPERIODICITYBOUNDED offset" << RJ(7, offset_) << " period " << RJ(7, periodicity_) << " ";
	parent::show_definition( out );
}

Size OffsetPeriodicBoundFunc::show_violations( std::ostream& out, Real x, Size verbose_level, Real threshold ) const
{
	return parent::show_violations( out, basic::periodic_range( x, periodicity_ ), verbose_level, threshold );
}



}
}
} //core

#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
core::scoring::constraints::OffsetPeriodicBoundFunc::OffsetPeriodicBoundFunc() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::constraints::OffsetPeriodicBoundFunc::save( Archive & arc ) const {
	arc( cereal::base_class< BoundFunc >( this ) );
	arc( CEREAL_NVP( periodicity_ ) ); // Real
	arc( CEREAL_NVP( offset_ ) ); // Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::constraints::OffsetPeriodicBoundFunc::load( Archive & arc ) {
	arc( cereal::base_class< BoundFunc >( this ) );
	arc( periodicity_ ); // Real
	arc( offset_ ); // Real
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::constraints::OffsetPeriodicBoundFunc );
CEREAL_REGISTER_TYPE( core::scoring::constraints::OffsetPeriodicBoundFunc )


/// @brief Default constructor required by cereal to deserialize this class
core::scoring::constraints::PeriodicBoundFunc::PeriodicBoundFunc() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::constraints::PeriodicBoundFunc::save( Archive & arc ) const {
	arc( cereal::base_class< BoundFunc >( this ) );
	arc( CEREAL_NVP( periodicity_ ) ); // Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::constraints::PeriodicBoundFunc::load( Archive & arc ) {
	arc( cereal::base_class< BoundFunc >( this ) );
	arc( periodicity_ ); // Real
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::constraints::PeriodicBoundFunc );
CEREAL_REGISTER_TYPE( core::scoring::constraints::PeriodicBoundFunc )


/// @brief Default constructor required by cereal to deserialize this class
core::scoring::constraints::BoundFunc::BoundFunc() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::constraints::BoundFunc::save( Archive & arc ) const {
	arc( cereal::base_class< func::Func >( this ) );
	arc( CEREAL_NVP( lb_ ) ); // Real
	arc( CEREAL_NVP( ub_ ) ); // Real
	arc( CEREAL_NVP( sd_ ) ); // Real
	arc( CEREAL_NVP( rswitch_ ) ); // Real
	arc( CEREAL_NVP( type_ ) ); // std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::constraints::BoundFunc::load( Archive & arc ) {
	arc( cereal::base_class< func::Func >( this ) );
	arc( lb_ ); // Real
	arc( ub_ ); // Real
	arc( sd_ ); // Real
	arc( rswitch_ ); // Real
	arc( type_ ); // std::string
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::constraints::BoundFunc );
CEREAL_REGISTER_TYPE( core::scoring::constraints::BoundFunc )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_constraints_BoundConstraint )
#endif // SERIALIZATION
