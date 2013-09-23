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

// Project Headers
// AUTO-REMOVED #include <core/scoring/constraints/ConstraintIO.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

// Utility headers

static basic::Tracer tr("core.constraints.BoundFunc",basic::t_info);

// C++ headers

namespace core {
namespace scoring {
namespace constraints {



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


Real delta;
if ( x > ub_ ) {
delta = x - ub_;
} else if ( lb_ <= x ) {
delta = 0;
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
	if (verbose_level > 100) {
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
	} else if (verbose_level > 70 ) {
		if ( x < lb_  && ( this->func(x) > threshold ) ) out << "-";
		else if ( x > ub_ && ( this->func(x) > threshold )) out << "+";
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
	basic::Tracer trInfo("core.io.constraints", basic::t_info );
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
		//		getline( in, line );
		type_ = tag;//+line;
	}
}


void
PeriodicBoundFunc::show_definition( std::ostream &out ) const
{
	using namespace ObjexxFCL::format;
	out << "PERIODICITYBOUNDED " << RJ(7, periodicity_) << " ";
	parent::show_definition( out );
}

void
PeriodicBoundFunc::read_data( std::istream& in )
{
	in >> periodicity_ ;
	parent::read_data( in );
}


void
OffsetPeriodicBoundFunc::show_definition( std::ostream &out ) const
{
	using namespace ObjexxFCL::format;
	out << "OFFSETPERIODICITYBOUNDED offset" << RJ(7, offset_) << " period " << RJ(7, periodicity_) << " ";
	parent::show_definition( out );
}

void
OffsetPeriodicBoundFunc::read_data( std::istream& in )
{
	in >> offset_ >> periodicity_ ;
	parent::read_data( in );
}


}
}
} //core
