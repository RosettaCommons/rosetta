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
/// @author Oliver Lange

#ifndef INCLUDED_core_scoring_constraints_BoundConstraint_hh
#define INCLUDED_core_scoring_constraints_BoundConstraint_hh

#include <core/scoring/constraints/BoundConstraint.fwd.hh>
//#include <core/scoring/constraints/Constraint.hh>
//#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/func/Func.hh>
#include <core/types.hh>
#include <basic/basic.hh>


// C++ Headers
//#include <cstdlib>
#include <iostream>

//#include <map>
//#include <utility>


namespace core {
namespace scoring {
namespace constraints {


class BoundFunc : public func::Func {
public:
	BoundFunc( Real const lb, Real const ub, Real sd, std::string type ): lb_( lb ), ub_( ub ), sd_ ( sd ), rswitch_( 0.5 ), type_( type ) {}
	BoundFunc( Real const lb, Real const ub, Real sd, Real rswitch, std::string type )
	  : lb_( lb ), ub_( ub ), sd_ ( sd ), rswitch_( rswitch ), type_( type ) {}

	virtual
	func::FuncOP clone() const { return func::FuncOP( new BoundFunc( *this ) ); };

	virtual
	void read_data( std::istream& );

	virtual
	Real
	func( Real const x ) const;

	virtual
	Real
	dfunc( Real const x ) const;

	virtual
	void show_definition( std::ostream &out ) const;

	virtual
	Size show_violations( std::ostream& out, Real x, Size verbose_level, Real threshold = 1 ) const;

	Real lb() const { return lb_; }
	Real ub() const { return ub_; }
	Real sd() const { return sd_; }
	Real rswitch() const {return rswitch_;}
	std::string type() const {return type_;}

private:
	Real lb_;
	Real ub_;
	Real sd_;
	Real rswitch_;
	std::string type_;
};


/// a variant of the bound func that is periodic
class PeriodicBoundFunc : public BoundFunc
{
public:
	typedef BoundFunc parent;

public:
	PeriodicBoundFunc(
		Real const lb, Real const ub, Real sd, std::string type, Real const periodicity_in
	) :
	BoundFunc(
		basic::periodic_range(lb, periodicity_in),
		basic::periodic_range(ub,periodicity_in),
		sd, type
	),
	periodicity_( periodicity_in )
	{}

	func::FuncOP clone() const { return func::FuncOP( new PeriodicBoundFunc( *this ) ); };

	void read_data( std::istream& );

	Real func(Real const x ) const
	{
		return parent::func( basic::periodic_range(x , periodicity_ ) );
	}

  Real dfunc( Real const x ) const
	{
		return parent::dfunc( basic::periodic_range(x , periodicity_ ) );
	}

	void show_definition( std::ostream& out ) const;

	virtual Size show_violations( std::ostream& out, Real x, Size verbose_level, Real threshold = 1) const
	{
		return parent::show_violations( out, basic::periodic_range( x, periodicity_ ), verbose_level, threshold );
	}

private:
	Real periodicity_;

};


/// a variant of the bound func that is periodic
class OffsetPeriodicBoundFunc : public BoundFunc
{
public:
	typedef BoundFunc parent;

public:
	OffsetPeriodicBoundFunc(
		Real const lb, Real const ub, Real sd, std::string type, Real const periodicity_in, Real const offset_in
	) :
	BoundFunc(
		lb,ub,
		sd, type
	),
	periodicity_( periodicity_in ),
	offset_( offset_in )
	{}

	func::FuncOP clone() const { return func::FuncOP( new OffsetPeriodicBoundFunc( *this ) ); };

	void read_data( std::istream& );

	Real func(Real const x ) const
	{
		return parent::func( basic::periodic_range(x - offset_, periodicity_ ) );
	}

	Real dfunc( Real const x ) const
	{
		return parent::dfunc( basic::periodic_range(x - offset_, periodicity_ ) );
	}

	void show_definition( std::ostream& out ) const;

	virtual Size show_violations( std::ostream& out, Real x, Size verbose_level, Real threshold = 1) const
	{
		return parent::show_violations( out, basic::periodic_range( x, periodicity_ ), verbose_level, threshold );
	}


private:
	Real periodicity_;
	Real offset_;

};


} //constraints
} //scoring
} //core

#endif
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Oliver Lange

#ifndef INCLUDED_core_scoring_constraints_BoundConstraint_hh
#define INCLUDED_core_scoring_constraints_BoundConstraint_hh

#include <core/scoring/constraints/BoundConstraint.fwd.hh>
//#include <core/scoring/constraints/Constraint.hh>
//#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/func/Func.hh>
#include <core/types.hh>
#include <basic/basic.hh>


// C++ Headers
//#include <cstdlib>
#include <iostream>

//#include <map>
//#include <utility>


namespace core {
namespace scoring {
namespace constraints {


class BoundFunc : public func::Func {
public:
 BoundFunc( Real const lb, Real const ub, Real sd, std::string type ): lb_( lb ), ub_( ub ), sd_ ( sd ), rswitch_( 0.5 ), type_( type ) {}
 BoundFunc( Real const lb, Real const ub, Real sd, Real rswitch, std::string type )
 : lb_( lb ), ub_( ub ), sd_ ( sd ), rswitch_( rswitch ), type_( type ) {}

 virtual
 func::FuncOP clone() const { return new BoundFunc( *this ); };

 virtual
 void read_data( std::istream& );

 virtual
 Real
 func( Real const x ) const;

 virtual
 Real
 dfunc( Real const x ) const;

 virtual
 void show_definition( std::ostream &out ) const;

 virtual
 Size show_violations( std::ostream& out, Real x, Size verbose_level, Real threshold = 1 ) const;

 Real lb() const { return lb_; }
 Real ub() const { return ub_; }
 Real sd() const { return sd_; }

private:
 Real lb_;
 Real ub_;
 Real sd_;
 Real rswitch_;
 std::string type_;
};


/// a variant of the bound func that is periodic
class PeriodicBoundFunc : public BoundFunc
{
public:
 typedef BoundFunc parent;

public:
 PeriodicBoundFunc(
 Real const lb, Real const ub, Real sd, std::string type, Real const periodicity_in
 ) :
 BoundFunc(
 basic::periodic_range(lb, periodicity_in),
 basic::periodic_range(ub,periodicity_in),
 sd, type
 ),
 periodicity_( periodicity_in )
 {}

 func::FuncOP clone() const { return new PeriodicBoundFunc( *this ); };

 void read_data( std::istream& );

 Real func(Real const x ) const
 {
 return parent::func( basic::periodic_range(x , periodicity_ ) );
 }

 Real dfunc( Real const x ) const
 {
 return parent::dfunc( basic::periodic_range(x , periodicity_ ) );
 }

 void show_definition( std::ostream& out ) const;

 virtual Size show_violations( std::ostream& out, Real x, Size verbose_level, Real threshold = 1) const
 {
 return parent::show_violations( out, basic::periodic_range( x, periodicity_ ), verbose_level, threshold );
 }

private:
 Real periodicity_;

};


/// a variant of the bound func that is periodic
class OffsetPeriodicBoundFunc : public BoundFunc
{
public:
 typedef BoundFunc parent;

public:
 OffsetPeriodicBoundFunc(
 Real const lb, Real const ub, Real sd, std::string type, Real const periodicity_in, Real const offset_in
 ) :
 BoundFunc(
 lb,ub,
 sd, type
 ),
 periodicity_( periodicity_in ),
 offset_( offset_in )
 {}

 func::FuncOP clone() const { return new OffsetPeriodicBoundFunc( *this ); };

 void read_data( std::istream& );

 Real func(Real const x ) const
 {
 return parent::func( basic::periodic_range(x - offset_, periodicity_ ) );
 }

 Real dfunc( Real const x ) const
 {
 return parent::dfunc( basic::periodic_range(x - offset_, periodicity_ ) );
 }

 void show_definition( std::ostream& out ) const;

 virtual Size show_violations( std::ostream& out, Real x, Size verbose_level, Real threshold = 1) const
 {
 return parent::show_violations( out, basic::periodic_range( x, periodicity_ ), verbose_level, threshold );
 }

private:
 Real periodicity_;
 Real offset_;

};


} //constraints
} //scoring
} //core

#endif
