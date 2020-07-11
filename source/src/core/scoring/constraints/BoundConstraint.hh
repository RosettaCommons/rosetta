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
#include <string>
#include <iosfwd>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace constraints {


class BoundFunc : public func::Func {
public:
	BoundFunc( Real const lb, Real const ub, Real sd, std::string const & type );
	BoundFunc( Real const lb, Real const ub, Real sd, Real rswitch, std::string const & type );

	func::FuncOP clone() const override;

	bool operator==( Func const & rhs ) const override;
	bool same_type_as_me( Func const & other ) const override;

	void read_data( std::istream& ) override;

	Real
	func( Real const x ) const override;

	Real
	dfunc( Real const x ) const override;

	void show_definition( std::ostream &out ) const override;

	Size show_violations( std::ostream& out, Real x, Size verbose_level, Real threshold = 1 ) const override;

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
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	BoundFunc();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


/// a variant of the bound func that is periodic
class PeriodicBoundFunc : public BoundFunc
{
public:
	typedef BoundFunc parent;

public:
	PeriodicBoundFunc(
		Real const lb, Real const ub, Real sd, std::string type, Real const periodicity_in
	);

	func::FuncOP clone() const override;

	bool operator == ( Func const & rhs ) const override;
	bool same_type_as_me( Func const & rhs ) const override;

	void read_data( std::istream& ) override;

	Real func(Real const x ) const override;
	Real dfunc( Real const x ) const override;

	void show_definition( std::ostream& out ) const override;
	Size show_violations( std::ostream& out, Real x, Size verbose_level, Real threshold = 1) const override;

private:
	Real periodicity_;

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	PeriodicBoundFunc();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


/// a variant of the bound func that is periodic
class OffsetPeriodicBoundFunc : public BoundFunc
{
public:
	typedef BoundFunc parent;

public:
	OffsetPeriodicBoundFunc(
		Real const lb, Real const ub, Real sd, std::string type, Real const periodicity_in, Real const offset_in
	);

	func::FuncOP clone() const override;

	bool operator == ( Func const & rhs ) const override;
	bool same_type_as_me( Func const & rhs ) const override;

	void read_data( std::istream& ) override;

	Real func(Real const x ) const override;
	Real dfunc( Real const x ) const override;

	void show_definition( std::ostream& out ) const override;
	Size show_violations( std::ostream& out, Real x, Size verbose_level, Real threshold = 1) const override;

private:
	Real periodicity_;
	Real offset_;

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	OffsetPeriodicBoundFunc();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} //constraints
} //scoring
} //core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_constraints_BoundConstraint )
#endif // SERIALIZATION


#endif
