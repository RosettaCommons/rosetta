// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/func/ScalarWeightedFunc.hh
/// @brief Weighted constraint function that encapsulates other constraints
/// @author James Thompson, Greg Taylor
/// @details The ScalarWeightedFunc is a class used to scale other
/// constraints by multiplying it by a constant value.
/// This is useful in cases where constraint templates are similar
/// to each other and should be down-weighted to avoid double counting
/// for the same structural information.


#ifndef INCLUDED_core_scoring_func_SkipViolFunc_hh
#define INCLUDED_core_scoring_func_SkipViolFunc_hh

#include <core/scoring/func/Func.hh>
#include <core/types.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace func {

class SkipViolFunc : public Func {
public:
	SkipViolFunc(
		Real const weight,
		FuncOP myfunc
	): weight_( weight ),
		count_viols_( 0 ),
		func_to_weight_( myfunc ) {}

	~SkipViolFunc() {};

	FuncOP
	clone() const { return FuncOP( new SkipViolFunc( *this ) ); }
	virtual bool operator == ( Func const & other ) const;
	virtual bool same_type_as_me( Func const & other ) const;

	Real func( Real const x ) const;
	Real dfunc( Real const x ) const;

	void read_data( std::istream& );

	/// @brief show some sort of stringified representation of the violations for this constraint.
	virtual core::Size show_violations( std::ostream& out, Real r, core::Size verbose_level, Real threshold = 1 ) const;

	/// @brief shows the definition of this function, usually the string type of function and the
	/// parameters passed in to the constructor.
	virtual void show_definition( std::ostream & out ) const;

	Size viols() const {
		return count_viols_;
	}

private:
	Real weight_;
	core::Size count_viols_;
	FuncOP func_to_weight_;

	//typedef std::map< std::string, scoring::func::FuncOP > FuncTypes;
	// static FuncTypes func_types_;

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	SkipViolFunc();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};
} // constraints
} // scoring
} // core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_func_SkipViolFunc )
#endif // SERIALIZATION


#endif
