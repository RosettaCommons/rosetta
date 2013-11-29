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


#ifndef INCLUDED_core_scoring_constraints_ScalarWeightedFunc_hh
#define INCLUDED_core_scoring_constraints_ScalarWeightedFunc_hh

#include <core/scoring/func/ScalarWeightedFunc.fwd.hh>

#include <core/scoring/func/Func.hh>

#include <core/types.hh>

#include <utility/pointer/ReferenceCount.hh>


// C++ Headers

namespace core {
namespace scoring {
namespace constraints {

	class ScalarWeightedFunc : public Func {
	public:
		ScalarWeightedFunc(
			Real const weight,
			FuncOP myfunc
		): weight_( weight ), func_to_weight_( myfunc ) {}

		virtual ~ScalarWeightedFunc();

		virtual
		FuncOP
		clone() const { return new ScalarWeightedFunc( *this ); }

		virtual
		Real func( Real const x ) const;
		virtual
		Real dfunc( Real const x ) const;

		virtual
		void read_data( std::istream& );

		virtual
		void show_definition( std::ostream &out ) const;
		virtual
		Size show_violations(std::ostream &out, Real x, Size verbose_level, Real threshold = 1) const;


	private:
		Real weight_;
		FuncOP func_to_weight_;
	};
} // constraints
} // scoring
} // core

#endif
