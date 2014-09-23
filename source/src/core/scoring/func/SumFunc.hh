// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/func/SumFunc.hh
/// @brief Sum of several individual constraints.
/// @author James Thompson, Greg Taylor


#ifndef INCLUDED_core_scoring_func_SumFunc_hh
#define INCLUDED_core_scoring_func_SumFunc_hh

#include <core/scoring/func/Func.hh>

#include <core/types.hh>

#include <utility/pointer/ReferenceCount.hh>

// AUTO-REMOVED #include <basic/Tracer.hh>

#include <utility/vector1_bool.hh>


// C++ Headers

namespace core {
namespace scoring {
namespace func {

	class SumFunc : public Func {
	public:
		// general func
		SumFunc() {}

		//SumFunc operator =( SumFunc const & other ):
		//	Func()
		//{
		//	funcs_ = other.funcs();
		//	return *this;
		//}

		//SumFunc( SumFunc const & src )
		//{
		//	*this = src;
		//}

		FuncOP
		clone() const { return FuncOP( new SumFunc( *this ) ); }

		Real func( Real const x ) const;
		Real dfunc( Real const x ) const;

		void read_data( std::istream& );

		// class-specific accessors
		utility::vector1< FuncOP > funcs() const {
			return funcs_;
		}

		void funcs( utility::vector1< FuncOP > funcs ) {
			funcs_ = funcs;
		}

		void add_func( FuncOP f ) {
			funcs_.push_back( f );
		}

	private:
  	utility::vector1< FuncOP > funcs_;
	};
} // constraints
} // scoring
} // core

#endif
