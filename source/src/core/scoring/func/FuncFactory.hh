// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/func/FuncFactory.hh
/// @brief Factory for creating various types of constraints.
/// @author Greg Taylor <gktaylor@u.washington.edu>

#ifndef INCLUDED_core_scoring_func_FuncFactory_hh
#define INCLUDED_core_scoring_func_FuncFactory_hh

// Unit headers
#include <core/scoring/func/FuncFactory.fwd.hh>

// Package headers
#include <core/scoring/func/Func.fwd.hh>

// Project headers

// C++ Headers
#include <string>
#include <map>

namespace core {
namespace scoring {
namespace func {

class FuncFactory {
public:
	typedef std::map< std::string, scoring::func::FuncOP > FuncTypes;

public:
	FuncFactory();
	~FuncFactory();
	void add_type( std::string type_name, FuncOP new_func );
	FuncOP new_func( std::string const& type ) const;

private:
	FuncTypes func_types_;
};

} //constraints
} //scoring
} //core

#endif
