// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   binder/test/T10.override.hpp
/// @brief  Binder self-test file. Bindings tests for overriding classes in Python.
/// @author Sergey Lyskov

#ifndef _INCLUDED_T10_override_hpp_
#define _INCLUDED_T10_override_hpp_

#include <memory>
#include <string>

class Base
{
public:
	std::string s_;

	virtual void foo_void() {}
	virtual void foo_void_pure() = 0;

	virtual int         foo_int() { return 0; }
	virtual std::string foo_sring() { return ""; }

	virtual std::string &foo_sring_r() { return s_; }
};

#endif // _INCLUDED_T10_override_hpp_
