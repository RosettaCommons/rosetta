// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   binder/test/T10.inheritance.hpp
/// @brief  Binder self-test file. Bindings tests for class hierarchy’s.
/// @author Sergey Lyskov

#ifndef _INCLUDED_T10_inheritance_hpp_
#define _INCLUDED_T10_inheritance_hpp_

#include <memory>
#include <string>

class Base
{
public:

	void foo() {}

	void maybe() {}
};


class Delete : public Base
{
public:
	Delete() = delete;

	void maybe() = delete;
};

class A : public std::enable_shared_from_this<A>
{};

class B : public std::string
{};

#endif // _INCLUDED_T10_inheritance_hpp_
