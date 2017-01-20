// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// Copyright (c) 2016 Sergey Lyskov <sergey.lyskov@jhu.edu>
//
// All rights reserved. Use of this source code is governed by a
// MIT license that can be found in the LICENSE file.

/// @file   binder/test/T10.inheritance.hpp
/// @brief  Binder self-test file. Bindings tests for class hierarchy’s.
/// @author Sergey Lyskov

#ifndef _INCLUDED_T10_inheritance_hpp_
#define _INCLUDED_T10_inheritance_hpp_

#include <memory>
#include <string>

class Base
{
protected:
	virtual void foo_protected() {}
	int data;

public:
	void foo() {}

	void maybe() {}
};

class Derived : public Base
{
public:
	using Base::foo_protected;
	using Base::data;
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
