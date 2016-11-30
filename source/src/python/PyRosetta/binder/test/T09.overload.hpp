// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// Copyright (c) 2016 Sergey Lyskov <sergey.lyskov@jhu.edu>
//
// All rights reserved. Use of this source code is governed by a
// MIT license that can be found in the LICENSE file.

/// @file   binder/test/T09.overload.hpp
/// @brief  Binder self-test file. Test for function overloads in Python
/// @author Sergey Lyskov

#pragma once

struct A1
{
	A1(int a=1) {}
};

struct A2
{
	A2(int a=1) {}

	virtual void f() {}
};

struct A3
{
	A3(int a=1) {}
	virtual void f() = 0;
};

struct A4
{
	A4(int a=1) {}
	virtual int operator=(int) = 0;
};
