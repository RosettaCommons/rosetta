// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// Copyright (c) 2016 Sergey Lyskov <sergey.lyskov@jhu.edu>
//
// All rights reserved. Use of this source code is governed by a
// MIT license that can be found in the LICENSE file.

/// @file   binder/test/T08.constructor.hpp
/// @brief  Binder self-test file. Bindings of classes constructors.
/// @author Sergey Lyskov

#ifndef _INCLUDED_T08_constructor_hpp_
#define _INCLUDED_T08_constructor_hpp_

class BaseA {};
class DerivedA : public BaseA {};


class BaseB
{
public:
	BaseB(int) {}
};
class DerivedB : public BaseB {};


class BaseC
{
public:
	BaseC(int) {}
private:
	BaseC() {}
};
class DerivedC : public BaseC {};


class BaseD
{
public:
	BaseD() = delete;
	BaseD(int) {}
};
class DerivedD : public BaseD {};

#endif // _INCLUDED_T08_constructor_hpp_
