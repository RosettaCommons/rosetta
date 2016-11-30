// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// Copyright (c) 2016 Sergey Lyskov <sergey.lyskov@jhu.edu>
//
// All rights reserved. Use of this source code is governed by a
// MIT license that can be found in the LICENSE file.

/// @file   binder/test/T02.function.hpp
/// @brief  Binder self-test file. Bindings for functions.
/// @author Sergey Lyskov

#ifndef _INCLUDED_T02_function_hpp_
#define _INCLUDED_T02_function_hpp_

#include <memory>

class A {};

void foo(bool)   {}
void foo(int)    {}
void foo(long)   {}
void foo(float)  {}
void foo(double) {}
void foo(A)      {}
void foo(std::shared_ptr<A> ) {}

void foo_r(bool &)   {}
void foo_r(int &)    {}
void foo_r(long &)   {}
void foo_r(float &)  {}
void foo_r(double &) {}
void foo_r(A &)      {}
void foo_r(std::shared_ptr<A> &) {}
void foo_r_not_binded(std::pair<int, int> &) {}

void foo_cr(bool const &)   {}
void foo_cr(int const &)    {}
void foo_cr(long const &)   {}
void foo_cr(float const &)  {}
void foo_cr(double const &) {}
void foo_cr(A const &)      {}
void foo_cr(std::shared_ptr<A> const &) {}
void foo_cr2(std::shared_ptr<A const> const &) {}
void foo_cr(std::pair<int, int> const &) {}

void foo_p(bool *)   {}
void foo_p(int *)    {}
void foo_p(long *)   {}
void foo_p(float *)  {}
void foo_p(double *) {}
void foo_p(A *)               {}
void foo_p(std::shared_ptr<A> *) {}
void foo_p_not_binded(std::pair<int, int> *) {}


#endif // _INCLUDED_T02_function_hpp_
