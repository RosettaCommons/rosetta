// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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

void foo_not_binded(bool *)   {}
void foo_not_binded(int *)    {}
void foo_not_binded(long *)   {}
void foo_not_binded(float *)  {}
void foo_not_binded(double *) {}
void foo_p(A *)               {}
void foo_p(std::shared_ptr<A> *) {}


#endif // _INCLUDED_T02_function_hpp_
