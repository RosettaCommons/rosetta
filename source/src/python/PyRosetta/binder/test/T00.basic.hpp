// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// Copyright (c) 2016 Sergey Lyskov <sergey.lyskov@jhu.edu>
//
// All rights reserved. Use of this source code is governed by a
// MIT license that can be found in the LICENSE file.

/// @file   binder/test/T00.basic.hpp
/// @brief  Binder self-test file. Basic functionality.
/// @author Sergey Lyskov

#ifndef _INCLUDED_T00_basic_hpp_
#define _INCLUDED_T00_basic_hpp_

int global_int;
long global_long;

unsigned int global_unsigned_int;
unsigned long global_unsigned_long;

float global_float;
double global_double;

void foo() {}

int foo_int(int)   { return 0; }
int foo_int(int &) { return 0; }
int foo_int(int *) { return 0; }

float foo_float(float)   { return 0; }
float foo_float(float &) { return 0; }
float foo_float(float *) { return 0; }

double foo_double(double)   { return 0; }
double foo_double(double &) { return 0; }
double foo_double(double *) { return 0; }

using VoidP = void *;

VoidP foo_void_p(VoidP) { return nullptr; }
VoidP foo_void_p(VoidP const &) { return nullptr; }
VoidP foo_void_p(VoidP &) { return nullptr; }

#endif // _INCLUDED_T00_basic_hpp_
