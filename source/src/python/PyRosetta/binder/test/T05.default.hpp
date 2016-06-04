// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   binder/test/T05.default.hpp
/// @brief  Binder self-test file. Bindings of functions with default arguments.
/// @author Sergey Lyskov

#ifndef _INCLUDED_T05_default_hpp_
#define _INCLUDED_T05_default_hpp_

int _static_int_ = 2;
double foo(int a, int i1=1, int &i2=_static_int_, float f=2.0f, double const &d=4.16, void *p=nullptr, int q=10) { return i1 + i2 + f + d; }

struct A {};

A _static_a_;

struct B
{
	B(int const a=9, float b=8, A& a2=_static_a_, A* a3=nullptr, int i=0) {}

	void foo(A a1=A(), A& a2=_static_a_, A* a3=nullptr, int i=0) {}

	double foo_member(int const i=1, float f=2.0f, double d=4.16) { return i + f + d; }

	double foo_member_const(int i=1, float f=2.0f, double d=4.16) const { return i + f + d; }

	static double foo_static(int i=1, float f=2.0f, double d=4.16) { return i + f + d; }
};


#endif // _INCLUDED_T05_default_hpp_
