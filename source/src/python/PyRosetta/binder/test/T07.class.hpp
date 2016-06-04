// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   binder/test/T07.class.hpp
/// @brief  Binder self-test file. Bindings of classes data members functionality.
/// @author Sergey Lyskov

#ifndef _INCLUDED_T07_class_hpp_
#define _INCLUDED_T07_class_hpp_

struct Integral
{
	int i;
	unsigned int ui;
	long l;
	unsigned long ul;
};

struct Floating_point
{
	float f;
	double d;
	long double ld;
};

struct Arrays
{
	int   i__not_binded[8];
	int   ii_not_binded[8][8];
	float f__not_binded[8];
	float ff_not_binded[8][8];
};

struct Pointers_and_References
{
	Pointers_and_References() = delete;

	int   * ip_not_binded;
	float * fp_not_binded;
	int   & i_r;
	float & f_r;
};


class Enums
{
public:
	enum E1 {E1V=2};
	enum class E2 {E2V=4};
	enum class E3 {};

protected:
	enum EP0_not_binded {EP0V};

private:
	enum EP1_not_binded {EP1V};
};


class Access
{
public:
	int i;
	void foo_public() {};

	static float foo(double) { return 1.0; }

	static double a_double;

protected:
	int i_protected_not_binded;
	void foo_protected_not_binded() {};

private:
	int i_not_binded;
	void foo_private_not_binded() {};
};


struct ConstOverload
{
	void foo(int) {}
	void foo(int) const {}
};

#endif // _INCLUDED_T07_class_hpp_
