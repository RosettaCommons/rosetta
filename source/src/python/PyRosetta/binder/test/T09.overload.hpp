// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
