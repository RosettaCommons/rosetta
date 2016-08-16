// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
