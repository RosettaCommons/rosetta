// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// Copyright (c) 2016 Sergey Lyskov <sergey.lyskov@jhu.edu>
//
// All rights reserved. Use of this source code is governed by a
// MIT license that can be found in the LICENSE file.

/// @file   binder/test/T40.stl.hpp
/// @brief  Binder self-test file. Tests for bindings of STL classes.
/// @author Sergey Lyskov

#ifndef _INCLUDED_T40_stl_hpp_
#define _INCLUDED_T40_stl_hpp_

#include <vector>

void foo(std::vector<int>) {}

template class std::vector<float>;

#endif // _INCLUDED_T40_stl_hpp_
