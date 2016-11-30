// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// Copyright (c) 2016 Sergey Lyskov <sergey.lyskov@jhu.edu>
//
// All rights reserved. Use of this source code is governed by a
// MIT license that can be found in the LICENSE file.

/// @file   binder/test/T30.include.hpp
/// @brief  Binder self-test file. Tests for include-request functionality.
/// @author Sergey Lyskov

#ifndef _INCLUDED_T30_include_hpp_
#define _INCLUDED_T30_include_hpp_

#include <T30.include.incl.a.include>
#include <T30.include.incl.b.include>

void foo(bbbb::B<aaaa::A>) {}

#endif // _INCLUDED_T30_include_hpp_
