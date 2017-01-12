// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// Copyright (c) 2016 Sergey Lyskov <sergey.lyskov@jhu.edu>
//
// All rights reserved. Use of this source code is governed by a
// MIT license that can be found in the LICENSE file.

/// @file   binder/binder.hpp
/// @brief  Options
/// @author Sergey Lyskov


#ifndef _INCLUDED_binder_hpp_
#define _INCLUDED_binder_hpp_


#include "llvm/Support/CommandLine.h"

extern llvm::cl::opt<bool> O_annotate_includes;
extern llvm::cl::opt<bool> O_single_file;
extern llvm::cl::opt<bool> O_trace;
extern llvm::cl::opt<bool> O_verbose;

namespace binder {
} // namespace binder

#endif // _INCLUDED_binder_hpp_
