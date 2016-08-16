// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   binder/binder.hpp
/// @brief  Options
/// @author Sergey Lyskov


#ifndef _INCLUDED_binder_hpp_
#define _INCLUDED_binder_hpp_


#include "llvm/Support/CommandLine.h"

extern llvm::cl::opt<bool> O_annotate_includes;
extern llvm::cl::opt<bool> O_single_file;

namespace binder {
} // namespace binder

#endif // _INCLUDED_binder_hpp_
