// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// Copyright (c) 2016 Sergey Lyskov <sergey.lyskov@jhu.edu>
//
// All rights reserved. Use of this source code is governed by a
// MIT license that can be found in the LICENSE file.

/// @file   binder/type.hpp
/// @brief  Various functionality for handline clang::QualType and clang::NamedDecl
/// @author Sergey Lyskov

#ifndef _INCLUDED_type_hpp_
#define _INCLUDED_type_hpp_

#include <context.hpp>
#include <config.hpp>

#include <clang/AST/Decl.h>

#include <string>

namespace binder {

/// check if user requested binding for the given QualType
bool is_binding_requested(clang::QualType const &qt, Config const &config);


/// check if user requested skipping for the given QualType
bool is_skipping_requested(clang::QualType const &qt, Config const &config);


// extract include path needed for declaration itself (without template dependency if any), return empty string if no include could be found
std::string relevant_include(clang::NamedDecl const *decl);


// extract include path needed for declaration itself (without template dependency if any), do nothing if include could not be found (ie for build-in's)
void add_relevant_include_for_decl(clang::NamedDecl const *decl, IncludeSet &includes);


/// extract include needed for this generator and add it to includes vector
void add_relevant_includes(clang::QualType const &qt, /*clang::ASTContext const &context,*/ IncludeSet &includes, int level);


// check if given QualType is bindable
bool is_bindable(clang::QualType const &qt);


/// extract type info from QualType if any and bind relative type
void request_bindings(clang::QualType const &qt, Context &context);


// transform give type name to standard form
std::string standard_name(std::string const &type);


/// check if given class/struct is builtin in Python and therefor should not be binded
bool is_python_builtin(clang::NamedDecl const *C);


} // namespace binder

#endif  // _INCLUDED_type_hpp_
