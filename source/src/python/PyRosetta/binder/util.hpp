// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   binder/util.hpp
/// @brief  Various helper functions
/// @author Sergey Lyskov


#ifndef _INCLUDED_util_hpp_
#define _INCLUDED_util_hpp_

#include <context.hpp>

#include <clang/AST/Decl.h>
#include <clang/AST/DeclTemplate.h>

#include <string>
#include <vector>

namespace binder {


/// Split string using given separator
std::vector<std::string> split(std::string const &buffer, std::string const & separator="\n");


/// Replace all occurrences of string inplace
void replace(std::string &r, std::string const & from, std::string const &to);

/// Replace all occurrences of string and return result as new string
std::string replace_(std::string const &s, std::string const & from, std::string const &to);

/// check if string begins with given prefix
bool begins_with(std::string const &source, std::string const &prefix);

/// check if string ends with given prefix
bool ends_with(std::string const &source, std::string const &prefix);

/// indent given code
std::string indent(std::string const &code, std::string const &indentation);


/// Try to read exisitng file and if content does not match to code - write a new version. Also create nested dirs starting from prefix if nessesary.
void update_source_file(std::string const &prefix, std::string const &file_name, std::string const &code);



/// calculate namespace path from given NamedDecl, like: std, core::pose
std::string namespace_from_named_decl(clang::NamedDecl const *decl);


std::string typename_from_type_decl(clang::TypeDecl *decl);


/// Calculate base (upper) namespace for given one: core::pose::motif --> core::pose
std::string base_namespace(std::string const & ns);


/// Calculate last namespace for given one: core::pose::motif --> motif
std::string last_namespace(std::string const & ns);

/// replace all _Bool types with bool
void fix_boolean_types(std::string &type);

/// Generate string representation of given expression
std::string expresion_to_string(clang::Expr *e);

/// Generate string representation of given TemplateArgument
std::string template_argument_to_string(clang::TemplateArgument const &);


/// calcualte line in source file for NamedDecl
std::string line_number(clang::NamedDecl const *decl);

// generate string represetiong class name that could be used in python
std::string mangle_type_name(std::string const &name, bool mark_template=true);

// generate C++ comment line for given declartion along with file path and line number: // core::scoring::vdwaals::VDWAtom file:core/scoring/vdwaals/VDWTrie.hh line:43
std::string generate_comment_for_declaration(clang::NamedDecl const *decl);

// extract doc string (Doxygen comments) for given declaration and convert it to C++ source code string
std::string generate_documentation_string_for_declaration(clang::Decl const*);


} // namespace binder

#endif // _INCLUDED_util_hpp_
