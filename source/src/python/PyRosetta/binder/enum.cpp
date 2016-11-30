// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// Copyright (c) 2016 Sergey Lyskov <sergey.lyskov@jhu.edu>
//
// All rights reserved. Use of this source code is governed by a
// MIT license that can be found in the LICENSE file.

/// @file   binder/enum.cpp
/// @brief  Binding generation for C++ enums
/// @author Sergey Lyskov


#include <enum.hpp>

#include <type.hpp>
#include <util.hpp>

#include <fmt/format.h>

#include <clang/AST/ASTContext.h>

using namespace llvm;
using namespace clang;

using std::string;
using std::vector;

using namespace fmt::literals;

namespace binder {


/// extract include needed for this generator and add it to includes vector
void add_relevant_includes(clang::EnumDecl const *E, IncludeSet &includes, int level)
{
	if( !includes.add_decl(E, level) ) return;
	add_relevant_include_for_decl(E, includes);
}


// Generate binding for given function: py::enum_<MyEnum>(module, "MyEnum")...
std::string bind_enum(std::string const & module, EnumDecl *E)
{
	string name { E->getNameAsString() };
	string qualified_name { E->getQualifiedNameAsString() };

	string r = "\tpybind11::enum_<{}>({}, \"{}\", \"{}\")\n"_format(qualified_name, module, name, generate_documentation_string_for_declaration(E));

	for(auto e = E->enumerator_begin(); e != E->enumerator_end(); ++e) {
		r += "\t\t.value(\"{}\", {})\n"_format(e->getNameAsString(), e->getQualifiedNameAsString());
	}
	r.pop_back();

	return r + ( E->isScopedUsingClassTag() ? ";\n\n" : "\n\t\t.export_values();\n\n" );
}


/// Generate string id that uniquly identify C++ binding object. For functions this is function prototype and for classes forward declaration.
string EnumBinder::id() const
{
	return E->getQualifiedNameAsString();
}


/// check if generator can create binding
bool EnumBinder::bindable() const
{
	if( E->isCXXInstanceMember()  or  E->isCXXClassMember() ) return false;
	else return true;
}


/// check if user requested binding for the given declaration
void EnumBinder::request_bindings_and_skipping(Config const &config)
{
	if( config.is_namespace_binding_requested( namespace_from_named_decl(E) ) ) Binder::request_bindings();
}


/// extract include needed for this generator and add it to includes vector
void EnumBinder::add_relevant_includes(IncludeSet &includes) const
{
	binder::add_relevant_includes(E, includes, 0);
}

/// generate binding code for this object and all its dependencies
void EnumBinder::bind(Context &context)
{
	if( is_binded() ) return;

	string const module_variable_name = context.module_variable_name( namespace_from_named_decl(E) );

	code()  = "\t" + generate_comment_for_declaration(E);
	code() += bind_enum(module_variable_name, E) + ";\n\n";
}



} // namespace binder
