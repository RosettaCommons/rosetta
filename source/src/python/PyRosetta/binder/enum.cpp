// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   binder/enum.cpp
/// @brief  Binding generation for C++ enums
/// @author Sergey Lyskov


#include <enum.hpp>

#include <type.hpp>
#include <util.hpp>

#include <fmt/format.h>

using namespace llvm;
using namespace clang;

using std::string;
using std::vector;

using namespace fmt::literals;

namespace binder {


/// extract include needed for this generator and add it to includes vector
void add_relevant_includes(clang::EnumDecl const *E, std::vector<std::string> &includes, std::set<clang::NamedDecl const *> &stack)
{
	add_relevant_include_for_decl(E, includes);
}


// Generate binding for given function: py::enum_<MyEnum>(module, "MyEnum")...
std::string bind_enum(std::string const & module, EnumDecl *E)
{
	string name { E->getNameAsString() };
	string qualified_name { E->getQualifiedNameAsString() };

	string r = "\tpybind11::enum_<{}>({}, \"{}\")\n"_format(qualified_name, module, name);

	// for(auto d = E->decls_begin(); d != E->decls_end(); ++d) {
	// 	if(EnumConstantDecl *e = dyn_cast<EnumConstantDecl>(*d) ) {
	// 		//outs() << "EnumConstant: " << e->getQualifiedNameAsString() << "\n";
	// 		r += "\t\t.value(\"{}\", {})\n"_format(e->getNameAsString(), e->getQualifiedNameAsString());
	// 	}
	// }

	for(auto e = E->enumerator_begin(); e != E->enumerator_end(); ++e) {
		//outs() << "EnumConstant: " << e->getQualifiedNameAsString() << "\n";
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
void EnumBinder::add_relevant_includes(std::vector<std::string> &includes, std::set<clang::NamedDecl const *> &stack) const
{
	binder::add_relevant_includes(E, includes, stack);
}

/// generate binding code for this object and all its dependencies
void EnumBinder::bind(Context &context)
{
	if( is_binded() ) return;

	string const module_variable_name = context.module_variable_name( namespace_from_named_decl(E) );
	string const include = relevant_include(E);

	code()  = "\t// " + E->getQualifiedNameAsString() + " file:" + include.substr(1, include.size()-2) + " line:" + line_number(E) + "\n";
	code() += bind_enum(module_variable_name, E) + ";\n\n";
}



} // namespace binder
