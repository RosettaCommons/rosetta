// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// Copyright (c) 2016 Sergey Lyskov <sergey.lyskov@jhu.edu>
//
// All rights reserved. Use of this source code is governed by a
// MIT license that can be found in the LICENSE file.

/// @file   binder/function.hpp
/// @brief  Binding generation for static and member functions
/// @author Sergey Lyskov

#ifndef _INCLUDED_function_hpp_
#define _INCLUDED_function_hpp_

#include <context.hpp>

#include <clang/AST/Decl.h>
#include <clang/AST/DeclCXX.h>

#include <string>
#include <tuple>

namespace binder {


// Generate function argument list separate by comma
std::string function_arguments(clang::FunctionDecl const *record);

// Generate function argument list separate by comma
// name_arguments - if arguments should be named: a1, a2, ...
// n - number of arguments to generate. If n > num_of_function_parameters - generate only list with num_of_function_parameters
std::pair<std::string, std::string> function_arguments_for_lambda(clang::FunctionDecl const *record, uint n);


// Generate function argument list with types separate by comma and with only arguments names
// name_arguments - if arguments should be named: a1, a2, ...
std::tuple<std::string, std::string, std::string> function_arguments_for_py_overload(clang::FunctionDecl const *record);


// generate string represetiong class name that could be used in python
std::string python_function_name(clang::FunctionDecl const *F);


// Generate function pointer type string for given function. Example void (*)(int, doule)_ or  void (ClassName::*)(int, doule)_ for memeber function
std::string function_pointer_type(clang::FunctionDecl const *record);


// generate qualified function name that could be used in bindings code indcluding template specialization if any
std::string function_qualified_name(clang::FunctionDecl const *F);

/// check if user requested binding for the given declaration
bool is_binding_requested(clang::FunctionDecl const *F, Config const &config);


/// check if user requested skipping for the given declaration
bool is_skipping_requested(clang::FunctionDecl const *F, Config const &config);


// Generate binding for given function: .def("foo", (std::string (aaaa::A::*)(int) ) &aaaa::A::foo, "doc")
// If function have default arguments generate set of bindings by creating separate bindings for each argument with default.
// if parent is not nullptr then bind function as-if it a member of that CXXRecordDecl (for handling visibility changes with 'using' directive)
std::string bind_function(std::string const & module, clang::FunctionDecl const *F, Context &, clang::CXXRecordDecl const *parent = nullptr);


/// extract include needed for this generator and add it to includes vector
void add_relevant_includes(clang::FunctionDecl const *F, IncludeSet &includes, int level/*, bool for_template_arg_only=false*/);


/// check if generator can create binding
bool is_bindable(clang::FunctionDecl const *F);

/// check if methods could be overload in Python
bool is_overloadable(clang::CXXMethodDecl const *M);


class FunctionBinder : public Binder
{
public:
	FunctionBinder(clang::FunctionDecl *f) : F(f) {}

	/// Generate string id that uniquly identify C++ binding object. For functions this is function prototype and for classes forward declaration.
	string id() const override;

	// return Clang AST NamedDecl pointer to original declaration used to create this Binder
	clang::NamedDecl * named_decl() const override { return F; };

	/// check if generator can create binding
	bool bindable() const override;

	/// check if user requested binding for the given declaration
	void request_bindings_and_skipping(Config const &) override;

	/// extract include needed for this generator and add it to includes vector
	void add_relevant_includes(IncludeSet &includes) const override;

	/// generate binding code for this object and all its dependencies
	void bind(Context &) override;

private:
	clang::FunctionDecl *F;
};



} // namespace binder

#endif // _INCLUDED_function_hpp_
