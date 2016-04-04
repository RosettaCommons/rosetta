// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   binder/function.hpp
/// @brief  Binding generation for static and member functions
/// @author Sergey Lyskov

#ifndef _INCLUDED_function_hpp_
#define _INCLUDED_function_hpp_

#include <context.hpp>

#include <clang/AST/Decl.h>

#include <string>

namespace binder {


// Generate function argument list separate by comma
std::string function_arguments(clang::FunctionDecl const *record);


// Generate function pointer type string for given function. Example void (*)(int, doule)_ or  void (ClassName::*)(int, doule)_ for memeber function
std::string function_pointer_type(clang::FunctionDecl const *record);


// generate qualified function name that could be used in bindings code indcluding template specialization if any
std::string function_qualified_name(clang::FunctionDecl const *F);

/// check if user requested binding for the given declaration
bool is_binding_requested(clang::FunctionDecl const *F, Config const &config);


/// check if user requested skipping for the given declaration
bool is_skipping_requested(clang::FunctionDecl const *F, Config const &config);


// Generate binding for given function: .def("foo", (std::string (aaaa::A::*)(int) ) &aaaa::A::foo, "doc")
std::string bind_function(clang::FunctionDecl *F, Context &);


/// extract include needed for this generator and add it to includes vector
void add_relevant_includes(clang::FunctionDecl const *F, std::vector<std::string> &includes, std::set<clang::NamedDecl const *> &stack, int level/*, bool for_template_arg_only=false*/);


/// check if generator can create binding
bool is_bindable(clang::FunctionDecl const *F);



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
	void add_relevant_includes(std::vector<std::string> &includes, std::set<clang::NamedDecl const *> &stack) const override;

	/// generate binding code for this object and all its dependencies
	void bind(Context &) override;

private:
	clang::FunctionDecl *F;
};



} // namespace binder

#endif // _INCLUDED_function_hpp_
