// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   binder/class.hpp
/// @brief  Binding generation for C++ struct and class objects
/// @author Sergey Lyskov

#ifndef _INCLUDED_class_hpp_
#define _INCLUDED_class_hpp_

#include <context.hpp>

#include <clang/AST/DeclCXX.h>

#include <string>
#include <vector>
#include <set>

namespace binder {

// generate classtemplate specialization for ClassTemplateSpecializationDecl or empty string otherwise
std::string template_specialization(clang::CXXRecordDecl const *C);


// generate class name that could be used in bindings code indcluding template specialization if any
std::string class_name(clang::CXXRecordDecl const *C);


// generate qualified class name that could be used in bindings code indcluding template specialization if any
std::string class_qualified_name(clang::CXXRecordDecl const *C);


// generate vector<QualType> with all types that class uses if it tempalated
std::vector<clang::QualType> get_type_dependencies(clang::CXXRecordDecl const *C /*, bool include_members=false*/);


// generate vector<NamedDecl const *> with all declarations related to class
std::vector<clang::NamedDecl const *> get_decl_dependencies(clang::CXXRecordDecl const *C);


/// check if generator can create binding
bool is_bindable(clang::CXXRecordDecl const *C);


// /// check if bindings for particular object was requested
// bool is_binding_requested(clang::CXXRecordDecl const *C, std::vector<std::string> const & namespaces_to_bind);


// // check if bindings for object should be skipped
// bool is_skipping_requested(clang::CXXRecordDecl const *C, std::vector<std::string> const & namespaces_to_skip);


/// check if user requested binding for the given declaration
bool is_binding_requested(clang::CXXRecordDecl const *C, Config const &config);


/// check if user requested skipping for the given declaration
bool is_skipping_requested(clang::CXXRecordDecl const *C, Config const &config);


// extract include needed for declaration and add it to includes
void add_relevant_includes(clang::CXXRecordDecl const *C, std::vector<std::string> &includes, std::set<clang::NamedDecl const *> &stack, int level);


class ClassBinder : public Binder
{
public:
	ClassBinder(clang::CXXRecordDecl *c) : C(c) {}

	/// Generate string id that uniquly identify C++ binding object. For functions this is function prototype and for classes forward declaration.
	string id() const override;

	// return Clang AST NamedDecl pointer to original declaration used to create this Binder
	clang::NamedDecl * named_decl() const override { return C; };

	/// check if generator can create binding
    bool bindable() const override;

	/// check if user requested binding for the given declaration
	virtual void request_bindings_and_skipping(Config const &) override;

	/// extract include needed for this generator and add it to includes vector
	void add_relevant_includes(std::vector<std::string> &includes, std::set<clang::NamedDecl const *> &stack) const override;

	/// generate binding code for this object and all its dependencies
	void bind(Context &) override;


private:
	clang::CXXRecordDecl *C;

	// set of members which user asked to exclude from bindings
	//std::set<clang::NamedDecl const *> members_to_skip;
};



} // namespace binder

#endif // _INCLUDED_class_hpp_
