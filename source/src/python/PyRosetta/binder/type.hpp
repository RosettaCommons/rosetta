// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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


/*
template<typename R, R default_value>
class QualType_functor
{
public:
	R apply(clang::QualType const &qt) {
		if( clang::PointerType const *pt = clang::dyn_cast<clang::PointerType>( qt.getTypePtr() ) ) return apply( pt->getPointeeType() );

		if( clang::ReferenceType const *rt = clang::dyn_cast<clang::ReferenceType>( qt.getTypePtr() ) ) return apply( rt->getPointeeType() );

		if( clang::CXXRecordDecl const *C = qt->getAsCXXRecordDecl() ) return cxx_record_decl(C);

		if( clang::EnumDecl const *E = clang::dyn_cast_or_null<clang::EnumDecl>( qt->getAsTagDecl() ) ) return enum_decl(E);

		return default_value;
	}

	R cxx_record_decl(clang::CXXRecordDecl const *C) { return default_value; }
	R enum_decl(clang::EnumDecl const *E) { return default_value; }
};


class QualType_void_functor
{
public:
	void apply(clang::QualType const &qt) {
		if( clang::PointerType const *pt = clang::dyn_cast<clang::PointerType>( qt.getTypePtr() ) ) apply( pt->getPointeeType() );

		if( clang::ReferenceType const *rt = clang::dyn_cast<clang::ReferenceType>( qt.getTypePtr() ) ) apply( rt->getPointeeType() );

		if( clang::CXXRecordDecl const *C = qt->getAsCXXRecordDecl() ) cxx_record_decl(C);

		if( clang::EnumDecl const *E = clang::dyn_cast_or_null<clang::EnumDecl>( qt->getAsTagDecl() ) ) enum_decl(E);
	}

	void cxx_record_decl(clang::CXXRecordDecl const *C) {}
	void enum_decl(clang::EnumDecl const *E) {}
};


/// check if user requested binding for the given QualType
class f_is_binding_requested : public QualType_functor<bool, false>
{
	Config const &config;
public:
	f_is_binding_requested(Config const &config_) : config(config_) {}

	bool cxx_record_decl(clang::CXXRecordDecl const *C);
};


/// check if user requested skipping for the given QualType
class f_is_skipping_requested : public QualType_functor<bool, false>
{
	Config const &config;
public:
	f_is_skipping_requested(Config const &config_) : config(config_) {}

	bool cxx_record_decl(clang::CXXRecordDecl const *C);
};



/// extract include needed for this generator and add it to includes vector
class f_add_relevant_includes : public QualType_void_functor
{
	std::vector<std::string> &includes;
	std::set<clang::NamedDecl const *> &stack;

public:
	f_add_relevant_includes(std::vector<std::string> &includes_, std::set<clang::NamedDecl const *> &stack_) : includes(includes_), stack(stack_) {}

	void cxx_record_decl(clang::CXXRecordDecl const *C);
	void enum_decl(clang::EnumDecl const *E);
};
*/

} // namespace binder

#endif  // _INCLUDED_type_hpp_
