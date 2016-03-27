// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   binder/context.hpp
/// @brief  Data structures to represent root context and modules
/// @author Sergey Lyskov

#ifndef _INCLUDED_context_hpp_
#define _INCLUDED_context_hpp_

#include <config.hpp>

#include <clang/AST/Decl.h>

#include <llvm/Support/raw_ostream.h>

#include <string>
#include <vector>
#include <set>
#include <unordered_map>


namespace binder {

class Context;

/// structure to hold include information and set of NamedDecl objects that was already queried
// struct Includes
// {
// 	std::set<std::string> &includes;
// 	std::set<clang::NamedDecl const *> stack;
// };

/// Bindings Generator - represent object that can generate binding info for function, class, enum or data variable
class Binder
{
public:
	typedef std::string string;

	virtual ~Binder() {}

	/// Generate string id that uniquly identify C++ binding object. For functions this is function prototype and for classes forward declaration.
	virtual string id() const = 0;

	// return Clang AST NamedDecl pointer to original declaration used to create this Binder
	virtual clang::NamedDecl * named_decl() const = 0;

	/// check if generator can create binding
	virtual bool bindable() const = 0;

	bool binding_requested() const { return binding_requested_; };
	bool skipping_requested() const { return skipping_requested_; };

	/// request bindings for this generator
	void request_bindings() { binding_requested_=true; }

	/// request skipping for this generator
	void request_skipping() { skipping_requested_=true; }

	/// check if user supplied config requested binding for the given declaration and if so request it
	virtual void request_bindings_and_skipping(Config const &) = 0;

	/// extract include needed for this generator and add it to includes vector
	virtual void add_relevant_includes(std::vector<std::string> &includes, std::set<clang::NamedDecl const *> &stack) const = 0;

	/// generate binding code for this object and all its dependencies
	virtual void bind(Context &) = 0;

	// return true if code was already generate for this object
	bool is_binded() const;

	// return binding code
	string & code() { return code_; }
	string const & code() const { return code_; }

	/// return true if object declared in system header
	bool is_in_system_header();

	explicit operator std::string() const { return id(); /*named_decl()->getQualifiedNameAsString();*/ }

private:
	bool binding_requested_=false, skipping_requested_=false;

	string code_;
};


typedef std::shared_ptr< Binder > BinderOP;

typedef std::vector<BinderOP> Binders;



llvm::raw_ostream & operator << (llvm::raw_ostream & os, Binder const &b);

/// Module - represent bindings of individual Python module
// struct Module
// {
// 	std::string path;
// 	std::vector<BinderOP> binders;
// };


/// Context - root, hold bindings info for whole TranslationUnit
class Context
{
	typedef std::string string;

public:

	void add(BinderOP &);

	void generate(Config const &config);

	// find binder related to given type name and bind it
	void request_bindings(std::string const &type);

	/// generate C++ expression for module variable for namespace_
	string module_variable_name(string const & namespace_);

	//void request_bindings(std::vector<string> const & namespaces);
	// is binding requested

private:

	/// array of all binderes from translation unit
	std::vector<BinderOP> binders;

	/// types → binder
	std::unordered_map<string, BinderOP> types;

	std::set<string> ids;

	// binder.id() → binder
	//std::unordered_map<string, Binders> binders_map;

	//std::unordered_map<string, Binders> modules;
	//std::unordered_map<string, BinderOP> system_binders;


	/// bind all objects residing in namespaces and it dependency
	void bind(Config const &config);

	std::set<string> create_all_nested_namespaces();

	/// create vector of all namespaces and sort it
	//std::vector<string> sorted_namespaces();
	//std::vector<string> bind_namespaces(string const &namespace_, size_t max_code_size);
};




} // namespace binder

#endif // _INCLUDED_context_hpp_
