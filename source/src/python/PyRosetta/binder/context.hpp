// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// Copyright (c) 2016 Sergey Lyskov <sergey.lyskov@jhu.edu>
//
// All rights reserved. Use of this source code is governed by a
// MIT license that can be found in the LICENSE file.

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
#include <unordered_set>
#include <unordered_map>

namespace binder {

class Context;

// structure to hold include set information and set of NamedDecl objects that was already queried for includes
class IncludeSet
{
	std::vector<std::string> includes_;

	std::unordered_map<clang::NamedDecl const *, int> stack_;

public:
	// add include to the set
	void add_include(std::string const &i) { includes_.push_back(i); }

	// check if declaration is already in stack with level at lease as 'level' or lower and add it if it is not - return true if declaration was added
	bool add_decl(clang::NamedDecl const *, int level);

	std::vector<std::string> const &includes() const { return includes_; }

	// remove all includes and clear up the stack
	void clear();
};


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
	virtual void add_relevant_includes(IncludeSet &) const = 0;

	/// generate binding code for this object and all its dependencies
	virtual void bind(Context &) = 0;

	// return true if code was already generate for this object
	bool is_binded() const;

	// return binding code
	string & code() { return code_; }
	string const & code() const { return code_; }

	/// return true if object declared in system header
	bool is_in_system_header();

	// return vector of declarations that need to be binded before this one could
	virtual std::vector<clang::CXXRecordDecl const *> dependencies() const { return std::vector<clang::CXXRecordDecl const *>(); }

	/// return prefix portion of bindings code
	virtual string prefix_code() const { return string(); }

	/// return unique strting ID for this binder
	explicit operator std::string() const { return id(); /*named_decl()->getQualifiedNameAsString();*/ }

private:
	bool binding_requested_=false, skipping_requested_=false;

	string code_;
};


typedef std::shared_ptr< Binder > BinderOP;

typedef std::vector<BinderOP> Binders;


llvm::raw_ostream & operator << (llvm::raw_ostream & os, Binder const &b);


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

	void add_insertion_operator(clang::FunctionDecl const *f);

	/// find gloval insertion operator for given type, return nullptr if not such operator find
	clang::FunctionDecl const * global_insertion_operator(clang::CXXRecordDecl const *);

	/// generate unique trace line containing `info` to insert into the code
	std::string trace_line(std::string const &info);

private:

	/// bind all objects residing in namespaces and it dependency
	void bind(Config const &config);

	std::set<string> create_all_nested_namespaces();

		/// check if forward declaration for CXXRecordDecl needed
	bool is_forward_needed(clang::CXXRecordDecl const *);

	/// add given class to 'aleady binded' set
	void add_to_binded(clang::CXXRecordDecl const *);

	/// sort vector of binders by dependecy so python imports could work
	void sort_binders();


	/// map of function-ptr -> Decl* for global instertion operators so we can determent for which types can bind __repr__
	std::map<std::string, clang::FunctionDecl const *> insertion_operators;

	/// array of all binderes from translation unit
	std::vector<BinderOP> binders;

	/// types → binder
	std::unordered_map<string, BinderOP> types;

	/// set of items unique id's to keep track of whats binders being added
	std::unordered_set<string> ids;

	/// set of items unique id's to keep track of whats binded and not
	std::set<string> binded;

	/// counter to generate unique trace lines for debug
	int trace_counter = -1;
};


} // namespace binder

#endif // _INCLUDED_context_hpp_
