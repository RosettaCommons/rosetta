// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   binder/enum.hpp
/// @brief  Binding generation for C++ enums
/// @author Sergey Lyskov


#ifndef _INCLUDED_enum_hpp_
#define _INCLUDED_enum_hpp_

#include <context.hpp>

#include <clang/AST/Decl.h>

#include <string>

namespace binder {

/// extract include needed for this generator and add it to includes vector
void add_relevant_includes(clang::EnumDecl const *E, IncludeSet &includes, int level);


// Generate binding for given function: py::enum_<MyEnum>(module, "MyEnum")...
std::string bind_enum(std::string const & module, clang::EnumDecl *E);


class EnumBinder : public Binder
{
public:
	EnumBinder(clang::EnumDecl *e) : E(e) {}

	/// Generate string id that uniquly identify C++ binding object. For functions this is function prototype and for classes forward declaration.
	string id() const override;

	// return Clang AST NamedDecl pointer to original declaration used to create this Binder
	clang::NamedDecl * named_decl() const override { return E; };

	/// check if generator can create binding
    bool bindable() const override;

	/// check if user requested binding for the given declaration
	virtual void request_bindings_and_skipping(Config const &) override;

	/// extract include needed for this generator and add it to includes vector
	void add_relevant_includes(IncludeSet &includes) const override;

	/// generate binding code for this object and all its dependencies
	void bind(Context &) override;


private:
	clang::EnumDecl *E;
};


} // namespace binder

#endif // _INCLUDED_enum_hpp_
