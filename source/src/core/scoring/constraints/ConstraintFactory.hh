// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/scoring/constraints/ConstraintFactory.hh
/// @brief Factory for creating various types of constraints.
/// @author Greg Taylor <gktaylor@u.washington.edu>

#ifndef INCLUDED_core_scoring_constraints_ConstraintFactory_hh
#define INCLUDED_core_scoring_constraints_ConstraintFactory_hh

// Unit headers

// Package headers
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/ConstraintCreator.fwd.hh>

// Utility headers
#include <utility/factory/WidgetRegistrator.hh>

// C++ Headers
#include <map>

#ifdef MULTI_THREADED
#ifdef CXX11
// C++11 Headers
#include <thread>
#endif
#endif

namespace core {
namespace scoring {
namespace constraints {

class ConstraintFactory {
private:
	ConstraintFactory();

	ConstraintFactory(ConstraintFactory const &); // unimplemented
	ConstraintFactory const & operator=( ConstraintFactory const & ); // unimplemented

	/// @brief private singleton creation function to be used with
	/// utility::thread::threadsafe_singleton
	static ConstraintFactory * create_singleton_instance();

public:
	static ConstraintFactory * get_instance();

	//void add_type( std::string const & type_name, ConstraintOP new_cst );
	//void add_type( ConstraintOP new_cst );

	void factory_register( ConstraintCreatorCOP creator );
	scoring::constraints::ConstraintOP newConstraint( std::string const & type_name );
	utility::vector1< std::string > get_cst_names() const;

	/// @brief Replace the load-time ConstraintCreator with another creator.
	void replace_creator( ConstraintCreatorCOP creator );

	ConstraintCreatorCOP
	get_creator( std::string const & type_name );

#ifdef MULTI_THREADED
#ifdef CXX11
public:

	/// @brief This public method is meant to be used only by the
	/// utility::thread::safely_create_singleton function and not meant
	/// for any other purpose.  Do not use.
	static std::mutex & singleton_mutex();

private:
	static std::mutex singleton_mutex_;
#endif
#endif

private:

	static ConstraintFactory * instance_;

	typedef std::map< std::string, scoring::constraints::ConstraintCreatorCOP > ConstraintCreatorMap;
	ConstraintCreatorMap cst_types_;
};

/// @brief This templated class will register an instance of an
/// ConstraintCreator (class T) with the ConstraintFactory.  It will ensure
/// that no ConstraintCreator is registered twice, and, centralizes
/// this registration logic so that thread safety issues can be handled in
/// one place
template < class T >
class ConstraintRegistrator : public utility::factory::WidgetRegistrator< ConstraintFactory, T >
{
public:
	typedef utility::factory::WidgetRegistrator< ConstraintFactory, T > parent;
public:
	ConstraintRegistrator() : parent() {}
};


} //constraints
} //scoring
} //core

#endif
