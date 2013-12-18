// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/loops/loops_definers/LoopsDefinerFactory.hh
/// @brief Factory for creating LoopsDefiner objects
/// @author Matthew O'Meara (mattjomeara@gmail.com)


#ifndef INCLUDED_protocols_loops_loops_definers_LoopsDefinerFactory_hh
#define INCLUDED_protocols_loops_loops_definers_LoopsDefinerFactory_hh

// Unit Headers
#include <protocols/loops/loops_definers/LoopsDefinerFactory.fwd.hh>
#include <protocols/loops/loops_definers/LoopsDefinerCreator.fwd.hh>

// Project Headers
#include <protocols/loops/loops_definers/LoopsDefiner.fwd.hh>
#include <protocols/loops/Loops.fwd.hh>

// Platform Headers
#include <core/pose/Pose.fwd.hh>

// Utility Headers
#include <utility/factory/WidgetRegistrator.hh>
#include <utility/vector1.hh>

// C++ Headers
#include <map>

#ifdef MULTI_THREADED
#ifdef CXX11
// C++11 Headers
#include <thread>
#endif
#endif

namespace protocols {
namespace loops {
namespace loops_definers {


/// Create LoopsDefiner Reporters
class LoopsDefinerFactory {

public: // types

	typedef std::map< std::string, LoopsDefinerCreatorCOP > LoopsDefinerCreatorMap;

private: // constructors

	// Private constructor to make it singleton managed
	LoopsDefinerFactory();
	LoopsDefinerFactory(const LoopsDefinerFactory & src); // unimplemented

	LoopsDefinerFactory const &
	operator=( LoopsDefinerFactory const & ); // unimplemented

	/// @brief private singleton creation function to be used with
	/// utility::thread::threadsafe_singleton
	static LoopsDefinerFactory * create_singleton_instance();

public:

	// Warning this is not called because of the singleton pattern
	virtual
	~LoopsDefinerFactory();

	static
	LoopsDefinerFactory *
	get_instance();

	void
	factory_register(
		LoopsDefinerCreatorOP);

	bool
	has_type(
		std::string const & ) const;

	utility::vector1< std::string >
	get_all_loops_definer_names() const;

	LoopsDefinerOP
	create_loops_definer(
		std::string const & type_name);

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

	static LoopsDefinerFactory * instance_;
	LoopsDefinerCreatorMap types_;
};

/// @brief This templated class will register an instance of an
/// LoopsDefinerCreator (class T) with the
/// LoopsDefinerFactory.  It will ensure that no
/// LoopsDefinerCreator is registered twice, and, centralizes this
/// registration logic so that thread safety issues can be handled in
/// one place
template < class T >
class LoopsDefinerRegistrator : public utility::factory::WidgetRegistrator< LoopsDefinerFactory, T >
{
public:
	typedef utility::factory::WidgetRegistrator< LoopsDefinerFactory, T > parent;
public:
	LoopsDefinerRegistrator() : parent() {}
};



} // namespace
} // namespace
} // namespace

#endif
