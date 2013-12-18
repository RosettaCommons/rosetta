// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/sequence/SequenceFactory.hh
/// @brief Factory for creating various types of sequence.
/// @author James Thompson

#ifndef INCLUDED_core_sequence_SequenceFactory_hh
#define INCLUDED_core_sequence_SequenceFactory_hh

// AUTO-REMOVED #include <core/pose/Pose.fwd.hh>
#include <core/sequence/Sequence.fwd.hh>
#include <core/sequence/SequenceCreator.fwd.hh>

#include <utility/vector1.hh>
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
namespace sequence {

class SequenceFactory {
private:
	SequenceFactory();

	SequenceFactory(SequenceFactory const &); // unimplemented
	SequenceFactory const & operator=( SequenceFactory const & ); // unimplemented
public:
	static SequenceFactory * get_instance();

	void factory_register( SequenceCreatorCOP creator );
	SequenceOP get_sequence( std::string const & type_name );
	SequenceOP seq_from_file(
		std::string const & fn,
		std::string const & type_name
	);
	utility::vector1< std::string > get_seq_names() const;

	/// @brief Replace the load-time SequenceCreator with another creator.
	void replace_creator( SequenceCreatorCOP creator );

	SequenceCreatorCOP
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

	/// @brief private singleton creation function to be used with
	/// utility::thread::threadsafe_singleton
	static SequenceFactory * create_singleton_instance();

private:

	static SequenceFactory * instance_;

	typedef std::map< std::string, SequenceCreatorCOP > SequenceCreatorMap;
	SequenceCreatorMap seq_types_;
};

/// @brief This templated class will register an instance of an
/// SequenceCreator (class T) with the SequenceFactory.  It will ensure
/// that no SequenceCreator is registered twice, and, centralizes
/// this registration logic so that thread safety issues can be handled in
/// one place
template < class T >
class SequenceRegistrator : public utility::factory::WidgetRegistrator< SequenceFactory, T >
{
public:
	typedef utility::factory::WidgetRegistrator< SequenceFactory, T > parent;
public:
	SequenceRegistrator() : parent() {}
};


} //sequence
} //core

#endif
