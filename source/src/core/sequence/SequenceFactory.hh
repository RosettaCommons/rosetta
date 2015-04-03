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

#include <core/sequence/Sequence.fwd.hh>
#include <core/sequence/SequenceCreator.fwd.hh>

// Utility headers
#include <utility/SingletonBase.hh>
#include <utility/vector1.hh>
#include <utility/factory/WidgetRegistrator.hh>

// C++ Headers
#include <map>

namespace core {
namespace sequence {

class SequenceFactory : public utility::SingletonBase< SequenceFactory >
{
public:
	friend class utility::SingletonBase< SequenceFactory >;

private:
	SequenceFactory();

	SequenceFactory(SequenceFactory const &); // unimplemented
	SequenceFactory const & operator=( SequenceFactory const & ); // unimplemented

public:
	void factory_register( SequenceCreatorCOP creator );
	SequenceOP get_sequence( std::string const & type_name );
	SequenceOP seq_from_file(
		std::string const & fn,
		std::string const & type_name
	);
	utility::vector1< std::string > get_seq_names() const;

	SequenceCreatorCOP
	get_creator( std::string const & type_name );

private:

	/// @brief private singleton creation function to be used with
	/// utility::thread::threadsafe_singleton
	static SequenceFactory * create_singleton_instance();

private:

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
