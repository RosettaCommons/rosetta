// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/parser/DataLoaderFactory.hh
/// @brief  Factory class for the parser's DataLoader classes
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_jd2_parser_DataLoaderFactory_hh
#define INCLUDED_protocols_jd2_parser_DataLoaderFactory_hh

// Package Headers
#include <protocols/jd2/parser/DataLoader.fwd.hh>
#include <protocols/jd2/parser/DataLoaderCreator.hh>

// Utility Headers
#include <utility/SingletonBase.hh>
#include <utility/factory/WidgetRegistrator.hh>
#include <utility/pointer/ReferenceCount.hh>

// c++ headers
#include <map>

#ifdef MULTI_THREADED
#ifdef CXX11
// C++11 Headers
#include <mutex>
#include <atomic>
#endif
#endif

namespace protocols {
namespace jd2 {
namespace parser {

/// @brief A factory for creating DataLoaders, which are able to load
/// arbitrary data into the basic::datacache::DataMap used in the XML-based parser.
/// This factory supports the load-time registration scheme allowing
/// DataLoaders to be defined in libraries outside of protocols.lib
class DataLoaderFactory
{
public:
	typedef std::map< std::string, DataLoaderCreatorOP > LoaderMap;

public:
	virtual ~DataLoaderFactory();

	static
	DataLoaderFactory * get_instance();

	void factory_register( DataLoaderCreatorOP creator );

	/// @brief Create a DataLoader given its identifying string
	DataLoaderOP newDataLoader( std::string const & ) const;

private:
	DataLoaderFactory();

	// Unimplemented -- uncopyable
	DataLoaderFactory( DataLoaderFactory const & );
	DataLoaderFactory const & operator = ( DataLoaderFactory const & );

	/// @brief private singleton creation function to be used with
	/// utility::thread::threadsafe_singleton
	static DataLoaderFactory * create_singleton_instance();

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
	/// Singleton instance pointer
#if defined MULTI_THREADED && defined CXX11
	static std::atomic< DataLoaderFactory * > instance_;
#else
	static DataLoaderFactory * instance_;
#endif

	LoaderMap dataloader_creator_map_;

};

/// @brief This templated class will register an instance of an
/// DataLoaderCreator (class T) with the DataLoaderFactory.  It will ensure
/// that no DataLoaderCreator is registered twice, and, centralizes
/// this registration logic so that thread safety issues can be handled in
/// one place
template < class T >
class DataLoaderRegistrator : public utility::factory::WidgetRegistrator< DataLoaderFactory, T >
{
public:
	typedef utility::factory::WidgetRegistrator< DataLoaderFactory, T > parent;
public:
	DataLoaderRegistrator() : parent() {}
};

} //namespace parser
} //namespace jd2
} //namespace protocols

#endif
