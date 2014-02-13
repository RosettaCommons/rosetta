// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/datacache/WriteableCacheableDataFactory.hh
/// @brief
/// @author Justin Porter

#ifndef INCLUDED_basic_datacache_WriteableCacheableDataFactory_hh
#define INCLUDED_basic_datacache_WriteableCacheableDataFactory_hh

// Unit Headers
#include <basic/datacache/WriteableCacheableDataFactory.fwd.hh>
#include <basic/datacache/WriteableCacheableDataCreator.hh>
#include <basic/datacache/WriteableCacheableData.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>

// Utility Headers
#include <utility/factory/WidgetRegistrator.hh>

// c++ headers
#include <map>
#include <set>
#include <utility/vector1.hh>

#include <boost/utility.hpp>

#ifdef MULTI_THREADED
#ifdef CXX11
// C++11 Headers
#include <thread>
#endif
#endif

namespace basic {
namespace datacache {

/// @brief This templated class will register an instance of an
/// WriteableCacheableDataCreator (class T) with the WriteableCacheableDataFactory.  It will ensure
/// that no WriteableCacheableDataCreator is registered twice, and, centralizes
/// this registration logic so that thread safety issues can be handled in
/// one place
template < class T >
class WriteableCacheableDataRegistrator : public utility::factory::WidgetRegistrator< WriteableCacheableDataFactory, T >
{
public:
  typedef utility::factory::WidgetRegistrator< WriteableCacheableDataFactory, T > parent;
public:
  WriteableCacheableDataRegistrator() : parent() {}
};

class WriteableCacheableDataFactory : boost::noncopyable { // Singleton: should not derive from ReferenceCount

public:
  typedef std::map< std::string, WriteableCacheableDataCreatorOP > WriteableCacheableDataMap;

public:
  virtual ~WriteableCacheableDataFactory() {};

  static
  WriteableCacheableDataFactory * get_instance();

  void factory_register( WriteableCacheableDataCreatorOP creator );

  /// @brief Create a data instance given its identifying string
  WriteableCacheableDataOP new_data_instance( std::string const& data_type_name,
                                              std::istream &in );

private:
  WriteableCacheableDataFactory() {};

  /// @brief private singleton creation function to be used with
  /// utility::thread::threadsafe_singleton
  static WriteableCacheableDataFactory* create_singleton_instance();

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
  static WriteableCacheableDataFactory* instance_;

  WriteableCacheableDataMap data_creator_map_;

};

} //namespace datacache
} //namespace basic

#endif
