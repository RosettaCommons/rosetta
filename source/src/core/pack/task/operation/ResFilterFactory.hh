// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/operation/ResFilterFactory.hh
/// @brief
/// @author ashworth

#ifndef INCLUDED_core_pack_task_operation_ResFilterFactory_hh
#define INCLUDED_core_pack_task_operation_ResFilterFactory_hh

// Unit Headers
#include <core/pack/task/operation/ResFilterFactory.fwd.hh>

// Package Headers
#include <core/pack/task/operation/ResFilter.fwd.hh>
#include <core/pack/task/operation/ResFilterCreator.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// c++ headers
#include <string>
#include <map>

#include <utility/vector0.hh>

#ifdef MULTI_THREADED
#include <atomic>
#include <mutex>
#endif

namespace core {
namespace pack {
namespace task {
namespace operation {

class ResFilterFactory // singleton; no need to derive from RefCount : public utility::pointer::ReferenceCount
{
public:
	//typedef utility::pointer::ReferenceCount parent;
	typedef std::map< std::string, ResFilterCreatorOP > ResFilterCreatorMap;
	typedef utility::tag::Tag Tag;
	typedef utility::tag::TagOP TagOP;
	typedef utility::tag::TagCOP TagCOP;

public:
	static ResFilterFactory * get_instance();
	void factory_register( ResFilterCreatorOP );

	/// @brief add a prototype, using its default type name as the map key
	void add_creator( ResFilterCreatorOP );
	bool has_type( std::string const & ) const;

	/// @brief return new ResFilter by key lookup in filter_map_; new ResFilter does not parse an input Tag
	ResFilterOP newResFilter( std::string const & ) const;

	/// @brief return new ResFilter by key lookup in filter_map_ new ResFilter parses Tag if provided
	ResFilterOP newResFilter( std::string const &, TagCOP tag ) const;

	/// @brief The %ResFilterFactory is the point of entry for the definition of the XML Schemas
	/// for every ResFilter that may be instantiated from a file.  It is  responsible for defining
	/// an xs:group named "res_filter" listing each of the residue-filter-complex types that
	/// may be initialized using the %ResFilterFactory and to iterate across each of the
	/// ResFilterCreator it contains asking them for the XML schema of the ResFilter they
	/// are responsible for creating.
	void define_res_filter_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const;

	/// @brief The name given to the XML schema group of all ResFilter s.
	static std::string res_filter_xml_schema_group_name();

#ifdef MULTI_THREADED
public:

	/// @brief This public method is meant to be used only by the
	/// utility::thread::safely_create_singleton function and not meant
	/// for any other purpose.  Do not use.
	static std::mutex & singleton_mutex();

private:
	static std::mutex singleton_mutex_;
#endif

private:
	ResFilterFactory();
	virtual ~ResFilterFactory();
	/// @brief private singleton creation function to be used with
	/// utility::thread::threadsafe_singleton
	static ResFilterFactory * create_singleton_instance();

private:
#if defined MULTI_THREADED
	static std::atomic< ResFilterFactory * > instance_;
#else
	static ResFilterFactory * instance_;
#endif

	ResFilterCreatorMap filter_creator_map_;

};

} //namespace operation
} //namespace task
} //namespace pack
} //namespace core

#endif
