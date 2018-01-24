// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/filters/FilterFactory.hh
/// @brief
/// @author ashworth

#ifndef INCLUDED_protocols_filters_FilterFactory_hh
#define INCLUDED_protocols_filters_FilterFactory_hh

// Unit Headers
#include <protocols/filters/FilterFactory.fwd.hh>
#include <protocols/filters/FilterCreator.hh>
#include <protocols/filters/Filter.fwd.hh>

#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
//#include <protocols/filters/Filter.fwd.hh> // Filters_map (typedef)

// Project headers
#include <core/pose/Pose.fwd.hh>

// Utility Headers
#include <utility/factory/WidgetRegistrator.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/SingletonBase.hh>

// c++ headers
#include <map>

#include <utility/vector1.hh>

#ifdef MULTI_THREADED
#include <atomic>
#include <mutex>
#endif

namespace protocols {
namespace filters {

/// @brief This templated class will register an instance of an
/// FilterCreator (class T) with the FilterFactory.  It will ensure
/// that no FilterCreator is registered twice, and, centralizes
/// this registration logic so that thread safety issues can be handled in
/// one place
template < class T >
class FilterRegistrator : public utility::factory::WidgetRegistrator< FilterFactory, T >
{
public:
	typedef utility::factory::WidgetRegistrator< FilterFactory, T > parent;
public:
	FilterRegistrator() : parent() {}
};


class FilterFactory : public utility::SingletonBase< FilterFactory >
{
public:
	friend class utility::SingletonBase< FilterFactory >;

	typedef std::map< std::string, FilterCreatorOP > FilterMap;
	typedef utility::tag::Tag Tag;
	typedef utility::tag::TagCOP TagCOP;
	typedef core::pose::Pose Pose;

public:
	virtual ~FilterFactory();

	void factory_register( FilterCreatorOP creator );

	/// @brief Is there a filter with the given name that's known to Rosetta?
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	bool filter_exists( std::string const &filter_name ) const;

	/// @brief Get the XML schema for a given filter.
	/// @details Throws an error if the filter is unknown to Rosetta.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void provide_xml_schema( std::string const &filter_name, utility::tag::XMLSchemaDefinition & xsd ) const;

	/// @brief Create a filter given its identifying string
	FilterOP newFilter( std::string const & ) const;

	/// @brief return new Filter by Tag parsing; the identifying string for the Filter is in the Tag
	FilterOP
	newFilter(
		TagCOP,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		Pose const &,
		std::string user_defined_name = ""
	) const;

	/// @brief Read access to the set of all FilterCreators; for unit testing purposes
	FilterMap const & filter_creator_map() const;

	void define_filter_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const;
	static std::string filter_xml_schema_group_name();

private:
	FilterFactory();

	// Unimplemented -- uncopyable
	FilterFactory( FilterFactory const & ) = delete;
	FilterFactory const & operator = ( FilterFactory const & ) = delete;

private:

	FilterMap filter_creator_map_;

};

} //namespace moves
} //namespace protocols

#endif
