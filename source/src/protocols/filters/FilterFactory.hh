// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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

// c++ headers
#include <map>

#include <utility/vector1.hh>


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


class FilterFactory // singletons need not derive from RefCount : public utility::pointer::ReferenceCount
{
public:
	typedef std::map< std::string, FilterCreatorOP > FilterMap;
	typedef utility::tag::Tag Tag;
	typedef utility::tag::TagCOP TagCOP;
	typedef core::pose::Pose Pose;

public:
	virtual ~FilterFactory();

	static
	FilterFactory * get_instance();

	void factory_register( FilterCreatorOP creator );

	/// @brief Create a mover given its identifying string
	FilterOP newFilter( std::string const & );

	///@brief return new Mover by Tag parsing; the identifying string for the Mover is in the Tag
	FilterOP
	newFilter(
		TagCOP const,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		Pose const &
	);

private:
	FilterFactory();

	// Unimplemented -- uncopyable
	FilterFactory( FilterFactory const & );
	FilterFactory const & operator = ( FilterFactory const & );

private:
	static FilterFactory * instance_;

	FilterMap filter_creator_map_;

};

} //namespace moves
} //namespace protocols

#endif
