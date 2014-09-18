// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/qsar/scoring_grid/GridFactory.hh
/// @author Sam DeLuca

#ifndef INCLUDED_protocols_qsar_scoring_grid_GridFactory_HH
#define INCLUDED_protocols_qsar_scoring_grid_GridFactory_HH

//Unit Headers
#include <protocols/qsar/scoring_grid/GridFactory.fwd.hh>
#include <protocols/qsar/scoring_grid/GridCreator.hh>
#include <utility/json_spirit/json_spirit_reader.h>

// Utility Headers
#include <utility/factory/WidgetRegistrator.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>
// c++ headers
#include <map>

#ifdef MULTI_THREADED
#ifdef CXX11
// C++11 Headers
#include <atomic>
#include <mutex>
#endif
#endif

namespace protocols {
namespace qsar {
namespace scoring_grid {

/// @brief this templated class will register an instance of a
/// GridCreator (class T) with the GridFactory.  It will ensure that
/// no GridCreator is registered twice and centralizes the registration logic
template< class T >
class GridRegistrator : public utility::factory::WidgetRegistrator<GridFactory,T>
{
public:
	typedef utility::factory::WidgetRegistrator<GridFactory,T> parent;
public:
	GridRegistrator() : parent() {}
};

class GridFactory
{
public:

	typedef std::map<std::string, GridCreatorOP > GridMap;
	typedef utility::tag::Tag Tag;
	typedef utility::tag::TagCOP TagCOP;

public:
	virtual ~GridFactory();

	static
	GridFactory * get_instance();

	void factory_register(GridCreatorOP creator);

	///@brief create Grid given grid tag
	GridBaseOP new_grid(utility::tag::TagCOP tag) const;

	///@brief create Grid given a serialized grid object
	GridBaseOP new_grid(utility::json_spirit::mObject data ) const;

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
	GridFactory();
	//unimplemented -- uncopyable
	GridFactory(GridFactory const & );
	GridFactory const & operator = (GridFactory const &);

	/// @brief private singleton creation function to be used with
	/// utility::thread::threadsafe_singleton
	static GridFactory * create_singleton_instance();

private:
#if defined MULTI_THREADED && defined CXX11
	static std::atomic< GridFactory * > instance_;
#else
	static GridFactory * instance_;
#endif

	GridMap grid_creator_map_;

};

}
}
}


#endif /* GRIDFACTORY_HH_ */
