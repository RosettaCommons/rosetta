// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/rosetta_scripts/PosePropertyReporterFactory.hh
/// @brief	Factory for PosePropertyReporters
/// @author Luki Goldschmidt <lugo@uw.edu>

#ifndef INCLUDED_protocols_moves_PosePropertyReporterFactory_hh
#define INCLUDED_protocols_moves_PosePropertyReporterFactory_hh

// Unit Headers
#include <protocols/rosetta_scripts/PosePropertyReporterFactory.fwd.hh>
#include <protocols/rosetta_scripts/PosePropertyReporter.fwd.hh>
#include <protocols/rosetta_scripts/PosePropertyReporter.hh>
#include <protocols/rosetta_scripts/PosePropertyReporterCreator.hh>

// Project Headers
#include <basic/datacache/DataMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// Utility Headers
#include <utility/factory/WidgetRegistrator.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>

// c++ headers
#include <map>
#include <set>


#ifdef MULTI_THREADED
#ifdef CXX11
// C++11 Headers
#include <atomic>
#include <mutex>
#endif
#endif

namespace protocols {
namespace rosetta_scripts {

/// @brief This templated class will register an instance of an
/// PosePropertyReporterCreator (class T) with the PosePropertyReporterFactory.  It will ensure
/// that no PosePropertyReporterCreator is registered twice, and, centralizes
/// this registration logic so that thread safety issues can be handled in
/// one place
template < class T >
class PosePropertyReporterRegistrator : public utility::factory::WidgetRegistrator< PosePropertyReporterFactory, T >
{
public:
	typedef utility::factory::WidgetRegistrator< PosePropertyReporterFactory, T > parent;
public:
	PosePropertyReporterRegistrator() : parent() {}
};


class PosePropertyReporterFactory // Singleton: should not derive from ReferenceCount
{
public:
	typedef std::map< std::string, PosePropertyReporterCreatorOP > PosePropertyReporterMap;
	typedef utility::tag::Tag Tag;
	typedef utility::tag::TagCOP TagCOP;
	typedef core::pose::Pose Pose;

public:
	virtual ~PosePropertyReporterFactory();

	static
	PosePropertyReporterFactory * get_instance();

	void factory_register( PosePropertyReporterCreatorOP creator );

	/// @brief Create a PosePropertyReporter given its identifying string
	PosePropertyReporterOP newPosePropertyReporter( std::string const & );

	///@brief return new PosePropertyReporter by Tag parsing; the identifying string for the PosePropertyReporter is in the Tag
	PosePropertyReporterOP
	newPosePropertyReporter(
		TagCOP,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	);

private:
	PosePropertyReporterFactory();

	// Unimplemented -- uncopyable
	PosePropertyReporterFactory( PosePropertyReporterFactory const & );
	PosePropertyReporterFactory const & operator = ( PosePropertyReporterFactory const & );

	/// @brief private singleton creation function to be used with
	/// utility::thread::threadsafe_singleton
	static PosePropertyReporterFactory * create_singleton_instance();

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
#if defined MULTI_THREADED && defined CXX11
	static std::atomic< PosePropertyReporterFactory * > instance_;
#else
	static PosePropertyReporterFactory * instance_;
#endif

	PosePropertyReporterMap reporter_creator_map_;
};

} //namespace rosetta_scripts
} //namespace protocols

#endif
