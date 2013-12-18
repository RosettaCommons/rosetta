// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/JobInputterFactory.hh
/// @brief  JobInputterFactory
/// @author Steven Lewis smlewi@gmail.com, Brian Weitzner brian.weitzner@gmail.com


#ifndef INCLUDED_protocols_jd2_JobInputterFactory_hh
#define INCLUDED_protocols_jd2_JobInputterFactory_hh

// Unit Headers
#include <protocols/jd2/JobInputterFactory.fwd.hh>
#include <protocols/jd2/JobInputterCreator.hh>
#include <protocols/jd2/JobInputter.fwd.hh>

// Utility Headers
#include <utility/factory/WidgetRegistrator.hh>
#include <utility/pointer/ReferenceCount.hh>

// c++ headers
#include <map>

#include <utility/vector1.hh>

#ifdef MULTI_THREADED
#ifdef CXX11
// C++11 Headers
#include <thread>
#endif
#endif

namespace protocols {
namespace jd2 {

/// @brief This templated class will register an instance of an
/// JobInputterCreator (class T) with the JobInputterFactory.  It will ensure
/// that no JobInputterCreator is registered twice, and, centralizes
/// this registration logic so that thread safety issues can be handled in
/// one place
template < class T >
class JobInputterRegistrator : public utility::factory::WidgetRegistrator< JobInputterFactory, T >
{
public:
	typedef utility::factory::WidgetRegistrator< JobInputterFactory, T > parent;
public:
	JobInputterRegistrator() : parent() {}
};


class JobInputterFactory
{
public:
	typedef std::map< std::string, JobInputterCreatorOP > JobInputterMap;

public:
	virtual ~JobInputterFactory();

	static
	JobInputterFactory * get_instance();

	void factory_register( JobInputterCreatorOP creator );

	///@brief return JobInputter defined by input parameters (contained in option system and #defines for MPI, etc)
	JobInputterOP get_new_JobInputter();

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
	JobInputterOP get_JobInputter_from_string( std::string const & job_inputter_type );

  JobInputterFactory();

	// Unimplemented -- uncopyable
	JobInputterFactory( JobInputterFactory const & );
	JobInputterFactory const & operator = ( JobInputterFactory const & );

	/// @brief private singleton creation function to be used with
	/// utility::thread::threadsafe_singleton
	static JobInputterFactory * create_singleton_instance();

private:
	static JobInputterFactory * instance_;

	JobInputterMap job_inputter_creator_map_;

};

} //namespace jd2
} //namespace protocols

#endif //INCLUDED_protocols_jd2_JobInputterFactory_hh
