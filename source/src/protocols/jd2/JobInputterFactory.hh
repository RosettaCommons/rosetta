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
#include <utility/SingletonBase.hh>
#include <utility/vector1.hh>
#include <utility/factory/WidgetRegistrator.hh>
#include <utility/pointer/ReferenceCount.hh>

// c++ headers
#include <map>

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

class JobInputterFactory : public utility::SingletonBase< JobInputterFactory >
{
public:
	typedef std::map< std::string, JobInputterCreatorOP > JobInputterMap;
	friend class utility::SingletonBase< JobInputterFactory >;

public:
	virtual ~JobInputterFactory();

	void factory_register( JobInputterCreatorOP creator );

	/// @brief return JobInputter defined by input parameters (contained in option system and #defines for MPI, etc)
	JobInputterOP get_new_JobInputter();

private:
	JobInputterOP get_JobInputter_from_string( std::string const & job_inputter_type );

	JobInputterFactory();

	// Unimplemented -- uncopyable
	JobInputterFactory( JobInputterFactory const & );
	JobInputterFactory const & operator = ( JobInputterFactory const & );

	/// @brief private singleton creation function to be used with
	/// utility::thread::threadsafe_singleton; used by the base class
	/// utility::SingletonBase
	static JobInputterFactory * create_singleton_instance();

private:

	JobInputterMap job_inputter_creator_map_;

};

} //namespace jd2
} //namespace protocols

#endif //INCLUDED_protocols_jd2_JobInputterFactory_hh
