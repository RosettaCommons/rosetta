// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/JobOutputterFactory.hh
/// @brief  JobOutputterFactory
/// @author Steven Lewis smlewi@gmail.com, Brian Weitzner brian.weitzner@gmail.com


#ifndef INCLUDED_protocols_jd2_JobOutputterFactory_hh
#define INCLUDED_protocols_jd2_JobOutputterFactory_hh

// Unit Headers
#include <protocols/jd2/JobOutputterFactory.fwd.hh>
#include <protocols/jd2/JobOutputterCreator.hh>
#include <protocols/jd2/JobOutputter.fwd.hh>

// Utility Headers
#include <utility/SingletonBase.hh>
#include <utility/factory/WidgetRegistrator.hh>
#include <utility/pointer/ReferenceCount.hh>

// c++ headers
#include <map>

#include <utility/vector1.hh>

namespace protocols {
namespace jd2 {

/// @brief This templated class will register an instance of an
/// JobOutputterCreator (class T) with the JobOutputterFactory.  It will ensure
/// that no JobOutputterCreator is registered twice, and, centralizes
/// this registration logic so that thread safety issues can be handled in
/// one place
template < class T >
class JobOutputterRegistrator : public utility::factory::WidgetRegistrator< JobOutputterFactory, T >
{
public:
	typedef utility::factory::WidgetRegistrator< JobOutputterFactory, T > parent;
public:
	JobOutputterRegistrator() : parent() {}
};


class JobOutputterFactory : public utility::SingletonBase< JobOutputterFactory >
{
public:
	friend class utility::SingletonBase< JobOutputterFactory >;
	typedef std::map< std::string, JobOutputterCreatorOP > JobOutputterMap;

public:
	virtual ~JobOutputterFactory();


	void factory_register( JobOutputterCreatorOP creator );

	/// @brief return JobOutputter defined by output parameters (contained in option system and #defines for MPI, etc)
	JobOutputterOP get_new_JobOutputter();

	/// @brief return JobOutputter defined by output parameters (contained in option system and #defines for MPI, etc).  The difference is that if the option system, etc, says nothing about output (which as of this writing defaults to PDBJobOutputter), this function leaves the input Outputter unchanged.  This allows overriding the default outputter choice in your executable (without abusing the mutability of the options system)
	JobOutputterOP get_new_JobOutputter( JobOutputterOP default_jobout );

private:
	JobOutputterOP get_JobOutputter_from_string( std::string const & job_outputter_type );

	JobOutputterFactory();

	// Unimplemented -- uncopyable
	JobOutputterFactory( JobOutputterFactory const & );
	JobOutputterFactory const & operator = ( JobOutputterFactory const & );

	/// @brief private singleton creation function to be used with
	/// utility::thread::threadsafe_singleton
	static JobOutputterFactory * create_singleton_instance();

private:

	JobOutputterMap job_outputter_creator_map_;

};

} //namespace jd2
} //namespace protocols

#endif //INCLUDED_protocols_jd2_JobOutputterFactory_hh
