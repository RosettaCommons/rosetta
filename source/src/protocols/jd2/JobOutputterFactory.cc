// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/JobOutputterFactory.cc
/// @brief  JobOutputterFactory
/// @author Steven Lewis smlewi@gmail.com, Brian Weitzner brian.weitzner@gmail.com

#include <protocols/jd2/JobOutputterFactory.hh>
#include <protocols/jd2/JobOutputterCreator.hh>
#include <protocols/jd2/JobOutputter.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/jd2.OptionKeys.gen.hh>
//#include <basic/options/keys/run.OptionKeys.gen.hh>
//#include <basic/options/keys/parser.OptionKeys.gen.hh>
//#include <basic/options/keys/enzdes.OptionKeys.gen.hh>

// Utility headers
#include <utility/exit.hh> // runtime_assert, utility_exit_with_message
#include <basic/Tracer.hh>

#include <utility/thread/threadsafe_creation.hh>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>

// Singleton instance and mutex static data members
namespace utility {

using protocols::jd2::JobOutputterFactory;

#if defined MULTI_THREADED
template <> std::mutex utility::SingletonBase< JobOutputterFactory >::singleton_mutex_{};
template <> std::atomic< JobOutputterFactory * > utility::SingletonBase< JobOutputterFactory >::instance_( 0 );
#else
template <> JobOutputterFactory * utility::SingletonBase< JobOutputterFactory >::instance_( nullptr );
#endif

}

namespace protocols {
namespace jd2 {

static THREAD_LOCAL basic::Tracer TR( "protocols.jd2.JobOutputterFactory" );

JobOutputterFactory *
JobOutputterFactory::create_singleton_instance()
{
	return new JobOutputterFactory;
}

JobOutputterFactory::JobOutputterFactory()
{}

JobOutputterFactory::~JobOutputterFactory()= default;

/// @brief add a JobOutputter prototype, using its default type name as the map key
void
JobOutputterFactory::factory_register( JobOutputterCreatorOP creator )
{
	runtime_assert( creator != nullptr );
	std::string const job_outputter_type( creator->keyname() );
	if ( job_outputter_creator_map_.find( job_outputter_type ) != job_outputter_creator_map_.end() ) {
		utility_exit_with_message("JobOutputterFactory::factory_register already has a mover creator with name \"" + job_outputter_type + "\".  Conflicting JobOutputter names" );
	}
	job_outputter_creator_map_[ job_outputter_type ] = creator;
}

/// @details return get_JobOutputter_from_string( "JobOutputter by key lookup in map
JobOutputterOP
JobOutputterFactory::get_JobOutputter_from_string( std::string const & job_outputter_type )
{
	//get pointer to Creator
	JobOutputterMap::const_iterator iter( job_outputter_creator_map_.find( job_outputter_type ) );
	if ( iter != job_outputter_creator_map_.end() ) { //if Creator has an entry
		if ( ! iter->second ) { //if Creator inexplicably fails to exist, crash
			utility_exit_with_message( "Error: JobOutputterCreatorOP for " + job_outputter_type + " is NULL, you should never have been able to get here!" );
		}

		//if creator exists, return a JobOutputter from it (this is good)
		return iter->second->create_JobOutputter();
	} else { //else, a non-existent JobOutputter has been requested.  Print existing ones and exit.
		TR << "Available : ";
		for ( JobOutputterMap::const_iterator mover_it = job_outputter_creator_map_.begin(); mover_it != job_outputter_creator_map_.end(); ++mover_it ) {
			TR << mover_it->first<<", ";
		}
		TR << std::endl;
		utility_exit_with_message( job_outputter_type + " is not known to the JobOutputterFactory. Was it registered via a JobOutputterRegistrator in one of the init.cc files (devel/init.cc or protocols/init.cc)?" );
		return nullptr;
	}
}

/// @brief return new JobOutputter from logic of option system plus compilation options.  All the logic for determining job output type lives here.
JobOutputterOP
JobOutputterFactory::get_new_JobOutputter()
{
	//initial copy of this code copied at XRW2 by SML+BDW from about SVN:46190 from JobDistributorFactory.cc
	if ( basic::options::option[ basic::options::OptionKeys::out::file::silent ].user() ) {
		return get_JobOutputter_from_string( "SilentFileJobOutputter" );
	} else if ( basic::options::option[basic::options::OptionKeys::out::file::atom_tree_diff].user() ) {
		return get_JobOutputter_from_string( "AtomTreeDiffJobOutputter" );
	} else if ( basic::options::option[basic::options::OptionKeys::out::file::score_only].user() ) {
		return get_JobOutputter_from_string( "ScoreOnlyJobOutputter" );
	} else if ( basic::options::option[ basic::options::OptionKeys::jd2::no_output ].value() || basic::options::option[ basic::options::OptionKeys::out::nooutput ] ) {
		return get_JobOutputter_from_string( "NoOutputJobOutputter" );
	} else if ( basic::options::option[ basic::options::OptionKeys::jd2::enzdes_out].user() ) {
		return get_JobOutputter_from_string( "EnzdesJobOutputter" );
	} else if ( basic::options::option[ basic::options::OptionKeys::out::use_database].user() ) {
		return get_JobOutputter_from_string( "DatabaseJobOutputter" );
	} else if ( basic::options::option[ basic::options::OptionKeys::out::mmCIF].value() ) {
		return get_JobOutputter_from_string( "mmCIFJobOutputter" );
	} else { //currently default; may need an if in the future
		return get_JobOutputter_from_string( "PDBJobOutputter" );
	}
	return get_JobOutputter_from_string( "PDBJobOutputter" ); //default case may change in the future

}

/// @brief return JobOutputter defined by output parameters (contained in option system and #defines for MPI, etc).  The difference is that if the option system, etc, says nothing about output (which as of this writing defaults to PDBJobOutputter), this function leaves the input Outputter unchanged.  This allows overriding the default outputter choice in your executable (without abusing the mutability of the options system)
JobOutputterOP
JobOutputterFactory::get_new_JobOutputter( JobOutputterOP default_jobout ) {

	//it would be really nice to figure out how to combine the logic in these two versions of the function in to one - perhaps with default NULL pointer?

	if ( basic::options::option[ basic::options::OptionKeys::out::file::silent ].user() ) {
		return get_JobOutputter_from_string( "SilentFileJobOutputter" );
	} else if ( basic::options::option[basic::options::OptionKeys::out::pdb].user() ) {
		return get_JobOutputter_from_string( "PDBJobOutputter" );
	} else if ( basic::options::option[basic::options::OptionKeys::out::file::atom_tree_diff].user() ) {
		return get_JobOutputter_from_string( "AtomTreeDiffJobOutputter" );
	} else if ( basic::options::option[basic::options::OptionKeys::out::file::score_only].user() ) {
		return get_JobOutputter_from_string( "ScoreOnlyJobOutputter" );
	} else if ( basic::options::option[ basic::options::OptionKeys::jd2::no_output ].value() || basic::options::option[ basic::options::OptionKeys::out::nooutput ] ) {
		return get_JobOutputter_from_string( "NoOutputJobOutputter" );
	} else if ( basic::options::option[ basic::options::OptionKeys::jd2::enzdes_out].user() ) {
		return get_JobOutputter_from_string( "EnzdesJobOutputter" );
	} else if ( basic::options::option[ basic::options::OptionKeys::out::use_database].user() ) {
		return get_JobOutputter_from_string( "DatabaseJobOutputter" );
	} else {
		return default_jobout;
	}
	return default_jobout;

}

} //namespace jd2
} //namespace protocols
