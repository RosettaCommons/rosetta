// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/tensorflow_manager/RosettaTensorflowManager.cc
/// @brief A manager class for loading Tensorflow sessions and controlling mapping to the CPU or GPU.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#include <basic/tensorflow_manager/RosettaTensorflowManager.hh>


// Unit headers

// Project header
#include <platform/types.hh>

// Utility headers
#include <utility/thread/threadsafe_creation.hh>
#include <utility/pointer/memory.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/tensorflow_manager/RosettaTensorflowSessionContainer.hh>
#include <basic/citation_manager/CitationManager.hh>
#include <basic/citation_manager/UnpublishedModuleInfo.hh>
#include <basic/citation_manager/CitationCollection.fwd.hh>

// C++ headers
#include <iostream>
#include <string>
#include <tuple>
#include <functional>

// Construct tracer.
static basic::Tracer TR( "basic.tensorflow_manager.RosettaTensorflowManager" );

namespace basic {
namespace tensorflow_manager {

// Public methods /////////////////////////////////////////////////////////////
// Static constant data access

// Private methods ////////////////////////////////////////////////////////////
// Empty constructor
RosettaTensorflowManager::RosettaTensorflowManager()
{
#ifdef USE_TENSORFLOW
    // Make sure that the authors are cited!
    register_unpublished_author_info();
#endif
}

/// @brief Destructor.
/// @details Checks for sessions that are still in use and spits out warnings (but does not throw exceptions).
RosettaTensorflowManager::~RosettaTensorflowManager() {
#ifdef USE_TENSORFLOW
    {
        //Scope 1: check for sessions still in use, indexed by tuple.
#ifdef MULTI_THREADED
        utility::thread::ReadLockGuard readlock( sessions_mutex_ );
#endif //MULTI_THREADED
        for( std::map< RosettaTensorflowSessionIdentifier, RosettaTensorflowSessionContainerOP >::const_iterator it = sessions_.begin(); it!=sessions_.end(); ++it ) {
            platform::Size const usecount( it->second.use_count() );
            if( usecount > 1 ) {
                TR.Error << "Error in RosettaTensorflowManager destructor: The Tensorflow session indexed with <\"" << std::get<0>( it->first )  << "\", \"" << std::get<1>( it->first ) << "\", " << std::get<2>(it->first) << "> is in use by " <<  usecount - 1 << " other users!" << std::endl;
            }
        }
    }

    {
        //Scope 2: check for sessions still in use, indexed by string.
#ifdef MULTI_THREADED
        utility::thread::ReadLockGuard readlock( sessions_by_string_mutex_ );
#endif //MULTI_THREADED
        for( std::map< std::string, RosettaTensorflowSessionContainerOP >::const_iterator it = sessions_by_string_.begin(); it!=sessions_by_string_.end(); ++it ) {
            platform::Size const usecount( it->second.use_count() );
            if( usecount > 1 ) {
                TR.Error << "Error in RosettaTensorflowManager destructor: The Tensorflow session indexed with <\"" << it->first << "\"> is in use by " <<  usecount - 1 << " other users!" << std::endl;
            }
        }
    }
#endif //USE_TENSORFLOW
}

#ifdef USE_TENSORFLOW
/// @brief Create a Tensorflow session.
/// @details Triggers disk access!
RosettaTensorflowSessionContainerOP
RosettaTensorflowManager::create_tensorflow_session(
    std::string const & filename,
    std::string const & tags,
    TF_SessionOptions * sess_options
) {
    // Note that this unfortunately must be new instead of make_shared, since the RosettaTensorflowSessionContainer
    // constructor is private:
    std::cout << std::flush;
    return RosettaTensorflowSessionContainerOP( new RosettaTensorflowSessionContainer( filename, tags, sess_options ) );
}

/// @brief Get a previously-loaded Tensorflow session, or create it, cache it, and return it if it has not been previously loaded.
/// @details Threadsafe and lazy-loading.
/// @note If the session was created using custom session options, and you want the same session back, you must provide a pointer to
/// the session options!
RosettaTensorflowSessionContainerCOP
RosettaTensorflowManager::get_session(
    std::string const & filename,
    std::string const & tags,
    TF_SessionOptions* sess_options/*=nullptr*/
) {
    RosettaTensorflowSessionIdentifier const key( std::make_tuple( filename, tags, sess_options ) );
    std::function< RosettaTensorflowSessionContainerOP () > creator( std::bind( &RosettaTensorflowManager::create_tensorflow_session, std::cref( filename ), std::cref( tags ), sess_options ) );
	return utility::thread::safely_check_map_for_key_and_insert_if_absent( creator, SAFELY_PASS_MUTEX( sessions_mutex_ ), key, sessions_ );
}

/// @brief Get a previously-loaded Tensorflow session, or create it, cache it, and return it if it has not been previously loaded.
/// @details Threadsafe and lazy-loading.
/// @note This version allows a cached session to be retrieved by key string, without necessarily having the original options object.
RosettaTensorflowSessionContainerCOP
RosettaTensorflowManager::get_session(
    std::string const & filename,
    std::string const & tags,
    std::string const & key,
    TF_SessionOptions* sess_options/*=nullptr*/
) {
    std::function< RosettaTensorflowSessionContainerOP () > creator( std::bind( &RosettaTensorflowManager::create_tensorflow_session, std::cref( filename ), std::cref( tags ), sess_options ) );
    return utility::thread::safely_check_map_for_key_and_insert_if_absent( creator, SAFELY_PASS_MUTEX( sessions_by_string_mutex_ ), key, sessions_by_string_ );
}

/// @brief Get a previously-loaded Tensorflow session by its key string alone.
/// @details Threadsafe.
/// @note This assumes that the session was alredy created, and will throw an error if it was not!
RosettaTensorflowSessionContainerCOP
RosettaTensorflowManager::get_session(
    std::string const & key
) const {
#ifdef MULTI_THREADED
    utility::thread::ReadLockGuard lockguard( sessions_by_string_mutex_ );
#endif //MULTI_THREADED
    runtime_assert_string_msg( sessions_by_string_.count( key ) != 0, "Error in RosettaTensorflowManager::get_session(): No Tensorflow session corresponding to the key \"" + key + "\" has been loaded!" );
    return sessions_by_string_.at( key );
}

/// @brief Want to give the manager your own session that you loaded yourself? Use this.
/// @details Threadsafe.
/// @details The manager takes lifetime ownership of this session
/// @author Jack Maguire, jackmaguire1444@gmail.com
void
RosettaTensorflowManager::set_session(
    std::string const & key,
		RosettaTensorflowSessionContainerOP session
) {
	std::function< RosettaTensorflowSessionContainerOP () > creator = [=](){
		return session;
	};
	RosettaTensorflowSessionContainerCOP ptr = utility::thread::safely_check_map_for_key_and_insert_if_absent( creator, SAFELY_PASS_MUTEX( sessions_by_string_mutex_ ), key, sessions_by_string_ );
	runtime_assert( ptr != nullptr );
}

/// @brief Returns true if a session exists, false otherwise.
/// @details Note that, if you use this, it is possible for another thread to create a session with this
/// key right after this function resturns "false".  For this reason, all of the get_session functions,
/// above, are threadsafe, and will return the already-created session if another thread created a session.
/// @note This version checks the sessions indexed by filename, tag, and (optionally) options.
bool
RosettaTensorflowManager::session_exists(
    std::string const & filename,
    std::string const & tags,
    TF_SessionOptions* sess_options/*=nullptr*/
) const {
    RosettaTensorflowSessionIdentifier const key( std::make_tuple( filename, tags, sess_options ) );
#ifdef MULTI_THREADED
    utility::thread::ReadLockGuard lockguard( sessions_mutex_ );
#endif //MULTI_THREADED
    return ( sessions_.count(key) != 0 );
}

/// @brief Returns true if a session exists, false otherwise.
/// @details Note that, if you use this, it is possible for another thread to create a session with this
/// key right after this function resturns "false".  For this reason, all of the get_session functions,
/// above, are threadsafe, and will return the already-created session if another thread created a session.
/// @note This version checks the sessions indexed by a single key string.
bool
RosettaTensorflowManager::session_exists(
    std::string const & key
) const {
#ifdef MULTI_THREADED
    utility::thread::ReadLockGuard lockguard( sessions_by_string_mutex_ );
#endif //MULTI_THREADED
    return ( sessions_by_string_.count(key) != 0 );
}

/// @brief Get the current Tensorflow version against which Rosetta is linked.
std::string
RosettaTensorflowManager::get_tensorflow_version() const {
    return std::string( TF_Version() );
}

#endif //USE_TENSORFLOW

////////////////// PRIVATE METHODS //////////////////

#ifdef USE_TENSORFLOW

/// @brief Register the author information with the Rosetta CitationManager.
/// @details This lists Vikram K. Mulligan, Jack Maguire, and Sergey Lyskov as
/// contributors to the RosettaTensorflowManager.
void
RosettaTensorflowManager::register_unpublished_author_info() const {
    using namespace basic::citation_manager;

    UnpublishedModuleInfoOP info( utility::pointer::make_shared< UnpublishedModuleInfo >(
        "RosettaTensorflowManager",
        CitedModuleType::Singleton,
        "Vikram K. Mulligan",
        "Systems Biology, Center for Computational Biology, Flatiron Institute",
        "vmulligan@flatironinstitute.org",
        "Created the RosettaTensorflowManager."
    ) );

    info->add_author(
        "Jack Magure",
        "Menten AI",
        "jack.maguire@menten.ai",
        "Expanded RosettaTensorflowManager capabilities for multi-head jobs and wrote tests."
    );

    info->add_author(
        "Sergey Lyskov",
        "Gray Lab, Department of Chemical & Biomolecular Engineering, Johns Hopkins University",
        "Sergey.Lyskov@jhu.edu",
        "Added testing infrastructure and helped to create the Rosetta-Tensorflow linked build."
    );

    utility::vector1< UnpublishedModuleInfoCOP > infovect{ info };
    CitationManager::get_instance()->add_unpublished_modules( infovect );
}

#endif //USE_TENSORFLOW

} //tensorflow_manager
} //basic
