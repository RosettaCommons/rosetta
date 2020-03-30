// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/tensorflow_manager/RosettaTensorflowManager
/// @brief A manager class for loading Tensorflow sessions and controlling mapping to the CPU or GPU.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)


#ifndef INCLUDED_basic_tensorflow_manager_RosettaTensorflowManager_hh
#define INCLUDED_basic_tensorflow_manager_RosettaTensorflowManager_hh

// Unit headers
#include <basic/tensorflow_manager/RosettaTensorflowManager.fwd.hh>

// Basic headers
#include <basic/tensorflow_manager/RosettaTensorflowSessionContainer.fwd.hh>

// Utility header
#include <utility/SingletonBase.hh>
#include <utility/vector1.hh>

// C++ header
#include <tuple>
#include <map>

#ifdef MULTI_THREADED
#include <utility/thread/ReadWriteMutex.hh>
#endif

#ifdef USE_TENSORFLOW
#include <tensorflow/c/c_api.h>
#endif

namespace basic {
namespace tensorflow_manager {

#ifdef USE_TENSORFLOW
typedef std::tuple< std::string, std::string, TF_SessionOptions* > RosettaTensorflowSessionIdentifier;
#endif

/// @brief A manager class for loading Tensorflow sessions and controlling mapping to the CPU or GPU.
class RosettaTensorflowManager : public utility::SingletonBase< RosettaTensorflowManager > {
	friend class utility::SingletonBase< RosettaTensorflowManager >;

private:  // Private methods ///////////////////////////////////////////////////
	/// @brief  Empty constructor
	RosettaTensorflowManager();

	/// @brief Copy constructor is explicitly deleted.
	RosettaTensorflowManager(RosettaTensorflowManager const & ) = delete;

	/// @brief Assignment operator is explicitly deleted.
	RosettaTensorflowManager operator=(RosettaTensorflowManager const & ) = delete;

	/// @brief Destructor.
	/// @details Checks for sessions that are still in use and spits out warnings (but does not throw exceptions).
	~RosettaTensorflowManager();

#ifdef USE_TENSORFLOW
	/// @brief Create a Tensorflow session.
	/// @details Triggers disk access!
	static RosettaTensorflowSessionContainerOP
	create_tensorflow_session(
		std::string const & filename,
		std::string const & tags,
		TF_SessionOptions * sess_options
	);
#endif //USE_TENSORFLOW

public:   // Public methods ///////////////////////////////////////////////////

#ifdef USE_TENSORFLOW
	/// @brief Get a previously-loaded Tensorflow session, or create it, cache it, and return it if it has not been previously loaded.
	/// @details Threadsafe and lazy-loading.
	/// @note If the session was created using custom session options, and you want the same session back, you must provide a pointer to
	/// the session options!
	RosettaTensorflowSessionContainerCOP get_session( std::string const & filename, std::string const & tags, TF_SessionOptions* sess_options=nullptr );

	/// @brief Get a previously-loaded Tensorflow session, or create it, cache it, and return it if it has not been previously loaded.
	/// @details Threadsafe and lazy-loading.
	/// @note This version allows a cached session to be retrieved by key string, without necessarily having the original options object.
	RosettaTensorflowSessionContainerCOP get_session( std::string const & filename, std::string const & tags, std::string const & key, TF_SessionOptions* sess_options=nullptr );

	/// @brief Get a previously-loaded Tensorflow session by its key string alone.
	/// @details Threadsafe.
	/// @note This assumes that the session was alredy created, and will throw an error if it was not!
	RosettaTensorflowSessionContainerCOP get_session( std::string const & key ) const;

	/// @brief Want to give the manager your own session that you loaded yourself? Use this.
	/// @details Threadsafe.
	/// @details The manager takes lifetime ownership of this session
	/// @author Jack Maguire, jackmaguire1444@gmail.com
	void set_session( std::string const & key, RosettaTensorflowSessionContainerOP session );

	/// @brief Returns true if a session exists, false otherwise.
	/// @details Note that, if you use this, it is possible for another thread to create a session with this
	/// key right after this function resturns "false".  For this reason, all of the get_session functions,
	/// above, are threadsafe, and will return the already-created session if another thread created a session.
	/// @note This version checks the sessions indexed by filename, tag, and (optionally) options.
	bool session_exists( std::string const & filename, std::string const & tags, TF_SessionOptions* sess_options=nullptr ) const;

	/// @brief Returns true if a session exists, false otherwise.
	/// @details Note that, if you use this, it is possible for another thread to create a session with this
	/// key right after this function resturns "false".  For this reason, all of the get_session functions,
	/// above, are threadsafe, and will return the already-created session if another thread created a session.
	/// @note This version checks the sessions indexed by a single key string.
	bool session_exists( std::string const & key ) const;

	/// @brief Get the current Tensorflow version against which Rosetta is linked.
	std::string get_tensorflow_version() const;

#endif //USE_TENSORFLOW

private:  // Private methods //////////////////////////////////////////////////

#ifdef USE_TENSORFLOW

	/// @brief Register the author information with the Rosetta CitationManager.
	/// @details This lists Vikram K. Mulligan, Jack Maguire, and Sergey Lyskov as
	/// contributors to the RosettaTensorflowManager.
	void register_unpublished_author_info() const;

#endif //USE_TENSORFLOW

private:  // Private data /////////////////////////////////////////////////////

#ifdef USE_TENSORFLOW
	/// @brief The Tensorflow sessions that we have loaded, stored in containers for safety.
	std::map< RosettaTensorflowSessionIdentifier, RosettaTensorflowSessionContainerOP > sessions_;

	/// @brief The Tensorflow sessions that we have loaded that were indexed by a string.
	std::map< std::string, RosettaTensorflowSessionContainerOP > sessions_by_string_;

#ifdef MULTI_THREADED
	mutable utility::thread::ReadWriteMutex sessions_mutex_;
	mutable utility::thread::ReadWriteMutex sessions_by_string_mutex_;
#endif //MULTI_THREADED

#endif //USE_TENSORFLOW

};

} //tensorflow_manager
} //basic

#endif //INCLUDED_basic/tensorflow_manager_RosettaTensorflowManager_fwd_hh



