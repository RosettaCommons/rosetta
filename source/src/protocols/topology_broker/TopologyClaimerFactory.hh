// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/topology_broker/TopologyClaimerFactory.hh
/// @author Oliver Lange
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef INCLUDED_protocols_topology_broker_TopologyClaimerFactory_hh
#define INCLUDED_protocols_topology_broker_TopologyClaimerFactory_hh

// Project headers
#include <protocols/topology_broker/TopologyClaimer.fwd.hh>

// Utility headers
#include <utility/SingletonBase.hh>

// External headers
#include <boost/utility.hpp>

// C/C++ headers
#include <map>

#ifdef MULTI_THREADED
#ifdef CXX11
// C++11 Headers
#include <atomic>
#include <mutex>
#endif
#endif

namespace protocols {
namespace topology_broker {

/// @brief A non-copyable factory for instantiating TopologyClaimers by name.
/// Commonly used TopologyClaimers are registered in the constructor. Additional
/// claimers can be registered after the fact using the add_type() method.
class TopologyClaimerFactory : boost::noncopyable {
public:
	/// @brief Returns an instance to the singleton
	static TopologyClaimerFactory const& get_instance();

	/// @brief Returns a new instance of the TopologyClaimer identified by <name>
	TopologyClaimerOP newTopologyClaimer(const std::string& name) const;

	/// @brief Registers the TopologyClaimer with name <name>. Instances of this
	/// class can be retrieved by calling newTopologyClaimer(<name>).
	void add_type(const std::string& name, TopologyClaimerOP claimer);

	/// @brief Registers the TopologyClaimer using the name returned by
	/// claimer->type(). Instances of this class can be retrieved by calling
	/// newTopologyClaimer(<name>).
	void add_type(TopologyClaimerOP claimer);

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
	/// @brief Constructs a new instance and initializes the lookup table
	/// <claimers_> with commonly used types
	TopologyClaimerFactory();

	/// @brief Frees resources associated with this object
	~TopologyClaimerFactory();

	/// @brief private singleton creation function to be used with
	/// utility::thread::threadsafe_singleton
	static TopologyClaimerFactory * create_singleton_instance();

private:

	/// @brief A map that associates claimer names with claimer types. Used by the
	/// newTopologyClaimer() method to instantiate new claimers by name.
	mutable std::map<std::string, TopologyClaimerOP> claimers_;

	/// @brief A pointer to the singleton instance of the factory object.
	/// Resources associated with the object are released on destruction.
	/// APL Question: Should this be one-instance-per-program (singleton) or one-instance-per-job?
#if defined MULTI_THREADED && defined CXX11
	static std::atomic< TopologyClaimerFactory * > instance_;
#else
	static TopologyClaimerFactory * instance_;
#endif

};

}
}

#endif  // INCLUDED_protocols_topology_broker_TopologyClaimerFactory_hh
