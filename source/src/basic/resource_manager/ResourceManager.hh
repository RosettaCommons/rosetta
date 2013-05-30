// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/resource_manager/ResourceManager.hh
/// @brief
/// @author

#ifndef INCLUDED_basic_resource_manager_ResourceManager_hh
#define INCLUDED_basic_resource_manager_ResourceManager_hh

//unit headers
#include <basic/resource_manager/ResourceManager.fwd.hh>
#include <basic/resource_manager/types.hh>

//project headers
#include <basic/resource_manager/ResourceOptions.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/options/keys/BooleanOptionKey.fwd.hh>
#include <utility/options/keys/BooleanVectorOptionKey.fwd.hh>
#include <utility/options/keys/FileOptionKey.fwd.hh>
#include <utility/options/keys/FileVectorOptionKey.fwd.hh>
#include <utility/options/keys/IntegerOptionKey.fwd.hh>
#include <utility/options/keys/IntegerVectorOptionKey.fwd.hh>
#include <utility/options/keys/OptionKey.fwd.hh>
#include <utility/options/keys/PathOptionKey.fwd.hh>
#include <utility/options/keys/PathVectorOptionKey.fwd.hh>
#include <utility/options/keys/RealOptionKey.fwd.hh>
#include <utility/options/keys/RealVectorOptionKey.fwd.hh>
#include <utility/options/keys/StringOptionKey.fwd.hh>
#include <utility/options/keys/StringVectorOptionKey.fwd.hh>
#include <utility/vector1.fwd.hh>
#include <utility/file/FileName.fwd.hh>
#include <utility/file/PathName.fwd.hh>

//C++ headers
#include <map>
#include <string>

namespace basic {
namespace resource_manager {

/// @brief The ResourceManager is a singleton class responsible for holding and eventually
/// deallocating resources which may be shared between multiple jobs.  A protocol
/// may communicate directly with the ResourceManager, requesting resources, but remaining
/// unaware of where those resources came from, or whether the same resource is
/// being used in multiple contexts.
///
/// A protocol should request a resource using a "resource description" (a string), which
/// is a generic way of referring to a piece of data.  For example, a protocol might request
/// the native Pose, by asking for a resource with the description "native".  The
/// ResourceManager's job is to return the resource matching that description; when
/// the "native" is requested, it is requested in some context (e.g. in a protocol running
/// under JD2) and the ResourceManager serves as a backbone in which to deliever
/// context-specific data to the protocol, while keeping the protocol ignorant and
/// independent of the surrounding context.
///
/// NOTE: Because the logic for deciding which of the derived ResourceManager classes
/// should be instantiated depends on the options system, the ResourceManager should
/// not be requested until after core::init::init() has been called (i.e. it should not
/// be requested at load time.)
class ResourceManager {
public:
	/// @brief Singleton accessor: request access to the resource manager.
	static ResourceManager * get_instance();

protected:
	/// @brief singleton, private constructor
	ResourceManager();

protected: // Derived class interface

	/// @brief Derived classes, when they instantiate a resource, should name it
	/// and give it to the base class for it to hold.  A resource may be later
	/// deallocated using the "free_resource" function.
	virtual
	void
	add_resource(
		ResourceTag const & resource_tag,
		ResourceOP resource);

	/// @brief Derived classes, when they are ready to deallocate a resource,
	/// may do so in this function.
	virtual
	void
	free_resource(
		ResourceTag const & resource_tag);

public:
	/// @brief Does a resource with a given name exist?
	virtual
	bool
	has_resource(
		ResourceTag const & resource_tag) const;

	/// @brief Get a resource with a given name.
	virtual
	ResourceOP
	find_resource(
		ResourceTag const & resource_tag);

	/// @brief wipe out all data; useful for unit testing, but hard to
	/// fathom how it would be useful otherwise.
	virtual
	void
	clear();

	virtual
	void
	show( std::ostream & out ) const;

	friend
	std::ostream &
	operator<< ( std::ostream & out, ResourceManager const & resource_manager );


public: // Protocol interface

	/// @brief A protocol may ask whether a resource with a given resource description
	/// has been provided.  It is possible that some resources are available sometimes
	/// when running a protocol, but not always (e.g. when running abinitio, sometimes
	/// you don't know what the native structure is).
	virtual
	bool
	has_resource_with_description(
		ResourceDescription const & resource_description) = 0;


	/// @brief Called by the protocol, returns the resource appropriate for the
	/// context in which it is requested (which the protocol should remain
	/// ignorant of).
	virtual
	ResourceOP
	get_resource(
		ResourceDescription const & resource_description) = 0;

public: // Options interface

	/// The following 12 functions allow protocols to request options that
	/// may have been set for the context in which the protocol is being run.

	virtual
	bool
	get_option(
		utility::options::BooleanOptionKey key ) const = 0;

	virtual
	utility::vector1< bool >  const &
	get_option(
		utility::options::BooleanVectorOptionKey key ) const = 0;

	virtual
	utility::file::FileName  const &
	get_option(
		utility::options::FileOptionKey key ) const = 0;

	virtual
	utility::vector1< utility::file::FileName >  const &
	get_option(
		utility::options::FileVectorOptionKey key ) const = 0;

	virtual
	int
	get_option(
		utility::options::IntegerOptionKey key ) const = 0;

	virtual
	utility::vector1< int >  const &
	get_option(
		utility::options::IntegerVectorOptionKey key ) const = 0;

	virtual
	utility::file::PathName  const &
	get_option(
		utility::options::PathOptionKey key ) const = 0;

	virtual
	utility::vector1< utility::file::PathName >  const &
	get_option(
		utility::options::PathVectorOptionKey key ) const = 0;

	virtual
	platform::Real
	get_option(
		utility::options::RealOptionKey key ) const = 0;

	virtual
	utility::vector1< platform::Real >  const &
	get_option(
		utility::options::RealVectorOptionKey key ) const = 0;

	virtual
	std::string  const &
	get_option(
		utility::options::StringOptionKey key ) const = 0;

	virtual
	utility::vector1< std::string > const &
	get_option(
		utility::options::StringVectorOptionKey key ) const = 0;

	/// The following 12 functions allow protocols to request if a particular
	/// option has been set (e.g. on the command line or for the current job).

	virtual
	bool
	has_option(
		utility::options::BooleanOptionKey key ) const = 0;

	virtual
	bool
	has_option(
		utility::options::BooleanVectorOptionKey key ) const = 0;

	virtual
	bool
	has_option(
		utility::options::FileOptionKey key ) const = 0;

	virtual
	bool
	has_option(
		utility::options::FileVectorOptionKey key ) const = 0;

	virtual
	bool
	has_option(
		utility::options::IntegerOptionKey key ) const = 0;

	virtual
	bool
	has_option(
		utility::options::IntegerVectorOptionKey key ) const = 0;

	virtual
	bool
	has_option(
		utility::options::PathOptionKey key ) const = 0;

	virtual
	bool
	has_option(
		utility::options::PathVectorOptionKey key ) const = 0;

	virtual
	bool
	has_option(
		utility::options::RealOptionKey key ) const = 0;

	virtual
	bool
	has_option(
		utility::options::RealVectorOptionKey key ) const = 0;

	virtual
	bool
	has_option(
		utility::options::StringOptionKey key ) const = 0;

	virtual
	bool
	has_option(
		utility::options::StringVectorOptionKey key ) const = 0;


private:
	static ResourceManager * instance_;

	typedef std::map< ResourceTag, ResourceOP > ResourcesMap;
	ResourcesMap resources_;

};

} // namespace resource_manager
} // namespace basic

#endif //INCLUDED_basic_resource_manager_ResourceManager_hh
