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

class ResourceManager {
public:
	static ResourceManager * get_instance();

protected:
	ResourceManager(); // singleton, private constructor

public: // ResourcePlanner interface

	virtual
	void
	add_resource(
		ResourceTag const & resource_tag,
		ResourceOP resource);

	virtual
	bool
	has_resource(
		ResourceTag const & resource_tag) const;

	virtual
	ResourceOP
	find_resource(
		ResourceTag const & resource_tag) const;

	virtual
	void
	free_resource(
		ResourceTag const & resource_tag);

	/// @brief wipe out all data; useful for unit testing, but hard to fathom how it would be useful
	/// otherwise.
	virtual
	void
	clear() = 0;

public: // Protocol interface

	virtual
	bool
	has_resource_with_description(
		ResourceDescription const & resource_description) = 0;


	///@brief called by the protocol, return initialized resource
	virtual
	ResourceOP
	get_resource(
		ResourceDescription const & resource_description) = 0;

public: // Options interface

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
