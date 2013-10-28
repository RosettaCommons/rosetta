// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/resource_manager/planner/JD2ResourceManager.hh
/// @brief  The ResourceManager that is compatible with the JD2 JobDistributor
/// @author Andrew Leaver-Fay
/// @author Brian Weitzner
/// @author Matthew O'Meara

#ifndef INCLUDED_protocols_jd2_JD2ResourceManager_hh
#define INCLUDED_protocols_jd2_JD2ResourceManager_hh

// Unit Headers

// Package headers
#include <protocols/jd2/JD2ResourceManager.fwd.hh>
#include <basic/resource_manager/LazyResourceManager.hh>
#include <basic/resource_manager/FallbackConfiguration.fwd.hh>
#include <basic/resource_manager/types.hh>

//C++ headers
#include <istream>
#include <string>

namespace protocols {
namespace jd2 {


/// @brief The %JD2ResourceManager is the ResourceManager that should be used when
/// running protocols under the jd2 JobDistributor.
///
/// @details The purpose of the ResourceManager is to disentangle the process of
/// feeding protocols with Resources they need from the job distribution system in
/// which they run.  When protocols request a resource under the ResourceManager,
/// they do so without knowing what kind of ResourceManager they are communicating
/// with.  Protocols that are designed to run under JD2 could conceivably be run
/// under any other job management scheme (e.g. a complicated MPI protocol).  It's
/// the job of the %JD2ResourceManager to determine the right resource to deliver
/// to a protocol, and the context that the %JD2ResourceManager uses to make that
/// decision is the job.  E.g., when a protocol asks for the "native" (the Pose of
/// the native structure), the %JD2ResourceManager looks up what job is currently
/// running and then delivers the (single) native Pose for that job.
///
/// The %JD2ResourceManager is meant to work with the JD2ResourceManagerJobInputter.
/// Jobs are defined in an input XML file and for each of the jobs that are defined,
/// resources can be mapped to them.  Resources are also declared in an XML file
/// and may be declared in the same XML file as the Jobs or in a different XML file.
/// The format of this XML file is specified in full between the documentation
/// in the JD2ResourceManagerJobInputter class and three functions that are
/// documented and implemented in this class.
///
/// The %JD2ResourceManager keeps track of which resources are used by which jobs
/// and keeps track of what jobs have completed.  It allocates a Resource the first
/// time that resource is requested by a protoocol, and it deallocates that Resource
/// when all jobs that require that resource have finished.  Thus a protocol that
/// needs an expensive-to-create resource (e.g. a set of 9-mers) can load that
/// resource once and then have that resource shared between all of the jobs that
/// require that resource (e.g. 10K abinitio trajectories) *without having to do
/// any legwork itself.*  If you have 100K abinitio trajectories for 10 different
/// targets, they could all be run in a single job. (Note: as of 9/2013, abinitio
/// does not work with the jd2 or the ResourceManager but it should!)
///
/// The main consequence of relying on the ResourceManager to hand resources to jobs
/// is that jobs which require different resources (usually specified on the command
/// line, and therefore requiring that each job run in a separate process) now can
/// be run together in a single process.  For MPI jobs, this is a major improvement,
/// especially in constructing benchmarks, where many separate, (possibly short)
/// jobs need to be run.
///
/// (Because Resources are not constructed until they are requested, resources can
/// be declared and never used; you can give a resource file that describes every
/// resource you use in every job, and a job declaration file that only uses one
/// of those resources.)
class JD2ResourceManager : public basic::resource_manager::LazyResourceManager
{
public:
	// this class is allowed to instantiate a the JD2 resource manager,
	// but no other class may do so.
	friend class JD2ResourceManagerCreator;

protected:
	JD2ResourceManager();

public:

	virtual
	void
	clear();

	static
	JD2ResourceManager *
	get_jd2_resource_manager_instance();

	virtual ~JD2ResourceManager();

	basic::resource_manager::ResourceOP
	get_resource(
		basic::resource_manager::ResourceDescription const & resource_description);

	bool
	has_resource_with_description(
		basic::resource_manager::ResourceDescription const & resource_description);

	/// @brief Read the portion of an XML file that declares ResourceLocator objects
	void
	read_resource_locators_tags(
		utility::tag::TagCOP tags );

	void
	read_resource_options_tags(
		utility::tag::TagCOP tags );


	void
	read_resources_tags(
		utility::tag::TagCOP tags );

private:
	// Functions to help parsing

	void
	read_resource_options_table_tag(
		utility::tag::TagCOP tag);

	void
	read_resource_option_item(
		utility::tag::TagCOP tag);

	void
	check_resource_loader_type(
		basic::resource_manager::LoaderType const & loader_type);

	basic::resource_manager::ResourceTag
	read_resource_tag_item(
		utility::tag::TagCOP tags,
		basic::resource_manager::LoaderType const & loader_type,
		basic::resource_manager::LocatorID const & locator_id);

	void
	read_resource_table_tag(
		utility::tag::TagCOP tags);


	basic::resource_manager::LocatorTag
	read_resource_locator_items(
		utility::tag::TagCOP tags,
		basic::resource_manager::LoaderType const & loader_type,
		basic::resource_manager::LocatorID & locator_id);

	basic::resource_manager::ResourceOptionsTag
	read_resource_options_tag_item(
		utility::tag::TagCOP tags,
		basic::resource_manager::LoaderType const & loader_type,
		basic::resource_manager::ResourceTag const & resource_tag);

public: // options access

	virtual
	bool
	get_option(
		utility::options::BooleanOptionKey key ) const;

	virtual
	utility::vector1< bool > const &
	get_option(
		utility::options::BooleanVectorOptionKey key ) const;

	virtual
	utility::file::FileName const &
	get_option(
		utility::options::FileOptionKey key ) const;

	virtual
	utility::vector1< utility::file::FileName > const &
	get_option(
		utility::options::FileVectorOptionKey key ) const;

	virtual
	int
	get_option(
		utility::options::IntegerOptionKey key ) const;

	virtual
	utility::vector1< int > const &
	get_option(
		utility::options::IntegerVectorOptionKey key ) const;

	virtual
	utility::file::PathName const &
	get_option(
		utility::options::PathOptionKey key ) const;

	virtual
	utility::vector1< utility::file::PathName > const &
	get_option(
		utility::options::PathVectorOptionKey key ) const;

	virtual
	platform::Real
	get_option(
		utility::options::RealOptionKey key ) const;

	virtual
	utility::vector1< platform::Real > const &
	get_option(
		utility::options::RealVectorOptionKey key ) const;

	virtual
	std::string const &
	get_option(
		utility::options::StringOptionKey key ) const;

	virtual
	utility::vector1< std::string > const &
	get_option(
		utility::options::StringVectorOptionKey key ) const;

	virtual
	bool
	has_option(
		utility::options::BooleanOptionKey key ) const;

	virtual
	bool
	has_option(
		utility::options::BooleanVectorOptionKey key ) const;

	virtual
	bool
	has_option(
		utility::options::FileOptionKey key ) const;

	virtual
	bool
	has_option(
		utility::options::FileVectorOptionKey key ) const;

	virtual
	bool
	has_option(
		utility::options::IntegerOptionKey key ) const;

	virtual
	bool
	has_option(
		utility::options::IntegerVectorOptionKey key ) const;

	virtual
	bool
	has_option(
		utility::options::PathOptionKey key ) const;

	virtual
	bool
	has_option(
		utility::options::PathVectorOptionKey key ) const;

	virtual
	bool
	has_option(
		utility::options::RealOptionKey key ) const;

	virtual
	bool
	has_option(
		utility::options::RealVectorOptionKey key ) const;

	virtual
	bool
	has_option(
		utility::options::StringOptionKey key ) const;

	virtual
	bool
	has_option(
		utility::options::StringVectorOptionKey key ) const;

private:

	basic::resource_manager::ResourceOP
	create_resource_from_fallback(
		basic::resource_manager::FallbackConfigurationCOP fallback,
		basic::resource_manager::ResourceDescription const & resource_description
	);

private:
	std::map< std::string, std::string > fallback_resource_descriptions_created_;

};



} // namespace resource_manager
} // namespace basic

#endif
