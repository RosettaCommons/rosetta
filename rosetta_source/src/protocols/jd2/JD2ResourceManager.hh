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

	void
	read_resource_locators_tags(
		utility::tag::TagPtr tags );

	void
	read_resource_options_tags(
		utility::tag::TagPtr tags );

	void
	read_resources_tags(
		utility::tag::TagPtr tags );


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
