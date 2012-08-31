// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/resource_manager/JobOptions.hh
/// @author Matthew O'Meara (mattjomeara@gmail.com)

#ifndef INCLUDED_basic_resource_manager_JobOptions_hh
#define INCLUDED_basic_resource_manager_JobOptions_hh

// Unit Headers
#include <basic/resource_manager/JobOptions.fwd.hh>

#include <utility/options/keys/BooleanOptionKey.hh>
#include <utility/options/keys/BooleanVectorOptionKey.hh>
#include <utility/options/keys/FileOptionKey.hh>
#include <utility/options/keys/FileVectorOptionKey.hh>
#include <utility/options/keys/IntegerOptionKey.hh>
#include <utility/options/keys/IntegerVectorOptionKey.hh>
#include <utility/options/keys/OptionKey.hh>
#include <utility/options/keys/PathOptionKey.hh>
#include <utility/options/keys/PathVectorOptionKey.hh>
#include <utility/options/keys/RealOptionKey.hh>
#include <utility/options/keys/RealVectorOptionKey.hh>
#include <utility/options/keys/StringOptionKey.hh>
#include <utility/options/keys/StringVectorOptionKey.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <utility/file/FileName.fwd.hh>
#include <utility/file/PathName.fwd.hh>

// Platform Headers
#include <platform/types.hh>

//C++ Headers
#include <string>
#include <map>

namespace basic {
namespace resource_manager {

/// @brief This class will hold job-specific options.  It can be used
/// by a ResourceManager to hold options for a particular job, so that
/// the ResourceManager can retrieve those options as needed.
/// It is basically a bag for 12 OptionKey/OptionKeyValue maps, one
/// for every kind of OptionKey.
class JobOptions : public utility::pointer::ReferenceCount {

public: // management methods
	///@brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~JobOptions();

	virtual
	void
	show( std::ostream & out ) const;

	friend
	std::ostream &
	operator<< (
		std::ostream & out,
		const JobOptions & job_options );



public: // accessor methods

	void
	add_option(
		utility::options::BooleanOptionKey key,
		bool val);

	bool
	has_option(
		utility::options::BooleanOptionKey key) const;

	bool
	get_option(
		utility::options::BooleanOptionKey key) const;


	void
	add_option(
		utility::options::BooleanVectorOptionKey key,
		utility::vector1< bool > const & val);

	bool
	has_option(
		utility::options::BooleanVectorOptionKey key) const;

	utility::vector1< bool > const &
	get_option(
		utility::options::BooleanVectorOptionKey key) const;


	void
	add_option(
		utility::options::FileOptionKey key,
		utility::file::FileName const & val);

	bool
	has_option(
		utility::options::FileOptionKey key) const;

	utility::file::FileName const &
	get_option(
		utility::options::FileOptionKey key) const;

	void
	add_option(
		utility::options::FileVectorOptionKey key,
		utility::vector1< utility::file::FileName > const & val);

	bool
	has_option(
		utility::options::FileVectorOptionKey key) const;

	utility::vector1< utility::file::FileName > const &
	get_option(
		utility::options::FileVectorOptionKey key) const;


	void
	add_option(
		utility::options::IntegerOptionKey key,
		int val);

	bool
	has_option(
		utility::options::IntegerOptionKey key) const;

	int
	get_option(
		utility::options::IntegerOptionKey key) const;


	void
	add_option(
		utility::options::IntegerVectorOptionKey key,
		utility::vector1< int > const & val);

	bool
	has_option(
		utility::options::IntegerVectorOptionKey key) const;

	utility::vector1< int > const &
	get_option(
		utility::options::IntegerVectorOptionKey key) const;


	void
	add_option(
		utility::options::PathOptionKey key,
		utility::file::PathName const & val);

	bool
	has_option(
		utility::options::PathOptionKey key) const;

	utility::file::PathName const &
	get_option(
		utility::options::PathOptionKey key) const;


	void
	add_option(
		utility::options::PathVectorOptionKey key,
		utility::vector1< utility::file::PathName > const & val);

	bool
	has_option(
		utility::options::PathVectorOptionKey key) const;

	utility::vector1< utility::file::PathName > const &
	get_option(
		utility::options::PathVectorOptionKey key) const;


	void
	add_option(
		utility::options::RealOptionKey key,
		platform::Real val);

	bool
	has_option(
		utility::options::RealOptionKey key) const;

	platform::Real
	get_option(
		utility::options::RealOptionKey key) const;


	void
	add_option(
		utility::options::RealVectorOptionKey key,
		utility::vector1< platform::Real > const & val);

	bool
	has_option(
		utility::options::RealVectorOptionKey key) const;

	utility::vector1< platform::Real > const &
	get_option(
		utility::options::RealVectorOptionKey key) const;


	void
	add_option(
		utility::options::StringOptionKey key,
		std::string val);

	bool
	has_option(
		utility::options::StringOptionKey key) const;

	std::string const &
	get_option(
		utility::options::StringOptionKey key) const;


	void
	add_option(
		utility::options::StringVectorOptionKey key,
		utility::vector1< std::string > const & val);

	bool
	has_option(
		utility::options::StringVectorOptionKey key) const;

	utility::vector1< std::string > const &
	get_option(
		utility::options::StringVectorOptionKey key) const;

private:
	std::map< utility::options::BooleanOptionKey, bool >
		boolean_options_;
	std::map< utility::options::BooleanVectorOptionKey, utility::vector1< bool > >
		boolean_vector_options_;

	std::map< utility::options::FileOptionKey, utility::file::FileName >
		file_options_;
	std::map< utility::options::FileVectorOptionKey, utility::vector1< utility::file::FileName > >
		file_vector_options_;

	std::map< utility::options::IntegerOptionKey, int >
		integer_options_;
	std::map< utility::options::IntegerVectorOptionKey, utility::vector1< int > >
		integer_vector_options_;

	std::map< utility::options::PathOptionKey, utility::file::PathName >
		path_options_;
	std::map< utility::options::PathVectorOptionKey, utility::vector1< utility::file::PathName > >
		path_vector_options_;

	std::map< utility::options::RealOptionKey, platform::Real >
		real_options_;
	std::map< utility::options::RealVectorOptionKey, utility::vector1< platform::Real > >
		real_vector_options_;

	std::map< utility::options::StringOptionKey, std::string >
		string_options_;
	std::map< utility::options::StringVectorOptionKey, utility::vector1< std::string > >
		string_vector_options_;

};

} // namespace
} // namespace

#endif
