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
#include <iosfwd>
#include <map>

namespace basic {
namespace resource_manager {

/// @brief The %JobOptions class holds job-specific options (i.e. command line flags).  It can be used
/// by the ResourceManager to hold options for a particular job, so that the ResourceManager can
/// retrieve those options as needed. It is basically a bag for 12 OptionKey/OptionKeyValue
/// maps, one for every kind of OptionKey.
class JobOptions : public utility::pointer::ReferenceCount {

public: // management methods
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~JobOptions();

	/// @brief Describe this %JobOptions to the given output stream
	virtual
	void
	show( std::ostream & out ) const;

	/// @brief This friend function output-operator function invokes the %JobOption's show() method
	friend
	std::ostream &
	operator<< (
		std::ostream & out,
		const JobOptions & job_options );

public: // accessor methods

	/// @brief Set the value for the indicated boolean option
	void
	add_option(
		utility::options::BooleanOptionKey key,
		bool val);

	/// @brief Return true if a value for the indicated boolean option has been set
	bool
	has_option(
		utility::options::BooleanOptionKey key) const;

	/// @brief Return the value of the indicated boolean option
	bool
	get_option(
		utility::options::BooleanOptionKey key) const;


	/// @brief Set the value for the indicated boolean-vector option
	void
	add_option(
		utility::options::BooleanVectorOptionKey key,
		utility::vector1< bool > const & val);

	/// @brief Return true if a value for the indicated boolean-vector option has been set
	bool
	has_option(
		utility::options::BooleanVectorOptionKey key) const;

	/// @brief Return the value of the indicated boolean-vector option
	utility::vector1< bool > const &
	get_option(
		utility::options::BooleanVectorOptionKey key) const;


	/// @brief Set the value for the indicated file option
	void
	add_option(
		utility::options::FileOptionKey key,
		utility::file::FileName const & val);

	/// @brief Return true if a value for the indicated file option has been set
	bool
	has_option(
		utility::options::FileOptionKey key) const;

	/// @brief Return the value of the indicated file option
	utility::file::FileName const &
	get_option(
		utility::options::FileOptionKey key) const;

	/// @brief Set the value for the indicated file-vector option
	void
	add_option(
		utility::options::FileVectorOptionKey key,
		utility::vector1< utility::file::FileName > const & val);

	/// @brief Return true if a value for the indicated file-vector option has been set
	bool
	has_option(
		utility::options::FileVectorOptionKey key) const;

	/// @brief Return the value of the indicated file-vector option
	utility::vector1< utility::file::FileName > const &
	get_option(
		utility::options::FileVectorOptionKey key) const;


	/// @brief Set the value for the indicated integer option
	void
	add_option(
		utility::options::IntegerOptionKey key,
		int val);

	/// @brief Return true if a value for the indicated integer option has been set
	bool
	has_option(
		utility::options::IntegerOptionKey key) const;

	/// @brief Return the value of the indicated integer option
	int
	get_option(
		utility::options::IntegerOptionKey key) const;


	/// @brief Set the value for the indicated integer-vector option
	void
	add_option(
		utility::options::IntegerVectorOptionKey key,
		utility::vector1< int > const & val);

	/// @brief Return true if a value for the indicated integer-vector option has been set
	bool
	has_option(
		utility::options::IntegerVectorOptionKey key) const;

	/// @brief Return the value of the indicated integer-vector option
	utility::vector1< int > const &
	get_option(
		utility::options::IntegerVectorOptionKey key) const;


	/// @brief Set the value for the indicated path option
	void
	add_option(
		utility::options::PathOptionKey key,
		utility::file::PathName const & val);

	/// @brief Return true if a value for the indicated path option has been set
	bool
	has_option(
		utility::options::PathOptionKey key) const;

	/// @brief Return the value of the indicated path option
	utility::file::PathName const &
	get_option(
		utility::options::PathOptionKey key) const;


	/// @brief Set the value for the indicated path-vector option
	void
	add_option(
		utility::options::PathVectorOptionKey key,
		utility::vector1< utility::file::PathName > const & val);

	/// @brief Return true if a value for the indicated path-vector option has been set
	bool
	has_option(
		utility::options::PathVectorOptionKey key) const;

	/// @brief Return the value of the indicated path-vector option
	utility::vector1< utility::file::PathName > const &
	get_option(
		utility::options::PathVectorOptionKey key) const;


	/// @brief Set the value for the indicated real option
	void
	add_option(
		utility::options::RealOptionKey key,
		platform::Real val);

	/// @brief Return true if a value for the indicated real option has been set
	bool
	has_option(
		utility::options::RealOptionKey key) const;

	/// @brief Return the value of the indicated real option
	platform::Real
	get_option(
		utility::options::RealOptionKey key) const;

	/// @brief Set the value for the indicated real-vector option
	void
	add_option(
		utility::options::RealVectorOptionKey key,
		utility::vector1< platform::Real > const & val);

	/// @brief Return true if a value for the indicated real-vector option has been set
	bool
	has_option(
		utility::options::RealVectorOptionKey key) const;

	/// @brief Return the value of the indicated real-vector option
	utility::vector1< platform::Real > const &
	get_option(
		utility::options::RealVectorOptionKey key) const;

	/// @brief Set the value for the indicated string option
	void
	add_option(
		utility::options::StringOptionKey key,
		std::string val);

	/// @brief Return true if a value for the indicated string option has been set
	bool
	has_option(
		utility::options::StringOptionKey key) const;

	/// @brief Return the value of the indicated string option
	std::string const &
	get_option(
		utility::options::StringOptionKey key) const;

	/// @brief Set the value for the indicated string-vector option
	void
	add_option(
		utility::options::StringVectorOptionKey key,
		utility::vector1< std::string > const & val);

	/// @brief Return true if a value for the indicated string-vector option has been set
	bool
	has_option(
		utility::options::StringVectorOptionKey key) const;

	/// @brief Return the value of the indicated string-vector option
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
