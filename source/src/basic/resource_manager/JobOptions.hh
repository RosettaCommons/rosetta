// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/resource_manager/JobOptions.hh
/// @author Matthew O'Meara (mattjomeara@gmail.com)

#ifndef INCLUDED_basic_resource_manager_JobOptions_hh
#define INCLUDED_basic_resource_manager_JobOptions_hh

// Unit Headers
#include <basic/resource_manager/JobOptions.fwd.hh>

#include <utility/options/OptionCollection.fwd.hh> // for the OptionTypes enum
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
#include <utility/file/FileName.hh>
#include <utility/file/PathName.hh>

// Platform Headers
#include <platform/types.hh>

//C++ Headers
#include <iosfwd>
#include <list>
#include <map>

namespace basic {
namespace resource_manager {

/// @brief The %JobOptions class holds job-specific options (i.e. command line flags).  It can be used
/// by the ResourceManager to hold options for a particular job, so that the ResourceManager can
/// retrieve those options as needed. It is basically a bag for 12 OptionKey/OptionKeyValue
/// maps, one for every kind of OptionKey.
class JobOptions : public utility::pointer::ReferenceCount {

public: // management methods

	JobOptions();

	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~JobOptions();

	/// @brief Describe this %JobOptions to the given output stream
	virtual
	void
	show( std::ostream & out ) const;

public: // accessor methods

	/// @brief Set the value for the indicated boolean option
	void
	add_option(
		utility::options::BooleanOptionKey key,
		bool val);

	/// @brief Unset the value for the indicated boolean option
	void
	remove_option( utility::options::BooleanOptionKey key );

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

	/// @brief Unset the value for the indicated boolean-vector option
	void
	remove_option( utility::options::BooleanVectorOptionKey key );

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

	/// @brief Unset the value for the indicated file option
	void
	remove_option( utility::options::FileOptionKey key );

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

	/// @brief Unset the value for the indicated file-vector option
	void
	remove_option( utility::options::FileVectorOptionKey key );

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

	/// @brief Unset the value for the indicated integer option
	void
	remove_option( utility::options::IntegerOptionKey key );

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

	/// @brief Unset the value for the indicated integer-vector option
	void
	remove_option( utility::options::IntegerVectorOptionKey key );

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

	/// @brief Unset the value for the indicated path option
	void
	remove_option( utility::options::PathOptionKey key );

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

	/// @brief Unset the value for the indicated path-vector option
	void
	remove_option( utility::options::PathVectorOptionKey key );

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

	/// @brief Unset the value for the indicated real option
	void
	remove_option( utility::options::RealOptionKey key );

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

	/// @brief Unset the value for the indicated real-vector option
	void
	remove_option( utility::options::RealVectorOptionKey key );

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

	/// @brief Unset the value for the indicated string option
	void
	remove_option( utility::options::StringOptionKey key );

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

	/// @brief Unset the value for the indicated string-vector option
	void
	remove_option( utility::options::StringVectorOptionKey key );

	/// @brief Return true if a value for the indicated string-vector option has been set
	bool
	has_option(
		utility::options::StringVectorOptionKey key) const;

	/// @brief Return the value of the indicated string-vector option
	utility::vector1< std::string > const &
	get_option(
		utility::options::StringVectorOptionKey key) const;

	bool operator == ( JobOptions const & rhs ) const;

	void track_insertion_order( bool setting );

	std::list< std::pair< utility::options::OptionTypes, std::string > > const &
	insertion_order() const;

private:

	std::map< utility::options::BooleanOptionKey, bool > const &
	map_for_key( utility::options::BooleanOptionKey const & ) const;
	std::map< utility::options::BooleanVectorOptionKey, utility::vector1< bool > > const &
	map_for_key( utility::options::BooleanVectorOptionKey const & ) const;
	std::map< utility::options::FileOptionKey, utility::file::FileName > const &
	map_for_key( utility::options::FileOptionKey const & ) const;
	std::map< utility::options::FileVectorOptionKey, utility::vector1< utility::file::FileName > > const &
	map_for_key( utility::options::FileVectorOptionKey const & ) const;
	std::map< utility::options::IntegerOptionKey, int > const &
	map_for_key( utility::options::IntegerOptionKey const & ) const;
	std::map< utility::options::IntegerVectorOptionKey, utility::vector1< int > > const &
	map_for_key( utility::options::IntegerVectorOptionKey const & ) const;
	std::map< utility::options::PathOptionKey, utility::file::PathName > const &
	map_for_key( utility::options::PathOptionKey const & ) const;
	std::map< utility::options::PathVectorOptionKey, utility::vector1< utility::file::PathName > > const &
	map_for_key( utility::options::PathVectorOptionKey const & ) const;
	std::map< utility::options::RealOptionKey, platform::Real > const &
	map_for_key( utility::options::RealOptionKey const & ) const;
	std::map< utility::options::RealVectorOptionKey, utility::vector1< platform::Real > > const &
	map_for_key( utility::options::RealVectorOptionKey const & ) const;
	std::map< utility::options::StringOptionKey, std::string > const &
	map_for_key( utility::options::StringOptionKey const & ) const;
	std::map< utility::options::StringVectorOptionKey, utility::vector1< std::string > > const &
	map_for_key( utility::options::StringVectorOptionKey const & ) const;

	std::map< utility::options::BooleanOptionKey, bool > &
	map_for_key( utility::options::BooleanOptionKey const & );
	std::map< utility::options::BooleanVectorOptionKey, utility::vector1< bool > > &
	map_for_key( utility::options::BooleanVectorOptionKey const & );
	std::map< utility::options::FileOptionKey, utility::file::FileName > &
	map_for_key( utility::options::FileOptionKey const & );
	std::map< utility::options::FileVectorOptionKey, utility::vector1< utility::file::FileName > > &
	map_for_key( utility::options::FileVectorOptionKey const & );
	std::map< utility::options::IntegerOptionKey, int > &
	map_for_key( utility::options::IntegerOptionKey const & );
	std::map< utility::options::IntegerVectorOptionKey, utility::vector1< int > > &
	map_for_key( utility::options::IntegerVectorOptionKey const & );
	std::map< utility::options::PathOptionKey, utility::file::PathName > &
	map_for_key( utility::options::PathOptionKey const & );
	std::map< utility::options::PathVectorOptionKey, utility::vector1< utility::file::PathName > > &
	map_for_key( utility::options::PathVectorOptionKey const & );
	std::map< utility::options::RealOptionKey, platform::Real > &
	map_for_key( utility::options::RealOptionKey const & );
	std::map< utility::options::RealVectorOptionKey, utility::vector1< platform::Real > > &
	map_for_key( utility::options::RealVectorOptionKey const & );
	std::map< utility::options::StringOptionKey, std::string > &
	map_for_key( utility::options::StringOptionKey const & );
	std::map< utility::options::StringVectorOptionKey, utility::vector1< std::string > > &
	map_for_key( utility::options::StringVectorOptionKey const & );

private:
	bool track_insertion_order_;

	// keep track of the order in which options are added to the JobOptions object
	// so that it can be used to document the options used by a system or application.
	std::list< std::pair< utility::options::OptionTypes, std::string > > insertion_order_;

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

/// @brief This output-operator function invokes the %JobOption's show() method
std::ostream &
operator << (
	std::ostream & out,
	JobOptions const & job_options
);

} // namespace
} // namespace

#endif
