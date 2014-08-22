// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/options/OptionCollection.hh
/// @brief  Program options collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_utility_options_OptionCollection_hh
#define INCLUDED_utility_options_OptionCollection_hh


// Unit headers
#include <utility/options/OptionCollection.fwd.hh>

// Package headers
#include <utility/options/Option.hh>
#include <utility/options/BooleanOption.hh>
#include <utility/options/IntegerOption.hh>
#include <utility/options/RealOption.hh>
#include <utility/options/StringOption.hh>
#include <utility/options/FileOption.hh>
#include <utility/options/PathOption.hh>
#include <utility/options/AnyOption.hh>
#include <utility/options/BooleanVectorOption.hh>
#include <utility/options/IntegerVectorOption.hh>
#include <utility/options/RealVectorOption.hh>
#include <utility/options/ResidueChainVectorOption.hh>
#include <utility/options/StringVectorOption.hh>
#include <utility/options/FileVectorOption.hh>
#include <utility/options/PathVectorOption.hh>
#include <utility/options/AnyVectorOption.hh>
#include <utility/options/VariantOption.hh>
#include <utility/options/keys/BooleanOptionKey.hh>
#include <utility/options/keys/IntegerOptionKey.hh>
#include <utility/options/keys/RealOptionKey.hh>
#include <utility/options/keys/StringOptionKey.hh>
#include <utility/options/keys/FileOptionKey.hh>
#include <utility/options/keys/PathOptionKey.hh>
#include <utility/options/keys/AnyOptionKey.hh>
#include <utility/options/keys/BooleanVectorOptionKey.hh>
#include <utility/options/keys/IntegerVectorOptionKey.hh>
#include <utility/options/keys/ResidueChainVectorOptionKey.hh>
#include <utility/options/keys/RealVectorOptionKey.hh>
#include <utility/options/keys/StringVectorOptionKey.hh>
#include <utility/options/keys/FileVectorOptionKey.hh>
#include <utility/options/keys/PathVectorOptionKey.hh>
#include <utility/options/keys/AnyVectorOptionKey.hh>

// Project headers
#include <utility/down_cast.hh>
#include <utility/keys/SmallKeyVector.hh>
#include <utility/exit.hh>

// C++ headers
#include <stdlib.h>
#include <cassert>
#include <iosfwd>
#include <list>
#include <string>
#include <vector>

namespace utility {
namespace options {


/// @brief Program options collection
class OptionCollection
{


private: // Types


	typedef  std::list< std::string >  ValueStrings;

	/// @brief Option types enumeration
	enum OptionTypes {
		UNKNOWN_OPTION,
		BOOLEAN_OPTION,
		INTEGER_OPTION,
		REAL_OPTION,
		STRING_OPTION,
		FILE_OPTION,
		PATH_OPTION,
		ANY_OPTION,
		BOOLEAN_VECTOR_OPTION,
		INTEGER_VECTOR_OPTION,
		REAL_VECTOR_OPTION,
		RESIDUE_CHAIN_VECTOR_OPTION,
		STRING_VECTOR_OPTION,
		FILE_VECTOR_OPTION,
		PATH_VECTOR_OPTION,
		ANY_VECTOR_OPTION
	};

	// Note: These waste a small amount of space for the sake of constant-time option lookup
	typedef  utility::keys::SmallKeyVector< BooleanOptionKey, BooleanOption >  Booleans;
	typedef  utility::keys::SmallKeyVector< IntegerOptionKey, IntegerOption >  Integers;
	typedef  utility::keys::SmallKeyVector< RealOptionKey, RealOption >  Reals;
	typedef  utility::keys::SmallKeyVector< StringOptionKey, StringOption >  Strings;
	typedef  utility::keys::SmallKeyVector< FileOptionKey, FileOption >  Files;
	typedef  utility::keys::SmallKeyVector< PathOptionKey, PathOption >  Paths;
	typedef  utility::keys::SmallKeyVector< AnyOptionKey, VariantOption< ScalarOption > >  Anys;
	typedef  utility::keys::SmallKeyVector< BooleanVectorOptionKey, BooleanVectorOption >  BooleanVectors;
	typedef  utility::keys::SmallKeyVector< IntegerVectorOptionKey, IntegerVectorOption >  IntegerVectors;
	typedef  utility::keys::SmallKeyVector< RealVectorOptionKey, RealVectorOption >  RealVectors;
	typedef  utility::keys::SmallKeyVector< ResidueChainVectorOptionKey, ResidueChainVectorOption >  ResidueChainVectors;
	typedef  utility::keys::SmallKeyVector< StringVectorOptionKey, StringVectorOption >  StringVectors;
	typedef  utility::keys::SmallKeyVector< FileVectorOptionKey, FileVectorOption >  FileVectors;
	typedef  utility::keys::SmallKeyVector< PathVectorOptionKey, PathVectorOption >  PathVectors;
	typedef  utility::keys::SmallKeyVector< AnyVectorOptionKey, VariantOption< VectorOption > >  AnyVectors;
	typedef  utility::keys::SmallKeyVector< OptionKey, OptionTypes >  All;


public: // Creation


	/// @brief Default constructor
	inline
	OptionCollection()
		: argv_copy_("")
	{}


	/// @brief Destructor
	~OptionCollection();


public: // Methods


	/// @brief Add the built-in options
	void
	add_built_in_options();

	/// @brief Add a BooleanOption
	inline
	BooleanOption &
	add( BooleanOptionKey const & key, std::string const & description )
	{
		check_key( key );
		all_( key ) = BOOLEAN_OPTION;
		return ( booleans_( key ) = BooleanOption( key, description ) );
	}


	/// @brief Add an IntegerOption
	inline
	IntegerOption &
	add( IntegerOptionKey const & key, std::string const & description )
	{
		check_key( key );
		all_( key ) = INTEGER_OPTION;
		return ( integers_( key ) = IntegerOption( key, description ) );
	}


	/// @brief Add a RealOption
	inline
	RealOption &
	add( RealOptionKey const & key, std::string const & description )
	{
		check_key( key );
		all_( key ) = REAL_OPTION;
		return ( reals_( key ) = RealOption( key, description ) );
	}


	/// @brief Add a StringOption
	inline
	StringOption &
	add( StringOptionKey const & key, std::string const & description )
	{
		check_key( key );
		all_( key ) = STRING_OPTION;
		return ( strings_( key ) = StringOption( key, description ) );
	}


	/// @brief Add a FileOption
	inline
	FileOption &
	add( FileOptionKey const & key, std::string const & description )
	{
		check_key( key );
		all_( key ) = FILE_OPTION;
		return ( files_( key ) = FileOption( key, description ) );
	}


	/// @brief Add a PathOption
	inline
	PathOption &
	add( PathOptionKey const & key, std::string const & description )
	{
		check_key( key );
		all_( key ) = PATH_OPTION;
		return ( paths_( key ) = PathOption( key, description ) );
	}


	/// @brief Add a BooleanVectorOption
	inline
	BooleanVectorOption &
	add( BooleanVectorOptionKey const & key, std::string const & description )
	{
		check_key( key );
		all_( key ) = BOOLEAN_VECTOR_OPTION;
		return ( boolean_vectors_( key ) = BooleanVectorOption( key, description ) );
	}


	/// @brief Add an IntegerVectorOption
	inline
	IntegerVectorOption &
	add( IntegerVectorOptionKey const & key, std::string const & description )
	{
		check_key( key );
		all_( key ) = INTEGER_VECTOR_OPTION;
		return ( integer_vectors_( key ) = IntegerVectorOption( key, description ) );
	}


	/// @brief Add a RealVectorOption
	inline
	RealVectorOption &
	add( RealVectorOptionKey const & key, std::string const & description )
	{
		check_key( key );
		all_( key ) = REAL_VECTOR_OPTION;
		return ( real_vectors_( key ) = RealVectorOption( key, description ) );
	}

	/// @brief Add an ResidueChainVectorOption
	inline
	ResidueChainVectorOption &
	add( ResidueChainVectorOptionKey const & key, std::string const & description )
	{
		check_key( key );
		all_( key ) = RESIDUE_CHAIN_VECTOR_OPTION;
		return ( residue_chain_vectors_( key ) = ResidueChainVectorOption( key, description ) );
	}


	/// @brief Add a StringVectorOption
	inline
	StringVectorOption &
	add( StringVectorOptionKey const & key, std::string const & description )
	{
		check_key( key );
		all_( key ) = STRING_VECTOR_OPTION;
		return ( string_vectors_( key ) = StringVectorOption( key, description ) );
	}


	/// @brief Add a FileVectorOption
	inline
	FileVectorOption &
	add( FileVectorOptionKey const & key, std::string const & description )
	{
		check_key( key );
		all_( key ) = FILE_VECTOR_OPTION;
		return ( file_vectors_( key ) = FileVectorOption( key, description ) );
	}


	/// @brief Add a PathVectorOption
	inline
	PathVectorOption &
	add( PathVectorOptionKey const & key, std::string const & description )
	{
		check_key( key );
		all_( key ) = PATH_VECTOR_OPTION;
		return ( path_vectors_( key ) = PathVectorOption( key, description ) );
	}


	/// @brief Add a BooleanOption
	inline
	BooleanOption &
	add( BooleanOption const & option )
	{
		check_key( option );
		all_( option.key() ) = BOOLEAN_OPTION;
		return ( booleans_( option.key() ) = option );
	}


	/// @brief Add an IntegerOption
	inline
	IntegerOption &
	add( IntegerOption const & option )
	{
		check_key( option );
		all_( option.key() ) = INTEGER_OPTION;
		return ( integers_( option.key() ) = option );
	}


	/// @brief Add a RealOption
	inline
	RealOption &
	add( RealOption const & option )
	{
		check_key( option );
		all_( option.key() ) = REAL_OPTION;
		return ( reals_( option.key() ) = option );
	}


	/// @brief Add a StringOption
	inline
	StringOption &
	add( StringOption const & option )
	{
		check_key( option );
		all_( option.key() ) = STRING_OPTION;
		return ( strings_( option.key() ) = option );
	}


	/// @brief Add a FileOption
	inline
	FileOption &
	add( FileOption const & option )
	{
		check_key( option );
		all_( option.key() ) = FILE_OPTION;
		return ( files_( option.key() ) = option );
	}


	/// @brief Add a PathOption
	inline
	PathOption &
	add( PathOption const & option )
	{
		check_key( option );
		all_( option.key() ) = PATH_OPTION;
		return ( paths_( option.key() ) = option );
	}


	/// @brief Add an AnyOption
	template< typename T >
	inline
	AnyOption< T > &
	add( AnyOption< T > const & option )
	{
		check_key( option );
		all_( option.key() ) = ANY_OPTION;
		return ( anys_( option.key() ) = option );
	}


	/// @brief Add a BooleanVectorOption
	inline
	BooleanVectorOption &
	add( BooleanVectorOption const & option )
	{
		check_key( option );
		all_( option.key() ) = BOOLEAN_VECTOR_OPTION;
		return ( boolean_vectors_( option.key() ) = option );
	}


	/// @brief Add an IntegerVectorOption
	inline
	IntegerVectorOption &
	add( IntegerVectorOption const & option )
	{
		check_key( option );
		all_( option.key() ) = INTEGER_VECTOR_OPTION;
		return ( integer_vectors_( option.key() ) = option );
	}


	/// @brief Add a RealVectorOption
	inline
	RealVectorOption &
	add( RealVectorOption const & option )
	{
		check_key( option );
		all_( option.key() ) = REAL_VECTOR_OPTION;
		return ( real_vectors_( option.key() ) = option );
	}

	/// @brief Add an ResidueChainVectorOption
	inline
	ResidueChainVectorOption &
	add( ResidueChainVectorOption const & option )
	{
		check_key( option );
		all_( option.key() ) = RESIDUE_CHAIN_VECTOR_OPTION;
		return ( residue_chain_vectors_( option.key() ) = option );
	}


	/// @brief Add a StringVectorOption
	inline
	StringVectorOption &
	add( StringVectorOption const & option )
	{
		check_key( option );
		all_( option.key() ) = STRING_VECTOR_OPTION;
		return ( string_vectors_( option.key() ) = option );
	}


	/// @brief Add a FileVectorOption
	inline
	FileVectorOption &
	add( FileVectorOption const & option )
	{
		check_key( option );
		all_( option.key() ) = FILE_VECTOR_OPTION;
		return ( file_vectors_( option.key() ) = option );
	}


	/// @brief Add a PathVectorOption
	inline
	PathVectorOption &
	add( PathVectorOption const & option )
	{
		check_key( option );
		all_( option.key() ) = PATH_VECTOR_OPTION;
		return ( path_vectors_( option.key() ) = option );
	}


	/// @brief Add an AnyVectorOption
	template< typename T >
	inline
	AnyVectorOption< T > &
	add( AnyVectorOption< T > const & option )
	{
		check_key( option );
		all_( option.key() ) = ANY_VECTOR_OPTION;
		return ( any_vectors_( option.key() ) = option );
	}


//	/// @brief Specify mutually exclusive options
//	void
//	exclusive(
//		OptionKey const & key1,
//		OptionKey const & key2
//	);


//	/// @brief Specify an option that requires a second option to also be specified
//	void
//	requires(
//		OptionKey const & key1,
//		OptionKey const & key2
//	);

	/// @brief Check for problems in the option specifications
	void
	check_specs() const;

	/// @brief Load the user-specified option values
	void
	load(
		const std::vector< std::string> & args,
		bool const free_args // Support free argument (without - prefix)?
	);

	/// @brief Load the user-specified option values
	void
	load(
		int const argc,
		char * const argv[],
		bool const free_args = false // Support free argument (without - prefix)?
	);

	/// @brief Load the user-specified option values
  void
  load(
		std::string executable_name, // usually argv[ 0 ]
		ValueStrings& arg_strings,
    bool const free_args // Support free argument (without - prefix)?
	);

	///@brief Load all options in a flags file
	void load_options_from_file(
		std::string const & file_string,
		std::string const & cid=""
	);

	///@brief same as load_options_from_file, but throws exception instead of call to std::exit
	void load_options_from_file_exception(
		std::string const & file_string,
		std::string const & cid=""
	);

	///@brief Load all options in a flags file
	void load_options_from_stream(
    std::istream& stream,
		std::string const & file_string="STREAM", //for error msg, set if stream is connected to file_string
		std::string const & cid=""
	);

	/// @brief Check for problems in the option values
	void
	check_values() const;


	/// @brief Show all the options and their descriptions
	void
	show_help( std::ostream & stream );


	/// @brief Show one option and it description
	void show_option_help(OptionKey const &, std::string &group, std::ostream & stream );


	/// @brief Show all the options and their descriptions in a hierarchy format
	void
	show_help_hier( std::ostream & stream );

	/// @brief Show one option and it description in a hierarchy format
	void show_option_help_heir(OptionKey const &, std::string &group, std::ostream & stream );


	/// @brief Show the user-specified options and their values
	void
	show_user( std::ostream & stream ) const;


	/// @brief Show all the options and their values
	void
	show_all( std::ostream & stream ) const;


	/// @brief Show all the options and their values in a hierarchy format
	void
	show_all_hier( std::ostream & stream ) const;


	/// @brief Show the options definitions table in text format
	void
	show_table_text( std::ostream & stream ) const;


	/// @brief Show the options definitions table in Wiki format
	void
	show_table_Wiki( std::ostream & stream ) const;


	/// @brief Show accessed list of options.
	void show_accessed_options(std::ostream & stream) const;

	/// @brief Show inaccessed user-specified options
	void show_inaccessed_user_options(std::ostream & stream) const;

	/// @brief Show options that were user-specified, but not yet accessed.
	/// Differs from show_inaccessed_user_options() mainly by formatting
	/// And how the options are accessed. (This one should be destructor safe.)
	void show_unused_options(std::ostream & stream) const;

	/// @brief Returns a copy of the concatenated argv strings that were initialized
	/// in load().
	std::string get_argv() const;

public: // Properties: predicate


	/// @brief Is there an option with a BooleanOptionKey?
	inline
	bool
	has( BooleanOptionKey const & key ) const
	{
		return booleans_.has( key );
	}


	/// @brief Is there an option with an IntegerOptionKey?
	inline
	bool
	has( IntegerOptionKey const & key ) const
	{
		return integers_.has( key );
	}


	/// @brief Is there an option with a RealOptionKey?
	inline
	bool
	has( RealOptionKey const & key ) const
	{
		return reals_.has( key );
	}


	/// @brief Is there an option with a StringOptionKey?
	inline
	bool
	has( StringOptionKey const & key ) const
	{
		return strings_.has( key );
	}


	/// @brief Is there an option with a FileOptionKey?
	inline
	bool
	has( FileOptionKey const & key ) const
	{
		return files_.has( key );
	}


	/// @brief Is there an option with a PathOptionKey?
	inline
	bool
	has( PathOptionKey const & key ) const
	{
		return paths_.has( key );
	}


	/// @brief Is there an option with an AnyOptionKey?
	inline
	bool
	has( AnyOptionKey const & key ) const
	{
		return anys_.has( key );
	}


	/// @brief Is there an option with a BooleanVectorOptionKey?
	inline
	bool
	has( BooleanVectorOptionKey const & key ) const
	{
		return boolean_vectors_.has( key );
	}


	/// @brief Is there an option with an IntegerVectorOptionKey?
	inline
	bool
	has( IntegerVectorOptionKey const & key ) const
	{
		return integer_vectors_.has( key );
	}


	/// @brief Is there an option with a RealVectorOptionKey?
	inline
	bool
	has( RealVectorOptionKey const & key ) const
	{
		return real_vectors_.has( key );
	}

	/// @brief Is there an option with an IntegerVectorOptionKey?
	inline
	bool
	has( ResidueChainVectorOptionKey const & key ) const
	{
		return residue_chain_vectors_.has( key );
	}


	/// @brief Is there an option with a StringVectorOptionKey?
	inline
	bool
	has( StringVectorOptionKey const & key ) const
	{
		return string_vectors_.has( key );
	}


	/// @brief Is there an option with a FileVectorOptionKey?
	inline
	bool
	has( FileVectorOptionKey const & key ) const
	{
		return file_vectors_.has( key );
	}


	/// @brief Is there an option with a PathVectorOptionKey?
	inline
	bool
	has( PathVectorOptionKey const & key ) const
	{
		return path_vectors_.has( key );
	}


	/// @brief Is there an option with an AnyVectorOptionKey?
	inline
	bool
	has( AnyVectorOptionKey const & key ) const
	{
		return any_vectors_.has( key );
	}


	/// @brief Is there an option with an OptionKey?
	inline
	bool
	has( OptionKey const & key ) const
	{
		return all_.has( key );
	}


public: // Properties


	/// @brief Option by BooleanOptionKey
	inline
	BooleanOption const &
	option( BooleanOptionKey const & key ) const
	{

		return operator []( key );
	}


	/// @brief Option by BooleanOptionKey
	inline
	BooleanOption &
	option( BooleanOptionKey const & key )
	{
		return operator []( key );
	}


	/// @brief Option by IntegerOptionKey
	inline
	IntegerOption const &
	option( IntegerOptionKey const & key ) const
	{
		return operator []( key );
	}


	/// @brief Option by IntegerOptionKey
	inline
	IntegerOption &
	option( IntegerOptionKey const & key )
	{
		return operator []( key );
	}


	/// @brief Option by RealOptionKey
	inline
	RealOption const &
	option( RealOptionKey const & key ) const
	{
		return operator []( key );
	}


	/// @brief Option by RealOptionKey
	inline
	RealOption &
	option( RealOptionKey const & key )
	{
		return operator []( key );
	}


	/// @brief Option by StringOptionKey
	inline
	StringOption const &
	option( StringOptionKey const & key ) const
	{
		return operator []( key );
	}


	/// @brief Option by StringOptionKey
	inline
	StringOption &
	option( StringOptionKey const & key )
	{
		return operator []( key );
	}


	/// @brief Option by FileOptionKey
	inline
	FileOption const &
	option( FileOptionKey const & key ) const
	{
		return operator []( key );
	}


	/// @brief Option by FileOptionKey
	inline
	FileOption &
	option( FileOptionKey const & key )
	{
		return operator []( key );
	}


	/// @brief Option by PathOptionKey
	inline
	PathOption const &
	option( PathOptionKey const & key ) const
	{
		return operator []( key );
	}


	/// @brief Option by PathOptionKey
	inline
	PathOption &
	option( PathOptionKey const & key )
	{
		return operator []( key );
	}


	/// @brief Option by AnyOptionKey
	inline
	Option const &
	option( AnyOptionKey const & key ) const
	{
		return operator []( key );
	}


	/// @brief Option by AnyOptionKey
	inline
	Option &
	option( AnyOptionKey const & key )
	{
		return operator []( key );
	}


	/// @brief Option by AnyOptionKey with option type template argument
	template< typename OptionType >
	inline
	OptionType const &
	option( AnyOptionKey const & key ) const
	{
		return operator []< OptionType >( key );
	}


	/// @brief Option by AnyOptionKey with option type template argument
	template< typename OptionType >
	inline
	OptionType &
	option( AnyOptionKey const & key )
	{
		return operator []< OptionType >( key );
	}


	/// @brief VectorOption by BooleanVectorOptionKey
	inline
	BooleanVectorOption const &
	option( BooleanVectorOptionKey const & key ) const
	{
		return operator []( key );
	}


	/// @brief VectorOption by BooleanVectorOptionKey
	inline
	BooleanVectorOption &
	option( BooleanVectorOptionKey const & key )
	{
		return operator []( key );
	}


	/// @brief VectorOption by IntegerVectorOptionKey
	inline
	IntegerVectorOption const &
	option( IntegerVectorOptionKey const & key ) const
	{
		return operator []( key );
	}


	/// @brief VectorOption by IntegerVectorOptionKey
	inline
	IntegerVectorOption &
	option( IntegerVectorOptionKey const & key )
	{
		return operator []( key );
	}


	/// @brief VectorOption by RealVectorOptionKey
	inline
	RealVectorOption const &
	option( RealVectorOptionKey const & key ) const
	{
		return operator []( key );
	}


	/// @brief VectorOption by RealVectorOptionKey
	inline
	RealVectorOption &
	option( RealVectorOptionKey const & key )
	{
		return operator []( key );
	}


	/// @brief VectorOption by StringVectorOptionKey
	inline
	StringVectorOption const &
	option( StringVectorOptionKey const & key ) const
	{
		return operator []( key );
	}


	/// @brief VectorOption by StringVectorOptionKey
	inline
	StringVectorOption &
	option( StringVectorOptionKey const & key )
	{
		return operator []( key );
	}


	/// @brief VectorOption by FileVectorOptionKey
	inline
	FileVectorOption const &
	option( FileVectorOptionKey const & key ) const
	{
		return operator []( key );
	}


	/// @brief VectorOption by FileVectorOptionKey
	inline
	FileVectorOption &
	option( FileVectorOptionKey const & key )
	{
		return operator []( key );
	}


	/// @brief VectorOption by PathVectorOptionKey
	inline
	PathVectorOption const &
	option( PathVectorOptionKey const & key ) const
	{
		return operator []( key );
	}


	/// @brief VectorOption by PathVectorOptionKey
	inline
	PathVectorOption &
	option( PathVectorOptionKey const & key )
	{
		return operator []( key );
	}


	/// @brief VectorOption by AnyVectorOptionKey
	inline
	VectorOption const &
	option( AnyVectorOptionKey const & key ) const
	{
		return operator []( key );
	}


	/// @brief VectorOption by AnyVectorOptionKey
	inline
	VectorOption &
	option( AnyVectorOptionKey const & key )
	{
		return operator []( key );
	}


	/// @brief VectorOption by AnyVectorOptionKey with option type template argument
	template< typename VectorOptionType >
	inline
	VectorOptionType const &
	option( AnyVectorOptionKey const & key ) const
	{
		return operator []< VectorOptionType >( key );
	}


	/// @brief VectorOption by AnyVectorOptionKey with option type template argument
	template< typename VectorOptionType >
	inline
	VectorOptionType &
	option( AnyVectorOptionKey const & key )
	{
		return operator []< VectorOptionType >( key );
	}


	/// @brief Option by OptionKey
	inline
	Option const &
	option( OptionKey const & key ) const
	{
		return operator []( key );
	}


	/// @brief Option by OptionKey
	inline
	Option &
	option( OptionKey const & key )
	{
		return operator []( key );
	}


	/// @brief Option by OptionKey with option type template argument
	template< typename OptionType >
	inline
	OptionType const &
	option( OptionKey const & key ) const
	{
		return operator []< OptionType >( key );
	}


	/// @brief Option by OptionKey with option type template argument
	template< typename OptionType >
	inline
	OptionType &
	option( OptionKey const & key )
	{
		return operator []< OptionType >( key );
	}


public: // Indexers


	/// @brief OptionCollection[ BooleanOptionKey ]
	inline
	BooleanOption const &
	operator []( BooleanOptionKey const & key ) const
	{
		BooleanOption const & option( booleans_[ key ] );
		option.check_restricted_access(true);
		return option;
	}


	/// @brief OptionCollection[ BooleanOptionKey ]
	inline
	BooleanOption &
	operator []( BooleanOptionKey const & key )
	{
		BooleanOption & option( booleans_[ key ] );
		option.check_restricted_access(true);
		return option;
	}


	/// @brief OptionCollection[ IntegerOptionKey ]
	inline
	IntegerOption const &
	operator []( IntegerOptionKey const & key ) const
	{
		IntegerOption const & option( integers_[ key ] );
		option.check_restricted_access(true);
		return option;
	}


	/// @brief OptionCollection[ IntegerOptionKey ]
	inline
	IntegerOption &
	operator []( IntegerOptionKey const & key )
	{
		IntegerOption & option( integers_[ key ] );
		option.check_restricted_access(true);
		return option;
	}


	/// @brief OptionCollection[ RealOptionKey ]
	inline
	RealOption const &
	operator []( RealOptionKey const & key ) const
	{
		RealOption const & option( reals_[ key ] );
		option.check_restricted_access(true);
		return option;
	}


	/// @brief OptionCollection[ RealOptionKey ]
	inline
	RealOption &
	operator []( RealOptionKey const & key )
	{
		RealOption & option( reals_[ key ] );
		option.check_restricted_access(true);
		return option;
	}


	/// @brief OptionCollection[ StringOptionKey ]
	inline
	StringOption const &
	operator []( StringOptionKey const & key ) const
	{
		StringOption const & option( strings_[ key ] );
		option.check_restricted_access(true);
		return option;
	}


	/// @brief OptionCollection[ StringOptionKey ]
	inline
	StringOption &
	operator []( StringOptionKey const & key )
	{
		StringOption & option( strings_[ key ] );
		option.check_restricted_access(true);
		return option;
	}


	/// @brief OptionCollection[ FileOptionKey ]
	inline
	FileOption const &
	operator []( FileOptionKey const & key ) const
	{
		FileOption const & option( files_[ key ] );
		option.check_restricted_access(true);
		return option;
	}


	/// @brief OptionCollection[ FileOptionKey ]
	inline
	FileOption &
	operator []( FileOptionKey const & key )
	{
		FileOption & option( files_[ key ] );
		option.check_restricted_access(true);
		return option;
	}


	/// @brief OptionCollection[ PathOptionKey ]
	inline
	PathOption const &
	operator []( PathOptionKey const & key ) const
	{
		PathOption const & option( paths_[ key ] );
		option.check_restricted_access(true);
		return option;
	}


	/// @brief OptionCollection[ PathOptionKey ]
	inline
	PathOption &
	operator []( PathOptionKey const & key )
	{
		PathOption & option( paths_[ key ] );
		option.check_restricted_access(true);
		return option;
	}


	/// @brief OptionCollection[ AnyOptionKey ]
	inline
	Option const &
	operator []( AnyOptionKey const & key ) const
	{
		Option const & option( anys_[ key ] );
		option.check_restricted_access(true);
		return option;
	}


	/// @brief OptionCollection[ AnyOptionKey ]
	inline
	Option &
	operator []( AnyOptionKey const & key )
	{
		Option & option( anys_[ key ] );
		option.check_restricted_access(true);
		return option;
	}


	/// @brief OptionCollection[ AnyOptionKey ] with option type template argument
	template< typename OptionType >
	inline
	OptionType const &
	operator []( AnyOptionKey const & key ) const
	{
		OptionType const & option(
			utility::down_cast< OptionType const & >( anys_[ key ] ));
		option.check_restricted_access(true);
		return option;
	}


	/// @brief OptionCollection[ AnyOptionKey ] with option type template argument
	template< typename OptionType >
	inline
	OptionType &
	operator []( AnyOptionKey const & key )
	{
		OptionType & option(
			utility::down_cast< OptionType const & >( anys_[ key ] ));
		option.check_restricted_access(true);
		return option;
	}


	/// @brief OptionCollection[ BooleanVectorOptionKey ]
	inline
	BooleanVectorOption const &
	operator []( BooleanVectorOptionKey const & key ) const
	{
		BooleanVectorOption const & option( boolean_vectors_[ key ] );
		option.check_restricted_access(true);
		return option;
	}


	/// @brief OptionCollection[ BooleanVectorOptionKey ]
	inline
	BooleanVectorOption &
	operator []( BooleanVectorOptionKey const & key )
	{
		BooleanVectorOption & option( boolean_vectors_[ key ] );
		option.check_restricted_access(true);
		return option;
	}


	/// @brief OptionCollection[ IntegerVectorOptionKey ]
	inline
	IntegerVectorOption const &
	operator []( IntegerVectorOptionKey const & key ) const
	{
		IntegerVectorOption const & option( integer_vectors_[ key ] );
		option.check_restricted_access(true);
		return option;
	}


	/// @brief OptionCollection[ IntegerVectorOptionKey ]
	inline
	IntegerVectorOption &
	operator []( IntegerVectorOptionKey const & key )
	{
		IntegerVectorOption & option( integer_vectors_[ key ] );
		option.check_restricted_access(true);
		return option;
	}


	/// @brief OptionCollection[ RealVectorOptionKey ]
	inline
	RealVectorOption const &
	operator []( RealVectorOptionKey const & key ) const
	{
		RealVectorOption const & option( real_vectors_[ key ] );
		option.check_restricted_access(true);
		return option;
	}


	/// @brief OptionCollection[ RealVectorOptionKey ]
	inline
	RealVectorOption &
	operator []( RealVectorOptionKey const & key )
	{
		RealVectorOption & option( real_vectors_[ key ] );
		option.check_restricted_access(true);
		return option;
	}

	/// @brief OptionCollection[ ResidueChainVectorOptionKey ]
	inline
	ResidueChainVectorOption const &
	operator []( ResidueChainVectorOptionKey const & key ) const
	{
		ResidueChainVectorOption const & option( residue_chain_vectors_[ key ] );
		option.check_restricted_access(true);
		return option;
	}


	/// @brief OptionCollection[ ResidueChainVectorOptionKey ]
	inline
	ResidueChainVectorOption &
	operator []( ResidueChainVectorOptionKey const & key )
	{
		ResidueChainVectorOption & option( residue_chain_vectors_[ key ] );
		option.check_restricted_access(true);
		return option;
	}


	/// @brief OptionCollection[ StringVectorOptionKey ]
	inline
	StringVectorOption const &
	operator []( StringVectorOptionKey const & key ) const
	{
		StringVectorOption const & option( string_vectors_[ key ] );
		option.check_restricted_access(true);
		return option;
	}


	/// @brief OptionCollection[ StringVectorOptionKey ]
	inline
	StringVectorOption &
	operator []( StringVectorOptionKey const & key )
	{
		StringVectorOption & option( string_vectors_[ key ] );
		option.check_restricted_access(true);
		return option;
	}


	/// @brief OptionCollection[ FileVectorOptionKey ]
	inline
	FileVectorOption const &
	operator []( FileVectorOptionKey const & key ) const
	{
		FileVectorOption const & option( file_vectors_[ key ] );
		option.check_restricted_access(true);
		return option;
	}


	/// @brief OptionCollection[ FileVectorOptionKey ]
	inline
	FileVectorOption &
	operator []( FileVectorOptionKey const & key )
	{
		FileVectorOption & option( file_vectors_[ key ] );
		option.check_restricted_access(true);
		return option;
	}


	/// @brief OptionCollection[ PathVectorOptionKey ]
	inline
	PathVectorOption const &
	operator []( PathVectorOptionKey const & key ) const
	{
		PathVectorOption const & option( path_vectors_[ key ] );
		option.check_restricted_access(true);
		return option;
	}


	/// @brief OptionCollection[ PathVectorOptionKey ]
	inline
	PathVectorOption &
	operator []( PathVectorOptionKey const & key )
	{
		PathVectorOption & option( path_vectors_[ key ] );
		option.check_restricted_access(true);
		return option;
	}


	/// @brief OptionCollection[ AnyVectorOptionKey ]
	inline
	VectorOption const &
	operator []( AnyVectorOptionKey const & key ) const
	{
		VectorOption const & option( any_vectors_[ key ] );
		option.check_restricted_access(true);
		return option;
	}


	/// @brief OptionCollection[ AnyVectorOptionKey ]
	inline
	VectorOption &
	operator []( AnyVectorOptionKey const & key )
	{
		VectorOption & option( any_vectors_[ key ] );
		option.check_restricted_access(true);
		return option;
	}


	/// @brief OptionCollection[ AnyVectorOptionKey ] with option type template argument
	template< typename VectorOptionType >
	inline
	VectorOptionType const &
	operator []( AnyVectorOptionKey const & key ) const
	{
		VectorOptionType const & option(
			utility::down_cast< VectorOptionType const & >( any_vectors_[ key ] ));
		option.check_restricted_access(true);
		return option;
	}


	/// @brief OptionCollection[ AnyVectorOptionKey ] with option type template argument
	template< typename VectorOptionType >
	inline
	VectorOptionType &
	operator []( AnyVectorOptionKey const & key )
	{
		VectorOptionType & option(
			utility::down_cast< VectorOptionType const & >( any_vectors_[ key ] ));
		option.check_restricted_access(true);
		return option;
	}


	/// @brief OptionCollection[ OptionKey ]
	inline
	Option const &
	operator []( OptionKey const & key ) const
	{
		using utility::down_cast;
		runtime_assert( all_.has( key ) );
		Option const * option;
		switch ( all_[ key ] ) {
		case BOOLEAN_OPTION:
			option = &booleans_[ down_cast< BooleanOptionKey const & >( key ) ];
			break;
		case INTEGER_OPTION:
			option = &integers_[ down_cast< IntegerOptionKey const & >( key ) ];
			break;
		case REAL_OPTION:
			option = &reals_[ down_cast< RealOptionKey const & >( key ) ];
			break;
		case STRING_OPTION:
			option = &strings_[ down_cast< StringOptionKey const & >( key ) ];
			break;
		case FILE_OPTION:
			option = &files_[ down_cast< FileOptionKey const & >( key ) ];
			break;
		case PATH_OPTION:
			option = &paths_[ down_cast< PathOptionKey const & >( key ) ];
			break;
		case ANY_OPTION:
			option = anys_[ down_cast< AnyOptionKey const & >( key ) ];
			break;
		case BOOLEAN_VECTOR_OPTION:
			option = &boolean_vectors_[ down_cast< BooleanVectorOptionKey const & >( key ) ];
			break;
		case INTEGER_VECTOR_OPTION:
			option = &integer_vectors_[ down_cast< IntegerVectorOptionKey const & >( key ) ];
			break;
		case REAL_VECTOR_OPTION:
			option = &real_vectors_[ down_cast< RealVectorOptionKey const & >( key ) ];
			break;
		case RESIDUE_CHAIN_VECTOR_OPTION:
			option = &residue_chain_vectors_[ down_cast< ResidueChainVectorOptionKey const & >( key ) ];
			break;
		case STRING_VECTOR_OPTION:
			option = &string_vectors_[ down_cast< StringVectorOptionKey const & >( key ) ];
			break;
		case FILE_VECTOR_OPTION:
			option = &file_vectors_[ down_cast< FileVectorOptionKey const & >( key ) ];
			break;
		case PATH_VECTOR_OPTION:
			option = &path_vectors_[ down_cast< PathVectorOptionKey const & >( key ) ];
			break;
		case ANY_VECTOR_OPTION:
			option = any_vectors_[ down_cast< AnyVectorOptionKey const & >( key ) ];
			break;
		case UNKNOWN_OPTION:
			std::cerr << "Offending key is: " << all_[ key ] << "  " << (int)all_[ key ] << std::endl;
			runtime_assert( false ); // Shouldn't get here
			return *booleans_.begin(); // Keep compiler happy
			break;
		default :
			std::cerr << "Offending key is: " << all_[ key ] << "  " << (int)all_[ key ] << std::endl;
			runtime_assert( false ); // Shouldn't get here
			return *booleans_.begin(); // Keep compiler happy
			break;
		}

		option->check_restricted_access(true);
		return *option;
	}


	/// @brief OptionCollection[ OptionKey ]
	inline
	Option &
	operator []( OptionKey const & key )
	{
		using utility::down_cast;
		runtime_assert( all_.has( key ) );
		Option * option;
		switch ( all_[ key ] ) {
		case BOOLEAN_OPTION:
			option = &booleans_[ down_cast< BooleanOptionKey const & >( key ) ];
			break;
		case INTEGER_OPTION:
			option = &integers_[ down_cast< IntegerOptionKey const & >( key ) ];
			break;
		case REAL_OPTION:
			option = &reals_[ down_cast< RealOptionKey const & >( key ) ];
			break;
		case STRING_OPTION:
			option = &strings_[ down_cast< StringOptionKey const & >( key ) ];
			break;
		case FILE_OPTION:
			option = &files_[ down_cast< FileOptionKey const & >( key ) ];
			break;
		case PATH_OPTION:
			option = &paths_[ down_cast< PathOptionKey const & >( key ) ];
			break;
		case ANY_OPTION:
			option = anys_[ down_cast< AnyOptionKey const & >( key ) ];
			break;
		case BOOLEAN_VECTOR_OPTION:
			option = &boolean_vectors_[ down_cast< BooleanVectorOptionKey const & >( key ) ];
			break;
		case INTEGER_VECTOR_OPTION:
			option = &integer_vectors_[ down_cast< IntegerVectorOptionKey const & >( key ) ];
			break;
		case REAL_VECTOR_OPTION:
			option = &real_vectors_[ down_cast< RealVectorOptionKey const & >( key ) ];
			break;
		case RESIDUE_CHAIN_VECTOR_OPTION:
			option = &residue_chain_vectors_[ down_cast< ResidueChainVectorOptionKey const & >( key ) ];
			break;
		case STRING_VECTOR_OPTION:
			option = &string_vectors_[ down_cast< StringVectorOptionKey const & >( key ) ];
			break;
		case FILE_VECTOR_OPTION:
			option = &file_vectors_[ down_cast< FileVectorOptionKey const & >( key ) ];
			break;
		case PATH_VECTOR_OPTION:
			option = &path_vectors_[ down_cast< PathVectorOptionKey const & >( key ) ];
			break;
		case ANY_VECTOR_OPTION:
			option = any_vectors_[ down_cast< AnyVectorOptionKey const & >( key ) ];
			break;
		case UNKNOWN_OPTION:
			std::cerr << "Offending key is: " << all_[ key ] << "  " << (int)all_[ key ] << std::endl;
			runtime_assert( false ); // Shouldn't get here
			return *booleans_.begin(); // Keep compiler happy
			break;
		default :
			runtime_assert( false ); // Shouldn't get here
			option = &(*booleans_.begin()); // Keep compiler happy
			break;
		}
		option->check_restricted_access(true);
		return *option;
	}


	/// @brief OptionCollection[ OptionKey ] with option type template argument
	template< typename OptionType >
	inline
	OptionType const &
	operator []( OptionKey const & key ) const
	{
		return utility::down_cast< OptionType const & >( operator []( key ) );
	}


	/// @brief OptionCollection[ OptionKey ] with option type template argument
	template< typename OptionType >
	inline
	OptionType &
	operator []( OptionKey const & key )
	{
		return utility::down_cast< OptionType & >( operator []( key ) );
	}


	/// @brief Option by BooleanOptionKey
	inline
	BooleanOption const &
	operator ()( BooleanOptionKey const & key, bool check_restricted_access=true ) const
	{
		BooleanOption const & option(booleans_[ key ]);
		option.check_restricted_access(check_restricted_access);
		return option;
	}


	/// @brief Option by BooleanOptionKey
	inline
	BooleanOption &
	operator ()( BooleanOptionKey const & key, bool check_restricted_access=true )
	{
		BooleanOption & option( booleans_[ key ]);
		option.check_restricted_access(check_restricted_access);
		return option;
	}


	/// @brief Option by IntegerOptionKey
	inline
	IntegerOption const &
	operator ()( IntegerOptionKey const & key, bool check_restricted_access=true ) const
	{
		IntegerOption const & option( integers_[ key ]);
		option.check_restricted_access(check_restricted_access);
		return option;
	}


	/// @brief Option by IntegerOptionKey
	inline
	IntegerOption &
	operator ()( IntegerOptionKey const & key, bool check_restricted_access=true )
	{
		IntegerOption & option( integers_[ key ]);
		option.check_restricted_access(check_restricted_access);
		return option;
	}


	/// @brief Option by RealOptionKey
	inline
	RealOption const &
	operator ()( RealOptionKey const & key, bool check_restricted_access=true ) const
	{
		RealOption const & option( reals_[ key ]);
		option.check_restricted_access(check_restricted_access);
		return option;
	}


	/// @brief Option by RealOptionKey
	inline
	RealOption &
	operator ()( RealOptionKey const & key, bool check_restricted_access=true )
	{
		RealOption & option( reals_[ key ]);
		option.check_restricted_access(check_restricted_access);
		return option;
	}


	/// @brief Option by StringOptionKey
	inline
	StringOption const &
	operator ()( StringOptionKey const & key, bool check_restricted_access=true ) const
	{
		StringOption const & option( strings_[ key ]);
		option.check_restricted_access(check_restricted_access);
		return option;
	}


	/// @brief Option by StringOptionKey
	inline
	StringOption &
	operator ()( StringOptionKey const & key, bool check_restricted_access=true )
	{
		StringOption & option( strings_[ key ]);
		option.check_restricted_access(check_restricted_access);
		return option;
	}


	/// @brief Option by FileOptionKey
	inline
	FileOption const &
	operator ()( FileOptionKey const & key, bool check_restricted_access=true ) const
	{
		FileOption const & option( files_[ key ]);
		option.check_restricted_access(check_restricted_access);
		return option;
	}


	/// @brief Option by FileOptionKey
	inline
	FileOption &
	operator ()( FileOptionKey const & key, bool check_restricted_access=true )
	{
		FileOption & option( files_[ key ]);
		option.check_restricted_access(check_restricted_access);
		return option;
	}


	/// @brief Option by PathOptionKey
	inline
	PathOption const &
	operator ()( PathOptionKey const & key, bool check_restricted_access=true ) const
	{
		PathOption const & option( paths_[ key ]);
		option.check_restricted_access(check_restricted_access);
		return option;
	}


	/// @brief Option by PathOptionKey
	inline
	PathOption &
	operator ()( PathOptionKey const & key, bool check_restricted_access=true )
	{
		PathOption & option( paths_[ key ]);
		option.check_restricted_access(check_restricted_access);
		return option;
	}


	/// @brief Option by AnyOptionKey
	inline
	Option const &
	operator ()( AnyOptionKey const & key, bool check_restricted_access=true ) const
	{
		Option const & option( anys_[ key ]);
		option.check_restricted_access(check_restricted_access);
		return option;
	}


	/// @brief Option by AnyOptionKey
	inline
	Option &
	operator ()( AnyOptionKey const & key, bool check_restricted_access=true )
	{
		Option & option( anys_[ key ]);
		option.check_restricted_access(check_restricted_access);
		return option;
	}


	/// @brief Option by AnyOptionKey with option type template argument
	template< typename OptionType >
	inline
	OptionType const &
	operator ()( AnyOptionKey const & key, bool check_restricted_access=true ) const
	{
		OptionType const & option(
			utility::down_cast< OptionType const & >( anys_[ key ] ));
		option.check_restricted_access(check_restricted_access);
		return option;
	}


	/// @brief Option by AnyOptionKey with option type template argument
	template< typename OptionType >
	inline
	OptionType &
	operator ()( AnyOptionKey const & key, bool check_restricted_access=true )
	{
		OptionType & option(
			utility::down_cast< OptionType const & >( anys_[ key ] ));
		option.check_restricted_access(check_restricted_access);
		return option;
	}


	/// @brief VectorOption by BooleanVectorOptionKey
	inline
	BooleanVectorOption const &
	operator ()( BooleanVectorOptionKey const & key, bool check_restricted_access=true ) const
	{
		BooleanVectorOption const & option( boolean_vectors_[ key ]);
		option.check_restricted_access(check_restricted_access);
		return option;
	}


	/// @brief VectorOption by BooleanVectorOptionKey
	inline
	BooleanVectorOption &
	operator ()( BooleanVectorOptionKey const & key, bool check_restricted_access=true )
	{
		BooleanVectorOption & option( boolean_vectors_[ key ]);
		option.check_restricted_access(check_restricted_access);
		return option;;
	}


	/// @brief VectorOption by IntegerVectorOptionKey
	inline
	IntegerVectorOption const &
	operator ()( IntegerVectorOptionKey const & key, bool check_restricted_access=true ) const
	{
		IntegerVectorOption const & option( integer_vectors_[ key ]);
		option.check_restricted_access(check_restricted_access);
		return option;
	}


	/// @brief VectorOption by IntegerVectorOptionKey
	inline
	IntegerVectorOption &
	operator ()( IntegerVectorOptionKey const & key, bool check_restricted_access=true )
	{
		IntegerVectorOption & option( integer_vectors_[ key ]);
		option.check_restricted_access(check_restricted_access);
		return option;
	}


	/// @brief VectorOption by RealVectorOptionKey
	inline
	RealVectorOption const &
	operator ()( RealVectorOptionKey const & key, bool check_restricted_access=true ) const
	{
		RealVectorOption const & option( real_vectors_[ key ]);
		option.check_restricted_access(check_restricted_access);
		return option;
	}


	/// @brief VectorOption by RealVectorOptionKey
	inline
	RealVectorOption &
	operator ()( RealVectorOptionKey const & key, bool check_restricted_access=true )
	{
		RealVectorOption & option( real_vectors_[ key ]);
		option.check_restricted_access(check_restricted_access);
		return option;
	}


	/// @brief VectorOption by ResidueChainVectorOptionKey
	inline
	ResidueChainVectorOption const &
	operator ()( ResidueChainVectorOptionKey const & key, bool check_restricted_access=true ) const
	{
		ResidueChainVectorOption const & option( residue_chain_vectors_[ key ]);
		option.check_restricted_access(check_restricted_access);
		return option;
	}


	/// @brief VectorOption by ResidueChainVectorOptionKey
	inline
	ResidueChainVectorOption &
	operator ()( ResidueChainVectorOptionKey const & key, bool check_restricted_access=true )
	{
		ResidueChainVectorOption & option( residue_chain_vectors_[ key ]);
		option.check_restricted_access(check_restricted_access);
		return option;
	}


	/// @brief VectorOption by StringVectorOptionKey
	inline
	StringVectorOption const &
	operator ()( StringVectorOptionKey const & key, bool check_restricted_access=true ) const
	{
		StringVectorOption const & option( string_vectors_[ key ]);
		option.check_restricted_access(check_restricted_access);
		return option;
	}


	/// @brief VectorOption by StringVectorOptionKey
	inline
	StringVectorOption &
	operator ()( StringVectorOptionKey const & key, bool check_restricted_access=true )
	{
		StringVectorOption & option( string_vectors_[ key ]);
		option.check_restricted_access(check_restricted_access);
		return option;
	}


	/// @brief VectorOption by FileVectorOptionKey
	inline
	FileVectorOption const &
	operator ()( FileVectorOptionKey const & key, bool check_restricted_access=true ) const
	{
		FileVectorOption const & option( file_vectors_[ key ]);
		option.check_restricted_access(check_restricted_access);
		return option;
	}


	/// @brief VectorOption by FileVectorOptionKey
	inline
	FileVectorOption &
	operator ()( FileVectorOptionKey const & key, bool check_restricted_access=true )
	{
		FileVectorOption & option( file_vectors_[ key ]);
		option.check_restricted_access(check_restricted_access);
		return option;
	}


	/// @brief VectorOption by PathVectorOptionKey
	inline
	PathVectorOption const &
	operator ()( PathVectorOptionKey const & key, bool check_restricted_access=true ) const
	{
		PathVectorOption const & option( path_vectors_[ key ]);
		option.check_restricted_access(check_restricted_access);
		return option;
	}


	/// @brief VectorOption by PathVectorOptionKey
	inline
	PathVectorOption &
	operator ()( PathVectorOptionKey const & key, bool check_restricted_access=true )
	{
		PathVectorOption & option( path_vectors_[ key ]);
		option.check_restricted_access(check_restricted_access);
		return option;
	}


	/// @brief VectorOption by AnyVectorOptionKey
	inline
	VectorOption const &
	operator ()( AnyVectorOptionKey const & key, bool check_restricted_access=true ) const
	{
		VectorOption const & option(any_vectors_[ key ] );
		option.check_restricted_access(check_restricted_access);
		return option;
	}


	/// @brief VectorOption by AnyVectorOptionKey
	inline
	VectorOption &
	operator ()( AnyVectorOptionKey const & key, bool check_restricted_access=true )
	{
		VectorOption & option( any_vectors_[ key ] );
		option.check_restricted_access(check_restricted_access);
		return option;
	}


	/// @brief VectorOption by AnyVectorOptionKey with option type template argument
	template< typename VectorOptionType >
	inline
	VectorOptionType const &
	operator ()( AnyVectorOptionKey const & key, bool check_restricted_access=true ) const
	{
		VectorOptionType const & option(
			utility::down_cast< VectorOptionType const & >( any_vectors_[ key ] ));
		option.check_restricted_access(check_restricted_access);
		return option;
	}


	/// @brief VectorOption by AnyVectorOptionKey with option type template argument
	template< typename VectorOptionType >
	inline
	VectorOptionType &
	operator ()( AnyVectorOptionKey const & key, bool check_restricted_access=true )
	{
		VectorOptionType & option(
			utility::down_cast< VectorOptionType & >( any_vectors_[ key ] ));
		option.check_restricted_access(check_restricted_access);
		return option;
	}


	/// @brief Option by OptionKey
	inline
	Option const &
	operator ()( OptionKey const & key, bool check_restricted_access=true ) const
	{
		using utility::down_cast;
		runtime_assert( all_.has( key ) );
		Option const * option;
		switch ( all_[ key ] ) {
		case BOOLEAN_OPTION:
			option = &booleans_[ down_cast< BooleanOptionKey const & >( key ) ];
			break;
		case INTEGER_OPTION:
			option = &integers_[ down_cast< IntegerOptionKey const & >( key ) ];
			break;
		case REAL_OPTION:
			option = &reals_[ down_cast< RealOptionKey const & >( key ) ];
			break;
		case STRING_OPTION:
			option = &strings_[ down_cast< StringOptionKey const & >( key ) ];
			break;
		case FILE_OPTION:
			option = &files_[ down_cast< FileOptionKey const & >( key ) ];
			break;
		case PATH_OPTION:
			option = &paths_[ down_cast< PathOptionKey const & >( key ) ];
			break;
		case ANY_OPTION:
			option = anys_[ down_cast< AnyOptionKey const & >( key ) ];
			break;
		case BOOLEAN_VECTOR_OPTION:
			option = &boolean_vectors_[ down_cast< BooleanVectorOptionKey const & >( key ) ];
			break;
		case INTEGER_VECTOR_OPTION:
			option = &integer_vectors_[ down_cast< IntegerVectorOptionKey const & >( key ) ];
			break;
		case REAL_VECTOR_OPTION:
			option = &real_vectors_[ down_cast< RealVectorOptionKey const & >( key ) ];
			break;
		case RESIDUE_CHAIN_VECTOR_OPTION:
			option = &residue_chain_vectors_[ down_cast< ResidueChainVectorOptionKey const & >( key ) ];
			break;
		case STRING_VECTOR_OPTION:
			option = &string_vectors_[ down_cast< StringVectorOptionKey const & >( key ) ];
			break;
		case FILE_VECTOR_OPTION:
			option = &file_vectors_[ down_cast< FileVectorOptionKey const & >( key ) ];
			break;
		case PATH_VECTOR_OPTION:
			option = &path_vectors_[ down_cast< PathVectorOptionKey const & >( key ) ];
			break;
		case ANY_VECTOR_OPTION:
			option = any_vectors_[ down_cast< AnyVectorOptionKey const & >( key ) ];
			break;
		case UNKNOWN_OPTION:
			std::cerr << "Offending key is: " << all_[ key ] << "  " << (int)all_[ key ] << std::endl;
			runtime_assert( false ); // Shouldn't get here
			return *booleans_.begin(); // Keep compiler happy
			break;
		default :
			std::cerr << "Offending key is: " << all_[ key ] << "  " << (int)all_[ key ] << std::endl;
			runtime_assert( false ); // Shouldn't get here
			return *booleans_.begin(); // Keep compiler happy
			break;
		}

		option->check_restricted_access(check_restricted_access);
		return *option;
	}


	/// @brief Option by OptionKey
	inline
	Option &
	operator ()( OptionKey const & key, bool check_restricted_access=true )
	{
		using utility::down_cast;
		runtime_assert( all_.has( key ) );
		Option * option;
		switch ( all_[ key ] ) {
		case BOOLEAN_OPTION:
			option = &booleans_[ down_cast< BooleanOptionKey const & >( key ) ];
			break;
		case INTEGER_OPTION:
			option = &integers_[ down_cast< IntegerOptionKey const & >( key ) ];
			break;
		case REAL_OPTION:
			option = &reals_[ down_cast< RealOptionKey const & >( key ) ];
			break;
		case STRING_OPTION:
			option = &strings_[ down_cast< StringOptionKey const & >( key ) ];
			break;
		case FILE_OPTION:
			option = &files_[ down_cast< FileOptionKey const & >( key ) ];
			break;
		case PATH_OPTION:
			option = &paths_[ down_cast< PathOptionKey const & >( key ) ];
			break;
		case ANY_OPTION:
			option = anys_[ down_cast< AnyOptionKey const & >( key ) ];
			break;
		case BOOLEAN_VECTOR_OPTION:
			option = &boolean_vectors_[ down_cast< BooleanVectorOptionKey const & >( key ) ];
			break;
		case INTEGER_VECTOR_OPTION:
			option = &integer_vectors_[ down_cast< IntegerVectorOptionKey const & >( key ) ];
			break;
		case REAL_VECTOR_OPTION:
			option = &real_vectors_[ down_cast< RealVectorOptionKey const & >( key ) ];
			break;
		case RESIDUE_CHAIN_VECTOR_OPTION:
			option = &residue_chain_vectors_[ down_cast< ResidueChainVectorOptionKey const & >( key ) ];
			break;
		case STRING_VECTOR_OPTION:
			option = &string_vectors_[ down_cast< StringVectorOptionKey const & >( key ) ];
			break;
		case FILE_VECTOR_OPTION:
			option = &file_vectors_[ down_cast< FileVectorOptionKey const & >( key ) ];
			break;
		case PATH_VECTOR_OPTION:
			option = &path_vectors_[ down_cast< PathVectorOptionKey const & >( key ) ];
			break;
		case ANY_VECTOR_OPTION:
			option = any_vectors_[ down_cast< AnyVectorOptionKey const & >( key ) ];
			break;
		case UNKNOWN_OPTION:
			std::cerr << "Offending key is: " << all_[ key ] << "  " << (int)all_[ key ] << std::endl;
			runtime_assert( false ); // Shouldn't get here
			return *booleans_.begin(); // Keep compiler happy
			break;
		default :
			runtime_assert( false ); // Shouldn't get here
			option = &(*booleans_.begin()); // Keep compiler happy
			break;
		}
		option->check_restricted_access(check_restricted_access);
		return *option;
	}


	/// @brief Option by OptionKey with option type template argument
	template< typename OptionType >
	inline
	OptionType const &
	operator ()( OptionKey const & key, bool check_restricted_access=true ) const
	{
		return utility::down_cast< OptionType const & >(
			operator ()< OptionType >( key, check_restricted_access ));
	}


	/// @brief Option by OptionKey with option type template argument
	template< typename OptionType >
	inline
	OptionType &
	operator ()( OptionKey const & key, bool check_restricted_access=true )
	{
		return utility::down_cast< OptionType & >(
			operator ()< OptionType >( key, check_restricted_access ));
	}


public: // I/O


	/// @brief Output to stream
	friend
	std::ostream &
	operator <<( std::ostream & stream, OptionCollection const & options );


private: // Methods


	/// @brief Load a user-specified option argument from a command line
	void
	load_option_cl(
		std::string const & arg_string, // Lead argument string
		ValueStrings & arg_strings, // Argument strings: Value string(s) in front
		std::string const & pid // Previous option id
	);

	/// @brief Load one option from user specified file
	void load_option_from_file(
		std::string const & arg_string, // Lead argument string
		ValueStrings & arg_strings, // Argument strings: Value string(s) in front
		std::string const & pid // Previous option id
	);

	/// @brief Load a user-specified option argument from an @file
	void
	load_option_file(
		std::string const & arg_string, // Argument string
		ValueStrings & val_strings, // Value strings
		std::string & cid, // Context option id
		bool const cl_context = false // Use command line key context?
	);


	/// @brief Set a user-specified option value from a command line
	void
	set_option_value_cl(
		std::string const & key_id, // Option key id
		ValueStrings & arg_strings // Argument strings: Value string(s) in front
	);


	/// @brief Set a user-specified option value from an @file
	void
	set_option_value_file(
		std::string const & key_id, // Option key id
		ValueStrings & val_strings // Value strings
	);


public: // Static functions

	/// @brief add OptionKey to list of application relevant options
	static
	void add_relevant(const OptionKey & key)
	{
		relevant_.push_back( &key );
	}

	/// @brief Checks if option has been registered as relevant.
	static
	bool
	is_relevant( OptionKey const & key );


	/// @brief Space-prefixed string except blank if string is empty
	inline
	static
	std::string
	space_prefixed(
		std::string const & s, // String
		int const n = 1 // Number of spaces
	)
	{
		return ( s.empty() ? "" : std::string( n, ' ' ) + s );
	}


	/// @brief Tab-prefixed string except blank if string is empty
	inline
	static
	std::string
	tab_prefixed(
		std::string const & s, // String
		int const n = 1 // Number of tabs
	)
	{
		return ( s.empty() ? "" : std::string( n, '\t' ) + s );
	}


	/// @brief Check that a key's identifiers are legal
	static
	void
	check_key( OptionKey const & key );


	/// @brief Check that an option's identifiers are legal
	static
	void
	check_key( Option const & option );


	/// @brief Find a user-specified option key in a command line context
	/// @note  Searches up the context to find a match
	static
	std::string
	find_key_cl(
		std::string const & key_string, // Option key string entered
		std::string const & cid, // Context option id
		bool const top // Top-level context?
	);


	/// @brief Find a user-specified option key in an indented @file context
	/// @note  Looks in the context to find a match
	static
	std::string
	find_key_file(
		std::string const & key_string, // Option key string entered
		std::string const & cid, // Context option id
		bool const top // Top-level context?
	);


	/// @brief Number of parts in an option id
	static
	std::string::size_type
	n_part( std::string const & s );


	/// @brief Number of prefix parts of two ids that match
	static
	std::string::size_type
	n_part_prefix_match(
		std::string const & s,
		std::string const & t
	);


	/// @brief Prefix of an option id with a specified number of parts
	static
	std::string
	prefix(
		std::string const & s, // String
		int const n = 1 // Number of prefix parts desired
	);


	/// @brief Suffix of an option id with a specified number of parts
	static
	std::string
	suffix(
		std::string const & s, // String
		int const n = 1 // Number of suffix parts desired
	);


	/// @brief Trim a specified number of parts from the suffix of an option id
	static
	std::string &
	trim(
		std::string & s, // String
		int const n = 1 // Number of suffix parts to trim
	);


	/// @brief Prefix of an option id with a specified number of suffix parts removed
	static
	std::string
	trimmed(
		std::string const & s, // String
		int const n = 1 // Number of suffix parts to trim
	);


	/// @brief Cleaned option id with repeat colons condensed
	static
	std::string
	cleaned( std::string const & s );


	/// @brief Merged option ids with the minimal suffix-prefix overlap, if any, removed
	static
	std::string
	merged(
		std::string const & s, // Lead id
		std::string const & t // Tail id
	);


	/// @brief String wrapped and indented
	static
	std::string
	wrapped(
		std::string const & s, // String to wrap
		std::string::size_type const indent = 10, // Width to indent continuation lines
		std::string::size_type const width = 80, // Column width to wrap at [80]
		std::string const header_for_extra_lines = ""
	);


	/// @brief modify 'show_accessed_options' flag;
	static void set_show_accessed_options_flag(bool v) { show_accessed_options_ = v; }

	/// @brief modify 'show_unused_options' flag;
	static void set_show_unused_options_flag(bool v) { show_unused_options_ = v; }

private: // Fields


	/// @brief Boolean options
	Booleans booleans_;

	/// @brief Integer options
	Integers integers_;

	/// @brief Real options
	Reals reals_;

	/// @brief String options
	Strings strings_;

	/// @brief File options
	Files files_;

	/// @brief Path options
	Paths paths_;

	/// @brief Any value type options
	Anys anys_;

	/// @brief Boolean vector options
	BooleanVectors boolean_vectors_;

	/// @brief Integer vector options
	IntegerVectors integer_vectors_;

	/// @brief Real vector options
	RealVectors real_vectors_;

	/// @brief Residue/chain vector options
	ResidueChainVectors residue_chain_vectors_;

	/// @brief String vector options
	StringVectors string_vectors_;

	/// @brief File vector options
	FileVectors file_vectors_;

	/// @brief Path vector options
	PathVectors path_vectors_;

	/// @brief Any value type vector options
	AnyVectors any_vectors_;

	/// @brief All options (non-owning pointers)
	All all_;

	/// @brief global list of application relevant options
	static std::vector< OptionKey const *> relevant_;

	/// @brief keep a copy of argv around in case people want to
	/// get at it elsewhere in the code.
	std::string argv_copy_;


	/// @brief Flag indicating that list of accessed option should be printed when destructor of
	///             OptionCollection is called.
	///        This flag is false by default.
	static bool show_accessed_options_;

	/// @brief Flag indicating that list of user specified but inaccessed options should be printed
	///             when destructor of OptionCollection is called.
	///        This flag is false by default.
	static bool show_unused_options_;

}; // OptionCollection


} // namespace options
} // namespace utility


#endif // INCLUDED_utility_options_OptionCollection_HH
