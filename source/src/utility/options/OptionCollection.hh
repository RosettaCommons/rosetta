// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/options/OptionCollection.hh
/// @brief  Program options collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author Modified by Vikram K. Mulligan (vmulligan@flatironinstitute.org) -- added thread safety.


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
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <list>
#include <string>
#include <vector>

// Multithreading headers
#ifdef MULTI_THREADED
#include <utility/thread/ReadWriteMutex.hh>
#endif

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace utility {
namespace options {


/// @brief Program options collection
class OptionCollection : public utility::pointer::ReferenceCount
{


private: // Types


	typedef  std::list< std::string >  ValueStrings;

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

	/// @brief Copy constructor.
	/// @details Needed because the mutex can't be copied.
	OptionCollection( OptionCollection const & src );


	/// @brief Destructor
	~OptionCollection() override;

	/// @brief Copy assignment operator.
	/// @details Needed because the mutext can't be copied.
	OptionCollection & operator=( OptionCollection const & src );

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
#ifdef MULTI_THREADED
		utility::thread::WriteLockGuard writelock( mutex_ );
#endif
		all_( key ) = BOOLEAN_OPTION;
		return ( booleans_( key ) = BooleanOption( key, description ) );
	}


	/// @brief Add an IntegerOption
	inline
	IntegerOption &
	add( IntegerOptionKey const & key, std::string const & description )
	{
		check_key( key );
#ifdef MULTI_THREADED
		utility::thread::WriteLockGuard writelock( mutex_ );
#endif
		all_( key ) = INTEGER_OPTION;
		return ( integers_( key ) = IntegerOption( key, description ) );
	}


	/// @brief Add a RealOption
	inline
	RealOption &
	add( RealOptionKey const & key, std::string const & description )
	{
		check_key( key );
#ifdef MULTI_THREADED
		utility::thread::WriteLockGuard writelock( mutex_ );
#endif
		all_( key ) = REAL_OPTION;
		return ( reals_( key ) = RealOption( key, description ) );
	}


	/// @brief Add a StringOption
	inline
	StringOption &
	add( StringOptionKey const & key, std::string const & description )
	{
		check_key( key );
#ifdef MULTI_THREADED
		utility::thread::WriteLockGuard writelock( mutex_ );
#endif
		all_( key ) = STRING_OPTION;
		return ( strings_( key ) = StringOption( key, description ) );
	}


	/// @brief Add a FileOption
	inline
	FileOption &
	add( FileOptionKey const & key, std::string const & description )
	{
		check_key( key );
#ifdef MULTI_THREADED
		utility::thread::WriteLockGuard writelock( mutex_ );
#endif
		all_( key ) = FILE_OPTION;
		return ( files_( key ) = FileOption( key, description ) );
	}


	/// @brief Add a PathOption
	inline
	PathOption &
	add( PathOptionKey const & key, std::string const & description )
	{
		check_key( key );
#ifdef MULTI_THREADED
		utility::thread::WriteLockGuard writelock( mutex_ );
#endif
		all_( key ) = PATH_OPTION;
		return ( paths_( key ) = PathOption( key, description ) );
	}


	/// @brief Add a BooleanVectorOption
	inline
	BooleanVectorOption &
	add( BooleanVectorOptionKey const & key, std::string const & description )
	{
		check_key( key );
#ifdef MULTI_THREADED
		utility::thread::WriteLockGuard writelock( mutex_ );
#endif
		all_( key ) = BOOLEAN_VECTOR_OPTION;
		return ( boolean_vectors_( key ) = BooleanVectorOption( key, description ) );
	}


	/// @brief Add an IntegerVectorOption
	inline
	IntegerVectorOption &
	add( IntegerVectorOptionKey const & key, std::string const & description )
	{
		check_key( key );
#ifdef MULTI_THREADED
		utility::thread::WriteLockGuard writelock( mutex_ );
#endif
		all_( key ) = INTEGER_VECTOR_OPTION;
		return ( integer_vectors_( key ) = IntegerVectorOption( key, description ) );
	}


	/// @brief Add a RealVectorOption
	inline
	RealVectorOption &
	add( RealVectorOptionKey const & key, std::string const & description )
	{
		check_key( key );
#ifdef MULTI_THREADED
		utility::thread::WriteLockGuard writelock( mutex_ );
#endif
		all_( key ) = REAL_VECTOR_OPTION;
		return ( real_vectors_( key ) = RealVectorOption( key, description ) );
	}

	/// @brief Add an ResidueChainVectorOption
	inline
	ResidueChainVectorOption &
	add( ResidueChainVectorOptionKey const & key, std::string const & description )
	{
		check_key( key );
#ifdef MULTI_THREADED
		utility::thread::WriteLockGuard writelock( mutex_ );
#endif
		all_( key ) = RESIDUE_CHAIN_VECTOR_OPTION;
		return ( residue_chain_vectors_( key ) = ResidueChainVectorOption( key, description ) );
	}


	/// @brief Add a StringVectorOption
	inline
	StringVectorOption &
	add( StringVectorOptionKey const & key, std::string const & description )
	{
		check_key( key );
#ifdef MULTI_THREADED
		utility::thread::WriteLockGuard writelock( mutex_ );
#endif
		all_( key ) = STRING_VECTOR_OPTION;
		return ( string_vectors_( key ) = StringVectorOption( key, description ) );
	}


	/// @brief Add a FileVectorOption
	inline
	FileVectorOption &
	add( FileVectorOptionKey const & key, std::string const & description )
	{
		check_key( key );
#ifdef MULTI_THREADED
		utility::thread::WriteLockGuard writelock( mutex_ );
#endif
		all_( key ) = FILE_VECTOR_OPTION;
		return ( file_vectors_( key ) = FileVectorOption( key, description ) );
	}


	/// @brief Add a PathVectorOption
	inline
	PathVectorOption &
	add( PathVectorOptionKey const & key, std::string const & description )
	{
		check_key( key );
#ifdef MULTI_THREADED
		utility::thread::WriteLockGuard writelock( mutex_ );
#endif
		all_( key ) = PATH_VECTOR_OPTION;
		return ( path_vectors_( key ) = PathVectorOption( key, description ) );
	}


	/// @brief Add a BooleanOption
	inline
	BooleanOption &
	add( BooleanOption const & opt )
	{
		check_key( opt );
#ifdef MULTI_THREADED
		utility::thread::WriteLockGuard writelock( mutex_ );
#endif
		all_( opt.key() ) = BOOLEAN_OPTION;
		return ( booleans_( opt.key() ) = opt );
	}


	/// @brief Add an IntegerOption
	inline
	IntegerOption &
	add( IntegerOption const & opt )
	{
		check_key( opt );
#ifdef MULTI_THREADED
		utility::thread::WriteLockGuard writelock( mutex_ );
#endif
		all_( opt.key() ) = INTEGER_OPTION;
		return ( integers_( opt.key() ) = opt );
	}


	/// @brief Add a RealOption
	inline
	RealOption &
	add( RealOption const & opt )
	{
		check_key( opt );
#ifdef MULTI_THREADED
		utility::thread::WriteLockGuard writelock( mutex_ );
#endif
		all_( opt.key() ) = REAL_OPTION;
		return ( reals_( opt.key() ) = opt );
	}


	/// @brief Add a StringOption
	inline
	StringOption &
	add( StringOption const & opt )
	{
		check_key( opt );
#ifdef MULTI_THREADED
		utility::thread::WriteLockGuard writelock( mutex_ );
#endif
		all_( opt.key() ) = STRING_OPTION;
		return ( strings_( opt.key() ) = opt );
	}


	/// @brief Add a FileOption
	inline
	FileOption &
	add( FileOption const & opt )
	{
		check_key( opt );
#ifdef MULTI_THREADED
		utility::thread::WriteLockGuard writelock( mutex_ );
#endif
		all_( opt.key() ) = FILE_OPTION;
		return ( files_( opt.key() ) = opt );
	}


	/// @brief Add a PathOption
	inline
	PathOption &
	add( PathOption const & opt )
	{
		check_key( opt );
#ifdef MULTI_THREADED
		utility::thread::WriteLockGuard writelock( mutex_ );
#endif
		all_( opt.key() ) = PATH_OPTION;
		return ( paths_( opt.key() ) = opt );
	}


	/// @brief Add an AnyOption
	template< typename T >
	inline
	AnyOption< T > &
	add( AnyOption< T > const & opt )
	{
		check_key( opt );
#ifdef MULTI_THREADED
		utility::thread::WriteLockGuard writelock( mutex_ );
#endif
		all_( opt.key() ) = ANY_OPTION;
		return ( anys_( opt.key() ) = opt );
	}


	/// @brief Add a BooleanVectorOption
	inline
	BooleanVectorOption &
	add( BooleanVectorOption const & opt )
	{
		check_key( opt );
#ifdef MULTI_THREADED
		utility::thread::WriteLockGuard writelock( mutex_ );
#endif
		all_( opt.key() ) = BOOLEAN_VECTOR_OPTION;
		return ( boolean_vectors_( opt.key() ) = opt );
	}


	/// @brief Add an IntegerVectorOption
	inline
	IntegerVectorOption &
	add( IntegerVectorOption const & opt )
	{
		check_key( opt );
#ifdef MULTI_THREADED
		utility::thread::WriteLockGuard writelock( mutex_ );
#endif
		all_( opt.key() ) = INTEGER_VECTOR_OPTION;
		return ( integer_vectors_( opt.key() ) = opt );
	}


	/// @brief Add a RealVectorOption
	inline
	RealVectorOption &
	add( RealVectorOption const & opt )
	{
		check_key( opt );
#ifdef MULTI_THREADED
		utility::thread::WriteLockGuard writelock( mutex_ );
#endif
		all_( opt.key() ) = REAL_VECTOR_OPTION;
		return ( real_vectors_( opt.key() ) = opt );
	}

	/// @brief Add an ResidueChainVectorOption
	inline
	ResidueChainVectorOption &
	add( ResidueChainVectorOption const & opt )
	{
		check_key( opt );
#ifdef MULTI_THREADED
		utility::thread::WriteLockGuard writelock( mutex_ );
#endif
		all_( opt.key() ) = RESIDUE_CHAIN_VECTOR_OPTION;
		return ( residue_chain_vectors_( opt.key() ) = opt );
	}


	/// @brief Add a StringVectorOption
	inline
	StringVectorOption &
	add( StringVectorOption const & opt )
	{
		check_key( opt );
#ifdef MULTI_THREADED
		utility::thread::WriteLockGuard writelock( mutex_ );
#endif
		all_( opt.key() ) = STRING_VECTOR_OPTION;
		return ( string_vectors_( opt.key() ) = opt );
	}


	/// @brief Add a FileVectorOption
	inline
	FileVectorOption &
	add( FileVectorOption const & opt )
	{
		check_key( opt );
#ifdef MULTI_THREADED
		utility::thread::WriteLockGuard writelock( mutex_ );
#endif
		all_( opt.key() ) = FILE_VECTOR_OPTION;
		return ( file_vectors_( opt.key() ) = opt );
	}


	/// @brief Add a PathVectorOption
	inline
	PathVectorOption &
	add( PathVectorOption const & opt )
	{
		check_key( opt );
#ifdef MULTI_THREADED
		utility::thread::WriteLockGuard writelock( mutex_ );
#endif
		all_( opt.key() ) = PATH_VECTOR_OPTION;
		return ( path_vectors_( opt.key() ) = opt );
	}


	/// @brief Add an AnyVectorOption
	template< typename T >
	inline
	AnyVectorOption< T > &
	add( AnyVectorOption< T > const & opt )
	{
		check_key( opt );
#ifdef MULTI_THREADED
		utility::thread::WriteLockGuard writelock( mutex_ );
#endif
		all_( opt.key() ) = ANY_VECTOR_OPTION;
		return ( any_vectors_( opt.key() ) = opt );
	}


	// /// @brief Specify mutually exclusive options
	// void
	// exclusive(
	//  OptionKey const & key1,
	//  OptionKey const & key2
	// );


	// /// @brief Specify an option that requires a second option to also be specified
	// void
	// requires(
	//  OptionKey const & key1,
	//  OptionKey const & key2
	// );

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

	/// @brief Load all options in a flags file
	void load_options_from_file(
		std::string const & file_string,
		std::string const & cid=""
	);

	/// @brief same as load_options_from_file, but throws exception instead of call to std::exit
	void load_options_from_file_exception(
		std::string const & file_string,
		std::string const & cid=""
	);

	/// @brief Load all options in a flags file
	void load_options_from_stream(
		std::istream& stream,
		std::string const & file_string="STREAM", //for error msg, set if stream is connected to file_string
		std::string const & cid="",
		bool print_lines = false
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
	void show_option_help_hier(OptionKey const &, std::string &group, std::ostream & stream );


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
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		return booleans_.has( key );
	}


	/// @brief Is there an option with an IntegerOptionKey?
	inline
	bool
	has( IntegerOptionKey const & key ) const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		return integers_.has( key );
	}


	/// @brief Is there an option with a RealOptionKey?
	inline
	bool
	has( RealOptionKey const & key ) const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		return reals_.has( key );
	}


	/// @brief Is there an option with a StringOptionKey?
	inline
	bool
	has( StringOptionKey const & key ) const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		return strings_.has( key );
	}


	/// @brief Is there an option with a FileOptionKey?
	inline
	bool
	has( FileOptionKey const & key ) const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		return files_.has( key );
	}


	/// @brief Is there an option with a PathOptionKey?
	inline
	bool
	has( PathOptionKey const & key ) const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		return paths_.has( key );
	}


	/// @brief Is there an option with an AnyOptionKey?
	inline
	bool
	has( AnyOptionKey const & key ) const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		return anys_.has( key );
	}


	/// @brief Is there an option with a BooleanVectorOptionKey?
	inline
	bool
	has( BooleanVectorOptionKey const & key ) const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		return boolean_vectors_.has( key );
	}


	/// @brief Is there an option with an IntegerVectorOptionKey?
	inline
	bool
	has( IntegerVectorOptionKey const & key ) const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		return integer_vectors_.has( key );
	}


	/// @brief Is there an option with a RealVectorOptionKey?
	inline
	bool
	has( RealVectorOptionKey const & key ) const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		return real_vectors_.has( key );
	}

	/// @brief Is there an option with an IntegerVectorOptionKey?
	inline
	bool
	has( ResidueChainVectorOptionKey const & key ) const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		return residue_chain_vectors_.has( key );
	}


	/// @brief Is there an option with a StringVectorOptionKey?
	inline
	bool
	has( StringVectorOptionKey const & key ) const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		return string_vectors_.has( key );
	}


	/// @brief Is there an option with a FileVectorOptionKey?
	inline
	bool
	has( FileVectorOptionKey const & key ) const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		return file_vectors_.has( key );
	}


	/// @brief Is there an option with a PathVectorOptionKey?
	inline
	bool
	has( PathVectorOptionKey const & key ) const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		return path_vectors_.has( key );
	}


	/// @brief Is there an option with an AnyVectorOptionKey?
	inline
	bool
	has( AnyVectorOptionKey const & key ) const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		return any_vectors_.has( key );
	}


	/// @brief Is there an option with an OptionKey?
	inline
	bool
	has( OptionKey const & key ) const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
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
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		BooleanOption const & opt( booleans_[ key ] );
		opt.check_restricted_access(true);
		return opt;
	}


	/// @brief OptionCollection[ BooleanOptionKey ]
	inline
	BooleanOption &
	operator []( BooleanOptionKey const & key )
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		BooleanOption & opt( booleans_[ key ] );
		opt.check_restricted_access(true);
		return opt;
	}


	/// @brief OptionCollection[ IntegerOptionKey ]
	inline
	IntegerOption const &
	operator []( IntegerOptionKey const & key ) const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		IntegerOption const & opt( integers_[ key ] );
		opt.check_restricted_access(true);
		return opt;
	}


	/// @brief OptionCollection[ IntegerOptionKey ]
	inline
	IntegerOption &
	operator []( IntegerOptionKey const & key )
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		IntegerOption & opt( integers_[ key ] );
		opt.check_restricted_access(true);
		return opt;
	}


	/// @brief OptionCollection[ RealOptionKey ]
	inline
	RealOption const &
	operator []( RealOptionKey const & key ) const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		RealOption const & opt( reals_[ key ] );
		opt.check_restricted_access(true);
		return opt;
	}


	/// @brief OptionCollection[ RealOptionKey ]
	inline
	RealOption &
	operator []( RealOptionKey const & key )
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		RealOption & opt( reals_[ key ] );
		opt.check_restricted_access(true);
		return opt;
	}


	/// @brief OptionCollection[ StringOptionKey ]
	inline
	StringOption const &
	operator []( StringOptionKey const & key ) const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		StringOption const & opt( strings_[ key ] );
		opt.check_restricted_access(true);
		return opt;
	}


	/// @brief OptionCollection[ StringOptionKey ]
	inline
	StringOption &
	operator []( StringOptionKey const & key )
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		StringOption & opt( strings_[ key ] );
		opt.check_restricted_access(true);
		return opt;
	}


	/// @brief OptionCollection[ FileOptionKey ]
	inline
	FileOption const &
	operator []( FileOptionKey const & key ) const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		FileOption const & opt( files_[ key ] );
		opt.check_restricted_access(true);
		return opt;
	}


	/// @brief OptionCollection[ FileOptionKey ]
	inline
	FileOption &
	operator []( FileOptionKey const & key )
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		FileOption & opt( files_[ key ] );
		opt.check_restricted_access(true);
		return opt;
	}


	/// @brief OptionCollection[ PathOptionKey ]
	inline
	PathOption const &
	operator []( PathOptionKey const & key ) const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		PathOption const & opt( paths_[ key ] );
		opt.check_restricted_access(true);
		return opt;
	}


	/// @brief OptionCollection[ PathOptionKey ]
	inline
	PathOption &
	operator []( PathOptionKey const & key )
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		PathOption & opt( paths_[ key ] );
		opt.check_restricted_access(true);
		return opt;
	}


	/// @brief OptionCollection[ AnyOptionKey ]
	inline
	Option const &
	operator []( AnyOptionKey const & key ) const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		Option const & opt( anys_[ key ] );
		opt.check_restricted_access(true);
		return opt;
	}


	/// @brief OptionCollection[ AnyOptionKey ]
	inline
	Option &
	operator []( AnyOptionKey const & key )
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		Option & opt( anys_[ key ] );
		opt.check_restricted_access(true);
		return opt;
	}


	/// @brief OptionCollection[ AnyOptionKey ] with option type template argument
	template< typename OptionType >
	inline
	OptionType const &
	operator []( AnyOptionKey const & key ) const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		OptionType const & opt(
			utility::down_cast< OptionType const & >( anys_[ key ] ));
		opt.check_restricted_access(true);
		return opt;
	}


	/// @brief OptionCollection[ AnyOptionKey ] with option type template argument
	template< typename OptionType >
	inline
	OptionType &
	operator []( AnyOptionKey const & key )
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		OptionType & opt(
			utility::down_cast< OptionType const & >( anys_[ key ] ));
		opt.check_restricted_access(true);
		return opt;
	}


	/// @brief OptionCollection[ BooleanVectorOptionKey ]
	inline
	BooleanVectorOption const &
	operator []( BooleanVectorOptionKey const & key ) const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		BooleanVectorOption const & opt( boolean_vectors_[ key ] );
		opt.check_restricted_access(true);
		return opt;
	}


	/// @brief OptionCollection[ BooleanVectorOptionKey ]
	inline
	BooleanVectorOption &
	operator []( BooleanVectorOptionKey const & key )
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		BooleanVectorOption & opt( boolean_vectors_[ key ] );
		opt.check_restricted_access(true);
		return opt;
	}


	/// @brief OptionCollection[ IntegerVectorOptionKey ]
	inline
	IntegerVectorOption const &
	operator []( IntegerVectorOptionKey const & key ) const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		IntegerVectorOption const & opt( integer_vectors_[ key ] );
		opt.check_restricted_access(true);
		return opt;
	}


	/// @brief OptionCollection[ IntegerVectorOptionKey ]
	inline
	IntegerVectorOption &
	operator []( IntegerVectorOptionKey const & key )
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		IntegerVectorOption & opt( integer_vectors_[ key ] );
		opt.check_restricted_access(true);
		return opt;
	}


	/// @brief OptionCollection[ RealVectorOptionKey ]
	inline
	RealVectorOption const &
	operator []( RealVectorOptionKey const & key ) const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		RealVectorOption const & opt( real_vectors_[ key ] );
		opt.check_restricted_access(true);
		return opt;
	}


	/// @brief OptionCollection[ RealVectorOptionKey ]
	inline
	RealVectorOption &
	operator []( RealVectorOptionKey const & key )
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		RealVectorOption & opt( real_vectors_[ key ] );
		opt.check_restricted_access(true);
		return opt;
	}

	/// @brief OptionCollection[ ResidueChainVectorOptionKey ]
	inline
	ResidueChainVectorOption const &
	operator []( ResidueChainVectorOptionKey const & key ) const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		ResidueChainVectorOption const & opt( residue_chain_vectors_[ key ] );
		opt.check_restricted_access(true);
		return opt;
	}


	/// @brief OptionCollection[ ResidueChainVectorOptionKey ]
	inline
	ResidueChainVectorOption &
	operator []( ResidueChainVectorOptionKey const & key )
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		ResidueChainVectorOption & opt( residue_chain_vectors_[ key ] );
		opt.check_restricted_access(true);
		return opt;
	}


	/// @brief OptionCollection[ StringVectorOptionKey ]
	inline
	StringVectorOption const &
	operator []( StringVectorOptionKey const & key ) const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		StringVectorOption const & opt( string_vectors_[ key ] );
		opt.check_restricted_access(true);
		return opt;
	}


	/// @brief OptionCollection[ StringVectorOptionKey ]
	inline
	StringVectorOption &
	operator []( StringVectorOptionKey const & key )
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		StringVectorOption & opt( string_vectors_[ key ] );
		opt.check_restricted_access(true);
		return opt;
	}


	/// @brief OptionCollection[ FileVectorOptionKey ]
	inline
	FileVectorOption const &
	operator []( FileVectorOptionKey const & key ) const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		FileVectorOption const & opt( file_vectors_[ key ] );
		opt.check_restricted_access(true);
		return opt;
	}


	/// @brief OptionCollection[ FileVectorOptionKey ]
	inline
	FileVectorOption &
	operator []( FileVectorOptionKey const & key )
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		FileVectorOption & opt( file_vectors_[ key ] );
		opt.check_restricted_access(true);
		return opt;
	}


	/// @brief OptionCollection[ PathVectorOptionKey ]
	inline
	PathVectorOption const &
	operator []( PathVectorOptionKey const & key ) const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		PathVectorOption const & opt( path_vectors_[ key ] );
		opt.check_restricted_access(true);
		return opt;
	}


	/// @brief OptionCollection[ PathVectorOptionKey ]
	inline
	PathVectorOption &
	operator []( PathVectorOptionKey const & key )
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		PathVectorOption & opt( path_vectors_[ key ] );
		opt.check_restricted_access(true);
		return opt;
	}


	/// @brief OptionCollection[ AnyVectorOptionKey ]
	inline
	VectorOption const &
	operator []( AnyVectorOptionKey const & key ) const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		VectorOption const & opt( any_vectors_[ key ] );
		opt.check_restricted_access(true);
		return opt;
	}


	/// @brief OptionCollection[ AnyVectorOptionKey ]
	inline
	VectorOption &
	operator []( AnyVectorOptionKey const & key )
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		VectorOption & opt( any_vectors_[ key ] );
		opt.check_restricted_access(true);
		return opt;
	}


	/// @brief OptionCollection[ AnyVectorOptionKey ] with option type template argument
	template< typename VectorOptionType >
	inline
	VectorOptionType const &
	operator []( AnyVectorOptionKey const & key ) const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		VectorOptionType const & opt(
			utility::down_cast< VectorOptionType const & >( any_vectors_[ key ] ));
		opt.check_restricted_access(true);
		return opt;
	}


	/// @brief OptionCollection[ AnyVectorOptionKey ] with option type template argument
	template< typename VectorOptionType >
	inline
	VectorOptionType &
	operator []( AnyVectorOptionKey const & key )
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		VectorOptionType & opt(
			utility::down_cast< VectorOptionType const & >( any_vectors_[ key ] ));
		opt.check_restricted_access(true);
		return opt;
	}


	/// @brief OptionCollection[ OptionKey ]
	inline
	Option const &
	operator []( OptionKey const & key ) const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		using utility::down_cast;
		runtime_assert( all_.has( key ) );
		Option const * opt;
		switch ( all_[ key ] ) {
		case BOOLEAN_OPTION :
			opt = &booleans_[ down_cast< BooleanOptionKey const & >( key ) ];
			break;
		case INTEGER_OPTION :
			opt = &integers_[ down_cast< IntegerOptionKey const & >( key ) ];
			break;
		case REAL_OPTION :
			opt = &reals_[ down_cast< RealOptionKey const & >( key ) ];
			break;
		case STRING_OPTION :
			opt = &strings_[ down_cast< StringOptionKey const & >( key ) ];
			break;
		case FILE_OPTION :
			opt = &files_[ down_cast< FileOptionKey const & >( key ) ];
			break;
		case PATH_OPTION :
			opt = &paths_[ down_cast< PathOptionKey const & >( key ) ];
			break;
		case ANY_OPTION :
			opt = anys_[ down_cast< AnyOptionKey const & >( key ) ];
			break;
		case BOOLEAN_VECTOR_OPTION :
			opt = &boolean_vectors_[ down_cast< BooleanVectorOptionKey const & >( key ) ];
			break;
		case INTEGER_VECTOR_OPTION :
			opt = &integer_vectors_[ down_cast< IntegerVectorOptionKey const & >( key ) ];
			break;
		case REAL_VECTOR_OPTION :
			opt = &real_vectors_[ down_cast< RealVectorOptionKey const & >( key ) ];
			break;
		case RESIDUE_CHAIN_VECTOR_OPTION :
			opt = &residue_chain_vectors_[ down_cast< ResidueChainVectorOptionKey const & >( key ) ];
			break;
		case STRING_VECTOR_OPTION :
			opt = &string_vectors_[ down_cast< StringVectorOptionKey const & >( key ) ];
			break;
		case FILE_VECTOR_OPTION :
			opt = &file_vectors_[ down_cast< FileVectorOptionKey const & >( key ) ];
			break;
		case PATH_VECTOR_OPTION :
			opt = &path_vectors_[ down_cast< PathVectorOptionKey const & >( key ) ];
			break;
		case ANY_VECTOR_OPTION :
			opt = any_vectors_[ down_cast< AnyVectorOptionKey const & >( key ) ];
			break;
		case UNKNOWN_OPTION : // Joined with default case
		default :
			std::cerr << "Offending key is: " << all_[ key ] << "  " << (int)all_[ key ] << std::endl;
			utility_exit_with_message("Cannot handle option key type."); // Shouldn't get here
			return *booleans_.begin(); // Keep compiler happy
		}

		opt->check_restricted_access(true);
		return *opt;
	}


	/// @brief OptionCollection[ OptionKey ]
	inline
	Option &
	operator []( OptionKey const & key )
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		using utility::down_cast;
		runtime_assert( all_.has( key ) );
		Option * opt;
		switch ( all_[ key ] ) {
		case BOOLEAN_OPTION :
			opt = &booleans_[ down_cast< BooleanOptionKey const & >( key ) ];
			break;
		case INTEGER_OPTION :
			opt = &integers_[ down_cast< IntegerOptionKey const & >( key ) ];
			break;
		case REAL_OPTION :
			opt = &reals_[ down_cast< RealOptionKey const & >( key ) ];
			break;
		case STRING_OPTION :
			opt = &strings_[ down_cast< StringOptionKey const & >( key ) ];
			break;
		case FILE_OPTION :
			opt = &files_[ down_cast< FileOptionKey const & >( key ) ];
			break;
		case PATH_OPTION :
			opt = &paths_[ down_cast< PathOptionKey const & >( key ) ];
			break;
		case ANY_OPTION :
			opt = anys_[ down_cast< AnyOptionKey const & >( key ) ];
			break;
		case BOOLEAN_VECTOR_OPTION :
			opt = &boolean_vectors_[ down_cast< BooleanVectorOptionKey const & >( key ) ];
			break;
		case INTEGER_VECTOR_OPTION :
			opt = &integer_vectors_[ down_cast< IntegerVectorOptionKey const & >( key ) ];
			break;
		case REAL_VECTOR_OPTION :
			opt = &real_vectors_[ down_cast< RealVectorOptionKey const & >( key ) ];
			break;
		case RESIDUE_CHAIN_VECTOR_OPTION :
			opt = &residue_chain_vectors_[ down_cast< ResidueChainVectorOptionKey const & >( key ) ];
			break;
		case STRING_VECTOR_OPTION :
			opt = &string_vectors_[ down_cast< StringVectorOptionKey const & >( key ) ];
			break;
		case FILE_VECTOR_OPTION :
			opt = &file_vectors_[ down_cast< FileVectorOptionKey const & >( key ) ];
			break;
		case PATH_VECTOR_OPTION :
			opt = &path_vectors_[ down_cast< PathVectorOptionKey const & >( key ) ];
			break;
		case ANY_VECTOR_OPTION :
			opt = any_vectors_[ down_cast< AnyVectorOptionKey const & >( key ) ];
			break;
		case UNKNOWN_OPTION : //Joined with default case
		default :
			std::cerr << "Offending key is: " << all_[ key ] << "  " << (int)all_[ key ] << std::endl;
			utility_exit_with_message("Cannot handle key type."); // Shouldn't get here
			return *booleans_.begin(); // Keep compiler happy
		}
		opt->check_restricted_access(true);
		return *opt;
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
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		BooleanOption const & opt(booleans_[ key ]);
		opt.check_restricted_access(check_restricted_access);
		return opt;
	}


	/// @brief Option by BooleanOptionKey
	inline
	BooleanOption &
	operator ()( BooleanOptionKey const & key, bool check_restricted_access=true )
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		BooleanOption & opt( booleans_[ key ]);
		opt.check_restricted_access(check_restricted_access);
		return opt;
	}


	/// @brief Option by IntegerOptionKey
	inline
	IntegerOption const &
	operator ()( IntegerOptionKey const & key, bool check_restricted_access=true ) const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		IntegerOption const & opt( integers_[ key ]);
		opt.check_restricted_access(check_restricted_access);
		return opt;
	}


	/// @brief Option by IntegerOptionKey
	inline
	IntegerOption &
	operator ()( IntegerOptionKey const & key, bool check_restricted_access=true )
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		IntegerOption & opt( integers_[ key ]);
		opt.check_restricted_access(check_restricted_access);
		return opt;
	}


	/// @brief Option by RealOptionKey
	inline
	RealOption const &
	operator ()( RealOptionKey const & key, bool check_restricted_access=true ) const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		RealOption const & opt( reals_[ key ]);
		opt.check_restricted_access(check_restricted_access);
		return opt;
	}


	/// @brief Option by RealOptionKey
	inline
	RealOption &
	operator ()( RealOptionKey const & key, bool check_restricted_access=true )
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		RealOption & opt( reals_[ key ]);
		opt.check_restricted_access(check_restricted_access);
		return opt;
	}


	/// @brief Option by StringOptionKey
	inline
	StringOption const &
	operator ()( StringOptionKey const & key, bool check_restricted_access=true ) const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		StringOption const & opt( strings_[ key ]);
		opt.check_restricted_access(check_restricted_access);
		return opt;
	}


	/// @brief Option by StringOptionKey
	inline
	StringOption &
	operator ()( StringOptionKey const & key, bool check_restricted_access=true )
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		StringOption & opt( strings_[ key ]);
		opt.check_restricted_access(check_restricted_access);
		return opt;
	}


	/// @brief Option by FileOptionKey
	inline
	FileOption const &
	operator ()( FileOptionKey const & key, bool check_restricted_access=true ) const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		FileOption const & opt( files_[ key ]);
		opt.check_restricted_access(check_restricted_access);
		return opt;
	}


	/// @brief Option by FileOptionKey
	inline
	FileOption &
	operator ()( FileOptionKey const & key, bool check_restricted_access=true )
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		FileOption & opt( files_[ key ]);
		opt.check_restricted_access(check_restricted_access);
		return opt;
	}


	/// @brief Option by PathOptionKey
	inline
	PathOption const &
	operator ()( PathOptionKey const & key, bool check_restricted_access=true ) const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		PathOption const & opt( paths_[ key ]);
		opt.check_restricted_access(check_restricted_access);
		return opt;
	}


	/// @brief Option by PathOptionKey
	inline
	PathOption &
	operator ()( PathOptionKey const & key, bool check_restricted_access=true )
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		PathOption & opt( paths_[ key ]);
		opt.check_restricted_access(check_restricted_access);
		return opt;
	}


	/// @brief Option by AnyOptionKey
	inline
	Option const &
	operator ()( AnyOptionKey const & key, bool check_restricted_access=true ) const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		Option const & opt( anys_[ key ]);
		opt.check_restricted_access(check_restricted_access);
		return opt;
	}


	/// @brief Option by AnyOptionKey
	inline
	Option &
	operator ()( AnyOptionKey const & key, bool check_restricted_access=true )
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		Option & opt( anys_[ key ]);
		opt.check_restricted_access(check_restricted_access);
		return opt;
	}


	/// @brief Option by AnyOptionKey with option type template argument
	template< typename OptionType >
	inline
	OptionType const &
	operator ()( AnyOptionKey const & key, bool check_restricted_access=true ) const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		OptionType const & opt(
			utility::down_cast< OptionType const & >( anys_[ key ] ));
		opt.check_restricted_access(check_restricted_access);
		return opt;
	}


	/// @brief Option by AnyOptionKey with option type template argument
	template< typename OptionType >
	inline
	OptionType &
	operator ()( AnyOptionKey const & key, bool check_restricted_access=true )
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		OptionType & opt(
			utility::down_cast< OptionType const & >( anys_[ key ] ));
		opt.check_restricted_access(check_restricted_access);
		return opt;
	}


	/// @brief VectorOption by BooleanVectorOptionKey
	inline
	BooleanVectorOption const &
	operator ()( BooleanVectorOptionKey const & key, bool check_restricted_access=true ) const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		BooleanVectorOption const & opt( boolean_vectors_[ key ]);
		opt.check_restricted_access(check_restricted_access);
		return opt;
	}


	/// @brief VectorOption by BooleanVectorOptionKey
	inline
	BooleanVectorOption &
	operator ()( BooleanVectorOptionKey const & key, bool check_restricted_access=true )
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		BooleanVectorOption & opt( boolean_vectors_[ key ]);
		opt.check_restricted_access(check_restricted_access);
		return opt;
	}


	/// @brief VectorOption by IntegerVectorOptionKey
	inline
	IntegerVectorOption const &
	operator ()( IntegerVectorOptionKey const & key, bool check_restricted_access=true ) const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		IntegerVectorOption const & opt( integer_vectors_[ key ]);
		opt.check_restricted_access(check_restricted_access);
		return opt;
	}


	/// @brief VectorOption by IntegerVectorOptionKey
	inline
	IntegerVectorOption &
	operator ()( IntegerVectorOptionKey const & key, bool check_restricted_access=true )
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		IntegerVectorOption & opt( integer_vectors_[ key ]);
		opt.check_restricted_access(check_restricted_access);
		return opt;
	}


	/// @brief VectorOption by RealVectorOptionKey
	inline
	RealVectorOption const &
	operator ()( RealVectorOptionKey const & key, bool check_restricted_access=true ) const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		RealVectorOption const & opt( real_vectors_[ key ]);
		opt.check_restricted_access(check_restricted_access);
		return opt;
	}


	/// @brief VectorOption by RealVectorOptionKey
	inline
	RealVectorOption &
	operator ()( RealVectorOptionKey const & key, bool check_restricted_access=true )
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		RealVectorOption & opt( real_vectors_[ key ]);
		opt.check_restricted_access(check_restricted_access);
		return opt;
	}


	/// @brief VectorOption by ResidueChainVectorOptionKey
	inline
	ResidueChainVectorOption const &
	operator ()( ResidueChainVectorOptionKey const & key, bool check_restricted_access=true ) const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		ResidueChainVectorOption const & opt( residue_chain_vectors_[ key ]);
		opt.check_restricted_access(check_restricted_access);
		return opt;
	}


	/// @brief VectorOption by ResidueChainVectorOptionKey
	inline
	ResidueChainVectorOption &
	operator ()( ResidueChainVectorOptionKey const & key, bool check_restricted_access=true )
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		ResidueChainVectorOption & opt( residue_chain_vectors_[ key ]);
		opt.check_restricted_access(check_restricted_access);
		return opt;
	}


	/// @brief VectorOption by StringVectorOptionKey
	inline
	StringVectorOption const &
	operator ()( StringVectorOptionKey const & key, bool check_restricted_access=true ) const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		StringVectorOption const & opt( string_vectors_[ key ]);
		opt.check_restricted_access(check_restricted_access);
		return opt;
	}


	/// @brief VectorOption by StringVectorOptionKey
	inline
	StringVectorOption &
	operator ()( StringVectorOptionKey const & key, bool check_restricted_access=true )
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		StringVectorOption & opt( string_vectors_[ key ]);
		opt.check_restricted_access(check_restricted_access);
		return opt;
	}


	/// @brief VectorOption by FileVectorOptionKey
	inline
	FileVectorOption const &
	operator ()( FileVectorOptionKey const & key, bool check_restricted_access=true ) const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		FileVectorOption const & opt( file_vectors_[ key ]);
		opt.check_restricted_access(check_restricted_access);
		return opt;
	}


	/// @brief VectorOption by FileVectorOptionKey
	inline
	FileVectorOption &
	operator ()( FileVectorOptionKey const & key, bool check_restricted_access=true )
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		FileVectorOption & opt( file_vectors_[ key ]);
		opt.check_restricted_access(check_restricted_access);
		return opt;
	}


	/// @brief VectorOption by PathVectorOptionKey
	inline
	PathVectorOption const &
	operator ()( PathVectorOptionKey const & key, bool check_restricted_access=true ) const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		PathVectorOption const & opt( path_vectors_[ key ]);
		opt.check_restricted_access(check_restricted_access);
		return opt;
	}


	/// @brief VectorOption by PathVectorOptionKey
	inline
	PathVectorOption &
	operator ()( PathVectorOptionKey const & key, bool check_restricted_access=true )
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		PathVectorOption & opt( path_vectors_[ key ]);
		opt.check_restricted_access(check_restricted_access);
		return opt;
	}


	/// @brief VectorOption by AnyVectorOptionKey
	inline
	VectorOption const &
	operator ()( AnyVectorOptionKey const & key, bool check_restricted_access=true ) const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		VectorOption const & opt(any_vectors_[ key ] );
		opt.check_restricted_access(check_restricted_access);
		return opt;
	}


	/// @brief VectorOption by AnyVectorOptionKey
	inline
	VectorOption &
	operator ()( AnyVectorOptionKey const & key, bool check_restricted_access=true )
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		VectorOption & opt( any_vectors_[ key ] );
		opt.check_restricted_access(check_restricted_access);
		return opt;
	}


	/// @brief VectorOption by AnyVectorOptionKey with option type template argument
	template< typename VectorOptionType >
	inline
	VectorOptionType const &
	operator ()( AnyVectorOptionKey const & key, bool check_restricted_access=true ) const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		VectorOptionType const & opt(
			utility::down_cast< VectorOptionType const & >( any_vectors_[ key ] ));
		opt.check_restricted_access(check_restricted_access);
		return opt;
	}


	/// @brief VectorOption by AnyVectorOptionKey with option type template argument
	template< typename VectorOptionType >
	inline
	VectorOptionType &
	operator ()( AnyVectorOptionKey const & key, bool check_restricted_access=true )
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		VectorOptionType & opt(
			utility::down_cast< VectorOptionType & >( any_vectors_[ key ] ));
		opt.check_restricted_access(check_restricted_access);
		return opt;
	}


	/// @brief Option by OptionKey
	inline
	Option const &
	operator ()( OptionKey const & key, bool check_restricted_access=true ) const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		using utility::down_cast;
		runtime_assert( all_.has( key ) );
		Option const * opt;
		switch ( all_[ key ] ) {
		case BOOLEAN_OPTION :
			opt = &booleans_[ down_cast< BooleanOptionKey const & >( key ) ];
			break;
		case INTEGER_OPTION :
			opt = &integers_[ down_cast< IntegerOptionKey const & >( key ) ];
			break;
		case REAL_OPTION :
			opt = &reals_[ down_cast< RealOptionKey const & >( key ) ];
			break;
		case STRING_OPTION :
			opt = &strings_[ down_cast< StringOptionKey const & >( key ) ];
			break;
		case FILE_OPTION :
			opt = &files_[ down_cast< FileOptionKey const & >( key ) ];
			break;
		case PATH_OPTION :
			opt = &paths_[ down_cast< PathOptionKey const & >( key ) ];
			break;
		case ANY_OPTION :
			opt = anys_[ down_cast< AnyOptionKey const & >( key ) ];
			break;
		case BOOLEAN_VECTOR_OPTION :
			opt = &boolean_vectors_[ down_cast< BooleanVectorOptionKey const & >( key ) ];
			break;
		case INTEGER_VECTOR_OPTION :
			opt = &integer_vectors_[ down_cast< IntegerVectorOptionKey const & >( key ) ];
			break;
		case REAL_VECTOR_OPTION :
			opt = &real_vectors_[ down_cast< RealVectorOptionKey const & >( key ) ];
			break;
		case RESIDUE_CHAIN_VECTOR_OPTION :
			opt = &residue_chain_vectors_[ down_cast< ResidueChainVectorOptionKey const & >( key ) ];
			break;
		case STRING_VECTOR_OPTION :
			opt = &string_vectors_[ down_cast< StringVectorOptionKey const & >( key ) ];
			break;
		case FILE_VECTOR_OPTION :
			opt = &file_vectors_[ down_cast< FileVectorOptionKey const & >( key ) ];
			break;
		case PATH_VECTOR_OPTION :
			opt = &path_vectors_[ down_cast< PathVectorOptionKey const & >( key ) ];
			break;
		case ANY_VECTOR_OPTION :
			opt = any_vectors_[ down_cast< AnyVectorOptionKey const & >( key ) ];
			break;
		case UNKNOWN_OPTION : //Joined with default case
		default :
			std::cerr << "Offending key is: " << all_[ key ] << "  " << (int)all_[ key ] << std::endl;
			utility_exit_with_message("Cannot handle key type."); // Shouldn't get here
			return *booleans_.begin(); // Keep compiler happy
		}

		opt->check_restricted_access(check_restricted_access);
		return *opt;
	}


	/// @brief Option by OptionKey
	inline
	Option &
	operator ()( OptionKey const & key, bool check_restricted_access=true )
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		using utility::down_cast;
		runtime_assert( all_.has( key ) );
		Option * opt;
		switch ( all_[ key ] ) {
		case BOOLEAN_OPTION :
			opt = &booleans_[ down_cast< BooleanOptionKey const & >( key ) ];
			break;
		case INTEGER_OPTION :
			opt = &integers_[ down_cast< IntegerOptionKey const & >( key ) ];
			break;
		case REAL_OPTION :
			opt = &reals_[ down_cast< RealOptionKey const & >( key ) ];
			break;
		case STRING_OPTION :
			opt = &strings_[ down_cast< StringOptionKey const & >( key ) ];
			break;
		case FILE_OPTION :
			opt = &files_[ down_cast< FileOptionKey const & >( key ) ];
			break;
		case PATH_OPTION :
			opt = &paths_[ down_cast< PathOptionKey const & >( key ) ];
			break;
		case ANY_OPTION :
			opt = anys_[ down_cast< AnyOptionKey const & >( key ) ];
			break;
		case BOOLEAN_VECTOR_OPTION :
			opt = &boolean_vectors_[ down_cast< BooleanVectorOptionKey const & >( key ) ];
			break;
		case INTEGER_VECTOR_OPTION :
			opt = &integer_vectors_[ down_cast< IntegerVectorOptionKey const & >( key ) ];
			break;
		case REAL_VECTOR_OPTION :
			opt = &real_vectors_[ down_cast< RealVectorOptionKey const & >( key ) ];
			break;
		case RESIDUE_CHAIN_VECTOR_OPTION :
			opt = &residue_chain_vectors_[ down_cast< ResidueChainVectorOptionKey const & >( key ) ];
			break;
		case STRING_VECTOR_OPTION :
			opt = &string_vectors_[ down_cast< StringVectorOptionKey const & >( key ) ];
			break;
		case FILE_VECTOR_OPTION :
			opt = &file_vectors_[ down_cast< FileVectorOptionKey const & >( key ) ];
			break;
		case PATH_VECTOR_OPTION :
			opt = &path_vectors_[ down_cast< PathVectorOptionKey const & >( key ) ];
			break;
		case ANY_VECTOR_OPTION :
			opt = any_vectors_[ down_cast< AnyVectorOptionKey const & >( key ) ];
			break;
		case UNKNOWN_OPTION : //Joined with default case
		default :
			std::cerr << "Offending key is: " << all_[ key ] << "  " << (int)all_[ key ] << std::endl;
			utility_exit_with_message("Cannot handle key type."); // Shouldn't get here
			return *booleans_.begin(); // Keep compiler happy
		}
		opt->check_restricted_access(check_restricted_access);
		return *opt;
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


public: //Methods, some static:

	/// @brief add OptionKey to list of application relevant options
	void
	add_relevant(const OptionKey & key)
	{
#ifdef MULTI_THREADED
		utility::thread::WriteLockGuard writelock( relevant_accessed_unused_mutex_ );
#endif
		relevant_.push_back( &key );
	}

	/// @brief Checks if option has been registered as relevant.
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

	/// @brief Remove underscores and lowercase string
	static
	std::string
	lower_no_under( std::string const & instring );

	/// @brief Add to the set possible typos for a given word
	/// @details These are strings with Levenshtein distance of 1
	/// - Based off of the approach from http://norvig.com/spell-correct.html
	static
	void
	add_edits( std::set<std::string> & items );

	/// @brief Find a user-specified option key in a command line context
	/// @note  Searches up the context to find a match
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
		std::string::size_type indent = 10, // Width to indent continuation lines
		std::string::size_type width = 80, // Column width to wrap at [80]
		std::string const & header_for_extra_lines = ""
	);


	/// @brief modify 'show_accessed_options' flag;
	void
	set_show_accessed_options_flag(bool v) {
#ifdef MULTI_THREADED
		utility::thread::WriteLockGuard writelock( relevant_accessed_unused_mutex_ );
#endif
		show_accessed_options_ = v;
	}

	/// @brief modify 'show_unused_options' flag;
	void
	set_show_unused_options_flag(bool v) {
#ifdef MULTI_THREADED
		utility::thread::WriteLockGuard writelock( relevant_accessed_unused_mutex_ );
#endif
		show_unused_options_ = v;
	}

private: // Fields


#ifdef MULTI_THREADED
	/// @brief Mutex for locking the OptionsCollection for read or write:
	mutable utility::thread::ReadWriteMutex mutex_;
#endif

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

	/// @brief List of application-relevant options
	std::vector< OptionKey const *> relevant_;

	/// @brief keep a copy of argv around in case people want to
	/// get at it elsewhere in the code.
	std::string argv_copy_;

#ifdef MULTI_THREADED
	/// @brief Mutex that protects relevant_, show_accessed_options_, and show_unused_options_
	mutable utility::thread::ReadWriteMutex relevant_accessed_unused_mutex_;
#endif

	/// @brief Flag indicating that list of accessed option should be printed when destructor of
	///             OptionCollection is called.
	///        This flag is false by default.
	bool show_accessed_options_ = false;

	/// @brief Flag indicating that list of user specified but inaccessed options should be printed
	///             when destructor of OptionCollection is called.
	///        This flag is false by default.
	bool show_unused_options_ = false;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; // OptionCollection


} // namespace options
} // namespace utility


#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( utility_options_OptionCollection )
#endif // SERIALIZATION


#endif // INCLUDED_utility_options_OptionCollection_HH
