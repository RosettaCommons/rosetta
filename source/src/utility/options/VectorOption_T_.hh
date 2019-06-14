// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/options/VectorOption_T_.hh
/// @brief  Program vector-valued option abstract base class
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author Modified by Sergey Lyskov (Sergey.Lyskov@jhu.edu)
/// @author Modified by Rhiju Das (rhiju@stanford.edu)
/// @author Modified by Vikram K. Mulligan (vmulligan@flatironinstitute.org) for thread-safety.


#ifndef INCLUDED_utility_options_VectorOption_T__HH
#define INCLUDED_utility_options_VectorOption_T__HH


// Unit headers
#include <utility/options/VectorOption_T_.fwd.hh>

// Package headers
#include <utility/options/VectorOption.hh>
#include <utility/Bound.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>

// ObjexxFCL headers
#include <ObjexxFCL/string.functions.hh>

// C++ headers
#include <utility/assert.hh>
#include <cstdlib>
#include <iostream>
#include <set>
#include <sstream>

#ifdef SERIALIZATION
#include <utility/vector1.srlz.hh>
#include <cereal/types/set.hpp>
#include <cereal/types/string.hpp>
#endif

#ifdef MULTI_THREADED
#include <utility/thread/ReadWriteMutex.hh>
#endif

namespace utility {
namespace options {


/// @brief Program vector-valued option abstract base class
template< typename K, typename T >
class VectorOption_T_ :
	public VectorOption
{


private: // Types


	typedef  VectorOption  Super;

	typedef  std::set< T >  LegalValues;
	typedef  Bound< T >  LegalBound;


public: // Types


	typedef  utility::vector1< T >  Values;

	// STL/boost style
	typedef  K  key_type;
	typedef  T  value_type;
	typedef  typename Values::const_iterator  const_iterator;
	typedef  typename Values::iterator  iterator;

	// Project style
	typedef  K  Key;
	typedef  T  Value;
	typedef  typename Values::const_iterator  ConstIterator;
	typedef  typename Values::iterator  Iterator;


protected: // Creation


	/// @brief Default constructor
	inline
	VectorOption_T_() :
		n_( 0 ),
		n_lower_( 0 ),
		n_upper_( 0 ),
		default_state_( INACTIVE ),
		state_( INACTIVE )
	{}


	/// @brief Copy constructor
	inline
	VectorOption_T_( VectorOption_T_ const & option ) :
		Super( option )
#ifdef MULTI_THREADED
		,
		mutex_()
#endif
	{
		(*this) = option;
	}


	/// @brief Key + description constructor
	inline
	VectorOption_T_(
		Key const & key_a,
		std::string const & description_a
	) :
#ifdef MULTI_THREADED
		mutex_(),
#endif
		key_( key_a ),
		description_( description_a ),
		short_description_( description_a ),
		n_( 0 ),
		n_lower_( 0 ),
		n_upper_( 0 ),
		default_state_( INACTIVE ),
		state_( INACTIVE )
	{}


public: // Creation


	/// @brief Clone this

	VectorOption_T_ *
	clone() const override = 0;


	/// @brief Destructor
	inline

	~VectorOption_T_() {}


protected: // Assignment


	/// @brief Copy assignment
	/// @details Threadsafe.
	inline
	VectorOption_T_ &
	operator =( VectorOption_T_ const & option )
	{
		Option::operator=(option);
		if ( this != &option ) {
#ifdef MULTI_THREADED
			utility::thread::PairedReadLockWriteLockGuard( option.mutex_ /*Gets read-lock.*/, mutex_ /*Gets write-lock.*/ );
#endif
			key_ = option.key_;
			description_ = option.description_;
			short_description_ = option.short_description_;
			legal_ = option.legal_;
			lower_ = option.lower_;
			upper_ = option.upper_;
			n_ = option.n_;
			n_lower_ = option.n_lower_;
			n_upper_ = option.n_upper_;
			default_state_ = option.default_state_;
			default_value_ = option.default_value_;
			state_ = option.state_;
			value_ = option.value_;
		}
		return *this;
	}

public: // copying

	/// @brief Copy operation
	/// #details Relies on assignment operator; does not lock mutexes itself.
	void copy_from( Option const & other ) override {

		debug_assert( dynamic_cast< VectorOption_T_ const * > ( & other ));

		VectorOption_T_ const & vect_opt_other =
			static_cast< VectorOption_T_ const & > ( other );

		*this = vect_opt_other; // rely on assignment operator
	}


public: // Conversion


	/// @brief Value conversion
	inline
	operator Values const &() const
	{
		been_accessed();
		if ( state_ == INACTIVE ) inactive_error();
		return value_;
	}

	/* Commenting this out for clear separation of read/write operators
	/// @brief Value conversion
	inline
	operator Values &()
	{
	if ( state_ == INACTIVE ) inactive_error();
	return value_;
	} */

	/// @brief Iterator access for range for loops
	/// @details Not threadsafe.
	typename Values::const_iterator
	begin() const {
		return (*this)().begin();
	}

	/// @brief Iterator access for range for loops
	/// @details Not threadsafe.
	typename Values::const_iterator
	end() const {
		return (*this)().end();
	}

public: // Methods


	/// @brief Activate
	inline
	VectorOption_T_ &
	activate() override
	{
#ifdef MULTI_THREADED
		utility::thread::WriteLockGuard writelock( mutex_ );
#endif
		state_ = USER;
		return *this;
	}


	/// @brief Deactivate
	inline
	VectorOption_T_ &
	deactivate() override
	{
#ifdef MULTI_THREADED
		utility::thread::WriteLockGuard writelock( mutex_ );
#endif
		state_ = INACTIVE;
		return *this;
	}


	/// @brief Set to default value, if any
	inline
	VectorOption_T_ &
	to_default() override
	{
#ifdef MULTI_THREADED
		utility::thread::WriteLockGuard writelock( mutex_ );
#endif
		if ( default_state_ == DEFAULT ) {
			state_ = DEFAULT;
			value_ = default_value_;
		}
		return *this;
	}


	/// @brief Clear
	inline
	VectorOption_T_ &
	clear() override
	{
#ifdef MULTI_THREADED
		utility::thread::WriteLockGuard writelock( mutex_ );
#endif
		n_ = 0;
		n_lower_ = 0;
		n_upper_ = 0;
		default_state_ = INACTIVE;
		default_value_.clear();
		state_ = INACTIVE;
		value_.clear();
		return *this;
	}


	/// @brief Add a legal value
	inline
	VectorOption_T_ &
	legal( Value const & value_a )
	{
#ifdef MULTI_THREADED
		utility::thread::WriteLockGuard writelock( mutex_ );
#endif
		legal_.insert( value_a );
		return *this;
	}

	/// @brief Set a short description
	inline
	VectorOption_T_ &
	shortd( std::string const & s)
	{
#ifdef MULTI_THREADED
		utility::thread::WriteLockGuard writelock( mutex_ );
#endif
		short_description_ = s ;
		return *this;
	}

	/// @brief Set a lower bound
	inline
	VectorOption_T_ &
	lower( Value const & value_a )
	{
#ifdef MULTI_THREADED
		utility::thread::WriteLockGuard writelock( mutex_ );
#endif
		lower_( value_a );
		return *this;
	}


	/// @brief Set a strict lower bound
	inline
	VectorOption_T_ &
	strict_lower( Value const & value_a )
	{
#ifdef MULTI_THREADED
		utility::thread::WriteLockGuard writelock( mutex_ );
#endif
		lower_( value_a, true );
		return *this;
	}


	/// @brief Set an upper bound
	inline
	VectorOption_T_ &
	upper( Value const & value_a )
	{
#ifdef MULTI_THREADED
		utility::thread::WriteLockGuard writelock( mutex_ );
#endif
		upper_( value_a );
		return *this;
	}


	/// @brief Set a strict upper bound
	inline
	VectorOption_T_ &
	strict_upper( Value const & value_a )
	{
#ifdef MULTI_THREADED
		utility::thread::WriteLockGuard writelock( mutex_ );
#endif
		upper_( value_a, true );
		return *this;
	}


	/// @brief Fixed number of values required assignment
	inline
	VectorOption_T_ &
	n( Size const n_a ) override
	{
#ifdef MULTI_THREADED
		utility::thread::WriteLockGuard writelock( mutex_ );
#endif
		n_ = n_a;
		return *this;
	}


	/// @brief Lower number of values allowed assignment
	inline
	VectorOption_T_ &
	n_lower( Size const n_a ) override
	{
#ifdef MULTI_THREADED
		utility::thread::WriteLockGuard writelock( mutex_ );
#endif
		n_lower_ = n_a;
		return *this;
	}


	/// @brief Upper number of values allowed assignment
	inline
	VectorOption_T_ &
	n_upper( Size const n_a ) override
	{
#ifdef MULTI_THREADED
		utility::thread::WriteLockGuard writelock( mutex_ );
#endif
		n_upper_ = n_a;
		return *this;
	}


	/// @brief Default value assignment
	inline
	virtual
	VectorOption_T_ &
	default_value( Value const & value_a )
	{
		{ //Scope for possible mutex lock
#ifdef MULTI_THREADED
			utility::thread::WriteLockGuard writelock( mutex_ );
#endif
			default_state_ = DEFAULT;
			default_value_.push_back( value_a );
			if ( ( state_ == INACTIVE ) || ( state_ == DEFAULT ) ) {
				state_ = DEFAULT;
				value_.push_back( value_a );
			}
		} //End mutex lock guard scope
		legal_default_check( value_a );
		return *this;
	}


	/// @brief Default value assignment
	inline
	virtual
	VectorOption_T_ &
	def( Value const & value_a )
	{
		{ //Scope for possible mutex lock
#ifdef MULTI_THREADED
			utility::thread::WriteLockGuard writelock( mutex_ );
#endif
			default_state_ = DEFAULT;
			default_value_.push_back( value_a );
			if ( ( state_ == INACTIVE ) || ( state_ == DEFAULT ) ) {
				state_ = DEFAULT;
				value_.push_back( value_a );
			}
		} //End mutex lock guard scope
		legal_default_check( value_a );
		return *this;
	}


	/// @brief Default value assignment to empty vector.
	inline
	virtual
	VectorOption_T_ &
	def()
	{
#ifdef MULTI_THREADED
		utility::thread::WriteLockGuard writelock( mutex_ );
#endif
		default_state_ = DEFAULT;
		if ( ( state_ == INACTIVE ) || ( state_ == DEFAULT ) ) {
			state_ = DEFAULT;
		}
		return *this;
	}


	/// @brief Default value vector assignment
	inline
	VectorOption_T_ &
	default_value( Values const & value_a )
	{
		{
#ifdef MULTI_THREADED
			utility::thread::WriteLockGuard writelock( mutex_ );
#endif
			default_state_ = DEFAULT;
			default_value_ = value_a;
			if ( ( state_ == INACTIVE ) || ( state_ == DEFAULT ) ) {
				state_ = DEFAULT;
				value_ = value_a;
			}
		}
		legal_default_check();
		return *this;
	}


	/// @brief Default value vector assignment
	inline
	VectorOption_T_ &
	def( Values const & value_a )
	{
		{
#ifdef MULTI_THREADED
			utility::thread::WriteLockGuard writelock( mutex_ );
#endif
			default_state_ = DEFAULT;
			default_value_ = value_a;
			if ( ( state_ == INACTIVE ) || ( state_ == DEFAULT ) ) {
				state_ = DEFAULT;
				value_ = value_a;
			}
		}
		legal_default_check();
		return *this;
	}


	/// @brief Value assignment from a command line string
	inline
	VectorOption_T_ &
	cl_value( std::string const & value_str ) override
	{
#ifdef MULTI_THREADED
		utility::thread::WriteLockGuard writelock( mutex_ );
#endif
		std::string const stripped_value_str( ObjexxFCL::stripped( value_str, "\"'" ) );
		if ( ! stripped_value_str.empty() ) {
			if ( state_ == DEFAULT ) value_.clear(); // Clear out the defaulted values
			state_ = USER;
			Values const vs( values_of( stripped_value_str ) );

			for ( Size i = 1; i <= vs.size(); i++ ) {
				Value v = vs[ i ];

				value_.push_back( v );
				if ( ! legal_value( v, true /*mutex_ already locked*/ ) ) {
					std::cerr << "ERROR: Illegal value specified for option -" << id()
						<< " : " << value_str << std::endl;
					std::exit( EXIT_FAILURE );
				}
			}

		}
		return *this;
	}


	/// @brief Value assignment
	inline
	virtual
	VectorOption_T_ &
	value( Value const & value_a )
	{
		{
#ifdef MULTI_THREADED
			utility::thread::WriteLockGuard writelock( mutex_ );
#endif
			if ( state_ == DEFAULT ) value_.clear(); // Clear out the defaults
			state_ = USER;
			value_.push_back( value_a );
		}
		legal_check( value_a );
		return *this;
	}

	/// @brief Add to values. Note that this just calls value() which has push_back functionality but confusing name.
	inline
	virtual
	VectorOption_T_ &
	push_back( Value const & value_a )
	{
		value( value_a ); //Threadsafe
		return *this;
	}


	/// @brief Value assignment
	inline
	virtual
	VectorOption_T_ &
	operator ()( Value const & value_a )
	{
		{
#ifdef MULTI_THREADED
			utility::thread::WriteLockGuard writelock( mutex_ );
#endif
			if ( state_ == DEFAULT ) value_.clear(); // Clear out the defaults
			state_ = USER;
			value_.push_back( value_a );
		}
		legal_check( value_a );
		return *this;
	}


	/// @brief Value vector assignment
	inline
	VectorOption_T_ &
	value( Values const & value_a )
	{
		{
#ifdef MULTI_THREADED
			utility::thread::WriteLockGuard writelock( mutex_ );
#endif
			if ( state_ == DEFAULT ) value_.clear(); // Clear out the defaults
			state_ = USER;
			value_ = value_a;
		}
		legal_check();
		return *this;
	}


	/// @brief Value vector assignment
	inline
	VectorOption_T_ &
	operator ()( Values const & value_a )
	{
		{
#ifdef MULTI_THREADED
			utility::thread::WriteLockGuard writelock( mutex_ );
#endif
			if ( state_ == DEFAULT ) value_.clear(); // Clear out the defaults
			state_ = USER;
			value_ = value_a;
		}
		legal_check();
		return *this;
	}


	/// @brief Default to another option's value
	inline
	VectorOption_T_ &
	default_to( VectorOption_T_ const & option )
	{
		bool do_default_value;
		{ //Scope for thread.
#ifdef MULTI_THREADED
			utility::thread::ReadLockGuard readlock( mutex_ );
#endif
			do_default_value = ( ( state_ == INACTIVE ) || ( state_ == DEFAULT ) );
		}
		if ( do_default_value ) {
			if ( option.active() /*Internally threadsafe*/ ) default_value( option.value() ) /*Also internally threasafe*/;
		}
		return *this;
	}


	/// @brief Legal specifications check: Report and return error state
	inline
	bool
	legal_specs_report() const override
	{
		return ( ( legal_limits_report() ) && ( legal_size_report() ) && ( legal_default_report() ) );
	}


	/// @brief Legal value limits check: Report and return error state
	inline
	bool
	legal_limits_report() const override
	{
		bool error( false );
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		if ( ( lower_.active() ) && ( upper_.active() ) ) {
			if ( ( lower_.strict() ) || ( upper_.strict() ) ) {
				if ( lower_() >= upper_() ) error = true;
			} else {
				if ( lower_() > upper_() ) error = true;
			}
			if ( error ) {
				std::cerr << "ERROR: Inconsistent lower and upper limits in option -" << id()
					<< " : " << legal_string(true) << std::endl;
			}
		}
		return ( ! error );
	}


	/// @brief Legal size limits check: Report and return error state
	inline
	bool
	legal_size_report() const override
	{
		bool error( false );
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		if ( ( n_upper_ > 0 ) && ( n_lower_ > n_upper_ ) ) error = true;
		if ( ( n_ > 0 ) && ( n_lower_ > 0 ) && ( n_ != n_lower_ ) ) error = true;
		if ( ( n_ > 0 ) && ( n_upper_ > 0 ) && ( n_ != n_upper_ ) ) error = true;
		if ( error ) {
			std::cerr << "ERROR: Inconsistent size limits in option -" << id() << std::endl;
		}
		return ( ! error );
	}


	/// @brief Legal default value check: Report and return error state
	inline
	bool
	legal_default_report() const override
	{
		bool error( false );
		if ( ! legal_default_size() ) {
			std::cerr << "ERROR: Illegal number of default values in option -" << id()
				<< " : " << default_string(true) << std::endl;
			error = true;
		}
		if ( ! legal_default_value() ) {
			std::cerr << "ERROR: Illegal default value in option -" << id()
				<< " : " << default_string(true) << std::endl;
			error = true;
		}
		return ( ! error );
	}


	/// @brief Legal default value check
	inline
	void
	legal_default_check() const override
	{
		bool error( false );
		if ( ! legal_default_size() ) {
			std::cerr << "ERROR: Illegal number of default values in option -" << id()
				<< " : " << default_string(true) << std::endl;
			error = true;
		}
		if ( ! legal_default_value() ) {
			std::cerr << "ERROR: Illegal default value in option -" << id()
				<< " : " << default_string(true) << std::endl;
			error = true;
		}
		if ( error ) std::exit( EXIT_FAILURE );
	}


	/// @brief Legal default value check
	inline
	void
	legal_default_check( Value const & value_a ) const
	{
		if ( ! legal_value( value_a ) ) {
			std::cerr << "ERROR: Illegal default value in option -" << id()
				<< " : " << value_string_of( value_a ) << std::endl;
			std::exit( EXIT_FAILURE );
		}
	}


	/// @brief Legal value check: Report and return error state
	inline
	bool
	legal_report() const override
	{
		bool error( false );
		if ( ! legal_size() ) {
			std::cerr << "ERROR: Illegal number of values specified for option -" << id()
				<< " : " << value_string() << std::endl;
			error = true;
		}
		if ( ! legal_value() ) {
			std::cerr << "ERROR: Illegal value specified for option -" << id()
				<< " : " << value_string() << std::endl;
			error = true;
		}
		return ( ! error );
	}


	/// @brief Legal value check
	inline
	void
	legal_check() const override
	{
		bool error( false );
		if ( ! legal_size() ) {
			std::cerr << "ERROR: Illegal number of values specified for option -" << id()
				<< " : " << value_string() << std::endl;
			error = true;
		}
		if ( ! legal_value() ) {
			std::cerr << "ERROR: Illegal value specified for option -" << id()
				<< " : " << value_string() << std::endl;
			error = true;
		}
		if ( error ) std::exit( EXIT_FAILURE );
	}


	/// @brief Legal value check
	inline
	void
	legal_check( Value const & value_a ) const
	{
		if ( ! legal_value( value_a ) ) {
			std::cerr << "ERROR: Illegal value specified for option -" << id()
				<< " : " << value_string_of( value_a ) << std::endl;
			std::exit( EXIT_FAILURE );
		}
	}


	/// @brief Required specified option check: Report and return error state
	inline
	bool
	specified_report() const override
	{
		bool error( false );
		if ( ! user() ) {
			std::cerr << "ERROR: Unspecified option -" << id() << " is required" << std::endl;
			error = true;
		}
		return ( ! error );
	}


	/// @brief Required specified option check
	inline
	void
	specified_check() const override
	{
		if ( ! user() ) {
			std::cerr << "ERROR: Unspecified option -" << id() << " is required" << std::endl;
			std::exit( EXIT_FAILURE );
		}
	}


public: // Properties


	/// @brief Key
	inline
	Key const &
	key() const override
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		return key_;
	}


	/// @brief ID
	inline
	std::string const &
	id() const override
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		return key_.id();
	}


	/// @brief Identifier
	inline
	std::string const &
	identifier() const override
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		return key_.identifier();
	}


	/// @brief Code
	inline
	std::string const &
	code() const override
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		return key_.code();
	}


	/// @brief Name
	inline
	std::string const &
	name() const override
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		return key_.id();
	}


	/// @brief Description
	inline
	std::string const &
	description() const override
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		return description_;
	}

	/// @brief Short Description
	inline
	std::string const &
	short_description() const override
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		return short_description_;
	}

	inline void short_description(std::string const & sd)
	{
#ifdef MULTI_THREADED
		utility::thread::WriteLockGuard writelock( mutex_ );
#endif
		short_description_ = sd;
	}


	/// @brief Legal or inactive default value and size?
	inline
	bool
	legal_default() const override
	{
		return ( ( legal_default_value() ) && ( default_size_ok() ) );
	}


	/// @brief Legal or inactive default value?
	inline
	bool
	legal_default_value() const override
	{
		return ( ( default_inactive() ) || ( unconstrained() ) || ( default_is_legal() ) || ( default_obeys_bounds() ) );
	}


	/// @brief Legal default value size?
	inline
	bool
	legal_default_size() const override
	{
		return default_size_ok();
	}


	/// @brief Legal value and size?
	inline
	bool
	legal() const override
	{
		return ( ( legal_value() ) && ( value_size_ok() ) );
	}


	/// @brief Legal value?
	inline
	bool
	legal_value() const override
	{
		return ( ( !active() ) || ( unconstrained() ) || ( value_is_legal() ) || ( value_obeys_bounds() ) );
	}


	/// @brief Legal value size?
	inline
	bool
	legal_size() const override
	{
		return value_size_ok();
	}

private:

	/// @brief Is the given value legal?
	/// @details If already_locked is true, there is already a write-lock on mutex_, so additional
	/// locking should not be done.
	inline
	bool
	legal_value( Value const & value_a, bool const already_locked ) const
	{
		return ( ( unconstrained(already_locked) ) || ( value_is_legal( value_a, already_locked ) ) || ( value_obeys_bounds( value_a, already_locked ) ) );
	}

public:

	/// @brief Overload of legal_value.  Assumes mutex_ is NOT locked, so that functions that this will call will obtain locks.
	inline bool legal_value( Value const & value_a ) const { return legal_value( value_a, false ); }


	/// @brief Has a default?
	inline
	bool
	has_default() const override
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		return ( default_state_ == DEFAULT );
	}


	/// @brief Default active?
	inline
	bool
	default_active() const override
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		return ( default_state_ == DEFAULT );
	}


	/// @brief Default inactive?
	inline
	bool
	default_inactive() const override
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		return ( default_state_ == INACTIVE );
	}


	/// @brief Active?  That is, the option has some value, either the default one or specified on the command line.
	inline
	bool
	active() const override
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		return ( state_ != INACTIVE );
	}


	/// @brief User-specified?  That is, the option value was specified on the command line.
	/// You should probably use active() instead in almost all cases!
	inline
	bool
	user() const override
	{
		been_accessed();
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		return ( state_ == USER );
	}


	/// @brief Can another value be added and stay within any size constraints?
	inline
	bool
	can_hold_another() const override
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		Size const s( size(true) );
		if ( ( n_ > 0 ) && ( s >= n_ ) ) return false;
		if ( ( n_upper_ > 0 ) && ( s >= n_upper_ ) ) return false;
		return true;
	}


	/// @brief Default size (number of default values)
	inline
	Size
	default_size() const override
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		return ( default_state_ == INACTIVE ? 0u : default_value_.size() );
	}


	/// @brief Number of default values (default size)
	inline
	Size
	n_default_value() const override
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		return ( default_state_ == INACTIVE ? 0u : default_value_.size() );
	}

private:

	/// @brief Size (number of values)
	/// @details Allows specification of whether mutex_ is already locked.
	inline
	Size
	size(
#ifdef MULTI_THREADED
		bool const already_locked
#else
		bool const
#endif
	) const {
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_, already_locked );
#endif
		return ( state_ == INACTIVE ? 0u : value_.size() );
	}

public:

	/// @brief Size( number of values).
	/// @details Assmues mutex_ is unlocked.
	inline Size size() const override { return size( false ); }


	/// @brief Number of values (size)
	inline
	Size
	n_value() const override
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		return ( state_ == INACTIVE ? 0u : value_.size() );
	}

	/// @brief Has Any Character of a std::string?
	/// non ambiguous vesrion for Python binding
	inline bool has_any_of_characters(std::string const & str1, std::string const & s ) const
	{
		size_type const s_len( s.length() );
		for ( char i : str1 ) {
			for ( size_type j = 0; j < s_len; ++j ) {
				if ( i == s[ j ] ) return true;
			}
		}
		return false; // No matches
	}

private:

	/// @brief Legal value string representation.
	/// @details Allows specification of whether mutex_ is locked.
	inline
	std::string
	legal_string(
#ifdef MULTI_THREADED
		bool const already_locked
#else
		bool const
#endif
	) const {
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_, already_locked );
#endif
		if ( ( legal_.empty() ) && ( lower_.inactive() ) && ( upper_.inactive() ) ) {
			return std::string();
		} else {
			std::ostringstream stream;
			stream_setup( stream );
			if ( ! legal_.empty() ) {
				//using ObjexxFCL::has_any_of;
				stream << '{';
				for ( typename LegalValues::const_iterator b = legal_.begin(), i = b, e = legal_.end(); i != e; ++i ) {
					if ( i != b ) stream << ',';
					std::string const s( value_string_of( *i ) );
					if ( has_any_of_characters( s, " ,\"" ) ) { // Quote-wrap the string
						stream << '"' << s << '"';
					} else {
						stream << s;
					}
				}
				stream << '}';
				if ( ( lower_.active() ) || ( upper_.active() ) ) stream << " or ";
			}
			if ( ( lower_.active() ) || ( upper_.active() ) ) {
				if ( lower_.active() ) {
					stream << ( lower_.strict() ? '(' : '[' ) << lower_();
				} else {
					stream << '<';
				}
				stream << '-';
				if ( upper_.active() ) {
					stream << upper_() << ( upper_.strict() ? ')' : ']' );
				} else {
					stream << '>';
				}
			}
			return stream.str();
		}
	}

public:

	/// @brief Legal value string representation.
	/// @details Assumes that mutex_ is unlocked.
	inline
	std::string
	legal_string() const override {
		return legal_string(false);
	}


	/// @brief Size constraint string representation
	inline
	std::string
	size_constraint_string() const override
	{
		std::ostringstream stream;
		if ( n_ > 0 ) {
			stream << '|' << n_ << '|';
		} else if ( ( n_lower_ > 0 ) || ( n_upper_ > 0 ) ) {
			stream << '|';
			if ( n_lower_ > 0 ) {
				stream << n_lower_;
			} else {
				stream << '<';
			}
			stream << '-';
			if ( n_upper_ > 0 ) {
				stream << n_upper_;
			} else {
				stream << '>';
			}
			stream << '|';
		}
		return stream.str();
	}

private:

	/// @brief Default value string representation
	inline
	std::string
	default_string(
#ifdef MULTI_THREADED
		bool const already_locked
#else
		bool const
#endif
	) const {
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard lock( mutex_, already_locked );
#endif
		if ( ( default_state_ == DEFAULT ) && ( ! default_value_.empty() ) ) {
			//using ObjexxFCL::has_any_of;
			std::ostringstream stream;
			stream_setup( stream );
			stream << '[';
			for ( ConstIterator i = default_value_.begin(), e = default_value_.end(); i != e; ++i ) {
				stream << ' ';
				std::string const s( value_string_of( *i ) );
				if ( has_any_of_characters( s, " \"" ) ) { // Quote-wrap the string
					stream << '"' << s << '"';
				} else {
					stream << s;
				}
			}
			stream << " ]";
			return stream.str();
		} else { // Default inactive or empty
			return std::string();
		}
	}

public:

	/// @brief Default value string representation
	inline
	std::string
	default_string() const override {
		return default_string(false);
	}

private:

	/// @brief Same as default_string, but without the "[" and "]"s wrapping the value list
	inline
	std::string
	raw_default_string(
#ifdef MULTI_THREADED
		bool const already_locked
#else
		bool const
#endif
	) const {
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard lock( mutex_, already_locked );
#endif
		if ( ( default_state_ == DEFAULT ) && ( ! default_value_.empty() ) ) {
			//using ObjexxFCL::has_any_of;
			std::ostringstream stream;
			stream_setup( stream );
			for ( ConstIterator i = default_value_.begin(), e = default_value_.end(); i != e; ++i ) {
				if ( i != default_value_.begin() ) { stream << ' '; }
				std::string const s( value_string_of( *i ) );
				if ( has_any_of_characters( s, " \"" ) ) { // Quote-wrap the string
					stream << '"' << s << '"';
				} else {
					stream << s;
				}
			}
			return stream.str();
		} else { // Default inactive or empty
			return std::string();
		}
	}

public:

	/// @brief Same as default_string, but without the "[" and "]"s wrapping the value list
	inline
	std::string
	raw_default_string() const override {
		return raw_default_string(false);
	}

	/// @brief Value string representation
	inline
	std::string
	value_string() const override
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		if ( ( state_ != INACTIVE ) && ( ! value_.empty() ) ) {
			//using ObjexxFCL::has_any_of;
			std::ostringstream stream;
			stream_setup( stream );
			for ( ConstIterator i = value_.begin(), e = value_.end(); i != e; ++i ) {
				if ( i != value_.begin() ) stream << ' ';
				std::string const s( value_string_of( *i ) );
				if ( has_any_of_characters( s, " \"" ) ) { // Quote-wrap the string
					stream << '"' << s << '"';
				} else {
					stream << s;
				}
			}
			return stream.str();
		} else { // Value inactive or empty
			return std::string();
		}
	}

	inline
	std::string
	raw_value_string() const override {
		return value_string();
	}


	/// @brief =Value string representation
	inline
	std::string
	equals_string() const override
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		if ( ( state_ != INACTIVE ) && ( ! value_.empty() ) ) {
			return '=' + value_string();
		} else { // Value inactive or empty
			return "=";
		}
	}


	/// @brief Lower bound
	inline
	LegalBound const &
	lower() const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		return lower_;
	}


	/// @brief Upper bound
	inline
	LegalBound const &
	upper()
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		return upper_;
	}


	/// @brief Fixed number of values required?
	inline
	bool
	fixed_size() const override
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		return ( n_ > 0 );
	}


	/// @brief Fixed number of values required (zero if none)
	inline
	Size
	n() const override
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		return n_;
	}


	/// @brief Lower number of values allowed (zero if none)
	inline
	Size
	n_lower() const override
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		return n_lower_;
	}


	/// @brief Upper number of values allowed (zero if none)
	inline
	Size
	n_upper() const override
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		return n_upper_;
	}


	/// @brief Default value
	inline
	Values const &
	default_value() const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		if ( default_state_ == INACTIVE ) default_inactive_error();
		return default_value_;
	}


	/// @brief Default value at a given index
	inline
	Value const &
	default_value( Size const i ) const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		if ( default_state_ == INACTIVE ) default_inactive_error();
		return default_value_[ i ];
	}


	/// @brief Value
	inline
	Values const &
	value() const
	{
		been_accessed();
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		if ( state_ == INACTIVE ) inactive_error();
		return value_;
	}


	/// @brief Value
	inline
	Values const &
	operator ()() const
	{
		been_accessed();
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		if ( state_ == INACTIVE ) inactive_error();
		return value_;
	}


	/// @brief Value or passed default if inactive
	inline
	Values // Have to return by value: Not efficient for many or large Value types
	value_or( Values const & value_a ) const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		if ( state_ != INACTIVE ) { // Return active value
			been_accessed(); //Threadsafe
			return value_;
		} else { // Return passed value
			return value_a;
		}
	}


	/// @brief Value or passed default if not user-specified
	inline
	Values // Have to return by value: Not efficient for many or large Value types
	user_or( Values const & value_a ) const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		if ( state_ == USER ) { // Return user-specified value
			been_accessed(); //Threadsafe
			return value_;
		} else { // Return passed value
			return value_a;
		}
	}


	/// @brief Value at a given index
	inline
	Value const &
	value( Size const i ) const
	{
		been_accessed(); //Threadsafe
#ifdef MULTI_THREADED
		utility::thread::WriteLockGuard writelock( mutex_ );
#endif
		if ( state_ == INACTIVE ) inactive_error();
		return value_[ i ];
	}


	/// @brief Does the VectorOption contain the given value?
	inline
	bool
	has_value( Value const & value ) {
		been_accessed();
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		if ( state_ == INACTIVE ) inactive_error();
		return value_.has_value( value );
	}

	/// @brief Value at a given index
	inline
	Value const &
	operator ()( Size const i ) const
	{
		been_accessed();
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		if ( state_ == INACTIVE ) inactive_error();
		return value_[ i ];
	}


	/// @brief Value at a given index
	inline
	Value const &
	operator []( Size const i ) const
	{
		been_accessed();
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		if ( state_ == INACTIVE ) inactive_error();
		return value_[ i ];
	}


	/// @brief Value at a given index or passed default if inactive
	inline
	Value // Have to return by value: Not efficient for large Value types
	value_or( Size const i, Value const & value_a ) const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		if ( state_ != INACTIVE ) { // Return active value
			return value_[ i ];
		} else { // Return passed value
			return value_a;
		}
	}


	/// @brief Value at a given index or passed default if not user-specified
	inline
	Value // Have to return by value: Not efficient for large Value types
	user_or( Size const i, Value const & value_a ) const
	{
		been_accessed();
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		if ( state_ == USER ) { // Return user-specified value
			return value_[ i ];
		} else { // Return passed value
			return value_a;
		}
	}


protected: // Methods


	/// @brief Value of a string
	virtual
	Value
	value_of( std::string const & value_str ) const = 0;

	/// @brief Value of a string
	virtual
	Values
	values_of( std::string const & value_str ) const
	{
		Values vs;
		vs.push_back( value_of( value_str ) );
		return vs;
	}


	/// @brief String representation of a given value
	inline
	virtual
	std::string
	value_string_of( Value const & v ) const
	{
		std::ostringstream stream;
		stream_setup( stream );
		stream << std::boolalpha << v;
		return stream.str();
	}


	/// @brief Error handler for using inactive option value
	/// @note  Not handled with assert because unspecified option can be a user error
	inline
	void
	default_inactive_error() const
	{
		debug_assert( default_state_ == INACTIVE ); // Or else why are we here
		std::cerr << "ERROR: Inactive default value of option accessed: -" << key_.id() << std::endl;
		std::exit( EXIT_FAILURE );
	}


	/// @brief Error handler for using inactive option value
	/// @note  Not handled with assert because unspecified option can be a user error
	inline
	void
	inactive_error() const
	{
		debug_assert( state_ == INACTIVE ); // Or else why are we here
		std::cerr << "ERROR: Value of inactive option accessed: -" << key_.id() << std::endl;
		std::exit( EXIT_FAILURE );
	}


	/// @brief Setup stream state for the Option value type
	virtual
	void
	stream_setup( std::ostream & ) const
	{}


private: // Properties -- most are protected; first one is private.


	/// @brief Value is unconstrained?
	/// @details Pass true if mutex_ is already locked.
	inline
	bool
	unconstrained(
#ifdef MULTI_THREADED
		bool const already_locked
#else
		bool const
#endif
	) const {
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_, already_locked );
#endif
		return ( ( legal_.empty() ) && ( lower_.inactive() ) && ( upper_.inactive() ) );
	}

protected:

	/// @brief Value is unconstrained?
	/// @details Assumes mutex_ is unlocked.
	inline bool unconstrained() const { return unconstrained( false ); }

private:

	/// @brief Default value is a specified legal value?
	/// @details Pass true if mutex_ is already locked.
	inline
	bool
	default_is_legal(
#ifdef MULTI_THREADED
		bool const already_locked
#else
		bool const
#endif
	) const {

#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_, already_locked );
#endif
		if ( ( default_state_ == INACTIVE ) || ( default_value_.empty() ) ) {
			return false;
		} else {
			for ( ConstIterator i = default_value_.begin(), e = default_value_.end(); i != e; ++i ) {
				if ( legal_.find( *i ) == legal_.end() ) return false;
			}
			return true;
		}
	}

protected:

	/// @brief Default value is a specified legal value?
	/// @details Assumes mutex_ is unlocked.
	inline
	bool
	default_is_legal() const {
		return default_is_legal(false);
	}

	/// @brief Value is a specified legal value?
	inline
	bool
	value_is_legal() const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		if ( ( state_ == INACTIVE ) || ( value_.empty() ) ) {
			return false;
		} else {
			for ( ConstIterator i = value_.begin(), e = value_.end(); i != e; ++i ) {
				if ( legal_.find( *i ) == legal_.end() ) return false;
			}
			return true;
		}
	}

private:

	/// @brief Value is legal?
	/// @details Pass true if mutex_ is already locked, false otherwise.
	inline
	bool
	value_is_legal(
		Value const & value_a,
#ifdef MULTI_THREADED
		bool const already_locked
#else
		bool const
#endif
	) const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_, already_locked );
#endif
		return ( legal_.find( value_a ) != legal_.end() );
	}

protected:

	/// @brief Value is legal?
	/// @details Assumes mutex_ is not locked.
	inline
	bool
	value_is_legal( Value const & value_a ) const {
		return value_is_legal( value_a, false );
	}


	/// @brief Default value obeys specified bounds?
	inline
	bool
	default_obeys_bounds() const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		if ( default_state_ == INACTIVE ) {
			return false;
		} else {
			if ( lower_.active() ) {
				for ( ConstIterator i = default_value_.begin(), e = default_value_.end(); i != e; ++i ) {
					Value const & v( *i );
					if ( lower_.strict() ) {
						if ( lower_() >= v ) return false;
					} else {
						if ( lower_() > v ) return false;
					}
				}
				if ( ! upper_.active() ) return true;
			}
			if ( upper_.active() ) {
				for ( ConstIterator i = default_value_.begin(), e = default_value_.end(); i != e; ++i ) {
					Value const & v( *i );
					if ( upper_.strict() ) {
						if ( v >= upper_() ) return false;
					} else {
						if ( v > upper_() ) return false;
					}
				}
				return true;
			}
			return false; // No bounds specified
		}
	}


	/// @brief Value obeys specified bounds?
	/// @details Assumes mutex_ is unlocked.
	inline
	bool
	value_obeys_bounds() const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		if ( state_ == INACTIVE ) {
			return false;
		} else {
			if ( lower_.active() ) {
				for ( ConstIterator i = value_.begin(), e = value_.end(); i != e; ++i ) {
					Value const & v( *i );
					if ( lower_.strict() ) {
						if ( lower_() >= v ) return false;
					} else {
						if ( lower_() > v ) return false;
					}
				}
				if ( ! upper_.active() ) return true;
			}
			if ( upper_.active() ) {
				for ( ConstIterator i = value_.begin(), e = value_.end(); i != e; ++i ) {
					Value const & v( *i );
					if ( upper_.strict() ) {
						if ( v >= upper_() ) return false;
					} else {
						if ( v > upper_() ) return false;
					}
				}
				return true;
			}
			return false; // No bounds specified
		}
	}

private:

	/// @brief Given value obeys specified bounds?
	/// @details Allows specification of whether mutex_ is already locked.
	inline
	bool
	value_obeys_bounds(
		Value const & value_a,
#ifdef MULTI_THREADED
		bool const already_locked
#else
		bool const
#endif
	) const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_, already_locked );
#endif
		if ( lower_.active() ) {
			if ( lower_.strict() ) {
				if ( lower_() >= value_a ) return false;
			} else {
				if ( lower_() > value_a ) return false;
			}
			if ( ! upper_.active() ) return true;
		}
		if ( upper_.active() ) {
			if ( upper_.strict() ) {
				return ( value_a < upper_() );
			} else {
				return ( value_a <= upper_() );
			}
		}
		return false; // No bounds specified
	}

protected:

	/// @brief Given value obeys specified bounds?
	/// @details Assumes mutex_ is unlocked.
	inline
	bool
	value_obeys_bounds( Value const & value_a ) const {
		return value_obeys_bounds( value_a, false );
	}


	/// @brief Default value size is OK?
	inline
	bool
	default_size_ok() const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		if ( default_state_ == INACTIVE ) {
			return true;
		} else {
			Size const nd( default_value_.size() );
			if ( ( n_ > 0 ) && ( nd != n_ ) ) return false;
			if ( ( n_lower_ > 0 ) && ( nd < n_lower_ ) ) return false;
			if ( ( n_upper_ > 0 ) && ( nd > n_upper_ ) ) return false;
			return true;
		}
	}


	/// @brief Value size is OK?
	inline
	bool
	value_size_ok() const
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		if ( state_ == INACTIVE ) {
			return true;
		} else {
			Size const nv( value_.size() );
			if ( ( n_ > 0 ) && ( nv != n_ ) ) return false;
			if ( ( n_lower_ > 0 ) && ( nv < n_lower_ ) ) return false;
			if ( ( n_upper_ > 0 ) && ( nv > n_upper_ ) ) return false;
			return true;
		}
	}


private: // Fields

#ifdef MULTI_THREADED
	/// @brief Mutex for controlling access.
	mutable utility::thread::ReadWriteMutex mutex_;
#endif

	/// @brief Key
	Key key_;

	/// @brief Description
	std::string description_;

	/// @brief Short Description
	std::string short_description_;

	/// @brief Legal values
	LegalValues legal_;

	/// @brief Bound values
	LegalBound lower_, upper_;

	/// @brief Size requirements
	Size n_, n_lower_, n_upper_;

	/// @brief Default state
	State default_state_;

	/// @brief Value
	Values default_value_;

	/// @brief State
	State state_;

	/// @brief Value
	Values value_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const {
		cereal::base_class< utility::options::VectorOption >( this );
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock( mutex_ );
#endif
		arc(key_);
		arc(description_);
		arc(short_description_);
		arc(legal_);
		arc(lower_);
		arc(upper_);
		arc(n_, n_lower_, n_upper_);
		arc(default_state_);
		arc(default_value_);
		arc(state_);
		arc(value_);
	}

	template< class Archive > void load( Archive & arc ) {
		cereal::base_class< utility::options::VectorOption >( this );
#ifdef MULTI_THREADED
		utility::thread::WriteLockGuard writelock( mutex_ );
#endif
		arc(key_);
		arc(description_);
		arc(short_description_);
		arc(legal_);
		arc(lower_);
		arc(upper_);
		arc(n_, n_lower_, n_upper_);
		arc(default_state_);
		arc(default_value_);
		arc(state_);
		arc(value_);
	}
#endif // SERIALIZATION


}; // VectorOption_T_


} // namespace options
} // namespace utility


#endif // INCLUDED_utility_options_VectorOption_T__HH
