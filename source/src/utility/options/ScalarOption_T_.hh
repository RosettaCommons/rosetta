// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/options/ScalarOption_T_.hh
/// @brief  Program scalar-valued option abstract base class
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author Modified by Sergey Lyskov


#ifndef INCLUDED_utility_options_ScalarOption_T__hh
#define INCLUDED_utility_options_ScalarOption_T__hh


// Unit headers
#include <utility/options/ScalarOption_T_.fwd.hh>
#include <utility/options/mpi_stderr.hh>

// Package headers
#include <utility/options/ScalarOption.hh>
#include <utility/Bound.hh>
// Auto-header: duplicate removed #include <utility/options/mpi_stderr.hh>
// ObjexxFCL headers
#include <ObjexxFCL/string.functions.hh>

// C++ headers
#include <utility/assert.hh>
#include <cstdlib>
#include <iostream>
#include <set>
#include <sstream>


namespace utility {
namespace options {


/// @brief Program scalar-valued option abstract base class
template< typename K, typename T >
class ScalarOption_T_ :
	public ScalarOption
{


private: // Types


	typedef  ScalarOption  Super;

	typedef  std::set< T >  LegalValues;
	typedef  Bound< T >  LegalBound;


public: // Types


	// STL/boost style
	typedef  K  key_type;
	typedef  T  value_type;

	// Project style
	typedef  K  Key;
	typedef  T  Value;


protected: // Creation


	/// @brief Default constructor
	inline
	ScalarOption_T_() :
		default_state_( INACTIVE ),
		default_value_( Value() ),
		state_( INACTIVE ),
		value_( Value() )
	{}


	/// @brief Copy constructor
	inline
	ScalarOption_T_( ScalarOption_T_ const & option ) :
		Super( option ),
		key_( option.key_ ),
		description_( option.description_ ),
		short_description_( option.short_description_ ),
		legal_( option.legal_ ),
		lower_( option.lower_ ),
		upper_( option.upper_ ),
		default_state_( option.default_state_ ),
		default_value_( option.default_value_ ),
		state_( option.state_ ),
		value_( option.value_ )
	{}


	/// @brief Key + description constructor
	inline
	ScalarOption_T_(
		Key const & key_a,
		std::string const & description_a
	) :
		key_( key_a ),
		description_( description_a ),
		short_description_( description_a ),
		default_state_( INACTIVE ),
		default_value_( Value() ),
		state_( INACTIVE ),
		value_( Value() )
	{}


public: // Creation


	/// @brief Clone this
	virtual
	ScalarOption_T_ *
	clone() const = 0;


	/// @brief Destructor
	inline
	virtual
	~ScalarOption_T_()
	{}


protected: // Assignment


	/// @brief Copy assignment
	inline
	ScalarOption_T_ &
	operator =( ScalarOption_T_ const & option )
	{
		Option::operator=(option);

		if ( this != &option ) {
			key_ = option.key_;
			description_ = option.description_;
			short_description_ = option.short_description_;
			legal_ = option.legal_;
			lower_ = option.lower_;
			upper_ = option.upper_;
			default_state_ = option.default_state_;
			default_value_ = option.default_value_;
			state_ = option.state_;
			value_ = option.value_;
		}
		return *this;
	}


public: // Conversion


	/// @brief Value conversion
	inline
	operator Value const &() const
	{
		been_accessed();
		if ( state_ == INACTIVE ) inactive_error();
		return value_;
	}

	/* Commenting this out for clear separation of read/write operators
	/// @brief Value conversion
	inline
	operator Value &()
	{
		//been_accessed();
		if ( state_ == INACTIVE ) inactive_error();
		return value_;
	}
	*/

public: // Methods


	/// @brief Activate
	inline
	ScalarOption_T_ &
	activate()
	{
		state_ = USER;
		return *this;
	}


	/// @brief Deactivate
	inline
	ScalarOption_T_ &
	deactivate()
	{
		state_ = INACTIVE;
		return *this;
	}


	/// @brief Set to default value, if any
	inline
	ScalarOption_T_ &
	to_default()
	{
		if ( default_state_ == DEFAULT ) {
			state_ = DEFAULT;
			value_ = default_value_;
		}
		return *this;
	}


	/// @brief Clear
	inline
	ScalarOption_T_ &
	clear()
	{
		default_state_ = INACTIVE;
		default_value_ = Value();
		state_ = INACTIVE;
		value_ = Value();
		return *this;
	}


	/// @brief Add a legal value
	inline
	ScalarOption_T_ &
	legal( Value const & value_a )
	{
		legal_.insert( value_a );
		return *this;
	}

	/// @brief Set a short description
	inline
	ScalarOption_T_ &
	shortd( std::string const & s)
	{
		short_description_ = s;
		return *this;
	}


	/// @brief Set a lower bound
	inline
	ScalarOption_T_ &
	lower( Value const & value_a )
	{
		lower_( value_a );
		return *this;
	}


	/// @brief Set a strict lower bound
	inline
	ScalarOption_T_ &
	strict_lower( Value const & value_a )
	{
		lower_( value_a, true );
		return *this;
	}


	/// @brief Set an upper bound
	inline
	ScalarOption_T_ &
	upper( Value const & value_a )
	{
		upper_( value_a );
		return *this;
	}


	/// @brief Set a strict upper bound
	inline
	ScalarOption_T_ &
	strict_upper( Value const & value_a )
	{
		upper_( value_a, true );
		return *this;
	}


	/// @brief Default value assignment
	inline
	virtual
	ScalarOption_T_ &
	default_value( Value const & value_a )
	{
		default_state_ = DEFAULT;
		default_value_ = value_a;
		if ( ( state_ == INACTIVE ) || ( state_ == DEFAULT ) ) {
			state_ = DEFAULT;
			value_ = value_a;
		}
		legal_default_check();
		return *this;
	}


	/// @brief Default value assignment
	inline
	virtual
	ScalarOption_T_ &
	def( Value const & value_a )
	{
		default_state_ = DEFAULT;
		default_value_ = value_a;
		if ( ( state_ == INACTIVE ) || ( state_ == DEFAULT ) ) {
			state_ = DEFAULT;
			value_ = value_a;
		}
		legal_default_check();
		return *this;
	}


	/// @brief Value assignment from a command line string
	inline
	ScalarOption_T_ &
	cl_value( std::string const & value_str )
	{
		Value const old_value( value_ );
		bool const check_override( state_ == USER );
		value( value_of( ObjexxFCL::stripped( value_str, "\"'" ) ) );
		if ( ( check_override ) && ( value_ != old_value ) ) {
#ifndef BOINC
			mpi_safe_std_err("WARNING: Override of option -"+id()+" sets a different value");
#endif
		}
		if ( ! legal_value( value_ ) ) {
			mpi_safe_std_err( "ERROR: Illegal value specified for option -" +id()+ " : " + value_str );
			std::exit( EXIT_FAILURE );
		}
		return *this;
	}


	/// @brief Value assignment
	inline
	virtual
	ScalarOption_T_ &
	value( Value const & value_a )
	{
		state_ = USER;
		value_ = value_a;
		legal_check();
		return *this;
	}


	/// @brief Value assignment
	inline
	virtual
	ScalarOption_T_ &
	operator ()( Value const & value_a )
	{
		state_ = USER;
		value_ = value_a;
		legal_check();
		return *this;
	}


	/// @brief Default to another option's value
	inline
	ScalarOption_T_ &
	default_to( ScalarOption_T_ const & option )
	{
		if ( ( state_ == INACTIVE ) || ( state_ == DEFAULT ) ) {
			if ( option.active() ) default_value( option.value() );
		}
		return *this;
	}


	/// @brief Legal specifications check: Report and return error state
	inline
	bool
	legal_specs_report() const
	{
		return ( ( legal_limits_report() ) && ( legal_default_report() ) );
	}


	/// @brief Legal value limits check: Report and return error state
	inline
	bool
	legal_limits_report() const
	{
		bool error( false );
		if ( ( lower_.active() ) && ( upper_.active() ) ) {
			if ( ( lower_.strict() ) || ( upper_.strict() ) ) {
				if ( lower_() >= upper_() ) error = true;
			} else {
				if ( lower_() > upper_() ) error = true;
			}
			if ( error ) {
				std::cerr << "ERROR: Inconsistent lower and upper limits in option -" << id()
				 << " : " << legal_string() << std::endl;
			}
		}
		return ( ! error );
	}


	/// @brief Legal size limits check: Report and return error state
	inline
	bool
	legal_size_report() const
	{
		return true; // Scalar options have size 1
	}


	/// @brief Legal default value check: Report and return error state
	inline
	bool
	legal_default_report() const
	{
		bool error( false );
		if ( ! legal_default() ) {
			std::cerr << "ERROR: Illegal default value in option -" << id()
			 << " : " << default_string() << std::endl;
			error = true;
		}
		return ( ! error );
	}


	/// @brief Legal default value check
	inline
	void
	legal_default_check() const
	{
		if ( ! legal_default() ) {
			std::cerr << "ERROR: Illegal default value in option -" << id()
			 << " : " << default_string() << std::endl;
			std::exit( EXIT_FAILURE );
		}
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
	legal_report() const
	{
		bool error( false );
		if ( ! legal() ) {
			mpi_safe_std_err( "ERROR: Illegal value specified for option -" +id()+ " : "+ value_string());
			error = true;
		}
		return ( ! error );
	}


	/// @brief Legal value check
	inline
	void
	legal_check() const
	{
		if ( ! legal() ) {
			mpi_safe_std_err("ERROR: Illegal value specified for option -"+id()+" : "+value_string());

			std::exit( EXIT_FAILURE );
		}
	}


	/// @brief Legal value check
	inline
	void
	legal_check( Value const & value_a ) const
	{
		if ( ! legal_value( value_a ) ) {
			mpi_safe_std_err("ERROR: Illegal value specified for option -"+id()+" : "+value_string_string_of(value_a));
			std::exit( EXIT_FAILURE );
		}
	}


	/// @brief Required specified option check: Report and return error state
	inline
	bool
	specified_report() const
	{
		bool error( false );
		if ( ! user() ) {
			mpi_safe_std_err("ERROR: Unspecified option -"+id()+ " is required");
			error = true;
		}
		return ( ! error );
	}


	/// @brief Required specified option check
	inline
	void
	specified_check() const
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
	key() const
	{
		return key_;
	}


	/// @brief ID
	inline
	std::string const &
	id() const
	{
		return key_.id();
	}


	/// @brief Identifier
	inline
	std::string const &
	identifier() const
	{
		return key_.identifier();
	}


	/// @brief Code
	inline
	std::string const &
	code() const
	{
		return key_.code();
	}


	/// @brief Name
	inline
	std::string const &
	name() const
	{
		return key_.id();
	}


	/// @brief Description
	inline
	std::string const &
	description() const
	{
		return description_;
	}

	/// @brief Short Description
	inline
	std::string const &
	short_description() const
	{
		return short_description_;
	}

	inline void	short_description(std::string const & sd)
	{
		short_description_ = sd;
	}

	inline void	description(std::string const & sd )
	{
		description_ = sd;
	}


	/// @brief Legal or inactive default value?
	inline
	bool
	legal_default() const
	{
		return ( ( default_inactive() ) || ( unconstrained() ) || ( default_is_legal() ) || ( default_obeys_bounds() ) );
	}


	/// @brief Legal value?
	inline
	bool
	legal() const
	{
		return ( ( !active() ) || ( unconstrained() ) || ( value_is_legal() ) || ( value_obeys_bounds() ) );
	}


	/// @brief Is the given value legal?
	inline
	bool
	legal_value( Value const & value_a ) const
	{
		return ( ( unconstrained() ) || ( value_is_legal( value_a ) ) || ( value_obeys_bounds( value_a ) ) );
	}


	/// @brief Has a default?
	inline
	bool
	has_default() const
	{
		return ( default_state_ == DEFAULT );
	}


	/// @brief Default active?
	inline
	bool
	default_active() const
	{
		return ( default_state_ == DEFAULT );
	}


	/// @brief Default inactive?
	inline
	bool
	default_inactive() const
	{
		return ( default_state_ == INACTIVE );
	}


	/// @brief Active?  That is, the option has some value, either the default one or specified on the command line.
	inline
	bool
	active() const
	{
		return ( state_ != INACTIVE );
	}


	/// @brief User-specified?  That is, the option value was specified on the command line.
	/// You should probably use active() instead in almost all cases!
	inline
	bool
	user() const
	{
		been_accessed();
		return ( state_ == USER );
	}


	/// @brief Can another value be added and stay within any size constraints?
	inline
	bool
	can_hold_another() const
	{
		return true; // Always hold another since new value replaces old in scalar
	}


	/// @brief Default size (number of default values)
	inline
	Size
	default_size() const
	{
		return ( default_state_ == INACTIVE ? 0u : 1u );
	}


	/// @brief Number of default values (default size)
	inline
	Size
	n_default_value() const
	{
		return ( default_state_ == INACTIVE ? 0u : 1u );
	}


	/// @brief Size (number of values)
	inline
	Size
	size() const
	{
		return ( state_ == INACTIVE ? 0u : 1u );
	}


	/// @brief Number of values (size)
	inline
	Size
	n_value() const
	{
		return ( state_ == INACTIVE ? 0u : 1u );
	}

	/// @brief Has Any Character of a std::string?
	/// non ambiguous vesrion for Python binding
	inline bool has_any_of_characters(std::string const str_, std::string const s ) const
	// AMW: cppcheck wants references.
	//I believe we cannot use references here due to the Python binding aspect
	{
		size_type const s_len( s.length() );
		for ( size_type i = 0; i < str_.size(); ++i ) {
			for ( size_type j = 0; j < s_len; ++j ) {
				if ( str_[ i ] == s[ j ] ) return true;
			}
		}
		return false; // No matches
	}


	/// @brief Legal value string representation
	inline
	std::string
	legal_string() const
	{
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


	/// @brief Size constraint string representation
	inline
	std::string
	size_constraint_string() const
	{
		return std::string(); // Scalar options have no size constraints
	}


	/// @brief Default value string representation
	inline
	std::string
	default_string() const
	{
		if ( default_state_ == DEFAULT ) {
			return '[' + value_string_of( default_value_ ) + ']';
		} else { // Default inactive
			return std::string();
		}
	}


	/// @brief Value string representation
	inline
	std::string
	value_string() const
	{
		if ( state_ != INACTIVE ) {
			return value_string_of( value_ );
		} else { // Value inactive
			return std::string();
		}
	}


	/// @brief =Value string representation
	inline
	std::string
	equals_string() const
	{
		if ( state_ != INACTIVE ) {
			return '=' + value_string();
		} else { // Value inactive
			return "=";
		}
	}


	/// @brief Lower bound
	inline
	LegalBound const &
	lower() const
	{
		return lower_;
	}


	/// @brief Upper bound
	inline
	LegalBound const &
	upper()
	{
		return upper_;
	}


	/// @brief Default value
	inline
	Value const &
	default_value() const
	{
		if ( default_state_ == INACTIVE ) default_inactive_error();
		return default_value_;
	}


	/// @brief Value
	inline
	Value const &
	value() const
	{
		been_accessed();
		if ( state_ == INACTIVE ) inactive_error();
		return value_;
	}


	/// @brief Value
	inline
	Value const &
	operator ()() const
	{
		been_accessed();
		if ( state_ == INACTIVE ) inactive_error();
		return value_;
	}


	/// @brief Value or passed default if inactive
	inline
	Value // Have to return by value: Not efficient for large Value types
	value_or( Value const & value_a ) const
	{
		//been_accessed();
		if ( state_ != INACTIVE ) { // Return active value
			return value_;
		} else { // Return passed value
			return value_a;
		}
	}


	/// @brief Value or passed default if not user-specified
	inline
	Value // Have to return by value: Not efficient for large Value types
	user_or( Value const & value_a ) const
	{
		been_accessed();
		if ( state_ == USER ) { // Return user-specified value
			return value_;
		} else { // Return passed value
			return value_a;
		}
	}


protected: // Methods


	/// @brief Value of a string
	virtual
	Value
	value_of( std::string const & value_str ) const = 0;


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


protected: // Properties


	/// @brief Value is unconstrained?
	inline
	bool
	unconstrained() const
	{
		return ( ( legal_.empty() ) && ( lower_.inactive() ) && ( upper_.inactive() ) );
	}


	/// @brief Default value is a specified legal value?
	inline
	bool
	default_is_legal() const
	{
		if ( default_state_ == INACTIVE ) {
			return false;
		} else {
			return ( legal_.find( default_value_ ) != legal_.end() );
		}
	}


	/// @brief Value is a specified legal value?
	inline
	bool
	value_is_legal() const
	{
		if ( state_ == INACTIVE ) {
			return false;
		} else {
			return ( legal_.find( value_ ) != legal_.end() );
		}
	}


	/// @brief Value is legal?
	inline
	bool
	value_is_legal( Value const & value_a ) const
	{
		return ( legal_.find( value_a ) != legal_.end() );
	}


	/// @brief Default value obeys specified bounds?
	inline
	bool
	default_obeys_bounds() const
	{
		if ( default_state_ == INACTIVE ) {
			return false;
		} else {
			if ( lower_.active() ) {
				if ( lower_.strict() ) {
					if ( lower_() >= default_value_ ) return false;
				} else {
					if ( lower_() > default_value_ ) return false;
				}
				if ( ! upper_.active() ) return true;
			}
			if ( upper_.active() ) {
				if ( upper_.strict() ) {
					return ( default_value_ < upper_() );
				} else {
					return ( default_value_ <= upper_() );
				}
			}
			return false; // No bounds specified
		}
	}


	/// @brief Value obeys specified bounds?
	inline
	bool
	value_obeys_bounds() const
	{
		if ( state_ == INACTIVE ) {
			return false;
		} else {
			if ( lower_.active() ) {
				if ( lower_.strict() ) {
					if ( lower_() >= value_ ) return false;
				} else {
					if ( lower_() > value_ ) return false;
				}
				if ( ! upper_.active() ) return true;
			}
			if ( upper_.active() ) {
				if ( upper_.strict() ) {
					return ( value_ < upper_() );
				} else {
					return ( value_ <= upper_() );
				}
			}
			return false; // No bounds specified
		}
	}


	/// @brief Given value obeys specified bounds?
	inline
	bool
	value_obeys_bounds( Value const & value_a ) const
	{
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


private: // Fields


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

	/// @brief Default state
	State default_state_;

	/// @brief Value
	Value default_value_;

	/// @brief State
	State state_;

	/// @brief Value
	Value value_;


}; // ScalarOption_T_


} // namespace options
} // namespace utility


#endif // INCLUDED_utility_options_ScalarOption_T__HH
