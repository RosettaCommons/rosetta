// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief An object for reading/writing a simple xml-like language.
/// @author Paul Murphy

#ifndef INCLUDED_utility_tag_tag_HH
#define INCLUDED_utility_tag_tag_HH

#if defined(__clang__) || defined(__llvm__)
#if __clang_major__ < 4
#define OLDER_CLANG
#endif
#endif

// Unit headers
#include <utility/tag/Tag.fwd.hh>

// C++ Headers
#include <sstream> // DO NOT AUTOREMOVE - needed for templated functions
#include <map>
#include <string>
#include <sstream>

// Utility headers
#include <utility/tag/AutoBool.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/VirtualBase.hh>
#include <utility/vector0.hh>

// Boost headers
#include <boost/lexical_cast.hpp>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

#ifdef PYROSETTA
#include <platform/types.hh>
#endif

namespace utility {
namespace tag {

class Tag : public utility::VirtualBase, public utility::pointer::enable_shared_from_this< Tag >
{
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from VirtualBase
	~Tag() override;

	typedef std::map<std::string, std::string > options_t;
	typedef utility::vector0< TagCOP > tags_t;

public:
	Tag();

	// BAD: Tag( std::string const & tag_string );
	// It's tempting to have a direct constructor instead of doing a create/read two-step
	// The problem is that we rely on enable_shared_from_this to add subtags.
	// This only works if the Tag is in an OP already (e.g. it doesn't work in the constructor.)

	/// self pointers
	inline TagCOP get_self_ptr() const { return shared_from_this(); }
	inline TagOP get_self_ptr() { return shared_from_this(); }
	inline TagCAP get_self_weak_ptr() const { return TagCAP( shared_from_this() ); }
	inline TagAP get_self_weak_ptr() { return TagAP( shared_from_this() ); }

	void clear();
	size_t size() const;

	void setName(std::string const& name);
	std::string const& getName() const { return name_; }

	TagCAP const & getParent() const { return parentTag_; }

	void addTag( TagOP tag );
	utility::vector0< TagCOP > const & getTags() const;
	utility::vector0< TagCOP > const & getTags( std::string const& name ) const;
	TagCOP const & getTag( std::string const& name ) const;

	///@brief Does the Tag have a specific sub-tag (branch) within (such as a MoveMap, etc.)?
	bool
	hasTag( std::string const& name ) const;

	///@brief Does the Tag have a specific option?
	bool
	hasOption( std::string const& key ) const;

	template< class T >
	void
	setOption(std::string const& key, T const& value) {
		if ( mOptions_.find(key) != mOptions_.end() ) {
			//runtime_assert( false );
		}
		std::ostringstream out;
		out << value;
		mOptions_[key] = out.str();
		accessed_options_.erase(key);
	} // setOption

	/// @brief Set the 'accessed' annotation of this option, without bothering to get the value.
	/// Used to prevent no-longer relevant options from crashing the XML parsing.
	void
	setAccessed(std::string const& key) const {
		accessed_options_[key]= key;
	}

	/// @brief Retrieve an option from the Tag with the given key name, using the
	/// provided default value (t_default) if the option is not present in the tag.
	/// @throws Throws a utility::excn::EXCN_Msg_Exception if the boost::lexical_cast
	/// fails to convert the input type as requested.
	template< class T >
	T
	getOption(std::string const& key, T const& t_default) const {
		auto i = mOptions_.find(key);
		if ( i == mOptions_.end() ) {
			accessed_options_[key]= key;
			return t_default;
		}
		accessed_options_[key]= i->second;
		//  T t = t_default; // not used?
		try{
			return boost::lexical_cast<T>(i->second);
		} catch(boost::bad_lexical_cast &) {
			std::stringstream error_message;
			error_message << "getOption: key= " << key << " stream extraction failed! Tried to parse '" << i->second << "'\n";
			throw CREATE_EXCEPTION(utility::excn::Exception,  error_message.str() );
		}
		return t_default;
	}

	/// @brief Variant that will only be defined for Sizes.
	/// @note Older clang compilers have trouble with the general case being
	/// deleted and then specialized cases being defined, even though this is supposed to be supported
	/// by the cxx11 standard.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
	template< class T >
	T
#ifdef OLDER_CLANG
	getOption( std::string const &, int const ) const;
#else
	getOption( std::string const &, int const ) const = delete;
#endif

	/// @brief Variant for the case in which the developer has mistakenly provided a string literal
	/// instead of a value of type T.  For anything but a boolean or a string type, this produces a
	/// compilation error.
	/// @note Older clang compilers have trouble with the general case being
	/// deleted and then specialized cases being defined, even though this is supposed to be supported
	/// by the cxx11 standard.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
	template< class T >
	T
#ifdef OLDER_CLANG
	getOption(std::string const& key, char const * default_as_string_literal ) const {
		utility_exit_with_message( "Program error: the developer has erroneously provided the string literal \"" + std::string(default_as_string_literal) + "\" as the default for option key \"" + key + "\".  This is effectively a compilation error, but it is only detectable at runtime.  The temporary workaround is to provide a value for the " + key + " option.  Please also inform a developer." );
		return T(0);
	}
#else
	getOption(std::string const& key, char const * default_as_string_literal ) const = delete;
#endif

	/// @brief Retrieve an option from the Tag with the given key name.
	/// @throws Throws a utility::excn::EXCN_Msg_Exception if the an option
	/// with the given key is not present, or if the boost::lexical_cast fails
	/// to convert the input type as requested.
	template< class T >
	T
	getOption(std::string const& key) const {
		auto i = mOptions_.find(key);
		if ( i == mOptions_.end() ) {
			std::stringstream error_message;
			error_message << "Option " << key << " not found.\n";
			throw CREATE_EXCEPTION(utility::excn::Exception,  error_message.str() );
		}
		accessed_options_[key]= i->second;
		try{
			return boost::lexical_cast<T>(i->second);
		} catch(boost::bad_lexical_cast &) {
			std::stringstream error_message;
			error_message << "getOption: key= " << key << " stream extraction failed! Tried to parse '" << i->second << "'\n";
			throw CREATE_EXCEPTION(utility::excn::Exception, error_message.str() );
		}
		return T(); // appease compiler
	}

#ifdef PYROSETTA
	bool get_option_bool(std::string const& key) const;
	bool get_option_bool(std::string const& key, char const * default_as_string_literal) const;
	bool get_option_bool(std::string const& key, bool const& t_default) const;
 	int get_option_int(std::string const& key) const;
	int get_option_int(std::string const& key, int const& t_default) const;
	platform::Real get_option_real(std::string const& key) const;
	platform::Real get_option_real(std::string const& key, platform::Real const& t_default) const;
	platform::Size get_option_size(std::string const& key) const;
	platform::Size get_option_size(std::string const& key, platform::Size const& t_default) const;
	std::string get_option_string(std::string const& key) const;
	std::string get_option_string(std::string const& key, std::string const& t_default) const;
#endif

	// See also the explicit specialization of getOption() below.

	options_t const& getOptions() const;
	void setOptions( options_t const& options );

	void read(std::istream& in);
	void read( std::string const & tag_string );
	void write(std::ostream& out, int num_tabs = 0 ) const;

	///@brief returns the string that would be written by write()
	std::string to_string( int num_tabs = 0 ) const;

	TagCOP const &
	operator[](std::string const& key) const;

	Tag &
	operator=(Tag const &other);

	/// @brief Recursively reset that accessed_options_ variable so that a re-parsing of the tag
	/// can identify options that have been given but that have not been read.
	void reset_accessed_options() const;

	void die_for_unaccessed_options() const;
	void die_for_unaccessed_options_recursively() const;

	static
	TagOP create(std::istream& in); // creates a new tag and reads into it

	static
	TagOP create(std::string const &instring); // creates a new tag and reads into it

	TagOP clone() const;

	/// @brief if true, options will be quoted when the tag is outputted
	///        if false, options will be left as-is (default)
	/// @param[in] quote_options_val Whether or not option values should be quoted.
	///                              Default=false
	void
	set_quote_options( bool const quote_options_val );

private:
	std::string name_;
	options_t mOptions_;
	mutable options_t accessed_options_;
	tags_t vTags_;
	std::map<std::string,tags_t> mvTags_;
	static utility::vector0<TagCOP> const vEmpty_; // need to return this from getTags
	TagCAP parentTag_;

	/// @brief if true, options will be quoted when the tag is outputted (default)
	///        if false, options will be left without quotes
	bool quote_options_;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; // class Tag

std::ostream& operator<<(std::ostream& out, Tag const & tag);
std::ostream& operator<<(std::ostream& out, TagCOP const & tag);

/// @brief This is explicit specialization for boolean values
/// to allow for use of "true" "false" etc. in addition to 1 and 0
template<>
bool
Tag::getOption<bool>(std::string const& key, bool const& t_default) const;

/// @brief This is for the variant in which someone has specified a default
/// using "true" instead of true, or "false" instead of false.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
template<>
bool
Tag::getOption<bool>(std::string const& key, char const * default_as_string_literal) const;

template<>
bool
Tag::getOption<bool>(std::string const& key) const;

template<>
AutoBool
Tag::getOption<AutoBool>(std::string const& key, AutoBool const& t_default) const;

/// @brief This is for the variant in which someone has specified a default
/// using "true" instead of true, or "false" instead of false.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
template<>
AutoBool
Tag::getOption<AutoBool>(std::string const& key, char const * default_as_string_literal) const;

template<>
AutoBool
Tag::getOption<AutoBool>(std::string const& key) const;

/// @brief Special-casing the string literal version for string options.  In this case,
/// there shouldn't be an error thrown.  A string literal should be allowed to set the
/// default value for a string.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
template<>
std::string
Tag::getOption<std::string>( std::string const & key, char const * default_as_string_literal ) const;

// @brief If this were uncommented, this would add a special-case treatment to ensure that integer defaults
// get interpreted as Reals when setting a Real option.  So, for example, a developer could write
// tag->getOption<core::Real>("myoption", 1) instead of tag->getOption<core::Real>("myoption", 1.0).
// @details This has been deliberately REMOVED because there's value in not allowing a Real's default to be set with a Size
// -- removing this revealed a number of actual errors.  I'm leaving this here, commented out, in case we someday
// decide to allow a Real option to be given a Size default.
// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
// template<>
// platform::Real
// Tag::getOption<platform::Real>( std::string const & key, int const default_int ) const;

#ifdef OLDER_CLANG
/// @brief Special-casing the string literal version for integer options, too.
/// @note Only needed for older clang compilers, which have trouble with the general case being deleted and
/// specializations being provided.
/// @details Needed to be deleted explicitly since 0 gets interpreted as a char const * and not a Size.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
template<>
platform::Size
Tag::getOption<platform::Size>( std::string const &, char const * ) const = delete;

/// @brief Special-casing the string literal version for Real options, too.
/// @note Only needed for older clang compilers, which have trouble with the general case being deleted and
/// specializations being provided.
/// @details Needed to be deleted explicitly since 0 gets interpreted as a char const * and not a Real.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
template<>
platform::Real
Tag::getOption<platform::Real>( std::string const &, char const * ) const = delete;

/// @brief Special-casing the string literal version for single-precision float options, too.
/// @note Only needed for older clang compilers, which have trouble with the general case being deleted and
/// specializations being provided.
/// @details Needed to be deleted explicitly since 0 gets interpreted as a char const * and not a float.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
template<>
float
Tag::getOption<float>( std::string const &, char const * ) const = delete;
#endif

/// @brief Special-casing to ensure that 0 gets interpreted as Size(0) rather than nullptr.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
template<>
platform::Size
Tag::getOption<platform::Size>( std::string const & key, int const default_int ) const;

/// @brief Special-casing for signed ints.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
template<>
int
Tag::getOption<int>( std::string const & key, int const default_int ) const;

/// @brief Special-casing for int64_t.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
template<>
int64_t
Tag::getOption<int64_t>( std::string const & key, int const default_int ) const;


} // namespace tag
} // namespace utility

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( utility_tag_Tag )
#endif // SERIALIZATION


#endif // INCLUDED_utility_tag_tag_HH
