// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief An object for reading/writing a simple xml-like language.
/// @author Paul Murphy

#ifndef INCLUDED_utility_tag_tag_HH
#define INCLUDED_utility_tag_tag_HH

// Unit headers
#include <utility/tag/Tag.fwd.hh>

// C++ Headers
#include <iosfwd>
#include <map>
#include <sstream>
#include <string>
#include <utility/assert.hh>
#include <iostream>

// Utility headers
#include <utility/exit.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector0.hh>

// Boost headers
#include <boost/lexical_cast.hpp>

namespace utility {
namespace tag {

class Tag : public utility::pointer::ReferenceCount, public utility::pointer::enable_shared_from_this< Tag >
{
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~Tag();

	typedef std::map<std::string, std::string > options_t;
	typedef utility::vector0< TagCOP > tags_t;

public:
	Tag();

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
	bool hasTag( std::string const& name ) const;
	bool hasOption( std::string const& key ) const;

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

	/// @brief Retrieve an option from the Tag with the given key name, using the
	/// provided default value (t_default) if the option is not present in the tag.
	/// @throws Throws a utility::excn::EXCN_Msg_Exception if the boost::lexical_cast
	/// fails to convert the input type as requested.
	template< class T >
	T
	getOption(std::string const& key, T const& t_default) const {
		options_t::const_iterator i = mOptions_.find(key);
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
			throw utility::excn::EXCN_Msg_Exception( error_message.str() );
		}
		return t_default;
	}

	/// @brief Retrieve an option from the Tag with the given key name.
	/// @throws Throws a utility::excn::EXCN_Msg_Exception if the an option
	/// with the given key is not present, or if the boost::lexical_cast fails
	/// to convert the input type as requested.
	template< class T >
	T
	getOption(std::string const& key) const {
		options_t::const_iterator i = mOptions_.find(key);
		if ( i == mOptions_.end() ) {
			std::stringstream error_message;
			error_message << "Option " << key << " not found.\n";
			throw utility::excn::EXCN_Msg_Exception( error_message.str() );
		}
		accessed_options_[key]= i->second;
		try{
			return boost::lexical_cast<T>(i->second);
		} catch(boost::bad_lexical_cast &) {
			std::stringstream error_message;
			error_message << "getOption: key= " << key << " stream extraction failed! Tried to parse '" << i->second << "'\n";
			throw utility::excn::EXCN_Msg_Exception( error_message.str() );
		}
		return T(); // appease compiler
	}

	// See also the explicit specialization of getOption() below.

	options_t const& getOptions() const;
	void setOptions( options_t const& options );

	void read(std::istream& in);
	void write(std::ostream& out, int num_tabs = 0 ) const;

	TagCOP const &
	operator[](std::string const& key) const;

	Tag &
	operator=(Tag const &other);

	void die_for_unaccessed_options() const;
	void die_for_unaccessed_options_recursively() const;

	static
	TagOP create(std::istream& in); // creates a new tag and reads into it

	static
	TagOP create(std::string instring); // creates a new tag and reads into it

	TagOP clone() const;

private:
	std::string name_;
	options_t mOptions_;
	mutable options_t accessed_options_;
	tags_t vTags_;
	std::map<std::string,tags_t> mvTags_;
	static utility::vector0<TagCOP> const vEmpty_; // need to return this from getTags
	TagCAP parentTag_;

}; // class Tag

std::ostream& operator<<(std::ostream& out, Tag const & tag);
std::ostream& operator<<(std::ostream& out, TagCOP const & tag);

//This is explicit specialization for boolean values
//to allow for use of "true" "false" etc. in addition to 1 and 0
template<>
bool
Tag::getOption<bool>(std::string const& key, bool const& t_default) const;

template<>
bool
Tag::getOption<bool>(std::string const& key) const;

} // namespace tag
} // namespace utility

#endif // INCLUDED_utility_tag_tag_HH
