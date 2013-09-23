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

#include <iosfwd>
#include <map>
#include <sstream>
#include <string>
//#include <vector>
#include <cassert>
#include <utility/tag/Tag.fwd.hh>
#include <utility/exit.hh>
#include <iostream>

#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector0.hh>
#include <boost/lexical_cast.hpp>

namespace utility {
namespace tag {

class Tag : public utility::pointer::ReferenceCount {
public:
	///@brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~Tag();

	typedef std::map<std::string, std::string > options_t;
	typedef utility::vector0< TagPtr > tags_t;

public:
	static size_t num_tags; // why?!

	Tag();

	void clear();
	size_t size() const;

	void setName(std::string const& name);
	std::string const& getName() const { return name_; }

	void addTag( TagPtr const& tag );
	utility::vector0< TagPtr > const & getTags() const;
	utility::vector0< TagPtr > const & getTags( std::string const& name ) const;
	TagPtr const & getTag( std::string const& name ) const;
	bool hasTag( std::string const& name ) const;
	bool hasOption( std::string const& key ) const;

	template< class T >
	void
	setOption(std::string const& key, T const& value) {
		if( mOptions_.find(key) != mOptions_.end() ) {
			//runtime_assert( false );
		}
		std::ostringstream out;
		out << value;
		mOptions_[key] = out.str();
		accessed_options_.erase(key);
	} // setOption

	template< class T >
	T
	getOption(std::string const& key, T const& t_default) const {
		options_t::const_iterator i = mOptions_.find(key);
		if( i == mOptions_.end() ) {
			accessed_options_[key]= key;
			return t_default;
		}
		accessed_options_[key]= i->second;
		//		T t = t_default; // not used?
		try{
			return boost::lexical_cast<T>(i->second);
		} catch(boost::bad_lexical_cast &) {
			std::cerr << "getOption: key= " << key << " stream extraction failed! Tried to parse '" << i->second <<
					"' returning default value: '" << t_default << "'" << std::endl;
		}
		return t_default;
	}

	template< class T >
	T
	getOption(std::string const& key) const {
		options_t::const_iterator i = mOptions_.find(key);
		if( i == mOptions_.end() ) {
			std::cerr << "Option " << key << " not found." << std::endl;
			runtime_assert( false );
		}
		accessed_options_[key]= i->second;
		try{
			return boost::lexical_cast<T>(i->second);
		} catch(boost::bad_lexical_cast &) {
			std::cerr << "getOption: key= " << key << " stream extraction failed! Tried to parse '" << i->second <<
				"' returning uninitialized value: '" << T() << "'" << std::endl;
		}
		return T();
	}

	// See also the explicit specialization of getOption() below.

	options_t const& getOptions() const { return mOptions_; }
	void setOptions( options_t const& options ) {
		mOptions_ = options;
		accessed_options_.clear();
	}

	void read(std::istream& in);
	void write(std::ostream& out, int num_tabs = 0 ) const;

	TagPtr const &
	operator[](std::string const& key) const {
		return getTag(key);
	} // operator[]

	void die_for_unaccessed_options();
	void die_for_unaccessed_options_recursively(){
		die_for_unaccessed_options();
		tags_t::const_iterator begin= vTags_.begin();
		for(; begin != vTags_.end(); ++begin){
			(*begin)->die_for_unaccessed_options_recursively();
		}
	}

	static
	TagPtr create(std::istream& in); // creates a new tag and reads into it

	static
	TagPtr create(std::string instring); // creates a new tag and reads into it

	TagPtr clone() const;

private:

	std::string name_;
	options_t mOptions_;
	mutable options_t accessed_options_;
	tags_t vTags_;
	std::map<std::string,tags_t> mvTags_;
	static utility::vector0<TagPtr> const vEmpty_; // need to return this from getTags

}; // class Tag

std::ostream& operator<<(std::ostream& out, Tag const& tag);
std::ostream& operator<<(std::ostream& out, TagPtr const& tag);

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
