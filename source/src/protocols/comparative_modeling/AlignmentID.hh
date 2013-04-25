// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/id/AlignmentID.hh
/// @author James Thompson


#ifndef INCLUDED_protocols_comparative_modeling_AlignmentID_hh
#define INCLUDED_protocols_comparative_modeling_AlignmentID_hh

#include <core/types.hh>

#include <ObjexxFCL/string.functions.hh>

#include <sstream>


namespace protocols {
namespace id {

/// @brief Unique template identifier class
class AlignmentID
{

public: // Creation

	/// @brief Default constructor
	inline
	AlignmentID() :
		align_idx_( 0 ),
		template_name_( "empty" )
	{}

	/// @brief Copy constructor
	inline
	AlignmentID( AlignmentID const & src ) :
		align_idx_( src.align_idx() ),
		template_name_( src.template_name() )
	{}

	/// @brief Property constructor
	inline
	AlignmentID(
		core::Size const align_idx_in,
		std::string const template_name_in
	) :
		align_idx_( align_idx_in ),
		template_name_( template_name_in )
	{}

public: // Properties

	inline
	std::string
	template_name() const {
	 	return template_name_;
	}

	inline
	std::string &
	template_name() {
		return template_name_;
	}

	inline
	core::Size
	align_idx() const {
		return align_idx_;
	}

	inline
	core::Size &
	align_idx() {
		return align_idx_;
	}

	void
	template_name( std::string const & name ) {
		template_name_ = name;
	}

	void
	align_idx( core::Size const & idx ) {
		align_idx_ = idx;
	}

	/// @brief Is this id valid?
	inline
	bool
	valid() const {
		return ( align_idx_ > 0 ) && ( template_name_ != "" );
	}

	inline
	std::string to_string() const {
		std::ostringstream out;
		out << template_name() << "_" << align_idx();
		return out.str();
	}

public: // Friends

	friend
	std::ostream &
	operator <<(
		std::ostream & os,
		AlignmentID const & a
	) {
		os << a.to_string();
		return os;
	}

	friend
	std::istream &
	operator >>(
		std::istream & is,
		AlignmentID & a
	) {
		using std::string;

		string tag;
		is >> tag;
		string::size_type pos = tag.find_first_of("_", 0);
		string t_name = tag.substr(0,pos);
		string index  = tag.substr(pos);

		a.template_name( t_name );
		a.align_idx( static_cast< core::Size > (ObjexxFCL::int_of( index )) );

		return is;
	}

	/// @brief a and b are the same atom
	friend
	inline
	bool
	operator ==(
		AlignmentID const & a,
		AlignmentID const & b
	) {
		return a.align_idx() == b.align_idx() &&
			a.template_name() == b.template_name();
	}

	/// @brief a and b are different atom
	friend
	inline
	bool
	operator !=(
		AlignmentID const & a,
		AlignmentID const & b
	) {
		return !( a == b );
	}

	/// @brief a is LOWER than b (e.g., first by smaller template name then
	/// by smaller alignment index)
	friend
	inline
	bool
	operator <(
		AlignmentID const & a,
		AlignmentID const & b
	) {
		return ( a.template_name_ <  b.template_name_ ||
					 ( a.template_name_ == b.template_name_ && a.align_idx_ < b.align_idx_ ) );
	}

private: // Fields
	/// @brief numerical identifier for Nth alignment to this template
	core::Size align_idx_;

	/// @brief template name
	std::string template_name_;
}; // AlignmentID

} // namespace id
} // namespace protocols


#endif // INCLUDED_protocols_comparative_modeling_AlignmentID_HH
