// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ./src/protocols/fldsgn/filters/HelixKinkFilter.hh
/// @brief header file for HelixKinkFilter class.
/// @details
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )


#ifndef INCLUDED_protocols_fldsgn_filters_HelixKinkFilter_hh
#define INCLUDED_protocols_fldsgn_filters_HelixKinkFilter_hh

// Unit Headers
#include <protocols/fldsgn/filters/HelixKinkFilter.fwd.hh>

// Package Headers
#include <protocols/filters/Filter.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>

// Utility headers

// Parser headers
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <utility/tag/Tag.fwd.hh>



//// C++ headers

namespace protocols {
namespace fldsgn {
namespace filters {

class HelixKinkFilter : public protocols::filters::Filter {
public:


	typedef protocols::filters::Filter Super;
	typedef protocols::filters::Filter Filter;
	typedef protocols::filters::FilterOP FilterOP;

	typedef std::string String;
	typedef core::Real Real;
	typedef core::Size Size;

	typedef core::pose::Pose Pose;

	typedef utility::tag::TagCOP TagCOP;
	typedef basic::datacache::DataMap DataMap;


public:// constructor/destructor


	// @brief default constructor
	HelixKinkFilter();

	~HelixKinkFilter() override{}


public:// virtual constructor


	// @brief make clone
	FilterOP clone() const override { return utility::pointer::make_shared< HelixKinkFilter >( *this ); }

	// @brief make fresh instance
	FilterOP fresh_instance() const override { return utility::pointer::make_shared< HelixKinkFilter >(); }


public:// accessors and mutators

	/// @brief Set the helix bend angle.
	void set_bend_angle( Real const bend_angle );

	/// @brief Set the secondary structure string.
	void set_secstruct( String const &secstruct );

	/// @brief Set the resnums (comma-separated string to be evaluated)
	/// @details This function sets string_resnums_ and select_resnums_.  select_resnums_ will be set to arguments'
	/// value.  If an empty string is passed (default), select_resnums_ will be set to false.  Otherwise,
	/// select_resnums_ will be set to true.
	void set_resnums( String const &resnums = "" );

	/// @brief Set helix_start_, helix_end_, and select_range_
	/// @details This function sets helix_start_, helix_end_, and sets select_range_ to true
	void set_helix_range( core::Size const helix_start, core::Size const helix_end );

	/// @brief Set select_range_
	void set_select_range( bool const select_range );

	/// @brief Get the helix bend angle.
	Real get_bend_angle() const;

	/// @brief Get the secondary structure string.
	String get_secstruct() const;

	/// @brief Get the resnum string.
	String get_string_resnums() const;

	/// @brief Get select_resnums_
	bool get_select_resnums() const;

	/// @brief Get helix_start_
	core::Size get_helix_start() const;

	/// @brief Get helix_end_
	core::Size get_helix_end() const;

	/// @brief Get select_range_
	bool get_select_range() const;

public:// parser


	void parse_my_tag( TagCOP tag,
		basic::datacache::DataMap &
	) override;


public:// virtual main operation


	// @brief returns true if the given pose passes the filter, false otherwise.
	// In this case, the test is whether the give pose is the topology we want.
	bool apply( Pose const & pose ) const override;

	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );



private:


	/// @brief
	Real bend_angle_;

	/// @brief
	String secstruct_;

	/// @brief
	String string_resnums_;
	bool select_resnums_;
	bool select_range_;
	core::Size helix_start_;
	core::Size helix_end_;

};

} // filters
} // fldsgn
} // protocols

#endif
