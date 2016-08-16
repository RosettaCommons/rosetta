// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/simple_filters/DomainInterfaceFilter.hh
/// @brief header file for DomainInterfaceFilter which checks whether a defined region of a protein
/// is buried in an interface by evaluating the number of non-water neighbors of all residues
///
/// @author Robert Lindner <rlindner@mpimf-heidelberg.mpg.de>

#ifndef INCLUDED_protocols_simple_filters_DomainInterfaceFilter_hh
#define INCLUDED_protocols_simple_filters_DomainInterfaceFilter_hh

// Unit Headers
#include <protocols/filters/Filter.hh>
#include <protocols/simple_filters/DomainInterfaceFilter.fwd.hh>

// Package Headers
#include <basic/datacache/DataMap.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/tag/Tag.fwd.hh>

//// C++ headers
#include <string>
#include <set>

#ifdef WIN32
#include <utility/tag/Tag.hh>
#endif

namespace protocols {
namespace simple_filters {

class DomainInterfaceFilter : public filters::Filter {
public:
	DomainInterfaceFilter();
	DomainInterfaceFilter( DomainInterfaceFilter const & );
	virtual ~DomainInterfaceFilter();

	std::string const & query_region_string();
	std::string const & target_region_string();
	std::set< core::Size > const & query_region();
	std::set< core::Size > const & target_region();

	void query_region( std::set< core::Size > & region );
	void target_region( std::set< core::Size > & region );
	void query_region_string( std::string const & region_str  );
	void target_region_string( std::string const & region_str  );

	core::Real cb_dist_cut() const;
	core::Real nearby_atom_cut() const;
	core::Real vector_angle_cut() const;
	core::Real vector_dist_cut() const;
	void cb_dist_cut( core::Real setting );
	void nearby_atom_cut( core::Real setting );
	void vector_angle_cut( core::Real setting );
	void vector_dist_cut( core::Real setting );

	std::string get_type() const{ return type_; }
	std::string get_user_defined_name() const { return user_defined_name_; }
	void set_user_defined_name( std::string const & name ) { user_defined_name_ = name; };

	/// @brief used to clear internal variables if needed. Using fresh_instance is preferred since it's a pure virtual
	//void clear() {};

	void parse_my_tag(
		utility::tag::TagCOP,
		basic::datacache::DataMap &,
		filters::Filters_map const &,
		moves::Movers_map const &,
		core::pose::Pose const &
	);

	virtual filters::FilterOP clone() const;
	virtual filters::FilterOP fresh_instance() const;

	/// @brief Returns true if the given pose passes the filter, false otherwise.
	virtual bool apply( core::pose::Pose const & pose ) const;

	virtual
	std::string name() const { return "DomainInterfaceFilter"; };

private:
	std::string const type_;
	std::string user_defined_name_;
	core::Real cb_dist_cut_;
	core::Real nearby_atom_cut_;
	core::Real vector_angle_cut_;
	core::Real vector_dist_cut_;
	std::set< core::Size > query_region_;
	std::set< core::Size > target_region_;
	std::string query_region_str_;
	std::string target_region_str_;
};

} // simple_filters
} // protocols

#endif
