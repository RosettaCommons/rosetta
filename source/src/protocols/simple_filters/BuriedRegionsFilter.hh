// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/simple_filters/BuriedRegionsFilter.hh
/// @brief header file for BuriedRegionsFilter which checks whether a defined region of a protein
/// is buried in an interface by evaluating the number of non-water neighbors of all residues
///
/// @author Robert Lindner <rlindner@mpimf-heidelberg.mpg.de>

#ifndef INCLUDED_protocols_simple_filters_BuriedRegionsFilter_hh
#define INCLUDED_protocols_simple_filters_BuriedRegionsFilter_hh

// Unit Headers
#include <protocols/filters/Filter.hh>
#include <protocols/simple_filters/BuriedRegionsFilter.fwd.hh>

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

class BuriedRegionsFilter : public filters::Filter {
public:
	BuriedRegionsFilter();
	BuriedRegionsFilter( core::Real distance_cutoff, core::Size neighbor_cutoff );
	BuriedRegionsFilter( BuriedRegionsFilter const & );
	virtual ~BuriedRegionsFilter();

	core::Real distance_cutoff();
	void distance_cutoff( core::Real distance_cutoff );
	core::Size neighbor_cutoff();
	void neighbor_cutoff( core::Size neighbor_cutoff );
	std::string const & region_string();
	void region_string( std::string const & region_str  );
	std::set< core::Size > const & region();
	void region( std::set< core::Size > & region );
  
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
	std::string name() const { return "BuriedRegionsFilter"; };

private:
	std::string const type_;
	std::string user_defined_name_;
	std::set< core::Size > region_;
	std::string region_str_;
	core::Real distance_cutoff_;
	core::Size neighbor_cutoff_;	
};

} // simple_filters
} // protocols

#endif
