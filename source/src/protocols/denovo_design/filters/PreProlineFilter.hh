// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/denovo_design/filters/PreProlineFilter.hh
/// @brief Tom's Denovo Protocol. This is freely mutable and used for playing around with stuff
/// @details
/// @author Tom Linsky (tlinsky@gmail.com)


#ifndef INCLUDED_protocols_denovo_design_filters_PreProlineFilter_hh
#define INCLUDED_protocols_denovo_design_filters_PreProlineFilter_hh

// Unit headers
#include <protocols/denovo_design/filters/PreProlineFilter.fwd.hh>

// Project headers
#include <protocols/filters/Filter.hh>

// Core headers
#include <core/pack/task/residue_selector/ResidueSelector.fwd.hh>

// C++ headers

namespace protocols {
namespace denovo_design {
namespace filters {

class PreProlineFilter : public protocols::filters::Filter {
public:

	/// @brief Initialize PreProlineFilter
	PreProlineFilter();

	/// @brief virtual constructor to allow derivation
	virtual ~PreProlineFilter();

	/// @brief Parses the PreProlineFilter tags
	void parse_my_tag(
			utility::tag::TagCOP tag,
			basic::datacache::DataMap & data,
			protocols::filters::Filters_map const &,
			protocols::moves::Movers_map const &,
			core::pose::Pose const & );

	/// @brief Return the name of this mover.
	virtual std::string get_name() const;

	/// @brief return a fresh instance of this class in an owning pointer
	virtual protocols::filters::FilterOP clone() const;
	virtual protocols::filters::FilterOP fresh_instance() const;

	/// @brief Apply the PreProlineFilter. Overloaded apply function from filter base class.
	virtual core::Real report_sm( core::pose::Pose const & pose ) const;
	virtual core::Real compute( core::pose::Pose const & pose ) const;
	void report( std::ostream & out, core::pose::Pose const & pose ) const;
	virtual bool apply( core::pose::Pose const & pose ) const;

	void set_selector( core::pack::task::residue_selector::ResidueSelectorCOP selector );
private:   // options
	/// @brief If total calculated filter score is <= theshold_, the filter passes, otherwise it fails.
	core::Real threshold_;

private:   // other data
	core::pack::task::residue_selector::ResidueSelectorCOP selector_;

};


} // filters
} // denovo_design
} // protocols

#endif
