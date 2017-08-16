// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file MotifScoreFilter.cc
/// @brief Filter for MotifScore
/// @author TJ Brunette

// Unit headers

#ifndef INCLUDED_protocols_simple_filters_MotifScoreFilter_hh
#define INCLUDED_protocols_simple_filters_MotifScoreFilter_hh

#include <protocols/filters/Filter.hh>
#include <protocols/simple_filters/MotifScoreFilter.fwd.hh>

#include <core/scoring/motif/motif_hash_stuff.hh>
#include <core/scoring/motif/util.hh>
#include <core/pose/motif/reference_frames.hh>

#include <core/pose/Pose.fwd.hh>

// C++ Headers
#include <string>
#include <utility/vector1.hh>


namespace protocols {
namespace simple_filters {

class MotifScoreFilter : public protocols::filters::Filter {
public:
	MotifScoreFilter();

	~MotifScoreFilter() override = default;

	filters::FilterOP clone() const override {
		return filters::FilterOP( new MotifScoreFilter( *this ) ); }

	filters::FilterOP fresh_instance() const override{
		return filters::FilterOP( new MotifScoreFilter() );
	}


	/// @brief Returns true if the given pose passes the filter, false otherwise.

	bool apply( core::pose::Pose const & pose ) const override;

	void report( std::ostream & out, core::pose::Pose const & pose ) const override;

	core::Real report_sm( core::pose::Pose const & pose ) const override;

	core::Real compute( core::pose::Pose const &pose ) const;

	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) override;

	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );



private:
	core::scoring::motif::MotifHashManager *mman_;
	core::Real score_threshold_;
};

} // filters
} // protocols

#endif
