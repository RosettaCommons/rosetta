// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/ReportFilter.hh
/// @brief Simple filter that tests whether a file exists. Useful to test whether we're recovering from a checkpoint
/// @author Sarel Fleishman

#ifndef INCLUDED_protocols_simple_filters_ReportFilter_hh
#define INCLUDED_protocols_simple_filters_ReportFilter_hh

//unit headers
#include <protocols/simple_filters/ReportFilter.fwd.hh>

// Project Headers
#include <core/scoring/ScoreFunction.hh>
#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <basic/datacache/DataMapObj.hh>

namespace protocols {
namespace simple_filters {

class ReportFilter : public filters::Filter
{
public:
	//default ctor
	ReportFilter();
	bool apply( core::pose::Pose const & pose ) const;
	filters::FilterOP clone() const {
		return filters::FilterOP( new ReportFilter( *this ) );
	}
	filters::FilterOP fresh_instance() const{
		return filters::FilterOP( new ReportFilter() );
	}

	void report( std::ostream & out, core::pose::Pose const & pose ) const;
	core::Real report_sm( core::pose::Pose const & pose ) const;
	core::Real compute( core::pose::Pose const &pose ) const;
	virtual ~ReportFilter();
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );

	void report_string( std::string const s );
	std::string report_string() const;

	void filter( protocols::filters::FilterOP f ){ filter_ = f; }
	protocols::filters::FilterOP filter() const{ return filter_; }
	std::string report_filter_name() const{ return report_filter_name_; }
	core::Real filter_val() const{ return filter_val_; }
	std::string checkpointing_file() const{ return checkpointing_file_; }
	void checkpointing_file( std::string const s ){ checkpointing_file_ = s ;}
private:
	void checkpoint_read() const; // read from a checkpoint
	void checkpoint_write() const; // write from a checkpoint
	utility::pointer::shared_ptr< basic::datacache::DataMapObj< std::string > > report_string_; //dflt ""
	protocols::filters::FilterOP filter_; //dflt NULL; either filter or report_string should be turned on
	std::string report_filter_name_; // the user defined filter name used in reporting
	mutable core::Real filter_val_; // stores filter_'s report_sm value at the time of apply call. Then, when report_sm for teh ReportFilter is called it reports the saved value
	std::string checkpointing_file_; //dflt ""; stores the filter_'s report_sm if the protocol is checkpointed. If a checkpointing file exists, the filter sticks to the value within the checkpointing file
};

}
}

#endif
