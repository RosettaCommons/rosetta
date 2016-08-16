// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/FileExistFilter.hh
/// @brief Simple filter that tests whether a file exists. Useful to test whether we're recovering from a checkpoint
/// @author Sarel Fleishman

#ifndef INCLUDED_protocols_simple_filters_FileExistFilter_hh
#define INCLUDED_protocols_simple_filters_FileExistFilter_hh

//unit headers
#include <protocols/simple_filters/FileExistFilter.fwd.hh>

// Project Headers
#include <core/scoring/ScoreFunction.hh>
#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>

namespace protocols {
namespace simple_filters {

class FileExistFilter : public filters::Filter
{
public:
	//default ctor
	FileExistFilter();
	bool apply( core::pose::Pose const & pose ) const;
	filters::FilterOP clone() const {
		return filters::FilterOP( new FileExistFilter( *this ) );
	}
	filters::FilterOP fresh_instance() const{
		return filters::FilterOP( new FileExistFilter() );
	}

	void report( std::ostream & out, core::pose::Pose const & pose ) const;
	core::Real report_sm( core::pose::Pose const & pose ) const;
	core::Real compute( core::pose::Pose const &pose ) const;
	virtual ~FileExistFilter();
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );

	void filename( std::string const & f );
	std::string filename() const;
	bool ignore_zero_byte() const{ return ignore_zero_byte_; }
	void ignore_zero_byte( bool const i ){ ignore_zero_byte_ = i; }
private:
	std::string filename_;
	bool ignore_zero_byte_; // dflt false; if the file is 0b then it counts as if the file doesn't exist
};

}
}

#endif
