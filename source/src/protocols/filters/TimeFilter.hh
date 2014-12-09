// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/filters/TimeFilter.hh
/// @author Sarel Fleishman (sarelf@uw.edu)

#ifndef INCLUDED_protocols_filters_TimeFilter_hh
#define INCLUDED_protocols_filters_TimeFilter_hh


// Project Headers
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>

#include <utility/vector1.hh>

#ifdef WIN32
	#include <ctime>
#endif

namespace protocols {
namespace filters {

class TimeFilter : public Filter
{
public:
	/// @brief default ctor
	TimeFilter();
	virtual bool apply( core::pose::Pose const & pose ) const;
	virtual core::Real compute( core::pose::Pose const & pose ) const;
	virtual void report( std::ostream & out, core::pose::Pose const & pose ) const;
	virtual core::Real report_sm( core::pose::Pose const & pose ) const;
	virtual FilterOP clone() const {
		return FilterOP( new TimeFilter( *this ) );
	}
	virtual FilterOP fresh_instance() const{
		return FilterOP( new TimeFilter() );
	}

	virtual ~TimeFilter();
	virtual void parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & );
private:
	core::Size mutable tic_, toc_;
};

} // filters
} // protocols

#endif
