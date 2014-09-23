// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/FoldUnitsFilter.hh
/// @brief Simple filter that tests whether a file exists. Useful to test whether we're recovering from a checkpoint
/// @author Sarel Fleishman

#ifndef INCLUDED_protocols_simple_filters_FoldUnitsFilter_hh
#define INCLUDED_protocols_simple_filters_FoldUnitsFilter_hh

//unit headers
#include <devel/splice/FoldUnitsFilter.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <string>

namespace devel {
namespace splice {

class FoldUnitsFilter : public protocols::filters::Filter
{
public:
	//default ctor
	FoldUnitsFilter();
	bool apply( core::pose::Pose const & pose ) const;
	protocols::filters::FilterOP clone() const {
		return protocols::filters::FilterOP( new FoldUnitsFilter( *this ) );
	}
	protocols::filters::FilterOP fresh_instance() const{
		return protocols::filters::FilterOP( new FoldUnitsFilter() );
	}

	core::Real compute( core::pose::Pose const & pose ) const;
	void report( std::ostream & out, core::pose::Pose const & pose ) const;
	core::Real report_sm( core::pose::Pose const & pose ) const;
	virtual ~FoldUnitsFilter();
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );

	void write_to_file( std::string const pdb, core::Size const from_res, core::Size const to_res, std::string const dssp, std::string const sequence, core::pose::Pose const & pose ) const;

	core::Size minimal_length() const{ return minimal_length_; }
	void minimal_length( core::Size const s ){ minimal_length_ = s; }
	core::Size maximal_length() const{ return maximal_length_; }
	void maximal_length( core::Size const s ){ maximal_length_ = s; }
	core::Real ends_distance() const{ return ends_distance_; }
	void ends_distance( core::Real const c ){ ends_distance_ = c; }
	void filename( std::string const s ){ filename_ = s; }
	std::string filename() const{ return filename_; }

private:
	core::Size minimal_length_, maximal_length_; //dflt 20,40; the minimal and maximal length of the sliding window
	core::Real ends_distance_; //dflt 4.5A; the distance between the ends of the segment
	std::string filename_; // dflt "FoldUnits.out"; the output filename
};

}
}

#endif
