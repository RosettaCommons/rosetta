// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/SigmoidFilter.hh

#ifndef INCLUDED_protocols_simple_filters_SigmoidFilter_hh
#define INCLUDED_protocols_simple_filters_SigmoidFilter_hh

//unit headers
#include <protocols/simple_filters/SigmoidFilter.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>

namespace protocols {
namespace simple_filters {

///@brief transform the output from a regular filter into a sigmoid ranging from 0-1 according to:
/// fx = 1/[1 + exp[ ( x - offset ) * steepness ]
/// The function asymptotically reaches 1 at negative values and 0 at positive values. It's 0.5 at the offset and steepness determines its slope at the offset
class Sigmoid : public filters::Filter
{
  public:
    Sigmoid();
    virtual ~Sigmoid();
		filters::FilterOP clone() const {
			return filters::FilterOP( new Sigmoid( *this ) );
		}
		filters::FilterOP fresh_instance() const{
			return filters::FilterOP( new Sigmoid() );
		}

		virtual bool apply( core::pose::Pose const & pose ) const;
		virtual void report( std::ostream & out, core::pose::Pose const & pose ) const;
		virtual core::Real report_sm( core::pose::Pose const & pose ) const;
		void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &filters, moves::Movers_map const &, core::pose::Pose const & );
		core::Real compute( core::pose::Pose const & pose ) const;
    core::Real steepness() const{ return steepness_; }
    void steepness( core::Real const s ){ steepness_ = s; }
    core::Real offset(){ return offset_; }
    void offset( core::Real const o ){ offset_ = o; }
    bool negate() const{ return negate_; }
    void negate( bool const b ){ negate_ = b; }
    protocols::filters::FilterOP filter() const;
    void filter( protocols::filters::FilterOP f );
		void reset_baseline( core::pose::Pose const & pose, bool const attempt_read_from_checkpoint ); /// allows within-trajectory resetting of the baseline. Notice this is nonconst, so can't be called from apply. attempt_read_from_checkpoint should be true for MC trials > 1, but false otherwise
  	core::Real threshold() const{ return threshold_; }
  	void threshold( core::Real const t ){ threshold_ = t; }
		void baseline_checkpointing_filename( std::string const s ){ baseline_checkpointing_filename_ = s; }
		std::string baseline_checkpointing_filename() const{ return baseline_checkpointing_filename_; }
  private:
    protocols::filters::FilterOP filter_; /// dflt NULL
    core::Real steepness_; //dflt 1
    core::Real offset_;/// dflt 0
		core::Real baseline_; /// dflt 0; this is tricky; used internally to keep track of where the pose started. It is only reset by reset_baseline, and cannot, due to constness, be changed by apply
    bool negate_; /// dflt false
		core::Real threshold_; /// dflt 0 (always accept)
		std::string baseline_checkpointing_filename_; // dflt ""; If this is set, the baseline value is saved to a checkpointing file to allow for recovery of the baseline after failure.
};
}
}

#endif
