// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/filters/OperatorFilter.hh

#ifndef INCLUDED_protocols_filters_OperatorFilter_hh
#define INCLUDED_protocols_filters_OperatorFilter_hh

//unit headers
#include <protocols/filters/OperatorFilter.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>

namespace protocols {
namespace filters {

enum Operation { SUM, PRODUCT, NORMALIZED_SUM, MAX, MIN };
///@brief simply take a list of filters and combine them using the operation above
class Operator : public filters::Filter
{
  public:
    Operator();
    virtual ~Operator();
		filters::FilterOP clone() const {
			return new Operator( *this );
		}
		filters::FilterOP fresh_instance() const{
			return new Operator();
		}

		virtual bool apply( core::pose::Pose const & pose ) const;
		virtual void report( std::ostream & out, core::pose::Pose const & pose ) const;
		virtual core::Real report_sm( core::pose::Pose const & pose ) const;
		void parse_my_tag( utility::tag::TagPtr const tag, moves::DataMap &, filters::Filters_map const &filters, moves::Movers_map const &, core::pose::Pose const & );
		core::Real compute( core::pose::Pose const & pose ) const;
    utility::vector1< protocols::filters::FilterOP > filters() const;
    void add_filter( protocols::filters::FilterOP f );
		void reset_baseline( core::pose::Pose const & pose ); /// goes over Sigmoid filters and resets them. Note that this is nonconst, and cannot be called from apply
  	core::Real threshold() const{ return threshold_; }
  	void threshold( core::Real const t ){ threshold_ = t; }
		Operation operation() const{ return operation_; }
		void operation( Operation const o ){ operation_ = o; }
  private:
    utility::vector1< protocols::filters::FilterOP > filters_;
		Operation operation_; // dflt PRODUCT
		core::Real threshold_; // dflt 0
};
}
}

#endif
