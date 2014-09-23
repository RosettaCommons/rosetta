// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/PoseCommentFilter.hh

#ifndef INCLUDED_protocols_simple_filters_PoseCommentFilter_hh
#define INCLUDED_protocols_simple_filters_PoseCommentFilter_hh

//unit headers
#include <protocols/simple_filters/PoseCommentFilter.fwd.hh>

// Project Headers
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>

namespace protocols {
namespace simple_filters {

///@brief test whether a pose contains a comment that evaluates to a predefined value. This is useful in controlling execution flow in RosettaScripts.
class PoseComment : public filters::Filter
{
  public:
    PoseComment();
    virtual ~PoseComment();
		filters::FilterOP clone() const {
			return filters::FilterOP( new PoseComment( *this ) );
		}
		filters::FilterOP fresh_instance() const{
			return filters::FilterOP( new PoseComment() );
		}

		virtual bool apply( core::pose::Pose const & pose ) const;
		virtual void report( std::ostream & out, core::pose::Pose const & pose ) const;
		virtual core::Real report_sm( core::pose::Pose const & pose ) const;
		void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &filters, moves::Movers_map const &, core::pose::Pose const & );
		core::Real compute( core::pose::Pose const & pose ) const;

		std::string comment_name() const{ return comment_name_; }
		void comment_name( std::string const s ){ comment_name_ = s; }

		std::string comment_value() const{ return comment_value_; }
		void comment_value( std::string const s ){ comment_value_ = s; }

		bool comment_exists() const { return comment_exists_; }
		void comment_exists( bool const c ){ comment_exists_ = c; }

  private:
		std::string comment_name_; //dflt ""; define the comment name
		std::string comment_value_; //dflt ""; define the comment value
	  bool comment_exists_; //dflt false; simply test whether or not the comment is there. If it is, return true, regardless of its value
};
}
}

#endif
