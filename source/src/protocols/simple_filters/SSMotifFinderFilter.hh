// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/SSMotifFinderFilter.hh

#ifndef INCLUDED_protocols_simple_filters_SSMotifFinderFilter_hh
#define INCLUDED_protocols_simple_filters_SSMotifFinderFilter_hh

//unit headers
#include <protocols/simple_filters/SSMotifFinderFilter.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/kinematics/Jump.hh>

namespace protocols {
namespace simple_filters {

class SSMotifFinder : public filters::Filter
{
  public:
    SSMotifFinder();
    virtual ~SSMotifFinder();
		filters::FilterOP clone() const {
			return new SSMotifFinder( *this );
		}
		filters::FilterOP fresh_instance() const{
			return new SSMotifFinder();
		}

		virtual bool apply( core::pose::Pose const & pose ) const;
		virtual void report( std::ostream & out, core::pose::Pose const & pose ) const;
		virtual core::Real report_sm( core::pose::Pose const & pose ) const;
		void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &filters, moves::Movers_map const &, core::pose::Pose const & );

		core::Size from_res() const{ return from_res_; }
		void from_res( core::Size const r ){ from_res_ = r; }
		core::Size to_res() const{ return to_res_; }
		void to_res( core::Size const r ){ to_res_ = r; }
		core::Size template_stem1() const{ return template_stem1_; }
		void template_stem1( core::Size const r ){ template_stem1_ = r; }
		core::Size template_stem2() const{ return template_stem2_; }
		void template_stem2( core::Size const r ){ template_stem2_ = r; }
		core::Real rmsd() const{ return rmsd_; }
		void rmsd( core::Real const c ){ rmsd_ = c ;}

		core::kinematics::Jump jump() const;
		core::kinematics::Jump compute_jump( core::pose::Pose const & pose, core::Size const start, core::Size const end ) const;

		// Undefined, commenting out to fix PyRosetta build  core::Real compute_rmsd() const;

		std::string filename() const{ return filename_; }
		void filename( std::string const s ){ filename_ = s;}
	private:
		core::Size from_res_, to_res_; // dflt 0,0; minimal and maximal sequence separation for motif search
		core::pose::PoseOP template_pose_; // dflt NULL
		core::Size template_stem1_, template_stem2_; //dflt 0,0; the start and end sites of the motif relative to the template
		core::kinematics::Jump jump_; //dflt NULL; template jump information. Computed once at parse time and used during apply
		core::Real rmsd_; //dflt 0; the maximal RMSd
		std::string filename_; //dflt ""; the file name into which to write the report
};
}
}

#endif
