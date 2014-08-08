// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/CircularPermutationFinderFilter.hh

#ifndef INCLUDED_protocols_simple_filters_CircularPermutationFinderFilter_hh
#define INCLUDED_protocols_simple_filters_CircularPermutationFinderFilter_hh

//unit headers
#include <protocols/simple_filters/CircularPermutationFinder.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/kinematics/Jump.hh>

namespace protocols {
namespace simple_filters {

class CircularPermutationFinder : public filters::Filter
{
  public:
    CircularPermutationFinder();
    virtual ~CircularPermutationFinder();
		filters::FilterOP clone() const {
			return new CircularPermutationFinder( *this );
		}
		filters::FilterOP fresh_instance() const{
			return new CircularPermutationFinder();
		}

		virtual bool apply( core::pose::Pose const & pose ) const;
		virtual void report( std::ostream & out, core::pose::Pose const & pose ) const;
		virtual core::Real report_sm( core::pose::Pose const & pose ) const;
		void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &filters, moves::Movers_map const &, core::pose::Pose const & );

		core::Real rmsd() const{ return rmsd_; }
		void rmsd( core::Real const c ){ rmsd_ = c ;}

		std::string filename() const{ return filename_; }
		void filename( std::string const s ){ filename_ = s;}

		void N_C_distance( core::Real const r ){ N_C_distance_ = r;}
		core::Real N_C_distance() const{ return N_C_distance_; }
    void circular_permutation( core::pose::Pose const & pose, core::Size const chainA, core::Size const chainB ) const;

		bool align_only_interface() const{ return align_only_interface_; }
		void align_only_interface( bool const c ){ align_only_interface_ = c; }
	private:
		core::Real rmsd_; //dflt 0; the maximal RMSd
		core::Real N_C_distance_; //dflt ; maximal distance between N and C termini
		std::string filename_; //dflt ""; the file name into which to write the report
		bool align_only_interface_; //dflt true; true->align only the binding residues; false->align all residues in the span from the start to the end of the interface
};
}
}

#endif
