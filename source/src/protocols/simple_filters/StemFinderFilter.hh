// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/StemFinderFilter.hh

#ifndef INCLUDED_protocols_simple_filters_StemFinderFilter_hh
#define INCLUDED_protocols_simple_filters_StemFinderFilter_hh

//unit headers
#include <protocols/simple_filters/StemFinderFilter.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <string>
#include <utility/vector1.hh>

namespace protocols {
namespace simple_filters {

class StemFinder : public filters::Filter
{
public:
	StemFinder();
	virtual ~StemFinder();
	filters::FilterOP clone() const {
		return filters::FilterOP( new StemFinder( *this ) );
	}
	filters::FilterOP fresh_instance() const{
		return filters::FilterOP( new StemFinder() );
	}

	virtual bool apply( core::pose::Pose const & pose ) const;
	virtual void report( std::ostream & out, core::pose::Pose const & pose ) const;
	virtual core::Real report_sm( core::pose::Pose const & pose ) const;
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &filters, moves::Movers_map const &, core::pose::Pose const & );

	core::Real from_res() const{ return from_res_; }
	void from_res( core::Real const r ){ from_res_ = r; }
	core::Real to_res() const{ return to_res_; }
	void to_res( core::Real const r ){ to_res_ = r; }
	core::Real rmsd() const{ return rmsd_; }
	void rmsd( core::Real const c ){ rmsd_ = c ;}

	utility::vector1< std::string > filenames() const{ return filenames_; }
	void add_filename( std::string const & s );

	bool stems_on_sse() const{ return stems_on_sse_; }
	void stems_on_sse( bool const b ){ stems_on_sse_ = b; }
	core::Real neighbor_distance() const{ return neighbor_distance_; }
	void neighbor_distance( core::Real const n ) { neighbor_distance_ = n; }
	void stems_are_neighbors( bool const b ){ stems_are_neighbors_ = b; }
	bool stems_are_neighbors() const{ return stems_are_neighbors_; }
	void neighbor_separation( core::Size const c ){ neighbor_separation_ = c; }
	core::Size neighbor_separation() const{ return neighbor_separation_; }
private:
	core::Size from_res_, to_res_; // dflt 0,0; If not specified goes across the entire range of the -s pdb
	core::Real rmsd_; //dflt 0.7; the maximal RMSd. If no positions fall under this, the filter fails
	utility::vector1< std::string > filenames_; //dflt empty; the PDB file names to search
	bool stems_on_sse_; //dflt false; if false look for stems on any bb
	bool stems_are_neighbors_; //dflt true; stems need to be within a certain 3D distance
	core::Real neighbor_distance_; //dflt 4.0A
	core::Size neighbor_separation_; //dflt 10; at least 10 aa separation between residues
};

/// potentially useful utility functions

/// @brief read dssp for a pose and return a string
std::string dssp( core::pose::Pose const & pose );

/// @brief find helix, sheet positions in dssp
utility::vector1< core::Size > positions_in_secstruct( core::pose::Pose const & pose );

/// @brief load PDBs into a vector
utility::vector1< core::pose::PoseOP > load_poses( utility::vector1< std::string > const & filenames );

/// @brief find the minimal atom-atom distance between two residues
core::Real
res_res_min_distance( core::pose::Pose const &p1, core::Size const r1,
	core::pose::Pose const &p2, core::Size const r2 );

}
}

#endif
