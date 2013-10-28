// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/NonSequentialNeighborsFilter.hh
/// @brief Simple filter that tests whether a file exists. Useful to test whether we're recovering from a checkpoint
/// @author Gabi Pszolla & Sarel Fleishman

#ifndef INCLUDED_protocols_simple_filters_NonSequentialNeighborsFilter_hh
#define INCLUDED_protocols_simple_filters_NonSequentialNeighborsFilter_hh

//unit headers
#include <protocols/simple_filters/NonSequentialNeighborsFilter.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>

namespace protocols {
namespace simple_filters {

class NonSequentialNeighborsFilter : public filters::Filter
{
public:
	//default ctor
	NonSequentialNeighborsFilter();
	bool apply( core::pose::Pose const & pose ) const;
	filters::FilterOP clone() const {
		return new NonSequentialNeighborsFilter( *this );
	}
	filters::FilterOP fresh_instance() const{
		return new NonSequentialNeighborsFilter();
	}

	void report( std::ostream & out, core::pose::Pose const & pose ) const;
	core::Real report_sm( core::pose::Pose const & pose ) const;
	core::Real compute( core::pose::Pose const &pose ) const;
	virtual ~NonSequentialNeighborsFilter();
	void parse_my_tag( utility::tag::TagCOP const tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );

	core::Real distance_threshold() const { return distance_threshold_; }
	void distance_threshold( core::Real const r ){ distance_threshold_ = r; }
	core::Size neighbor_cutoff() const {return neighbor_cutoff_; }
	void neighbor_cutoff( core::Size const n ) { neighbor_cutoff_ = n; }
	bool bound() const{ return bound_; }
	void bound( bool const b ) { bound_ = b;}
	core::Size resnum() const{ return resnum_; }
	void resnum( core::Size const s ){ resnum_ = s; }
	void jump( core::Size const j ){ jump_ = j; }
	core::Size jump() const { return jump_; }
	core::Size residue_neighbors( core::pose::Pose const & pose, core::Size const resi ) const; // number of residue neighbours outside the sequence neighbours.
private:
	core::Real distance_threshold_; // dflt 8.0A; sphere around the residue
	core::Size neighbor_cutoff_; // dflt 10; how many residues in the sequence around the target residue to ignore in computing neighbours
	bool bound_; // dflt false; treat the bound or unbound pose
	core::Size resnum_; // dflt 0; 0 - means all residues in the pose; any other number targets a particular residue
	core::Size jump_; // dflt 1;
};

}
}

#endif
