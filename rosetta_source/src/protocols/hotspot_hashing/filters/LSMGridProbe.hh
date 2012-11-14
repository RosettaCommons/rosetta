// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 sw=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/hotspot_hashing/filters/LSMGridProbe.hh
/// @brief
/// @author Alex Ford fordas@uw.edu

#ifndef INCLUDED_protocols_hotspot_hashing_filters_LSMGridProbe_hh
#define INCLUDED_protocols_hotspot_hashing_filters_LSMGridProbe_hh


// Project Headers
#include <string>

#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>

#include <protocols/hotspot_hashing/filters/LSMGridProbe.fwd.hh>
#include <protocols/hotspot_hashing/PlaceMinimizeSearch.hh>

// Unit headers

namespace protocols
{
namespace hotspot_hashing
{
namespace filters
{

class LSMGridProbe : public protocols::filters::Filter
{
  private:
    typedef protocols::filters::Filter parent;
  public:
    /// @brief default ctor
    LSMGridProbe();
    virtual ~LSMGridProbe();

    virtual bool apply( core::pose::Pose const & pose ) const;

    virtual protocols::filters::FilterOP clone() const;
    virtual protocols::filters::FilterOP fresh_instance() const;

    void parse_my_tag( utility::tag::TagPtr const tag,
                       protocols::moves::DataMap &,
                       protocols::filters::Filters_map const &,
                       protocols::moves::Movers_map const &,
                       core::pose::Pose const & );

		static bool parse_lsm_spec(std::string lsmspec, VectorPair & lsmspec);
  private:
		std::string residue_name_;

    VectorPair lsm_spec_;
		core::Size lsm_id_;

		core::Real angle_sampling_;
		core::Real translocation_sampling_;
		core::Real max_radius_;
		core::Real distance_sampling_;
		core::Real max_distance_;

    protocols::moves::MoverOP relax_mover_;
    protocols::filters::FilterOP triage_filter_;
};

} // filters
} // hotspot_hashing
} // protocols

#endif 
