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

#ifndef INCLUDED_protocols_hotspot_hashing_movers_PlaceLSMProbe_hh
#define INCLUDED_protocols_hotspot_hashing_movers_PlaceLSMProbe_hh


// Project Headers
#include <string>

#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/moves/DataMap.fwd.hh>

#include <protocols/hotspot_hashing/movers/PlaceLSMProbe.fwd.hh>
#include <protocols/hotspot_hashing/SearchPattern.hh>

// Unit headers

namespace protocols
{
namespace hotspot_hashing
{
namespace movers
{

class PlaceLSMProbe : public protocols::moves::Mover
{
  private:
    typedef protocols::moves::Mover parent;
  public:
    PlaceLSMProbe();

    PlaceLSMProbe(
			std::string residue_name,
			VectorPair lsm_spec,
			core::Size lsm_id,
			core::Real angle_sampling,
			core::Real translocation_sampling,
			core::Real max_radius,
			core::Real distance_sampling,
			core::Real max_distance,
			bool log_pose,
			core::conformation::ResidueCOP target_residue);

    virtual void apply( Pose & );

    virtual std::string get_name() const { return "PlaceLSMProbe"; }

		virtual protocols::moves::MoverOP clone() const;

    void parse_my_tag(
         utility::tag::TagPtr const tag,
         protocols::moves::DataMap &,
         protocols::filters::Filters_map const &,
         protocols::moves::Movers_map const &,
         core::pose::Pose const &);

    virtual bool reinitialize_for_new_input() const { return true; }

  protected:
    void check_and_initialize();

  private:
		std::string residue_name_;

    VectorPair lsm_spec_;
		core::Size lsm_id_;

		core::Real angle_sampling_;
		core::Real translocation_sampling_;
		core::Real max_radius_;
		core::Real distance_sampling_;
		core::Real max_distance_;

    bool log_pose_;

    core::conformation::ResidueCOP target_residue_;

    bool initialized_pattern_;
    utility::vector1<TransformPair> search_points_;
};

} // movers
} // hotspot_hashing
} // protocols

#endif 
