// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/hotspot_hashing/movers/PlaceSurfaceProbe.hh
/// @brief
/// @author Alex Ford fordas@uw.edu

#ifndef INCLUDED_protocols_hotspot_hashing_movers_PlaceSurfaceProbe_hh
#define INCLUDED_protocols_hotspot_hashing_movers_PlaceSurfaceProbe_hh


// Project Headers
#include <string>

#include <utility/tag/Tag.fwd.hh>

#include <core/conformation/Residue.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.hh>

#include <protocols/hotspot_hashing/SearchPattern.fwd.hh>
#include <protocols/hotspot_hashing/movers/PlaceProbeMover.hh>
#include <protocols/hotspot_hashing/movers/PlaceSurfaceProbe.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.hh>


// Unit headers

namespace protocols
{
namespace hotspot_hashing
{
namespace movers
{

class PlaceSurfaceProbe : public protocols::hotspot_hashing::movers::PlaceProbeMover
{
public:
	PlaceSurfaceProbe();

	PlaceSurfaceProbe(
		std::string residue_name,
		core::Real search_density,
		core::Real x_angle_sampling,
		core::Real y_angle_sampling,
		core::Real refinement_distance_sampling,
		core::Real refinement_distance,
		core::Real refinement_translation_sampling,
		core::conformation::ResidueCOP target_residue,
		core::pack::task::TaskFactoryOP surface_selection = NULL,
		core::Size search_partition = 1,
		core::Size total_search_partition = 1);

	PlaceSurfaceProbeOP shared_from_this() { return utility::pointer::dynamic_pointer_cast<PlaceSurfaceProbe>( PlaceProbeMover::shared_from_this() ); }


	virtual std::string get_name() const { return "PlaceSurfaceProbe"; }

	virtual protocols::moves::MoverOP clone() const;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &);

	virtual bool reinitialize_for_new_input() const { return false; }

protected:
	virtual SearchPatternOP create_search_pattern(core::pose::Pose const & target_pose);
	virtual SearchPatternOP create_partitioned_search_pattern(core::pose::Pose const & target_pose);
	virtual SearchPatternOP create_refinement_pattern(core::pose::Pose const & target_pose, core::Size target_residue);

private:

	core::pack::task::TaskFactoryOP surface_selection_;

	core::Real search_density_;
	core::Real coarse_angle_sampling_;
	core::Real coarse_sampling_;

	core::Real refinement_distance_;
	core::Real refinement_angle_sampling_;
	core::Real refinement_sampling_;

	SearchPatternOP initialize_refinement_pattern();
	SearchPatternOP refinement_pattern_;
};

}
}
}

#endif
