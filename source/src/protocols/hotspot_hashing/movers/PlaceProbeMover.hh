// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 sw=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/hotspot_hashing/movers/PlaceProbeMover.hh
/// @brief
/// @author Alex Ford fordas@uw.edu

#ifndef INCLUDED_protocols_hotspot_hashing_movers_PlaceProbeMover_hh
#define INCLUDED_protocols_hotspot_hashing_movers_PlaceProbeMover_hh


// Project Headers
#include <string>

#include <utility/tag/Tag.fwd.hh>

#include <core/conformation/Residue.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.hh>

#include <protocols/hotspot_hashing/SearchPattern.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.hh>

// Unit headers

namespace protocols
{
namespace hotspot_hashing
{
namespace movers
{

class PlaceProbeMover : virtual public protocols::moves::Mover
{
  public:
		enum ExecutionMode
		{
			RunAll, OnePerStruct
		};

		enum StructureOutputMode
		{
			None, Probe, Full
		};

    PlaceProbeMover();

    PlaceProbeMover(
			std::string residue_name,
			core::conformation::ResidueCOP target_residue,
      core::Size search_partition = 1,
      core::Size total_search_partition = 1);

    virtual void apply( Pose & );

    virtual bool reinitialize_for_new_input() const { return false; }

  protected:
		///@brief Parses tag compoments for PlaceProbeMover
    void parse_place_probe_tag(
         utility::tag::TagCOP const tag,
         basic::datacache::DataMap &,
         protocols::filters::Filters_map const &,
         protocols::moves::Movers_map const &,
         core::pose::Pose const &);

	  StructureOutputMode parse_output_mode(std::string name);

    void check_and_initialize(core::pose::Pose const & target_pose);

    void execute_one_search(core::pose::Pose & target_pose, core::Size search_index);

		void perform_local_refinement(core::pose::Pose & target_pose, core::Size target_residue);

    virtual SearchPatternOP create_search_pattern(core::pose::Pose const & target_pose) = 0;
    virtual SearchPatternOP create_partitioned_search_pattern(core::pose::Pose const & target_pose);

		virtual SearchPatternOP create_refinement_pattern(core::pose::Pose const & target_pose, core::Size target_residue);
		virtual core::pack::task::PackerTaskOP create_refinement_packing_task(core::pose::Pose const & target_pose, core::Size target_residue);

		std::string residue_name_;
		std::vector<std::string> residue_variants_;

    core::conformation::ResidueCOP target_residue_;
		core::scoring::ScoreFunctionCOP refinement_scorefxn_;

		ExecutionMode current_mode_;

    core::Size search_partition_;
    core::Size total_search_partition_;

    bool initialized_pattern_;
    utility::vector1<core::kinematics::Stub> search_points_;
};

}
}
}

#endif 
