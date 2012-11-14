// Project Headers
// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 sw=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/hotspot_hashing/movers/PlaceProbeMover.cc
/// @brief
/// @author Alex Ford fordas@uw.edu
//

#include <basic/Tracer.hh>

#include <utility/tag/Tag.hh>
#include <utility/exit.hh>

#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pack/task/TaskFactory.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/moves/DataMap.fwd.hh>

#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

#include <protocols/rosetta_scripts/util.hh>

#include <protocols/hotspot_hashing/SearchPattern.hh>
#include <protocols/hotspot_hashing/SurfaceSearchPattern.hh>

#include <protocols/hotspot_hashing/StubGenerator.hh>

#include <protocols/hotspot_hashing/movers/PlaceProbeMover.hh>

namespace protocols
{
namespace hotspot_hashing
{
namespace movers
{

static basic::Tracer TR( "protocols.hotspot_hashing.movers.PlaceProbeMover" ); 

PlaceProbeMover::PlaceProbeMover() :
    residue_name_(""),
    target_residue_(NULL),
    search_partition_(0),
    total_search_partition_(1),
    initialized_pattern_(false),
    search_points_()
{} 

PlaceProbeMover::PlaceProbeMover(
		std::string residue_name,
    core::conformation::ResidueCOP target_residue,
		core::Size search_partition,
		core::Size total_search_partition) :
		residue_name_(residue_name),
    target_residue_(target_residue),
    search_partition_(search_partition),
    total_search_partition_(total_search_partition),
    initialized_pattern_(false),
    search_points_()
{}

void PlaceProbeMover::apply(core::pose::Pose & pose)
{
	check_and_initialize(pose);

	core::Size nstruct = jd2::JobDistributor::get_instance()->current_job()->nstruct_index();

	core::Size search_index = nstruct % search_points_.size();
	if( search_index == 0 )
	{
		search_index = search_points_.size();
	}
	
	core::kinematics::Stub transform = search_points_[search_index];

	core::Size residuejumpindex;
	core::Size residueindex;

	StubGenerator::placeResidueAtTransform(pose, target_residue_, transform, residuejumpindex, residueindex);

	core::pose::add_variant_type_to_pose_residue( pose, "SHOVE_BB", residueindex );
}

void PlaceProbeMover::check_and_initialize(core::pose::Pose const & target_pose)
{
	using namespace protocols::hotspot_hashing;

  if (initialized_pattern_)
  {
    return;
  }

	TR.Debug << "Initializing search pattern." << std::endl;

  SearchPatternOP search_pattern = create_partitioned_search_pattern(target_pose);

	TR.Info << "Initialized search pattern. Size: " << search_points_.size() << std::endl;

  protocols::jd2::JobOP current_job(jd2::JobDistributor::get_instance()->current_job());

  if (current_job->nstruct_max() < search_points_.size())
	{
		TR.Error << "Current job nstruct_max: " << current_job->nstruct_max() << " less than search pattern size: " << search_points_.size() << std::endl;
	}

  if (current_job->nstruct_max() < search_points_.size())
	{
		TR.Warning << "Current job nstruct_max: " << current_job->nstruct_max() << " greater than search pattern size: " << search_points_.size() << " Search points will be repeated" << std::endl;
	}
}

SearchPatternOP PlaceProbeMover::create_partitioned_search_pattern(core::pose::Pose const & target_pose)
{
	return new PartitionedSearchPattern(create_search_pattern(target_pose), search_partition_, total_search_partition_);
}

void
PlaceProbeMover::parse_place_probe_tag( utility::tag::TagPtr const tag,
                                protocols::moves::DataMap & data,
                                protocols::filters::Filters_map const &,
                                protocols::moves::Movers_map const &,
                                core::pose::Pose const &)
{
	// Residue spec
	if(tag->hasOption("residue_name"))
	{
		residue_name_ = tag->getOption<std::string>("residue_name");
	}
	else
	{
    utility_exit_with_message( "residue_name not specified" );
	}

  // Partition Spec
  search_partition_ = tag->getOption< core::Size >( "search_partition", 0 );
  total_search_partition_ = tag->getOption< core::Size >( "total_search_partition", 1 );

	if (!(search_partition_ >= 0 && search_partition_ < total_search_partition_ && total_search_partition_ > 0))
	{
		TR.Error << "Invalid search partition specficition. Partition: " << search_partition_ << " Total partitions: " << total_search_partition_ << std::endl;

		utility_exit_with_message("Invalid search partition specification.");
	}

	// Initialize residue representation
	target_residue_ = StubGenerator::getStubByName( residue_name_ );
}

}
}
}
