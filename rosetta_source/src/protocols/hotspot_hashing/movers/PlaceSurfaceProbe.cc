// Project Headers
// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 sw=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/hotspot_hashing/movers/PlaceSurfaceProbe.cc
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
#include <protocols/hotspot_hashing/SICSearchPattern.hh>
#include <protocols/hotspot_hashing/PlaceMinimizeSearch.hh>

#include <protocols/hotspot_hashing/StubGenerator.hh>

#include <protocols/hotspot_hashing/movers/PlaceSurfaceProbe.hh>
#include <protocols/hotspot_hashing/movers/PlaceSurfaceProbeCreator.hh>


namespace protocols
{
namespace hotspot_hashing
{
namespace movers
{

static basic::Tracer TR( "protocols.hotspot_hashing.movers.PlaceSurfaceProbe" ); 

PlaceSurfaceProbe::PlaceSurfaceProbe() :
	  protocols::moves::Mover("PlaceSurfaceProbe"),
	  protocols::hotspot_hashing::movers::PlaceProbeMover(),
		search_density_(0),
		surface_selection_(NULL),
		x_angle_sampling_(0),
		y_angle_sampling_(0),
		refinement_distance_sampling_(0),
		refinement_distance_(0),
		refinement_translation_sampling_(0),
		refinement_pattern_(NULL)
{
 
} 


PlaceSurfaceProbe::PlaceSurfaceProbe(
		std::string residue_name,
		core::Real search_density,
		core::Real x_angle_sampling,
		core::Real y_angle_sampling,
		core::Real refinement_distance_sampling,
		core::Real refinement_distance,
		core::Real refinement_translation_sampling,
    core::conformation::ResidueCOP target_residue,
		core::pack::task::TaskFactoryOP surface_selection,
    core::Size search_partition,
    core::Size total_search_partition) :
    protocols::moves::Mover( "PlaceSurfaceProbe" ),
	  protocols::hotspot_hashing::movers::PlaceProbeMover(
			residue_name,
			target_residue,
			search_partition,
			total_search_partition),
		search_density_(search_density),
		surface_selection_(surface_selection),
		x_angle_sampling_(x_angle_sampling),
		y_angle_sampling_(y_angle_sampling),
		refinement_distance_sampling_(refinement_distance_sampling),
		refinement_distance_(refinement_distance),
		refinement_translation_sampling_(refinement_translation_sampling),
		refinement_pattern_(initialize_refinement_pattern())
{
}

protocols::moves::MoverOP PlaceSurfaceProbe::clone() const
{
	return new PlaceSurfaceProbe(*this);
}

SearchPatternOP PlaceSurfaceProbe::create_search_pattern(core::pose::Pose const & target_pose)
{
	using core::conformation::ResidueOP;

	SearchPatternOP surface_pattern(
			new SurfaceSearchPattern(
							target_pose,
							surface_selection_,
							search_density_));

	SearchPatternOP residue_sampling_pattern(
			new RotationSearchPattern(x_angle_sampling_, y_angle_sampling_));

	core::pose::Pose residue_pose;

	ResidueOP virtual_bb_residue = core::pose::add_variant_type_to_residue(*target_residue_, "VIRTUAL_BB", target_pose);
	StubGenerator::placeResidueOnPose(residue_pose, virtual_bb_residue);

	SearchPatternOP sampled_surface_pattern(
			new SICPatternAtTransform(
				target_pose,
				residue_pose,
				surface_pattern,
				residue_sampling_pattern));

	return sampled_surface_pattern;
}

SearchPatternOP PlaceSurfaceProbe::create_partitioned_search_pattern(core::pose::Pose const & target_pose)
{
	using core::conformation::ResidueOP;

	SearchPatternOP surface_pattern(
			new SurfaceSearchPattern(
							target_pose,
							surface_selection_,
							search_density_));

	SearchPatternOP partitioned_surface_pattern(
			new PartitionedSearchPattern(surface_pattern, search_partition_, total_search_partition_));

	SearchPatternOP residue_sampling_pattern(
			new RotationSearchPattern(x_angle_sampling_, y_angle_sampling_));

	core::pose::Pose residue_pose;

	ResidueOP virtual_bb_residue = core::pose::add_variant_type_to_residue(*target_residue_, "VIRTUAL_BB", target_pose);
	StubGenerator::placeResidueOnPose(residue_pose, virtual_bb_residue);

	SearchPatternOP sampled_surface_pattern(
			new SICPatternAtTransform(
				target_pose,
				residue_pose,
				partitioned_surface_pattern,
				residue_sampling_pattern));

	return sampled_surface_pattern;
}

SearchPatternOP PlaceSurfaceProbe::create_refinement_pattern(core::pose::Pose const & target_pose, core::Size target_residue)
{
	return refinement_pattern_;
}

SearchPatternOP PlaceSurfaceProbe::initialize_refinement_pattern()
{
	core::Real expected_bound = std::sqrt(search_density_);

	return new CartesianSearchPattern(
				refinement_distance_sampling_,
				refinement_translation_sampling_,
				refinement_translation_sampling_,
				-refinement_distance_,
				refinement_distance_,
				-(expected_bound / 2),
				expected_bound / 2,
				-(expected_bound / 2),
				expected_bound / 2
			);
}

void
PlaceSurfaceProbe::parse_my_tag( utility::tag::TagPtr const tag,
                                protocols::moves::DataMap & data,
                                protocols::filters::Filters_map const & filters_map,
                                protocols::moves::Movers_map const & movers_map,
                                core::pose::Pose const & target_pose)
{
	parse_place_probe_tag(
		tag,
		data,
		filters_map,
		movers_map,
		target_pose
	);

	// Surface Spec
	search_density_ = tag->getOption< core::Real >( "search_density", 1);
	surface_selection_ = protocols::rosetta_scripts::parse_task_operations( tag, data );

	// Sampling Spec
	x_angle_sampling_ = tag->getOption< core::Real >( "x_angle_sampling", 5 );
	y_angle_sampling_ = tag->getOption< core::Real >( "y_angle_sampling", 5 );

	refinement_distance_sampling_ = tag->getOption< core::Real >( "refinement_distance_sampling", .05 );
	refinement_distance_ = tag->getOption< core::Real >( "refinement_distance", 1 );
	refinement_translation_sampling_ = tag->getOption< core::Real >( "refinement_translation_sampling", .05 );

	refinement_pattern_ = initialize_refinement_pattern();
}

protocols::moves::MoverOP
PlaceSurfaceProbeCreator::create_mover() const
{
  return new PlaceSurfaceProbe;
}

std::string
PlaceSurfaceProbeCreator::keyname() const
{
  return "PlaceSurfaceProbe";
}

}
}
}
