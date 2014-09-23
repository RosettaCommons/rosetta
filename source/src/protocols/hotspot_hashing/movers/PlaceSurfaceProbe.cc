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
/// @author Alex Ford fordas@uw.edu
//

#include <basic/Tracer.hh>

#include <utility/tag/Tag.hh>
#include <utility/exit.hh>

#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueProperty.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pack/task/TaskFactory.hh>

#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.fwd.hh>

#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

#include <protocols/rosetta_scripts/util.hh>

#include <protocols/hotspot_hashing/SearchPattern.hh>
#include <protocols/hotspot_hashing/SurfaceSearchPattern.hh>
#include <protocols/hotspot_hashing/SICSearchPattern.hh>

#include <protocols/hotspot_hashing/StubGenerator.hh>

#include <protocols/hotspot_hashing/movers/PlaceSurfaceProbe.hh>
#include <protocols/hotspot_hashing/movers/PlaceSurfaceProbeCreator.hh>


namespace protocols
{
namespace hotspot_hashing
{
namespace movers
{

static thread_local basic::Tracer TR( "protocols.hotspot_hashing.movers.PlaceSurfaceProbe" );

PlaceSurfaceProbe::PlaceSurfaceProbe() :
	  protocols::moves::Mover("PlaceSurfaceProbe"),
	  protocols::hotspot_hashing::movers::PlaceProbeMover(),
		surface_selection_(/* NULL */),
		search_density_(0),
		coarse_angle_sampling_(0),
		coarse_sampling_(0),
		refinement_distance_(0),
		refinement_angle_sampling_(0),
		refinement_sampling_(0),
		refinement_pattern_(/* NULL */)
{

}


PlaceSurfaceProbe::PlaceSurfaceProbe(
		std::string residue_name,
		core::Real search_density,
		core::Real coarse_angle_sampling,
		core::Real coarse_sampling,
		core::Real refinement_distance,
		core::Real refinement_angle_sampling,
		core::Real refinement_sampling,
		core::conformation::ResidueCOP target_residue,
		core::pack::task::TaskFactoryOP /*surface_selection*/,
		core::Size search_partition,
		core::Size total_search_partition) :
		protocols::moves::Mover( "PlaceSurfaceProbe" ),
		protocols::hotspot_hashing::movers::PlaceProbeMover(
			residue_name,
			target_residue,
			search_partition,
			total_search_partition),
		search_density_(search_density),
		coarse_angle_sampling_(coarse_angle_sampling),
		coarse_sampling_(coarse_sampling),
		refinement_distance_(refinement_distance),
		refinement_angle_sampling_(refinement_angle_sampling),
		refinement_sampling_(refinement_sampling),
		refinement_pattern_(initialize_refinement_pattern())
{
}

protocols::moves::MoverOP PlaceSurfaceProbe::clone() const
{
	return protocols::moves::MoverOP( new PlaceSurfaceProbe(*this) );
}

SearchPatternOP PlaceSurfaceProbe::create_search_pattern(core::pose::Pose const & target_pose)
{
	using core::conformation::ResidueOP;

	SearchPatternOP surface_pattern( new SurfaceSearchPattern(
							target_pose,
							surface_selection_,
							search_density_) );

	core::Real expected_course_search_bound = std::sqrt(search_density_);

	SearchPatternOP spherical_rotation_pattern( new SphericalRotationSearchPattern(
		coarse_angle_sampling_,
		coarse_angle_sampling_,
		coarse_angle_sampling_,
		0, 360,
		0, 90,
		0, 360) );

	SearchPatternOP cartesian_pattern( new CartesianSearchPattern(
		0,
		coarse_sampling_,
		coarse_sampling_,
		0, 0,
		-(expected_course_search_bound / 2),
		expected_course_search_bound / 2,
		-(expected_course_search_bound / 2),
		expected_course_search_bound / 2) );
		
	SearchPatternOP residue_sampling_pattern( new ComposeSearchPatterns(spherical_rotation_pattern, cartesian_pattern) );
	core::pose::Pose residue_pose;

	ResidueOP virtual_bb_residue = core::pose::add_variant_type_to_residue(*target_residue_, core::chemical::VIRTUAL_BB, target_pose);
	StubGenerator::placeResidueOnPose(residue_pose, virtual_bb_residue);

	SearchPatternOP sampled_surface_pattern( new SICPatternAtTransform(
				target_pose,
				residue_pose,
				surface_pattern,
				residue_sampling_pattern) );

	return sampled_surface_pattern;
}

SearchPatternOP PlaceSurfaceProbe::create_partitioned_search_pattern(core::pose::Pose const & target_pose)
{
	using core::conformation::ResidueOP;

	SearchPatternOP surface_pattern( new SurfaceSearchPattern(
							target_pose,
							surface_selection_,
							search_density_) );

	SearchPatternOP partitioned_surface_pattern( new PartitionedSearchPattern(surface_pattern, search_partition_, total_search_partition_) );

	SearchPatternOP residue_sampling_pattern( new SphericalRotationSearchPattern(
				coarse_angle_sampling_,
				coarse_angle_sampling_,
				coarse_angle_sampling_,
				0, 360,
				0, 180,
				0, 360) );

	core::pose::Pose residue_pose;

	ResidueOP virtual_bb_residue = core::pose::add_variant_type_to_residue(*target_residue_, core::chemical::VIRTUAL_BB, target_pose);
	StubGenerator::placeResidueOnPose(residue_pose, virtual_bb_residue);

	SearchPatternOP sampled_surface_pattern( new SICPatternAtTransform(
				target_pose,
				residue_pose,
				partitioned_surface_pattern,
				residue_sampling_pattern) );

	return sampled_surface_pattern;
}

SearchPatternOP PlaceSurfaceProbe::create_refinement_pattern(core::pose::Pose const & /*target_pose*/, core::Size /*target_residue*/)
{
	return refinement_pattern_;
}

SearchPatternOP PlaceSurfaceProbe::initialize_refinement_pattern()
{
	SearchPatternOP spherical_rotation_pattern( new SphericalRotationSearchPattern(
		refinement_angle_sampling_,
		refinement_angle_sampling_,
		refinement_angle_sampling_,
		0, coarse_angle_sampling_,
		0, coarse_angle_sampling_,
		0, coarse_angle_sampling_) );
		
	SearchPatternOP cartesian_pattern( new CartesianSearchPattern(
		refinement_sampling_,
		refinement_sampling_,
		refinement_sampling_,
		-refinement_distance_, refinement_distance_,
		-(coarse_sampling_ / 2),
		coarse_sampling_ / 2,
		-(coarse_sampling_ / 2),
		coarse_sampling_ / 2) );
		
	return SearchPatternOP( new ComposeSearchPatterns(spherical_rotation_pattern, cartesian_pattern) );
}

void
PlaceSurfaceProbe::parse_my_tag( utility::tag::TagCOP tag,
                                basic::datacache::DataMap & data,
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

	// Coarse grid spec
	coarse_angle_sampling_ = tag->getOption< core::Real >( "coarse_angle_sampling");
	coarse_sampling_ = tag->getOption< core::Real >( "coarse_sampling");

	refinement_angle_sampling_ = tag->getOption< core::Real >( "refinement_angle_sampling");
	refinement_sampling_ = tag->getOption< core::Real >( "refinement_sampling");
	refinement_distance_ = tag->getOption< core::Real >( "refinement_distance");

	refinement_pattern_ = initialize_refinement_pattern();
}

protocols::moves::MoverOP
PlaceSurfaceProbeCreator::create_mover() const
{
  return protocols::moves::MoverOP( new PlaceSurfaceProbe );
}

std::string
PlaceSurfaceProbeCreator::keyname() const
{
  return "PlaceSurfaceProbe";
}

std::string PlaceSurfaceProbeCreator::mover_name()
{
	return "PlaceSurfaceProbeMover";
}

}
}
}
