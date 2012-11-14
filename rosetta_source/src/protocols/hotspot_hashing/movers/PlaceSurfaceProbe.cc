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
		angle_sampling_(0),
		translocation_sampling_(0),
		max_radius_(0),
		distance_sampling_(0),
		max_distance_(0)
{
 
} 


PlaceSurfaceProbe::PlaceSurfaceProbe(
		std::string residue_name,
		core::Real search_density,
		core::Real angle_sampling,
		core::Real translocation_sampling,
		core::Real max_radius,
		core::Real distance_sampling,
		core::Real max_distance,
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
		angle_sampling_(angle_sampling),
		translocation_sampling_(translocation_sampling),
		max_radius_(max_radius),
		distance_sampling_(distance_sampling),
		max_distance_(max_distance)
{
}

protocols::moves::MoverOP PlaceSurfaceProbe::clone() const
{
	return new PlaceSurfaceProbe(*this);
}

SearchPatternOP PlaceSurfaceProbe::create_search_pattern(core::pose::Pose const & target_pose)
{
	return new SurfaceSearchPattern(
							target_pose,
							surface_selection_,
							search_density_,
							angle_sampling_,
							translocation_sampling_,
							max_radius_,
							distance_sampling_,
							max_distance_);
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
	angle_sampling_ = tag->getOption< core::Real >( "angle_sampling", 15 );
	translocation_sampling_ = tag->getOption< core::Real >( "translocation_sampling", 1 );
	max_radius_ = tag->getOption< core::Real >( "max_radius", 0 );
	distance_sampling_ = tag->getOption< core::Real >( "distance_sampling", 1 );
	max_distance_ = tag->getOption< core::Real >( "max_distance", 0 );

  // Partition Spec
  search_partition_ = tag->getOption< core::Size >( "search_partition", 0 );
  total_search_partition_ = tag->getOption< core::Size >( "total_search_partition", 1 );

  TR<<"<PlaceSurfaceProbe " <<
		"residue_name=\"" << residue_name_ << "\" " <<
		"search_partition=\"" << search_partition_ << "\" " <<
		"total_search_partition=\"" << total_search_partition_ << "\" " <<
 	  "search_density=" << search_density_ << " "<< 
 	  "angle_sampling=" << angle_sampling_ << " "<< 
 	  "translocation_sampling=" << translocation_sampling_ << " "<< 
 	  "max_radius=" << max_radius_ << " "<< 
 	  "distance_sampling=" << distance_sampling_ << " "<< 
 	  "max_distance=" << max_distance_ << " "<< 
		">"<< std::endl;
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
