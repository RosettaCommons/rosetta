// Project Headers
// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 sw=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/hotspot_hashing/movers/PlaceLSMProbe.cc
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

#include <protocols/moves/Mover.hh>
#include <protocols/moves/DataMap.fwd.hh>

#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

#include <protocols/hotspot_hashing/SearchPattern.hh>
#include <protocols/hotspot_hashing/PlaceMinimizeSearch.hh>

#include <protocols/hotspot_hashing/movers/PlaceLSMProbe.hh>
#include <protocols/hotspot_hashing/movers/PlaceLSMProbeCreator.hh>


namespace protocols
{
namespace hotspot_hashing
{
namespace movers
{

static basic::Tracer TR( "protocols.hotspot_hashing.movers.PlaceLSMProbe.cc" ); 

PlaceLSMProbe::PlaceLSMProbe() :
    parent( "PlaceLSMProbe" ),
    residue_name_(""),
    lsm_spec_(),
    lsm_id_(0),
		angle_sampling_(0),
		translocation_sampling_(0),
		max_radius_(0),
		distance_sampling_(0),
		max_distance_(0),
		log_pose_(false),
    target_residue_(NULL),
    initialized_pattern_(false),
    search_points_()
{
 
} 


PlaceLSMProbe::PlaceLSMProbe(
		std::string residue_name,
    VectorPair lsm_spec,
		core::Size lsm_id,
		core::Real angle_sampling,
		core::Real translocation_sampling,
		core::Real max_radius,
		core::Real distance_sampling,
		core::Real max_distance,
    bool log_pose,
    core::conformation::ResidueCOP target_residue) :
    parent( "PlaceLSMProbe" ),
    residue_name_(residue_name),
    lsm_spec_(lsm_spec),
    lsm_id_(lsm_id),
		angle_sampling_(angle_sampling),
		translocation_sampling_(translocation_sampling),
		max_radius_(max_radius),
		distance_sampling_(distance_sampling),
		max_distance_(max_distance),
		log_pose_(log_pose),
    target_residue_(target_residue),
    initialized_pattern_(false),
    search_points_()
{
}

protocols::moves::MoverOP PlaceLSMProbe::clone() const
{
	return new PlaceLSMProbe(*this);
}

void PlaceLSMProbe::apply(core::pose::Pose & pose)
{
	check_and_initialize();

	core::Size nstruct = jd2::JobDistributor::get_instance()->current_job()->nstruct_index();
	
  TransformPair transform = search_points_[nstruct];

	core::Size residuejumpindex;
	core::Size residueindex;

	PlaceMinimizeSearch::placeResidueAtTransform(pose, *target_residue_, transform, residuejumpindex, residueindex);

	core::pose::add_variant_type_to_pose_residue( pose, "SHOVE_BB", residueindex );
}

void PlaceLSMProbe::check_and_initialize()
{
  if (initialized_pattern_)
  {
    return;
  }

  LSMSearchPattern search_pattern(
						lsm_spec_, 
						angle_sampling_,
						translocation_sampling_,
						max_radius_,
						distance_sampling_,
						max_distance_);

	search_points_ = search_pattern.Searchpoints();

	TR.Info << "Initialized search pattern. Size: " << search_points_.size() << "\n";

  protocols::jd2::JobOP current_job(jd2::JobDistributor::get_instance()->current_job());

  if (current_job->nstruct_max() != search_points_.size())
	{
		TR.Warning << "Current job nstruct_max: " << current_job->nstruct_max() << " not equal to search pattern size: " << search_points_.size() << std::endl;
	}
}

void
PlaceLSMProbe::parse_my_tag( utility::tag::TagPtr const tag,
                                protocols::moves::DataMap &,
                                protocols::filters::Filters_map const &,
                                protocols::moves::Movers_map const &,
                                core::pose::Pose const & pose)
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
	
	// LSM Spec
	if(tag->hasOption("lsm_id"))
	{
		lsm_id_ = tag->getOption<core::Size>("lsm_id");
	}
	else
	{
    utility_exit_with_message( "lsm_id not specified" );
	}

	if(tag->hasOption("lsm_spec"))
	{
		// Parses spec string into ref.
		if(!LSMSearchPattern::parse_lsm_spec(tag->getOption<std::string>("lsm_spec"), lsm_spec_))
		{
			utility_exit_with_message( "unable to parse lsm_spec \"" + tag->getOption<std::string>("lsm_spec") + "\"" );
		}
	}
	else
	{
    utility_exit_with_message( "lsm_spec not specified" );
	}

	angle_sampling_ = tag->getOption< core::Real >( "angle_sampling", 90 );
	translocation_sampling_ = tag->getOption< core::Real >( "translocation_sampling", 1 );
	max_radius_ = tag->getOption< core::Real >( "max_radius", 2 );
	distance_sampling_ = tag->getOption< core::Real >( "distance_sampling", 1 );
	max_distance_ = tag->getOption< core::Real >( "max_distance", 2 );

	log_pose_ = tag->getOption<bool>("log_pose", false);

	// Initialize residue representation
	core::chemical::ResidueTypeSet const & residue_set ( pose.residue(1).residue_type_set() );
	core::chemical::ResidueType const & restype( residue_set.name_map( residue_name_ ) );
	target_residue_ = core::conformation::ResidueFactory::create_residue( restype );

  TR<<"<PlaceLSMProbe " <<
		"residue_name=\"" << residue_name_ << "\" " <<
		"lsm_spec=\"" << lsm_spec_ << "\" " <<
 	  "lsm_id=" << lsm_id_ << " "<< 
 	  "angle_sampling=" << angle_sampling_ << " "<< 
 	  "translocation_sampling=" << translocation_sampling_ << " "<< 
 	  "max_radius=" << max_radius_ << " "<< 
 	  "distance_sampling=" << distance_sampling_ << " "<< 
 	  "max_distance=" << max_distance_ << " "<< 
 	  "log_pose=" << log_pose_ << " "<< 
		">"<< std::endl;
}

protocols::moves::MoverOP
PlaceLSMProbeCreator::create_mover() const
{
  return new PlaceLSMProbe;
}

std::string
PlaceLSMProbeCreator::keyname() const
{
  return "PlaceLSMProbe";
}

}
}
}
