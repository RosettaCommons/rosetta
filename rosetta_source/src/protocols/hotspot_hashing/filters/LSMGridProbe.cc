// Project Headers
// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 sw=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/hotspot_hashing/filters/LSMGridProbe.cc
/// @brief
/// @author Alex Ford fordas@uw.edu

#include <sstream>

#include <basic/Tracer.hh>

#include <utility/tag/Tag.hh>
#include <utility/exit.hh>

#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/pose/Pose.fwd.hh>

#include <protocols/filters/Filter.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/DataMap.fwd.hh>

#include <protocols/hotspot_hashing/filters/LSMGridProbe.hh>
#include <protocols/hotspot_hashing/filters/LSMGridProbeFilterCreator.hh>
#include <protocols/hotspot_hashing/PlaceMinimizeSearch.hh>
	

namespace protocols
{
namespace hotspot_hashing
{
namespace filters
{

static basic::Tracer TR( "protocols.hotspot_hashing.filters.LSMGridProbe" );

LSMGridProbe::LSMGridProbe() :
    parent( "LSMGridProbe" ),
    residue_name_(""),
    lsm_spec_(),
		angle_sampling_(0),
		translocation_sampling_(0),
		max_radius_(0),
		distance_sampling_(0),
		max_distance_(0),
    relax_mover_(NULL),
    triage_filter_(NULL)
{
}

void
LSMGridProbe::parse_my_tag( utility::tag::TagPtr const tag,
                                protocols::moves::DataMap &,
                                protocols::filters::Filters_map const &filters,
                                protocols::moves::Movers_map const & movers,
                                core::pose::Pose const & )
{
	// Triage filter
  std::string const triage_filter_name( tag->getOption< std::string >( "triage_filter", "true_filter" ) );

	if (filters.count(triage_filter_name) == 0)
	{
    utility_exit_with_message( "Triage filter " + triage_filter_name + " not found" );
	}
	else
	{
		triage_filter_ = filters.find( triage_filter_name )->second;
	}

	// Relax mover
  std::string relax_mover_name( tag->getOption< std::string >( "relax_mover", "null" ) );

	if (movers.count(relax_mover_name) == 0)
	{
    utility_exit_with_message( "relax mover " + relax_mover_name + " not found" );
	}
	else
	{
		relax_mover_ = movers.find( relax_mover_name )->second;
	}

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

  TR<<"<LSMGridProbe " <<
		"triage_filter=\"" << triage_filter_name << "\" " <<
		"relax_mover=\"" << relax_mover_name << "\" " <<
		"residue_name=\"" << residue_name_ << "\" " <<
		"lsm_spec=\"" << lsm_spec_ << "\" " <<
 	  "lsm_id=" << lsm_id_ << " "<< 
 	  "angle_sampling=" << angle_sampling_ << " "<< 
 	  "translocation_sampling=" << translocation_sampling_ << " "<< 
 	  "max_radius=" << max_radius_ << " "<< 
 	  "distance_sampling=" << distance_sampling_ << " "<< 
 	  "max_distance=" << max_distance_ << " "<< 
		">"<< std::endl;
}

bool LSMGridProbe::apply( core::pose::Pose const & pose) const
{
	SearchPatternOP pattern_ref(new LSMSearchPattern(
			lsm_spec_, 
			angle_sampling_,
			translocation_sampling_,
			max_radius_,
			distance_sampling_,
			max_distance_));

	// Initialize residue representation
	core::chemical::ResidueTypeSet const & residue_set ( pose.residue(1).residue_type_set() );
	core::chemical::ResidueType const & restype( residue_set.name_map( residue_name_ ) );
	core::conformation::ResidueCOP residue(core::conformation::ResidueFactory::create_residue( restype ));

	std::stringstream output_tag (std::stringstream::in | std::stringstream::out);
	output_tag << residue_name_ << "." << lsm_id_;

  PlaceMinimizeSearch search(pose, residue, pattern_ref, relax_mover_, triage_filter_, output_tag.str());

  search.execute();

  return true;
}

protocols::filters::FilterOP
LSMGridProbe::fresh_instance() const
{
  return new LSMGridProbe();
}

LSMGridProbe::~LSMGridProbe() {}

protocols::filters::FilterOP
LSMGridProbe::clone() const
{
  return new LSMGridProbe( *this );
}

protocols::filters::FilterOP
LSMGridProbeFilterCreator::create_filter() const
{
  return new LSMGridProbe;
}

std::string
LSMGridProbeFilterCreator::keyname() const
{
  return "LSMGridProbe";
}


}
}
}
