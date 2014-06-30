// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/environment/claims/EnvLabelSelector.cc
/// @brief  The EnvLabelSelector holds a selection that another object can set inside of it.
/// @author Justin R. Porter

// Unit headers
#include <protocols/environment/claims/EnvLabelSelector.hh>

// Package headers
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>

#include <core/environment/LocalPosition.hh>
#include <protocols/environment/ProtectedConformation.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>

#include <boost/foreach.hpp>

// C++ headers
#include <cassert>

namespace protocols {
namespace environment {
namespace claims {

EnvLabelSelector::EnvLabelSelector() {}

EnvLabelSelector::EnvLabelSelector( LocalPositions const& positions_in ) {
  this->set_local_positions( positions_in );
}

EnvLabelSelector::~EnvLabelSelector() {}

void EnvLabelSelector::apply(
  core::pose::Pose const & pose,
  ResidueSubset& subset
) const
{
  using core::environment::LocalPositionOP;
  assert( subset.size() == pose.total_residue() );

  for( Size i = 1; i <= subset.size(); ++i ){
    subset[i] = false;
  }

  ProtectedConformationCOP conf = dynamic_cast< ProtectedConformation const* >( &( pose.conformation() ) );
  core::environment::SequenceAnnotationCOP ann = conf->annotations();

  BOOST_FOREACH( LocalPositionOP pos, positions ){
    subset[ ann->resolve_seq( *pos ) ] = true;
  }
}

void EnvLabelSelector::parse_my_tag(
  utility::tag::TagCOP,
  basic::datacache::DataMap & )
{
  throw utility::excn::EXCN_RosettaScriptsOption( "Not to be used in RosettaScripts. For legacy compatibility only." );
}

void EnvLabelSelector::set_local_positions( LocalPositions const& positions_in ){
  using namespace core::environment;

  BOOST_FOREACH( LocalPositionOP pos, positions_in ){
    positions.push_back( new LocalPosition( *pos ) );
  }
}

std::string EnvLabelSelector::get_name() const {
  return EnvLabelSelector::class_name();
}

std::string EnvLabelSelector::class_name() {
  return "EnvLabel";
}

} //namespace claims
} //namespace environment
} //namespace protocols

