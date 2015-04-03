// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file edge reweighting
/// @brief collection of routines to assign different weights to IG edges
/// @author Florian Richter, floric@u.washington.edu, june 08

// Unit headers
#include <protocols/toolbox/IGEdgeReweighters.hh>

// Package headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/pack/task/PackerTask.hh>

#include <utility/vector1.hh>


namespace protocols{
namespace toolbox{

using namespace core;

core::Real
IGLigandDesignEdgeUpweighter::get_edge_reweight(
  pose::Pose const & pose,
  pack::task::PackerTask const & task,
  Size res1,
  Size res2
) const
{

  if( ( pose.residue( res1 ).is_ligand() && task.design_residue( res2 ) )
    ||( pose.residue( res2 ).is_ligand() && task.design_residue( res1 ) ) ){
    return weight_factor_;
  }

  else return default_weight_;

}

/*
template <class T>
core::Real
ResidueGroupIGEdgeUpweighter::get_edge_reweight(
  pose::Pose const & pose,
  pack::task::PackerTask const & task,
  Size res1,
  Size res2
) const
{

	if( (group1_.find(res1) != group1_.end()) && (group2_.find(res2) != group2_.end() )
		||(group2_.find(res1) != group2_.end()) && (group1_.find(res2) != group1_.end() ) ){
		return weight_factor_;
	}
	else return default_weight_;

}
*/

} //namespace toolbox
} //namespace protocols
