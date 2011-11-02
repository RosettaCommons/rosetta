// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    devel/metal_interface/AddMetalSiteConstraints.hh
/// @brief   Adds metal binding site geometry constraints to pose.
/// @details I typically work with zinc binding sites having tetrahedral coordination geometry, in which case there are 4 distance constraints, 4 angle constraints, 6 tetrahedral constraints, and 4 dihedral constraints.  The code is flexibile enough to accommodate 2, 3, or 4-residue metal binding sites.
/// @author Bryan Der

#ifndef INCLUDED_devel_metal_interface_AddMetalSiteConstraints_HH
#define INCLUDED_devel_metal_interface_AddMetalSiteConstraints_HH

#include <devel/metal_interface/AddMetalSiteConstraints.fwd.hh>
#include <devel/metal_interface/MetalSiteResidue.fwd.hh>
// AUTO-REMOVED #include <core/pose/Pose.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>



using namespace core;

namespace devel {
namespace metal_interface {

///@brief Add metalsite geometry constraints to pose
class AddMetalSiteConstraints : public utility::pointer::ReferenceCount {

public:

  AddMetalSiteConstraints( utility::vector1< devel::metal_interface::MetalSiteResidueOP > msr );
  virtual ~AddMetalSiteConstraints();

	virtual void add_constraints  ( pose::Pose & pose );

private:
	utility::vector1< devel::metal_interface::MetalSiteResidueOP > msr_;

};//end AddMetalSiteConstraints


}//namespace metal_interface
}//namespace devel

#endif // INCLUDED_devel_metal_interface_AddMetalSiteConstraints_HH

