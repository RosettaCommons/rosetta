// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_sic_dock_RigidScore_FWD_hh
#define INCLUDED_protocols_sic_dock_RigidScore_FWD_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace sic_dock {

class RigidScore;
typedef utility::pointer::shared_ptr< RigidScore > RigidScoreOP;
typedef utility::pointer::shared_ptr< RigidScore const > RigidScoreCOP;

class CBScore;
typedef utility::pointer::shared_ptr< CBScore > CBScoreOP;
typedef utility::pointer::shared_ptr< CBScore const > CBScoreCOP;

class LinkerScore;
typedef utility::pointer::shared_ptr< LinkerScore > LinkerScoreOP;
typedef utility::pointer::shared_ptr< LinkerScore const > LinkerScoreCOP;

class ConstraintSetScore;
typedef utility::pointer::shared_ptr< ConstraintSetScore > ConstraintSetScoreOP;
typedef utility::pointer::shared_ptr< ConstraintSetScore const > ConstraintSetScoreCOP;

// class EdgeStandScore;
// typedef utility::pointer::owning_ptr< EdgeStandScore > EdgeStandScoreOP;
// typedef utility::pointer::owning_ptr< EdgeStandScore const > EdgeStandScoreCOP;

// class HelixScore;
// typedef utility::pointer::owning_ptr< HelixScore > HelixScoreOP;
// typedef utility::pointer::owning_ptr< HelixScore const > HelixScoreCOP;

// class BuriedPolarScore;
// typedef utility::pointer::owning_ptr< BuriedPolarScore > BuriedPolarScoreOP;
// typedef utility::pointer::owning_ptr< BuriedPolarScore const > BuriedPolarScoreCOP;

class JointScore;
typedef utility::pointer::shared_ptr< JointScore > JointScoreOP;
typedef utility::pointer::shared_ptr< JointScore const > JointScoreCOP;

// class CachedScore;
// typedef utility::pointer::owning_ptr< CachedScore > CachedScoreOP;
// typedef utility::pointer::owning_ptr< CachedScore const > CachedScoreCOP;

}
}

#endif
