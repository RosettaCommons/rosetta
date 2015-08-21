// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/metal_interface/MatchGrafter.hh
/// @brief  Takes a scaffold protein and a match pdb from RosettaMatch, grafts the match onto the protein.  For zinc homodimer design, it can then combine two grafted poses by overlaying the zinc atoms.
/// @author Bryan Der

#ifndef INCLUDED_protocols_metal_interface_MatchGrafter_HH
#define INCLUDED_protocols_metal_interface_MatchGrafter_HH

#include <protocols/metal_interface/MatchGrafter.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>


namespace protocols {
namespace metal_interface {


/// @details
class MatchGrafter : public utility::pointer::ReferenceCount {

public:

	typedef core::pose::Pose Pose;

	/// @brief
	MatchGrafter();

	virtual ~MatchGrafter();

	virtual Pose graft( Pose & match, Pose & partner_ungrafted );
	virtual Pose build_combined_pose_with_zinc_overlay( Pose & partner1, Pose & partner2 );
	virtual void ensure_proper_his_tautomers( Pose & combined, utility::vector1< core::Size > metalsite_seqpos );

private:


};//end MatchGrafter


}//namespace metal_interface
}//namespace protocols

#endif // INCLUDED_protocols_metal_interface_MatchGrafter_HH
