// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /devel/moves/MoverStub.hh
/// @brief
/// @author

#ifndef INCLUDED_devel_denovo_protein_design_FragmentSequenceMover_hh
#define INCLUDED_devel_denovo_protein_design_FragmentSequenceMover_hh

// Unit Headers
#include <devel/denovo_protein_design/FragmentSequenceMover.fwd.hh>

// Package Headers
#include <core/fragment/FragSet.hh>
#include <core/kinematics/MoveMap.hh>
// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/Mover.hh>

// Utility Headers
#include <core/types.hh>

namespace devel {
namespace denovo_protein_design {

///@details
class FragmentSequenceMover : public protocols::moves::Mover {

public:

	///@brief
  FragmentSequenceMover();

	///@brief
	FragmentSequenceMover( core::fragment::FragSetCOP fragset, core::kinematics::MoveMapCOP movemap );

	virtual ~FragmentSequenceMover();

	virtual void apply( core::pose::Pose & pose );

private:
	core::fragment::FragSetCOP fragset_;
	core::kinematics::MoveMapCOP movemap_;

};//end FragmentSequenceMover

}//namespace denovo_protein_design
}//namespace devel

#endif // INCLUDED_devel_DenovoProteinDesign_FragmentSequenceMover_HH
