// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

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
