// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /protocols/moves/ReturnSidechainMover.hh
/// @brief Mover to "return" proper sidechains to a pose that was temporarily in centroid mode (can be used for any sidechain copying)
/// @author Steven Lewis

#ifndef INCLUDED_protocols_moves_ReturnSidechainMover_hh
#define INCLUDED_protocols_moves_ReturnSidechainMover_hh

// Unit Headers
#include <protocols/moves/ReturnSidechainMover.fwd.hh>

// Project Headers
#include <core/pose/Pose.hh> //we're going to contain a pose
#include <protocols/moves/Mover.hh>

// Utility Headers
#include <core/types.hh>
#include <utility/vector1.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace moves {

///@details This class takes two poses, one at instantiation and one at apply.  It copies the residue type set and chi information from its kept pose to the apply pose.  The intended purpose is for "returning" known sidechains to a pose that was temporarily in centroid mode, although it may work with other residue type sets.
class ReturnSidechainMover : public protocols::moves::Mover {

public:

	///@brief default constructor
	ReturnSidechainMover();

	///@brief constructor with pose
	ReturnSidechainMover(
		core::pose::Pose const & pose_in,
		core::Size start_res = 0,
		core::Size end_res = 0 );

	///@brief constructor with pose and copiable residue array
	ReturnSidechainMover(
		core::pose::Pose const & pose_in,
		utility::vector1<bool> allow_chi_in,
		core::Size start_res = 0,
		core::Size end_res = 0 );

	virtual ~ReturnSidechainMover();

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

	bool copy_all_chi_;
	utility::vector1<bool> allow_chi_copy_;

private:
	///@brief remembered old pose
	core::pose::Pose const remembered_pose_;

	///@brief residue numbers for which residues to loop over for recovery
	core::Size start_res_, end_res_;

};//end ReturnSidechainMover

}//namespace moves
}//namespace protocols

#endif // INCLUDED_protocols_moves_ReturnSidechainMover_HH
