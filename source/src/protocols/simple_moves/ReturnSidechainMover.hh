// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /protocols/simple_moves/ReturnSidechainMover.hh
/// @brief protocols::moves::Mover to "return" proper sidechains to a pose that was temporarily in centroid mode (can be
/// used for any sidechain copying)
/// @author Steven Lewis

#ifndef INCLUDED_protocols_simple_moves_ReturnSidechainMover_hh
#define INCLUDED_protocols_simple_moves_ReturnSidechainMover_hh

// Unit Headers
#include <protocols/simple_moves/ReturnSidechainMover.fwd.hh>

// Project Headers
#include <core/pose/Pose.hh> //we're going to contain a pose
#include <protocols/moves/Mover.hh>

// Utility Headers
#include <core/types.hh>
#include <utility/vector1.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace simple_moves {

/// @details This class takes two poses, one at instantiation and one at apply.
/// It copies the residue type set and chi information from its kept pose to the apply pose.
/// The intended purpose is for "returning" known sidechains to a pose that was temporarily in centroid mode, although
/// it may work with other residue type sets.
class ReturnSidechainMover : public protocols::moves::Mover {

public:

	/// @brief default constructor
	ReturnSidechainMover();

	/// @brief constructor with pose
	ReturnSidechainMover(
		core::pose::Pose const & pose_in,
		core::Size start_res = 0,
		core::Size end_res = 0 );

	/// @brief constructor with pose and copiable residue array
	ReturnSidechainMover(
		core::pose::Pose const & pose_in,
		utility::vector1<bool> allow_chi_in,
		core::Size start_res = 0,
		core::Size end_res = 0 );

	/// @brief copy constructor
	ReturnSidechainMover(ReturnSidechainMover const & object_to_copy);

	virtual ~ReturnSidechainMover();

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;
	virtual void show(std::ostream & output=std::cout) const;
	virtual protocols::moves::MoverOP clone() const;
	virtual protocols::moves::MoverOP fresh_instance() const;

	bool copy_all_chi_;
	utility::vector1<bool> allow_chi_copy_;

	core::Size get_start_res() const;
	core::Size get_end_res() const;

private:
	/// @brief remembered old pose
	core::pose::Pose const remembered_pose_;

	/// @brief residue numbers for which residues to loop over for recovery
	core::Size start_res_, end_res_;

};//end ReturnSidechainMover

std::ostream &operator<< (std::ostream &os, ReturnSidechainMover const &mover);

}//namespace moves
}//namespace protocols

#endif // INCLUDED_protocols_simple_moves_ReturnSidechainMover_HH
