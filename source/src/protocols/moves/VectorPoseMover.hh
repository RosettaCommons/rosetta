// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/moves/VectorPoseMover.hh
/// @brief Designates a mover that can be passed multiple poses by the VectorPoseJobDistributor
/// Any movers deriving from this subclass can then act on all of the input poses simultaneously
/// Only accessible through recon application.
/// @author Alex Sevy (alex.sevy@gmail.com)

#ifndef INCLUDED_PROTOCOLS_MOVES_VectorPoseMover_HH
#define INCLUDED_PROTOCOLS_MOVES_VectorPoseMover_HH

#include <protocols/moves/VectorPoseMover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace moves {

/// @brief Designates a mover that can be passed multiple poses by the VectorPoseJobDistributor
/// Any movers deriving from this subclass can then act on all of the input poses simultaneously
/// Only accessible through recon application.
class VectorPoseMover : public Mover {

public:
	/// @brief Constructor
	VectorPoseMover();

	~VectorPoseMover() override;

	VectorPoseMover( std::string const& name );

	VectorPoseMover( VectorPoseMover const & other );

	std::string get_name() const override;

	void apply ( core::pose::Pose& pose ) override = 0;

	/// @brief Pure virtual method to apply this mover under MPI
	virtual void apply_mpi ( core::pose::Pose& pose ) = 0;

	/// @brief Set the vector of poses for the mover to act upon
	void set_poses( utility::vector1< core::pose::PoseOP > const & poses );

protected:
	/// @brief Vector of all the input poses
	utility::vector1< core::pose::PoseOP > poses_;
};

}  // namespace moves
}  // namespace protocols


#endif  // PROTOCOLS_MOVES_VectorPoseMover_HH_
