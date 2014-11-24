// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

///@file protocols/simple_moves/VectorPoseMover.hh
///@brief Designates a mover that can be passed multiple poses by the MSDJobDistributor
///@brief Any movers deriving from this subclass can then act on all of the input poses simultaneously
///@author Alex Sevy (alex.sevy@gmail.com)

#ifndef INCLUDED_PROTOCOLS_MOVES_VectorPoseMover_HH
#define INCLUDED_PROTOCOLS_MOVES_VectorPoseMover_HH

#include <protocols/moves/VectorPoseMover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace moves {

/// @detail A simple class used for a mover that acts on a vector of poses
class VectorPoseMover : public Mover {

public:
	/// @brief Constructor
	VectorPoseMover();

	virtual ~VectorPoseMover();

	VectorPoseMover( std::string const& name );

	VectorPoseMover( VectorPoseMover const & other );

	virtual std::string get_name() const;

	virtual void apply ( core::pose::Pose& pose ) = 0;

	/// @brief Set the vector of poses for the mover to act upon
	void set_poses( utility::vector1< core::pose::PoseOP > const & poses );

	/// @brief Sets current pose in case you want to act on only one of the present poses
//	void set_current_pose( core::Size current );

protected:
	utility::vector1< core::pose::PoseOP > poses_;
//	core::Size current_pose_;
};

}  // namespace moves
}  // namespace protocols



#endif  // PROTOCOLS_MOVES_VectorPoseMover_HH_
