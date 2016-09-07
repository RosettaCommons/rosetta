// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file PoseMembraneRigidBodyMover.hh
/// @brief
/// @author Yifan Song

#ifndef INCLUDED_protocols_rigid_PoseMembraneRigidBodyMover_hh
#define INCLUDED_protocols_rigid_PoseMembraneRigidBodyMover_hh

// Package headers
#include <protocols/moves/Mover.hh>


// Utility Headers
#include <core/types.hh>
#include <numeric/xyzVector.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace rigid {

/// move pose into a membrane
class MovePoseToMembraneCenterMover : public moves::Mover {

public:
	//default constructor
	MovePoseToMembraneCenterMover();

	void apply( core::pose::Pose & pose ) override;

	core::Vector estimate_membrane_center( core::pose::Pose & pose );
	std::string  get_name() const override;
};

/// perturb the pose along membrane normal
class MembraneCenterPerturbationMover : public moves::Mover {
public:
	MembraneCenterPerturbationMover (); // default constructor
	MembraneCenterPerturbationMover (core::Real const & trans_mag_in); // constructor

	void apply( core::pose::Pose & pose ) override;

	std::string  get_name() const override;
private:
	///data
	core::Real trans_mag_; //maximum translation magnitude

};

/// rotation pose around membrane center, perturb the membrane normal vector relative to the pose
class MembraneNormalPerturbationMover : public moves::Mover {
public:
	MembraneNormalPerturbationMover (); // default constructor
	MembraneNormalPerturbationMover ( core::Real const & rotation_mag_in ); // constructor

	void apply( core::pose::Pose & pose ) override;

	std::string  get_name() const override;
private:
	///data
	core::Real rotation_mag_; //maximum translation magnitude

};

/// translate the whole pose
class WholeBodyTranslationMover :  public moves::Mover {
public:
	WholeBodyTranslationMover ( core::Vector const & trans_in );

	void apply( core::pose::Pose & pose ) override;
	std::string  get_name() const override;

private:
	/// data
	core::Vector trans_;
};

/// rotate the whole pose
class WholeBodyRotationMover :  public moves::Mover {
public:
	WholeBodyRotationMover ( core::Vector const & axis, core::Vector const & center, core::Real const & alpha /* degrees */ );

	void apply( core::pose::Pose & pose ) override;
	std::string  get_name() const override;

private:
	/// data
	//numeric::xyzMatrix< core::Real > rotation_m_;
	core::Vector axis_;
	core::Vector center_;
	core::Real alpha_; /* degrees */
};

}
}

#endif //INCLUDED_protocols_rigid_PoseMembraneRigidBodyMover_HH
