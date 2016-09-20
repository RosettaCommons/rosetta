// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file DockingPrepackProtocol
/// @brief Prepacking of the bound structure before
///        docking
/// @author Robin A Thottungal (raugust1@jhu.edu)
///   added to: JKLeman (julia.koehler1982@gmail.com)
#ifndef INCLUDED_protocols_docking_DockingPrepackProtocol_hh
#define INCLUDED_protocols_docking_DockingPrepackProtocol_hh

// Unit Headers
#include <protocols/docking/DockingPrepackProtocol.fwd.hh>
#include <protocols/docking/DockingHighRes.hh>
#include <protocols/docking/DockingInitialPerturbation.fwd.hh>

// Package headers
#include <protocols/docking/SidechainMinMover.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

#include <protocols/moves/MoverContainer.fwd.hh>
#include <protocols/rigid/RigidBodyMover.fwd.hh>
#include <protocols/simple_moves/PackRotamersMover.fwd.hh>
#include <protocols/simple_moves/RotamerTrialsMinMover.fwd.hh>

// Utility headers
#include <core/types.hh>
#include <numeric/xyzVector.fwd.hh>
#include <utility/vector1.fwd.hh>


namespace protocols {
namespace docking {

class DockingPrepackProtocol : public DockingHighRes {

public:
	/// @brief Default constructor
	DockingPrepackProtocol();

	~DockingPrepackProtocol() override;

	/// @brief Assigns default values to members
	void setup_defaults();

	/// @brief Instantiates and configures movers used by DockingPrepackProtocol
	void setup_pack_operation_movers();


	void apply( core::pose::Pose & ) override;

	std::string get_name() const override;

	/// @brief Scores and outputs the pose - jd2 compatible.
	void score_and_output(std::string filename,core::pose::Pose &);
	void set_dock_ppk(bool dock_ppk);

private:
	// add @brief for members
	utility::vector1< rigid::RigidBodyTransMoverOP > trans_away_vec_;
	utility::vector1< rigid::RigidBodyTransMoverOP > trans_back_vec_;

	core::Real trans_magnitude_;

	/// @brief membrane for translating in the membrane plane
	bool membrane_;

	protocols::simple_moves::RotamerTrialsMinMoverOP rtmin_mover_;
	protocols::simple_moves::PackRotamersMoverOP prepack_full_repack_;
	SidechainMinMoverOP scmin_mover_;
	protocols::moves::SequenceMoverOP pack_operations_;
	bool dock_ppk_;
	FaDockingSlideIntoContactOP fa_dock_slide_into_contact_;
	bool movers_setup_;

	/// @brief Performs setup that requires a pose
	void finalize_setup( core::pose::Pose & );
	void register_options();
	void init_from_options();

};

}
}
#endif
