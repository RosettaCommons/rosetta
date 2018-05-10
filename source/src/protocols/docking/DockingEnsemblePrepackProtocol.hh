// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file DockingEnsemblePrepackProtocol
/// @brief Prepacking of the bound structure before docking with ensembles
/// @author Monica Berrondo
/// @author Modified by Jeliazko Jeliazkov -- check_ensemble_member_compatibility()

#ifndef INCLUDED_protocols_docking_DockingEnsemblePrepackProtocol_HH
#define INCLUDED_protocols_docking_DockingEnsemblePrepackProtocol_HH

// Unit Headers
#include <protocols/docking/DockingEnsemblePrepackProtocol.fwd.hh>
#include <protocols/docking/DockingHighRes.hh>

// Package headers
#include <protocols/docking/SidechainMinMover.fwd.hh>
#include <protocols/docking/DockingEnsemble.fwd.hh>

// Project headers
#include <protocols/moves/MoverContainer.fwd.hh>
#include <protocols/rigid/RigidBodyMover.fwd.hh>
#include <protocols/minimization_packing/PackRotamersMover.fwd.hh>
#include <protocols/minimization_packing/RotamerTrialsMinMover.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace docking {


class DockingEnsemblePrepackProtocol : public DockingHighRes {

public:
	/// @brief Default constructor
	DockingEnsemblePrepackProtocol();

	~DockingEnsemblePrepackProtocol() override;

	/// @brief Assigns default values to members
	void setup_defaults();

	/// @brief Instantiates and configures movers used by DockingEnsemblePrepackProtocol
	void setup_pack_operation_movers();

	/// @brief Returns a vector1 of the pdb chains in a pose, in order
	utility::vector1< char > get_pose_chains( core::pose::Pose & );

	/// @brief Ensures all members of either ensemble are compatible with one another
	void check_ensemble_member_compatibility();

	void move_away( core::pose::Pose & pose );
	void move_back( core::pose::Pose & pose );

	void apply( core::pose::Pose & ) override;

	std::string get_name() const override;

	void set_ensemble1( std::string ensemble1 ) { ensemble1_filename_ = ensemble1; }
	void set_ensemble2( std::string ensemble2 ) { ensemble2_filename_ = ensemble2; }

private:
	// add @brief for members
	utility::vector1< rigid::RigidBodyTransMoverOP > trans_away_vec_;
	utility::vector1< rigid::RigidBodyTransMoverOP > trans_back_vec_;

	core::Real trans_magnitude_;
	bool movers_setup_; //only append sequence mover once

	protocols::minimization_packing::RotamerTrialsMinMoverOP rtmin_mover_;
	protocols::minimization_packing::PackRotamersMoverOP prepack_full_repack_;
	SidechainMinMoverOP scmin_mover_;
	protocols::moves::SequenceMoverOP pack_operations_;

	// ensemble objects
	DockingEnsembleOP ensemble1_;
	DockingEnsembleOP ensemble2_;
	std::string ensemble1_filename_, ensemble2_filename_;

	/// @brief Performs setup that requires a pose
	void finalize_setup( core::pose::Pose & );
	void register_options();
	void init_from_options();
};

}
}
#endif
