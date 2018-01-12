// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file SurfaceDocking.hh
/// @brief <add a description of the class>
/// @author Robin A Thottungal (rathottungal@gmail.com)
/// @author Michael Pacella (mpacella88@gmail.com)

#ifndef INCLUDED_protocols_surface_docking_SurfaceDockingProtocol_hh
#define INCLUDED_protocols_surface_docking_SurfaceDockingProtocol_hh

// Unit Headers
#include <protocols/moves/Mover.hh>
#include <protocols/surface_docking/SurfaceDockingProtocol.fwd.hh>

// Package headers
#include <protocols/surface_docking/SurfaceParameters.fwd.hh>
#include <protocols/surface_docking/SurfaceOrientMover.fwd.hh>
#include <protocols/surface_docking/CentroidRelaxMover.fwd.hh>
#include <protocols/surface_docking/FullatomRelaxMover.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.fwd.hh>
#include <protocols/abinitio/ClassicAbinitio.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <protocols/minimization_packing/PackRotamersMover.fwd.hh>
#include <protocols/docking/DockingInitialPerturbation.fwd.hh>
#include <protocols/rigid/RigidBodyMover.fwd.hh>

// C++ Headers
#include <string>
#include <sstream>

namespace protocols {
namespace surface_docking {

class SurfaceDockingProtocol : public moves::Mover {

public:
	//Standard Methods /////////////////////////////////////////////
	/// @brief Default constructor
	SurfaceDockingProtocol();

	/// @brief Copy constructor
	SurfaceDockingProtocol(SurfaceDockingProtocol const & src);

	//// @brief Assignment operator
	// Undefined, commenting out to fix PyRosetta build  SurfaceDockingProtocol & operator=(SurfaceDockingProtocol const & src);

	//destructor
	~SurfaceDockingProtocol() override;

	//Standard Rosetta methods /////////////////////////////////////
	//General methods
	/// @brief Register options with the option system.
	// Undefined, commenting out to fix PyRosetta build  static void register_options();

	/// @brief Generate string representation of SurfaceDockingProtocol for debugging purposes
	void show(std::ostream & ouput=std::cout) const override;

	/// Insertion operator (overloaded so that SurfaceDockingProtocol can be "printed" in Pyrosetta).
	//friend std::ostream & operator<<(std::ostream & output, SurfaceDockingProtocol const & object_to_output);

	/// Assignment operator


	// Mover methods
	/// @brief Return the name of the Mover.
	std::string get_name() const override;

	protocols::moves::MoverOP clone() const override;

	protocols::moves::MoverOP fresh_instance() const override;

	/// @brief Apply the corresponding move to the pose
	void apply( core::pose::Pose & pose) override;

private:
	// Private methods //////////////////////////////////////////

	//Initialize data members
	void init();

	// Copy all data members src to destinatioin
	void copy_data(SurfaceDockingProtocol object_to_copy_to, SurfaceDockingProtocol object_to_copy_from);

	bool valid_surface_pose(core::pose::Pose const & pose);

	void calc_secondary_structure(core::pose::Pose & pose);

	void calc_secondary_structure_with_surface(core::pose::Pose const & pose);

	void set_secondary_structure(core::pose::Pose & pose);

	void initialize_surface_energies (core::pose::Pose & pose, Size first_protein_residue);

	void set_surface_parameters ( core::pose::Pose & pose );

	void setup_movers ( core::pose::Pose const & pose, Size const first_protein_residue );

	void setup_abinitio ();

	void setup_slide_movers( core::pose::Pose const & pose );

	void split_protein_surface_poses (core::pose::Pose const & pose, core::pose::Pose & surface, core::pose::Pose & protein );

	void merge_protein_surface_poses (core::pose::Pose & pose, core::pose::Pose const & surface, core::pose::Pose const & protein );

	core::pack::task::PackerTaskOP create_surface_packer_task ( core::pose::Pose const & pose, Size const first_protein_residue );

	// Private data /////////////////////////////////////////////

	core::scoring::ScoreFunctionOP score_sidechain_pack_;
	std::string sec_struct_;
	protocols::surface_docking::SurfaceParametersOP surface_parameters_;
	simple_moves::SwitchResidueTypeSetMoverOP to_centroid_;
	simple_moves::SwitchResidueTypeSetMoverOP to_full_atom_;
	SurfaceOrientMoverOP surface_orient_;
	protocols::abinitio::ClassicAbinitioOP abinitio_;
	protocols::surface_docking::CentroidRelaxMoverOP centroid_relax_;
	protocols::minimization_packing::PackRotamersMoverOP pack_rotamers_fullatom_;
	protocols::rigid::RigidBodyTransMoverOP slide_away_from_surface_;
	protocols::docking::FaDockingSlideIntoContactOP slide_into_surface_;
	protocols::surface_docking::FullatomRelaxMoverOP fullatom_relax_;
	protocols::rigid::RigidBodyTransMoverOP position_above_surface_;


};


} // surface_docking
} // protocols

#endif  // INCLUDED_protocols_SurfaceDockingProtocol_hh
