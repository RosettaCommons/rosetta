// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/topology_broker/TMHTopologySamplerClaimer.hh
/// @brief header file for TMHTopologySamplerClaimer, to be used by the TopologyBroker for de novo folding of TM proteins
/// @details header for sampling protocol that treats transmembrane helices as rigid bodies and moves them around to improve
///  sampling of membrane protein topologies
///
/// @author Stephanie H. DeLuca (stephanie.h.deluca@vanderbilt.edu)

#ifndef INCLUDED_protocols_topology_broker_TMHTopologySamplerClaimer_hh
#define INCLUDED_protocols_topology_broker_TMHTopologySamplerClaimer_hh

// Unit Headers
#include <protocols/topology_broker/TMHTopologySamplerClaimer.fwd.hh>

// Package Headers
#include <protocols/topology_broker/claims/DofClaim.fwd.hh>
#include <protocols/topology_broker/TopologyClaimer.hh>
#include <protocols/topology_broker/TopologyBroker.hh>

// Project Headers
#include <protocols/topology_broker/MembraneTopologyClaimer.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/scoring/MembraneTopology.hh>
#include <core/scoring/MembranePotential.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>

#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// C++ headers
#include <string>
#include <vector>

namespace protocols {
namespace topology_broker {

class TMHTopologySamplerClaimer : public TopologyClaimer{
	typedef TopologyClaimer Parent;

public:
	TMHTopologySamplerClaimer();

	/// @brief constructor: supply mover classes for Fragment Moves
	TMHTopologySamplerClaimer(topology_broker::TopologyBrokerOP broker);

	/// @brief destructor
	
	~TMHTopologySamplerClaimer() override;

	/// @brief TMHTopologySamplerClaimer has virtual functions... use this to obtain a new instance
	TopologyClaimerOP clone() const override {
		return TopologyClaimerOP( new TMHTopologySamplerClaimer( *this ) );
	}

	/// @brief type() is specifying the output name of the TopologyClaimer
	std::string type() const override {
		return _static_type_name();
	}

	static std::string _static_type_name() {
		return "TMHTopologySamplerClaimer";
	}

	/// @brief register cmd-line options in option system ( call before core::init )
	static void
	register_options();

	/// @brief get the pose from the boker and set it as this object's pose
	void set_pose_from_broker(core::pose::Pose& pose) override;

	/// @brief read in the pose's spans via the MembraneTopology store din the pose, determine jumps, TMHs, loops, etc.
	void pre_process(core::pose::Pose& pose) override;

	/// @brief the broker checks if the claimer builds its own fold tree to figure out if needs to build one itself
	bool claimer_builds_own_fold_tree() override;

	//getter function to retrieve this claimer's fold tree for the broker
	core::kinematics::FoldTreeOP get_fold_tree(core::pose::Pose& pose) override;

	//getter function to retrieve this claimer's current_pose (to get after moving spans)
	core::pose::PoseOP get_pose_from_claimer() override;

	/// @brief generate DofClaims
	void generate_claims( claims::DofClaims& dof_claims) override; //add to list ( never call clear() on list )

	/// @brief make move_map and add the DoF claims from generate_claims() to the movemap.  Now we can move certain parts with certain DOFs.
	void initialize_dofs( core::pose::Pose& pose, claims::DofClaims const& init_dofs, claims::DofClaims& failed_to_init) override;

	/// @brief claimers can add movers to the RandomMover (Container).
	/// add your moves, make it dependent on stage if you want to. So far this is called only by abinitio...
	/// if you don't want to do anything special --- don't overload this method!
	/// default: adds mover given by virtual call get_mover()  with stage-dependent weight given by abinitio_mover_weight_
	void add_mover(
		moves::RandomMover& /* random_mover */,
		core::pose::Pose const& /*pose*/,
		abinitio::StageID /*stageID*/, /* abinitio sampler stage */
		core::scoring::ScoreFunction const& /*scorefxn*/,
		core::Real /* progress */ /* progress within stage */
	) override;

	/// @brief this claimer builds its own radial fold tree based on read-in spanfile
	void build_fold_tree(core::pose::Pose& pose, core::kinematics::FoldTree& fold_tree) override;

	/// @brief read tag from topology broker file (setup.tpb)
	bool read_tag( std::string tag, std::istream& is ) override;

protected:
	/// @brief called by constructor ---  calls all set_default_XXX methods
	void set_defaults() override;

	/// @brief Make extended chain an idealized helix
	void set_pose_torsions(core::pose::Pose& pose);

	/// @brief move helices closer together
	void move_spans(core::pose::Pose& pose);

	/// @brief output membrane center and normal as virtual atoms
	core::Vector output_membrane_vector(core::pose::Pose& pose);

	/// @brief pre-compute grid points where you want to move the helices
	utility::vector1<core::Vector> pre_compute_grid_points(core::pose::Pose& pose);

	/// @brief get membrane topology information
	core::scoring::MembraneTopologyOP get_membrane_topology(core::pose::Pose& pose);

	/// @brief get membrane embedding information
	core::scoring::MembraneEmbed get_membrane_embed(core::pose::Pose& pose);

private:
	core::Real rotation_mag_;
	core::Real translation_mag_;
	core::Real rb_mover_stage1_weight_;
	core::pose::Pose current_pose_;
	int nres_;
	core::Size njumps_;
	ObjexxFCL::FArray2D_int jump_array_;
	ObjexxFCL::FArray1D_int cuts_;
	int topology_root_res_;
	core::Size tmhelix_;
};
} // namespace topology_broker
} // namespace protocols

#endif // TMHTopologySamplerClaimer.hh
