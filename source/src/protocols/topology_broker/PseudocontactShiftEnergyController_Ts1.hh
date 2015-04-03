// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

 //////////////////////////////////////////////
 /// @file PseudocontactShiftEnergyController.hh
 ///
 /// @authorv Christophe Schmitz //kalabharath & Oliver Lange
 ///
 ////////////////////////////////////////////////

#ifndef INCLUDED_protocols_topology_broker_PseudocontactShiftEnergyController_Ts1_hh
#define INCLUDED_protocols_topology_broker_PseudocontactShiftEnergyController_Ts1_hh

// Unit Headers
#include <protocols/topology_broker/PseudocontactShiftEnergyController_Ts1.fwd.hh>

// Package Headers
#include <protocols/topology_broker/TopologyClaimer.hh>
#include <protocols/topology_broker/claims/DofClaim.fwd.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/moves/MoverContainer.fwd.hh>


// ObjexxFCL Headers

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

//Auto Headers
#include <core/scoring/constraints/Constraint.hh>


//// C++ headers

namespace protocols {
namespace topology_broker {

class PseudocontactShiftEnergyController_Ts1 : public TopologyClaimer {

	typedef TopologyClaimer Parent;
public:
                 	PseudocontactShiftEnergyController_Ts1(); //for factory

	virtual TopologyClaimerOP clone() const {
		return TopologyClaimerOP( new PseudocontactShiftEnergyController_Ts1( *this ) );
	}

	/// @brief type() is specifying the output name of the TopologyClaimer
	virtual std::string type() const {
		return _static_type_name();
	}

	static std::string _static_type_name() {
		return "PseudocontactShiftEnergyController_Ts1";
	}

	virtual bool read_tag( std::string tag, std::istream & );

	virtual void set_defaults(); //eg before reading starts.


	virtual void add_mover(
    moves::RandomMover& /* random_mover */,
		core::pose::Pose const& /*pose*/,
		abinitio::StageID /*stageID*/, /* abinitio sampler stage */
		core::scoring::ScoreFunction const& /*scorefxn*/,
		core::Real /* progress */ /* progress within stage */
	);

	virtual void init_after_reading();

private:
	void
	control_grid_param();

	core::Real grid_edge_stage1_;
	core::Real grid_edge_stage2_;
	core::Real grid_edge_stage3_;
	core::Real grid_edge_stage4_;

	core::Real grid_step_stage1_;
	core::Real grid_step_stage2_;
	core::Real grid_step_stage3_;
	core::Real grid_step_stage4_;

	core::Real grid_small_cutoff_stage1_;
	core::Real grid_small_cutoff_stage2_;
	core::Real grid_small_cutoff_stage3_;
	core::Real grid_small_cutoff_stage4_;

	core::Real grid_large_cutoff_stage1_;
	core::Real grid_large_cutoff_stage2_;
	core::Real grid_large_cutoff_stage3_;
	core::Real grid_large_cutoff_stage4_;

	core::Real grid_cone_angle_cutoff_stage1_;
	core::Real grid_cone_angle_cutoff_stage2_;
	core::Real grid_cone_angle_cutoff_stage3_;
	core::Real grid_cone_angle_cutoff_stage4_;

	std::string grid_atom_name_1_stage1_;
	std::string grid_atom_name_1_stage2_;
	std::string grid_atom_name_1_stage3_;
	std::string grid_atom_name_1_stage4_;

	std::string grid_atom_name_2_stage1_;
	std::string grid_atom_name_2_stage2_;
	std::string grid_atom_name_2_stage3_;
	std::string grid_atom_name_2_stage4_;

	core::SSize grid_residue_num_1_stage1_;
	core::SSize grid_residue_num_1_stage2_;
	core::SSize grid_residue_num_1_stage3_;
	core::SSize grid_residue_num_1_stage4_;

	core::SSize grid_residue_num_2_stage1_;
	core::SSize grid_residue_num_2_stage2_;
	core::SSize grid_residue_num_2_stage3_;
	core::SSize grid_residue_num_2_stage4_;

	core::Real grid_k_vector_stage1_;
	core::Real grid_k_vector_stage2_;
	core::Real grid_k_vector_stage3_;
	core::Real grid_k_vector_stage4_;

	bool minimize_best_tensor_stage1_;
	bool minimize_best_tensor_stage2_;
	bool minimize_best_tensor_stage3_;
	bool minimize_best_tensor_stage4_;

	core::Real pcs_weight_stage1_;
	core::Real pcs_weight_stage2_;
	core::Real pcs_weight_stage3_;
	core::Real pcs_weight_stage4_;

	utility::vector1<std::string> filenames_;
	utility::vector1<core::Real> individual_weights_;

}; //class PseudocontactShiftEnergyController_Ts1

}
}

#endif
