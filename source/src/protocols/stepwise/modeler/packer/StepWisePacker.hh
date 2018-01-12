// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file SWA_Packer.hh
/// @brief
/// @details
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_stepwise_protein_StepWisePacker_HH
#define INCLUDED_protocols_stepwise_protein_StepWisePacker_HH

#include <protocols/stepwise/modeler/packer/StepWisePacker.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/minimization_packing/GreenPacker.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <utility/vector1.hh>
#include <string>

namespace protocols {
namespace stepwise {
namespace modeler {
namespace packer {

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
class StepWisePacker: public protocols::moves::Mover {
public:

	//constructor!
	StepWisePacker( utility::vector1< core::Size > const & working_moving_residues );

	//destructor -- necessary?
	~StepWisePacker();

	/// @brief Apply the minimizer to one pose
	virtual void apply( core::pose::Pose & pose );

	virtual std::string get_name() const;

	void
	do_packing( core::pose::Pose & pose );

	void
	reset( core::pose::Pose const & pose );

	void
	set_scorefxn( core::scoring::ScoreFunctionCOP scorefxn );

	void
	set_use_packer_instead_of_rotamer_trials( bool const & setting ){ use_packer_instead_of_rotamer_trials_ = setting; }

	void
	do_prepack( core::pose::Pose & pose );

	void set_allow_virtual_side_chains( bool const & setting ){ allow_virtual_side_chains_ = setting; }
	bool allow_virtual_side_chains() const{ return allow_virtual_side_chains_; }

	void set_allow_virtual_o2prime_hydrogens( bool const & setting ){ allow_virtual_o2prime_hydrogens_ = setting; }
	bool allow_virtual_o2prime_hydrogens() const{ return allow_virtual_o2prime_hydrogens_; }

	void set_pack_all_side_chains( bool const & setting ){ pack_all_side_chains_ = setting; }
	bool pack_all_side_chains() const{ return pack_all_side_chains_; }

	void set_pack_o2prime_hydrogens( bool const & setting ){ pack_o2prime_hydrogens_ = setting; }
	bool pack_o2prime_hydrogens() const { return pack_o2prime_hydrogens_; }

	bool working_pack_res_was_inputted() const { return working_pack_res_was_inputted_; }

	void set_working_pack_res( utility::vector1< core::Size > const & setting ) { working_pack_res_ = setting; }
	utility::vector1< core::Size > const & working_pack_res() const { return working_pack_res_; }
	utility::vector1< core::Size > const & previous_working_pack_res() const { return previous_working_pack_res_; }

	void
	figure_out_neighbors( core::pose::Pose & pose );

private:

	void
	reinstate_side_chain_angles( core::pose::Pose & pose, core::pose::Pose const & src_pose );

	void
	setup_pack_task( core::pose::Pose const & pose );

private:

	utility::vector1< core::Size > const working_moving_res_list_;
	core::scoring::ScoreFunctionCOP scorefxn_;
	bool use_packer_instead_of_rotamer_trials_;
	bool allow_virtual_side_chains_;
	bool allow_virtual_o2prime_hydrogens_;
	bool pack_o2prime_hydrogens_;
	core::pack::task::PackerTaskOP pack_task_;
	core::pose::PoseOP pose_init_;

	utility::vector1< core::Size > working_pack_res_; // set this from outside, or packer will figure it out from working_moving_res
	utility::vector1< core::Size > previous_working_pack_res_; // cached of residues packed last time, used for reinstating side chains.
	bool working_pack_res_was_inputted_;
	bool pack_all_side_chains_;

};

} //packer
} //modeler
} //stepwise
} //protocols

#endif
