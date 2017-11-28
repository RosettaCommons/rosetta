// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/rna/denovo/movers/RNA_DeNovoMasterMover.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_rna_denovo_movers_RNA_DeNovoMasterMover_HH
#define INCLUDED_protocols_rna_denovo_movers_RNA_DeNovoMasterMover_HH

#include <protocols/moves/Mover.hh>
#include <protocols/rna/denovo/options/RNA_FragmentMonteCarloOptions.fwd.hh>
#include <protocols/rna/denovo/movers/RNA_DeNovoMasterMover.fwd.hh>
#include <protocols/rna/denovo/movers/RNA_FragmentMover.fwd.hh>
#include <protocols/rna/denovo/fragments/RNA_Fragments.fwd.hh>
#include <protocols/rna/denovo/base_pairs/RNA_BasePairHandler.fwd.hh>
#include <protocols/rna/denovo/libraries/RNA_ChunkLibrary.fwd.hh>
#include <protocols/rna/denovo/movers/RNA_JumpMover.fwd.hh>
#include <protocols/rna/denovo/movers/RNA_HelixMover.fwd.hh>
#include <protocols/rna/movers/RNA_LoopCloser.fwd.hh>
#include <protocols/rigid/RigidBodyMover.fwd.hh>
#include <protocols/toolbox/AtomLevelDomainMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/types.hh>

namespace protocols {
namespace rna {
namespace denovo {
namespace movers {

class RNA_DeNovoMasterMover: public protocols::moves::Mover {

public:

	//constructor
	RNA_DeNovoMasterMover( options::RNA_FragmentMonteCarloOptionsCOP options,
		protocols::toolbox::AtomLevelDomainMapCOP atom_level_domain_map,
		base_pairs::RNA_BasePairHandlerCOP rna_base_pair_handler,
		protocols::rna::movers::RNA_LoopCloserOP rna_loop_closer,
		libraries::RNA_ChunkLibraryOP rna_chunk_library  );

	//destructor
	~RNA_DeNovoMasterMover();

public:

	virtual void apply( core::pose::Pose & pose );

	virtual std::string get_name() const{ return "RNA_DeNovoMasterMover"; }

	void
	apply( core::pose::Pose & pose,
		core::Size const & cycle_number );

	void
	setup_rnp_fold_tree( core::pose::Pose & pose, bool const & refine_pose, bool const & randomize );

	void
	setup_rigid_body_mover( core::pose::Pose const & pose,
		core::Real const & rot_mag,
		core::Real const & trans_mag );

	void
	setup_dock_into_density_mover( core::pose::Pose const & pose,
		core::Real const & rot_mag,
		core::Real const & trans_mag );

	void
	setup_rna_protein_docking_mover( core::pose::Pose const & pose,
		core::Real const & rot_mag,
		core::Real const & trans_mag );

	movers::RNA_FragmentMoverOP rna_fragment_mover(){ return rna_fragment_mover_; }
	movers::RNA_JumpMoverOP rna_jump_mover(){ return rna_jump_mover_; }

	void
	do_random_moves( core::pose::Pose & pose, Size const & monte_carlo_cycles, bool const & check_num_rna_res = false );

	void set_close_loops( bool const & setting ){ close_loops_ = setting; }
	bool close_loops() const { return close_loops_; }

	void set_frag_size( core::Size const & setting ){ frag_size_ = setting; }
	core::Size frag_size() const { return frag_size_; }

	std::string move_type() const { return move_type_; }
	bool success() const { return ( move_type_.size() > 0 ); }

	void
	search_rigid_body_orientation( core::pose::Pose & pose );

	void
	set_rna_helix_mover( movers::RNA_HelixMoverOP const rna_helix_mover ){ rna_helix_mover_ = rna_helix_mover; }

	movers::RNA_HelixMoverOP rna_helix_mover(){ return rna_helix_mover_; }

	void
	set_helix_mover_magnitude( core::Real const & rot_mag, core::Real const & trans_mag );

private:

	void
	do_move_trial( core::Size const & i, core::pose::Pose & pose );

	void
	RNA_move_trial( core::pose::Pose & pose );

	void
	dock_into_density_trial( core::pose::Pose & pose );

	void
	random_jump_trial( core::pose::Pose & pose );

	void
	rnp_docking_trial( core::pose::Pose & pose );

	void
	random_fragment_trial( core::pose::Pose & pose );

	bool
	random_chunk_trial( core::pose::Pose & pose );

	void
	randomize_rigid_body_orientations( core::pose::Pose & pose );

	void
	randomize_rnp_rigid_body_orientations( core::pose::Pose & pose );


	// core::kinematics::FoldTree
	// get_rnp_docking_fold_tree( core::pose::Pose const & pose );

private:

	options::RNA_FragmentMonteCarloOptionsCOP options_;

	core::Size frag_size_;
	core::Real jump_change_frequency_;
	core::Real dock_into_density_freq_;
	bool close_loops_;
	bool do_rnp_docking_;
	std::string move_type_;

	movers::RNA_FragmentMoverOP rna_fragment_mover_;
	movers::RNA_JumpMoverOP rna_jump_mover_;
	libraries::RNA_ChunkLibraryOP rna_chunk_library_;
	protocols::rna::movers::RNA_LoopCloserOP rna_loop_closer_;
	protocols::rigid::RigidBodyPerturbMoverOP rigid_body_mover_, rnp_docking_mover_,
		dock_into_density_mover_;
	core::kinematics::FoldTree rnp_docking_ft_;
	movers::RNA_HelixMoverOP rna_helix_mover_;

};

} //movers
} //denovo
} //rna
} //protocols

#endif
