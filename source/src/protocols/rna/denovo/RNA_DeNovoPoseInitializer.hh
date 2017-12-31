// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @author Rhiju Das

#ifndef INCLUDED_protocols_rna_denovo_RNA_DeNovoPoseInitializer_HH
#define INCLUDED_protocols_rna_denovo_RNA_DeNovoPoseInitializer_HH

#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/import_pose/RNA_JumpMover.fwd.hh>
#include <core/pose/rna/BasePairStep.hh>
#include <core/import_pose/RNA_DeNovoParameters.hh>
#include <core/pose/toolbox/AtomLevelDomainMap.fwd.hh>
#include <core/import_pose/libraries/RNA_ChunkLibrary.fwd.hh>
#include <core/types.hh>
#include <utility/pointer/ReferenceCount.hh>
// ObjexxFCL Headers

// C++ Headers
#include <string>
#include <map>
#include <iostream>
#include <list>

#include <utility/vector1.hh>

namespace protocols {
namespace rna {
namespace denovo {

//////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Object used in RNA_FragmentMonteCarlo to handle setup of pose & fold-tree, and
///   jump changes.
class RNA_DeNovoPoseInitializer : public utility::pointer::ReferenceCount {
public:

	//constructor
	RNA_DeNovoPoseInitializer( core::import_pose::RNA_DeNovoParameters const & rna_params_file_info );
	virtual ~RNA_DeNovoPoseInitializer();

	/// @brief "classic" setup, used in denovo protocol. Note: changes pose (virtualizes phosphate).
	///    and requires later call of setup_fold_tree_and_jumps_and_variants
	void
	initialize_for_de_novo_protocol(
		core::pose::Pose & pose,
		bool const ignore_secstruct = false );

	void
	setup_fold_tree_and_jumps_and_variants( core::pose::Pose & pose,
		core::import_pose::RNA_JumpMover const & rna_jump_mover,
		core::pose::toolbox::AtomLevelDomainMapOP atom_level_domain_map,
		core::import_pose::libraries::RNA_ChunkLibrary const & rna_chunk_library,
		bool const & enumerate = false ) const;

	void
	setup_fold_tree_and_jumps_and_variants( core::pose::Pose & pose ) const;

	void
	set_root_at_first_rigid_body( bool const setting ){ root_at_first_rigid_body_ = setting; }

	void
	set_dock_each_chunk( bool const & setting ){ dock_each_chunk_ = setting; }

	void
	set_dock_each_chunk_per_chain( bool const & setting ){ dock_each_chunk_per_chain_ = setting; }

	void
	set_center_jumps_in_single_stranded( bool const & setting ){ center_jumps_in_single_stranded_ = setting; }

	void
	set_new_fold_tree_initializer( bool const & setting ){ new_fold_tree_initializer_ = setting; }

	bool new_fold_tree_initializer() const { return new_fold_tree_initializer_; }

	void
	set_model_with_density( bool const & setting ){ model_with_density_ = setting; }

	bool model_with_density() const { return model_with_density_; }

	void set_bps_moves( bool const & setting ){ bps_moves_ = setting; }

	core::import_pose::RNA_DeNovoParameters const & rna_params() const { return rna_params_; }

	void
	setup_chainbreak_variants( core::pose::Pose & pose,
		core::pose::toolbox::AtomLevelDomainMapOP atom_level_domain_map ) const;


private:

	void
	initialize_secstruct( core::pose::Pose & pose );

	void
	override_secstruct( core::pose::Pose & pose );

	void setup_virtual_phosphate_variants( core::pose::Pose & pose ) const;

	void
	append_virtual_anchor( core::pose::Pose & pose );

	void
	setup_jumps( core::pose::Pose & pose,
		core::import_pose::RNA_JumpMover const & rna_jump_mover,
		core::import_pose::libraries::RNA_ChunkLibrary const & rna_chunk_library,
		bool const enumerate = false ) const;

	void
	setup_chainbreak_variants( core::pose::Pose & pose,
		core::pose::toolbox::AtomLevelDomainMapOP atom_level_domain_map,
		core::import_pose::RNA_JumpMover const & rna_jump_mover,
		core::import_pose::libraries::RNA_ChunkLibrary const & rna_chunk_library,
		bool const & enumerate = false ) const;


	void
	setup_fold_tree_through_build_full_model_info(
		core::pose::Pose & pose,
		core::import_pose::libraries::RNA_ChunkLibrary const & chunk_library,
		bool const & enumerate = false ) const;

	core::kinematics::FoldTree
	setup_fold_tree_legacy( core::pose::Pose & pose,
		core::import_pose::RNA_JumpMover const & rna_jump_mover ) const;

	void
	setup_block_stack_variants(
		core::pose::Pose & pose,
		core::pose::toolbox::AtomLevelDomainMapOP atom_level_domain_map ) const;

	void
	insert_base_pair_jumps( core::pose::Pose & pose, core::import_pose::RNA_JumpMover const & jump_mover,  bool & success ) const;

private:

	core::import_pose::RNA_DeNovoParameters rna_params_;

	bool const assume_non_stem_is_loop; // legacy parameter.
	bool bps_moves_;
	bool root_at_first_rigid_body_;
	bool dock_each_chunk_;
	bool dock_each_chunk_per_chain_;
	bool center_jumps_in_single_stranded_;
	bool new_fold_tree_initializer_;
	bool model_with_density_;

};


} //denovo
} //rna
} //protocols

#endif
