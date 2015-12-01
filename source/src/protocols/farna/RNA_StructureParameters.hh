// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @author Rhiju Das

#ifndef INCLUDED_protocols_rna_RNA_StructureParameters_HH
#define INCLUDED_protocols_rna_RNA_StructureParameters_HH

#include <core/pose/Pose.fwd.hh>
#include <core/pose/rna/BasePair.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <protocols/farna/RNA_JumpLibrary.fwd.hh>
#include <protocols/farna/BasePairStep.hh>
#include <protocols/toolbox/AllowInsert.fwd.hh>
#include <core/kinematics/Jump.hh>
#include <core/types.hh>
#include <core/id/AtomID.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
// ObjexxFCL Headers

// C++ Headers
#include <string>
#include <map>
#include <iostream>
#include <list>

#include <utility/vector1.hh>


namespace protocols {
namespace farna {

typedef utility::vector1< core::pose::rna::BasePair > RNA_BasePairList;

//////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Object used in RNA_FragmentMonteCarlo to handle setup of pose & fold-tree, and
///   jump changes.
class RNA_StructureParameters : public utility::pointer::ReferenceCount {
public:

	//constructor
	RNA_StructureParameters();
	virtual ~RNA_StructureParameters();

	/// @brief if pose already has desired fold-tree, cutpoints, etc., use this!
	void
	initialize_from_pose( core::pose::Pose const & pose );

	/// @brief "classic" setup, used in denovo protocol. Note: changes pose (virtualizes phosphate).
	///    and requires later call of setup_fold_tree_and_jumps_and_variants
	void
	initialize_for_de_novo_protocol(
		core::pose::Pose & pose,
		std::string const rna_params_file,
		std::string const jump_library_file,
		bool const ignore_secstruct );

	void
	setup_fold_tree_and_jumps_and_variants( core::pose::Pose & pose ) const;

	bool
	random_jump_change( core::pose::Pose & pose ) const;

	toolbox::AllowInsertOP
	allow_insert() const;

	void
	set_allow_insert(toolbox::AllowInsertOP allow_insert );

	void
	set_root_at_first_rigid_body( bool const setting ){ root_at_first_rigid_body_ = setting; }

	void
	set_suppress_bp_constraint( core::Real const setting ){ suppress_bp_constraint_ = setting; }

	void
	set_bps_moves( Size const setting ){ bps_moves_ = setting; }

	bool
	check_base_pairs( core::pose::Pose & pose ) const;

	std::map< Size, Size >
	connections() const;

	std::list< Size >
	get_stem_residues(  core::pose::Pose const & pose ) const;

	void
	setup_base_pair_constraints( core::pose::Pose & pose ) const;

	void
	setup_virtual_phosphate_variants( core::pose::Pose & pose );

	utility::vector1 < utility::vector1 <core::Size > >
	get_obligate_pairing_sets() {
		return obligate_pairing_sets_;
	}

	RNA_BasePairList
	get_rna_pairing_list() {
		return rna_pairing_list_;
	}

	utility::vector1< BasePairStep >
	get_canonical_base_pair_steps() const;

	utility::vector1< BasePairStep >
	get_noncanonical_base_pair_steps() const;

	void
	set_jump_library( RNA_JumpLibraryCOP rna_jump_library );

private:

	void
	initialize_secstruct( core::pose::Pose & pose );

	void
	override_secstruct( core::pose::Pose & pose );

	void update_allow_insert_to_move_internal_phosphates( core::pose::Pose const & pose );

	void add_virtual_phosphate_variants( core::pose::Pose & pose ) const;

	void update_allow_insert_to_not_move_virtual_phosphates( core::pose::Pose const & pose );

	void
	append_virtual_anchor( core::pose::Pose & pose );

	void
	initialize_allow_insert( core::pose::Pose const & pose  );

	void
	get_pairings_from_line(
		std::istringstream & line_stream,
		bool const in_stem );

	void
	save_res_lists_to_chain_connections_and_clear( utility::vector1< Size > & res_list1,
		utility::vector1< Size > & res_list2 );

	void
	read_chain_connection( std::istringstream & line_stream );

	void
	read_parameters_from_file( std::string const & pairing_file );

	void
	setup_jumps( core::pose::Pose & pose ) const;

	void
	setup_chainbreak_variants( core::pose::Pose & pose ) const;

	std::string const
	read_secstruct_from_file( std::string const & rna_secstruct_file );

	core::Size
	check_in_chain_connections( core::Size const & pos1, core::Size const & pos2 ) const;

	bool
	check_forward_backward(
		core::pose::Pose & pose,
		core::Size const jump_pos ) const;

	void
	add_new_RNA_jump(
		core::pose::Pose & pose,
		core::Size const & which_jump,
		bool & success ) const;

	void
	sample_alternative_chain_connection( core::pose::Pose & pose, core::Size const & which_jump ) const;

	void
	insert_base_pair_jumps( core::pose::Pose & pose, bool & success ) const;

	void
	fill_in_default_jump_atoms( core::kinematics::FoldTree & f, core::pose::Pose const & pose ) const;

	void
	figure_out_partner( std::map< Size, Size > & partner, bool const force_canonical ) const;

	Size
	check_in_pairing_sets( utility::vector1 < utility::vector1 <core::Size > > pairing_sets,
												 core::pose::rna::BasePair const & rna_pairing_check ) const;

	utility::vector1< BasePairStep >
	get_base_pair_steps( bool const just_canonical ) const;

private:

	RNA_JumpLibraryCOP rna_jump_library_;
	RNA_BasePairList rna_pairing_list_;

	utility::vector1 < utility::vector1 <core::Size > > obligate_pairing_sets_;
	utility::vector1 < utility::vector1 <core::Size > > stem_pairing_sets_;

	utility::vector1 < std::pair< utility::vector1 <core::Size >, utility::vector1 <core::Size > > > chain_connections_;

	//  int force_stems_;// deprecated.

	utility::vector1 <core::Size > cutpoints_open_;
	utility::vector1 <core::Size > cutpoints_closed_;
	utility::vector1 <core::Size > virtual_anchor_attachment_points_;

	bool secstruct_defined_;
	std::string rna_secstruct_;
	bool assume_non_stem_is_loop;
	bool bps_moves_;

	bool add_virtual_anchor_;
	bool root_at_first_rigid_body_;
	core::Real suppress_bp_constraint_;

	utility::vector1 < core::Size  > allow_insert_res_;
	toolbox::AllowInsertOP allow_insert_;


};

typedef utility::pointer::shared_ptr< RNA_StructureParameters > RNA_StructureParametersOP;


} //farna
} //protocols

#endif
