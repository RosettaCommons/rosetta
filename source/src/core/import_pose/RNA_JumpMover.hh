// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/import_pose/RNA_JumpMover.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_core_import_pose_RNA_JumpMover_HH
#define INCLUDED_core_import_pose_RNA_JumpMover_HH

#include <utility/pointer/ReferenceCount.hh>
#include <core/import_pose/RNA_JumpMover.fwd.hh>
#include <core/import_pose/libraries/RNA_JumpLibrary.fwd.hh>
#include <core/pose/toolbox/AtomLevelDomainMap.hh>
#include <core/types.hh>
#include <core/pose/rna/BasePair.hh>
#include <utility/vector1.hh>

namespace core {
namespace import_pose {

class RNA_JumpMover: public utility::pointer::ReferenceCount {

public:

	//constructor
	RNA_JumpMover( core::import_pose::libraries::RNA_JumpLibraryCOP rna_jump_library,
		core::pose::toolbox::AtomLevelDomainMapCOP atom_level_domain_map );

	//destructor
	~RNA_JumpMover();

public:

	bool
	random_jump_change( core::pose::Pose & pose ) const;


	void
	add_new_RNA_jump(
		core::pose::Pose & pose,
		core::Size const & which_jump,
		bool & success ) const;


	void
	set_chain_connections(  utility::vector1 < std::pair< utility::vector1 <core::Size >, utility::vector1 <core::Size > > > const & setting ) {
		chain_connections_ = setting;
	}

	core::pose::rna::RNA_BasePairList
	rna_pairing_list() const { return rna_pairing_list_; }

	void
	set_rna_pairing_list( core::pose::rna::RNA_BasePairList const & setting ) { rna_pairing_list_ = setting; }

	core::pose::toolbox::AtomLevelDomainMapCOP atom_level_domain_map() const { return atom_level_domain_map_; }

private:

	bool
	check_forward_backward(
		core::pose::Pose & pose,
		Size const jump_pos ) const;

	Size
	check_in_chain_connections( Size const & pos1, Size const & pos2 ) const;

	void
	sample_alternative_chain_connection( core::pose::Pose & pose, Size const & which_jump ) const;

private:

	core::import_pose::libraries::RNA_JumpLibraryCOP rna_jump_library_;
	core::pose::toolbox::AtomLevelDomainMapCOP atom_level_domain_map_;
	core::pose::rna::RNA_BasePairList rna_pairing_list_;
	utility::vector1 < std::pair< utility::vector1 <core::Size >, utility::vector1 <core::Size > > > chain_connections_;

};

} //movers
} //protocols

#endif
