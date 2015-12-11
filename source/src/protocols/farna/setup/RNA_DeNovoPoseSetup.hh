// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @author Rhiju Das

#ifndef INCLUDED_protocols_rna_RNA_DeNovoPoseSetup_HH
#define INCLUDED_protocols_rna_RNA_DeNovoPoseSetup_HH

#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <protocols/farna/movers/RNA_JumpMover.fwd.hh>
#include <protocols/farna/base_pairs/BasePairStep.hh>
#include <protocols/farna/setup/RNA_DeNovoParameters.hh>
#include <protocols/toolbox/AtomLevelDomainMap.fwd.hh>
#include <core/types.hh>
#include <utility/pointer/ReferenceCount.hh>
// ObjexxFCL Headers

// C++ Headers
#include <string>
#include <map>
#include <iostream>
#include <list>

#include <utility/vector1.hh>

using namespace protocols::farna::movers;
using namespace protocols::farna::base_pairs;

namespace protocols {
namespace farna {
namespace setup {

//////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Object used in RNA_FragmentMonteCarlo to handle setup of pose & fold-tree, and
///   jump changes.
class RNA_DeNovoPoseSetup : public utility::pointer::ReferenceCount {
public:

	//constructor
	RNA_DeNovoPoseSetup( RNA_DeNovoParameters const & rna_params_file_info );
	virtual ~RNA_DeNovoPoseSetup();

	/// @brief "classic" setup, used in denovo protocol. Note: changes pose (virtualizes phosphate).
	///    and requires later call of setup_fold_tree_and_jumps_and_variants
	void
	initialize_for_de_novo_protocol(
		core::pose::Pose & pose,
		bool const ignore_secstruct = false );

	void
	setup_fold_tree_and_jumps_and_variants( core::pose::Pose & pose,
		RNA_JumpMover const & rna_jump_mover,
		protocols::toolbox::AtomLevelDomainMapOP atom_level_domain_map ) const;

	void
	setup_fold_tree_and_jumps_and_variants( core::pose::Pose & pose ) const;

	void
	set_root_at_first_rigid_body( bool const setting ){ root_at_first_rigid_body_ = setting; }

	void set_bps_moves( bool const & setting ){ bps_moves_ = setting; }

	RNA_DeNovoParameters const & rna_params() const { return rna_params_; }

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
		RNA_JumpMover const & rna_jump_mover ) const;

	void
	setup_chainbreak_variants( core::pose::Pose & pose,
		protocols::toolbox::AtomLevelDomainMapOP atom_level_domain_map ) const;

	void
	insert_base_pair_jumps( core::pose::Pose & pose, RNA_JumpMover const & jump_mover,  bool & success ) const;

private:

	RNA_DeNovoParameters rna_params_;

	bool const assume_non_stem_is_loop; // legacy parameter.
	bool bps_moves_;
	bool root_at_first_rigid_body_;

};


} //setup
} //farna
} //protocols

#endif
