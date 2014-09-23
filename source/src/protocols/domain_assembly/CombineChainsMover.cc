/// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/domain_assembly/CombineChainsMover.cc
/// @brief  Takes a multichain pose and returns a single-chain one
/// @author Dominik Gront

#include <protocols/domain_assembly/CombineChainsMover.hh>

#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/ChemicalManager.hh>


#include <core/sequence/Sequence.hh>
#include <core/sequence/SequenceAlignment.hh>

#include <protocols/comparative_modeling/ThreadingMover.hh>

// AUTO-REMOVED #include <basic/options/option_macros.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <core/types.hh>

#include <core/pose/annotated_sequence.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace domain_assembly {

void CombineChainsMover::apply( core::pose::Pose & pose ) {

  using namespace basic::options;
  using namespace basic::options::OptionKeys;

	core::sequence::SequenceAlignment aln;
	core::sequence::SequenceOP seq( new core::sequence::Sequence( pose.sequence(), "this comment here is really irrelevant..." ) );

	core::pose::Pose extended_pose;
        core::chemical::ResidueTypeSetCOP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set(
	                option[ in::file::residue_type_set ]() );
	core::pose::make_pose_from_sequence(extended_pose, pose.sequence(), *rsd_set);

	aln.add_sequence(seq);
	aln.add_sequence(seq);
	protocols::comparative_modeling::ThreadingMover mover( aln, pose );
	mover.randomize_loop_coords( false );
	mover.build_loops( false );
	mover.apply( extended_pose );


	extended_pose.dump_pdb("out.pdb");

	pose = extended_pose;
} // apply


}//namespace domain_assembly
}//namespace protocols

