// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/match_enzdes_util/InvrotTreeNodeBase.hh
/// @brief  Forward declaration for inverse rotamer tree node base
/// @author Florian Richter, flosopher@gmail.com, mar 2012


/// unit headeers
#include <protocols/toolbox/match_enzdes_util/AlignPoseToInvrotTreeMover.hh>

// package headers
#include <protocols/toolbox/match_enzdes_util/InvrotTree.hh>
#include <protocols/toolbox/match_enzdes_util/InvrotTarget.hh>
#include <protocols/toolbox/match_enzdes_util/AllowedSeqposForGeomCst.hh>


///project headers
#include <basic/options/option.hh>
#include <basic/options/keys/match.OptionKeys.gen.hh>

#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/rms_util.hh>
#include <core/types.hh>

#include <basic/Tracer.hh>

#include <numeric/random/random.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.fwd.hh>

// C++ headers

namespace protocols {
namespace toolbox {
namespace match_enzdes_util {

static basic::Tracer TR( "protocols.toolbox.match_enzdes_util.AlignPoseToInvrotTreeMover" );


AlignPoseToInvrotTreeMover::AlignPoseToInvrotTreeMover(
  InvrotTreeCOP invrot_tree,
  AllowedSeqposForGeomCstCOP seqpos
  ) : Mover(),
      invrot_tree_(invrot_tree), seqpos_(seqpos)
{}

AlignPoseToInvrotTreeMover::~AlignPoseToInvrotTreeMover(){}

std::string
AlignPoseToInvrotTreeMover::get_name() const{
  return "AlignPoseToInvrotTreeMover";
}


/// @details this function does two things
/// 1. align one of the seqpos in the pose onto one of the invrots
///    in the tree
/// 2. grab the InvrotTarget residues and add them to the pose
///    this entails setting up the foldtree such that the target
///    residues are upstream of the rest of the pose
void
AlignPoseToInvrotTreeMover::apply( core::pose::Pose & pose ){

  //for now, this only works for trees that have one target only
  //changing it won't be hard, but need to think about how to communicate
  //chosen state between this mover and other things using it
  //prolly just use pose datacache
  runtime_assert( invrot_tree_->num_target_states() == 1 );

  //1a. get the invrots and pick a random one from the first list
  utility::vector1< InvrotCollectorCOP > all_invrots( invrot_tree_->collect_all_inverse_rotamers() );
  Size picked_collector( numeric::random::random_range( 1, all_invrots.size() ) );
  Size picked_geomcst( numeric::random::random_range( 1, all_invrots[picked_collector]->invrots().size() - 1 ) );
  Size picked_rotamer( numeric::random::random_range(1, all_invrots[ picked_collector ]->invrots()[picked_geomcst].size() ) );

  //temp debug
  picked_collector = 1;
  picked_geomcst = 1;
  picked_rotamer = 1;
  //temp debug over

  std::list<core::conformation::ResidueCOP>::const_iterator list_it( all_invrots[ picked_collector ]->invrots()[picked_geomcst].begin() );
  for( Size i =1; i < picked_rotamer; ++i ) list_it++; //not ideal, but a list is what we have
  core::conformation::ResidueCOP ranrot( *list_it );

  //1b. superimpose pose onto the ranrot,
  //need to create a pose to use existing functionality
  core::pose::Pose temp_pose;
  temp_pose.append_residue_by_jump( *ranrot, (Size) 0 );
  //temp_pose.dump_pdb("align_to_invrot_template.pdb");

  //pick a random residue
  Size picked_seqpos( seqpos_->seqpos_for_geomcst(picked_geomcst)[ numeric::random::random_range(1,seqpos_->seqpos_for_geomcst(picked_geomcst).size() ) ] );
  //temp debug
  picked_seqpos = seqpos_->seqpos_for_geomcst(picked_geomcst)[1];
  //temp debug over

  core::id::AtomID_Map< core::id::AtomID > atom_map( core::id::BOGUS_ATOM_ID );
  core::pose::initialize_atomid_map( atom_map, pose, core::id::BOGUS_ATOM_ID );
  atom_map.set( core::id::AtomID(pose.residue( picked_seqpos ).atom_index("CA"), picked_seqpos ), core::id::AtomID( temp_pose.residue(1).atom_index("CA"), 1 ) );
  atom_map.set( core::id::AtomID(pose.residue( picked_seqpos ).atom_index("N"), picked_seqpos ), core::id::AtomID( temp_pose.residue(1).atom_index("N"), 1 ) );
  atom_map.set( core::id::AtomID(pose.residue( picked_seqpos ).atom_index("C"), picked_seqpos ), core::id::AtomID( temp_pose.residue(1).atom_index("C"), 1 ) );
  atom_map.set( core::id::AtomID(pose.residue( picked_seqpos ).atom_index("CB"), picked_seqpos ), core::id::AtomID( temp_pose.residue(1).atom_index("CB"), 1 ) );

  core::scoring::superimpose_pose( pose, temp_pose, atom_map );

  //pose.dump_pdb("align_to_invrot_after_align.pdb");

  //2. now we need to add the target residues to the aligned pose,
  //and try to setup the fold tree the right way
  //2a
  //std::list<core::conformation::ResidueCOP>::const_iterator target_it( all_invrots[ picked_collector ]->invrots()[0].begin() );
  //Size first_target_seqpos( pose.total_residue() + 1 );
  //pose.append_residue_by_jump( **target_it, pose.total_residue() );
  //target_it++;
  //Size jump_num = pose.num_jump();
  //for( ; target_it != all_invrots[ picked_collector ]->invrots()[0].end(); ++target_it){
  //  pose.append_residue_by_jump( **target_it, first_target_seqpos );
  //}

  //2b. fold tree setup

}



}
}
}
