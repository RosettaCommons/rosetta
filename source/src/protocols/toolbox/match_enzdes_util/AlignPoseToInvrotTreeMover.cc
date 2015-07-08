// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
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
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh>
#include <protocols/toolbox/match_enzdes_util/MatchConstraintFileInfo.hh>


///project headers
#include <basic/options/option.hh>
#include <basic/options/keys/match.OptionKeys.gen.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/rms_util.hh>
#include <core/types.hh>
#include <core/util/SwitchResidueTypeSet.hh>

#include <basic/Tracer.hh>

#include <numeric/random/random.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.fwd.hh>

// C++ headers

namespace protocols {
namespace toolbox {
namespace match_enzdes_util {

static thread_local basic::Tracer TR( "protocols.toolbox.match_enzdes_util.AlignPoseToInvrotTreeMover" );


AlignPoseToInvrotTreeMover::AlignPoseToInvrotTreeMover(
  InvrotTreeCOP invrot_tree,
  AllowedSeqposForGeomCstCOP seqpos
  ) : Mover(),
      add_target_to_pose_(false), invrot_tree_(invrot_tree),
      seqpos_(seqpos), all_invrots_( invrot_tree->collect_all_inverse_rotamers() )
{
	geomcsts_for_superposition_.clear();
	//this is assuming that every InvrotCollector has the same number of geomcsts.
	//i can't imagine a situation where this wouldn't be the case...
	for( core::Size i =1; i <= all_invrots_[1]->invrots().size() - 1; ++i) geomcsts_for_superposition_.push_back(i);
}

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
/// WARNING: right now this absolutely only works for InvrotTrees
///          where the target is one ligand
void
AlignPoseToInvrotTreeMover::apply( core::pose::Pose & pose ){

  //TR << "pose fold tree at apply start: " << pose.fold_tree() << std::endl; //debug
  //for now, this only works for trees that have one target only
  //changing it won't be hard, but need to think about how to communicate
  //chosen state between this mover and other things using it
  //prolly just use pose datacache
  runtime_assert( invrot_tree_->num_target_states() == 1 );

  //1a. get the invrots and pick a random one from the first list
  //utility::vector1< InvrotCollectorCOP > all_invrots(  );
  Size picked_collector( numeric::random::random_range( 1, all_invrots_.size() ) );
  Size picked_geomcst( geomcsts_for_superposition_[ numeric::random::random_range( 1, geomcsts_for_superposition_.size() ) ] );
  Size picked_rotamer( numeric::random::random_range(1, all_invrots_[ picked_collector ]->invrots()[picked_geomcst].size() ) );

  //temp debug
  //picked_collector = 1;
  //picked_geomcst = 1;
  //picked_rotamer = 1;
  //temp debug over

  std::list<core::conformation::ResidueCOP>::const_iterator list_it( all_invrots_[ picked_collector ]->invrots()[picked_geomcst].begin() );
  for( Size i =1; i < picked_rotamer; ++i ) ++list_it; //not ideal, but a list is what we have
  core::conformation::ResidueCOP ranrot( *list_it );

  //1b. superimpose pose onto the ranrot,
  //need to create a pose to use existing functionality
  core::pose::Pose temp_pose;
  temp_pose.append_residue_by_jump( *ranrot, (Size) 0 );
  //temp_pose.dump_pdb("align_to_invrot_template.pdb");

  //pick a random residue
  Size picked_seqpos( seqpos_->seqpos_for_geomcst(picked_geomcst)[ numeric::random::random_range(1,seqpos_->seqpos_for_geomcst(picked_geomcst).size() ) ] );
  //temp debug
  //picked_seqpos = seqpos_->seqpos_for_geomcst(picked_geomcst)[1];
  //TR << "there are " << all_invrots.size() << " invrot collectors, picked collector has " << all_invrots[ picked_collector ]->invrots()[0].size() << " rotamers in 0th element." << std::endl;
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
  std::list<core::conformation::ResidueCOP>::const_iterator target_it( all_invrots_[ picked_collector ]->invrots()[0].begin() );
  Size first_target_seqpos( pose.total_residue() );
  if( add_target_to_pose_ ){
    first_target_seqpos++;
    //have to add something to switch target res to centroid here
    core::conformation::ResidueCOP ligres( this->switch_residue_type_set( *target_it,  pose.residue(1).residue_type_set().name()) );

    pose.append_residue_by_jump( *ligres, pose.total_residue() );
    ++target_it;
    //Size jump_num = pose.num_jump();

    //below commented out for now. need to think about how to best approach a
    //case where the ligand can have different rotameric states
    //for( ; target_it != all_invrots[ picked_collector ]->invrots()[0].end(); ++target_it){
    //  core::conformation::ResidueCOP ligres( this->switch_residue_type_set( *target_it,  pose.residue(1).residue_type_set().name()) );
    //  pose.append_residue_by_jump( *ligres, first_target_seqpos );
    //}
  }
  // if the target already was in the pose, that means its
  // position got fucked up during the above superimpose call,
  // so we need to reset it to the position in the invrot tree
  // current implemenation absolutely only works for one ligand case
  else{
    pose.replace_residue( first_target_seqpos, *(this->switch_residue_type_set( *target_it,  pose.residue(1).residue_type_set().name())), false );
  }

  //2b. fold tree setup
  //TR << "pose fold tree before mod: " << pose.fold_tree() << std::endl; //debug
  this->setup_foldtree_around_anchor_invrot( pose, picked_seqpos, first_target_seqpos );
  //TR << "pose fold tree after mod: " << pose.fold_tree() << std::endl; //debug

}

void
AlignPoseToInvrotTreeMover::set_add_target_to_pose(
  bool const setting
)
{
  add_target_to_pose_ = setting;
}

void
AlignPoseToInvrotTreeMover::set_geomcst_for_superposition_from_enz_io(
	EnzConstraintIOCOP enzcst_io )
{

	utility::vector1< core::Size > found_geomcsts;
	for( core::Size i = 1; i <= enzcst_io->num_mcfi_lists(); ++i ){
		std::map< std::string, utility::vector1< std::string > > const &
			alg_info(  enzcst_io->mcfi_list( i )->mcfi( 1 )->algorithm_inputs() );

		if( alg_info.find("invrot_tree") != alg_info.end() ){
			utility::vector1< std::string > const & info( alg_info.find( "invrot_tree" )->second );
			for( core::Size line = 1; line <= info.size(); ++line ){
				std::istringstream infostr( info[ line ] );
				std::string first;
				infostr >> first;
				std::cout << "'" << first << "'" << std::endl;
				if( first == "superpose" ){
					found_geomcsts.push_back( i );
					break;
				}
			}
		}
	}
	if( found_geomcsts.size() != 0 ){
		geomcsts_for_superposition_ = found_geomcsts;
		TR << "Superposition of the pose allowed on invrot/seqpos pairs from the following constraint blocks: ";
		for( core::Size i =1; i <= geomcsts_for_superposition_.size(); ++i ) TR << geomcsts_for_superposition_[i] << " ";
		TR << std::endl;
	}
}

/// @details the simplest possible implementation for now
/// assumes the pose only has one chain
void
AlignPoseToInvrotTreeMover::setup_foldtree_around_anchor_invrot(
  core::pose::Pose & pose,
   Size const anchor_seqpos,
  Size const first_target_seqpos ) const
{

  using namespace core::kinematics;
  FoldTree new_fold_tree;
  new_fold_tree.add_edge( anchor_seqpos, 1, Edge::PEPTIDE );
  new_fold_tree.add_edge( anchor_seqpos, first_target_seqpos - 1, Edge::PEPTIDE );
  Size num_jumps_to_add( pose.total_residue() - first_target_seqpos + 1 );
  for( Size i =0; i < num_jumps_to_add; ++i ){
    //TR << "URZ adding jump between res " << anchor_seqpos << " of restype " << pose.residue_type( anchor_seqpos ).name() << " and " << first_target_seqpos + i << ", which is of restype " << pose.residue_type(  first_target_seqpos + i ).name() << std::endl;
    new_fold_tree.add_edge( anchor_seqpos, first_target_seqpos +i, i + 1 );
  }
  if( !new_fold_tree.check_fold_tree() ) {
    utility_exit_with_message("Invalid fold tree after trying to set up around invrot anchor residue");
  }
  pose.fold_tree( new_fold_tree );
}

core::conformation::ResidueCOP
AlignPoseToInvrotTreeMover::switch_residue_type_set(
  core::conformation::ResidueCOP residue,
  std::string const desired_restype_set_name
) const{

  if( desired_restype_set_name != residue->residue_type_set().name() ){
    core::pose::PoseOP temp_pose( new core::pose::Pose() );
    temp_pose->append_residue_by_jump( *residue, (Size) 0 );
    core::util::switch_to_residue_type_set( *temp_pose, desired_restype_set_name );
    residue = temp_pose->residue(1).get_self_ptr();
  }
  return residue;
}

}
}
}
