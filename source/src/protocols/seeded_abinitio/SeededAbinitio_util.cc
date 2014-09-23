// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/seeded_abinitio/SeededAbinitio_util.cc
/// @brief
/// @author Eva-Maria Strauch ( evas01@u.washington.edu )

// Unit Headers
#include <protocols/seeded_abinitio/SeededAbinitio_util.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBPoseMap.hh>
#include <core/pose/PDBInfo.hh>
#include <core/conformation/Conformation.hh>

#include <core/conformation/Residue.hh>
//#include <core/chemical/ResidueSelector.hh>
//#include <core/chemical/ResidueTypeSet.hh>
//#include <core/chemical/ResidueTypeSet.fwd.hh>
//#include <core/conformation/ResidueFactory.hh>
#include <core/kinematics/FoldTree.hh>


#include <core/id/SequenceMapping.hh>
#include <core/pose/datacache/CacheableObserverType.hh>
#include <core/pose/datacache/cacheable_observers.hh>
#include <core/pose/datacache/ObserverCache.hh>

#include <core/types.hh>
#include <basic/Tracer.hh>

#include <core/pose/Pose.hh>
#include <core/pose/selection.hh>
#include <utility/vector1.hh>

#include <core/pose/util.hh>
#include <core/pose/PDBInfo.hh>

//#include <basic/options/keys/in.OptionKeys.gen.hh>
//#include <basic/options/option.hh>
//#include <core/id/AtomID_Map.hh>
//#include <core/import_pose/import_pose.hh>
//#include <core/init/init.hh>
//#include <core/io/pdb/pose_io.hh>
//#include <core/scoring/rms_util.hh>

#include <core/conformation/ResidueFactory.hh>
#include <core/kinematics/MoveMap.hh>


#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/Loops.fwd.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/util.hh>

//Auto Headers
#include <core/chemical/AtomType.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>

// C++ headers
#include <map>
#include <algorithm>

static thread_local basic::Tracer TR( "seeded_abinitio.SeededAbinitio_util" );

namespace protocols {
	namespace seeded_abinitio{
		
		using namespace core::scoring;
		//using namespace protocols::moves;
		using namespace core;
		using namespace std;
		using utility::vector1;

		
/// @brief method to update numbering of a movemap if the pose length was changed previously
///	only works if a LengthEventCollector is set in the pose's observer cache 
///default behavior is to set everything new to be movable and adjust the numbering of the residues that were previously
/// not allowed to move to the appropriate new numbering
		
void
adjust_mm_to_length( core::pose::Pose const & pose, core::kinematics::MoveMapOP & mm ){
	
	TR<<"adjusting the movemap to current numbering of the decoy" <<std::endl;
	
	if( !pose.observer_cache().has( pose::datacache::CacheableObserverType::LENGTH_EVENT_COLLECTOR ) ){
		TR<<"WARNING there is no length observer attached to the pose! no adjustments can be made" <<std::endl;
		return;
	}
	
	//get data from the pose observer
	pose::datacache::CacheableObserverCOP len_obs = pose.observer_cache().get_const_ptr( pose::datacache::CacheableObserverType::LENGTH_EVENT_COLLECTOR );
	pose::datacache::LengthEventCollectorCOP lencollect( utility::pointer::static_pointer_cast< pose::datacache::LengthEventCollector const >( len_obs ) );
	
	utility::vector1< core::conformation::signals::LengthEvent > const & events( lencollect->events() );
	utility::vector1< core::id::SequenceMapping > smaps;
	
	for( Size i =1; i <= events.size(); ++i ){
		smaps.push_back( core::id::SequenceMapping( events[i] ) );
	}
	core::id::SequenceMappingOP fullsmap( core::id::combine_sequence_mappings( smaps ) );

	//make a new movemap in which everything is set to true (set_bb..)
	//iterate through the movemap and if get_bb (false) adjust that position
	//replace the old movemap with the new one, assert total residue size with the new mm size and all is good!
	
	core::kinematics::MoveMapOP new_mm( new core::kinematics::MoveMap() );
	new_mm->set_bb( false ); 
 	new_mm->set_chi( false ); 
// 	new_mm->set_jump( false ); 
	
	for( Size resi = 1; resi <= pose.total_residue(); ++resi ){
		Size previous_pos ( (*fullsmap)[resi] );
		TR.Debug<<"previous position "<<previous_pos <<", current residue: " << resi << std::endl;
		
		//if the residue didnt exist before
		if( previous_pos == 0 || !previous_pos) continue;
		
		//if the residue was set to NOT move, update the position
		TR.Debug<<"mm set to: "<< mm->get_bb( previous_pos) << std::endl;
		if( !mm->get_bb( previous_pos ) ){
			TR.Debug<<"adjusting NOT movable position to: " << resi <<std::endl; 
		   	new_mm->set_bb( resi, false );
		}
	   if( !mm->get_chi( previous_pos ) ){
			TR.Debug<<"adjusting NOT movable position to: " << resi <<std::endl; 
			new_mm->set_chi( resi, false );		 
		}   
	}
}	  		

//runtime parsing of seeds
protocols::loops::Loops
parse_seeds( core::pose::Pose const & pose, utility::vector1 < std::pair < std::string, std::string > > seed_vector){

  protocols::loops::Loops tmpseed;

    for( Size iter = 1 ; iter <= seed_vector.size() ; ++ iter ){
      TR.Debug<<"sanity check, seed_vector[iter].first " <<seed_vector[iter].first <<std::endl; ////////////
      core::Size const begin = core::pose::parse_resnum( seed_vector[iter].first, pose ) ;
      core::Size const end   = core::pose::parse_resnum( seed_vector[iter].second, pose );
      //runtime_assert( end > begin );
      //runtime_assert( begin>=1);
      //runtime_assert( end<=pose.total_residue() );
      tmpseed.add_loop( begin , end , 0, 0, false );
  }
  TR.Debug<<"runtime parsed: "<< tmpseed <<std::endl;
  return tmpseed;
}//end parse seeds
		
void
combine_two_poses( core::pose::Pose design_pose , core::pose::PoseOP target_chain ){

	core::pose::PoseOP combo_pose( new core::pose::Pose );
	combo_pose = target_chain;
	
	TR<<"new poseOP total number should contain the additional target chain number: " << combo_pose->total_residue();

	core::pose::PDBInfoOP pdb_info_design( new core::pose::PDBInfo( design_pose ) );
	core::pose::PDBInfoOP pdb_info_target( new core::pose::PDBInfo( *target_chain ) );
	pdb_info_target->set_chains( 'A');
	pdb_info_design->set_chains('B');

	// take the target chain pose and append through a jump the newly folded protein pose

	TR<< "folded proteins has " << design_pose.total_residue() << " total residues \n";
	TR<< "the target protein (chain A) has " << target_chain->total_residue() << " total residues "<<std::endl;

	for ( core::Size i=1; i<=design_pose.total_residue(); ++i ) {
		core::conformation::ResidueCOP new_rsd = design_pose.residue(i).clone();
		if ( i == 1 ) {
			combo_pose->append_residue_by_jump( *new_rsd, combo_pose->total_residue(), "", "", true /*new chain*/ );//anchor and start atomes are not defined, wonder whether the constructor can deal wtih taht
		} 
		else {
			combo_pose->append_residue_by_bond( *new_rsd );
		}
	}//done adding chain by residues

	combo_pose->dump_pdb( "target_plus_folded.pdb" );
	TR.Debug<<" total residues of new pdb: " <<combo_pose->total_residue()<<std::endl;
	
	design_pose = *combo_pose;
	
}//end combining 
		
	/*
void
dump_pymol( std::string fn ) const {
		utility::io::ozstream out( fn );
		if ( out ) {
			out << "from pymol import cmd\n";
			for ( Size i=1; i<=size(); i++ ) {
				out << "cmd.select( \"jumps"<<i<<"\", \"resi "<<jumps_(1,i)<<"+"<<jumps_(2,i)<<"\");" << std::endl;
				out << "cmd.show( \"sticks\", \"jumps"<<i<<"\");\n";
				out << "cmd.select( \"cut"<<i<<"\", \"resi " <<cuts_(i)<<"\");\n";
				out << "cmd.show( \"sphere\",\"cut"<<i<<"\");" << std::endl;
			}
	*/
	}
}
