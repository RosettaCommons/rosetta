// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file util
/// @brief protocols that are specific to docking low resolution
/// @detailed
/// @author Brian Weitzner

#include <protocols/docking/util.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <utility/vector1.hh>
#include <core/conformation/Conformation.hh>
#include <protocols/scoring/Interface.hh>
#include <core/types.hh>
#include <core/pose/PDBInfo.hh>
#include <protocols/geometry/RB_geometry.hh>
#include <basic/Tracer.hh>

#include <protocols/scoring/InterfaceInfo.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>

//#include <basic/datacache/DiagnosticData.hh>


static basic::Tracer TR("protocols.docking.util");

namespace protocols {
namespace docking {

///@brief setup the docking fold tree
///@details ensure that the fold tree is set up such that the jump points
///		are at the center of masses of the two partners
//

void
setup_foldtree( core::pose::Pose & pose, std::string const partner_chainID, DockJumps movable_jumps )
{

	using namespace core;
	//using namespace basic::options;
	runtime_assert( movable_jumps.size() > 0 );

	// identify the chainIDs for partner1 and partner2
	//TR.Debug << "Get the jump point for jump " << rb_jump_ << "............." << std::endl;

	core::pose::PDBInfoCOP pdb_info = pose.pdb_info();
	runtime_assert( pdb_info );
	char second_chain = '_';
	Size cutpoint = 0;

	using namespace kinematics;
	FoldTree f( pose.fold_tree() );

	// identify cutpoint for first movable jump
	if ( partner_chainID == "_"){//chains not specified
		if (f.num_jump()){	//if a jump exists, use first jump for cutpoint
			Size const min_moving_jump = *(std::min_element( movable_jumps.begin(), movable_jumps.end() )); //dereference iterator returned by min_element
			cutpoint = f.cutpoint_by_jump( min_moving_jump );
		} else {	//otherwise use second chain. this should never happen! docking with <1 jump makes no sense.
			second_chain = pdb_info->chain( pose.total_residue() );
		}
	}
	else {//chains specified only works with two-body docking
		for (Size i=1; i<=partner_chainID.length()-1; i++){ //identify second chain from input partner_chainID
			if (partner_chainID[i-1] == '_') second_chain = partner_chainID[i];
		}
		for ( Size i=2; i<= pose.total_residue(); ++i ) {
			if(pdb_info->chain( i ) == second_chain){ //identify cutpoint corresponding to second chain in partner_chainID
				cutpoint = i-1;
				break;
			}
		}
		// Specifying chains overrides movable_jumps
		for (Size i=1; i<=f.num_jump(); ++i){
			Size const current_cutpoint = f.cutpoint_by_jump(i);
			if( current_cutpoint  == cutpoint){
				movable_jumps.clear();
				movable_jumps.push_back( i );
			}
		}
	}

	// calculate jump points from partners
	runtime_assert( cutpoint );
	runtime_assert( movable_jumps.size() > 0 );

	// build docking fold tree one of two ways: two-body docking or multibody docking

	// first case: two-body docking, one jump movable
	// SJF The default foldtree (when the pose is read from disk) sets all the jumps from the first residue to the first residue of each chain.
	// We want a rigid body jump to be between centres of mass of the two partners (jump_pos1, 2) and all of the rest of the jumps to be
	// sequential so that all of the chains are traversed from the jump position onwards
	// default tree. pas simple...
	if( movable_jumps.size() == 1 ) {

		//identify center of masses for jump points
		Size jump_pos1 ( geometry::residue_center_of_mass( pose, 1, cutpoint ) );
		Size jump_pos2 ( geometry::residue_center_of_mass( pose, cutpoint+1, pose.total_residue() ) );
		TR.Debug << "cutpoint: " << cutpoint << std::endl;
		TR.Debug << "jump1: " << jump_pos1 << std::endl;
		TR.Debug << "jump2: " << jump_pos2 << std::endl;

		//setup fold tree based on cutpoints and jump points
		f.clear();
		f.simple_tree( pose.total_residue() );
		f.new_jump( jump_pos1, jump_pos2, cutpoint);
		movable_jumps.clear();
		movable_jumps.push_back( 1 );

		Size chain_begin(0), chain_end(0);

		//rebuild jumps between chains N-terminal to the docking cutpoint
		chain_end = cutpoint;
		chain_begin = pose.conformation().chain_begin( pose.chain(chain_end) );
		while (chain_begin != 1){
			chain_end = chain_begin-1;
			f.new_jump( chain_end, chain_begin, chain_end);
			chain_begin = pose.conformation().chain_begin( pose.chain(chain_end) );
		}

		//rebuild jumps between chains C-terminal to the docking cutpoint
		chain_begin = cutpoint+1;
		chain_end = pose.conformation().chain_end( pose.chain(chain_begin) );
		while (chain_end != pose.total_residue()){
			chain_begin = chain_end+1;
			f.new_jump( chain_end, chain_begin, chain_end);
			chain_end = pose.conformation().chain_end( pose.chain(chain_begin) );
		}
	}

	// second case: multibody docking, more than one jump movable
	// anchor all jumps relative to the CoM of the "downstream" chains, which are defined as first sequential chains that are not movable
	// this will always be at least chain 1, but could be chains 1+2+...
	// the jumps for all nonmoving chains are left alone -- they anchor to N terminal residue
	else {
		std::sort( movable_jumps.begin(), movable_jumps.end() ); // sort jumps to be sure they're ordered

		//Size const min_moving_jump = *(std::min_element( movable_jumps_.begin(), movable_jumps_.end() )); //dereference iterator returned by min_element
		Size const base_cutpoint = cutpoint;
		Size const base_jump_pos( geometry::residue_center_of_mass( pose, 1, base_cutpoint ) );
		for( utility::vector1_int::const_iterator it = movable_jumps.begin(); it != movable_jumps.end(); ++it ) {
			Size const curr_jump = *it;
			Size const curr_cutpoint = f.cutpoint_by_jump( curr_jump ); // used to get the index of a residue in the moving chain (curr_cutpoint+1)
			Size const chain_begin = pose.conformation().chain_begin( pose.chain(curr_cutpoint+1) );
			Size const chain_end = pose.conformation().chain_end( pose.chain(curr_cutpoint+1) );
			Size const moving_jump_pos( geometry::residue_center_of_mass( pose, chain_begin, chain_end ) );
			TR.Debug << "Adjusting Jump (cut) for #" << curr_jump << "("<<curr_cutpoint<<")" << ": begin " << chain_begin << "    end " << chain_end << "      base_cutpoint " << base_cutpoint<< "         base_jump_pos " << base_jump_pos << "      moving_jump_pos " << moving_jump_pos << std::endl;
			f.slide_jump( curr_jump, base_jump_pos, moving_jump_pos );
		}
	}

	// set docking fold tree to the pose
	f.reorder( 1 );
	f.check_fold_tree();
	runtime_assert( f.check_fold_tree() );
	pose.fold_tree( f );

	//set up InterfaceInfo object in pose to specify which interface(s) to calculate docking centroid mode scoring components from
	using namespace core::scoring;
	using namespace protocols::scoring;
	//using core::pose::datacache::CacheableDataType::INTERFACE_INFO;

	InterfaceInfoOP docking_interface = new InterfaceInfo( movable_jumps );
	pose.data().set( core::pose::datacache::CacheableDataType::INTERFACE_INFO, docking_interface );
}

} //docking
} //protocols
