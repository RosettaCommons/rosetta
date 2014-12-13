// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/docking/util.cc

/// @brief Construct a foldtree from a description of the chains
/// involved in the interface

/// @details
/// @author Brian Weitzner
/// @author Matthew O'Meara


// Unit Headers
#include <protocols/docking/util.hh>
#include <protocols/scoring/InterfaceInfo.hh>

// Platform headers
#include <protocols/rigid/RB_geometry.hh>
#include <basic/Tracer.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/database/sql_utils.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/conformation/Conformation.hh>
#include <core/types.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/util.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/string_util.hh>

// External Headers
#include <cppdb/frontend.h>

// C++ Headers
#include <map>
#include <sstream>

static thread_local basic::Tracer TR( "protocols.docking.util" );

namespace protocols {
namespace docking {

using std::endl;
using std::map;
using std::min_element;
using std::pair;
using std::string;
using std::stringstream;
using cppdb::statement;
using cppdb::result;
using utility::vector1;
using utility::sql_database::sessionOP;
using basic::database::safely_prepare_statement;
using basic::database::safely_read_from_database;
using core::Size;
using core::conformation::Conformation;
using core::kinematics::FoldTree;
using core::pose::Pose;
using core::pose::PDBInfoCOP;
using protocols::scoring::InterfaceInfoOP;

/// @brief Setup foldtree for across an interface specified by a
/// string for the partner chains (using pdb_chain
/// identification). The foldtree is set up such that the jump points
/// are at the center of masses of the two partners
///
/// @detailed
///  if partner_chainID == '_' -> use the first movable jump as the cutpoint
///  if partner_chainID is of the form (pdb_chain_id)+_(pdb_chain_id)+
///      then use the jump between the last chain of the first partner
///      and the first chain of the second partner as the cutpoint.
///      For example 'ABC_DEF' has chains 'ABC' as the first partner and 'DEF' as
///      the second partner.
///
/// With the current implementation, all of the chains of the first
/// partner must be upstream of the chains of the second partner.  If
/// this behavior is too limiting, then the behavior can be extended.
///
///
/// @return The foldtree in the pose will have a jump from the center
/// of mass of the first partner to the center of mass of the second
/// partner and no other jumps between residues in different partners
///
/// @return the movable_jumps vector contains as it's only entry the
/// number of the jump across the interface.
void
setup_foldtree(
	Pose & pose,
	string const & partner_chainID,
	DockJumps & movable_jumps) 
{
	core::kinematics::FoldTree f;
	setup_foldtree(pose, partner_chainID, movable_jumps, f);
	pose.fold_tree(f);

}

void
setup_foldtree(
	core::pose::Pose const & pose,
	std::string const & partner_chainID,
	DockJumps & movable_jumps,
	core::kinematics::FoldTree & ft)
{
	PDBInfoCOP pdb_info = pose.pdb_info();
	movable_jumps.clear(); //Why is this required and then cleared?

	if(!pdb_info){
		utility_exit_with_message("Attempting to setup foldtree between interface partners, however, the pdb_info object associated with the pose does not exits.");
	}

	// identify cutpoint for first movable jump
	Size cutpoint = 0;
	FoldTree const & f( pose.fold_tree() );
	if ( partner_chainID == "_") {
		// By default, set the first jump as the cutpoint

		if(f.num_jump() == 0){
			utility_exit_with_message("Attempting to auto-detect interface partner chains, however the pose contains no jumps.");
		}
		cutpoint = f.cutpoint_by_jump(1);
		movable_jumps.push_back(1);
	} else {
		char first_chain_second_partner = char();
		for ( Size i=1; i<=partner_chainID.length()-1; i++ ) {
			
			if (partner_chainID[i-1] == '_') {
				first_chain_second_partner = partner_chainID[i];
			}
		}
		
		cutpoint = pose.conformation().chain_begin(core::pose::get_chain_id_from_chain(first_chain_second_partner, pose)) - 1;
		
		if (cutpoint == 0){
			utility_exit_with_message(
				"Attempting to setup foldtree between interface partners, "
				"however, the cutpoint could not be identified because "
				"the first chain second partner cannot currently be the first chain in the PDB: " +utility::to_string(first_chain_second_partner));
		}
		
		/*
		for ( Size i=2; i<= pose.total_residue(); ++i ) {
			if ( pdb_info->chain( i ) == first_chain_second_partner ) {
				cutpoint = i-1;
				break;
			}
		}
		if(!cutpoint){
			utility_exit_with_message(
				"Attempting to setup foldtree between interface partners, "
				"however, the cutpoint could not be identified because "
				"the pdb chain identifier could not be found for "
				"the first chain of the second partner");
		}
		*/

		// Specifying chains overrides movable_jumps
		for (Size i=1; i<=f.num_jump(); ++i) {
			Size const current_cutpoint = f.cutpoint_by_jump(i);
			if( current_cutpoint  == cutpoint ) {
				movable_jumps.push_back( i );
			}
		}


	}

	setup_foldtree(pose, cutpoint, movable_jumps, ft);
}

/// Here is the protocol that the database should have:
///
///
///CREATE TABLE IF NOT EXISTS interfaces (
///	struct_id BLOB,
///	interface_id INTEGER,
///	FOREIGN KEY (struct_id) REFERENCES structures (struct_id)	DEFERRABLE INITIALLY DEFERRED,
///	PRIMARY KEY(struct_id, interface_id));
///
///CREATE TABLE IF NOT EXISTS interface_partners (
///	interface_id INTEGER,
///	partner_id INTEGER,
///	FOREIGN KEY (interface_id) REFERENCES interfaces (interface_id)	DEFERRABLE INITIALLY DEFERRED,
///	PRIMARY KEY(interface_id, partner_id));
///
///CREATE TABLE IF NOT EXISTS interface_partner_chains (
///	partner_id INTEGER,
///	chain_id INTEGER,
///	FOREIGN KEY (interface_interface_id) REFERENCES interface_partners (interface_partner_id)	DEFERRABLE INITIALLY DEFERRED,
///	PRIMARY KEY(partner_id, chain_id));
void
setup_foldtree(
	Pose & pose,
	Size const interface_id,
	sessionOP db_session,
	DockJumps & movable_jumps) 
{
	core::kinematics::FoldTree f;
	setup_foldtree(pose, interface_id, db_session, movable_jumps, f);
	pose.fold_tree(f);
}

void
setup_foldtree(
	core::pose::Pose const & pose,
	core::Size interface_id,
	utility::sql_database::sessionOP db_session,
	DockJumps & movable_jumps,
	core::kinematics::FoldTree & ft)
{
	string const sele(
		"SELECT\n"
		"	partner.partner_id AS partner,\n"
		"	chain.chain_id AS chain\n"
		"FROM\n"
		"	interface_partners AS partner,\n"
		"	interface_partner_chains AS chain\n"
		"WHERE\n"
		"	partner.interface_id = ? AND\n"
		"	chain.partner_id = partner.partner_id;");
	statement stmt(safely_prepare_statement(sele, db_session));
	stmt.bind(1, interface_id);
	result res(safely_read_from_database(stmt));

	map<Size, vector1< Size > > partner_to_chains;

	while(res.next()){
		Size partner, chain;
		res >> partner >> chain;
		partner_to_chains[partner].push_back(chain);
	}

	setup_foldtree(pose, partner_to_chains, movable_jumps, ft);
}


void
setup_foldtree(
	Pose & pose,
	map< Size, vector1< Size > > const & partner_to_chains,
	DockJumps & movable_jumps) 
{
	core::kinematics::FoldTree f;
	setup_foldtree(pose, partner_to_chains, movable_jumps, f);
	pose.fold_tree(f);
}

void
setup_foldtree(
	core::pose::Pose const & pose,
	std::map< core::Size, utility::vector1< core::Size > > const & partner_to_chains,
	DockJumps & movable_jumps,
	core::kinematics::FoldTree & ft)
{
	if(partner_to_chains.size() != 2){
		stringstream error_msg;
		error_msg
			<< "The current implementation of setup_foldtree requires "
			<< "that there be exactly two chains at an interface, however, "
			<< "'" << partner_to_chains.size()  << "' partners were given.";
		utility_exit_with_message(error_msg.str());
	}

	// The current implementation of the setup_foltree assumes all the
	// chains of one partner come before all the chains of the second
	// partner.
	vector1< pair< Size, Size > > first_and_last_chains;
	for(map< Size, vector1 < Size > >::const_iterator
				p_cs = partner_to_chains.begin(),
				p_cs_end = partner_to_chains.end();
				p_cs != p_cs_end; ++p_cs){

		first_and_last_chains.push_back(
			pair< Size, Size>(
				*min_element(p_cs->second.begin(), p_cs->second.end()),
				*max_element(p_cs->second.begin(), p_cs->second.end())));
	}

	Size first_chain_second_partner;
	if(first_and_last_chains[1].second < first_and_last_chains[2].first){
		first_chain_second_partner = first_and_last_chains[2].first;
	} else if(first_and_last_chains[2].second < first_and_last_chains[1].first){
		first_chain_second_partner = first_and_last_chains[1].first;
	} else {
		utility_exit_with_message(
			"The current implementation of setup_foldtree requires "
			"that the chains of one partner come before all "
			"the chains of the second partner.");
	}
	Size cutpoint = 0;
	for ( Size i=2; i<= pose.total_residue(); ++i ) {
		if ( (Size)pose.chain(i) == first_chain_second_partner ) {
			cutpoint = i-1;
			break;
		}
	}
	if(!cutpoint){
		utility_exit_with_message(
			"Attempting to setup foldtree between interface partners, "
			"however, the cutpoint could not be identified because "
			"the pdb chain identifier could not be found for "
			"the first chain of the second partner");
	}

	movable_jumps.clear();
	FoldTree const & f(pose.fold_tree());

	// Specifying chains overrides movable_jumps
	for (Size i=1; i<=f.num_jump(); ++i) {
		Size const current_cutpoint = f.cutpoint_by_jump(i);
		if( current_cutpoint  == cutpoint ) {
			movable_jumps.push_back( i );
		}
	}

	setup_foldtree(pose, cutpoint, movable_jumps, ft);
}

void
setup_foldtree(
	Pose & pose,
	Size const cutpoint,
	DockJumps & movable_jumps
){
	core::kinematics::FoldTree f;
	setup_foldtree(pose, cutpoint, movable_jumps);
	pose.fold_tree(f);
}

void
setup_foldtree(
	core::pose::Pose const & pose,
	core::Size cutpoint,
	DockJumps & movable_jumps,
	core::kinematics::FoldTree & ft)
{
	runtime_assert_string_msg(
		movable_jumps.size() > 0,
		"Unable to set up interface foldtree because there are no movable jumps");

	Conformation const & conformation(pose.conformation());

	// first case: two-body docking, one jump movable

	// SJF The default foldtree (when the pose is read from disk) sets
	// all the jumps from the first residue to the first residue of each
	// chain.  We want a rigid body jump to be between centres of mass
	// of the two partners (jump_pos1, 2) and all of the rest of the
	// jumps to be sequential so that all of the chains are traversed
	// from the jump position onwards default tree. pas simple...
	if( movable_jumps.size() == 1 ) {

		//identify center of masses for jump points
		Size jump_pos1 ( core::pose::residue_center_of_mass( pose, 1, cutpoint ) );
		Size jump_pos2 ( core::pose::residue_center_of_mass( pose, cutpoint+1, pose.total_residue() ) );
		TR.Debug << "cutpoint: " << cutpoint << std::endl;
		TR.Debug << "jump1: " << jump_pos1 << std::endl;
		TR.Debug << "jump2: " << jump_pos2 << std::endl;

		//setup fold tree based on cutpoints and jump points
		ft.clear();
		ft.simple_tree( pose.total_residue() );
		ft.new_jump( jump_pos1, jump_pos2, cutpoint);
		movable_jumps.clear();
		movable_jumps.push_back( 1 );

		Size chain_begin(0), chain_end(0);

		//rebuild jumps between chains N-terminal to the docking cutpoint
		chain_end = cutpoint;
		chain_begin = conformation.chain_begin( pose.chain(chain_end) );
		while (chain_begin != 1){
			chain_end = chain_begin-1;
			ft.new_jump( chain_end, chain_begin, chain_end);
			chain_begin = conformation.chain_begin( pose.chain(chain_end) );
		}

		//rebuild jumps between chains C-terminal to the docking cutpoint
		chain_begin = cutpoint+1;
		chain_end = conformation.chain_end( pose.chain(chain_begin) );
		while (chain_end != pose.total_residue()){
			chain_begin = chain_end+1;
			ft.new_jump( chain_end, chain_begin, chain_end);
			chain_end = conformation.chain_end( pose.chain(chain_begin) );
		}
	}

	// second case: multibody docking, more than one jump movable

	// anchor all jumps relative to the CoM of the "downstream" chains,
	// which are defined as first sequential chains that are not movable
	// this will always be at least chain 1, but could be chains 1+2+...
	// the jumps for all nonmoving chains are left alone -- they anchor
	// to N terminal residue
	else {
		std::sort( movable_jumps.begin(), movable_jumps.end() );

		Size const base_cutpoint = cutpoint;
		Size const base_jump_pos( core::pose::residue_center_of_mass( pose, 1, base_cutpoint ) );
		for(DockJumps::const_iterator
					curr_jump = movable_jumps.begin(),
					last_movable_jump = movable_jumps.end();
				curr_jump != last_movable_jump; ++curr_jump ) {
			Size const curr_cutpoint = ft.cutpoint_by_jump( *curr_jump ); // used to get the index of a residue in the moving chain (curr_cutpoint+1)
			Size const chain_begin = conformation.chain_begin( pose.chain(curr_cutpoint+1) );
			Size const chain_end = conformation.chain_end( pose.chain(curr_cutpoint+1) );
			Size const moving_jump_pos( core::pose::residue_center_of_mass( pose, chain_begin, chain_end ) );
			TR.Debug
				<< "Adjusting Jump (cut) for #" << *curr_jump
				<< "(" << curr_cutpoint << "): "
				<< "begin " << chain_begin
				<< "    end " << chain_end
				<< "      base_cutpoint " << base_cutpoint
				<< "         base_jump_pos " << base_jump_pos
				<< "      moving_jump_pos " << moving_jump_pos << endl;
			ft.slide_jump( *curr_jump, base_jump_pos, moving_jump_pos );
		}
	}

	// set docking fold tree to the pose
	ft.reorder( 1 );
	runtime_assert( ft.check_fold_tree() );
	
	/* Moved to DockLowRes!!! If it is 2015+, delete me!!

	//set up InterfaceInfo object in pose to specify which interface(s)
	//to calculate docking centroid mode scoring components from
	using namespace core::scoring;
	using namespace protocols::scoring;
	//using core::pose::datacache::CacheableDataType::INTERFACE_INFO;

	InterfaceInfoOP docking_interface( new InterfaceInfo( movable_jumps ) );
	pose.data().set(
		core::pose::datacache::CacheableDataType::INTERFACE_INFO,
		docking_interface);
	 */
}

} //docking
} //protocols
