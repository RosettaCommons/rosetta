// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pose/util.cc
/// @brief  Pose class utilities
/// @author Phil Bradley
/// @author Modified by Sergey Lyskov, Rhiju Das, Steven Lewis, Vikram K. Mulligan


// Unit header
#include <core/pose/chains_util.hh>

// Package headers
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>

// Project headers
#include <core/conformation/Conformation.hh>
#include <core/conformation/util.hh>
#include <core/id/DOF_ID_Map.hh>
#include <core/id/Exceptions.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/conformation/Residue.hh>

// Basic headers
#include <basic/Tracer.hh>

// Numeric headers
#include <numeric/xyzVector.string.hh>

// Utility headers
#include <utility/io/izstream.hh>
#include <utility/exit.hh>
#include <utility/string_constants.hh>
#include <utility/string_util.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>
#include <utility/vector1.functions.hh>

// C/C++ headers
#include <cmath>
#include <iostream>
#include <algorithm>
#include <numeric>

// External headers
#include <boost/functional/hash.hpp>

namespace core {
namespace pose {

static basic::Tracer TR( "core.pose.util" );

void jumps_from_pose(core::pose::Pose const & pose, Jumps & jumps) {
	for ( Size i = 1; i <= pose.num_jump(); ++i ) {
		jumps.insert(i);
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void conf2pdb_chain_default_map( core::pose::Pose const & pose, std::map<core::Size,char> & chainmap ) {
	chainmap.clear();
	char letter = 'A';
	for ( core::Size i = 1; i <= pose.conformation().num_chains(); ++i ) {
		chainmap[i] = letter;
		if ( 'Z'==letter ) utility_exit_with_message("too many chains to map to letters!!!");
		letter = static_cast<char>(letter + 1);
	}
}

std::map< core::Size, char > conf2pdb_chain( core::pose::Pose const & pose ) {
	using core::Size;
	using core::pose::PDBInfo;
	typedef std::map< core::Size, char > Conf2PDB;

	// TODO: Wiping out the mapping completely is rather silly
	// Ideally we should either get the "default" PDBInfo (if not present),
	// or if there isn't an inconsistency, we should give as much info as we can, marking the entries which aren't
	Conf2PDB conf2pdb;

	if ( !pose.pdb_info().get() ) {
		TR.Warning << "conf2pdb_chain(): PDBInfo does not exist, returning default map 1=A, 2=B, ..." << std::endl;
		conf2pdb_chain_default_map(pose,conf2pdb);
		return conf2pdb;
	}

	for ( Size i = 1, ie = pose.size(); i <= ie; ++i ) {
		core::Size const conf = pose.chain( i );
		char const pdb = pose.pdb_info()->chain( i );

		auto c2p = conf2pdb.find( conf );
		if ( c2p != conf2pdb.end() ) { // must check if existing record inconsistent
			if ( c2p->second != pdb ) {

				// three cases:
				//  (1) replace an existing empty record
				//  (2) it's an unneeded empty record, so continue
				//  (3) there is an actual problem
				if ( pdb != PDBInfo::empty_record() && c2p->second == PDBInfo::empty_record() ) {
					// replace the record
					c2p->second = pdb;
				} else if ( pdb == PDBInfo::empty_record() ) {
					continue; // skip empty record
				} else {
					// something is inconsistent
					//fd: tone this message down a bit
					TR.Warning << "conf2pdb_chain(): chain mapping inconsistent, returning default map p 1=A, 2=B, ... " << std::endl;
					TR.Warning << "existing " << c2p->first << " -> " << c2p->second << "  |  ";
					TR.Warning << "new " << conf << " -> " << pdb << std::endl;
					conf2pdb_chain_default_map(pose,conf2pdb);
					return conf2pdb;
				}
			}
		} else { // record doesn't exist yet
			conf2pdb[ conf ] = pdb;
		}

	} // foreach residue

	debug_assert( conf2pdb.size() == pose.conformation().num_chains() );
	return conf2pdb;
}

utility::vector1< core::Size > get_chains( core::pose::Pose const & pose ) {

	utility::vector1< core::Size > chains;
	for ( core::Size ii(1); ii <= pose.num_chains(); ++ii ) {
		chains.push_back( ii );
	}
	return chains;
}

core::Size chain_end_res( Pose const & pose, core::Size const chain ) {

	// check whether chain exists in pose
	if ( ! has_chain( chain, pose ) ) {
		TR << "??? Chain should be " << chain << "???" << std::endl;
		utility_exit_with_message("Cannot get chain end residue for chain that doesn't exist. Quitting");
	}

	return pose.conformation().chain_end( chain );
}

/// @brief compute last residue numbers of all chains
utility::vector1< core::Size > chain_end_res( Pose const & pose ) {

	utility::vector1< core::Size > endings( pose.conformation().chain_endings() );
	endings.push_back( pose.size() ); // Last residue is implicit in chain_endings()
	return endings;
}


/// @brief Compute uniq chains in a complex
/// @detail Returns a vector of pose length with true/false of uniq chain
utility::vector1< bool > compute_unique_chains( Pose & pose ) {

	// initilize vector of pose length with false
	utility::vector1< bool > uniq( pose.size(), true );

	// get chains, initialize uniq chains and vector of chain sequences
	utility::vector1< core::Size > chains( get_chains( pose ) );
	utility::vector1< bool > uniq_chains( chains.size(), true );
	utility::vector1< std::string > sequences;

	// get sequences of all chains into a vector of strings
	for ( core::Size i = 1; i <= chains.size(); ++i ) {
		std::string seq( pose.chain_sequence( i ) );
		sequences.push_back( seq );
	}

	// compare sequences with respect to each other
	// this is the simplest and dirties way to compute whether the chains
	//  are similar, it does NOT do a proper sequence alignment
	// the assumptions are:
	// - if the sequences have different lengths, they are different
	//   (this obviously goes wrong if there are single residue insertions
	//   or deletions while the rest of the sequences are the same!)
	// - if the sequences have the same length and the sequence identity is
	//   above 95%, then they are "not unique"

	// go through vectors of chain sequences to compare them
	for ( core::Size i = 1; i <= sequences.size(); ++i ) {

		std::string const & seq1 = sequences[ i ];

		// no double counting
		for ( core::Size j = i+1; j <= sequences.size(); ++j ) {

			std::string const & seq2 = sequences[ j ];

			// make sure that the sequences are of same length
			core::Size num_ident_res( 0 );
			if ( seq1.size() == seq2.size() ) {

				TR << "sequences " << i << " and " << j << " have same length" << std::endl;

				// go through sequence, std::string indexing from 0
				for ( core::Size k = 0; k <= seq1.size(); ++k ) {
					if ( seq1[ k ] == seq2[ k ] ) {
						++num_ident_res;
					}
				}

				// compute sequence identity
				core::Real seqid = num_ident_res / seq1.size();
				TR << "seqid " << seqid << std::endl;

				// if sequence identity >= 95%, add a true to the uniq_chains
				// vector
				if ( seqid >= 0.95 ) {
					TR << "adding a false to chain " << i << std::endl;
					uniq_chains[ i ] = false;
				}

			} // same length sequences
		} // iterate over chains vector
	} // iterate over chains vector

	// go through uniq chains vector
	for ( core::Size i = 1; i <= uniq_chains.size(); ++i ) {

		// go through residues
		for ( core::Size j = 1; j <= pose.size(); ++j ) {

			// if residue belongs to uniq chain, set residue to true in uniq seq vector
			if ( uniq_chains[ i ] == false && chains[ i ] == static_cast< core::Size >( pose.chain( j ) ) ) {
				uniq[ j ] = false;
			}
		}
	}

	return uniq;

} // compute unique chains

/// @brief renumber PDBInfo based on Conformation chains; each chain starts from 1
/// @param[in,out] pose The Pose to modify.
/// @param[in] fix_chains If true, the procedure will attempt to fix any empty record
///  characters it finds in the PDBInfo. (default true)
/// @param[in] start_from_existing_numbering If true, will attempt to start each
///  chain from the existing numbering in the PDBInfo.  E.g. if the first residue
///  of chain 2 in the Conformation is 27, then the renumbering of the chain in
///  PDBInfo will start from 27. (default true)
/// @param[in] keep_insertion_codes If true, will maintain insertion codes and
///  will not increment the pdb residue numbering for those residues.  This means
///  new numbering with insertion codes will only reflect properly if the
///  old numbering included the base numbering of the insertion code residues,
///  i.e. 100 100A 100B and not just 100A 100B (with 100 never appearing).
///  (default false)
/// @param[in] rotate_chain_ids If true, allows support for more than 26 pdb chains
///  by rotating [A,Z] continuously.  WARNING: This will break the assumption
///  made by the PDBPoseMap that each pdb chain id is unique, so make sure you
///  are not using the PDBPoseMap feature downstream in your code path without
///  corrections! (default false)
/// @remarks If fixing chains and there is only one chain and the PDBInfo exists
///  but all records are marked as empty, will renumber and set the PDBInfo chain
///  to 'A'.
/// @return true if renumbering successful, false otherwise
bool renumber_pdbinfo_based_on_conf_chains(
	core::pose::Pose & pose,
	bool fix_chains,
	bool const start_from_existing_numbering,
	bool const keep_insertion_codes,
	bool const rotate_chain_ids
) {
	using core::Size;
	using core::pose::PDBInfo;
	typedef std::map< core::Size, char > Conf2PDB;

	if ( !pose.pdb_info().get() ) {
		TR.Warning << "renumber_pdbinfo_based_on_conf_chains(): no PDBInfo, returning" << std::endl;
		return false;
	}

	Conf2PDB conf2pdb = conf2pdb_chain( pose );

	if ( fix_chains ) {
		if ( conf2pdb.empty() ) { // something is wrong with chain consistency
			TR.Warning << "renumber_pdbinfo_based_on_conf_chains(): Request to fix PDBInfo chains, but ";
			TR.Warning << "chain mapping is inconsistent, so that step will be skipped." << std::endl;
			fix_chains = false;
		} else { // Try to fill in any empty record characters.

			// two different schemes: rotating and fixed length
			// WARNING: Rotating will break assumption of unique chain ids
			// inside PDBPoseMap, so make sure you are not using the PDBPoseMap
			// feature after calling this function without correcting!
			// First either remove or rotate any existing chains to the end of
			// the list.
			std::string letters( utility::UPPERCASE_LETTERS );
			for ( auto & i : conf2pdb ) {
				if ( i.second == PDBInfo::empty_record() )  continue;

				std::string::size_type const j = letters.find( i.second );
				if ( j == std::string::npos )  continue;

				if ( rotate_chain_ids ) { // rotating
					letters.push_back( letters.at( j ) );
				}

				letters.erase( j, 1 );
			}

			// Now fill in empty records.
			Size lidx = 0;
			for ( auto & i : conf2pdb ) {
				if ( i.second != PDBInfo::empty_record() )  continue;

				if ( rotate_chain_ids ) { // rotating
					i.second = letters.at( lidx % letters.size() );
				} else { // fixed length
					runtime_assert( lidx < letters.size() );
					i.second = letters.at( lidx );
				}
				++lidx;
			}

		} // if conf2pdb.empty()
	} // if fix_chains

	PDBInfo & pdbinfo = *pose.pdb_info();

	// grab all the chain endings
	utility::vector1< Size > chain_endings = pose.conformation().chain_endings();
	chain_endings.push_back( pose.size() ); // add the last res, which is not in the list

	Size res = 1;
	for ( Size const chain_end : chain_endings ) {
		int pdb_res = 0; // new chain, so reset pdb_res counter
		if ( start_from_existing_numbering && pdbinfo.chain( res ) != PDBInfo::empty_record() ) {
			pdb_res = pdbinfo.number( res ) - 1;
		}

		// find the corresponding pdb chain
		Conf2PDB::const_iterator c2p = conf2pdb.find( pose.chain( chain_end ) );
		debug_assert( ( fix_chains && c2p != conf2pdb.end() ) || !fix_chains ); // otherwise something's very wrong

		for ( ; res <= chain_end; ++res ) {
			// handle the pdb chain only if necessary
			char chain = pdbinfo.chain( res );
			if ( pdbinfo.chain( res ) == PDBInfo::empty_record() && fix_chains && c2p != conf2pdb.end() ) {
				chain = c2p->second;
			}

			// If keeping insertion codes, increment pdb_res counter only if
			// no insertion code or we're at position 1, in case there's an
			// insertion code at 1.
			char icode = pdbinfo.icode( res );
			if ( keep_insertion_codes && ( pdbinfo.icode( res ) == ' ' || res == 1 ) ) {
				++pdb_res;
			} else if ( !keep_insertion_codes ) { // always increment and clear insertion code
				icode = ' ';
				++pdb_res;
			}

			// The new pdb info for this residue must be setup in one shot.
			// The way we're redoing the info in this function can cause the
			// pdb2pose map to become out-of-sync if we attempt to make the
			// changes by separately calling chain(), icode(), and number().
			pdbinfo.set_resinfo( res, chain, pdb_res, icode );
		}

	} // foreach chain ending

	// no point updating pdb_info if it's just thrown away
	pose.pdb_info()->obsolete( false );

	debug_assert( res == pose.size() + 1 );

	return true;
}

bool
has_chain(std::string const & chain, core::pose::Pose const & pose){
	debug_assert(chain.size()==1);// chain is one char
	char chain_char= chain[0];
	return has_chain(chain_char, pose);
}

bool
has_chain(char const & chain, core::pose::Pose const & pose){
	PDBInfoCOP pdb_info( pose.pdb_info() );
	if ( ! pdb_info ) {
		TR.Warning << "In has_chain(): Pose doesn't have a PDBInfo object - making a temporary default one." << std::endl;
		pdb_info = PDBInfoCOP( new PDBInfo( pose ) );
	}
	for ( core::Size ii=1; ii <= pose.size(); ++ii ) {
		if ( pdb_info->chain( ii ) == chain ) {
			return true;
		}
	}
	return false;
}

bool
has_chain(core::Size chain_id, core::pose::Pose const & pose){
	return chain_id >= 1 && chain_id <= pose.conformation().num_chains();
}

utility::vector1< core::Size >
get_chain_ids_from_chains(utility::vector1< std::string > const & chains, core::pose::Pose const & pose ) {
	utility::vector1< char > chains_char;
	for ( std::string const & chain: chains ) {
		debug_assert( chain.size() >= 1 );
		chains_char.push_back( chain[0] );
	}
	return get_chain_ids_from_chains( chains_char, pose );
}

utility::vector1< core::Size >
get_chain_ids_from_chains(utility::vector1< char > const & chains, core::pose::Pose const & pose ) {
	std::set< core::Size > chain_ids;

	PDBInfoCOP pdb_info( pose.pdb_info() );
	if ( ! pdb_info ) {
		TR.Warning << "In get_chain_ids_from_chains(): Pose doesn't have a PDBInfo object - making a temporary default one." << std::endl;
		pdb_info = PDBInfoCOP( new PDBInfo( pose ) );
	}

	for ( core::Size ii(1); ii <= pose.size(); ++ii ) {
		char chain( pdb_info->chain( ii ) );
		if ( std::find( chains.begin(), chains.end(), chain ) != chains.end() ) {
			core::Size chain_id( pose.chain(ii) );
			chain_ids.insert( chain_id );
		}
	}
	// Set is sorted, so the result will be sorted.
	return utility::vector1< core::Size >( chain_ids.begin(), chain_ids.end() );
}

utility::vector1<core::Size>
get_chain_ids_from_chain(std::string const & chain, core::pose::Pose const & pose){
	debug_assert(chain.size()==1);// chain is one char
	char chain_char= chain[0];
	return get_chain_ids_from_chain(chain_char, pose);
}

utility::vector1<core::Size>
get_chain_ids_from_chain(char const & chain, core::pose::Pose const & pose){
	utility::vector1< char > chains;
	chains.push_back( chain );
	return get_chain_ids_from_chains( chains, pose );
}

core::Size
get_chain_id_from_chain(std::string const & chain, core::pose::Pose const & pose){
	debug_assert(chain.size()==1);// chain is one char
	return get_chain_id_from_chain( chain[0], pose );
}

core::Size
get_chain_id_from_chain(char const & chain, core::pose::Pose const & pose){
	utility::vector1<core::Size> chain_ids = get_chain_ids_from_chain(chain, pose);
	if ( chain_ids.size() == 0 ) {
		throw utility::excn::EXCN_RangeError(" chain_id "+utility::to_string(chain)+" does not exist");
	} else if ( chain_ids.size() > 1 ) {
		throw utility::excn::EXCN_RangeError("chain_id "+utility::to_string(chain)+" represents more than one chain!");
	}
	return chain_ids[1];
}

char
get_chain_from_chain_id( core::Size const & chain_id, core::pose::Pose const & pose ) {
	debug_assert( chain_id <= pose.num_chains() );
	PDBInfoCOP pdb_info( pose.pdb_info() );
	if ( ! pdb_info ) {
		TR.Warning << "In res_in_chain(): Pose doesn't have a PDBInfo object - making a temporary default one." << std::endl;
		pdb_info = PDBInfoCOP( new PDBInfo( pose ) );
	}
	core::Size first_chaisize = pose.conformation().chain_begin( chain_id );
	return pose.pdb_info()->chain(first_chaisize);
}

std::set<core::Size>
get_jump_ids_from_chain_ids(std::set<core::Size> const & chain_ids, core::pose::Pose const & pose){
	std::set<core::Size> jump_ids;

	for ( core::Size jump_id=1; jump_id <= pose.num_jump(); jump_id++ ) {
		core::Size downstream_residue_id( pose.fold_tree().downstream_jump_residue(jump_id) );
		core::Size downstream_chain_id( pose.chain(downstream_residue_id) );
		if ( chain_ids.count( downstream_chain_id ) ) {
			jump_ids.insert( jump_id );
		}
	}

	return jump_ids;
}

core::Size
get_jump_id_from_chain_id( core::Size const & chain_id, const core::pose::Pose & pose ){
	if ( ! has_chain( chain_id, pose ) ) {
		utility_exit_with_message( "Pose does not have a chain " + utility::to_string(chain_id) );
	}
	std::set<core::Size> chain_ids{ chain_id }; // initializaiton list with one element
	std::set<core::Size> jumps( get_jump_ids_from_chain_ids( chain_ids, pose ) );

	if ( jumps.empty() ) {
		utility_exit_with_message("Chain '" + utility::to_string(chain_id) + "' is not directly built by any jump.");
	} else if ( jumps.size() > 1 ) {
		TR.Warning << "get_jump_id_from_chain_id(): Chain " << chain_id << " has multiple jumps, returning the lowest number one." << std::endl;
	}

	// std::set is sorted, so the lowest number jump will be the first entry.
	return *(jumps.begin());
}

utility::vector1<core::Size>
get_jump_ids_from_chain( char const & chain, core::pose::Pose const & pose ) {
	utility::vector1<core::Size> jump_ids;

	PDBInfoCOP pdb_info( pose.pdb_info() );
	if ( ! pdb_info ) {
		TR.Warning << "In get_jump_ids_from_chain(): Pose doesn't have a PDBInfo object - making a temporary default one." << std::endl;
		pdb_info = PDBInfoCOP( new PDBInfo( pose ) );
	}

	for ( core::Size jump_id=1; jump_id <= pose.num_jump(); jump_id++ ) {
		core::Size downstream_residue_id( pose.fold_tree().downstream_jump_residue(jump_id) );
		char downstream_chain( pdb_info->chain(downstream_residue_id) );
		if ( chain == downstream_chain ) {
			jump_ids.push_back( jump_id );
		}
	}

	return jump_ids; // Will be returned in sorted order
}

utility::vector1<core::Size>
get_jump_ids_from_chain(std::string const & chain, core::pose::Pose const & pose){
	debug_assert(chain.size()==1);// chain is one char
	char chain_char= chain[0];
	return get_jump_ids_from_chain(chain_char, pose);
}

core::Size
get_jump_id_from_chain( std::string const & chain, core::pose::Pose const & pose ) {
	debug_assert( chain.size() == 1 );
	return get_jump_id_from_chain( chain[0], pose );
}

core::Size
get_jump_id_from_chain( char const & chain, core::pose::Pose const & pose ) {
	utility::vector1<core::Size> jumps( get_jump_ids_from_chain( chain, pose ) );
	if ( jumps.empty() ) {
		// Explicit cast here so we can use operator+
		utility_exit_with_message(std::string("Chain '") + chain + "' is not directly built by any jump.");
	} else if ( jumps.size() > 1 ) {
		TR.Warning << "get_jump_id_from_chain(): Chain " << chain << " has multiple jumps, returning the lowest number one." << std::endl;
	}
	// Because of the way the return value from get_jump_ids_from_chain is filled, it should be sorted order.
	return jumps[1];
}

core::Size
get_chain_id_from_jump_id(core::Size const & jump_id, core::pose::Pose const & pose){
	debug_assert( jump_id <= pose.num_jump() );
	core::Size downstream_residue_id( pose.fold_tree().downstream_jump_residue(jump_id) );
	return pose.chain(downstream_residue_id);
}

char
get_chain_from_jump_id(core::Size const & jump_id, core::pose::Pose const & pose){
	debug_assert( jump_id <= pose.num_jump() );
	core::Size downstream_residue_id( pose.fold_tree().downstream_jump_residue(jump_id) );

	PDBInfoCOP pdb_info( pose.pdb_info() );
	if ( ! pdb_info ) {
		TR.Warning << "In get_chain_from_jump_id(): Pose doesn't have a PDBInfo object - making a temporary default one." << std::endl;
		pdb_info = PDBInfoCOP( new PDBInfo( pose ) );
	}
	return pdb_info->chain( downstream_residue_id );
}

utility::vector1<core::Size>
get_resnums_for_chain( core::pose::Pose const & pose, char chain ) {
	PDBInfoCOP pdb_info( pose.pdb_info() );
	if ( ! pdb_info ) {
		TR.Warning << "Attempted to find chain for Pose without annotated chain letters - making default chain letters." << std::endl;
		pdb_info = PDBInfoCOP( new PDBInfo( pose ) );
	}

	utility::vector1<core::Size> resnums_in_chain;
	for ( core::Size ii(1); ii <= pose.size(); ++ii ) {
		if ( pdb_info->chain( ii ) == chain ) {
			resnums_in_chain.push_back( ii );
		}
	}
	return resnums_in_chain;
}

utility::vector1<core::Size>
get_resnums_for_chain_id( core::pose::Pose const & pose, core::Size chain_id ) {
	debug_assert( chain_id <= pose.num_chains() );

	utility::vector1<core::Size> resnums_in_chain;
	core::Size begin( pose.conformation().chain_begin(chain_id) );
	core::Size const end( pose.conformation().chain_end(chain_id) );
	for ( ; begin <= end; ++begin ) {
		resnums_in_chain.push_back( begin );
	}
	return resnums_in_chain;
}

core::conformation::ResidueCOPs
get_chain_residues( core::pose::Pose const & pose, core::Size const chain_id ) {
	debug_assert( chain_id <= pose.num_chains() );
	core::Size begin( pose.conformation().chain_begin(chain_id) );
	core::Size const end( pose.conformation().chain_end(chain_id) );
	core::conformation::ResidueCOPs residues;
	for ( ; begin <= end; ++begin ) {
		residues.push_back( core::conformation::ResidueOP( new core::conformation::Residue(pose.residue(begin)) ) );
	}
	return residues;
}

core::conformation::ResidueCOPs
get_residues_from_chains(core::pose::Pose const & pose, utility::vector1<core::Size> const & chain_ids)
{
	core::conformation::ResidueCOPs chains_residue_pointers;
	for ( core::Size chain_id : chain_ids ) {
		core::conformation::ResidueCOPs chain_residue_pointers = core::pose::get_chain_residues(pose, chain_id);
		for ( core::conformation::ResidueCOP chain_residue_pointer : chain_residue_pointers ) {
			chains_residue_pointers.push_back(chain_residue_pointer);
		}
	}
	return chains_residue_pointers;
}

bool res_in_chain( core::pose::Pose const & pose, core::Size resnum, std::string const & chain ) {
	debug_assert( chain.size() == 1 );
	PDBInfoCOP pdb_info( pose.pdb_info() );
	if ( ! pdb_info ) {
		TR.Warning << "In res_in_chain(): Pose doesn't have a PDBInfo object - making a temporary default one." << std::endl;
		pdb_info = PDBInfoCOP( new PDBInfo( pose ) );
	}
	return pdb_info->chain( resnum ) == chain[0];
}

core::Size get_hash_from_chain(char const & chain, core::pose::Pose const & pose, std::string const & extra_label)
{
	PDBInfoCOP pdb_info( pose.pdb_info() );
	if ( ! pdb_info ) {
		TR.Warning << "In get_hash_from_chain: Pose doesn't have a PDBInfo object - making a temporary default one." << std::endl;
		pdb_info = PDBInfoCOP( new PDBInfo( pose ) );
	}

	std::size_t hash = 0;

	for ( core::Size res_num = 1; res_num <= pose.size(); ++res_num ) {
		if ( pdb_info->chain( res_num ) != chain ) {
			continue;
		}
		core::Size natoms = pose.conformation().residue(res_num).natoms();
		for ( core::Size atom_num = 1; atom_num <= natoms; ++atom_num ) {
			id::AtomID atom_id(atom_num,res_num);
			PointPosition current_xyz = pose.conformation().xyz(atom_id);
			boost::hash_combine(hash,current_xyz);
		}
	}
	if ( ! extra_label.empty() ) {
		boost::hash_combine(hash,extra_label);
	}

	return hash;
}

core::Size get_hash_excluding_chain(char const & chain, core::pose::Pose const & pose, std::string const & extra_label)
{
	PDBInfoCOP pdb_info( pose.pdb_info() );
	if ( ! pdb_info ) {
		TR.Warning << "In get_hash_excluding_chain: Pose doesn't have a PDBInfo object - making a temporary default one." << std::endl;
		pdb_info = PDBInfoCOP( new PDBInfo( pose ) );
	}

	std::size_t hash = 0;

	for ( core::Size res_num = 1; res_num <= pose.size(); ++res_num ) {
		if ( pdb_info->chain( res_num ) == chain ) {
			continue;
		}
		core::Size natoms = pose.conformation().residue(res_num).natoms();
		for ( core::Size atom_num = 1; atom_num <= natoms; ++atom_num ) {
			id::AtomID atom_id(atom_num,res_num);
			PointPosition current_xyz = pose.conformation().xyz(atom_id);
			boost::hash_combine(hash,current_xyz);
		}
	}
	if ( ! extra_label.empty() ) {
		boost::hash_combine(hash,extra_label);
	}

	return hash;
}

std::string get_sha1_hash_from_chain(char const & chain, core::pose::Pose const & pose, std::string const & extra_label)
{
	utility::vector1< std::string > chains;
	chains.push_back( utility::to_string( chain ) );
	return get_sha1_hash_from_chains(chains, pose, extra_label);
}

std::string
get_sha1_hash_from_chains(utility::vector1< std::string > const & chains, core::pose::Pose const & pose, std::string const & extra_label) {
	PDBInfoCOP pdb_info( pose.pdb_info() );
	if ( ! pdb_info ) {
		TR.Warning << "In get_sha1_hash_from_chain: Pose doesn't have a PDBInfo object - making a temporary default one." << std::endl;
		pdb_info = PDBInfoCOP( new PDBInfo( pose ) );
	}

	std::stringstream coord_stream;

	for ( core::Size res_num = 1; res_num <= pose.size(); ++res_num ) {
		if ( ! chains.contains( utility::to_string( pdb_info->chain( res_num ) ) ) ) {
			continue;
		}
		core::Size natoms = pose.conformation().residue(res_num).natoms();
		for ( core::Size atom_num = 1; atom_num <= natoms; ++atom_num ) {
			id::AtomID atom_id(atom_num,res_num);
			PointPosition current_xyz = pose.conformation().xyz(atom_id);
			coord_stream << numeric::truncate_and_serialize_xyz_vector(current_xyz,3);
		}
	}
	if ( ! extra_label.empty() ) {
		coord_stream << extra_label;
	}
	return utility::string_to_sha1(coord_stream.str());

}

std::string get_sha1_hash_excluding_chain(char const & chain, core::pose::Pose const & pose, std::string const & extra_label)
{
	utility::vector1< std::string > chains;
	chains.push_back( utility::to_string( chain ) );
	return get_sha1_hash_excluding_chains(chains, pose, extra_label);
}

std::string
get_sha1_hash_excluding_chains(utility::vector1< std::string > const & chains, core::pose::Pose const & pose, std::string const & extra_label)
{
	PDBInfoCOP pdb_info( pose.pdb_info() );
	if ( ! pdb_info ) {
		TR.Warning << "In get_sha1_hash_excluding_chain: Pose doesn't have a PDBInfo object - making a temporary default one." << std::endl;
		pdb_info = PDBInfoCOP( new PDBInfo( pose ) );
	}

	std::stringstream coord_stream;

	for ( core::Size res_num = 1; res_num <= pose.size(); ++res_num ) {
		if ( chains.contains( utility::to_string( pdb_info->chain( res_num ) ) ) ) {
			continue;
		}
		core::Size natoms = pose.conformation().residue(res_num).natoms();
		for ( core::Size atom_num = 1; atom_num <= natoms; ++atom_num ) {
			id::AtomID atom_id(atom_num,res_num);
			PointPosition current_xyz = pose.conformation().xyz(atom_id);
			coord_stream << numeric::truncate_and_serialize_xyz_vector(current_xyz,5);
		}
	}
	if ( ! extra_label.empty() ) {
		coord_stream << extra_label;
	}
	return utility::string_to_sha1(coord_stream.str());
}

} // pose
} // core
