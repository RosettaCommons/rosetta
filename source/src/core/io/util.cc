// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/io/util.cc
/// @brief Util functions for Input and Output.  Very general IO should go to utility/io.
///   These should be related to core in a deep way or not able to be called from utility.
///
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com), XRW 2016 Team

// Unit headers
#include <core/io/util.hh>

// Package headers
#include <core/io/StructFileRep.hh>
#include <core/io/StructFileRepOptions.hh>
#include <core/io/pose_to_sfr/PoseToStructFileRepConverter.hh>
#include <core/io/pose_from_sfr/PoseFromSFRBuilder.hh>
#include <core/io/ResidueInformation.hh>

// Project headers
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeFinder.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfoManager.hh>
#include <core/scoring/Energies.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>


// Basic headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableString.hh>
#include <basic/datacache/CacheableStringFloatMap.hh>
#include <basic/datacache/CacheableStringMap.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/io/ozstream.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <utility/string_util.hh>

// External headers
#include <ObjexxFCL/format.hh>

#include <map>

namespace core {
namespace io {

// Tracer instance for this file
static basic::Tracer TR( "core.io.util" );


using namespace ObjexxFCL::format; // AUTO USING NS



/// @brief Write Pose energies information into a string and return it.
/// @details Added during the 2016 Chemical XRW.
/// @author Vikram K. Mulligan (vmullig@uw.edu) Jared Adolf-Bryfogle (jadolfbr@gmail.com)
std::string pose_energies_from_sfr(
	StructFileRep const & sfr
) {
	std::stringstream out;
	pose_energies_from_sfr(sfr, out);
	return out.str();
}

void pose_energies_from_sfr(
	StructFileRep const & sfr,
	std::stringstream & out
)
{
	using namespace core::io::pose_to_sfr;

	// This version is formatted for easy parsing by R, Excel, etc.

	utility::vector1< std::string > const & score_names = sfr.score_table_labels();
	utility::vector1< std::vector< std::string > > const & score_lines = sfr.score_table_lines();

	if ( score_names.size() == 0 || score_lines.size() == 0 ) return; //Was not extracted!

	out << "# All scores below are weighted scores, not raw scores.\n";

	if ( ! sfr.score_table_filename().empty() ) {
		out << "#BEGIN_POSE_ENERGIES_TABLE " << sfr.score_table_filename() << std::endl;
	} else {
		out << "#BEGIN_POSE_ENERGIES_TABLE " << std::endl;
	}

	out << "label";

	for ( std::string const & score_name : score_names ) {
		out << " " << score_name;
	}
	out << "\n";
	out << "weights";
	utility::vector1< core::Real > const & score_weights = sfr.score_table_weights();

	for ( core::Real const weight : score_weights ) {
		out << " " << weight;
	}
	out << " NA\n";


	for ( std::vector<std::string> const & score_line : score_lines ) {
		std::string line = "";
		for ( std::string const & column : score_line ) {
			line = line+" "+column;
		}
		line = utility::strip(line);
		out << line << "\n";
	}
	if ( ! sfr.score_table_filename().empty() ) {
		out << "#END_POSE_ENERGIES_TABLE " << sfr.score_table_filename() << std::endl;
	} else {
		out << "#END_POSE_ENERGIES_TABLE " << std::endl;
	}
}


/// @brief Write Pose energies information into a string and return it.
/// @details Added during the 2016 Chemical XRW.
/// @author Vikram K. Mulligan (vmullig@uw.edu) + Jared Adolf-Bryfogle (jadolfbr@gmail.com)
std::string pose_data_cache_from_sfr(
	StructFileRep const & sfr
) {
	std::stringstream out;
	pose_data_cache_from_sfr(sfr, out);
	return out.str();
}

void pose_data_cache_from_sfr(
	StructFileRep const & sfr,
	std::stringstream & out
)
{

	//If either of these are empty, will not do anything.
	std::map< std::string, std::string > const & string_data = sfr.pose_cache_string_data();
	std::map< std::string,     float   > const & float_data =  sfr.pose_cache_float_data();

	// ARBITRARY_STRING_DATA
	for ( auto const & it : string_data ) {
		//TR << it->first << " " << it->second << std::endl;
		out << it.first << " " << it.second << std::endl;
	}

	// ARBITRARY_FLOAT_DATA
	for ( auto const & it : float_data ) {
		//TR << it->first << " " << it->second << std::endl;
		out << it.first << " " << it.second << std::endl;
	}
}

void
pose_from_pose(
	pose::Pose & new_pose,
	pose::Pose const & old_pose,
	utility::vector1< core::Size > const & residue_indices
){
	StructFileRepOptions options;
	pose_from_pose( new_pose, old_pose, residue_indices, options );
}


void
pose_from_pose(
	pose::Pose & new_pose,
	pose::Pose const & old_pose,
	utility::vector1< core::Size > const & residue_indices,
	StructFileRepOptions const & options
){
	using namespace chemical;
	ResidueTypeSetCOP residue_set(
		ChemicalManager::get_instance()->residue_type_set( FA_STANDARD )
	);
	pose_from_pose( new_pose, old_pose, *residue_set,  residue_indices, options);
}


void
pose_from_pose(
	pose::Pose & new_pose,
	pose::Pose const & old_pose,
	chemical::ResidueTypeSet const & residue_set,
	utility::vector1< core::Size > const & residue_indices
){
	StructFileRepOptions options;
	pose_from_pose( new_pose, old_pose, residue_set, residue_indices, options );
}


/// Creates a subpose from a pose, to include only certain
/// residues, using StructFileRep::init_from_pose() to construct the
/// pose, and build_pose_as_is1() to construct the pose
/// with the given options.
void
pose_from_pose(
	pose::Pose & new_pose,
	pose::Pose const & old_pose,
	chemical::ResidueTypeSet const & residue_set,
	utility::vector1< core::Size > const & residue_indices,
	StructFileRepOptions const & options
){
	core::io::pose_to_sfr::PoseToStructFileRepConverter converter;
	converter.init_from_pose( old_pose, residue_indices );
	pose_from_sfr::PoseFromSFRBuilder builder( residue_set.get_self_ptr(), options );
	builder.build_pose( *converter.sfr(), new_pose );

}


////////////////////////////////////////////////////////////////////
///  Glycan IO

utility::vector1< core::Size >
fix_glycan_order( utility::vector1< core::io::ResidueInformation > & rinfos,
	utility::vector1< core::Size > const & glycan_positions,
	StructFileRepOptions const & options,
	std::map< std::string, std::map< std::string, std::pair< std::string, std::string > > > const & known_links )

{
	utility::vector1< core::Size > chain_ends;
	// Note that we shouldn't re-use the found glycan link map, as the ordering can/will change
	utility::vector1< core::Size > correct_order( find_carbohydrate_order( rinfos, glycan_positions, chain_ends,
		determine_glycan_links( rinfos, options ), known_links ) );
	reorder_glycan_residues( rinfos, correct_order, glycan_positions );

	utility::vector1< core::Size > reorganized_ends;
	for ( core::Size orig_end: chain_ends ) {
		reorganized_ends.push_back( glycan_positions[ correct_order.index_of(orig_end) ] );
	}
	return reorganized_ends;
}

utility::vector1< core::Size >
find_carbohydrate_order( utility::vector1< core::io::ResidueInformation > const & rinfos,
	utility::vector1< core::Size > const & glycan_positions,
	utility::vector1< core::Size > & chain_ends, // return-by-reference for (non-reducing) end sugars
	// map of anomeric positions to where they are connected to
	std::map< std::pair< core::Size, std::string >, std::pair< core::Size, std::string > > const & link_map,
	std::map< std::string, std::map< std::string, std::pair< std::string, std::string > > > const & known_links )
{
	std::set< core::Size > roots; // Sorted!
	for ( core::Size resi: glycan_positions ) {
		bool is_root = true;
		for ( auto const & link: link_map ) {
			// Not a root if we have an anomeric link to a glycan residue
			if ( link.first.first == resi &&
					core::chemical::carbohydrates::CarbohydrateInfoManager::is_valid_sugar_code( rinfos[ link.second.first ].rosetta_resName() )
					) {
				is_root = false;
				break;
			}
		}
		//Add all the known link records
		std::string rinfo1 = rinfos[resi].resid();
		for ( auto const & link1: known_links ) {
			for ( auto const & link2: link1.second ) {
				if ( (link1.first == rinfo1 && link2.first == " C1 " ) || (link2.second.first == rinfo1 && link2.second.second == " C1 ") ) {
					is_root = false;
					break;
				}
			}
		}
		if ( is_root ) { roots.insert( resi ); }
	}

	std::set< core::Size > addressed;
	utility::vector1< core::Size > full_order;
	// Unfortunately we need to redo the data structure a bit
	std::map< core::Size, std::map< std::string, std::pair< core::Size, std::string > > > connectivity;
	for ( auto const & link: link_map ) {
		connectivity[ link.second.first ][ link.second.second ] = link.first;
	}
	//this will override any auto detected links with the links for the input
	for ( auto const & link: known_links ) {
		core::Size res1 = 1;
		for ( core::Size i=1; i<=rinfos.size(); i++ ) {
			if ( link.first == rinfos[i].resid() ) {
				res1 = i;
			}
		}
		for ( auto const & link2: link.second ) {
			core::Size res2 = 1;
			for ( core::Size i=1; i<=rinfos.size(); i++ ) {
				if ( link2.second.first == rinfos[i].resid() ) {
					res2 = i;
				}
			}
			//Only add the connection lower->upper
			if ( res1 < res2 ) {
				connectivity[ res1 ][ link2.first ] = std::make_pair(res2,link2.second.second);
			}
		}
	}

	// Now find the chains, prefering the roots with lower numbers
	for ( core::Size root: roots ) {
		full_order.append( find_carbohydrate_subbranch_order( root, chain_ends, connectivity, addressed ) );
	}

	// Okay, now address residues which aren't findable from a root - e.g. cycles, loners
	for ( core::Size resi: glycan_positions ) {
		if ( addressed.count( resi ) == 0 ) {
			full_order.append( find_carbohydrate_subbranch_order( resi, chain_ends, connectivity, addressed ) );
		}
	}

	return full_order;
}

utility::vector1< core::Size >
find_carbohydrate_subbranch_order( core::Size current_res,
	utility::vector1< core::Size > & chain_ends, // return-by-reference for (non-reducing) end sugars
	// Nested map positions & (non-anomeric) positions to attached anomeric position
	std::map< core::Size, std::map< std::string, std::pair< core::Size, std::string > > > const & connectivity,
	std::set< core::Size > & addressed )
{

	utility::vector1< core::Size > subtree_order;
	if ( addressed.count( current_res ) ) {
		// Cycle - we don't add ourself to the order, we're in another branch.
		return subtree_order; // empty - not even self
	}
	subtree_order.push_back( current_res );
	addressed.insert( current_res );

	if ( connectivity.count( current_res ) == 0 ) {
		// No children
		chain_ends.push_back( current_res );
		return subtree_order;
	}

	// Currently going in atom sorted order, rather than branch size order
	for ( auto const & per_atom_name_pair: connectivity.at( current_res ) ) { // std::map sorts keys
		core::Size child = per_atom_name_pair.second.first;
		utility::vector1< core::Size > const & subbranch( find_carbohydrate_subbranch_order( child, chain_ends, connectivity, addressed ) );
		subtree_order.append( subbranch );
	}

	return subtree_order;
}

/// @brief Bring glycans into the correct order, which corresponds to connectivity of ech glycan tree
/// This requires reordering rinfos
void reorder_glycan_residues( utility::vector1< core::io::ResidueInformation >& rinfos,
	utility::vector1< core::Size > const & correct_order,
	utility::vector1< core::Size > const & glycan_positions  )
{
	utility::vector1< core::io::ResidueInformation > rinfos_orig = rinfos;
	using namespace core::chemical;
	TR << "Automatic glycan connection is activated." << std::endl;
	TR << "Start reordering residues." << std::endl;
	TR << "Corrected glycan residue order (internal numbering): " << correct_order << std::endl;
	TR << std::endl;
	if ( correct_order.size() != glycan_positions.size() ) {
		std::stringstream err_msg;
		err_msg << "Not all glycans have been detected by the neighbor search. Glycans in pdb: " << glycan_positions.size() <<  ", detected in neighbor search: " <<  correct_order.size() << std::endl;
		utility_exit_with_message( err_msg.str() );
	}
	size_t indx = 0;
	for ( core::Size i : glycan_positions ) {
		++indx;
		// rinfos
		Size corrected_pos = correct_order[ indx ];
		ResidueInformation new_res = rinfos_orig[ corrected_pos ];
		ResidueInformation current_res = rinfos[ i ];
		// swap in the correct ResidueInfo objects in rinfos
		TR.Trace << "Swapping " << current_res.resSeq() << "( " << i << " ) against " << new_res.resSeq() << " ( " << corrected_pos << " ) " << std::endl;
		rinfos[ i ] = new_res;
		// name3s
	}
}

std::map< std::pair< core::Size, std::string >, std::pair< core::Size, std::string > >
determine_glycan_links( utility::vector1< core::io::ResidueInformation > const & rinfos,
	StructFileRepOptions const & options
) {

	// Linkage map of anomeric carbons to the carbon that it is (nominally) attached to.
	std::map< std::pair< core::Size, std::string >, std::pair< core::Size, std::string > > linkage_map;

	core::Real max_cutoff = options.max_bond_length();
	core::Real min_cutoff = options.min_bond_length();

	for ( Size ii(1); ii <= rinfos.size(); ++ii ) {
		std::string const & res_name( rinfos[ii].rosetta_resName() );
		if ( ! core::chemical::carbohydrates::CarbohydrateInfoManager::is_valid_sugar_code( res_name ) ) { continue; }
		char anomeric = core::chemical::carbohydrates::CarbohydrateInfoManager::anomeric_position_from_code( res_name );
		std::string anomeric_name = " C"; // Need to break this out as we're dealing with C-string & char
		anomeric_name += anomeric;
		anomeric_name += " ";
		if ( ! rinfos[ ii ].xyz().count( anomeric_name ) ) {
			TR << "Sugar residue " << rinfos[ii].resid() << " doesn't have coordinates for the anomeric carbon " << anomeric_name << std::endl;
			continue;
		}
		Vector const & anomeric_coords( rinfos[ ii ].xyz().at( anomeric_name ) );
		std::pair< core::Size, std::string > anomeric_pair( ii, anomeric_name );

		// There's an intrinsic directionality here - ii is the anomeric position, so we need to do the full matrix
		for ( Size jj(1); jj <= rinfos.size(); ++jj ) { // There's an intrinsic directionality
			if ( ii == jj ) { continue; } // Don't do self links
			if ( linkage_map.count( anomeric_pair ) ) { break; } // Only one link is needed.
			std::map< std::string, Vector > XYZs(rinfos[ jj ].xyz());
			for ( auto const &pair: XYZs ) {
				if ( rinfos[jj].rosetta_resName() == "ASN" && pair.first != " ND2" ) continue;
				// TODO: only check oxygens ??
				Vector const & pos( pair.second );
				core::Real distance = pos.distance( anomeric_coords );
				//TR.Trace << " Distance to " << pair.first << ": " << distance << std::endl;
				if ( (min_cutoff <= distance) && (distance <= max_cutoff) ) {
					linkage_map[ anomeric_pair ] = make_pair( jj, pair.first );
				}
			}
		}
	}
	return linkage_map;
}

std::map< std::string, std::map< std::string, std::pair< std::string, std::string > > >
explicit_links_from_sfr_linkage( std::map< std::string, utility::vector1< LinkInformation > > const & link_map,
	utility::vector1< core::io::ResidueInformation > const & rinfos )
{
	std::set< std::string > known_resid;
	for ( core::io::ResidueInformation const & resinfo : rinfos ) {
		known_resid.insert( resinfo.resid() );
	}

	std::map< std::string, std::map< std::string, std::pair< std::string, std::string > > > bi_map;
	for ( auto const & lm_pair : link_map ) {
		for ( LinkInformation const & link_info : lm_pair.second ) {
			if ( known_resid.count( link_info.resID1 ) == 0 || known_resid.count( link_info.resID2 ) == 0 ) {
				TR << "Link between " << link_info.resID1 << " and " << link_info.resID2 << " is ill-formed - one/both residues don't exist!." << std::endl;
				continue;
			}
			bi_map[ link_info.resID1 ][ link_info.name1 ] = make_pair( link_info.resID2, link_info.name2 );
			bi_map[ link_info.resID2 ][ link_info.name2 ] = make_pair( link_info.resID1, link_info.name1 );
		}
	}
	return bi_map;
}

void
add_glycan_links_to_map(
	std::map< std::string, std::map< std::string, std::pair< std::string, std::string > > > & known_links,
	std::map< std::pair< core::Size, std::string >, std::pair< core::Size, std::string > > const & link_map,
	utility::vector1< core::io::ResidueInformation > const & rinfos )
{

	if ( !basic::options::option[ basic::options::OptionKeys::in::maintain_links ].value() ) {
		known_links.clear(); // Ignore existing links, as they may be garbage -- we may want to reconsider this
	}

	std::map< std::string, core::Size > resid_to_pos;
	for ( core::Size ii(1); ii <= rinfos.size(); ++ii ) {
		resid_to_pos[ rinfos[ii].resid() ] = ii;
	}

	for ( auto const & pair: link_map ) {
		std::string rinfo1 = rinfos[pair.first.first].resid();
		std::string atom1 = pair.first.second;
		std::string rinfo2 = rinfos[pair.second.first].resid();
		std::string atom2 = pair.second.second;

		// Don't overwrite explicit links
		if ( known_links[ rinfo1 ].count( atom1 ) == 0 ) {
			known_links[ rinfo1 ][ atom1 ] = make_pair( rinfo2, atom2 );
		}
		if ( known_links[ rinfo2 ].count( atom2 ) == 0 ) {
			known_links[ rinfo2 ][ atom2 ] = make_pair( rinfo1, atom1 );
		}
	}
}


} //core
} //io
