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
#include <core/io/NomenclatureManager.hh>
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

namespace core {
namespace io {

// Tracer instance for this file
static THREAD_LOCAL basic::Tracer TR( "core.io.util" );


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

/// @details When LINK-records in a PDB file are missing, the connections of glycan residues have to be traced
/// using the xyz coordines directly. find_downstream neighbor looks for a) glycans connected to ansparagine and
/// b) glycans connected to other glycans. The fuction is used by "find_mainchain_connection()" and "find_branch_points()"
/// to get the connections right. Specifically, branching information and the two booleans "same_chain_prev" and
/// the current residue.

/// @brief Identify residues that are branch points (have more than one downstream neighbor)
/// One eneighbpr will be considered as the mainchain continuations. All others are branches.
/// @author Sebastian Rämisch, raemisch@scripps.edu
void
find_branch_points( Size const & seqpos, chemical::ResidueTypeCOP & RT, bool & is_branch_point, utility::vector1< std::string > & branch_points_on_this_residue, utility::vector1< std::string > const & rosetta_residue_name3s, Size mainchain_neighbor, utility::vector1< core::io::ResidueInformation > & rinfos, core::io::StructFileRep::Strings & branch_lower_termini, utility::vector1< Size >& glycan_positions, StructFileRepOptions const & options ) {
	std::string const name3 = rosetta_residue_name3s[ seqpos ];
	TR.Trace << "Find branch point at this residue" << std::endl;
	TR.Trace << "It's name3 is >" << name3 << "<" << std::endl;
	bool is_carbohydrate = RT->is_carbohydrate();
	if ( name3 == "ASN" || name3 == "SER" || name3 == "THR" || is_carbohydrate ) {
		TR.Trace << "It's Asn/Ser/Thr or carbohydrate" << std::endl;
		Size n_branches = 0;
		std::pair< core::Size,std::string  > neighbor = std::make_pair(0,"C0"); // <residue number,atom name>
		std::map< std::string, Vector > const& XYZs(rinfos[ seqpos ].xyz());
		//for ( Size atom = 1; atom <= RT->natoms(); ++atom ) {
		for ( std::map< std::string, Vector >::const_iterator XYZs_it = XYZs.begin(); XYZs_it != XYZs.end(); ++XYZs_it ) {
			//std::string const & atom_name = RT->atom_name(atom);
			std::string const & atom_name = XYZs_it->first;
			if ( ( is_carbohydrate && atom_name.find(" O") != std::string::npos ) || ( atom_name == " ND2" || atom_name == " OG1" || atom_name == " OG ") ) {
				TR.Trace << "Check for branch point at atom: " << atom_name << std::endl;
				Vector const xyz = XYZs_it->second;
				/////////////   Do it!  //////////////////////////////
				find_downstream_neighbor( seqpos,xyz,neighbor,rinfos, glycan_positions, rosetta_residue_name3s, options );
				/////////////////////////////////////////////////////
				if ( neighbor.first != 0 || neighbor.second != "C0" ) {
					TR.Trace << "New branch found" << std::endl;
					++n_branches;
					if ( is_carbohydrate == false || neighbor.first != mainchain_neighbor ) {
						std::string const & child_resid = rinfos[ neighbor.first ].resid();
						if ( std::find( branch_lower_termini.begin(), branch_lower_termini.end(), child_resid ) == branch_lower_termini.end() ) {
							TR.Trace << "This residue is a branch point with corresponding branch lower terminus: " << neighbor.first << std::endl;
							is_branch_point = true;
							if ( branch_lower_termini.contains( child_resid ) == false ) {
								branch_lower_termini.push_back( child_resid );
							}
							branch_points_on_this_residue.push_back( atom_name );
						}
						if ( is_carbohydrate == false ) { break; }
					}
					neighbor = std::make_pair(0,"C0");
				}
			}
		}// Atom loop
	}// Asn/Ser/Thr or sugar?
}// function

/// @brief Helper function to find connected residues.
/// @detail Given an atom (seqpos, xyz-coordinates ...), find a neighboring residue.
/// This residue can be the mainchain coninuation, or a branched-off residue
/// Note: Only one neighbor will be returned. If more that one residue branches off from
/// this atom, only the first to be found will be returned.
/// @author Sebastian Rämisch, raemisch@scripps.edu
void
find_downstream_neighbor( core::Size const seqpos, Vector const & upstream_atom_xyz, std::pair<core::Size, std::string> & neighbor, utility::vector1< core::io::ResidueInformation > const & rinfos, utility::vector1< Size > const & glycan_positions, utility::vector1< std::string > const & rosetta_residue_name3s, StructFileRepOptions const & options )
{
	using namespace core::chemical;

	// NST: Asn/Ser/Thr
	bool const is_NST = ( rinfos[ seqpos ].resName() == "ASN" || rinfos[ seqpos ].resName() == "SER" || rinfos[ seqpos ].resName() == "THR" );
	core::Real max_cutoff = options.max_bond_length();
	core::Real min_cutoff = options.min_bond_length();
	Vector downstream_atom_xyz(0.0, 0.0, 0.0);
	for ( core::Size i : glycan_positions ) {
		if ( !is_NST && ( i > ( seqpos + 20 ) ) ) { return; } //For carbohydrates, the 20 next residue should be enough to check
		if ( !is_NST && i == seqpos ) { continue; } // Only check residues with higher numbers than seqpos
		TR.Trace << "check residue " << i << " ( " << rinfos[ i ].resid() << " )" << std::endl;
		// Loop through atoms of potential neighbor residue
		std::map< std::string, Vector > const & child_atom_coords = rinfos[ i ].xyz();
		for ( std::map< std::string, Vector >::const_iterator child_atom = child_atom_coords.begin(); child_atom != child_atom_coords.end(); ++child_atom ) {
			if ( ( child_atom->first == " C1 " ) || ( ( child_atom->first == " C2 " ) && ( rosetta_residue_name3s[ i ] == "Neu" )  )  ) {
				downstream_atom_xyz = child_atom->second; // xyz coordinates
				core::Real distance = downstream_atom_xyz.distance(upstream_atom_xyz);
				TR.Trace << " Distance to " << child_atom->first << ": " << distance << std::endl;
				if ( (distance <= max_cutoff) && (distance >= min_cutoff) /*&& (i > seqpos)*/ ) {
					TR.Trace << "One neighbor was found to be " << i << ", at atom " << child_atom->first << std::endl;
					neighbor = std::make_pair(i,child_atom->first);
					return;
				}
			}
		}
	}
	// Apparently no neighbor was found
}

/// @detail Loop through a residue's atoms to see if there are othe residues connected to it.
/// Then the residue position that contiues the main chain (if any) is returned.
Size
find_mainchain_connection( utility::vector1< core::io::ResidueInformation >& rinfos, core::io::StructFileRep& sfr,  std::string const & resid, Size const & seqpos, utility::vector1< std::string > const & rosetta_residue_name3s , bool & same_chain_next, bool & is_upper_terminus, int const CARB_MAINCHAIN_CONN_POS, utility::vector1< Size >& glycan_positions, StructFileRepOptions const & options ) {
	TR.Trace << "No LINK record found or -auto_detect_glycan_connections set to true. Checking for neighbors using coordinates." << std::endl;
	TR.Trace << "Residue Nr: " << resid  << ", at sequence position: " << seqpos <<  std::endl;
	std::pair< core::Size,std::string  > neighbor = std::make_pair(0,"C0"); // <residue number,atom name>
	std::map< std::string, std::pair< core::Size,std::string > > neighbor_map;
	std::map< Size,std::string > connection_map;
	std::map< std::string, Vector > const & XYZs = rinfos[ seqpos ].xyz();
	for ( std::map< std::string, Vector >::const_iterator XYZs_it = XYZs.begin(); XYZs_it != XYZs.end(); ++XYZs_it ) {
		// Find the first exocyclic oxigen
		std::string const & atom_name = XYZs_it->first;
		TR.Trace << "Check main chain connection at atom " << atom_name << std::endl;
		if ( atom_name.find(" O") != std::string::npos ) {  // is oxygen
			//TR.Trace << "Name: " << atom_name << " " << xyz.at(0) << " "<< xyz.at(1) << " " << xyz.at(2) << std::endl;
			/////////////// Do it! //////////////////////////////
			Vector const xyz = XYZs_it->second;
			find_downstream_neighbor( seqpos,xyz,neighbor,rinfos, glycan_positions, rosetta_residue_name3s, options);
			/////////////////////////////////////////////////////
			if ( neighbor.first != 0 && neighbor.second != "C0" ) {
				connection_map[ neighbor.first ] = atom_name; // <residue number,upstream atom name>
				//std::string const & neighbor_resid = rinfos[ neighbor.first ].resid(); JAB - comment out unused variable.
				// Main chain continues at residue neighbor.first
				TR.Trace << "Main chain continues at residue " << neighbor.first << "( " << rinfos[ neighbor.first ].resid() << " )"  << std::endl;
				TR.Trace << "Assigning main-chain connectivity to position " << connection_map[ neighbor.first ];
				TR.Trace << " of this residue." << std::endl;
				sfr.residue_type_base_names()[ resid ].second[ CARB_MAINCHAIN_CONN_POS ] = connection_map[ neighbor.first ][2];
				same_chain_next = true;
				is_upper_terminus = false;
				return neighbor.first;
			} else {
				continue;
			}
		}
	}// atom loop
	// If the function arrives here, there was no neighbor found
	TR.Trace << "No downstream mainchain connection. is_upper_terminus = true" << std::endl;
	same_chain_next = false;
	return 0;
}


/// @brief Recursively find a child residue and it's children and it's children ....
/// This function figures out in which order glycan residues are connected
/// @author Sebastian Rämisch, raemisch@scripps.edu
void
find_children( Size const seqpos, utility::vector1< core::io::ResidueInformation > const & rinfos, chemical::ResidueTypeSetCOP residue_type_set, utility::vector1< std::string > const & rosetta_residue_name3s, utility::vector1< core::Size > & correct_order, utility::vector1< Size > const & glycan_positions, utility::vector1< core::Size > & glycan_positions_temp, StructFileRepOptions const & options ) {
	using namespace core::chemical;
	TR.Trace << "Find children for residue " << seqpos << std::endl;
	correct_order.push_back( seqpos );
	TR.Trace << "Correct order: " <<  correct_order << std::endl;
	std::pair< core::Size,std::string  > neighbor = std::make_pair(0,"C0"); // <residue number,atom name>
	std::string const & name3 = rosetta_residue_name3s[ seqpos ];
	ResidueTypeCOP RT = ResidueTypeFinder( *residue_type_set ).name3( name3 ).get_representative_type();
	//TODO: remove the neighbor_list and the call to it
	utility::vector1< Size > neighbor_list;
	for ( Size atom = 1; atom <= RT->natoms(); ++atom ) {
		neighbor = std::make_pair(0,"C0");
		// Find the first exocyclic oxigen
		std::string const & atom_name = RT->atom_name(atom);
		TR.Trace << "Find child residue at atom " << atom_name << std::endl;
		std::map< std::string, Vector > XYZs(rinfos[ seqpos ].xyz());
		if ( atom_name.find(" O") != std::string::npos ) {  // is oxygen
			/////////////// Do it! //////////////////////////////
			if ( XYZs.count( atom_name ) ) {
				Vector const xyz = XYZs.at( atom_name );
				find_downstream_neighbor( seqpos,xyz,neighbor,rinfos, glycan_positions, rosetta_residue_name3s, options);
			} else {
				//TR.Trace << atom_name << " for residue " << seqpos << " not found " << std::endl;
				continue;
			}
			/////////////////////////////////////////////////////
			TR.Trace << "Neighbor pair for residue " << seqpos << ": " << neighbor << std::endl;
			if ( (neighbor.first != 0) && (neighbor.second != "C0") && (neighbor.first != seqpos) && ( std::find(correct_order.begin(), correct_order.end(), neighbor.first) == correct_order.end())  )  {
				Size neighbor_resnum = neighbor.first;
				neighbor_list.push_back( neighbor_resnum );
				// recursive function call
				find_children( neighbor_resnum, rinfos, residue_type_set, rosetta_residue_name3s, correct_order, glycan_positions , glycan_positions_temp, options );
			}
		}
	}// atom loop
	if ( neighbor_list.size() == 0 ) { TR.Trace << "No neighbors found for residue: " << seqpos << std::endl; }
	// remove this residue from the temp_list while reordering
	if ( glycan_positions_temp.size() > 0 ) {
		glycan_positions_temp.erase(std::remove(glycan_positions_temp.begin(), glycan_positions_temp.end(), seqpos), glycan_positions_temp.end());
	}

}//find_children

/// @brief Test if a glycan residue is the first of it's tree
/// This is done by checking if there is anything conncted at the C1 position
/// NOTE: if a parent redue has a sequence position that is greater by more than 20 ( e.g. Glc1021 -> Man1000 )
/// this fail. Can be changed easily, though. It's set like that to avoid checking too many residues
/// @author Sebastian Rämisch, raemisch@scripps.edu
bool is_root( core::Size const seqpos, utility::vector1< core::io::ResidueInformation >const & rinfos, utility::vector1< Size > const & glycan_positions, StructFileRep::Strings & branch_lower_termini_extra, StructFileRepOptions const & options )
{
	using namespace core::chemical;
	TR.Trace << "Find upstream residue for residue " << seqpos << "(" << rinfos[ seqpos ].resid() << ")" <<  std::endl;
	// Find the first exocyclic oxigen
	std::map< std::string, Vector > XYZs(rinfos[ seqpos ].xyz());
	Vector C1_atom_coords = XYZs[ " C1 " ];
	core::Real max_cutoff = options.max_bond_length();
	core::Real min_cutoff = options.min_bond_length();
	for ( core::Size i : glycan_positions ) {
		//std::map< std::string, Vector > xyzs(rinfos[ *it ].xyz());
		if ( i > ( seqpos + 20 ) ) { break; } //For carbohydrates, the 20 next residue should be enough to check
		if ( i == seqpos ) { continue; } // Exclude distance checking to self
		TR.Trace << "check if residue " << i << " is upstream of residue " << seqpos << std::endl;
		// Loop through atoms of potential neighbor residue
		std::map< std::string, Vector > const & parent_atom_coords = rinfos[ i ].xyz();
		for ( std::map< std::string, Vector >::const_iterator parent_atom = parent_atom_coords.begin(); parent_atom != parent_atom_coords.end(); ++parent_atom ) {
			// TODO: only check oxygens
			//if ( parent_atom->first != " O_something " ) { continue; }
			Vector parent_atom_xyz = parent_atom->second; // xyz coordinates
			core::Real distance = parent_atom_xyz.distance( C1_atom_coords );
			TR.Trace << " Distance to " << parent_atom->first << ": " << distance << std::endl;
			if ( (distance <= max_cutoff) && (distance >= min_cutoff) /*&& (i > seqpos)*/ ) {
				TR.Trace << "Residue " << seqpos << " is not a glycan tree root" << std::endl;
				TR.Trace << "Upstream neighbor was found to be " << i << ", at atom " << parent_atom->first << std::endl;
				return 0;
			}
		}
	}//glycan_positions loop
	// Is it attached to an Asn?

	//JAB - This should really go in it's own function due to changing Branch lower termini.
	size_t i = 0;
	for ( utility::vector1< core::io::ResidueInformation >::const_iterator it = rinfos.begin() ; it != rinfos.end(); ++it ) {
		++i;
		if ( it->resName() == "ASN" ) {
			std::map< std::string, Vector > const & parent_atom_coords = it->xyz();
			for ( std::map< std::string, Vector >::const_iterator parent_atom = parent_atom_coords.begin(); parent_atom != parent_atom_coords.end(); ++parent_atom ) {
				Vector const & parent_atom_xyz  = parent_atom->second; // xyz coordinates
				core::Real distance = parent_atom_xyz.distance( C1_atom_coords );
				if ( ( distance <= max_cutoff ) && ( distance >= min_cutoff ) ) {
					TR.Trace << "Sugar residue " << seqpos << " is bound to Asn " << it->resid() << std::endl;
					// If glycan is listed before the protein-branch point, it will create problems
					if ( seqpos < i ) {
						std::stringstream err_msg;
						err_msg << "Glycan " << rinfos[ seqpos ].resid() << " comes before Asn " << it->resid()
							<< ". Please fix the input file, so that the connected glycans appear after this Asn." <<    std::endl;
						utility_exit_with_message( err_msg.str() );
					}
					branch_lower_termini_extra.push_back( rinfos[ seqpos ].resid() );
					return 1;
				}
			}
		}
	}
	return 1;
}

/// @brief Function to identify glycans and fix their rosetta names
/// This function determines the correct name3s and fills in the glycan_positions vector
/// @author Sebastian Rämisch, raemisch@scripps.edu
void fix_residue_info_and_order(utility::vector1< core::io::ResidueInformation >& rinfos, core::io::StructFileRep& sfr, chemical::ResidueTypeSetCOP residue_type_set, utility::vector1< std::string >& rosetta_residue_name3s, StructFileRep::Strings & branch_lower_termini_extra, utility::vector1< std::string >& glycan_tree_roots, utility::vector1< core::Size >& glycan_positions, StructFileRepOptions const & options )
{
	TR.Trace << "Detecting glycans and residue order! " << std::endl;
	using namespace core::chemical;
	for ( Size ii = 1; ii <= rinfos.size(); ++ii ) {
		// Convert PDB 3-letter code to Rosetta 3-letter code, if a list of alternative codes has been provided.
		std::pair< std::string, std::string > const & rosetta_names(
			NomenclatureManager::get_instance()->rosetta_names_from_pdb_code( rinfos[ ii ].resName() ) );
		std::string const & name3( rosetta_names.first );
		rosetta_residue_name3s[ ii ] = name3;
		std::string const & resid = rinfos[ ii ].resid();
		if ( carbohydrates::CarbohydrateInfoManager::is_valid_sugar_code( name3 ) ) {
			TR.Trace << "Identified glycan at position " << ii << std::endl;
			glycan_positions.push_back( ii );
			sfr.residue_type_base_names()[ resid ] = std::make_pair( name3, rosetta_names.second );
		}
	}
	// Sugar residues have to be in a certain order. E.g.
	// 1--2--3--4--5--6
	//          |--7--8
	TR.Trace << "Glycan vector:" << std::endl;
	TR.Trace << glycan_positions << std::endl;

	// the glycan_positions_temp will shrink with every identified sugar until it's size = 0
	// it's needed to be able to find the next yet unchecked glycan tree root. Once the next
	// root is picked, all it's children are examined by the find_children function and those
	// will be removed from the glycan_positions_temp vector.
	utility::vector1< core::Size > correct_order;
	utility::vector1< core::Size > glycan_positions_temp( glycan_positions );
	size_t indx = 1;
	while ( glycan_positions_temp.size() > 0 ) {
		bool root = is_root( glycan_positions_temp[ indx ], rinfos, glycan_positions, branch_lower_termini_extra, options  );
		if ( root ) {
			TR.Trace << "Residue " << glycan_positions_temp[ indx ] << " (" << rinfos[ glycan_positions_temp[ indx ] ].resid()  << ") is a root." << std::endl;
			glycan_tree_roots.push_back( rinfos[ glycan_positions_temp[ indx ] ].resid() );
			//////////////////////////////////////////////
			find_children( glycan_positions_temp[ indx ], rinfos, residue_type_set, rosetta_residue_name3s, correct_order, glycan_positions, glycan_positions_temp, options );
			/////////////////////////////////////////////
			indx = 1;
		} else {
			++indx;
		}
	}
	// At this point, the correct sugar order is determined and the ResidueInfo objects in rinfos
	// need to be re-ordered to match the connectivity
	reorder_glycan_residues( rinfos, rosetta_residue_name3s, correct_order, glycan_positions );
}

/// @brief Bring glycans into the correct order, which corresponds to connectivity of ech glycan tree
/// This requires reordering rinfos and rosetta_residue_name3s.
void reorder_glycan_residues( utility::vector1< core::io::ResidueInformation >& rinfos, utility::vector1< std::string >& rosetta_residue_name3s, utility::vector1< core::Size >& correct_order, utility::vector1< core::Size > const & glycan_positions  )
{
	utility::vector1< core::io::ResidueInformation > rinfos_orig = rinfos;
	utility::vector1< std::string > rosetta_residue_name3sorig = rosetta_residue_name3s;
	using namespace core::chemical;
	TR << "Automatic glycan connection is activated." << std::endl;
	TR << "Start reordering residues." << std::endl;
	TR << "Corrected glycan residue order (internal numbering): " << correct_order << std::endl;
	TR.Trace << std::endl;
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
		std::string new_name3 = rosetta_residue_name3sorig[ corrected_pos ];
		rosetta_residue_name3s[ i ] = new_name3;
	}


}



/// @details Helper function for sorting LinkInformation records
/// Bubble sorting the connections obtained from LINK records to
/// make sure that the mainchain continues with the next residue
/// according to PDB-numbering. Building a pose would fail otherwise.
/// This automatically fixes LINK records that are in the wrong order.
/// - Sebastian Rämisch
void
sort_link_records( utility::vector1<LinkInformation> &link_records ) {
	// sorting
	Size indx = 1;
	bool was_swapped = true;
	while ( was_swapped ) {
		was_swapped = false;
		for ( utility::vector1<LinkInformation>::iterator it = link_records.begin(); indx < link_records.size(); ++it ) {
			LinkInformation holder = *it;
			if ( it->resSeq2 > (it+1)->resSeq2 ) {
				// swap their positions
				link_records[indx] = *(it+1);
				link_records[indx+1] = holder;
				was_swapped = true;
			}
			++indx;
		}
	}//sorting
}//sort_link_records


} //core
} //io
