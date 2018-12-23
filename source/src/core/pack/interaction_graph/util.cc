// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/interaction_graph/util.cc
/// @brief  Simple utilities for interaction graphs.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

// Package Headers
#include <core/pack/interaction_graph/util.hh>

// Core Headers
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/pack/interaction_graph/AnnealableGraphBase.hh>
#include <core/pack/interaction_graph/InteractionGraphBase.hh>
#include <core/pack/interaction_graph/DensePDInteractionGraph.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileOptions.hh>

// Utility Headers
#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace interaction_graph {

/// @brief Given an interaction graph, get a summary describing it fully.
/// @details This is intended to be a machine-readable summary that facilitates a Rosetta interface to
/// external annealers or optimizers.
/// @param[out] igstream The output stream for the summary.
/// @param[in] pose The pose, for reference.
/// @param[in] rot_sets The RotamerSets object containing the set of sets of rotamers for each position.
/// @param[in] anngraph The pre-computed interaction graph.
/// @param[in] short_version  If true, only the interaction graph summary is returned.  If false, information needed to reconstruct the pose
/// is also included.  False by default.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
void
get_annealable_graph_summary(
	std::stringstream & igstream,
	core::pose::Pose const & pose,
	core::pack::rotamer_set::RotamerSets const & rot_sets,
	AnnealableGraphBaseCOP const & anngraph,
	bool const short_version
) {
	std::string const errmsg( "Error in core::pack::interaction_graph::get_annealable_graph_summary(): ");

	// Some checks up front:
	runtime_assert_string_msg( pose.total_residue() == rot_sets.total_residue(), errmsg + "The number of residues in the pose and the rotamer sets object do not match." );
	runtime_assert_string_msg( pose.total_residue() >= rot_sets.nmoltenres(), errmsg + "Thge number of molten residues is somehow greater than the number of residues in the pose.  This should not be possible." );
	InteractionGraphBaseCOP pdgraph( utility::pointer::dynamic_pointer_cast< InteractionGraphBase const >( anngraph ) );
	runtime_assert_string_msg( pdgraph != nullptr, "The annealable graph passed to this function is not an InteractionGraph.  Note that annealable graph summaries cannot be produced if the scoring function uses non-pairwise decomposible score terms." );

	if ( !short_version ) {
		// First, print out the pose to allow reconstruction:
		igstream << "[BEGIN POSE BINARY]" << "\n";
		core::io::silent::SilentFileOptions opts; // initialized from the command line
		core::io::silent::BinarySilentStruct silent_struct(opts, pose, "" );
		silent_struct.print_conformation(igstream);
		igstream << "[END POSE BINARY]" << "\n";

		// Next, print out the identities and chi values of each rotamer.
		igstream << "[BEGIN ROTAMER NAME/SEQPOS/ROTINDEX/CHIVALS]" << "\n";
		for ( core::Size i(1), imax( rot_sets.nmoltenres() ); i<=imax; ++i ) { //Iterate over the "molten" (packable) residues.
			core::pack::rotamer_set::RotamerSetCOP cur_rot_set( rot_sets.rotamer_set_for_moltenresidue( i ) );
			core::Size const seqpos( rot_sets.moltenres_2_resid( i ) );
			for ( core::Size j(1), jmax( cur_rot_set->num_rotamers() ); j<=jmax; ++j ) {
				core::conformation::ResidueCOP curres( cur_rot_set->rotamer(j) );
				igstream << "\t[BEGIN_ROT]" << "\n";
				igstream << "\t" << curres->name() << "\t" << seqpos << "\t" << j;
				for ( core::Size k(1), kmax( curres->nchi() ); k<=kmax; ++k ) {
					igstream << "\t" << curres->chi( k );
				}
				igstream << "\n";
				core::pose::Pose temppose;
				temppose.append_residue_by_jump(*curres, 0 );
				core::io::silent::BinarySilentStruct res_ss( opts, temppose, "" );
				res_ss.print_conformation(igstream);
				igstream << "\t[END_ROT]\n";
			}
		}
		igstream << "[END ROTAMER NAME/SEQPOS/INDEX/CHIVALS]" << "\n";
	}

	// Next, print out the onebody energies of each rotamer.
	igstream << "[BEGIN ONEBODY SEQPOS/ROTINDEX/ENERGY]" << "\n";
	for ( core::Size i(1), imax(pdgraph->get_num_nodes()); i<=imax; ++i ) {
		PrecomputedPairEnergiesNode const * curnode( dynamic_cast< PrecomputedPairEnergiesNode const * >( pdgraph->get_node(i) ) );
		core::Size const seqpos( rot_sets.moltenres_2_resid(i) );
		runtime_assert_string_msg( curnode != nullptr, errmsg + "Node " + std::to_string( i ) + " is not a PrecomputedPairEnergiesNode!" );
		for ( core::Size j(1), jmax(curnode->get_num_states()); j<=jmax; ++j ) {
			igstream << seqpos << "\t" << j << "\t" << curnode->get_one_body_energy(j) << "\n";
		}
	}
	igstream << "[END ONEBODY SEQPOS/ROTINDEX/ENERGY]" << "\n";

	// Next, print out the twobody energies of each rotamer pair.
	igstream << "[BEGIN TWOBODY SEQPOS1/ROTINDEX1/SEQPOS2/ROTINDEX2/ENERGY]" << "\n";
	for ( core::Size i(1), imax(pdgraph->get_num_nodes()); i<imax; ++i ) {
		core::Size const seqpos1( rot_sets.moltenres_2_resid(i) );
		for ( core::Size j(i+1); j<=imax; ++j ) {
			core::Size const seqpos2( rot_sets.moltenres_2_resid(j) );
			if ( !pdgraph->get_edge_exists(i,j) ) continue; //Skip if no edge between these nodes.
			PrecomputedPairEnergiesEdge const * curedge( dynamic_cast< PrecomputedPairEnergiesEdge const * >( pdgraph->find_edge(i, j) ) );
			runtime_assert_string_msg( curedge != nullptr, errmsg + "Edge between nodes " + std::to_string(i) + " and " + std::to_string(j) + " is not a PrecomputedPairEnergiesEdge!" );
			//Iterate over all states in nodes i and j:
			for ( core::Size ii(1), iimax(pdgraph->get_num_states_for_node(i)); ii<=iimax; ++ii ) {
				for ( core::Size jj(1), jjmax( pdgraph->get_num_states_for_node(j)); jj<=jjmax; ++jj ) {
					igstream << seqpos1 << "\t" << ii << "\t" << seqpos2 << "\t" << jj << "\t" << curedge->get_two_body_energy( ii, jj ) << "\n";
				}
			}
		}
	}
	igstream << "[END TWOBODY SEQPOS1/ROTINDEX1/SEQPOS2/ROTINDEX2/ENERGY]" << "\n";
}

/// @brief Given a description of a packer problem, a description of the solution, and a pose to put a result in, generate
/// a pose representing the solution.
/// @details The contents of the pose are destroyed by this operation, and replaced by the pose described in the packer problem
/// and the rotamers described in the solutions!
/// @param[out] pose The pose, cleared and rebuilt by this operation.
/// @param[in] packer_problem_description A description of the packer problem, including blocks for the original pose, the rotamers,
/// and the energies, as generated by get_annealable_graph_summary().
/// @param[out] packer_solution_description A solution to the packer problem, in the format of a series of lines, one for each packable position,
/// consisting of <seqpos #> <rotamer #>.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
void
set_externally_generated_packer_solution(
	core::pose::Pose & pose,
	std::string const & packer_problem_description,
	std::string const & packer_solution_description
) {
	//First, generate the pose:
	generate_pose_from_packer_description( pose, packer_problem_description );

	//Next, add the appropriate rotamers:
	utility::vector1< std::pair< core::Size, core::conformation::ResidueCOP > > rotamers;
	generate_rotamers_from_packer_solution_description( rotamers, packer_problem_description, packer_solution_description );
	add_rotamers_to_pose( pose, rotamers );
}

/// @brief Given the description of a packer problem as generated by get_annealable_graph_summary(), pull out
/// and re-generate the pose.
/// @details The contents of the pose are destroyed and replaced by the pose from the description in this operation.
void
generate_pose_from_packer_description(
	core::pose::Pose & pose,
	std::string const & packer_problem_description
) {
	static std::string const errmsg( "Error in core::pack::interaction_graph::generate_pose_from_packer_description(): " );

	pose.clear();
	bool in_block(false);
	std::stringstream problem_ss( packer_problem_description );
	std::string linestring;
	utility::vector1< std::string > lines;
	while ( getline(problem_ss, linestring ) ) {
		if ( !in_block ) {
			utility::strip_whitespace( linestring );
			if ( linestring == "[BEGIN POSE BINARY]" ) {
				in_block = true;
				continue;
			}
		} else {
			if ( linestring == "[END POSE BINARY]" ) {
				break;
			}
			lines.push_back( linestring );
		}
	}

	runtime_assert_string_msg( in_block && !lines.empty(), errmsg + "File did not contain a \"[BEGIN POSE BINARY]...[END POSE BINARY]\" block." );
	core::io::silent::SilentFileOptions opts;
	core::io::silent::BinarySilentStruct silentstruct(opts);
	core::io::silent::SilentFileData dummy_container(opts);
	silentstruct.init_from_lines(lines, dummy_container, true);
	silentstruct.fill_pose(pose);
}

/// @brief Given a packer problem and its solution, generate rotamers from it.
/// @details The rotamer vector (a vector of pairs of <seqpos, residueCOP>) is cleared and populated by this operation.
void
generate_rotamers_from_packer_solution_description(
	utility::vector1< std::pair< core::Size, core::conformation::ResidueCOP > > & rotamers,
	std::string const & packer_problem_description,
	std::string const & packer_solution_description
) {
	static std::string const errmsg( "Error in core::pack::interaction_graph::generate_rotamers_from_packer_solution_description(): " );

	core::io::silent::SilentFileOptions opts;

	std::stringstream packer_solution_stream( packer_solution_description );
	std::string curline;
	utility::vector1< core::Size > seqpos_seen;
	std::map< std::pair< core::Size, core::Size >, bool > rotamers_needed; // Key pair is <seqpos, rotindex>; bool indicates whether it has been found.

	// Figure out which rotamers I need:
	while ( getline(packer_solution_stream, curline) ) {
		std::stringstream curline_stream( curline );
		core::Size seqpos, rotindex;
		curline_stream >> seqpos;
		runtime_assert_string_msg( !( curline_stream.fail() || curline_stream.bad() ), errmsg + "Could not parse sequence position from line \"" + curline + "\"." );
		runtime_assert_string_msg( !seqpos_seen.has_value(seqpos), "Position " + std::to_string(seqpos) + " occurs multiple times in the solution description." );
		seqpos_seen.push_back(seqpos);
		curline_stream >> rotindex;
		runtime_assert_string_msg( !( curline_stream.fail() || curline_stream.bad() ), errmsg + "Could not parse rotamer index from line \"" + curline + "\"." );
		debug_assert( rotamers_needed.count(std::make_pair(seqpos, rotindex)) == 0 ); //There should not already exist this key in the map.
		rotamers_needed[ std::make_pair( seqpos, rotindex ) ] = false;
	}

	//Find those rotamers:
	rotamers.clear();
	bool block_found(false);
	std::stringstream problem_description_stream( packer_problem_description );
	while ( getline( problem_description_stream, curline) ) {
		utility::strip_whitespace(curline);
		if ( curline == "[BEGIN ROTAMER NAME/SEQPOS/ROTINDEX/CHIVALS]" ) {
			block_found = true;
			break;
		}
	}
	runtime_assert_string_msg( block_found, errmsg + "The packer problem description file does not have a \"[BEGIN ROTAMER NAME/SEQPOS/ROTINDEX/CHIVALS]...[END ROTAMER NAME/SEQPOS/ROTINDEX/CHIVALS]\" block." );
	//Now I'm in the rotamers block:
	{
		bool in_rotamer_block(false), parsed_first_line(false), parsing_this_rotamer(false);
		utility::vector1< std::string > currot_lines;
		core::Size seqpos(0), rotno(0);
		while ( getline( problem_description_stream, curline ) ) {
			if ( !in_rotamer_block ) {
				utility::strip_whitespace(curline);
				if ( curline == "[BEGIN_ROT]" ) {
					in_rotamer_block = true;
					parsed_first_line = false;
					parsing_this_rotamer = false;
					continue;
				} else if ( curline == "[END ROTAMER NAME/SEQPOS/ROTINDEX/CHIVALS]" ) {
					break;
				}
			} else {
				std::string curline2 = curline;
				utility::strip_whitespace(curline2);
				if ( curline2 == "[END_ROT]" ) {
					if ( parsing_this_rotamer ) {
						core::io::silent::BinarySilentStruct silentstruct(opts);
						core::io::silent::SilentFileData dummy_container(opts);
						silentstruct.init_from_lines(currot_lines, dummy_container, true);
						core::pose::Pose temppose;
						silentstruct.fill_pose(temppose);
						runtime_assert_string_msg(temppose.total_residue() == 1, errmsg + "Rotamer pose was improperly stored.  Could not regenerate rotamer for sequence position " + std::to_string(seqpos) );
						rotamers.push_back( std::pair< core::Size, core::conformation::ResidueCOP >(seqpos, temppose.conformation().residue_cop(1)) );
						currot_lines.clear();
					}
					in_rotamer_block = false;
					parsed_first_line = false;
					parsing_this_rotamer = false;
					continue;
				}
				if ( !parsed_first_line ) {
					utility::strip_whitespace(curline);
					parsed_first_line = true;
					std::stringstream curline_ss( curline );
					std::string rotname;
					curline_ss >> rotname >> seqpos >> rotno;
					runtime_assert_string_msg( !( curline_ss.fail() || curline_ss.bad() ), errmsg + "Could not parse line \"" + curline + "\"." );
					std::pair< core::Size, core::Size > currot_key( std::make_pair(seqpos, rotno) );
					parsing_this_rotamer = rotamers_needed.count(currot_key);
					if ( parsing_this_rotamer ) {
						runtime_assert_string_msg( rotamers_needed.at(currot_key) == false, errmsg + "Multiple entries found in the packer problem definition for sequence position " + std::to_string(seqpos) + ", rotamer index " + std::to_string(rotno) + "." );
						currot_lines.clear();
						rotamers_needed[currot_key] = true; //Found this rotamer.
					}
					continue; //Not strictly necessary, but safer.
				} else { //We have parsed the first line
					if ( parsing_this_rotamer ) {
						currot_lines.push_back(curline);
					}
				}
			}
		}
	}

	//Check that I found all rotamers:
	for ( auto const & rotamer_needed : rotamers_needed ) { //rotamer_needed is a std::pair< core::Size, core::Size > const &
		runtime_assert_string_msg(rotamer_needed.second, errmsg + "Rotamer " + std::to_string(rotamer_needed.first.second) + " at sequence position " + std::to_string(rotamer_needed.first.first) + " does not exist in the packer problem description, but was specified for the solution." );
	}
}

/// @brief Given a vector of solution rotamers and a pose, put the rotamers on the pose.
/// @details The vector is a vector of pairs of the form <seqpos, rotamerCOP>.
void
add_rotamers_to_pose(
	core::pose::Pose & pose,
	utility::vector1< std::pair< core::Size, core::conformation::ResidueCOP > > & rotamers
) {
	for ( core::Size i(1), imax(rotamers.size()); i<=imax; ++i ) {
		runtime_assert_string_msg( rotamers[i].first > 0 && rotamers[i].first <= pose.total_residue(), "Error in core::pack::interaction_graph::add_rotamers_to_pose(): Sequence position " + std::to_string(rotamers[i].first) + " does not exist in the pose!  (The pose is " + std::to_string(pose.total_residue()) + " residues long." );
		pose.replace_residue(rotamers[i].first, *(rotamers[i].second), false);
	}
}

} // namespace interaction_graph
} // namespace pack
} // namespace core
