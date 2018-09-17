// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/select/util.cc
/// @brief Utilities to help in selecting residues.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)
/// @author Mike Tyka
/// @author Chu Wang
/// @author Daniel J. Mandell


#include <core/select/util.hh>

#include <basic/Tracer.hh>

#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <core/select/residue_selector/SymmetricalResidueSelector.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/MirrorSymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/scoring/Energies.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/symmetry/util.hh>
#include <core/id/AtomID.hh>
#include <core/chemical/AtomType.hh>



#include <numeric/xyzVector.hh>

#include <utility/string_util.hh>

static basic::Tracer TR( "core.select.util" );


namespace core {
namespace select {

utility::vector1< Size >
get_residues_from_subset( utility::vector1< bool > const & subset, bool select){
	utility::vector1< Size > residues;
	for ( core::Size i = 1; i <= subset.size(); ++i ) {
		if ( subset[i] == select ) {
			residues.push_back( i );
		}
	}
	return residues;
}

std::set< core::Size >
get_residue_set_from_subset( utility::vector1< bool > const & subset, bool select){
	std::set< Size > residues;
	for ( core::Size i = 1; i <= subset.size(); ++i ) {
		if ( subset[i] == select ) {
			residues.insert( i );
		}
	}
	return residues;

}

utility::vector1< bool >
get_subset_from_residues( utility::vector1< core::Size > const & residues, core::Size total_nres, bool invert) {
	// If invert is false, we default to false and re-set residues in the passed vector to not-false (true).
	// If invert is true, we default to true and re-set residues in the passed vectorr to not-true (false).
	utility::vector1< bool > subset( total_nres, invert );

	for ( core::Size resnum: residues ) {
		debug_assert( resnum != 0 && resnum <= total_nres );
		subset[ resnum ] = !invert; // We set it to the opposite of the default.
	}

	return subset;
}

core::select::residue_selector::ResidueIndexSelectorOP
get_residue_selector_from_subset(
	core::select::residue_selector::ResidueSubset subset) {

	utility::vector1<Size> seq_pos = get_residues_from_subset( subset );
	std::string res_string = utility::join( seq_pos, "," );

	return residue_selector::ResidueIndexSelectorOP( new residue_selector::ResidueIndexSelector( res_string ) );
}

utility::vector1< bool >
get_neighbor_residues(
	core::pose::Pose const & pose,
	utility::vector1< bool > const & residue_positions,
	core::Real neighbor_dis
)
{
	utility::vector1< bool > selection_and_neighbors = residue_positions;

	if ( neighbor_dis <= 10.0 ) {
		fill_neighbor_residues(pose, selection_and_neighbors, neighbor_dis);

		//Make sure to turn off subset residues!
		for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
			if ( residue_positions[i] ) {
				selection_and_neighbors[i] = false;
			}
		}

		return selection_and_neighbors;

	} else {
		utility_exit_with_message("get_neighbor_residues only currently works for neighbor distances of 10A or less!  Please use NeighborhoodResidueSelector instead");
	}


}


void
fill_neighbor_residues(
	core::pose::Pose const & pose,
	utility::vector1< bool > & residue_positions,
	core::Real neighbor_dis
)
{
	utility::vector1< bool > selection = residue_positions;

	fill_tenA_neighbor_residues(pose, residue_positions);
	filter_neighbors_by_distance(pose, selection, residue_positions, neighbor_dis);
}


utility::vector1< bool >
trim_neighbors_by_distance(
	core::pose::Pose const & pose,
	utility::vector1<bool> const & selection,
	utility::vector1<bool> const & all_neighbors,
	core::Real & dist_cutoff
)
{
	utility::vector1< bool > trimmed_neighbors = all_neighbors;
	filter_neighbors_by_distance(pose, selection, trimmed_neighbors, dist_cutoff);
	return trimmed_neighbors;

}


void
filter_neighbors_by_distance(
	core::pose::Pose const & pose,
	utility::vector1<bool> const & selection,
	utility::vector1<bool> & selection_and_neighbors,
	core::Real & dist_cutoff
)
{

	for ( Size i = 1; i <= selection_and_neighbors.size(); ++i ) {
		if ( selection_and_neighbors[ i ] == false ) continue;

		selection_and_neighbors[ i ] = false; //Get ready to change this.
		for ( Size x = 1; x <= selection.size(); ++x ) {
			if ( ! selection[x] ) continue;

			// Get the atom vectors for loop and scaffold CB, or CA if GLY
			numeric::xyzVector< Real > neighbor_vec;
			numeric::xyzVector< Real > select_vec;
			neighbor_vec = pose.residue( i ).xyz( pose.residue( i ).nbr_atom() );
			select_vec = pose.residue( x ).xyz( pose.residue( x ).nbr_atom() );
			// only keep as neighbor if dist within cutoff
			Real dist = neighbor_vec.distance( select_vec );
			if ( dist <= dist_cutoff ) {
				selection_and_neighbors[ i ] = true;
			}
		}
	}


}

utility::vector1< bool >
get_tenA_neighbor_residues(
	pose::Pose const & pose,
	utility::vector1<bool> const & residue_positions
)
{
	utility::vector1< bool > positions_and_neighbors = residue_positions;
	fill_tenA_neighbor_residues(pose, positions_and_neighbors);
	return positions_and_neighbors;
}


void fill_tenA_neighbor_residues(
	pose::Pose const & pose,
	utility::vector1<bool> & residue_positions
)
{
	//make a local copy first because we will change content in residue_positions
	utility::vector1<bool> local_residue_positions = residue_positions;
	core::scoring::TenANeighborGraph const & tenA_neighbor_graph( pose.energies().tenA_neighbor_graph() );
	for ( Size i=1; i <= local_residue_positions.size(); ++i ) {
		if ( ! local_residue_positions[i] ) continue;
		utility::graph::Node const * current_node( tenA_neighbor_graph.get_node(i)); // find neighbors for this node
		for ( utility::graph::Node::EdgeListConstIter it = current_node->const_edge_list_begin();
				it != current_node->const_edge_list_end(); ++it ) {
			Size pos = (*it)->get_other_ind(i);
			residue_positions[ pos ] = true;
		}
	}
}

//Author: Jared Adolf-Bryfogle (jadolfbr@gmail.com)
std::string
get_pymol_selection_for_atoms(pose::Pose const & pose, utility::vector1< id::AtomID > const & atoms, std::string const & sele_name, bool skip_virts /*true*/ ){

	std::string selection = "select "+sele_name+", ";
	std::map< std::string, utility::vector1<std::string >> chain_residues;
	std::map< std::string, utility::vector1<std::string >> residue_to_atoms;
	utility::vector1< std::string > chains;
	for ( auto atom : atoms ) {

		core::Size i = atom.rsd();

		if ( pose.residue(i).atom_type( atom.atomno() ).is_virtual() && skip_virts ) continue;

		std::string chain = utility::to_string( pose.pdb_info()->chain(i) );
		std::string num = utility::to_string( pose.pdb_info()->number(i) );
		std::string icode = utility::to_string( pose.pdb_info()->icode(i) );

		std::string res = "";
		if ( icode == utility::to_string(' ') ) {
			res = num;
		} else {
			res = num+icode;
		}

		std::string atom_name = pose.residue( i ).atom_name( atom.atomno());
		if ( ! residue_to_atoms.count(res) ) {
			utility::vector1< std::string > atom_names;
			residue_to_atoms[res] = atom_names;
			residue_to_atoms[res].push_back(atom_name);

		} else {
			residue_to_atoms[res].push_back(atom_name);
		}

		if ( ! chain_residues.count( chain ) ) {
			utility::vector1< std::string > residues;
			chain_residues[chain] = residues;
			chain_residues[chain].push_back( res );
			chains.push_back( chain );
		} else if ( !chain_residues[chain].contains(res) ) {
			chain_residues[chain].push_back( res );
		}
	}

	for ( core::Size i = 1; i <= chains.size(); ++i ) {

		std::string chain = chains[i];
		std::string subselection = "";

		if ( i == 1 ) {
			subselection = "(chain "+chain+" and resid ";
		} else {
			subselection = subselection + " or "+"(chain "+chain+" and resid ";
		}
		for ( core::Size x = 1; x <= chain_residues[ chain ].size(); ++x ) {
			std::string res = chain_residues[chain][x];

			if ( x == 1 ) {
				subselection += res;
			} else {
				subselection += ("resid "+res);
			}

			//Atom Subselection
			std::string atom_subselection = "";
			//std::cout << x <<" " << residue_to_atoms[res] << std::endl;
			for ( core::Size A = 1; A <= residue_to_atoms[res].size(); ++A ) {
				std::string const & atom_name = residue_to_atoms[res][A];
				if ( A == 1 ) {
					atom_subselection = " and name ";
				}
				atom_subselection += utility::strip(atom_name);

				if ( A != residue_to_atoms[res].size() ) {
					atom_subselection += "+";
				}
			}
			subselection += atom_subselection;

			if ( x != chain_residues[chain].size() ) {
				subselection += " or ";
			}

		}
		subselection+=")";
		//TR << chain <<" " <<subselection << std::endl;
		selection += subselection;
	}

	return selection;
}

///@brief This is for a 'symmetrical' selection (iE Normal selection for symmetrical poses!).
///  Turns off all residues that are not part of the master subunit.
utility::vector1< bool >
get_master_subunit_selection(pose::Pose const & pose, utility::vector1<bool> const & full_subset){
	using namespace core::conformation::symmetry;

	utility::vector1< bool > result = full_subset;
	if ( ! core::pose::symmetry::is_symmetric( pose) ) {
		return result;
	}

	auto const & symm_conf (
		dynamic_cast<SymmetricConformation const & > ( pose.conformation()) );
	SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );


	for ( core::Size i = 1; i <= full_subset.size(); ++i ) {
		if ( ! full_subset[i] ) continue;

		if ( symm_info->bb_is_independent( i ) ) {
			continue;
		} else {
			result[i] = false;
		}
	}
	return result;
}



} //core
} //select


