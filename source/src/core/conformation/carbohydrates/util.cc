// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/conformation/carbohydrates/util.cc
/// @brief Utility functions that DO NOT require a pose.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <core/conformation/carbohydrates/util.hh>

// Package Headers
#include <core/chemical/carbohydrates/CarbohydrateInfo.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfoManager.hh>
#include <core/chemical/carbohydrates/database_io.hh>

// Project Headers
#include <core/id/types.hh>
#include <core/id/AtomID.hh>
#include <core/id/TorsionID.hh>
#include <core/types.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/AtomType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/util.hh>

// Utility Headers
#include <utility/string_constants.hh>
#include <utility/vector1.hh>
#include <utility/vector1.functions.hh>
#include <utility/exit.hh>
#include <utility/string_util.hh>



// External Headers
#include <boost/lexical_cast.hpp>
#include <basic/basic.hh>

// C++ Header
#include <list>
#include <utility>

#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "core.conformation.carbohydrates.util" );


namespace core {
namespace conformation {
namespace carbohydrates {

// Helper Functions ///////////////////////////////////////////////////////////
// Use a saccharide residue's connections to find the residue from which it follows or branches.
/// @return  The sequence position of the residue before this one (n-1) or the residue in the parent chain from which
/// the branch occurs or zero if N/A, i.e., if this is the lower terminus.
core::uint
find_seqpos_of_saccharides_parent_residue( conformation::Residue const & residue ) {
	debug_assert( residue.is_carbohydrate() );


	if ( ! residue.is_lower_terminus() ) {
		uint const id_of_connection_to_parent(
			residue.type().residue_connection_id_for_atom( residue.carbohydrate_info()->anomeric_carbon_index() ) );
		return residue.residue_connection_partner( id_of_connection_to_parent );
	} else  {
		TR.Debug << "This residue is a lower terminus! Returning 0." << std::endl;
		return 0;
	}


	//JAB - this fails during pose loading, even though the residue types should be finalized.
	// Not sure exactly why this would fail.  So, for now, we use the original code.
	core::Size anomeric_carbon = residue.carbohydrate_info()->anomeric_carbon_index();
	uint const id_of_connection_to_parent(
		residue.type().residue_connection_id_for_atom( anomeric_carbon ) );


	return residue.residue_connection_partner( id_of_connection_to_parent );

}


core::uint
find_seqpos_of_saccharides_mainchain_child( conformation::Residue const & residue ){
	core::uint linkage_position = residue.carbohydrate_info()->mainchain_glycosidic_bond_acceptor();
	return find_seqpos_of_saccharides_child_residue_at( residue, linkage_position);
}


// TODO: What if this is is a sialic acid as reducing end?
// Use a saccharide residue's connections to find the residue following it from a given linkage position.
/// @param   <linkage_position>: an integer n of (1->n) of polysaccharide nomenclature, where n specifies the attachment
/// point on the parent monosaccharide residue; e.g., 4 specifies O4
/// @return  The sequence position of the residue downstream of this one attached to the given linkage positions.  This
/// is n+1, if the linkage position is the same as the main chain connectivity, or zero if N/A, i.e., if this is the
/// upper terminus.
core::uint
find_seqpos_of_saccharides_child_residue_at( conformation::Residue const & residue, core::uint linkage_position )
{
	using namespace std;
	using namespace chemical::carbohydrates;

	debug_assert( residue.is_carbohydrate() );
	debug_assert( linkage_position <= CarbohydrateInfo::MAX_C_SIZE_LIMIT );

	CarbohydrateInfoCOP info( residue.carbohydrate_info() );

	if ( residue.is_upper_terminus() ) {
		TR.Debug << "Residue " << residue.seqpos() << " is an upper terminus! Returning 0." << endl;
		return 0;
	} else if ( linkage_position > info->n_carbons() ) {
		TR.Debug << "Residue " << residue.seqpos() << " does not have a position " << linkage_position << '!' << endl;
		return 0;
	} else {
		string atom_name( "O" + string( 1, ( '0' + linkage_position ) ) );  // to get "O1", "O2", etc.

		if ( ! residue.has( atom_name ) ) {
			TR.Warning << "Residue " << residue.seqpos() << " does not have an " << atom_name << endl;
			return 0;
		}

		TR.Debug << "Finding seqpos of child residue of " << residue.seqpos() << " at " << atom_name;

		uint index( residue.atom_index( atom_name ) );
		if ( info->cyclic_oxygen_index() == index ) {  // This means that this is an exocyclic linkage.
			atom_name = "O" + string( 1, ( '0' + linkage_position + 1 ) );  // Try the next oxygen instead.
			index = residue.atom_index( atom_name );
			TR.Debug << "; cyclic oxygen, finding " << atom_name << " instead.";
		}

		TR.Debug << endl;

		chemical::ResidueType const & type( residue.type() );
		if ( type.atom_forms_residue_connection( index ) ) {
			uint const id_of_connection_to_child( type.residue_connection_id_for_atom( index ) );
			return residue.residue_connection_partner( id_of_connection_to_child );
		} else {
			return 0;
		}
	}
}


core::uint
get_linkage_position_of_saccharide_residue( conformation::Residue const & rsd, conformation::Residue const & parent_rsd )
{
	using namespace std;

	if ( rsd.seqpos() == parent_rsd.seqpos() ) {  // This occurs when there is no parent residue.
		TR.Debug << "This residue is a lower terminus! Returning 0." << endl;
		return 0;
	}
	if ( ! parent_rsd.is_carbohydrate() ) {
		TR.Debug << "Parent is a non-saccharide. Returning 0." << endl;
		return 0;
	}
	if ( rsd.seqpos() == parent_rsd.seqpos() + 1 ) {
		// We have a main chain connection. ( JAB - can we trust this numbering )?

		return parent_rsd.carbohydrate_info()->mainchain_glycosidic_bond_acceptor();
	}

	// If we get this far, we have a branch and need more information.

	//JAB - this looks OK to me - Jason - please check this over!
	// This will be replaced with a graph-search of the ResidueType once I know how to do that.
	core::Size parent_resnum = parent_rsd.seqpos();
	core::Size parent_connect= 0;
	//core::Size res_connect = 0;
	for ( core::Size con = 1; con <= rsd.n_possible_residue_connections(); ++con ) {
		if ( rsd.connected_residue_at_resconn(con) == parent_resnum ) {

			parent_connect = rsd.residue_connection_conn_id(con);
			//res_connect = con;
			break;
		}
	}
	core::Size connect_atom = parent_rsd.residue_connect_atom_index(parent_connect);

	core::Size c_atom_index = 0;

	//Get all bonded neighbors of this atom from the ResidueType.
	//This is important.  Our upper neighbor Carbon is not present, and the only other bonded atom is the H that would be a monosacharide.
	// So, we are good to go here as long as other carbons wouldn't be coming off of the O in any monosacharide (which i don't think is possible and still have a glycan linkage?
	utility::vector1< core::Size > bonded_atoms = parent_rsd.bonded_neighbor( connect_atom );
	for ( core::Size i  = 1; i <= bonded_atoms.size(); ++i ) {
		c_atom_index = bonded_atoms[ i ];
		std::string element = parent_rsd.atom_type( c_atom_index ).element();
		if ( element == "C" ) {
			break;
		}
	}

	std::string connect_atom_name = parent_rsd.atom_name(c_atom_index);
	return utility::string2Size(utility::to_string(connect_atom_name[2]) );
}



// Get whether the glycosidic linkage between the residue and previous residue (parent residue) has an exocyclic carbon.
/// @details  Does not currently work for aa->glycan.  Returns false if previous residue is not carbohydrate.
/// @author   Jared Adolf-Bryfogle (jadolfbr@gmail.com)
bool
has_exocyclic_glycosidic_linkage( conformation::Conformation const & conf, uint const seqpos )
{
	conformation::Residue const & rsd = conf.residue( seqpos );
	core::Size lower_resnum = find_seqpos_of_saccharides_parent_residue( rsd );
	//Lowest of saccharide chain.  Return false.
	if ( lower_resnum == 0 ) {
		TR.Debug << "has_exocyclic_glycosidic_linkage: This residue is a lower terminus! Returning false." << std::endl;
		return false;
	}
	//TR << "lower resnum: " << lower_resnum << std::endl;

	conformation::Residue const & prev_rsd = conf.residue( lower_resnum );
	return has_exocyclic_glycosidic_linkage( rsd, prev_rsd );
}

// Get whether the glycosidic linkage between the residue and previous residue (parent residue) has an exocyclic carbon.
/// @details  Does not currently work for aa->glycan.  Returns false if previous residue is not carbohydrate.
/// @author   Jared Adolf-Bryfogle (jadolfbr@gmail.com)
bool
has_exocyclic_glycosidic_linkage( conformation::Residue const & rsd, conformation::Residue const & parent_rsd )
{
	// What does this mean for ASN-glycan connections?? Technically, it won't be an exocyclic atom -
	// but it WILL have omega and omega2 if ASN - so be careful here!
	if ( ! parent_rsd.is_carbohydrate() ) {
		TR.Debug << "has_exocyclic_glycosidic_linkage: Previous residue is not a carbohydrate! Returning false. " << std::endl;
		return false;
	}

	core::Size const n_carbons = parent_rsd.carbohydrate_info()->n_carbons();
	core::Size linkage_position = get_linkage_position_of_saccharide_residue( rsd, parent_rsd );
	core::Size last_carbon = parent_rsd.carbohydrate_info()->last_carbon_in_ring();

	if ( ( n_carbons == linkage_position ) && ( last_carbon != linkage_position ) ) {
		return true;
	} else {
		return false;
	}
}

/// @brief  Recursive function to get branches of a set of residues, etc.
///  list_of_residues and tips are arrays are non-const references and modified by this function.
///
///  Children Residues:  Residue nums of parent residue connected that we are interested in finding connected branchs.
///  List Of  Residues:  All the residue nums of the branching from children residues
///  Tips:  All 'ends' of all the branches found using this function.
///
///  See Also: get_carbohydrate_residues_and_tips_of_branch
///            trim_carbohydrate_branch_from_X
void
get_branching_residues( conformation::Conformation const & conf,
	Size parent_residue,
	utility::vector1< Size > & children_residues,
	utility::vector1< Size > & list_of_residues,
	utility::vector1< Size > & tips )
{
	for ( core::Size i =1; i <= children_residues.size(); ++i ) {
		Size res = children_residues[ i ];
		utility::vector1< Size > children;
		fill_upstream_children_res_and_tips( conf, res, parent_residue, children, list_of_residues, tips );

		if ( children.size() != 0 ) {
			get_branching_residues( conf, res, children, list_of_residues, tips);
		}
	}
}

/// @brief  Find all children residues, list of residues, and any found tips from a given residue not including parent
///
///  Children Residues:  Filled in list of children residues found if not tips.
///  List Of  Residues:  All the residue nums found.
///  Tips:  All 'ends' of of children found.
///
///  See Also: get_carbohydrate_residues_and_tips_of_branch
///            trim_carbohydrate_branch_from_X
void
fill_upstream_children_res_and_tips( conformation::Conformation const & conf,
	Size res,
	Size parent_residue,
	utility::vector1< Size > & children_residues,
	utility::vector1< Size > & list_of_residues,
	utility::vector1< Size > & tips )
{
	Size connections = conf.residue( res ).n_possible_residue_connections(); //Want the index to match here.
	for ( core::Size connection = 1; connection <= connections; ++connection ) {

		Size connecting_res = conf.residue( res ).connected_residue_at_resconn( connection );

		if ( connecting_res != 0 && connecting_res != parent_residue && conf.residue( connecting_res ).is_carbohydrate() ) {
			list_of_residues.push_back( connecting_res );
			if ( conf.residue( connecting_res ).n_current_residue_connections() == 1 ) {
				tips.push_back( connecting_res );
			} else {
				children_residues.push_back( connecting_res );
			}
		}
	}
}

core::Size
get_glycan_tree_size( conformation::Conformation const & conf, core::Size const first_glycan_resnum ){
	utility::vector1< core::Size > glycan_resnums = get_carbohydrate_residues_of_branch(conf, first_glycan_resnum);
	return glycan_resnums.size() + 1;
}

core::Size
get_largest_glycan_tree_size( conformation::Conformation const & conf ){

	utility::vector1< core::Size > tree_sizes;
	utility::vector1< bool > glycan_start_points = get_glycan_start_points( conf );
	for ( core::Size glycan_start = 1; glycan_start <= glycan_start_points.size(); ++glycan_start ) {
		if ( glycan_start_points[ glycan_start ] ) {
			core::Size glycan_length = get_glycan_tree_size( conf, glycan_start );
			tree_sizes.push_back( glycan_length );
		}
	}
	return utility::max( tree_sizes );
}

core::Size
get_distance_to_start( conformation::Conformation const & conf, core::Size const position){
	core::Size res_distance = 0;
	core::Size parent = find_seqpos_of_saccharides_parent_residue(conf.residue(position));
	if ( parent == 0 ) {
		return 0;
	}
	while ( ( parent != 0 ) && ( conf.residue(parent).is_carbohydrate() ) ) {
		res_distance+=1;
		parent = find_seqpos_of_saccharides_parent_residue(conf.residue(parent));
	}
	return res_distance;
}

utility::vector1< bool >
get_glycan_start_points(conformation::Conformation const & conf){

	utility::vector1< bool > glycan_start_points(conf.size(), false);

	for ( core::Size i = 1; i <= conf.size(); ++i ) {

		//Branch point - definitely a start the glycan
		if ( conf.residue( i ).is_branch_point() && ! conf.residue(i).is_carbohydrate() ) {
			Size glycan_plus_one = get_glycan_connecting_protein_branch_point(conf, i);
			if ( glycan_plus_one != 0 ) {
				glycan_start_points[ glycan_plus_one ] = true;
			}
		} else if ( conf.residue( i ).is_carbohydrate() ) {
			//Indicates a glycan that is not part of a protein chain. - SO we make sure the parent residue is 0 to indicate the start point.
			if ( find_seqpos_of_saccharides_parent_residue( conf.residue(i)) == 0 ) {
				glycan_start_points[ i ] = true;
			}
		}
	}
	return glycan_start_points;
}

/// @brief Get residues further down the branch from this residue.  starting_position ->
/// @details Calls get_carbohydrate_residues_and_tips_of_branch
utility::vector1< core::Size >
get_carbohydrate_residues_of_branch(
	conformation::Conformation const & conf,
	uint const starting_position)
{
	return get_carbohydrate_residues_and_tips_of_branch(conf, starting_position).first;
}

/// @brief Get tips (end residue of linear components of branches) further down the branch from this residue.  starting_position ->
/// @details Calls get_carbohydrate_residues_and_tips_of_branch
utility::vector1< core::Size >
get_carbohydrate_tips_of_branch(
	conformation::Conformation const & conf,
	uint const starting_position)
{
	return get_carbohydrate_residues_and_tips_of_branch(conf, starting_position).second;
}

std::pair< utility::vector1< core::Size >, utility::vector1< core::Size > >
get_carbohydrate_residues_and_tips_of_branch(
	conformation::Conformation const & conf,
	uint const starting_position,
	bool include_starting_position /*false*/)
{
	using namespace core::chemical::carbohydrates;

	utility::vector1< Size > tips;
	utility::vector1< Size > list_of_residues;
	utility::vector1< Size > children_residues;

	if ( include_starting_position ) {
		list_of_residues.push_back(starting_position);
	}


	if ( ! conf.residue( starting_position ).is_carbohydrate() && ! conf.residue( starting_position ).is_branch_point() ) {
		TR << "Delete to residue is not carbohydrate and not a branch point.  Nothing to be done." << std::endl;
		return std::make_pair( list_of_residues, tips);
	}


	Size parent_residue;
	if ( conf.residue( starting_position).is_carbohydrate() ) {
		parent_residue = find_seqpos_of_saccharides_parent_residue( conf.residue( starting_position ) );
	} else {
		parent_residue = starting_position - 1;
	}

	fill_upstream_children_res_and_tips( conf, starting_position, parent_residue, children_residues, list_of_residues, tips );

	//TR << "Children: " << utility::to_string( children_residues) << std::endl;
	get_branching_residues( conf, starting_position, children_residues, list_of_residues, tips );

	return std::make_pair( list_of_residues, tips );
}

///@brief Get the carbohydrate residue connecting the protein branch point.
core::Size
get_glycan_connecting_protein_branch_point(conformation::Conformation const & conf, core::Size const protein_branch_point_resnum){

	debug_assert(conf.residue(protein_branch_point_resnum).is_branch_point());

	core::Size parent_residue = protein_branch_point_resnum - 1;

	Size connections = conf.residue( protein_branch_point_resnum ).n_possible_residue_connections();
	for ( core::Size connection = 1; connection <= connections; ++connection ) {

		Size connecting_res = conf.residue( protein_branch_point_resnum ).connected_residue_at_resconn( connection );

		if ( connecting_res != 0 && connecting_res != parent_residue && conf.residue( connecting_res ).is_carbohydrate() ) {
			return connecting_res;
		}
	}
	return 0;

}

///@brief Get the particular resnum from a glycan position, givin the protein branch point.
/// The glycan_position is numbered 1 -> length of glycan. This is useful for easily identifying a particular glycan position.
///
core::Size
get_resnum_from_glycan_position(conformation::Conformation const & conf, core::Size const glycan_one, core::Size const glycan_position){
	using namespace utility;

	if ( ( glycan_position == 1 )| ( glycan_one == 0 ) ) {
		return glycan_one;
	} else {
		Size glycan_length = get_glycan_tree_size(conf, glycan_one);
		if ( glycan_position <= glycan_length ) {
			core::Size glycan_residue  = glycan_one + glycan_position - 1;
			return glycan_residue;
		} else {
			return 0;
		}
	}

}

core::Size
get_glycan_position_from_resnum(conformation::Conformation const & conf, core::Size const glycan_one, core::Size const glycan_residue ){

	if ( glycan_one == glycan_residue ) {
		return 1;
	} else {
		Size glycan_length = get_glycan_tree_size(conf, glycan_one);
		core::Size glycan_position = glycan_residue + 1 - glycan_one;

		if ( glycan_position > glycan_length ) {
			return 0;
		} else {
			return glycan_position;
		}
	}



}


} //core
} //conformation
} //carbohydrates


