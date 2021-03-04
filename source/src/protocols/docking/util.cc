// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/docking/util.cc

/// @brief Construct a foldtree from a description of the chains
/// involved in the interface

/// @details
/// @author Brian Weitzner
/// @author Matthew O'Meara


// Unit Headers
#include <protocols/docking/util.hh>

// Project headers
#include <core/kinematics/Edge.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/select/residue_selector/ChainSelector.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/types.hh>

// Basic headers
#include <basic/Tracer.hh>

// Numeric headers
#include <numeric/random/random.hh>

// Utility headers
#include <utility/string_util.hh>
#include <utility/vector1.hh>

// C++ Headers
#include <sstream>

static basic::Tracer TR( "protocols.docking.util" );

namespace protocols {
namespace docking {

/// @details Converts a string representation of the partner from the dock_partners flag to a comma-separated list of
/// chains that can be passed into a ChainSeletor.
/// For example, the dock_partners flags "ABC_D" will split into two partners: "ABC" and "D".
/// This function will convert "ABC" to "A,B,C".
/// @return A comma-separated list of chains (e.g. "A,B,C")
std::string comma_separated_partner_chains( std::string const & chains )
{
	using std::endl;
	using std::string;

	// abort if chains is an empty string
	if ( ! chains.size() ) { return chains; }
	if ( TR.Debug.visible() ) { TR.Debug << "Chain group: " << chains << endl; }

	string r;
	r.reserve( ( chains.size() * 2 ) );
	for ( char chain : chains ) {
		r.push_back( chain );
		r.push_back( ',' );
	}

	// remove the last comma
	r.resize( r.size() - 1 );
	if ( TR.Debug.visible() ) { TR.Debug << "Chains to pass to ChainResidueSelector ctor: " << r << endl; }
	return r;
}

/// @details Creates Edges for each of the continuous stretches of residues belonging to a particular partner.
/// If one partner spans multiple chains, the chains will be connected by a Jump.
/// No Edge will span the center_of_mass_residue of the partner.
/// Instead, two Edges will be created: the first will span the beginning of the chain and end at the CoM residue,
/// and the second will begin at the CoM residue and stop at the end of the chain.
/// @return The appropriate FoldTree is configured and accessible by the caller of this function.
void setup_edges_for_partner(
	core::pose::Pose const & pose,
	utility::vector1< bool > const & partner,
	core::Size const center_of_mass_residue,
	core::kinematics::FoldTree & ft )
{
	using utility::vector1;
	using core::Size;
	using core::kinematics::Edge;

	// a separate list of edges is stored to keep track of which edges belong with this partner
	vector1< Edge > edge_list;

	// create an `Edge` for each contiguous stretch of residues in the partner
	bool prev_residue_is_in_partner( false );
	core::Size edge_start( 0 );
	core::Size partner_start( 0 );

	for ( core::Size i = partner.l(); i <= partner.u(); ++i ) {
		bool const residue_is_in_partner( partner[ i ] );

		if ( ! residue_is_in_partner && prev_residue_is_in_partner ) {
			// 'false' after seeing a 'true' - create an `Edge`
			Edge const edge( edge_start, i - 1, Edge::PEPTIDE );
			edge_list.push_back( edge );
			ft.add_edge( edge );

			edge_start = 0;
		} else if ( residue_is_in_partner && prev_residue_is_in_partner && pose.chain( i - 1 ) != pose.chain( i ) ) {
			// new chain - create an `Edge` that ends where the chain ends
			Edge const edge( edge_start, i - 1, Edge::PEPTIDE );
			edge_list.push_back( edge );
			ft.add_edge( edge );

			edge_start = i;
		} else if ( residue_is_in_partner && ! edge_start ) {
			// first occurrence of 'true' - set `edge_start` to this position
			edge_start = i;
			partner_start = i;
		}

		prev_residue_is_in_partner = residue_is_in_partner;
	}

	// account for edges that end with the last residue
	if ( prev_residue_is_in_partner && edge_start ) {
		Edge const edge( edge_start, partner.size(), Edge::PEPTIDE );
		edge_list.push_back( edge );
		ft.add_edge( edge );
	}

	// connect the edges associated with the partner with a jump
	for ( core::Size i = edge_list.l() + 1; i <= edge_list.u(); ++i ) {
		ft.add_edge( edge_list[ i - 1 ].stop(), edge_list[ i ].start(), ft.num_jump() + 1 );
	}

	// adjust the FoldTree to have no edges that span the CoM residue.
	if ( partner_start != center_of_mass_residue ) {
		ft.split_existing_edge_at_residue( center_of_mass_residue );
	}
}

/// @details Creates a Jump connecting the center of mass residues in the FoldTree that is passed through.
/// By default, the dock Jump will be `ft.num_jump() + 1`.
/// This function can optionally ensure that the dock jump has the label "1".
/// @return The label of the dock Jump as a core::Size.
/// @return The appropriate FoldTree is configured and accessible by the caller of this function.
core::Size setup_dock_jump(
	core::Size const partner1_CoM,
	core::Size const partner2_CoM,
	core::kinematics::FoldTree & ft,
	bool const make_dock_jump_label_1 = false )
{
	using core::kinematics::Edge;
	using core::Size;

	core::Size const last_jump_number( ft.num_jump() + 1 );
	core::Size const dock_jump_number( make_dock_jump_label_1 ?  1 : last_jump_number );

	if ( make_dock_jump_label_1 && ft.num_jump() ) {
		// change the label of Jump 1 to be the last Jump so the rigid body docking jump can be "1"
		Edge const & jump_to_update( ft.jump_edge( dock_jump_number ) );
		ft.update_edge_label( jump_to_update.start(), jump_to_update.stop(), jump_to_update.label(), last_jump_number );
	}

	// make the dock jump with the correct label
	ft.add_edge( partner1_CoM, partner2_CoM, dock_jump_number );
	return dock_jump_number;
}

/// @details If partner_chainID is "_", the first chain will be docked to the rest of the complex.
/// If partner_chainID is of the form (pdb_chain_id)+_(pdb_chain_id)+, a jump will be created between the residue
/// nearest to the center of mass of the first partner and the residue nearest to the center of mass of the second
/// partner.
/// For example, "ABC_DEF" has chains "ABC" as the first partner and "DEF" as the second partner.
/// Chains can be listed in any order.
/// @return The constructed FoldTree is set to the pose.
/// @return The movable_jumps vector contains the number of the jump across the interface.
void
setup_foldtree(
	core::pose::Pose & pose,
	std::string const & partner_chainID,
	DockJumps & movable_jumps,
	bool rand_jump_res_partner2 ) // default = false
{
	using std::string;
	using std::stringstream;
	using utility::string_split;
	using utility::vector1;
	using core::Size;
	using core::kinematics::FoldTree;
	using core::select::residue_selector::ChainSelector;

	FoldTree f;
	vector1< bool > partner1( pose.size(), false );
	if ( partner_chainID == "_" ) {
		debug_assert( pose.chain( pose.size() ) > 1 );

		core::Size const last_res_of_first_chain( pose.conformation().chain_end( 1 ) );
		for ( core::Size i = partner1.l(); i <= last_res_of_first_chain; ++i ) { partner1[ i ] = true; }
	} else {
		vector1< string > const partners = string_split( partner_chainID, '_' );

		if ( partners.size() != 2 ) {
			stringstream error_msg;
			error_msg << "Automatic FoldTree setup only works for two-body docking. The value of the partners flag \"";
			error_msg << partner_chainID << "\" implies there are " << partners.size() << " independently movable chains.";
			utility_exit_with_message( error_msg.str() );
		}

		for ( auto const & partner : partners ) {
			if ( partner == "" ) {
				stringstream error_msg;
				error_msg << "Cannot create FoldTree using the provided partners flag \"" << partner_chainID;
				error_msg << "\". At least one of the partner chains is empty.";
				utility_exit_with_message( error_msg.str() );
			}
		}

		ChainSelector const partner1_selector( comma_separated_partner_chains( partners[ 1 ] ) );
		partner1 = partner1_selector.apply( pose );
	}

	setup_foldtree( pose, partner1, movable_jumps, f, rand_jump_res_partner2 );
	pose.fold_tree( f );
}

/// @details This function only supports two-body docking.
/// The vector of boolean values is used to differentiate residues that belong to one partner (a true value) or the
/// other (a false value).
/// Edges are constructed for continuous stretches of residues belonging to a particular partner.
/// If one partner spans multiple chains, the chains will be connected by a Jump.
/// The docking Jump is always set to be Jump number 1.
/// WARNING: This function clears the incoming FoldTree.
/// CHEMICAL Edges are then restored.
/// @return The appropriate FoldTree is configured and accessible by the caller of this function.
/// @return The movable_jumps vector contains the number of the jump across the interface.
void
setup_foldtree(
	core::pose::Pose const & pose,
	utility::vector1< bool > const & partner1,
	DockJumps & movable_jumps,
	core::kinematics::FoldTree & ft,
	bool rand_jump_res_partner2 ) // default = false
{
	using std::endl;
	using utility::vector1;
	using core::Size;
	using core::kinematics::Edge;
	using core::pose::residue_center_of_mass;

	debug_assert( pose.size() );
	debug_assert( partner1.size() == pose.size() );
	TR.Debug << "Setting up the FoldTree for two-body docking" << std::endl;

	// compute which residues belong in partner 2
	vector1< bool > const partner2( partner1.invert() );

	// Save CHEMICAL Edges so that branch connections are not destroyed.
	// Chemical Edge stops will be used to ensure a carbohydrate residue that
	//  is the child of a branch point is not chosen as that breaks the FoldTree
	vector1< Edge > const & chemical_edges( pose.fold_tree().get_chemical_edges() );
	vector1< bool > is_upstream_branch_connect( pose.size() );
	for ( Edge const & chem_edge : chemical_edges ) {
		is_upstream_branch_connect[ chem_edge.stop() ] = true;
	}

	// identify the residue closest to the center of masses of partner1
	Size const jump_pos1( residue_center_of_mass( pose, partner1 ) );

	// identify the residue in partner2 to serve as the Jump residue
	Size jump_pos2;
	// DEFAULT behavior
	// identify the residue closest to the center of mass of partner2
	if ( ! rand_jump_res_partner2 ) {
		jump_pos2 = residue_center_of_mass( pose, partner2 );
	} else {
		// NOT default behavior - must be triggered with flag
		// otherwise, choose a random residue in partner2 for the Jump residue
		// most use cases should not use this behavior
		// This FoldTree setup behavior is used in GlycanDock (@mlnance)
		// glycans (chains of carbohydrates) are ~short, flexible oligomers.
		// by allowing the kinematic propagation to flow in random directions
		// this strategy serves to enhance sampling coverage for glycoligand docking
		// added by @mlnance March 2021
		TR.Debug << "Identifying a random residue from partner2 for the Jump " << std::endl;
		vector1< Size > partner2_resnums;
		for ( Size p2_resnum( 1 ); p2_resnum <= partner2.size(); ++p2_resnum ) {
			// If this residue is part of partner2
			if ( partner2[ p2_resnum ] ) {
				// Don't allow upper connects of a branch point to be
				// a potential Jump residue as they are not functional Jump points
				// (Currently one cannot be a Jump residue AND a chemical Edge)
				if ( ! is_upstream_branch_connect[ p2_resnum ] ) {
					// This is a potential residue to set as the partner2 Jump
					partner2_resnums.push_back( p2_resnum );
				}
			}
		} // done identifying all potential partner2 residues
		// We expect there to be one or more options for the Jump residue
		// so, pick a random residue in partner2 as the Jump point
		if ( ! partner2_resnums.empty() ) {
			jump_pos2 = numeric::random::rg().random_element( partner2_resnums );
		} else {
			// however, if no residues in partner2 were deemed allowable Jump points,
			// (for whatever reason it might be. this probably shouldn't ever happen)
			// fall back to default behavior and identify the residue closest to the CoM
			jump_pos2 = residue_center_of_mass( pose, partner2 );
			TR.Debug << "Partner2 residue selection came up empty "
				"when selecting a random residue as the Jump residue." << std::endl;
			TR.Debug << "Falling back to default behavior of identifying the "
				"residue closest to the center-of-mass of partner2" << std::endl;
		}
	} // END identifying random residue from partner2 as the Jump

	if ( TR.Debug.visible() ) {
		TR.Debug << "jump1: " << jump_pos1 << endl;
		TR.Debug << "jump2: " << jump_pos2 << endl;
	}

	// setup the FoldTree appropriately with the new Jump residues
	// the provided ft and movable_jumps from input are cleared
	ft.clear();
	movable_jumps.clear();
	setup_edges_for_partner( pose, partner1, jump_pos1, ft );
	setup_edges_for_partner( pose, partner2, jump_pos2, ft );
	movable_jumps.push_back( setup_dock_jump( jump_pos1, jump_pos2, ft, true ) );

	ft.delete_self_edges();  // "self-edges" can happen in cases with branches.
	ft.reorder( 1 );

	// Now add back any CHEMICAL Edges, which will have been converted to JUMPs.
	for ( auto const & chemical_edge : chemical_edges ) {
		vector1< Edge > const jump_edges( ft.get_jump_edges() );
		for ( auto const & jump_edge : jump_edges ) {
			if ( jump_edge.stop() == chemical_edge.stop() ) {
				ft.replace_edge( jump_edge, chemical_edge );
				break;
			}
		}
	}
	if ( chemical_edges.size() ) {
		ft.renumber_jumps();  // ...since some were replaced.
		core::uint const current_movable_jump_label
			( ft.edge_label( jump_pos1, jump_pos2 ) );
		if ( current_movable_jump_label != 1 ) {
			// Fix if this is now the case because of renumbering.
			Edge const & current_first_jump( ft.jump_edge( 1 ) );
			ft.update_edge_label( current_first_jump.start(),
				current_first_jump.stop(), 1,
				current_movable_jump_label );
			ft.update_edge_label( jump_pos1, jump_pos2, current_movable_jump_label, 1 );
		}
	}

	debug_assert( ft.check_fold_tree() );
}

} //docking
} //protocols
