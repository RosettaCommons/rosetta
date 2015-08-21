// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/magnesium/minimize_util.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/magnesium/minimize_util.hh>
#include <protocols/magnesium/util.hh>
#include <protocols/magnesium/MgMinimizer.hh>
#include <protocols/stepwise/setup/FullModelInfoSetupFromCommandLine.hh>
#include <core/chemical/rna/util.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.magnesium.minimize_util" );

using namespace core;
using utility::vector1;

//////////////////////////////////////////////////////////////////////////////////
//
// Unclear if this code is going to be useful in long term
// May want to encapsulate into a MgMinimizer class.
//
//                 -- rhiju, june 2015
//
//////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace magnesium {

////////////////////////////////////////////////////////////////////////////////
// @brief useful function recording choices for how to minimize Mg(2+)
////////////////////////////////////////////////////////////////////////////////
void
minimize_magnesium_and_hydration_shell( pose::Pose & pose,
	utility::vector1< Size > const mg_res,
	core::scoring::ScoreFunctionCOP minimize_scorefxn /* = 0 */,
	core::Distance const mg_coord_cst_dist /* = 0.2 */ ) {
	MgMinimizer mg_minimizer;
	mg_minimizer.set_mg_res( mg_res );
	mg_minimizer.set_minimize_scorefxn( minimize_scorefxn );
	mg_minimizer.set_mg_coord_cst_dist( mg_coord_cst_dist );
	mg_minimizer.apply( pose );
}

//////////////////////////////////////////
void
minimize_magnesium_and_hydration_shell( pose::Pose & pose /*for viewing*/,
	utility::vector1< pose::PoseOP > & pose_list,
	utility::vector1< Size > const mg_res,
	core::scoring::ScoreFunctionCOP minimize_scorefxn /* = 0 */,
	core::Distance const mg_coord_cst_dist /* = 0.2 */ ) {
	MgMinimizer mg_minimizer;
	mg_minimizer.set_mg_res( mg_res );
	mg_minimizer.set_minimize_scorefxn( minimize_scorefxn );
	mg_minimizer.set_mg_coord_cst_dist( mg_coord_cst_dist );
	mg_minimizer.apply( pose_list, pose );
}


//////////////////////////////////////////////////////////////////////////////////
// @brief sets up all jumps for Mg(2+) and HOH in a pose.
//
// @details This is kind of extreme -- total rewiring of fold tree at all Mg(2+) and
//  all HOH. Would be better to break down into fold_tree remodeling for each Mg(2+),
//  perhaps with a final cleanup for HOH. Will do that if a good use case
//  becomes apparent.
//////////////////////////////////////////////////////////////////////////////////
void
update_mg_hoh_fold_tree( pose::Pose & pose ){
	using namespace core::kinematics;
	using namespace core::id;
	using namespace protocols::magnesium;
	FoldTree const & f( pose.fold_tree() );

	// get jumps, except for mg(2+) & hoh
	// get cuts
	vector1< Size > jump_partners1, jump_partners2, cuts, water_res, mg_res;
	vector1< std::string > jump_atoms1, jump_atoms2;
	vector1< std::pair< Size, std::string > > is_connected( pose.total_residue(), std::make_pair(0,"") );
	for ( Size n = 1; n <= f.num_jump(); n++ ) {
		Size const i = f.upstream_jump_residue( n );
		Size const j = f.downstream_jump_residue( n );
		bool i_is_mg_or_hoh_res = ( pose.residue( i ).name3() == "HOH" || pose.residue( i ).name3() == " MG" );
		bool j_is_mg_or_hoh_res = ( pose.residue( j ).name3() == "HOH" || pose.residue( j ).name3() == " MG" );
		if ( i_is_mg_or_hoh_res && j_is_mg_or_hoh_res ) continue;
		if ( i_is_mg_or_hoh_res && !j_is_mg_or_hoh_res ) {
			std::pair< Size, std::string > const & connection( is_connected[ i ] );
			if (  connection.first > 0 ) {
				jump_partners1.push_back( connection.first );
				jump_partners2.push_back( j );
				jump_atoms1.push_back( connection.second );
				jump_atoms2.push_back( chemical::rna::default_jump_atom( pose.residue( j ) ) );
			} else {
				is_connected[ i ] = std::make_pair( j, chemical::rna::default_jump_atom( pose.residue( j ) )  );
			}
			continue;
		}
		if ( j_is_mg_or_hoh_res && !i_is_mg_or_hoh_res ) {
			std::pair< Size, std::string > const & connection( is_connected[ j ] );
			if ( connection.first > 0 ) {
				jump_partners1.push_back( i );
				jump_partners2.push_back( connection.first );
				jump_atoms1.push_back( chemical::rna::default_jump_atom( pose.residue( i ) ) );
				jump_atoms2.push_back( connection.second );
			} else {
				is_connected[ j ] = std::make_pair( i, chemical::rna::default_jump_atom( pose.residue( i ) ) );
			}
			continue;
		}
		jump_partners1.push_back( i );
		jump_partners2.push_back( j );
		jump_atoms1.push_back( f.upstream_atom( n ) );
		jump_atoms2.push_back( f.downstream_atom( n ) );
	}
	for ( Size n = 1; n <= pose.total_residue(); n++ ) {
		if ( f.is_cutpoint( n ) && n < pose.total_residue() ) cuts.push_back( n );
		if ( pose.residue( n ).name3() == "HOH" ) water_res.push_back( n );
		if ( pose.residue( n ).name3() == " MG" ) mg_res.push_back( n );
	}

	// setup jumps for all mg(2+), and for all hoh near mg(2+)
	vector1< bool > assigned_jump_to_water( pose.total_residue(), false );
	for ( Size n = 1; n <= mg_res.size(); n++ ) {
		Size const i = mg_res[ n ];
		vector1< AtomID > mg_ligands = get_mg_ligands( pose, i );
		bool found_non_hoh_ligand( false );
		AtomID mg_partner_id;
		for ( Size k = 1; k <= mg_ligands.size(); k++ ) {
			Size const j = mg_ligands[ k ].rsd();
			if ( pose.residue( j ).name3() == "HOH" ) {
				if ( !assigned_jump_to_water[ j ] ) {
					assigned_jump_to_water[ j ] = true;
					jump_partners1.push_back( i );
					jump_partners2.push_back( j );
					jump_atoms1.push_back( "MG  " );
					jump_atoms2.push_back( " O  " );
				}
			} else {
				if ( !found_non_hoh_ligand ) mg_partner_id = AtomID( mg_ligands[ k ].atomno(), j );
			}
		}
		if ( !found_non_hoh_ligand ) mg_partner_id = get_closest_non_hoh_contact( pose, i, " MG" );
		if ( mg_partner_id.rsd() == 0 ) continue; // special case, just Mg2+ and HOH.

		jump_partners1.push_back( i );
		jump_partners2.push_back( mg_partner_id.rsd() );
		jump_atoms1.push_back( "MG  " );
		jump_atoms2.push_back( pose.residue( mg_partner_id.rsd() ).atom_name( mg_partner_id.atomno() ) );
	}

	// setup jumps for all remaining hoh -- make sure they do not stick to mg2+ or mg-bound hoh, so
	// that they can move separately.
	for ( Size n = 1; n <= water_res.size(); n++ ) {
		Size const i = water_res[ n ];
		if ( assigned_jump_to_water[ i ] ) continue;
		// find closest hbond partner
		AtomID best_partner = get_closest_non_hoh_contact( pose, i );
		jump_partners1.push_back( i );
		jump_partners2.push_back( best_partner.rsd() );
		jump_atoms1.push_back( "MG  " );
		jump_atoms2.push_back( pose.residue( best_partner.rsd() ).atom_name( best_partner.atomno() ) );
	}

	// choose fold tree and put into pose.
	runtime_assert( cuts.size() == jump_partners1.size() );
	FoldTree f_new = protocols::stepwise::setup::get_tree( pose.total_residue(), cuts, jump_partners1, jump_partners2, jump_atoms1, jump_atoms2 );
	pose.fold_tree( f_new );
}

////////////////////////////////
// @brief helper function for update_mg_hoh_fold_tree();
core::id::AtomID
get_closest_non_hoh_contact( pose::Pose const & pose, Size const i, std::string const & exclude_rsd /* = "" */ ) {
	using namespace core::id;
	using namespace core::conformation;
	Vector const & i_xyz( pose.residue( i ).xyz( 1 ) );
	AtomID best_partner;
	Distance best_d = ( pose.residue( 1 ).xyz( 1 ) - i_xyz ).length(); // arbitrary init
	for ( Size j = 1; j <= pose.total_residue(); j++ ) {
		Residue const & rsd = pose.residue( j );
		if ( rsd.name3() == "HOH" ) continue;
		if ( rsd.name3() == exclude_rsd ) continue;
		for ( Size jj = 1; jj <= rsd.nheavyatoms(); jj++ ) {
			if ( rsd.heavyatom_is_an_acceptor( jj ) || rsd.heavyatom_has_polar_hydrogens( jj ) ) {
				Distance const d = ( rsd.xyz( jj ) - i_xyz ).length();
				if ( d < best_d ) {
					best_partner = AtomID( jj, j );
					best_d = d;
				}
			}
		}
	}
	return best_partner;
}

} //magnesium
} //protocols
