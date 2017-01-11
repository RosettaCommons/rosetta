// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pose/metalloproteins/util.cc
/// @brief  Utilities for working with metalloproteins.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
/// @author Andy Watkins (amw579@nyu.edu)


// Unit header
#include <core/util/metalloproteins_util.hh>

// Package headers
#include <core/pose/copydofs/CopyDofs.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

// Project headers
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
// AMW XRW TODO remove dependence on cache
#include <core/chemical/ResidueTypeSetCache.hh>
#include <core/chemical/Patch.hh>
#include <core/chemical/PatchOperation.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/id/AtomID.hh>

#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <core/scoring/func/ScalarWeightedFunc.hh>
#include <core/scoring/constraints/NamedAtomPairConstraint.hh>
#include <core/scoring/constraints/NamedAngleConstraint.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// Numeric headers
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>

// C/C++ headers
#include <cmath>
#include <iostream>

// External headers
#include <ObjexxFCL/string.functions.hh>

namespace core {
namespace util {

static THREAD_LOCAL basic::Tracer TR( "core.util.metalloproteins_util" );

// removed from the header file so that it cannot be called directly, only from this file ~Labonte
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief: This is a helper function for the add_covalent_linkage function.
/// @details:  This is useful for adding covalent linkages between metal-binding
/// side-chains and metal atoms.  This code was shamelessly stolen from
/// Florian's EnzConstraintParameters.cc but now does not modify the RTS
/// and would be threadsafe.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
/// @author Florian Richter (flosopher@gmail.com)
/// @author Andy Watkins (andy.watkins2@gmail.com)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
add_covalent_linkage_helper(
	core::pose::Pose & pose,
	core::Size const res_pos,
	core::Size const Atpos,
	numeric::xyzVector< core::Real > const & /*partner_xyz*/, //Coordinates of the atom to which this will be bonded
	bool const remove_hydrogens
)
{
	using namespace core;
	using namespace chemical;

	conformation::Residue const & res = pose.residue( res_pos );

	std::string res_atom = res.atom_name( Atpos );
	ObjexxFCL::strip_whitespace( res_atom );

	std::string current_residue_type_basename( residue_type_base_name( res.type() ) );
	std::string current_residue_type_patches_name( residue_type_all_patches_name( res.type() ) );

	std::string res_patchname( "MP-" + res_atom + "-connect" );
	if ( remove_hydrogens ) { res_patchname += ":MP-" + res_atom + "-pruneH"; }
	if ( res.type().is_metal() ) { res_patchname = "MP-" + res_atom + "-metal_connect"; }

	std::string res_type_mod_name( current_residue_type_basename + ':' + res_patchname + current_residue_type_patches_name );

	conformation::Residue new_res( core::pose::get_restype_for_pose( pose, res_type_mod_name, res.type().mode() ), true);

	// Temporarily make a copy of the old residue:
	conformation::Residue old_res = ( *res.clone() );

	// replacing the residue, using first three atoms alone for metals.
	if ( res.is_metal() ) {
		utility::vector1< std::pair< std::string, std::string > > atom_pairs;
		for ( core::Size ia=1; ia<=3; ia++ ) atom_pairs.emplace_back( res.atom_name(ia),new_res.atom_name(ia) );
		pose.replace_residue( res_pos, new_res, atom_pairs);
	} else {
		pose.replace_residue( res_pos, new_res, true);
	}

	// and resetting the xyz positions
	for ( core::Size at_ct = 1, at_ctmax = old_res.natoms(); at_ct <= at_ctmax; ++at_ct ) {
		// BOTH residues must have it.
		if ( new_res.has( old_res.atom_name( at_ct ) ) ) {
			pose.set_xyz( id::AtomID( at_ct, res_pos ), old_res.xyz( old_res.atom_name( at_ct ) ) );
		}
	}

} //add_covalent_linkage_helper


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief: Adds an arbitrary covalent linkage between two atoms (resA_At and resB_At) in two residues (at positions resA_pos and resB_pos).
/// @details:  This is useful for adding covalent linkages between metal-binding side-chains and metal atoms.  This code was shamelessly
/// stolen from Florian's EnzConstraintParameters.cc in protocols/toolbox/match_enzdes_utils, and was modified to permit deletion of
/// unnecessary protons.
/// @author:  Vikram K. Mulligan (vmullig@uw.edu), Florian Richter (flosopher@gmail.com)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
add_covalent_linkage(
	core::pose::Pose & pose,
	Size const resA_pos,
	Size const resB_pos,
	Size const resA_At,
	Size const resB_At,
	bool const remove_hydrogens //Should extraneous hydrogens on the bonding atoms be removed?
)
{
	using namespace core::chemical;

	TR << "Adding covalent linkage between residue " << resA_pos << "'s " << pose.residue(resA_pos).atom_name(resA_At).c_str() <<
		" atom and residue " << resB_pos << "'s " << pose.residue(resB_pos).atom_name(resB_At).c_str() << " atom." << std::endl;

	add_covalent_linkage_helper( pose, resA_pos, resA_At, pose.residue(resB_pos).xyz(resB_At), remove_hydrogens);
	add_covalent_linkage_helper( pose, resB_pos, resB_At, pose.residue(resA_pos).xyz(resA_At), remove_hydrogens);

	std::string resA_atomname = pose.residue( resA_pos ).atom_name( resA_At );
	std::string resB_atomname = pose.residue( resB_pos ).atom_name( resB_At );

	pose.conformation().declare_chemical_bond(
		resA_pos, resA_atomname,
		resB_pos, resB_atomname
	);

} //add_covalent_linkage


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Function that generates a list of metal-binding atoms that coordinate a metal in a protein.
/// @details This function generates the list by looping through all residues and checking all metal-binding atoms of all
/// metal-binding residues, so it's not super speedy.
/// Inputs:
///  pose (The pose that we'll operate on, unchanged by operation)
///  metal_postion (The residue number of the metal)
///  dist_cutoff_multiplier (A float for the distance cutoff multiplier; the cutoff is the sum of the Lennard-Jones radii times the multiplier)
/// @author Vikram K. Mulligan (vmulligan@uw.edu)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< core::id::AtomID >
find_metalbinding_atoms (
	core::pose::Pose const &pose,
	core::Size const metal_position,
	core::Real const dist_cutoff_multiplier
) {
	using namespace core;

	if ( !pose.residue( metal_position ).is_metal() ) {
		std::string message = "Error!  Asked to find metal-binding atoms coordinating something that isn't a metal.";
		utility_exit_with_message(message);
	}

	// Metal residues have a single atom
	numeric::xyzVector< core::Real > const & metal_xyz = pose.residue( metal_position ).xyz( 1 );

	Size const nres = pose.size();

	if ( ( metal_position < 1 ) || ( metal_position > nres ) ) {
		utility_exit_with_message( "Error!  Asked to find metal-binding atoms coordinating a residue that's not in the pose (metal_position < 1 or > n_residue." );
	}

	utility::vector1< id::AtomID > coordinating_atoms;

	// Loop through all metalbinding residues, skipping the metal itself
	// AMW: does not necessarily avoid ALL metal positions, just the one in question
	for ( Size ir = 1; ir <= nres; ++ir ) {
		if ( ir == metal_position ) continue;
		if ( ! pose.residue( ir ).is_metalbinding() ) continue;

		utility::vector1< Size > binding_atom_list;
		pose.residue( ir ).get_metal_binding_atoms( binding_atom_list );

		for ( Size const binding_atom : binding_atom_list ) {
			numeric::xyzVector< Real > binding_atom_xyz = pose.residue( ir ).xyz( binding_atom );
			Real distsq = metal_xyz.distance_squared( binding_atom_xyz );
			Real metal_rad = pose.residue( metal_position ).atom_type( 1 ).lj_radius();
			Real lig_rad = pose.residue( ir ).atom_type( binding_atom ).lj_radius();
			Real distcutoffsq = dist_cutoff_multiplier * ( metal_rad + lig_rad )
				* dist_cutoff_multiplier * ( metal_rad + lig_rad );

			if ( distsq > distcutoffsq ) continue;

			TR.Debug << "Residue " << ir << " atom " << pose.residue( ir ).atom_name( binding_atom ) << " binds the residue " << metal_position << " metal." << std::endl;
			id::AtomID curatom( binding_atom, ir );
			coordinating_atoms.push_back( curatom );

		}
	}

	return coordinating_atoms;

} //find_metalbinding_atoms
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Function to add covalent linkages between a metal atom and all the liganding atoms provided in a vector of AtomIDs.
/// @details The inputs are:
///    pose (The pose to be modified)
///    metal_position (The residue number of the metal in the pose)
///    liganding_atomids (A list of AtomIDs on other residues that will be covalently linked to the metal)
///    remove_hydrogens (Should hydrogens on the liganding atoms be removed automatically?  Default true.)
/// This function uses core::pose::add_covalent_linkage, which can strip off extraneous hydrogens.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
add_covalent_linkages_to_metal (
	core::pose::Pose &pose,
	core::Size const metal_position,
	utility::vector1< core::id::AtomID > & liganding_atomids,
	bool const remove_hydrogens
) {
	if ( !pose.residue(metal_position).is_metal() ) {
		std::string message = "Error!  Asked to add covalent bonds between metal-binding atoms and something that isn't a metal.";
		utility_exit_with_message(message);
	}

	for ( auto const & liganding_atomid : liganding_atomids ) {
		add_covalent_linkage( pose, liganding_atomid.rsd(), metal_position, liganding_atomid.atomno(), 1, remove_hydrogens);
	}

	return;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Function to auto-detect and add covalent connections to ALL metal ions in a pose.
/// @details This function iteratively calls auto_setup_metal_bonds.
/// Inputs:
///  pose (The pose that we'll operate on, unchanged by operation)
///  dist_cutoff_multiplier (A float for the distance cutoff multiplier; the cutoff is the sum of the Lennard-Jones radii times the multiplier)
///   remove_hydrogesn (Should hydrogens on the liganding atoms be auto-removed?  Default true.)
/// @author Vikram K. Mulligan (vmulligan@uw.edu)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
auto_setup_all_metal_bonds(
	core::pose::Pose &pose,
	core::Real const dist_cutoff_multiplier,
	bool const remove_hydrogens
) {

	TR << "Automatically setting covalent bonds between metal ions and metal-binding residues." << std::endl ;

	for ( core::Size ir=1; ir<=pose.size(); ++ir ) { //Loop through all residues.
		if ( pose.residue(ir).is_metal() ) {
			utility::vector1< core::id::AtomID > metalbinding_atomids = find_metalbinding_atoms( pose, ir, dist_cutoff_multiplier );
			add_covalent_linkages_to_metal( pose, ir, metalbinding_atomids, remove_hydrogens);
		}
	}

	return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Function to set up distance and angle constraints between metals and
/// the residues that bind them.
/// @details This function constrains the distances to be whatever they are in
/// the input pose.  This version does not set the weights for the constraints
/// terms in the scorefunction.
/// Inputs:
///  pose (The pose that we'll operate on, changed by operation)
///  distance_constraint_multiplier (A float for the strength of the metal - binding atom distance constraint.  A value of 2.0 doubles it, for example.)
///  angle_constraint_multiplier (A float for the strength of the metal - binding atom - binding atom parent angle constraint.)
/// @author Vikram K. Mulligan (vmulligan@uw.edu)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
auto_setup_all_metal_constraints(
	core::pose::Pose &pose,
	core::Real const distance_constraint_multiplier,
	core::Real const angle_constraint_multiplier
	//core::scoring::ScoreFunctionOP sfxn
) {
	using namespace core;
	using namespace scoring;
	using namespace scoring::func;
	using namespace scoring::constraints;
	using namespace id;

	TR << "Automatically setting up constraints between metal ions and metal-binding residues." << std::endl ;

	for ( Size ir = 1; ir <= pose.size(); ++ir ) {
		conformation::Residue const & ir_res = pose.residue( ir );
		if ( ! ir_res.is_metal() ) continue;
		Size const ir_nconn = ir_res.n_possible_residue_connections();

		// Make a vector of indices of virtual atoms in this residue:
		utility::vector1< std::string > ir_virt_names;
		utility::vector1< Size > ir_virt_atomnos;
		for ( Size ia = 1, ia_max = ir_res.natoms(); ia <= ia_max; ++ia ) {
			if ( ir_res.is_virtual( ia ) ) {
				ir_virt_atomnos.push_back( ia );
				ir_virt_names.push_back( ir_res.atom_name( ia ) );
			}
		}

		// We must have at least as many virt atoms as connections
		if ( ir_nconn > ir_virt_names.size() ) {
			std::stringstream ss;
			ss << "Error!  The number of connections to a metal, " << ir_nconn << ", is greater than the number of virtual atoms , " << ir_virt_names.size() << ", in that metal's residuetype, " << ir_res.name() << ".  Unable to continue.";
			std::string const & message = ss.str();
			utility_exit_with_message(message);
		}

		// Constrain all the metal's residue connections
		// Note: since allowing design we must use Named*Constraints
		// because atom indices may change.
		for ( Size jr = 1; jr <= ir_nconn; ++jr ) {
			chemical::ResConnID const & jr_conn = ir_res.connect_map( jr );
			Size const conn_res = jr_conn.resid();
			conformation::Residue const & jr_res = pose.residue( conn_res );
			Size const conn_res_conid = jr_conn.connid(); //The other residue's connection id for the bond to the metal.
			Size const conn_res_atomno = jr_res.residue_connection(conn_res_conid).atomno(); //The atom index of the atom to which this metal is bonded.
			std::string const conn_res_atom = jr_res.atom_name( conn_res_atomno ); //The atom index of the atom to which this metal is bonded.
			std::string const conn_res_atom_parent = jr_res.atom_name( jr_res.type().atom_base( conn_res_atomno ) ); //The atom index of the parent of the atom to which this metal is bonded in the other residue.

			// Note: assumes that metal residue types have the first atom as the metal.
			NamedAtomID metalID( ir_res.atom_name( 1 ), ir );
			NamedAtomID virtID(ir_virt_names[jr], ir);
			NamedAtomID otherID(conn_res_atom, conn_res);
			NamedAtomID otherparentID(conn_res_atom_parent, conn_res);

			pose.set_xyz(virtID, pose.xyz(otherID)); //Move this residue's virt to the other residue's metal-binding atom position.

			// Setting up distance constraints:
			if ( distance_constraint_multiplier > 1.0e-10 ) {
				FuncOP hfunc( new ScalarWeightedFunc( distance_constraint_multiplier, FuncOP(new HarmonicFunc( 0.0, 0.1 )) ));
				//Atom pair constraint holding the virt at the position of the metal-binding atom.
				NamedAtomPairConstraintOP pairconst(
					new NamedAtomPairConstraint(virtID, otherID, hfunc, core::scoring::metalbinding_constraint) );
				pose.add_constraint(pairconst); //Add the constraint to the pose, and carry on.
			}

			// Setting up angle constraints
			if ( angle_constraint_multiplier > 1.0e-10 ) {
				core::Real const ang1 = numeric::angle_radians( pose.residue(ir).xyz(1), pose.residue(conn_res).xyz(conn_res_atom),  pose.residue(conn_res).xyz(conn_res_atom_parent) ); //Angle between metal-bonding atom-bonding atom's parent.
				//Circular harmonic function for constraining angles (works in RADIANS).
				FuncOP circfunc1( new ScalarWeightedFunc( angle_constraint_multiplier, FuncOP(new CircularHarmonicFunc( ang1, 0.05 )) ));
				NamedAngleConstraintOP angleconst1(
					new NamedAngleConstraint( metalID, otherID, otherparentID, circfunc1, core::scoring::metalbinding_constraint ) );
				pose.add_constraint(angleconst1);
			}
		} //Loop through all of the metal's connections

		// set up coordinate constraints to the generated virtuals (for cartrefine)
		if ( distance_constraint_multiplier > 1.0e-10 ) {
			for ( Size jr = 1; jr <= ir_nconn; ++jr ) {
				NamedAtomID metalID( ir_res.atom_name( 1 ), ir );
				NamedAtomID virtID(ir_virt_names[jr], ir);
				core::Real const mvdist =
					(pose.residue(ir).xyz(ir_virt_atomnos[jr]) - pose.residue(ir).xyz(1)).length();
				FuncOP mvfunc( new ScalarWeightedFunc( distance_constraint_multiplier, FuncOP(new HarmonicFunc( mvdist, 0.1 )) ));
				NamedAtomPairConstraintOP pairconst_mv(
					new NamedAtomPairConstraint(metalID, virtID, mvfunc, core::scoring::metalbinding_constraint) );
				pose.add_constraint(pairconst_mv);

				for ( Size kr = jr+1; kr <= ir_nconn; ++kr ) {
					NamedAtomID virtID_j(ir_virt_names[kr], ir);
					core::Real const vvdist =
						(pose.residue(ir).xyz(ir_virt_atomnos[jr]) - pose.residue(ir).xyz(ir_virt_atomnos[kr])).length();
					FuncOP vvfunc( new ScalarWeightedFunc( distance_constraint_multiplier, FuncOP(new HarmonicFunc( vvdist, 0.1 )) ));
					NamedAtomPairConstraintOP pairconst_vv(
						new NamedAtomPairConstraint(virtID, virtID_j, vvfunc, core::scoring::metalbinding_constraint) );
					pose.add_constraint(pairconst_vv);
				}
			}
		}
	} //Loop through all residues
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Function to set up distance and angle constraints between metals and the residues that bind them.
/// @details This function constrains the distances to be whatever they are in the input pose.  This version
/// sets the weights for the constraints terms in the scorefunction to 1.0 if they're off, or scales the
/// constraints themselves appropriately if they're already on.
/// Inputs:
///  pose (The pose that we'll operate on, changed by operation)
///   sfxn (An owning pointer to the scorefunction, changed by operation)
///   distance_constraint_multiplier (A float for the strength of the metal - binding atom distance constraint.  A value of 2.0 doubles
///   it, for example.)
///   angle_constraint_multiplier (A float for the strength of the metal - binding atom - binding atom parent angle constraint.)
/// @author Vikram K. Mulligan (vmulligan@uw.edu)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
auto_setup_all_metal_constraints(
	core::pose::Pose &pose,
	core::scoring::ScoreFunctionOP sfxn,
	core::Real const distance_constraint_multiplier,
	core::Real const angle_constraint_multiplier
) {
	using namespace core::scoring;

	if ( std::fabs( sfxn->get_weight(metalbinding_constraint) ) < 1.0e-10 ) {
		sfxn->set_weight(metalbinding_constraint, 1.0);
	}

	auto_setup_all_metal_constraints( pose, distance_constraint_multiplier, angle_constraint_multiplier );

	return;
}

} // util
} // core
