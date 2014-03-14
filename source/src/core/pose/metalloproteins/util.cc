// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pose/metalloproteins/util.cc
/// @brief  Pose class utilities for working with metalloproteins.
/// @author Vikram K. Mulligan (vmullig@uw.edu)


// Unit header
#include <core/pose/metalloproteins/util.hh>

// Package headers
#include <core/pose/copydofs/CopyDofs.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/MiniPose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/datacache/PositionConservedResiduesStore.hh>
#include <core/pose/util.tmpl.hh>

// Project headers

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfo.hh>
#include <core/chemical/Patch.hh>
#include <core/chemical/PatchOperation.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/util.hh>
#include <core/id/DOF_ID_Map.hh>
#include <core/id/Exceptions.hh>
#include <core/id/NamedStubID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/TorsionID.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/io/raw_data/DisulfideFile.hh>
#include <core/io/silent/BinaryProteinSilentStruct.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/conformation/Residue.hh>
#include <core/id/SequenceMapping.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/util.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/SingleLigandRotamerLibrary.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableString.hh>
#include <basic/datacache/CacheableStringFloatMap.hh>
#include <basic/datacache/CacheableStringMap.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// Numeric headers
#include <numeric/constants.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.string.hh>

// Utility headers
#include <utility/io/izstream.hh>
#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>

// C/C++ headers
#include <cmath>
#include <iostream>
#include <stdio.h>

// External headers
#include <ObjexxFCL/string.functions.hh>
//#include <boost/functional/hash.hpp>
//#include <boost/foreach.hpp>


#define foreach BOOST_FOREACH


namespace core {
namespace pose {
namespace metalloproteins {

	static basic::Tracer TR("core.pose.util.metalloproteins");

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Function that generates a list of metal-binding atoms that coordinate a metal in a protein.
/// @details This function generates the list by looping through all residues and checking all metal-binding atoms of all
/// metal-binding residues, so it's not super speedy.
/// Inputs:
///		pose (The pose that we'll operate on, unchanged by operation)
///		metal_postion (The residue number of the metal)
/// 	dist_cutoff_multiplier (A float for the distance cutoff multiplier; the cutoff is the sum of the Lennard-Jones radii times the multiplier)
/// @author Vikram K. Mulligan (vmulligan@uw.edu)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1 < core::id::AtomID >
find_metalbinding_atoms (
	core::pose::Pose const &pose,
	core::Size const metal_position,
	core::Real const dist_cutoff_multiplier
) {
	if(!pose.residue(metal_position).is_metal()) {
		std::string message = "Error!  Asked to find metal-binding atoms coordinating something that isn't a metal.";
		utility_exit_with_message(message);
	}

	numeric::xyzVector < core::Real > const &metalatomxyz = pose.residue(metal_position).xyz(1); //Atom 1 should be the metal in any metal residue; this gets its xyz position.

	core::Size const nres = pose.n_residue(); //Number of residues in the pose

	if((metal_position < 1) || (metal_position > nres)) {
		std::string message = "Error!  Asked to find metal-binding atoms coordinating a residue that's not in the pose (metal_position < 1 or > n_residue.";
		utility_exit_with_message(message);
	}

	utility::vector1 < core::id::AtomID > outputlist; //The list of AtomIDs that we will output.

	for(core::Size ir=1; ir<=nres; ++ir) { //Loop through all residues
		if (ir==metal_position) continue; //The metal doesn't bind itself.
		if (pose.residue(ir).is_metalbinding()) { //If this is a metal-binding residue
			utility::vector1 < core::Size > bindingatomlist;
			pose.residue(ir).get_metal_binding_atoms( bindingatomlist );
			if(bindingatomlist.size()>0) {
				for(core::Size ia=1, iamax=bindingatomlist.size(); ia<=iamax; ia++) { //Loop through all the metal-binding atoms in the residue.
					numeric::xyzVector < core::Real > bindingatomxyz = pose.residue(ir).xyz(bindingatomlist[ia]);
					core::Real distsq = metalatomxyz.distance_squared(bindingatomxyz);
					core::Real distcutoffsq = pow(dist_cutoff_multiplier*(pose.residue(metal_position).atom_type(1).lj_radius()+pose.residue(ir).atom_type(bindingatomlist[ia]).lj_radius()), 2);
					if ( distsq <= distcutoffsq ) {
						TR.Debug << "Residue " << ir << " atom " << pose.residue(ir).atom_name( bindingatomlist[ia] ) << " binds the residue " << metal_position << " metal." << std::endl;
						core::id::AtomID curatom( bindingatomlist[ia], ir);
						outputlist.push_back(curatom);
					} //if below distance cutoff
				}//for loop through all metal-binding atoms
			}//if bindingatomlist.size()>0
		}//if is_metalbinding
	}//for loop through all residues

	return outputlist;

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
	utility::vector1 < core::id::AtomID > &liganding_atomids,
	bool const remove_hydrogens
) {
	if(!pose.residue(metal_position).is_metal()) {
		std::string message = "Error!  Asked to add covalent bonds between metal-binding atoms and something that isn't a metal.";
		utility_exit_with_message(message);
	}

	for(core::Size ir=1; ir<=liganding_atomids.size(); ++ir) {
			core::pose::add_covalent_linkage( pose, liganding_atomids[ir].rsd(), metal_position, liganding_atomids[ir].atomno(), 1, remove_hydrogens);
	}

	return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Function to auto-detect and add covalent connections to metal ions.
/// @details This function calls find_metalbinding_atoms and then passes the list of metal binding atoms to add_covalent_linkages_to_metal.
/// Inputs:
///		pose (The pose that we'll operate on, unchanged by operation)
///		metal_postion (The residue number of the metal)
/// 	dist_cutoff_multiplier (A float for the distance cutoff multiplier; the cutoff is the sum of the Lennard-Jones radii times the multiplier)
///   remove_hydrogesn (Should hydrogens on the liganding atoms be auto-removed?  Default true.)
/// @author Vikram K. Mulligan (vmulligan@uw.edu)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
auto_setup_metal_bonds (
	core::pose::Pose &pose,
	core::Size const metal_position,
	core::Real const dist_cutoff_multiplier,
	bool const remove_hydrogens
) {
	utility::vector1 < core::id::AtomID > metalbinding_atomids = find_metalbinding_atoms(pose, metal_position, dist_cutoff_multiplier);
	add_covalent_linkages_to_metal (pose, metal_position, metalbinding_atomids, remove_hydrogens);

	return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Function to auto-detect and add covalent connections to ALL metal ions in a pose.
/// @details This function iteratively calls auto_setup_metal_bonds.
/// Inputs:
///		pose (The pose that we'll operate on, unchanged by operation)
/// 	dist_cutoff_multiplier (A float for the distance cutoff multiplier; the cutoff is the sum of the Lennard-Jones radii times the multiplier)
///   remove_hydrogesn (Should hydrogens on the liganding atoms be auto-removed?  Default true.)
/// @author Vikram K. Mulligan (vmulligan@uw.edu)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
auto_setup_all_metal_bonds (
	core::pose::Pose &pose,
	core::Real const dist_cutoff_multiplier,
	bool const remove_hydrogens
) {

	TR << "Automatically setting covalent bonds between metal ions and metal-binding residues." << std::endl ;

	core::Size nres = pose.n_residue(); //Residue count
	if(nres > 0) {
		for(core::Size ir=1; ir<=nres; ++ir) { //Loop through all residues.
			if(pose.residue(ir).is_metal()) auto_setup_metal_bonds (pose, ir, dist_cutoff_multiplier, remove_hydrogens); //If this is a metal, detect bonding atoms and add covalent connections to them.
		}
	}

	return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Function to set up distance and angle constraints between metals and the residues that bind them.
/// @details This function constrains the distances to be whatever they are in the input pose.  This version
/// does not set the weights for the constraints terms in the scorefunction.
/// Inputs:
///		pose (The pose that we'll operate on, changed by operation)
///   distance_constraint_multiplier (A float for the strength of the metal - binding atom distance constraint.  A value of 2.0 doubles
///   it, for example.)
///   angle_constraint_multiplier (A float for the strength of the metal - binding atom - binding atom parent angle constraint.)
/// @author Vikram K. Mulligan (vmulligan@uw.edu)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
auto_setup_all_metal_constraints (
	core::pose::Pose &pose,
	core::Real const distance_constraint_multiplier,
	core::Real const angle_constraint_multiplier
	//core::scoring::ScoreFunctionOP sfxn
) {
	using namespace core::scoring;
	using namespace core::scoring::func;
	using namespace core::scoring::constraints;
	using namespace core::id;

	TR << "Automatically setting up constraints between metal ions and metal-binding residues." << std::endl ;

	core::Size const nres = pose.n_residue();

	for(core::Size ir=1; ir<=nres; ++ir) { //Loop through all residues
		if(pose.residue(ir).is_metal()) { //When we find a metal, check what it's bound to and set up distance constraints accordingly.
			core::Size const n_resconn = pose.residue(ir).n_residue_connections(); //Number of residue connections

			//Make a vector of indices of virtual atoms in this residue:
			utility::vector1 < core::Size > virtindices;
			for(core::Size ia=1, ia_max=pose.residue(ir).natoms(); ia<=ia_max; ++ia) {
				if(pose.residue(ir).is_virtual(ia)) virtindices.push_back(ia);
			}

			for (core::Size jr=1; jr<=n_resconn; ++jr) { //Loop through all of the metal's residue connections
				if(jr > virtindices.size()) { //We'll use virtual atoms for the constraints.  This means that if we have more connections than virtual residues, we can't do this.
					std::string message = "Error!  The number of connections to a metal is greater than the number of virtual atoms in that metal's residuetype.  Unable to continue.";
					utility_exit_with_message(message);
				}

				core::chemical::ResConnID const & curconnection = pose.residue(ir).connect_map(jr); //The current residue connection
				core::Size const otherres = curconnection.resid(); //The other residue to which this one is connected
				core::Size const otherres_conid = curconnection.connid(); //The other residue's connection id for the bond to the metal.
				core::Size const otherres_atom = pose.residue(otherres).residue_connection(otherres_conid).atomno(); //The atom index of the atom to which this metal is bonded.
				core::Size const otherres_atom_parent = pose.residue(otherres).type().atom_base(otherres_atom); //The atom index of the parent of the atom to which this metal is bonded in the other residue.

				AtomID metalID(1, ir);
				AtomID virtID(virtindices[jr], ir);
				AtomID otherID(otherres_atom, otherres);
				AtomID otherparentID(otherres_atom_parent, otherres);

				pose.set_xyz(virtID, pose.xyz(otherID)); //Move this residue's virt to the other residue's metal-binding atom position.

				//Setting up distance constraints:
				if(distance_constraint_multiplier > 1.0e-10) {
					HarmonicFuncOP hfunc = new HarmonicFunc( 0.0, 0.1 / sqrt(distance_constraint_multiplier)); //Harmonic function for constraining the position.
					AtomPairConstraintOP pairconst = new AtomPairConstraint(virtID, otherID, hfunc); //Atom pair constraint holding the virt at the position of the metal-binding atom.
					pose.add_constraint(pairconst); //Add the constraint to the pose, and carry on.
				}

				//Setting up angle constraints:
				if(angle_constraint_multiplier > 1.0e-10) {
					core::Real const ang1 = numeric::angle_radians( pose.residue(ir).xyz(1), pose.residue(otherres).xyz(otherres_atom),  pose.residue(otherres).xyz(otherres_atom_parent) ); //Angle between metal-bonding atom-bonding atom's parent.
					CircularHarmonicFuncOP circfunc1 = new CircularHarmonicFunc( ang1, 0.05/sqrt(angle_constraint_multiplier) ); //Circular harmonic function for constraining angles (works in RADIANS).
					AngleConstraintOP angleconst1 = new AngleConstraint( metalID, otherID, otherparentID, circfunc1 ); //Angle constraint holding the metal
					pose.add_constraint(angleconst1);
				}
				
			} //Loop through all of the metal's connections
		} //If this is a metal
	} //Loop through all residues

	return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Function to set up distance and angle constraints between metals and the residues that bind them.
/// @details This function constrains the distances to be whatever they are in the input pose.  This version
/// sets the weights for the constraints terms in the scorefunction to 1.0 if they're off, or scales the
/// constraints themselves appropriately if they're already on.
/// Inputs:
///		pose (The pose that we'll operate on, changed by operation)
///   sfxn (An owning pointer to the scorefunction, changed by operation)
///   distance_constraint_multiplier (A float for the strength of the metal - binding atom distance constraint.  A value of 2.0 doubles
///   it, for example.)
///   angle_constraint_multiplier (A float for the strength of the metal - binding atom - binding atom parent angle constraint.)
/// @author Vikram K. Mulligan (vmulligan@uw.edu)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
auto_setup_all_metal_constraints (
	core::pose::Pose &pose,
	core::scoring::ScoreFunctionOP sfxn,
	core::Real const distance_constraint_multiplier,
	core::Real const angle_constraint_multiplier
) {
	using namespace core::scoring;

	core::Real distmult = 1.0; //Constraint weight multiplier, for cases in which the atom_pair_constraint has already been turned on but is not equal to 1.0.
	if(sfxn->get_weight(atom_pair_constraint) < 1.0e-10) { //If this score term is turned off, turn it on.
		sfxn->set_weight(atom_pair_constraint, 1.0);
	} else { //Otherwise, if this score term is already turned on, set the weight multiplier appropriately.
		distmult = 1.0 / sfxn->get_weight(atom_pair_constraint); //Note -- we will take the square root of this later on.
	}

	core::Real angmult=1.0; //Constraint weight multiplier, for cases in which the atom_pair_constraint has already been turned on but is not equal to 1.0.
	if(sfxn->get_weight(angle_constraint) < 1.0e-10) { //If this score term is turned off, turn it on.
		sfxn->set_weight(angle_constraint, 1.0);
	} else { //Otherwise, if this score term is already turned on, set the weight multiplier appropriately.
		angmult = 1.0 / sfxn->get_weight(angle_constraint); //Note -- we will take the square root of this later on.
	}

	auto_setup_all_metal_constraints ( pose, distmult*distance_constraint_multiplier, angmult*angle_constraint_multiplier );

	return;
}
	
} //metalloproteins
} // pose
} // core
