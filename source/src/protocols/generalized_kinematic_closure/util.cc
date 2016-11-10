// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/generalized_kinematic_closure/util.cc
/// @brief  Utility functions for generalized kinematic closure.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Unit Headers
#include <protocols/generalized_kinematic_closure/GeneralizedKIC.hh>
#include <protocols/generalized_kinematic_closure/GeneralizedKICCreator.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>

#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <utility/tag/Tag.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>

#include <core/pose/util.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <numeric/random/random.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzVector.io.hh>
#include <numeric/kinematic_closure/bridgeObjects.hh>
#include <numeric/kinematic_closure/kinematic_closure_helpers.hh>
#include <numeric/conversions.hh>

#include <boost/foreach.hpp>

//Auto Headers
#include <utility/excn/Exceptions.hh>
#include <core/pose/Pose.hh>


using basic::T;
using basic::Error;
using basic::Warning;

namespace protocols {
namespace generalized_kinematic_closure {

static THREAD_LOCAL basic::Tracer TR( "protocols.generalized_kinematic_closure.util" );


/// @brief Function to determine whether a value is in a list.
bool is_in_list (core::Size const val, utility::vector1 < core::Size > const &list)
{
	core::Size const listsize = list.size();
	if ( listsize==0 ) return false;
	for ( core::Size i=1; i<=listsize; ++i ) {
		if ( list[i]==val ) return true;
	}

	return false;
}


/// @brief function to determine whether a residue index from the original pose is in the residue_map list.
/// @details Returns true if the residue index in the original pose has been mapped to a residue in the loop, false otherwise.
/// @param[in] residue_index - The index of a residue in the original pose.
/// @param[in] residue_map - A vector of pairs of indices.  Each pair is (index in loop, index in original pose).
bool original_pose_residue_is_in_residue_map (core::Size const residue_index, utility::vector1 < std::pair<core::Size, core::Size> > const &residue_map) {
	core::Size const residue_map_size=residue_map.size();
	if ( residue_map_size==0 ) return false;
	for ( core::Size i=1; i<=residue_map_size; ++i ) { //Loop through the residue map.
		if ( residue_map[i].second == residue_index ) return true;
	}
	return false;
}


/// @brief Set the loop pose conformation based on a set of results from kinematic closure.
/// @details
/// @param[in,out] pose -- A pose consisting of the loop to be closed only.
/// @param[in] atomlist -- A list of AtomIDs and corresponding xyz coordinates (though we don't use the latter) of the chain of atoms closed by KIC.  Note that the residue indices refer to the loop pose, not the original pose.
/// @param[in] t_ang -- The torsion angles values to set.
/// @param[in] b_ang -- The bond angle values to set.
/// @param[in] b_len -- The bond length values to set.
void set_loop_pose (
	core::pose::Pose &pose,
	utility::vector1 < std::pair < core::id::AtomID, numeric::xyzVector<core::Real> > > const &atomlist, //I want this to be const access
	utility::vector1 < core::Real > const &t_ang,
	utility::vector1 < core::Real > const &b_ang,
	utility::vector1 < core::Real > const &b_len
) {
	using namespace numeric::conversions;
	using namespace core::id;

	if ( TR.Debug.visible() ) { TR.Debug << "Setting loop conformation based on solution from generalized kinematic closure." << std::endl; TR.Debug.flush(); }

	//pose.dump_pdb("before.pdb"); //DELETE ME -- for debugging only.

	for ( core::Size ii=1, iimax=(b_len.size()-1); ii<=iimax; ++ii ) { //Set bond lengths
		pose.conformation().set_bond_length( atomlist[ii].first, atomlist[ii+1].first, b_len[ii] );
	}
	for ( core::Size ii=2, iimax=(b_ang.size()-1); ii<=iimax; ++ii ) { //Set bond angles
		pose.conformation().set_bond_angle( atomlist[ii-1].first, atomlist[ii].first, atomlist[ii+1].first, radians(b_ang[ii]) );
	}
	for ( core::Size ii=2, iimax=(t_ang.size()-2); ii<=iimax; ++ii ) { //Set torsion angles
		//TR << "Setting at1=" << atomlist[ii-1].first.atomno() << " rsd1=" << atomlist[ii-1].first.rsd() << " at2=" << atomlist[ii].first.atomno() << " rsd2=" << atomlist[ii].first.rsd() << " at3=" << atomlist[ii+1].first.atomno() << " rsd3=" << atomlist[ii+1].first.rsd() << " at4=" << atomlist[ii+2].first.atomno() << " rsd4=" << atomlist[ii+2].first.rsd() << " angle=" << t_ang[ii] << std::endl ; //DELETE ME
		pose.conformation().set_torsion_angle(atomlist[ii-1].first, atomlist[ii].first, atomlist[ii+1].first, atomlist[ii+2].first, radians(t_ang[ii]));
	}

	pose.update_residue_neighbors();

	//Rebuild the H and O atoms on peptide bonds and other atoms that are dependent on the connection to another atom.
	for ( core::Size ir=2, irmax=pose.size()-1; ir<=irmax; ++ir ) {
		core::Size const nresconn = pose.residue(ir).n_possible_residue_connections();
		if ( nresconn>0 ) {
			for ( core::Size ic=1; ic<=nresconn; ++ic ) {
				if ( !pose.residue(ir).connection_incomplete(ic) ) {
					core::Size const conn_at_index = pose.residue(ir).residue_connect_atom_index(ic); //The index of the connection atom
					for ( core::Size ia=1, iamax=pose.residue(ir).natoms(); ia<=iamax; ++ia ) {
						if ( pose.residue(ir).icoor(ia).stub_atom(1).atomno()==conn_at_index /*||
								pose.residue(ir).icoor(ia).stub_atom(2).atomno()==conn_at_index ||
								pose.residue(ir).icoor(ia).stub_atom(3).atomno()==conn_at_index */ ) {
							//TR << "Rebuilding rsd " << ir << " atom " << ia << " (" << pose.residue(ir).atom_name(ia) << ")" << std::endl; TR.flush();
							pose.conformation().set_xyz( AtomID( ia, ir ), pose.residue(ir).icoor(ia).build(pose.residue(ir), pose.conformation()) );
							//pose.conformation().set_xyz( AtomID( ia, ir ), pose.residue(ir).icoor(ia).build(pose.residue(ir) ) );
						}
					}
					//pose.conformation().rebuild_residue_connection_dependent_atoms(ir, ic);
				}
			}
		}
		//pose.conformation().rebuild_polymer_bond_dependent_atoms(ir);
	}

	pose.update_residue_neighbors();

	//pose.dump_pdb("after.pdb"); //DELETE ME -- for debugging only.

	if ( TR.visible() ) TR.flush();

	return;
} //set_loop_pose

/// @brief Set the loop pose conformation based on a set of results from kinematic closure.
/// @details  This version ONLY sets mainchain torsions, and does no rebuilding of mainchain O or H atoms.
/// @param[in,out] pose -- A pose consisting of the loop to be closed only.
/// @param[in] atomlist -- A list of AtomIDs and corresponding xyz coordinates (though we don't use the latter) of the chain of atoms closed by KIC.  Note that the residue indices refer to the loop pose, not the original pose.
/// @param[in] torsions -- The torsion angles values to set.
void set_loop_pose (
	core::pose::Pose &pose,
	utility::vector1 < std::pair < core::id::AtomID, numeric::xyzVector<core::Real> > > const &atomlist, //I want this to be const access
	utility::vector1 < core::Real > const &t_ang
) {
	using namespace numeric::conversions;
	using namespace core::id;

	for ( core::Size ii=2, iimax=(t_ang.size()-2); ii<=iimax; ++ii ) { //Set torsion angles
		//TR << "Setting at1=" << atomlist[ii-1].first.atomno() << " rsd1=" << atomlist[ii-1].first.rsd() << " at2=" << atomlist[ii].first.atomno() << " rsd2=" << atomlist[ii].first.rsd() << " at3=" << atomlist[ii+1].first.atomno() << " rsd3=" << atomlist[ii+1].first.rsd() << " at4=" << atomlist[ii+2].first.atomno() << " rsd4=" << atomlist[ii+2].first.rsd() << " angle=" << t_ang[ii] << std::endl ; //DELETE ME
		pose.conformation().set_torsion_angle(atomlist[ii-1].first, atomlist[ii].first, atomlist[ii+1].first, atomlist[ii+2].first, radians(t_ang[ii]));
	}
	pose.update_residue_neighbors();

	return;
} //set_loop_pose

/// @brief Copy the atom positions of the residues in the loop pose to the original pose.
///
void copy_loop_pose_to_original (
	core::pose::Pose &original_pose,
	core::pose::Pose const &loop_pose,
	utility::vector1 < std::pair < core::Size, core::Size > > const &residue_map,
	utility::vector1 < std::pair < core::Size, core::Size > > const &tail_residue_map
) {
	core::Size nres = residue_map.size();
	if ( nres==0 ) return;

	for ( core::Size repeats=1; repeats<=2; ++repeats ) {
		for ( core::Size ir=1; ir<=nres; ++ir ) { //Loop through all residues in the residue map (repeats==1), then in the tail_residue_map (repeats==2)
			core::Size const loop_index = ((repeats==1)?residue_map[ir].first:tail_residue_map[ir].first);
			core::Size const original_pose_index = ((repeats==1)?residue_map[ir].second:tail_residue_map[ir].second);
			core::Size const natom = loop_pose.residue( loop_index ).natoms();
			for ( core::Size ia=1; ia<=natom; ++ia ) { //Loop through all atoms in the current residue.
				std::string const curatom = loop_pose.residue(loop_index).atom_name(ia);
				assert( original_pose.residue(original_pose_index).has( curatom ) ); //It shouldn't be possible for the original pose to lack this atom, but let's check in any case
				original_pose.set_xyz(core::id::NamedAtomID( curatom, original_pose_index ), loop_pose.residue(loop_index).xyz(ia)); //Copy the xyz coordinates of each atom.
			}
		}
		nres=tail_residue_map.size();
		if ( nres==0 ) break;
	}

	original_pose.update_residue_neighbors();

	if ( TR.Debug.visible() ) { TR.Debug << "Copied coordinates of the loop pose to the original pose." << std::endl;  TR.Debug.flush(); }

	return;
}

/// @brief Update the polymer bond-dependent atom positions for every residue in the
/// residue map.
void
correct_polymer_dependent_atoms_in_pose_segment (
	core::pose::Pose &original_pose,
	utility::vector1 < std::pair < core::Size, core::Size > > const &residue_map
) {
	for(core::Size i=1, imax=residue_map.size(); i<=imax; ++i ) {
		original_pose.conformation().rebuild_polymer_bond_dependent_atoms_this_residue_only( residue_map[i].second );
	}
}

/// @brief Sets phi for an L-alpha amino acid, even if the nitrogen has a nonstandard connection.
///
void general_set_phi (
	core::pose::Pose &pose,
	core::Size const residue_index,
	core::Real const &phi_value
) {
	using namespace core::id;
	using namespace numeric::conversions;

	core::Size const N_connection_index = 1;

	if ( !pose.residue(residue_index).type().is_alpha_aa() ) {
		if ( TR.Warning.visible() ) { TR.Warning << "Residue " << residue_index << " was passed to general_set_phi, but this isn't an alpha-amino acid.  Skipping." << std::endl;  TR.Warning.flush(); }
		return;
	}
	if ( pose.residue(residue_index).is_lower_terminus() ) {
		if ( TR.Debug.visible() ) { TR.Debug << "Residue " << residue_index << " was passed to general_set_phi, but this is a lower terminus.  Skipping." << std::endl;  TR.Debug.flush(); }
		return;
	}
	if ( pose.residue(residue_index).connection_incomplete(N_connection_index) ) {
		if ( TR.Debug.visible() ) { TR.Debug << "Residue " << residue_index << " was passed to general_set_phi, but this isn't connected to anything at its N-terminus.  Skipping." << std::endl;  TR.Debug.flush(); }
		return; //If nothing is attached here, we can't set this phi.
	}

	core::Size const other_res = pose.residue(residue_index).connected_residue_at_resconn(N_connection_index); //The residue connected to the nitrogen.
	//This is convoluted, but we want the OTHER residue's atom index corresponding to the OTHER residue's connection to THIS residue:
	core::Size const other_atmno = pose.residue(other_res).residue_connect_atom_index( pose.residue( residue_index ).residue_connection_conn_id(N_connection_index) );

	pose.conformation().set_torsion_angle(
		AtomID(other_atmno, other_res),
		AtomID(1, residue_index),
		AtomID(2, residue_index),
		AtomID(3, residue_index),
		radians(phi_value)
	);

	pose.update_residue_neighbors();

	return;
}

/// @brief Sets psi for an L-alpha amino acid, even if the carbonyl carbon has a nonstandard connection.
///
void general_set_psi (
	core::pose::Pose &pose,
	core::Size const residue_index,
	core::Real const &psi_value
) {
	using namespace core::id;
	using namespace numeric::conversions;

	core::Size C_connection_index = 2;

	if ( !pose.residue(residue_index).type().is_alpha_aa() ) {
		if ( TR.Warning.visible() ) { TR.Warning << "Residue " << residue_index << " was passed to general_set_psi, but this isn't an alpha-amino acid.  Skipping." << std::endl;  TR.Warning.flush(); }
		return;
	}
	if ( pose.residue(residue_index).is_lower_terminus() ) {
		--C_connection_index;
	}
	if ( pose.residue(residue_index).is_upper_terminus() ) {
		if ( TR.Debug.visible() ) { TR.Debug << "Residue " << residue_index << " was passed to general_set_psi, but this is an upper terminus.  Skipping." << std::endl;  TR.Debug.flush(); }
		return;
	}
	if ( pose.residue(residue_index).connection_incomplete(C_connection_index) ) {
		if ( TR.Debug.visible() ) { TR.Debug << "Residue " << residue_index << " was passed to general_set_psi, but this isn't connected to anything at its N-terminus.  Skipping." << std::endl;  TR.Debug.flush(); }
		return; //If nothing is attached here, we can't set this phi.
	}

	core::Size const other_res = pose.residue(residue_index).connected_residue_at_resconn(C_connection_index); //The residue connected to the carbonyl carbon.
	//This is convoluted, but we want the OTHER residue's atom index corresponding to the OTHER residue's connection to THIS residue:
	core::Size const other_atmno = pose.residue(other_res).residue_connect_atom_index( pose.residue( residue_index ).residue_connection_conn_id(C_connection_index) );

	pose.conformation().set_torsion_angle(
		AtomID(1, residue_index),
		AtomID(2, residue_index),
		AtomID(3, residue_index),
		AtomID(other_atmno, other_res),
		radians(psi_value)
	);

	pose.update_residue_neighbors();

	return;
}

/// @brief Given a residue_map vector of pairs, where each pair is < residue_index_in_perturbedloop_pose, residue_index_in_original_pose >,
/// and a residue index in the perturbed loop, return the corresponding residue index in the original pose.
core::Size get_original_pose_rsd ( core::Size const perturbedloop_rsd, utility::vector1 < std::pair < core::Size, core::Size > > const &residue_map )
{
	for ( core::Size i=1, imax=residue_map.size(); i<=imax; ++i ) {
		if ( residue_map[i].first == perturbedloop_rsd ) return residue_map[i].second;
	}

	utility_exit_with_message("Error in generalized_kinematic_closure::get_original_pose_rsd.");
	return 0;
}

} //namespace generalized_kinematic_closure
} //namespace protocols
