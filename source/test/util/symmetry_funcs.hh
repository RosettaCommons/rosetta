// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/util/symmetry_funcs.hh
/// @brief  Functions to help test symmetry
/// @author Vikram K. Mulligan (vmullig@uw.edu)
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_util_symmetry_funcs_HH
#define INCLUDED_util_symmetry_funcs_HH

// Test headers

// Project headers
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/util.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

/// @brief Given a pose, remove all disulfides in it.
///
inline
void remove_disulfides( core::pose::PoseOP pose ) {
	utility::vector1< std::pair < core::Size, core::Size > > disulfides;
	core::conformation::disulfide_bonds( pose->conformation(), disulfides );

	for ( core::Size i=1, imax=disulfides.size(); i<=imax; ++i ) {
		core::conformation::break_disulfide( pose->conformation(), disulfides[i].first, disulfides[i].second );
	}
}

/// @brief Given a pose, add disulfides between the first cysteine and the next.
///
inline
void form_disulfides( core::pose::PoseOP pose) {
	bool breaknow(false);
	for ( core::Size ir=1, nres=pose->size(); ir<nres; ++ir ) { //Loop through all residues except the last
		if ( pose->residue(ir).name3() == "CYS" || pose->residue(ir).name3() == "DCS" ) {
			for ( core::Size jr=ir+1; jr<=nres; ++jr ) { //Loop through rest of residues
				if ( pose->residue(jr).name3() == "CYS" || pose->residue(jr).name3() == "DCS" ) {
					core::conformation::form_disulfide( pose->conformation(), ir, jr, true, false );
					breaknow=true;
					break;
				}
			}
		}
		if ( breaknow ) break;
	}
}

/// @brief Given a residue type, get its mirror-image type.
///
inline
core::chemical::ResidueType const & get_mirror_type( core::chemical::ResidueType const &master_type, core::chemical::ResidueTypeSet const & residue_type_set ) {
	if ( !master_type.is_l_aa() && !master_type.is_d_aa() ) return master_type;
	if ( master_type.is_l_aa() ) {
		return residue_type_set.name_map( "D"+master_type.name() );
	}
	return residue_type_set.name_map( master_type.name().substr(1) );
}

/// @brief Given a pose, flip the L-residues to D-residues.
///
inline
void flip_residues( core::pose::Pose &pose ) {
	// Create the new residue and replace it
	for ( core::Size ir=1, irmax=pose.size(); ir<=irmax; ++ir ) {
		core::chemical::ResidueType const & mirror_type( get_mirror_type( pose.residue(ir).type(), *pose.residue_type_set_for_pose( pose.residue(ir).type().mode() ) ) );
		core::conformation::ResidueOP new_res = core::conformation::ResidueFactory::create_residue( mirror_type, pose.residue( ir ), pose.conformation());
		core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms( pose.residue( ir ), *new_res, pose.conformation(), true );
		pose.replace_residue( ir, *new_res, false );
	}
}

/// @brief Given a pose, construct its mirror image.
///
inline
void mirror_pose( core::pose::Pose const &master, core::pose::Pose &mirror) {

	for ( core::Size ir=1, irmax=mirror.size(); ir<=irmax; ++ir ) {
		for ( core::Size ia=1, iamax=mirror.residue(ir).type().natoms(); ia<=iamax; ++ia ) {
			core::id::AtomID const curatom(ia, ir);
			numeric::xyzVector<core::Real> xyztemp( master.xyz( curatom ) );
			xyztemp.z( -1.0*xyztemp.z() );
			mirror.set_xyz( curatom, xyztemp );
		}
	}

	mirror.update_residue_neighbors();

	return;
}

/// @brief Given a pose, construct its mirror image.
///
inline
core::pose::PoseOP
mirror_pose( core::pose::PoseCOP master) {
	core::pose::PoseOP output_pose( master->clone() );
	flip_residues(*output_pose);

	for ( core::Size ir=1, irmax=output_pose->size(); ir<=irmax; ++ir ) {
		for ( core::Size ia=1, iamax=output_pose->residue(ir).type().natoms(); ia<=iamax; ++ia ) {
			core::id::AtomID const curatom(ia, ir);
			numeric::xyzVector<core::Real> xyztemp( master->xyz( curatom ) );
			xyztemp.z( -1.0*xyztemp.z() );
			output_pose->set_xyz( curatom, xyztemp );
		}
	}

	return output_pose;
}

/// @brief Given a pose, construct its mirror image.
///
inline
core::pose::PoseOP
mirror_pose_with_disulfides( core::pose::PoseCOP master) {
	core::pose::PoseOP refpose( master->clone() );
	remove_disulfides(refpose);
	core::pose::PoseOP output_pose( refpose->clone() );
	flip_residues(*output_pose);

	for ( core::Size ir=1, irmax=output_pose->size(); ir<=irmax; ++ir ) {
		for ( core::Size ia=1, iamax=output_pose->residue(ir).type().natoms(); ia<=iamax; ++ia ) {
			core::id::AtomID const curatom(ia, ir);
			numeric::xyzVector<core::Real> xyztemp( refpose->xyz( curatom ) );
			xyztemp.z( -1.0*xyztemp.z() );
			output_pose->set_xyz( curatom, xyztemp );
		}
	}

	form_disulfides(output_pose);

	return output_pose;
}

/// @brief Given an input cyclic pose, circularly permute by 1 residue.
inline
core::pose::PoseOP permute(
	core::pose::PoseOP input_pose
) {
	core::pose::PoseOP new_pose( new core::pose::Pose );

	for ( core::Size i=2, imax=input_pose->size(); i<=imax; ++i ) {
		core::conformation::ResidueOP new_rsd( input_pose->residue(i).clone() );
		if ( i == 2 ) {
			new_pose->append_residue_by_jump(*new_rsd, 1);
		} else {
			new_pose->append_residue_by_bond(*new_rsd, false);
		}
	}
	new_pose->append_residue_by_bond( *(input_pose->residue(1).clone()), false);
	new_pose->conformation().declare_chemical_bond(1, "N", 9, "C");
	return new_pose;
}

/// @brief Given an input cyclic pose, circularly permute by 1 residue.
inline
core::pose::PoseOP permute_with_disulfides(
	core::pose::PoseOP input_pose
) {
	core::pose::PoseOP ref_pose( input_pose->clone() );
	remove_disulfides(ref_pose);
	core::pose::PoseOP new_pose( new core::pose::Pose );

	for ( core::Size i=2, imax=ref_pose->size(); i<=imax; ++i ) {
		core::conformation::ResidueOP new_rsd( ref_pose->residue(i).clone() );
		if ( i == 2 ) {
			new_pose->append_residue_by_jump(*new_rsd, 1);
		} else {
			new_pose->append_residue_by_bond(*new_rsd, false);
		}
	}
	new_pose->append_residue_by_bond( *(ref_pose->residue(1).clone()), false);
	new_pose->conformation().declare_chemical_bond(1, "N", new_pose->size(), "C");
	form_disulfides(new_pose);
	return new_pose;
}

#endif
