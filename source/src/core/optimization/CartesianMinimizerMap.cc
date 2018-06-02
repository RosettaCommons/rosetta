// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   CartesianMinimizerMap
/// @brief  Kinematics
/// @author Frank Dimaio
/// @author Jared Adolf-Bryfogle (per-atom min control)

// Unit headers
#include <core/optimization/CartesianMinimizerMap.hh>


// Project headers
#include <core/pose/util.hh>
#include <core/pose/Pose.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/PDBInfo.hh> //JAB - temp

// Numeric headers

// C++ headers
#include <cstdlib>

#include <core/id/DOF_ID_Map.hh>
#include <core/id/AtomID_Map.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/pose/util.tmpl.hh>

// Options - TEMP
#include <core/select/util.hh>

namespace core {
namespace optimization {


/////////////////////////////////////////////////////////////////////////////
CartesianMinimizerMap::~CartesianMinimizerMap() {
	moving_atoms_.clear();
	moving_dofids_.clear();
	moving_torsionids_.clear();
}


/////////////////////////////////////////////////////////////////////////////
void
CartesianMinimizerMap::reset( pose::Pose const & pose ) {
	moving_atoms_.clear();
	moving_dofids_.clear();
	moving_torsionids_.clear();

	atom_derivatives_.resize( pose.size() );
	for ( Size ii = 1; ii <= pose.size(); ++ii ) {
		atom_derivatives_[ ii ].resize( pose.residue( ii ).natoms() );
	}
}

/////////////////////////////////////////////////////////////////////////////
void
CartesianMinimizerMap::copy_dofs_from_pose(
	pose::Pose const & pose,
	Multivec & dofs
) const {
	runtime_assert(dofs.size() == 3*moving_atoms_.size());
	Size natoms=moving_atoms_.size();
	utility::vector1< numeric::xyzVector<core::Real> > dofs_v(natoms);
	pose.batch_get_xyz( moving_atoms_, dofs_v );

	for ( Size i=1; i<=natoms; ++i ) {
		dofs[i*3-2] = dofs_v[i].x();
		dofs[i*3-1] = dofs_v[i].y();
		dofs[i*3]   = dofs_v[i].z();
	}
}


/////////////////////////////////////////////////////////////////////////////
void
CartesianMinimizerMap::copy_dofs_to_pose(
	pose::Pose & pose,
	Multivec const & dofs
) const
{
	runtime_assert(dofs.size() == 3*moving_atoms_.size());

	Size natoms=moving_atoms_.size();
	utility::vector1< numeric::xyzVector<core::Real> > dofs_v(natoms);
	for ( Size i=1; i<=natoms; ++i ) {
		dofs_v[i].x() = dofs[i*3-2];
		dofs_v[i].y() = dofs[i*3-1];
		dofs_v[i].z() = dofs[i*3];
	}
	pose.batch_set_xyz( moving_atoms_, dofs_v );
}

////////////////////////////////////////////////////////////////////////////

void
CartesianMinimizerMap::add_torsion(
	DOF_ID const & dof_id,
	DOF_ID const &
) {
	// only care about torsional DOFs
	if ( dof_id.type() == id::PHI ) {
		moving_dofids_.push_back( dof_id );
	}
}

/////////////////////////////////////////////////////////////////////////////
void
CartesianMinimizerMap::add_atom(
	AtomID const &,
	DOF_ID const &
)
{
	; // do nothing
}

/////////////////////////////////////////////////////////////////////////////
void
CartesianMinimizerMap::zero_stored_derivs()
{
	for ( Size ii = 1; ii <= atom_derivatives_.size(); ++ii ) {
		for ( Size jj = 1; jj <= atom_derivatives_[ ii ].size(); ++jj ) {
			atom_derivatives_[ ii ][ jj ].f1() = 0.0;
			atom_derivatives_[ ii ][ jj ].f2() = 0.0;
		}
	}
}


/////////////////////////////////////////////////////////////////////////////
// private
void
CartesianMinimizerMap::assign_rosetta_torsions_and_trim( pose::Pose const & pose )
{
	utility::vector1<id::DOF_ID> new_moving_dofids;

	// mapping from AtomTree DOF ID's to bb/chi torsion angle ids
	id::DOF_ID_Map< id::TorsionID > dof_map ( id::TorsionID::BOGUS_TORSION_ID() );
	pose::setup_dof_to_torsion_map( pose, dof_map );

	Size ndofs = moving_dofids_.size();
	for ( Size i=1; i<=ndofs; ++i ) {
		DOF_ID const& dof_id( moving_dofids_[i] );
		id::TorsionID const & tor_id( dof_map[ dof_id ] );

		if ( tor_id.valid() ) {  // we don't care otherwise
			// check atoms
			id::AtomID id1,id2,id3,id4;
			pose.conformation().get_torsion_angle_atom_ids( tor_id, id1,id2,id3,id4 );
			if ( atom_indices_.has(id1)
					|| atom_indices_.has(id2)
					|| atom_indices_.has(id3)
					|| atom_indices_.has(id4) ) {
				new_moving_dofids.push_back(dof_id);
				moving_torsionids_.push_back(tor_id);
			}
		}
	}

	moving_dofids_ = new_moving_dofids;
}


/////////////////////////////////////////////////////////////////////////////

void
CartesianMinimizerMap::setup(
	pose::Pose & pose,
	kinematics::MoveMap const & mm
)
{

	// this clears moving_atoms_
	reset( pose );

	/////////////////////
	// convert the allow_bb,allow_chi,allow_jump information
	//   in the MoveMap into a simple boolean mask over movable atoms
	// because the movemap is designed with torsion-space refinement in mind,
	//   interpret the meaning as best as possible
	// implicit DOFs and jumps don't make a lot of sense so ignore them
	Size const n_res( pose.size() );
	core::pose::initialize_atomid_map( atom_indices_, pose );

	// special case for symmetry
	core::conformation::symmetry::SymmetryInfoOP symm_info;
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		auto & SymmConf (
			dynamic_cast<core::conformation::symmetry::SymmetricConformation &> ( pose.conformation()) );
		symm_info = SymmConf.Symmetry_Info();
	}

	/////////////////////
	// setup the domain_map which indicates what rsd pairs are fixed/moving
	id::AtomID_Mask moving_dof, moving_xyz;
	core::pose::initialize_atomid_map( moving_xyz, pose, false );
	core::pose::initialize_atomid_map( moving_dof, pose, false );

	id::AtomID_Map< bool > atom_mask;
	core::pose::initialize_atomid_map( atom_mask, pose, false );

	for ( Size i = 1; i <= n_res; ++i ) {

		conformation::Residue const & rsd( pose.residue(i) );

		bool const bb_move( mm.get_bb(i) );
		bool const chi_move( mm.get_chi(i) );

		//fpd  do not let aa_vrt move
		if ( (rsd.aa() == chemical::aa_vrt )|| rsd.is_virtual_residue() ) continue;

		// cartesian logic ...
		//    if (chi_move && !bb_move) sc atoms only
		//    if (bb_move) all atoms
		//
		//JAB - specific torsions do not make sense as multiple atoms can potentially define a specific torsion,
		//  So enabling a torsion does not nessessarily enable the atoms that you are expecting.
		//  Instead, I have implemented the ability for the MoveMap to set specific atoms.


		if ( chi_move || bb_move ) {
			Size start1 = rsd.first_sidechain_atom();
			Size stop1  = rsd.nheavyatoms();
			Size start2 = rsd.first_sidechain_hydrogen();
			Size stop2  = rsd.natoms();

			//std::cout << "original enable: rsd=" << i << " : " << pose.pdb_info()->pose2pdb(i) << std::endl;
			//std::cout << start1 << " " << stop1 << " " << start2 <<" "<< stop2 << std::endl;
			for ( Size j=1; j<=stop2; ++j ) {
				if ( !bb_move && j<start1 ) continue;
				if ( !bb_move && j>stop1 && j<start2 ) continue;

				id::AtomID atm_ij( j,i );

				atom_mask.set( atm_ij, true);
			}
		}
	} // i=1,n_res

	//JAB - go through specific atom OVERRIDEs to the default behavior.
	std::map< AtomID, bool > const & atoms = mm.get_atoms();
	for ( auto atom_pair : atoms ) {
		atom_mask.set( atom_pair.first, atom_pair.second);
	}

	//Now go through all atoms and add to the list of movable atoms.
	for ( core::Size i= 1; i <= pose.size(); ++i ) {
		bool const is_independent = (!symm_info || symm_info->bb_is_independent(i));
		for ( core::Size x = 1; x <= atom_mask.n_atom( i ); ++x ) {
			AtomID atom = AtomID( x, i ); //Think about adding a Size,Size accessor to AtomID to increase speed here.
			if ( atom_mask.get(atom) ) {
				moving_xyz[ atom ] = true;

				if ( is_independent ) {
					moving_atoms_.push_back( atom );
					atom_indices_[ atom ] = moving_atoms_.size();
				}
			}
		}
	}

	//std::cout << "Total moving atoms: " << moving_atoms_.size() << std::endl;

	//JAB - enable this for debugging to get all atoms moving during cartesian min.
	//std::cout << select::get_pymol_selection_for_atoms(pose, moving_atoms_, "moving") << std::endl;

	/////////////////////
	// get a list of torsional DOFs which are implicitly moved by these xyzs
	kinematics::MoveMap move_map_torsional;
	move_map_torsional.set_bb(true);
	move_map_torsional.set_chi(true);
	move_map_torsional.set_nu(true);
	move_map_torsional.set_branches(true);
	id::DOF_ID_Mask dof_mask(false);
	pose::setup_dof_mask_from_move_map( move_map_torsional, pose, dof_mask );

	// fill the torsion list
	DOF_ID tmp( id::DOF_ID::BOGUS_DOF_ID() );
	pose.atom_tree().root()->setup_min_map( tmp, dof_mask, *this );

	// now trim this list ensuring that at least one of the
	//   four atoms that define the torsion are in our movable set
	assign_rosetta_torsions_and_trim( pose );

	domain_map_.dimension( pose.size() );
	pose.conformation().atom_tree().update_domain_map( domain_map_, moving_dof, moving_xyz );
}


} // namespace kinematics
} // namespace core
