// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/backbone_moves/local_backbone_mover/FreePeptide.cc
/// @brief FreePeptide represents a free peptide
/// @author xingjiepan (xingjiepan@gmail.com)

#include <protocols/backbone_moves/local_backbone_mover/FreePeptide.hh>
#include <protocols/backbone_moves/local_backbone_mover/util.hh>
#include <basic/Tracer.hh>

// Numeric headers
#include <numeric/conversions.hh>
#include <numeric/alignment/QCP_Kernel.hh>

// Core headers
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>

static basic::Tracer TR( "protocols.backbone_moves.local_backbone_mover.FreePeptide" );


namespace protocols {
namespace backbone_moves {
namespace local_backbone_mover {

Residue::Residue(Size seqpos, core::pose::Pose const& pose){
	using core::id::AtomID;

	seqpos_ = seqpos;

	// Get the relevant atom ids

	AtomID n_id(1, seqpos_);
	AtomID ca_id(2, seqpos_);
	AtomID c_id(3, seqpos_);
	AtomID c_id_prev(3, seqpos_ - 1);
	AtomID n_id_next(1, seqpos_ + 1);

	// Get the reference to the conformation

	core::conformation::Conformation const& conformation(pose.conformation());

	// Set bond_lengths

	n_ca_bond_ = conformation.bond_length(n_id, ca_id);
	ca_c_bond_ = conformation.bond_length(ca_id, c_id);
	c_n_bond_ = conformation.bond_length(c_id, n_id_next);

	// Set bond angles

	n_ca_c_angle_ = conformation.bond_angle(n_id, ca_id, c_id);
	ca_c_n_angle_ = conformation.bond_angle(ca_id, c_id, n_id_next);
	c_n_ca_angle_ = conformation.bond_angle(c_id_prev, n_id, ca_id);

	// Set torsions

	phi_ = pose.phi(seqpos_);
	psi_ = pose.psi(seqpos_);
	omega_ = pose.omega(seqpos_);

	// Get the xyz coordinates

	n_xyz_ = pose.xyz(n_id);
	ca_xyz_ = pose.xyz(ca_id);
	c_xyz_ = pose.xyz(c_id);
}

void
Residue::apply_to_pose(core::pose::Pose &pose){
	using core::id::AtomID;

	// Get the relevant atom ids

	AtomID n_id(1, seqpos_);
	AtomID ca_id(2, seqpos_);
	AtomID c_id(3, seqpos_);
	AtomID c_id_prev(3, seqpos_ - 1);
	AtomID n_id_next(1, seqpos_ + 1);

	// Get the reference to the conformation

	core::conformation::Conformation &conformation(pose.conformation());

	// Set bond_lengths

	conformation.set_bond_length(n_id, ca_id, n_ca_bond_);
	conformation.set_bond_length(ca_id, c_id, ca_c_bond_);
	conformation.set_bond_length(c_id, n_id_next, c_n_bond_);

	// Set bond angles

	conformation.set_bond_angle(n_id, ca_id, c_id, n_ca_c_angle_);
	conformation.set_bond_angle(ca_id, c_id, n_id_next, ca_c_n_angle_);
	conformation.set_bond_angle(c_id_prev, n_id, ca_id, c_n_ca_angle_);

	// Set torsions

	pose.set_phi(seqpos_, phi_);
	pose.set_psi(seqpos_, psi_);
	pose.set_omega(seqpos_, omega_);
}

FreePeptide::FreePeptide(core::pose::Pose const& pose, Size pivot1, Size pivot2):
	utility::pointer::ReferenceCount(),
	pivot1_(pivot1),
	pivot2_(pivot2)
{
	using core::id::AtomID;

	// The free peptide must have at least 2 residues

	runtime_assert(pivot2 > pivot1 + 2);

	// Initialize the residues

	for ( Size i = pivot1 - 2; i <= pivot2 + 2; ++i ) {
		residues_.push_back( Residue(i, pose) );
	}

	// Initialize old stubs

	old_stubs_.push_back( pose.xyz( AtomID(2, pivot1 + 1) ) ); // fisrt CA
	old_stubs_.push_back( pose.xyz( AtomID(3, pivot1 + 1) ) ); // fisrt C
	old_stubs_.push_back( pose.xyz( AtomID(1, pivot2 - 1) ) ); // fisrt N
	old_stubs_.push_back( pose.xyz( AtomID(2, pivot2 - 1) ) ); // fisrt CA

	old_stub_centor_of_mass_.zero();
	for ( Size i=1; i<=4; ++i ) {
		old_stub_centor_of_mass_ += old_stubs_[i];
	}
	old_stub_centor_of_mass_ *= 0.25;


	// Initialize reference coordinates

	reference_atom_ids_.push_back( AtomID(3, pivot1) ); // C
	reference_atom_ids_.push_back( AtomID(1, pivot1 + 1) ); // N
	reference_atom_ids_.push_back( AtomID(2, pivot1 + 1) ); // CA

	for ( AtomID &aid : reference_atom_ids_ ) {
		reference_coordinates_.push_back( pose.xyz(aid) );
	}

}

FreePeptide::~FreePeptide()= default;

FreePeptideOP
FreePeptide::clone() const {
	return FreePeptideOP( new FreePeptide( *this ) );
}

void
FreePeptide::apply_to_pose(core::pose::Pose &pose){
	for ( Size i = pivot1_ - 1; i <= pivot2_ + 1; ++i ) {
		residues_[res_id(i)].apply_to_pose(pose);
	}
}

void
FreePeptide::update_xyz_coords(){

	using numeric::conversions::radians;

	// Assign the first 2 positions from reference coordinates

	residues_[4].n_xyz(reference_coordinates_[2]);
	residues_[4].ca_xyz(reference_coordinates_[3]);

	// Calculate the rest of xyz coordinates from internal coordinates

	residues_[4].c_xyz( xyz_from_internal_coords(residues_[4].ca_xyz(),
		residues_[4].n_xyz(), reference_coordinates_[1],
		radians(residues_[4].phi()), residues_[4].n_ca_c_angle(), residues_[4].ca_c_bond() ) );

	for ( Size i=5; i<res_id(pivot2_); ++i ) {
		residues_[i].n_xyz( xyz_from_internal_coords(residues_[i - 1].c_xyz(),
			residues_[i - 1].ca_xyz(), residues_[i - 1].n_xyz(),
			radians(residues_[i - 1].psi()), residues_[i - 1].ca_c_n_angle(), residues_[i - 1].c_n_bond() ) );

		residues_[i].ca_xyz( xyz_from_internal_coords(residues_[i].n_xyz(),
			residues_[i - 1].c_xyz(), residues_[i - 1].ca_xyz(),
			radians(residues_[i - 1].omega()), residues_[i].c_n_ca_angle(), residues_[i].n_ca_bond() ) );

		residues_[i].c_xyz( xyz_from_internal_coords(residues_[i].ca_xyz(),
			residues_[i].n_xyz(), residues_[i - 1].c_xyz(),
			radians(residues_[i].phi()), residues_[i].n_ca_c_angle(), residues_[i].ca_c_bond() ) );
	}

	xyz_updated_ = true;
}

void
FreePeptide::translate(xyzVector <Real> v){
	for ( xyzVector <Real> &coord : reference_coordinates_ ) {
		coord += v;
	}
	xyz_updated_ = false;
}

void
FreePeptide::rotate(xyzMatrix <Real> M){
	xyzVector <Real> center = reference_coordinates_[3];

	for ( Size i=1; i<=2; ++i ) {
		xyzVector <Real> delta = reference_coordinates_[i] - center;
		reference_coordinates_[i] = center + M * delta;
	}
	xyz_updated_ = false;
}

void
FreePeptide::align(){
	numeric::alignment::QCP_Kernel<Real> kernel;
	Real old_stubs[12];
	Real current_stubs[12];
	vector1 < xyzVector<Real> > current_stubs_vec;
	Real rot_matrix[9];

	// Initialize the coordinates for alignment

	current_stubs_vec.push_back( ca_xyz(pivot1() + 1) );
	current_stubs_vec.push_back( c_xyz(pivot1() + 1) );
	current_stubs_vec.push_back( n_xyz(pivot2() - 1) );
	current_stubs_vec.push_back( ca_xyz(pivot2() - 1) );

	for ( Size i=1; i<=4; ++i ) {
		for ( Size j=1; j<=3; ++j ) {
			old_stubs[3 * (i-1) + j - 1] = old_stubs_[i](j);
			current_stubs[3 * (i-1) + j - 1] = current_stubs_vec[i](j);
		}
	}

	// Align stubs

	kernel.calc_coordinate_rmsd(old_stubs, current_stubs, 4, rot_matrix);

	// Translate the free peptide to overlap the centers of mass

	xyzVector <Real> current_stub_center_of_mass(0, 0, 0);
	for ( Size i=1; i<=4; ++i ) {
		current_stub_center_of_mass += current_stubs_vec[i];
	}
	current_stub_center_of_mass *= 0.25;

	for ( Size i=1; i<=3; ++i ) {
		reference_coordinates_[i] += old_stub_centor_of_mass_ - current_stub_center_of_mass;
	}

	// Rotate the free peptide

	xyzMatrix <Real> rot_M( xyzMatrix <Real>::rows( rot_matrix ) );

	for ( Size i=1; i<=3; ++i ) {
		reference_coordinates_[i] = old_stub_centor_of_mass_ \
			+ rot_M * (reference_coordinates_[i] - old_stub_centor_of_mass_);
	}

	xyz_updated_ = false;
}

} //protocols
} //backbone_moves
} //local_backbone_mover






