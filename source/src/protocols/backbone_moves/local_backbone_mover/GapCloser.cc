// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/backbone_moves/local_backbone_mover/GapCloser.cc
/// @brief GapCloser closes the gaps after moving the free peptide.
/// @author xingjiepan (xingjiepan@gmail.com)

#include <protocols/backbone_moves/local_backbone_mover/util.hh>
#include <protocols/backbone_moves/local_backbone_mover/FreePeptide.hh>
#include <protocols/backbone_moves/local_backbone_mover/GapCloser.hh>
#include <protocols/backbone_moves/local_backbone_mover/gap_solution_pickers/RandomGapSolutionPicker.hh>
#include <basic/Tracer.hh>

// Numeric headers
#include <numeric/kinematic_closure/bridgeObjects_nonredundant.hh>
#include <numeric/conversions.hh>

// Core headers
#include <core/kinematics/Edge.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <utility/fixedsizearray1.hh>


static THREAD_LOCAL basic::Tracer TR( "protocols.backbone_moves.local_backbone_mover.GapCloser" );


namespace protocols {
namespace backbone_moves {
namespace local_backbone_mover {

GapCloser::GapCloser():
	utility::pointer::ReferenceCount()
{
	set_solution_picker(gap_solution_pickers::GapSolutionPickerOP(new gap_solution_pickers::RandomGapSolutionPicker));
}

GapCloser::~GapCloser(){}

GapCloser::GapCloser( GapCloser const & ) {

}

GapCloserOP
GapCloser::clone() const {
	return GapCloserOP( new GapCloser( *this ) );
}

void
GapCloser::solve_gaps(FreePeptide &free_peptide){
	Size pivot1 = free_peptide.pivot1();
	Size pivot2 = free_peptide.pivot2();

	gap_solved_ = solve_a_gap(free_peptide, pivot1, pivot_torsions1_);

	if ( gap_solved_ ) {
		gap_solved_ = solve_a_gap(free_peptide, pivot2, pivot_torsions2_);
	}
}

void
GapCloser::apply_closure(core::pose::Pose &pose, FreePeptide &free_peptide){
	assert(gap_solved_);

	using core::kinematics::Edge;
	using core::kinematics::FoldTree;

	Size pivot1 = free_peptide.pivot1();
	Size pivot2 = free_peptide.pivot2();

	// Temporarily change the fold tree such that numerical errors
	// don't propogate

	FoldTree original_ft(pose.fold_tree());
	FoldTree tmp_ft;

	tmp_ft.add_edge(1, pivot2 + 2, Edge::PEPTIDE);
	tmp_ft.add_edge(pivot1 - 2, pivot2 + 3, 1);
	tmp_ft.add_edge(pivot2 + 3, pose.size(), Edge::PEPTIDE);

	pose.fold_tree(tmp_ft);

	// Apply the internal coordinates of the free peptide

	free_peptide.apply_to_pose(pose);

	// Pick solutions

	Size sol1, sol2;
	pick_solutions(sol1, sol2, pose, free_peptide);

	// Apply the torsions at gaps

	pose.set_phi(pivot1 - 1, pivot_torsions1_[sol1][1]);
	pose.set_psi(pivot1 - 1, pivot_torsions1_[sol1][2]);
	pose.set_phi(pivot1, pivot_torsions1_[sol1][3]);
	pose.set_psi(pivot1, pivot_torsions1_[sol1][4]);
	pose.set_phi(pivot1 + 1, pivot_torsions1_[sol1][5]);
	pose.set_psi(pivot1 + 1, pivot_torsions1_[sol1][6]);

	pose.set_phi(pivot2 - 1, pivot_torsions2_[sol2][1]);
	pose.set_psi(pivot2 - 1, pivot_torsions2_[sol2][2]);
	pose.set_phi(pivot2, pivot_torsions2_[sol2][3]);
	pose.set_psi(pivot2, pivot_torsions2_[sol2][4]);
	pose.set_phi(pivot2 + 1, pivot_torsions2_[sol2][5]);
	pose.set_psi(pivot2 + 1, pivot_torsions2_[sol2][6]);

	// Reset the fold tree

	pose.fold_tree(original_ft);
}

bool
GapCloser::solve_a_gap(FreePeptide &free_peptide, Size pivot, vector1<vector1<Real> > &pivot_torsions){

	using numeric::conversions::degrees;
	using utility::fixedsizearray1;
	
	vector1<fixedsizearray1 <Real,3> > stub1(3);
	vector1<fixedsizearray1 <Real,3> > stub2(3);
	vector1<Real> torsions_chain1(1);
	vector1<Real> torsions_chain2(1);
	vector1<Real> angles(7);
	vector1<Real> bonds(6);
	int nsol = 0;

	xyz_to_vec1(free_peptide.c_xyz(pivot - 2), stub1[1]);
	xyz_to_vec1(free_peptide.n_xyz(pivot - 1), stub1[2]);
	xyz_to_vec1(free_peptide.ca_xyz(pivot - 1), stub1[3]);

	xyz_to_vec1(free_peptide.ca_xyz(pivot + 1), stub2[1]);
	xyz_to_vec1(free_peptide.c_xyz(pivot + 1), stub2[2]);
	xyz_to_vec1(free_peptide.n_xyz(pivot + 2), stub2[3]);

	torsions_chain1[1] = free_peptide.omega(pivot - 1);
	torsions_chain2[1] = free_peptide.omega(pivot);

	angles[1] = degrees( free_peptide.n_ca_c_angle(pivot - 1) );
	angles[2] = degrees( free_peptide.ca_c_n_angle(pivot - 1) );
	angles[3] = degrees( free_peptide.c_n_ca_angle(pivot) );
	angles[4] = degrees( free_peptide.n_ca_c_angle(pivot) );
	angles[5] = degrees( free_peptide.ca_c_n_angle(pivot) );
	angles[6] = degrees( free_peptide.c_n_ca_angle(pivot + 1) );
	angles[7] = degrees( free_peptide.n_ca_c_angle(pivot + 1) );

	bonds[1] = free_peptide.ca_c_bond(pivot - 1);
	bonds[2] = free_peptide.c_n_bond(pivot - 1);
	bonds[3] = free_peptide.n_ca_bond(pivot);
	bonds[4] = free_peptide.ca_c_bond(pivot);
	bonds[5] = free_peptide.c_n_bond(pivot);
	bonds[6] = free_peptide.n_ca_bond(pivot - 1);

	numeric::kinematic_closure::bridgeObjects_nonredundant(stub1, stub2,
		torsions_chain1, torsions_chain2, angles, bonds, pivot_torsions, nsol);

	return nsol > 0;
}

void
GapCloser::pick_solutions(Size &index1, Size &index2, core::pose::Pose &pose, FreePeptide &free_peptide){
	index1 = solution_picker_->pick(pose, free_peptide, pivot_torsions1_, free_peptide.pivot1());
	index2 = solution_picker_->pick(pose, free_peptide, pivot_torsions2_, free_peptide.pivot2());
}

} //protocols
} //backbone_moves
} //local_backbone_mover






