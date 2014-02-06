// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   ClosureProblem.cc
/// @brief  Source file for ClosureProblem.
/// @author Kale Kundert (kale.kundert@ucsf.edu)

// Unit headers
#include <protocols/kinematic_closure/types.hh>
#include <protocols/kinematic_closure/internal.hh>
#include <protocols/kinematic_closure/ClosureProblem.hh>
#include <protocols/kinematic_closure/ClosureSolution.hh>
#include <protocols/kinematic_closure/perturbers/Perturber.hh>
#include <protocols/kinematic_closure/pivot_pickers/PivotPicker.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Edge.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/util.hh>
#include <core/id/AtomID.hh>
#include <core/id/TorsionID.hh>

// Protocols headers
#include <protocols/loops/Loop.hh>
#include <protocols/loop_modeling/loggers/Logger.hh>
#include <protocols/loop_modeling/loggers/NullLogger.hh>

// Numeric headers
#include <numeric/kinematic_closure/closure.hh>
#include <numeric/kinematic_closure/vector.hh>
#include <numeric/kinematic_closure/kinematic_closure_helpers.hh>
#include <numeric/conversions.hh>

// Utility headers
#include <utility/exit.hh>
#include <boost/foreach.hpp>
#include <iostream>

#define foreach BOOST_FOREACH
using namespace std;

namespace protocols {
namespace kinematic_closure {

using core::id::AtomID;
using numeric::kinematic_closure::operator <<;
using numeric::conversions::radians;
using numeric::conversions::from_radians;
using protocols::loop_modeling::loggers::LoggerOP;
using protocols::loop_modeling::loggers::NullLogger;

ClosureProblem::ClosureProblem() { // {{{1
	logger_ = new NullLogger();
	frame_called_ = false;
}

ClosureProblem::ClosureProblem(LoggerOP logger) { // {{{1
	logger_ = logger;
	frame_called_ = false;
}
// }}}1

void ClosureProblem::frame( // {{{1
		Pose const & pose, Loop const & loop,
		pivot_pickers::PivotPickerOP pivot_picker) {

	pivots_ = pivot_picker->pick(pose, loop);
	pivots_.set_cut(pivots_.cut() ? pivots_.cut() : pivots_.midpoint());

	frame_called_ = true;

	extract_cartesian_coordinates(pose,
			unperturbed_xyzs_);

	extract_internal_coordinates(pose,
			unperturbed_lengths_,
			unperturbed_angles_,
			unperturbed_torsions_);

	perturbed_torsions_ = ParameterList(unperturbed_torsions_);
	perturbed_angles_ = ParameterList(unperturbed_angles_);
	perturbed_lengths_ = ParameterList(unperturbed_lengths_);
}

void ClosureProblem::solve(SolutionList & solutions) const { // {{{1
	ParameterMatrix solution_torsions;
	ParameterMatrix solution_angles;
	ParameterMatrix solution_lengths;
	int num_solutions;

	IndexList order(3);
	order[1] = 1;
	order[2] = 2;
	order[3] = 3;

	runtime_assert(frame_called_);

	// Solve closure problem.
	numeric::kinematic_closure::radians::bridge_objects(
			unperturbed_xyzs_,
			perturbed_torsions_,
			perturbed_angles_,
			perturbed_lengths_,
			pivot_atoms(), order,
			solution_torsions,
			solution_angles,
			solution_lengths,
			num_solutions);

	// Cast the results into an easy-to-use format.
	solutions.clear();
	solutions.resize(num_solutions);

	for (Size index = 1; index <= num_solutions; index++) {
		solutions[index] = new ClosureSolution(
				this, index,
				solution_torsions[index],
				solution_angles[index],
				solution_lengths[index]);
	}
}

// {{{1

/// @details This method is typically used to restore the pose if no solutions 
/// to the closure problem could be found.  It should only be called after 
/// frame() and solve().

void ClosureProblem::restore(Pose & pose) const {
	runtime_assert(frame_called_);
	apply_internal_coordinates(
			unperturbed_lengths_, unperturbed_angles_, unperturbed_torsions_, pose);
}
// }}}1

void ClosureProblem::perturb_phi( // {{{1
		Size residue, Real value, AngleUnit unit) {

	Size loop_residue = residue - first_residue() + 1;
	perturbed_torsions_[3 * loop_residue + 1] = radians(value, unit);
}

void ClosureProblem::perturb_psi( // {{{1
		Size residue, Real value, AngleUnit unit) {

	Size loop_residue = residue - first_residue() + 1;
	perturbed_torsions_[3 * loop_residue + 2] = radians(value, unit);
}

void ClosureProblem::perturb_omega( // {{{1
		Size residue, Real value, AngleUnit unit) {

	Size loop_residue = residue - first_residue() + 1;
	perturbed_torsions_[3 * loop_residue + 0] = radians(value, unit);
}

void ClosureProblem::perturb_n_ca_c( // {{{1
		Size residue, Real value, AngleUnit unit) {

	Size loop_residue = residue - first_residue() + 1;
	perturbed_angles_[3 * loop_residue + 2] = radians(value, unit);
}

void ClosureProblem::perturb_ca_c_n( // {{{1
		Size residue, Real value, AngleUnit unit) {

	Size loop_residue = residue - first_residue() + 1;
	perturbed_angles_[3 * loop_residue + 0] = radians(value, unit);
}

void ClosureProblem::perturb_c_n_ca( // {{{1
		Size residue, Real value, AngleUnit unit) {

	Size loop_residue = residue - first_residue() + 1;
	perturbed_angles_[3 * loop_residue + 1] = radians(value, unit);
}

void ClosureProblem::perturb_n_ca(Size residue, Real value) { // {{{1
	Size loop_residue = residue - first_residue() + 1;
	perturbed_lengths_[3 * loop_residue + 1] = value;
}

void ClosureProblem::perturb_ca_c(Size residue, Real value) { // {{{1
	Size loop_residue = residue - first_residue() + 1;
	perturbed_lengths_[3 * loop_residue + 2] = value;
}

void ClosureProblem::perturb_c_n(Size residue, Real value) { // {{{1
	Size loop_residue = residue - first_residue() + 1;
	perturbed_lengths_[3 * loop_residue + 0] = value;
}

ClosureProblem::Memento::Memento(ClosureProblemOP problem) // {{{1
	: problem_(problem),
	  perturbed_torsions_(problem->perturbed_torsions_),
	  perturbed_angles_(problem->perturbed_angles_),
	  perturbed_lengths_(problem->perturbed_lengths_) {

	runtime_assert(problem->frame_called_);
}

void ClosureProblem::Memento::restore() const { // {{{1
	  problem_->perturbed_torsions_ = perturbed_torsions_;
	  problem_->perturbed_angles_ = perturbed_angles_;
	  problem_->perturbed_lengths_ = perturbed_lengths_;
}
// }}}1

Size ClosureProblem::first_residue() const { // {{{1
	runtime_assert(frame_called_);
	return pivots_.start();
}

Size ClosureProblem::middle_residue() const { // {{{1
	runtime_assert(frame_called_);
	return pivots_.cut();
}

Size ClosureProblem::last_residue() const { // {{{1
	runtime_assert(frame_called_);
	return pivots_.stop();
}

Size ClosureProblem::num_residues() const { // {{{1
	runtime_assert(frame_called_);
	return pivots_.length();
}

Size ClosureProblem::num_atoms() const { // {{{1
	return 3 * (num_residues() + 2);
}

Real ClosureProblem::phi(Size residue, AngleUnit unit) const { // {{{1
	Size loop_residue = residue - first_residue() + 1;
	return from_radians(perturbed_torsions_[3 * loop_residue + 1], unit);
}

Real ClosureProblem::psi(Size residue, AngleUnit unit) const { // {{{1
	Size loop_residue = residue - first_residue() + 1;
	return from_radians(perturbed_torsions_[3 * loop_residue + 2], unit);
}

Real ClosureProblem::omega(Size residue, AngleUnit unit) const { // {{{1
	Size loop_residue = residue - first_residue() + 1;
	return from_radians(perturbed_torsions_[3 * loop_residue + 0], unit);
}

Real ClosureProblem::n_ca_c(Size residue, AngleUnit unit) const { // {{{1
	Size loop_residue = residue - first_residue() + 1;
	return from_radians(perturbed_angles_[3 * loop_residue + 2], unit);
}

Real ClosureProblem::ca_c_n(Size residue, AngleUnit unit) const { // {{{1
	Size loop_residue = residue - first_residue() + 1;
	return from_radians(perturbed_angles_[3 * loop_residue + 0], unit);
}

Real ClosureProblem::c_n_ca(Size residue, AngleUnit unit) const { // {{{1
	Size loop_residue = residue - first_residue() + 1;
	return from_radians(perturbed_angles_[3 * loop_residue + 1], unit);
}

Real ClosureProblem::n_ca(Size residue) const { // {{{1
	Size loop_residue = residue - first_residue() + 1;
	return perturbed_lengths_[3 * loop_residue + 1];
}

Real ClosureProblem::ca_c(Size residue) const { // {{{1
	Size loop_residue = residue - first_residue() + 1;
	return perturbed_lengths_[3 * loop_residue + 2];
}

Real ClosureProblem::c_n(Size residue) const { // {{{1
	Size loop_residue = residue - first_residue() + 1;
	return perturbed_lengths_[3 * loop_residue + 0];
}

bool ClosureProblem::is_pivot_residue(Size residue) const { // {{{1
	return residue == first_residue() or
	       residue == middle_residue() or
	       residue == last_residue();
}

bool ClosureProblem::is_nonpivot_residue(Size residue) const { // {{{1
	return not is_pivot_residue(residue);
}

Loop ClosureProblem::pivot_loop() const { // {{{1
	return pivots_;
}

IndexList ClosureProblem::residues() const { // {{{1
	IndexList residues;

	for (Size i = first_residue(); i <= last_residue(); i++) {
		residues.push_back(i);
	}

	return residues;
}

IndexList ClosureProblem::pivot_residues() const { // {{{1
	IndexList pivot_residues(3);

	pivot_residues[1] = first_residue();
	pivot_residues[2] = middle_residue();
	pivot_residues[3] = last_residue();

	return pivot_residues;
}

IndexList ClosureProblem::nonpivot_residues() const { // {{{1
	IndexList nonpivot_residues;

	for (Size i = first_residue() + 1; i < last_residue(); i++) {
		if (i == middle_residue()) continue;
		nonpivot_residues.push_back(i);
	}

	return nonpivot_residues;
}

IndexList ClosureProblem::pivot_atoms() const { // {{{1
	IndexList pivot_atoms(3);

	Size middle_offset = middle_residue() - first_residue();
	Size end_offset = last_residue() - first_residue();

	pivot_atoms[1] = 5;
	pivot_atoms[2] = 5 + (3 * middle_offset);
	pivot_atoms[3] = 5 + (3 * end_offset);

	return pivot_atoms;
}
// }}}1

// {{{1
/// @details This method does not change anything about the problem object 
/// itself.  All changes are made to the given `atom_xyzs` pseudo-matrix.

void ClosureProblem::extract_cartesian_coordinates (
		Pose const & pose, CoordinateList & atom_xyzs) const {

	atom_xyzs.resize(num_atoms());

	// The closure algorithm requires coordinates for the two residues which 
	// frame the region being closed.  If the region coincides with the end of a 
	// chain, this is becomes a nontrivial problem.  For that reason, this step 
	// has been delegated out to the two functions called below.
	
	frame_lower_pivot(pose, atom_xyzs);
	frame_upper_pivot(pose, atom_xyzs);

	// The remaining code simply copies Cartesian coordinates from the given pose 
	// into the given matrix structure.  The copying starts with the 4th position 
	// in the matrix, to prevent overwriting the work done above.

	Size index = 4;

	for (Size i = first_residue(); i <= last_residue(); i++) {
		core::conformation::Residue const &residue = pose.residue(i);

		// This inner loop assumes that atoms 1, 2, & 3 are the backbone atoms.  
		// This will break if the residue in question is, for example, a metal ion.  
		// Metal ions only have one atom, so attempts to access the second and 
		// third ones will crash the program.

		for (Size j = 1; j <= 3; j++) {
			atom_xyzs[index].resize(3);
			atom_xyzs[index][1] = residue.xyz(j).x();
			atom_xyzs[index][2] = residue.xyz(j).y();
			atom_xyzs[index][3] = residue.xyz(j).z();
			index++;
		}
	}
}

void ClosureProblem::extract_internal_coordinates( // {{{1
		Pose const & pose,
		ParameterList & bond_lengths,
		ParameterList & bond_angles,
		ParameterList & torsion_angles) const {

	using core::conformation::Conformation;
	using core::id::AtomID;

	Conformation const & conformation = pose.conformation();
	AtomID ids[4];

	Size num_bond_lengths = num_atoms() - 1;
	Size num_bond_angles = num_atoms() - 2;
	Size num_torsion_angles = num_atoms() - 3;

	bond_lengths.resize(num_atoms());
	bond_angles.resize(num_atoms());
	torsion_angles.resize(num_atoms());

	for (Size index = 1; index <= num_atoms(); index++) {
		ids[0] = id_from_index(index + 0);
		ids[1] = id_from_index(index + 1);
		ids[2] = id_from_index(index + 2);
		ids[3] = id_from_index(index + 3);

		bond_lengths[index] = (index > num_bond_lengths) ?
			0 : conformation.bond_length(ids[0], ids[1]);

		bond_angles[(index % num_atoms()) + 1] = (index > num_bond_angles) ?
			0 : conformation.bond_angle(ids[0], ids[1], ids[2]);

		torsion_angles[(index % num_atoms()) + 1] = (index > num_torsion_angles) ?
			0 : conformation.torsion_angle(ids[0], ids[1], ids[2], ids[3]);
	}
}

void ClosureProblem::apply_internal_coordinates( // {{{1
		ParameterList const & bond_lengths,
		ParameterList const & bond_angles,
		ParameterList const & torsion_angles,
		Pose & pose) const {

	using core::conformation::Conformation;
	using core::id::AtomID;

	Conformation & conformation = pose.conformation();
	AtomID ids[4];

	Size num_bond_lengths = num_atoms() - 1;
	Size num_bond_angles = num_atoms() - 2;
	Size num_torsion_angles = num_atoms() - 3;

	// Set the bond lengths in the solution pose.
	for (Size index = 1; index <= num_bond_lengths; index++) {
		ids[0] = id_from_index(index + 0);
		ids[1] = id_from_index(index + 1);

		conformation.set_bond_length(
				ids[0], ids[1], bond_lengths[index]);
	}

	// Set the bond angles in the solution pose.  Note that in bond_angles, index 
	// x refers to the angle between atoms x-1 and x+1.
	for (Size index = 1; index <= num_bond_angles; index++) {
		ids[0] = id_from_index(index + 0);
		ids[1] = id_from_index(index + 1);
		ids[2] = id_from_index(index + 2);

		conformation.set_bond_angle(
				ids[0], ids[1], ids[2], bond_angles[index + 1]);
	}

	// Set the torsion angles in the solution pose.  Note that in torsion_angles,
	// index x refers to the angle between atoms x-1 and x+2.
	for (Size index = 1; index <= num_torsion_angles; index++) {
		ids[0] = id_from_index(index + 0);
		ids[1] = id_from_index(index + 1);
		ids[2] = id_from_index(index + 2);
		ids[3] = id_from_index(index + 3);

		conformation.set_torsion_angle(
					ids[0], ids[1], ids[2], ids[3], torsion_angles[index + 1]);
	}
}

void ClosureProblem::frame_lower_pivot( // {{{1
		Pose const & pose, CoordinateList & atom_xyzs) const {

	using core::conformation::Residue;
	using core::conformation::ResidueOP;

	ResidueOP frame_residue;
	Residue pivot_residue = pose.residue(first_residue());

	// The closure algorithm requires coordinates for the two residues framing 
	// the region being sampled.  This method finds those coordinates for the 
	// residue framing the lower (N-terminal) side of the region.  If the lower 
	// pivot is also the N-terminus, this residue doesn't exist and the needed 
	// coordinates are from a "pseudo-residue" that is artificially build with 
	// ideal coordinates.

	if (pivot_residue.is_lower_terminus()) {
		core::conformation::Conformation buffer;

		// Make a copy of the N-terminal residue.
		Residue pivot_residue_copy = core::conformation::Residue(
				pivot_residue.type(), pivot_residue,
				pose.conformation(), false);

		// Create a pseudo-residue of the same type.
		Residue pseudo_residue = Residue(
				pivot_residue.type(), true /* dummy argument to pick constructor. */);

		// Store the N-terminus in a temporary conformation.
		buffer.append_residue_by_bond(
				pivot_residue_copy);

		// Attach the pseudo-residue before the N-terminus in the copy.
		buffer.safely_prepend_polymer_residue_before_seqpos(
				pseudo_residue, 1, true);

		// Set the junction omega angle to its ideal value.
		buffer.set_torsion(
				core::id::TorsionID(1, core::id::BB, 3), 
				IdealParameters::omega_dihedral);

		// Get a pointer to the pseudo-residue.
		frame_residue = buffer.residue(1).clone();
	}

	// If the lower pivot isn't the N-terminus, then simply identify the frame 
	// residue as the one right below it.  This is the standard case.
	
	else {
		frame_residue = pose.residue(first_residue() - 1).clone();
	}

	// Append the coordinates of the frame residue into the list of atoms for the 
	// kinematic closure algorithm to consider.  The coordinates are added in the 
	// order: N, CA, C.

	for (Size i = 1; i <= 3; i++) {
		atom_xyzs[i].resize(3);
		atom_xyzs[i][1] = frame_residue->xyz(i).x();
		atom_xyzs[i][2] = frame_residue->xyz(i).y();
		atom_xyzs[i][3] = frame_residue->xyz(i).z();
	}
}

void ClosureProblem::frame_upper_pivot( // {{{1
		Pose const & pose, CoordinateList & atom_xyzs) const {

	using core::conformation::Residue;
	using core::conformation::ResidueOP;

	ResidueOP frame_residue;
	Residue pivot_residue = pose.residue(last_residue());

	// The closure algorithm requires coordinates for the two residues framing 
	// the region being sampled.  This method finds those coordinates for the 
	// residue framing the upper (C-terminal) side of the region.  This residue 
	// won't exist if the lower pivot is the C-terminus.  In this case, the 
	// needed coordinates are extracted from a "pseudo-residue" that is 
	// artificially build with ideal coordinates.

	if (pivot_residue.is_upper_terminus()) {
		core::conformation::Conformation buffer;

		// Make a copy of the N-terminal residue.
		Residue pivot_residue_copy = core::conformation::Residue(
				pivot_residue.type(), pivot_residue,
				pose.conformation(), false);

		// Create a pseudo-residue of the same type.
		Residue pseudo_residue = core::conformation::Residue(
				pivot_residue.type(), true /* dummy argument to pick constructor. */);

		// Store the N-terminus in a temporary conformation.
		buffer.append_residue_by_bond(
				pivot_residue_copy);

		// Attach the pseudo-residue before the N-terminus in the copy.
		buffer.safely_append_polymer_residue_after_seqpos(
				pseudo_residue, 1, true);

		// Set the junction omega angle to its ideal value.
		buffer.set_torsion(
				core::id::TorsionID(2, core::id::BB, 3), 
				IdealParameters::omega_dihedral);

		// Get a pointer to the pseudo-residue.
		frame_residue = buffer.residue(2).clone();
	}

	// If the upper pivot isn't the C-terminus, then simply identify the frame 
	// residue as the one right above it.  This is the standard case.
	
	else {
		frame_residue = pose.residue(last_residue() + 1).clone();
	}

	// Append the coordinates of the frame residue into the list of atoms for the 
	// kinematic closure algorithm to consider.  The coordinates are added in the 
	// order: N, CA, C.
	
	Size offset = atom_xyzs.size() - 3;

	for (Size i = 1; i <= 3; i++) {
		atom_xyzs[i + offset].resize(3);
		atom_xyzs[i + offset][1] = frame_residue->xyz(i).x();
		atom_xyzs[i + offset][2] = frame_residue->xyz(i).y();
		atom_xyzs[i + offset][3] = frame_residue->xyz(i).z();
	}
}

// {{{1
/// @details More specifically, the first atom in the loop is taken to be the 
/// backbone nitrogen of the residue preceding the first pivot.  Thus the index 
/// of the first pivot is 5 (pivots must be alpha carbons).

AtomID ClosureProblem::id_from_index(Size index) const {
	return AtomID(
			((index - 1) % 3) + 1,
			((index - 1) / 3) + first_residue() - 1);
}
// }}}1

} // end namespace kinematic_closure
} // end namespace protocols

