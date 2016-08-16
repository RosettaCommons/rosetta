// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   ClosureProblem.hh
/// @brief  Header file for ClosureProblem.
/// @author Kale Kundert (kale.kundert@ucsf.edu)

#ifndef INCLUDED_protocols_kinematic_closure_ClosureProblem_HH
#define INCLUDED_protocols_kinematic_closure_ClosureProblem_HH

// Unit headers
#include <protocols/kinematic_closure/types.hh>
#include <protocols/kinematic_closure/ClosureProblem.fwd.hh>
#include <protocols/kinematic_closure/ClosureSolution.fwd.hh>
#include <protocols/kinematic_closure/perturbers/Perturber.fwd.hh>
#include <protocols/kinematic_closure/pivot_pickers/PivotPicker.fwd.hh>

// Protocols headers
#include <protocols/loops/Loop.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/id/AtomID.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <boost/noncopyable.hpp>

namespace protocols {
namespace kinematic_closure {

/// @brief Represent and solve a kinematic closure problem.
///
/// @details A number of parameters are needed to define a kinematic closure
/// problem.  These parameters include a region of the protein backbone, a set
/// of pivot residues, and a new set of torsion angles for the nonpivot
/// residues.  All of this can be specified using the frame() and \b perturb_*
/// methods.  Once that has been done, the solve() method can be called to
/// return the complete set of backbone conformations consistent with the given
/// parameters.  These conformations are represented by the ClosureSolution
/// class.  The restore() method can be used to undo a solution once it has
/// been applied, which is useful in some Monte Carlo schemes.

class ClosureProblem
	: public utility::pointer::ReferenceCount, private boost::noncopyable {

	friend class ClosureSolution;

	// Constructors {{{1
public:

	/// @brief Default constructor.
	ClosureProblem();

	/// @brief Default destructor.
	~ClosureProblem();

	// Methods to solve the problem {{{1
public:

	/// @brief Choose pivots and copy coordinate information into the problem.
	void frame(
		Pose const & pose, Loop const & loop,
		pivot_pickers::PivotPickerOP pivot_picker);

	/// @brief Return every possible solution to this problem.
	SolutionList solve() const;

	/// @brief Undo any changes made to the pose during the last move.
	void restore(Pose & pose) const;

	// Methods to perturb the problem {{{1
public:

	/// @brief Perturb the given phi angle.
	void perturb_phi(Size residue, Real value, AngleUnit unit);

	/// @brief Perturb the given psi angle.
	void perturb_psi(Size residue, Real value, AngleUnit unit);

	/// @brief Perturb the given omega angle.
	void perturb_omega(Size residue, Real value, AngleUnit unit);

	/// @brief Perturb the N-CA-C bond angle in the given residue.
	void perturb_n_ca_c(Size residue, Real value, AngleUnit unit);

	/// @brief Perturb the CA-C-N bond angle in the given residue.
	void perturb_ca_c_n(Size residue, Real value, AngleUnit unit);

	/// @brief Perturb the C-N-CA bond angle in the given residue.
	void perturb_c_n_ca(Size residue, Real value, AngleUnit unit);

	/// @brief Perturb the N-CA bond length in the given residue.
	void perturb_n_ca(Size residue, Real value);

	/// @brief Perturb the CA-C bond length in the given residue.
	void perturb_ca_c(Size residue, Real value);

	/// @brief Perturb the C-N bond length in the given residue.
	void perturb_c_n(Size residue, Real value);

	/// @brief Provide non-const access to the raw torsion angle list.
	ParameterList & perturb_torsions();

	/// @brief Provide non-const access to the raw bond angle list.
	ParameterList & perturb_angles();

	/// @brief Provide non-const access to the raw bond length list.
	ParameterList & perturb_lengths();

	/// @brief Save the current state of the closure problem.  This is meant to
	/// facilitate undoing rejected perturbations when necessary.
	class Memento {

	public:
		Memento(ClosureProblemOP problem);
		void restore() const;

	private:
		ClosureProblemOP problem_;
		ParameterList const perturbed_torsions_;
		ParameterList const perturbed_angles_;
		ParameterList const perturbed_lengths_;
	};

	friend class Memento;

	// Attribute Access {{{1
public:

	/// @brief Return the index of the first pivot residue.
	Size first_residue() const;

	/// @brief Return the index of the second pivot residue.
	Size cut_residue() const;

	/// @brief Return the index of the third pivot residue.
	Size last_residue() const;

	/// @brief Return the number of residues in the loop.
	Size num_residues() const;

	/// @brief Return the number of atoms in the loop.
	Size num_atoms() const;

	/// @brief Return the current value of the given phi angle.
	Real phi(Size residue, AngleUnit unit) const;

	/// @brief Return the current value of the given psi angle.
	Real psi(Size residue, AngleUnit unit) const;

	/// @brief Return the current value of the given omega angle.
	Real omega(Size residue, AngleUnit unit) const;

	/// @brief Return the current value of the given N-CA-C bond angle.
	Real n_ca_c(Size residue, AngleUnit unit) const;

	/// @brief Return the current value of the given CA-C-N bond angle.
	Real ca_c_n(Size residue, AngleUnit unit) const;

	/// @brief Return the current value of the given C-N-CA bond angle.
	Real c_n_ca(Size residue, AngleUnit unit) const;

	/// @brief Return the current value of the given N-CA bond length.
	Real n_ca(Size residue) const;

	/// @brief Return the current value of the given CA-C bond length.
	Real ca_c(Size residue) const;

	/// @brief Return the current value of the given C-N bond length.
	Real c_n(Size residue) const;

	/// @brief Return true if the given residue is a pivot.
	bool is_pivot_residue(Size residue) const;

	/// @brief Return true if the given residue is not a pivot.
	bool is_nonpivot_residue(Size residue) const;

	/// @brief Return a loop specifying the three pivot residues.
	Loop pivot_loop() const;

	/// @brief Return the indices of every residue in the loop.
	IndexList residues() const;

	/// @brief Return the indices of all three pivot residues.
	IndexList pivot_residues() const;

	/// @brief Return the indices of every nonpivot residue in the loop.
	IndexList nonpivot_residues() const;

	/// @brief Return the indices of all three pivot atoms.
	IndexList pivot_atoms() const;

	// Private Helpers {{{1
private:

	/// @brief Extract cartesian coordinates for those backbone atoms which are
	/// relevant to the closure algorithm.  This includes all atoms from within
	/// the loop and from the two residues on either side of the loop.
	void extract_cartesian_coordinates(
		Pose const & pose,
		CoordinateList & atom_xyzs) const;

	/// @brief Extract torsion angles, bond angles, and bond lengths for those
	/// backbone atoms which are relevant to the closure algorithm.
	void extract_internal_coordinates(
		CoordinateList const & atom_xyzs,
		ParameterList & bond_lengths,
		ParameterList & bond_angles,
		ParameterList & torsion_angles) const;

	void extract_internal_coordinates(
		Pose const & pose,
		ParameterList & bond_lengths,
		ParameterList & bond_angles,
		ParameterList & torsion_angles) const;

	/// @brief Apply the given internal coordinates to the given pose.
	void apply_internal_coordinates(
		ParameterList const & bond_lengths,
		ParameterList const & bond_angles,
		ParameterList const & torsion_angles,
		Pose & pose) const;

	/// @brief Extract coordinates for the residue just beyond the N-terminal
	/// pivot.  These coordinates help anchor the closure problem.
	void frame_lower_pivot(
		Pose const & pose,
		CoordinateList &atom_xyzs) const;

	/// @brief Extract coordinates for the residue just beyond the C-terminal
	/// pivot.  These coordinates help anchor the closure problem.
	void frame_upper_pivot(
		Pose const & pose,
		CoordinateList &atom_xyzs) const;

	/// @brief Return true if the two given atoms are on opposite sides of the
	/// cutpoint specified by the loop that was used to create this problem.
	bool ids_span_cut(core::id::AtomID left, core::id::AtomID right) const;

	/// @brief Return an AtomID referring to the atom at the given index
	/// relative to the start of the loop.
	core::id::AtomID id_from_index(Size index) const;

	// Data Members {{{1
private:

	Loop pivots_;
	bool frame_called_;

	ParameterList perturbed_lengths_;
	ParameterList perturbed_angles_;
	ParameterList perturbed_torsions_;

	ParameterList unperturbed_lengths_;
	ParameterList unperturbed_angles_;
	ParameterList unperturbed_torsions_;
	CoordinateList unperturbed_xyzs_;

	// }}}1

};

} // end namespace kinematic_closure
} // end namespace protocols

#endif


