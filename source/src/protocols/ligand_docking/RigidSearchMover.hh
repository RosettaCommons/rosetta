// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/ligand_docking/RigidSearchMover.hh
///
/// @brief
/// @author Ian W. Davis


#ifndef INCLUDED_protocols_ligand_docking_RigidSearchMover_hh
#define INCLUDED_protocols_ligand_docking_RigidSearchMover_hh

#include <protocols/ligand_docking/RigidSearchMover.fwd.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/moves/Mover.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace ligand_docking {


/// @brief An optimized mover for Monte Carlo trial of rigid body perturbations.
///
/// @details
///
class RigidSearchMover : public protocols::moves::Mover
{
public:

	RigidSearchMover(int jump_id, int num_trials, core::scoring::ScoreFunctionCOP scorefxn);
	RigidSearchMover(RigidSearchMover const & that);
	virtual ~RigidSearchMover();

	virtual void apply(core::pose::Pose & pose);
	virtual std::string get_name() const;

	/// @brief Will the absolute lowest-energy pose be recovered at the end of apply()?
	bool recover_low() const { return recover_low_; }
	/// @brief Should the absolute lowest-energy pose be recovered at the end of apply()?
	void recover_low(bool b) { recover_low_ = b; }

	/// @brief Amount of random (Gaussian) translation, in Angstroms
	core::Real translation() const { return translate_Ang_; }
	/// @brief Amount of random (Gaussian) translation, in Angstroms
	void translation(core::Real angstroms) { translate_Ang_ = angstroms; }

	/// @brief Amount of random (Gaussian) rotation, in degrees
	core::Real rotation() const { return rotate_deg_; }
	/// @brief Amount of random (Gaussian) rotation, in degrees
	void rotation(core::Real degrees) { rotate_deg_ = degrees; }

	/// @brief Rotation occurs around centroid of downstream half of the jump (default)
	void rotate_around_downstream_centroid() { rotate_rsd_ = 0; rotate_atom_ = 0; }
	/// @brief Rotation occurs around the specified atom
	void rotate_around_atom(core::Size rsd, core::Size atom) { rotate_rsd_ = rsd; rotate_atom_ = atom; }

private:
	int jump_id_;
	int num_trials_;
	core::scoring::ScoreFunctionCOP scorefxn_;
	core::Real temperature_;
	core::Real rotate_deg_;
	core::Real translate_Ang_;
	// If these are 0, rotate around the downstream centroid, else rotate around this atom.
	core::Size rotate_rsd_;
	core::Size rotate_atom_;
	bool recover_low_;

}; // RigidSearchMover


} // namespace ligand_docking
} // namespace protocols

#endif // INCLUDED_protocols_ligand_docking_RigidSearchMover_HH
