// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/ResfileReader.cc
/// @brief  implementation of resfile reader and its command classes
/// @author Gordon Lemmon (glemmon@gmail.com), adapted from the ResfileReader code
/// by Steven Lewis (smlewi@unc.edu) and Andrew Leaver-Fay

// Project Headers
#include <core/pose/Pose.hh>

// // Unit Headers
#include <protocols/ligand_docking/TetherLigand.hh>
#include <core/pose/util.hh>

#include <core/conformation/Conformation.hh>

#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>

// Utility Headers
#include <utility/exit.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>
// AUTO-REMOVED #include <utility/tag/Tag.hh>

//Auto Headers
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/tree/Atom.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


using basic::T;
using basic::Error;
using basic::Warning;

namespace protocols {
namespace ligand_docking {

static basic::Tracer tether_ligand_tracer("protocols.ligand_docking.ligand_options.Tether_ligand", basic::t_debug);

TetherLigand::TetherLigand(){}

TetherLigand::TetherLigand(const char & chain, const core::Real & angstroms):
		protocols::moves::Mover(),
		chain_(chain),
		angstroms_(angstroms)
{}

TetherLigand::TetherLigand(TetherLigand const & that):
		protocols::moves::Mover( that ),
		chain_(that.chain_),
		angstroms_(that.angstroms_), //size of one stdev for ligand restraint
		ligand_tether_(that.ligand_tether_)

{}

TetherLigand::~TetherLigand() {}

std::string TetherLigand::get_name() const{
	return "TetherLigand";
}

void
TetherLigand::apply( core::pose::Pose & pose ){
	core::Size chain_id= core::pose::get_chain_id_from_chain(chain_, pose);
	core::Size const residue_id = pose.conformation().chain_begin(chain_id);
	///TODO find the centroid positioned residue rather than just taking the first (above)
	ligand_tether_= restrain_ligand_nbr_atom(residue_id, angstroms_, pose);
}

void TetherLigand::release(core::pose::Pose & pose) {
	///TODO make this a const_iterator
	pose.remove_constraint(ligand_tether_);
}

core::scoring::constraints::ConstraintCOP const &
TetherLigand::get_ligand_tether() const {
	return ligand_tether_;
}

core::scoring::constraints::ConstraintCOP
restrain_ligand_nbr_atom(
		core::Size const lig_id,
		core::Real const stddev_Angstroms,
		core::pose::Pose & pose
){
	tether_ligand_tracer.Debug<< "stddev: " << stddev_Angstroms << std::endl;
	core::scoring::func::FuncOP const restraint_function = new core::scoring::func::HarmonicFunc(0, stddev_Angstroms);

	core::id::AtomID const fixed_pt(pose.atom_tree().root()->atom_id());

	core::conformation::Residue const & residue = pose.residue(lig_id);
	core::scoring::constraints::ConstraintCOP constraint = new core::scoring::constraints::CoordinateConstraint(
			core::id::AtomID( residue.nbr_atom(), lig_id),
			fixed_pt,
			residue.nbr_atom_xyz(),
			restraint_function
	);
	constraint = pose.add_constraint(constraint);

	return constraint;
}


} //namespace ligand_docking
} //namespace protocols
