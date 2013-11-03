// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/enzdes/EnzdesMovers.hh
/// @brief a collection of movers that are used at different stages in enzyme design
/// @author Sinisa Bjelic, Florian Richter, floric@u.washington.edu


#ifndef INCLUDED_protocols_enzdes_EnzdesMovers_hh
#define INCLUDED_protocols_enzdes_EnzdesMovers_hh

#include <protocols/enzdes/EnzdesMovers.fwd.hh>
// AUTO-REMOVED #include <protocols/enzdes/EnzdesBaseProtocol.fwd.hh>

// Unit headers
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>

// Package headers
#include <core/scoring/constraints/AtomPairConstraint.fwd.hh>
// AUTO-REMOVED #include <core/scoring/constraints/Func.fwd.hh>
#include <core/scoring/constraints/MultiConstraint.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Project headers
// AUTO-REMOVED #include <core/conformation/Residue.fwd.hh>
// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/io/silent/SilentEnergy.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
// AUTO-REMOVED #include <core/scoring/constraints/ConstraintSet.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>
// AUTO-REMOVED #include <numeric/xyzVector.io.hh>

#include <utility/vector1.hh>


//Utility Headers

// C++ Headers

namespace protocols {
namespace enzdes {

class EnzdesConstraintReporter : public utility::pointer::ReferenceCount
{
public:
	EnzdesConstraintReporter();
	virtual ~EnzdesConstraintReporter();

	EnzdesConstraintReporter( EnzdesConstraintReporter const & );

	/// @brief Recurse through all the constraints in the pose to the ligand,
	/// through all the constraint-container constraints (e.g. Ambiguous constraints
	/// and multi-constraints) to find all the atoms that participate in various constraints
	/// to ligand atoms in the input Pose.
	void
	find_constraints_to_ligand(
		core::pose::Pose const &pose
	);

	/// @brief Read access to the set of atoms that participate in distance constraints
	/// to ligand atoms.
	utility::vector1< core::Size > const & constrained_lig_atoms() const {
		return constrained_lig_atoms_;
	}


	utility::vector1< core::id::AtomID > const & constrained_nonligand_atoms() const {
		return constrained_nonligand_atoms_;
	}

	/// @brief Set the (one) ligand residue index
	void ligand_resno(core::Size res_no) { ligand_seqpos_ = res_no; }

	/// @brief Get the (one) ligand residue index
	core::Size ligand_resno() const { return ligand_seqpos_; }

protected:

	void
	add_constrained_atoms_from_multiconstraint(
		core::scoring::constraints::MultiConstraintCOP real_multi_constraint
	);

	void
	add_constrained_atoms_from_atom_pair_constraint(
		core::scoring::constraints::AtomPairConstraintCOP atom_pair_constraint
	);

	void
	add_constrained_lig_atom(
		core::Size atom_no
	);

	void
	add_constrained_nonligand_atom(
		core::id::AtomID const & atid
	);


private:
	utility::vector1< core::Size   > constrained_lig_atoms_;
	utility::vector1< core::id::AtomID > constrained_nonligand_atoms_;
	core::Size ligand_seqpos_; //Ligand's sequence position


};

class PredesignPerturbMover : public protocols::rigid::RigidBodyPerturbMover
{

public:
  //Deafault constructor
  PredesignPerturbMover();

  //Deafault constructor
  ~PredesignPerturbMover();

	void
	set_docking_pose(
		core::pose::Pose &pose,
		core::pack::task::PackerTaskCOP task );

	void
	reinstate_pose(
		core::pose::Pose &pose,
		core::pose::Pose const &old_Pose );

	//void
	//add_constrained_lig_atoms_from_multiconstraint(
	//	core::scoring::constraints::MultiConstraintCOP real_multi_constraint );

	//void
	//add_constrained_lig_atom(
	//	core::Size atom_no );

	void
	find_constraints_to_ligand(
		core::pose::Pose const &pose );

	void
	apply(
		core::pose::Pose &pose );

	virtual std::string get_name() const;

	void
	set_ligand(core::Size res_no);

	/// @brief Allow an outside interpretter to the set of atoms for which there are AtomPair constraints to
	/// the ligand.
	utility::vector1< core::Size > const &
	get_constrained_lig_atoms() const
	{
		return constraint_reporter_.constrained_lig_atoms();
	}

	core::Vector
	find_rotation_center( core::pose::Pose const &pose );

	core::Vector
	find_geometric_center_for_constrained_lig_atoms(
		core::pose::Pose const &pose );

//parser functions

	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datamap ,
		protocols::filters::Filters_map const & ,
		protocols::moves::Movers_map const & ,
		core::pose::Pose const & );

	protocols::moves::MoverOP clone() const;

	protocols::moves::MoverOP fresh_instance() const;

private:
	EnzdesConstraintReporter constraint_reporter_;
	utility::vector1< core::Size >positions_to_replace_;
	//utility::vector1< core::Size >constrained_lig_atoms_;
	core::Size dock_trials_;
	core::pack::task::TaskFactoryOP task_factory_;
	//core::Size ligand_seqpos_;//Ligand's sequence position
};

/// @brief class that will identify the region around the ligand,
/// remove it, and then do a repack. It can also calculate the following
/// parameters: E diff after the repack, (in essence a crude delta G calc)
/// rmsd of the repacked site after the repack and rmsd of catalytic residues
class RepackLigandSiteWithoutLigandMover : public protocols::moves::Mover
{

public:
	RepackLigandSiteWithoutLigandMover();

	RepackLigandSiteWithoutLigandMover(
		core::scoring::ScoreFunctionCOP sfxn,
		bool calculate_silent_Es
	);

	~RepackLigandSiteWithoutLigandMover();

	void
	apply( core::pose::Pose & pose );

	virtual std::string get_name() const;

	void
	set_cstio(
		toolbox::match_enzdes_util::EnzConstraintIOCOP enzcst_io );

	void
	set_sfxn(
		core::scoring::ScoreFunctionCOP sfxn );

	void
	set_calculate_silent_Es(
		bool calculate );

    void
    set_separate_prt_ligand( bool separate_prt_ligand );


	utility::vector1< core::io::silent::SilentEnergy > const &
	silent_Es(){
		return silent_Es_; }

	core::pack::task::PackerTaskCOP get_ptask() const;

private:

	void
	separate_protein_and_ligand( core::pose::Pose & pose ) const;

	core::scoring::ScoreFunctionCOP sfxn_;
	core::Size lig_seqpos_;
	toolbox::match_enzdes_util::EnzConstraintIOCOP enzcst_io_;
	bool calculate_silent_Es_;
	core::pack::task::PackerTaskOP ptask_;
	utility::vector1< core::io::silent::SilentEnergy > silent_Es_;
    // PG 21-05-2013
    bool separate_prt_ligand_;
};


} //enzdes
} //protocols


#endif
