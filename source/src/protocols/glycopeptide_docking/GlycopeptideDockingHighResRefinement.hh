// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/glycopeptide_docking/GlycopeptideDockingHighResRefinement.hh
/// @brief High resolution refinement of glycosyltransferase peptide/glycopeptide substrate complex.
/// This protocol applies MCM sampling in torsional and rigid body space with sidechain
/// repacking of the substrate and substrate-enzyme interface. The attractive and repulsive terms are
/// ramped down and up respectively to allow softer potential at the beginning of the run.
/// Protocol similar to high resolution refinement of flexpepdock.
/// @author Yashes Srinivasan (yashess@gmail.com)

#ifndef INCLUDED_protocols_glycopeptide_docking_GlycopeptideDockingHighResRefinement_HH
#define INCLUDED_protocols_glycopeptide_docking_GlycopeptideDockingHighResRefinement_HH

// Unit headers
#include <protocols/glycopeptide_docking/GlycopeptideDockingHighResRefinement.fwd.hh>
#include <protocols/moves/Mover.hh>
// Protocol headers
#include <protocols/glycopeptide_docking/GlycopeptideDockingFlags.fwd.hh>
#include <protocols/glycopeptide_docking/utils.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>

#include <protocols/constraint_movers/ConstraintSetMover.fwd.hh>

#include <protocols/minimization_packing/MinMover.fwd.hh>
#include <protocols/simple_moves/BackboneMover.fwd.hh>
#include <protocols/rigid/RigidBodyMover.fwd.hh>
#include <protocols/minimization_packing/PackRotamersMover.fwd.hh>
#include <protocols/simple_task_operations/RestrictToInterface.fwd.hh>
#include <protocols/task_operations/RestrictResiduesToRepackingOperation.fwd.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pack/task/operation/TaskOperations.fwd.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.fwd.hh>
#include <core/pack/task/operation/OperateOnResidueSubset.fwd.hh>
#include <core/select/residue_selector/ResidueSpanSelector.fwd.hh>
#include <core/select/residue_selector/NotResidueSelector.fwd.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/types.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreType.hh>


// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.fwd.hh>
#include <utility/excn/Exceptions.fwd.hh>

// Numeric headers
#include <numeric/random/random.fwd.hh>
#include <numeric/xyzVector.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
//#include <utility/tag/XMLSchemaGeneration.fwd.hh> //transcluded from Mover

#include <string>

namespace protocols {
namespace glycopeptide_docking {

///@brief High resolution refinement of glycosyltransferase substrate.
///@details Applies a set of samplers to the pose in high resolution mode.
/// At the beginning, the fa_atr term is upweighted and the fa_rep term is downweighted.
/// The terms are ramped down and up respectively to their score function weights over
/// n_cycles_. For  each cycles, the pose is subjected to: \n
/// 1. N cycles of rigid body moves
/// 1.1. Rigid body perturbation
/// 1.2 Sidechain packing of substrate
/// 1.3 Sidechain packing of interface
/// 1.4 Minimization across the substrate-enzyme jump
/// 1.5 Metropolis criterion (MCM with 1.4)
/// 2. N cycles of torsion moves "induced fit"
/// 2.1 Apply small or shear mover
/// 2.2 Sidechain packing of substrate
/// 2.3 Sidechain packing of interface
/// 2.4 Minimize and apply Metropolis criterion (MCM)

class GlycopeptideDockingHighResRefinement : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	/// @details Default constructors does not setup anything. Various getters and setters are provided
	/// to setup all options.
	GlycopeptideDockingHighResRefinement();

	/// @brief Constructor requiring flags, scorefunction sf and folding trees
	/// @details Constructor initialized with the GlycopeptideDockingFlagsOP flags \n
	/// (all options for setting up glycopeptide_docking), scorefunction sf and fold trees for peptide-enzyme
	/// docking and sampling the substrate. Currently, the use of ft_substrate is disabled as it
	/// worsened predictions.
	GlycopeptideDockingHighResRefinement( protocols::glycopeptide_docking::GlycopeptideDockingFlagsOP flags, core::scoring::ScoreFunctionOP sf, core::kinematics::FoldTreeOP ft_docking, core::kinematics::FoldTreeOP ft_substrate);

	/// @brief Constructor that only requires the GlycopeptideDockingFlagsOP flags
	/// and scorefunction sf.
	/// @details The foldtree is automatically setup with a util function based on options \n
	/// in the flags.
	GlycopeptideDockingHighResRefinement( protocols::glycopeptide_docking::GlycopeptideDockingFlagsOP flags, core::scoring::ScoreFunctionOP sf);

	/// @brief Copy constructor (not needed unless you need deep copies)
	GlycopeptideDockingHighResRefinement( GlycopeptideDockingHighResRefinement const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~GlycopeptideDockingHighResRefinement() override;

	/////////////////////
	/// Mover Methods ///
	/////////////////////

public:
	/// @brief setup default for the mover
	void
	setup();
	/// @brief Apply the mover
	void
	apply( core::pose::Pose & pose ) override;

	/// @brief show mover details
	void
	show( std::ostream & output = std::cout ) const override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	clone() const override;

	/// @brief get name
	std::string
	get_name() const override;

	/// @brief get mover name
	static
	std::string
	mover_name();

	/// Getters and setters for sf, foldtrees, flags
	void
	set_options(protocols::glycopeptide_docking::GlycopeptideDockingFlagsOP flags){
		flags_ = flags;
	}

	protocols::glycopeptide_docking::GlycopeptideDockingFlagsOP
	get_options() {
		return flags_;
	}

	void
	set_scorefunction(core::scoring::ScoreFunctionOP sf){
		sf_ = sf;
	}

	core::scoring::ScoreFunctionOP
	get_scorefunction() {
		return sf_;
	}

	void
	set_docking_foldtree(core::kinematics::FoldTreeOP ft_docking){
		ft_docking_ = ft_docking;
	}

	core::kinematics::FoldTreeCOP
	get_docking_foldtree() {
		return ft_docking_;
	}

	void
	set_substrate_foldtree(core::kinematics::FoldTreeOP ft_substrate){
		ft_substrate_ = ft_substrate;
	}

	core::kinematics::FoldTreeCOP
	get_substrate_foldtree() {
		return ft_substrate_;
	}

	void
	set_target_attractive_weight(core::Real w){
		target_atr_ = w;
	}

	core::Real
	get_target_attractive_weight(){
		return target_atr_;
	}

	void
	set_target_repulsive_weight(core::Real w){
		target_rep_ = w;
	}

	core::Real
	get_target_repulsive_weight(){
		return target_rep_;
	}

private: // methods

	/// @brief Ramps up weight for repulsive term and ramps down attractive LJ term.
	void
	ramp_score_weight( core::scoring::ScoreType const method,
		core::Real const target,
		core::Real const fraction_completion );


private: // data

	protocols::glycopeptide_docking::GlycopeptideDockingFlagsOP flags_;

	//////////////////////////
	//////// MC and Sf ramping parameters ///////
	//////////////////////////

	/// @brief Boltmann sampling temperature
	core::Real kt_;
	/// @brief Final weight of the attractive force term
	/// @details Set to be the same as ref2015 weight
	core::Real target_atr_;
	/// @brief Final weight of the repulsive force term
	/// @details Set to be the same as ref2015 weight
	core::Real target_rep_;
	/// @brief Fractional factor by which the attractive weight is ramped down
	core::Real starting_ramp_down_factor_;
	/// @brief Fractional factor by which the repulsive weight is ramped up
	core::Real starting_ramp_up_factor_;


	//////////////////////////
	//////// Scoring /////////
	//////////////////////////

	core::scoring::ScoreFunctionOP sf_;
	core::kinematics::FoldTreeOP ft_docking_;
	core::kinematics::FoldTreeOP ft_substrate_;
	protocols::moves::MonteCarloOP mc_;

	//////////////////////////
	//////// Movers //////////
	//////////////////////////

	protocols::rigid::RigidBodyPerturbMoverOP perturber_;
	// protocols::docking::FaDockingSlideIntoContactOP slider_;

	protocols::minimization_packing::MinMoverOP jump_minimizer_;
	protocols::minimization_packing::MinMoverOP torsion_minimizer_;

	protocols::simple_moves::SmallMoverOP small_mover_;
	protocols::simple_moves::ShearMoverOP shear_mover_;

	//protocols::LinkageConformerMoverOP linkage_conformer_mover_;

	protocols::minimization_packing::PackRotamersMoverOP packer_substrate_;
	protocols::minimization_packing::PackRotamersMoverOP packer_interface_;
};

std::ostream &
operator<<( std::ostream & os, GlycopeptideDockingHighResRefinement const & mover );

} //glycopeptide_docking
} //protocols

#endif //protocols_glycopeptide_docking_GlycopeptideDockingHighResRefinement_HH
