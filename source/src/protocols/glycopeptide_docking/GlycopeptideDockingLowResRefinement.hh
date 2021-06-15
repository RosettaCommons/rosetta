// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/glycopeptide_docking/GlycopeptideDockingLowResRefinement.hh
/// @brief Low resolution refinement of substrate
/// @details Low resolution refinement of glycosyltransferase peptide/glycopeptide substrate complex.
/// This protocol allows enhanced sampling of the peptide/glycopeptide dofs in the low resolution mode with
/// simulated annealing. Current protocol was only tested for monosaccharides.
/// It does not treat the glycan-branch in a special manner. This behavior will change as the protocol is
/// developed in future versions for higher glycans.
/// @author Yashes Srinivasan (yashess@gmail.com)

#ifndef INCLUDED_protocols_glycopeptide_docking_GlycopeptideDockingLowResRefinement_HH
#define INCLUDED_protocols_glycopeptide_docking_GlycopeptideDockingLowResRefinement_HH

// Unit headers
#include <protocols/glycopeptide_docking/GlycopeptideDockingLowResRefinement.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/PyMOLMover.fwd.hh>
// Protocol headers
#include <protocols/glycopeptide_docking/GlycopeptideDockingFlags.fwd.hh>

#include <protocols/filters/Filter.fwd.hh>
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
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pack/task/operation/TaskOperations.fwd.hh>
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

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
//#include <utility/tag/XMLSchemaGeneration.fwd.hh> //transcluded from Mover

#include <string>

namespace protocols {
namespace glycopeptide_docking {

///@brief Low resolution refinement of substrate
/// @details Low resolution refinement primarly for sampling
/// peptide conformational space. We use simulated annealing
/// varying mc temperature from 2.0 to 0.6 (final) temperature.
/// The protocol applies the following steps over 30 cycles of simulated annealing:
/// 1. rigidbody perturber mover\n
/// 2. minimization across the jump\n
/// 3. Samples torsions of the peptide with small and shear mover\n
/// 4. Minimizes peptide\n

class GlycopeptideDockingLowResRefinement : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	GlycopeptideDockingLowResRefinement();

	/// @brief Constructor requiring flags, scorefunction sf and folding trees
	/// @details Constructor initialized with the GlycopeptideDockingFlagsOP flags \n
	/// (all options for setting up glycopeptide_docking), scorefunction sf and fold trees for peptide-enzyme
	/// docking and sampling the substrate. Currently, the use of ft_substrate is disabled as it
	/// worsened predictions.
	GlycopeptideDockingLowResRefinement( protocols::glycopeptide_docking::GlycopeptideDockingFlagsOP flags,
		core::scoring::ScoreFunctionOP sf , core::kinematics::FoldTreeOP ft_docking,
		core::kinematics::FoldTreeOP ft_substrate);

	/// @brief Constructor that only requires the GlycopeptideDockingFlagsOP flags
	/// and scorefunction sf.
	/// @details The foldtree is automatically setup with a util function based on options \n
	/// in the flags.
	GlycopeptideDockingLowResRefinement( protocols::glycopeptide_docking::GlycopeptideDockingFlagsOP flags,
		core::scoring::ScoreFunctionOP sf);


	/// @brief Copy constructor (not needed unless you need deep copies)
	GlycopeptideDockingLowResRefinement( GlycopeptideDockingLowResRefinement const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~GlycopeptideDockingLowResRefinement() override;

	/////////////////////
	/// Mover Methods ///
	/////////////////////

public:
	/// @brief Apply the mover
	void
	apply( core::pose::Pose & pose ) override;

	void
	show( std::ostream & output = std::cout ) const override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	clone() const override;

	std::string
	get_name() const override;

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

private: // methods

	void
	setup_backbone_and_jump_movemap();

	void
	sample_torsions( core::pose::Pose & pose, core::Size cycles, core::Real & acceptance_rate, core::kinematics::MoveMapOP mm );

private: // data

	/// @brief Flags object with options specifying glycopeptide_docking options
	protocols::glycopeptide_docking::GlycopeptideDockingFlagsOP flags_;
	/// @brief Scorefunction
	core::scoring::ScoreFunctionOP sf_;
	/// @brief Foldtree for peptide-enzyme docking
	core::kinematics::FoldTreeOP ft_docking_=nullptr;
	/// @brief Foldtree for sampling peptide degrees of freedon
	core::kinematics::FoldTreeOP ft_substrate_=nullptr;
	/// @brief Monte Carlo object for sampling
	protocols::moves::MonteCarloOP mc_;
	/// @brief Min mover for minimizing torsional dofs
	protocols::minimization_packing::MinMoverOP torsion_minimizer_;
	/// @brief Min mover for minmizing peptide-enzyme jump
	protocols::minimization_packing::MinMoverOP jump_minimizer_;
};

std::ostream &
operator<<( std::ostream & os, GlycopeptideDockingLowResRefinement const & mover );

} //protocols
} //glycopeptide_docking

#endif //protocols_glycopeptide_docking_GlycopeptideDockingLowResRefinement_HH
