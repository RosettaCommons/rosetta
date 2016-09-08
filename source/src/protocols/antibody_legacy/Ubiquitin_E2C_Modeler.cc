// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @author Aroop Sircar ( aroopsircar@yahoo.com )


#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/util.hh>
#include <core/pose/util.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <basic/options/option.hh>
#include <core/pack/rotamer_set/UnboundRotamersOperation.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <core/pack/task/operation/OperateOnCertainResidues.hh>
#include <core/pack/task/operation/ResFilters.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <protocols/toolbox/task_operations/RestrictToInterface.hh>
#include <protocols/scoring/Interface.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pack/dunbrack/RotamerConstraint.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/DiagnosticData.hh>
#include <basic/Tracer.hh>

#include <numeric/random/random.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
using namespace ObjexxFCL::format;

#include <protocols/antibody_legacy/Ubiquitin_E2C_Modeler.hh>
#include <protocols/docking/DockingInitialPerturbation.hh>
#include <protocols/rigid/RB_geometry.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/ConstraintSetMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/simple_moves/ReturnSidechainMover.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/simple_moves/RotamerTrialsMinMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/jd2/ScoreMap.hh>

//#include <utility/vector1.hh>

#include <protocols/scoring/InterfaceInfo.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/chemical/ChemicalManager.fwd.hh>

using basic::T;
using basic::Error;
using basic::Warning;

static THREAD_LOCAL basic::Tracer TR( "protocols.ub_e2c.ubi_e2c_modeler" );

namespace protocols {
namespace ub_e2c {

#ifndef WIN32
// Aroop's Magic number, do not change it
// (and dont try and use it anywhere else)

using namespace core;
using namespace core::scoring;

// default constructor
ubi_e2c_modeler::ubi_e2c_modeler() : Mover(),
	e2_k48r_jump_( 1 ),
	e2_d77_jump_( 2 ),
	max_k48_cter_dist_( 15.00 ),
	temperature_( 0.8 ),
	flex_cter_( 3 ),
	CA( 2 ),
	max_repeats_( 100000 ),
	centroid_allowed_CSP_fraction_( 0.7 ),
	fullatom_allowed_CSP_fraction_( 0.7 ),
	centroid_CSP_weight_( 5.0 ),
	centroid_non_CSP_weight_( 5.0 ),
	fullatom_CSP_weight_( 1000.0 ),
	fullatom_non_CSP_weight_( 1000.0 ),
	cen_vdw_( 1.00 ),
	cen_constraint_( 10.00 ),
	// full_vdw_( 1.00 ),
	refinement_mode_( false ), // false for initial run
	cov_bond_only_flag_( true ), // false for initial run
	monoub_mode_( false ), // false for our model system ::: SET THIS TO FALSE
	higher_d77_pert_mode_( false ) { // false for regular runs ::: SET THIS TO TRUE
	Mover::type( "ubi_e2c_modeler" );
	set_default();
}

// default destructor
ubi_e2c_modeler::~ubi_e2c_modeler() = default;

//clone
protocols::moves::MoverOP ubi_e2c_modeler::clone() const {
	return( protocols::moves::MoverOP( new ubi_e2c_modeler() ) );
}

void ubi_e2c_modeler::set_default() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	TR <<  "UBI Setting up default settings" << std::endl;

	// README
	// The 76th residue (C-terminal) of Ub1 (UbK48R) is linked to the 48th
	// Lysine of Ub2 (UbD77)
	// The 48th residue of Ub1 (UbK48R) has a spin labeled tag
	// The 75th residue of Ub2 (UbD77) has a spin labeled tag
	e2_ctr_of_mass_ = 0;
	k48r_ctr_of_mass_ = 0;
	d77_ctr_of_mass_ = 0;
	e2_end_ = 0;
	k48r_end_ = 0;
	d77_end_ = 0;
	k48r_48_lys_ = 0;
	d77_48_lys_ = 0;
	k48r_trim_end_ = 0 ;
	d77_trim_end_ = 0;
	d77_trim_ctr_mass_ = 0;

	k48r_swap_ = false;

	k48r_48k_mtsl_ = 0;
	d77_75g_mtsl_ = 0;

	passed_centroid_filter_ = true;
	passed_fullatom_filter_ = true;

	applied_fullatom_pert_ = false;

	min_tolerance_ = 0.1;
	min_type_ = std::string( "lbfgs_armijo_nonmonotone" );
	nb_list_ = true;

	if ( cov_bond_only_flag_ ) {
		full_constraint_ = 100000.00;
	} else {
		full_constraint_ = 100.00;
	}

	if ( monoub_mode_ ) {
		fullatom_constraint_cutoff_ = 3 * 1.5;
	} else {
		fullatom_constraint_cutoff_ = 1.5;
	}

	dock_lowres_scorefxn_ = core::scoring::ScoreFunctionFactory::
		create_score_function("interchain_cen");
	dock_lowres_scorefxn_->set_weight( core::scoring::interchain_vdw, cen_vdw_ );
	dock_lowres_scorefxn_->set_weight( core::scoring::interchain_contact,
		cen_constraint_ );

	lowres_scorefxn_ = core::scoring::ScoreFunctionFactory::
		create_score_function( "cen_std" );
	lowres_scorefxn_->set_weight( core::scoring::vdw, cen_vdw_ );
	lowres_scorefxn_->set_weight( core::scoring::cbeta, cen_vdw_ );

	dock_lowres_cst_scorefxn_ = core::scoring::ScoreFunctionFactory::
		create_score_function("interchain_cen");
	// adding constraints
	dock_lowres_cst_scorefxn_->set_weight( core::scoring::atom_pair_constraint,
		cen_constraint_ );
	dock_lowres_cst_scorefxn_->set_weight( core::scoring::interchain_vdw,cen_vdw_);
	dock_lowres_cst_scorefxn_->set_weight( core::scoring::interchain_contact,
		cen_vdw_ );

	lowres_cst_scorefxn_ = core::scoring::ScoreFunctionFactory::
		create_score_function( "cen_std" );
	// adding constraints
	lowres_cst_scorefxn_->set_weight( core::scoring::atom_pair_constraint,
		cen_constraint_ );
	lowres_cst_scorefxn_->set_weight( core::scoring::vdw, cen_vdw_ );
	lowres_cst_scorefxn_->set_weight( core::scoring::cbeta, cen_vdw_ );

	pack_scorefxn_ = core::scoring::get_score_function_legacy("pre_talaris_2013_standard.wts");

	dockfa_scorefxn_ = core::scoring::ScoreFunctionFactory::
		create_score_function( "docking" );
	//dockfa_scorefxn_->set_weight( core::scoring::fa_sol, full_vdw_ );
	//dockfa_scorefxn_->set_weight( core::scoring::fa_rep, full_vdw_ );

	dockfa_min_scorefxn_ = core::scoring::ScoreFunctionFactory::
		create_score_function( "docking", "docking_min" );
	//dockfa_min_scorefxn_->set_weight( core::scoring::fa_sol, full_vdw_ );
	//dockfa_min_scorefxn_->set_weight( core::scoring::fa_rep, full_vdw_ );

	dockfa_cst_scorefxn_ = core::scoring::ScoreFunctionFactory::
		create_score_function( "docking" );
	//dockfa_cst_scorefxn_->set_weight( core::scoring::fa_sol, full_vdw_ );
	//dockfa_cst_scorefxn_->set_weight( core::scoring::fa_rep, full_vdw_ );
	dockfa_cst_scorefxn_->set_weight( core::scoring::atom_pair_constraint,
		full_constraint_ );

	dockfa_cst_min_scorefxn_ = core::scoring::ScoreFunctionFactory::
		create_score_function( "docking", "docking_min" );
	//dockfa_cst_min_scorefxn_->set_weight( core::scoring::fa_sol, full_vdw_ );
	//dockfa_cst_min_scorefxn_->set_weight( core::scoring::fa_rep, full_vdw_ );
	dockfa_cst_min_scorefxn_->set_weight( core::scoring::atom_pair_constraint,
		full_constraint_ );

	pack_cst_scorefxn_ = core::scoring::get_score_function("pre_talaris_2013_standard.wts");
	// adding constraints
	pack_cst_scorefxn_->set_weight( core::scoring::atom_pair_constraint,
		full_constraint_ );
	// pack_cst_scorefxn_->set_weight( core::scoring::fa_sol, full_vdw_ );

	output_cen_scorefxn_ = core::scoring::ScoreFunctionFactory::
		create_score_function( "cen_std" );

	output_full_scorefxn_ = core::scoring::ScoreFunctionFactory::
		create_score_function( "docking" );

	TR <<  "UBI Done:Setting up default settings" << std::endl;

	return;
} // set_default

void ubi_e2c_modeler::compute_trim_CSPs() {

	TR << "UBI Assigning Trimmed CSPs" << std::endl;

	for ( Size i = 1; i <= CSP_.size(); i++ ) {
		if ( CSP_[i] <= ( k48r_end_ - ( flex_cter_ + 1)) ) {
			CSP_trim_.push_back( CSP_[i] );
		} else if ( CSP_[i] <= ( e2_end_ - ( flex_cter_ + 1)) ) {
			CSP_trim_.push_back( CSP_[i] - (flex_cter_ + 1));
		}
	}

	for ( Size i = 1; i <= non_CSP_.size(); i++ ) {
		if ( non_CSP_[i] <= ( k48r_end_ - ( flex_cter_ + 1)) ) {
			non_CSP_trim_.push_back( non_CSP_[i] );
		} else if ( non_CSP_[i] <= ( e2_end_ - ( flex_cter_ + 1)) ) {
			non_CSP_trim_.push_back( CSP_[i] - (flex_cter_ + 1));
		}
	}

	TR << "UBI Done: Assigning Trimmed CSPs" << std::endl;

	return;
} // compute_trim_CSPs

void ubi_e2c_modeler::compute_swap_trim_CSPs() {

	TR << "UBI Assigning Swap Trimmed CSPs" << std::endl;

	Size ub_trim = d77_trim_end_ - k48r_trim_end_;

	for ( Size i = 1; i <= CSP_.size(); i++ ) {
		if ( CSP_[i] <= e2_end_ ) {
			CSP_swap_trim_.push_back( CSP_[i] );
		} else if ( CSP_[i] <= k48r_end_ - (flex_cter_ + 1 ) ) {
			CSP_swap_trim_.push_back( CSP_[i] + ub_trim );
		} else if ( CSP_[i] <= d77_end_ - (flex_cter_ + 1 ) ) {
			CSP_swap_trim_.push_back( CSP_[i] - ( ub_trim  + flex_cter_ + 1));
		}
	}

	for ( Size i = 1; i <= non_CSP_.size(); i++ ) {
		if ( non_CSP_[i] <= e2_end_ ) {
			non_CSP_swap_trim_.push_back( non_CSP_[i] );
		} else if ( non_CSP_[i] <= k48r_end_ - (flex_cter_ + 1 ) ) {
			non_CSP_swap_trim_.push_back( non_CSP_[i] + ub_trim );
		} else if ( non_CSP_[i] <= d77_end_ - (flex_cter_ + 1 ) ) {
			non_CSP_swap_trim_.push_back( non_CSP_[i]-(ub_trim + flex_cter_ +1));
		}
	}

	TR << "UBI Done: Assigning Swap Trimmed CSPs" << std::endl;

	return;
} // compute_swap_trim_CSPs

void ubi_e2c_modeler::assign_CSPs(
	const pose::Pose & pose_in ) {

	TR << "UBI Assigning CSPs" << std::endl;

	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 23 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 24 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 25 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 26 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 29 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 32 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 33 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 34 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 36 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 40 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 41 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 42 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 52 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 53 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 54 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 57 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 58 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 74 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 160 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 164 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 14 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 17 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 18 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 27 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 35 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 38 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 39 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 43 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 46 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 50 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 51 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 56 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 59 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 60 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 90 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 163 ) );

	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 7 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 13 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 45 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 47 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 48 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 49 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 65 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 68 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 69 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 70 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 71 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 73 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 74 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 6 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 8 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 11 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 14 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 23 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 32 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 42 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 43 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 51 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 54 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 66 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 67 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 72 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 76 ) );

	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 7 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 43 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 48 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 70 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 71 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 23 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 29 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 45 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 51 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 69 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 73 ) );

	assign_non_CSPs( pose_in );

	compute_trim_CSPs();
	compute_swap_trim_CSPs();

	TR << "UBI Done: Assigning CSPs" << std::endl;

	return;
} // assign_CSPs

void ubi_e2c_modeler::assign_non_CSPs(
	const pose::Pose & pose_in ) {

	TR << "UBI Assigning non CSPs" << std::endl;

	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A',   4 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A',   5 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A',   6 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A',   7 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A',   8 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A',   9 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A',  10 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A',  11 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A',  12 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A',  13 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A',  15 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A',  16 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A',  19 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A',  22 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A',  30 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A',  31 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A',  37 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A',  45 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A',  47 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A',  48 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A',  49 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A',  62 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A',  63 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A',  64 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A',  66 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A',  67 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A',  70 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A',  71 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A',  72 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A',  73 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A',  75 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A',  76 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A',  77 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A',  78 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A',  79 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A',  81 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A',  82 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A',  83 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A',  85 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A',  86 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A',  87 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A',  88 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A',  89 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A',  91 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A',  92 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A',  93 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A',  94 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A',  95 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A',  97 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A',  98 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A',  99 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 101 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 102 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 103 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 104 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 105 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 106 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 107 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 108 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 109 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 110 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 111 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 113 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 114 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 115 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 116 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 117 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 118 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 119 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 120 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 121 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 122 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 123 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 124 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 125 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 126 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 127 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 128 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 129 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 131 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 132 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 133 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 134 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 135 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 136 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 137 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 138 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 139 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 140 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 141 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 142 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 143 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 144 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 145 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 146 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 147 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 148 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 149 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 150 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 151 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 152 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 153 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 154 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 155 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 156 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 157 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 158 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 159 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 161 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 162 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 165 ) );

	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B',  2 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B',  3 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B',  4 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B',  5 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 10 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 12 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 15 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 16 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 17 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 18 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 20 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 21 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 22 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 25 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 26 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 27 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 28 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 29 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 30 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 31 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 33 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 34 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 35 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 36 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 39 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 40 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 41 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 44 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 50 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 52 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 55 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 56 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 57 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 58 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 59 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 60 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 61 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 62 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 63 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 64 ) );

	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C',  2 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C',  3 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C',  4 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C',  5 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C',  6 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C',  8 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 10 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 11 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 12 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 13 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 14 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 15 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 16 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 17 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 18 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 20 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 21 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 22 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 25 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 26 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 27 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 28 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 30 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 31 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 32 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 33 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 34 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 35 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 36 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 39 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 40 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 41 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 42 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 44 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 47 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 49 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 50 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 52 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 54 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 55 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 56 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 57 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 58 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 59 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 60 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 61 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 62 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 63 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 64 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 65 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 66 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 67 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 68 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 72 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 74 ) );
	non_CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'C', 76 ) );

	TR << "UBI Done: Assigning non CSPs" << std::endl;

	return;
} // assign_non_CSPs


void ubi_e2c_modeler::apply( pose::Pose & pose_in ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace chemical;
	using namespace id;
	//using namespace jobdist;
	//using pose::datacache::CacheableDataType::SCORE_MAP;
	using utility::file::FileName;

	// control monoubiquitin monomer
	if ( monoub_mode_ ) {
		monoub_apply( pose_in );
		return;
	}

	/*
	{
	evaluate_native( pose_in );
	return;
	}
	*/

	protocols::simple_moves::ConstraintSetMoverOP mtsl_constraint( new protocols::simple_moves::ConstraintSetMover() );
	mtsl_constraint->apply( pose_in );

	const pose::Pose start_pose( pose_in );

	// setup key residues
	setup_key_residues( start_pose );

	// assign Chemical Shift Perturbations
	assign_CSPs( start_pose );

	// setup Move Maps
	setup_move_maps();

	// Residue movers
	protocols::simple_moves::SwitchResidueTypeSetMover to_centroid( chemical::CENTROID );
	protocols::simple_moves::SwitchResidueTypeSetMover to_all_atom( chemical::FA_STANDARD );
	protocols::simple_moves::ReturnSidechainMover recover_sidechains( start_pose);

	// centroid mode
	// Filter distance for K48-Cterminal
	Real ubi_cov_bond_dist( 100.00 );
	Size tries = 0;

	if ( !cov_bond_only_flag_ ) {
		//start loop of decoy creation until filters are all passed
		for ( Size r = 1; r <= max_repeats_; r++ ) {
			tries = r;
			pose_in = start_pose;

			if ( !refinement_mode_ ) {
				// convert to centroid mode
				to_centroid.apply( pose_in );
				ubi_cov_bond_dist = centroid_mode_perturbation( pose_in );

				// add scores to map for output
				score_map_[ "AD_k48_CA_CA" ] = ubi_cov_bond_dist;
				protocols::jd2::ScoreMap::nonzero_energies( score_map_, lowres_cst_scorefxn_, pose_in);

				// check low-res docking filter
				( *dock_lowres_cst_scorefxn_ )( pose_in );
				passed_centroid_filter_ = centroid_filter( pose_in );
			} // if( !refinement_mode_ )

			applied_fullatom_pert_ = false;

			// fullatom mode
			if ( passed_centroid_filter_ || r == max_repeats_ || refinement_mode_ ) {
				//pose_in.dump_pdb( "main_passed_centroid.pdb" );
				// convert to full atom
				if ( !refinement_mode_ ) {
					to_all_atom.apply( pose_in );
					recover_sidechains.apply( pose_in );
				} // if( !refinement_mode_ )
				ubi_cov_bond_dist = fullatom_mode_perturbation( pose_in );

				// add scores to map for output
				score_map_[ "AD_k48_CA_CA" ] = ubi_cov_bond_dist;
				protocols::jd2::ScoreMap::nonzero_energies( score_map_, pack_cst_scorefxn_, pose_in);

				// check highres docking filter
				( *dockfa_cst_scorefxn_ )( pose_in );

				passed_fullatom_filter_ = fullatom_filter( pose_in );
			} // if fullatom mode
			// pose_in.dump_pdb( "fullatom.pdb" );

			if ( passed_centroid_filter_ && passed_fullatom_filter_ ) break;
			else  TR <<"UBI Repeating structure " << r << " times" << std::endl;

		} // for max_repeats_
	}

	if ( !refinement_mode_ ) {
		optimize_cov_bond( pose_in );
		( *dockfa_cst_scorefxn_ )( pose_in );
		fullatom_filter( pose_in );
	} // if ( !refinement_mode_ )

	// add scores to map for output
	if ( applied_fullatom_pert_ || cov_bond_only_flag_ ) {
		( *output_full_scorefxn_ )( pose_in );
		protocols::jd2::ScoreMap::nonzero_energies(score_map_, output_full_scorefxn_, pose_in);
	} else {
		( *output_cen_scorefxn_ )( pose_in );
		protocols::jd2::ScoreMap::nonzero_energies( score_map_, output_cen_scorefxn_, pose_in);
	}

	if ( !cov_bond_only_flag_ ) {
		score_map_["AJ_k48r_rms"] = calc_Lrmsd( pose_in, start_pose,
			e2_k48r_jump_ );
		score_map_["AK_d77_rms"] = calc_Lrmsd( pose_in, start_pose,
			e2_d77_jump_ );
		score_map_["AI_rms"] = score_map_["AJ_k48r_rms"] +
			score_map_["AK_d77_rms"];
	} else {
		score_map_["AL_k48r_tail_rmsg"] = calc_Lrmsd( pose_in, start_pose,
			flex_cter_ ); // option
	}

	using namespace basic::datacache;
	pose_in.data().set( core::pose::datacache::CacheableDataType::SCORE_MAP, DataCache_CacheableData::DataOP( new basic::datacache::DiagnosticData(score_map_) ) );

	TR << "UBI Outputing structure after " << tries << " times" << std::endl;

	//constraints::ConstraintFactory cst_factory;
	//constraints::ConstraintOP cst_op;
	//cst_op = cst_factory.newConstraint( "AtomPair" );
	//cst_op->read_def( "AtomPair CA 213 CA 74 BOUNDED 0.0 15.0 4.0", pose_in, func_factory );
}// end apply

std::string
ubi_e2c_modeler::get_name() const {
	return "ubi_e2c_modeler";
}

void
ubi_e2c_modeler::setup_key_residues(
	const pose::Pose & pose_in ) {

	TR << "UBI Setting Up Key Residues" << std::endl;

	pose::PDBInfoCOP pdb_info = pose_in.pdb_info();
	d77_end_ = pose_in.size();

	//char chain = '_';
	char old_chain = '_';
	Size next_chain( 0 );
	for ( Size i = 1; i <= d77_end_; i++ ) {
		char chain = pdb_info->chain( i );

		// if initial condition
		if ( i == 1 ) {
			old_chain = pdb_info->chain( i );
		}

		if ( chain != old_chain ) {
			if ( next_chain == 0  ) {
				next_chain = 1;
				e2_end_ = i - 1;
			} else {
				k48r_end_ = i -1;
			}
		} // if new chain
		old_chain = chain;
	}// for i <= d77_end_

	e2_ctr_of_mass_ = core::pose::residue_center_of_mass( pose_in, 1,
		e2_end_ );
	k48r_ctr_of_mass_ = core::pose::residue_center_of_mass( pose_in,
		e2_end_ + 1,
		k48r_end_ );
	d77_ctr_of_mass_ = core::pose::residue_center_of_mass( pose_in,
		k48r_end_ + 1,
		d77_end_ );

	// C-terminal trimmed parameters default
	k48r_trim_end_ = k48r_end_;
	d77_trim_end_ = d77_end_;
	d77_trim_ctr_mass_ = d77_ctr_of_mass_;

	k48r_48_lys_ = pose_in.pdb_info()->pdb2pose( 'B', 48 );
	d77_48_lys_ = pose_in.pdb_info()->pdb2pose( 'C', 48 );
	k48r_48k_mtsl_ = pose_in.pdb_info()->pdb2pose( 'B', 48 );
	d77_75g_mtsl_ = pose_in.pdb_info()->pdb2pose( 'C', 75 );

	TR << "UBI Done: Setting Up Key Residues" << std::endl;

	//utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
} // setup_key_residues

void
ubi_e2c_modeler::setup_move_maps() {
	using namespace core;
	using namespace kinematics;

	if ( init_all_dof_map_ ) {
		all_dof_map_ = core::kinematics::MoveMapOP( new MoveMap( *init_all_dof_map_ ) );
		k48r_docking_map_ = core::kinematics::MoveMapOP( new MoveMap( *init_k48r_docking_map_ ) );
		d77_docking_map_ = core::kinematics::MoveMapOP( new MoveMap( *init_d77_docking_map_ ) );
		docking_map_ = core::kinematics::MoveMapOP( new MoveMap( *init_docking_map_ ) );
		flex_cter_map_ = core::kinematics::MoveMapOP( new MoveMap( *init_flex_cter_map_ ) );

		TR << "UBI Reinitializing Move Maps" << std::endl;
		return;
	}

	TR << "UBI Setting Up Move Maps" << std::endl;

	bool bb = false;
	bool chi = true;

	k48r_docking_map_ = core::kinematics::MoveMapOP( new MoveMap() );
	k48r_docking_map_->clear();
	k48r_docking_map_->set_chi( chi );
	k48r_docking_map_->set_bb( bb );
	k48r_docking_map_->set_jump( e2_k48r_jump_, true );
	k48r_docking_map_->set_jump( e2_d77_jump_, false );

	d77_docking_map_ = core::kinematics::MoveMapOP( new MoveMap() );
	d77_docking_map_->clear();
	d77_docking_map_->set_chi( chi );
	d77_docking_map_->set_bb( bb );
	d77_docking_map_->set_jump( e2_k48r_jump_, false );
	d77_docking_map_->set_jump( e2_d77_jump_, true );

	docking_map_ = core::kinematics::MoveMapOP( new MoveMap() );
	docking_map_->clear();
	docking_map_->set_chi( chi );
	docking_map_->set_bb( bb );
	docking_map_->set_jump( e2_k48r_jump_, true );
	docking_map_->set_jump( e2_d77_jump_, true );

	flex_cter_map_ = core::kinematics::MoveMapOP( new MoveMap() );
	flex_cter_map_->clear();
	flex_cter_map_->set_chi( chi );
	flex_cter_map_->set_bb( bb );
	flex_cter_map_->set_jump( e2_k48r_jump_, false );
	flex_cter_map_->set_jump( e2_d77_jump_, false );
	for ( Size i = 0; i <= flex_cter_; i++ ) {
		flex_cter_map_->set_bb( k48r_end_ - i, true );
		flex_cter_map_->set_bb( d77_end_ - i, true );
	}

	if ( cov_bond_only_flag_ ) {
		// allow 2 additional residues of d77 C-terminal tail to move
		flex_cter_map_->set_bb( d77_end_ - ( flex_cter_ + 1 ), true );
		flex_cter_map_->set_bb( d77_end_ - ( flex_cter_ + 2 ), true );
	}

	all_dof_map_ = core::kinematics::MoveMapOP( new MoveMap() );
	all_dof_map_->clear();
	all_dof_map_->set_chi( chi );
	all_dof_map_->set_bb( bb );
	all_dof_map_->set_jump( e2_k48r_jump_, true );
	all_dof_map_->set_jump( e2_d77_jump_, true );
	for ( Size i = 0; i <= flex_cter_; i++ ) {
		all_dof_map_->set_bb( k48r_end_ - i, true );
		all_dof_map_->set_bb( d77_end_ - i, true );
	}

	init_k48r_docking_map_ = k48r_docking_map_;
	init_d77_docking_map_ = d77_docking_map_;
	init_docking_map_ = docking_map_;
	init_flex_cter_map_ = flex_cter_map_;
	init_all_dof_map_ = all_dof_map_;

	TR << "UBI Done: Setting Up Move Maps" << std::endl;

	return;
} // setup_move_maps

void
ubi_e2c_modeler::setup_complex_fold_tree(
	pose::Pose & pose_in,
	bool trim ) {

	using namespace kinematics;

	TR << "UBI Setting up complex fold tree" << std::endl;

	FoldTree f;
	f.clear();
	Size nres = d77_end_;
	Size jumppoint1 = e2_ctr_of_mass_;
	Size jumppoint2 = k48r_ctr_of_mass_;
	Size jumppoint3 = d77_ctr_of_mass_;
	Size cutpoint1 = e2_end_;
	Size cutpoint2 = k48r_end_;

	if ( trim ) {
		nres = d77_end_ - ( 2 * ( flex_cter_ + 1 ));
		jumppoint1 = e2_ctr_of_mass_;
		jumppoint2 = k48r_ctr_of_mass_;
		jumppoint3 = d77_ctr_of_mass_ - ( flex_cter_ + 1 );
		cutpoint1 = e2_end_;
		cutpoint2 = k48r_end_ - ( flex_cter_ + 1 );
	}

	f.add_edge( 1, jumppoint1, Edge::PEPTIDE );
	f.add_edge( jumppoint1, cutpoint1, Edge::PEPTIDE );
	f.add_edge( cutpoint1 + 1, jumppoint2, Edge::PEPTIDE );
	f.add_edge( jumppoint2, cutpoint2, Edge::PEPTIDE );
	f.add_edge( cutpoint2 + 1, jumppoint3, Edge::PEPTIDE );
	f.add_edge( jumppoint3, nres, Edge::PEPTIDE );
	f.add_edge( jumppoint1, jumppoint2, e2_k48r_jump_ );
	f.add_edge( jumppoint1, jumppoint3, e2_d77_jump_ );
	f.reorder( 1 );

	pose_in.fold_tree( f );

	TR << "UBI Done: Setting up complex fold tree" << std::endl;

} // setup_complex_fold_tree

void
ubi_e2c_modeler::initial_cter_perturbation(
	pose::Pose & pose_in ) {
	using namespace protocols::moves;

	TR << "UBI Initial C-Terminal Perturbation" << std::endl;

	setup_complex_fold_tree( pose_in );

	// idealize c-terminals
	for ( Size i = 0; i <= flex_cter_; i++ ) {
		conformation::idealize_position( k48r_end_ - ( flex_cter_ - i ),
			pose_in.conformation());
		conformation::idealize_position( d77_end_ - ( flex_cter_ - i),
			pose_in.conformation());
	}

	// pose_in.dump_pdb( "idealized.pdb" );

	// extend c-terminals
	Real const init_phi  ( -150.0 );
	Real const init_psi  (  150.0 );
	Real const init_omega(  180.0 );

	for ( Size i = 0; i <= flex_cter_; i++ ) {
		if ( i != flex_cter_ ) {
			pose_in.set_phi( k48r_end_ - i, init_phi );
			pose_in.set_phi( d77_end_ - i, init_phi );
		}
		if ( i != 0 ) {
			pose_in.set_psi( k48r_end_ - i, init_psi );
			pose_in.set_psi( d77_end_ - i, init_psi );
		}
		if ( ( i != flex_cter_ ) && ( i != 0 ) ) {
			pose_in.set_omega( k48r_end_ - i, init_omega );
			pose_in.set_omega( d77_end_ - i, init_omega );
		}
	}

	// pose_in.dump_pdb( "extended.pdb" );


	kinematics::MoveMapOP ub_cter_map( new kinematics::MoveMap() );
	ub_cter_map->clear();
	ub_cter_map->set_chi( true );
	ub_cter_map->set_bb( false );
	for ( Size i = 0; i <= flex_cter_; i++ ) {
		ub_cter_map->set_bb( k48r_end_ - i, true );
		ub_cter_map->set_bb( d77_end_ - i, true );
	}

	SequenceMoverOP perturb_min_cter( new SequenceMover() );

	protocols::simple_moves::BackboneMoverOP small_mover( new protocols::simple_moves::SmallMover(ub_cter_map,
		temperature_, 8 ) );
	small_mover->angle_max( 90.0 );
	protocols::simple_moves::BackboneMoverOP shear_mover( new protocols::simple_moves::ShearMover(ub_cter_map,
		temperature_, 8 ) );

	shear_mover->angle_max( 90.0 );

	protocols::simple_moves::MinMoverOP min_mover( new protocols::simple_moves::MinMover ( ub_cter_map, lowres_cst_scorefxn_,
		"linmin", min_tolerance_, nb_list_, false, false ) );

	perturb_min_cter->add_mover( small_mover );
	perturb_min_cter->add_mover( shear_mover );
	perturb_min_cter->add_mover( min_mover );

	MonteCarloOP mc;
	mc = MonteCarloOP( new moves::MonteCarlo( pose_in, *lowres_cst_scorefxn_,temperature_) );
	TrialMoverOP cter_pert_trial( new TrialMover( perturb_min_cter, mc ) );
	RepeatMoverOP cter_cycle;
	cter_cycle = RepeatMoverOP( new RepeatMover( cter_pert_trial, 40 ) ); // cycles
	cter_cycle->apply( pose_in );
	mc->recover_low( pose_in );


	TR << "UBI Done: Initial C-Terminal Perturbation" << std::endl;

} // initial_cter_perturbation

void
ubi_e2c_modeler::setup_simple_fold_tree(
	Size jumppoint1,
	Size cutpoint,
	Size jumppoint2,
	Size nres,
	pose::Pose & pose_in ) {

	using namespace kinematics;

	TR << "UBI Setting up simple fold tree" << std::endl;

	FoldTree f;
	f.clear();

	f.add_edge( 1, jumppoint1, Edge::PEPTIDE );
	f.add_edge( jumppoint1, cutpoint, Edge::PEPTIDE );
	f.add_edge( cutpoint + 1, jumppoint2, Edge::PEPTIDE );
	f.add_edge( jumppoint2, nres, Edge::PEPTIDE );
	f.add_edge( jumppoint1, jumppoint2, e2_k48r_jump_ );
	f.reorder( 1 );

	pose_in.fold_tree( f );

	TR << "UBI Done: Setting up simple fold tree" << std::endl;

} // setup_simple_fold_tree

void
ubi_e2c_modeler::trim_cter(
	pose::Pose & pose_in ) {

	TR << "UBI Trimming C-terminal" << std::endl;

	for ( Size i = 0; i <= flex_cter_; i++ ) {
		pose_in.delete_polymer_residue( d77_end_ - i );
	}

	for ( Size i = 0; i <= flex_cter_; i++ ) {
		pose_in.delete_polymer_residue( k48r_end_ - i );
	}

	k48r_trim_end_ = k48r_end_ - ( flex_cter_ + 1 );
	d77_trim_end_ = pose_in.size();
	d77_trim_ctr_mass_ = d77_ctr_of_mass_ - ( flex_cter_ + 1 );

	TR << "UBI Done: Trimming C-terminal" << std::endl;

	return;
} // trim_cter

void
ubi_e2c_modeler::restore_cter(
	pose::Pose & pose_in,
	pose::Pose without_cter ) {

	TR << "UBI Restoring C-terminal" << std::endl;

	setup_complex_fold_tree( pose_in );
	Size ub_trim_size = d77_trim_end_ - k48r_trim_end_;
	Size ub_size = d77_end_ - k48r_end_;
	pose::Pose k48r( without_cter, e2_end_ + 1, k48r_trim_end_ );
	pose::Pose d77( without_cter, k48r_trim_end_ + 1, d77_trim_end_ );

	for ( Size i = 0, j = 0; i <= flex_cter_; i++, j++ ) {
		k48r.conformation().safely_append_polymer_residue_after_seqpos(
			pose_in.residue( k48r_end_  - (flex_cter_ - i ) ), ub_trim_size +
			j, true);
		d77.conformation().safely_append_polymer_residue_after_seqpos(
			pose_in.residue( d77_end_  - (flex_cter_ - i ) ), ub_trim_size +
			j, true);
	}

	pose_in.copy_segment( ub_size, k48r, e2_end_ + 1, 1 );
	protocols::simple_moves::SwitchResidueTypeSetMover to_centroid( chemical::CENTROID );
	to_centroid.apply( pose_in );
	pose_in.copy_segment( ub_size, d77, k48r_end_ + 1, 1 );
	to_centroid.apply( pose_in );

	//pose_in.dump_pdb( "post_append.pdb" );
	setup_complex_fold_tree( pose_in );

	TR << "UBI Done: Restoring C-terminal" << std::endl;

	return;
} // restore_cter

void
ubi_e2c_modeler::init_k48r_perturbation(
	pose::Pose & pose_in ) {
	//using pose::datacache::CacheableDataType::INTERFACE_INFO;
	using namespace core::scoring;
	using namespace protocols::docking;
	using namespace protocols::moves;

	TR << "UBI Initial K48R Perturbation" << std::endl;

	pose::Pose e2_k48r( pose_in, 1, k48r_trim_end_ );

	Size ub_trim_size = d77_trim_end_ - k48r_trim_end_;

	Size jumppoint1 = e2_ctr_of_mass_;
	Size cutpoint = e2_end_;
	Size jumppoint2 = k48r_ctr_of_mass_;
	Size nres = k48r_trim_end_;
	setup_simple_fold_tree( jumppoint1, cutpoint, jumppoint2, nres, e2_k48r);

	protocols::simple_moves::SwitchResidueTypeSetMover to_centroid( chemical::CENTROID );
	to_centroid.apply( e2_k48r );

	// make starting perturbations based on command-line flags
	DockingInitialPerturbationOP init_e2_mono_ub_dock( new
		DockingInitialPerturbation( e2_k48r_jump_,true/*slide into contact*/) );
	// pose_in.dump_pdb( "pre.pdb" );
	if ( higher_d77_pert_mode_ && k48r_swap_ ) {
		rigid::RigidBodyPerturbMover mover( e2_k48r_jump_,
			15, // rot magnitude
			5 ); // trans magnitude
		mover.apply( e2_k48r );
		rigid::RigidBodySpinMover spin( e2_k48r_jump_ );
		spin.apply( e2_k48r );
		DockingSlideIntoContact slide( e2_k48r_jump_ );
		slide.apply( e2_k48r );
	} else {
		init_e2_mono_ub_dock->apply( e2_k48r );
	}

	{
		// dock movers
		rigid::RigidBodyPerturbNoCenterMoverOP dock_e2_mono_ub( new
			rigid::RigidBodyPerturbNoCenterMover( e2_k48r_jump_, 10.0, // rot_magnitude
			1.0 ) ); // trans_magnitude_
		MonteCarloOP mc;
		mc = MonteCarloOP( new moves::MonteCarlo( e2_k48r, *dock_lowres_scorefxn_,
			temperature_ ) );
		TrialMoverOP dock_trial( new TrialMover( dock_e2_mono_ub, mc ) );
		Size const cycles( 50 );
		RepeatMoverOP rbp_cycle( new RepeatMover( dock_trial, cycles ) );
		rbp_cycle->apply( e2_k48r );
		mc->recover_low( e2_k48r );
	}
	pose_in.copy_segment( ub_trim_size, e2_k48r, e2_end_ + 1, e2_end_ + 1 );
	to_centroid.apply( pose_in );

	Size e2_k48r_ctr_of_mass = core::pose::residue_center_of_mass( pose_in,
		1, k48r_trim_end_ );
	// pose_in.dump_pdb( "mid.pdb" );

	jumppoint1 = e2_k48r_ctr_of_mass;
	cutpoint = k48r_trim_end_;
	jumppoint2 = d77_trim_ctr_mass_;
	nres = d77_trim_end_;
	setup_simple_fold_tree( jumppoint1, cutpoint, jumppoint2, nres, pose_in);
	if ( higher_d77_pert_mode_ && !k48r_swap_ ) {
		rigid::RigidBodyPerturbMover mover( e2_k48r_jump_,
			15, // rot magnitude
			5 ); // trans magnitude
		mover.apply( pose_in );
		rigid::RigidBodySpinMover spin( e2_k48r_jump_ );
		spin.apply( pose_in );
		DockingSlideIntoContact slide( e2_k48r_jump_ );
		slide.apply( pose_in );
	} else {
		init_e2_mono_ub_dock->apply( pose_in );
	}
	{
		// dock movers
		rigid::RigidBodyPerturbNoCenterMoverOP dock_e2_mono_ub( new rigid::RigidBodyPerturbNoCenterMover( e2_k48r_jump_, 10.0, // rot_magnitude
			1.0 ) ); // trans_magnitude_

		setup_complex_fold_tree( pose_in, true );
		protocols::scoring::InterfaceInfoOP tri_interface( new protocols::scoring::InterfaceInfo( e2_k48r_jump_ ) );
		tri_interface->add_jump( e2_d77_jump_ );
		pose_in.data().set( core::pose::datacache::CacheableDataType::INTERFACE_INFO, tri_interface );
		Real score = ( *dock_lowres_scorefxn_ )( pose_in );
		Real current_CSP_fraction = CSP_fraction( pose_in,
			false, // non_CSP
			true, // trim
			k48r_swap_);
		Real current_non_CSP_fraction = CSP_fraction( pose_in,
			true, // non_CSP
			true, // trim
			k48r_swap_);
		Real last_accepted_CSP_score = score
			+ ( centroid_non_CSP_weight_ * current_non_CSP_fraction)
			- ( centroid_CSP_weight_ * current_CSP_fraction );
		Real lowest_CSP_score = last_accepted_CSP_score;
		setup_simple_fold_tree( jumppoint1, cutpoint, jumppoint2, nres,
			pose_in);
		protocols::scoring::InterfaceInfoOP one_interface( new protocols::scoring::InterfaceInfo( e2_k48r_jump_ ) );
		pose_in.data().set( core::pose::datacache::CacheableDataType::INTERFACE_INFO, one_interface );
		pose::Pose best_CSP_pose( pose_in );
		pose::Pose lowest_CSP_score_pose( pose_in );

		for ( Size i = 1; i <= 50 ; i++ ) {
			dock_e2_mono_ub->apply( pose_in );
			( *dock_lowres_scorefxn_ )( pose_in );
			setup_complex_fold_tree( pose_in, true );
			pose_in.data().set( core::pose::datacache::CacheableDataType::INTERFACE_INFO, tri_interface );
			score = ( *dock_lowres_scorefxn_ )( pose_in );
			current_CSP_fraction = CSP_fraction(pose_in, false, true,k48r_swap_);
			current_non_CSP_fraction=CSP_fraction(pose_in,true,true, k48r_swap_);
			Real current_CSP_score = score
				+ ( centroid_non_CSP_weight_ * current_non_CSP_fraction )
				- ( centroid_CSP_weight_ * current_CSP_fraction );
			Real boltz_factor = ( last_accepted_CSP_score - current_CSP_score ) /
				temperature_;
			Real probability = std::exp( std::min (40.0, std::max(-40.0,
				boltz_factor)));
			Real random = numeric::random::rg().uniform();
			if ( (current_CSP_score < last_accepted_CSP_score) ||
					( random < probability ) ) {
				last_accepted_CSP_score = current_CSP_score;
				best_CSP_pose = pose_in;
				if ( current_CSP_score < lowest_CSP_score ) {
					lowest_CSP_score = current_CSP_score;
					lowest_CSP_score_pose = pose_in;
				}
			}
			pose_in = best_CSP_pose;
			setup_simple_fold_tree( jumppoint1, cutpoint, jumppoint2, nres,
				pose_in);
			pose_in.data().set( core::pose::datacache::CacheableDataType::INTERFACE_INFO, one_interface );
		}
		pose_in = lowest_CSP_score_pose;
		setup_simple_fold_tree( jumppoint1, cutpoint, jumppoint2, nres,
			pose_in);
		pose_in.data().set( core::pose::datacache::CacheableDataType::INTERFACE_INFO, one_interface );
	}

	// setup_complex_fold_tree( pose_in );
	// pose_in.dump_pdb( "final.pdb" );

	TR << "UBI Done: Initial K48R Perturbation" << std::endl;

	return;
} // init_k48r_perturbation

void
ubi_e2c_modeler::init_d77_perturbation(
	pose::Pose & pose_in ) {

	TR << "UBI Initial D77 Perturbation" << std::endl;

	pose::Pose e2_k48r( pose_in, 1, k48r_trim_end_ );
	pose::Pose e2_d77( e2_k48r );
	Size ub_trim_size = d77_trim_end_ - k48r_trim_end_;
	e2_d77.copy_segment( ub_trim_size, pose_in, e2_end_ + 1,
		k48r_trim_end_ + 1 );

	protocols::simple_moves::SwitchResidueTypeSetMover to_centroid( chemical::CENTROID );
	to_centroid.apply( e2_d77 );

	pose_in.copy_segment( ub_trim_size,e2_k48r, k48r_trim_end_+1, e2_end_+1);
	to_centroid.apply( pose_in );
	pose_in.copy_segment( ub_trim_size, e2_d77, e2_end_ + 1, e2_end_+1);
	to_centroid.apply( pose_in );

	init_k48r_perturbation( pose_in );

	pose::Pose pseudo_e2_k48r( pose_in, 1, k48r_trim_end_ );
	pose::Pose pseudo_e2_d77( e2_k48r );
	pseudo_e2_d77.copy_segment( ub_trim_size, pose_in,e2_end_ + 1,
		k48r_trim_end_ + 1 );
	to_centroid.apply( pseudo_e2_d77 );

	pose_in.copy_segment( ub_trim_size, pseudo_e2_k48r, k48r_trim_end_ + 1,
		e2_end_ + 1 );
	to_centroid.apply( pose_in );
	pose_in.copy_segment( ub_trim_size,pseudo_e2_d77,e2_end_+1, e2_end_ + 1);
	to_centroid.apply( pose_in );

	TR << "UBI Done: Initial D77 Perturbation" << std::endl;

	return;
} // init_d77_perturbation

Real
ubi_e2c_modeler::initial_perturbation(
	pose::Pose & pose_in ) {
	//using pose::datacache::CacheableDataType::INTERFACE_INFO;
	using namespace core::scoring;
	using namespace protocols::docking;
	using namespace protocols::moves;

	TR << "UBI Initial Perturbation" << std::endl;

	protocols::scoring::InterfaceInfoOP one_interface( new protocols::scoring::InterfaceInfo( e2_k48r_jump_ ) );
	pose_in.data().set( core::pose::datacache::CacheableDataType::INTERFACE_INFO, one_interface );

	// initial_cter_perturbation( pose_in );

	pose::Pose without_cter( pose_in );

	trim_cter( without_cter );

	Real k48r_vs_d77 = numeric::random::rg().uniform();
	if ( k48r_vs_d77 < 0.5 ) {
		k48r_swap_ = false;
		init_k48r_perturbation( without_cter );
	} else {
		k48r_swap_ = true;
		init_d77_perturbation( without_cter );
	}

	restore_cter( pose_in, without_cter );
	initial_cter_perturbation( pose_in );

	/*
	// make starting perturbations based on command-line flags
	DockingInitialPerturbationOP init_e2_k48r_dock( new
	DockingInitialPerturbation( e2_k48r_jump_,
	true)); // slide into contact
	DockingInitialPerturbationOP init_e2_d77_dock( new
	DockingInitialPerturbation( e2_d77_jump_,
	true)); // slide into contact

	SequenceMoverOP init_e2_k48r_d77_docker( new SequenceMover() );
	init_e2_k48r_d77_docker->add_mover( init_e2_k48r_dock );
	init_e2_k48r_d77_docker->add_mover( init_e2_d77_dock );
	SequenceMoverOP init_e2_d77_k48r_docker( new SequenceMover() );
	init_e2_d77_k48r_docker->add_mover( init_e2_d77_dock );
	init_e2_d77_k48r_docker->add_mover( init_e2_k48r_dock );

	RandomMoverOP init_docker( new RandomMover() );
	init_docker->add_mover( init_e2_k48r_d77_docker, 0.5 );
	init_docker->add_mover( init_e2_d77_k48r_docker, 0.5 );

	init_docker->apply( pose_in );
	//pose_in.dump_pdb( "init_perturb.pdb" );
	// checking distance on randomizing
	*/

	Real current_ubi_cov_bond_dist( 100.00 );
	numeric::xyzVector_float d77_lys_CA(
		pose_in.residue( d77_48_lys_ ).xyz( CA ) ),
		k48r_gly_CA( pose_in.residue( k48r_end_ ).xyz( CA ) );
	current_ubi_cov_bond_dist = d77_lys_CA.distance( k48r_gly_CA );

	// Resetting Interfaces to three
	protocols::scoring::InterfaceInfoOP tri_interface( new protocols::scoring::InterfaceInfo( e2_k48r_jump_ ) );
	tri_interface->add_jump( e2_d77_jump_ );
	pose_in.data().set( core::pose::datacache::CacheableDataType::INTERFACE_INFO, tri_interface );

	TR << "UBI Done: Initial Perturbation" << std::endl;

	return( current_ubi_cov_bond_dist );
} // initial_perturbation

Real
ubi_e2c_modeler::centroid_mode_perturbation(
	pose::Pose & pose_in ) {

	//using pose::datacache::CacheableDataType::INTERFACE_INFO;
	using namespace core::scoring;
	using namespace protocols::docking;
	using namespace protocols::moves;

	TR << "UBI Centroid Mode Perturbation" << std::endl;

	setup_complex_fold_tree( pose_in );

	protocols::scoring::InterfaceInfoOP tri_interface( new protocols::scoring::InterfaceInfo( e2_k48r_jump_ ) );
	tri_interface->add_jump( e2_d77_jump_ );
	pose_in.data().set( core::pose::datacache::CacheableDataType::INTERFACE_INFO, tri_interface );

	pose::Pose start_pose;
	start_pose = pose_in;
	Real max_di_ubi_dist =  1.25 * max_k48_cter_dist_;
	Real current_ubi_cov_bond_dist( 100.00 );

	Size init_trials(0);
	while ( current_ubi_cov_bond_dist > max_di_ubi_dist ) {
		init_trials++;
		pose_in = start_pose;
		current_ubi_cov_bond_dist = initial_perturbation( pose_in );
		TR << "UBI Init Trial: " << init_trials
			<< "\t" << "Cov Bond Distance: "
			<< current_ubi_cov_bond_dist << std::endl;
	}

	//pose_in.dump_pdb( "final_init.pdb" );

	// movers
	// dock movers
	rigid::RigidBodyPerturbNoCenterMoverOP dock_e2_k48r( new rigid::RigidBodyPerturbNoCenterMover( e2_k48r_jump_, 5.0, // rot_magnitude
		0.7 ) ); // trans_magnitude_
	rigid::RigidBodyPerturbNoCenterMoverOP dock_k48r_d77( new rigid::RigidBodyPerturbNoCenterMover( e2_d77_jump_, 5.0, // rot_magnitude
		0.7 ) ); // trans_magnitude_

	RandomMoverOP docker( new RandomMover() );
	//docker->add_mover( dock_e2_k48r, 0.5 );
	//docker->add_mover( dock_k48r_d77, 0.5 );

	RepeatMoverOP cter_cycle;

	{
		//////////Small/ShearMovers//////////////////////////////////////
		SequenceMoverOP perturb_min_cter( new SequenceMover() );

		protocols::simple_moves::BackboneMoverOP small_mover( new protocols::simple_moves::SmallMover( flex_cter_map_,
			temperature_, 5 /*n_moves*/ ) );
		small_mover->angle_max( 90.0 );

		protocols::simple_moves::BackboneMoverOP shear_mover( new protocols::simple_moves::ShearMover( flex_cter_map_,
			temperature_, 5 /*n_moves*/ ) );
		shear_mover->angle_max( 90.0 );

		protocols::simple_moves::MinMoverOP min_mover( new protocols::simple_moves::MinMover(flex_cter_map_, lowres_cst_scorefxn_,
			"linmin", min_tolerance_, nb_list_, false /*deriv_check*/,
			false /* non verbose-deriv-check, default*/ ) );

		perturb_min_cter->add_mover( small_mover );
		perturb_min_cter->add_mover( shear_mover );
		perturb_min_cter->add_mover( min_mover );

		MonteCarloOP c_ter_mc;
		c_ter_mc = MonteCarloOP( new moves::MonteCarlo( pose_in, *lowres_cst_scorefxn_,
			temperature_ ) );
		TrialMoverOP cter_pert_trial( new TrialMover( perturb_min_cter,
			c_ter_mc ) );
		cter_cycle = RepeatMoverOP( new RepeatMover( cter_pert_trial, 10 ) ); // cycles
	}

	//SequenceMoverOP dock_n_perturb( new SequenceMover() );
	//dock_n_perturb->add_mover( docker );
	//dock_n_perturb->add_mover( cter_cycle );

	docker->add_mover( dock_e2_k48r, 0.375 );
	docker->add_mover( dock_k48r_d77, 0.375 );
	docker->add_mover( cter_cycle, 0.25 );

	//MonteCarloOP mc;
	//mc = new moves::MonteCarlo( pose_in, *dock_lowres_cst_scorefxn_,
	//              temperature_ );
	//
	//TrialMoverOP dock_pert_trial = new TrialMover( dock_n_perturb, mc );
	//
	//Size const cycles( 50 );
	//RepeatMoverOP rbp_cycle = new RepeatMover( dock_pert_trial, cycles );

	//rbp_cycle->apply( pose_in );
	//mc->recover_low( pose_in );

	Real score = ( *dock_lowres_cst_scorefxn_ )( pose_in );
	Real last_accepted_CSP_score = score
		+ ( centroid_non_CSP_weight_ * CSP_fraction( pose_in, true ) )
		- ( centroid_CSP_weight_ * CSP_fraction( pose_in ) );
	Real lowest_CSP_score = last_accepted_CSP_score;
	pose::Pose best_CSP_pose( pose_in );
	pose::Pose lowest_CSP_score_pose( pose_in );

	for ( Size i = 1; i <= 50 ; i++ ) {
		docker->apply( pose_in );
		score = ( *dock_lowres_cst_scorefxn_ )( pose_in );
		Real current_CSP_fraction = CSP_fraction( pose_in );
		Real current_non_CSP_fraction = CSP_fraction( pose_in, true );
		Real current_CSP_score = score
			+ ( centroid_non_CSP_weight_ * current_non_CSP_fraction )
			- ( centroid_CSP_weight_ * current_CSP_fraction );
		Real boltz_factor = ( last_accepted_CSP_score - current_CSP_score ) /
			temperature_;
		Real probability = std::exp( std::min (40.0, std::max(-40.0,
			boltz_factor)));
		Real random = numeric::random::rg().uniform();
		if ( (current_CSP_score < last_accepted_CSP_score) ||
				( random < probability ) ) {
			last_accepted_CSP_score = current_CSP_score;
			best_CSP_pose = pose_in;
			if ( current_CSP_score < lowest_CSP_score ) {
				lowest_CSP_score = current_CSP_score;
				lowest_CSP_score_pose = pose_in;
			}
		}
		pose_in = best_CSP_pose;
	}
	pose_in = lowest_CSP_score_pose;

	// pose_in.dump_pdb( "post_centroid.pdb" );
	numeric::xyzVector_float d77_lys_CA(
		pose_in.residue( d77_48_lys_ ).xyz( CA ) ),
		k48r_gly_CA( pose_in.residue( k48r_end_ ).xyz( CA ) );
	current_ubi_cov_bond_dist = d77_lys_CA.distance( k48r_gly_CA );
	TR << "UBI Done: Centroid Mode Perturbation" << std::endl;

	return( current_ubi_cov_bond_dist );
} // centroid_mode_perturbation

Real
ubi_e2c_modeler::fullatom_mode_perturbation(
	pose::Pose & pose_in ) {
	using namespace moves;
	using namespace pack::task;
	using namespace pack::task::operation;

	TR << "UBI Fullatom Mode Perturbation" << std::endl;

	initial_repack( pose_in ); // **REAL** REPACKING

	SequenceMoverOP fullatom_optimizer( new SequenceMover() ); // **MAIN**

	protocols::simple_moves::PackRotamersMoverOP pack_interface_repack( new protocols::simple_moves::PackRotamersMover(
		pack_scorefxn_ ) );
	pack_interface_repack->task_factory(tf_);

	//fullatom_optimizer->add_mover( pack_interface_repack ); // **MAIN**
	//pack_interface_repack->apply( pose_in );

	//set up minimizer movers
	protocols::simple_moves::MinMoverOP k48r_dock_min_mover( new protocols::simple_moves::MinMover( k48r_docking_map_,
		dockfa_cst_min_scorefxn_, min_type_, min_tolerance_, nb_list_ ) );
	protocols::simple_moves::MinMoverOP d77_dock_min_mover( new protocols::simple_moves::MinMover( d77_docking_map_,
		dockfa_cst_min_scorefxn_, min_type_, min_tolerance_, nb_list_ ) );
	protocols::simple_moves::MinMoverOP dock_min_mover( new protocols::simple_moves::MinMover( docking_map_,
		dockfa_cst_min_scorefxn_, min_type_, min_tolerance_, nb_list_ ) );
	protocols::simple_moves::MinMoverOP flex_cter_min_mover( new protocols::simple_moves::MinMover( flex_cter_map_,
		dockfa_cst_min_scorefxn_, min_type_, min_tolerance_, nb_list_ ) );
	protocols::simple_moves::MinMoverOP all_dof_min_mover( new protocols::simple_moves::MinMover( all_dof_map_,
		dockfa_cst_min_scorefxn_, min_type_, min_tolerance_, nb_list_ ) );

	protocols::simple_moves::RotamerTrialsMinMoverOP rtmin( new protocols::simple_moves::RotamerTrialsMinMover(
		pack_scorefxn_, tf_ ) );

	// set up rigid body movers
	Real trans_magnitude = 0.1; // default high-res docking values
	Real rot_magnitude = 5.0; // default high-res docking values
	rigid::RigidBodyPerturbMoverOP k48r_perturb( new rigid::RigidBodyPerturbMover(
		e2_k48r_jump_, rot_magnitude, trans_magnitude , rigid::partner_downstream,
		true ) );
	rigid::RigidBodyPerturbMoverOP d77_perturb( new rigid::RigidBodyPerturbMover(
		e2_d77_jump_, rot_magnitude, trans_magnitude, rigid::partner_downstream,
		true ) );

	protocols::simple_moves::RotamerTrialsMoverOP pack_rottrial( new protocols::simple_moves::RotamerTrialsMover(
		pack_scorefxn_, tf_ ) );

	MonteCarloOP mc;
	mc = MonteCarloOP( new moves::MonteCarlo( pose_in, *dockfa_cst_scorefxn_,temperature_) );

	// cter movers

	RepeatMoverOP cter_mover_rot;
	//////////Small/ShearMovers//////////////////////////////////////
	SequenceMoverOP perturb_min_cter( new SequenceMover() );

	protocols::simple_moves::BackboneMoverOP small_mover( new protocols::simple_moves::SmallMover( flex_cter_map_,
		temperature_, 5 ) );
	small_mover->angle_max( 5.0 );
	protocols::simple_moves::BackboneMoverOP shear_mover( new protocols::simple_moves::ShearMover( flex_cter_map_,
		temperature_, 5 ) );
	shear_mover->angle_max( 5.0 );

	perturb_min_cter->add_mover( small_mover );
	perturb_min_cter->add_mover( shear_mover );
	perturb_min_cter->add_mover( pack_rottrial );
	perturb_min_cter->add_mover( flex_cter_min_mover );

	TrialMoverOP cter_pert_trial( new TrialMover( perturb_min_cter, mc ) );
	cter_mover_rot = RepeatMoverOP( new RepeatMover( cter_pert_trial, 10 ) ); // cycles

	SequenceMoverOP cter_mover_repack( new SequenceMover() );
	cter_mover_repack->add_mover( cter_mover_rot );
	cter_mover_repack->add_mover( pack_interface_repack );

	////////////////////////////////////////////////////////////////////

	cter_mover_repack->apply( pose_in ); // **REAL** FIRST C-TER RELAX

	SequenceMoverOP k48r_perturb_rot( new SequenceMover() );
	SequenceMoverOP d77_perturb_rot( new SequenceMover() );
	SequenceMoverOP cter_perturb_rot( new SequenceMover() );


	k48r_perturb_rot->add_mover( k48r_perturb );
	k48r_perturb_rot->add_mover( pack_rottrial );

	d77_perturb_rot->add_mover( d77_perturb );
	d77_perturb_rot->add_mover( pack_rottrial );

	RandomMoverOP random_moves( new RandomMover() );
	random_moves->add_mover( k48r_perturb_rot, 0.375 );
	random_moves->add_mover( d77_perturb_rot, 0.375 );
	random_moves->add_mover( perturb_min_cter, 0.25 );

	//TrialMoverOP mover_trial = new TrialMover( random_moves, mc );

	SequenceMoverOP min_n_repack( new SequenceMover() );
	min_n_repack->add_mover( pack_interface_repack );

	//CycleMoverOP mover_min_trial  = new CycleMover;
	//for ( Size i=1; i < 8; ++i )
	// mover_min_trial->add_mover( mover_trial );
	//mover_min_trial->add_mover( min_n_repack );

	//RepeatMoverOP multi_cycles = new RepeatMover( mover_min_trial, 50 );

	Real score = ( *dockfa_cst_scorefxn_ )( pose_in );
	Real last_accepted_CSP_score = score
		+ ( fullatom_non_CSP_weight_ * CSP_fraction( pose_in, true ) )
		- ( fullatom_CSP_weight_ * CSP_fraction( pose_in ) );
	Real lowest_CSP_score = last_accepted_CSP_score;
	pose::Pose best_CSP_pose( pose_in );
	pose::Pose lowest_CSP_score_pose( pose_in );

	for ( Size i = 1; i <= 100 ; i++ ) {
		if ( i == 1 ) {
			perturb_min_cter->apply( pose_in ); // **REAL** FIRST MIN CTER
		} else if ( i%8 == 0 ) {
			min_n_repack->apply( pose_in ); // **REAL** OCCASIONAL REPACK
		} else {
			random_moves->apply( pose_in ); // **REAL** RANDOM MOVE SELECTION
		}
		score = ( *dockfa_cst_scorefxn_ )( pose_in );
		Real current_CSP_fraction = CSP_fraction( pose_in );
		Real current_non_CSP_fraction = CSP_fraction( pose_in, true );
		Real current_CSP_score = score
			+ ( fullatom_non_CSP_weight_ * current_non_CSP_fraction )
			- ( fullatom_CSP_weight_ * current_CSP_fraction );
		Real boltz_factor = ( last_accepted_CSP_score - current_CSP_score ) /
			temperature_;
		Real probability = std::exp( std::min (40.0, std::max(-40.0,
			boltz_factor)));
		Real random = numeric::random::rg().uniform();
		if ( ( current_CSP_score < last_accepted_CSP_score ) ||
				( random < probability ) ) {
			last_accepted_CSP_score = current_CSP_score;
			best_CSP_pose = pose_in;
			if ( current_CSP_score < lowest_CSP_score ) {
				lowest_CSP_score = current_CSP_score;
				lowest_CSP_score_pose = pose_in;
			}
		}
		pose_in = best_CSP_pose;
	}
	pose_in = lowest_CSP_score_pose;

	//pack_interface_repack->apply( pose_in );

	//fullatom_optimizer->add_mover( multi_cycles ); // **MAIN**
	fullatom_optimizer->add_mover( pack_interface_repack ); // **MAIN**

	fullatom_optimizer->apply( pose_in ); // **MAIN**
	rtmin->apply( pose_in );
	//mc->recover_low( pose_in );
	//pose_in.dump_pdb( "post_fullatom.pdb" );

	Real current_ubi_cov_bond_dist( 100.00 );

	numeric::xyzVector_float d77_lys_CA(
		pose_in.residue( d77_48_lys_ ).xyz( CA ) ),
		k48r_gly_CA( pose_in.residue( k48r_end_ ).xyz( CA ) );
	current_ubi_cov_bond_dist = d77_lys_CA.distance( k48r_gly_CA );

	applied_fullatom_pert_ = true;

	TR << "UBI Done: Fullatom Mode Perturbation" << std::endl;

	return( current_ubi_cov_bond_dist );

} // fullatom_mode_perturbation

void
ubi_e2c_modeler::initial_repack(
	pose::Pose & pose_in ) {
	using namespace moves;
	using namespace pack::task;
	using namespace pack::task::operation;

	TR << "UBI Initial Repack" << std::endl;

	moves::MonteCarloOP mc;
	mc = moves::MonteCarloOP( new moves::MonteCarlo( pose_in, *pack_scorefxn_, temperature_ ) );

	setup_packer_task( pose_in );
	// restrict_to_interfacial_loop_packing( pose_in );

	protocols::simple_moves::PackRotamersMoverOP pack_interface_repack( new protocols::simple_moves::PackRotamersMover(
		pack_scorefxn_ ) );
	pack_interface_repack->task_factory(tf_);
	TrialMoverOP pack_interface_and_loops_trial( new TrialMover(
		pack_interface_repack, mc ) );

	if ( !refinement_mode_ ) {
		pack_interface_and_loops_trial->apply( pose_in );
	}
	//pose_in.dump_pdb( "post_initial_repack" );
	TR << "UBI Done: Initial Repack" << std::endl;

	return;
} // initial_repack

void
ubi_e2c_modeler::setup_packer_task(
	pose::Pose & pose_in ) {
	// using namespace basic::options;
	using namespace pack::task;
	using namespace pack::task::operation;


	if ( init_task_factory_ ) {
		tf_ = core::pack::task::TaskFactoryOP( new TaskFactory( *init_task_factory_ ) );
		TR << "UBI Reinitializing Packer Task" << std::endl;
		return;
	} else {
		tf_ = core::pack::task::TaskFactoryOP( new TaskFactory );
	}

	TR << "UBI Setting Up Packer Task" << std::endl;

	tf_->push_back( TaskOperationCOP( new OperateOnCertainResidues( ResLvlTaskOperationOP( new PreventRepackingRLT ), ResFilterOP( new ResidueLacksProperty("PROTEIN") ) ) ) );
	tf_->push_back( TaskOperationCOP( new InitializeFromCommandline ) );
	tf_->push_back( TaskOperationCOP( new IncludeCurrent ) );
	tf_->push_back( TaskOperationCOP( new RestrictToRepacking ) );
	tf_->push_back( TaskOperationCOP( new NoRepackDisulfides ) );

	// incorporating Ian's UnboundRotamer operation.
	// note that nothing happens if unboundrot option is inactive!
	pack::rotamer_set::UnboundRotamersOperationOP unboundrot( new pack::rotamer_set::UnboundRotamersOperation() );
	unboundrot->initialize_from_command_line();
	operation::AppendRotamerSetOP unboundrot_operation( new operation::AppendRotamerSet( unboundrot ) );
	tf_->push_back( unboundrot_operation );
	// adds scoring bonuses for the "unbound" rotamers, if any
	core::pack::dunbrack::load_unboundrot( pose_in );

	init_task_factory_ = tf_;

	TR << "UBI Done: Setting Up Packer Task" << std::endl;

} // setup_packer_task

void
ubi_e2c_modeler::restrict_to_interfacial_loop_packing(
	pose::Pose & pose_in ) {
	using namespace pack::task;
	using namespace pack::task::operation;

	TR << "UBI Restricting To Interface" << std::endl;

	( *pack_scorefxn_ )( pose_in );

	// selecting movable c-terminal residues
	ObjexxFCL::FArray1D_bool loop_residues( pose_in.size(), false );
	for ( Size i = 0; i <= flex_cter_; i++ ) {
		loop_residues( k48r_end_ - i) = true;
		loop_residues( d77_end_ - i) = true;
	}

	// setting up triple interfaces
	utility::vector1_size rb_jumps;
	rb_jumps.push_back( e2_k48r_jump_ );
	rb_jumps.push_back( e2_d77_jump_ );
	tf_->push_back( TaskOperationCOP( new protocols::toolbox::task_operations::RestrictToInterface( rb_jumps, loop_residues) ) );

	TR << "UBI Done: Restricting To Interface" << std::endl;
	return;
} // restrict_to_interfacial_loop_packing

void
ubi_e2c_modeler::set_e2g2_diubi_fold_tree(
	pose::Pose & pose_in ) {

	using namespace kinematics;

	TR << "UBI Setting up E2G2-DiUbiqutin Fold Tree" << std::endl;

	FoldTree f;
	f.clear();
	Size nres = d77_end_;
	Size jumppoint1 = e2_ctr_of_mass_;
	Size diubi_center_of_mass = core::pose::residue_center_of_mass( pose_in,
		e2_end_ + 1, nres );
	Size jumppoint2 = diubi_center_of_mass;
	Size cutpoint = e2_end_;

	f.add_edge( 1, jumppoint1, Edge::PEPTIDE );
	f.add_edge( jumppoint1, cutpoint, Edge::PEPTIDE );
	f.add_edge( cutpoint + 1, jumppoint2, Edge::PEPTIDE );
	f.add_edge( jumppoint2, nres, Edge::PEPTIDE );
	f.add_edge( jumppoint1, jumppoint2, 1 );
	f.reorder( 1 );

	pose_in.fold_tree( f );

	TR << "UBI Done: Setting up E2G2-DiUbiqutin Fold Tree" << std::endl;

	return;
} //set_e2g2_diubi_fold_tree

Real
ubi_e2c_modeler::calc_interaction_energy(
	const pose::Pose & pose_in,
	bool dimer ) {
	using namespace core::scoring;

	TR << "UBI Calculating Interaction Energy" << std::endl;

	ScoreFunctionOP docking_scorefxn;
	pose::Pose complex_pose = pose_in;

	moves::SequenceMoverOP separate( new moves::SequenceMover() );
	Real trans_magnitude = 1000;
	if ( dimer ) {
		set_e2g2_diubi_fold_tree( complex_pose );
		rigid::RigidBodyTransMoverOP translate_away_diubi( new rigid::RigidBodyTransMover( complex_pose, 1 ) );
		translate_away_diubi->step_size( trans_magnitude );
		separate->add_mover( translate_away_diubi );
	} else {
		setup_complex_fold_tree( complex_pose );
		rigid::RigidBodyTransMoverOP translate_away_k48r( new rigid::RigidBodyTransMover( complex_pose, e2_k48r_jump_ ) );
		translate_away_k48r->step_size( trans_magnitude );
		rigid::RigidBodyTransMoverOP translate_away_d77( new rigid::RigidBodyTransMover( complex_pose, e2_d77_jump_ ) );
		translate_away_d77->step_size( 0.00 - trans_magnitude );
		separate->add_mover( translate_away_k48r );
		separate->add_mover( translate_away_d77 );
	}

	if ( pose_in.is_fullatom() ) {
		docking_scorefxn = ScoreFunctionFactory::create_score_function(
			"docking" );
	} else {
		docking_scorefxn = ScoreFunctionFactory::create_score_function(
			"interchain_cen" );
	}

	Real bound_energy = ( *docking_scorefxn )( complex_pose );
	separate->apply( complex_pose );
	Real unbound_energy = ( *docking_scorefxn )( complex_pose );

	if ( dimer ) {
		TR << "UBI Done: Calculating Dimer Interaction Energy" << std::endl;
	} else {
		TR << "UBI Done: Calculating Trimer Interaction Energy" << std::endl;
	}

	return (bound_energy - unbound_energy);
} // calc_interaction_energy

core::Real
ubi_e2c_modeler::CSP_fraction(
	const core::pose::Pose & pose_in,
	bool non_CSP,
	bool trim,
	bool swap ) {

	if ( non_CSP && ( ( fullatom_non_CSP_weight_ == 0.00 ) ||
			( centroid_non_CSP_weight_ == 0.00 ) ) ) {
		return( 0.00 );
	}

	/*
	bool fullatom = pose_in.is_fullatom();

	if( fullatom )
	TR << "UBI Calculating Fullatom CSP Fraction" << std::endl;
	else
	TR << "UBI Calculating Centroid CSP Fraction" << std::endl;
	*/

	utility::vector1<Size> current_CSP;
	if ( !non_CSP ) {
		if ( trim ) {
			if ( swap ) {
				current_CSP = CSP_swap_trim_;
			} else {
				current_CSP = CSP_trim_;
			}
		} else {
			current_CSP = CSP_;
		}
	} else  {
		if ( trim ) {
			if ( swap ) {
				current_CSP = non_CSP_swap_trim_;
			} else {
				current_CSP = non_CSP_trim_;
			}
		} else {
			current_CSP = non_CSP_;
		}
	}

	Size total_CSPs = current_CSP.size();

	Size nres = pose_in.size();
	utility::vector1<bool> is_interface( nres, false );

	Size num_jump = 2;
	for ( Size jj = 1; jj <= num_jump; jj++ ) {
		// "interface" is a reserved c++ keyword (on VC++ )
		protocols::scoring::Interface _interface( jj );
		// distance decided by inspection of PDB 2FUH
		_interface.distance( 8.0 );
		_interface.calculate( pose_in );

		for ( Size ii=1; ii <= nres; ++ii ) {
			if ( _interface.is_interface(ii) ) {
				is_interface[ii] = true;
			}
		}
	}

	Size satisfied_CSPs = 0;
	for ( Size ii = 1; ii <= total_CSPs; ++ii ) {
		if ( is_interface[ current_CSP[ii] ] ) {
			satisfied_CSPs++;
		}
	}

	Real CSP_fraction = ( (Real) satisfied_CSPs ) / ( (Real) total_CSPs);

	/*
	if( fullatom )
	TR << "UBI Done: Calculating Fullatom CSP Fraction " << CSP_fraction
	<< std::endl;
	else
	TR << "UBI Done: Calculating Centroid CSP Fraction " << CSP_fraction
	<< std::endl;
	*/

	return( CSP_fraction );

} // CSP_fraction

bool
ubi_e2c_modeler::centroid_filter( pose::Pose & pose_in ) {
	using namespace core::scoring;

	TR << "UBI Checking Centroid Filter" << std::endl;

	bool passed_filter = true;

	if ( score_map_[ "AD_k48_CA_CA" ] >= max_k48_cter_dist_ ) {
		passed_filter = false;
		TR << "UBI Failed Centroid Filter K48_CA_CA: "
			<< score_map_[ "AD_k48_CA_CA" ] << std::endl;
	} else {
		TR << "UBI Success Centroid Filter K48_CA_CA: "
			<< score_map_[ "AD_k48_CA_CA" ] << std::endl;
	}

	dock_lowres_cst_scorefxn_->set_weight( core::scoring::atom_pair_constraint, 1);
	( *dock_lowres_cst_scorefxn_ )( pose_in );
	Real constraint_score = pose_in.energies().total_energies()[
		atom_pair_constraint ];
	score_map_[ "AC_constraint" ] = constraint_score;
	if ( passed_filter && ( constraint_score >= fullatom_constraint_cutoff_) ) {
		passed_filter = false;
		TR << "UBI Failed Centroid Filter Constraint Score: "
			<< constraint_score
			<< std::endl;
	} else {
		if ( constraint_score < fullatom_constraint_cutoff_ ) {
			TR << "UBI Success Centroid Filter Constraint Score: "
				<< constraint_score
				<< std::endl;
		}
	}
	dock_lowres_cst_scorefxn_->set_weight( core::scoring::atom_pair_constraint,
		cen_constraint_ );
	( *dock_lowres_cst_scorefxn_ )( pose_in );

	Real CSP_fraction_ratio = CSP_fraction( pose_in );
	score_map_[ "AA_CSP_fraction" ] = CSP_fraction_ratio;
	Real non_CSP_fraction_ratio = CSP_fraction( pose_in, true );
	score_map_[ "AB_non_CSP_fraction" ] = non_CSP_fraction_ratio;

	if ( passed_filter && ( CSP_fraction_ratio <
			centroid_allowed_CSP_fraction_ ) ) {
		passed_filter = false;
		TR << "UBI Failed Centroid Filter CSP Fraction: "
			<< CSP_fraction_ratio << std::endl;
	} else {
		if ( CSP_fraction_ratio >= centroid_allowed_CSP_fraction_ ) {
			TR << "UBI Success Centroid Filter CSP Fraction: "
				<< CSP_fraction_ratio << std::endl;
		}
	}

	Real cen_score( 10000.00 );
	if ( passed_filter ) {
		cen_score = ( *output_cen_scorefxn_ )( pose_in );
	}
	if ( passed_filter && ( cen_score > 0.00 ) ) {
		passed_filter = false;
		TR << "UBI Failed Centroid Score: " << cen_score << std::endl;
	} else {
		if ( cen_score <= 0.00 ) {
			TR << "UBI Success Centroid Score: " << cen_score << std::endl;
		}
	}

	if ( !passed_filter ) {
		TR << "UBI STRUCTURE FAILED LOW-RES FILTER" << std::endl;
	} else {
		TR << "UBI Done: Successfully Passed Centroid Filter" << std::endl;
	}

	return( passed_filter );
} // centroid_filter


bool
ubi_e2c_modeler::fullatom_filter( pose::Pose & pose_in ) {
	using namespace core::scoring;

	TR << "UBI Checking Fullatom Filter" << std::endl;

	bool passed_filter = true;

	if ( passed_filter && (score_map_["AD_k48_CA_CA"] >= max_k48_cter_dist_) ) {
		passed_filter = false;
		TR << "UBI Failed Fullatom Filter K48_CA_CA: "
			<< score_map_[ "AD_k48_CA_CA" ] << std::endl;
	} else {
		TR << "UBI Success Fullatom Filter K48_CA_CA: "
			<< score_map_[ "AD_k48_CA_CA" ] << std::endl;
	}

	dockfa_cst_scorefxn_->set_weight( core::scoring::atom_pair_constraint, 1 );
	( *dockfa_cst_scorefxn_ )( pose_in );
	Real constraint_score = pose_in.energies().total_energies()[
		atom_pair_constraint ];
	if ( refinement_mode_ ) {
		constraint_score = constraint_score / 1000000;
	}
	score_map_[ "AC_constraint" ] = constraint_score;
	if ( passed_filter &&  ( constraint_score >= fullatom_constraint_cutoff_
			&& !refinement_mode_ ) ) {
		passed_filter = false;
		TR << "UBI Failed Fullatom Filter Constraint Score: "
			<< constraint_score
			<< std::endl;
	} else {
		if ( constraint_score < fullatom_constraint_cutoff_ ) {
			TR << "UBI Success Fullatom Filter Constraint Score: "
				<< constraint_score
				<< std::endl;
		}
	}
	dockfa_cst_scorefxn_->set_weight( core::scoring::atom_pair_constraint,
		full_constraint_ );
	( *dockfa_cst_scorefxn_ )( pose_in );

	Real CSP_fraction_ratio = CSP_fraction( pose_in );
	score_map_[ "AA_CSP_fraction" ] = CSP_fraction_ratio;
	Real non_CSP_fraction_ratio = CSP_fraction( pose_in, true  );
	score_map_[ "AB_non_CSP_fraction" ] = non_CSP_fraction_ratio;
	if ( passed_filter && ( CSP_fraction_ratio <
			fullatom_allowed_CSP_fraction_ ) ) {
		passed_filter = false;
		TR << "UBI Failed Fullatom Filter CSP Fraction: "
			<< CSP_fraction_ratio << std::endl;
	} else {
		if ( CSP_fraction_ratio >= fullatom_allowed_CSP_fraction_ ) {
			TR << "UBI Success Fullatom Filter CSP Fraction: "
				<< CSP_fraction_ratio << std::endl;
		}
	}

	score_map_[ "AH_I_sc_trimer" ] = 10.00;
	//if( passed_filter )
	score_map_[ "AH_I_sc_trimer" ] = calc_interaction_energy( pose_in,
		false ); //dimer

	score_map_[ "AG_I_sc_e2_ubi" ] = calc_interaction_energy( pose_in );

	/*
	if ( passed_filter && score_map_[ "A_I_sc" ] >= 0.0) {
	passed_filter = false;
	TR << "UBI Failed Fullatom Filter I_sc: "
	<< score_map_[ "A_I_sc" ] << std::endl;
	}
	else {
	if( score_map_[ "A_I_sc" ] < 0.0)
	TR << "UBI Success Fullatom Filter I_sc: "
	<< score_map_ [ "A_I_sc" ] << std::endl;
	}
	*/

	if ( !passed_filter ) {
		TR << "UBI STRUCTURE FAILED HIGH-RES FILTER" << std::endl;
	} else {
		TR << "UBI Done: Successfully Passed Fullatom Filter" << std::endl;
	}

	return passed_filter;
} // fullatom_filter

core::Real
ubi_e2c_modeler::calc_Lrmsd (
	const pose::Pose & pose_in,
	const pose::Pose & native_pose,
	Size ubiquitin ) {

	using namespace core::scoring;

	ObjexxFCL::FArray1D_bool superpos_partner ( pose_in.size(), false );

	Size compute_rmsd_start(0), compute_rmsd_end(0);

	if ( ubiquitin == e2_k48r_jump_ ) {
		compute_rmsd_start = e2_end_ + 1;
		compute_rmsd_end = k48r_end_;
	} else if ( ubiquitin == e2_d77_jump_ ) {
		compute_rmsd_start = k48r_end_ + 1;
		compute_rmsd_end = d77_end_;
	} else {
		compute_rmsd_start = k48r_end_ - flex_cter_;
		compute_rmsd_end = k48r_end_;
	}

	for ( Size i = compute_rmsd_start; i <= compute_rmsd_end; ++i ) {
		superpos_partner(i) = true;
	}

	Real Lrmsd = rmsd_no_super_subset( native_pose, pose_in,
		superpos_partner, is_protein_CA );
	return ( Lrmsd );
} // calc_Lrmsd


void
ubi_e2c_modeler::evaluate_native( pose::Pose & pose_in ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace chemical;
	using namespace id;
	//using namespace jobdist;
	//using pose::datacache::CacheableDataType::INTERFACE_INFO;
	using utility::file::FileName;
	using namespace moves;
	using namespace pack::task;
	using namespace pack::task::operation;
	using namespace kinematics;

	// score functions
	core::scoring::ScoreFunctionOP dock_lowres_score;
	core::scoring::ScoreFunctionOP lowres_score;
	core::scoring::ScoreFunctionOP pack_score;
	core::scoring::ScoreFunctionOP dockfa_score;

	dock_lowres_score = core::scoring::ScoreFunctionFactory::
		create_score_function("interchain_cen");
	lowres_score = core::scoring::ScoreFunctionFactory::
		create_score_function( "cen_std" );
	dockfa_score = core::scoring::ScoreFunctionFactory::
		create_score_function( "docking"  );
	pack_score = core::scoring::get_score_function("pre_talaris_2013_standard.wts");

	setup_key_residues( pose_in );
	setup_complex_fold_tree( pose_in );

	protocols::scoring::InterfaceInfoOP tri_interface( new protocols::scoring::InterfaceInfo( e2_k48r_jump_ ) );
	tri_interface->add_jump( e2_d77_jump_ );
	pose_in.data().set( core::pose::datacache::CacheableDataType::INTERFACE_INFO, tri_interface );

	pose::Pose start_pose;
	start_pose = pose_in;

	protocols::simple_moves::SwitchResidueTypeSetMover to_centroid( chemical::CENTROID );
	protocols::simple_moves::SwitchResidueTypeSetMover to_all_atom( chemical::FA_STANDARD );
	protocols::simple_moves::ReturnSidechainMover recover_sidechains( start_pose);

	Real score( 0.00 );
	to_centroid.apply( pose_in );

	score = ( *lowres_score )( pose_in );
	TR << "Native Cen Score     : " << score << std::endl;

	score = ( *dock_lowres_score )( pose_in );
	TR << "Native DockLow Score : " << score << std::endl;

	to_all_atom.apply( pose_in );
	recover_sidechains.apply( pose_in );

	score = ( *dockfa_score )( pose_in );
	TR << "Native DockFA Score  : " << score << std::endl;

	score = ( *pack_score )( pose_in );
	TR << "Native Pack Score    : " << score << std::endl;

	setup_packer_task( pose_in );
	protocols::simple_moves::PackRotamersMoverOP repack( new protocols::simple_moves::PackRotamersMover( pack_score ) );
	repack->task_factory( tf_ );
	repack->apply( pose_in );

	score = ( *dockfa_score )( pose_in );
	TR << "Packed DockFA Score  : " << score << std::endl;

	score = ( *pack_score )( pose_in );
	TR << "Packed Pack Score    : " << score << std::endl;

	set_e2g2_diubi_fold_tree( start_pose );

	Real trans_magnitude = 1000;
	rigid::RigidBodyTransMoverOP translate_away_diubi( new rigid::RigidBodyTransMover( start_pose, 1 ) );
	translate_away_diubi->step_size( trans_magnitude );

	Real bound_energy = ( *dockfa_score )( start_pose );
	translate_away_diubi->apply( start_pose );
	Real unbound_energy = ( *dockfa_score )( start_pose );

	TR << "Interaction Score    : " << bound_energy - unbound_energy
		<< std::endl;

	return;
} // evaluate_native

void
ubi_e2c_modeler::optimize_cov_bond(
	pose::Pose & pose_in ) {
	using namespace kinematics;
	//using pose::datacache::CacheableDataType::INTERFACE_INFO;
	using namespace core::scoring;
	using namespace protocols::docking;
	using namespace protocols::moves;

	TR << "UBI Optimizing Covalent Bond" << std::endl;

	Real const local_full_constraint( 100000.00 );
	dockfa_cst_scorefxn_->set_weight( atom_pair_constraint,
		local_full_constraint );
	dockfa_cst_min_scorefxn_->set_weight( atom_pair_constraint,
		local_full_constraint );
	pack_cst_scorefxn_->set_weight( atom_pair_constraint,
		local_full_constraint );


	setup_complex_fold_tree( pose_in );

	protocols::scoring::InterfaceInfoOP tri_interface( new protocols::scoring::InterfaceInfo( e2_k48r_jump_ ) );
	tri_interface->add_jump( e2_d77_jump_ );
	pose_in.data().set( core::pose::datacache::CacheableDataType::INTERFACE_INFO, tri_interface );

	setup_packer_task( pose_in );

	protocols::simple_moves::PackRotamersMoverOP packer( new protocols::simple_moves::PackRotamersMover( pack_cst_scorefxn_ ) );
	packer->task_factory(tf_);

	//set up minimizer
	protocols::simple_moves::MinMoverOP flex_cter_min_mover( new protocols::simple_moves::MinMover( flex_cter_map_,
		dockfa_cst_min_scorefxn_, min_type_, min_tolerance_, nb_list_ ) );

	protocols::simple_moves::RotamerTrialsMoverOP pack_rottrial( new protocols::simple_moves::RotamerTrialsMover(
		pack_cst_scorefxn_, tf_ ) );

	flex_cter_min_mover->apply( pose_in ); // **REAL** MINIMIZE C TER

	// cter movers
	SequenceMoverOP perturb_min_cter( new SequenceMover() );

	protocols::simple_moves::BackboneMoverOP small_mover( new protocols::simple_moves::SmallMover( flex_cter_map_,
		temperature_, 5 ) );
	small_mover->angle_max( 5.0 );
	protocols::simple_moves::BackboneMoverOP shear_mover( new protocols::simple_moves::ShearMover( flex_cter_map_,
		temperature_, 5 ) );
	shear_mover->angle_max( 5.0 );

	perturb_min_cter->add_mover( small_mover );
	perturb_min_cter->add_mover( shear_mover );
	perturb_min_cter->add_mover( pack_rottrial );
	perturb_min_cter->add_mover( flex_cter_min_mover );

	Real score = ( *dockfa_cst_scorefxn_ )( pose_in );
	Real last_accepted_CSP_score = score
		+ ( fullatom_non_CSP_weight_ * CSP_fraction( pose_in, true ) )
		- ( fullatom_CSP_weight_ * CSP_fraction( pose_in ) );
	Real lowest_CSP_score = last_accepted_CSP_score;
	pose::Pose best_CSP_pose( pose_in );
	pose::Pose lowest_CSP_score_pose( pose_in );

	for ( Size i = 1; i <= 100 ; i++ ) {
		if ( i%8 == 0 ) {
			packer->apply( pose_in );
		} else {
			perturb_min_cter->apply( pose_in ); // **REAL** C-TER MIN
		}
		score = ( *dockfa_cst_scorefxn_ )( pose_in );
		Real current_CSP_fraction = CSP_fraction( pose_in );
		Real current_non_CSP_fraction = CSP_fraction( pose_in, true );
		Real current_CSP_score = score +
			( fullatom_non_CSP_weight_ * current_non_CSP_fraction ) -
			( fullatom_CSP_weight_ * current_CSP_fraction );
		Real boltz_factor = ( last_accepted_CSP_score - current_CSP_score ) /
			temperature_;
		Real probability = std::exp( std::min (40.0, std::max(-40.0,
			boltz_factor)));
		Real random = numeric::random::rg().uniform();
		if ( ( current_CSP_score < last_accepted_CSP_score ) ||
				( random < probability ) ) {
			last_accepted_CSP_score = current_CSP_score;
			best_CSP_pose = pose_in;
			if ( current_CSP_score < lowest_CSP_score ) {
				lowest_CSP_score = current_CSP_score;
				lowest_CSP_score_pose = pose_in;
			}
		}
		pose_in = best_CSP_pose;
	}
	pose_in = lowest_CSP_score_pose;

	packer->apply( pose_in ); // **MAIN**

	Real ubi_cov_bond_dist( 100.00 );

	numeric::xyzVector_float d77_lys_CA(
		pose_in.residue( d77_48_lys_ ).xyz( CA ) ),
		k48r_gly_CA( pose_in.residue( k48r_end_ ).xyz( CA ) );

	ubi_cov_bond_dist = d77_lys_CA.distance( k48r_gly_CA );

	score_map_[ "AD_k48_CA_CA" ] = ubi_cov_bond_dist;

	Size const NZ_atom (9);
	Size const O_atom (4);

	numeric::xyzVector_float d77_lys_NZ(
		pose_in.residue( d77_48_lys_ ).xyz( NZ_atom ) ),
		k48r_gly_O( pose_in.residue( k48r_end_ ).xyz( O_atom ) );
	ubi_cov_bond_dist = d77_lys_NZ.distance( k48r_gly_O );

	score_map_[ "AE_k48_O_NZ" ] = ubi_cov_bond_dist;

	applied_fullatom_pert_ = true;

	// Restoring full atom constraint weights
	dockfa_cst_scorefxn_->set_weight( atom_pair_constraint,
		full_constraint_ );
	dockfa_cst_min_scorefxn_->set_weight( atom_pair_constraint,
		full_constraint_ );
	pack_cst_scorefxn_->set_weight( atom_pair_constraint,
		full_constraint_ );

	TR << "UBI Done: Optimizing Covalent Bond" << std::endl;

	return;

} // optimize_cov_bond

void ubi_e2c_modeler::monoub_assign_CSPs(
	const pose::Pose & pose_in ) {

	TR << "UBI Mono Ubi Assigning CSPs" << std::endl;

	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 23 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 24 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 25 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 26 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 29 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 32 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 33 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 34 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 36 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 40 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 41 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 43 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 52 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 53 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 54 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 57 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 74 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 17 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 18 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 27 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 38 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 39 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 42 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 46 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 50 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 58 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 90 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 160 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 162 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 163 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'A', 164 ) );

	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 7 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 13 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 14 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 32 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 42 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 45 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 47 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 48 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 49 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 68 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 69 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 70 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 71 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 75 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 6 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 8 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 11 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 15 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 29 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 34 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 43 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 51 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 65 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 66 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 67 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 72 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 73 ) );
	CSP_.push_back( pose_in.pdb_info()->pdb2pose( 'B', 74 ) );


	TR << "UBI Done: Mono Ubi Assigning CSPs" << std::endl;

	return;
} // monoub_assign_CSPs

void ubi_e2c_modeler::monoub_apply( pose::Pose & pose_in ) {
	using namespace basic::datacache;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace chemical;
	using namespace id;
	//using namespace jobdist;
	//using pose::datacache::CacheableDataType::SCORE_MAP;
	using utility::file::FileName;

	TR << "UBI Mono Ubi Apply Start" << std::endl;

	protocols::simple_moves::ConstraintSetMoverOP mtsl_constraint( new protocols::simple_moves::ConstraintSetMover() );
	mtsl_constraint->apply( pose_in );

	const pose::Pose start_pose( pose_in );

	// setup key residues
	monoub_setup_key_residues( start_pose );

	// assign Chemical Shift Perturbations
	monoub_assign_CSPs( start_pose );

	// setup Move Maps
	monoub_setup_move_maps();

	// Residue movers
	protocols::simple_moves::SwitchResidueTypeSetMover to_centroid( chemical::CENTROID );
	protocols::simple_moves::SwitchResidueTypeSetMover to_all_atom( chemical::FA_STANDARD );
	protocols::simple_moves::ReturnSidechainMover recover_sidechains( start_pose);

	// centroid mode
	Size tries = 0;

	//start loop of decoy creation until filters are all passed
	for ( Size r = 1; r <= max_repeats_; r++ ) {
		tries = r;
		pose_in = start_pose;

		// convert to centroid mode
		to_centroid.apply( pose_in );
		monoub_centroid_mode_perturbation( pose_in );
		protocols::jd2::ScoreMap::nonzero_energies( score_map_, lowres_cst_scorefxn_, pose_in);

		// check low-res docking filter
		( *dock_lowres_cst_scorefxn_ )( pose_in );
		passed_centroid_filter_ = monoub_centroid_filter( pose_in );

		applied_fullatom_pert_ = false;

		// fullatom mode
		if ( passed_centroid_filter_ || r == max_repeats_ ) {
			// convert to full atom
			to_all_atom.apply( pose_in );
			recover_sidechains.apply( pose_in );

			monoub_fullatom_mode_perturbation( pose_in );

			// add scores to map for output

			protocols::jd2::ScoreMap::nonzero_energies( score_map_, pack_cst_scorefxn_, pose_in);

			// check highres docking filter
			( *dockfa_cst_scorefxn_ )( pose_in );
			passed_fullatom_filter_ = monoub_fullatom_filter( pose_in );
		} // if fullatom mode

		if ( passed_centroid_filter_ && passed_fullatom_filter_ ) break;
		else {
			TR <<"UBI Mono Ubi Repeating structure " << r << " times"
				<< std::endl;
		}
	} // for max_repeats_

	// add scores to map for output
	if ( applied_fullatom_pert_ ) {
		( *output_full_scorefxn_ )( pose_in );
		protocols::jd2::ScoreMap::nonzero_energies(score_map_, output_full_scorefxn_, pose_in);
	} else {
		( *output_cen_scorefxn_ )( pose_in );
		protocols::jd2::ScoreMap::nonzero_energies( score_map_, output_cen_scorefxn_, pose_in);
	}

	score_map_["AJ_monoub_rms"] = monoub_calc_Lrmsd( pose_in, start_pose );

	pose_in.data().set( core::pose::datacache::CacheableDataType::SCORE_MAP,
		DataCache_CacheableData::DataOP( new basic::datacache::DiagnosticData(score_map_) ));

	TR << "UBI Mono Ubi Outputing structure after " << tries << " times"
		<< std::endl;

}// end monoub_apply

void
ubi_e2c_modeler::monoub_setup_key_residues(
	const pose::Pose & pose_in ) {

	TR << "UBI Mono Ubi Setting Up Key Residues" << std::endl;

	pose::PDBInfoCOP pdb_info = pose_in.pdb_info();
	monoub_end_ = pose_in.size();

	//char chain = '_';
	char old_chain = '_';
	for ( Size i = 1; i <= monoub_end_; i++ ) {
		char chain = pdb_info->chain( i );

		// if initial condition
		if ( i == 1 ) {
			old_chain = pdb_info->chain( i );
		}

		if ( chain != old_chain ) {
			e2_end_ = i - 1;
			break;
		}
		old_chain = chain;
	}// for i <= monoub_end_

	e2_ctr_of_mass_ = core::pose::residue_center_of_mass( pose_in, 1,
		e2_end_ );
	monoub_ctr_of_mass_ = core::pose::residue_center_of_mass( pose_in,
		e2_end_ + 1,
		monoub_end_ );

	TR << "UBI Done: Mono Ubi Setting Up Key Residues" << std::endl;

} // monoub_setup_key_residues

void
ubi_e2c_modeler::monoub_setup_move_maps() {
	using namespace core;
	using namespace kinematics;

	if ( init_all_dof_map_ ) {
		monoub_all_dof_map_ = core::kinematics::MoveMapOP( new MoveMap( *init_monoub_all_dof_map_ ) );
		monoub_docking_map_ = core::kinematics::MoveMapOP( new MoveMap( *init_monoub_docking_map_ ) );
		monoub_flex_cter_map_ = core::kinematics::MoveMapOP( new MoveMap( *init_monoub_flex_cter_map_ ) );

		TR << "UBI Mono Ubi Reinitializing Move Maps" << std::endl;
		return;
	}

	TR << "UBI Mono Ubi Setting Up Move Maps" << std::endl;

	bool bb = false;
	bool chi = true;

	monoub_docking_map_ = core::kinematics::MoveMapOP( new MoveMap() );
	monoub_docking_map_->clear();
	monoub_docking_map_->set_chi( chi );
	monoub_docking_map_->set_bb( bb );
	monoub_docking_map_->set_jump( 1, true );

	monoub_flex_cter_map_ = core::kinematics::MoveMapOP( new MoveMap() );
	monoub_flex_cter_map_->clear();
	monoub_flex_cter_map_->set_chi( chi );
	monoub_flex_cter_map_->set_bb( bb );
	monoub_flex_cter_map_->set_jump( 1, false );
	for ( Size i = 0; i <= flex_cter_; i++ ) {
		monoub_flex_cter_map_->set_bb( monoub_end_ - i, true );
	}

	monoub_all_dof_map_ = core::kinematics::MoveMapOP( new MoveMap() );
	monoub_all_dof_map_->clear();
	monoub_all_dof_map_->set_chi( chi );
	monoub_all_dof_map_->set_bb( bb );
	monoub_all_dof_map_->set_jump( 1, true );
	for ( Size i = 0; i <= flex_cter_; i++ ) {
		monoub_all_dof_map_->set_bb( monoub_end_ - i, true );
	}

	init_monoub_docking_map_ = monoub_docking_map_;
	init_monoub_flex_cter_map_ = monoub_flex_cter_map_;
	init_monoub_all_dof_map_ = monoub_all_dof_map_;

	TR << "UBI Done: Mono Ubi Setting Up Move Maps" << std::endl;

	return;
} // monoub_setup_move_maps

void
ubi_e2c_modeler::monoub_fold_tree(
	pose::Pose & pose_in ) {

	using namespace kinematics;

	TR << "UBI Setting up mono ubi fold tree" << std::endl;

	FoldTree f;
	f.clear();
	Size nres = monoub_end_;
	Size jumppoint1 = e2_ctr_of_mass_;
	Size jumppoint2 = monoub_ctr_of_mass_;
	Size cutpoint1 = e2_end_;

	f.add_edge( 1, jumppoint1, Edge::PEPTIDE );
	f.add_edge( jumppoint1, cutpoint1, Edge::PEPTIDE );
	f.add_edge( cutpoint1 + 1, jumppoint2, Edge::PEPTIDE );
	f.add_edge( jumppoint2, nres, Edge::PEPTIDE );
	f.add_edge( jumppoint1, jumppoint2, 1 );
	f.reorder( 1 );

	pose_in.fold_tree( f );

	TR << "UBI Done: Setting up mono ubi fold tree" << std::endl;

} // monoub_fold_tree

void
ubi_e2c_modeler::monoub_initial_cter_perturbation(
	pose::Pose & pose_in ) {
	using namespace protocols::moves;

	TR << "UBI Mono Ubi Initial C-Terminal Perturbation" << std::endl;

	// idealize c-terminals
	for ( Size i = 0; i <= flex_cter_; i++ ) {
		conformation::idealize_position( monoub_end_ - ( flex_cter_ - i ),
			pose_in.conformation());
	}

	// pose_in.dump_pdb( "idealized.pdb" );

	// extend c-terminals
	Real const init_phi  ( -150.0 );
	Real const init_psi  (  150.0 );
	Real const init_omega(  180.0 );

	for ( Size i = 0; i <= flex_cter_; i++ ) {
		if ( i != flex_cter_ ) {
			pose_in.set_phi( monoub_end_ - i, init_phi );
		}
		if ( i != 0 ) {
			pose_in.set_psi( monoub_end_ - i, init_psi );
		}
		if ( ( i != flex_cter_ ) && ( i != 0 ) ) {
			pose_in.set_omega( monoub_end_ - i, init_omega );
		}
	}

	// pose_in.dump_pdb( "extended.pdb" );


	kinematics::MoveMapOP ub_cter_map( new kinematics::MoveMap() );
	ub_cter_map->clear();
	ub_cter_map->set_chi( true );
	ub_cter_map->set_bb( false );
	for ( Size i = 0; i <= flex_cter_; i++ ) {
		ub_cter_map->set_bb( monoub_end_ - i, true );
	}

	SequenceMoverOP perturb_min_cter( new SequenceMover() );

	protocols::simple_moves::BackboneMoverOP small_mover( new protocols::simple_moves::SmallMover( ub_cter_map,
		temperature_, 8 ) );
	small_mover->angle_max( 90.0 );
	protocols::simple_moves::BackboneMoverOP shear_mover( new protocols::simple_moves::ShearMover( ub_cter_map,
		temperature_, 8 ) );

	shear_mover->angle_max( 90.0 );

	protocols::simple_moves::MinMoverOP min_mover( new protocols::simple_moves::MinMover ( ub_cter_map, lowres_cst_scorefxn_,
		"linmin", min_tolerance_, nb_list_, false, false ) );

	perturb_min_cter->add_mover( small_mover );
	perturb_min_cter->add_mover( shear_mover );
	perturb_min_cter->add_mover( min_mover );

	MonteCarloOP mc;
	mc = MonteCarloOP( new moves::MonteCarlo( pose_in, *lowres_cst_scorefxn_,temperature_) );
	TrialMoverOP cter_pert_trial( new TrialMover( perturb_min_cter, mc ) );
	RepeatMoverOP cter_cycle;
	cter_cycle = RepeatMoverOP( new RepeatMover( cter_pert_trial, 40 ) ); // cycles
	cter_cycle->apply( pose_in );
	mc->recover_low( pose_in );

	TR << "UBI Done: Mono Ubi Initial C-Terminal Perturbation" << std::endl;

} // monoub_initial_cter_perturbation

void
ubi_e2c_modeler::monoub_first_perturbation(
	pose::Pose & pose_in ) {
	//using pose::datacache::CacheableDataType::INTERFACE_INFO;
	using namespace core::scoring;
	using namespace protocols::docking;
	using namespace protocols::moves;

	TR << "UBI Mono Ubi First Perturbation" << std::endl;

	protocols::simple_moves::SwitchResidueTypeSetMover to_centroid( chemical::CENTROID );
	to_centroid.apply( pose_in );

	// make starting perturbations based on command-line flags
	DockingInitialPerturbationOP init_e2_mono_ub_dock( new
		DockingInitialPerturbation( 1 ,true/*slide into contact*/) );
	// pose_in.dump_pdb( "pre.pdb" );
	init_e2_mono_ub_dock->apply( pose_in );

	// dock movers
	rigid::RigidBodyPerturbNoCenterMoverOP dock_e2_mono_ub( new
		rigid::RigidBodyPerturbNoCenterMover( 1 , 10.0, // rot_magnitude
		1.0 ) ); // trans_magnitude_

	protocols::scoring::InterfaceInfoOP one_interface( new protocols::scoring::InterfaceInfo( 1 ) );
	pose_in.data().set( core::pose::datacache::CacheableDataType::INTERFACE_INFO, one_interface );
	Real score = ( *dock_lowres_scorefxn_ )( pose_in );
	Real current_CSP_fraction = monoub_CSP_fraction( pose_in );
	Real last_accepted_CSP_score = score
		- ( centroid_CSP_weight_ * current_CSP_fraction );
	Real lowest_CSP_score = last_accepted_CSP_score;
	pose::Pose best_CSP_pose( pose_in );
	pose::Pose lowest_CSP_score_pose( pose_in );

	for ( Size i = 1; i <= 50 ; i++ ) {
		dock_e2_mono_ub->apply( pose_in );
		( *dock_lowres_scorefxn_ )( pose_in );
		score = ( *dock_lowres_scorefxn_ )( pose_in );
		current_CSP_fraction = monoub_CSP_fraction( pose_in );
		Real current_CSP_score = score
			- ( centroid_CSP_weight_ * current_CSP_fraction );
		Real boltz_factor = ( last_accepted_CSP_score - current_CSP_score ) /
			temperature_;
		Real probability = std::exp( std::min (40.0, std::max(-40.0,
			boltz_factor)));
		Real random = numeric::random::rg().uniform();
		if ( (current_CSP_score < last_accepted_CSP_score) ||
				( random < probability ) ) {
			last_accepted_CSP_score = current_CSP_score;
			best_CSP_pose = pose_in;
			if ( current_CSP_score < lowest_CSP_score ) {
				lowest_CSP_score = current_CSP_score;
				lowest_CSP_score_pose = pose_in;
			}
		}
		pose_in = best_CSP_pose;
	}
	pose_in = lowest_CSP_score_pose;

	TR << "UBI Done: Mono Ubi First Perturbation" << std::endl;
	return;
} // monoub_first_perturbation

void
ubi_e2c_modeler::monoub_initial_perturbation(
	pose::Pose & pose_in ) {
	//using pose::datacache::CacheableDataType::INTERFACE_INFO;
	using namespace core::scoring;
	using namespace protocols::docking;
	using namespace protocols::moves;

	TR << "UBI Mono Ubi Initial Perturbation" << std::endl;

	protocols::scoring::InterfaceInfoOP one_interface( new protocols::scoring::InterfaceInfo( 1 ) );
	pose_in.data().set( core::pose::datacache::CacheableDataType::INTERFACE_INFO, one_interface );

	monoub_first_perturbation( pose_in );
	monoub_initial_cter_perturbation( pose_in );

	TR << "UBI Done: Mono Ubi Initial Perturbation" << std::endl;

	return;
} // monoub_initial_perturbation

void
ubi_e2c_modeler::monoub_centroid_mode_perturbation(
	pose::Pose & pose_in ) {

	//using pose::datacache::CacheableDataType::INTERFACE_INFO;
	using namespace core::scoring;
	using namespace protocols::docking;
	using namespace protocols::moves;

	TR << "UBI Centroid Mode Perturbation" << std::endl;

	monoub_fold_tree( pose_in );

	protocols::scoring::InterfaceInfoOP interface( new protocols::scoring::InterfaceInfo( 1 ) );
	pose_in.data().set( core::pose::datacache::CacheableDataType::INTERFACE_INFO, interface );

	monoub_initial_perturbation( pose_in );

	//pose_in.dump_pdb( "final_init.pdb" );

	// movers
	// dock movers
	rigid::RigidBodyPerturbNoCenterMoverOP dock_e2_monoub( new
		rigid::RigidBodyPerturbNoCenterMover( 1, 5.0, // rot_magnitude
		0.7 ) ); // trans_magnitude_

	RandomMoverOP docker( new RandomMover() );

	RepeatMoverOP cter_cycle;

	{
		//////////Small/ShearMovers//////////////////////////////////////
		SequenceMoverOP perturb_min_cter( new SequenceMover() );

		protocols::simple_moves::BackboneMoverOP small_mover( new protocols::simple_moves::SmallMover( monoub_flex_cter_map_,
			temperature_, 5 /*n_moves*/ ) );
		small_mover->angle_max( 90.0 );

		protocols::simple_moves::BackboneMoverOP shear_mover( new protocols::simple_moves::ShearMover( monoub_flex_cter_map_,
			temperature_, 5 /*n_moves*/ ) );
		shear_mover->angle_max( 90.0 );

		protocols::simple_moves::MinMoverOP min_mover( new protocols::simple_moves::MinMover( monoub_flex_cter_map_,
			lowres_cst_scorefxn_, "linmin", min_tolerance_, nb_list_,
			false /*deriv_check*/,false /* non verbose-deriv-check,default*/) );

		perturb_min_cter->add_mover( small_mover );
		perturb_min_cter->add_mover( shear_mover );
		perturb_min_cter->add_mover( min_mover );

		MonteCarloOP c_ter_mc;
		c_ter_mc = MonteCarloOP( new moves::MonteCarlo( pose_in, *lowres_cst_scorefxn_,
			temperature_ ) );
		TrialMoverOP cter_pert_trial( new TrialMover( perturb_min_cter,
			c_ter_mc ) );
		cter_cycle = RepeatMoverOP( new RepeatMover( cter_pert_trial, 10 ) ); // cycles
	}

	docker->add_mover( dock_e2_monoub, 0.75 );
	docker->add_mover( cter_cycle, 0.25 );

	Real score = ( *dock_lowres_cst_scorefxn_ )( pose_in );
	Real last_accepted_CSP_score = score
		- ( centroid_CSP_weight_ * monoub_CSP_fraction( pose_in ) );
	Real lowest_CSP_score = last_accepted_CSP_score;
	pose::Pose best_CSP_pose( pose_in );
	pose::Pose lowest_CSP_score_pose( pose_in );

	for ( Size i = 1; i <= 50 ; i++ ) {
		docker->apply( pose_in );
		score = ( *dock_lowres_cst_scorefxn_ )( pose_in );
		Real current_CSP_fraction = monoub_CSP_fraction( pose_in );
		Real current_CSP_score = score
			- ( centroid_CSP_weight_ * current_CSP_fraction );
		Real boltz_factor = ( last_accepted_CSP_score - current_CSP_score ) /
			temperature_;
		Real probability = std::exp( std::min (40.0, std::max(-40.0,
			boltz_factor)));
		Real random = numeric::random::rg().uniform();
		if ( (current_CSP_score < last_accepted_CSP_score) ||
				( random < probability ) ) {
			last_accepted_CSP_score = current_CSP_score;
			best_CSP_pose = pose_in;
			if ( current_CSP_score < lowest_CSP_score ) {
				lowest_CSP_score = current_CSP_score;
				lowest_CSP_score_pose = pose_in;
			}
		}
		pose_in = best_CSP_pose;
	}
	pose_in = lowest_CSP_score_pose;

	TR << "UBI Done: Mono UbiCentroid Mode Perturbation" << std::endl;

	return;
} // monoub_centroid_mode_perturbation

void
ubi_e2c_modeler::monoub_fullatom_mode_perturbation(
	pose::Pose & pose_in ) {
	using namespace moves;
	using namespace pack::task;
	using namespace pack::task::operation;

	TR << "UBI Mono Ubi Fullatom Mode Perturbation" << std::endl;

	initial_repack( pose_in ); // **REAL** REPACKING

	SequenceMoverOP fullatom_optimizer( new SequenceMover() ); // **MAIN**

	protocols::simple_moves::PackRotamersMoverOP pack_interface_repack( new protocols::simple_moves::PackRotamersMover(
		pack_scorefxn_ ) );
	pack_interface_repack->task_factory(tf_);

	//set up minimizer movers
	protocols::simple_moves::MinMoverOP monoub_dock_min_mover( new protocols::simple_moves::MinMover( monoub_docking_map_,
		dockfa_cst_min_scorefxn_, min_type_, min_tolerance_, nb_list_ ) );
	protocols::simple_moves::MinMoverOP flex_cter_min_mover( new protocols::simple_moves::MinMover( monoub_flex_cter_map_,
		dockfa_cst_min_scorefxn_, min_type_, min_tolerance_, nb_list_ ) );
	protocols::simple_moves::MinMoverOP all_dof_min_mover( new protocols::simple_moves::MinMover( monoub_all_dof_map_,
		dockfa_cst_min_scorefxn_, min_type_, min_tolerance_, nb_list_ ) );

	protocols::simple_moves::RotamerTrialsMinMoverOP rtmin( new protocols::simple_moves::RotamerTrialsMinMover(
		pack_scorefxn_, tf_ ) );

	// set up rigid body movers
	Real trans_magnitude = 0.1; // default high-res docking values
	Real rot_magnitude = 5.0; // default high-res docking values
	rigid::RigidBodyPerturbMoverOP monoub_perturb( new rigid::RigidBodyPerturbMover(
		1, rot_magnitude, trans_magnitude , rigid::partner_downstream,
		true ) );
	protocols::simple_moves::RotamerTrialsMoverOP pack_rottrial( new protocols::simple_moves::RotamerTrialsMover(
		pack_scorefxn_, tf_ ) );

	MonteCarloOP mc;
	mc = MonteCarloOP( new moves::MonteCarlo( pose_in, *dockfa_cst_scorefxn_,temperature_) );

	// cter movers

	RepeatMoverOP cter_mover_rot;
	//////////Small/ShearMovers//////////////////////////////////////
	SequenceMoverOP perturb_min_cter( new SequenceMover() );

	protocols::simple_moves::BackboneMoverOP small_mover( new protocols::simple_moves::SmallMover( monoub_flex_cter_map_,
		temperature_, 5 ) );
	small_mover->angle_max( 5.0 );
	protocols::simple_moves::BackboneMoverOP shear_mover( new protocols::simple_moves::ShearMover( monoub_flex_cter_map_,
		temperature_, 5 ) );
	shear_mover->angle_max( 5.0 );

	perturb_min_cter->add_mover( small_mover );
	perturb_min_cter->add_mover( shear_mover );
	perturb_min_cter->add_mover( pack_rottrial );
	perturb_min_cter->add_mover( flex_cter_min_mover );

	TrialMoverOP cter_pert_trial( new TrialMover( perturb_min_cter, mc ) );
	cter_mover_rot = RepeatMoverOP( new RepeatMover( cter_pert_trial, 10 ) ); // cycles

	SequenceMoverOP cter_mover_repack( new SequenceMover() );
	cter_mover_repack->add_mover( cter_mover_rot );
	cter_mover_repack->add_mover( pack_interface_repack );

	////////////////////////////////////////////////////////////////////

	cter_mover_repack->apply( pose_in ); // **REAL** FIRST C-TER RELAX

	SequenceMoverOP monoub_perturb_rot( new SequenceMover() );
	monoub_perturb_rot->add_mover( monoub_perturb );
	monoub_perturb_rot->add_mover( pack_rottrial );

	RandomMoverOP random_moves( new RandomMover() );
	random_moves->add_mover( monoub_perturb_rot, 0.75 );
	random_moves->add_mover( perturb_min_cter, 0.25 );

	SequenceMoverOP min_n_repack( new SequenceMover() );
	min_n_repack->add_mover( pack_interface_repack );

	Real score = ( *dockfa_cst_scorefxn_ )( pose_in );
	Real last_accepted_CSP_score = score
		- ( fullatom_CSP_weight_ * monoub_CSP_fraction( pose_in ) );
	Real lowest_CSP_score = last_accepted_CSP_score;
	pose::Pose best_CSP_pose( pose_in );
	pose::Pose lowest_CSP_score_pose( pose_in );

	for ( Size i = 1; i <= 100 ; i++ ) {
		if ( i == 1 ) {
			perturb_min_cter->apply( pose_in ); // **REAL** FIRST MIN CTER
		} else if ( i%8 == 0 ) {
			min_n_repack->apply( pose_in ); // **REAL** OCCASIONAL REPACK
		} else {
			random_moves->apply( pose_in ); // **REAL** RANDOM MOVE SELECTION
		}
		score = ( *dockfa_cst_scorefxn_ )( pose_in );
		Real current_CSP_fraction = monoub_CSP_fraction( pose_in );
		Real current_CSP_score = score
			- ( fullatom_CSP_weight_ * current_CSP_fraction );
		Real boltz_factor = ( last_accepted_CSP_score - current_CSP_score ) /
			temperature_;
		Real probability = std::exp( std::min (40.0, std::max(-40.0,
			boltz_factor)));
		Real random = numeric::random::rg().uniform();
		if ( ( current_CSP_score < last_accepted_CSP_score ) ||
				( random < probability ) ) {
			last_accepted_CSP_score = current_CSP_score;
			best_CSP_pose = pose_in;
			if ( current_CSP_score < lowest_CSP_score ) {
				lowest_CSP_score = current_CSP_score;
				lowest_CSP_score_pose = pose_in;
			}
		}
		pose_in = best_CSP_pose;
	}
	pose_in = lowest_CSP_score_pose;

	fullatom_optimizer->add_mover( pack_interface_repack ); // **MAIN**
	fullatom_optimizer->apply( pose_in ); // **MAIN**
	rtmin->apply( pose_in );

	applied_fullatom_pert_ = true;

	TR << "UBI Done: Mono Ubi Fullatom Mode Perturbation" << std::endl;

	return;
} // monoub_fullatom_mode_perturbation

Real
ubi_e2c_modeler::monoub_calc_interaction_energy(
	const pose::Pose & pose_in ) {
	using namespace core::scoring;

	TR << "UBI Mono Ubi Calculating Interaction Energy" << std::endl;

	ScoreFunctionOP docking_scorefxn;
	pose::Pose complex_pose = pose_in;

	moves::SequenceMoverOP separate( new moves::SequenceMover() );
	Real trans_magnitude = 1000;

	monoub_fold_tree( complex_pose );
	rigid::RigidBodyTransMoverOP translate_away_ubi( new rigid::RigidBodyTransMover( complex_pose, 1 ) );
	translate_away_ubi->step_size( trans_magnitude );
	separate->add_mover( translate_away_ubi );

	if ( pose_in.is_fullatom() ) {
		docking_scorefxn = ScoreFunctionFactory::create_score_function(
			"docking" );
	} else {
		docking_scorefxn = ScoreFunctionFactory::create_score_function(
			"interchain_cen" );
	}

	Real bound_energy = ( *docking_scorefxn )( complex_pose );
	separate->apply( complex_pose );
	Real unbound_energy = ( *docking_scorefxn )( complex_pose );

	TR << "UBI Done: Calculating Mono Ubi Interaction Energy" << std::endl;

	return (bound_energy - unbound_energy);
} // monoub_calc_interaction_energy

core::Real
ubi_e2c_modeler::monoub_CSP_fraction(
	const core::pose::Pose & pose_in ) {

	/*
	if( non_CSP && ( ( fullatom_non_CSP_weight_ == 0.00 ) ||
	( centroid_non_CSP_weight_ == 0.00 ) ) )
	return( 0.00 );
	*/

	utility::vector1<Size> current_CSP;
	current_CSP = CSP_;
	Size total_CSPs = current_CSP.size();

	Size nres = pose_in.size();
	utility::vector1<bool> is_interface( nres, false );

	protocols::scoring::Interface _interface( 1 );
	// distance decided by inspection of PDB 2FUH
	_interface.distance( 8.0 );
	_interface.calculate( pose_in );

	for ( Size ii=1; ii <= nres; ++ii ) {
		if ( _interface.is_interface(ii) ) {
			is_interface[ii] = true;
		}
	}

	Size satisfied_CSPs = 0;
	for ( Size ii = 1; ii <= total_CSPs; ++ii ) {
		if ( is_interface[ current_CSP[ii] ] ) {
			satisfied_CSPs++;
		}
	}

	Real CSP_fraction = ( (Real) satisfied_CSPs ) / ( (Real) total_CSPs);

	return( CSP_fraction );

} // monoub_CSP_fraction

bool
ubi_e2c_modeler::monoub_centroid_filter( pose::Pose & pose_in ) {
	using namespace core::scoring;

	TR << "UBI Mono Ubi Checking Centroid Filter" << std::endl;

	bool passed_filter = true;

	dock_lowres_cst_scorefxn_->set_weight( core::scoring::atom_pair_constraint, 1);
	( *dock_lowres_cst_scorefxn_ )( pose_in );
	Real constraint_score = pose_in.energies().total_energies()[
		atom_pair_constraint ];
	score_map_[ "AC_constraint" ] = constraint_score;
	if ( passed_filter && ( constraint_score >= fullatom_constraint_cutoff_) ) {
		passed_filter = false;
		TR << "UBI Failed Centroid Filter Constraint Score: "
			<< constraint_score
			<< std::endl;
	} else {
		if ( constraint_score < fullatom_constraint_cutoff_ ) {
			TR << "UBI Success Centroid Filter Constraint Score: "
				<< constraint_score
				<< std::endl;
		}
	}
	dock_lowres_cst_scorefxn_->set_weight( core::scoring::atom_pair_constraint,
		cen_constraint_ );
	( *dock_lowres_cst_scorefxn_ )( pose_in );

	Real CSP_fraction_ratio = monoub_CSP_fraction( pose_in );
	score_map_[ "AA_CSP_fraction" ] = CSP_fraction_ratio;

	if ( passed_filter && ( CSP_fraction_ratio <
			centroid_allowed_CSP_fraction_ ) ) {
		passed_filter = false;
		TR << "UBI Failed Centroid Filter CSP Fraction: "
			<< CSP_fraction_ratio << std::endl;
	} else {
		if ( CSP_fraction_ratio >= centroid_allowed_CSP_fraction_ ) {
			TR << "UBI Success Centroid Filter CSP Fraction: "
				<< CSP_fraction_ratio << std::endl;
		}
	}

	Real cen_score( 10000.00 );
	if ( passed_filter ) {
		cen_score = ( *output_cen_scorefxn_ )( pose_in );
	}
	if ( passed_filter && ( cen_score > 0.00 ) ) {
		passed_filter = false;
		TR << "UBI Failed Centroid Score: " << cen_score << std::endl;
	} else {
		if ( cen_score <= 0.00 ) {
			TR << "UBI Success Centroid Score: " << cen_score << std::endl;
		}
	}

	if ( !passed_filter ) {
		TR << "UBI MONO UBI STRUCTURE FAILED LOW-RES FILTER" << std::endl;
	} else {
		TR << "UBI Done: Mono Ubi Successfully Passed Centroid Filter" << std::endl;
	}

	return( passed_filter );
} // monoub_centroid_filter


bool
ubi_e2c_modeler::monoub_fullatom_filter( pose::Pose & pose_in ) {
	using namespace core::scoring;

	TR << "UBI Mono Ubi Checking Fullatom Filter" << std::endl;

	bool passed_filter = true;

	dockfa_cst_scorefxn_->set_weight( core::scoring::atom_pair_constraint, 1 );
	( *dockfa_cst_scorefxn_ )( pose_in );
	Real constraint_score = pose_in.energies().total_energies()[
		atom_pair_constraint ];
	score_map_[ "AC_constraint" ] = constraint_score;
	if ( passed_filter &&  ( constraint_score >= fullatom_constraint_cutoff_
			&& !refinement_mode_ ) ) {
		passed_filter = false;
		TR << "UBI Failed Fullatom Filter Constraint Score: "
			<< constraint_score
			<< std::endl;
	} else {
		if ( constraint_score < fullatom_constraint_cutoff_ ) {
			TR << "UBI Success Fullatom Filter Constraint Score: "
				<< constraint_score
				<< std::endl;
		}
	}
	dockfa_cst_scorefxn_->set_weight( core::scoring::atom_pair_constraint,
		full_constraint_ );
	( *dockfa_cst_scorefxn_ )( pose_in );

	Real CSP_fraction_ratio = monoub_CSP_fraction( pose_in );
	score_map_[ "AA_CSP_fraction" ] = CSP_fraction_ratio;
	if ( passed_filter && ( CSP_fraction_ratio <
			fullatom_allowed_CSP_fraction_ ) ) {
		passed_filter = false;
		TR << "UBI Failed Fullatom Filter CSP Fraction: "
			<< CSP_fraction_ratio << std::endl;
	} else {
		if ( CSP_fraction_ratio >= fullatom_allowed_CSP_fraction_ ) {
			TR << "UBI Success Fullatom Filter CSP Fraction: "
				<< CSP_fraction_ratio << std::endl;
		}
	}

	score_map_[ "AG_I_sc" ] = monoub_calc_interaction_energy( pose_in );

	if ( !passed_filter ) {
		TR << "UBI MONO UBI STRUCTURE FAILED HIGH-RES FILTER" << std::endl;
	} else {
		TR << "UBI Done: Mono Ubi Successfully Passed Fullatom Filter"
			<< std::endl;
	}

	return( passed_filter );
} // monoub_fullatom_filter

core::Real
ubi_e2c_modeler::monoub_calc_Lrmsd (
	const pose::Pose & pose_in,
	const pose::Pose & native_pose ) {

	using namespace core::scoring;

	ObjexxFCL::FArray1D_bool superpos_partner ( monoub_end_, false );

	Size compute_rmsd_start(0), compute_rmsd_end(0);

	compute_rmsd_start = e2_end_ + 1;
	compute_rmsd_end = monoub_end_;
	for ( Size i = compute_rmsd_start; i <= compute_rmsd_end; ++i ) {
		superpos_partner(i) = true;
	}

	Real Lrmsd = rmsd_no_super_subset( native_pose, pose_in,
		superpos_partner, is_protein_CA );
	return ( Lrmsd );
} // monoub_calc_Lrmsd
#endif
} // end ub_e2c

} // end protocols
