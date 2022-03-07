// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/glycan_docking/GlycanDockProtocol.cc
/// @brief Local docking refinement of an oligomeric glycoligand in a binding site
/// @author Morgan Nance (morganlnance@gmail.com)
/// Replaces dock_glycans app by Jason Labonte (JWLabonte@jhu.edu)


// Package Headers
#include <protocols/glycan_docking/GlycanDockProtocol.hh>
#include <protocols/glycan_docking/GlycanDockProtocolCreator.hh>
#include <protocols/glycan_docking/util.hh>

// Project Headers
#include <protocols/constraint_movers/ConstraintSetMover.hh>
#include <protocols/minimization_packing/MinMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/analysis/InterfaceAnalyzerMover.hh>
#include <protocols/carbohydrates/RingPlaneFlipMover.hh>
#include <protocols/docking/DockingInitialPerturbation.hh>
#include <protocols/docking/metrics.hh>
#include <protocols/docking/util.hh>
#include <protocols/rigid/RigidBodyMover.hh>

#include <core/simple_metrics/metrics/SelectedResiduesPyMOLMetric.hh>

#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/BBDihedralSamplerMover.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/minimization_packing/RotamerTrialsMover.hh>
//#include <protocols/moves/PyMOLMover.hh>

#include <core/scoring/constraints/util.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>

#include <core/pack/task/TaskFactory.hh>


#include <core/kinematics/MoveMap.fwd.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/types.hh>

#include <core/select/util.hh>
#include <core/select/residue_selector/ChainSelector.hh>
#include <core/select/residue_selector/FalseResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <protocols/residue_selectors/HBondSelector.hh>

// Basic Headers
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/carbohydrates.OptionKeys.gen.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/rings.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>

// Numeric Headers
#include <numeric/random/random.hh>

// Utility Headers
#include <utility/excn/Exceptions.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

// XSD Includes
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

#include <core/kinematics/FoldTree.hh> // AUTO IWYU For operator<<, FoldTree
#include <core/pack/task/PackerTask.hh> // AUTO IWYU For PackerTask


static basic::Tracer TR( "protocols.glycan_docking.GlycanDockProtocol" );

namespace protocols {
namespace glycan_docking {


////////////////////////////////////////////////////////////////////////////////
/// @brief Default constructor
GlycanDockProtocol::GlycanDockProtocol():
	protocols::moves::Mover( GlycanDockProtocol::mover_name() )
{
	type("GlycanDockProtocol");

	register_options();
	set_defaults();
}


////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
GlycanDockProtocol::GlycanDockProtocol( GlycanDockProtocol const & ) = default;


////////////////////////////////////////////////////////////////////////////////
// Default deconstructor
// (important for properly forward-declaring smart-pointer members)
GlycanDockProtocol::~GlycanDockProtocol(){}


////////////////////////////////////////////////////////////////////////////////
/// @brief Register options for Movers relevant to GlycanDockProtocol
void
GlycanDockProtocol::register_options()
{
	using namespace basic::options;

	// Register relevant options for this GlycanDock protocol
	option.add_relevant( OptionKeys::constraints::cst_fa_file ); // no carb centroids yet
	option.add_relevant( OptionKeys::docking::partners ); // required for docking
	option.add_relevant( OptionKeys::in::file::s );
	option.add_relevant( OptionKeys::in::file::native );
	option.add_relevant( OptionKeys::rings::lock_rings );
	option.add_relevant( OptionKeys::run::n_cycles );
	option.add_relevant( OptionKeys::run::min_type );
	option.add_relevant( OptionKeys::run::min_tolerance );

	// Options that are GlycanDockProtocol specific
	option.add_relevant( OptionKeys::carbohydrates::glycan_dock::n_repeats );
	option.add_relevant( OptionKeys::carbohydrates::glycan_dock::refine_only );
	option.add_relevant( OptionKeys::carbohydrates::glycan_dock::prepack_only );
	option.add_relevant( OptionKeys::carbohydrates::glycan_dock::stage1_perturb_glycan_com_rot_mag );
	option.add_relevant( OptionKeys::carbohydrates::glycan_dock::stage1_perturb_glycan_com_trans_mag );
	option.add_relevant( OptionKeys::carbohydrates::glycan_dock::stage1_rotate_glycan_about_com );
	option.add_relevant( OptionKeys::carbohydrates::glycan_dock::stage1_torsion_uniform_pert_mag );
	option.add_relevant( OptionKeys::carbohydrates::glycan_dock::n_rigid_body_rounds );
	option.add_relevant( OptionKeys::carbohydrates::glycan_dock::slide_glycan_into_contact );
	option.add_relevant( OptionKeys::carbohydrates::glycan_dock::stage2_rot_mag );
	option.add_relevant( OptionKeys::carbohydrates::glycan_dock::stage2_trans_mag );
	option.add_relevant( OptionKeys::carbohydrates::glycan_dock::n_torsion_rounds );
	option.add_relevant( OptionKeys::carbohydrates::glycan_dock::mc_kt );
	option.add_relevant( OptionKeys::carbohydrates::glycan_dock::interface_packing_distance );
	option.add_relevant( OptionKeys::carbohydrates::glycan_dock::full_packing_frequency );
	option.add_relevant( OptionKeys::carbohydrates::glycan_dock::ramp_scorefxn );
	option.add_relevant( OptionKeys::carbohydrates::glycan_dock::watch_in_pymol );

	// Call register_options() on other Movers used by this protocol
	protocols::carbohydrates::RingPlaneFlipMover::register_options();
	protocols::constraint_movers::ConstraintSetMover::register_options();
	protocols::docking::FaDockingSlideIntoContact::register_options();
	protocols::minimization_packing::MinMover::register_options();
	protocols::minimization_packing::PackRotamersMover::register_options();
	//protocols::moves::PyMOLMover::register_options();
	protocols::moves::RandomMover::register_options();
	protocols::rigid::RigidBodyRandomizeMover::register_options();
	protocols::rigid::RigidBodyPerturbMover::register_options();
	protocols::simple_moves::BackboneMover::register_options();
	protocols::simple_moves::BBDihedralSamplerMover::register_options();
	protocols::simple_moves::ShearMover::register_options();

} // END register_options


////////////////////////////////////////////////////////////////////////////////
void
GlycanDockProtocol::set_defaults()
{
	using namespace basic::options;

	// First set protocol-level default behaviors set by cmd-line
	// Will also set any behaviors passed by user's flags
	set_cmd_line_defaults();

	// dock_jump_num_ = 1 will define the Jump between the
	// upstream protein and the downstream glycoligand
	movable_jump_.push_back( dock_jump_num_ );

	clear_counters();

} // END set_defaults


////////////////////////////////////////////////////////////////////////////////
void
GlycanDockProtocol::set_cmd_line_defaults()
{
	using namespace basic::options;

	// Default behavior of a Rosetta protocol and specifically GlycanDockProtocol
	// These options already have defaults and must be grabbed for the
	// code to work, but the user could have also set an explicit value
	// Either way, grab the inputs for each protocol-necessary flag default
	n_cycles_ = option[ OptionKeys::run::n_cycles ]; // default = 10 in .hh
	if ( option[ OptionKeys::rings::lock_rings ].user() ) {
		lock_rings_ = option[ OptionKeys::rings::lock_rings ]; // default = true in .hh
	}
	min_type_ = option[ OptionKeys::run::min_type ]; // default = lbfgs_armijo_nonmonotone
	min_tolerance_ = option[ OptionKeys::run::min_tolerance ]; // default here = 0.01

	// GlycanDock specific flags
	// defaults are set in options_rosetta.py
	refine_only_ =
		option[OptionKeys::carbohydrates::glycan_dock::refine_only](); // default = false
	prepack_only_ =
		option[OptionKeys::carbohydrates::glycan_dock::prepack_only](); // default = false

	partners_provided_ = option[ OptionKeys::docking::partners ].user();
	partners_ = option[ OptionKeys::docking::partners ]();

	// For Stage 1 GlycanDock trajectory initialization
	stage1_rot_mag_ =
		option[OptionKeys::carbohydrates::glycan_dock::stage1_perturb_glycan_com_rot_mag](); // default = 7.5 degrees
	stage1_trans_mag_ =
		option[OptionKeys::carbohydrates::glycan_dock::stage1_perturb_glycan_com_trans_mag](); // default = 0.5 Angstrom
	stage1_rotate_glyc_about_com_ =
		option[OptionKeys::carbohydrates::glycan_dock::stage1_rotate_glycan_about_com](); // default = false
	stage1_torsion_uniform_pert_mag_ =
		option[OptionKeys::carbohydrates::glycan_dock::stage1_torsion_uniform_pert_mag](); // default = 12.5 degrees

	// For GlycanDock Stage 2 docking and optimization
	mc_kt_ =
		option[OptionKeys::carbohydrates::glycan_dock::mc_kt](); // default = 0.6
	n_repeats_ =
		option[OptionKeys::carbohydrates::glycan_dock::n_repeats](); // default = 1
	n_rb_rounds_ =
		option[OptionKeys::carbohydrates::glycan_dock::n_rigid_body_rounds]; // default = 8
	slide_glyc_into_contact_ =
		option[OptionKeys::carbohydrates::glycan_dock::slide_glycan_into_contact]; // default = True
	stage2_rot_mag_ =
		option[OptionKeys::carbohydrates::glycan_dock::stage2_rot_mag](); // default = 7.5
	stage2_trans_mag_ =
		option[OptionKeys::carbohydrates::glycan_dock::stage2_trans_mag](); // default = 0.5
	rand_glyc_jump_res_ =
		option[OptionKeys::carbohydrates::glycan_dock::rand_glycan_jump_res](); // default = true
	n_tor_rounds_ =
		option[OptionKeys::carbohydrates::glycan_dock::n_torsion_rounds]; // default = 8

	// Misc protocol-changing behaviors
	intf_pack_dist_ =
		option[OptionKeys::carbohydrates::glycan_dock::interface_packing_distance](); // default = 16
	full_pack_freq_ =
		option[OptionKeys::carbohydrates::glycan_dock::full_packing_frequency](); // default = 8
	ramp_scorefxn_ =
		option[OptionKeys::carbohydrates::glycan_dock::ramp_scorefxn](); // default = true
	watch_in_pymol_ =
		option[OptionKeys::carbohydrates::glycan_dock::watch_in_pymol](); // default = false

} // END set_cmd_line_defaults


////////////////////////////////////////////////////////////////////////////////
// PyRosetta-friendly version
void
GlycanDockProtocol::show() const
{
	show( TR );
}


////////////////////////////////////////////////////////////////////////////////
void
GlycanDockProtocol::show( std::ostream & output ) const
{
	protocols::moves::Mover::show( output );  // name, type, tag

	// ( bool ) ? <if true> : <if false>
	std::string const locked( ( lock_rings_ ) ? "locked" : "flexible" );
	// If refine_only_ is true, then NOT performing initial randomization
	std::string const refine_only( ( ! refine_only_ ) ? "yes" : "no" );
	std::string const prepack_only( ( prepack_only_ ) ? "yes" : "no" );
	std::string const stage1_rotate_glyc_about_com
		( ( stage1_rotate_glyc_about_com_ ) ? "yes" : "no" );
	std::string const slide_glycan_into_contact
		( ( slide_glyc_into_contact_ ) ? "yes" : "no" );
	std::string const full_packing_frequency
		( utility::to_string( full_pack_freq_ ) );
	std::string const ramp_scorefxn( ( ramp_scorefxn_ ) ? "yes" : "no" );
	std::string const watch_in_pymol( ( watch_in_pymol_ ) ? "yes" : "no" );

	// Flags relevant to any usage of the GlycanDock app
	output <<
		"Rings: " << locked << std::endl;

	// If user requests to only prepack the input structure,
	// inform them and return as other options now not relevant
	if ( prepack_only_ ) {
		output <<
			"Performing pre-packing of input structure ONLY" << std::endl <<
			"-Separating complex for packing" << std::endl;
		return;
	}

	// Options relevant to GlycanDock (with modification for refine_only_ flag)
	// Skip information here if not performing Stage 1 conformation initialization
	if ( ! refine_only_ ) {
		// Stage 1 glycosidic torsion angles
		output <<
			"--Perturbation magnitude for uniform torsion randomization "
			"(degrees): " << stage1_torsion_uniform_pert_mag_ << std::endl;
		// Stage 1 rigid-body
		output <<
			"-Rotate glycoligand about CoM during Stage 1?: " <<
			stage1_rotate_glyc_about_com << std::endl <<
			"-Stage 1 perturb CoM rotation mag. (degrees): " <<
			stage1_rot_mag_ << std::endl <<
			"-Stage 1 perturb CoM translation mag. (Ang.): " <<
			stage1_trans_mag_ << std::endl;
	}

	// The rest of the information is about Stage 2 of GlycanDock
	// n_repeats_ == 0 is for generating a conformation in Stage 1
	// only and not refining it in Stage 2
	// useful for generating benchmark structures quickly
	if ( n_repeats_ > 0 ) {
		output <<
			"During Stage 2 docking refinement" << std::endl <<
			"-Docking protein-glycoligand complex using docking partners: " <<
			partners_ << std::endl <<
			"-Cycles: " << n_cycles_ << std::endl <<
			"-Repeats (if decoy fails filter): " << n_repeats_ << std::endl <<
			"-Number of rigid-body moves per cycle: " << n_rb_rounds_ << std::endl <<
			"--Occasionally slide glycoligand into contact with protein?: " <<
			slide_glycan_into_contact << std::endl <<
			"-Number of torsion moves per cycle: " << n_tor_rounds_ << std::endl <<
			"-Distance in Ang. that defines interface residues: " <<
			intf_pack_dist_ << std::endl <<
			"-Full packing frequency: " << full_packing_frequency << std::endl <<
			"-Ramp ScoreFunction: " << ramp_scorefxn << std::endl <<
			"-Watch in PyMOL: " << watch_in_pymol << std::endl;
	}

	output << std::endl;

} // END show


////////////////////////////////////////////////////////////////////////////////
void
GlycanDockProtocol::clear_counters()
{
	TR.Info << "Clearing counters" << std::endl;

	n_rb_moves_made_ = 0;
	n_rb_moves_accepted_ = 0;
	n_rb_cycles_performed_ = 0;
	n_tor_moves_made_ = 0;
	n_tor_moves_accepted_ = 0;
	n_tor_cycles_performed_ = 0;
	n_stage2_cycles_performed_ = 0;

} // END clear_counters


////////////////////////////////////////////////////////////////////////////////
void
GlycanDockProtocol::setup_scorefxn( core::pose::Pose & pose )
{
	using namespace core::scoring;
	using namespace basic::options;

	if ( ! scorefxn_ ) {
		TR.Debug << "Setting up ScoreFunction" << std::endl;

		// Create a standard full-atom ScoreFunction
		// GlycanDock was benchmarked using ref2015 as the base sf
		scorefxn_ = get_score_function();

		// Set some terms used by carbohydrates to default weights for GlycanDock
		// sugar_bb weight as suggested by JAB
		// These terms/weights were used during benchmarking of GlycanDock
		scorefxn_->set_weight( sugar_bb, 0.5 );
		scorefxn_->set_weight( fa_intra_rep_nonprotein,
			scorefxn_->get_weight( fa_rep ) );

		TR.Debug << "Finished setting up ScoreFunction" << std::endl;
	}

	// If user set weights via the command-line, respect them
	// could override sugar_bb and fa_intra_rep_nonprotein set above
	if ( option[ OptionKeys::score::set_weights ].user() ) {
		apply_set_weights( scorefxn_, option[ OptionKeys::score::set_weights ]() );
	}

	// For score ramping purposes during docking
	// Keep track of the original weights for fa_atr and fa_rep
	target_atr_ = scorefxn_->get_weight( fa_atr );
	target_rep_ = scorefxn_->get_weight( fa_rep );

	// Setup any cmd-line full-atom constraints
	// Apply any cmd-line constraints (full-atom respected only)
	if ( option[ OptionKeys::constraints::cst_fa_file ].user() ) {
		TR << "Setting constraint(s) to Pose and ScoreFunction" << std::endl;
		// Initialization begins with selecting the first cst_fa_file given
		protocols::constraint_movers::ConstraintSetMoverOP constraint_setter =
			utility::pointer::make_shared
			< protocols::constraint_movers::ConstraintSetMover >();
		// If given multiple -cst_fa_files, pick one at random to use
		// Calling this here (instead at construction) ensures
		// each decoy uses a random cst file per apply
		constraint_setter->constraint_file
			( core::scoring::constraints::get_cst_fa_file_option() );
		// Read chosen constraint file and apply fa constraints to the pose
		constraint_setter->apply( pose );
		// Add all constraint terms to sf using weights set via -cst_fa_weight
		// Default value for -cst_fa_weight is 1.0
		core::scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn
			( *scorefxn_ );
	}

} // END setup_scorefxn


////////////////////////////////////////////////////////////////////////////////
void
GlycanDockProtocol::setup_docking_foldtree
( core::pose::Pose & pose,
	core::pose::Pose & ref_pose )
{
	TR << "Setting up protein-glycoligand docking FoldTree" << std::endl;

	// setup an appropriate protein-glycoligand docking FoldTree
	TR << "Original FoldTree for Pose:" << std::endl;
	TR << pose.fold_tree() << std::endl;
	protocols::docking::setup_foldtree
		( pose, partners_,
		movable_jump_ /* 1 */, rand_glyc_jump_res_ /* true */ );
	TR << "New docking FoldTree for Pose:" << std::endl;
	TR << pose.fold_tree() << std::endl;

	// Setup an equivalent FoldTree for the ref_pose
	// Note, it may not be the exact same FoldTree! but that is okay
	// just needs to have Jump 1 define the correct docking partners
	TR.Info << "Setting up equivalent FoldTree for reference Pose" << std::endl;
	protocols::docking::setup_foldtree
		( ref_pose, partners_,
		movable_jump_, rand_glyc_jump_res_ );

	TR.Debug << "Finished setting up protein-glycoligand docking FoldTree" << std::endl;

} // END setup_docking_foldtree


////////////////////////////////////////////////////////////////////////////////
// If ScoreTypes other than fa_atr and fa_rep are desired for ramping,
// this function must be altered accordingly
void
GlycanDockProtocol::ramp_score_weight
( core::scoring::ScoreType const score_type,
	core::Real const target_weight,
	core::Size const current_cycle,
	core::Size const n_cycles )
{
	using namespace core;
	using namespace core::scoring;

	// Set score_type weight according to the fraction completion of protocol

	if ( ( current_cycle == 1 ) && ( n_cycles != 1 ) ) {
		////////////////////
		// IF FIRST CYCLE //
		////////////////////
		// Manually checking ScoreTypes because it is assumed
		// that only fa_atr and fa_rep are being ramped
		if ( score_type == fa_atr ) {
			scorefxn_->set_weight( fa_atr, target_atr_ * starting_ramp_down_atr_factor_ );
		} else if ( score_type == fa_rep ) {
			scorefxn_->set_weight( fa_rep, target_rep_ * starting_ramp_up_rep_factor_ );
		} else {
			// catcher for if an unexpected score_type was given
			utility_exit_with_message("During score ramping, an unexpected "
				"ScoreType (" + utility::to_string(score_type)
				+ ") was given. GlycanDock was written "
				"to only ramp fa_atr and fa_rep.");
		}
	} else if ( current_cycle == n_cycles ) {
		///////////////////
		// IF LAST CYCLE //
		///////////////////
		// Sanity check, ensure score weights are set to their target values
		scorefxn_->set_weight( fa_atr, target_atr_ );
		scorefxn_->set_weight( fa_rep, target_rep_ );
	} else {
		////////////////////////
		// IF DURING PROTOCOL //
		////////////////////////
		// Main ramping functionality
		// Adjust the weight for the specified score_type accordingly
		// STILL assumes that there is only fa_atr and fa_rep to ramp
		// Calculate the fraction complete of the protocol
		Real fraction_completion( Real( current_cycle ) / n_cycles );

		Real adjustment_factor;
		Real const current_weight( scorefxn_->get_weight( score_type ) );

		// Calculate the adjustment_factor needed to adjust the score weight
		if ( current_weight < target_weight ) {
			// if current_weight is lower than target_weight, it's fa_rep ScoreType
			// (fa_rep starts low and gets ramped up higher to its target_weight)
			Real const factor_range_size( 1 - starting_ramp_up_rep_factor_ );
			adjustment_factor =
				starting_ramp_up_rep_factor_ +
				fraction_completion * factor_range_size;
		} else if ( current_weight > target_weight ) {
			// if current_weight is higher than target_weight, it's fa_atr ScoreType
			// (fa_atr starts high and gets ramped up lower to its target_weight)
			Real const factor_range_size( starting_ramp_down_atr_factor_ - 1 );
			adjustment_factor =
				starting_ramp_down_atr_factor_ -
				fraction_completion * factor_range_size;
		} else {
			// Otherwise, if the current_weight is the target_weight, change nothing
			adjustment_factor = 1;
		}

		// Set the ramped weight based on the adjustment_factor calculated
		Real const new_weight( target_weight * adjustment_factor );
		scorefxn_->set_weight( score_type, new_weight );

	} // end score ramping

	// Relay the score_type's new weight to the user
	std::string const cycle_status
		( ( current_cycle == n_cycles ) ? "Final " : "Current " );
	TR.Info << cycle_status + utility::to_string(score_type) +
		" weight: " << scorefxn_->get_weight( score_type ) << std::endl;

} // END ramp_score_weight


void
GlycanDockProtocol::perform_stage2_docking_and_optimization
( core::pose::Pose & pose,
	protocols::moves::RandomMoverOP stage2_rb_randmover,
	protocols::moves::RandomMoverOP stage2_tor_randmover,
	minimization_packing::PackRotamersMoverOP full_packer,
	minimization_packing::EnergyCutRotamerTrialsMoverOP ecut_packer,
	minimization_packing::MinMoverOP minimizer,
	moves::MonteCarloOP mc )
{
	// Pick the Mover type(s) to perform this cycle
	// Doing this manually to be able to control order of sampling type
	core::Real const mover_choice = numeric::random::rg().uniform();

	if ( mover_choice <= 0.5 ) {
		/////
		// Perform rigid-body moves then torsion moves //
		/////
		n_rb_moves_accepted_ +=
			do_stage2_sample_pack_min_cycle
			( pose, stage2_rb_randmover, n_rb_rounds_, // rb
			full_packer, ecut_packer, full_pack_freq_,
			minimizer, mc );
		n_rb_moves_made_ += n_rb_rounds_;
		n_rb_cycles_performed_ += 1;
		// Follow up RB sampling with torsion angle sampling
		if ( glycoligand_has_dihedrals_ ) {
			n_tor_moves_accepted_ +=
				do_stage2_sample_pack_min_cycle
				( pose, stage2_tor_randmover, n_tor_rounds_, // tor
				full_packer, ecut_packer, full_pack_freq_,
				minimizer, mc );
			n_tor_moves_made_ += n_tor_rounds_;
			n_tor_cycles_performed_ += 1;
		} else {
			// If this is a monosaccharide, there are no torsions,
			// so do more rigid-body sampling instead
			n_rb_moves_accepted_ +=
				do_stage2_sample_pack_min_cycle
				( pose, stage2_rb_randmover, n_rb_rounds_, // rb
				full_packer, ecut_packer, full_pack_freq_,
				minimizer, mc );
			n_rb_moves_made_ += n_rb_rounds_;
			n_rb_cycles_performed_ += 1;
		} // mover_choice <= 0.5

	} else {
		/////
		// Perform torsion moves then rigid-body moves //
		/////
		if ( glycoligand_has_dihedrals_ ) {
			n_tor_moves_accepted_ +=
				do_stage2_sample_pack_min_cycle
				( pose, stage2_tor_randmover, n_tor_rounds_, // tor
				full_packer, ecut_packer, full_pack_freq_,
				minimizer, mc );
			n_tor_moves_made_ += n_tor_rounds_;
			n_tor_cycles_performed_ += 1;
		} else {
			// If this is a monosaccharide, there are no torsions,
			// so do rigid-body sampling instead
			n_rb_moves_accepted_ +=
				do_stage2_sample_pack_min_cycle
				( pose, stage2_rb_randmover, n_rb_rounds_, // rb
				full_packer, ecut_packer, full_pack_freq_,
				minimizer, mc );
			n_rb_moves_made_ += n_rb_rounds_;
			n_rb_cycles_performed_ += 1;
		}
		// And then perform rigid-body moves
		n_rb_moves_accepted_ +=
			do_stage2_sample_pack_min_cycle
			( pose, stage2_rb_randmover, n_rb_rounds_, // rb
			full_packer, ecut_packer, full_pack_freq_,
			minimizer, mc );
		n_rb_moves_made_ += n_rb_rounds_;
		n_rb_cycles_performed_ += 1;
	}

	n_stage2_cycles_performed_ += 1;

} // END perform_stage2_docking_and_optimization


////////////////////////////////////////////////////////////////////////////////
bool
GlycanDockProtocol::docking_filter
( core::pose::Pose const & pose,
	core::Size const repeat )
{
	bool passed_filter = true;
	core::pose::Pose working_pose( pose );

	// Interface score should be below 0
	// TODO future development
	// Make a setable value for intE from the command line
	// and/or allow user to define which filter
	// although this would likely be controlled via RosettaScripts filters
	core::Real const interaction_energy_cutoff = 0.0;

	// Interface score calculations
	core::Real const interaction_energy =
		protocols::docking::calc_interaction_energy
		( working_pose, scorefxn_, movable_jump_ );

	// Criterion for failure
	if ( interaction_energy >= interaction_energy_cutoff ) {
		passed_filter = false;
	}

	std::string const filter_status( ( passed_filter ) ? "PASSED" : "FAILED" );
	TR << "STRUCTURE " << filter_status <<
		" DOCKING FILTER: " << interaction_energy << std::endl;
	// If it didn't pass the docking filter, inform user based
	// on which repeat of n_repeats_ this is
	if ( ! passed_filter ) {
		// If it didn't pass and there are more n_repeats_,
		// repeat docking cycles and make new decoy
		if ( repeat != n_repeats_ ) {
			TR << "TRAJECTORY REPEAT " << repeat <<
				" of " << n_repeats_ << std::endl;
		} else {
			// Otherwise, this was our last repeat and it failed, dump last made decoy
			TR << "ALL " << n_repeats_ <<
				" REPEATS FAILED. DUMPING LAST DECOY" << std::endl;
		}
	}

	return passed_filter;

} // END docking_filter


////////////////////////////////////////////////////////////////////////////////
void
GlycanDockProtocol::record_pose_metrics
( core::pose::Pose & pose, // metrics written to pose
	//std::string const & partners,
	utility::vector1< bool > const & glycolig_subset,
	core::pose::Pose const & ref_pose, // should be already scored!
	utility::vector1< bool > const & ref_glycolig_subset )
{
	using namespace core::scoring; // carb ring atom and rmsd
	using namespace protocols::docking; // interaction energy

	// Convert the ChainSelector --> glycolig_subset to --> glycolig_resnums
	utility::vector1< core::Size > const glycolig_resnums =
		core::select::get_residues_from_subset( glycolig_subset, true );
	utility::vector1< core::Size > const ref_glycolig_resnums =
		core::select::get_residues_from_subset( ref_glycolig_subset, true );

	// Calculate metrics relevant to output from standard GlycanDock
	// this also ensures the energies neighbor graph is up-to-datee
	core::Real const total_score = ( *scorefxn_ )( pose );
	core::Real const ring_Srmsd
		( rmsd_with_super( ref_pose, pose, // const
		ref_glycolig_resnums, glycolig_resnums,
		is_carbohydrate_ring_atom) );
	core::Real const heavy_Srmsd
		( rmsd_with_super( ref_pose, pose, // const
		ref_glycolig_resnums, glycolig_resnums,
		is_heavyatom) );
	// ring_Srmsd is the rmsd of the shape of the glycan
	// (aligns the pose glycoligand onto the ref glycoligand)
	core::pose::setPoseExtraScore( pose, "ring_Srmsd", ring_Srmsd );
	core::pose::setPoseExtraScore( pose, "heavy_Srmsd", heavy_Srmsd );

	// Get HBond information at the interface
	/*
	protocols::residue_selectors::HBondSelectorOP intf_hbond_stor
	( utility::pointer::make_shared
	< protocols::residue_selectors::HBondSelector >() );
	intf_hbond_stor->set_hbond_energy_cutoff( -0.5 );
	intf_hbond_stor->set_include_bb_bb( true );
	intf_hbond_stor->set_scorefxn( scorefxn_ );
	core::simple_metrics::metrics::SelectedResiduesPyMOLMetricOP SRPyMOL
	( utility::pointer::make_shared
	< core::simple_metrics::metrics::SelectedResiduesPyMOLMetric >() );
	for ( core::Size glycolig_resnum : glycolig_resnums ) {
	utility::vector1< bool > single_glycolig_subset =
	( utility::pointer::make_shared
	< core::select::residue_selector::FalseResidueSelector >() )->apply( pose );
	single_glycolig_subset[glycolig_resnum] = true;
	intf_hbond_stor->set_input_set_selector
	( core::select::get_residue_selector_from_subset( single_glycolig_subset ) );
	SRPyMOL->set_residue_selector
	( core::select::get_residue_selector_from_subset
	( intf_hbond_stor->apply( pose ) ) );
	core::pose::setPoseExtraScore( pose,
	"pymol_hbonders_to-" +
	pose.pdb_info()->pose2pdb
	( glycolig_resnum, "" ),
	SRPyMOL->calculate( pose ) );
	}
	*/

	// Calculate interaction energy
	core::Real const interaction_energy
		( calc_interaction_energy( pose, // const
		scorefxn_, movable_jump_ ) );

	// Calculate the standard docking decoy metrics using the reference pose
	// If no native pose given on input, the starting structure is used
	// Note the calc functions use the FoldTree definition from pose
	core::Real const ring_Lrmsd
		( rmsd_no_super( ref_pose, pose, // const
		ref_glycolig_resnums, glycolig_resnums,
		is_carbohydrate_ring_atom) );
	core::Real const heavy_Lrmsd
		( rmsd_no_super( ref_pose, pose, // const
		ref_glycolig_resnums, glycolig_resnums,
		is_heavyatom) );
	utility::vector1< core::Real > GD_intf_metrics =
		calc_GlycanDock_intf_metrics( pose, ref_pose,
		glycolig_subset, 5 );
	core::Real const n_intf_residues = GD_intf_metrics[1];
	core::Real const n_nat_intf_residues = GD_intf_metrics[2];
	core::Real const Fnat_intf_res = GD_intf_metrics[3];
	core::Real const n_intf_res_contacts = GD_intf_metrics[4];
	core::Real const n_nat_intf_res_contacts = GD_intf_metrics[5];
	core::Real const Fnat = GD_intf_metrics[6];

	// Record and output decoy pose metrics
	TR << "START Metrics for decoy:" << std::endl;
	TR << "-total_score:                     " << total_score << std::endl;
	TR << "-interaction_energy:              " << interaction_energy << std::endl;
	TR << "-ring_Lrmsd:                      " << ring_Lrmsd << std::endl;
	TR << "-heavy_Lrmsd:                     " << heavy_Lrmsd << std::endl;
	TR << "-ring_Srmsd:                      " << ring_Srmsd << std::endl;
	TR << "-heavy_Srmsd:                     " << heavy_Srmsd << std::endl;
	TR << "-N interface residues:            " << n_intf_residues << std::endl;
	TR << "-N interface residue contacts:    " << n_intf_res_contacts << std::endl;
	TR << "-Fraction of native contacts:     " << Fnat << std::endl;
	TR << "-Fnat interface residues:         " << Fnat_intf_res << std::endl;
	TR << "END Metrics for decoy" << std::endl;

	// Add the extra decoy metrics to the scorefile and the decoy
	core::pose::setPoseExtraScore( pose, "interaction_energy",
		interaction_energy );
	core::pose::setPoseExtraScore( pose, "n_intf_residues",
		n_intf_residues );
	core::pose::setPoseExtraScore( pose, "n_nat_intf_residues",
		n_nat_intf_residues );
	core::pose::setPoseExtraScore( pose, "Fnat_intf_residues", Fnat_intf_res );
	core::pose::setPoseExtraScore( pose, "n_intf_res_contacts",
		n_intf_res_contacts );
	core::pose::setPoseExtraScore( pose, "n_nat_intf_res_contacts",
		n_nat_intf_res_contacts );
	core::pose::setPoseExtraScore( pose, "Fnat", Fnat );
	core::pose::setPoseExtraScore( pose, "ring_Lrmsd", ring_Lrmsd );
	core::pose::setPoseExtraScore( pose, "heavy_Lrmsd", heavy_Lrmsd );
	// And relay which carbohydrate residue was the docking Jump res
	core::pose::setPoseExtraScore( pose, "glycan_Jump_res",
		pose.fold_tree().jump_edge(dock_jump_num_).stop() );

	// After manually getting above metrics, use the InterfaceAnalyzerMover
	// the apply method calls apply_const first, so the pose is const
	// until the data are written to the pose object (output file) itself
	//protocols::analysis::InterfaceAnalyzerMoverOP intf_analyzer
	// ( get_GlycanDock_IAM( partners, scorefxn_ ) );
	//intf_analyzer->apply( pose );

} // END record_pose_metrics


////////////////////////////////////////////////////////////////////////////////
void
GlycanDockProtocol::record_protocol_info( core::pose::Pose & pose )
{
	using namespace core::pose; //::setPoseExtraScore

	TR << "Number of rigid-body cycles performed: " <<
		n_rb_cycles_performed_ << std::endl;
	setPoseExtraScore( pose, "n_rb_cycles", n_rb_cycles_performed_ );
	TR << "Number of rigid-body moves made: " <<
		n_rb_moves_made_ << std::endl;
	setPoseExtraScore( pose, "n_rb_moves_made", n_rb_moves_made_ );
	TR << "Number of rigid-body moves accepted: " <<
		n_rb_moves_accepted_ << std::endl;
	setPoseExtraScore( pose, "n_rb_moves_accepted", n_rb_moves_accepted_ );

	TR << "Number of torsion cycles performed: " <<
		n_tor_cycles_performed_ << std::endl;
	setPoseExtraScore( pose, "n_tor_cycles", n_tor_cycles_performed_ );
	TR << "Number of torsion moves made: " <<
		n_tor_moves_made_ << std::endl;
	setPoseExtraScore( pose, "n_tor_moves_made", n_tor_moves_made_ );
	TR << "Number of torsion moves accepted: " <<
		n_tor_moves_accepted_ << std::endl;
	setPoseExtraScore( pose, "n_tor_moves_accepted", n_tor_moves_accepted_ );

	TR << "Total number of Stage 2 cycles performed: " <<
		n_stage2_cycles_performed_ << std::endl;
	TR << "Total number of moves made: " <<
		n_rb_moves_made_ + n_tor_moves_made_ << std::endl;
	TR << "Total number of moves accepted: " <<
		n_rb_moves_accepted_ + n_tor_moves_accepted_ << std::endl;
	core::Real mc_acceptance =
		( 1. * n_rb_moves_accepted_ + n_tor_moves_accepted_ ) /
		( n_rb_moves_made_ + n_tor_moves_made_ );
	TR << "MonteCarlo acceptance rate: " << mc_acceptance << std::endl;
	setPoseExtraScore( pose, "mc_acceptance", mc_acceptance );

} // END record_protocol_info


////////////////////////////////////////////////////////////////////////////////
void
GlycanDockProtocol::prepack_only
( core::pose::Pose & pose,
	protocols::minimization_packing::PackRotamersMoverOP prepack_packer )
{
	// Prepare the sliders to move the glycoligand away from the protein
	core::Real const trans_magnitude = 1000;
	protocols::rigid::RigidBodyTransMover translate_away =
		protocols::rigid::RigidBodyTransMover( pose, dock_jump_num_ );
	translate_away.step_size( trans_magnitude );
	// translating back to original position is using the same
	// translation mover as above but negated
	protocols::rigid::RigidBodyTransMover translate_back =
		protocols::rigid::RigidBodyTransMover( pose, dock_jump_num_ );
	translate_back.step_size( trans_magnitude );
	translate_back.trans_axis().negate();

	// Slide the glycoligand away
	TR.Debug << "Separating glycoligand from protein receptor" << std::endl;
	translate_away.apply( pose );
	//if ( watch_in_pymol_ ) { pmm_->apply( pose ); }

	// Pack as specified by the given prepack_task
	TR.Debug << "Packing separated complex" << std::endl;
	prepack_packer->apply( pose );
	//if ( watch_in_pymol_ ) { pmm_->apply( pose ); }

	// TODO
	// Any reason not to minimize the packed sidechains?

	// Return the glycoligand to its starting position
	TR.Debug << "Returning glycoligand to starting position" << std::endl;
	translate_back.apply( pose );
	//if ( watch_in_pymol_ ) { pmm_->apply( pose ); }

}


////////////////////////////////////////////////////////////////////////////////
//////////////////////////////
// APPLY GlycanDockProtocol //
//////////////////////////////
void
GlycanDockProtocol::apply( core::pose::Pose & pose )
{
	using namespace protocols::glycan_docking;
	using namespace basic::options;
	using namespace core::scoring; // fa_atr/fa_rep

	TR.Debug << "Start GlycanDockProtocol apply" << std::endl;

	// Only in:file:s is respected - not in:file:l
	if ( ! option[ OptionKeys::in::file::s ].user() ) {
		utility_exit_with_message("No input structure was provided "
			"via the -in:file:s flag! "
			"Nothing to do here. Exiting");
	}

	// Check more flags that could immediately stop GlycanDock if not right
	// We expect a two-body system (either standard GlycanDock,
	// GlycanDock prepack_only, or GlycanDock refine_only, all of which
	// require a protein-glycoligand system as input and their chain IDs)
	if ( ! partners_provided_ ) {
		utility_exit_with_message("Please provide the chain IDs that "
			"define the protein-glycoligand docking "
			"partners via the '-docking:partners' flag. "
			"E.g -docking:partners AB_X, "
			"where chain A and B define the upstream "
			"protein receptor chains and chain X defines the "
			"downstream glycoligand chain. Note, GlycanDock "
			"has only been tested using glycoligands "
			"defined by a single chain ID.");
	} else {
		// will need pose to ensure chain IDs passed exist and seem correct
		if ( partners_ == "" ) {
			utility_exit_with_message("During GlycanDock, the -docking:partners "
				"flag was passed, but nothing was provided! "
				"Check your input. For example, passing "
				"-docking:partners AB_X tells GlycanDock "
				"that chain A and B define the upstream "
				"protein receptor chains and chain X "
				"defines the downstream glycoligand chain. "
				"Note, GlycanDock has only been tested using "
				"glycoligands defined by a single chain ID.");
		}
	}

	//////////////////////////////
	// Setup Step 0 - Adjust (cmd-line) options depending on
	//                user's selected GlycanDockProtocol usage
	//  (1) - If prepack_only, do... nothing
	//  (2) - If refine_only, do.... nothing
	//--------------------
	// At current implementation, nothing needed for Step 0


	//////////////////////////////
	// Setup Step 1 - Setup an appropriate ScoreFunction with any cmd-line
	//                score weights or constraints provided
	//                 * If run through RosettaScripts, the user can provide
	//                   a scorefxn_ and this method will keep that sf
	//--------------------
	setup_scorefxn( pose ); // fills scorefxn_


	//////////////////////////////
	// Setup Step 2 - Get local copy of reference pose for metric calculations
	//                Looks for -in:file:native or uses input pose as reference
	//  (1) - Create ref_pose
	//--------------------
	core::pose::PoseOP ref_pose = nullptr;
	if ( option[ OptionKeys::in::file::native ].user() ) {
		// User provided a reference structure, grab it
		ref_pose = utility::pointer::make_shared< core::pose::Pose >();
		core::import_pose::pose_from_file( *ref_pose,
			option[ OptionKeys::in::file::native ](),
			core::import_pose::PDB_file);
	} else {
		// Input pose will be the reference structure
		ref_pose = pose.clone();
	}
	( *scorefxn_ )( *ref_pose );


	//////////////////////////////
	// Setup Step 3 - Setup an appropriate docking FoldTree for the pose
	//                 * Requires knowing partners_ (e.g. AB_X)
	//                 * ref_pose needs docking FoldTree too for metric calcs
	//--------------------
	setup_docking_foldtree( pose, *ref_pose );


	//////////////////////////////
	// Setup Step 4 - Setup the MonteCarlo object
	//--------------------
	protocols::moves::MonteCarloOP mc =
		utility::pointer::make_shared< protocols::moves::MonteCarlo >
		( pose, *scorefxn_, // affected by ramp_scorefxn
		mc_kt_ ); // 0.6


	//////////////////////////////
	// Setup Step 5 - Prepare ResidueSelectors and ResidueSubsets for pose
	//   (1) - glycolig_stor: selects the specific glycoligand chain
	//         by user through (-docking:)partners flag/xml option;
	//          * util setup function also performs some quality control
	//   (2) - glycolig_subset: subset for all glycoligand residues
	//          * is used to prepare sampling Movers and for metrics
	//          * n_glyc_residues == size glycoligand
	//   (3) - glycolig_subset_with_torsions: selects the glycan
	//         residues selected by glycolig_stor that have (movable)
	//         torsions. Reasoning is that reducing-end sugars have no parent
	//         residue and therefore have no glycosidic torsion angles,
	//         so there is no relevant torsion sampling to do.
	//         Therefore, don't give that residue to a sampling Mover
	//         as that would just be wasting computational resources/time
	//          * if there are no residues with torsions, then it should
	//            be a monosaccharide (warn user if n_glyc_residues != 1)
	//   (4) - Same as 1 but for the reference pose
	//--------------------
	// (1)
	core::select::residue_selector::ChainSelectorCOP glycolig_stor =
		setup_glycoligand_selector( partners_ );
	// (2)
	utility::vector1< bool > const glycolig_subset = glycolig_stor->apply( pose );
	core::Size const n_glyc_residues =
		core::select::residue_selector::count_selected( glycolig_subset );
	if ( n_glyc_residues == 0 ) {
		utility_exit_with_message("There are no carbohydrate residues "
			"found in the glycoligand! Did you pass the "
			"right glycan chain ID to GlycanDock? "
			"Please check your input.");
	}
	if ( TR.Debug.visible() ) {
		TR.Debug << "The following " << n_glyc_residues << " residues "
			"have been identified as components of the glycoligand "
			"given the provided chain ID" << std::endl;
		for ( core::Size ii = 1; ii <= pose.size(); ++ii ) {
			if ( glycolig_subset[ii] ) {
				TR.Debug << pose.pdb_info()->pose2pdb(ii) << std::endl;
			}
		}
	}
	// (3)
	utility::vector1< bool > const glycolig_subset_with_torsions =
		get_glycolig_subset_with_torsions( pose, glycolig_stor );
	if ( core::select::residue_selector::all_false_selection
			( glycolig_subset_with_torsions ) ) {
		glycoligand_has_dihedrals_ = false;
		TR << "Glycoligand has no dihedrals! Not performing "
			"any glycosidic torsion angle sampling" << std::endl;
		if ( n_glyc_residues != 1 ) {
			TR.Warning << "If your glycoligand is NOT a monosaccharide, "
				"check your input!" << std::endl;
		}
	}
	// (4)
	// the same glycolig_stor for the pose should work on the ref_pose,
	// but going through the setup function for its quality control
	core::select::residue_selector::ChainSelectorCOP ref_glycolig_stor =
		setup_glycoligand_selector( partners_ );
	utility::vector1< bool > const ref_glycolig_subset =
		ref_glycolig_stor->apply( *ref_pose );
	if ( core::select::residue_selector::count_selected
			( ref_glycolig_subset ) == 0 ) {
		utility_exit_with_message("There are no carbohydrate residues "
			"found in the glycoligand in the reference pose!! Did you pass "
			"the right glycan chain ID to GlycanDock? "
			"Please check your input.");
	}


	//////////////////////////////
	// Setup Step 6 - create packing objects
	//  (1) - packer_tf: the TaskFactory for all packers
	//  (2) - pack_rot: 'full' packing i.e. PackRotamersMover
	//  (3) - ecut_rot_trials: when not doing 'full' packing, use this
	//--------------------
	// (1)
	// if prepack_only_ == true, the tf will not include RestrictToInterface
	core::pack::task::TaskFactoryOP packer_tf =
		setup_GlycanDock_taskfactory( dock_jump_num_, intf_pack_dist_,
		prepack_only_ );
	// Debug the initial packer task
	// (Helpful during benchmarking to ensure interface defined correctly)
	if ( TR.Debug.visible() ) {
		TR.Debug << "Initial packer task created for GlycanDockProtocol:" << std::endl;
		core::pack::task::PackerTaskOP task =
			packer_tf->create_task_and_apply_taskoperations( pose );
		task->show( TR.Debug );
	}
	// (2)
	protocols::minimization_packing::PackRotamersMoverOP pack_rot =
		utility::pointer::make_shared
		< protocols::minimization_packing::PackRotamersMover >();
	pack_rot->score_function( scorefxn_ ); // affected by ramp_scorefxn
	// TaskFactory always overrides/regenerates task in PackRotamersMover
	pack_rot->task_factory( packer_tf );
	// (3)
	protocols::minimization_packing::EnergyCutRotamerTrialsMoverOP ecut_rot_trials =
		utility::pointer::make_shared
		< protocols::minimization_packing::EnergyCutRotamerTrialsMover >
		( scorefxn_, packer_tf, mc, // affected by ramp_scorefxn
		rt_energycut_ ); // 0.05


	//////////////////////////////
	// Setup Step 7 - create movemap and minimization objects
	//  (1) - minimizer_mm: the MoveMap for the minimizer
	//  (2) - minimizer: the MinMover for minimization
	//--------------------
	// (1)
	core::kinematics::MoveMapOP minimizer_mm =
		setup_GlycanDock_movemap( pose, glycolig_stor, lock_rings_, dock_jump_num_ );
	// (2)
	protocols::minimization_packing::MinMoverOP minimizer =
		utility::pointer::make_shared< protocols::minimization_packing::MinMover >
		( minimizer_mm, scorefxn_, // affected by ramp_scorefxn
		min_type_, // "lbfgs_armijo_nonmonotone"
		min_tolerance_, // 0.01
		true );
	if ( TR.Debug.visible() ) {
		TR.Debug << "Minimizer for GlycanDockProtocol" << std::endl;
		minimizer->show( TR.Debug );
	}


	//////////////////////////////
	// Setup Step 8 - setup Stage 1 conformation initialization movers
	//  (1) - stage1_rb_seqmover: gaussian perturb on rb CoM (trans + rot)
	//         * Performs sequence of rigid-body perturbations
	//            1. Rotational + Translational perturbation of glycolig CoM
	//            2. (If specified) Complete, random perturbation of
	//                              glycolig CoM in 360-degree space
	//                               * (redundant with 1. CoM rot perturbation)
	//  ~I don't know of a uniform torsion angle pertub mover~
	//    Need something like protocols::dihedral::UniformPerturbMover
	//    TODO future development (stage1_tor_mover)
	//  (2) - ~~~: uniform +/- perturb on current torsions
	//         * For now, perform stage1 tor perturbation in util file
	//--------------------
	protocols::moves::SequenceMoverOP stage1_rb_seqmover = nullptr;
	if ( ! refine_only_ ) {
		// (1)
		stage1_rb_seqmover = setup_GlycanDock_stage1_rb_seqmover
			( pose, dock_jump_num_,
			stage1_rot_mag_, stage1_trans_mag_,
			stage1_rotate_glyc_about_com_ );
		// ( 2)
		// ~~~ stage1_tor_mover: defined as a method of this protocol in util
	}


	//////////////////////////////
	// Setup Step 9 - setup Stage 2 docking and optimization movers
	//  (1) - stage2_rb_randmover: weighted RandomMover for rigid-body
	//  (2) - stage2_tor_randmover: weighted RandomMover for torsion angle
	//--------------------
	// (1)
	protocols::moves::RandomMoverOP stage2_rb_randmover =
		setup_GlycanDock_stage2_rb_randmover
		( dock_jump_num_, stage2_rot_mag_, stage2_trans_mag_,
		slide_glyc_into_contact_ );
	// (2)
	protocols::moves::RandomMoverOP stage2_tor_randmover =
		setup_GlycanDock_stage2_tor_randmover
		( pose, glycolig_subset_with_torsions,
		scorefxn_, mc_kt_, n_shear_moves_, refine_only_ );


	//------------------------------------------------------------
	// DONE WITH SETUP STEPS
	// Inform the user of the settings and behavior for this protocol
	show( TR );

	///// Mutually exclusive protocol usages

	if ( prepack_only_ ) {
		//////////////////////////
		// PREPACK COMPLEX ONLY //
		//////////////////////////
		TR << "Begin GlycanDock prepack_only protocol" << std::endl;
		prepack_only( pose, pack_rot );
		TR << "End GlycanDock prepack_only protocol" << std::endl;


	} else {
		//////////////////////////////
		// MAIN GlycanDock PROTOCOL //
		//////////////////////////////
		TR << "Begin GlycanDock protein-glycoligand docking protocol" << std::endl;

		////////////////////////////////////////////////////
		// STAGE 1 GlycanDock CONFORMATION INITIALIZATION //
		////////////////////////////////////////////////////
		// Start GlycanDock with an initial randomized trajectory
		// If refine_only_ is true, will NOT perform Stage 1
		//  * use refine_only if you have high(ish) confidence in your
		//    input structure and just want to get the glycoligand
		//    conformation optimized without moving much in its pocket
		if ( ! refine_only_ ) {
			TR << "Begin GlycanDock Stage 1 conformation initialization" << std::endl;

			// Conformation initialization technique developed for GlycanDock
			// Really, just performing a random rigid-body perturbation
			// and randomly perturbing glycosidic torsion angles a bit
			// no MonteCarlo check or minimization - just apply it and go
			// --Note, stage1_tor_mover does not currently exist--
			//    * Applying uniform torsion angle perturbation manually
			do_stage1_conformation_initialization
				( pose, stage1_rb_seqmover,
				/* stage1_tor_mover when that exists; done manually */
				glycolig_subset,
				glycoligand_has_dihedrals_,
				stage1_torsion_uniform_pert_mag_ );

		} else { // END Stage 1 conformation initialization
			TR << "The 'refine_only' behavior was set." << std::endl <<
				"Skipping GlycanDock Stage 1 conformation initialization" << std::endl;
		}

		// Store a copy of the starting structure (after Stage 1, if applied)
		// If n_repeats_ > 1, this starting pose will be used to replace
		// the working decoy getting docked if it fails the docking filter
		// TODO QUESTION
		// How to copy pose correctly? Is the cloning required? Use .assign()?
		core::pose::Pose const starting_pose( *pose.clone() );
		// Use starting_pose for RMSD metrics for any benchmarking data

		// Output starting information for the trajectory
		TR << "Starting score for complex:" << std::endl;
		scorefxn_->show( TR, pose );


		/////////////////////////////////////////////////
		// STAGE 2 GlycanDock DOCKING AND OPTIMIZATION //
		/////////////////////////////////////////////////
		TR << "Begin GlycanDock Stage 2 docking and optimization" << std::endl;

		// Wrap docking protocol with a docking filter repeat
		// If the decoy at the end of the protocol has an
		// unfavorable interaction energy, redo the trajectory
		// starting from the stored starting_pose
		core::Size repeat;
		for ( repeat = 1; repeat <= n_repeats_; ++repeat ) {
			TR << "Repeat: " << repeat << " of " << n_repeats_ << std::endl;

			// If this is a repeat docking cycle (i.e. repeat > 1),
			// Get a fresh copy of the starting pose and
			// Reset counters used for tracking protocol stats
			if ( repeat > 1 ) {
				// TODO QUESTION
				// This feels wrong since I'm not using an assign
				pose = starting_pose;
				clear_counters();
			}

			// START Stage 2 outer cycles
			// outer cycles control ramping the scorefxn
			for ( core::Size cycle( 1 ); cycle <= n_cycles_; ++cycle ) {
				TR << "Stage 2 inner sampling and optimization cycle: " <<
					cycle << " of " << n_cycles_ << std::endl;

				// ramping scorefunction weights is true by default
				if ( ramp_scorefxn_ ) {
					ramp_score_weight( fa_atr, target_atr_, cycle, n_cycles_ );
					ramp_score_weight( fa_rep, target_rep_, cycle, n_cycles_ );
					// make double sure that objects that use the scorefxn_
					// are aware that it got changed
					// The MonteCarlo object relies on the last pose scored
					mc->reset_scorefxn( pose, *scorefxn_ );
					pack_rot->score_function( scorefxn_ );
					ecut_rot_trials->score_function( scorefxn_ );
					minimizer->score_function( scorefxn_ );
				}

				// START Stage 2 inner sampling cycles
				// inner cycles control rigid-body and torsion sampling and optimization
				perform_stage2_docking_and_optimization
					( pose,
					stage2_rb_randmover, stage2_tor_randmover,
					pack_rot, ecut_rot_trials,
					minimizer, mc );

			} // END Stage 2 outer cycles

			///////////////////////////////////////////
			// Check if decoy passed docking filters //
			///////////////////////////////////////////
			// scorefxn_ should be returned to normal weights by/during the last cycle
			const bool passed_filter = docking_filter( pose, repeat );
			// if passed filter, don't repeat protocol, dump this pose
			if ( passed_filter ) { break; }

		} // END docking repeat filter

		//--------------------------------------------------
		TR << "End GlycanDock protein-glycoligand docking protocol" << std::endl;

		/////////////////////////////////
		// REPORT PROTOCOL INFORMATION //
		/////////////////////////////////
		// Report number of rigid-body and torsion moves made and accepted
		record_protocol_info( pose );

	} // END main GlycanDockProtocol


	//--------------------------------------------------
	/////////////////////////
	// RECORD POSE METRICS //
	/////////////////////////
	// Output final scoring information to user
	// Doing final score here to ensure scorefile contains information
	// when the user employs any non-docking version of the protocol
	// e.g. prepack_only
	TR << "Final score for decoy:" << std::endl;
	scorefxn_->show( TR, pose );

	// Record output decoy pose metrics to be stored in the score file
	// but for prepacking, we are not interested in these data
	if ( ! prepack_only_ ) {
		record_pose_metrics( pose, glycolig_subset,
			*ref_pose, ref_glycolig_subset );
	}

	TR.Debug << "Finished GlycanDockProtocol apply" << std::endl;

} // END apply


//////////////////////////////////////////////////////////////////////
///////////////////////////////
/// Rosetta Scripts Support ///
///////////////////////////////

/// @brief parse XML tag (to use GlycanDockProtocol in Rosetta Scripts)
void
GlycanDockProtocol::parse_my_tag( utility::tag::TagCOP tag,
	basic::datacache::DataMap & )
{
	// If options or defaults change, change in basic/options/options_rosetta.py
	// AND below in the XML documentation (in provide_xml_schema)
	// The tag part is getting the user-provided value for the flag,
	// but if the user did not specify the flag in the XML,
	// we are passing a default value to use in its place
	// tag->getOption<type>( "option/flag name", default used if it wasn't set)

	// mutually-exclusive GlycanDock protocols
	set_refine_only( (tag->getOption<bool>( "refine_only", false )) );
	set_prepack_only( (tag->getOption<bool>( "prepack_only", false )) );

	if ( tag->hasOption( "partners" ) && tag->hasOption( "glycan_chain" ) ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,
			"'partners' and 'glycan_chain' are mutually "
			"exclusive parameters. Use either 'partners' "
			"(for performing protein-glycoligand docking "
			"or pre-packing) or 'glycan_chain' (for "
			"performing glycan sampling only)");
	}

	set_partners( tag->getOption<std::string>( "partners", "_") );

	// Stage 1 - conformation initialization
	set_stage1_rotate_glycan_about_com( (tag->getOption<bool>( "stage1_rotate_glycan_about_com", false )) );
	set_stage1_perturb_glycan_com_trans_mag( (tag->getOption<core::Real>( "stage1_perturb_glycan_com_trans_mag", 0.5 )) );
	set_stage1_perturb_glycan_com_rot_mag( (tag->getOption<core::Real>( "stage1_perturb_glycan_com_rot_mag", 7.5 )) );
	set_stage1_torsion_uniform_pert_mag( (tag->getOption<core::Real>( "stage1_torsion_uniform_pert_mag", 12.5 )) );

	// Stage 2 - local docking and refinement
	set_n_repeats( (tag->getOption<core::Size>( "n_repeats", 1 )) );
	set_mc_kt( (tag->getOption<core::Real>( "mc_kt", 0.6 )) );
	set_n_rigid_body_rounds( (tag->getOption<core::Size>( "n_rigid_body_rounds", 8 )) );
	set_n_torsion_rounds( (tag->getOption<core::Size>( "n_torsion_rounds", 8 )) );
	set_stage2_trans_mag( (tag->getOption<core::Real>( "stage2_trans_mag", 0.5 )) );
	set_stage2_rot_mag( (tag->getOption<core::Real>( "stage2_rot_mag", 7.5 )) );
	set_full_packing_frequency( (tag->getOption<core::Size>( "full_packing_frequency", 8 )) );
	set_interface_packing_distance( (tag->getOption<core::Real>( "interface_packing_distance", 16.0 )) );
	set_ramp_scorefxn( (tag->getOption<bool>( "ramp_scorefxn", true )) );

	// misc
	set_watch_in_pymol( (tag->getOption<bool>( "watch_in_pymol", false )) );

} // END parse_my_tag


void
GlycanDockProtocol::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	// If options or defaults change, change in basic/options/options_rosetta.py too
	attlist
		// use 'XMLSchemaAttribute::required_attribute' a required XML attribute
		+ XMLSchemaAttribute::attribute_w_default
		("refine_only", xsct_rosetta_bool,
		"Perform refinement of the input putative complex only. Skips Stage 1 "
		"(conformational initialization via a random perturbation) and, "
		"during Stage 2, do not perform large perturbations in glycosidic "
		"torsion angle space. Default = false",
		"false")
		+ XMLSchemaAttribute::attribute_w_default
		("prepack_only", xsct_rosetta_bool,
		"Perform Stage 0 pre-packing of the input putative complex only. "
		"Separates the glycoligand from its protein receptor and optimizes "
		"all sidechain rotamer conformations. Default = false",
		"false")
		+ XMLSchemaAttribute::attribute_w_default
		("partners", xs_string,
		"Chain IDs of protein-glycoligand docking partneers. Synonymous "
		"with the -docking:partners flag. Required when performing "
		"docking or the pre-pack protocol. E.g. A_X or HL_X. Default = '_'",
		"_")
		+ XMLSchemaAttribute::attribute_w_default
		("stage1_rotate_glycan_about_com", xsct_rosetta_bool,
		"During Stage 1 conformation initialization, rotate the glycoligand "
		"about its center-of-mass in uniform 3D space. Default = false. "
		"Recommended to set to true if confidence of the glycoligand's "
		"rigid-body orientation in the putative binding site is low.",
		"false")
		+ XMLSchemaAttribute::attribute_w_default
		("stage1_perturb_glycan_com_trans_mag", xsct_real,
		"During Stage 1 conformation initialization, this is the magnitude "
		"for performing a random translational Gaussian perturbation on the "
		"glycoligand's center-of-mass. Default = 0.5 (Angstroms)",
		"0.5")
		+ XMLSchemaAttribute::attribute_w_default
		("stage1_perturb_glycan_com_rot_mag", xsct_real,
		"During Stage 1 conformation initialization, this is the magnitude "
		"for performing a random rotational Gaussian perturbation on the "
		"glycoligand's center-of-mass. Default = 7.5 (degrees)",
		"7.5")
		+ XMLSchemaAttribute::attribute_w_default
		("stage1_torsion_uniform_pert_mag", xsct_real,
		"During Stage 1 conformation initialization, the magnitude used to "
		"perform a uniform perturbation on each phi, psi, and omega glycosidic "
		"torsion angle. Default = 12.5 (degrees).",
		"12.5")
		+ XMLSchemaAttribute::attribute_w_default
		("n_repeats", xsct_non_negative_integer,
		"Number of times to run the GlycanDock protocol on a "
		"protein-glycoligand system if the final docked structure does not "
		"pass the quality filter (negative Rosetta interaction energy). Default = 1",
		"1")
		+ XMLSchemaAttribute::attribute_w_default
		("mc_kt", xsct_real,
		"During Stage 2 docking and refinement, the value of kT used to "
		"accept or reject moves based on the Metropolis criterion. Default = 0.6",
		"0.6")
		+ XMLSchemaAttribute::attribute_w_default
		("n_rigid_body_rounds", xsct_non_negative_integer,
		"During Stage 2 docking and refinement, the number of rigid-body "
		"sampling rounds to perform each cycle. Default = 8",
		"8")
		+ XMLSchemaAttribute::attribute_w_default
		("n_torsion_rounds", xsct_non_negative_integer,
		"During Stage 2 docking and refinement, the number of glycosidic "
		"torsion angle sampling rounds to perform each cycle. Default = 8",
		"8")
		+ XMLSchemaAttribute::attribute_w_default
		("stage2_trans_mag", xsct_real,
		"During Stage 2 docking and refinement, this is the magnitude for "
		"performing a translational Gaussian perturbation on the glycoligand's "
		"center-of-mass. The perturbation is accepted or rejected based on "
		"the Metropolis criterion. Default = 0.5 (Angstroms)",
		"0.5")
		+ XMLSchemaAttribute::attribute_w_default
		("stage2_rot_mag", xsct_real,
		"During Stage 2 docking and refinement, this is the magnitude for "
		"performing a rotational Gaussian perturbation on the glycoligand's "
		"center-of-mass. The perturbation is accepted or rejected based on "
		"the Metropolis criterion. Default = 7.5 (degrees)",
		"7.5")
		+ XMLSchemaAttribute::attribute_w_default
		("full_packing_frequency", xsct_non_negative_integer,
		"During Stage 2 docking and refinement, the frequency at which a "
		"full packing operation via the PackRotamersMover should be applied. "
		"Default = 8. Should be a factor of n_rigid_body_rounds and "
		"n_torsion_rounds. When not performing a full packing operation, "
		"the faster EnergyCutRotamerTrialsMover is applied.",
		"8")
		+ XMLSchemaAttribute::attribute_w_default
		("interface_packing_distance", xsct_real,
		"During Stage 2 docking and refinement, the distance used to define "
		"protein-glycoligand interface residues for packing. Default = 16 "
		"(Angstroms). Used for the RestrictToInterface task operation.",
		"16.0")
		+ XMLSchemaAttribute::attribute_w_default
		("ramp_scorefxn", xsct_rosetta_bool,
		"During Stage 2 docking and refinement, ramp the fa_atr and fa_rep "
		"score terms. fa_atr is set high and fa_rep is set low, and then "
		"ramped to their starting weights incrementally over the course of "
		"n_cycles. Default = true. Used to promote sampling by not strictly "
		"enforcing rigid sterics in the early stages of the protocol.",
		"true")
		+ XMLSchemaAttribute::attribute_w_default
		("watch_in_pymol", xsct_rosetta_bool,
		"Watch the GlycanDock protocol in PyMOL? Sends the Pose at specific, "
		"hard-coded steps to PyMOL. Default = false. Used as an alternative "
		"to -show_simulation_in_pymol.",
		"false");

	protocols::moves::xsd_type_definition_w_attributes
		( xsd, mover_name(),
		"GlycanDock performs protein-glycoligand local docking given "
		"a putative complex as input. All GlycanDockProtocol attributes "
		"are set to their default, and suggested values based on the "
		"results from the GlycanDock benchmark assessment.",
		attlist );

}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
GlycanDockProtocol::fresh_instance() const
{
	return protocols::moves::MoverOP
		( utility::pointer::make_shared< GlycanDockProtocol >() );
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
GlycanDockProtocol::clone() const
{
	return utility::pointer::make_shared< GlycanDockProtocol >( *this );
}

std::string
GlycanDockProtocol::get_name() const
{
	return mover_name();
}

std::string
GlycanDockProtocol::mover_name()
{
	return "GlycanDockProtocol";
}

///////////////////////////////
//////////  Creator  //////////
///////////////////////////////
protocols::moves::MoverOP
GlycanDockProtocolCreator::create_mover() const
{
	return protocols::moves::MoverOP
		( utility::pointer::make_shared< GlycanDockProtocol >() );
}

std::string
GlycanDockProtocolCreator::keyname() const
{
	return GlycanDockProtocol::mover_name();
}

void
GlycanDockProtocolCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	GlycanDockProtocol::provide_xml_schema( xsd );
}


//////////////////////////////////////////////////////////////////////
///////////////////////////////
/////// Setters/Getters ///////
///////////////////////////////

// all defined in .hh file


//////////////////////////////////////////////////////////////////////
///////////////////////////////
/////// Private Methods ///////
///////////////////////////////

std::ostream &
operator<<( std::ostream & os, GlycanDockProtocol const & mover )
{
	mover.show(os);
	return os;
}


} // glycan_docking
} // protocols
