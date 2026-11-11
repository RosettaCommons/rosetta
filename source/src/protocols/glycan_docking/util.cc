// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/glycan_docking/util.cc
/// @brief Utility methods for protocols/glycan_docking/GlycanDockProtocol
/// @author Morgan Nance (morganlnance@gmail.com)


// Project Headers
#include <protocols/minimization_packing/MinMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/carbohydrates/RingPlaneFlipMover.hh>
#include <protocols/docking/DockingInitialPerturbation.hh>
#include <protocols/docking/util.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/BBDihedralSamplerMover.hh>
#include <protocols/simple_moves/bb_sampler/SugarBBSampler.hh>
#include <protocols/simple_moves/bb_sampler/SmallBBSampler.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/minimization_packing/RotamerTrialsMover.hh>
#include <protocols/analysis/InterfaceAnalyzerMover.hh>

#include <core/scoring/ScoreFunction.hh>

#include <core/simple_metrics/metrics/SelectedResiduesPyMOLMetric.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <protocols/simple_task_operations/RestrictToInterface.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/carbohydrates/GlycanTreeSet.hh>
#include <core/conformation/carbohydrates/util.hh>

#include <core/kinematics/MoveMap.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/carbohydrates/util.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <core/pose/DockingPartners.hh>
#include <core/types.hh>

#include <core/select/util.hh>
#include <core/select/residue_selector/ChainSelector.hh>
#include <core/select/residue_selector/FalseResidueSelector.hh>
#include <core/select/residue_selector/NeighborhoodResidueSelector.hh>
#include <core/select/residue_selector/ReturnResidueSubsetSelector.hh>
#include <core/select/residue_selector/util.hh>

// Basic Headers
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

// Numeric Headers
#include <numeric/angle.functions.hh>
#include <numeric/random/random.hh>

// Utility Headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>
#include <ObjexxFCL/FArray2D.hh>


static basic::Tracer TR( "protocols.glycan_docking.util" );

namespace protocols {
namespace glycan_docking {


////////////////////////////////////////////////////////////////////////////////

//--------------------------------------------------
core::select::residue_selector::ChainSelectorCOP
setup_glycoligand_selector( core::pose::DockingPartners const & docking_partners )
{
	TR.Debug << "Setting up glycoligand ChainSelector" << std::endl;

	// Here, we are using GlycanDock for a two-body system
	// (an upstream protein receptor and a downstream glycoligand)
	// Turn the downstream component of docking_partners

	// Assert that this is a two-body docking problem
	if ( docking_partners.partner1.empty() || docking_partners.partner2.empty() ) {
		utility_exit_with_message("Provided docking partners '" +
			docking_partners.str() + "' did not result in "
			"a two-body docking problem! E.g. A_X or AB_X. "
			"Check your input and retry.");
	}

	// Systems like A_X and AB_X are fine in GlycanDock,
	// but have not tested systems like A_XY. Likely would break app
	if ( docking_partners.partner2.size() > 1 ) {
		TR.Warning << "GlycanDock has not been tested on "
			"systems in which the glycoligand has more than one chain ID! "
			"Unexpected behavior and results may occur." << std::endl;
	}

	// Create the ChainSelector object
	core::select::residue_selector::ChainSelectorCOP out_glycolig_selector
		( utility::pointer::make_shared
		< core::select::residue_selector::ChainSelector >( docking_partners.partner2 ) );

	TR.Debug << "Finished setting up glycoligand ChainSelector" << std::endl;

	return out_glycolig_selector;

} // END setup_glycoligand_selector


//--------------------------------------------------
utility::vector1< bool >
get_glycolig_subset_with_torsions
( core::pose::Pose const & pose,
	core::select::residue_selector::ChainSelectorCOP glycolig_selector )
{

	TR.Debug << "Collecting glycoligand subset with "
		"movable glycosidic torsions " << std::endl;

	// The full glycoligand subset is needed and won't be changed
	utility::vector1< bool > const glycolig_subset
		( glycolig_selector->apply( pose ) );
	// The glycoligand subset with movable torsion angles first
	// starts as the full glycolig subset and gets pruned as needed
	utility::vector1< bool > out_glycolig_subset_with_movable_torsions
		( glycolig_subset );

	// For every residue in the pose that is in the glycolig_subset,
	//  1 - ensure it is a carbohydrate residue
	//     * if not, exit
	//  2 - determine if it has a carbohydrate parent residue
	//     * if not, exit
	//  3 - if it does NOT have movable torsion angles,
	//      update out_glycolig_subset_with_movable_torsions
	for ( core::uint resnum( 1 ); resnum <= pose.total_residue(); ++resnum ) {
		if ( ! glycolig_subset[resnum] ) { continue; }
		// Exit criteria 1 - selected residue is not a carbohydrate
		if ( ! pose.residue( resnum ).is_carbohydrate() ) {
			utility_exit_with_message("Residue " + pose.pdb_info()->pose2pdb(resnum)
				+ " is set by the provided glycoligand "
				"chain ID but is not a carbohdyrate residue! "
				"Is your glycoligand at the bottom of "
				"the PDB file? Does your glycoligand "
				"have any unrecognized residue types?");
		}

		// Determine presence of glycosidic torsions
		core::Size parent = pose.glycan_tree_set()->get_parent( resnum );
		if ( parent == 0 ) {
			// A carbohydrate residue with no parent has no
			// dihedrals available for torsion sampling
			// Turn this residue off in out_glycolig_subset_with_movable_torsions
			// so we do not choose it and then do nothing with it
			// when performing docking and sampling refinement (wasteful)
			out_glycolig_subset_with_movable_torsions[resnum] = false;
		} else if ( ! pose.residue( parent ).is_carbohydrate() ) {
			// Exit criteria 2 - selected residue has non-carb parent
			// likely a protein-conjugated glycan, which is out-of-scope here
			// NOTE
			// This would need to be changed if the GlycanDockProtocol
			// is ever altered to perform glycoprotein docking
			utility_exit_with_message("Residue " + utility::to_string(resnum) +
				" is set by the provided glycoligand chain ID, "
				"but it has a non-carbohydrate parent residue! "
				"GlycanDock is not setup to function "
				"with protein-conjugated glycoligands.");
		} else {
			// a glycosidic linkage has two or more torsion angles,
			// enforce that n glycosidic torsions is at least not 0
			if ( core::pose::carbohydrates::get_n_glycosidic_torsions_in_res
					( pose, resnum ) == 0 ) {
				utility_exit_with_message("Residue " + utility::to_string(resnum) +
					" is set by the provided glycoligand "
					"chain ID, is a carbohydrate residue, "
					"and has a carbohydrate parent, but it "
					"has 0 or glycosidic torsions! "
					"Is this type of sugar biologically "
					"realistic? Is it in the database? Is it "
					"a monosaccharide that escaped detection?");
			}
			// Otherwise, do nothing to out_glycolig_subset_with_movable_torsions
			// This residue is (1) a carbohydrate, (2) has a carbohydrate parent,
			// and (3) has a non-zero number of glycosidic torsion angles
			// and is already true in out_glycolig_subset_with_movable_torsions
		}
	}

	TR.Debug << "Finished collecting glycoligand subset with "
		"movable glycosidic torsions " << std::endl;

	return out_glycolig_subset_with_movable_torsions;

} // END get_glycolig_subset_with_movable_torsions


//--------------------------------------------------
core::pack::task::TaskFactoryOP
setup_GlycanDock_taskfactory
( core::Size const docking_jump_num,
	core::Real const interface_packing_distance,
	bool const prepack_only )
{
	using namespace core::pack;
	using namespace core::pack::task;
	using namespace core::pack::task::operation;

	// The TaskFactory WILL be overwritten by any user-inputted resfile
	TR.Debug << "Setting up TaskFactory for GlycanDock packing" << std::endl;

	core::pack::task::TaskFactoryOP out_packer_tf
		( utility::pointer::make_shared< TaskFactory >() );

	// Call some standard TaskOperations
	out_packer_tf->push_back
		(TaskOperationCOP(utility::pointer::make_shared
		< InitializeFromCommandline >() ) ); // Extra flags
	out_packer_tf->push_back
		(TaskOperationCOP(utility::pointer::make_shared
		< RestrictToRepacking >() ) ); // No design in GlycanDockProtocol
	out_packer_tf->push_back
		(TaskOperationCOP(utility::pointer::make_shared
		< NoRepackDisulfides >() ) ); // Don't pack disulfides
	// User determines if IncludeCurrent is used based on -use_input_sc flag

	// NOTE
	// Default packing radius is 8 Ang, but upped to 16 Ang default
	// for GlycanDockProtocol during benchmarking because
	// the neighbor atom (nbr atom) of carbohydrates wasn't giving me
	// expected/good behavior during packing
	// --------
	// There is technically NO interface if doing prepack_only since
	//  the protein-glycoligand system gets split apart then packed
	if ( ! prepack_only ) {
		// Restrict to the interface (neighboratoms w/in X Ang)
		out_packer_tf->push_back
			(TaskOperationCOP(utility::pointer::make_shared
			<protocols::simple_task_operations::RestrictToInterface>
			(docking_jump_num, interface_packing_distance)));
	}

	// Give/overwrite the packer task any user-provided resfile
	if ( basic::options::option[basic::options::OptionKeys::packing::resfile].user() ) {
		TR.Debug << "Using user-provided resfile" << std::endl;
		out_packer_tf->push_back
			( utility::pointer::make_shared< operation::ReadResfile >() );
	}

	TR.Debug << "Finished setting up TaskFactory for GlycanDock packing" << std::endl;

	return out_packer_tf;

} // END setup_GlycanDock_taskfactory


//--------------------------------------------------
core::kinematics::MoveMapOP
setup_GlycanDock_movemap
( core::pose::Pose const & pose,
	core::select::residue_selector::ChainSelectorCOP glycolig_selector,
	bool const lock_rings,
	core::Size const dock_jump_num )
{
	TR.Debug << "Setting up MoveMap for GlycanDock minimization" << std::endl;

	///////////////////////
	// MINIMIZER MOVEMAP //
	// Glycan bb+sc      //
	// Protein sc        //
	// Interface Jump    //
	///////////////////////

	//// Glycan bb+sc
	// First, set the MoveMap to respect the backbone and IUPAC chi of the glycan
	// TODO future development
	// Default is True for glycan branch torsions. Should this be settable?
	core::kinematics::MoveMapOP out_minimizer_mm =
		core::pose::carbohydrates::create_glycan_movemap_from_residue_selector
		( pose,
		glycolig_selector,
		true, // include_iupac_chi
		// if lock_rings == true, do NOT include ring torsions
		! lock_rings, // include_glycan_ring_torsions
		true, // include_bb_torsions
		false ); // cartesian

	//// Receptor sc
	// Set chi for all residues not specified to be the glycoligand
	// Note that the protein partner may have conjugated glycans attached to it
	// and the blanket 'set_chi' functionality in MoveMap may NOT
	// respect carbohydrate-specific chis
	utility::vector1< bool > glycolig_subset = glycolig_selector->apply( pose );
	for ( core::Size ii( 1 ); ii <= pose.size(); ++ii ) {
		if ( ! glycolig_subset[ii] ) {
			if ( pose.residue( ii ).is_carbohydrate() ) {
				// If this is a carb residue that is not part of the glycolig_subset
				// (like a glycan conjugated to the protein receptor partner)
				// Turn on its chis that are NOT also backbone torsions
				core::pose::carbohydrates::set_glycan_iupac_chi_torsions
					( pose,
					*out_minimizer_mm,
					ii,
					true, // set chi
					false); // cartesian
			} else {
				// Otherwise, this should be a protein residue. Set chis as normal
				// But do NOT allow the chis of the protein residues (like ASN)
				// that are conjugated (to a glycan) to minimize (moves entire glycan)
				if ( pose.residue(ii).is_branch_point() ) { continue; }
				out_minimizer_mm->set_chi( ii, true );
			}
		}
	}

	//// Interface Jump
	// Allow for the dock_jump_num_ to be flexible as well
	out_minimizer_mm->set_jump( false ); // ensure all other Jumps are false
	out_minimizer_mm->set_jump( dock_jump_num, true );

	// Override MoveMap given user's custom MoveMap
	if ( basic::options::option[ basic::options::OptionKeys::in::file::movemap ].user() ) {
		TR.Debug << "Using user-provided MoveMap file" << std::endl;
		out_minimizer_mm->init_from_file
			(basic::options::option[ basic::options::OptionKeys::in::file::movemap ]);
	}

	// Don't need to debug movemap here because debug output of the
	// minimizer will also show the movemap

	return out_minimizer_mm;

} // END setup_GlycanDock_movemap


//--------------------------------------------------
protocols::moves::SequenceMoverOP
setup_GlycanDock_stage1_rb_seqmover
( core::pose::Pose const & pose,
	core::Size const docking_jump_num,
	core::Real const stage1_rotation_magnitude,
	core::Real const stage1_translation_magnitude,
	bool const stage1_rotate_glycolig_about_com )
{

	TR.Debug << "Setting up rigid-body perturbation Movers for "
		"Stage 1 GlycanDock conformation initialization" << std::endl;

	protocols::moves::SequenceMoverOP out_stage1_rb_seqmover =
		utility::pointer::make_shared< protocols::moves::SequenceMover >();

	// For perturbing the glycoligand's center-of-mass using both
	// a translational and rotational perturbation
	// Default functionality is to use gaussian moves for CoM perturbation
	out_stage1_rb_seqmover->add_mover
		( utility::pointer::make_shared
		< protocols::rigid::RigidBodyPerturbMover >
		( docking_jump_num,
		stage1_rotation_magnitude, stage1_translation_magnitude,
		protocols::rigid::partner_downstream, false /*interface_in*/) );

	// For randomizing the glycoligand's rigid-body position in 360 degree space
	// Much more aggressive rotational perturbation of glycoligand CoM
	// Should only use this if you are not at all confident in the initial
	// rigid-body 360 degree orientation of the glycoligand in the binding pocket
	if ( stage1_rotate_glycolig_about_com ) {
		TR.Debug << "Adding a RigidBodyRandomizeMover to the Stage 1 "
			"rigid-body conformation initialization mover list" << std::endl;
		out_stage1_rb_seqmover->add_mover
			( utility::pointer::make_shared
			< protocols::rigid::RigidBodyRandomizeMover >
			// pose, docking Jump, partner up or downstream (down = glycan = 2),
			// phi angle, psi angle, update_center_after_move
			( pose, docking_jump_num,
			protocols::rigid::partner_downstream, 360, 360, false) );
	}

	TR.Debug << "Finished setting up rigid-body perturbation Movers for "
		"Stage 1 GlycanDock conformation initialization" << std::endl;

	return out_stage1_rb_seqmover;

} // END setup_GlycanDock_stage1_rb_seqmover


//--------------------------------------------------
protocols::moves::RandomMoverOP
setup_GlycanDock_stage2_rb_randmover
( core::Size const docking_jump_num,
	core::Real const stage2_rotation_magnitude,
	core::Real const stage2_translation_magnitude,
	bool const slide_glycolig_into_contact )
{

	TR.Debug << "Setting up rigid-body perturbation Movers for "
		"Stage 2 GlycanDock sampling/docking moves" << std::endl;

	protocols::moves::RandomMoverOP out_stage2_rb_randmover =
		utility::pointer::make_shared< protocols::moves::RandomMover >();

	//////////////////////////////
	// RIGID BODY MOVER WEIGHTS //
	//////////////////////////////
	// Default behavior; occasionally slide glycan into contact with protein
	core::Real rb_perturber_weight = 0.67; // Default 2/3 of the time
	core::Real rb_perturb_with_slide_weight = 0.33; // Default 1/3 of the time
	if ( ! slide_glycolig_into_contact ) {
		// If set to false, we are only doing RB perturbations
		// Perhaps the sliding gives the user weird behavior?
		rb_perturber_weight = 1.0;
		rb_perturb_with_slide_weight = 0.0;
	}

	///////////////////////
	// RIGID BODY MOVERS //
	///////////////////////
	// Given the docking_jump_num (1), perturb the downstream partner
	// (the glycoligand) a max of stage2_rot_mag and stage2_trans_mag in space
	protocols::rigid::RigidBodyPerturbMoverOP rb_perturber =
		utility::pointer::make_shared
		< protocols::rigid::RigidBodyPerturbMover >
		( docking_jump_num,
		stage2_rotation_magnitude, stage2_translation_magnitude,
		protocols::rigid::partner_downstream, false /*interface_in*/ );
	// Add it to the weighted random mover at the given weight
	out_stage2_rb_randmover->add_mover( rb_perturber, rb_perturber_weight );

	if ( slide_glycolig_into_contact ) {
		// Perform the same perturbation as above, but follow the perturbation
		// up with a slide into contact Mover
		// This sequence is to ensure that, occassionaly,
		// the glycoligand is slide back toward its protein partner
		// in the cases that it drifts too far away and any
		// glycoligand constraints to the protein aren't enough
		// Create the SequenceMover using the above rb_perturber and slide_into_contact
		// Instantiating with 'true' as this sets the SequenceMover
		// to use the Mover status to control moves
		protocols::moves::SequenceMoverOP
			rb_perturb_with_slide_into_contact =
			utility::pointer::make_shared
			< protocols::moves::SequenceMover >( true );
		rb_perturb_with_slide_into_contact->add_mover( rb_perturber );
		rb_perturb_with_slide_into_contact->add_mover
			( utility::pointer::make_shared
			< protocols::docking::FaDockingSlideIntoContact >
			( docking_jump_num ) );
		// Add the SequenceMover to the weighted random mover at the given weight
		// if slide_glycolig_into_contact = false, the weight is 0.0
		out_stage2_rb_randmover->add_mover
			( rb_perturb_with_slide_into_contact, rb_perturb_with_slide_weight );
	}

	// TODO future development
	// any other RB moves to employ? diff sequences, amounts, weights?

	TR.Debug << "Finished setting up rigid-body perturbation Movers for "
		"Stage 2 GlycanDock sampling/docking moves" << std::endl;

	return out_stage2_rb_randmover;

} // END setup_GlycanDock_stage2_rb_randmover


//--------------------------------------------------
protocols::moves::RandomMoverOP
setup_GlycanDock_stage2_tor_randmover
( core::pose::Pose const & pose,
	utility::vector1< bool > const & glycolig_subset_with_torsions,
	core::scoring::ScoreFunctionOP sf,
	core::Real const mc_kt,
	core::Size const n_shear_moves,
	bool const refine_only )
{

	using namespace protocols::carbohydrates;
	using namespace protocols::simple_moves;
	using namespace protocols::simple_moves::bb_sampler;

	TR.Debug << "Setting up glycosidic torsion angle perturbation Movers "
		"for Stage 2 GlycanDock sampling/docking moves" << std::endl;

	// Emulate the setup of GlycanSampler
	// NOTE
	// Not using GlycanSampler because:
	// 1) it uses packing and minimization as potential Movers, which
	// I do not want because I want to switch between these Movers
	// 2) it uses LinkageConformerMover, which I don't think should
	// be applied during a docking protocol
	protocols::moves::RandomMoverOP out_stage2_tor_randmover =
		utility::pointer::make_shared< protocols::moves::RandomMover >();

	// Given the filtered glycolig_subset_with_torsions,
	// make a local selector that returns that subset
	core::select::residue_selector::ReturnResidueSubsetSelectorOP
		glycolig_stor_with_torsions =
		utility::pointer::make_shared
		< core::select::residue_selector::ReturnResidueSubsetSelector >
		( glycolig_subset_with_torsions );

	// NOTE dihedral mask setup ONLY uses glycan residues that
	// are part of the glycolig_subset_with_torsions
	// Dihedral masks tell samplers which bb torsions each residue actually has
	// < pose residue number, < list of torsion IDs > >
	std::map< core::Size, utility::vector1< core::Size >> dihedral_mask;
	// Also take note of the max glycosidic dihedrals in the glycan
	// At minimum should have phi (1) and psi (2), but may have an omega (3)
	core::Size max_torsion_id = 0;
	for ( core::Size resnum = 1; resnum <= pose.total_residue(); ++resnum ) {
		if ( ! glycolig_subset_with_torsions[resnum] ) { continue; }
		dihedral_mask[ resnum ];
		// 1 : phi, 2 : psi, 3 : omega, 4 : omega2
		core::Size n_glycosidic_torsions
			( core::conformation::carbohydrates::get_n_glycosidic_torsions_in_res
			( pose.conformation(), resnum ) );
		for ( core::Size torsion_id = 1;
				torsion_id <= n_glycosidic_torsions; ++torsion_id ) {
			dihedral_mask[ resnum ].push_back( torsion_id );
			// Update the max_torsion_id seen among all glycosidic torsions
			if ( torsion_id > max_torsion_id ) { max_torsion_id = torsion_id; }
		}
	} // end dihedral mask setup

	///////////////////////////
	// TORSION MOVER WEIGHTS //
	///////////////////////////
	// Again, based on GlycanSampler
	core::Real const sugarbb_sampler_weight( 0.45 );
	// SmallBBSampler is the only Mover that will never fail, of these choices
	// IT IS THEREFORE NOT RECOMMENDED TO REMOVE THE SmallBBSampler
	core::Real const total_smallbb_weight( 0.3 );
	// SmallBBSamplerMover splits weights by size of sampler
	// small size moves (+/- 30 degrees): 0.171429
	// medium size moves (+/- 90 degrees): 0.0857143
	// large size moves (+/- 180 degrees): 0.0428571
	core::Real const shear_mover_weight( 0.2 );
	core::Real const rpf_mover_weight( 0.05 );

	////////////////////
	// SugarBBSampler //
	////////////////////
	// SugarBBSampler uses sugar_bb data to set the specified torsion
	// If set in dihedral mask, want to sample all possible torsions per residue
	// Default args are to use prob for BBSampleType with sampling_step_size of 0.1
	SugarBBSamplerOP phi_sugarbb_sampler =
		utility::pointer::make_shared< SugarBBSampler >
		( core::id::phi_dihedral );
	SugarBBSamplerOP psi_sugarbb_sampler =
		utility::pointer::make_shared< SugarBBSampler >
		( core::id::psi_dihedral );
	SugarBBSamplerOP omega_sugarbb_sampler =
		utility::pointer::make_shared< SugarBBSampler >
		( core::id::omega_dihedral );

	// Add the SugarBBSamplers to a unified sugarbb_sampler_mover
	BBDihedralSamplerMoverOP sugarbb_sampler_mover =
		utility::pointer::make_shared< BBDihedralSamplerMover >();
	sugarbb_sampler_mover->add_sampler( phi_sugarbb_sampler );
	sugarbb_sampler_mover->add_sampler( psi_sugarbb_sampler );
	sugarbb_sampler_mover->add_sampler( omega_sugarbb_sampler );
	// And set the local glycolig_stor_with_torsions that has
	// been modified to ignore glycans with no dihedrals (sugars with no parent)
	sugarbb_sampler_mover->set_residue_selector( glycolig_stor_with_torsions );
	// And set the dihedral mask that clarifies which carbohydrate
	// residues have which torsions available to it
	// BBDihedralSampler has phi, psi, and omega set, but not every
	// carbohydrate will have all three of those torsions available
	sugarbb_sampler_mover->set_dihedral_mask( dihedral_mask );

	// Add SugarBBSampler to the Mover container
	out_stage2_tor_randmover->add_mover
		( sugarbb_sampler_mover , sugarbb_sampler_weight );

	////////////////////
	// SmallBBSampler //
	////////////////////
	// Note that this is different than the SmallMover of BackboneMover
	BBDihedralSamplerMoverOP small_size_bb_mover =
		utility::pointer::make_shared< BBDihedralSamplerMover >();
	BBDihedralSamplerMoverOP medium_size_bb_mover =
		utility::pointer::make_shared< BBDihedralSamplerMover >();
	BBDihedralSamplerMoverOP large_size_bb_mover =
		utility::pointer::make_shared< BBDihedralSamplerMover >();
	// Set the residue selectors
	small_size_bb_mover->set_residue_selector( glycolig_stor_with_torsions );
	medium_size_bb_mover->set_residue_selector( glycolig_stor_with_torsions );
	large_size_bb_mover->set_residue_selector( glycolig_stor_with_torsions );
	// And add the dihedral masks
	small_size_bb_mover->set_dihedral_mask( dihedral_mask );
	medium_size_bb_mover->set_dihedral_mask( dihedral_mask );
	large_size_bb_mover->set_dihedral_mask( dihedral_mask );
	// Make small, medium, and large SmallMovers per available torsion
	// TODO QUESTION
	// wouldn't the dihedral masks handle this?
	for ( core::Size ii(1); ii <= max_torsion_id; ++ii ) {
		auto dih_type = static_cast< core::id::MainchainTorsionType >( ii );

		SmallBBSamplerOP small_size_bb_sampler =
			utility::pointer::make_shared< SmallBBSampler >
			( dih_type, 30 ); // +/- 15
		SmallBBSamplerOP medium_size_bb_sampler =
			utility::pointer::make_shared< SmallBBSampler >
			( dih_type, 90 ); // +/- 45
		SmallBBSamplerOP large_size_bb_sampler =
			utility::pointer::make_shared< SmallBBSampler >
			( dih_type, 180); // +/- 90

		small_size_bb_mover->add_sampler( small_size_bb_sampler );
		medium_size_bb_mover->add_sampler( medium_size_bb_sampler );
		large_size_bb_mover->add_sampler( large_size_bb_sampler );
	}

	// Taken from @jadolfbr GlycanSampler
	// Add Small bb samplers to the Mover container based on a given ratio
	//  4:2:1 ratio of sampling Small,Medium,Large
	//  .3 = 4X (small) + 2X (medium) + X(large)
	constexpr core::Real const expr = 0.3/(1+2+4);
	core::Real const small_bb_weight = 4*expr;
	core::Real const medium_bb_weight = 2*expr;
	core::Real const large_bb_weight = expr;

	// If doing refine_only, do not make huge random movements
	if ( ! refine_only ) {
		out_stage2_tor_randmover->add_mover
			( small_size_bb_mover, small_bb_weight );
		out_stage2_tor_randmover->add_mover
			( medium_size_bb_mover, medium_bb_weight );
		out_stage2_tor_randmover->add_mover
			( large_size_bb_mover, large_bb_weight );
	} else {
		// Doing refine_only, make small backbone moves only
		out_stage2_tor_randmover->add_mover
			( small_size_bb_mover, total_smallbb_weight );
	}

	////////////////
	// ShearMover //
	////////////////
	// temperature is needed for ShearMover to make a move with
	// rama if no sf given, but will be given
	// n_shear_moves is default to 1 in ShearMover code and here
	ShearMoverOP shear_mover =
		utility::pointer::make_shared< ShearMover >();
	shear_mover->set_residue_selector( glycolig_stor_with_torsions );
	shear_mover->nmoves( n_shear_moves ); // default = 1
	shear_mover->temperature( mc_kt ); // default = 0.6
	shear_mover->scorefxn( sf );
	// Add ShearMover to the Mover container
	out_stage2_tor_randmover->add_mover( shear_mover, shear_mover_weight );

	////////////////////////
	// RingPlaneFlipMover //
	////////////////////////
	RingPlaneFlipMoverOP rpf_mover
		( utility::pointer::make_shared< RingPlaneFlipMover >() );
	rpf_mover->selector( glycolig_stor_with_torsions );
	// Add RingPlaneFlipMover to the Mover container
	out_stage2_tor_randmover->add_mover( rpf_mover, rpf_mover_weight );

	TR.Debug << "Finished setting up glycosidic torsion angle perturbation Movers "
		"for Stage 2 GlycanDock sampling/docking moves" << std::endl;

	return out_stage2_tor_randmover;

} // END setup_GlycanDock_stage2_tor_randmover


//--------------------------------------------------
void
do_stage1_tor_uniform_perturbation
( core::pose::Pose & pose,
	utility::vector1< bool > const & glycolig_subset,
	core::Real const stage1_tor_perturb_mag )
{

	using namespace numeric::random; // uniform()

	for ( core::Size resnum( 1 ); resnum <= glycolig_subset.size(); ++resnum ) {
		if ( ! glycolig_subset[ resnum ] ) { continue; }
		// Uniformally perturb phi, psi, (and omega) glycosidic torsion angles
		// First, create the perturbed values for each dihedral
		// Perturbation =  initial torsion angle +/- (uniform([0, 1]) * magnitude)
		// NOTE omega2 not available to set directly via the Pose
		// Would need the util function set_glycosidic_torsion to set omega2
		core::Real const phi_pert =
			uniform() * stage1_tor_perturb_mag * 2 - stage1_tor_perturb_mag;
		core::Real const psi_pert =
			uniform() * stage1_tor_perturb_mag * 2 - stage1_tor_perturb_mag;
		core::Real const omega_pert =
			uniform() * stage1_tor_perturb_mag * 2 - stage1_tor_perturb_mag;
		// Now, set each dihedral to its perturbed torsion angle
		// If residue does not have that torsion (i.e. omega), nothing happens
		pose.set_phi
			( resnum,
			numeric::principal_angle_degrees( pose.phi(resnum) + phi_pert ) );
		pose.set_psi
			( resnum,
			numeric::principal_angle_degrees( pose.psi(resnum) + psi_pert ) );
		pose.set_omega
			( resnum,
			numeric::principal_angle_degrees( pose.omega(resnum) + omega_pert ) );
	}

} // END do_stage1_tor_uniform_perturbation


//--------------------------------------------------
void
do_stage1_conformation_initialization
( core::pose::Pose & pose,
	protocols::moves::SequenceMoverOP stage1_rb_seqmover,
	/* protocols::moves::RandomMoverOP stage1_tor_mover TODO */
	utility::vector1< bool > const & glycolig_subset,
	bool const glycolig_has_dihedrals,
	core::Real const stage1_tor_perturb_mag )
{
	// TODO future development
	// this method will use input stage1_tor_mover when that can exist
	// Until then, perform uniform torsion angle perturbation manually

	debug_assert( stage1_rb_seqmover );

	TR << "Performing Stage 1 GlycanDock "
		"conformation initialization" << std::endl;

	// If glycoligand has glycosidic torsion angles to perturb, perturb them
	if ( glycolig_has_dihedrals ) {
		TR.Debug << "Manually performing Stage 1 glycosidic "
			"torsion angle perturbation" << std::endl;
		do_stage1_tor_uniform_perturbation( pose,
			glycolig_subset, stage1_tor_perturb_mag );
	} else {
		TR << "Glycoligand did not have any detected glycosidic torsion angles! "
			"Cannot perform any Stage 1 torsion angle perturbation" << std::endl;
	}

	TR << "Finished performing Stage 1 GlycanDock "
		"conformation initialization" << std::endl;

} // END do_stage1_conformation_initialization


//--------------------------------------------------
core::Size
do_stage2_sample_pack_min_cycle
( core::pose::Pose & pose,
	protocols::moves::RandomMoverOP stage2_randmover,
	core::Size const n_rounds,
	minimization_packing::PackRotamersMoverOP full_packer,
	minimization_packing::EnergyCutRotamerTrialsMoverOP ecut_packer,
	core::Size const full_pack_every_x_rounds,
	minimization_packing::MinMoverOP minimizer,
	moves::MonteCarloOP mc )
{

	// if doing anything other than using the full PackRotamersMover
	// every round, assert that the ecut_packer is not a nullptr
	if ( full_pack_every_x_rounds != 1 ) {
		debug_assert( ecut_packer );
	}

	core::Size n_moves_accepted( 0 ); // returned
	for ( core::Size rnd( 1 ); rnd <= n_rounds; ++rnd ) {
		// APPLY RANDOM SAMPLING MOVER
		stage2_randmover->apply( pose );

		// CHECK MOVER APPLY STATUS (AND REDO IF NEEDED)
		// If the randomly chosen Mover did not work
		// try another Mover at least 100 times until one does
		// all rigid-body Movers will be successfully applied,
		// but some torsion Movers may not be accepted on a glycoligand
		// (e.g. RingPlaneFlipMover and ShearMover require special circumstances)
		if ( stage2_randmover->get_last_move_status() !=
				protocols::moves::MS_SUCCESS ) {
			TR.Info << "Mover type " << stage2_randmover->type() <<
				" didn't work, trying another" << std::endl;
			for ( core::Size ii( 1 ); ii <= 100; ++ii ) {
				stage2_randmover->apply( pose );
				if ( stage2_randmover->get_last_move_status() ==
						protocols::moves::MS_SUCCESS ) {
					TR.Info << stage2_randmover->type() << " worked" << std::endl;
					break;
				}
			}
		}

		// APPLY PACKER
		// Pack at interface after glycoligand has been perturbed
		// Default is to use full PackRotamersMover every 8 inner cycles
		// Can use only EnergyCutRotamerTrialsMover if full_pack_every_x_rounds == 0
		if ( ( full_pack_every_x_rounds != 0 ) &&
				(( rnd % full_pack_every_x_rounds == 0 ) ||
				( rnd == full_pack_every_x_rounds )) ) {
			full_packer->apply( pose );
			//if ( TR.Debug.visible() ) {
			//TR.Debug << "PackRotamersMover task used" << std::endl;
			//packer_tf_->create_task_and_apply_taskoperations(pose)->show( TR );
			//}
		} else {
			ecut_packer->apply( pose );
		}

		// APPLY MINIMIZER
		// TODO future development and benchmarking
		// is the minimization quick? Should it happen more or less frequently?
		if ( rnd % 2 == 0 || rnd == n_rounds ) {
			minimizer->apply( pose );
		}

		// APPLY MONTECARLO METROPOLIS CRITERION
		const bool accepted = mc->boltzmann( pose );
		if ( accepted ) {
			++n_moves_accepted;
		}
	}

	return n_moves_accepted;

} // END do_stage2_sample_pack_min_cycle


//--------------------------------------------------
bool
calc_res_contact
( core::conformation::ResidueCOP rsd_ii,
	core::conformation::ResidueCOP rsd_jj,
	core::Real const dist_cutoff )
{
	core::Real dist_cutoff_sq = dist_cutoff * dist_cutoff;

	for ( core::Size at_ii = 1; at_ii <= rsd_ii->natoms(); ++at_ii ) {
		if ( rsd_ii->atom_is_hydrogen( at_ii ) ) { continue; }
		for ( core::Size at_jj = 1; at_jj <= rsd_jj->natoms(); ++at_jj ) {
			if ( rsd_jj->atom_is_hydrogen( at_jj ) ) { continue; }
			// if here, atom is not a hydrogen, but it could be a virtual atom
			// since virtual atoms in carbs track real atoms, this should be okay
			core::Real dist_sq =
				rsd_ii->xyz( at_ii ).distance_squared( rsd_jj->xyz( at_jj ) );
			if ( dist_sq <= dist_cutoff_sq ) { return true; }
		}
	}

	return false;

} // END calc_res_contact


//--------------------------------------------------
ObjexxFCL::FArray2D_bool
gen_intf_contact_map
( core::pose::Pose const & pose,
	utility::vector1< bool > const & receptor_nbrhood_subset,
	utility::vector1< bool > const & glycolig_subset,
	core::Real const dist_cutoff )
{
	debug_assert( receptor_nbrhood_subset.size() == glycolig_subset.size() );

	ObjexxFCL::FArray2D_bool contact_map
		( receptor_nbrhood_subset.size(), glycolig_subset.size(), false );

	for ( core::Size ii = 1; ii <= receptor_nbrhood_subset.size(); ++ii ) {
		if ( ! receptor_nbrhood_subset[ii] ) { continue; }
		for ( core::Size jj = 1; jj <= glycolig_subset.size(); ++jj ) {
			if ( ! glycolig_subset[jj] ) { continue; }
			// (protein) receptor residue
			core::conformation::ResidueCOP rsd_ii
				( utility::pointer::make_shared
				< core::conformation::Residue >( pose.residue( ii ) ) );
			// glycoligand residue
			core::conformation::ResidueCOP rsd_jj
				( utility::pointer::make_shared
				< core::conformation::Residue >( pose.residue( jj ) ) );
			// check each heavy atom distance to each glycolig heavy atom
			contact_map( ii, jj ) =
				calc_res_contact( rsd_ii, rsd_jj, dist_cutoff );
		}
	}

	return contact_map;

} // END gen_intf_contact_map


//--------------------------------------------------
utility::vector1< core::Real >
calc_GlycanDock_intf_metrics
( core::pose::Pose & pose, // writing data to pose object
	core::pose::Pose const & ref_pose,
	utility::vector1< bool > const & glycolig_subset,
	core::Real const dist_cutoff /* = 5 Ang */)
{
	using namespace core::select::residue_selector;

	// not pretty, but it works
	bool skip_ref_comparison = false;
	if ( pose.size() != ref_pose.size() ) {
		TR.Warning << "Pose size does not match the reference Pose size! "
			"Cannot calculate interface metrics with confidence. "
			"Returning bogus values for native interface residue and contact "
			"recoveries as well as for Fnat" << std::endl;
		skip_ref_comparison = true;
	}

	// Get the residues within 10A of the glycoligand
	// using 10A as this selector can use the 10A neighbor graph
	NeighborhoodResidueSelector neighborhood_stor =
		NeighborhoodResidueSelector();
	neighborhood_stor.set_focus( glycolig_subset );
	neighborhood_stor.set_distance( 10.0 ); // tenA neighbor graph
	// use the glycoligand to get neighbor residues, but
	// do not include the glycoligand in the neighborhood subset
	neighborhood_stor.set_include_focus_in_subset( false );

	////////////////////
	// Get contact map info for the decoy pose
	// filters down 10A energy graph down to 5A heavy-atom contacts
	utility::vector1< bool > const decoy_nbrhood_subset =
		neighborhood_stor.apply( pose );
	ObjexxFCL::FArray2D_bool const decoy_contact_map =
		gen_intf_contact_map( pose, decoy_nbrhood_subset,
		glycolig_subset, dist_cutoff );

	////////////////////
	// setup counters and subsets for interface metrics
	core::Real decoy_intf_contacts = 0.0; // ALL unique res-res contact pairs
	core::Real decoy_nat_intf_contacts = 0.0; // NATIVE unique res-res contact pairs
	utility::vector1< bool > decoy_intf_res_subset = // ALL unique intf res
		( utility::pointer::make_shared
		< FalseResidueSelector >() )->apply( pose );
	utility::vector1< bool > decoy_nat_intf_res_subset = // NATIVE unique intf res
		( utility::pointer::make_shared
		< FalseResidueSelector >() )->apply( pose );

	// Count interfacial residue contacts and the number of interface residues
	for ( core::Size ii = 1; ii <= pose.size(); ++ii ) {
		if ( ! decoy_nbrhood_subset[ii] ) { continue; }
		for ( core::Size jj = 1; jj <= pose.size(); ++jj ) {
			if ( ! glycolig_subset[jj] ) { continue; }
			// Check the decoy contact pair
			if ( decoy_contact_map( ii, jj ) == true ) {
				decoy_intf_contacts += 1.0;
				decoy_intf_res_subset[ii] = true;
				decoy_intf_res_subset[jj] = true;
			} // if this is a contact
		} // for every glycolig resnum
	} // for every neighbor resnum (does not include any glycolig resnums)

	// count the number of interface residues in the decoy
	// and return a pymol selection of those decoy interface residues
	core::Real const n_intf_residues = // yes, change Size to Real
		core::select::residue_selector::count_selected( decoy_intf_res_subset );
	//core::simple_metrics::metrics::SelectedResiduesPyMOLMetric SRPyMOL =
	// core::simple_metrics::metrics::SelectedResiduesPyMOLMetric();
	//SRPyMOL.set_residue_selector( core::select::get_residue_selector_from_subset
	// ( decoy_intf_res_subset ) );
	//core::pose::setPoseExtraScore( pose, "pymol_intf_res_5A_cutoff",
	// SRPyMOL.calculate( pose ) );

	// If we are not comparing the decoy to the reference pose,
	// then we are done getting interface stats
	// The rest (like Fnat) will be returned as -1
	if ( skip_ref_comparison ) {
		// return data like so
		// < n_intf_res, n_nat_intf_res, Fnat_intf_res,
		//   n_intf_contacts, n_nat_intf_contacts, Fnat >
		utility::vector1< core::Real > intf_data;
		intf_data.push_back( n_intf_residues );
		intf_data.push_back( -1.0 ); // n_nat_intf_res
		intf_data.push_back( -1.0 ); // Fnat_intf_res
		intf_data.push_back( decoy_intf_contacts );
		intf_data.push_back( -1.0 ); // decoy_nat_intf_contacts
		intf_data.push_back( -1.0 ); // Fnat

		return intf_data;

	} else {
		// We can compare the decoy to the reference pose
		// i.e. pose.size() == ref_pose.size()
		////////////////////
		// First, get the info needed from the reference pose
		// get the contact map for the reference pose
		utility::vector1< bool > const ref_nbrhood_subset =
			neighborhood_stor.apply( ref_pose );
		ObjexxFCL::FArray2D_bool const ref_contact_map =
			gen_intf_contact_map( ref_pose,
			ref_nbrhood_subset, // like partner1, but trimmed a bit
			glycolig_subset, dist_cutoff );
		core::Real ref_intf_contacts = 0.0; // unique res-res contact pairs
		utility::vector1< bool > ref_intf_res_subset = // unique intf residues
			( utility::pointer::make_shared
			< FalseResidueSelector >() )->apply( pose );

		for ( core::Size ii = 1; ii <= pose.size(); ++ii ) {
			for ( core::Size jj = 1; jj <= pose.size(); ++jj ) {
				if ( ! glycolig_subset[jj] ) { continue; }
				// if this contact pair is in the ref_contact_map,
				// it is a native intf contact and are both native intf residues
				if ( ref_contact_map( ii, jj ) == true ) {
					ref_intf_contacts += 1.0;
					ref_intf_res_subset[ii] = true;
					ref_intf_res_subset[jj] = true;
					// now we can check if the decoy recovered this native contact
					if ( decoy_contact_map( ii, jj ) == true ) {
						if ( TR.Debug.visible() ) {
							TR.Debug << pose.pdb_info()->pose2pdb( ii ) <<
								" to " << pose.pdb_info()->pose2pdb( jj ) <<
								" is a recovered native contact " << std::endl;
						}
						decoy_nat_intf_contacts += 1.0;
						decoy_nat_intf_res_subset[ii] = true;
						decoy_nat_intf_res_subset[jj] = true;
					} // endif decoy recovered native contact
				} // endif this is a native contact
			} // for every glycolig resnum
		} // for every resnum (does not include any glycolig resnums)

		// these metrics we can only calculate if we can compare to the ref_pose
		core::Real const Fnat = decoy_nat_intf_contacts / ref_intf_contacts;
		core::Real const n_nat_intf_res =
			core::select::residue_selector::count_selected( decoy_nat_intf_res_subset );
		// internet says float / int is okay (i.e. returns a float)
		core::Real const Fnat_intf_res = n_nat_intf_res /
			core::select::residue_selector::count_selected( ref_intf_res_subset );

		// return data like so
		// < n_intf_res, n_nat_intf_res, Fnat_intf_res,
		//   n_intf_contacts, n_nat_intf_contacts, Fnat >
		utility::vector1< core::Real > intf_data;
		intf_data.push_back( n_intf_residues );
		intf_data.push_back( n_nat_intf_res );
		intf_data.push_back( Fnat_intf_res );
		intf_data.push_back( decoy_intf_contacts );
		intf_data.push_back( decoy_nat_intf_contacts );
		intf_data.push_back( Fnat );

		return intf_data;

	} // end if ( skip_ref_comparison )

	// Should not get here
	// hence, no return statement right here

} // END calc_GlycanDock_Fnat


//--------------------------------------------------
protocols::analysis::InterfaceAnalyzerMoverOP
get_GlycanDock_IAM(
	core::pose::DockingPartners const & docking_partners,
	core::scoring::ScoreFunctionOP sf
) {
	protocols::analysis::InterfaceAnalyzerMoverOP IAM
		( utility::pointer::make_shared
		< protocols::analysis::InterfaceAnalyzerMover >() );

	// By setting the interface like LH_X instead of using dock_jump_num_
	// we instantiate InterfaceAnalyzerMover using an explicit constructor
	// which will work correctly for Poses with > 1 upstream chains
	// compared to the 1 downstream glycan chain
	IAM->set_interface( docking_partners );

	// By now sf should be adjusted as needed for glycoligand docking
	IAM->set_scorefunction( sf->clone() );

	// Turn off un-needed computations
	// no centroid capability with glycans at the moment
	IAM->set_use_centroid_dG( false );
	// do not change the input pose
	IAM->set_pack_input( false );

	// Turn on desired computations
	// we calc intE separately (w/ no packing), so get dG (w/ packing)
	IAM->set_pack_separated( true );
	// Lawrence and Coleman Shape Complementarity score
	IAM->set_compute_interface_sc( true );
	// costly calculation getting packing efficiency at interface
	IAM->set_compute_packstat( true );

	// We want to use the tracer to write osstreams
	// like the pymol info to the pose
	IAM->set_use_tracer( false );
	// but since we are using the apply (vs apply_const) method,
	// we do not want to write the same data twice-so skip reporting
	IAM->set_skip_reporting( true );

	return IAM;

} // END get_GlycanDock_IAM

} // glycan_docking
} // protocols
