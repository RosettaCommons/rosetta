// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/ligand_docking/LigandDockProtocol.cc
///
/// @brief
/// @author Ian W. Davis


#include <protocols/ligand_docking/LigandDockProtocol.hh>

#include <core/types.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/graph/Graph.hh>
#include <core/grid/CartGrid.hh>
#include <core/kinematics/Edge.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <basic/options/option.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/sasa.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/dunbrack/RotamerConstraint.hh>
#include <basic/prof.hh>
#include <basic/Tracer.hh>
#include <basic/MetricValue.hh>

//#include <protocols/relax_protocols.hh>
#include <protocols/docking/DockingInitialPerturbation.hh>
#include <protocols/rigid/RB_geometry.hh>
#include <protocols/ligand_docking/grid_functions.hh>
#include <protocols/ligand_docking/ligand_functions.hh>
#include <protocols/ligand_docking/LigandBaseProtocol.hh>
#include <protocols/ligand_docking/RandomConformerMover.hh>
#include <protocols/ligand_docking/ResidueTorsionRestraints.hh>
#include <protocols/ligand_docking/UnconstrainedTorsionsMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/JumpOutMover.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>


#include <numeric/conversions.hh>
#include <numeric/random/random.hh>
#include <numeric/xyzVector.io.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/string.functions.hh>
#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/file/FileName.hh>

#include <utility/pointer/owning_ptr.hh>
#include <utility/tools/make_vector1.hh>

#include <fstream>
#include <sstream>


// option key includes

#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <basic/options/keys/enzdes.OptionKeys.gen.hh>

#include <core/pose/Pose.hh>
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh>
#include <utility/vector0.hh>


namespace protocols {
namespace ligand_docking {


static THREAD_LOCAL basic::Tracer TR( "protocols.ligand_docking.LigandDockProtocol" );


LigandDockProtocol::LigandDockProtocol():
	//utility::pointer::ReferenceCount(),
	protocols::moves::Mover(),
	LigandBaseProtocol(),
	start_from_pts_(),
	ligand_torsion_restraints_()
{
	Mover::type( "LigandDockProtocol" );

	using basic::options::option;
	using namespace basic::options;

	// options
	protocol_           = option[ OptionKeys::docking::ligand::protocol ];
	minimize_ligand_    = option[ OptionKeys::docking::ligand::minimize_ligand ];
	minimize_backbone_  = option[ OptionKeys::docking::ligand::minimize_backbone ];
	tether_ligand_      = option[ OptionKeys::docking::ligand::tether_ligand ].user();
	minimize_all_rsds_  = false;
	repack_all_rsds_    = false;
	rottrials_all_rsds_ = false;
	ligand_protonation_ = option[ OptionKeys::docking::ligand::mutate_same_name3 ];
	minimize_water_     = true;
	sc_interface_padding_  = 0.0;//5.0; // 5A ends up repacking half the protein! (literally)
	bb_interface_cutoff_   = 7.0; //5.0;
	ligand_chi_stddev_deg_ = option[ OptionKeys::docking::ligand::harmonic_torsions ];
	protein_CA_stddev_Ang_ = option[ OptionKeys::docking::ligand::harmonic_Calphas ];
	ligand_shear_moves_    = (core::Size) option[ OptionKeys::docking::ligand::shear_moves ];
	ligand_tether_stddev_Ang_ = (tether_ligand_ ? option[ OptionKeys::docking::ligand::tether_ligand ]() : -1.0);

	if ( option[ OptionKeys::docking::ligand::start_from ].user() ) {
		utility::vector1< core::Real > start_from = option[ OptionKeys::docking::ligand::start_from ]();
		// Make sure a whole number of triples were supplied
		if ( start_from.size() % 3 != 0 ) {
			utility_exit_with_message("-start_from requires one or more X,Y,Z triples -- you didn't provide enough numbers!");
		}
		for ( Size ii = 1; ii <= start_from.size(); ii += 3 ) {
			start_from_pts_.push_back( core::Vector(start_from[ii], start_from[ii+1], start_from[ii+2]) );
		}
	}
	if ( ligand_shear_moves_ && !basic::options::option[ basic::options::OptionKeys::docking::ligand::use_ambig_constraints ]() ) {
		utility_exit_with_message("Must use ambiguous torsion constraints with ligand shear moves!");
	}
}

LigandDockProtocol::LigandDockProtocol(
	std::string const & protocol,
	bool const minimize_ligand,
	bool const minimize_backbone,
	bool const tether_ligand,
	bool const mutate_same_name3,
	core::Real const ligand_chi_stddev_deg,
	core::Real const protein_CA_stddev_Ang,
	core::Real const ligand_tether_stddev_Ang,
	core::Size const ligand_shear_moves
):
	LigandBaseProtocol(),
	protocol_(protocol),
	minimize_ligand_(minimize_ligand),
	minimize_backbone_(minimize_backbone),
	tether_ligand_(tether_ligand),
	ligand_protonation_(mutate_same_name3),
	ligand_chi_stddev_deg_(ligand_chi_stddev_deg),
	protein_CA_stddev_Ang_(protein_CA_stddev_Ang),
	ligand_tether_stddev_Ang_(ligand_tether_stddev_Ang),
	ligand_shear_moves_(ligand_shear_moves),
	minimize_all_rsds_(false),
	repack_all_rsds_(false),
	rottrials_all_rsds_(false),
	minimize_water_(true),
	start_from_pts_(),
	ligand_torsion_restraints_()
{
	Mover::type( "LigandDockProtocol" );
	sc_interface_padding_= 0.0;//5.0; // 5A ends up repacking half the protein! (literally)
	bb_interface_cutoff_ = 7.0; //5.0;

	using basic::options::option;
	using namespace basic::options;

	if ( ligand_shear_moves_ && !option[ OptionKeys::docking::ligand::use_ambig_constraints ]() ) {
		utility_exit_with_message("Must use ambiguous torsion constraints with ligand shear moves!");
	}
}

void
LigandDockProtocol::add_start_from(core::Real x, core::Real y, core::Real z){
	start_from_pts_.push_back( core::Vector(x, y, z) );
}

LigandDockProtocol::LigandDockProtocol(LigandDockProtocol const & /*that*/):
	//utility::pointer::ReferenceCount(),
	protocols::moves::Mover(),
	LigandBaseProtocol()
{
	utility_exit_with_message("copy c-tor not allowed!");
}


LigandDockProtocol::~LigandDockProtocol() {}


/// @brief Creates a new hierarchy of Movers for each Pose passed in.
/// @details Some Movers (e.g. repack) require knowledge of the Pose to create,
///  and are only valid for that Pose and other conformations of it
///  (i.e. poses with the same number of residues, jumps, etc).
///  In general, we expect each Pose coming in here to be from a different PDB file.
///  The cost of creating these objects anew should be minimal compared
///  to the processing work they do.
void
LigandDockProtocol::apply( core::pose::Pose & pose )
{
	basic::prof_reset();
	using namespace protocols::moves;
	using utility::file::FileName;

	core::pack::dunbrack::load_unboundrot(pose); // adds scoring bonuses for the "unbound" rotamers, if any

	if ( protocol_ == "rescore" ) return;

	// Run most of search with soft-rep, but do final minimization and scoring with hard-rep.
	scorefxn_ = ( use_soft_rep_ ? soft_scorefxn_ : hard_scorefxn_ );

	// If we change the fold tree during the course of this run,
	// we should restore it afterwards for scoring and output!
	core::kinematics::FoldTree fold_tree_copy = pose.fold_tree();

	core::Size const jump_id= get_ligand_jump_id(pose);
	core::Size const lig_id = get_ligand_id(pose, jump_id);

	// Scoring function already set up by superclass
	// BUG:  currently need to score pose explicitly to get everything initialized properly
	(*scorefxn_)( pose );

	if ( protocol_ != "unbound" ) random_conformer(pose); // now only "pose" parameter is needed, actually

	move_ligand_to_desired_centroid(pose, jump_id, start_from_pts_);

	// Includes initial random perturbation (-dock_pert, -randomize2, etc)
	if ( protocol_ != "unbound" ) optimize_orientation3(pose, jump_id, lig_id);

	// Safer to add these after the initial rotamer juggling, especially in grid mode.
	// However, it means there should be no minimization done before this point!!
	// Have to do this very early, before MC or anyone else makes copies!
	// Only turn these on if needed? used to interfere with ligand packing (e.g. proton rotamers)
	ligand_torsion_restraints_.clear();
	if ( minimize_ligand_ ) {
		restrain_ligand_chis( pose);
	}

	// Make sure ligand doesn't start out in outer space.
	// We DO NOT start with a slide apart step -- for enclosed binding pockets,
	// it's seriously unlikely you'll ever find it again, especially if there are small bumps!
	if ( protocol_ != "unbound" ) {
		protocols::docking::FaDockingSlideIntoContact slideTogether(jump_id);
		slideTogether.apply( pose );
	}

	// Modifies pose (foldtree) and jump_id!
	if ( minimize_backbone_ ) {
		setup_bbmin_foldtree(pose, jump_id, bb_interface_cutoff_, protein_CA_stddev_Ang_);
	}
	// Put the move-map here so the interface matches up with the backbone constraints!
	core::kinematics::MoveMapOP movemap = make_movemap(pose, jump_id, sc_interface_padding_, minimize_all_rsds_, minimize_backbone_, minimize_ligand_, minimize_water_);

	// Only want to do this once the ligand is in its "final" starting place.
	core::scoring::constraints::ConstraintOP ligand_tether( NULL );
	if ( protocol_ != "unbound" ) {
		if ( tether_ligand_ ) ligand_tether = restrain_ligand_nbr_atom(pose, lig_id, ligand_tether_stddev_Ang_);
	}

	// Create a MonteCarlo object
	// Want to do this after perturb so we don't reset to pre-perturb state
	MonteCarloOP monteCarlo( new MonteCarlo(pose, *scorefxn_, 2.0 /* temperature, from RosettaLigand paper */) );

	if ( protocol_ == "meiler2006" ) classic_protocol(pose, jump_id, scorefxn_, monteCarlo, 50, 8); // Meiler and Baker 2006
	// pack - rottrials - rottrials - rottrials - pack
	else if ( protocol_ == "abbreviated" ) classic_protocol(pose, jump_id, scorefxn_, monteCarlo, 5, 4); // Davis ca. 2007
	// pack - RT - RT - pack - RT - RT (avoids ending on pack to avoid noise?)
	else if ( protocol_ == "abbrev2" ) classic_protocol(pose, jump_id, scorefxn_, monteCarlo, 6, 3); // Davis ca. April 2008
	else if ( protocol_ == "shear_min" ) shear_min_protocol(pose, jump_id, scorefxn_, monteCarlo, 20);
	else if ( protocol_ == "min_only" || protocol_ == "unbound" ) {} // no docking steps, just minimize (mostly for debugging/testing)
	else utility_exit_with_message("Unknown protocol '"+protocol_+"'");

	// Remove the ligand tether.  Could wait until after the final minimization,
	// but have to do it *some* time, or it will interfere with value of interface_delta.
	if ( tether_ligand_ ) pose.remove_constraint( ligand_tether );

	// keep the best structure we found, not the current one
	monteCarlo->show_scores();
	monteCarlo->recover_low(pose);

	// Run most of search with soft-rep, but do final minimization and scoring with hard-rep.
	scorefxn_ = hard_scorefxn_;


	//movemap->show(TR, pose.n_residue());
	// Do the final, "tight" minimization of all DOFs
	// Set up move map for minimizing.
	// Have to do this after initial perturb so we get the "right" interface defn.
	// Putting it here, we will get a slightly different interface than is used during classic_protocol() ...
	protocols::simple_moves::MinMoverOP dfpMinTightTol( new protocols::simple_moves::MinMover( movemap, scorefxn_, "lbfgs_armijo_nonmonotone_atol", 0.02, true /*use_nblist*/ ) );
	dfpMinTightTol->min_options()->nblist_auto_update(true);
	dfpMinTightTol->apply(pose);

	// Fast full-atom relax to eliminate backbone clashes
	// This moves the whole protein around, a lot
	//{
	// protocols::relax::FastRelax relax( scorefxn_ );
	// core::pack::task::PackerTaskOP relax_task;
	// relax_task = pack::task::TaskFactory::create_packer_task( pose );
	// //relax_task->initialize_from_command_line().restrict_to_repacking();
	// relax_task->restrict_to_repacking(); // -ex1, -ex2, etc make this take too long for whole protein!
	// relax_task->or_include_current( true );
	// protocols::simple_moves::PackRotamersMoverOP relax_full_repack = new protocols::simple_moves::PackRotamersMover( scorefxn_, relax_task );
	// relax.set_full_repack(relax_full_repack); // ligand torsions are still constrained during this...
	// relax.apply( pose );
	//}

	// If we change the fold tree during the course of this run,
	// we should restore it afterwards for scoring and output!
	pose.fold_tree( fold_tree_copy );

	basic::prof_show();
}

std::string
LigandDockProtocol::get_name() const {
	return "LigandDockProtocol";
}


void
LigandDockProtocol::classic_protocol(
	core::pose::Pose & pose,
	int jump_id,
	core::scoring::ScoreFunctionOP scorefxn,
	protocols::moves::MonteCarloOP monteCarlo,
	core::Size num_cycles,
	core::Size repack_every_Nth
) const
{
	using namespace protocols::moves;
	using core::pack::task::PackerTaskOP;
	runtime_assert(repack_every_Nth >= 2);

	// Set up move map for minimizing.
	// Have to do this after initial perturb so we get the "right" interface defn.
	//core::kinematics::MoveMapOP movemap = make_movemap(pose, jump_id, sc_interface_padding_, minimize_all_rsds_, minimize_backbone_, minimize_ligand_);
	core::kinematics::MoveMapOP movemap = make_movemap(pose, jump_id, sc_interface_padding_, minimize_all_rsds_, false /*minimize_backbone_*/, minimize_ligand_, minimize_water_);

	// Set up the packer task
	PackerTaskOP repack_task = make_packer_task(pose, jump_id, sc_interface_padding_, repack_all_rsds_, ligand_protonation_);
	PackerTaskOP rottrials_task = make_packer_task(pose, jump_id, sc_interface_padding_, rottrials_all_rsds_, ligand_protonation_);

	// Rigid body exploration
	MoverOP simple_rigbod( new rigid::RigidBodyPerturbMover( jump_id, numeric::conversions::degrees(0.05), 0.1) );

	for ( core::Size cycle = 1; cycle <= num_cycles; ++cycle ) {
		// RotamerTrialsMover actually asks for a non-const OP to scorefxn, sadly.
		// this problem did not manifest until I fixed the ScoreFunctionCOP definition in ScoreFunction.fwd.hh
		MoverOP pack_mover(
			(cycle % repack_every_Nth == 1) ?
			(Mover *) new protocols::simple_moves::PackRotamersMover(scorefxn, repack_task) :
			(Mover *) new protocols::simple_moves::RotamerTrialsMover(scorefxn, *rottrials_task)
		);
		// Wrap it in something to disable the torsion constraints before packing!
		pack_mover = MoverOP( new protocols::ligand_docking::UnconstrainedTorsionsMover( pack_mover, ligand_torsion_restraints_ ) );

		//MoverOP dockmcm_mover = make_dockmcm_mover(pose, jump_id, pack_mover, simple_rigbod, movemap, scorefxn, monteCarlo);
		//dockmcm_mover->apply(pose);
		{
			protocols::simple_moves::MinMoverOP min_mover( new protocols::simple_moves::MinMover( movemap, scorefxn, "lbfgs_armijo_nonmonotone_atol", 1.0, true /*use_nblist*/ ) );
			min_mover->min_options()->nblist_auto_update(true); // does this cost us lots of time in practice?

			core::Real const score1 = (*scorefxn)( pose );
			simple_rigbod->apply(pose);

			pack_mover->apply(pose);

			core::Real const score2 = (*scorefxn)( pose );
			if ( score2 - score1 < 15.0 ) {
				//std::cout << "YES, minimizing" << std::endl;
				min_mover->apply(pose);
			} else {
				//std::cout << "NO, not minimizing" << std::endl;
			}
			(*scorefxn)( pose ); // no effect at all

			monteCarlo->boltzmann( pose );

			if ( ligand_shear_moves_ ) shear_min_protocol(pose, jump_id, scorefxn, monteCarlo, ligand_shear_moves_);
		}

		// We always want the option (after the initial unbiased pack)
		// of sticking with our current nicely minimized conformation.
		repack_task->or_include_current(true);
		rottrials_task->or_include_current(true); // no effect -- rottrials always includes current.
	}
}


void
LigandDockProtocol::shear_min_protocol(
	core::pose::Pose & pose,
	int jump_id,
	core::scoring::ScoreFunctionOP scorefxn,
	protocols::moves::MonteCarloOP monteCarlo,
	core::Size num_cycles
) const
{
	using namespace protocols::moves;
	using core::pack::task::PackerTaskOP;

	core::Size const lig_id = get_ligand_id(pose, jump_id);
	core::conformation::Residue const & lig_rsd = pose.residue(lig_id);
	if ( lig_rsd.nchi() == 0 ) {
		TR << "Warning! Shear minimization protocol attempted on ligand without movable chis. (Lig id " << lig_id << ")" << std::endl;
		// Protocol is effectively a no-op for rigid ligands.
		return;
	}

	// Set up move map for minimizing.
	// Have to do this after initial perturb so we get the "right" interface defn.
	core::kinematics::MoveMapOP movemap = make_movemap(pose, jump_id, sc_interface_padding_, minimize_all_rsds_, /*minimize_backbone_=*/ false, minimize_ligand_, minimize_water_);
	//// Really simple movemap -- just ligand DOFs
	//// This works very poorly!
	//core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap();
	//movemap->set_jump(jump_id, true);
	//movemap->set_chi(lig_id, true);

	//monteCarlo->reset_counters();
	for ( core::Size cycle = 1; cycle <= num_cycles; ++cycle ) {
		//TR << "shear_min_protocol(), cycle " << cycle << std::endl;
		protocols::simple_moves::MinMoverOP min_mover( new protocols::simple_moves::MinMover( movemap, scorefxn, "lbfgs_armijo_nonmonotone_atol", 1.0, true /*use_nblist*/ ) );
		min_mover->min_options()->nblist_auto_update(true); // does this cost us lots of time in practice?
		//core::Real const score1 = (*scorefxn)( pose );

		// First, dumb version:  compensating changes in two dihedrals, but without checking that they're 1 bond apart!
		// For one test case (1GWX), this is significantly better than dumb version 2, below.
		core::Real const max_angle = 90.0; // 10 is too small, 180 is too large
		core::Real const angle_delta = max_angle * 2. * (numeric::random::rg().uniform() - 0.5);
		core::Size const chi1 = numeric::random::rg().random_range(1, lig_rsd.nchi());
		core::Size const chi2 = numeric::random::rg().random_range(1, lig_rsd.nchi());
		pose.set_chi(chi1, lig_id, lig_rsd.chi(chi1) + angle_delta);
		pose.set_chi(chi2, lig_id, lig_rsd.chi(chi2) - angle_delta);

		//// Second dumb version:  uncorrelated changes in 1-3 angles
		//core::Real const max_angle = 90.0;
		//core::Size const num_pert( numeric::random::rg().random_range(1, 3) );
		//core::conformation::Residue const & lig_rsd = pose.residue(lig_id);
		//for(core::Size i = 1; i <= num_pert; ++i) {
		// core::Real const angle_delta = max_angle * 2. * (numeric::random::rg().uniform() - 0.5);
		// core::Size const chi1 = numeric::random::rg().random_range(1, lig_rsd.nchi());
		// pose.set_chi(chi1, lig_id, lig_rsd.chi(chi1) + angle_delta);
		//}

		//core::Real const score2 = (*scorefxn)( pose );
		//if(score2 - score1 < 15.0) {
		//std::cout << "YES, minimizing" << std::endl;
		min_mover->apply(pose);
		//} else {
		// //std::cout << "NO, not minimizing" << std::endl;
		//}
		(*scorefxn)( pose ); // no effect at all

		monteCarlo->boltzmann( pose );
	}
	//monteCarlo->show_state();
}


/// @brief Replace current ligand(s) with conformers picked randomly from the library.
void
LigandDockProtocol::random_conformer(
	core::pose::Pose & pose
) const
{
	using namespace basic::options;
	using namespace protocols::moves;
	using core::conformation::ResidueOP;

	if ( option[ OptionKeys::docking::ligand::random_conformer ]() ) {
		for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
			if ( pose.residue(i).is_polymer() ) continue;
			//(*scorefxn_)( pose ); scorefxn_->accumulate_residue_total_energies( pose ); std::cout << "Constraints before: " << pose.energies().total_energies()[ core::scoring::dihedral_constraint ] << std::endl;
			RandomConformerMoverOP rcm( new RandomConformerMover(i) );
			UnconstrainedTorsionsMoverOP utm( new UnconstrainedTorsionsMover( rcm, ligand_torsion_restraints_ ) );
			utm->apply(pose);
			//(*scorefxn_)( pose ); scorefxn_->accumulate_residue_total_energies( pose ); std::cout << "Constraints after: " << pose.energies().total_energies()[ core::scoring::dihedral_constraint ] << std::endl;
		}
	}

	// Rotamer trials superimposes around the "nbr atom", which can shift the
	// computed centroid by significant distances.
	// Thus, after substituting, we re-center the new conformer.
	// (Re-centering now occurs in main protocol.)
}


/// @brief Repeatedly randomize ligand orientation in a given location
/// to find where *some* ligand rotamer has minimal clashes with the backbone.
void
LigandDockProtocol::optimize_orientation3(
	core::pose::Pose & pose,
	int jump_id,
	core::Size lig_id
)
{
	//clock_t start_time = clock();
	using namespace protocols::moves;
	using namespace core::scoring;
	using namespace core::pack::task;
	using namespace basic::options;
	using utility::vector1;

	// Takes about 3.5 sec to make grid for 250 residue protein,
	// compared to 15 sec for 1000 cycles with large ligand.
	// Given overall time cost of docking protocol and difficulty of telling when
	// we've moved to a new starting (input) structure, it's not worth caching the grid.

	TR << "Making ligand grid ..." << std::endl;
	//clock_t start_time = clock();
	utility::pointer::shared_ptr<core::grid::CartGrid<int> >  grid = make_atr_rep_grid(pose, pose.residue(lig_id).nbr_atom_xyz());
	//clock_t end_time = clock();
	//TR << "Elapsed time: " << double(end_time - start_time) / CLOCKS_PER_SEC << " seconds" << std::endl;

	int atr(0), rep(0);
	grid_score_atr_rep(*grid, pose.residue(lig_id), atr, rep);
	TR << "Input score atr = " << atr << " ; rep = " << rep << " ; ha = " << pose.residue(lig_id).nheavyatoms() << std::endl;

	// Kinemage output for debugging -- much less useful than ED map format
	if ( option[ OptionKeys::docking::ligand::grid::grid_kin ].user() ) {
		std::ofstream of( option[ OptionKeys::docking::ligand::grid::grid_kin ]().name().c_str() );
		of << "@kinemage\n";
		of << "@group {grid}\n";
		of << "@dotlist {atr} color= greentint\n";
		grid_to_kin<int>(of, *grid, -1, -1, 2);
		of << "@dotlist {rep} color= pinktint off\n";
		grid_to_kin<int>(of, *grid, 1, 1, 2);
		of.close();
	}

	// OMap output for debugging
	if ( option[ OptionKeys::docking::ligand::grid::grid_map ].user() ) {
		grid->write_to_BRIX( option[ OptionKeys::docking::ligand::grid::grid_map ]().name() );
	}

	MoverOP initialPerturb( new protocols::docking::DockingInitialPerturbation(jump_id, false /*don't slide*/) );
	// Make sure the initial perturb lands our nbr_atom in an empty space.
	{
		core::pose::Pose const orig_pose( pose );
		// With no initial perturbation, this can get trapped in an endless loop otherwise!
		//while(true) {
		for ( Size cnt = 0; cnt < 50; ++cnt ) {
			initialPerturb->apply(pose);
			core::Vector c = pose.residue(lig_id).nbr_atom_xyz();
			// did our nbr_atom land in an empty space on the grid?
			// Don't want to insist the nbr_atom is in the attractive region, because it may not be!
			if ( grid->is_in_grid( c.x(), c.y(), c.z() ) && grid->getValue( c.x(), c.y(), c.z() ) <= 0 ) {
				TR << "Accepting ligand position with nbr_atom at " << c << std::endl;
				break;
			}
			TR << "Rejecting ligand position with nbr_atom at " << c << std::endl;
			pose = orig_pose; // reset and try again
		}
	}

	if ( ! option[ OptionKeys::docking::ligand::improve_orientation ].active() ) {
		return;
	}

	int const num_cycles = std::max(0, option[ OptionKeys::docking::ligand::improve_orientation ]());
	runtime_assert( num_cycles >= 0 );
	TR << "Doing " << num_cycles << " trials of improvement for ligand orientation..." << std::endl;

	// Want to keep the same center of rotation even after rotamer trials,
	// so the ligand can't "walk" away from the protein!
	// RigidBodyRandomizeMover can't do that, so we didn't use it here.
	//core::Vector dummy_up, rot_center;
	//protocols::geometry::centroids_by_jump(pose, jump_id, dummy_up, rot_center);
	// Better to rotate around the nbr_atom, because that's where rotamer trials are centered.
	// Yes, I know this code will ONLY work for single residue ligands.  Tough.
	core::Vector const rot_center = pose.residue(lig_id).nbr_atom_xyz();

	// Can't refer to starting position (e.g. in case we're not using -randomize2, makes it too easy)
	// As long as we're doing at least one cycle, these will be overwritten by the first random.
	// These are kept here so we get the expected behavior if the user asks for 0 cycles.
	// Scale back the "perfect" values a little to avoid over-convergence.
	int const perfect_rep = 0, perfect_atr = -(int)(0.85 * pose.residue(lig_id).nheavyatoms());
	int best_rep(0), best_atr(0); // dummy initial values for compiler
	grid_score_atr_rep(*grid, pose.residue(lig_id), best_atr, best_rep);
	TR << "Starting orientation energy atr = " << best_atr << " ; rep = " << best_rep << std::endl;
	core::kinematics::Jump best_jump = pose.jump( jump_id );

	// Many poses will fall into broad energy wells, and a few will fall into narrow ones.
	// By collecting "perfect" poses and enforcing a minimum level of diversity,
	// we should be able to pick from them randomly and enrich for rare poses.
	// These next two parameters are wild-ass heuristic guesses that seem OK for the Meiler x-dock set.
	core::Real const diverse_rms = 0.65 * std::sqrt( (double)pose.residue_type(lig_id).nheavyatoms() );
	// Setting this too low leads to little enrichment in rare poses,
	// but setting it too high wastes a lot of time in rms calculation:
	core::Size const max_diversity = 5 * (pose.residue_type(lig_id).nchi()+1);
	vector1< core::kinematics::Jump > perfect_jumps;
	vector1< core::conformation::ResidueOP > perfect_rsds;
	core::Size perfect_count(0), diverse_count(0);

	// If we have constraints, we may want to use them as a filter.
	// Since we're not storing poses, though, we have to score as we go.
	bool const score_csts = option[ OptionKeys::enzdes::cstfile ].user();
	core::scoring::ScoreFunction cst_scorefxn;
	cst_scorefxn.set_weight(core::scoring::coordinate_constraint, 1.0 );
	cst_scorefxn.set_weight(core::scoring::atom_pair_constraint, 1.0 );
	cst_scorefxn.set_weight(core::scoring::angle_constraint, 1.0 );
	cst_scorefxn.set_weight(core::scoring::dihedral_constraint, 1.0 );
	vector1< core::Real > perfect_cstscores;

	for ( int i = 0; i < num_cycles; ++i ) {
		//randomize2->apply(pose);
		// Randomize orientation:  copied from RigidBodyRandomizeMover.apply()
		core::kinematics::Jump flexible_jump = pose.jump( jump_id );
		core::kinematics::Stub upstream_stub = pose.conformation().upstream_jump_stub( jump_id );
		core::kinematics::Stub downstream_stub = pose.conformation().downstream_jump_stub( jump_id );
		// comments for set_rb_center() explain which stub to use when!
		flexible_jump.set_rb_center( 1 /* n2c -- isn't defd in any .hh file :( */, downstream_stub, rot_center );
		flexible_jump.rotation_by_matrix( upstream_stub, rot_center, protocols::geometry::random_reorientation_matrix() );
		pose.set_jump( jump_id, flexible_jump );

		grid_rotamer_trials_atr_rep(*grid, pose, lig_id);
		int curr_rep(0), curr_atr(0); // dummy initial values for compiler
		grid_score_atr_rep(*grid, pose.residue(lig_id), curr_atr, curr_rep);
		// Always accept first try as best so far ... sort of a do-while loop.
		if ( i == 0 || curr_rep < best_rep || (curr_rep == best_rep && curr_atr < best_atr) ) {
			best_rep = curr_rep;
			best_atr = curr_atr;
			best_jump = pose.jump( jump_id );
			//if(best_rep <= perfect_rep && best_atr <= perfect_atr) {
			// TR << "Aborting after " << i+1 << " cycles because score is perfect!" << std::endl;
			// break; // can't do any better than this!
			//}
		}
		// Accumulate poses with a certain minimum level of diversity
		if ( curr_rep <= perfect_rep && curr_atr <= perfect_atr ) {
			perfect_count += 1;
			core::Real min_rms = 9999;
			for ( Size j = 1, j_end = perfect_rsds.size(); j <= j_end; ++j ) {
				core::Real rms = automorphic_rmsd(pose.residue(lig_id), *perfect_rsds[j], false /*don't superimpose*/);
				if ( rms < min_rms ) min_rms = rms;
			}
			if ( min_rms >= diverse_rms ) {
				diverse_count += 1;
				perfect_rsds.push_back( pose.residue(lig_id).clone() );
				perfect_jumps.push_back( pose.jump( jump_id ) );
				perfect_cstscores.push_back( cst_scorefxn(pose) );
				if ( diverse_count >= max_diversity ) {
					TR << "Aborting after " << i+1 << " cycles with " << diverse_count << " diverse 'perfect' poses" << std::endl;
					break; // can't do any better than this!
				}
			}
		}
		TR.Debug << "  Randomize energy atr = " << curr_atr << " ; rep = " << curr_rep << std::endl;
		TR.Debug << "    NBR_ATOM at " << pose.residue(lig_id).nbr_atom_xyz() << std::endl;
	}
	TR << "Best random orientation energy atr = " << best_atr << " ; rep = " << best_rep << std::endl;
	TR << "Found " << perfect_count << " 'perfect' poses, kept " << diverse_count << " with rms >= " << diverse_rms << std::endl;

	if ( !perfect_rsds.empty() ) {
		// Found multiple diverse and high-quality poses.  Choose one at random.
		core::Size which_perfect = (core::Size) numeric::random::rg().random_range(1, perfect_rsds.size());
		// If we have constraints, take the perfect pose with the best constraint score
		if ( score_csts ) {
			core::Real min_cstscore = 1e99;
			for ( core::Size j = 1, j_end = perfect_rsds.size(); j <= j_end; ++j ) {
				if ( perfect_cstscores[j] < min_cstscore ) {
					min_cstscore = perfect_cstscores[j];
					which_perfect = j;
					//std::cout << "*** choosing perfect pose " << j << " with cstscore = " << min_cstscore << std::endl;
				}
			}
		}
		pose.set_jump( jump_id, perfect_jumps[which_perfect] );
		// Just copy over XYZ coordinates
		for ( core::Size j = 1, j_end = pose.residue_type(lig_id).natoms(); j <= j_end; ++j ) { // no refolds required
			core::id::AtomID atom_id(j, lig_id); // atom, residue
			pose.set_xyz(atom_id, perfect_rsds[which_perfect]->xyz(j));
		}
	} else {
		// Found only one, half-decent pose
		// recover best jump
		pose.set_jump( jump_id, best_jump );
		// recover best ligand pose
		// Not *guaranteed* to be the same one if multiple confs are iso-energetic, actually.
		// Not that that matters for our purposes.
		grid_rotamer_trials_atr_rep(*grid, pose, lig_id);
	}

	//clock_t end_time = clock();
	//TR << "Elapsed time: " << double(end_time - start_time) / CLOCKS_PER_SEC << " seconds" << std::endl;
}


/// @brief Build Mover for Monte Carlo Minimization protocol
protocols::moves::MoverOP
LigandDockProtocol::make_dockmcm_mover(
	core::pose::Pose const & /*pose*/,
	int /*jump_id*/,
	protocols::moves::MoverOP repack_mover,
	protocols::moves::MoverOP rigbod_mover,
	core::kinematics::MoveMapOP movemap, //< would be COP but MinMover wants OP
	core::scoring::ScoreFunctionOP scorefxn,
	protocols::moves::MonteCarloOP monteCarlo
) const
{
	using namespace protocols::moves;

	protocols::simple_moves::MinMoverOP min_mover( new protocols::simple_moves::MinMover( movemap, scorefxn, "lbfgs_armijo_nonmonotone_atol", 1.0, true /*use_nblist*/ ) );
	min_mover->min_options()->nblist_auto_update(true); // does this cost us lots of time in practice?

	MoverOP sequence_mover( new SequenceMover(
		rigbod_mover,
		repack_mover
		) );

	MoverOP jump_out_mover( new JumpOutMover(
		sequence_mover,
		min_mover,
		scorefxn,
		15.0 // energy units, taken from Rosetta++ docking_minimize.cc
		) );

	TrialMoverOP mctrial( new TrialMover( jump_out_mover, monteCarlo ) );

	//mctrial->set_keep_stats(false); // big time sink for MC
	mctrial->keep_stats_type( no_stats );
	return mctrial;
}


/// @brief Helper function to tame the PoseMetricCalculator madness.
core::Size count_buried_unsat_Hbonds(core::pose::Pose const & pose)
{
	using namespace protocols::toolbox::pose_metric_calculators;
	// Despite the name, it's counting H-bonders, not any old polars.
	BuriedUnsatisfiedPolarsCalculator calc("default", "default");
	basic::MetricValue< core::Size > val;
	calc.get("all_bur_unsat_polars", val, pose);
	return val.value();
}


/// Because inquiring minds want to know: what Hbonds are precluded by this docking?
void print_buried_unsat_Hbonds(core::pose::Pose const & bound, core::pose::Pose const & unbound)
{
	using namespace protocols::toolbox::pose_metric_calculators;
	using core::id::AtomID_Map;
	// Despite the name, it's counting H-bonders, not any old polars.
	BuriedUnsatisfiedPolarsCalculator calc_bound("default", "default"), calc_unbound("default", "default");
	basic::MetricValue< AtomID_Map< bool > > map_bound, map_unbound;
	calc_bound.get("atom_bur_unsat", map_bound, bound);
	calc_unbound.get("atom_bur_unsat", map_unbound, unbound);
	// This is ridiculously inefficient but I don't see a better way...
	for ( core::Size r = 1; r <= map_bound.value().n_residue(); ++r ) {
		for ( core::Size a = 1; a <= map_bound.value().n_atom(r); ++a ) {
			if ( map_bound.value()(r,a) && !map_unbound.value()(r,a) ) {
				TR << "Unsatisfied interface H-bond: " << bound.residue_type(r).name3() << " " << r << " " << bound.residue_type(r).atom_name(a) << std::endl;
			}
		}
	}
}


/// @brief Scores to be output that aren't normal scorefunction terms.
void
LigandDockProtocol::append_ligand_docking_scores(
	core::pose::Pose const & before,
	core::pose::Pose const & after,
	core::scoring::ScoreFunctionCOP scorefxn,
	std::map< std::string, core::Real > & scores, //< appended to for output
	protocols::toolbox::match_enzdes_util::EnzConstraintIOCOP constraint_io /*= NULL*/
) const
{
	using namespace core::scoring;
	Size const jump_id = get_ligand_jump_id(after);

	// Figure out energy across the interface.
	// A truer "binding energy" would allow the components to relax (repack)
	// once separated, but Chu's JMB paper found this doesn't really help,
	// at least for protein-protein docking.
	// Also, it's very slow, requiring 10+ independent repacks.
	core::pose::PoseOP after_unbound( new core::pose::Pose( after ) );
	// If constraints aren't removed, pulling the components apart gives big penalties.
	if ( constraint_io ) constraint_io->remove_constraints_from_pose(*after_unbound, /*keep_covalent=*/ false, /*fail_on_constraints_missing=*/ true);

	// A very hacky way of guessing whether the components are touching:
	// if pushed together by 1A, does fa_rep change at all?
	// (The docking rb_* score terms aren't implemented as of this writing.)
	core::Real const together_score = (*scorefxn)( *after_unbound );
	EnergyMap const together_energies = after_unbound->energies().total_energies();
	core::Real const initial_fa_rep = after_unbound->energies().total_energies()[ fa_rep ];
	protocols::rigid::RigidBodyTransMover trans_mover( *after_unbound, jump_id );
	trans_mover.trans_axis( trans_mover.trans_axis().negate() ); // now move together
	trans_mover.step_size(1);
	trans_mover.apply( *after_unbound );
	(*scorefxn)( *after_unbound );
	core::Real const push_together_fa_rep = after_unbound->energies().total_energies()[ fa_rep ];
	bool const are_touching = (std::abs(initial_fa_rep - push_together_fa_rep) > 1e-4);
	scores["ligand_is_touching"] = are_touching;

	// Now pull apart by 500 A to determine the reference E for calculating interface E.
	trans_mover.trans_axis( trans_mover.trans_axis().negate() ); // now move apart
	trans_mover.step_size(500); // make sure they're fully separated!
	trans_mover.apply( *after_unbound );
	core::Real const separated_score = (*scorefxn)( *after_unbound );
	EnergyMap const separated_energies = after_unbound->energies().total_energies();
	scores["interface_delta"] = together_score - separated_score;

	// Interface delta, broken down by component
	for ( int i = 1; i <= n_score_types; ++i ) {
		ScoreType ii = ScoreType(i);
		if ( !scorefxn->has_nonzero_weight(ii) ) continue;
		scores[ "if_"+name_from_score_type(ii) ] = ( scorefxn->get_weight(ii) * (together_energies[ii] - separated_energies[ii]) );
	}

	// Fun stuff from the pose metrics calculators:
	core::Real const sasa_radius = 1.4;
	// Explicit cast to Real, as you can get negative numbers when the ligand is making buried hydrogen bonds to the protein
	scores["if_buried_unsat_hbonds"] = core::Real(count_buried_unsat_Hbonds(after)) - core::Real(count_buried_unsat_Hbonds(*after_unbound));
	scores["if_buried_sasa"] = -( core::scoring::calc_total_sasa(after, sasa_radius) - core::scoring::calc_total_sasa(*after_unbound, sasa_radius) );
	print_buried_unsat_Hbonds(after, *after_unbound);

	// Another interesting metric -- how far does the ligand centroid move?
	// Large values indicate we're outside of the intended binding site.
	core::Vector upstream_dummy, downstream_before, downstream_after;
	protocols::geometry::centroids_by_jump(before, jump_id, upstream_dummy, downstream_before);
	protocols::geometry::centroids_by_jump(after,  jump_id, upstream_dummy, downstream_after);
	if ( start_from_pts_.empty() ) {
		// Compare to starting position or...
		core::Real const ligand_centroid_travel = downstream_before.distance( downstream_after );
		scores["ligand_centroid_travel"] = ligand_centroid_travel;
	} else {
		// Compare to *nearest*  -start_from  point
		core::Real min_travel = 1e99;
		for ( Size ii = 1; ii <= start_from_pts_.size(); ++ii ) {
			core::Real travel = start_from_pts_[ii].distance( downstream_after );
			min_travel = std::min( min_travel, travel );
		}
		scores["ligand_centroid_travel"] = min_travel;
	}

	// Calculate radius of gyration for downstream non-H atoms
	// Ligands tend to bind in outstretched conformations...
	core::Real lig_rg = 0;
	int lig_rg_natoms = 0;
	ObjexxFCL::FArray1D_bool is_upstream ( before.total_residue(), false );
	before.fold_tree().partition_by_jump( jump_id, is_upstream );
	for ( core::Size i = 1, i_end = before.total_residue(); i <= i_end; ++i ) {
		if ( is_upstream(i) ) continue; // only downstream residues
		core::conformation::Residue const & rsd = before.residue(i);
		for ( core::Size j = 1, j_end = rsd.nheavyatoms(); j <= j_end; ++j ) {
			lig_rg += downstream_before.distance_squared( rsd.xyz(j) );
			lig_rg_natoms += 1;
		}
	}
	lig_rg = std::sqrt( lig_rg / lig_rg_natoms );
	scores["ligand_radius_of_gyration"] = lig_rg;

	// These RMSD values don't account for symmetry (e.g. phenyl rings)
	//scores["ligand_rms_no_super"] = rmsd_no_super(before, after, is_ligand_heavyatom);
	//scores["ligand_rms_with_super"] = rmsd_with_super(before, after, is_ligand_heavyatom);

	core::Size const lig_id = get_ligand_id(before, jump_id);
	runtime_assert( !before.residue(lig_id).is_polymer() );
	scores["ligand_auto_rms_with_super"] = automorphic_rmsd(before.residue(lig_id), after.residue(lig_id), true /*superimpose*/);
	scores["ligand_auto_rms_no_super"] = automorphic_rmsd(before.residue(lig_id), after.residue(lig_id), false /*don't superimpose*/);

	// These might be interesting for cases where we get part of the ligand
	// in the right place, but not all of it.
	utility::vector1< core::Real > cuts = utility::tools::make_vector1(0.5, 1.0, 2.0);
	utility::vector1< core::Real > fracs;
	protocols::ligand_docking::frac_atoms_within(before.residue(lig_id), after.residue(lig_id), cuts, fracs);
	scores["frac_atoms_within_0.5"] = fracs[1];
	scores["frac_atoms_within_1.0"] = fracs[2];
	scores["frac_atoms_within_2.0"] = fracs[3];

	// This might be useful in understanding entropic effects:
	scores["ligand_num_chi"] = before.residue(lig_id).nchi();
}

void
LigandDockProtocol::restrain_ligand_chis(
	core::pose::Pose & pose
){
	using basic::options::option;
	using namespace basic::options;

	if ( option[ OptionKeys::docking::ligand::use_ambig_constraints ] ) {
		constrain_ligand_torsions(pose, ligand_chi_stddev_deg_);
	} else {
		for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
			if ( pose.residue(i).is_polymer() ) continue;
			ligand_torsion_restraints_.push_back( protocols::ligand_docking::ResidueTorsionRestraintsOP( new protocols::ligand_docking::ResidueTorsionRestraints(pose, i, ligand_chi_stddev_deg_) ) );
		}
	}
}


} // namespace ligand_docking
} // namespace protocols
