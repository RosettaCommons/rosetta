// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    dock_glycans.cc
/// @brief   This application performs a bound-bound dock between an oligosaccharide and a carbohydrate-binding protein.
/// @note    The ligand must be the final chain (or chains, if it has branches).
/// @example TODO: fill in later
/// @author  Labonte <JWLabonte@jhu.edu>
/// @remarks This is an attempt to port dock_glycans.py to C++.


// Project headers
#include <devel/init.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/RingConformationMover.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/docking/DockingInitialPerturbation.hh>
#include <protocols/docking/util.hh>
#include <protocols/docking/metrics.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.hh>

#include <core/types.hh>
#include <core/id/AtomID.hh>
#include <core/chemical/rings/RingConformerSet.hh>
#include <core/kinematics/Edge.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/rms_util.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/import_pose/import_pose.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>

// Numeric headers
#include <numeric/random/random.hh>

// Basic headers
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/rings.OptionKeys.gen.hh>

// C++ headers
#include <string>
#include <iostream>


using namespace std;
using namespace utility;
using namespace core;
using namespace protocols;


// Class Definitions //////////////////////////////////////////////////////////

/// @brief  The protocol ported from dock_glycans.py.
class DockGlycansProtocol : public moves::Mover {
public:  // Standard methods
	/// @brief  Default constructor.
	DockGlycansProtocol() : moves::Mover()
	{
		init();
	}

	/// @brief  Copy constructor.
	DockGlycansProtocol( DockGlycansProtocol const & object_to_copy ) : Mover( object_to_copy )
	{
		copy_data( *this, object_to_copy );
	}

	// Assignment operator
	DockGlycansProtocol &
	operator=( DockGlycansProtocol const & object_to_copy )
	{
		// Abort self-assignment.
		if ( this == &object_to_copy ) {
			return *this;
		}

		moves::Mover::operator=( object_to_copy );
		copy_data( *this, object_to_copy );
		return *this;
	}

	// Destructor
	virtual ~DockGlycansProtocol() {}


public:  // Standard Rosetta methods
	// General methods
	/// @brief  Register options with the option system.
	static void
	register_options()
	{
		using namespace basic::options;
		using namespace protocols;

		option.add_relevant( OptionKeys::rings::idealize_rings );
		option.add_relevant( OptionKeys::rings::lock_rings );

		// Call register_options() on all other Movers used by this class.
		rigid::RigidBodyRandomizeMover::register_options();
		rigid::RigidBodyRandomizeMover::register_options();
		docking::FaDockingSlideIntoContact::register_options();
		rigid::RigidBodyPerturbMover::register_options();
		simple_moves::MinMover::register_options();
		simple_moves::RingConformationMover::register_options();
		simple_moves::SmallMover::register_options();
		simple_moves::ShearMover::register_options();
		simple_moves::MinMover::register_options();
	}

	/// @brief  Generate string representation of DockGlycansProtocol for debugging purposes.
	virtual
	void
	show( std::ostream & output=std::cout ) const
	{
		moves::Mover::show( output );  // name, type, tag
	}


	// Mover methods
	/// @brief  Return the name of the Mover.
	virtual
	std::string
	get_name() const
	{
		return type();
	}

	virtual
	moves::MoverOP
	clone() const
	{
		return protocols::moves::MoverOP( new DockGlycansProtocol( *this ) );
	}

	virtual
	moves::MoverOP
	fresh_instance() const
	{
		return protocols::moves::MoverOP( new DockGlycansProtocol );
	}

	/// @brief  Apply the corresponding protocol to <pose>.
	virtual
	void
	apply( core::pose::Pose & pose )
	{
		using namespace id;
		using namespace scoring;
		using namespace conformation;
		using namespace docking;
		using namespace rigid;

		show( cout );

		// Prepare the FoldTree.
		determine_docking_partners( pose );

		string const partners( upstream_chains_ + "_" + downstream_chains_ );
		vector1< int > movable_jumps( 1, JUMP_NUM );
		setup_foldtree( pose, partners, movable_jumps );

		// Prepare the ligand constraint, which holds it locally in the active site.
		kinematics::Edge const & jump_edge( pose.fold_tree().jump_edge( JUMP_NUM ) );
		core::uint const start_resnum( jump_edge.start() );
		core::uint const stop_resnum( jump_edge.stop() );
		Residue const & start_res( pose.residue( start_resnum ) );
		Residue const & stop_res( pose.residue( stop_resnum ) );
		core::uint const start_atomnum( start_res.nbr_atom() );
		core::uint const stop_atomnum( stop_res.nbr_atom() );
		AtomID const start_atom( start_atomnum, start_resnum );
		AtomID const stop_atom( stop_atomnum, stop_resnum );
		Distance const distance( ( stop_res.xyz( stop_atomnum ) - start_res.xyz( start_atomnum ) ).norm() );
		Distance const deviation( 3.00 );

		func::HarmonicFuncOP function( new func::HarmonicFunc( distance, deviation ) );
		constraints::ConstraintOP ligand_constraint(
			new constraints::AtomPairConstraint( start_atom, stop_atom, function ) );
		pose.add_constraint( ligand_constraint );

		// Print some information about the starting pose.
		cout << endl << "Starting pose:" << endl;
		cout << ' ' << pose.fold_tree();
		cout << " Ligand [chain(s) " << downstream_chains_ << "] center is ";
		cout << distance;
		cout << " angstroms from protein center [chain(s) " << upstream_chains_ << "]." << endl;

		prepare_scoring_function();

		cout << " Starting Score:" << endl;
		sf_->show( pose );
		cout << " Interface score:                          ";
		cout << calc_interaction_energy( pose, sf_, movable_jumps ) << endl;

		// Prepare MoveMaps and Movers that require Pose information.
		Size const n_residues( pose.total_residue() );
		jump_mm_->set_jump( JUMP_NUM, true );
		ring_mm_->set_nu_true_range( first_ligand_residue_, n_residues );
		torsion_mm_->set_jump( JUMP_NUM, true );
		torsion_mm_->set_bb_true_range( first_ligand_residue_ + 1, n_residues );
		torsion_mm_->set_chi_true_range( first_ligand_residue_, n_residues );
		torsion_mm_->set_nu( false );  // TEMP... until rings are treated properly by the MinMover

		randomizerA_ = RigidBodyRandomizeMoverOP(
			new RigidBodyRandomizeMover( pose, 1, partner_downstream, 360, 360, false ) );
		randomizerB_ = RigidBodyRandomizeMoverOP(
			new RigidBodyRandomizeMover( pose, 1, partner_upstream, 360, 360, false ) );


		cout << "Randomizing ligand conformation..." << endl;
		for ( core::uint residue( first_ligand_residue_ ); residue <= n_residues; ++residue ) {
			pose.set_phi( residue, numeric::random::rg().uniform() * 360 );
			pose.set_psi( residue, numeric::random::rg().uniform() * 360 );
			//pose.set_ring_conformation( residue, 1,
			//pose.residue( residue ).type().ring_conformer_set( 1 )->get_random_conformer() );
		}


		// Begin docking protocol.
		cout << endl << "Docking..." << endl;
		if ( idealize_rings_ ) {
			cout << " Idealizing rings..." << endl;
			for ( core::uint residue( first_ligand_residue_ ); residue <= n_residues; ++residue ) {
				pose.set_ring_conformation( residue, 1,
					pose.residue( residue ).type().ring_conformer_set( 1 )->get_lowest_energy_conformer() );
			}
		}

		cout << " Randomizing positions..." << endl;
		randomizerA_->apply( pose );
		slider_->apply( pose );

		cout << " Ligand center is ";
		cout << ( stop_res.xyz( stop_atomnum ) - start_res.xyz( start_atomnum ) ).norm();
		cout << " angstroms from protein center." << endl;

		cout << " Refining..." << endl;

		// Set beginning values for weights to ramp.
		sf_->set_weight( fa_atr, target_atr_ * STARTING_RAMP_DOWN_FACTOR );
		sf_->set_weight( fa_rep, target_atr_ * STARTING_RAMP_UP_FACTOR );

		mc_ = moves::MonteCarloOP( new moves::MonteCarlo( pose, *sf_, kt_ ) );

		for ( core::uint cycle( 1 ); cycle <= n_cycles_; ++cycle ) {
			if ( cycle % ( n_cycles_ / 10 ) == 0 ) {  // Ramp every ~10% of n_cycles.
				Real fraction = Real( cycle ) / n_cycles_;
				ramp_score_weight( fa_atr, target_atr_, fraction );
				ramp_score_weight( fa_rep, target_rep_, fraction );
				mc_->reset( pose );
			}

			// Rigid-body moves
			perturber_->apply( pose );
			slider_->apply( pose );
			jump_minimizer_->apply( pose );

			// Ligand torsion moves
			ring_mover_->apply( pose );
			small_mover_->apply( pose );
			shear_mover_->apply( pose );
			slider_->apply( pose );
			torsion_minimizer_->apply( pose );

			// Metropolis criterion.
			mc_->boltzmann( pose );

			if ( cycle % 5 == 0 ) {
				cout << "  Cycle " << cycle << "\tCurrent Score:" << ( *sf_ )( pose ) << endl;
			}
		}

		cout << " Final score for decoy: " << ( *sf_ )( pose ) << endl;

		if ( ref_pose_ ) {
			jd2::JobOP job( jd2::JobDistributor::get_instance()->current_job() );
			record_pose_metrics( pose, *job, movable_jumps );
		}
	}


	// Accessors/Mutators
	void
	score_function( core::scoring::ScoreFunctionOP input_sf )
	{
		sf_ = input_sf;
	}

	void
	set_ref_pose_from_filename( std::string const & filename )
	{
		ref_pose_ = import_pose::pose_from_file( filename, core::import_pose::PDB_file );
	}


private:  // Private methods
	// Set command-line options.  (Called by init())
	void
	set_commandline_options()
	{
		using namespace basic::options;

		idealize_rings_ = option[ OptionKeys::rings::idealize_rings ].active();
		lock_rings_ = option[ OptionKeys::rings::lock_rings ].active();
	}


	// Initialize data members from arguments.
	void
	init()
	{
		using namespace simple_moves;

		type( "DockGlycansProtocol" );

		Hbond_mult_ = 1.0;
		elec_mult_ = 1.0;
		sol_mult_ = 1.0;
		atr_mult_ = 1.0;
		rep_mult_ = 1.0;

		rot_ = 2.0;
		trans_ = 0.5;

		idealize_rings_ = true;
		lock_rings_ = true;

		// Instantiate the Movers.
		// Note: The randomizer objects require a Pose and so cannot be initialized here.

		slider_ = docking::FaDockingSlideIntoContactOP ( new docking::FaDockingSlideIntoContact( JUMP_NUM ) );
		perturber_ = rigid::RigidBodyPerturbMoverOP ( new rigid::RigidBodyPerturbMover( JUMP_NUM, rot_, trans_ ) );

		jump_mm_ = kinematics::MoveMapOP( new kinematics::MoveMap );
		ring_mm_ = kinematics::MoveMapOP( new kinematics::MoveMap );
		torsion_mm_ = kinematics::MoveMapOP( new kinematics::MoveMap );

		jump_minimizer_ = MinMoverOP( new MinMover( jump_mm_, sf_, "lbfgs_armijo_nonmonotone", 0.01, true ) );

		ring_mover_ = RingConformationMoverOP( new RingConformationMover( ring_mm_ ) );
		small_mover_ = SmallMoverOP( new SmallMover( torsion_mm_, kt_, 3 ) );
		shear_mover_ = ShearMoverOP( new ShearMover( torsion_mm_, kt_, 3 ) );

		torsion_minimizer_ = MinMoverOP( new MinMover( torsion_mm_, sf_, "lbfgs_armijo_nonmonotone", 0.01, true ) );

		n_cycles_ = 100;

		kt_ = 0.8;

		set_commandline_options();
	}

	// Copy all data members from <object_to_copy_from> to <object_to_copy_to>.
	void
	copy_data( DockGlycansProtocol & object_to_copy_to, DockGlycansProtocol const & object_to_copy_from )
	{
		object_to_copy_to.upstream_chains_ = object_to_copy_from.upstream_chains_;
		object_to_copy_to.downstream_chains_ = object_to_copy_from.downstream_chains_;
		object_to_copy_to.first_ligand_residue_ = object_to_copy_from.first_ligand_residue_;

		object_to_copy_to.sf_ = object_to_copy_from.sf_;
		object_to_copy_to.ref_pose_ = object_to_copy_from.ref_pose_;

		object_to_copy_to.Hbond_mult_ = object_to_copy_from.Hbond_mult_;
		object_to_copy_to.elec_mult_ = object_to_copy_from.elec_mult_;
		object_to_copy_to.sol_mult_ = object_to_copy_from.sol_mult_;
		object_to_copy_to.atr_mult_ = object_to_copy_from.atr_mult_;
		object_to_copy_to.rep_mult_ = object_to_copy_from.rep_mult_;

		object_to_copy_to.target_atr_ = object_to_copy_from.target_atr_;
		object_to_copy_to.target_rep_ = object_to_copy_from.target_rep_;

		object_to_copy_to.rot_ = object_to_copy_from.rot_;
		object_to_copy_to.trans_ = object_to_copy_from.trans_;

		object_to_copy_to.idealize_rings_ = object_to_copy_from.idealize_rings_;
		object_to_copy_to.lock_rings_ = object_to_copy_from.lock_rings_;

		// MoveMaps
		object_to_copy_to.jump_mm_ = object_to_copy_from.jump_mm_;
		object_to_copy_to.ring_mm_ = object_to_copy_from.ring_mm_;
		object_to_copy_to.torsion_mm_ = object_to_copy_from.torsion_mm_;

		// Movers
		object_to_copy_to.randomizerA_ = object_to_copy_from.randomizerA_;
		object_to_copy_to.randomizerB_ = object_to_copy_from.randomizerB_;
		object_to_copy_to.slider_ = object_to_copy_from.slider_;
		object_to_copy_to.perturber_ = object_to_copy_from.perturber_;
		object_to_copy_to.jump_minimizer_ = object_to_copy_from.jump_minimizer_;
		object_to_copy_to.ring_mover_ = object_to_copy_from.ring_mover_;
		object_to_copy_to.small_mover_ = object_to_copy_from.small_mover_;
		object_to_copy_to.shear_mover_ = object_to_copy_from.shear_mover_;
		object_to_copy_to.torsion_minimizer_ = object_to_copy_from.torsion_minimizer_;

		object_to_copy_to.n_cycles_ = object_to_copy_from.n_cycles_;

		object_to_copy_to.kt_ = object_to_copy_from.kt_;
		object_to_copy_to.mc_ = object_to_copy_from.mc_;
	}


	// Set the strings designating the upstream and downstream docking partners for this pose.
	// The strings determined are of the pdb file chain designations.
	// Currently, this assumes that the ligand is the final chain(s).
	// In addition, this subroutine sets the variable storing the first residue of the ligand.
	void
	determine_docking_partners( core::pose::Pose const & pose )
	{
		// Clear data.
		upstream_chains_ = "";
		downstream_chains_ = "";

		pose::PDBInfoCOP pdb_info( pose.pdb_info() );
		vector1< core::uint > const & chain_endings( pose.conformation().chain_endings() );
		Size const n_chain_endings( chain_endings.size() );
		//cout << chain_endings[ 1 ] << ", ..., " << chain_endings[ n_chain_endings ] << endl;
		core::uint chain_num_of_last_cut_point( n_chain_endings );
		//cout << chain_num_of_last_cut_point << endl;
		core::uint last_cut_point( chain_endings[ chain_num_of_last_cut_point ] );
		//cout << last_cut_point << endl;

		// Check the last chain to see if it is a branch.
		while ( pose.residue( last_cut_point + 1 ).is_branch_lower_terminus() ) {
			--chain_num_of_last_cut_point;
			last_cut_point = chain_endings[ chain_num_of_last_cut_point ];
			//cout << last_cut_point << endl;
		}

		core::uint cut_point;
		for ( core::uint chain_number( 1 ); chain_number <= n_chain_endings; ++chain_number ) {
			cut_point = chain_endings[ chain_number ];
			//cout << cut_point << endl;
			if ( cut_point < last_cut_point ) {
				upstream_chains_ += pdb_info->chain( cut_point );
			} else if ( cut_point == last_cut_point ) {
				upstream_chains_ += pdb_info->chain( cut_point );
				downstream_chains_ += pdb_info->chain( cut_point + 1 );
			} else /* cut_point > last_cut_point */ {
				downstream_chains_ += pdb_info->chain( cut_point + 1 );
			}
		}

		first_ligand_residue_ = last_cut_point + 1;
	}

	// Adjust the scoring function weights according to supplied multipliers and prepare for ramping.
	void
	prepare_scoring_function()
	{
		using namespace scoring;

		//sf_ = get_score_function();
		vector1< string > const patches( 1, "docking" );
		sf_ = ScoreFunctionFactory::create_score_function( "talaris2014", patches );

		sf_->set_weight( sugar_bb, 1.0 );
		sf_->set_weight( atom_pair_constraint, 1.0 );

		sf_->set_weight( hbond_sr_bb, sf_->get_weight( hbond_sr_bb ) * Hbond_mult_ );
		sf_->set_weight( hbond_lr_bb, sf_->get_weight( hbond_lr_bb ) * Hbond_mult_ );
		sf_->set_weight( hbond_bb_sc, sf_->get_weight( hbond_bb_sc ) * Hbond_mult_ );
		sf_->set_weight( hbond_sc, sf_->get_weight( hbond_sc ) * Hbond_mult_ );
		sf_->set_weight( fa_elec, sf_->get_weight( fa_elec ) * elec_mult_ );
		sf_->set_weight( fa_sol, sf_->get_weight( fa_sol ) * sol_mult_ );

		target_atr_ = sf_->get_weight( fa_atr ) * atr_mult_;
		target_rep_ = sf_->get_weight( fa_rep ) * rep_mult_;

		sf_->set_weight( fa_atr, target_atr_ );
		sf_->set_weight( fa_rep, target_rep_ );;
	}

	// Set the weight for the given method within this Mover's ScoreFunction using the appropriate ramping factor and
	// the fraction complete.
	void
	ramp_score_weight( core::scoring::ScoreType const method,
		core::Real const target,
		core::Real const fraction_completion )
	{
		Real factor;
		Real const current_weight( sf_->get_weight( method ) );

		if ( current_weight < target ) {
			Real const factor_range_size( 1 - STARTING_RAMP_UP_FACTOR );
			factor = STARTING_RAMP_UP_FACTOR + fraction_completion * factor_range_size;
		} else if ( current_weight > target ) {
			Real const factor_range_size( STARTING_RAMP_DOWN_FACTOR - 1 );
			factor = STARTING_RAMP_DOWN_FACTOR - fraction_completion * factor_range_size;
		} else {
			factor = 1;
		}

		Real const new_weight( target * factor );
		sf_->set_weight( method, new_weight );
	}

	// Record a collection of decoy metrics.
	void
	record_pose_metrics( core::pose::Pose const & pose,
		protocols::jd2::Job & job,
		utility::vector1< int > const jumps )
	{
		using namespace scoring;
		using namespace docking;

		Real const ligand_rmsd( non_peptide_heavy_atom_RMSD( pose, *ref_pose_ ) );
		Real const interaction_energy( calc_interaction_energy( pose, sf_, jumps ) );
		Real const fraction_native_contacts( calc_Fnat( pose, *ref_pose_, sf_, jumps ) );

		cout << " Metrics for this decoy:" << endl;
		cout << "  Ligand RMSD:                 " << ligand_rmsd << endl;
		cout << "  Interaction Energy:          " << interaction_energy << endl;
		cout << "  Fraction of native contacts: " << fraction_native_contacts << endl;

		job.add_string_real_pair( "ligand_rmsd", ligand_rmsd );
		job.add_string_real_pair( "interaction_energy", interaction_energy );
		job.add_string_real_pair( "Fnat", fraction_native_contacts );
	}

private:  // Private data
	std::string upstream_chains_;  // e.g., "AB"
	std::string downstream_chains_; // e.g., "XY"
	core::uint first_ligand_residue_;

	core::scoring::ScoreFunctionOP sf_;
	core::pose::PoseOP ref_pose_;

	core::Real Hbond_mult_;
	core::Real elec_mult_;
	core::Real sol_mult_;
	core::Real atr_mult_;
	core::Real rep_mult_;

	core::Real target_atr_;
	core::Real target_rep_;

	core::Angle rot_;
	core::Distance trans_;

	bool idealize_rings_;
	bool lock_rings_;


	// MoveMaps
	core::kinematics::MoveMapOP jump_mm_;
	core::kinematics::MoveMapOP ring_mm_;
	core::kinematics::MoveMapOP torsion_mm_;

	// Movers
	protocols::rigid::RigidBodyRandomizeMoverOP randomizerA_;
	protocols::rigid::RigidBodyRandomizeMoverOP randomizerB_;
	protocols::docking::FaDockingSlideIntoContactOP slider_;
	protocols::rigid::RigidBodyPerturbMoverOP perturber_;
	protocols::simple_moves::MinMoverOP jump_minimizer_;
	protocols::simple_moves::RingConformationMoverOP ring_mover_;
	protocols::simple_moves::SmallMoverOP small_mover_;
	protocols::simple_moves::ShearMoverOP shear_mover_;
	protocols::simple_moves::MinMoverOP torsion_minimizer_;

	core::Size n_cycles_;

	core::Real kt_;
	protocols::moves::MonteCarloOP mc_;


private:  // Constants
	static core::uint const JUMP_NUM;

	static core::Real const STARTING_RAMP_DOWN_FACTOR;
	static core::Real const STARTING_RAMP_UP_FACTOR;
};


// Constants & Type Definitions ///////////////////////////////////////////////
int const SUCCESS( 0 );
int const FAILURE( -1 );

core::uint const DockGlycansProtocol::JUMP_NUM( 1 );
core::Real const DockGlycansProtocol::STARTING_RAMP_DOWN_FACTOR( 3.25 );  // values used by Krishna
core::Real const DockGlycansProtocol::STARTING_RAMP_UP_FACTOR( 0.45455 );  // values used by Krishna

typedef utility::pointer::shared_ptr< DockGlycansProtocol > DockGlycansProtocolOP;


// Main ///////////////////////////////////////////////////////////////////////
int
main( int argc, char *argv[] )
{
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		// Initialize Rosetta.
		cout << "Initializing Rosetta..." << endl;
		devel::init( argc, argv );

		// Construct the protocol.
		DockGlycansProtocolOP protocol( new DockGlycansProtocol );

		// Set user options.
		if ( option[ in::file::native ].active() ) {
			protocol->set_ref_pose_from_filename( option[ in::file::native ] );
		}

		// Distribute the mover.
		protocols::jd2::JobDistributor::get_instance()->go( protocol );
	} catch ( utility::excn::EXCN_Base const & e ) {
		cerr << "Caught exception: " << e.msg() << endl;
		return FAILURE;
	}
	return SUCCESS;
}
