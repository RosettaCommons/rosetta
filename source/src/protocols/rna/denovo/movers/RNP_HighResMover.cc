// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file RNP_HighResMover.cc
/// @brief protocols that are specific to RNP_HighResMover
/// @details High-resolution sampling of an RNA/protein structure -- RNA frag insertions and
///          protein sidechain packing
/// @author Kalli Kappel


#include <protocols/rna/denovo/movers/RNP_HighResMover.hh>
#include <protocols/rna/denovo/movers/RNA_FragmentMover.hh>
#include <protocols/rna/denovo/movers/RNA_Minimizer.hh>
#include <core/import_pose/options/RNA_FragmentMonteCarloOptions.fwd.hh>
#include <protocols/rna/denovo/util.hh>
#include <protocols/rna/movers/RNA_LoopCloser.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/pose/Pose.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/rotamer_trials.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/conformation/Residue.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/Energies.hh>


#include <core/scoring/constraints/ConstraintSet.hh>

#include <basic/options/option.hh>


//Relaxer stuff

// ObjexxFCL Headers

#include <core/types.hh>
#include <basic/Tracer.hh>

#include <numeric/random/random.hh>

// External library headers


//C++ headers
#include <utility>
#include <vector>
#include <string>
#include <sstream>

#if defined(WIN32) || defined(__CYGWIN__)
#include <ctime>
#endif

// option key includes

#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <core/kinematics/Jump.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/format.hh>

//Auto using namespaces
namespace ObjexxFCL { namespace format { } } using namespace ObjexxFCL::format; // AUTO USING NS
//Auto using namespaces end


using namespace core;


static basic::Tracer TR( "protocols.rna.denovo.movers.RNP_HighResMover" );

namespace protocols {
namespace rna {
namespace denovo {
namespace movers {

//////////////////////////////////////////////////////////////////////////////////////////
RNP_HighResMover::RNP_HighResMover( RNA_FragmentMoverOP rna_fragment_mover, protocols::rna::movers::RNA_LoopCloserOP rna_loop_closer, core::import_pose::options::RNA_FragmentMonteCarloOptionsCOP options ):
	Mover(),
	rna_fragment_mover_(std::move( rna_fragment_mover )),
	rna_loop_closer_(std::move( rna_loop_closer )),
	options_(std::move( options )),
	fragment_cycles_( 200 ),
	frag_size_( 1 ),
	is_init_( false ),
	docking_( false )
{
	protein_sc_pack_cutoff_ = 20.0;
	Mover::type("RNP_HighResMover");
}

RNP_HighResMover::~RNP_HighResMover() = default;

// initialize, this is going to set up the rnp docking mover and the rnp docking ft
// this class can then store the mover and the fold tree
// hopefully this can all be removed as soon as we update the fold tree setup in rna_denovo
void
RNP_HighResMover::initialize( core::pose::Pose const & pose ) {

	// New method of setting up rigid body docking mover
	core::kinematics::MoveMap movemap;
	movemap.set_jump( false );

	core::Real rot_mag( 0.2 );
	core::Real trans_mag( 0.1 );

	bool move_first_rigid_body = false;
	if ( options_->move_first_rigid_body() || options_->dock_into_density() ) {
		move_first_rigid_body = true;
	}

	bool allow_single_rigid_body = false;
	if ( options_->dock_into_density() ) allow_single_rigid_body = true;

	// this updates the movemap
	bool rigid_body_moves = let_rigid_body_jumps_move( movemap, pose, move_first_rigid_body, allow_single_rigid_body /* for moving absolute coordinates*/ );

	if ( !rigid_body_moves ) {
		docking_ = false;
	} else {
		docking_ = true;
		rnp_docking_mover_ = protocols::rigid::RigidBodyPerturbMoverOP( new protocols::rigid::RigidBodyPerturbMover( pose, movemap, rot_mag, trans_mag, protocols::rigid::partner_upstream /*because virtual anchor should be root*/ ) );
	}

	is_init_ = true;
}

/// @details  Apply the RNP high-res mover
///
///////////////////////////////////////////////////////////////////////////////////////////
void RNP_HighResMover::apply( core::pose::Pose & pose )
{
	using namespace core::scoring;
	using namespace core::kinematics;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	//check that we've actually initialized stuff
	if ( !is_init_ ) utility_exit_with_message( "RNP_HighResMover has not been initialized before apply!" );

	time_t pdb_start_time = time(nullptr);

	ScoreFunctionOP scorefxn;

	if ( option[ score::weights ].user() ) {
		scorefxn = get_score_function();
	} else {
		// this score function definitely needs to be updated!!
		scorefxn = ScoreFunctionFactory::create_score_function( "rna/denovo/rna_hires_with_protein.wts" );
	}
	// turn on constraint weights
	if ( pose.constraint_set()->has_constraints() ) {
		if ( !scorefxn->has_nonzero_weight( atom_pair_constraint ) )  scorefxn->set_weight( atom_pair_constraint, 1.0 );
		if ( !scorefxn->has_nonzero_weight( base_pair_constraint ) )  scorefxn->set_weight( base_pair_constraint, 1.0 );
		if ( !scorefxn->has_nonzero_weight( coordinate_constraint ) ) scorefxn->set_weight( coordinate_constraint, 1.0 );
	}

	// slowly introduce fa_rep (we just added in protein sidechains, they might really be clashing with the RNA)

	// make sure that the score function has density term, if we have density
	if ( options_->model_with_density() ) {
		if ( scorefxn->get_weight( core::scoring::elec_dens_fast ) <= 0 ) {
			scorefxn->set_weight( core::scoring::elec_dens_fast, 10.0 );
		}
	}


	protocols::moves::MonteCarloOP monte_carlo_( new protocols::moves::MonteCarlo( pose, *scorefxn, 2.0 ) );

	// also want to do small RNP docking moves
	core::Real const final_fa_rep_weight = scorefxn->get_weight( fa_rep );

	for ( Size r = 1; r <= options_->rnp_high_res_cycles(); r++ ) { // default is 10 cycles

		TR << "RNP high-res mover ROUND " << r << " of " << options_->rnp_high_res_cycles() << std::endl;

		if ( options_->rnp_ramp_rep() ) {
			// scale the fa_rep weight based on the round
			core::Real const scaling_factor = static_cast<Real>( r ) / options_->rnp_high_res_cycles();
			scorefxn->set_weight( fa_rep, scaling_factor*final_fa_rep_weight );
		}

		// we should start with protein sidechain packing!
		if ( options_->rnp_pack_first() ) {
			protein_sidechain_packing( pose, scorefxn );
			monte_carlo_->boltzmann( pose, "initial_pack" );
		}

		// do RNA fragment insertions
		for ( Size j = 1; j<=fragment_cycles_; ++j ) {
			if ( docking_ && j % 10 == 0 ) {
				rnp_docking( pose );
				monte_carlo_->boltzmann( pose, "rnp_dock" );
			} else {
				rna_fragment_mover_->random_fragment_insertion( pose, frag_size_ );
				// should we do this after every move?
				rna_loop_closer_->apply( pose );
				monte_carlo_->boltzmann( pose, "rnp_frag" );
			}
		}

		// more small docking moves, hard code number of tries for now
		if ( docking_ ) {
			for ( Size j = 1; j<=10; ++j ) {
				rnp_docking( pose );
				monte_carlo_->boltzmann( pose, "rnp_dock" );
			}
		}

		// get the lowest score pose before packing protein sidechains
		pose = monte_carlo_->lowest_score_pose();

		// do protein sidechain repacking
		// (just residues that are near the RNA)
		protein_sidechain_packing( pose, scorefxn );

		monte_carlo_->boltzmann( pose, "rnp_high_res" );
		monte_carlo_->show_counters();
	}

	monte_carlo_->show_counters();

	pose = monte_carlo_->lowest_score_pose();

	time_t pdb_end_time = time(nullptr);

	scorefxn->show( TR, pose );

	TR << "RNP high-res mover finished in " << (long)(pdb_end_time - pdb_start_time) << " seconds." << std::endl;

}

///////////////////////////////////////////////////////////////////////////////////////////
std::string
RNP_HighResMover::get_name() const {
	return "RNP_HighResMover";
}
///////////////////////////////////////////////////////////////////////////////////////////
void
RNP_HighResMover::protein_sidechain_packing( core::pose::Pose & pose, core::scoring::ScoreFunctionCOP const & packer_scorefxn ) const {

	pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));

	// get the list of protein residues that are near the RNA
	utility::vector1< Size > residues_to_pack = get_residues_within_dist_of_RNA( pose, protein_sc_pack_cutoff_);

	for ( Size i=1; i<= pose.size(); ++i ) {
		// Check whether this is in the list of residues to pack
		if ( std::find(residues_to_pack.begin(), residues_to_pack.end(), i ) != residues_to_pack.end() ) {
			// if it's RNA
			if ( pose.residue( i ).is_RNA() ) {
				task->nonconst_residue_task(i).and_extrachi_cutoff( 0 );
				task->nonconst_residue_task(i).or_include_current( true );
			} else if ( pose.residue( i ).is_protein() ) {
				task->nonconst_residue_task(i).restrict_to_repacking();
				task->nonconst_residue_task(i).or_include_current( true );
			}
		} else { // don't pack this
			task->nonconst_residue_task(i).prevent_repacking();
		}
	}

	std::cout << residues_to_pack << std::endl;
	TR << "packing sidechains" << std::endl;
	pack::rotamer_trials( pose, *packer_scorefxn, task );

}

void
RNP_HighResMover::rnp_docking( core::pose::Pose & pose ) const {

	rnp_docking_mover_->apply( pose );

}


} //movers
} //denovo
} //rna
} //protocols
