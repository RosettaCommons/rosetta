// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/UnfoldedStateEnergyCalculator/UnfoldedStateEnergyCalculatorMover.cc
/// @brief UnfoldedStateEnergyCalculatorMover class definitions
/// @author P. Douglas Renfrew (renfrew@nyu.edu)

#include <protocols/simple_moves/MutateResidue.hh>

// Unit Headers
#include <protocols/unfolded_state_energy_calculator/UnfoldedStateEnergyCalculatorMover.hh>

// Protocol headers
#ifdef USEMPI
#include <protocols/unfolded_state_energy_calculator/UnfoldedStateEnergyCalculatorMPIWorkPoolJobDistributor.hh>
#else
#include <protocols/unfolded_state_energy_calculator/UnfoldedStateEnergyCalculatorJobDistributor.hh>
#endif



// Core headers
#include <core/types.hh>

#include <core/pose/Pose.hh>

#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>

// Utility eaders
#include <utility/vector1.hh>

// Numeric headers
#include <numeric/random/random.hh>



using namespace core;
using namespace core::pose;
using namespace core::chemical;
using namespace core::conformation;
using namespace utility;

namespace protocols {
namespace unfolded_state_energy_calculator {

///@brief dtor
UnfoldedStateEnergyCalculatorMover::~UnfoldedStateEnergyCalculatorMover()
{}

///@brief cctor
UnfoldedStateEnergyCalculatorMover::UnfoldedStateEnergyCalculatorMover( UnfoldedStateEnergyCalculatorMover const & usecm ) :
	//utility::pointer::ReferenceCount(),
	Mover( "UnfoldedStateEnergyCalculatorMover" ),
	job_dist_( usecm.job_dist_ ),
	pack_scrfxn_( usecm.pack_scrfxn_ ),
	score_scrfxn_( usecm.score_scrfxn_ ),
	frag_length_( usecm.frag_length_ ),
	mut_aa_( usecm.mut_aa_ ),
	repack_fragments_( usecm.repack_fragments_ ),
	native_sequence_( usecm.native_sequence_ ),
	sequence_match_sequence_( usecm.sequence_match_sequence_ ),
	sequence_match_position_( usecm.sequence_match_position_ ),
	sequence_matched_fragments_( usecm.sequence_matched_fragments_ )
{}

///@brief alternate ctor
UnfoldedStateEnergyCalculatorMover::UnfoldedStateEnergyCalculatorMover(
#ifdef USEMPI
	UnfoldedStateEnergyCalculatorMPIWorkPoolJobDistributor & job_dist,
#else
	UnfoldedStateEnergyCalculatorJobDistributor & job_dist,
#endif
  core::scoring::ScoreFunctionCOP pack_scrfxn,
	core::scoring::ScoreFunctionCOP score_scrfxn,
	core::Size frag_length,
	std::string mut_aa,
	bool repack_fragments,
	bool native_sequence,
	std::string sequence_match_sequence,
	Size sequence_match_position,
	bool sequence_matched_fragments
):
	Mover( "UnfoldedStateEnergyCalculatorMover" ),
	job_dist_( job_dist ),
	pack_scrfxn_( pack_scrfxn ),
	score_scrfxn_( score_scrfxn ),
	frag_length_( frag_length ),
	mut_aa_( mut_aa ),
	repack_fragments_( repack_fragments ),
	native_sequence_( native_sequence ),
	sequence_match_sequence_( sequence_match_sequence ),
	sequence_match_position_( sequence_match_position ),
	sequence_matched_fragments_( sequence_matched_fragments )
{}

protocols::moves::MoverOP
UnfoldedStateEnergyCalculatorMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new UnfoldedStateEnergyCalculatorMover( *this ) );
}

void
UnfoldedStateEnergyCalculatorMover::create_random_fragments( Pose & pose, vector1< Pose > & fragments )
{
	// get number of protein residues
	Size num_protein_res( 0 );
	for ( Size i = 1; i < pose.total_residue(); ++i ) {
		if ( pose.residue( i ).type().is_protein() )
			num_protein_res++;
	}

	// get number of fragments
	Size frag_number( num_protein_res / frag_length_ ); // this seems to give a reasonable number of fragments

	// check frag number
	if ( frag_number == 0 ) {
		return;
	}

	// create fragments
	do {

		// generate a random start location within this structure; use a random number whose range does not include the termini.
		Size frag_start( numeric::random::random_range( 1,  num_protein_res - frag_length_ + 1 ) );

		// if the fragment is bad
		if ( !fragment_check( pose, frag_start ) ) continue;

		// correct varient types (termini, disulfides, etc.)

		// create fragment pose
		Pose frag_temp;
		frag_temp.append_residue_by_jump( pose.residue( frag_start ), 1 );
		for ( Size i(frag_start + 1); i < frag_start + frag_length_; ++i ) {
			frag_temp.append_residue_by_bond( pose.residue( i ) );
		}

		// add to vector of fragments
		fragments.push_back( frag_temp );
	} while ( fragments.size() < frag_number );

}

void
UnfoldedStateEnergyCalculatorMover::create_sequence_match_fragments( Pose & pose, vector1< Pose > & fragments )
{

	// find all the places in the sequence that match our input sequence match
	utility::vector1< Size > match_start_pos;
	Size temp( std::string::npos );
	temp = pose.sequence().find( sequence_match_sequence_ );
	while ( temp != std::string::npos ) {
		match_start_pos.push_back( temp );
		temp = pose.sequence().find( sequence_match_sequence_, temp+1 );
	}

	// use frag length and central residue number to determine if fragment is good
	for ( Size i( 1 ); i <= match_start_pos.size(); ++i ) {

		// calculate start position
		Size frag_start( match_start_pos[ i ] + sequence_match_position_ - frag_length_ / 2 );

		// if the fragment is bad
		if ( !fragment_check( pose, frag_start) ) continue;

		// correct varient types (termini, disulfides, etc.)

		// create fragment pose
		Pose frag_temp;
		frag_temp.append_residue_by_jump( pose.residue( frag_start ), 1 );
		for ( Size i(frag_start + 1); i < frag_start + frag_length_; ++i ) {
			frag_temp.append_residue_by_bond( pose.residue( i ) );
		}

		// add to vector of fragments
		fragments.push_back( frag_temp );
	}
}

/// @brief Check fragment residues should be polymeric and on same chain
bool
UnfoldedStateEnergyCalculatorMover::fragment_check( Pose & pose, Size frag_start )
{
	// do not go past the begining or end of the pose
	if ( frag_start <= 0 ) return false;
	if ( frag_start + frag_length_ > pose.total_residue() )	return false;

	// all residues should be on the same chain
	for ( Size i( frag_start ); i <= frag_start + frag_length_; ++i ) {
		if ( pose.chain( frag_start ) != pose.chain( i ) ) return false;
	}

	// all residues should be polymeric
	for ( Size i( frag_start ); i <= frag_start + frag_length_; ++i ) {
		if ( !pose.residue( i ).type().is_polymer() ) return false;
	}

	// all residues should not be nucleic acid
	for ( Size i( frag_start ); i <= frag_start + frag_length_; ++i ) {
		if ( pose.residue( i ).type().is_NA() ) return false;
	}

	return true;
}

void
UnfoldedStateEnergyCalculatorMover::apply( Pose & pose )
{

	// create fragments vector
	vector1< Pose > fragments;

	// get central residue size
	Size frag_cent_res_number( ( frag_length_ / 2 ) + 1 );

	// create sequence matched frgments or random fragments
	if ( sequence_matched_fragments_ ) {
		create_sequence_match_fragments( pose, fragments );
	} else {
		create_random_fragments( pose, fragments );
	}

	// mutate central residue if native_sequence flag is false
	if ( !native_sequence_ ) {
		for ( Size i( 1 ); i <= fragments.size(); ++i ) {
			protocols::simple_moves::MutateResidueOP mut_res_mvr( new protocols::simple_moves::MutateResidue( frag_cent_res_number, mut_aa_ ) );
			mut_res_mvr->apply( fragments[i] );
		}
	}

	// check repack fragment central residues
	if ( repack_fragments_ ) {
		for ( Size i(1); i <= fragments.size(); ++i ) {

			// create packer task
			pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( fragments[ i ] ) );
			task->initialize_from_command_line();

			// set all residues in fragment to only repack
			for ( Size j(1); j <= frag_length_; ++j ) {
				task->nonconst_residue_task( j ).restrict_to_repacking();
			}

			// repack fragment
			pack::pack_rotamers( fragments[ i ], *pack_scrfxn_, task);
		}
	}

	/*
	// DEBUG output pdbs of fragments
	for ( Size i( 1 ); i <= fragments.size(); ++i ) {
	 	std::stringstream outputfilename;
	 	outputfilename << "pack_frag_" << i;
	 	for ( Size j( 1 ); j <= frag_length_; ++j ) {
	 		outputfilename << "_" << fragments[ i ].residue( j ).type().name3();
	 	}
	 	outputfilename << ".pdb";
	 	fragments[i].dump_scored_pdb( outputfilename.str(), *score_scrfxn_ );
	}
	*/

	// check fragment scores
	for ( Size i(1); i <= fragments.size(); ++i ) {

		// clear energies
		( fragments[ i ] ).energies().clear();

		// score with scoring score function
		( *score_scrfxn_ )( fragments[ i ] );

		// send central residue energies to job distributor
		job_dist_.add_unfolded_energy_data( fragments[ i ].residue(frag_cent_res_number).type().name3(), fragments[ i ].energies().residue_total_energies( frag_cent_res_number ) );

	}

	// send weights set to job distributor
	job_dist_.set_energy_terms( ( *score_scrfxn_ ).weights() );

}

std::string
UnfoldedStateEnergyCalculatorMover::get_name() const {
	return "UnfoldedStateEnergyCalculatorMover";
}

bool
UnfoldedStateEnergyCalculatorMover::reinitialize_for_each_job() const
{
	return false;
}

bool
UnfoldedStateEnergyCalculatorMover::reinitialize_for_new_input() const
{
	return false;
}

} // UnfoldedStateEnergyCalculator
} // protocols
