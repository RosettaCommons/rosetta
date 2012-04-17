// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is protocolsoped by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/UnfoldedStateEnergyCalculator/UnfoldedStateEnergyCalculatorMover.cc
/// @brief UnfoldedStateEnergyCalculatorMover class definitions
/// @author P. douglas Renfrew (renfrew@unc.edu)

// Unit Headers
#include <protocols/unfolded_state_energy_calculator/UnfoldedStateEnergyCalculatorMover.hh>

// Package headers

// Project headers
#include <core/types.hh>

#include <core/pose/Pose.hh>

#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>

#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>

// Utility Headers
#include <utility/vector1.hh>

// Numeric headers
#include <numeric/random/random.hh>

#include <protocols/unfolded_state_energy_calculator/UnfoldedStateEnergyCalculatorJobDistributor.hh>
#include <utility/vector0.hh>


// C++ Headers

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
	native_sequence_( usecm.native_sequence_ )
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
	bool native_sequence ):
	Mover( "UnfoldedStateEnergyCalculatorMover" ),
	job_dist_( job_dist ),
	pack_scrfxn_( pack_scrfxn ),
	score_scrfxn_( score_scrfxn ),
	frag_length_( frag_length ),
	mut_aa_( mut_aa ),
	repack_fragments_( repack_fragments ),
	native_sequence_( native_sequence )
{}

protocols::moves::MoverOP
UnfoldedStateEnergyCalculatorMover::fresh_instance() const
{
	return new UnfoldedStateEnergyCalculatorMover( *this );
}

void
UnfoldedStateEnergyCalculatorMover::apply( Pose & pose )
{

	// get number of protein residues
	Size num_protein_res( 0 );
	for ( Size i = 1; i < pose.total_residue(); ++i ) {
		if ( pose.residue( i ).type().is_protein() )
			num_protein_res++;
	}
	//std::cout << "NUM PROT RES: " << num_protein_res  << std::endl; //debug

	// get number of fragments
	Size frag_number( num_protein_res / frag_length_ );
	//std::cout << "FRAG NUM: " << frag_number  << std::endl; //debug

	// get central residue size
	Size frag_cent_res_number( ( frag_length_ / 2 ) + 1 );
	//std::cout << "FRAG NUM CENT RES: " << frag_cent_res_number << std::endl; //debug

	// check frag number
	if ( frag_number == 0 ) {
		return;
	}

	// create fragments vector
	vector1< Pose > fragments;

	// get residue type set and a copy of the residue to mutate to
	ResidueTypeSetCAP rts( chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" ) );
	ResidueType mut_aa_type( rts->name_map( mut_aa_ ) );
	ResidueOP mut_aa_res( ResidueFactory::create_residue( mut_aa_type ) );

	// create fragments
	do {

		// generate a random start location within this structure; use a random number whose range does not include the termini.
		Size frag_start( numeric::random::random_range( 1,  num_protein_res - frag_length_ + 1 ) );

		// check fragment residues should be polymeric and on same chain (and central residue should be correct type if replace typr is turned on)
		bool frag_check(true);

		// all residues should be on the same chain
		for ( Size i( frag_start ); i <= frag_start + frag_length_; ++i ) {
			if ( pose.chain( frag_start ) != pose.chain( i ) ) frag_check = false;
		}

		// all residues should be polymeric
		for ( Size i( frag_start ); i <= frag_start + frag_length_; ++i ) {
			if ( !pose.residue( i ).type().is_polymer() ) frag_check = false;
		}

		// all residues should not be nucleic acid
		for ( Size i( frag_start ); i <= frag_start + frag_length_; ++i ) {
			if ( pose.residue( i ).type().is_NA() ) frag_check = false;
		}

		// if replace_res is on central residue should be the same type as the residue being replaced
		// hardcoded for proteins, peptoids and beta peptides currently, don't know of a good way compare
		// compatible properties

		if ( !native_sequence_ ) {
			if ( pose.residue( frag_cent_res_number ).type().has_property("PROTEIN") && !mut_aa_type.has_property("PROTEIN") ) frag_check = false;
			if ( pose.residue( frag_cent_res_number ).type().has_property("BETA_PEPTIDE") && !mut_aa_type.has_property("BETA_PEPTIDE") ) frag_check = false;
			if ( pose.residue( frag_cent_res_number ).type().has_property("PEPTOID") && !mut_aa_type.has_property("PEPTOID") ) frag_check = false;
		}

		// // DEBUG
		// std::cout << "FRAG:\t" << fragments.size() << "\t" << frag_check << "\t" << std::flush;
		// for ( Size i( frag_start ); i <= frag_start + frag_length_; ++i ) {
		// 	std::cout << i << ":" << pose.residue( i ).type().name() << ":" << pose.residue( i ).chain() << "\t" << std::flush;
		// }
		// std::cout << std::endl;

		// if the fragment is bad
		if ( !frag_check ) continue;

		// correct varient types (termini, disulfides, etc.)

		// create fragment pose
		Pose frag_temp;
		//std::cout << "BUILD FRAG START" << std::endl;
		//std::cout << "DEBUG:\t" << frag_start << "\t" << pose.residue( frag_start ).type().name() << "\t" << pose.residue( frag_start ).chain() << std::endl;
		frag_temp.append_residue_by_jump( pose.residue( frag_start ), 1 );
		for ( Size i(frag_start + 1); i < frag_start + frag_length_; ++i ) {
			//std::cout << "DEBUG:\t" << i << "\t" <<pose.residue( i ).type().name() << "\t" << pose.residue( i ).chain() << std::endl;
			frag_temp.append_residue_by_bond( pose.residue( i ) );
		}

		// mutate central residue if native_sequence flag is false
		if ( !native_sequence_ ) {
			frag_temp.replace_residue( frag_cent_res_number, *mut_aa_res, true );
		}

		// add to vector of fragments
		fragments.push_back( frag_temp);
	} while ( fragments.size() < frag_number );

	// // DEBUG output pdbs of fragments
	// for ( Size i( 1 ); i <= fragments.size(); ++i ) {
	// 	std::stringstream outputfilename;
	// 	outputfilename << "frag_" << i;
	// 	for ( Size j( 1 ); j <= frag_length_; ++j ) {
	// 		outputfilename << "_" << fragments[ i ].residue( j ).type().name3();
 	// 	}
	// 	outputfilename << ".pdb";
	// 	fragments[i].dump_pdb( outputfilename.str() );
	// }


	// check repack fragment central residues
	if ( repack_fragments_ ) {
		for ( Size i(1); i <= frag_number; ++i ) {

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

	// // DEBUG output pdbs of fragments
	// for ( Size i( 1 ); i <= fragments.size(); ++i ) {
	// 	std::stringstream outputfilename;
	// 	outputfilename << "pack_frag_" << i;
	// 	for ( Size j( 1 ); j <= frag_length_; ++j ) {
	// 		outputfilename << "_" << fragments[ i ].residue( j ).type().name3();
	// 	}
	// 	outputfilename << ".pdb";
	// 	fragments[i].dump_scored_pdb( outputfilename.str(), *score_scrfxn_ );
	// }

	// check fragment scores
	for ( Size i(1); i <= frag_number; ++i ) {

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
