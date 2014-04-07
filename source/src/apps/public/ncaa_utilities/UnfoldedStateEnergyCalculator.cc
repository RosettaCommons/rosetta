// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/apps/pilot/doug/UnfoldedStateEnergyCalculator.cc
/// @brief UnfoldedStateEnergyCalculator application
/// @author P. douglas Renfrew (renfrew@unc.edu)

/* WARNING WARNING WARNING
 * This code in intended to run over a large set of pdb files. Not all PDBs are well behaved with mini.
 * This code will not run well without some other changes to mini. To use it, you are strongly
 * encouraged to "robustify" mini.  This will cause bad PDBs to be ignored rather than causing crashes!
 *
 * To robustify mini:
 * replace all assert statements in the vectorL (vector1) class with runtime_assert statements
 * replace all assert statements in the Conformation class with runtime_assert statements
 */

// MPI headers
#ifdef USEMPI
#include <mpi.h> //keep this first
#endif

// Unit Headers
#include <protocols/unfolded_state_energy_calculator/UnfoldedStateEnergyCalculatorMover.hh>

// Package headers

// Project headers
#include <core/types.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <devel/init.hh>
#include <protocols/unfolded_state_energy_calculator/UnfoldedStateEnergyCalculatorJobDistributor.hh>
// AUTO-REMOVED #include <protocols/unfolded_state_energy_calculator/UnfoldedStateEnergyCalculatorMPIWorkPoolJobDistributor.hh>

#include <basic/options/option.hh>

// Utility Headers
#include <devel/init.hh>
#include <basic/Tracer.hh>

// C++ Headers
#include <iostream>
#include <string>

#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>

using namespace core;
using namespace core::scoring;
using namespace basic::options;
using namespace basic::options::OptionKeys;

using namespace protocols::unfolded_state_energy_calculator;

using basic::T;

// Local options
namespace usec {
	IntegerOptionKey const frag_size( "usec::frag_size" );
	StringOptionKey const residue_name( "usec::residue_name" );
	BooleanOptionKey const repack_fragments( "usec::repack_fragments" );
	BooleanOptionKey const native_sequence( "usec::native_sequence" );
}

int main( int argc, char* argv[] )
{
	try {
 	// add application specific options to options system
	option.add( usec::frag_size, "Sets the number of residues in each fragment, should be an odd number" ).def( 5 );
	option.add( usec::residue_name, "Sets the three letter code of the residue typpe which the central residue will be mutated to" ).def( "LEU" );
	option.add( usec::repack_fragments, "Controls if the fragments will be repacked before scoring" ).def( true );
	option.add( usec::native_sequence, "Controls if the central residue will be mutated before scoring" ).def( false );

	// init
	devel::init(argc, argv);

  // setup score function for packing fragments
	ScoreFunctionOP packing_scrfxn( ScoreFunctionFactory::create_score_function( MM_STD_WTS ) );
	packing_scrfxn->set_weight( unfolded, 0.0 );

  core::scoring::methods::EnergyMethodOptions pemo( packing_scrfxn->energy_method_options() );
  pemo.hbond_options().decompose_bb_hb_into_pair_energies( true);
  packing_scrfxn->set_energy_method_options( pemo );

  // setup score function for scoring fragments
  ScoreFunctionOP scoring_scrfxn( ScoreFunctionFactory::create_score_function( MM_STD_WTS ) );
	scoring_scrfxn->set_weight( unfolded, 0.0 );

  core::scoring::methods::EnergyMethodOptions semo( scoring_scrfxn->energy_method_options() );
  semo.hbond_options().decompose_bb_hb_into_pair_energies( true );
  scoring_scrfxn->set_energy_method_options( semo );

	// get job distributor based on whether or not we are using MPI
#ifdef USEMPI
  UnfoldedStateEnergyCalculatorMPIWorkPoolJobDistributor jd;
#else
  UnfoldedStateEnergyCalculatorJobDistributor jd;
#endif

  // setup UnfoldedStateEnergyCalculatorMover
	Size frag_length( option[ usec::frag_size ].value() );
	std::string mut_aa( option[ usec::residue_name ].value() );
	bool repack_fragments( option[ usec::repack_fragments ].value() );
	bool native_sequence( option[ usec::native_sequence ].value() );
  UnfoldedStateEnergyCalculatorMoverOP usecm( new UnfoldedStateEnergyCalculatorMover( jd, packing_scrfxn, scoring_scrfxn, frag_length, mut_aa, repack_fragments, native_sequence ) );

	// call job distributor with mover
	jd.go( usecm );

  T("UnfoldedStateEnergyCalculator") << "\n+-----------------------------------------------------------------+\n"
																		 <<   "|                              DONE                               |\n"
																		 <<   "+-----------------------------------------------------------------+" << std::endl;
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
  return 0;
}
