// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/scoring/IsopeptideBondTests.cxxtest.hh
/// @brief  Unit tests confirming that mainchain potential scoring works properly with isopeptide-bonded residues.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers


// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>

// Protocols Headers (for convenience)
#include <protocols/cyclic_peptide/PeptideStubMover.hh> //To build peptides easily
#include <protocols/cyclic_peptide/DeclareBond.hh> //Quick and dirty way to correct termini variants.

// Utility, etc Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("IsopeptideBondTests");


class IsopeptideBondTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init_with_additional_options( "-write_all_connect_info" );
		phipsi_.push_back( std::pair<core::Real, core::Real>( -61.0, -41.0 ) );
		phipsi_.push_back( std::pair<core::Real, core::Real>( -63.0, -38.0 ) );
		phipsi_.push_back( std::pair<core::Real, core::Real>( -135, -135 ) );
		phipsi_.push_back( std::pair<core::Real, core::Real>( 61, 41 ) );
		phipsi_.push_back( std::pair<core::Real, core::Real>( 63, 43 ) );
		phipsi_.push_back( std::pair<core::Real, core::Real>( 135, 135 ) );
		phipsi_.push_back( std::pair<core::Real, core::Real>( -65, 120 ) );
		phipsi_.push_back( std::pair<core::Real, core::Real>( 65, -120 ) );
	}

	void tearDown() {}


	void common_test_fxn_nterm( core::scoring::ScoreFunctionOP sfxn, core::scoring::ScoreType const scoretype ) {
		core::pose::Pose pose1, pose2;

		{
			protocols::cyclic_peptide::PeptideStubMover builder1;
			builder1.add_residue( "Append", "ALA", 1, true, "", 0, 1, "" );
			builder1.add_residue( "Append", "ALA", 2, false, "", 0, 1, "" );
			builder1.add_residue( "Append", "ALA", 3, false, "", 0, 2, "" );
			builder1.apply(pose1);
			protocols::cyclic_peptide::DeclareBond fix_termini;
			fix_termini.set( 1, "C", 2, "N", true );
			fix_termini.apply(pose1);
		}
		{
			protocols::cyclic_peptide::PeptideStubMover builder2;
			builder2.add_residue( "Append", "ALA", 1, true, "", 0, 1, "" );
			builder2.add_residue( "Append", "ALA", 2, false, "", 0, 1, "" );
			builder2.add_residue( "Append", "GLU:SidechainConjugation", 3, false, "CD", 0, 1, "N" );
			builder2.apply(pose2);
			protocols::cyclic_peptide::DeclareBond fix_termini;
			fix_termini.set( 1, "C", 2, "N", true );
			fix_termini.apply(pose2);
		}

		TR << "\nN-terminal connection\n";
		TR << "PHI\tPSI\t" << core::scoring::name_from_score_type( scoretype ) << "\tCOMPARISON\n";
		for ( core::Size i(1), imax(phipsi_.size()); i<=imax; ++i ) {
			pose1.set_phi( 2, phipsi_[i].first );
			pose1.set_psi( 2, phipsi_[i].second );
			pose1.set_omega( 2, 180 );
			pose2.set_phi( 1, phipsi_[i].first );
			pose2.set_psi( 1, phipsi_[i].second );
			pose2.set_omega( 1, 180 );
			(*sfxn)(pose1);
			(*sfxn)(pose2);
			TR << phipsi_[i].first << "\t" << phipsi_[i].second << "\t" << pose1.energies().total_energy() << "\t" << pose2.energies().total_energy() << "\n";
			TS_ASSERT_DELTA( pose1.energies().total_energy(), pose2.energies().total_energy(), 1e-5 );
			//DELETE THE FOLLOWING -- FOR DEBUGGING ONLY
			/*char outfile [256];
			sprintf( outfile, "RAMACTRL_N%04lu.pdb", i );
			pose1.dump_pdb( std::string(outfile) );
			sprintf( outfile, "RAMATEST_N%04lu.pdb", i );
			pose2.dump_pdb( std::string(outfile) );*/
		}
		TR << std::endl;
	}

	void common_test_fxn_cterm( core::scoring::ScoreFunctionOP sfxn, core::scoring::ScoreType const scoretype ) {
		core::pose::Pose pose1, pose2;

		{
			protocols::cyclic_peptide::PeptideStubMover builder1;
			builder1.add_residue( "Append", "ALA", 1, true, "", 0, 1, "" );
			builder1.add_residue( "Append", "ALA", 2, false, "", 0, 1, "" );
			builder1.add_residue( "Append", "ALA", 3, false, "", 0, 2, "" );
			builder1.apply(pose1);
			protocols::cyclic_peptide::DeclareBond fix_termini;
			fix_termini.set( 1, "C", 2, "N", true );
			fix_termini.apply(pose1);
		}
		{
			protocols::cyclic_peptide::PeptideStubMover builder2;
			builder2.add_residue( "Append", "ALA", 1, true, "", 0, 1, "" );
			builder2.add_residue( "Append", "ALA", 2, false, "", 0, 1, "" );
			builder2.add_residue( "Append", "DPP:SidechainConjugation", 3, false, "NG", 0, 2, "C" );
			builder2.apply(pose2);
			protocols::cyclic_peptide::DeclareBond fix_termini;
			fix_termini.set( 1, "C", 2, "N", true );
			fix_termini.apply(pose2);
		}

		TR << "\nC-terminal connection\n";
		TR << "PHI\tPSI\t" << core::scoring::name_from_score_type( scoretype ) << "\tCOMPARISON\n";
		for ( core::Size i(1), imax(phipsi_.size()); i<=imax; ++i ) {
			pose1.set_phi( 2, phipsi_[i].first );
			pose1.set_psi( 2, phipsi_[i].second );
			pose1.set_omega( 2, 180 );
			pose2.set_phi( 2, phipsi_[i].first );
			pose2.set_psi( 2, phipsi_[i].second );
			pose2.set_omega( 2, 180 );
			(*sfxn)(pose1);
			(*sfxn)(pose2);
			TR << phipsi_[i].first << "\t" << phipsi_[i].second << "\t" << pose1.energies().total_energy() << "\t" << pose2.energies().total_energy() << "\n";
			TS_ASSERT_DELTA( pose1.energies().total_energy(), pose2.energies().total_energy(), 1e-5 );
			//DELETE THE FOLLOWING -- FOR DEBUGGING ONLY
			/*char outfile [256];
			sprintf( outfile, "RAMACTRL_C%04lu.pdb", i );
			pose1.dump_pdb( std::string(outfile) );
			sprintf( outfile, "RAMATEST_C%04lu.pdb", i );
			pose2.dump_pdb( std::string(outfile) );*/
		}
		TR << std::endl;
	}

	void common_test_fxn( core::scoring::ScoreType const scoretype ) {
		core::scoring::ScoreFunctionOP sfxn( new core::scoring::ScoreFunction );
		sfxn->set_weight( scoretype, 1.0 );
		common_test_fxn_nterm(sfxn, scoretype);
		common_test_fxn_cterm(sfxn, scoretype);
	}

	void test_rama_prepro(){
		common_test_fxn( core::scoring::rama_prepro );
	}

	void test_p_aa_pp(){
		common_test_fxn( core::scoring::p_aa_pp );
	}

private:

	/// @brief The phi, psi pairs to test.
	utility::vector1< std::pair< core::Real, core::Real > > phipsi_;

};
