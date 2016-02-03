// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/optimize_weights/DDGBindOptEData.cxxtest.hh
/// @brief  test suite for the optE data class which holds data related to ddG bind optimization
/// @author Ron Jacak

// Test framework headers
#include <cxxtest/TestSuite.h>

// Core Headers
//#include <core/conformation/Residue.hh>
//#include <core/scoring/TenANeighborGraph.hh>
//#include <core/pack/packer_neighbors.hh>
//#include <basic/options/util.hh>

#include <protocols/optimize_weights/OptEData.hh>
#include <protocols/optimize_weights/DDGBindOptEData.hh>

#include <core/chemical/AA.hh> // to get aa_from_oneletter_code()
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreType.hh>

#include <utility/vector1.hh>

#include <ObjexxFCL/format.hh>

// Utility Headers
#include <basic/Tracer.hh>

// Numeric headers

// Test headers
#include <test/core/init_util.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>


static basic::Tracer TR("test.protocols.optimize_weights.DDGBindOptEData");

using namespace protocols;
using namespace protocols::optimize_weights;


// --------------- Test Class --------------- //

class DDGBindOptEDataTests : public CxxTest::TestSuite {

public:

	// Shared data elements go here.
	protocols::optimize_weights::DDGBindOptEDataOP ddg_bind_position_data;
	core::pose::Pose wt_complex, mut_complex, wt_unbounded, mut_unbounded;
	core::scoring::ScoreFunctionOP scorefxn;

	SingleStructureDataOP ssd;
	utility::vector1< core::Real > component_weights;

	core::optimization::Multivec dofs;
	core::optimization::Multivec dE_dvars;

	core::scoring::EnergyMap free_parameters, fixed_parameters; // energy maps
	core::scoring::ScoreTypes free_score_list, fixed_score_list; // list of ScoreType objects


	// --------------- Suite-level Fixture --------------- //

	DDGBindOptEDataTests() {

		// if the tests are run manually (or one suite at a time), that doesn't mute all of the tracer output by default.  Place
		// a mute here because the interaction graphs generate tons of debugging output (in DEBUG mode anyway).
		core_init_with_additional_options( "-no_optH -mute core.io core.init core.scoring core.mm core.pack.task" );


		// To create a DDGBindOptEDataTests object, we need to create a few other objects like a Pose, a ScoreFunction, etc
		// Create all of these objects here in the suite-level fixture since they'll get reused throughout the suite.


		// --- ScoreFunction ---
		// create a score function using the standard packer weights
		TR << "creating sfxn for calculating ddGs of binding" << std::endl;
		scorefxn = core::scoring::get_score_function();
		scorefxn->set_weight( core::scoring::surface, 0.5 );

		// pretend the input file looks as follows
		// 1UAD.wt_complex.pdb   1UADA.E38A.mut_complex.pdb   1UAD.wt_unbounded.pdb   1UADA.E38A.mut_unbounded.pdb   1.99
		// this isn't what the input file actually looks like, but it's close enough
		TR << "creating OptEData object" << std::endl;
		ddg_bind_position_data = protocols::optimize_weights::DDGBindOptEDataOP( new DDGBindOptEData );

		// save the experimental ddg for this mutant
		ddg_bind_position_data->set_experimental_ddg_bind( 1.99 );
		ddg_bind_position_data->tag( "1UAD.E38A.test_case" );

		TR << "reading in test PDBs" << std::endl;
		core::import_pose::pose_from_file( wt_complex, "protocols/optimize_weights/1UAD.wt_complex.pdb" , core::import_pose::PDB_file);
		core::import_pose::pose_from_file( mut_complex, "protocols/optimize_weights/1UAD.E38A.mut_complex.pdb" , core::import_pose::PDB_file);
		core::import_pose::pose_from_file( wt_unbounded, "protocols/optimize_weights/1UAD.wt_unbounded.pdb" , core::import_pose::PDB_file);
		core::import_pose::pose_from_file( mut_unbounded, "protocols/optimize_weights/1UAD.E38A.mut_unbounded.pdb" , core::import_pose::PDB_file);

		// sequences will be used to set the list of mutated amino acids held in the optE data object
		std::string wt_complex_seq = wt_complex.sequence();
		std::string mut_complex_seq = mut_complex.sequence();

		TR << "identifying mutated residues and setting mutation field in DDGBindOptEData object" << std::endl;
		// set the mutations field in the position data object
		for ( core::Size jj = 0; jj < wt_complex_seq.size(); ++jj ) {
			if ( wt_complex_seq[ jj ] != mut_complex_seq[ jj ] ) {
				ddg_bind_position_data->add_mutation(
					std::make_pair( jj, std::make_pair( core::chemical::aa_from_oneletter_code( wt_complex_seq[ jj ] ), core::chemical::aa_from_oneletter_code( mut_complex_seq[ jj ] ) ) ) );
			}
		}

		TR << "creating EMaps to hold free and fixed weights" << std::endl;

		free_parameters[ core::scoring::fa_atr ] = 0.8;
		free_parameters[ core::scoring::fa_rep ] = 0.44;
		free_parameters[ core::scoring::fa_sol ] = 0.65;
		free_parameters[ core::scoring::fa_intra_rep ] = 0.004;
		free_parameters[ core::scoring::fa_pair ] = 0.49;
		free_parameters[ core::scoring::fa_dun ] = 0.56;
		free_parameters[ core::scoring::hbond_lr_bb ] = 1.17;
		free_parameters[ core::scoring::hbond_sr_bb ] = 0.585;
		free_parameters[ core::scoring::hbond_bb_sc ] = 1.17;
		free_parameters[ core::scoring::hbond_sc ] = 1.1;
		free_parameters[ core::scoring::rama ] = 0.2;
		free_parameters[ core::scoring::p_aa_pp ] = 0.32;
		free_parameters[ core::scoring::pro_close ] = 1.0;
		free_parameters[ core::scoring::surface ] = 1.0;

		fixed_parameters[ core::scoring::omega ] = 0.5;
		fixed_parameters[ core::scoring::dslf_ss_dst ] = 1.0;
		fixed_parameters[ core::scoring::dslf_cs_ang ] = 1.0;
		fixed_parameters[ core::scoring::dslf_ss_dih ] = 1.0;
		fixed_parameters[ core::scoring::dslf_ca_dih ] = 1.0;

		core::Real const refEs[20] = { 0.16, 1.7, -0.67, -0.81, 0.63, -0.17, 0.56, 0.24, -0.65, -0.1, -0.34, -0.89, 0.02, -0.97, -0.98, -0.37, -0.27, 0.29, 0.91, 0.51 };

		utility::vector1< core::Real > reference_energies( core::chemical::num_canonical_aas );
		for ( core::Size ii = 1; ii <= core::chemical::num_canonical_aas; ++ii ) {
			reference_energies[ ii ] = refEs[ ii - 1 ];
		}

		core::scoring::EnergyMap include_terms;
		include_terms = free_parameters;
		include_terms += fixed_parameters;

		TR << "setting up ScoreTypes variables, list of ScoreType objects that denote what terms are in use" << std::endl;
		for ( int ii = 1 ; ii <= core::scoring::n_score_types ; ++ii ) {
			if ( include_terms[ ( core::scoring::ScoreType )ii ] != 0.0 ) {
				if ( fixed_parameters[ core::scoring::ScoreType( ii ) ] == 0.0 ) {
					free_score_list.push_back(  ( core::scoring::ScoreType )ii );
				} else {
					fixed_score_list.push_back( ( core::scoring::ScoreType )ii );
				}
			}
		}

		TR << "setting up dofs/vars and dE_dvars arrays" << std::endl;
		dofs.resize( 34 /* total number of dofs */, 0.0 );
		dE_dvars.resize( 34 /* total number of dofs */, 0.0 );

		// num_total_dofs is the number of free energy term weights (num_energy_dofs) plus the reference energy dofs
		core::Size dof_index = 1;
		for ( core::scoring::ScoreTypes::const_iterator itr = free_score_list.begin(), end_itr = free_score_list.end(); itr != end_itr; ++itr ) {
			dofs[ dof_index++ ] = include_terms[ *itr ];
		}
		for ( core::Size ii = 1; ii <= reference_energies.size(); ++ii ) {
			dofs[ dof_index++ ] = reference_energies[ ii ];
		}

		TR << "vars: ";
		for ( core::Size ii = 1; ii <= dofs.size(); ++ii ) { TR << ObjexxFCL::format::F( 6,3,dofs[ii] ) << ", "; }
		TR << "]" << std::endl;

		TR << "making component weights file" << std::endl;
		component_weights.resize( n_optE_data_types );
		std::fill( component_weights.begin(), component_weights.end(), 1.0 );

		TR << "filling OptEData object with unweighted energies" << std::endl;

		utility::vector1< core::Real > free_data( free_score_list.size(), 0.0 );
		utility::vector1< core::Real > fixed_data( fixed_score_list.size(), 0.0 );
		(*scorefxn)( wt_complex );
		for ( core::Size kk = 1; kk <= free_score_list.size(); ++kk ) {
			free_data[ kk ] = wt_complex.energies().total_energies()[ free_score_list[ kk ] ];
		}
		for ( core::Size kk = 1; kk <= fixed_score_list.size(); ++kk ) {
			fixed_data[ kk ] = wt_complex.energies().total_energies()[ fixed_score_list[ kk ] ];
		}

		ssd = SingleStructureDataOP( new SingleStructureData( free_data, fixed_data ) );
		ddg_bind_position_data->add_wt_complex( ssd );
		free_data.clear(); fixed_data.clear(); ssd = NULL;

		free_data.resize( free_score_list.size(), 0.0 );
		fixed_data.resize( fixed_score_list.size(), 0.0 );
		(*scorefxn)( mut_complex );
		for ( core::Size kk = 1; kk <= free_score_list.size(); ++kk ) {
			free_data[ kk ] = mut_complex.energies().total_energies()[ free_score_list[ kk ] ];
		}
		for ( core::Size kk = 1; kk <= fixed_score_list.size(); ++kk ) {
			fixed_data[ kk ] = mut_complex.energies().total_energies()[ fixed_score_list[ kk ] ];
		}

		ssd = SingleStructureDataOP( new SingleStructureData( free_data, fixed_data ) );
		ddg_bind_position_data->add_mutant_complex( ssd );
		free_data.clear(); fixed_data.clear(); ssd = NULL;


		free_data.resize( free_score_list.size(), 0.0 );
		fixed_data.resize( fixed_score_list.size(), 0.0 );
		(*scorefxn)( wt_unbounded );
		for ( core::Size kk = 1; kk <= free_score_list.size(); ++kk ) {
			free_data[ kk ] = wt_unbounded.energies().total_energies()[ free_score_list[ kk ] ];
		}
		for ( core::Size kk = 1; kk <= fixed_score_list.size(); ++kk ) {
			fixed_data[ kk ] = wt_unbounded.energies().total_energies()[ fixed_score_list[ kk ] ];
		}

		ssd = SingleStructureDataOP( new SingleStructureData( free_data, fixed_data ) );
		ddg_bind_position_data->add_wt_unbounds( ssd );
		free_data.clear(); fixed_data.clear(); ssd = NULL;


		free_data.resize( free_score_list.size(), 0.0 );
		fixed_data.resize( fixed_score_list.size(), 0.0 );
		(*scorefxn)( mut_unbounded );
		for ( core::Size kk = 1; kk <= free_score_list.size(); ++kk ) {
			free_data[ kk ] = mut_unbounded.energies().total_energies()[ free_score_list[ kk ] ];
		}
		for ( core::Size kk = 1; kk <= fixed_score_list.size(); ++kk ) {
			fixed_data[ kk ] = mut_unbounded.energies().total_energies()[ fixed_score_list[ kk ] ];
		}

		ssd = SingleStructureDataOP( new SingleStructureData( free_data, fixed_data ) );
		ddg_bind_position_data->add_mutant_unbounds( ssd );
		free_data.clear(); fixed_data.clear(); ssd = NULL;

		TR << "---" << std::endl;
	}

	virtual ~DDGBindOptEDataTests() {}

	static DDGBindOptEDataTests *createSuite() {
		return new DDGBindOptEDataTests();
	}

	static void destroySuite( DDGBindOptEDataTests *suite ) {
		delete suite;
	}


	// --------------- Test Fixture --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case.

	void setUp() {}

	// Shared finalization goes here.
	// All memory allocated via OPs; objects should destroy themselves so nothing else to do here.
	void tearDown() {}


public:

	// --------------- Test Cases --------------- //


	/// @details
	/// Tests the function get_score() in the DDGBindOptEData class.
	/// get_score() calls process_score() which takes the input array of weights and applies them to the unweighted
	/// energies for all of the structures contained within.
	/// it calculates a predicted ddG of binding and returns the error between that and the experimental one, squared
	///
	void test_get_score() {

		TR << "Running test_get_score..." << std::endl;

		// Test a couple of different possible weight sets
		core::Real score = 0.0;
		std::ofstream outlog( "dummy_logfile" );
		ddg_bind_position_data->print_score( outlog, component_weights, dofs, dE_dvars, 13 /* free_score_list.size() */,
			20 /* chemical::num_canonical_aas */,  33 /* free_score_list.size() + chemical::num_canonical_aas */, fixed_parameters, free_score_list, fixed_score_list );
		score = ddg_bind_position_data->get_score( component_weights, dofs, dE_dvars, 13 /* free_score_list.size() */,
			20 /* chemical::num_canonical_aas */,  33 /* free_score_list.size() + chemical::num_canonical_aas */, fixed_parameters, free_score_list, fixed_score_list );
		TS_ASSERT_DELTA( score, 4.1140, 0.1 );

		// also check to make sure the change in the reference energy for all of the mutated residues was taken into account by checking the dE_dvars array
		// The refE's are in order by one letter code, and a MultiVec is a vector1 of Reals. So E38A would touch indices (14 and 17)
		TS_ASSERT_DELTA( dE_dvars[ 14 ], -4.0566, 0.1 );
		TS_ASSERT_DELTA( dE_dvars[ 17 ],  4.0566, 0.1 );

		//score = ddg_bind_position_data->get_score( component_weights, vars, dE_dvars, free_score_list.size(),
		// chemical::num_canonical_aas, free_score_list.size() + chemical::num_canonical_aas, fixed_parameters, free_score_list, fixed_score_list );
		//TS_ASSERT_DELTA( score, 0.7064, 0.05 );


	}

};

