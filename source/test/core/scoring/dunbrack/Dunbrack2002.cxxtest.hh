// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/Dunbrack.cxxtest.hh
/// @brief  Dunbrack Unit Test
/// @author Oliver Lange

// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers

// Project Headers
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

//#include <core/chemical/ResidueTypeSet.hh>
// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>

#include <core/conformation/Residue.hh>

#include <core/graph/Graph.hh>

//amw
#include <core/pose/PDBInfo.hh>

// AUTO-REMOVED #include <core/pose/Pose.hh>

#include <core/scoring/ScoreFunction.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunctionFactory.hh>
// AUTO-REMOVED #include <core/scoring/ScoringManager.hh>
#include <core/pack/dunbrack/DunbrackRotamer.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>

#include <core/types.hh>

// Package headers
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
#include <utility/io/izstream.hh>

// AUTO-REMOVED #include <basic/Tracer.hh>
#include <test/UTracer.hh>

#include <ObjexxFCL/string.functions.hh>


// Project headers
#include <core/kinematics/MoveMap.hh>
//#include <core/optimization/MinimizerMap.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/NumericalDerivCheckResult.hh>

//Auto Headers
#include <utility/vector1.hh>


using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("core.scoring.Dunbrack.cxxtest");

// using declarations
using namespace core;
using namespace scoring;
using namespace conformation;
using namespace chemical;
using namespace scoring;
using namespace pose;
using namespace ObjexxFCL;

///////////////////////////////////////////////////////////////////////////
/// @name Dunbrack2002Test
/// @brief: unified tests for difference score functions/methods
///////////////////////////////////////////////////////////////////////////
class Dunbrack2002Test : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init_with_additional_options( "-restore_pre_talaris_2013_behavior -override_rsd_type_limit" );
	}

	void tearDown() {}

	void test_dunbrack_rotamers()
	{
		using namespace core;
		using namespace scoring;
		using namespace pack::dunbrack;

		pose::Pose pose( create_test_in_pdb_pose() );
		//core::import_pose::pose_from_pdb( pose, "core/scoring/test_in.pdb" );

		pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
		// just use standard task no command-line parsing to make test robust against commandline changes
		// bla

		utility::io::izstream data( "core/scoring/dunbrack/test_in_rotamers.dat" );
		std::string line;
		TS_ASSERT( data.good() );
		utility::vector1< utility::vector1< Real > > chi_gold;
		utility::vector1< Real > rot_ener_gold, pose_ener_gold;
		while ( getline( data, line ) ) {
			std::istringstream in( line );
			Size pos, nchi;
			std::string temp1, temp2;
			in >> temp1 >> temp2;

			assert( is_int(temp1) && is_int(temp2) );
			pos  = int_of(temp1);
			nchi = int_of(temp2);
            
            TS_ASSERT( pose.residue( pos ).nchi() == nchi );
			utility::vector1< Real > chis;
			chi_gold.push_back( chis );
			for ( Size ii = 1; ii <= nchi; ++ii ) {
				Real x;
				in >> x;
				chi_gold[ chi_gold.size() ].push_back( x );
			}
			Real rot_ener;
			Real pose_ener;
			in >> rot_ener >> pose_ener;
			rot_ener_gold.push_back( rot_ener );
			while ( pose_ener_gold.size() < pos ) pose_ener_gold.push_back( -1.0 );
			pose_ener_gold[ pos ] = pose_ener;
		}
		graph::GraphOP dummy_graph( new graph::Graph() );
		scoring::ScoreFunction dummy_scorefxn;
		Size ct( 1 );
		for (Size pos = 1; pos <= pose.total_residue(); pos++ ) {
			Residue const & residue( pose.residue( pos ) );

			utility::vector1< ResidueOP > suggested_rotamers;

			// generate empty list of extra_chi_steps
			utility::vector1< utility::vector1< Real > > extra_chi_steps( residue.nchi() );

			SingleResidueRotamerLibraryCOP rotlib = RotamerLibrary::get_instance()->get_rsd_library( residue.type() );
			if (rotlib) {
				rotlib->fill_rotamer_vector( pose, dummy_scorefxn, *task, dummy_graph, residue.type().get_self_ptr(), residue, extra_chi_steps, false /*buried*/, suggested_rotamers);
			}
            
            bool bOut ( false  );//switch to true to produce a new test_input file
			for ( utility::vector1< ResidueOP >::const_iterator it = suggested_rotamers.begin(),
							eit = suggested_rotamers.end();
						it!=eit;
						++it ) {

				// if the number of rotamers built does not match the number in the "gold" files, quit
				if ( ct > chi_gold.size() ) {
                    TS_ASSERT( ct <= chi_gold.size() ); // tell the user that something is wrong
					break;
				}

				if ( bOut )	std::cout << pos << " " << residue.nchi() << " ";
				else {
					if ( chi_gold[ct].size() != residue.nchi() ) {
                        
						TS_ASSERT( chi_gold[ ct ].size() == residue.nchi() );
						++ct;
						break; // something is very wrong
					}
				}

				for ( Size n = 1; n <= residue.nchi(); ++n ) {
					if ( bOut ) {
						std::cout << (*it)->chi()[ n ] << " " ;
					}
					else {
						TS_ASSERT_DELTA( chi_gold[ ct ][ n ],  (*it)->chi()[ n ] , 0.001);
					}
				}
				RotamerLibraryScratchSpace scratch;
				if ( bOut ) std::cout << rotlib->rotamer_energy( (**it), scratch ) << " " << rotlib->rotamer_energy( residue, scratch);
				else {
					TS_ASSERT_DELTA( rot_ener_gold[ ct ], rotlib->rotamer_energy( **it, scratch ), 0.001);
					TS_ASSERT_DELTA( pose_ener_gold[ pos ], rotlib->rotamer_energy( residue, scratch ), 0.001);
				}
				if ( residue.nchi() ) {
					ct++;
				}
				if ( bOut ) std::cout << std::endl;
			} // read-out
		} // residues

	}// test

	void test_dunbrack_best_rotamer_energy()
	{
		using namespace core;
		using namespace scoring;
		using namespace pack::dunbrack;
		test::UTracer UT("core/scoring/dunbrack/best_rotamer_energy.u"); // the name "UT" matters -- it goes with the UTRACE macro

		pose::Pose pose(create_test_in_pdb_pose() );
		//core::import_pose::pose_from_pdb( pose, "core/scoring/test_in.pdb" );

		for (Size pos = 1; pos <= pose.total_residue(); pos++ ) {
			Residue const & residue( pose.residue( pos ) );
			SingleResidueRotamerLibraryCOP rotlib = RotamerLibrary::get_instance()->get_rsd_library( residue.type() );
			if( rotlib.get() == NULL ) continue;

			RotamerLibraryScratchSpace scratch;
			Real const this_rotamerE = rotlib->best_rotamer_energy(residue, true /*current well only*/, scratch);
			Real const best_rotamerE = rotlib->best_rotamer_energy(residue, false /*global best*/, scratch);
			UTRACE << "this (ideal) = " << this_rotamerE << "; best (this phi,psi) = " << best_rotamerE << std::endl;
			TS_ASSERT( best_rotamerE <= this_rotamerE );
		}
	}

	void test_dunbrack_get_rotamer_probability()
	{
		using namespace core;
		using namespace core::chemical;
		using namespace core::scoring;
		using namespace core::pack::dunbrack;

		//mjo commenting out 'sc_manager' because it is unused and causes a warning
		//ScoringManager * sc_manager = ScoringManager::get_instance();
		RotamerLibrary & rotlib = *RotamerLibrary::get_instance();

		Real const phi_example = -59;
		Real const psi_example =  61;
        
		for ( Size ii = 1; ii <= num_canonical_aas; ++ii ) {

			if ( AA( ii ) == aa_ala || AA( ii ) == aa_gly ) continue;
			SingleResidueRotamerLibraryCOP aa_rotlib = rotlib.get_library_by_aa( (AA) ii );
			SingleResidueDunbrackLibraryCOP aa_dunlib( utility::pointer::dynamic_pointer_cast< core::pack::dunbrack::SingleResidueDunbrackLibrary const > ( aa_rotlib ) );
			TS_ASSERT( aa_dunlib );
			if ( ! aa_dunlib ) {
				std::cerr << "Failed to find dunbrack library for aa " << (AA) ii << std::endl;
				continue;
			}

			Size const ii_nrots = aa_dunlib->n_rotamer_bins();
            utility::vector1< DunbrackRotamerSampleData > aa_samples = aa_dunlib->get_all_rotamer_samples( phi_example, psi_example );

			TS_ASSERT( aa_samples.size() <= ii_nrots );
			for ( Size jj = 1; jj <= aa_samples.size(); ++jj ) {
                Real const jj_prob = aa_dunlib->get_probability_for_rotamer( phi_example, psi_example, jj );
				TS_ASSERT_DELTA( aa_samples[ jj ].probability(), jj_prob, 1e-10 );
			}
		}
	}

	void test_chi_derivatives_w_dunE()
	{
		using namespace core;
		using namespace core::id;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::constraints;
		using namespace core::optimization;

		Pose pose = create_trpcage_ideal_pose();
		ScoreFunction sfxn;
		sfxn.set_weight( fa_dun, 0.75 );


		AtomTreeMinimizer minimizer;
		//std::cout.precision( 16 );
		//std::cout << "start score: " << sfxn(pose) << std::endl;
		//Real start_score =
		sfxn(pose);
		//TS_ASSERT_DELTA( 38.86927045441701, start_score, 1e-12 );


		kinematics::MoveMap movemap;
		movemap.set_chi( true ); // apparently, there is a bug in the neighborlist-autoupdate code and I trip that bug in this test case

		MinimizerOptions min_options( "dfpmin_armijo", 0.01, true, true, false );
		min_options.deriv_check_to_stdout( false );
		//min_options.nblist_auto_update( true );

		minimizer.run( pose, movemap, sfxn, min_options );
		NumericalDerivCheckResultOP deriv_check_result = minimizer.deriv_check_result();
		for ( Size ii = 1, iiend = deriv_check_result->n_deriv_check_results(); ii <= iiend; ++ii ) {
			NumDerivCheckData const & iidata( deriv_check_result->deriv_check_result( ii ) );
			TS_ASSERT( iidata.nsteps() >= 1 );
			for ( Size jj = 1; jj <= iidata.nangles(); ++jj ) {
				TS_ASSERT_DELTA( iidata.dof_step_data( jj, 1 ).num_deriv(), iidata.dof_step_data( jj, 1 ).ana_deriv(), 1e-6 );
			}
		}

		//Real end_score = sfxn(pose);
		//std::cout << "end score: " << sfxn(pose) << std::endl;
		//TS_ASSERT_DELTA( 38.57005418761457, end_score, 1e-12 );
		//pose.dump_pdb( "cstetest1.pdb" );

	}

	void test_phipsi_derivatives_w_dunE()
	{
		using namespace core;
		using namespace core::id;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::constraints;
		using namespace core::optimization;

		Pose pose = create_trpcage_ideal_pose();
		ScoreFunction sfxn;
		sfxn.set_weight( fa_dun, 0.75 );


		AtomTreeMinimizer minimizer;
		std::cout.precision( 16 );
		std::cout << "start score: " << sfxn(pose) << std::endl;
		//Real start_score =
		sfxn(pose);
		//TS_ASSERT_DELTA( 38.86927045441701, start_score, 1e-12 );


		kinematics::MoveMap movemap;
		movemap.set_bb( true );

		MinimizerOptions min_options( "dfpmin_armijo", 0.01, true, true, false );
		min_options.deriv_check_to_stdout( false );
		//min_options.nblist_auto_update( true );

		minimizer.run( pose, movemap, sfxn, min_options );
		NumericalDerivCheckResultOP deriv_check_result = minimizer.deriv_check_result();
		for ( Size ii = 1, iiend = deriv_check_result->n_deriv_check_results(); ii <= iiend; ++ii ) {
			NumDerivCheckData const & iidata( deriv_check_result->deriv_check_result( ii ) );
			TS_ASSERT( iidata.nsteps() >= 1 );
			for ( Size jj = 1; jj <= iidata.nangles(); ++jj ) {
                if ( jj % 3 == 1 ) // new res
                    std::cout << "Now looking at residue " << ( ( jj+2 )/3 ) << " which is a " << pose.residue_type( (jj+2)/3 ).name() << std::endl;
                else
                    std::cout << "Angle jj " << jj << " or for this residue, specifically " << (jj%3 ) << std::endl;
                
				TS_ASSERT_DELTA( iidata.dof_step_data( jj, 1 ).num_deriv(), iidata.dof_step_data( jj, 1 ).ana_deriv(), 1e-6 );
			}
		}

		//Real end_score = sfxn(pose);
		std::cout << "end score: " << sfxn(pose) << std::endl;
		//TS_ASSERT_DELTA( 38.57005418761457, end_score, 1e-12 );
		pose.dump_pdb( "dunmin1.pdb" );

	}

	void dont_test_atomtree_minimization_with_etable_and_dunE()
	{
		using namespace core;
		using namespace core::id;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::constraints;
		using namespace core::optimization;

		Pose pose = create_trpcage_ideal_pose();
		ScoreFunction sfxn;
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.5 );
		sfxn.set_weight( fa_dun, 0.75 );


		AtomTreeMinimizer minimizer;
		std::cout.precision( 16 );
		std::cout << "start score: " << sfxn(pose) << std::endl;
		Real start_score = sfxn(pose);
		TS_ASSERT_DELTA( 38.86927045441701, start_score, 1e-12 );


		kinematics::MoveMap movemap;
		movemap.set_bb( true );
		//movemap.set_chi( true ); // apparently, there is a bug in the neighborlist-autoupdate code and I trip that bug in this test case

		MinimizerOptions min_options( "dfpmin_armijo", 0.01, true, false, false );
		//min_options.nblist_auto_update( true );

		minimizer.run( pose, movemap, sfxn, min_options );

		Real end_score = sfxn(pose);
		std::cout << "end score: " << end_score << std::endl;
		//TS_ASSERT_DELTA( 38.57005418761457, end_score, 1e-12 );
		pose.dump_pdb( "dunmintest1.pdb" );

	}

};
