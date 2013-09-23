// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief


// libRosetta headers
// AUTO-REMOVED #include <devel/dna/protocols.hh>
// AUTO-REMOVED #include <devel/dna/ProteinDNA_Relax.hh>
#include <devel/dna/relax_util.hh>
#include <devel/dna/base_movers.hh>
#include <devel/cartesian_frags/DNA_FragLib.hh>
#include <utility/excn/Exceptions.hh>
// #include <devel/dna/util.hh>
// //#include <devel/dna/util.hh>
// AUTO-REMOVED #include <protocols/loops/ccd_closure.hh>
// #include <protocols/loops/loops_main.hh>
// #include <protocols/loops/Loops.hh>
// #include <protocols/frags/TorsionFragment.hh>

#include <protocols/viewer/viewers.hh>
// #include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/MonteCarlo.hh>
// #include <protocols/moves/rigid_body_moves.hh>

#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/dna/setup.hh>
// AUTO-REMOVED #include <core/scoring/dna/base_geometry.hh>
// AUTO-REMOVED #include <core/scoring/dna/BasePartner.hh>
// #include <core/scoring/dna/DNA_BasePotential.hh>
// #include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/scoring/constraints/Func.hh>
// AUTO-REMOVED #include <core/scoring/constraints/AtomPairConstraint.hh>
// AUTO-REMOVED #include <core/scoring/constraints/AngleConstraint.hh>
// AUTO-REMOVED #include <core/scoring/constraints/CoordinateConstraint.hh>
// AUTO-REMOVED #include <core/scoring/constraints/ConstraintSet.hh>
// #include <core/scoring/etable/Etable.hh>
// #include <core/scoring/ScoringManager.hh>
// #include <core/scoring/AtomVDW.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunctionFactory.hh>
// #include <core/scoring/rms_util.hh>
// #include <core/scoring/hbonds/hbonds.hh>
// #include <core/scoring/hbonds/HBondSet.hh>
// #include <core/scoring/elec/FA_ElecEnergy.hh>
// #include <core/scoring/etable/EtableEnergy.hh>
// #include <core/scoring/etable/count_pair/CountPairAll.hh>
// #include <core/scoring/etable/count_pair/CountPairFunction.hh>
// #include <core/scoring/etable/count_pair/CountPairFactory.hh>
// //#include <core/scoring/etable/count_pair/CountPair1BC4.hh>

// #include <core/types.hh>

// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.hh>
// #include <core/chemical/ResidueSelector.hh>
// AUTO-REMOVED #include <core/chemical/VariantType.hh>
// AUTO-REMOVED
// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>
// #include <core/chemical/AtomTypeSet.hh>
// #include <core/chemical/MMAtomTypeSet.hh>
// #include <core/chemical/AA.hh>

// #include <core/conformation/util.hh>
// AUTO-REMOVED #include <core/conformation/ResidueFactory.hh>
// AUTO-REMOVED #include <core/conformation/Residue.hh>
// #include <core/conformation/ResidueMatcher.hh>

// #include <core/pack/rotamer_trials.hh>
// #include <core/pack/pack_rotamers.hh>
// AUTO-REMOVED #include <core/pack/task/PackerTask.hh>
// AUTO-REMOVED #include <core/pack/task/TaskFactory.hh>
// #include <core/pack/rotamer_set/RotamerCouplings.hh>
// #include <core/pack/rotamer_set/WaterPackingInfo.hh>

// AUTO-REMOVED #include <core/kinematics/util.hh>
// AUTO-REMOVED #include <core/kinematics/visualize.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/id/AtomID_Map.hh>
// AUTO-REMOVED #include <core/id/AtomID_Map.Pose.hh>

// AUTO-REMOVED #include <core/optimization/AtomTreeMinimizer.hh>
// AUTO-REMOVED #include <core/optimization/MinimizerOptions.hh>

#include <core/pose/Pose.hh>

#include <basic/options/util.hh>
// //#include <basic/options/after_opts.hh>

#include <basic/prof.hh> // profiling
// #include <basic/basic.hh>
// #include <core/id/SequenceMapping.hh>

#include <devel/init.hh>

#include <core/io/pdb/pose_io.hh>

// #include <utility/vector1.hh>

// #include <numeric/xyzVector.hh>
// #include <numeric/xyzMatrix.hh>
// #include <numeric/xyz.functions.hh>

#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>

#include <ObjexxFCL/string.functions.hh>


// // C++ headers
// //#include <cstdlib>
#include <fstream>
// #include <iostream>
// #include <string>
// #include <set>
// #include <cstdlib>
// #include <sstream>

//silly using/typedef


#include <basic/Tracer.hh>

// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/dna.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <ObjexxFCL/format.hh>



using basic::T;
using basic::Error;
using basic::Warning;

static numeric::random::RandomGenerator RG(54323); // <- Magic number, do not change it!!!

using namespace core;
//using namespace protocols;

using utility::vector1;
using std::string;


static basic::Tracer tt( "demo.phil.dna_frag_test", basic::t_trace );
static basic::Tracer td( "demo.phil.dna_frag_test", basic::t_debug );
static basic::Tracer ti( "demo.phil.dna_frag_test", basic::t_info );
static basic::Tracer tw( "demo.phil.dna_frag_test", basic::t_warning );

///////////////////////////////////////////////////////////////////////////////
std::string
filebase( std::string const & file )
{
	size_t found = file.find_last_of("/\\");
	if ( found == std::string::npos ) return file;
	else return file.substr(found+1);
}


///////////////////////////////////////////////////////////////////////////////
void
dna_stats()
{
	using namespace pose;
	using namespace kinematics;
	using namespace id;
	using namespace conformation;
	using namespace chemical;
	using namespace scoring;
	using namespace scoring::dna;
	using namespace optimization;


	using namespace devel::dna;
	using namespace devel::cartesian_frags;

	//// retrieve commandline options /////////////////
	std::string output_tag, score_function_file;
	Size nstruct( 10 );
	Size ninner( 20 );
	Size nouter( 10 );
	vector1< string > fraglib_files, input_files;
	core::Real const min_tol( 0.001 );

	{ // scope
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace basic::options::OptionKeys::dna;
		output_tag = option[ out::output_tag ];
		score_function_file = option[ specificity::score_function ];
		fraglib_files = option[ specificity::frag_files ]();
		input_files = start_files();
		if ( option[ out::nstruct              ].user() ) nstruct = option[ out::nstruct ];
		if ( option[ specificity::n_inner ].user() ) ninner  = option[ specificity::n_inner ];
		if ( option[ specificity::n_outer ].user() ) nouter  = option[ specificity::n_outer ];
	}

	// randomize the order of the input files:
	numeric::random::random_permutation( input_files, RG );

	// the simulation pose
	Pose pose;
	protocols::viewer::add_conformation_viewer( pose.conformation(), "jump_pose" );

	// collect the fragments
	DNA_FragLibOP fraglib( new DNA_FragLib() );
	build_frag_libraries( fraglib_files, *fraglib );


	for ( Size nn=1; nn<= input_files.size(); ++nn ) {

		////// read starting pdb, setup the simulation pose
		std::string const & input_file( input_files[nn] );
		Pose pdb_pose;
		std::cout << "Reading input_file= " << input_file << std::endl;
		core::import_pose::pose_from_pdb( pdb_pose, input_file );
		set_base_partner( pdb_pose ); // fills base partner inf

		//pdb_pose.dump_pdb("test"+filebase( input_file ));

		// takes just the dna, and sets basepair foldtree. clears current data in pose
		setup_dna_only_jump_pose( pdb_pose, pose );

		// setup the atompair and angle constraints to penalize chainbreaks
		setup_dna_chainbreak_constraints( pose );

		std::cout << "dna_jump_pose for " << input_file << " nres= " << pose.total_residue() << std::endl;

		Pose const start_pose( pose );

		// scorefunction
		scoring::ScoreFunctionOP cst_scorefxn( new ScoreFunction() );
		cst_scorefxn->set_weight( atom_pair_constraint, 1.0 );
		cst_scorefxn->set_weight(     angle_constraint, 1.0 );

		scoring::ScoreFunctionOP scorefxn( new ScoreFunction() );

		//scorefxn.energy_method_options().exclude_DNA_DNA( false );
		// Safer:
		methods::EnergyMethodOptions options( scorefxn->energy_method_options() );
		options.exclude_DNA_DNA( false );
		scorefxn->set_energy_method_options( options );

		scorefxn->add_weights_from_file( score_function_file );
		scorefxn->set_weight( atom_pair_constraint, 1.0 );
		scorefxn->set_weight(     angle_constraint, 1.0 );

		Real const start_score( (*scorefxn)( pose ) );
		ti << "start_scores:\n";
		scorefxn->show( ti, pose );
		ti << '\n'; // flush??

		for ( Size n=1; n<= nstruct; ++n ) {
			using namespace protocols;
			using namespace protocols::moves;

			pose.clear();
			pose = start_pose;

			// montecarlo
			MonteCarloOP mc( new MonteCarlo( pose, *scorefxn, 0.8 ) );

			// movemap
			MoveMapOP mm( new MoveMap() );
			mm->set_jump( true );
			mm->set_chi( true );
			mm->set_bb( true );

			// movers
			BasePairMoverOP bp_mover( new BasePairMover( fraglib, 2.0, 20, 2.0, cst_scorefxn ) );
			BaseStepMoverOP bs_mover( new BaseStepMover( fraglib, 2.0, 20, 2.0, cst_scorefxn ) );
			protocols::simple_moves::MinMoverOP minmover( new protocols::simple_moves::MinMover( mm, scorefxn, "dfpmin", min_tol, true ) );
			TrialMoverOP bp_trial, bs_trial;
			{
				SequenceMoverOP bpseq( new SequenceMover() ), bsseq( new SequenceMover() );
				bpseq->add_mover( bp_mover  );
				bpseq->add_mover( minmover );
				bp_trial = new TrialMover( bpseq, mc );
				bsseq->add_mover( bs_mover  );
				bsseq->add_mover( minmover );
				bs_trial = new TrialMover( bsseq, mc );
			}

			// initial min-trial
			minmover->apply( pose );
			core::Real const aftermin_rms  ( scoring::all_atom_rmsd( pose, start_pose ) );
			core::Real const aftermin_score( (*scorefxn)( pose ) );
			ti << "aftermin_scores:\n";
			scorefxn->show( ti, pose );
			ti << '\n';
			mc->boltzmann( pose );

			basic::prof_reset();
			for ( Size outer=1; outer<= nouter; ++outer ) {
				for ( Size inner=1; inner<= ninner; ++inner ) {
					bp_trial->apply( pose );
					bs_trial->apply( pose );
				}
				mc->recover_low( pose );

				// flip the foldtree
				core::Real const score_before( (*scorefxn)( pose ) );
				setup_dna_only_fold_tree( pose, outer%2 == 1 );
				core::Real const score_after( (*scorefxn)( pose ) );
				assert( std::abs( score_before - score_after ) < 1e-1 );
				mc->reset( pose );

				basic::prof_show();
			}

			// write the pdb
			string const outfilename( output_tag + filebase( input_file ) + "_final" + lead_zero_string_of( n, 4 )+".pdb" );
			std::ofstream out( outfilename.c_str() );
			pose.dump_pdb( out );
			scorefxn->show( out, pose );
			out.close();
			basic::prof_show();

			// status output
			mc->show_counters();
			ti << "final_scores:\n";
			scorefxn->show( ti, pose );
			ti << '\n';
			core::Real const final_rms( scoring::all_atom_rmsd( pose, start_pose ) );
			core::Real const final_score( (*scorefxn)(pose) );
			{
				using namespace ObjexxFCL::format;
				std::cout << "final_rmsd_and_score: " << outfilename << F(9,3,aftermin_rms) << F(9,3,final_rms) <<
					F(9,3,start_score) << F(9,3,aftermin_score) << F(9,3,final_score) << std::endl;
			}

		} // nstruct loop
	}
}

////////////////////////////////////////////////////////////////////
void*
my_main( void* )
{
	dna_stats();
	exit(0);
}

///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try{
	// initialize option and random number system
	devel::init( argc, argv );

	protocols::viewer::viewer_main( my_main );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl; 
	} 
}
