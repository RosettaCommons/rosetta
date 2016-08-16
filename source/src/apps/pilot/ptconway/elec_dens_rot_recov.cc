// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief


// libRosetta headers
#include <core/types.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <utility/excn/Exceptions.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <core/scoring/electron_density/util.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/Energies.hh>
#include <core/conformation/Residue.hh>

#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>

#include <basic/basic.hh>
#include <basic/database/open.hh>
#include <devel/init.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/electron_density/SetupForDensityScoringMover.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/simple_moves/GenericMonteCarloMover.hh>
#include <protocols/simple_moves/MissingDensityToJumpMover.hh>
#include <protocols/viewer/viewers.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/filters/BasicFilters.hh>


#include <core/io/pdb/pdb_writer.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>
#include <protocols/simple_moves/ConstraintSetMover.hh>

#include <core/scoring/constraints/util.hh>


#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>

// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

#ifdef USEMPI
#include <mpi.h>
#endif


//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <basic/Tracer.hh>
using basic::T;
using basic::Error;
using basic::Warning;
using namespace core;
using namespace protocols;
using utility::vector1;



OPT_1GRP_KEY(Boolean, edrr, dump)
OPT_1GRP_KEY(Boolean, edrr, dump_asymm)
OPT_1GRP_KEY(Boolean, edrr, debug)
OPT_1GRP_KEY(Boolean, edrr, print_all)
OPT_1GRP_KEY(Real, edrr, chi_diff_cutoff)
OPT_1GRP_KEY(Real, edrr, auto_rmsd_cutoff)
OPT_1GRP_KEY(Real, edrr, dens_diff_cutoff)
OPT_1GRP_KEY(Integer, edrr, gmc_trials)
OPT_1GRP_KEY(String, edrr, prefix)

// return true if line1 < line2
bool output_cmp( std::string line1, std::string line2 )
{
	Real key1, key2;
/*
	std::istringstream l1(line1);
	std::istringstream l2(line2);

	// sort column = 3 for dens_diff
	for ( int col = 1; col <= 3; ++col ) {
		l1 >> key1;
		l2 >> key2;
	}	*/

	key1 = atof(line1.substr(18, 14).c_str());
	key2 = atof(line2.substr(18, 14).c_str());

	//std::cout << key1 << " " << key2 << std::endl;
	return (key1 < key2);
}

// class definitions

class ElecDensMinPackMinReporter {
private:
	vector1<std::string> output;
	Real dens_cut, auto_cut, chi_cut;

public:
	ElecDensMinPackMinReporter(){}

	void find_flips( core::pose::Pose & first_pass, core::pose::Pose & second_pass, std::string output_id, int nres ) 		// use nres to ignore extra residues in symmetric pose
	{
		using namespace std;
		using namespace basic::options;

		float dens_diff_cutoff = option[ OptionKeys::edrr::dens_diff_cutoff ];			// above: high second pass - low first pass = positive
		float auto_rmsd_cutoff = option[ OptionKeys::edrr::auto_rmsd_cutoff ];			// above
		float chi_diff_cutoff = option[ OptionKeys::edrr::chi_diff_cutoff ];			// above

		float dens_diff;
		float auto_rmsd;		//other rms as well?
		float chi1_diff;
		float chi2_diff;
		float chi3_diff;
		float chi4_diff;

		for (int ii = 1; ii <= nres; ++ii)
		{
			stringstream output_line;
			dens_diff = second_pass.energies().residue_total_energies(ii)[ core::scoring::elec_dens_fast ] - first_pass.energies().residue_total_energies(ii)[ core::scoring::elec_dens_fast ];
			auto_rmsd = core::scoring::automorphic_rmsd( first_pass.residue( ii ), second_pass.residue( ii ), false );		//dont superimpose - fixed bb is sufficient

			switch(first_pass.residue(ii).nchi())
			{
				case 4:	chi4_diff = second_pass.residue(ii).chi()[4] - first_pass.residue(ii).chi()[4];
				case 3:	chi3_diff = second_pass.residue(ii).chi()[3] - first_pass.residue(ii).chi()[3];
				case 2:	chi2_diff = second_pass.residue(ii).chi()[2] - first_pass.residue(ii).chi()[2];
				case 1:	chi1_diff = second_pass.residue(ii).chi()[1] - first_pass.residue(ii).chi()[1];
				case 0:	break;
				default:	cout << "wtf? this shouldn't happen..." << endl;
			}

			if ( ! option[ OptionKeys::edrr::print_all ]() ) {
				if ( dens_diff < dens_diff_cutoff )	continue;
				if ( auto_rmsd < auto_rmsd_cutoff )	continue;
//			if ( chi1_diff[ res ] < chi_diff_cutoff && chi2_diff[ res ] < chi_diff_cutoff && chi3_diff[ res ] < chi_diff_cutoff && chi4_diff[ res ] < chi_diff_cutoff )	continue;
			}

			// print intersections
			// save to string and add string to list for later sorting
			output_line << fixed << setprecision(2) << setw(4) << output_id << setw(8) << first_pass.residue( ii ).name3() << setw(6) << ii << setw(14) <<
								setw(14) << dens_diff << setw(14) << auto_rmsd <<
								setw(14) << chi1_diff << setw(14) << chi2_diff <<
								setw(14) << chi3_diff << setw(14) << chi4_diff;// << endl;
			output.push_back( output_line.str() );
		}
	}

	void print_flips( std::ostream & out )
	{
		using namespace std;

		/*#ifdef USEMPI
		int myrank;
		MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
		if (myrank == 0)
			cout << "master" << endl;
			int num_nodes;
			MPI_Comm_size(MPI_COMM_WORLD, &num_nodes);
			MPI_Gather( NULL, NULL, NULL, output, 1, MPI_STRING, 0, MPI_COMM_WORLD );
			MPI_Barrier(MPI_COMM_WORLD);
		else
			cout << "slave" << endl;
			//send lists here
		#endif*/

	  std::sort( output.begin(), output.end(), output_cmp ); 		// output_cmp defined outside of class

		out << setw(4) << "id" << setw(8) << "type" << setw(6) << "res" << setw(14) << "dens_diff" << setw(14) << "auto_rmsd" << setw(14) <<
        "chi1_diff" << setw(14) << "chi2_diff" << setw(14) << "chi3_diff" << setw(14) << "chi4_diff" << endl;
    out << "-----------------------------------------------------------------------------------------------------" << endl;

		for ( vector1<string>::iterator ii = output.begin(); ii != output.end(); ++ii )
			out << *(ii) << endl;

	}

};

// global reporter for jd2 access
namespace globals {
	ElecDensMinPackMinReporter reporter;
}

class ElecDensMinPackMinMover : public protocols::moves::Mover {
public:
	ElecDensMinPackMinMover(){}
	void apply( core::pose::Pose & pose) {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		using namespace protocols::moves;
		using namespace scoring;

		// steal relax flags
		scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();

		kinematics::MoveMapOP mm = new kinematics::MoveMap();
		mm->set_bb  ( false );
		mm->set_chi ( true );
		mm->set_jump( true );
		mm->set( core::id::THETA, option[ OptionKeys::relax::minimize_bond_angles ]() );
		mm->set( core::id::D, option[ OptionKeys::relax::minimize_bond_lengths ]() );

		int nres = pose.total_residue();	//only consider the residues in the original pose, not symmetric copies

		if ( option[ OptionKeys::symmetry::symmetry_definition ].user() )  {
			protocols::simple_moves::symmetry::SetupForSymmetryMoverOP symm( new protocols::simple_moves::symmetry::SetupForSymmetryMover );
			symm->apply( pose );
			core::pose::symmetry::make_symmetric_movemap( pose, *mm );
		}

		// now add density scores from cmd line
		if ( option[ edensity::mapfile ].user() ) {
			core::scoring::electron_density::add_dens_scores_from_cmdline_to_scorefxn( *scorefxn );
		}

		// set pose for density scoring if a map was input
		//   + (potentially) dock map into density
		if ( option[ edensity::mapfile ].user() ) {
			protocols::electron_density::SetupForDensityScoringMoverOP edens
											 ( new protocols::electron_density::SetupForDensityScoringMover );
			edens->apply( pose );
		}

		// setup the minimizer mover
		protocols::simple_moves::MinMoverOP minimizer = new protocols::simple_moves::symmetry::SymMinMover( mm, scorefxn, option[ OptionKeys::relax::min_type ]() , 0.0001, true, false, false );
		minimizer->cartesian( true );

		// setup the packer mover
		core::pack::task::TaskFactoryOP task_factory( new core::pack::task::TaskFactory );
		task_factory->push_back( new core::pack::task::operation::InitializeFromCommandline );
		task_factory->push_back( new core::pack::task::operation::RestrictToRepacking );
		task_factory->push_back( new core::pack::task::operation::IncludeCurrent );
		core::pack::task::PackerTaskOP packer_task( task_factory->create_task_and_apply_taskoperations( pose ) );
		protocols::moves::MoverOP packer = new protocols::simple_moves::symmetry::SymPackRotamersMover( scorefxn, packer_task );					//extend this to handle asymmetry, etc

		// setup the sequence mover
		protocols::moves::SequenceMoverOP seqmov = new protocols::moves::SequenceMover;
    seqmov->add_mover( new protocols::simple_moves::MissingDensityToJumpMover );
		seqmov->add_mover( minimizer );
		seqmov->add_mover( packer );
		seqmov->add_mover( minimizer );

		// setup the GenericMonteCarlo mover
		protocols::simple_moves::GenericMonteCarloMoverOP gmc = new protocols::simple_moves::GenericMonteCarloMover();
		protocols::filters::FilterOP falseFilter = new protocols::filters::FalseFilter;
		gmc->stopping_condition( falseFilter );
		gmc->set_drift( false );
		gmc->set_preapply( false );
		gmc->set_scorefxn( scorefxn );
		gmc->set_mover( seqmov );
		gmc->set_maxtrials( option[ OptionKeys::edrr::gmc_trials ]() );

		// testing against Rosetta scripts output
		/*Pose test_pose_hi = pose;
		//scorefxn->set_weight( core::scoring::elec_dens_fast, 100 );
		minimizer->score_function()->show(std::cout);
		(*scorefxn)(test_pose_hi);
    scorefxn->show(std::cout, test_pose_hi);
		//minimizer->score_function( scorefxn );
		minimizer->apply( test_pose_hi );
		//scoring::ScoreFunction scorefxn_hi = *(minimizer->score_function());
		//scorefxn_hi.show(std::cout, test_pose_hi);
		(*scorefxn)(test_pose_hi);
    scorefxn->show(std::cout, test_pose_hi);
		test_pose_hi.dump_pdb("test_pose_hi.pdb");
		std::cout << minimizer;
		std::cout << "dumped test_pose_hi.pdb" << std::endl;

		Pose test_pose_lo = pose;
		scorefxn->set_weight( core::scoring::elec_dens_fast, 30 );
		minimizer->apply( test_pose_lo );
		test_pose_lo.dump_pdb("test_pose_lo.pdb");
		std::cout << "dumped test_pose_lo.pdb" << std::endl;*/


		// run first pass
		Pose first_pass = pose;
		scorefxn->set_weight( core::scoring::elec_dens_fast, 100 );
		gmc->apply( first_pass );

		// run second pass
		Pose second_pass = first_pass;
		scorefxn->set_weight( core::scoring::elec_dens_fast, 30 );
		gmc->apply( second_pass );
			/*Pose recov_low;
			gmc->recover_low( recov_low );
			std::cout << "last pose vs recover_low all atom rms " << core::scoring::all_atom_rmsd( second_pass, recov_low ) << std::endl;*/


		// output results
		std::string cur_output_name = protocols::jd2::JobDistributor::get_instance()->current_output_name();
		cur_output_name = cur_output_name.substr(0,4);

		if ( option[ OptionKeys::edrr::dump ]() ) {
			first_pass.dump_pdb( cur_output_name + ".high.pdb" );
			second_pass.dump_pdb( cur_output_name + ".low.pdb" );
		}

		if ( option[ OptionKeys::edrr::dump_asymm ]() ) {
			Pose asymm_first, asymm_second;
			pose::symmetry::extract_asymmetric_unit( first_pass, asymm_first, true );
			pose::symmetry::extract_asymmetric_unit( second_pass, asymm_second, true );
			asymm_first.dump_pdb( cur_output_name + ".asymm.high.pdb" );
			asymm_second.dump_pdb( cur_output_name + ".asymm.low.pdb" );
		}
		globals::reporter.find_flips( first_pass, second_pass, cur_output_name, nres );

		/*(if ( option[ OptionKeys::edrr::print_scores ]() )
    {
						(*scorefxn)(pose);
						scorefxn->show(cout, pose);

						(*scorefxn)(first_pass);
						scorefxn->show(cout, first_pass);

						(*scorefxn)(second_pass);
						scorefxn->show(cout, second_pass);
		}*/


	}

	virtual std::string get_name() const {
		return "ElecDensMinPackMinMover";
	}
};

void*
my_main( void* ) {
	using namespace protocols::moves;
	using namespace protocols::simple_moves::symmetry;

	try{
		protocols::jd2::JobDistributor::get_instance()->go( new ElecDensMinPackMinMover() );
	} catch ( utility::excn::EXCN_Base& excn ) {
		std::cerr << "Exception: " << std::endl;
		excn.show( std::cerr );
	}

	return 0;
}


///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
	try{
	NEW_OPT(edrr::dump, "dump high.pdb and low.pdb", false);
	NEW_OPT(edrr::dump_asymm, "dump asymm high.pdb and low.pdb", false);
	NEW_OPT(edrr::debug, "debug - trials=1", false);
	NEW_OPT(edrr::print_all, "ignore cutoffs", false);
	NEW_OPT(edrr::dens_diff_cutoff, "electron density score difference cutoff to qualify as rotamer flip", 1);
	NEW_OPT(edrr::auto_rmsd_cutoff, "automorphic rmsd cutoff to qualify as rotamer flip", 0.2);
	NEW_OPT(edrr::chi_diff_cutoff, "chi difference cutoff to qualify as a rotamer flip", 5);
	NEW_OPT(edrr::gmc_trials, "how many times to repeat repacking", 10);
	NEW_OPT(edrr::prefix, "prefix to use before .flips.txt in output name", "rotamer");

	devel::init(argc, argv);

	protocols::viewer::viewer_main( my_main );

	std::string outfile = basic::options::option[ basic::options::OptionKeys::edrr::prefix ]() + ".flips.txt";
	std::ofstream out( outfile.c_str() );
	globals::reporter.print_flips( out );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
