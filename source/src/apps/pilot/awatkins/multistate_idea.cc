// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


//awatkins: based heavily on code from kdrew/oop_dock_design.cc

// Project Headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Conformation.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/ScoringManager.hh>

#include <core/chemical/VariantType.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>

#include <utility/pointer/owning_ptr.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

// Mover headers
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/PyMOLMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/simple_moves/TaskAwareMinMover.hh>
#include <protocols/simple_moves/BackboneMover.fwd.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/hbs/HbsPatcher.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/rigid/RB_geometry.hh>

// Filter headers
#include <basic/MetricValue.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
//#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
#include <core/pose/metrics/simple_calculators/SasaCalculatorLegacy.hh>

#include <protocols/toolbox/pose_metric_calculators/NumberHBondsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/PackstatCalculator.hh>

// Utility Headers
#include <devel/init.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>

// C++ headers
#include <string>
#include <sstream>
#include <fstream>


//The original author used a lot of using declarations here.  This is a stylistic choice.
// Namespaces
using namespace core;
using namespace conformation;
using namespace chemical;
using namespace scoring;
using namespace pose;
using namespace protocols;
using namespace protocols::moves;
using namespace core::pack::task;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::id;
using basic::Error;
using basic::Warning;
using utility::file::FileName;

// tracer - used to replace cout
static basic::Tracer TR( "MESM" );

// application specific options
namespace mesm {
// pert options
RealOptionKey const boltz_kT( "mesm::boltz_kT" );
IntegerOptionKey const samples_per_input( "mesm::samples_per_input" );
}

int
main( int argc, char* argv[] )
{
	try {

		option.add( mesm::boltz_kT, "The kT to use for generating Boltzmann factors." ).def( 0.593 );
		option.add( mesm::samples_per_input, "The number of states to sample per input structure" ).def( 1 );

		// init command line options
		//you MUST HAVE THIS CALL near the top of your main function, or your code will crash when you first access the command line options
		devel::init(argc, argv);

		// get list of input files
		utility::vector1< std::string > infiles;
		if ( basic::options::option[in::file::l].user() ) {
			infiles = basic::options::option[in::file::l]();
		}

		core::Size nstruct = basic::options::option[ mesm::samples_per_input ]();
		core::Real kT = basic::options::option[ mesm::boltz_kT ]();
		core::Size nres = 0;

		scoring::ScoreFunctionOP score_fxn = get_score_function();
		scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn(*score_fxn);

		utility::vector1< core::pose::Pose > poses;
		utility::vector1< core::Real > scores;

		protocols::relax::FastRelaxOP fr( new protocols::relax::FastRelax( score_fxn, 5 ) );

		// now we have a vector of pdbs in infiles
		for ( Size i = 1; i <= infiles.size(); ++i ) {
			// Do something to each file
			core::pose::Pose pose;
			core::import_pose::pose_from_file( pose, infiles[ i ] , core::import_pose::PDB_file);

			if ( nres == 0 ) {
				nres = pose.pdb_info()->nres();
			}

			if ( false ) {
				// just score
				poses.push_back( pose );
				scores.push_back( ( *score_fxn ) ( pose ) );
			} else {
				core::pose::Pose stored_pose = pose;
				for ( Size j = 1; j <= nstruct; ++j ) {
					// 1 relax trajectory per nstruct
					//protocols::jd2::JobDistributor::get_instance()->go( fr );
					fr->apply( pose );
					poses.push_back( pose );
					scores.push_back( ( *score_fxn ) ( pose ) );
					pose = stored_pose;
				}
			}
		}

		// now we have our pose, score ensemble
		// really we could just score the poses here of course, and we'd not have to do it mid-loop, but w/e
		std::ofstream statefile( "statefile" );
		std::ofstream fitnessfile( "fitnessfile" );

		utility::vector1< core::Real > boltz;
		TR << "Scores: " << std::endl;
		for ( Size i = 1; i <= infiles.size() * nstruct; ++i ) {
			statefile << "pose_000" << i << ".pdb correspondence_" << i << " " << i << "_design.res" << std::endl;
			fitnessfile << "STATE_VECTOR sv statefile" << std::endl;

			TR << scores[ i ] << std::endl;
			std::stringstream fn;
			fn << "pose_000" << i << ".pdb";
			poses[ i ].dump_pdb( fn.str().c_str() );
			boltz.push_back( std::exp( -kT * scores[ i ] ) );

			std::stringstream corr;
			corr << "correspondence_" << i;
			std::ofstream correspondence( corr.str().c_str() );
			for ( Size resi = 1; resi <= nres; ++resi ) {
				//std::string pdbnum = poses[i].pdb_info()->pose2pdb( resi );
				if ( poses[i].pdb_info()->chain( resi ) == 'B' ) { // only designing chain B
					correspondence << resi << " " << poses[i].pdb_info()->chain( resi ) << " " << poses[i].pdb_info()->number( resi ) << std::endl;
				}
			}

			std::stringstream res;
			res << i << "_design.res";
			std::ofstream resfile( res.str().c_str() );
			resfile << "NATRO" << std::endl << "start" << std::endl;
			for ( Size resi = 1; resi <= nres; ++resi ) {
				std::string pdbnum = poses[i].pdb_info()->pose2pdb( resi );
				//TR << pdbnum << std::endl;
				if ( poses[i].pdb_info()->chain( resi ) == 'A' ) { // repacking chain A
					resfile << pdbnum << " NATAA EX 1 EX 2" << std::endl;
				}
			}
		}
		core::Real boltz_sum = 0;
		for ( Size i = 1; i <= infiles.size() * nstruct; ++i ) {
			boltz_sum += boltz[i];
		}
		for ( Size i = 1; i <= infiles.size() * nstruct; ++i ) {
			boltz[i] = boltz[i]/boltz_sum;
		}

		// write out entity resfile
		std::ofstream ent("ent.res");
		ent << poses[1].pdb_info()->nres() << std::endl;
		ent << "ALLAAxc EX1 EX ARO 2\nstart" << std::endl;




	} catch (utility::excn::Exception const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}//main
