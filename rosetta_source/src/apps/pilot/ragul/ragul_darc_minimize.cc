// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @brief
/// @author jk
#include <iostream>
#include <iomanip>
#include <fstream>
#include <ostream>
#include <string>
#include <sstream>
#include <cmath>
#include <map>

#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <core/scoring/Energies.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/after_opts.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/rigid/RB_geometry.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/BoundConstraint.hh>
#include <core/scoring/constraints/HarmonicFunc.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
// Utility Headers
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <utility/io/mpistream.hh>
#include <core/kinematics/MoveMap.hh>

using namespace core;
using namespace basic::options;
using namespace std;
using namespace core::scoring;
using namespace core::optimization;
using namespace basic::options::OptionKeys;
using namespace conformation;
using namespace core::pose::datacache;
using namespace core::scoring::constraints;
using namespace core::id;
using namespace protocols::moves;


OPT_KEY( Real, cst_force_constant )

static basic::Tracer TR( "apps.pilot.ragul_darc_minimize.main" );

/// General testing code
int main( int argc, char * argv [] ) {

	NEW_OPT( cst_force_constant, "coordinate constraint force constant", 0.5 );

	devel::init(argc, argv);

	TR << "Starting minimization and repacking" << std::endl;

	// scoring function
	scoring::ScoreFunctionOP scorefxn( ScoreFunctionFactory::create_score_function(STANDARD_WTS, SCORE12_PATCH) );
	scoring::ScoreFunctionOP repack_scorefxn( ScoreFunctionFactory::create_score_function(STANDARD_WTS, SCORE12_PATCH) );

	repack_scorefxn->set_weight( core::scoring::fa_dun, 0.1 ); // this works for BAFF (1oqe)

	// create pose for native pose from pdb
	pose::Pose native_pose;
	std::string const input_pdb_name ( basic::options::start_file() );
	core::import_pose::pose_from_pdb( native_pose, input_pdb_name );

	//calculate no. of heavy atoms in ligand
	core::Size lig_res_num = 0;
	for ( int j = 1, resnum = native_pose.total_residue(); j <= resnum; ++j ) {
		if (!native_pose.residue(j).is_protein()){
      lig_res_num = j;
      break;
    }
	}
		if (lig_res_num == 0){
			std::cout<<"Error, no ligand for PlaidFingerprint" << std::endl;
			exit(1);
		}
		conformation::Residue const & curr_rsd = native_pose.conformation().residue(lig_res_num);
		core::Real num_heavy_atoms = curr_rsd.nheavyatoms();

	//create tag for output filename
	int dot_index1 = input_pdb_name.rfind(".", input_pdb_name.size());
  assert(dot_index1 != -1 && "No dot found in filename");
	std::string tag = input_pdb_name.substr(0,dot_index1);
	std::string init_pdb = "init_" + tag + ".pdb";
	std::string mini_pdb = "mini_" + tag + ".pdb";
	std::string unbo_pdb = "unbo_" + tag + ".pdb";

	//Apply constraint
	if ( native_pose.residue( native_pose.fold_tree().root() ).aa() != core::chemical::aa_vrt ) {
		native_pose.append_residue_by_jump
			( *ResidueFactory::create_residue( native_pose.residue(1).residue_type_set().name_map( "VRT" ) ),
				native_pose.total_residue()/2 );
	}

	Size nres = native_pose.total_residue();
	Real const coord_sdev( option[ OptionKeys::cst_force_constant ] );
	// default is 0.5 (from idealize) -- maybe too small
	for ( Size i = 1; i<= nres; ++i ) {
		if ( (Size)i==(Size)native_pose.fold_tree().root() ) continue;
		Residue const & nat_i_rsd( native_pose.residue(i) );
		for ( Size ii = 1; ii<= nat_i_rsd.nheavyatoms(); ++ii ) {
			native_pose.add_constraint( new CoordinateConstraint(
				AtomID(ii,i), AtomID(1,nres), nat_i_rsd.xyz( ii ),
				new HarmonicFunc( 0.0, coord_sdev ) ) );
		}
	}

	scorefxn->set_weight( coordinate_constraint, 0.5 );

	(*scorefxn)(native_pose);
	native_pose.dump_scored_pdb( init_pdb, *scorefxn );
	TR << "Initial score: " << native_pose.energies().total_energies()[ total_score ] << std::endl;

	// setting degrees of freedom which can move during minimization - everything
	kinematics::MoveMap mm_all;
	mm_all.set_chi( true );
	mm_all.set_bb( true );
	mm_all.set_jump( true );

	// minimize protein
	TR << "Starting minimization...." << std::endl;
	AtomTreeMinimizer minimizer;
	MinimizerOptions min_options( "dfpmin", 0.00001, true, false );

	minimizer.run( native_pose, mm_all, *scorefxn, min_options );
	minimizer.run( native_pose, mm_all, *scorefxn, min_options );
	(*scorefxn)(native_pose);

	TR << "Post minimization 1 constrained score: " << native_pose.energies().total_energies()[ total_score ] << std::endl;

	native_pose.remove_constraints();
	(*scorefxn)(native_pose);

	TR << "Post minimization 1 UNconstrained score: " << native_pose.energies().total_energies()[ total_score ] << std::endl;

	// Setup packer task for repacking
  pack::task::PackerTaskOP base_packer_task( pack::task::TaskFactory::create_packer_task( native_pose ));
	base_packer_task->set_bump_check( false );
	base_packer_task->initialize_from_command_line();
	base_packer_task->or_include_current( true ); // jk absolutely critical for BAFF case, Tyr65 is unusual

	for ( Size ii = 1; ii <= native_pose.total_residue(); ++ii ) {
		base_packer_task->nonconst_residue_task(ii).restrict_to_repacking();
	}

	// First repack
	pack::pack_rotamers( native_pose, *repack_scorefxn, base_packer_task );

	// Report Scores
	(*scorefxn)(native_pose);
	native_pose.dump_scored_pdb( "repacked_once.pdb", *scorefxn );
	TR << "Score after repacking once: " << native_pose.energies().total_energies()[ total_score ] << std::endl << std::endl;

	// iterate over minimizing and repacking
	for ( Size iter = 1; iter <= 5; ++iter ) {
		minimizer.run( native_pose, mm_all, *scorefxn, min_options );
		minimizer.run( native_pose, mm_all, *scorefxn, min_options );
		minimizer.run( native_pose, mm_all, *scorefxn, min_options );
		(*scorefxn)(native_pose);
		TR << "Current score after minimizing: " << native_pose.energies().total_energies()[ total_score ] << std::endl << std::endl;
		pack::pack_rotamers( native_pose, *repack_scorefxn, base_packer_task );
		(*scorefxn)(native_pose);
		TR << "Current score after repacking: " << native_pose.energies().total_energies()[ total_score ] << std::endl << std::endl;
	}

	// final minimization
		minimizer.run( native_pose, mm_all, *scorefxn, min_options );
		minimizer.run( native_pose, mm_all, *scorefxn, min_options );
		minimizer.run( native_pose, mm_all, *scorefxn, min_options );
		minimizer.run( native_pose, mm_all, *scorefxn, min_options );
	(*scorefxn)(native_pose);
	native_pose.dump_scored_pdb( mini_pdb, *scorefxn );
	TR << "Final score: " << native_pose.energies().total_energies()[ total_score ] << std::endl << std::endl;

	TR << "Successfully finished minimizing input." << std::endl;


	//setup the unbound pose
	pose::Pose unbound_pose;
	unbound_pose = native_pose;
	core::Real const unbound_dist = 40.;
	Size const rb_jump = 1; // use the first jump as the one between partners
	protocols::rigid::RigidBodyTransMover trans_mover( unbound_pose, rb_jump );
  trans_mover.trans_axis( trans_mover.trans_axis() );
  trans_mover.step_size(unbound_dist);
  trans_mover.apply( unbound_pose );
	(*scorefxn)(unbound_pose);
	unbound_pose.dump_pdb( unbo_pdb );

	core::Real const bound_score = native_pose.energies().total_energy();
	core::Real const unbound_score = unbound_pose.energies().total_energy();
	std::cout << "bound score is: " << bound_score << std::endl;
	std::cout << "unbound score is: " << unbound_score << std::endl;
	core::Real diff = bound_score - unbound_score;
	std::cout << input_pdb_name << " Difference is:  " << diff << std::endl;
	core::Real norm_score = (diff/sqrt(num_heavy_atoms));
	std::cout << input_pdb_name << " num_heavy_atoms:  " << num_heavy_atoms <<" Norm_score is:  " << norm_score << std::endl;
	return 0;

}


