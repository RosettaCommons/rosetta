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
/// @author Liz Kellogg ekellogg@u.washington.edu

// libRosetta headers

#include <core/types.hh>

#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueMatcher.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/kinematics/MoveMap.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/pose/Pose.hh>

#include <basic/options/util.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/keys/OptionKeys.hh>

#include <devel/init.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>
#include <core/pack/task/ResfileReader.hh>

#include <fstream>
#include <iostream>
#include <sstream>
#include <ios>
#include <utility/io/izstream.hh>
#include <ObjexxFCL/format.hh>

// C++ headers
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

#include <basic/Tracer.hh>

//numeric::random::RandomGenerator RG(154313929); // <- Magic number, do not change it!!!

#include <core/types.hh>

#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/scoring/sasa.hh>
#include <core/scoring/rms_util.hh>
//#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
//#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/BoundConstraint.hh>
#include <core/scoring/func/Func.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/PDBInfo.hh>

#include <devel/init.hh>

#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/ddg.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <basic/basic.hh>
#include <basic/Tracer.hh>
#include <basic/database/open.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/BinarySilentStruct.hh>

#include <core/pack/pack_rotamers.hh>
#include <core/pack/rotamer_trials.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <basic/database/open.hh>

//protocols
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/moves/Mover.hh>
//#include <protocols/simple_filters/RmsdEvaluator.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
//#include <protocols/looprelax/looprelax_main.hh>
//#include <protocols/comparative_modeling/ConstraintRemodelMover.hh>

#include <utility/file/FileName.hh>
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/excn/Exceptions.hh>
//#include "james_util.hh" //for calculation of burial
#include <basic/Tracer.hh>
#include <time.h>
using basic::T;
using basic::Warning;
using basic::Error;

// C++ headers
#include <fstream>
#include <iostream>
#include <string>

//C++ filechek
#include <sys/stat.h>

//Auto Headers
#include <core/import_pose/import_pose.hh>


using namespace core;
using namespace std;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::scoring;
using namespace core::scoring::constraints;
using namespace ddg;
using namespace id;
using namespace protocols;
using namespace moves;


ScoreFunction&
reduce_fa_rep(float fraction_fa_rep, ScoreFunction & s){
	s.set_weight( core::scoring::score_type_from_name("fa_rep"),
		s.get_weight(core::scoring::score_type_from_name("fa_rep"))*fraction_fa_rep);
	return s;
}

void
minimize_with_constraints(pose::Pose & p, ScoreFunction & s){
	core::optimization::AtomTreeMinimizer min_struc;
	float const minimizer_tol = 0.0000001;
	core::optimization::MinimizerOptions options( "lbfgs_armijo_nonmonotone", minimizer_tol, true /*use_nb_list*/,
		false /*deriv_check_in*/, false /*deriv_check_verbose_in*/);
	options.nblist_auto_update( true );
	options.max_iter(5000); //otherwise, they don't seem to converge
	core::kinematics::MoveMap mm;
	if ( !basic::options::option[sc_min_only]() ) {
		mm.set_bb(true);
	} else {
		mm.set_bb(false);
	}
	mm.set_chi(true);
	if ( basic::options::option[OptionKeys::ddg::ramp_repulsive]() ) {
		ScoreFunctionOP one_tenth_orig(s.clone());
		reduce_fa_rep(0.1,*one_tenth_orig);
		min_struc.run(p,mm,*one_tenth_orig,options);
		ScoreFunctionOP one_third_orig(s.clone());
		reduce_fa_rep(0.33,*one_third_orig);
		min_struc.run(p,mm,*one_third_orig,options);
	}
	min_struc.run(p,mm,s,options);
}


void
setup_ca_constraints(pose::Pose & pose, ScoreFunction & s, float const CA_cutoff, float const /*cst_tol */){
	int nres = pose.size();
	for ( int i = 1; i <= nres; i++ ) {
		Vector const CA_i( pose.residue(i).xyz(" CA "));
		for ( int j = 1; j <=nres; j++ ) {
			Vector const CA_j(pose.residue(j).xyz(" CA "));
			Real const CA_dist = (CA_i - CA_j).length();
			if ( CA_dist < CA_cutoff && i != j ) {
				ConstraintCOP cst( ConstraintOP( new AtomPairConstraint( AtomID(pose.residue(i).atom_index(" CA "),i),AtomID(pose.residue(j).atom_index(" CA "),j),core::scoring::func::FuncOP(new BoundFunc(CA_dist-0.5,CA_dist+0.5,0.1,10,"BoundFunc"))) ) ); //bounded func
				pose.add_constraint(cst);
			}
		}
	}
	s.set_weight(atom_pair_constraint, basic::options::option[OptionKeys::ddg::constraint_weight]());
}


void
setup_movers(simple_moves::SmallMoverOP small, simple_moves::ShearMoverOP shear,
	core::Real small_H_angle_max, core::Real small_E_angle_max, core::Real small_L_angle_max,
	core::Real shear_H_angle_max, core::Real shear_E_angle_max, core::Real shear_L_angle_max){
	//smaller moves to boost acceptance rate
	small->angle_max( 'H', small_H_angle_max ); //def 0.5
	small->angle_max( 'E', small_E_angle_max ); //def 0.5
	small->angle_max( 'L', small_L_angle_max ); //def 0.75

	shear->angle_max( 'H', shear_H_angle_max); //def 2.0
	shear->angle_max( 'E', shear_E_angle_max); //def 2.0
	shear->angle_max( 'L', shear_L_angle_max); //def 3.0
}

double
run_mc(pose::Pose & p, ScoreFunctionOP s,
	core::Real temperature, int numstruct,
	std::string output_tag, bool output_pdbs){

	using namespace protocols;
	using namespace moves;
	using namespace core::scoring::constraints;

	core::Size nmoves = (core::Size)p.size()/4; //number of moves for each move type

	bool BACKBONE_MOVEMENT = !basic::options::option[OptionKeys::ddg::no_bb_movement]();
	core::Size nrounds = 1000;
	if ( !BACKBONE_MOVEMENT ) {
		nrounds = 100;
	}

	core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap() );
	mm->set_bb(true);

	//small mover
	simple_moves::SmallMoverOP small_mover( new simple_moves::SmallMover( mm, temperature, nmoves) ) ;
	simple_moves::ShearMoverOP shear_mover( new simple_moves::ShearMover(mm, temperature, nmoves) );

	setup_movers(small_mover,shear_mover,
		0.1,0.1,0.2,
		0.25,1,1.5);

	core::pack::task::PackerTaskOP pt(pack::task::TaskFactory::create_packer_task(p));
	pt->restrict_to_repacking();
	pt->or_include_current(true);
	protocols::simple_moves::RotamerTrialsMoverOP rottrial_mover( new protocols::simple_moves::RotamerTrialsMover(s,(*pt)) );
	moves::RandomMoverOP apply_random_move( new moves::RandomMover() );

	if ( BACKBONE_MOVEMENT ) {
		apply_random_move->add_mover( small_mover, .45);
		apply_random_move->add_mover( shear_mover, .45);
		apply_random_move->add_mover(rottrial_mover, 0.1);
	} else {
		apply_random_move->add_mover(rottrial_mover, 1);
	}

	//separate set of movers for low temperature (reduce angle-max)
	simple_moves::SmallMoverOP small_mover_low( new simple_moves::SmallMover( mm, (temperature*0.25), nmoves) ) ;
	simple_moves::ShearMoverOP shear_mover_low( new simple_moves::ShearMover(mm, (temperature*0.25), nmoves) );
	setup_movers(small_mover_low,shear_mover_low,
		0.05,0.05,0.1,
		0.25,1.0,1.5);

	moves::RandomMoverOP apply_random_move_low( new moves::RandomMover() );
	if ( BACKBONE_MOVEMENT ) {
		apply_random_move_low->add_mover( small_mover_low, .45);
		apply_random_move_low->add_mover( shear_mover_low, .45);
		apply_random_move_low->add_mover( rottrial_mover, 0.1);
	} else {
		apply_random_move_low->add_mover( rottrial_mover, 1);
	}

	pose::Pose init_pose(p);
	(*s)(init_pose);
	//silent file stuff
	bool write_silent_file = false;
	std::string silent_file_name="";
	core::io::silent::SilentFileData out_sfd;
	core::io::silent::SilentFileData out_sfd_min;
	utility::vector1<std::string> tags_done;

	if ( option[out::file::silent].user() ) {
		write_silent_file = true;
		silent_file_name = option[out::file::silent]();
		out_sfd= core::io::silent::SilentFileData(silent_file_name,false,true,"binary");
		out_sfd_min= core::io::silent::SilentFileData("min."+silent_file_name,false,true,"binary");
		struct stat stFileInfo;
		int file_out_stat = stat(silent_file_name.c_str(),
			&stFileInfo);

		if ( file_out_stat == 0 ) { //file already exists
			out_sfd.read_file(silent_file_name);
			//read everything back in and store
			//std::cout << "[WARNING]: changing the name of your output silent file from " << silent_file_name << " to " << silent_file_name << ".new" << std::endl;
			//silent_file_name += ".new";
			out_sfd.set_filename(silent_file_name);
			tags_done = out_sfd.tags();
		}
	}

	MonteCarloOP mc( new moves::MonteCarlo(init_pose,*s,temperature) );
	for ( int ns =1; ns<=numstruct; ns++ ) {
		time_t time_per_decoy = time(NULL);
		std::ostringstream curr;
		curr << ns;

		struct stat stFileInfo;
		int file_stat = stat((basic::options::option[OptionKeys::ddg::last_accepted_pose_dir]()+"lowest."+curr.str()+".pdb").c_str(),
			&stFileInfo);
		if ( file_stat != 0 || !output_pdbs || ((silent_file_name.compare("") != 0) && !out_sfd.has_tag("lowest."+curr.str())) ) { //file doesn't exist or file can exist and output_pdbs set to false

			std::cout << "nrounds is " << nrounds <<  " and backbone_movement status is " << BACKBONE_MOVEMENT << std::endl;
			mc->reset_counters();
			mc->reset(init_pose);
			mc->set_temperature(temperature);

			p = init_pose;
			moves::TrialMoverOP tm( new moves::TrialMover(apply_random_move,mc) );
			RepeatMoverOP full_cycle( new moves::RepeatMover( tm, nrounds ) );
			full_cycle->apply( p );

			//std::cout << "ended trial at high temp" << std::endl;
			//mc->show_scores();
			//mc->show_counters();
			//std::cout << "monte carlo scores at high temp" << std::endl;

			pose::Pose last_accepted = mc->last_accepted_pose();

			mc->reset_counters();
			mc->set_lowest_score_pose(last_accepted); //high energy mark for annealing back to native

			core::Real low_temp = temperature * 0.25;
			mc->set_temperature(low_temp);

			//mc->show_scores();
			//mc->show_counters();
			//std::cout << "end double-checking" << std::endl;

			moves::TrialMoverOP ltm( new moves::TrialMover(apply_random_move_low,mc) );
			RepeatMoverOP full_cycle_2( new moves::RepeatMover( ltm, nrounds ) );
			full_cycle_2->apply( p );

			//std::cout << "end low temp trial, showing monte carlo stats" << std::endl;
			//mc->show_scores();
			//mc->show_counters();

			if ( !write_silent_file ) {
				core::io::pdb::dump_pdb(mc->lowest_score_pose(), basic::options::option[OptionKeys::ddg::last_accepted_pose_dir]()+"lowest."+curr.str()+".pdb");
			} else {
				core::io::silent::BinarySilentStruct bss(mc->lowest_score_pose(),"lowest."+curr.str());
				out_sfd.write_silent_struct(bss,silent_file_name,false);
			}

			pose::Pose lowest_score_pose = mc->lowest_score_pose();

			if ( basic::options::option[OptionKeys::ddg::min_with_cst]() ) {
				ScoreFunctionOP s = get_score_function();

				lowest_score_pose.remove_constraints((lowest_score_pose.constraint_set())->get_all_constraints());
				setup_ca_constraints(lowest_score_pose, (*s), 9.0, 0.5);
				minimize_with_constraints(lowest_score_pose, (*s));
				if ( !write_silent_file ) {
					lowest_score_pose.dump_pdb(output_tag+"."+curr.str()+"_0001.pdb");
				} else {
					//std::cout << "C-alpha rmsd from start: " << core::scoring::CA_rmsd(lowest_score_pose,init_pose) << std::endl;
					core::io::silent::BinarySilentStruct bss(lowest_score_pose,"min_cst.lowest."+curr.str()+"_0001");
					out_sfd_min.write_silent_struct(bss,"min."+silent_file_name,false);
				}
			}

		} else { //if(file_stat != 0){
			std::cout << "file:  " << (basic::options::option[OptionKeys::ddg::last_accepted_pose_dir]()+"lowest."+curr.str()+".pdb") << " already exists, skipping to next iteration" << std::endl;
		}
		time_t time_per_decoy_finish = time(NULL);
		std::cout << "time to finish decoy " << (time_per_decoy_finish-time_per_decoy) << std::endl;
	}

	return 1.0;
}


void
create_ensemble(pose::Pose & p, ScoreFunctionOP s, std::string output_tag){

	using namespace protocols;
	using namespace moves;

	core::Real temperature=0;
	temperature = basic::options::option[OptionKeys::ddg::temperature]();
	/*double avg_ca_rmsd =*/ run_mc(p,s,temperature,basic::options::option[out::nstruct](),output_tag,true);
}


int
main( int argc, char* argv [] )
{
	try {
		using namespace core;
		using namespace core::pose;
		using namespace utility;
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace core::scoring;
		using namespace core::scoring::constraints;
		// options, random initialization. MAKE SURE THIS COMES BEFORE OPTIONS
		devel::init( argc, argv );
		Real cst_tol = basic::options::option[ OptionKeys::ddg::harmonic_ca_tether ]();
		pose::Pose pose;

		ScoreFunctionOP scorefxn( get_score_function());

		vector1<file::FileName> files;
		if ( basic::options::option[in::file::s].user() ) {
			std::cout << "using -s option" << std::endl;
			files=option[in::file::s]();
		} else if ( option[in::file::l].user() ) {
			std::cout << "using -l option " << std::endl;

			utility::vector1<file::FileName> list = basic::options::option[ in::file::l ]();
			for ( unsigned int h=1; h<=list.size(); h++ ) {
				utility::io::izstream pdbs(list[h]);
				std::string fname;
				while ( pdbs >> fname ) {
					files.push_back(fname);
				}
			}
		}


		for ( unsigned int f=1; f<=files.size(); f++ ) {
			core::import_pose::pose_from_file(pose, files[f], core::import_pose::PDB_file);
			setup_ca_constraints(pose,(*scorefxn),9.0,cst_tol);
			ConstraintSetCOP cs = pose.constraint_set();
			std::string output = pose.pdb_info()->name();
			std::string no_pdb = output.erase(output.find(".pdb",0));
			create_ensemble(pose, scorefxn , (no_pdb.erase(0,(no_pdb.find_last_of("/")+1))));
		}
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}
