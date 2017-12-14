// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// awatkins: based heavily on kdrew/oop_creator.cc

// Project Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/subpose_manipulation_util.hh>
#include <core/pose/ncbb/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/util.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/Patch.hh>
#include <core/chemical/VariantType.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>

#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

// Mover headers
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/PyMOLMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/RandomTorsionMover.hh>
#include <protocols/simple_moves/hbs/HbsPatcher.hh>
#include <protocols/simple_moves/a3b_hbs/A3BHbsPatcher.hh>
#include <protocols/simple_moves/chiral/ChiralMover.hh>
#include <protocols/rigid/RB_geometry.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/simple_moves/TaskAwareMinMover.hh>

#include <numeric/conversions.hh>
#include <numeric/random/random.hh>
#include <numeric/xyzVector.hh>


//Basic headers
#include <basic/resource_manager/ResourceManager.hh>

// Utility Headers
#include <devel/init.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
//#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
//#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/tools/make_vector1.hh>

// C++ headers
#include <string>
#include <sstream>

//The original author used a lot of using declarations here.  This is a stylistic choice.
// Namespaces
using namespace core;
using namespace conformation;
using namespace chemical;
using namespace scoring;
using namespace pose;
using namespace protocols;
using namespace protocols::moves;
using namespace protocols::simple_moves;
using namespace protocols::simple_moves::a3b_hbs;
using namespace protocols::simple_moves::chiral;
using namespace core::pack::task;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::id;
using basic::Error;
using basic::Warning;
using utility::file::FileName;

//kdrew: this app adds hbs patches to the given pdb strucure

// tracer - used to replace cout
static basic::Tracer TR("A3BHBS_Creator");

void
setup_pert_foldtree(
	core::pose::Pose & pose
);

// application specific options
namespace a3b_hbs_creator {
// pert options
StringOptionKey const hbs_chain ( "a3b_hbs_creator::hbs_chain" );
IntegerOptionKey const hbs_final_res ( "a3b_hbs_creator::hbs_final_res" );
IntegerOptionKey const hbs_length ( "a3b_hbs_creator::hbs_length" );
IntegerOptionKey const hbs_offset( "a3b_hbs_creator::hbs_offset" );
IntegerOptionKey const change_resnums( "a3b_hbs_creator::change_resnums" );
BooleanOptionKey const bridgehead_gly( "a3b_hbs_creator::bridgehead_gly" );
BooleanOptionKey const final_repack( "a3b_hbs_creator::final_repack" );
BooleanOptionKey const final_minimize( "a3b_hbs_creator::final_minimize" );
BooleanOptionKey const final_mc ( "a3b_hbs_creator::final_mc" );
}

class A3BHbsCreatorMover : public Mover {

public:

	//default ctor
	A3BHbsCreatorMover();

	//default dtor
	~A3BHbsCreatorMover() override= default;

	void apply( core::pose::Pose & pose ) override;
	void repack( core::pose::Pose & pose );
	void scan( core::pose::Pose & pose, char const hbs_chn );

	void make_a3b_pose(
		pose::Pose & pose,
		pose::Pose & a3bpose );

	void
	delete_extra_residues(
		core::pose::Pose & pose );

	void
	add_hbond_and_omega_constraints_starting_at_seqpos(
		Pose & pose,
		Size const seqpos );

	void
	do_mc(
		core::pose::Pose & pose );

	std::string get_name() const override { return "A3BHbsCreatorMover"; }

private:
	ScoreFunctionOP score_fxn_;
	ScoreFunctionOP score_fxn_cart_;
	Size hbs_length_;
	bool mut_gly_;
	std::map< Size, Vector > former_CAs_;
	Size offset_;
	char hbs_chain_;
	utility::vector1< Size > oldnums_;
};

using A3BHbsCreatorMoverOP = utility::pointer::shared_ptr<A3BHbsCreatorMover>;
using A3BHbsCreatorMoverCOP = utility::pointer::shared_ptr<const A3BHbsCreatorMover>;

A3BHbsCreatorMover::A3BHbsCreatorMover():
	Mover("A3BHbsCreatorMover"),
	hbs_length_( option[a3b_hbs_creator::hbs_length].value() ),
	mut_gly_( option[a3b_hbs_creator::bridgehead_gly].value() ),
	former_CAs_( std::map< Size, Vector >() ),
	offset_( ( option[a3b_hbs_creator::hbs_offset].value() == 1000 ) ? Size( numeric::random::uniform() * 4 ) : option[a3b_hbs_creator::hbs_offset].value() ),
	hbs_chain_( option[a3b_hbs_creator::hbs_chain].value()[0] )
{
	using namespace core::scoring::constraints;
	// create score function
	score_fxn_ = scoring::get_score_function();
	add_fa_constraints_from_cmdline_to_scorefxn(*score_fxn_);
	score_fxn_->set_weight_if_zero( atom_pair_constraint, 0.1 );
	score_fxn_->set_weight_if_zero( dihedral_constraint, 0.1 );
	score_fxn_->set_weight_if_zero( angle_constraint, 0.1 );

	// soften rep
	score_fxn_->set_weight_if_zero( fa_rep, score_fxn_->get_weight( fa_rep ) * 0.5 );

	score_fxn_cart_ = score_fxn_->clone();
	score_fxn_cart_->set_weight( pro_close, 0 );
	score_fxn_cart_->set_weight( cart_bonded, 1 );
}

std::string alpha_to_beta( std::string const & alpha ) {
	TR << "Converting name for " << alpha << std::endl;

	std::string postfix = alpha;
	postfix.erase( 0, postfix.find( ":" ) );

	std::string prefix = alpha;
	if ( prefix.find(":") != std::string::npos ) {
		prefix.erase( prefix.find(":"), std::string::npos-1 );
	}

	if ( prefix == "ALA" ) {
		prefix = "B3A";
	} else if ( prefix == "CYS" ) {
		prefix = "B3C";
	} else if ( prefix == "ASP" ) {
		prefix = "B3D";
	} else if ( prefix == "GLU" ) {
		prefix = "B3E";
	} else if ( prefix == "PHE" ) {
		prefix = "B3F";
	} else if ( prefix == "GLY" ) {
		prefix = "B3G";
	} else if ( prefix == "HIS" ) {
		prefix = "B3H";
	} else if ( prefix == "HIS_D" ) {
		prefix = "B3H";
	} else if ( prefix == "ILE" ) {
		prefix = "B3I";
	} else if ( prefix == "LYS" ) {
		prefix = "B3K";
	} else if ( prefix == "LEU" ) {
		prefix = "B3L";
	} else if ( prefix == "MET" ) {
		prefix = "B3M";
	} else if ( prefix == "ASN" ) {
		prefix = "B3N";
	} else if ( prefix == "PRO" ) {
		prefix = "B3P";
	} else if ( prefix == "GLN" ) {
		prefix = "B3Q";
	} else if ( prefix == "ARG" ) {
		prefix = "B3R";
	} else if ( prefix == "SER" ) {
		prefix = "B3S";
	} else if ( prefix == "THR" ) {
		prefix = "B3T";
	} else if ( prefix == "VAL" ) {
		prefix = "B3V";
	} else if ( prefix == "TRP" ) {
		prefix = "B3W";
	} else if ( prefix == "TYR" ) {
		prefix = "B3Y";
	}
	return prefix + postfix;
}

void A3BHbsCreatorMover::repack(
	core::pose::Pose & pose
) {
	// create a task factory and task operations
	TaskFactoryOP tf( new TaskFactory );
	tf->push_back( operation::TaskOperationCOP( new operation::InitializeFromCommandline ) );
	tf->push_back( operation::TaskOperationCOP( new operation::IncludeCurrent ) );

	using namespace basic::resource_manager;
	if ( ResourceManager::get_instance()->has_option( packing::resfile ) ||  option[ packing::resfile ].user() ) {
		operation::ReadResfileOP rrop( new operation::ReadResfile );
		rrop->default_filename();
		tf->push_back( rrop );
	} else {
		//kdrew: do not do design, makes NATAA if res file is not specified
		operation::RestrictToRepackingOP rtrp( new operation::RestrictToRepacking() );
		tf->push_back( rtrp );
	}

	// create a pack rotamers mover
	simple_moves::PackRotamersMoverOP packer( new protocols::simple_moves::PackRotamersMover );
	packer->task_factory( tf );
	packer->score_function( score_fxn_ );
	packer->apply(pose);
}

void A3BHbsCreatorMover::do_mc(
	Pose & pose
) {
	// create a monte carlo object for the full cycle
	moves::MonteCarloOP mc( new moves::MonteCarlo( pose, *score_fxn_, 1 ) );

	// create a monte carlo object for the pertubation phase
	moves::MonteCarloOP pert_mc( new moves::MonteCarlo( pose, *score_fxn_, 1 ) );

	// create a rigid body mover to move the peptide around in the pocket
	rigid::RigidBodyPerturbMoverOP pert_dock_rbpm( new rigid::RigidBodyPerturbMover(1, 1, 1 ) );

	// get peptide start and end positions
	Size pep_start( pose.conformation().chain_begin( 2 ) ); Size pep_end( pose.size() );
	TR << "pep_start: " << pep_start << " pep_end: " << pep_end << std::endl;

	// create movemap for peptide
	kinematics::MoveMapOP pert_alpha_mm( new kinematics::MoveMap() );
	kinematics::MoveMapOP pert_beta_mm( new kinematics::MoveMap() );

	core::Size hbs_seq_position = 0;
	core::Size hbs_length = 0;
	for ( core::Size i = 1; i <= pose.size(); ++i ) {

		if ( i >= pep_start+3 && i <= pep_end ) {
			// movemap settings
			if ( pose.residue(i).type().is_alpha_aa() ) {
				pert_alpha_mm->set_bb( i, true );
				pert_beta_mm->set_bb( i, false );
			} else if ( pose.residue(i).type().is_beta_aa() ) {
				pert_beta_mm->set_bb( i, true );
				pert_alpha_mm->set_bb( i, false );
			}
		}

		if ( hbs_seq_position>0 && hbs_seq_position <= i ) {
			hbs_length++;
		}
	}
	assert(hbs_seq_position != 0);

	// create small and shear movers
	simple_moves::SmallMoverOP pert_pep_alpha( new simple_moves::SmallMover( pert_alpha_mm, 1, 1 ) );
	pert_pep_alpha->angle_max( 'H', 2.0 );
	pert_pep_alpha->angle_max( 'L', 2.0 );
	pert_pep_alpha->angle_max( 'E', 2.0 );
	simple_moves::RandomTorsionMoverOP pert_pep_beta( new simple_moves::RandomTorsionMover( pert_beta_mm, 1, 1 ) );

	// create random mover
	moves::RandomMoverOP pert_pep_random( new moves::RandomMover() );
	pert_pep_random->add_mover( pert_pep_alpha, .75 );
	pert_pep_random->add_mover( pert_pep_beta, .25 );
	moves::RepeatMoverOP pert_pep_repeat( new moves::RepeatMover( pert_pep_random, 10 ) );

	/******************************************************************************
	Rotamer Trials Setup
	*******************************************************************************/

	// create a task factory and task operations
	TaskFactoryOP pert_tf(new TaskFactory());
	pert_tf->push_back( operation::TaskOperationCOP( new operation::InitializeFromCommandline ) );

	operation::ReadResfileOP pert_rrop( new operation::ReadResfile() );
	pert_rrop->default_filename();
	pert_tf->push_back( pert_rrop );

	operation::RestrictToRepackingOP pert_rtrp( new operation::RestrictToRepacking() );
	pert_tf->push_back( pert_rtrp );

	// create a rotamer trials mover
	simple_moves::RotamerTrialsMoverOP pert_rt(new simple_moves::EnergyCutRotamerTrialsMover( score_fxn_, pert_tf, pert_mc, 0.1 /*energycut*/ ) );

	/*********************************************************
	Common Setup
	**********************************************************/

	// create a random mover to hold the docking, and peptide pertubation movers
	moves::RandomMoverOP pert_random( new moves::RandomMover() );
	pert_random->add_mover( pert_dock_rbpm, 1 );
	pert_random->add_mover( pert_pep_repeat, 0.5 );

	// create a sequence move to hold random and rotamer trials movers
	moves::SequenceMoverOP pert_sequence( new moves::SequenceMover() );
	pert_sequence->add_mover( pert_random );
	pert_sequence->add_mover( pert_rt );

	// create a TrialMover for the pertubation
	moves::TrialMoverOP pert_trial( new moves::TrialMover( pert_sequence, pert_mc ) );

	/*********************************************************
	Design Min Phase
	**********************************************************/
	/*********************************************************
	Design Setup
	**********************************************************/

	// create a task factory and task operations
	TaskFactoryOP desn_tf( new TaskFactory() );
	desn_tf->push_back( operation::TaskOperationCOP( new operation::InitializeFromCommandline ) );

	operation::ReadResfileOP desn_rrop( new operation::ReadResfile() );
	desn_rrop->default_filename();
	desn_tf->push_back( desn_rrop );

	// create a pack rotamers mover
	simple_moves::PackRotamersMoverOP desn_pr( new simple_moves::PackRotamersMover() );
	desn_pr->task_factory( desn_tf );
	desn_pr->score_function( score_fxn_ );

	/*********************************************************
	Minimize Setup
	**********************************************************/

	// create move map for minimization
	kinematics::MoveMapOP desn_mm( new kinematics::MoveMap() );
	desn_mm->set_bb( false );
	desn_mm->set_bb_true_range( pep_start, pep_end );
	desn_mm->set_chi( true );
	desn_mm->set_jump( 1, true );

	// create minimization mover
	simple_moves::MinMoverOP desn_min( new simple_moves::MinMover( desn_mm, score_fxn_, option[ OptionKeys::run::min_type ].value(), 0.01, true ) );

	//definitely want sidechain minimization here
	using protocols::simple_moves::TaskAwareMinMoverOP;
	using protocols::simple_moves::TaskAwareMinMover;
	TaskAwareMinMoverOP desn_ta_min = TaskAwareMinMoverOP( new TaskAwareMinMover( desn_min, desn_tf ) );

	/*********************************************************
	Common Setup
	**********************************************************/

	// create a sequence mover to hold pack rotamers and minimization movers
	moves::SequenceMoverOP desn_sequence( new moves::SequenceMover() );
	desn_sequence->add_mover( desn_pr );
	desn_sequence->add_mover( desn_ta_min );

	TR << "Main loop..." << std::endl;

	protocols::jd2::JobOP curr_job( protocols::jd2::JobDistributor::get_instance()->current_job() );

	//pose.dump_pdb("pre_main_loop.pdb");
	for ( Size k = 1; k <= 5; ++k ) {
		pert_mc->reset(pose);

		// pert loop
		for ( Size j = 1; j <= 10; ++j ) {
			TR << "PERTURB: " << k << " / "  << j << std::endl;
			pert_trial->apply( pose );
			curr_job->add_string_real_pair( "ENERGY_PERT (pert score)", (*score_fxn_)(pose) );
		}
		pert_mc->recover_low( pose );
		curr_job->add_string_real_pair( "ENERGY_PERT (pert score) recovered low", (*score_fxn_)(pose) );

		// design
		TR << "DESIGN: " << k << std::endl;
		desn_sequence->apply( pose );
		curr_job->add_string_real_pair( "ENERGY_DESN (hard score)", (*score_fxn_)(pose) );

		TR<< "pre mc->boltzmann" << std::endl;
		mc->show_state();
		mc->boltzmann( pose );
		TR<< "post mc->boltzmann" << std::endl;
		mc->show_state();

	}//dock_design for loop

	mc->recover_low( pose );
}

void A3BHbsCreatorMover::scan(
	core::pose::Pose & pose,
	char const hbs_chn
) {
	core::pose::Pose test_pose = pose;
	core::pose::PDBInfoCOP pdb_info( test_pose.pdb_info() );

	kinematics::MoveMapOP littlemm( new kinematics::MoveMap() );
	littlemm->set_bb( false );
	for ( core::Size i = 1; i <= test_pose.size(); ++i ) {
		if ( pdb_info->chain(i) != hbs_chn ) {
			continue;
		}

		littlemm->set_bb( true );
	}

	littlemm->set_chi( true );
	littlemm->set_jump( 1, true );
	simple_moves::MinMoverOP littlemin( new protocols::simple_moves::MinMover( littlemm, score_fxn_, option[ OptionKeys::run::min_type ].value(), 1, true ) );

	core::Real bin = 15;

	Real best_score = 999999;
	core::pose::Pose best_pose;
	for ( Real alphaphi = -84; alphaphi <= -34; alphaphi += bin ) {
		for ( Real alphapsi = -80; alphapsi <= -30; alphapsi += bin ) {
			for ( Real betaphi = -131; betaphi <= -91; betaphi += bin ) {
				for ( Real betatht = 40; betatht <= 90; betatht += bin ) {
					for ( Real betapsi = -133; betapsi <= -83; betapsi += bin ) {

						TR << "Eval (" << alphaphi << ", " << alphapsi << ") (" << betaphi << ", " << betatht << ", " << betapsi << ")" << std::endl;

						bool first = true;
						for ( Size i = 1; i <= test_pose.size(); ++i ) {
							if ( pdb_info->chain(i) != hbs_chn ) {
								continue;
							}
							if ( pdb_info->chain(i) == hbs_chn && first ) {
								first = false;
								if ( test_pose.residue(i).type().is_beta_aa() ) {
									core::id::TorsionID tht( i, id::BB, 2);
									core::id::TorsionID psi( i, id::BB, 3);
									test_pose.conformation().set_torsion( tht, betatht);
									test_pose.conformation().set_torsion( psi, betapsi);

								} else {
									core::id::TorsionID psi( i, id::BB, 2);
									test_pose.conformation().set_torsion( psi, alphapsi);
								}

							} else {
								if ( test_pose.residue(i).type().is_beta_aa() ) {
									core::id::TorsionID phi( i, id::BB, 1);
									core::id::TorsionID tht( i, id::BB, 2);

									test_pose.conformation().set_torsion( phi, betaphi);
									test_pose.conformation().set_torsion( tht, betatht);
									if ( i != pose.size() ) {
										core::id::TorsionID psi( i, id::BB, 3);
										test_pose.conformation().set_torsion( psi, betapsi);
									}
								} else {
									core::id::TorsionID phi( i, id::BB, 1);

									test_pose.conformation().set_torsion( phi, alphaphi);
									if ( i != pose.size() ) {
										core::id::TorsionID psi( i, id::BB, 2);
										test_pose.conformation().set_torsion( psi, alphapsi);
									}
								}
							}
						}


						// create minimization mover
						//TR << "Minimizing away any major clashes (score start: " << ( *score_fxn )( test_pose );
						//littlemin->apply( test_pose );
						//TR << "; score end: " << ( *score_fxn )( test_pose ) << std::endl;
						//repack( test_pose, score_fxn );
						//TR << "; after repack: " << ( *score_fxn )( test_pose ) << ")" << std::endl;


						Real score = ( *score_fxn_ )( test_pose );
						if ( score <= best_score ) {
							TR << "Improvement! to "<< score << std::endl;
							best_pose = test_pose;
							best_score = score;
						}

					}
				}
			}
		}
	}

	pose = best_pose;
}


int
main( int argc, char* argv[] )
{
	try {
		utility::vector1< core::Size > empty_vector(0);

		option.add( a3b_hbs_creator::hbs_chain, "Chain from PDB to be mimicked. Default 'A'. Use letters." ).def("A");
		option.add( a3b_hbs_creator::hbs_final_res, "Residue number of the final residue for mimicry. Default 1." ).def(1);
		option.add( a3b_hbs_creator::hbs_length, "Number of residues to mimic. Default 12." ).def(12);
		option.add( a3b_hbs_creator::hbs_offset, "Size of offset, 0-3. Default random." ).def(1000);
		option.add( a3b_hbs_creator::change_resnums, "You may want to change the residue numbering of the HBS chain. Default 0" ).def(0);
		option.add( a3b_hbs_creator::bridgehead_gly, "Mutate bridgehead to gly?" ).def(true);
		option.add( a3b_hbs_creator::final_repack, "Do a final repack. Default false" ).def(false);
		option.add( a3b_hbs_creator::final_minimize, "Do a final minimization. Default false" ).def(false);
		option.add( a3b_hbs_creator::final_mc, "Do a final monte carlo on hbs. Default false" ).def(false);

		// init command line options
		//you MUST HAVE THIS CALL near the top of your main function, or your code will crash when you first access the command line options
		devel::init(argc, argv);

		//create mover instance
		A3BHbsCreatorMoverOP HC_mover( new A3BHbsCreatorMover() );

		//call job distributor
		protocols::jd2::JobDistributor::get_instance()->go( HC_mover );
	} catch (utility::excn::Exception const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;

}//main

void
A3BHbsCreatorMover::delete_extra_residues(
	core::pose::Pose & pose
) {
	core::Size final_res = option[a3b_hbs_creator::hbs_final_res].value();

	for ( Size i = 1; i <= pose.size(); ++i ) {
		char chn = pose.pdb_info()->chain(i);
		if ( chn != hbs_chain_ ) continue;
		core::Size pdb_res_num = pose.pdb_info()->number(i);

		// hbs pre is the smallest number of what we want to preserve
		// AMW: this actually doesn't work at all.
		if ( pdb_res_num < final_res ) {
			while ( pdb_res_num < final_res ) {
				//TR << "deleting residue " << pdb_res_num  << " which was " << core::chemical::oneletter_code_from_aa(pose.aa(i)) << std::endl;
				pose.delete_polymer_residue(i);
				//TR << "deleted residue " << pdb_res_num << std::endl;
				pdb_res_num = pose.pdb_info()->number(i);
				//TR << "now residue " << i << " refers to " << pdb_res_num << std::endl;
			}

		} else if ( pdb_res_num > final_res + hbs_length_ ) {
			//TR << "deleting residue " << pdb_res_num << std::endl;
			while ( chn == hbs_chain_ && i <= pose.size() ) {
				chn = pose.pdb_info()->chain(i);
				pose.delete_polymer_residue(i);
			}
		}
	}
}

void A3BHbsCreatorMover::make_a3b_pose(
	pose::Pose & pose,
	pose::Pose & a3bpose
) {
	chemical::ResidueTypeSetCOP restype_set = chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
	pose::PDBInfoCOP old_pdb_info = pose.pdb_info();
	Size first_hbs_chain_res = 0;
	for ( Size i = 1; i <= pose.size(); ++i ) {

		char chn = pose.pdb_info()->chain(i);
		//TR << "evaluating residue " << chn  << " " << pdb_info->number(i) << std::endl;

		if ( chn != hbs_chain_ ) continue;
		if ( first_hbs_chain_res == 0 ) first_hbs_chain_res = i;

		// replace each HBS residue with an a3b
		TR << "Recording CA for " << i << std::endl;
		former_CAs_[ pose.pdb_info()->number(i) ] = pose.residue( i ).atom( "CA" ).xyz();

		if ( (i-first_hbs_chain_res) == 2 && mut_gly_ ) {
			conformation::Residue ala( restype_set->name_map( "GLY"+pose.residue_type( i ).name().substr(3) ), true );
			pose.replace_residue( i, ala, true );
			conformation::idealize_position( i, pose.conformation() );
		}

		std::string name = (i-first_hbs_chain_res) % 4 == offset_ ?
			alpha_to_beta( pose.residue(i).name() ) : pose.residue(i).name();

		Residue r = *new Residue( restype_set->name_map( name ), true );
		r.set_all_chi(pose.residue(i).chi());

		if ( a3bpose.size() == 0 ) {
			a3bpose.append_residue_by_jump( r , 1 );
		} else {
			a3bpose.append_residue_by_bond( r, true );
		}
	}
}

void A3BHbsCreatorMover::add_hbond_and_omega_constraints_starting_at_seqpos(
	Pose & pose,
	Size const seqpos
) {
	using namespace core::scoring::constraints;
	AtomID fixed_atom( pose.residue(1).atom_index( "CA" ), 1 );

	pose.add_constraint(
		ConstraintOP(
		new AtomPairConstraint(
		AtomID( pose.residue( seqpos ).atom_index( "OY" ), seqpos ),
		AtomID( pose.residue( seqpos+3 ).atom_index( "H"  ), seqpos+3 ),
		core::scoring::func::HarmonicFuncOP( new core::scoring::func::HarmonicFunc( 1.9, 0.04 ) )
		) ) );

	for ( Size i = seqpos; i <= pose.size(); ++i ) {
		if ( oldnums_.size() > 0 ) {
			TR << "Setting chain for " << i << " to " << hbs_chain_ << std::endl;
			pose.pdb_info()->chain( i, hbs_chain_ );
			pose.pdb_info()->number( i, oldnums_[ i - seqpos + 1 ] );
		}
		if ( i <= pose.size() - 4 ) {
			pose.add_constraint(
				ConstraintOP(
				new AtomPairConstraint(
				AtomID( pose.residue( i   ).atom_index( "O" ), i ),
				AtomID( pose.residue( i+4 ).atom_index( "H" ), i+4 ),
				core::scoring::func::HarmonicFuncOP( new core::scoring::func::HarmonicFunc( 1.9, 0.04 ) )
				) ) );
		}

		std::string i_omg_atom = pose.residue_type( i ).is_beta_aa() ? "CM" : "CA";
		if ( i < pose.size() ) {
			pose.add_constraint(
				ConstraintOP(
				new DihedralConstraint(
				AtomID( pose.residue( i ).atom_index( i_omg_atom ), i ),
				AtomID( pose.residue( i ).atom_index( "C"  ), i ),
				AtomID( pose.residue( i+1 ).atom_index( "N"  ), i+1 ),
				AtomID( pose.residue( i+1 ).atom_index( "CA"  ), i+1 ),
				core::scoring::func::CircularHarmonicFuncOP( new core::scoring::func::CircularHarmonicFunc( 3.14159, 0.04 ) )
				) ) );
		}

		// No coordinate constraints within the macrocycle, I think.
		// AMW TODO: specify "helical portion" over which these are defined.
		/*if ( former_CAs_.find( pose.pdb_info()->number(i) ) != former_CAs_.end()
		&& i > seqpos+2 ) {

		pose.add_constraint(
		ConstraintOP(
		new CoordinateConstraint(
		AtomID( pose.residue( i ).atom_index( "CA" ), i ),
		fixed_atom,
		former_CAs_[ pose.pdb_info()->number(i) ],
		core::scoring::func::HarmonicFuncOP( new core::scoring::func::HarmonicFunc( 0.0, 0.1 ) )
		) ) );
		}*/
	}
	pose.add_constraint(
		ConstraintOP(
		new AtomPairConstraint(
		AtomID( pose.residue( pose.size()-3   ).atom_index( "O" ), pose.size()-3 ),
		AtomID( pose.residue( pose.size() ).atom_index( "HM" ), pose.size() ),
		core::scoring::func::HarmonicFuncOP( new core::scoring::func::HarmonicFunc( 1.9, 0.04 ) )
		) ) );

	Size const nres = pose.size();
	std::string const i_omg_atom = pose.residue_type( nres ).is_beta_aa() ? "CM" : "CA";
	pose.add_constraint(
		ConstraintOP(
		new DihedralConstraint(
		AtomID( pose.residue_type( nres ).atom_index( i_omg_atom ), nres ),
		AtomID( pose.residue_type( nres ).atom_index( "C"   ), nres ),
		AtomID( pose.residue_type( nres ).atom_index( "NM"  ), nres ),
		AtomID( pose.residue_type( nres ).atom_index( "CN"  ), nres ),
		core::scoring::func::CircularHarmonicFuncOP( new core::scoring::func::CircularHarmonicFunc( 3.14159, 0.04 ) )
		) ) );
}

void
A3BHbsCreatorMover::apply(
	core::pose::Pose & pose
) {
	using namespace core::scoring::constraints;

	add_fa_constraints_from_cmdline_to_pose(pose);

	chemical::ResidueTypeSetCOP restype_set = chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

	/******************************************\
	*  FOLD TREE MANIPULATION                *
	\******************************************/
	// shift the jump
	// account for the position of the offset so we never delete the jump res
	/*kinematics::FoldTree f = pose.fold_tree();
	f.slide_jump( 1, 1, pose.pdb_info()->pdb2pose( hbs_chain_, final_res+1+offset ) );
	pose.fold_tree( f );
	*/

	delete_extra_residues( pose );

	TR << "making an a3b peptide" << std::endl;

	// CONVERT TO A3B
	pose::Pose a3bpose;
	make_a3b_pose( pose, a3bpose );

	// NOW set dihedrals because you've added all the residues
	// otherwise you'll fail at setting psi EVERY TIME you doofus
	for ( Size i = 1; i <= a3bpose.size(); ++i ) {
		TR << "Setting dihedrals for " << i << std::endl;
		if ( a3bpose.residue_type( i ).is_beta_aa() ) {
			a3bpose.conformation().set_torsion( TorsionID( i, id::BB, 1), -105);
			a3bpose.conformation().set_torsion( TorsionID( i, id::BB, 2),   65);
			a3bpose.conformation().set_torsion( TorsionID( i, id::BB, 3), -115);
			a3bpose.conformation().set_torsion( TorsionID( i, id::BB, 4),  180);
		} else {
			a3bpose.conformation().set_torsion( TorsionID( i, id::BB, 1),  -60);
			a3bpose.conformation().set_torsion( TorsionID( i, id::BB, 2),  -47);
			a3bpose.conformation().set_torsion( TorsionID( i, id::BB, 3),  180);
		}
	}

	// Mutate any prolines
	for ( Size i = 1; i <= a3bpose.size(); ++i ) {
		if ( a3bpose.residue_type( i ).is_beta_aa() ) {
			if ( a3bpose.residue_type( i ).name3() != "B3P" ) continue; //&& pose.residue_type( i ).name3() != "B3P" ) continue;

			conformation::Residue ala( restype_set->name_map( "B3A"+a3bpose.residue_type( i ).name().substr(3) ), true );
			a3bpose.replace_residue( i, ala, true );
			//conformation::idealize_position( i, a3bpose.conformation() );
		} else {
			if ( a3bpose.residue_type( i ).name3() != "PRO" ) continue; //&& pose.residue_type( i ).name3() != "B3P" ) continue;

			conformation::Residue ala( restype_set->name_map( "ALA"+a3bpose.residue_type( i ).name().substr(3) ), true );
			a3bpose.replace_residue( i, ala, true );
			//conformation::idealize_position( i, a3bpose.conformation() );
		}
	}

	kinematics::MoveMapOP promutmm( new kinematics::MoveMap() );
	promutmm->set_bb( true );
	promutmm->set_chi( true );
	simple_moves::MinMoverOP promutmin( new protocols::simple_moves::MinMover( promutmm, score_fxn_/*cart_*/, "lbfgs_armijo_nonmonotone", 1, true ) );
	promutmin->cartesian( false );//true );
	promutmin->apply( a3bpose );

	for ( Size i = 1; i <= a3bpose.size() - 4; ++i ) {
		TR << "constraining " << i << std::endl;

		a3bpose.add_constraint(
			ConstraintOP(
			new AtomPairConstraint(
			AtomID( a3bpose.residue( i   ).atom_index( "O" ), i ),
			AtomID( a3bpose.residue( i+4 ).atom_index( "H" ), i+4 ),
			core::scoring::func::HarmonicFuncOP( new core::scoring::func::HarmonicFunc( 1.9, 0.04 ) )
			) ) );
	}

	// presently the final residue in the pose is the terminal residue of the hbs
	// replace with terminal variant
	conformation::Residue term( restype_set->get_residue_type_with_variant_added(a3bpose.residue(a3bpose.size()).type(), chemical::METHYLATED_CTERMINUS_VARIANT), true );
	term.set_all_chi(a3bpose.residue(a3bpose.size()).chi());
	a3bpose.replace_residue( a3bpose.size(), term, true );
	conformation::idealize_position( a3bpose.size(), a3bpose.conformation() );
	a3bpose.add_constraint(
		ConstraintOP(
		new AtomPairConstraint(
		AtomID( a3bpose.residue( a3bpose.size()-3   ).atom_index( "O" ), a3bpose.size()-3 ),
		AtomID( a3bpose.residue( a3bpose.size() ).atom_index( "HM" ), a3bpose.size() ),
		core::scoring::func::HarmonicFuncOP( new core::scoring::func::HarmonicFunc( 1.9, 0.04 ) )
		) ) );

	// PATCH
	// If offset is three then we need a special HBS patch for this pattern
	if ( offset_ == 3 ) {
		a3b_hbs::A3BHbsPatcherOP hbs_patcher( new a3b_hbs::A3BHbsPatcher( 1 ) );
		hbs_patcher->apply( a3bpose );
	} else {
		hbs::HbsPatcherOP hbs_patcher( new hbs::HbsPatcher( 1 ) );
		hbs_patcher->apply( a3bpose );
	}

	a3bpose.conformation().declare_chemical_bond( 1, "CYH", 3, "CZH" );

	a3bpose.conformation().detect_bonds();
	a3bpose.conformation().detect_pseudobonds();
	for ( core::Size i = 1; i <= a3bpose.size(); ++i ) {
		a3bpose.conformation().update_polymeric_connection(i);
	}

	add_hbond_and_omega_constraints_starting_at_seqpos( a3bpose, 1 );

	std::string cyn = offset_ == 0 ? "CY3" : "CY2";
	a3bpose.add_constraint(
		ConstraintOP(
		new DihedralConstraint(
		AtomID( a3bpose.residue( 1 ).atom_index( cyn ), 1 ),
		AtomID( a3bpose.residue( 1 ).atom_index( "CYH"  ), 1 ),
		AtomID( a3bpose.residue( 3 ).atom_index( "CZH"  ), 3 ),
		AtomID( a3bpose.residue( 3 ).atom_index( "N"  ), 3 ),
		core::scoring::func::CircularHarmonicFuncOP( new core::scoring::func::CircularHarmonicFunc( 3.14159, 0.1 ) )
		) ) );

	a3bpose.add_constraint(
		ConstraintOP(
		new DihedralConstraint(
		AtomID( a3bpose.residue( 1 ).atom_index( "H" ), 1 ),
		AtomID( a3bpose.residue( 1 ).atom_index( "N"  ), 1 ),
		AtomID( a3bpose.residue( 1 ).atom_index( "CY"  ), 1 ),
		AtomID( a3bpose.residue( 1 ).atom_index( "OY"  ), 1 ),
		core::scoring::func::CircularHarmonicFuncOP( new core::scoring::func::CircularHarmonicFunc( 3.14159, 0.1 ) )
		) ) );

	a3bpose.dump_pdb( "postpseudobonds.pdb");

	// minimize new pose a bit
	kinematics::MoveMapOP a3blittlemm( new kinematics::MoveMap() );
	a3blittlemm->set_bb( true );
	a3blittlemm->set_chi( true );
	simple_moves::MinMoverOP a3blittlemin( new protocols::simple_moves::MinMover( a3blittlemm, score_fxn_/*cart_*/, "lbfgs_armijo_nonmonotone", 1, true ) );
	a3blittlemin->cartesian( false );
	for ( Size ii = 0; ii <= 10; ++ii ) {
		score_fxn_->set_weight_if_zero( atom_pair_constraint, 0.01*ii );
		score_fxn_->set_weight_if_zero( dihedral_constraint, 0.01*ii );
		score_fxn_->set_weight_if_zero( angle_constraint, 0.01*ii );
		a3blittlemin->apply( a3bpose );
	}
	a3bpose.dump_pdb( "posta3blittlemin.pdb");

	PDBInfoCOP pdb_info = pose.pdb_info();
	for ( Size i = 1; i <= pose.size(); ++i ) {
		char chn = pose.pdb_info()->chain(i);
		//TR << "evaluating residue " << chn  << " " << pose.pdb_info()->number(i) << std::endl;
		if ( chn != hbs_chain_ ) continue;
		oldnums_.push_back( pose.pdb_info()->number(i)+option[a3b_hbs_creator::change_resnums].value() );
	}

	// delete res
	for ( Size i = 1; i <= pose.size(); ++i ) {
		char chn = pose.pdb_info()->chain(i);
		//TR << "evaluating residue " << chn  << " " << pose.pdb_info()->number(i) << std::endl;
		if ( chn != hbs_chain_ ) continue;
		pose.delete_residue_range_slow( i, pose.size() );
		break;
	}

	Size last_of_first_chain = pose.size();
	pose::append_pose_to_pose( pose, a3bpose, true );
	pose::ncbb::initialize_ncbbs( pose );

	// Constraint that we might end up applying in initialize_ncbbs, too
	pose.add_constraint(
		ConstraintOP(
		new DihedralConstraint(
		AtomID( pose.residue( last_of_first_chain + 1 ).atom_index( cyn ), last_of_first_chain + 1 ),
		AtomID( pose.residue( last_of_first_chain + 1 ).atom_index( "CYH"  ), last_of_first_chain + 1 ),
		AtomID( pose.residue( last_of_first_chain + 3 ).atom_index( "CZH"  ), last_of_first_chain + 3 ),
		AtomID( pose.residue( last_of_first_chain + 3 ).atom_index( "N"  ), last_of_first_chain + 3 ),
		core::scoring::func::CircularHarmonicFuncOP( new core::scoring::func::CircularHarmonicFunc( 3.14159, 0.1 ) )
		) ) );

	pose.add_constraint(
		ConstraintOP(
		new DihedralConstraint(
		AtomID( pose.residue( last_of_first_chain + 1 ).atom_index( "H" ), last_of_first_chain + 1 ),
		AtomID( pose.residue( last_of_first_chain + 1 ).atom_index( "N"  ), last_of_first_chain + 1 ),
		AtomID( pose.residue( last_of_first_chain + 1 ).atom_index( "CY"  ), last_of_first_chain + 1 ),
		AtomID( pose.residue( last_of_first_chain + 1 ).atom_index( "OY"  ), last_of_first_chain + 1 ),
		core::scoring::func::CircularHarmonicFuncOP( new core::scoring::func::CircularHarmonicFunc( 3.14159, 0.1 ) )
		) ) );

	// Add the same constraints to hbonds
	add_hbond_and_omega_constraints_starting_at_seqpos( pose, last_of_first_chain+1 );

	AtomID fixed_atom( pose.residue(1).atom_index( "CA" ), 1 );

	for ( Size i = last_of_first_chain+3; i <= pose.size(); ++i ) {
		if ( former_CAs_.find( pose.pdb_info()->number(i) ) == former_CAs_.end() ) continue;

		pose.add_constraint(
			ConstraintOP(
			new CoordinateConstraint(
			AtomID( pose.residue( i ).atom_index( "CA" ), i ),
			fixed_atom,
			former_CAs_[ pose.pdb_info()->number(i) ],
			core::scoring::func::HarmonicFuncOP( new core::scoring::func::HarmonicFunc( 0.0, 0.25 ) )
			) ) );
	}
	pose.dump_pdb( "postcombination.pdb");

	setup_pert_foldtree(pose);


	// TINYMIN
	// create move map for minimization
	kinematics::MoveMapOP littlemm( new kinematics::MoveMap );
	littlemm->set_bb( false );
	littlemm->set_chi( false );
	littlemm->set_jump( 1, true );

	TR << "Creating littlemin minmover " << std::endl;
	simple_moves::MinMoverOP littlemin( new protocols::simple_moves::MinMover( littlemm, score_fxn_/*cart_*/, option[ OptionKeys::run::min_type ].value(), 1, true ) );
	littlemin->cartesian( false );//true );
	// Ramp constraint weights?
	for ( Size ii = 1; ii <= 10; ++ii ) {
		score_fxn_->set_weight( coordinate_constraint, 0.01 * ii );
		littlemin->apply( pose );
	}
	pose.dump_pdb( "postfirstlittlemin.pdb");

	for ( core::Size i = 1; i <= pose.size(); ++i ) {
		if ( pose.pdb_info()->chain(i) != hbs_chain_ ) continue;

		if ( i > last_of_first_chain + 3 ) littlemm->set_bb( i, true );
	}
	littlemm->set_chi( true );

	for ( Size ii = 11; ii <= 20; ++ii ) {
		score_fxn_->set_weight( coordinate_constraint, 0.01 * ii );
		littlemin->apply( pose );
	}
	pose.dump_pdb( "postlittlemin.pdb");

	// SCAN DIHEDRALS
	//scan( pose, score_fxn_, hbs_chain_ );

	// TINYMIN
	TR << "Minimizing away any major clashes." << std::endl;
	TR << "score start: " << pose.energies().total_energy() << std::endl;//( *score_fxn_ )( pose ) << std::endl;
	littlemin->apply( pose );
	TR << "score end: " << ( *score_fxn_ )( pose ) << std::endl;
	TR << (*pose.pdb_info());
	repack( pose );
	TR << "after repack: " << ( *score_fxn_ )( pose ) << std::endl;

	if ( option[ a3b_hbs_creator::final_mc ].value() ) {

		do_mc( pose );

		/*moves::SequenceMoverOP pert_sequence( new moves::SequenceMover() );
		moves::MonteCarloOP pert_mc( new moves::MonteCarlo( pose, *score_fxn_, 0.2 ) );

		kinematics::MoveMapOP pert_pep_mm( new kinematics::MoveMap() );

		// core::Size hbs_position = 1;

		for ( Size i = 1; i <= pose.size(); ++i ) {
		if ( pose.pdb_info()->chain(i) != hbs_chain_ ) continue;
		if ( pose.pdb_info()->number(i) < int(last_of_first_chain + 3) ) continue;

		TR << "setting small movable resid: "<< i<<std::endl;
		//kdrew: commenting out because small mover fails randomly
		pert_pep_mm->set_bb( i, true );
		}

		simple_moves::RandomTorsionMoverOP pert_pep_small( new simple_moves::RandomTorsionMover( pert_pep_mm, 2, 1 ) );

		pert_sequence->add_mover( pert_pep_small );

		//awatkins: add all hbs_pre positions to random small mover
		//TODO: I would PAY for understanding as to why this is so broken.
		//hbs::HbsRandomSmallMoverOP hpm( new hbs::HbsRandomSmallMover ( hbs_position, 2.0));//option[hbs_creator::hbs_length].value(), 2.0 ) );
		//moves::RepeatMoverOP pert_pep_repeat( new moves::RepeatMover( hpm, 1000 ) );
		//pert_sequence->add_mover( pert_pep_repeat );

		moves::TrialMoverOP pert_trial( new moves::TrialMover( pert_sequence, pert_mc ) );

		pert_trial->apply( pose );
		pert_mc->recover_low( pose );
		*/
	}
	pose.dump_pdb( "postmc.pdb");

	if ( option[ a3b_hbs_creator::final_repack ].value() ) {
		repack( pose );
	}
	pose.dump_pdb( "postrepack.pdb");

	if ( option[ a3b_hbs_creator::final_minimize ].value() ) {
		using namespace core::id;
		using namespace core::scoring;
		using namespace core::scoring::constraints;

		// create move map for minimization
		kinematics::MoveMapOP mm( new kinematics::MoveMap() );
		for ( Size i = 1; i <= pose.size(); ++i ) {
			if ( pose.pdb_info()->chain(i) != hbs_chain_ ) continue;
			//if ( Size(pose.pdb_info()->number(i)) < last_of_first_chain + 3 ) continue;
			mm->set_bb( i, true );
		}
		mm->set_chi( true );
		//mm->set_jump( 1, true );
		simple_moves::MinMoverOP minM( new protocols::simple_moves::MinMover( mm, score_fxn_, "lbfgs_armijo_nonmonotone", 0.01, true ) );
		minM->cartesian( false );

		for ( Size ii = 5; ii >= 1; --ii ) {
			score_fxn_->set_weight( coordinate_constraint, 0.04 * (ii-1) );

			minM->apply( pose );
		}

		// Ramp rep back up. Had been halved, so add 1/5 of the weight 5x
		/*Real increment = score_fxn_->get_weight( fa_rep ) * 0.2;
		for ( Size ii = 1; ii <= 5; ++ii ) {
		score_fxn_->set_weight( fa_rep, score_fxn_->get_weight( fa_rep ) + increment );

		minM->apply( pose );
		}*/
	}
}

// this only works for two chains and assumes the protein is first and the peptide is second
// inspired by protocols/docking/DockingProtocol.cc
void
setup_pert_foldtree(
	core::pose::Pose & pose
) {
	using namespace kinematics;

	// get current fold tree
	FoldTree f( pose.fold_tree() );
	f.clear();

	// get the start and end for both chains
	Size pro_start( pose.conformation().chain_begin( 1 ) );
	Size pro_end( pose.conformation().chain_end( 1 ) );
	Size pep_start( pose.conformation().chain_begin( 2 ) );
	Size pep_end( pose.conformation().chain_end( 2 ) );

	// get jump positions based on the center of mass of the chains
	Size dock_jump_pos_pro( core::pose::residue_center_of_mass( pose, pro_start, pro_end ) );
	Size dock_jump_pos_pep( core::pose::residue_center_of_mass( pose, pep_start, pep_end ) );

	// build fold tree
	Size jump_index( f.num_jump() + 1 );
	//-1 is a magic number for PEPTIDE EDGE.  There is a constant defined with the fold tree that should have been used here.
	f.add_edge( pro_start, dock_jump_pos_pro, -1 );
	f.add_edge( dock_jump_pos_pro, pro_end, -1 );
	f.add_edge( pep_start, dock_jump_pos_pep, -1);
	f.add_edge( dock_jump_pos_pep, pep_end, -1 );
	f.add_edge( dock_jump_pos_pro, dock_jump_pos_pep, jump_index );

	// set pose foldtree to foldtree we just created
	f.reorder(1);
	f.check_fold_tree();
	assert( f.check_fold_tree() );

	std::cout << "AFTER: " << f << std::endl;

	pose.fold_tree( f );
}
