// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license.
// (c) The Rosetta software is developed by the contributing members of the
// (c) Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org.
// (c) Questions about this can be addressed to University of Washington UW
// (c) TechTransfer, email: license@u.washington.edu.

/// @file   build_a3b.cc
/// @brief  Miscellany with beta aas
/// @author Andy Watkins (amw579@nyu.edu)

// includes
#include <iostream>
#include <fstream>
#include <string>

#include <devel/init.hh>

#include <core/types.hh>

#include <core/pose/Pose.hh>
#include <core/pose/ncbb/util.hh>

#include <core/import_pose/import_pose.hh>

#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/conformation/Residue.hh>

#include <core/kinematics/MoveMap.hh>

#include <core/id/TorsionID.hh>
#include <core/id/types.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>

#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/FadeFunc.hh>
#include <core/scoring/func/SumFunc.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

// Mover headers
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/PyMolMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/simple_moves/TaskAwareMinMover.hh>
#include <protocols/simple_moves/BackboneMover.fwd.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/RandomTorsionMover.hh>
#include <protocols/simple_moves/a3b_hbs/A3BHbsPatcher.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/rigid/RB_geometry.hh>

#include <utility/vector1.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/excn/Exceptions.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <numeric/random/random.hh>

using namespace protocols;
using namespace protocols::moves;
using namespace protocols::simple_moves;
using namespace protocols::simple_moves::a3b_hbs;

using namespace basic;
using namespace basic::options;
using namespace basic::options::OptionKeys;

class A3BPeptideBuilder : public Mover {

public:

	//default ctor
	A3BPeptideBuilder(): Mover("A3BPeptideBuilder"){}

	//default dtor
	virtual ~A3BPeptideBuilder(){}

	//methods

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const { return "A3BPeptideBuilder"; }

};

typedef utility::pointer::shared_ptr< A3BPeptideBuilder > A3BPeptideBuilderOP;
typedef utility::pointer::shared_ptr< A3BPeptideBuilder const > A3BPeptideBuilderCOP;

int main ( int argc, char* argv[] )
{
	try {
		//option[ chemical::patch_selectors ].push_back( "CTERM_AMIDATION" );

		devel::init(argc, argv);

		A3BPeptideBuilderOP builder( new A3BPeptideBuilder() );
		protocols::jd2::JobDistributor::get_instance()->go( builder );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}

void
A3BPeptideBuilder::apply(
	core::pose::Pose & pose
) {

	using namespace core;
	using namespace utility;
	using namespace scoring;
	using namespace pose;
	using namespace core::chemical;
	using namespace conformation;
	using namespace func;
	using namespace constraints;

	using namespace core::id;
	using namespace core::pack;
	using namespace core::pack::task;

	//first, load the file of residue types to get min energies for.

	//now do initialization stuff.
	TaskFactoryOP task_factory = TaskFactoryOP( new TaskFactory );
	task_factory->push_back( operation::TaskOperationCOP( new operation::InitializeFromCommandline ) );
	//need these to keep pack_rotamers from redesigning the residue.
	operation::RestrictToRepackingOP rtrop = operation::RestrictToRepackingOP( new operation::RestrictToRepacking );
	task_factory->push_back( rtrop );

	ScoreFunctionOP scorefxn = get_score_function();

	//Get the residue set we are drawing from.
	core::chemical::ResidueTypeSetCOP residue_set_cap = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

	pose.clear();

	ResidueType const & restype_first = residue_set_cap->name_map( "PHE:NtermProteinFull" );
	ResidueType const & internal_B3A = residue_set_cap->name_map( "B3L" );
	ResidueType const & internal_ALA = residue_set_cap->name_map( "ALA" );
	ResidueType const & internal_GLY = residue_set_cap->name_map( "GLY" );
	ResidueType const & internal_ASN = residue_set_cap->name_map( "ASN" );
	//ResidueType const & restype_last = residue_set_cap->name_map( "ASN:MethylatedCtermProteinFull" );
	ResidueType const & restype_last = residue_set_cap->name_map( "ALA:MethylatedCtermProteinFull" );
	Residue res_first( restype_first, true );
	Residue res_int_B3A( internal_B3A, true );
	Residue res_int_ALA( internal_ALA, true );
	Residue res_int_GLY( internal_GLY, true );
	Residue res_int_ASN( internal_ASN, true );
	Residue res_last( restype_last, true );
	pose.append_residue_by_jump( res_first, 1 );
	pose.append_residue_by_bond( res_int_ALA, true );
	pose.append_residue_by_bond( res_int_GLY, true );
	pose.append_residue_by_bond( res_int_B3A, true );
	pose.append_residue_by_bond( res_int_ALA, true );
	pose.append_residue_by_bond( res_int_ALA, true );
	pose.append_residue_by_bond( res_int_ALA, true );
	pose.append_residue_by_bond( res_int_B3A, true );
	pose.append_residue_by_bond( res_int_ALA, true );
	pose.append_residue_by_bond( res_int_ALA, true );
	pose.append_residue_by_bond( res_int_ALA, true );
	pose.append_residue_by_bond( res_int_B3A, true );
	pose.append_residue_by_bond( res_int_ALA, true );
	pose.append_residue_by_bond( res_int_ALA, true );
	pose.append_residue_by_bond( res_int_ALA, true );
	pose.append_residue_by_bond( res_int_ASN, true );
	pose.append_residue_by_bond( res_last, true );

	// create movemap for peptide
	kinematics::MoveMapOP pert_alpha_mm( new kinematics::MoveMap() );
	kinematics::MoveMapOP pert_beta_mm( new kinematics::MoveMap() );

	for ( Size i = 1; i <= pose.n_residue(); ++i ) {

		if ( pose.residue(i).type().is_beta_aa() ) { //( i == 4 || i == 8 || i == 12 ) {

			id::TorsionID bb1( i, id::BB, 1 ); //phi
			id::TorsionID bb2( i, id::BB, 2 ); //theta
			id::TorsionID bb3( i, id::BB, 3 ); //psi
			id::TorsionID bb4( i, id::BB, 4 ); //omg

			pert_beta_mm->set_bb( i, true );
			pert_alpha_mm->set_bb( i, false );

			pert_beta_mm->set( bb4, false );

			pose.set_torsion( bb1, -104 );
			pose.set_torsion( bb2, 64 );
			pose.set_torsion( bb3, -116 );
			pose.set_torsion( bb4, 180 );
		} else {

			pert_alpha_mm->set_bb( i, true );
			pert_beta_mm->set_bb( i, false );

			pose.set_phi( i, -60);
			pose.set_omega( i, 180);
			pose.set_psi( i, -40);
		}
	}

	/*if ( pose.residue(pose.n_residue()).type().is_beta_aa() ) { //( i == 4 || i == 8 || i == 12 ) {

	id::TorsionID bb1( pose.n_residue(), id::BB, 1 ); //phi
	id::TorsionID bb2( pose.n_residue(), id::BB, 2 ); //theta
	id::TorsionID bb3( pose.n_residue(), id::BB, 3 ); //psi
	id::TorsionID bb4( pose.n_residue(), id::BB, 4 ); //omg

	pert_beta_mm->set_bb( pose.n_residue(), true );
	pert_alpha_mm->set_bb( pose.n_residue(), false );

	pert_beta_mm->set( bb4, false );

	pose.set_torsion( bb1, -104 );
	pose.set_torsion( bb2, 64 );
	pose.set_torsion( bb3, -116 );
	pose.set_torsion( bb4, 180 );
	} else {

	pert_alpha_mm->set_bb( pose.n_residue(), true );
	pert_beta_mm->set_bb( pose.n_residue(), false );
	id::TorsionID bb1( pose.n_residue(), id::BB, 1 ); //phi
	id::TorsionID bb2( pose.n_residue(), id::BB, 2 ); //psi
	id::TorsionID bb3( pose.n_residue(), id::BB, 3 ); //omg

	pose.set_torsion( bb1, -60 );
	pose.set_torsion( bb2, -40 );
	pose.set_torsion( bb3, 180  );
	}*/


	if ( scorefxn->has_zero_weight( core::scoring::atom_pair_constraint ) ) {
		scorefxn->set_weight( core::scoring::atom_pair_constraint, 5.0 );
	}

	HarmonicFuncOP hbond_func( new HarmonicFunc( 2.0, 0.2 ) );
	//FadeFuncOP hbond_func2( new FadeFunc( 0, 1.6, 0.1, 1000 ) );
	//SumFuncOP hbond_func;
	//hbond_func->add_func( hbond_func1 );
	//hbond_func->add_func( hbond_func2 );
	for ( Size i = 1; i <= pose.n_residue() - 4; ++i ) {
		AtomID aidO( pose.residue( i ).atom_index( "O" ), i );
		AtomID aidH( pose.residue( i+4 ).atom_index( "H" ), i+4 );

		ConstraintCOP hbond( new AtomPairConstraint( aidO, aidH, hbond_func ) );
		pose.add_constraint( hbond );
	}
	AtomID aidO( pose.residue( pose.n_residue() - 3 ).atom_index( "O" ), pose.n_residue() - 3 );
	AtomID aidH( pose.residue( pose.n_residue() ).atom_index( "HM" ), pose.n_residue() );

	ConstraintCOP hbond( new AtomPairConstraint( aidO, aidH, hbond_func ) );
	pose.add_constraint( hbond );

	A3BHbsPatcherOP hbs_patcher (new A3BHbsPatcher( 1 ) );
	hbs_patcher->apply( pose );

	pose.conformation().detect_bonds();
	pose.conformation().detect_pseudobonds();
	for ( core::Size i=1; i<=pose.total_residue(); ++i ) {
		pose.conformation().update_polymeric_connection(i);
	}

	kinematics::MoveMapOP desn_mm( new kinematics::MoveMap );
	desn_mm->set_bb( true );
	for ( Size i = 1; i <= pose.n_residue(); ++i ) {
		id::TorsionID bb3( i, id::BB, pose.residue(i).type().is_beta_aa() ? 4 : 3 ); //omg
		desn_mm->set( bb3, false ); // needed for mm_std, which has no omega tether
	}
	desn_mm->set_chi( true );
	protocols::simple_moves::MinMoverOP desn_min( new simple_moves::MinMover( desn_mm, scorefxn, "lbfgs_armijo_nonmonotone", 0.0001, true ) );

	std::cout << "Dump initial" << std::endl;
	pose.dump_scored_pdb( "B3A_initial.pdb", *scorefxn);
	desn_min->apply( pose );
	std::cout << "Dump initial min" << std::endl;
	pose.dump_scored_pdb( "B3A_min.pdb", *scorefxn);

	// create random torsion mover
	kinematics::MoveMapOP pert_pep_mm( new kinematics::MoveMap() );
	pert_pep_mm->set_bb( true );
	moves::MonteCarloOP mc( new moves::MonteCarlo( pose, *scorefxn, 1.0 ) );

	// create small and shear movers
	simple_moves::SmallMoverOP pert_pep_alpha( new simple_moves::SmallMover( pert_alpha_mm, .8, 1 ) );
	pert_pep_alpha->angle_max( 'H', 2 );
	pert_pep_alpha->angle_max( 'L', 2 );
	pert_pep_alpha->angle_max( 'E', 2 );

	simple_moves::RandomTorsionMoverOP pert_pep_beta( new simple_moves::RandomTorsionMover( pert_beta_mm, .8, 1 ) );


	moves::RandomMoverOP pert_pep_random( new moves::RandomMover() );
	pert_pep_random->add_mover( pert_pep_alpha, .75 );
	pert_pep_random->add_mover( pert_pep_beta, .25 );
	moves::RepeatMoverOP pert_repeat( new moves::RepeatMover( pert_pep_random, 100 ) );

	moves::TrialMoverOP pert_trial( new moves::TrialMover( pert_repeat, mc ) );

	//definitely want sidechain minimization here
	using core::pack::task::operation::TaskOperationCOP;
	TaskFactoryOP desn_tf( new TaskFactory() );
	desn_tf->push_back( TaskOperationCOP( new core::pack::task::operation::InitializeFromCommandline ) );

	using protocols::simple_moves::TaskAwareMinMoverOP;
	using protocols::simple_moves::TaskAwareMinMover;
	TaskAwareMinMoverOP desn_ta_min( new TaskAwareMinMover( desn_min, desn_tf ) );

	protocols::jd2::JobOP curr_job( protocols::jd2::JobDistributor::get_instance()->current_job() );

	for ( Size k = 1; k <= 10; ++k ) {

		mc->reset(pose);

		// pert loop
		for ( Size j = 1; j <= 200; ++j ) {
			std::cout <<  "PERTURB: " << k << " / "  << j << std::endl;
			pert_trial->apply( pose );
		}
		mc->recover_low( pose );

		PackerTaskOP task( TaskFactory::create_packer_task( pose ) );
		for ( Size i = 1; i <= pose.n_residue(); i++ ) {
			task->nonconst_residue_task(i).restrict_to_repacking();
			task->nonconst_residue_task(i).initialize_from_command_line();
		}
		// create a pack rotamers mover
		simple_moves::PackRotamersMoverOP desn_pr( new simple_moves::PackRotamersMover(scorefxn, task) );
		moves::SequenceMoverOP desn_sequence( new moves::SequenceMover() );
		desn_sequence->add_mover( desn_pr );
		desn_sequence->add_mover( desn_ta_min );
		desn_sequence->apply( pose );

		std::cout << "pre mc->boltzmann" << std::endl;
		mc->show_state();
		mc->boltzmann( pose );
		std::cout << "post mc->boltzmann" << std::endl;
		mc->show_state();
	}

	mc->recover_low( pose );

	std::cout <<  "Ending main loop..." << std::endl;

	std::cout <<  "Checking pose energy..." << std::endl;
	std::cout << "Final min" << std::endl;
	desn_min->apply( pose );

}
