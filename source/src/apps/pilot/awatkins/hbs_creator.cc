// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// awatkins: based heavily on kdrew/oop_creator.cc

// Project Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueTypeSet.hh>
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
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

// Mover headers
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/PyMolMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/hbs/HbsRandomSmallMover.hh>
#include <protocols/simple_moves/hbs/HbsPatcher.hh>
#include <protocols/simple_moves/chiral/ChiralMover.hh>
#include <protocols/rigid/RB_geometry.hh>

#include <numeric/conversions.hh>

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
using namespace protocols::simple_moves::hbs;
using namespace protocols::simple_moves::chiral;
using namespace core::pack::task;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::id;
using basic::T;
using basic::Error;
using basic::Warning;
using utility::file::FileName;

//kdrew: this app adds hbs patches to the given pdb strucure

// tracer - used to replace cout
static thread_local basic::Tracer TR( "HBS_Creator" );

void
setup_pert_foldtree(
					core::pose::Pose & pose
					);

// application specific options
namespace hbs_creator{
	// pert options
	StringOptionKey const hbs_chain ( "hbs_creator::hbs_chain" );
	IntegerOptionKey const hbs_final_res ( "hbs_creator::hbs_final_res" );
	IntegerOptionKey const hbs_length ( "hbs_creator::hbs_length" );
	BooleanOptionKey const final_repack( "hbs_creator::final_repack" );
	BooleanOptionKey const final_minimize( "hbs_creator::final_minimize" );
	BooleanOptionKey const final_mc ( "hbs_creator::final_mc" );
	// BooleanOptionKey const correct_hbs_dihedrals ( "hbs_creator::correct_hbs_dihedrals" ); to be implemented if possible

}

class HbsCreatorMover : public Mover {

	public:

		//default ctor
		HbsCreatorMover(): Mover("HbsCreatorMover"){}

		//default dtor
		virtual ~HbsCreatorMover(){}

		virtual void apply( core::pose::Pose & pose );
		virtual std::string get_name() const { return "HbsCreatorMover"; }

};

typedef utility::pointer::shared_ptr< HbsCreatorMover > HbsCreatorMoverOP;
typedef utility::pointer::shared_ptr< HbsCreatorMover const > HbsCreatorMoverCOP;


int
main( int argc, char* argv[] )
{
	try {
	utility::vector1< core::Size > empty_vector(0);

	option.add( hbs_creator::hbs_chain, "Chain from PDB to be mimicked. Default 'A'. Use letters." ).def("A");
	option.add( hbs_creator::hbs_final_res, "Residue number of the final residue for mimicry. Default 1." ).def(1);
	option.add( hbs_creator::hbs_length, "Number of residues to mimic. Default 12." ).def(12);
	option.add( hbs_creator::final_repack, "Do a final repack. Default false" ).def(false);
	option.add( hbs_creator::final_minimize, "Do a final minimization. Default false" ).def(false);
	option.add( hbs_creator::final_mc, "Do a final monte carlo on hbs. Default false" ).def(false);
	//option.add( hbs_creator::correct_hbs_dihedrals, "Correct hbs dihedral to low energy well. Default false" ).def(false);

	// init command line options
	//you MUST HAVE THIS CALL near the top of your main function, or your code will crash when you first access the command line options
	devel::init(argc, argv);

	//create mover instance
	HbsCreatorMoverOP HC_mover( new HbsCreatorMover() );

	//call job distributor
	protocols::jd2::JobDistributor::get_instance()->go( HC_mover );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;

}//main

void
HbsCreatorMover::apply(
	core::pose::Pose & pose
)
{

	// create score function
	//kdrew: old standard scoring function, using MM scoring function now because of NCAAs
	//scoring::ScoreFunctionOP score_fxn( ScoreFunctionFactory::create_score_function( scoring::STANDARD_WTS, scoring::SCORE12_PATCH ) );
	scoring::ScoreFunctionOP score_fxn( ScoreFunctionFactory::create_score_function( scoring::MM_STD_WTS) );
	scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn(*score_fxn);

	scoring::constraints::add_fa_constraints_from_cmdline_to_pose(pose);

	 score_fxn->set_weight( core::scoring::atom_pair_constraint, 1.0 );
	/*

	 TODO:

	 We probably aren't going to need this. We are probably going to:
	 1) Delete all residues in the chain hbs_creator::hbs_chain that are outside of hbs_creator::hbs_final_res to that plus hbs_creator::hbs_length
	 2) There aren't two puckers. So we'll probably just push back the three residues we care about to all_positions.
	*/
	utility::vector1< core::Size > all_positions;
	char hbs_chain = option[hbs_creator::hbs_chain].value()[0];
	core::Size final_res = option[hbs_creator::hbs_final_res].value();
	core::pose::PDBInfoCOP pdb_info( pose.pdb_info() );

	for(Size i=1; i<=pose.total_residue(); ++i) {
		char chn = pdb_info->chain(i);
		if (chn == hbs_chain) { // correct chain to be truncated and prepped

			core::Size pdb_res_num = pdb_info->number(i);

			// hbs pre is the smallest number of what we want to preserve
			if (pdb_res_num < final_res) {
				TR << "deleting residue " << pdb_res_num  << " which was " << core::chemical::oneletter_code_from_aa(pose.aa(i)) << std::endl;

				while (pdb_res_num < final_res) {
          pose.delete_polymer_residue(i);
				  pdb_res_num = pdb_info->number(i);
        }

        setup_pert_foldtree(pose);


			}
			else if (pdb_res_num == final_res) { // also applies the post patch
				hbs::HbsPatcherOP hbs_patcher( new hbs::HbsPatcher( i ) );
				//pose.dump_pdb( "prepatch.pdb");
				hbs_patcher->apply( pose );
				//pose.dump_pdb( "postpatch.pdb");
			}
			else if (pdb_res_num > final_res + option[hbs_creator::hbs_length].value()) {
				TR << "deleting residue " << pdb_res_num << std::endl;

				while (pdb_res_num > final_res) {
          pose.delete_polymer_residue(i);
				  pdb_res_num = pdb_info->number(i);
        }

				setup_pert_foldtree(pose);

			}

			if (pdb_res_num >= final_res && pdb_res_num <= final_res+2) {
				all_positions.push_back(i);

				if (pdb_res_num == final_res) {
					pose.set_phi(i, -78);
					pose.set_psi(i, -48);
				}
				else if (pdb_res_num == final_res + 1) {
					pose.set_phi(i+1, -72);
					pose.set_psi(i+1, -47);
				}
				//else if (pdb_res_num == final_res + 2) {
				//	pose.set_phi(i+2, -55);
				//	pose.set_psi(i+2, -45);
				//}

				/*
				pose.set_psi( i, -57) ;
				pose.set_phi( i, -48) ;
				*/
			}
		}
	}

	//pose.set_phi(i+3, -40);
	//pose.set_psi(i+3, -58);

	//pose.dump_pdb( "postdihedrals.pdb");

	//kdrew: since we added new connection types (i.e. hbs atoms) above, need to reset connections
	pose.conformation().detect_bonds();
	pose.conformation().detect_pseudobonds();
	for(core::Size i=1; i<=pose.total_residue(); ++i){
		pose.conformation().update_polymeric_connection(i);
	}

	//pose.dump_pdb( "postpseudobonds.pdb");

	//kdrew: monte carlo phi/psi of hbs to find low energy
	if( option[ hbs_creator::final_mc ].value() )
	{
		moves::SequenceMoverOP pert_sequence( new moves::SequenceMover() );
		moves::MonteCarloOP pert_mc( new moves::MonteCarlo( pose, *score_fxn, 0.2 ) );

		kinematics::MoveMapOP pert_pep_mm( new kinematics::MoveMap() );
		simple_moves::SmallMoverOP pert_pep_small( new simple_moves::SmallMover( pert_pep_mm, 0.2, 1 ) );
		pert_pep_small->angle_max( 'H', 2.0 );
		pert_pep_small->angle_max( 'L', 2.0 );
		pert_pep_small->angle_max( 'E', 2.0 );

		// core::Size hbs_position = 1;

		//kdrew: load all hbs-affected positions into vector and make all other positions movable by small mover
		for( Size i = 1; i <= pose.total_residue(); ++i )
		{
			//TR << "resid: " << i << " is last in the HBS macrocycle." << std::endl;
			//if( i < option[hbs_creator::hbs_final_res] || i>option[hbs_creator::hbs_final_res+2] ) // maybe 3 someday?
			if (i != ((unsigned) (option[hbs_creator::hbs_final_res].value())))
			{
				if( is_l_chiral( pose.residue_type( i ) ) )
				{
					TR << "setting small movable resid: "<< i<<std::endl;
					//kdrew: commenting out because small mover fails randomly
					//pert_pep_mm->set_bb( i );
				}
			}
			//else
			//{ hbs_position = i; }
		}

		pert_sequence->add_mover( pert_pep_small );

		//awatkins: add all hbs_pre positions to random small mover
		//TODO: I would PAY for understanding as to why this is so broken.
		//hbs::HbsRandomSmallMoverOP hpm( new hbs::HbsRandomSmallMover ( hbs_position, 2.0));//option[hbs_creator::hbs_length].value(), 2.0 ) );
		//moves::RepeatMoverOP pert_pep_repeat( new moves::RepeatMover( hpm, 1000 ) );
		//pert_sequence->add_mover( pert_pep_repeat );

		moves::TrialMoverOP pert_trial( new moves::TrialMover( pert_sequence, pert_mc ) );

		pert_trial->apply( pose );
    	pert_mc->recover_low( pose );

	}

	if( option[ hbs_creator::final_repack ].value() )
	{

		// create a task factory and task operations
		using core::pack::task::operation::TaskOperationCOP;
		TaskFactoryOP tf( new TaskFactory() );
		tf->push_back( TaskOperationCOP( new core::pack::task::operation::InitializeFromCommandline ) );

		using namespace basic::resource_manager;
		if ( ResourceManager::get_instance()->has_option( packing::resfile ) ||  option[ packing::resfile ].user() )
		{
			operation::ReadResfileOP rrop( new operation::ReadResfile() );
			rrop->default_filename();
			tf->push_back( rrop );
		}
		else
		{
			//kdrew: do not do design, makes NATAA if res file is not specified
			operation::RestrictToRepackingOP rtrp( new operation::RestrictToRepacking() );
			tf->push_back( rtrp );
		}


		// create a pack rotamers mover
		simple_moves::PackRotamersMoverOP packer( new protocols::simple_moves::PackRotamersMover() );
		packer->task_factory( tf );
		packer->score_function( score_fxn );

		packer->apply(pose);
	}


	if( option[ hbs_creator::final_minimize ].value() )
	{
		using namespace core::id;
		using namespace core::scoring;
		using namespace core::scoring::constraints;

		//kdrew: add constraints to omega angle, (this problem might have been fixed and these constraints are unnecessary)
		//iwatkins: this is general stuff, not hbs or oop specific
		for( Size i = 1; i < pose.conformation().chain_end( 1 ); ++i )
		{
			id::AtomID id1,id2,id3,id4;
			core::id::TorsionID torsion_id = TorsionID( i, id::BB, 3 ); //kdrew: 3 is omega angle

			//kdrew: put constraint on omega angle
			pose.conformation().get_torsion_angle_atom_ids( torsion_id, id1, id2, id3, id4 );

			Real torsion_value( pose.torsion( torsion_id ) );

			core::scoring::func::CircularHarmonicFuncOP circularharm_func( new core::scoring::func::CircularHarmonicFunc( numeric::conversions::radians( torsion_value ), numeric::conversions::radians( 10.0 ) ) );

			ConstraintCOP dihedral1( ConstraintOP( new DihedralConstraint( id1, id2, id3, id4, circularharm_func ) ) );

			pose.add_constraint( dihedral1 );
		}
		//kdrew: if constraint weight is not set on commandline or elsewhere, set to 1.0
		if( score_fxn->has_zero_weight( dihedral_constraint ) )
		{
        	score_fxn->set_weight( dihedral_constraint, 1.0 );
		}
		if( score_fxn->has_zero_weight( atom_pair_constraint ) )
		{
        	score_fxn->set_weight( atom_pair_constraint, 1.0 );
		}


		// create move map for minimization
		kinematics::MoveMapOP mm( new kinematics::MoveMap() );
		mm->set_bb( true );
		mm->set_chi( true );
		mm->set_jump( 1, true );

		// create minimization mover
		simple_moves::MinMoverOP minM( new protocols::simple_moves::MinMover( mm, score_fxn, option[ OptionKeys::run::min_type ].value(), 0.01,	true ) );

	//kdrew: only turn on pymol observer in debug mode
	//#ifndef NDEBUG
	 //   protocols::moves::PyMolObserverOP pymover = protocols::moves::AddPyMolObserver(pose);
	//#endif

		//kdrew: minimizer not working after appending/prepending residues, not sure why
		// final min (okay to use ta min here)
		minM->apply( pose );
	}
}

// this only works for two chains and assumes the protein is first and the peptide is second
// inspired by protocols/docking/DockingProtocol.cc
void
setup_pert_foldtree(
												core::pose::Pose & pose
												)
{
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
