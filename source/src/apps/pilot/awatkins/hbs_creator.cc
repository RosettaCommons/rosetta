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
#include <core/import_pose/import_pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/util.hh>
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
#include <protocols/ncbb/util.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/PyMOLMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/BackboneMover.hh>
//#include <protocols/simple_moves/hbs/HbsRandomSmallMover.hh>
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
using namespace protocols::ncbb;
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
static basic::Tracer TR("HBS_Creator");

// application specific options
namespace hbs_creator {
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
		HbsCreatorMoverOP HC_mover( new HbsCreatorMover );

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
	scoring::ScoreFunctionOP score_fxn( scoring::get_score_function() );
	scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn(*score_fxn);

	scoring::constraints::add_fa_constraints_from_cmdline_to_pose(pose);

	if ( score_fxn->has_zero_weight( atom_pair_constraint ) ) {
		score_fxn->set_weight( atom_pair_constraint, 0.1 );
	}

	if ( score_fxn->has_zero_weight( dihedral_constraint ) ) {
		score_fxn->set_weight( dihedral_constraint, 0.1 );
	}

	if ( score_fxn->has_zero_weight( angle_constraint ) ) {
		score_fxn->set_weight( angle_constraint, 0.1 );
	}

	chemical::ResidueTypeSetCOP restype_set = chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

	char hbs_chain = option[hbs_creator::hbs_chain].value()[0];
	core::Size final_res = option[hbs_creator::hbs_final_res].value();
	core::pose::PDBInfoCOP pdb_info( pose.pdb_info() );

	// shift the jump
	//core::kinematics::FoldTree f = pose.fold_tree();
	//f.slide_jump( 1, 1, pose.pdb_info()->pdb2pose( hbs_chain, final_res+1 ) );
	//pose.fold_tree( f );

	Size patch_ros_num = 0;

	for ( Size i = 1; i <= pose.size(); ++i ) {
		char chn = pdb_info->chain(i);
		core::Size pdb_res_num = pdb_info->number(i);
		TR << "evaluating residue " << chn  << " " << pdb_info->number(i) << std::endl;

		if ( chn != hbs_chain ) continue;
		// correct chain to be truncated and prepped

		// hbs pre is the smallest number of what we want to preserve
		while ( pdb_res_num < final_res ) {
			TR << "deleting residue " << pdb_res_num  << " which was " << core::chemical::oneletter_code_from_aa(pose.aa(i)) << std::endl;
			pose.delete_polymer_residue(i);
			pdb_res_num = pdb_info->number(i);
		}

		if ( pdb_res_num > final_res + option[hbs_creator::hbs_length].value() ) {
			//TR << "deleting residue " << pdb_res_num << std::endl;
			while ( chn == hbs_chain && i <= pose.size() ) {
				chn = pdb_info->chain(i);
				pose.delete_polymer_residue(i);
			}
		}

		if ( pdb_res_num == final_res ) patch_ros_num = i;
	}

	hbs::HbsPatcherOP hbs_patcher( new hbs::HbsPatcher( patch_ros_num ) );
	hbs_patcher->apply( pose );


	setup_pert_foldtree(pose);

	// presently the final residue in the pose is the terminal residue of the hbs
	// replace with terminal variant
	// AMW: This means that you have to clear other possible Cterm variant assignments
	// from this residue.
	conformation::Residue term( restype_set->get_residue_type_with_variant_added(
		restype_set->get_residue_type_with_variant_removed(
		pose.residue(pose.size()).type(),
		chemical::UPPER_TERMINUS_VARIANT ),
		chemical::METHYLATED_CTERMINUS_VARIANT), true );

	term.set_all_chi(pose.residue(pose.size()).chi());
	//replace_res_post.mainchain_torsions(pose.residue(oop_post_pos_).mainchain_torsions());

	pose.replace_residue( pose.size(), term, true );
	conformation::idealize_position( pose.size(), pose.conformation() );


	//pose.set_phi(i+3, -40);
	//pose.set_psi(i+3, -58);

	//pose.dump_pdb( "postdihedrals.pdb");

	pose.conformation().detect_bonds();
	pose.conformation().detect_pseudobonds();
	for ( core::Size i = 1; i <= pose.size(); ++i ) {
		pose.conformation().update_polymeric_connection(i);
	}

	pose.dump_pdb( "postpseudobonds.pdb");

	// TINYMIN
	// create move map for minimization
	kinematics::MoveMapOP littlemm( new kinematics::MoveMap() );
	littlemm->set_bb( false );
	littlemm->set_chi( true );
	//mm->set_jump( 1, true );

	// create minimization mover
	simple_moves::MinMoverOP littlemin( new protocols::simple_moves::MinMover( littlemm, score_fxn, option[ OptionKeys::run::min_type ].value(), 1, true ) );

	littlemin->apply( pose );

	if ( option[ hbs_creator::final_mc ].value() ) {
		moves::SequenceMoverOP pert_sequence( new moves::SequenceMover() );
		moves::MonteCarloOP pert_mc( new moves::MonteCarlo( pose, *score_fxn, 0.2 ) );

		kinematics::MoveMapOP pert_pep_mm( new kinematics::MoveMap() );

		// core::Size hbs_position = 1;

		for ( Size i = 1; i <= pose.size(); ++i ) {
			if ( pdb_info->chain(i) == hbs_chain ) { //i != ( ( unsigned ) ( option[ hbs_creator::hbs_final_res ].value() ) ) ) {
				if ( pose.residue_type( i ).is_l_aa() ) {
					TR << "setting small movable resid: "<< i<<std::endl;
					//kdrew: commenting out because small mover fails randomly
					pert_pep_mm->set_bb( i );
				}
			}
		}

		simple_moves::SmallMoverOP pert_pep_small( new simple_moves::SmallMover( pert_pep_mm, 0.2, 1 ) );
		pert_pep_small->angle_max( 'H', 2.0 );
		pert_pep_small->angle_max( 'L', 2.0 );
		pert_pep_small->angle_max( 'E', 2.0 );

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
	pose.dump_pdb( "postmc.pdb");

	if ( option[ hbs_creator::final_repack ].value() ) {

		// create a task factory and task operations
		TaskFactoryOP tf(new TaskFactory());
		tf->push_back( operation::TaskOperationCOP( new operation::InitializeFromCommandline ) );

		using namespace basic::resource_manager;
		if ( ResourceManager::get_instance()->has_option( packing::resfile ) ||  option[ packing::resfile ].user() ) {
			operation::ReadResfileOP rrop( new operation::ReadResfile() );
			rrop->default_filename();
			tf->push_back( rrop );
		} else {
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
	pose.dump_pdb( "postrepack.pdb");

	if ( option[ hbs_creator::final_minimize ].value() ) {
		using namespace core::id;
		using namespace core::scoring;
		using namespace core::scoring::constraints;

		// create move map for minimization
		kinematics::MoveMapOP mm( new kinematics::MoveMap() );
		mm->set_chi( true );
		//mm->set_jump( 1, true );
		for ( Size i = 1; i <= pose.size(); ++i ) {
			if ( pose.pdb_info()->chain(i) == hbs_chain ) {
				TR << "Can move " << i << " bb " << std::endl;
				mm->set_bb( i, true);
			} else {
				mm->set_bb( i, false );
			}
		}

		TR << *mm << std::endl;
		score_fxn->set_weight( atom_pair_constraint, 0.5 );
		score_fxn->set_weight( dihedral_constraint, 0.5 );
		score_fxn->set_weight( angle_constraint, 0.5 );

		// create minimization mover
		simple_moves::MinMoverOP minM( new protocols::simple_moves::MinMover( mm, score_fxn, "lbfgs_armijo_nonmonotone", 0.01, true ) );
		//minM->cartesian( true );
		minM->apply( pose );

		score_fxn->set_weight( atom_pair_constraint, 1 );
		score_fxn->set_weight( dihedral_constraint, 1 );
		score_fxn->set_weight( angle_constraint, 1 );
		minM->apply( pose );

	}
}

// this only works for two chains and assumes the protein is first and the peptide is second
// inspired by protocols/docking/DockingProtocol.cc
void
setup_pert_foldtree( core::pose::Pose & pose ) {
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
