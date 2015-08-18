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
#include <protocols/moves/PyMolMover.hh>
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
using basic::T;
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
namespace a3b_hbs_creator{
	// pert options
	StringOptionKey const hbs_chain ( "a3b_hbs_creator::hbs_chain" );
	IntegerOptionKey const hbs_final_res ( "a3b_hbs_creator::hbs_final_res" );
	IntegerOptionKey const hbs_length ( "a3b_hbs_creator::hbs_length" );
	BooleanOptionKey const final_repack( "a3b_hbs_creator::final_repack" );
	BooleanOptionKey const final_minimize( "a3b_hbs_creator::final_minimize" );
	BooleanOptionKey const final_mc ( "a3b_hbs_creator::final_mc" );
}

class A3BHbsCreatorMover : public Mover {

	public:

		//default ctor
		A3BHbsCreatorMover(): Mover("A3BHbsCreatorMover"){}

		//default dtor
		virtual ~A3BHbsCreatorMover(){}

		virtual void apply( core::pose::Pose & pose );
		virtual std::string get_name() const { return "A3BHbsCreatorMover"; }

};

typedef utility::pointer::shared_ptr< A3BHbsCreatorMover > A3BHbsCreatorMoverOP;
typedef utility::pointer::shared_ptr< A3BHbsCreatorMover const > A3BHbsCreatorMoverCOP;

std::string alpha_to_beta( std::string alpha ) {
	std::string beta = alpha;
	beta.erase( 0, beta.find( ":" ) );
	std::string prefix = "XXX";
	if ( alpha == "ALA") {
		prefix = "B3A";
	} else if ( alpha == "CYS") {
		prefix = "B3C";
	} else if ( alpha == "ASP") {
		prefix = "B3D";
	} else if ( alpha == "GLU") {
		prefix = "B3E";
	} else if ( alpha == "PHE") {
		prefix = "B3F";
	} else if ( alpha == "GLY") {
		prefix = "B3G";
	} else if ( alpha == "HIS") {
		prefix = "B3H";
	} else if ( alpha == "HIS_D") {
		prefix = "B3H";
	} else if ( alpha == "ILE") {
		prefix = "B3I";
	} else if ( alpha == "LYS") {
		prefix = "B3K";
	} else if ( alpha == "LEU") {
		prefix = "B3L";
	} else if ( alpha == "MET") {
		prefix = "B3M";
	} else if ( alpha == "ASN") {
		prefix = "B3N";
	} else if ( alpha == "PRO") {
		prefix = "B3P";
	} else if ( alpha == "GLN") {
		prefix = "B3Q";
	} else if ( alpha == "ARG") {
		prefix = "B3R";
	} else if ( alpha == "SER") {
		prefix = "B3S";
	} else if ( alpha == "THR") {
		prefix = "B3T";
	} else if ( alpha == "VAL") {
		prefix = "B3V";
	} else if ( alpha == "TRP") {
		prefix = "B3W";
	} else if ( alpha == "TYR") {
		prefix = "B3Y";
	}
	return prefix + beta;
}

void repack(
	core::pose::Pose & pose,
	core::scoring::ScoreFunctionOP score_fxn
) {
	// create a task factory and task operations
	TaskFactoryOP tf(new TaskFactory());
	tf->push_back( operation::TaskOperationCOP( new operation::InitializeFromCommandline ) );
	
	using namespace basic::resource_manager;
	if ( ResourceManager::get_instance()->has_option( packing::resfile ) ||  option[ packing::resfile ].user() ) {
		operation::ReadResfileOP rrop( new operation::ReadResfile() );
		rrop->default_filename();
		tf->push_back( rrop );
	}
	else {
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

void scan(
	core::pose::Pose & pose,
	core::scoring::ScoreFunctionOP score_fxn,
	char hbs_chn
) {
	core::pose::Pose test_pose = pose;
	core::pose::PDBInfoCOP pdb_info( test_pose.pdb_info() );
	
	
	
	kinematics::MoveMapOP littlemm( new kinematics::MoveMap() );
	littlemm->set_bb( false );
	for ( core::Size i = 1; i <= test_pose.total_residue(); ++i ) {
		if ( pdb_info->chain(i) != hbs_chn ) {
			continue;
		}
		
		littlemm->set_bb( true );
	}
	
	littlemm->set_chi( true );
	littlemm->set_jump( 1, true );
	simple_moves::MinMoverOP littlemin( new protocols::simple_moves::MinMover( littlemm, score_fxn, option[ OptionKeys::run::min_type ].value(), 1, true ) );
	
	
	
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
		for ( Size i = 1; i <= test_pose.total_residue(); ++i ) {
			if ( pdb_info->chain(i) != hbs_chn ) {
				continue;
			}
			if ( pdb_info->chain(i) == hbs_chn && first) {
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
				
			} else{
			if ( test_pose.residue(i).type().is_beta_aa() ) {
				core::id::TorsionID phi( i, id::BB, 1);
				core::id::TorsionID tht( i, id::BB, 2);
				
				test_pose.conformation().set_torsion( phi, betaphi);
				test_pose.conformation().set_torsion( tht, betatht);
				if ( i != pose.total_residue() ) {
					core::id::TorsionID psi( i, id::BB, 3);
					test_pose.conformation().set_torsion( psi, betapsi);
				}
			} else {
				core::id::TorsionID phi( i, id::BB, 1);
				
				test_pose.conformation().set_torsion( phi, alphaphi);
				if ( i != pose.total_residue() ) {
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

		
		Real score = ( *score_fxn )( test_pose );
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
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;

}//main

void
A3BHbsCreatorMover::apply(
	core::pose::Pose & pose
)
{
	using namespace core::scoring::constraints;

	// create score function
	scoring::ScoreFunctionOP score_fxn( scoring::get_score_function() );
	add_fa_constraints_from_cmdline_to_scorefxn(*score_fxn);
	add_fa_constraints_from_cmdline_to_pose(pose);
	
	if( score_fxn->has_zero_weight( atom_pair_constraint ) )
		score_fxn->set_weight( atom_pair_constraint, 0.1 );
	
	if( score_fxn->has_zero_weight( dihedral_constraint ) )
		score_fxn->set_weight( dihedral_constraint, 0.1 );
	
	if( score_fxn->has_zero_weight( angle_constraint ) )
		score_fxn->set_weight( angle_constraint, 0.1 );
	
	chemical::ResidueTypeSetCOP restype_set = chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
	
	char hbs_chain = option[a3b_hbs_creator::hbs_chain].value()[0];
	core::Size final_res = option[a3b_hbs_creator::hbs_final_res].value();
	core::pose::PDBInfoCOP pdb_info( pose.pdb_info() );

	// shift the jump
	// account for the position of the offset so we never delete the jump res
	Size offset = Size( numeric::random::uniform() * 4 );
	kinematics::FoldTree f = pose.fold_tree();
	f.slide_jump( 1, 1, pose.pdb_info()->pdb2pose( hbs_chain, final_res+1+offset ) );
	pose.fold_tree( f );
	
	// DELETE EXTRANEOUS
	for ( Size i = 1; i <= pose.total_residue(); ++i) {
		char chn = pdb_info->chain(i);
		//TR << "evaluating residue " << chn  << " " << pdb_info->number(i) << std::endl;
		
		if ( chn != hbs_chain) {
			continue;
		}
		// correct chain to be truncated and prepped
		
		core::Size pdb_res_num = pdb_info->number(i);
		
		// hbs pre is the smallest number of what we want to preserve
		if ( pdb_res_num < final_res ) {
			//TR << "deleting residue " << pdb_res_num  << " which was " << core::chemical::oneletter_code_from_aa(pose.aa(i)) << std::endl;
			while ( pdb_res_num < final_res ) {
				pose.delete_polymer_residue(i);
				pdb_res_num = pdb_info->number(i);
			}
			
		} else if ( pdb_res_num > final_res + option[a3b_hbs_creator::hbs_length].value() ) {
			//TR << "deleting residue " << pdb_res_num << std::endl;
			while ( chn == hbs_chain && i <= pose.total_residue() ) {
				chn = pdb_info->chain(i);
				pose.delete_polymer_residue(i);
			}
		}
	}
	TR << "making an a3b peptide" << std::endl;
	id::AtomID fixed_atom( pose.residue(1).atom_index( "CA" ), 1 );
	
	// CONVERT TO A3B
	pose::Pose a3bpose;
	utility::vector1< Vector > former_CAs;
	pose::PDBInfoCOP old_pdb_info = pose.pdb_info();
	for ( Size i = 1; i <= pose.total_residue(); ++i) {

		char chn = pdb_info->chain(i);
		//TR << "evaluating residue " << chn  << " " << pdb_info->number(i) << std::endl;
		
		if ( chn != hbs_chain) {
			continue;
		}
		// replace each HBS residue with an a3b
		TR << "Recording CA for " << i << std::endl;
		former_CAs.push_back( pose.residue( i ).atom( "CA" ).xyz() );
		std::string name;
		if ( i % 4 == offset ) {
			name = alpha_to_beta( pose.residue(i).name() );
		} else {
			name = pose.residue(i).name();
		}
		
		Residue r = *new Residue( restype_set->name_map( name ), true );
		r.set_all_chi(pose.residue(i).chi());
		
		if ( a3bpose.total_residue() == 0 ) {
			a3bpose.append_residue_by_jump( r , 1 );
		} else {
			a3bpose.append_residue_by_bond( r, true );
		}

		//a3bpose.pdb_info()->append_res( )
		// port residue info
		//a3bpose.pdb_info()->set_resinfo(
		//		a3bpose.total_residue(),
		//		hbs_chain,
		//		pose.pdb_info()->number(i) );
	}
	
	// NOW set dihedrals because you've added all the residues
	// otherwise you'll fail at setting psi EVERY TIME you doofus
	for ( Size i = 1; i <= a3bpose.total_residue(); ++i ) {
		TR << "Setting dihedrals for " << i << std::endl;
		if ( a3bpose.residue(i).type().is_beta_aa() ) {
			core::id::TorsionID phi( i, id::BB, 1);
			core::id::TorsionID tht( i, id::BB, 2);
			core::id::TorsionID psi( i, id::BB, 3);
			core::id::TorsionID omg( i, id::BB, 4);
			
			a3bpose.conformation().set_torsion( phi, -105);
			a3bpose.conformation().set_torsion( tht,   65);
			a3bpose.conformation().set_torsion( psi, -115);
			a3bpose.conformation().set_torsion( omg,  180);
		} else {
			
			core::id::TorsionID phi( i, id::BB, 1);
			core::id::TorsionID psi( i, id::BB, 2);
			core::id::TorsionID omg( i, id::BB, 3);
			
			a3bpose.conformation().set_torsion( phi,  -60);
			a3bpose.conformation().set_torsion( psi,  -47);
			a3bpose.conformation().set_torsion( omg,  180);
			
		}
	}
	
	for ( Size i = 1; i <= a3bpose.total_residue() - 4; ++i ) {
		TR << "constraining " << i << std::endl;
		
		a3bpose.add_constraint(
				ConstraintOP(
						new AtomPairConstraint(
								*new id::AtomID( a3bpose.residue( i   ).atom_index( "O" ), i ),
								*new id::AtomID( a3bpose.residue( i+4 ).atom_index( "H" ), i+4 ),
								core::scoring::func::HarmonicFuncOP( new core::scoring::func::HarmonicFunc( 1.9, 0.04 ) )
		) ) );
	}
	
	// presently the final residue in the pose is the terminal residue of the hbs
	// replace with terminal variant
	conformation::Residue term( restype_set->get_residue_type_with_variant_added(a3bpose.residue(a3bpose.total_residue()).type(), chemical::METHYLATED_CTERMINUS_VARIANT), true );
	term.set_all_chi(a3bpose.residue(a3bpose.total_residue()).chi());
	//replace_res_post.mainchain_torsions(pose.residue(oop_post_pos_).mainchain_torsions());
	a3bpose.replace_residue( a3bpose.total_residue(), term, true );
	conformation::idealize_position( a3bpose.total_residue(), a3bpose.conformation() );
	a3bpose.add_constraint(
			ConstraintOP(
						new AtomPairConstraint(
								*new id::AtomID( a3bpose.residue( a3bpose.total_residue()-4   ).atom_index( "O" ), a3bpose.total_residue()-4 ),
								*new id::AtomID( a3bpose.residue( a3bpose.total_residue() ).atom_index( "HM" ), a3bpose.total_residue() ),
								core::scoring::func::HarmonicFuncOP( new core::scoring::func::HarmonicFunc( 1.9, 0.04 ) )
	) ) );
	
	// PATCH
	// If offset is two then we need a normal HBS patch for this pattern
	
	if ( offset == 2 ) {
		hbs::HbsPatcherOP hbs_patcher( new hbs::HbsPatcher( 1 ) );
		hbs_patcher->apply( a3bpose );
	} else {
		a3b_hbs::A3BHbsPatcherOP hbs_patcher( new a3b_hbs::A3BHbsPatcher( 1 ) );
		hbs_patcher->apply( a3bpose );
	}
	
	//pose.set_phi(i+3, -40);
	//pose.set_psi(i+3, -58);
	
	//pose.dump_pdb( "postdihedrals.pdb");
	a3bpose.conformation().declare_chemical_bond( 1, "CYH", 3, "CZH" );

	a3bpose.conformation().detect_bonds();
	a3bpose.conformation().detect_pseudobonds();
	for ( core::Size i = 1; i <= a3bpose.total_residue(); ++i ) {
		a3bpose.conformation().update_polymeric_connection(i);
	}
	
	a3bpose.dump_pdb( "postpseudobonds.pdb");
	


	
	// minimize new pose a bit
	kinematics::MoveMapOP a3blittlemm( new kinematics::MoveMap() );
	a3blittlemm->set_bb( true );
	a3blittlemm->set_chi( true );
	simple_moves::MinMoverOP a3blittlemin( new protocols::simple_moves::MinMover( a3blittlemm, score_fxn, "lbfgs_armijo_nonmonotone", 1, true ) );
	a3blittlemin->cartesian( true );
	a3blittlemin->apply( a3bpose );
	a3bpose.dump_pdb( "posta3blittlemin.pdb");

	// delete res
	for ( Size i = 1; i <= pose.total_residue(); ++i) {
		char chn = pdb_info->chain(i);
		TR << "evaluating residue " << chn  << " " << pdb_info->number(i) << std::endl;
		
		if ( chn != hbs_chain) {
			continue;
		}
		
		pose.delete_residue_range_slow( i, pose.total_residue() );
		break;
	}
	
	Size last_of_first_chain = pose.total_residue();
	pose::append_pose_to_pose( pose, a3bpose, true );
	for ( Size i = last_of_first_chain + 1; i <= pose.total_residue(); ++i) {
		TR << "Setting chain for " << i << " to " << hbs_chain << std::endl;
		pose.pdb_info()->chain( i, hbs_chain );
		pose.pdb_info()->number( i, pdb_info->number( i ) );
		
		pose.add_constraint(
				ConstraintOP(
						new CoordinateConstraint(
								*new id::AtomID( pose.residue( i ).atom_index( "CA" ), i ),
								fixed_atom,
								former_CAs[ i - last_of_first_chain ],
								core::scoring::func::HarmonicFuncOP( new core::scoring::func::HarmonicFunc( 0, 0 ) )
		) ) );
		
		
	}
	pose.dump_pdb( "postcombination.pdb");

	setup_pert_foldtree(pose);
	
	// TINYMIN
	// create move map for minimization
	kinematics::MoveMapOP littlemm( new kinematics::MoveMap() );
	littlemm->set_bb( false );
	for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
		if ( pdb_info->chain(i) != hbs_chain ) {
			continue;
		}
		
		littlemm->set_bb( true );
	}
	//score_fxn->set_weight( atom_pair_constraint, 50 );
	littlemm->set_chi( true );
	littlemm->set_jump( 1, true );
	simple_moves::MinMoverOP littlemin( new protocols::simple_moves::MinMover( littlemm, score_fxn, option[ OptionKeys::run::min_type ].value(), 1, true ) );
	littlemin->cartesian( true );
	littlemin->apply( pose );
	pose.dump_pdb( "postlittlemin.pdb");

	// SCAN DIHEDRALS
	//scan( pose, score_fxn, hbs_chain );
	
	// TINYMIN
	TR << "Minimizing away any major clashes (score start: " << ( *score_fxn )( pose );
	littlemin->apply( pose );
	TR << "; score end: " << ( *score_fxn )( pose );
	repack( pose, score_fxn );
	TR << "; after repack: " << ( *score_fxn )( pose ) << ")" << std::endl;
	
	if( option[ a3b_hbs_creator::final_mc ].value() ) {
		moves::SequenceMoverOP pert_sequence( new moves::SequenceMover() );
		moves::MonteCarloOP pert_mc( new moves::MonteCarlo( pose, *score_fxn, 0.2 ) );
		
		kinematics::MoveMapOP pert_pep_mm( new kinematics::MoveMap() );
		
		// core::Size hbs_position = 1;
		
		for( Size i = 1; i <= pose.total_residue(); ++i ) {
			if ( pdb_info->chain(i) == hbs_chain ) { //i != ( ( unsigned ) ( option[ hbs_creator::hbs_final_res ].value() ) ) ) {
				//if ( pose.residue_type( i ).is_l_aa() ) {
				TR << "setting small movable resid: "<< i<<std::endl;
				//kdrew: commenting out because small mover fails randomly
				pert_pep_mm->set_bb( i );
				//}
			}
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
		
	}
	pose.dump_pdb( "postmc.pdb");
	
	if( option[ a3b_hbs_creator::final_repack ].value() ) {
		repack( pose, score_fxn );
	}
	pose.dump_pdb( "postrepack.pdb");
	
	if ( option[ a3b_hbs_creator::final_minimize ].value() ) {
		using namespace core::id;
		using namespace core::scoring;
		using namespace core::scoring::constraints;
		
		// create move map for minimization
		kinematics::MoveMapOP mm( new kinematics::MoveMap() );
		mm->set_bb( true );
		mm->set_chi( true );
		//mm->set_jump( 1, true );
		
		score_fxn->set_weight( atom_pair_constraint, 0.5 );
		score_fxn->set_weight( dihedral_constraint, 0.5 );
		score_fxn->set_weight( angle_constraint, 0.5 );
		
		// create minimization mover
		simple_moves::MinMoverOP minM( new protocols::simple_moves::MinMover( mm, score_fxn, "lbfgs_armijo_nonmonotone", 0.01,	true ) );
		minM->cartesian( true );
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
