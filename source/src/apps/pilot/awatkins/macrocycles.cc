// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Project Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/Patch.hh>
#include <core/chemical/VariantType.hh>

#include <numeric/xyz.functions.hh>

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
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <core/scoring/func/TopOutFunc.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

// Mover headers
#include <protocols/ncbb/util.hh>
#include <protocols/loops/loop_closure/kinematic_closure/KinematicMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/PyMolMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/CyclizationMover.hh>
#include <protocols/simple_moves/chiral/ChiralMover.hh>
#include <protocols/rigid/RB_geometry.hh>

#include <numeric/conversions.hh>
#include <numeric/random/random.hh>

//Basic headers
#include <basic/resource_manager/ResourceManager.hh>

// Utility Headers
#include <devel/init.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
//#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
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
using namespace core::chemical;
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
static basic::Tracer TR("MacrocycleCreator");

// application specific options
namespace macrocycles{
	// pert options
	IntegerOptionKey const macrocycle_size ( "macrocycles::macrocycle_size" );
	StringOptionKey const macrocycle_sequence ( "macrocycles::macrocycle_sequence" );
	BooleanOptionKey const final_repack( "macrocycles::final_repack" );
	BooleanOptionKey const final_minimize( "macrocycles::final_minimize" );
	BooleanOptionKey const final_mc ( "macrocycles::final_mc" );
}

class MacrocycleMover : public Mover {

	public:

		//default ctor
		MacrocycleMover(): Mover("MacrocycleMover"){}

		//default dtor
		virtual ~MacrocycleMover(){}

		virtual void apply( core::pose::Pose & pose );
		virtual std::string get_name() const { return "MacrocycleMover"; }

};

typedef utility::pointer::shared_ptr< MacrocycleMover > MacrocycleMoverOP;
typedef utility::pointer::shared_ptr< MacrocycleMover const > MacrocycleMoverCOP;


int
main( int argc, char* argv[] )
{
try {
	utility::vector1< core::Size > empty_vector(0);

	option.add( macrocycles::macrocycle_size, "Size of macrocycle. Default 6." ).def(6);
	option.add( macrocycles::macrocycle_sequence, "Sequence of macrocycle. Default ""." ).def("");
	option.add( macrocycles::final_repack, "Do a final repack. Default false" ).def(false);
	option.add( macrocycles::final_minimize, "Do a final minimization. Default false" ).def(false);
	option.add( macrocycles::final_mc, "Do a final monte carlo on macrocycle. Default false" ).def(false);

	// init command line options
	//you MUST HAVE THIS CALL near the top of your main function, or your code will crash when you first access the command line options
	devel::init(argc, argv);
	option[ OptionKeys::chemical::patch_selectors ].push_back( "NTERM_CONNECT" );
	option[ OptionKeys::chemical::patch_selectors ].push_back( "CTERM_CONNECT" );
	
	//create mover instance
	MacrocycleMoverOP Mac_mover( new MacrocycleMover() );

	//call job distributor
	protocols::jd2::JobDistributor::get_instance()->go( Mac_mover );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;

}//main

void
MacrocycleMover::apply(
	core::pose::Pose & pose
)
{
	scoring::ScoreFunctionOP score_fxn = scoring::get_score_function();
	scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn(*score_fxn);
	scoring::constraints::add_fa_constraints_from_cmdline_to_pose(pose);

	score_fxn->set_weight( core::scoring::atom_pair_constraint, 1.0 );
	score_fxn->set_weight( core::scoring::angle_constraint, 1.0 );
	score_fxn->set_weight( core::scoring::dihedral_constraint, 1.0 );//10.0 );

	pose.clear();
	
	bool closed = false;
	
	while ( !closed ) {
	core::chemical::ResidueTypeSetCOP residue_set_cap = ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
	
	Size nres;

	if ( option[ macrocycles::macrocycle_sequence ].user() ) {
		make_pose_from_sequence( pose, option[ macrocycles::macrocycle_sequence ].value(), *residue_set_cap, false );
		nres = pose.total_residue();
	} else {
		nres = option[ macrocycles::macrocycle_size ].value();
		ResidueType const & ala_type = residue_set_cap->name_map( "ALA" );
	
		Residue ala( ala_type, true );
		pose.conformation().append_residue_by_jump( ala, 1 );
		for ( Size i = 2; i <= nres; ++i ) {
			pose.conformation().append_residue_by_bond( ala, true );
		}
	}
	
	Size const chain = 1;
	Size const minimization_rounds = 15;
	
	/*
	// Add omega constraints
	for ( Size i = 1; i <= nres; ++i ) {
		
		using namespace core::id;
		using namespace core::scoring;
		using namespace core::scoring::func;
		using namespace core::scoring::constraints;
		
		AtomID aidO  ( pose.residue( i   ).atom_index("O"), i   );
		AtomID aidCA1( pose.residue( i   ).atom_index("CA"), i   );
		AtomID aidC  ( pose.residue( i   ).atom_index("C"),  i   );
		AtomID aidN  ( pose.residue( i%nres+1 ).atom_index("N"),  i%nres+1 );
		AtomID aidCA2( pose.residue( i%nres+1 ).atom_index("CA"), i%nres+1 );
		AtomID aidH( pose.residue( i%nres+1 ).atom_index("H"),  i%nres+1 );
		
		//pose.add_constraint( DihedralConstraintOP( new DihedralConstraint( aidCA1, aidC, aidN, aidCA2, CircularHarmonicFuncOP( new CircularHarmonicFunc( numeric::NumericTraits<float>::pi(), 0.02 ) ) ) ) );
		pose.add_constraint( DihedralConstraintOP( new DihedralConstraint( aidO, aidC, aidN, aidH, CircularHarmonicFuncOP( new CircularHarmonicFunc( numeric::NumericTraits<float>::pi(), 0.02 ) ) ) ) );
		
	}*/
	
	using namespace core::id;
	using namespace core::scoring;
	using namespace core::scoring::func;
	using namespace core::scoring::constraints;
	pose.add_constraint(
			DihedralConstraintOP( new DihedralConstraint(
					*new AtomID( pose.residue( nres ).atom_index("CA"), nres ),
					*new AtomID( pose.residue( nres ).atom_index("C" ), nres ),
					*new AtomID( pose.residue( 1    ).atom_index("N" ), 1    ),
					*new AtomID( pose.residue( 1    ).atom_index("CA"), 1    ),
					CircularHarmonicFuncOP( new CircularHarmonicFunc(
							numeric::NumericTraits<float>::pi(), 0.02 ) ) ) ) );
	
	pose.add_constraint(
			DihedralConstraintOP( new DihedralConstraint(
					*new AtomID( pose.residue( nres ).atom_index("O" ), nres ),
					*new AtomID( pose.residue( nres ).atom_index("C" ), nres ),
					*new AtomID( pose.residue( 1    ).atom_index("N" ), 1    ),
					*new AtomID( pose.residue( 1    ).atom_index("H" ), 1    ),
					CircularHarmonicFuncOP( new CircularHarmonicFunc(
							numeric::NumericTraits<float>::pi(), 0.02 ) ) ) ) );
	
	for ( Size ii = 1; ii <= nres; ++ii ) {
		pose.set_phi( ii, numeric::random::rg().uniform() * 360.0 );
		pose.set_psi( ii, numeric::random::rg().uniform() * 360.0 );
		pose.set_omega( ii, 180.0 );
	}
	
	// create move map for minimization
	kinematics::MoveMapOP mm( new kinematics::MoveMap() );
	mm->set_bb( true );
	mm->set_chi( true );
	mm->set_jump( 1, true );
	
	// create minimization mover
	simple_moves::MinMoverOP minM( new protocols::simple_moves::MinMover( mm, score_fxn, option[ OptionKeys::run::min_type ].value(), 0.01,	true ) );

	AtomPairConstraintOP topout( new AtomPairConstraint(
			*new AtomID( pose.residue( nres ).atom_index("C" ), nres ),
			*new AtomID( pose.residue( 1    ).atom_index("N" ), 1    ),
			TopOutFuncOP( new TopOutFunc( 100, 1.33, 5 ) ) ) );
	
	pose.add_constraint( topout );
	
	/********************************************************\
	 Initial MC
	\********************************************************/
	moves::SequenceMoverOP pert_sequence( new moves::SequenceMover() );
	moves::MonteCarloOP pert_mc( new moves::MonteCarlo( pose, *score_fxn, 0.2 ) );
	
	kinematics::MoveMapOP pert_pep_mm( new kinematics::MoveMap() );
	
	for( Size i = 1; i <= pose.total_residue(); ++i ) {
		pert_pep_mm->set_bb( i );
	}
	
	simple_moves::SmallMoverOP pert_pep_small( new simple_moves::SmallMover( pert_pep_mm, 0.2, 1 ) );
	pert_pep_small->angle_max( 'H', 2.0 );
	pert_pep_small->angle_max( 'L', 2.0 );
	pert_pep_small->angle_max( 'E', 2.0 );
	simple_moves::ShearMoverOP pert_pep_shear( new simple_moves::ShearMover( pert_pep_mm, 0.2, 1 ) );
	pert_pep_shear->angle_max( 'H', 2.0 );
	pert_pep_shear->angle_max( 'L', 2.0 );
	pert_pep_shear->angle_max( 'E', 2.0 );
	
	moves::RandomMoverOP pert_pep_random( new moves::RandomMover() );
	pert_pep_random->add_mover( pert_pep_small, 1 );
	pert_pep_random->add_mover( pert_pep_shear, 1 );
	moves::RepeatMoverOP pert_pep_repeat( new moves::RepeatMover( pert_pep_random, 10 ) );
	
	pert_sequence->add_mover( pert_pep_repeat );
	moves::TrialMoverOP pert_trial( new moves::TrialMover( pert_sequence, pert_mc ) );
	
	for ( Size j = 1; j <= 1000; ++j ) {
		pert_trial->apply( pose );
		
		pert_mc->recover_low( pose );
		//TR<< "pre mc->boltzmann" << std::endl;
		//pert_mc->show_state();
		pert_mc->boltzmann( pose );
		//TR<< "post mc->boltzmann" << std::endl;
		//pert_mc->show_state();
		if ( j % 100 == 0 ) {
			TR << "After round " << j << " of initial MC, score is " << ( *score_fxn )( pose ) << std::endl;
		}
		//minM->apply( pose );
		//TR << "After minimization, score is " << ( *score_fxn )( pose ) << std::endl;
	}
	pert_mc->recover_low( pose );
	minM->apply( pose );
	TR << "After minimization, score is " << ( *score_fxn )( pose ) << std::endl;
	
	pose.remove_constraint( topout );
	
	
	/********************************************************\
	 Cyclization
	\********************************************************/
	CyclizationMoverOP cm( new CyclizationMover( chain, true, true, minimization_rounds ) );
	cm->apply( pose );
	pose.update_residue_neighbors();
	
	if ( option[ macrocycles::final_mc ].value() ) {

		for ( Size j = 1; j <= 10; ++j ) {
		
			pert_trial->apply( pose );
			
			pert_mc->recover_low( pose );
			//TR<< "pre mc->boltzmann" << std::endl;
			//pert_mc->show_state();
			pert_mc->boltzmann( pose );
			//TR<< "post mc->boltzmann" << std::endl;
			//pert_mc->show_state();
			TR << "After round " << j << " of final MC, score is " << ( *score_fxn )( pose ) << std::endl;
		}
	}

	if ( option[ macrocycles::final_repack ].value() ) {

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
			operation::RestrictToRepackingOP rtrp( new operation::RestrictToRepacking() );
			tf->push_back( rtrp );
		}

		// create a pack rotamers mover
		simple_moves::PackRotamersMoverOP packer( new protocols::simple_moves::PackRotamersMover() );
		packer->task_factory( tf );
		packer->score_function( score_fxn );

		packer->apply(pose);
	}


	if ( option[ macrocycles::final_minimize ].value() ) {
		using namespace core::id;
		using namespace core::scoring;
		using namespace core::scoring::constraints;

		//kdrew: if constraint weight is not set on commandline or elsewhere, set to 1.0
		if( score_fxn->has_zero_weight( dihedral_constraint ) )
			score_fxn->set_weight( dihedral_constraint, 1.0 );

		if( score_fxn->has_zero_weight( atom_pair_constraint ) )
			score_fxn->set_weight( atom_pair_constraint, 1.0 );

		minM->apply( pose );
	}
	
	Real omg = numeric::dihedral_degrees(
			pose.residue( nres).xyz("CA"),
			pose.residue( nres).xyz("C"),
			pose.residue( 1 ).xyz("N"),
			pose.residue( 1 ).xyz("CA")
	);
		
		Real imp1 = numeric::dihedral_degrees(
											  pose.residue( nres).xyz("CA"),
											  pose.residue( nres).xyz("C"),
											  pose.residue( 1 ).xyz("N"),
											  pose.residue( nres ).xyz("O")
											  );
		
		Real imp2 = numeric::dihedral_degrees(
											  pose.residue( 1).xyz("CA"),
											  pose.residue( 1).xyz("H"),
											  pose.residue( 1 ).xyz("N"),
											  pose.residue( nres ).xyz("C")
											  );
		
	bool dihedral_good = ( omg > 170 && omg <= 180 ) || ( omg >= -180 && omg < 170 );
		dihedral_good &= ( imp1 > 170 && imp1 <= 180 ) || ( imp1 >= -180 && imp1 < 170 );
		dihedral_good &= ( imp2 > 170 && imp2 <= 180 ) || ( imp2 >= -180 && imp2 < 170 );
	bool dist_good = ( pose.residue( nres ).xyz("C").distance( pose.residue( 1 ).xyz( "N" ) ) < 1.4 );
		
	if ( dihedral_good && dist_good ) {
		closed = true;
	} else {
		TR << "Closure failed; restarting" << std::endl;
	}
	}
}
