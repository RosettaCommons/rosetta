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
#include <core/chemical/AtomType.hh>
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
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <core/scoring/func/HarmonicFunc.hh>
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
#include <protocols/simple_moves/RandomTorsionMover.hh>
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
namespace macrocycles {
// pert options
IntegerOptionKey const macrocycle_size ( "macrocycles::macrocycle_size" );
StringOptionKey const macrocycle_sequence ( "macrocycles::macrocycle_sequence" );
BooleanOptionKey const final_repack( "macrocycles::final_repack" );
BooleanOptionKey const final_minimize( "macrocycles::final_minimize" );
BooleanOptionKey const final_mc ( "macrocycles::final_mc" );
IntegerOptionKey const max_closure_cycles( "macrocycles::max_closure_cycles" );
}



bool
bump_check( Pose const & pose, Real const multiplier  ) {

	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		//TR << "Checking residue " << ii << std::endl;
		for ( Size ai = 1; ai <= pose.residue( ii ).natoms(); ++ai ) {
			//TR << "Checking atom " << pose.residue(ii).atom_name(ai) << std::endl;

			core::Real const ai_rad = pose.residue( ii ).type().atom_type( ai ).lj_radius();

			// intra
			for ( Size ai2 = 1; ai2 <= pose.residue( ii ).natoms(); ++ai2 ) {
				//TR << "vs atom " << pose.residue(ii).atom_name(ai2) << std::endl;

				if ( ai == ai2 ) continue;
				if ( pose.residue(ai).type().atoms_are_bonded(ai, ai2) ) continue;
				core::Real const ai2_rad = pose.residue( ii ).type().atom_type( ai2 ).lj_radius();

				if ( pose.residue( ii ).xyz( ai ).distance_squared(
						pose.residue( ii ).xyz( ai2 ) )  < pow((ai_rad+ai2_rad)*multiplier, 2) ) {
					TR << "Residue " << ii << " atom " << pose.residue(ii).atom_name(ai) << "  " << ii << " atom " << pose.residue(ai2).atom_name(ai2) << pose.residue( ii ).xyz( ai ).distance_squared( pose.residue( ii ).xyz( ai2 ) ) <<std::endl;

					return false;
				}
			}

			// inter
			for ( Size jj = 1; jj <= pose.total_residue(); ++jj ) {
				//TR << "vs residue " << jj << std::endl;

				if ( ii == jj ) continue;

				for ( Size aj = 1; aj <= pose.residue( jj ).natoms(); ++aj ) {
					//TR << "vs atom " << pose.residue(jj).atom_name(aj) << std::endl;

					// are ai and aj bonded?
					//bool are_bonded = false;
					if ( pose.residue(ii).is_bonded( pose.residue( jj ) ) ) {
						continue;
						// Get the list of connection ids in residue 1 that connect to residue 2:
						/*utility::vector1< core::Size > connlist = pose.residue(ii).connections_to_residue(jj);

						for ( core::Size c = 1; c<=connlist.size(); ++c ) {
						if ( pose.residue(ii).residue_connect_atom_index(connlist[c]) != ai ) continue;
						if ( pose.residue(jj).residue_connect_atom_index( pose.residue(ii).residue_connection_conn_id( connlist[c]) ) == aj) are_bonded = true;
						}*/
					}
					//if ( are_bonded ) continue;

					core::Real const aj_rad = pose.residue( jj ).type().atom_type( aj ).lj_radius();

					if ( pose.residue( ii ).xyz( ai ).distance_squared(
							pose.residue(jj).xyz( aj ) )  < pow((ai_rad+aj_rad)*multiplier, 2) ) {
						TR << "Residue " << ii << " atom " << pose.residue(ii).atom_name(ai) << "  " << jj << " atom " << pose.residue(jj).atom_name(aj) << pose.residue( jj ).xyz( aj ).distance_squared( pose.residue( jj ).xyz( aj ) ) << std::endl;
						return false;
					}
				}
			}
		}
	}
	return true;
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
		option.add( macrocycles::max_closure_cycles, "Maximum number of iterations through the closure loop.  Default 25." ).def(25);

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

	core::chemical::ResidueTypeSetCOP residue_set_cap = ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

	Size nres;

	if ( option[ macrocycles::macrocycle_sequence ].user() ) {
		// make_pose_from_sequence( pose, option[ macrocycles::macrocycle_sequence ].value(), *residue_set_cap, false );
		// nres = pose.total_residue();

		// Avoid bad behavior from make_pose_from_sequence for now...

		core::chemical::ResidueTypeCOPs requested_types = residue_types_from_sequence( option[ macrocycles::macrocycle_sequence ].value(), *residue_set_cap, false );
		nres = requested_types.size();
		Residue first( *requested_types[1], true );
		pose.conformation().append_residue_by_jump( first, 1 );
		for ( Size i = 2; i <= nres; ++i ) {
			Residue ala( *requested_types[i], true );
			pose.conformation().append_residue_by_bond( ala, true );
		}
	} else {
		nres = option[ macrocycles::macrocycle_size ].value();
		ResidueType const & ala_type = residue_set_cap->name_map( "ALA" );

		Residue ala( ala_type, true );
		pose.conformation().append_residue_by_jump( ala, 1 );
		for ( Size i = 2; i <= nres; ++i ) {
			pose.conformation().append_residue_by_bond( ala, true );
		}
	}
	// Removed - no longer using CyclizationMover
	//Size const chain = 1;
	//Size const minimization_rounds = 15;

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

	std::string next_to_C = pose.residue( nres ).type().is_beta_aa() ? "CM" : "CA";
	std::string next_to_N = pose.residue(  1   ).type().is_beta_aa() ? "CA" : "CA";

	pose.add_constraint(
		DihedralConstraintOP( new DihedralConstraint(
		*new AtomID( pose.residue( nres ).atom_index(next_to_C), nres ),
		*new AtomID( pose.residue( nres ).atom_index("C" ), nres ),
		*new AtomID( pose.residue( 1    ).atom_index("N" ), 1    ),
		*new AtomID( pose.residue( 1    ).atom_index(next_to_N), 1    ),
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

	pose.add_constraint(
		DihedralConstraintOP( new DihedralConstraint(
		*new AtomID( pose.residue( nres ).atom_index(next_to_C), nres ),
		*new AtomID( pose.residue( nres ).atom_index("C" ), nres ),
		*new AtomID( pose.residue( 1    ).atom_index("N" ), 1    ),
		*new AtomID( pose.residue( nres ).atom_index("O" ), nres ),
		CircularHarmonicFuncOP( new CircularHarmonicFunc(
		numeric::NumericTraits<float>::pi(), 0.02 ) ) ) ) );

	pose.add_constraint(
		DihedralConstraintOP( new DihedralConstraint(
		*new AtomID( pose.residue( 1    ).atom_index(next_to_N ), 1 ),
		*new AtomID( pose.residue( 1    ).atom_index("H" ), 1 ),
		*new AtomID( pose.residue( 1    ).atom_index("N" ), 1    ),
		*new AtomID( pose.residue( nres ).atom_index("C" ), nres ),
		CircularHarmonicFuncOP( new CircularHarmonicFunc(
		numeric::NumericTraits<float>::pi(), 0.02 ) ) ) ) );

	pose.add_constraint(
		AngleConstraintOP( new AngleConstraint(
		*new AtomID( pose.residue( 1    ).atom_index("H" ), 1 ),
		*new AtomID( pose.residue( 1    ).atom_index("N" ), 1 ),
		*new AtomID( pose.residue( nres ).atom_index("C" ), nres ),
		CircularHarmonicFuncOP( new CircularHarmonicFunc(
		numeric::NumericTraits<float>::pi()*2.0/3.0, 0.02 ) ) ) ) );

	pose.add_constraint(
		AngleConstraintOP( new AngleConstraint(
		*new AtomID( pose.residue( nres ).atom_index("O" ), nres ),
		*new AtomID( pose.residue( nres ).atom_index("C" ), nres ),
		*new AtomID( pose.residue( 1    ).atom_index("N" ), 1    ),
		CircularHarmonicFuncOP( new CircularHarmonicFunc(
		numeric::NumericTraits<float>::pi()*2.0/3.0, 0.02 ) ) ) ) );

	for ( Size ii = 1; ii <= nres; ++ii ) {
		if ( pose.residue(ii).type().is_beta_aa() ) {
			if ( ii > 1 ) {
				pose.set_torsion( core::id::TorsionID( ii, id::BB, 1 ), numeric::random::rg().uniform() * 360.0 );
			}
			if ( ii < nres ) {
				pose.set_torsion( core::id::TorsionID( ii, id::BB, 2 ), numeric::random::rg().uniform() * 360.0 );
				pose.set_torsion( core::id::TorsionID( ii, id::BB, 3 ), numeric::random::rg().uniform() * 360.0 );
				pose.set_torsion( core::id::TorsionID( ii, id::BB, 4 ), 180.0 );
			}
		} else {
			if ( ii > 1 ) {
				pose.set_phi( ii, numeric::random::rg().uniform() * 360.0 );
			}
			if ( ii < nres ) {
				pose.set_psi( ii, numeric::random::rg().uniform() * 360.0 );
				pose.set_omega( ii, 180.0 );
			}
		}
	}

	// create move map for minimization
	kinematics::MoveMapOP mm( new kinematics::MoveMap() );
	mm->set_bb( true );
	//mm->set_bb( 1, false );
	//mm->set_bb( nres, false );
	mm->set_chi( true );
	mm->set_jump( 1, true );

	// create minimization mover
	simple_moves::MinMoverOP minM( new protocols::simple_moves::MinMover( mm, score_fxn, option[ OptionKeys::run::min_type ].value(), 0.01, true ) );

	ScoreFunctionOP cart_sfxn = score_fxn->clone();
	cart_sfxn->set_weight( cart_bonded, 1.0 );
	cart_sfxn->set_weight( pro_close, 0.0 );

	simple_moves::MinMoverOP cart_min( new protocols::simple_moves::MinMover( mm, cart_sfxn, "lbfgs_armijo_nonmonotone", 0.01, true ) );
	cart_min->cartesian( true );

	AtomPairConstraintOP topout( new AtomPairConstraint(
		*new AtomID( pose.residue( nres ).atom_index("C" ), nres ),
		*new AtomID( pose.residue( 1    ).atom_index("N" ), 1    ),
		//TopOutFuncOP( new TopOutFunc( 100, 1.33, 5 ) ) ) );
		HarmonicFuncOP( new HarmonicFunc( 1.33, 0.03 ) ) ) );

	pose.add_constraint( topout );
	pose.conformation().declare_chemical_bond( nres, "C", 1, "N" );


	Pose oldpose = pose;

	moves::SequenceMoverOP pert_sequence( new moves::SequenceMover() );
	moves::MonteCarloOP pert_mc( new moves::MonteCarlo( pose, *score_fxn, 0.2 ) );

	kinematics::MoveMapOP pert_pep_smallshear_mm( new kinematics::MoveMap() );
	kinematics::MoveMapOP pert_pep_beta_mm( new kinematics::MoveMap() );

	bool is_all_beta = true;
	bool is_all_alpha = true;
	for ( Size i = 1; i <= pose.total_residue(); ++i ) {
		if ( pose.residue( i ).type().is_alpha_aa() ) {
			is_all_beta = false;
			pert_pep_smallshear_mm->set_bb( i );
		} else {
			is_all_alpha = false;
			pert_pep_beta_mm->set_bb( i );
		}
	}

	simple_moves::SmallMoverOP pert_pep_small( new simple_moves::SmallMover( pert_pep_smallshear_mm, 0.2, 1 ) );
	pert_pep_small->angle_max( 'H', 2.0 );
	pert_pep_small->angle_max( 'L', 2.0 );
	pert_pep_small->angle_max( 'E', 2.0 );
	simple_moves::ShearMoverOP pert_pep_shear( new simple_moves::ShearMover( pert_pep_smallshear_mm, 0.2, 1 ) );
	pert_pep_shear->angle_max( 'H', 2.0 );
	pert_pep_shear->angle_max( 'L', 2.0 );
	pert_pep_shear->angle_max( 'E', 2.0 );
	simple_moves::RandomTorsionMoverOP pert_beta( new simple_moves::RandomTorsionMover( pert_pep_beta_mm, 0.2, 1 ) );

	moves::RandomMoverOP pert_pep_random( new moves::RandomMover() );
	if ( !is_all_beta ) {
		pert_pep_random->add_mover( pert_pep_small, 1 );
		pert_pep_random->add_mover( pert_pep_shear, 1 );
	}
	if ( !is_all_alpha ) {
		pert_pep_random->add_mover( pert_beta, 2 );
	}
	moves::RepeatMoverOP pert_pep_repeat( new moves::RepeatMover( pert_pep_random, 10 ) );

	pert_sequence->add_mover( pert_pep_repeat );
	moves::TrialMoverOP pert_trial( new moves::TrialMover( pert_sequence, pert_mc ) );

	core::Size iteration_counter(0);
	core::Size const max_iterations( static_cast<core::Size>( option[ macrocycles::max_closure_cycles ] ) );
	while ( !closed ) {
		++iteration_counter;
		if ( iteration_counter > max_iterations ) {
			TR << "Max iterations (" << max_iterations << ") exceeded.  Exiting from closure cycle." << std::endl;
			break;
		}

		pose = oldpose;
		/********************************************************\
		Initial MC
		\********************************************************/
		pert_mc->reset( pose );
		for ( Size j = 1; j <= 300; ++j ) {
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

		for ( Real new_apc_wt = 0.1; new_apc_wt <= 1.0; new_apc_wt += 0.1 ) {
			minM->apply( pose );
			score_fxn->set_weight( atom_pair_constraint, new_apc_wt );
		}
		TR << "After minimization, score is " << ( *score_fxn )( pose ) << std::endl;

		cart_min->apply( pose );
		TR << "After cartesian minimization, score is " << ( *score_fxn )( pose ) << std::endl;
		minM->apply( pose );
		TR << "After minimization, score is " << ( *score_fxn )( pose ) << std::endl;

		//pose.remove_constraint( topout );


		/********************************************************\
		Cyclization
		\********************************************************/
		//CyclizationMoverOP cm( new CyclizationMover( chain, true, true, minimization_rounds ) );
		//cm->apply( pose );
		//pose.update_residue_neighbors();

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
			} else {
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

			if ( score_fxn->has_zero_weight( dihedral_constraint ) ) {
				score_fxn->set_weight( dihedral_constraint, 1.0 );
			}

			if ( score_fxn->has_zero_weight( atom_pair_constraint ) ) {
				score_fxn->set_weight( atom_pair_constraint, 1.0 );
			}

			minM->apply( pose );
		}

		Real omg = numeric::dihedral_degrees(
			pose.residue( nres).xyz(next_to_C),
			pose.residue( nres).xyz("C"),
			pose.residue( 1 ).xyz("N"),
			pose.residue( 1 ).xyz(next_to_N)
		);

		Real imp1 = numeric::dihedral_degrees(
			pose.residue( nres).xyz(next_to_C),
			pose.residue( nres).xyz("C"),
			pose.residue( 1 ).xyz("N"),
			pose.residue( nres ).xyz("O")
		);

		Real imp2 = numeric::dihedral_degrees(
			pose.residue( 1).xyz(next_to_N),
			pose.residue( 1).xyz("H"),
			pose.residue( 1 ).xyz("N"),
			pose.residue( nres ).xyz("C")
		);

		Real ang1 = numeric::angle_degrees(
			pose.residue( 1).xyz("H"),
			pose.residue( 1 ).xyz("N"),
			pose.residue( nres ).xyz("C")
		);

		Real ang2 = numeric::angle_degrees(
			pose.residue( 1).xyz("O"),
			pose.residue( 1 ).xyz("C"),
			pose.residue( nres ).xyz("N")
		);

		bool bump_good = true;//bump_check( pose, 0.3 );

		bool dihedral_good = ( omg > 170 && omg <= 180 ) || ( omg >= -180 && omg < 170 );
		dihedral_good &= ( imp1 > 170 && imp1 <= 180 ) || ( imp1 >= -180 && imp1 < 170 );
		dihedral_good &= ( imp2 > 170 && imp2 <= 180 ) || ( imp2 >= -180 && imp2 < 170 );
		bool dist_good = ( pose.residue( nres ).xyz("C").distance( pose.residue( 1 ).xyz( "N" ) ) < 1.4 );

		bool ang_good = ( ang1 > 110 && ang1 < 130 );
		ang_good &= ( ang2 > 110 && ang2 < 130 );

		TR << "Omega: " << omg  << std::endl;
		TR << "Imp1: " << imp1 << std::endl;
		TR << "Imp2: " << imp2 << std::endl;
		TR << "Ang1: " << ang1 << std::endl;
		TR << "Ang2: " << ang2 << std::endl;
		TR << "Dist: " << pose.residue( nres ).xyz("C").distance( pose.residue( 1 ).xyz( "N" ) ) << std::endl;
		TR << "Bump: " << bump_good << std::endl;

		if ( dihedral_good && dist_good && bump_good && ang_good ) {
			closed = true;
		} else {
			TR << "Closure failed; restarting" << std::endl;
		}
	}
}
