// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


//Headers are generally organized by either what they do or where they come from.  This organization is first core library headers, then protocols library, then utility stuff.


// Project Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueConnection.hh>
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
#include <protocols/simple_moves/oop/OopPuckMover.hh>
#include <protocols/simple_moves/oop/OopRandomSmallMover.hh>
#include <protocols/simple_moves/oop/OopPatcher.hh>
#include <protocols/simple_moves/chiral/ChiralMover.hh>

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
#include <utility/excn/Exceptions.hh>

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
using namespace protocols::simple_moves::oop;
using namespace protocols::simple_moves::chiral;
using namespace core::pack::task;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::id;
using basic::T;
using basic::Error;
using basic::Warning;
using utility::file::FileName;

//kdrew: this app adds oop patches to the given pdb strucure

// tracer - used to replace cout
static basic::Tracer TR("OOP_Creator");

// application specific options
namespace oop_creator{
	// pert options
	IntegerVectorOptionKey const oop_plus_positions( "oop_creator::oop_plus_positions" );
	IntegerVectorOptionKey const oop_minus_positions( "oop_creator::oop_minus_positions" );
	IntegerVectorOptionKey const oop_d_plus_positions( "oop_creator::oop_d_plus_positions" );
	IntegerVectorOptionKey const oop_d_minus_positions( "oop_creator::oop_d_minus_positions" );
	IntegerVectorOptionKey const oop_low_e_puck_positions( "oop_creator::oop_low_e_puck_positions" );
	IntegerOptionKey const prepend_n_residues( "oop_creator::prepend_n_residues" );
	IntegerOptionKey const append_n_residues( "oop_creator::append_n_residues" );
	BooleanOptionKey const final_repack( "oop_creator::final_repack" );
	BooleanOptionKey const final_minimize( "oop_creator::final_minimize" );
	BooleanOptionKey const final_mc ( "oop_creator::final_mc" );
	BooleanOptionKey const correct_oop_post ( "oop_creator::correct_oop_post" );

}

class OopCreatorMover : public Mover {

	public:

		//default ctor
		OopCreatorMover(): Mover("OopCreatorMover"){}

		//default dtor
		virtual ~OopCreatorMover(){}

		virtual void apply( core::pose::Pose & pose );
		virtual std::string get_name() const { return "OopCreatorMover"; }

};

typedef utility::pointer::owning_ptr< OopCreatorMover > OopCreatorMoverOP;
typedef utility::pointer::owning_ptr< OopCreatorMover const > OopCreatorMoverCOP;


int
main( int argc, char* argv[] )
{
    try {
	utility::vector1< core::Size > empty_vector(0);
	option.add( oop_creator::oop_plus_positions, "The positions of the first residues of plus oop rings" ).def( empty_vector );
	option.add( oop_creator::oop_minus_positions, "The positions of the first residues of minus oop rings" ).def( empty_vector );
	option.add( oop_creator::oop_d_plus_positions, "The positions of the first residues of chiral d plus oop rings" ).def( empty_vector );
	option.add( oop_creator::oop_d_minus_positions, "The positions of the first residues of chiral d minus oop rings" ).def( empty_vector );
	option.add( oop_creator::oop_low_e_puck_positions, "The positions of the first oop residues, pucker will change to low e conformation" ).def( empty_vector );
	option.add( oop_creator::prepend_n_residues, "Number of residues to prepend" ).def( 0 );
	option.add( oop_creator::append_n_residues, "Number of residues to append" ).def( 0 );
	option.add( oop_creator::final_repack, "Do a final repack. Default false" ).def(false);
	option.add( oop_creator::final_minimize, "Do a final minimization. Default false" ).def(false);
	option.add( oop_creator::final_mc, "Do a final monte carlo on oop. Default false" ).def(false);
	option.add( oop_creator::correct_oop_post, "Correct oop post phi/psi to low energy well. Default false" ).def(false);

	// init command line options
	//you MUST HAVE THIS CALL near the top of your main function, or your code will crash when you first access the command line options
	devel::init(argc, argv);
	//basic::options::option[ basic::options::OptionKeys::chemical::include_patches](utility::tools::make_vector1( std::string("patches/oop_pre.txt"), std::string("patches/oop_post.txt") ) );

	//create mover instance
	OopCreatorMoverOP OC_mover( new OopCreatorMover() );

	//call job distributor
	protocols::jd2::JobDistributor::get_instance()->go( OC_mover );

    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
    }
    return 0;
}//main

void
OopCreatorMover::apply(
	core::pose::Pose & pose
)
{

	// create score function
	//kdrew: old standard scoring function, using MM scoring function now because of NCAAs
	//scoring::ScoreFunctionOP score_fxn( getScoreFunction() );
	scoring::ScoreFunctionOP score_fxn( ScoreFunctionFactory::create_score_function( scoring::MM_STD_WTS) );
	scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn(*score_fxn);

	scoring::constraints::add_fa_constraints_from_cmdline_to_pose(pose);

	//kdrew: positions, plus, minus, random, patch
	//oop::OopPuckMoverOP opm_plus( new oop::OopPuckMover( option[ oop_creator::oop_plus_positions].value() , true, false, false, true) );
	//opm_plus->apply(pose);
	//oop::OopPuckMoverOP opm_minus( new oop::OopPuckMover( option[ oop_creator::oop_minus_positions].value(), false, true, false, true ) );
	//opm_minus->apply(pose);

	utility::vector1< core::Size > all_positions;

	utility::vector1< core::Size > const plus_positions = option[ oop_creator::oop_plus_positions].value();
	all_positions.insert( all_positions.end(), plus_positions.begin(), plus_positions.end() );
	for(Size i = 1; i <= plus_positions.size(); i++)
	{
		oop::OopPatcherOP oop_patcher (new oop::OopPatcher( plus_positions[i] ) );
		oop_patcher->apply( pose );

		oop::OopPuckPlusMoverOP opm_plus( new oop::OopPuckPlusMover( plus_positions[i] ) );
		opm_plus->apply( pose );
	}

	utility::vector1< core::Size > const minus_positions = option[ oop_creator::oop_minus_positions].value();
	all_positions.insert( all_positions.end(), minus_positions.begin(), minus_positions.end() );
	for(Size i = 1; i <= minus_positions.size(); i++)
	{
		oop::OopPatcherOP oop_patcher (new oop::OopPatcher( minus_positions[i] ) );
		oop_patcher->apply( pose );

		oop::OopPuckMinusMoverOP opm_minus( new oop::OopPuckMinusMover( minus_positions[i] ) );
		opm_minus->apply( pose );
	}

	utility::vector1< core::Size > const d_plus_positions = option[ oop_creator::oop_d_plus_positions].value();
	all_positions.insert( all_positions.end(), d_plus_positions.begin(), d_plus_positions.end() );
	for(Size i = 1; i <= d_plus_positions.size(); i++)
	{
		oop::OopPatcherOP oop_patcher (new oop::OopPatcher( d_plus_positions[i] ) );
		oop_patcher->apply( pose );

		oop::OopDPuckPlusMoverOP opm_plus( new oop::OopDPuckPlusMover( d_plus_positions[i] ) );
		opm_plus->apply( pose );
	}

	utility::vector1< core::Size > const d_minus_positions = option[ oop_creator::oop_d_minus_positions].value();
	all_positions.insert( all_positions.end(), d_minus_positions.begin(), d_minus_positions.end() );
	for(Size i = 1; i <= d_minus_positions.size(); i++)
	{
		oop::OopPatcherOP oop_patcher (new oop::OopPatcher( d_minus_positions[i] ) );
		oop_patcher->apply( pose );

		oop::OopDPuckMinusMoverOP opm_minus( new oop::OopDPuckMinusMover( d_minus_positions[i] ) );
		opm_minus->apply( pose );
	}

	utility::vector1< core::Size > const low_e_puck_positions = option[ oop_creator::oop_low_e_puck_positions ].value();
	all_positions.insert( all_positions.end(), low_e_puck_positions.begin(), low_e_puck_positions.end() );
	for(Size i = 1; i <= low_e_puck_positions.size(); ++i)
	{
		oop::OopPatcherOP oop_patcher (new oop::OopPatcher( low_e_puck_positions[i] ) );
		oop_patcher->apply( pose );

		//kdrew: this is currently not implemented, do some mc random puck and small phi/psi moves to find lowest energy pucker, TODO?
		//oop::OopLowEMoverOP opm( new oop::OopLowEMover( low_e_puck_positions[i] ) );
		//opm->apply( pose );

		//kdrew: poor man's version for now, if L use PuckPlus, if D use DPuckPlus
		if( is_l_chiral( pose.residue_type( low_e_puck_positions[i] ) ) )
		{
        	//kdrew: use PuckPlus
			oop::OopPuckPlusMoverOP opm_plus( new oop::OopPuckPlusMover( low_e_puck_positions[i] ) );
			opm_plus->apply( pose );
		}
		else if( is_d_chiral( pose.residue_type( low_e_puck_positions[i] ) ) )
		{
        	//kdrew: use DPuckPlus
			oop::OopDPuckPlusMoverOP opm_plus( new oop::OopDPuckPlusMover( low_e_puck_positions[i] ) );
			opm_plus->apply( pose );

		}
		else
		{
			TR << " residue: " << pose.residue_type( low_e_puck_positions[i] ).name() << " not found in chiral map" <<  std::endl;
			TR << " possibly achiral (ex GLY) or not listed in map" <<  std::endl;
			TR << " not changing" << std::endl;
		}

	}

	//kdrew: sets oop_post phi/psi near low energy well
	if( option[ oop_creator::correct_oop_post ].value() )
	{
		for( Size i = 1; i <= all_positions.size(); ++i )
		{
			//kdrew: the +1 is to get the oop_post position
			if( is_d_chiral( pose.residue_type( all_positions[i] +1 ) ) )
			{
				pose.set_phi( all_positions[i] +1, 135.0 ) ;
				pose.set_psi( all_positions[i] +1, -70.0 ) ;
			}
			//kdrew: defaults to L conformation
			else
			{
				pose.set_phi( all_positions[i] +1, -135.0 ) ;
				pose.set_psi( all_positions[i] +1, 70.0 ) ;
			}
		}

	}


	//kdrew: create glycine residue
	ResidueTypeSet const & rsd_set( pose.residue(1).residue_type_set() );
	ResidueOP gly( ResidueFactory::create_residue( rsd_set.name_map( "GLY" ) ) );

	Size pep_begin( pose.conformation().chain_begin( 1 ) );
	Size pep_end( pose.conformation().chain_end( 1 ) );

	//kdrew: since we probably added new connection types (i.e. oop CYP and CZP atoms) above, need to reset connections
	pose.conformation().detect_bonds();
	pose.conformation().detect_pseudobonds();
	for(core::Size i=1; i<=pose.total_residue(); ++i){
		pose.conformation().update_polymeric_connection(i);
	}


	//kdrew: grabbed code from chrisk pep_spec
	//kdrew: append residues , hard coded to glycine
	for ( Size i = 1; i <= Size( option[ oop_creator::append_n_residues ].value() ) ; ++i ) {
		TR << "in append: " << pep_end << std::endl;
		pose.conformation().safely_append_polymer_residue_after_seqpos( *gly, pep_end, true );
		pep_end = pep_end + 1;
		pose.set_omega( pep_end - 1, 180.0 );
		pose.conformation().update_polymeric_connection( pep_end );
		pose.conformation().update_polymeric_connection( pep_end - 1 );
	}
	//kdrew: prepend residues , hard coded to glycine
	for ( Size i = 1; i <= Size( option[ oop_creator::prepend_n_residues ].value() ) ; ++i ) {
		TR << "in prepend: " << pep_begin << std::endl;
		pose.conformation().safely_prepend_polymer_residue_before_seqpos( *gly, pep_begin, true );
		pep_end = pep_end + 1;
		pep_begin =  pose.conformation().chain_begin( 1 ) ; //reset pep beginning
		pose.set_omega( pep_begin, 180.0 );
		pose.conformation().update_polymeric_connection( pep_begin );
		pose.conformation().update_polymeric_connection( pep_begin + 1 );
	}


	//kdrew: monte carlo phi/psi of oop to find low energy
	if( option[ oop_creator::final_mc ].value() )
	{
		moves::SequenceMoverOP pert_sequence( new moves::SequenceMover() );
		moves::MonteCarloOP pert_mc( new moves::MonteCarlo( pose, *score_fxn, 0.2 ) );

		kinematics::MoveMapOP pert_pep_mm( new kinematics::MoveMap() );
		simple_moves::SmallMoverOP pert_pep_small( new simple_moves::SmallMover( pert_pep_mm, 0.2, 1 ) );
		pert_pep_small->angle_max( 'H', 2.0 );
		pert_pep_small->angle_max( 'L', 2.0 );
		pert_pep_small->angle_max( 'E', 2.0 );

		utility::vector1< core::Size > oop_pre_positions;
		//kdrew: load all oop_pre positions into vector and make all non-oop_pre positions movable by small mover
		for( Size i = 1; i <= pose.total_residue(); ++i )
		{
			TR << "resid: " << i << " is OOP_PRE: " << pose.residue(i).has_variant_type(chemical::OOP_PRE) << std::endl;
			if( pose.residue(i).has_variant_type(chemical::OOP_PRE) != 1 )
			{
				if( is_l_chiral( pose.residue_type( i ) ) )
				{
					TR << "setting small movable resid: "<< i<<std::endl;
					//kdrew: commenting out because small mover fails randomly
					//pert_pep_mm->set_bb( i );
				}
			}
			else
			{ oop_pre_positions.push_back(i); }
		}

		pert_sequence->add_mover( pert_pep_small );

		//kdrew: add all oop_pre positions to random small mover
		if( oop_pre_positions.size() > 0 )
		{
			oop::OopRandomSmallMoverOP opm( new oop::OopRandomSmallMover ( oop_pre_positions, 2.0 ) );
			moves::RepeatMoverOP pert_pep_repeat( new moves::RepeatMover( opm, oop_pre_positions.size() * 1000 ) );
			pert_sequence->add_mover( pert_pep_repeat );
		}

		moves::TrialMoverOP pert_trial( new moves::TrialMover( pert_sequence, pert_mc ) );

		pert_trial->apply( pose );
    	pert_mc->recover_low( pose );

	}

	if( option[ oop_creator::final_repack ].value() )
	{

		// create a task factory and task operations
		TaskFactoryOP tf(new TaskFactory());
		tf->push_back( new core::pack::task::operation::InitializeFromCommandline );

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

	//kdrew: monte carlo phi/psi of oop to find low energy
	if( option[ oop_creator::final_mc ].value() )
	{
		moves::SequenceMoverOP pert_sequence( new moves::SequenceMover() );
		moves::MonteCarloOP pert_mc( new moves::MonteCarlo( pose, *score_fxn, 0.2 ) );

		kinematics::MoveMapOP pert_pep_mm( new kinematics::MoveMap() );
		simple_moves::SmallMoverOP pert_pep_small( new simple_moves::SmallMover( pert_pep_mm, 0.2, 1 ) );
		pert_pep_small->angle_max( 'H', 2.0 );
		pert_pep_small->angle_max( 'L', 2.0 );
		pert_pep_small->angle_max( 'E', 2.0 );

		utility::vector1< core::Size > oop_pre_positions;
		//kdrew: load all oop_pre positions into vector and make all non-oop_pre positions movable by small mover
		for( Size i = 1; i <= pose.total_residue(); ++i )
		{
			if( pose.residue(i).has_variant_type(chemical::OOP_PRE) != 1)
			{ pert_pep_mm->set_bb( i ); }
			else
			{ oop_pre_positions.push_back(i); }
		}

		pert_sequence->add_mover( pert_pep_small );

		//kdrew: add all oop_pre positions to random small mover
		if( oop_pre_positions.size() > 0 )
		{
			oop::OopRandomSmallMoverOP opm( new oop::OopRandomSmallMover ( oop_pre_positions, 2.0 ) );
			moves::RepeatMoverOP pert_pep_repeat( new moves::RepeatMover( opm, oop_pre_positions.size() * 1000 ) );
			pert_sequence->add_mover( pert_pep_repeat );
		}

		moves::TrialMoverOP pert_trial( new moves::TrialMover( pert_sequence, pert_mc ) );

		pert_trial->apply( pose );
    	pert_mc->recover_low( pose );

	}

	if( option[ oop_creator::final_repack ].value() )
	{

		// create a task factory and task operations
		TaskFactoryOP tf(new TaskFactory());
		tf->push_back( new core::pack::task::operation::InitializeFromCommandline );

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


	if( option[ oop_creator::final_minimize ].value() )
	{
		using namespace core::id;
		using namespace core::scoring;
		using namespace core::scoring::constraints;

		//kdrew: add constraints to omega angle, (this problem might have been fixed and these constraints are unnecessary)
		for( Size i = 1; i < pose.conformation().chain_end( 1 ); ++i )
		{
			id::AtomID id1,id2,id3,id4;
			core::id::TorsionID torsion_id = TorsionID( i, id::BB, 3 ); //kdrew: 3 is omega angle

			//kdrew: put constraint on omega angle
			pose.conformation().get_torsion_angle_atom_ids( torsion_id, id1, id2, id3, id4 );

			Real torsion_value( pose.torsion( torsion_id ) );

			CircularHarmonicFuncOP circularharm_func  (new CircularHarmonicFunc( numeric::conversions::radians( torsion_value ), numeric::conversions::radians( 10.0 ) ) );

			ConstraintCOP dihedral1 = new DihedralConstraint( id1, id2, id3, id4, circularharm_func );

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


