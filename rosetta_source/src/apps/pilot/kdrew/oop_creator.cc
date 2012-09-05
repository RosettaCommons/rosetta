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
#include <core/chemical/ResidueTypeSet.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/constraints/util.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

// Mover headers
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/PyMolMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/oop/OopPuckMover.hh>
#include <protocols/simple_moves/oop/OopPatcher.hh>

// Utility Headers
#include <devel/init.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
//#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
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
using namespace protocols::simple_moves::oop;
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
	IntegerOptionKey const prepend_n_residues( "oop_creator::prepend_n_residues" );
	IntegerOptionKey const append_n_residues( "oop_creator::append_n_residues" );

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
	utility::vector1< core::Size > empty_vector(0);
	option.add( oop_creator::oop_plus_positions, "The positions of the first residues of plus oop rings" ).def( empty_vector );
	option.add( oop_creator::oop_minus_positions, "The positions of the first residues of minus oop rings" ).def( empty_vector );
	option.add( oop_creator::prepend_n_residues, "Number of residues to prepend" ).def( 0 );
	option.add( oop_creator::append_n_residues, "Number of residues to append" ).def( 0 );

	// init command line options
	//you MUST HAVE THIS CALL near the top of your main function, or your code will crash when you first access the command line options
	devel::init(argc, argv);
	//basic::options::option[ basic::options::OptionKeys::chemical::include_patches](utility::tools::make_vector1( std::string("patches/oop_pre.txt"), std::string("patches/oop_post.txt") ) );

	//create mover instance
	OopCreatorMoverOP OC_mover( new OopCreatorMover() );

	//call job distributor
	protocols::jd2::JobDistributor::get_instance()->go( OC_mover );

}//main

void
OopCreatorMover::apply(
	core::pose::Pose & pose
)
{

	// create score function
	//kdrew: old standard scoring function, using MM scoring function now because of NCAAs
	//scoring::ScoreFunctionOP score_fxn( ScoreFunctionFactory::create_score_function( scoring::STANDARD_WTS, scoring::SCORE12_PATCH ) );
	scoring::ScoreFunctionOP score_fxn( ScoreFunctionFactory::create_score_function( scoring::MM_STD_WTS) );
	scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn(*score_fxn);

	scoring::constraints::add_fa_constraints_from_cmdline_to_pose(pose);

	//kdrew: positions, plus, minus, random, patch
	//oop::OopPuckMoverOP opm_plus( new oop::OopPuckMover( option[ oop_creator::oop_plus_positions].value() , true, false, false, true) );
	//opm_plus->apply(pose);
	//oop::OopPuckMoverOP opm_minus( new oop::OopPuckMover( option[ oop_creator::oop_minus_positions].value(), false, true, false, true ) );
	//opm_minus->apply(pose);

	utility::vector1< core::Size > const plus_positions = option[ oop_creator::oop_plus_positions].value();
	for(Size i = 1; i <= plus_positions.size(); i++)
	{                                     
		oop::OopPatcherOP oop_patcher (new oop::OopPatcher( plus_positions[i] ) );
		oop_patcher->apply( pose );

		oop::OopPuckPlusMoverOP opm_plus( new oop::OopPuckPlusMover( plus_positions[i] ) );
		opm_plus->apply( pose );
	}

	utility::vector1< core::Size > const minus_positions = option[ oop_creator::oop_minus_positions].value();
	for(Size i = 1; i <= minus_positions.size(); i++)
	{
		oop::OopPatcherOP oop_patcher (new oop::OopPatcher( minus_positions[i] ) );
		oop_patcher->apply( pose );

		oop::OopPuckMinusMoverOP opm_minus( new oop::OopPuckMinusMover( minus_positions[i] ) );
		opm_minus->apply( pose );
	}


	// create a task factory and task operations
	TaskFactoryOP tf(new TaskFactory());
	tf->push_back( new core::pack::task::operation::InitializeFromCommandline );

	operation::ReadResfileOP rrop( new operation::ReadResfile() );
	rrop->default_filename();
	tf->push_back( rrop );
 
	// create a pack rotamers mover
	simple_moves::PackRotamersMoverOP packer( new protocols::simple_moves::PackRotamersMover() );
	packer->task_factory( tf );
	packer->score_function( score_fxn );

	packer->apply(pose);


	//kdrew: create glycine residue
	ResidueTypeSet const & rsd_set( pose.residue(1).residue_type_set() );
	ResidueOP gly( ResidueFactory::create_residue( rsd_set.name_map( "GLY" ) ) );

	Size pep_begin( pose.conformation().chain_begin( 1 ) );
	Size pep_end( pose.conformation().chain_end( 1 ) );

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

	// create move map for minimization
	kinematics::MoveMapOP mm( new kinematics::MoveMap() );
	mm->set_bb( true );
	mm->set_chi( true );
	mm->set_jump( 1, true );

	// create minimization mover
	simple_moves::MinMoverOP minM( new protocols::simple_moves::MinMover( mm, score_fxn, option[ OptionKeys::run::min_type ].value(), 0.01,	true ) );

//kdrew: only turn on pymol observer in debug mode
#ifndef NDEBUG
	protocols::moves::PyMolObserverOP pymover = protocols::moves::AddPyMolObserver(pose);
#endif

	//kdrew: minimizer not working after appending/prepending residues, not sure why
	// final min (okay to use ta min here)
	//minM->apply( pose );

}


