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
#include <core/pose/PDBInfo.hh>
#include <core/pose/PDBPoseMap.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueTypeSet.hh>
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
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/chiral/ChiralMover.hh>
#include <protocols/simple_moves/oop/OopMover.hh>

//Basic headers
#include <basic/resource_manager/ResourceManager.hh>

#include <numeric/conversions.hh>

// Utility Headers
#include <devel/init.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
//#include <basic/options/keys/OptionKeys.hh>
#include <utility/options/keys/StringVectorOptionKey.fwd.hh>
#include <utility/options/keys/StringVectorOptionKey.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/excn/Exceptions.hh>

// C++ headers
#include <string>
#include <iostream>
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
using namespace protocols::simple_moves::chiral;
using namespace protocols::simple_moves::oop;
using namespace core::pack::task;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::id;
using basic::T;
using basic::Error;
using basic::Warning;
using utility::file::FileName;

//kdrew: this app alters chiral configuration of the given pdb strucure

// tracer - used to replace cout
static thread_local basic::Tracer TR( "LD_CONVERTER" );

// application specific options
namespace ld_converter{
	// pert options
	IntegerVectorOptionKey const convert_L_positions( "ld_converter::convert_L_positions" );
	IntegerVectorOptionKey const convert_D_positions( "ld_converter::convert_D_positions" );
	StringVectorOptionKey const convert_L_pdb_positions( "ld_converter::convert_L_pdb_positions" );
	StringVectorOptionKey const convert_D_pdb_positions( "ld_converter::convert_D_pdb_positions" );
	BooleanOptionKey const final_repack( "ld_converter::final_repack" );
	BooleanOptionKey const final_minimize( "ld_converter::final_minimize" );
	BooleanOptionKey const fix_functional_group( "ld_converter::fix_functional_group" );

}

class LDConverterMover : public Mover {

	public:

		//default ctor
		LDConverterMover(): Mover("LDConverterMover"){}

		//default dtor
		virtual ~LDConverterMover(){}

		virtual void apply( core::pose::Pose & pose );
		virtual std::string get_name() const { return "LDConverterMover"; }

};

typedef utility::pointer::owning_ptr< LDConverterMover > LDConverterMoverOP;
typedef utility::pointer::owning_ptr< LDConverterMover const > LDConverterMoverCOP;


int
main( int argc, char* argv[] )
{
    try {
	utility::vector1< core::Size > empty_vector(0);
	utility::vector1< std::string > empty_vector_string(0);
	option.add( ld_converter::convert_L_positions, "The residue positions to convert to L configuration" ).def( empty_vector );
	option.add( ld_converter::convert_D_positions, "The residue positions to convert to D configuration" ).def( empty_vector );

	option.add( ld_converter::convert_D_pdb_positions, "The residue positions to convert to D configuration in pdb numbering (chain resnum)" ).def( empty_vector_string );
	option.add( ld_converter::convert_L_pdb_positions, "The residue positions to convert to L configuration in pdb numbering (chain resnum)" ).def( empty_vector_string );

	option.add( ld_converter::final_repack, "Do a final repack. Default false" ).def( false );
	option.add( ld_converter::final_minimize, "Do a final minimization. Default false" ).def( false );
	option.add( ld_converter::fix_functional_group, "Keep the side chain fixed rather than the backbone. Default false" ).def( false );

	// init command line options
	//you MUST HAVE THIS CALL near the top of your main function, or your code will crash when you first access the command line options
	devel::init(argc, argv);

	//basic::options::option[ basic::options::OptionKeys::chemical::include_patches](utility::tools::make_vector1( std::string("patches/oop_pre.txt"), std::string("patches/oop_post.txt") ) );

	//create mover instance
	LDConverterMoverOP LDC_mover( new LDConverterMover() );

	//call job distributor
	protocols::jd2::JobDistributor::get_instance()->go( LDC_mover );

    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
				return -1;
    }
    return 0;
}//main

void
LDConverterMover::apply(
	core::pose::Pose & pose
)
{
	// create score function
	//kdrew: using MM scoring function because of NCAAs
	scoring::ScoreFunctionOP score_fxn( ScoreFunctionFactory::create_score_function( scoring::MM_STD_WTS) );
	scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn(*score_fxn);

	scoring::constraints::add_fa_constraints_from_cmdline_to_pose(pose);

	utility::vector1< core::Size > l_positions = option[ ld_converter::convert_L_positions ].value();
	//kdrew: parse pdb numbering in convert_pdb_position options
	utility::vector1< std::string > const l_pdb_positions = option[ ld_converter::convert_L_pdb_positions].value();
	core::pose::PDBPoseMap const & pose_map(pose.pdb_info()->pdb2pose());
	//kdrew: chain and res numbers are a pair so increment by 2
	for( Size i = 1; i <= l_pdb_positions.size(); i=i+2 )
	{
		core::Size respos = pose_map.find( *l_pdb_positions[i].c_str(), atoi(l_pdb_positions[i+1].c_str()) );
        l_positions.push_back(respos);
	}

	for(Size i = 1; i <= l_positions.size(); ++i )
	{
		TR << "residue id: " << l_positions[i] << std::endl;
		protocols::simple_moves::chiral::ChiralMoverOP cm( new protocols::simple_moves::chiral::ChiralMover( l_positions[i], chiral::L_CHIRALITY, option[ ld_converter::fix_functional_group].value() ) );
		cm->apply( pose );

		//kdrew: if an oop residue, update hydrogens
		if ( pose.residue( l_positions[i] ).has_variant_type( core::chemical::OOP_PRE ) == 1 )
		{
			oop::OopMoverOP opm( new oop::OopMover( l_positions[i] ) );
			opm->update_hydrogens( pose );
		}
	}

	utility::vector1< core::Size > d_positions = option[ ld_converter::convert_D_positions ].value();
	utility::vector1< std::string > const d_pdb_positions = option[ ld_converter::convert_D_pdb_positions].value();
	//kdrew: parse pdb numbering in convert_pdb_position options
	//kdrew: chain and res numbers are a pair so increment by 2
	for( Size i = 1; i <= d_pdb_positions.size(); i=i+2 )
	{
		core::Size respos = pose_map.find( *d_pdb_positions[i].c_str(), std::atoi(d_pdb_positions[i+1].c_str()) );
        d_positions.push_back(respos);

	}
	for(Size i = 1; i <= d_positions.size(); ++i )
	{
		TR << "residue id: " << d_positions[i] << std::endl;
		protocols::simple_moves::chiral::ChiralMoverOP cm( new protocols::simple_moves::chiral::ChiralMover( d_positions[i], chiral::D_CHIRALITY , option[ ld_converter::fix_functional_group].value()) );
		cm->apply( pose );

		//kdrew: if an oop residue, update hydrogens
		if ( pose.residue( d_positions[i] ).has_variant_type( core::chemical::OOP_PRE ) == 1 )
		{
			TR << "updating oop hydrogens" << std::endl;
			oop::OopMoverOP opm( new oop::OopMover( d_positions[i] ) );
			opm->update_hydrogens( pose );
		}
	}

	//kdrew: if no positions were passed in, make all positions L configuration
	if ( l_positions.size() == 0 && d_positions.size() == 0 )
	{
		TR << "making all positions L configuration" << std::endl;
		for( Size i = 1; i <= pose.conformation().chain_end( 1 ); ++i )
		{
			TR << "residue id: " << i << std::endl;
			protocols::simple_moves::chiral::ChiralMoverOP cm( new protocols::simple_moves::chiral::ChiralMover( i, option[ ld_converter::fix_functional_group].value() ) );
			cm->apply( pose );

			//kdrew: if an oop residue, update hydrogens
			if ( pose.residue( i ).has_variant_type( core::chemical::OOP_PRE ) == 1 )
			{
				oop::OopMoverOP opm( new oop::OopMover( i ) );
				opm->update_hydrogens( pose );
			}
		}

	}


	if( option[ ld_converter::final_repack ].value() )
	{

		// create a task factory and task operations
		using core::pack::task::operation::TaskOperationCOP;
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

	if( option[ ld_converter::final_minimize].value() )
	{
		using namespace core::id;
		using namespace core::scoring;
		using namespace core::scoring::constraints;

		for( Size i = 1; i < pose.conformation().chain_end( 1 ); ++i )
		{
			id::AtomID id1,id2,id3,id4;
			core::id::TorsionID torsion_id = TorsionID( i, id::BB, 3 );

			//kdrew: put constraint on omega angle
			pose.conformation().get_torsion_angle_atom_ids( torsion_id, id1, id2, id3, id4 );

			Real torsion_value( pose.torsion( torsion_id ) );

			core::scoring::func::CircularHarmonicFuncOP circularharm_func  (new core::scoring::func::CircularHarmonicFunc( numeric::conversions::radians( torsion_value ), numeric::conversions::radians( 3.0 ) ) );

			ConstraintCOP dihedral1 = new DihedralConstraint( id1, id2, id3, id4, circularharm_func );

			pose.add_constraint( dihedral1 );
		}
		if( score_fxn->has_zero_weight( dihedral_constraint ) )
		{
        	score_fxn->set_weight( dihedral_constraint, 1.0 );
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
	//	protocols::moves::PyMolObserverOP pymover = protocols::moves::AddPyMolObserver(pose);
	//#endif

		//kdrew: minimizer not working after appending/prepending residues, not sure why
		// final min (okay to use ta min here)
		minM->apply( pose );
	}

}
