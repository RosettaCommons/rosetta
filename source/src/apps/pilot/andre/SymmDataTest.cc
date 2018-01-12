/// @file
/// @brief


// libRosetta headers

#include <protocols/viewer/viewers.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/util.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/conformation/symmetry/SymmData.fwd.hh>
//#include <protocols/symmetric_docking/SymDockProtocol.fwd.hh>
//#include <protocols/rigid/RigidBodyMover.hh>

#include <basic/options/util.hh>//option.hh>
#include <devel/init.hh>

#include <core/io/pdb/pdb_writer.hh>

// packing

//#include <protocols/symmetrical_docking/SymRestrictTaskForDocking.hh>
//#include <protocols/minimization_packing/symmetry/SymPackRotamersMover.hh>


#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <utility/excn/Exceptions.hh>


////////////////////////////////////////////////
// danger USING ////////////////////////////////
using namespace core;
using namespace protocols;
using namespace id;
using namespace pose;
using namespace conformation;
using namespace chemical;
using namespace scoring;
using namespace basic::options;
using namespace optimization;
using utility::vector1;
using std::string;
using core::import_pose::pose_from_file;

static basic::Tracer TR( "apps.pilot.ingemar.symm_test" );


///////////////////////////////////////////////////////////////////////////////
void
SymmDataTest()
{


	using namespace conformation::symmetry;

	Pose pdb_pose;
	core::import_pose::pose_from_file( pdb_pose, start_file() , core::import_pose::PDB_file);

	Pose pose = pdb_pose;
	make_symmetric_pose( pose );

	ScoreFunctionOP scorefxn_sym( ScoreFunctionFactory::create_score_function( core::scoring::PRE_TALARIS_2013_STANDARD_WTS, core::scoring::DOCK_PATCH ) );
	scorefxn_sym->show( std::cout, pose );
	pose.dump_pdb("tmp.pdb");
	Pose pose_asym;
	core::import_pose::pose_from_file( pose_asym, "tmp.pdb" , core::import_pose::PDB_file);
	scorefxn_sym->show( std::cout, pose_asym );
	/*
	SymmetricConformationOP symm_conf (
	dynamic_cast<SymmetricConformation &> ( pose.conformation()) );

	assert( !is_symmetric( pdb_pose ) );
	assert( is_symmetric( *symm_conf ) );
	assert( is_symmetric( pose ) );

	protocols::viewer::add_conformation_viewer( *symm_conf, "symm_conf" );

	ScoreFunctionOP scorefxn_sym( get_score_function_legacy( "score13" ) );

	pack::task::PackerTaskOP packer_task( pack::task::TaskFactory::create_packer_task( pose ));
	utility::vector1<bool> allow_repacked( pose.size(), false );
	for (Size res=1; res <= pose.size(); ++res )
	{
	if ( symm_conf->Symmetry_Info().fa_is_independent(res) ) allow_repacked.at(res) = true;
	}
	packer_task->restrict_to_residues( allow_repacked );
	packer_task->restrict_to_repacking();

	std::map< Size, SymDof > dofs ( symm_conf->Symmetry_Info().get_dofs() );
	pack::symmetric_pack_rotamers( pose, *scorefxn_sym, packer_task);
	//  protocols::rigid::RigidBodySymDofRandomTransMover mover( dofs );

	pose.dump_pdb("before_dofmover.pdb");
	for ( int i = 1; i <= 250; ++i ) {
	//  mover.apply(pose);
	}
	Pose min_pose;
	pose.dump_pdb("after_dofmover.pdb");
	pose_from_file( min_pose, "after_dofmover.pdb" , core::import_pose::PDB_file);
	//core::scoring::ScoreFunctionOP scorefxn_asym ( get_score_function_legacy( "score13" ) );
	std::cout << "score asym, sym " << std::endl;
	//scorefxn_asym->show( std::cout, min_pose );
	scorefxn_sym->show( std::cout, pose );
	kinematics::MoveMap mm;

	// setup moving dofs
	for ( Size i=1; i<= pose.conformation().size(); ++i ) {
	if ( symm_conf->Symmetry_Info().bb_is_independent(i) ) {
	mm.set_bb ( i, false );
	mm.set_chi( i, true );
	} else {
	mm.set_bb ( i, false );
	mm.set_chi( i, false );
	}
	}
	for ( int j=1; j<= 6; ++j ) {
	DOF_ID const & id
	( pose.conformation().dof_id_from_torsion_id(TorsionID(1,JUMP,j)));
	if ( j == 1 || j == 4 || j == 5 || j == 6 ) {
	mm.set(id, true );
	} else {
	mm.set( id, false );
	}
	}
	// mm.set_jump(1,true);
	// setup the options
	// MinimizerOptions options( "lbfgs_armijo_nonmonotone", 10.0, true, true, true );
	//  MinimizerOptions options( "linmin", 10.0, true use_nblist,
	//    true deriv_check, true );

	// optimization::AtomTreeMinimizer minimizer_asym;
	// minimizer_asym.run( min_pose, mm, *scorefxn_asym, options );
	//  optimization::symmetry::SymAtomTreeMinimizer minimizer;

	// minimizer.run( pose, mm, *scorefxn_sym, options );


	//  (*scorefxn_sym)( pose );

	// pose.dump_pdb("after_min.pdb");
	// std::exit(0);
	//  protocols::minimization_packing::PackRotamersMoverOP pack_interface_repack = new minimization_packing::symmetry::SymPackRotamersMover( scfx_cop, packer_task );
	// pack_interface_repack->apply(pose);
	protocols::symmetric_docking::SymDockProtocolOP dock_mover = new protocols::symmetric_docking::SymDockProtocol();
	dock_mover->apply(pose);

	scorefxn_sym->show(std::cout, pose );

	//  protocols::rigid::RigidBodySymDofPerturbMoverOP rb_mover = new protocols::rigid::RigidBodySymDofPerturbMover( 1 , dof, 3, 8 );
	// protocols::symmetry::SymmTranslateMoverOP trnsMoverOP( new protocols::symmetry::SymmTranslateMover() );
	//  rb_mover->apply( pose );

	pose.dump_pdb("test2.pdb");
	*/
}


///////////////////////////////////////////////////////////////////////////////

void*
my_main( void* )
{
	using namespace basic::options;

	SymmDataTest();
	//std::string const mode( option[ OK::dna::specificity::mode ].value() );

	return 0;
}

///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {
		// initialize option and random number system
		devel::init( argc, argv );

		protocols::viewer::viewer_main( my_main );
	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
