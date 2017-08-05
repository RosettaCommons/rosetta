#include <devel/init.hh>

#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/annotated_sequence.hh>  //make_pose_from_sequence
#include <core/pose/util.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/util/SwitchResidueTypeSet.hh>

#include <core/import_pose/pose_stream/util.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>
#include <core/import_pose/import_pose.hh>

#include <protocols/electron_density/util.hh>
#include <protocols/electron_density/SetupForDensityScoringMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/rigid/RB_geometry.hh>

#include <protocols/simple_moves/PackRotamersMover.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>

#include <core/scoring/electron_density/ElectronDensity.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/Energies.hh>

#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <ObjexxFCL/format.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzVector.io.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

#include <basic/Tracer.hh>

#include <iostream>
#include <string>
#include <algorithm>

#include <apps/pilot/rayyrw/util.hh>

#ifndef apps_pilot_rayyrw_min_pack_min_hh
#define apps_pilot_rayyrw_min_pack_min_hh


class MinPackMin;

class MinPackMin{
public:

	MinPackMin();

	void rigid_body_minimization( core::pose::Pose & pose );

	void backbone_minimization( core::pose::Pose & pose );

	void pack_sidechains( core::pose::Pose & pose );

	inline void set_density_wt( core::Real wt ){ density_wt_ = wt; }


private:
	//bool min_bb_;
	//core::Real cycles_;
	core::Real density_wt_;

}; // class declaration

static THREAD_LOCAL basic::Tracer tr("MinPackMin");
static THREAD_LOCAL basic::Tracer rbmin("MinPackMin.rigid_body_minimization");
static THREAD_LOCAL basic::Tracer bbmin("MinPackMin.backbone_minimization");
static THREAD_LOCAL basic::Tracer pack_sc("MinPackMin.pack_sidechains");


MinPackMin::MinPackMin(){
	set_density_wt( 10 );
}


// minimize jump against density
void
MinPackMin::
rigid_body_minimization(
	core::pose::Pose & pose
){
	remove_all_virtual_residues( pose );  // in case adding more than one Virt residue
	core::pose::addVirtualResAsRoot( pose );

	core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction() );

	rbmin << "rigid body minimization" << std::endl;
	rbmin << "use elec_dens_fast to minimize" << std::endl;
	scorefxn->set_weight( core::scoring::elec_dens_fast, density_wt_ );

	// movemap
	core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap );
	movemap->set_bb( false );
	movemap->set_chi( false );

	rbmin << "set up jumps for minimizing pose into density alongs jump(s)";
	// get jump index of root jump
	int root = pose.fold_tree().root();
	rbmin << "root: " << root << std::endl;
	utility::vector1< core::kinematics::Edge > root_edges = pose.fold_tree().get_outgoing_edges(root);
	for ( core::Size i=1; i<=root_edges.size(); ++i ) {
		rbmin << "  " << root_edges[i].label();
		movemap->set_jump ( root_edges[i].label() , true );
	}
	rbmin << std::endl;

	rbmin << "before applying minmover\n" << std::endl;
	scorefxn->show( rbmin, pose );

	protocols::moves::MoverOP min_mover( new protocols::simple_moves::MinMover( movemap, scorefxn, "dfpmin_armijo_nonmonotone", 1e-5, true ) );

	bool densInMinimizer = core::scoring::electron_density::getDensityMap().getUseDensityInMinimizer();
	core::scoring::electron_density::getDensityMap().setUseDensityInMinimizer( true );

	min_mover->apply( pose );
	core::scoring::electron_density::getDensityMap().setUseDensityInMinimizer( densInMinimizer ); // this seems to be unncessary

	rbmin << "after applying minmover\n" << std::endl;
	scorefxn->show( rbmin, pose );
	rbmin << "--------------- done rigid body minimization against density only ---------------" << std::endl;

}


// let everything move!
void
MinPackMin::
backbone_minimization(
	core::pose::Pose & pose
){
	bbmin << "backbone minimization" << std::endl;
	remove_all_virtual_residues( pose );  // in case adding more than one Virt residue
	core::pose::addVirtualResAsRoot( pose );

	if ( ! pose.is_fullatom() ) {
		bbmin << "switch to full atom: " << std::endl;
		core::util::switch_to_residue_type_set( pose, core::chemical::FA_STANDARD );
	}

	core::scoring::ScoreFunctionOP scorefxn =  core::scoring::ScoreFunctionFactory::create_score_function( "soft_rep" );

	scorefxn->set_weight( core::scoring::fa_sol, 0 );
	scorefxn->set_weight( core::scoring::elec_dens_fast, density_wt_ );

	// movemap
	core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap );
	movemap->set_bb( true );
	movemap->set_chi( true );
	movemap->set_jump( true );

	rbmin << "set up jumps for minimizing pose into density alongs jump(s)";
	// get jump index of root jump
	int root = pose.fold_tree().root();
	rbmin << "root: " << root << std::endl;
	utility::vector1< core::kinematics::Edge > root_edges = pose.fold_tree().get_outgoing_edges(root);
	for ( core::Size i=1; i<=root_edges.size(); ++i ) {
		rbmin << "  " << root_edges[i].label();
		movemap->set_jump ( root_edges[i].label() , true );
	}
	bbmin << std::endl;

	bbmin << "before applying minmover\n" << std::endl;
	scorefxn->show( bbmin, pose );

	protocols::moves::MoverOP min_mover( new protocols::simple_moves::MinMover( movemap, scorefxn, "dfpmin_armijo_nonmonotone", 0.01, true ) );

	bool densInMinimizer = core::scoring::electron_density::getDensityMap().getUseDensityInMinimizer();
	core::scoring::electron_density::getDensityMap().setUseDensityInMinimizer( true );

	min_mover->apply( pose );
	core::scoring::electron_density::getDensityMap().setUseDensityInMinimizer( densInMinimizer ); // this seems to be unncessary

	bbmin << "after applying minmover\n" << std::endl;
	scorefxn->show( bbmin, pose );
	bbmin << "--------------- done backbone minimization against density only ---------------" << std::endl;

}


void
MinPackMin::
pack_sidechains(
	core::pose::Pose & pose
){
	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	pack_sc << "------------------------------" << std::endl;
	pack_sc << "packing side chains!" << std::endl;

	remove_all_virtual_residues( pose );  // in case adding more than one Virt residue
	core::pose::addVirtualResAsRoot( pose );

	if ( ! pose.is_fullatom() ) {
		pack_sc << "switch to full atom!!" << std::endl;
		core::util::switch_to_residue_type_set( pose, core::chemical::FA_STANDARD );
	}
	//core::scoring::ScoreFunctionOP scorefxn = core::scoring::ScoreFunctionFactory::create_score_function( core::scoring::STANDARD_WTS, core::scoring::SCORE12_PATCH );
	core::scoring::ScoreFunctionOP scorefxn = core::scoring::ScoreFunctionFactory::create_score_function( "talaris2013" );
	scorefxn->set_weight( core::scoring::fa_sol, 0.0 );
	scorefxn->set_weight( core::scoring::elec_dens_fast, density_wt_ );
	scorefxn->show( pack_sc, pose );

	// Set up packer tasks and move map
	pack_sc << "restrict to packing!!!  no design!" << std::endl;
	TaskFactoryOP tf( new TaskFactory() );
	tf->push_back( TaskOperationCOP( new InitializeFromCommandline() )); // get extra rotamer flags from command line
	tf->push_back( TaskOperationCOP( new operation::IncludeCurrent )); // include current rotamer by default
	tf->push_back( TaskOperationCOP( new RestrictToRepacking() )); // do not design

	// set move map to allow sc minimization
	core::kinematics::MoveMapOP repack_movemap( new core::kinematics::MoveMap() );
	repack_movemap->set_bb( false );
	repack_movemap->set_chi( true );
	repack_movemap->set_jump( false );

	// Make PackRots movers
	protocols::simple_moves::PackRotamersMoverOP packer( new protocols::simple_moves::PackRotamersMover() );
	packer->task_factory( tf );
	packer->score_function( scorefxn );

	packer -> apply(pose);
	pack_sc << "applied packrotamer mover" << std::endl;
	scorefxn->show( pack_sc, pose );

	pack_sc << "--------------- done packing sidechains ---------------" << std::endl;
}


#endif
